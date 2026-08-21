!----------------------------------------------------------------------
! setup_grsett
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setup_grsett.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! Prepares a model for the gravitational-settling solve done by
! grsett.f90: locates the diffusion zone boundaries (excluding any
! convective core/envelope and any hydrogen- or helium-exhausted
! region), converts the model to Bahcall & Loeb (1990) natural units,
! and computes the helium (and, if enabled, iron) settling diffusion
! coefficients diffusion_coeff1/diffusion_coeff2 and their partial
! derivatives with respect to X (diffusion_coeff1_dx/diffusion_
! coeff2_dx) at every zone, using either the simple analytic fit
! (LTHOUL false) or the full Thoul et al. (1994) coefficients (via
! thdiff, LTHOUL true) with a choice of Coulomb-logarithm prescription
! selected by coulomb_log_choice.
subroutine setup_grsett(timestep_seconds, dlnp_dr, log_radius, &
     log_density, mass_grams, log_temperature, convective_flag, &
     num_zones, total_mass, diffusion_coeff1, diffusion_coeff2, &
     composition, radius_bl, temperature_bl, zone_begin, zone_end, &
     fully_convective_flag, diffusion_coeff1_dx, diffusion_coeff2_dx)

      use rotdiff_lib
      use scrtch_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: timestep_seconds
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json), &
           log_temperature(json)
      double precision, intent(inout) :: mass_grams(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      double precision, intent(out) :: diffusion_coeff1(json), &
           diffusion_coeff2(json)
      double precision, intent(in) :: composition(15,json)
      double precision, intent(out) :: radius_bl(json), temperature_bl(json)
      integer, intent(out) :: zone_begin, zone_end
      logical, intent(out) :: fully_convective_flag
      double precision, intent(out) :: diffusion_coeff1_dx(json), &
           diffusion_coeff2_dx(json)










! MHP 3/94 ADDED METAL DIFFUSION
! common/gravsz/: src_grid_metal_diffusion_coeff1/coeff2/coeff1_dz/
! coeff2_dz, all used here (originally COD1Z/COD2Z/QCOD1Z/QCOD2Z, the
! iron-settling analogs of diffusion_coeff1/diffusion_coeff2 on the
! original, non-equally-spaced model grid). Naming matches
! model_to_equal.f90.
      double precision :: src_grid_metal_diffusion_coeff1(json), &
           src_grid_metal_diffusion_coeff2(json), &
           src_grid_metal_diffusion_coeff1_dz(json), &
           src_grid_metal_diffusion_coeff2_dz(json)
      common/gravsz/ src_grid_metal_diffusion_coeff1, &
           src_grid_metal_diffusion_coeff2, &
           src_grid_metal_diffusion_coeff1_dz, &
           src_grid_metal_diffusion_coeff2_dz



! MHP 8/94 ADDED ATOMIC WEIGHTS AND CHARGES FOR H,HE,FE,ELECTRONS -
! NEEDED FOR FULL THOUL COEFFICIENTS
      double precision :: atomic_weight(4), atomic_charge(4), &
           settling_ap(4), settling_at(4), settling_ac(4,4), &
           coulomb_log(4,4), species_mass_fraction(4)
!     additional variables to be declared for variable ln lambda
!     intermediate variables:
! ac_scratch (originally AC) is reused for two unrelated quantities,
! as in the original: the mean charge-weighted atomic weight sum used
! to get the electron density (coulomb_log_choice==4 branch below),
! and, later per zone, the empirical Thoul-fit helium-settling
! amplitude coefficient feeding diffusion_coeff2.
      double precision :: charge_to_mass_sum, ac_scratch, &
           ion_number_density, mean_charge_sq, coulomb_xij, &
           electron_number_density, interion_distance, debye_length, &
           coulomb_lambda, concentration(4)
      double precision :: ln_lambda
      data atomic_weight/1.008d0,4.004d0,55.86d0,5.486d-4/
      data atomic_charge/1.0d0,2.0d0,26.0d0,-1.0d0/
      integer :: num_species
      data num_species/4/
      save

! --- locals ---
      double precision :: solar_radius_bl, seconds_per_year_bl
      integer :: zone_idx
      double precision :: hydrogen_fraction, metal_fraction_total, &
           iron_fraction, hydrogen_fraction_sq, hydrogen_fraction_cubed, &
           hydrogen_metal_product, thoul_denominator
      double precision :: settling_prefactor, iron_settling_ah
      double precision :: settling_coeff_p, settling_coeff_t, dap_dx, dat_dx, dac_dx
      double precision :: rho_local, temp_local
      integer :: species_a_idx, species_b_idx

!     OUTPUT VARIABLES :
!
!     radius_bl - VECTOR OF UNLOGGED RADII IN BAHCALL AND LOEB UNITS.
!     temperature_bl - VECTOR OF UNLOGGED TEMPERATURES IN BAHCALL AND
!          LOEB UNITS.
!     zone_begin - FIRST ZONE FOR DIFFUSION PURPOSES (EITHER THE FIRST
!          MODEL POINT OR THE OUTERMOST POINT OF A CENTRAL CONVECTION
!          ZONE).
!     zone_end - LAST ZONE FOR DIFFUSION PURPOSES (EITHER THE LAST
!          MODEL POINT OR THE INNERMOST POINT OF A SURFACE CONVECTION
!          ZONE).
!     THE VECTORS dlnp_dr AND mass_grams, AND THE SCALARS
!     timestep_seconds AND total_mass, ARE ALSO CONVERTED TO BAHCALL
!     AND LOEB UNITS.
!     CONSTANTS DEFINED :
!     ln10 = CONVERSION FACTOR FROM LN TO LOG10
!     solar_radius_bl = SOLAR RADIUS (CM)
!     seconds_per_year_bl = NUMBER OF SECONDS IN A YEAR.
      solar_radius_bl=6.9598d10
      seconds_per_year_bl=3.1558d7
!     fully_convective_flag=T FOR FULLY CONVECTIVE MODEL(AND IF TRUE,
!     DIFFUSION IS SKIPPED).
      fully_convective_flag=.false.
!     CHECK FOR CONVECTIVE CORE.
      if(convective_flag(1))then
         do 10 zone_idx=2,num_zones
            if(.not.convective_flag(zone_idx))goto 20
   10    continue
!        DIFFUSION NOT COMPUTED FOR FULLY CONVECTIVE MODELS.
         fully_convective_flag=.true.
         write(short_file_unit,15)
   15    format(1x,' FULLY CONVECTIVE MODEL - NO SETTLING')
         goto 9999
   20    continue
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_begin = zone_idx-1
      else
         zone_begin = 1
      endif
! MHP 6/90 CHECK FOR HYDROGEN-EXHAUSTED CORE.
      do 23 zone_idx = zone_begin,num_zones
         if(composition(1,zone_idx).gt.hydrogen_diffusion_floor)goto 25
   23 continue
!     HYDROGEN-FREE MODEL - EXIT.
      write(short_file_unit,16)hydrogen_diffusion_floor
   16 format(1x,'X BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
      fully_convective_flag = .true.
      goto 9999
   25 continue
      zone_begin = zone_idx
!     CHECK FOR CONVECTIVE ENVELOPE.
      if(convective_flag(num_zones))then
         do 30 zone_idx=num_zones-1,2,-1
            if(.not.convective_flag(zone_idx))goto 40
   30    continue
   40    continue
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_end = zone_idx+1
      else
         zone_end = num_zones
      endif
!     CHECK FOR HELIUM-EXHAUSTED SURFACE.
!     OUTER POINT IS SET WHEREVER Y>YMIN.
      do 45 zone_idx=zone_end,1,-1
         if(composition(2,zone_idx).gt.helium_diffusion_min) goto 47
   45 continue
!     HYDROGEN-FREE MODEL - EXIT.
      write(short_file_unit,17)helium_diffusion_min
   17 format(1x,'Y BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
      fully_convective_flag = .true.
      goto 9999
   47 continue
      zone_end = zone_idx
!     rot_diff%bl_mass_scale=CONVERSION FACTOR FOR MASS.
!     CON_RADIUS=CONVERSION FACTOR FOR RADIUS.
!     rot_diff%bl_temp_scale=CONVERSION FACOTR FOR TEMPERATURE.
!     rot_diff%bl_time_scale=CONVERSION FACTOR FOR TIME.
      rot_diff%bl_radius_scale=1.0d0/solar_radius_bl
      rot_diff%bl_mass_scale=1.0d-2*rot_diff%bl_radius_scale**3
      rot_diff%bl_temp_scale=1.0d-7
!     INCLUDES FACTOR OF 2.2 FROM LN LAMBDA
      rot_diff%bl_time_scale=2.7d13*seconds_per_year_bl
!     CONVERT LOG(RADIUS) AND LOG(TEMPERATURE) TO NATURAL UNITS.
!     ALSO CONVERT NATURAL UNITS TO BAHCALL AND LOEB UNITS.
      do 50 zone_idx=1,num_zones

         radius_bl(zone_idx)=exp(ln10*log_radius(zone_idx))*rot_diff%bl_radius_scale
         temperature_bl(zone_idx)=exp(ln10*log_temperature(zone_idx))*rot_diff%bl_temp_scale
         mass_grams(zone_idx)=mass_grams(zone_idx)*rot_diff%bl_mass_scale
         dlnp_dr(zone_idx)=dlnp_dr(zone_idx)/rot_diff%bl_radius_scale
!        SDEL(2,I)=0.4D0   !COMMENT OUT IN REAL CODE
   50 continue
      timestep_seconds=timestep_seconds/rot_diff%bl_time_scale
      total_mass=total_mass*rot_diff%bl_mass_scale
!     SET UP DIFFUSION COEFFICIENTS.
!     MODIFIED BY BC MAY/90 -- VALID FOR ALL X WITH VARIABLE LN LAMBDA
!     GENERAL EQUATION IS
!     DX/DT = (D/DR(D1)+D/DR(D2 DX/DR))/(RHO*R**2), WHERE
!     D1 = R**2/LN LAMBDA  * X  * T**5/2 * (DLNP/DR) * (1-X) *
!          [5/4 + DEL*6*(X-0.32)/(5.4+6.3X-4.5X**2)]
!     D2 = R**2/LN LAMBDA * T**5/2 * (X+3)/(X+1)/(3+5X)
      do 60 zone_idx = 1,num_zones

         hydrogen_fraction = composition(1,zone_idx)
! MHP 10/02 INITIALIZED X - WAS NOT DONE PRIOR TO USAGE IN SHELL 1
         if(coulomb_log_choice.eq.2)then
!           Noerdlinger's formula (1977 A&A 57,407) for LN LAMBDA:
!           Ln Lambda = -19.7 - ln[4/(3*X + 1)]/2 - ln(rho)/2 + 1.5*ln(T)
            ln_lambda=-1.97d1 - 0.5d0*log(1.0d0/(0.75d0*hydrogen_fraction+0.25d0)) &
                      -0.5d0*log_density(zone_idx)*ln10 + 1.5d0*log_temperature(zone_idx)*ln10
         else if(coulomb_log_choice.eq.3)then
!           Loeb's formula (1989 Phys. Rev. D 39, 1009) for LN LAMBDA
!           ln lambda = -19.105747 - ln(rho)/2 + 1.5*ln(T)
            ln_lambda=-1.9105747d1-0.50d0*log_density(zone_idx)*ln10+1.5d0*log_temperature(zone_idx)*ln10
         else
            ln_lambda = 2.2d0
         end if
!
         settling_prefactor=fgry*radius_bl(zone_idx)**2*temperature_bl(zone_idx)**2.5d0/ln_lambda
!         X = HCOMP(1,I)
         metal_fraction_total = composition(3,zone_idx) + composition(4,zone_idx)
         iron_fraction = composition(3,zone_idx)
         hydrogen_fraction_sq = hydrogen_fraction*hydrogen_fraction
         if(.not.lthoul)then
            hydrogen_metal_product = hydrogen_fraction*metal_fraction_total
            hydrogen_fraction_cubed = hydrogen_fraction_sq*hydrogen_fraction
            thoul_denominator=5.4d0+6.3d0*hydrogen_fraction-4.5d0*hydrogen_fraction_sq
            diffusion_coeff1(zone_idx)=settling_prefactor*dlnp_dr(zone_idx)* &
                 (hydrogen_fraction - hydrogen_fraction_sq - hydrogen_metal_product)*(1.25d0+ &
                 shell_diag%del_grad(2,zone_idx)*6.0d0*(hydrogen_fraction+0.32d0)/thoul_denominator)
            diffusion_coeff2(zone_idx)=settling_prefactor*(hydrogen_fraction+3.0d0)/ &
                 (5.0d0*hydrogen_fraction_sq + 8.0d0*hydrogen_fraction + 3.0d0)
            diffusion_coeff1_dx(zone_idx)=settling_prefactor*dlnp_dr(zone_idx)* &
                 ( (1.0d0-2.0d0*hydrogen_fraction-metal_fraction_total)*(1.25d0+ &
                 (6.0d0*shell_diag%del_grad(2,zone_idx)*(hydrogen_fraction+0.32d0))/thoul_denominator)+ &
                 (hydrogen_fraction-hydrogen_fraction_sq-hydrogen_metal_product)*6.0d0* &
                 shell_diag%del_grad(2,zone_idx)*(3.384d0+2.88d0*hydrogen_fraction+4.5d0*hydrogen_fraction_sq)/ &
                 thoul_denominator**2 )
            diffusion_coeff2_dx(zone_idx)=-settling_prefactor*(5.0d0*hydrogen_fraction_sq + &
                 3.0d1*hydrogen_fraction + 2.1d1)/ &
                 (5.0d0*hydrogen_fraction_sq + 8.0d0*hydrogen_fraction + 3.0d0)**2
         else
            species_mass_fraction(1) = composition(1,zone_idx)
            species_mass_fraction(2) = composition(2,zone_idx)
            species_mass_fraction(3) = composition(3,zone_idx)
            if(.not.use_thoul_fit)then
               if(coulomb_log_choice.eq.4)then
                  rho_local = exp(ln10*log_density(zone_idx))
                  temp_local = exp(ln10*log_temperature(zone_idx))
!                 calculate concentrations from mass fractions:
                  charge_to_mass_sum=0.d0
                  do species_a_idx=1,num_species-1
                   charge_to_mass_sum=charge_to_mass_sum+atomic_charge(species_a_idx)* &
                        species_mass_fraction(species_a_idx)/atomic_weight(species_a_idx)
                  enddo
                  do species_a_idx=1,num_species-1
                     concentration(species_a_idx)=species_mass_fraction(species_a_idx)/ &
                          (atomic_weight(species_a_idx)*charge_to_mass_sum)
                  enddo
                  concentration(num_species)=1.d0
!                 calculate density of electrons (NE) from mass density (RHO):
                  ac_scratch=0.d0
                  do species_a_idx=1,num_species
                   ac_scratch=ac_scratch+atomic_weight(species_a_idx)*concentration(species_a_idx)
                  enddo
                  electron_number_density=rho_local/(1.6726d-24*ac_scratch)
!                 calculate interionic distance (AO):
                  ion_number_density=0.d0
                  do species_a_idx=1,num_species-1
                     ion_number_density=ion_number_density+concentration(species_a_idx)*electron_number_density
                  enddo
                  interion_distance=(0.23873d0/ion_number_density)**cc13
!                 calculate Debye length (LAMBDAD):
                  mean_charge_sq=0.d0
                  do species_a_idx=1,num_species
                   mean_charge_sq=mean_charge_sq+concentration(species_a_idx)*atomic_charge(species_a_idx)**2
                  enddo
                  debye_length=6.9010d0*sqrt(temp_local/(electron_number_density*mean_charge_sq))
!                 calculate LAMBDA to use in Coulomb logarithm:
                  coulomb_lambda=max(debye_length,interion_distance)
!                 calculate Coulomb logarithms:
                  do species_a_idx=1,num_species
                   do species_b_idx=1,num_species
                      coulomb_xij=2.3939d3*temp_local*coulomb_lambda/ &
                           abs(atomic_charge(species_a_idx)*atomic_charge(species_b_idx))
                      coulomb_log(species_a_idx,species_b_idx)=0.81245d0 &
                             *log(1.d0+0.18769d0*coulomb_xij**1.2d0)
                   enddo
                  enddo
               else
                  do species_a_idx = 1,num_species
                     do species_b_idx = 1,num_species
                        coulomb_log(species_a_idx,species_b_idx) = ln_lambda
                     end do
                  end do
               endif
               call thdiff(num_species,atomic_weight,atomic_charge, &
                    species_mass_fraction,coulomb_log,settling_ap,settling_at,settling_ac)
               settling_coeff_p = -settling_ap(1)
               settling_coeff_t = -shell_diag%del_grad(2,zone_idx)*settling_at(1)
            else
               settling_coeff_p = 1.58d0 - 2.42d0*hydrogen_fraction + 0.844d0*hydrogen_fraction_sq
               settling_coeff_t = shell_diag%del_grad(2,zone_idx)*(1.90d0 - 2.69d0*hydrogen_fraction + 0.805d0*hydrogen_fraction_sq)
            endif
            ac_scratch = 1.15d0 - 1.42d0*hydrogen_fraction + 0.647d0*hydrogen_fraction_sq
            dap_dx = -2.42d0 + 1.688d0*hydrogen_fraction
            dat_dx = shell_diag%del_grad(2,zone_idx)*(-2.69d0 + 1.61d0*hydrogen_fraction)
            dac_dx = -1.42d0 + 1.294d0*hydrogen_fraction
!CFD 10/09 Mimic Mixing to reduce settling.
!            COD1(I) = FAC*HQPR(I)*X*(AP+AT)
            diffusion_coeff1(zone_idx) = cstmixing*settling_prefactor*dlnp_dr(zone_idx)* &
                 hydrogen_fraction*(settling_coeff_p+settling_coeff_t)
            diffusion_coeff2(zone_idx) = cstmixing*settling_prefactor*ac_scratch
            diffusion_coeff1_dx(zone_idx) = settling_prefactor*dlnp_dr(zone_idx)* &
                 (settling_coeff_p+settling_coeff_t+hydrogen_fraction*(dap_dx+dat_dx))
            diffusion_coeff2_dx(zone_idx) = settling_prefactor*dac_dx
         endif
!        METAL DIFFUSION, USING THE THOUL ET AL. COEFFICIENTS FOR FULLY
!        IONIZED IRON.
         if(use_diffusion_z)then

            settling_prefactor=fgrz*radius_bl(zone_idx)**2*temperature_bl(zone_idx)**2.5d0/ln_lambda
            if(lthoul)then
               if(use_thoul_fit)then
                  settling_coeff_p = -0.157d0 -0.511d0*hydrogen_fraction + 0.389d0*hydrogen_fraction_sq
                  settling_coeff_t = shell_diag%del_grad(2,zone_idx)*(-1.36d0 - 1.42d0*hydrogen_fraction + &
                       0.549d0*hydrogen_fraction_sq)
               else
                  settling_coeff_p = -settling_ap(3)
                  settling_coeff_t = -shell_diag%del_grad(2,zone_idx)*settling_at(3)
               endif
               iron_settling_ah = -0.0375d0 -0.193d0*hydrogen_fraction + 0.107d0*hydrogen_fraction_sq
!CFD 10/09 Mimic Mixing to reduce settling (cstmixing)
!         and add the uncertainties of differential mixing (cstdiffmix).
!
! old ver      COD1Z(I) = FAC*HQPR(I)*ZZ*(AP+AT)
               src_grid_metal_diffusion_coeff1(zone_idx) = cstdiffmix*cstmixing* &
                    settling_prefactor*dlnp_dr(zone_idx)*iron_fraction*(settling_coeff_p+settling_coeff_t)
!              POSITIVE DIFFUSION COEFFICIENTS NEEDED!
! old ver.     COD2Z(I) = ABS(FAC*AH)
! old ver.     QCOD1Z(I) = FAC*HQPR(I)*(AP+AT)
               src_grid_metal_diffusion_coeff2(zone_idx) = cstmixing*abs(settling_prefactor*iron_settling_ah)
               src_grid_metal_diffusion_coeff1_dz(zone_idx) = cstmixing*settling_prefactor* &
                    dlnp_dr(zone_idx)*(settling_coeff_p+settling_coeff_t)
               src_grid_metal_diffusion_coeff2_dz(zone_idx) = 0.0d0
            endif
         endif
!
! If using Noerdlinger's formula for LN LAMBDA, have a new term
! in D(D1)/DX
         if(coulomb_log_choice.eq.2)  then
           diffusion_coeff1_dx(zone_idx)=diffusion_coeff1_dx(zone_idx) + diffusion_coeff1(zone_idx)*1.5d0/ &
           (ln_lambda*(3.0d0*hydrogen_fraction+1.0d0))
       end if
   60 continue
 9999 continue
      return
end subroutine setup_grsett
