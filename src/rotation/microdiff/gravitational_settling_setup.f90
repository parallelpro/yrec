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
! gravitational_settling.f90: locates the diffusion zone boundaries (excluding any
! convective core/envelope and any hydrogen- or helium-exhausted
! region), converts the model to Bahcall & Loeb (1990) natural units,
! and computes the helium (and, if enabled, iron) settling diffusion
! coefficients diffusion_coeff1/diffusion_coeff2 and their partial
! derivatives with respect to X (diffusion_coeff1_dx/diffusion_
! coeff2_dx) at every zone, using either the simple analytic fit
! (LTHOUL false) or the full Thoul et al. (1994) coefficients (via
! thoul_diffusion, LTHOUL true) with a choice of Coulomb-logarithm prescription
! selected by coulomb_log_choice.
subroutine gravitational_settling_setup(timestep_seconds, dlnp_dr, log_radius, &
     log_density, mass_grams, log_temperature, del_grad, convective_flag, &
     num_zones, total_mass, diffusion_coeff1, diffusion_coeff2, &
     composition, radius_bl, temperature_bl, zone_begin, zone_end, &
     settling_skipped_flag, diffusion_coeff1_dx, diffusion_coeff2_dx, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star, json, i_h1, i_he4, i_metals, i_he3
      use species_table_lib, only: thoul_h1, thoul_he4, thoul_fe, thoul_electron, &
           thoul_col_h, thoul_col_he, thoul_col_metal
      use bahcall_loeb_units_lib, only: set_bahcall_loeb_scales
      use luout_lib
      use run_log_lib, only: solver_diagnostics
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(inout) :: timestep_seconds
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json), &
           log_temperature(json)
! del_grad: DEL (=DLNT/DLNP) per zone, from the caller (was star%gradT)
      double precision, intent(in) :: del_grad(json)
      double precision, intent(inout) :: mass_grams(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      double precision, intent(out) :: diffusion_coeff1(json), &
           diffusion_coeff2(json)
      double precision, intent(in) :: composition(15,json)
      double precision, intent(out) :: radius_bl(json), temperature_bl(json)
      integer, intent(out) :: zone_begin, zone_end
      logical, intent(out) :: settling_skipped_flag
      integer, intent(out) :: ierr
      double precision, intent(out) :: diffusion_coeff1_dx(json), &
           diffusion_coeff2_dx(json)
! MHP 8/94 ADDED ATOMIC WEIGHTS AND CHARGES FOR H,HE,FE,ELECTRONS -
! NEEDED FOR FULL THOUL COEFFICIENTS (VALUES FROM species_table_lib)
      double precision :: atomic_weight(4), atomic_charge(4), &
           settling_ap(4), settling_at(4), settling_ac(4,4), &
           coulomb_log(4,4), species_mass_fraction(4)
!     additional variables to be declared for variable ln lambda
!     intermediate variables:
! mass_weighted_conc (originally AC): concentration-weighted atomic
! weight sum used to get the electron density (coulomb_log_choice==4).
      double precision :: charge_to_mass_sum, mass_weighted_conc, &
           ion_number_density, mean_charge_sq, coulomb_xij, &
           electron_number_density, interion_distance, debye_length, &
           coulomb_lambda, concentration(4)
      double precision :: ln_lambda
      data atomic_weight/thoul_h1%weight,thoul_he4%weight,thoul_fe%weight, &
           thoul_electron%weight/
      data atomic_charge/thoul_h1%charge,thoul_he4%charge,thoul_fe%charge, &
           thoul_electron%charge/
      integer :: num_species
      data num_species/4/
! --- locals ---
      integer :: zone_idx
      double precision :: hydrogen_fraction, z_plus_he3_fraction, &
           metal_fraction, hydrogen_fraction_sq, &
           hydrogen_metal_product, thoul_denominator
! helium_ah_coeff (also originally AC): empirical Thoul-fit helium-
! settling amplitude coefficient feeding diffusion_coeff2.
      double precision :: settling_prefactor, iron_settling_ah, helium_ah_coeff
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
!     settling_skipped_flag=T IF SETTLING IS SKIPPED THIS MODEL (FULLY
!     CONVECTIVE, HYDROGEN-EXHAUSTED OR HELIUM-EXHAUSTED).
      ierr = 0
      settling_skipped_flag=.false.
!     CHECK FOR CONVECTIVE CORE.
      if(convective_flag(1))then
         do zone_idx=2,num_zones
            if(.not.convective_flag(zone_idx))exit
         end do
         if (zone_idx > num_zones) then
!        DIFFUSION NOT COMPUTED FOR FULLY CONVECTIVE MODELS.
         settling_skipped_flag=.true.
! print once per suspension; every model only under
! report_solver_diagnostics (2026 run-log verbosity sweep)
         if (solver_diagnostics() .or. &
              .not. star%settling_suspended_reported) then
            write(run_log_unit,15)
   15       format(1x,' FULLY CONVECTIVE MODEL - NO SETTLING')
            star%settling_suspended_reported = .true.
         end if
         return
         end if
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_begin = zone_idx-1
      else
         zone_begin = 1
      endif
! MHP 6/90 CHECK FOR HYDROGEN-EXHAUSTED CORE.
      do zone_idx = zone_begin,num_zones
         if(composition(i_h1,zone_idx).gt.star%ctrl%hydrogen_diffusion_floor)exit
      end do
      if (zone_idx > num_zones) then
!     HYDROGEN-FREE MODEL - EXIT.
! print once per suspension; every model only under
! report_solver_diagnostics (2026 run-log verbosity sweep)
      if (solver_diagnostics() .or. &
           .not. star%settling_suspended_reported) then
         write(run_log_unit,16)star%ctrl%hydrogen_diffusion_floor
   16    format(1x,'X BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
         star%settling_suspended_reported = .true.
      end if
      settling_skipped_flag = .true.
      return
      end if
      zone_begin = zone_idx
!     CHECK FOR CONVECTIVE ENVELOPE.
      if(convective_flag(num_zones))then
         do zone_idx=num_zones-1,2,-1
            if(.not.convective_flag(zone_idx))exit
         end do
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_end = zone_idx+1
      else
         zone_end = num_zones
      endif
!     CHECK FOR HELIUM-EXHAUSTED SURFACE.
!     OUTER POINT IS SET WHEREVER Y>YMIN.
      do zone_idx=zone_end,1,-1
         if(composition(i_he4,zone_idx).gt.star%ctrl%helium_diffusion_min) exit
      end do
      if (zone_idx < 1) then
!     HELIUM-EXHAUSTED MODEL - EXIT.
! print once per suspension; every model only under
! report_solver_diagnostics (2026 run-log verbosity sweep)
      if (solver_diagnostics() .or. &
           .not. star%settling_suspended_reported) then
         write(run_log_unit,17)star%ctrl%helium_diffusion_min
   17    format(1x,'Y BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
         star%settling_suspended_reported = .true.
      end if
      settling_skipped_flag = .true.
      return
      end if
      zone_end = zone_idx
! all suspension checks passed: settling proceeds this model
      if (star%settling_suspended_reported) then
         write(run_log_unit,916)
  916    format(1x,' SETTLING RESUMED')
         star%settling_suspended_reported = .false.
      end if
!     SET THE star%bl_*_scale CONVERSION FACTORS.
      call set_bahcall_loeb_scales()
!     CONVERT LOG(RADIUS) AND LOG(TEMPERATURE) TO NATURAL UNITS.
!     ALSO CONVERT NATURAL UNITS TO BAHCALL AND LOEB UNITS.
      do zone_idx=1,num_zones
         radius_bl(zone_idx)=exp(ln10*log_radius(zone_idx))*star%bl_radius_scale
         temperature_bl(zone_idx)=exp(ln10*log_temperature(zone_idx))*star%bl_temp_scale
         mass_grams(zone_idx)=mass_grams(zone_idx)*star%bl_mass_scale
         dlnp_dr(zone_idx)=dlnp_dr(zone_idx)/star%bl_radius_scale
      end do
      timestep_seconds=timestep_seconds/star%bl_time_scale
      total_mass=total_mass*star%bl_mass_scale
!     SET UP DIFFUSION COEFFICIENTS.
!     MODIFIED BY BC MAY/90 -- VALID FOR ALL X WITH VARIABLE LN LAMBDA
!     GENERAL EQUATION IS
!     DX/DT = (D/DR(D1)+D/DR(D2 DX/DR))/(RHO*R**2), WHERE
!     D1 = R**2/LN LAMBDA  * X  * T**5/2 * (DLNP/DR) * (1-X) *
!          [5/4 + DEL*6*(X-0.32)/(5.4+6.3X-4.5X**2)]
!     D2 = R**2/LN LAMBDA * T**5/2 * (X+3)/(X+1)/(3+5X)
      do zone_idx = 1,num_zones
         hydrogen_fraction = composition(i_h1,zone_idx)
! MHP 10/02 INITIALIZED X - WAS NOT DONE PRIOR TO USAGE IN SHELL 1
         if(star%ctrl%coulomb_log_choice.eq.2)then
!           Noerdlinger's formula (1977 A&A 57,407) for LN LAMBDA:
!           Ln Lambda = -19.7 - ln[4/(3*X + 1)]/2 - ln(rho)/2 + 1.5*ln(T)
            ln_lambda=-1.97d1 - 0.5d0*log(1.0d0/(0.75d0*hydrogen_fraction+0.25d0)) &
                      -0.5d0*log_density(zone_idx)*ln10 + 1.5d0*log_temperature(zone_idx)*ln10
         else if(star%ctrl%coulomb_log_choice.eq.3)then
!           Loeb's formula (1989 Phys. Rev. D 39, 1009) for LN LAMBDA
!           ln lambda = -19.105747 - ln(rho)/2 + 1.5*ln(T)
            ln_lambda=-1.9105747d1-0.50d0*log_density(zone_idx)*ln10+1.5d0*log_temperature(zone_idx)*ln10
         else
            ln_lambda = 2.2d0
         end if
!
         settling_prefactor=star%ctrl%fgry*radius_bl(zone_idx)**2*pow(temperature_bl(zone_idx), 2.5d0)/ln_lambda
         z_plus_he3_fraction = composition(i_metals,zone_idx) + composition(i_he3,zone_idx)
         metal_fraction = composition(i_metals,zone_idx)
         hydrogen_fraction_sq = hydrogen_fraction*hydrogen_fraction
         if(.not.star%ctrl%use_thoul_diffusion)then
            hydrogen_metal_product = hydrogen_fraction*z_plus_he3_fraction
            thoul_denominator=5.4d0+6.3d0*hydrogen_fraction-4.5d0*hydrogen_fraction_sq
            diffusion_coeff1(zone_idx)=settling_prefactor*dlnp_dr(zone_idx)* &
                 (hydrogen_fraction - hydrogen_fraction_sq - hydrogen_metal_product)*(1.25d0+ &
                 del_grad(zone_idx)*6.0d0*(hydrogen_fraction+0.32d0)/thoul_denominator)
            diffusion_coeff2(zone_idx)=settling_prefactor*(hydrogen_fraction+3.0d0)/ &
                 (5.0d0*hydrogen_fraction_sq + 8.0d0*hydrogen_fraction + 3.0d0)
            diffusion_coeff1_dx(zone_idx)=settling_prefactor*dlnp_dr(zone_idx)* &
                 ( (1.0d0-2.0d0*hydrogen_fraction-z_plus_he3_fraction)*(1.25d0+ &
                 (6.0d0*del_grad(zone_idx)*(hydrogen_fraction+0.32d0))/thoul_denominator)+ &
                 (hydrogen_fraction-hydrogen_fraction_sq-hydrogen_metal_product)*6.0d0* &
                 del_grad(zone_idx)*(3.384d0+2.88d0*hydrogen_fraction+4.5d0*hydrogen_fraction_sq)/ &
                 thoul_denominator**2 )
            diffusion_coeff2_dx(zone_idx)=-settling_prefactor*(5.0d0*hydrogen_fraction_sq + &
                 3.0d1*hydrogen_fraction + 2.1d1)/ &
                 (5.0d0*hydrogen_fraction_sq + 8.0d0*hydrogen_fraction + 3.0d0)**2
         else
            species_mass_fraction(thoul_col_h) = composition(i_h1,zone_idx)
            species_mass_fraction(thoul_col_he) = composition(i_he4,zone_idx)
            species_mass_fraction(thoul_col_metal) = composition(i_metals,zone_idx)
            if(.not.star%ctrl%use_thoul_fit)then
               if(star%ctrl%coulomb_log_choice.eq.4)then
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
                  mass_weighted_conc=0.d0
                  do species_a_idx=1,num_species
                   mass_weighted_conc=mass_weighted_conc+atomic_weight(species_a_idx)*concentration(species_a_idx)
                  enddo
                  electron_number_density=rho_local/(m_proton_cgs*mass_weighted_conc)
!                 calculate interionic distance (AO):
                  ion_number_density=0.d0
                  do species_a_idx=1,num_species-1
                     ion_number_density=ion_number_density+concentration(species_a_idx)*electron_number_density
                  enddo
                  interion_distance=pow((0.23873d0/ion_number_density), cc13)
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
                             *log(1.d0+0.18769d0*pow(coulomb_xij, 1.2d0))
                   enddo
                  enddo
               else
                  do species_a_idx = 1,num_species
                     do species_b_idx = 1,num_species
                        coulomb_log(species_a_idx,species_b_idx) = ln_lambda
                     end do
                  end do
               endif
               call thoul_diffusion(num_species,atomic_weight,atomic_charge, &
                    species_mass_fraction,coulomb_log,settling_ap,settling_at,settling_ac, &
                    ierr)
               if (ierr /= 0) return
               settling_coeff_p = -settling_ap(thoul_col_h)
               settling_coeff_t = -del_grad(zone_idx)*settling_at(thoul_col_h)
            else
               settling_coeff_p = 1.58d0 - 2.42d0*hydrogen_fraction + 0.844d0*hydrogen_fraction_sq
               settling_coeff_t = del_grad(zone_idx)*(1.90d0 - 2.69d0*hydrogen_fraction + 0.805d0*hydrogen_fraction_sq)
            endif
            helium_ah_coeff = 1.15d0 - 1.42d0*hydrogen_fraction + 0.647d0*hydrogen_fraction_sq
            dap_dx = -2.42d0 + 1.688d0*hydrogen_fraction
            dat_dx = del_grad(zone_idx)*(-2.69d0 + 1.61d0*hydrogen_fraction)
            dac_dx = -1.42d0 + 1.294d0*hydrogen_fraction
!CFD 10/09 Mimic Mixing to reduce settling.
            diffusion_coeff1(zone_idx) = star%ctrl%constant_mixing_coeff*settling_prefactor*dlnp_dr(zone_idx)* &
                 hydrogen_fraction*(settling_coeff_p+settling_coeff_t)
            diffusion_coeff2(zone_idx) = star%ctrl%constant_mixing_coeff*settling_prefactor*helium_ah_coeff
            diffusion_coeff1_dx(zone_idx) = settling_prefactor*dlnp_dr(zone_idx)* &
                 (settling_coeff_p+settling_coeff_t+hydrogen_fraction*(dap_dx+dat_dx))
            diffusion_coeff2_dx(zone_idx) = settling_prefactor*dac_dx
         endif
!        METAL DIFFUSION, USING THE THOUL ET AL. COEFFICIENTS FOR FULLY
!        IONIZED IRON.
         if(star%job%use_diffusion_z)then
            settling_prefactor=star%job%fgrz*radius_bl(zone_idx)**2*pow(temperature_bl(zone_idx), 2.5d0)/ln_lambda
            if(star%ctrl%use_thoul_diffusion)then
               if(star%ctrl%use_thoul_fit)then
                  settling_coeff_p = -0.157d0 -0.511d0*hydrogen_fraction + 0.389d0*hydrogen_fraction_sq
                  settling_coeff_t = del_grad(zone_idx)*(-1.36d0 - 1.42d0*hydrogen_fraction + &
                       0.549d0*hydrogen_fraction_sq)
               else
                  settling_coeff_p = -settling_ap(thoul_col_metal)
                  settling_coeff_t = -del_grad(zone_idx)*settling_at(thoul_col_metal)
               endif
               iron_settling_ah = -0.0375d0 -0.193d0*hydrogen_fraction + 0.107d0*hydrogen_fraction_sq
!CFD 10/09 Mimic Mixing to reduce settling (constant_mixing_coeff)
!         and add the uncertainties of differential mixing (constant_settling_reduction).
               rot_scr%src_grid_metal_diffusion_coeff1(zone_idx) = star%ctrl%constant_settling_reduction*star%ctrl%constant_mixing_coeff* &
                    settling_prefactor*dlnp_dr(zone_idx)*metal_fraction*(settling_coeff_p+settling_coeff_t)
!              POSITIVE DIFFUSION COEFFICIENTS NEEDED!
               rot_scr%src_grid_metal_diffusion_coeff2(zone_idx) = star%ctrl%constant_mixing_coeff*abs(settling_prefactor*iron_settling_ah)
               rot_scr%src_grid_metal_diffusion_coeff1_dz(zone_idx) = star%ctrl%constant_mixing_coeff*settling_prefactor* &
                    dlnp_dr(zone_idx)*(settling_coeff_p+settling_coeff_t)
               rot_scr%src_grid_metal_diffusion_coeff2_dz(zone_idx) = 0.0d0
            endif
         endif
!
! If using Noerdlinger's formula for LN LAMBDA, have a new term
! in D(D1)/DX
         if(star%ctrl%coulomb_log_choice.eq.2)  then
           diffusion_coeff1_dx(zone_idx)=diffusion_coeff1_dx(zone_idx) + diffusion_coeff1(zone_idx)*1.5d0/ &
           (ln_lambda*(3.0d0*hydrogen_fraction+1.0d0))
       end if
      end do
      return
end subroutine gravitational_settling_setup
