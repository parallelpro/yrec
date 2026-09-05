!----------------------------------------------------------------------
! viscos
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original viscos.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Computes the kinematic microscopic viscosity (radiative + molecular)
! and the radiative thermal diffusivity at every zone, for use by the
! rotational-mixing/instability diffusion routines.
subroutine viscos(composition, log_density, log_temperature, num_zones)
      use star_info_lib, only: star, json, i_metals
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_temperature(json)
      integer, intent(in) :: num_zones

      double precision :: amu
! Kinetic-theory weights and charges for composition rows 1..11
! (i_h1..i_o18 of star_info_lib); row i_metals is the lumped Z, skipped
! in the molecular-viscosity sums below.  These values differ from the
! Thoul-solver descriptors in species_table_lib and stay local.
      double precision :: weight(11), z(11)
      data amu/amu_cgs_legacy/
      data weight/1.007825d0,4.0026d0,1.0d0,3.01603d0,12.0d0, &
           13.00335d0,14.00307d0,15.0011d0,15.99491d0,16.99491d0, &
           17.99491d0/
      data z/1.0d0,2.0d0,0.0d0,2.0d0,6.0d0,6.0d0,7.0d0,7.0d0,8.0d0, &
           8.0d0,8.0d0/
! --- locals ---
      integer :: shell_idx, species_idx, species_idx2
      double precision :: opacity_local, mean_charge, number_density_sum
      double precision :: temperature_cgs, temperature_sq, density_cgs, &
           electron_number_density
      double precision :: viscosity_radiative
      double precision :: number_density(11)
      double precision :: viscosity_molecular_species(11)
      double precision :: coulomb_log_factor
      double precision :: mfp_temperature_factor, molecular_coeff, &
           viscosity_molecular
      double precision :: species_coeff, species_sum

!  LOOP OVER ALL ZONES (CONVECTIVE ZONES INCLUDED)
      do shell_idx = 1,num_zones
!  COMPUTE THE KINEMATIC MICROSCOPIC VISCOSITY DUE TO RADIATION AND IONS
!  CONVERT TO NUMBER DENSITIES AND FIND MEAN CHARGE PER ION(ZF) AND NE.
         opacity_local = star%opacity_zone(shell_idx)
         mean_charge = 0.0d0
         number_density_sum = 0.0d0
         do species_idx = 1,11
            number_density(species_idx) = composition(species_idx,shell_idx)/ &
                 weight(species_idx)
            number_density_sum = number_density_sum+number_density(species_idx)
            mean_charge = mean_charge+number_density(species_idx)*z(species_idx)
         end do
         mean_charge = mean_charge/number_density_sum
         temperature_cgs = exp(ln10*log_temperature(shell_idx))
         temperature_sq = temperature_cgs**2
         density_cgs = exp(ln10*log_density(shell_idx))
         electron_number_density = mean_charge*density_cgs/amu
!  RADIATIVE DYNAMIC VISCOSITY (ELECTRON SCATTERING ONLY):
!  METHOD USED IS FROM LEDOUX,1958,HANDBUCH DER PHYSIK VOL.LI,P.445
!  (THE ORIGINAL ALSO EVALUATED THE ENDAL-SOFIA / THOMAS 1930 FORMULAE
!  FOR COMPARISON; THOSE VALUES WERE NEVER USED AND ARE NOT COMPUTED.)
         viscosity_radiative = rad_viscosity_coeff_cgs*temperature_sq*temperature_sq/ &
              (opacity_local*density_cgs**2)
!  MOLECULAR DYNAMIC VISCOSITY
!  REF. SPITZER,1962,PHYSICS OF FULLY IONIZED GASES
         mfp_temperature_factor = dmin1(1.0d0,4.2d5/temperature_cgs)
         coulomb_log_factor = 9.424536845d0+0.5d0*log(temperature_sq* &
              temperature_cgs*mfp_temperature_factor/ &
              (electron_number_density*mean_charge**2))
         molecular_coeff = 3.125d-15*dsqrt(temperature_cgs)*temperature_sq
         viscosity_molecular = 0.0d0
!  VISCX(I) IS THE MOLECULAR VISCOSITY OF SPECIES I.
         do species_idx = 1,11
            if(species_idx.eq.i_metals) cycle
            species_coeff = molecular_coeff*number_density(species_idx)* &
                 dsqrt(weight(species_idx))/ &
                 ((coulomb_log_factor-log(z(species_idx)))*z(species_idx)**2)
            species_sum = 0.0d0
            do species_idx2 = 1,11
               if(species_idx2.eq.i_metals) cycle
               species_sum = species_sum+number_density(species_idx2)* &
                    z(species_idx2)**2* &
                    dsqrt((weight(species_idx)+weight(species_idx2))/ &
                    weight(species_idx2))
            end do
            if(dabs(species_sum).lt.1.0d-38) cycle
            viscosity_molecular_species(species_idx) = species_coeff/species_sum
            if(viscosity_molecular_species(species_idx).gt.0.0d0) &
                 viscosity_molecular = viscosity_molecular+ &
                 viscosity_molecular_species(species_idx)
         end do
         viscosity_molecular = viscosity_molecular/density_cgs
         star%visc(shell_idx) = viscosity_radiative+viscosity_molecular
!  THERMAL DIFFUSIVITY(THDIF) DUE TO RADIATION IS CALCULATED
!  COMPONENT DUE TO THERMAL CONDUCTION OF MATTER IS NEGLECTED
!  RADIATIVE DIFFUSIVITY = K*T**3/(O*RHO**2*CP)
         star%thdif(shell_idx) = 1.6d1*cc13*sigma_sb_cgs_legacy*temperature_cgs* &
              temperature_sq/(opacity_local*density_cgs**2*star%cp(shell_idx))
      end do

      return
end subroutine viscos
