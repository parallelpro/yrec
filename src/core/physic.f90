!----------------------------------------------------------------------
! physic
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original physic.f; only variable names, source form, and comment
! style were updated.
!
! For every shell, calls the equation of state (MEQOS or EQSTAT),
! GETOPAC, and TPGRAD to find the actual/adiabatic temperature
! gradients, opacity, and mean molecular weight, then calls VISCOS to
! find the thermometric diffusivity and kinematic viscosity. Finally,
! for every interface, finds the maximum allowable and actual
! d(omega)/d(ln r) for the shear instability of Endal & Sofia, ApJ
! 220:279 (1978): if the gradient of omega exceeds the critical value,
! a short-timescale instability occurs and adjacent unstable zones are
! mixed to a marginally stable state (elsewhere, not in this routine).
! KC 2025-05-31 removed the unused LCZ dummy argument (see the
! commented-out original signature below).
!       SUBROUTINE PHYSIC(FP,FT,HCOMP,HD,HG,HL,HP,HR,HS,HT,LC,LCZ,M,TEFFL)  ! KC 2025-05-31
subroutine physic(fp, ft, composition, log_density, hg, log_luminosity, &
     log_pressure, log_radius, log_mass, log_temperature, convective_flag, &
     num_zones, log_teff, ierr)

      use star_info_lib, only: star
      use star_info_lib, only: star
      use temp_lib
      use envelope_comp_lib
      use star_info_lib, only: star
      use const_lib
      use eos_lib
      use kap_lib
      use numerics_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: fp(json), ft(json)
      double precision, intent(in) :: composition(15,json), log_density(json)
      double precision, intent(in) :: hg(json)
      double precision, intent(in) :: log_luminosity(json), &
           log_pressure(json), log_radius(json), log_mass(json), &
           log_temperature(json)
      logical, intent(out) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: log_teff













      double precision :: atomic_weight(4)
      double precision :: log_mass_nodes(4), interp_weights(4)
      data atomic_weight /1.007825d0, 4.002603d0, 12.0d0, 3.01603d0/
      save

! --- locals ---
      logical :: want_derivatives, local_conductive_opacity_flag, &
           in_atmosphere
      integer :: i, im, k
      double precision :: log10_mass, log10_temperature, log10_pressure, &
           log10_radius, luminosity_lsun, hydrogen_fraction, &
           metal_fraction, log10_density, pressure_rotation_factor, &
           temperature_rotation_factor
      double precision :: temperature, pressure, density, beta, &
           beta_inverse, beta14, ion_fraction(3), specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      integer :: saha_state
      double precision :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
           convective_velocity
      logical :: is_convective
      double precision :: dfx1, dfx2, dfx3, dfx4, temp_scratch, amu2, emu2
      double precision :: log_mass_mid, pressure_mid, density_mid, &
           actual_grad_mid, adiabatic_grad_mid, gravity_mid

!  FIND ACTUAL AND ADIABATIC TEMPERATURE GRADIENTS,OPACITY,AND
!  MEAN MOLECULAR WEIGHT FOR ALL RADIATIVE SHELLS.
      integer, intent(out) :: ierr

      ierr = 0

      want_derivatives = .false.
      local_conductive_opacity_flag = .false.
      in_atmosphere = .false.
      idt = 15
      do 25 i = 1,4
         idd(i) = 5
   25 continue
      do 30 im = 1,num_zones
         log10_mass = log_mass(im)
         log10_temperature = log_temperature(im)
         log10_pressure = log_pressure(im)
         log10_radius = log_radius(im)
         luminosity_lsun = log_luminosity(im)
         hydrogen_fraction = composition(1,im)
         metal_fraction = composition(3,im)
         log10_density = log_density(im)
         pressure_rotation_factor = fp(im)
         temperature_rotation_factor = ft(im)

         call eos_get(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, &
              electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, &
              adiabatic_gradient_dp, specific_heat_cp_dt, &
              specific_heat_cp_dp, want_derivatives, in_atmosphere, &
              saha_state, composition_at_zone=composition(:,im))
         call kap_get(log10_density, log10_temperature, hydrogen_fraction, &
              metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, &
              dlnkap_dlnt, ion_fraction)
         iovim = im
         call tpgrad(log10_temperature, temperature, log10_pressure, &
              pressure, density, log10_radius, log10_mass, &
              luminosity_lsun, opacity, dlnrho_dlnt, dlnrho_dlnp, &
              dlnkap_dlnt, dlnkap_dlnrho, specific_heat_cp, &
              actual_gradient, radiative_gradient, adiabatic_gradient, &
              dlnrho_dlnt_dt, dlnrho_dlnp_dt, adiabatic_gradient_dt, &
              adiabatic_gradient_dp, dgrad_dt_component, &
              dgrad_dp_component, dgrad_dr_component, specific_heat_cp_dt, &
              specific_heat_cp_dp, convective_velocity, &
              want_derivatives, is_convective, pressure_rotation_factor, &
              temperature_rotation_factor, log_teff, ierr)
         if (ierr /= 0) return
         convective_flag(im) = is_convective
         star%diag%del_grad(1,im) = radiative_gradient
         star%diag%del_grad(2,im) = actual_gradient
         star%diag%del_grad(3,im) = adiabatic_gradient
!  FIND NEW RUN OF MEAN MOLECULAR WEIGHT ASSUMING FULLY IONIZED GAS.
!  AMUENV IS(1/MEAN MOLECULAR WEIGHT PER ION OF THE SURFACE MIXTURE.)
         dfx1 = composition(1,im) - env_comp%envelope_hydrogen_fraction
         dfx2 = composition(2,im) - star%run%envelope_helium_fraction
         dfx3 = composition(3,im) - env_comp%envelope_metal_fraction
         dfx4 = composition(4,im) - star%run%envelope_he3_fraction
         temp_scratch = env_comp%amuenv + dfx1/atomic_weight(1) + &
              dfx2/atomic_weight(2) + dfx3/atomic_weight(3) + &
              dfx4/atomic_weight(4)
         amu2 = 1.0d0/temp_scratch
         temp_scratch = composition(1,im)/atomic_weight(1) + &
              2.0d0*(composition(4,im)/atomic_weight(4) + &
              composition(2,im)/atomic_weight(2)) + 0.5d0*composition(3,im)
         emu2 = 1.0d0/temp_scratch
         shell_temp%mean_molecular_weight(im) = amu2*emu2/(amu2+emu2)
         star%diag%so(im) = opacity
         shell_temp%cp(im) = specific_heat_cp
         shell_temp%qdt(im) = dlnrho_dlnt
! JVS 10/13 Always want SVEL
         star%diag%svel(im) = convective_velocity
!         IF(LC(IM))THEN
!            SVEL(IM) = VEL
!         ELSE
!            SVEL(IM) = 0.0D0
!         ENDIF
   30 continue
!  FIND THE THERMOMETRIC DIFFUSIVITY AND KINEMATIC VISCOSITY.
!       CALL VISCOS(HCOMP,HD,HT,LC,M)  ! KC 2025-05-31
      call viscos(composition, log_density, log_temperature, num_zones)
!  FIND THE RUN OF MAXIMUM ALLOWABLE D OMEGA/D LNR AND THE RUN OF ACTUAL
!  D OMEGA/D LN R FOR ALL RADIATIVE ZONES.
!  IF THE GRADIENT OF OMEGA EXCEEDS THE CRITICAL VALUE FOR THE SHEAR
!  INSTABILITY, A SHORT TIMESCALE INSTABILITY OCCURS.  MIX ALL ADJACENT
!  UNSTABLE ZONES TO A MARGINALLY STABLE STATE.
      do 100 im = 2,num_zones
!  SKIP CONVECTIVE REGIONS
         if (convective_flag(im).and.convective_flag(im-1)) then
            star%rot%max_domega_dr(im) = 0.0d0
            goto 100
         end if
!  NOW CHECK FOR SHEAR INSTABILITY -REF.ENDAL&SOFIA APJ 220:279(1978)
!  THERMODYNAMIC QUANTITIES ARE CALCULATED AT THE SHELL MIDPOINT BY
!  4-POINT LAGRANGIAN INTERPOLATION.
         if (im.le.2) then
            k = 1
         else if (im.eq.num_zones) then
            k = num_zones - 3
         else
            k = im - 2
         end if
         do 90 i = 1, 4
            log_mass_nodes(i) = log_mass(i+k-1)
   90    continue
!  USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND PHYSICAL VARIABLES
!  AT THE INTERFACE BEING TESTED.
         log_mass_mid = 0.5d0*(log_mass(im) + log_mass(im-1))
         call intrp2(log_mass_nodes, interp_weights, log_mass_mid)
         pressure_mid = log_pressure(k)*interp_weights(1) + &
              log_pressure(k+1)*interp_weights(2) + &
              log_pressure(k+2)*interp_weights(3) + &
              log_pressure(k+3)*interp_weights(4)
         density_mid = log_density(k)*interp_weights(1) + &
              log_density(k+1)*interp_weights(2) + &
              log_density(k+2)*interp_weights(3) + &
              log_density(k+3)*interp_weights(4)
         actual_grad_mid = star%diag%del_grad(2,k)*interp_weights(1) + &
              star%diag%del_grad(2,k+1)*interp_weights(2) + &
              star%diag%del_grad(2,k+2)*interp_weights(3) + &
              star%diag%del_grad(2,k+3)*interp_weights(4)
         adiabatic_grad_mid = star%diag%del_grad(3,k)*interp_weights(1) + &
              star%diag%del_grad(3,k+1)*interp_weights(2) + &
              star%diag%del_grad(3,k+2)*interp_weights(3) + &
              star%diag%del_grad(3,k+3)*interp_weights(4)
         gravity_mid = hg(k)*interp_weights(1) + hg(k+1)*interp_weights(2) + &
              hg(k+2)*interp_weights(3) + hg(k+3)*interp_weights(4)
         temp_scratch = dexp(ln10*(density_mid - pressure_mid))* &
              (adiabatic_grad_mid - actual_grad_mid)*gravity_mid**2
         if (temp_scratch.gt.0.0d0) then
            star%rot%max_domega_dr(im) = 2.0d0*dsqrt(temp_scratch)
         else
            star%rot%max_domega_dr(im) = 0.0d0
         end if
  100 continue

      return
end subroutine physic
