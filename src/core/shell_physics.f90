!----------------------------------------------------------------------
! shell_physics (legacy physic)
!----------------------------------------------------------------------
! The local physics of every zone, for the rotation code and the
! convection-zone edges: the equation of state (eos_get), the opacity
! (kap_get) and the temperature gradients (compute_temperature_gradients,
! with the rotation factors fp and ft and the adiabatic-overshoot region
! star%iov1..iov2). Sets convective_flag and stores in star% the
! gradients grad_rad, grad and grad_ad (gradr, gradT, grada), the mean
! molecular weight mu, the opacity, c_P, dlnrho/dlnT (qdt) and the
! convective velocity; then viscos gives the kinematic viscosity and the
! thermal diffusivity. log_luminosity holds L/Lsun, not its log (the
! callers pass star%luminosity_lsun).
subroutine shell_physics(fp, ft, composition, log_density, log_luminosity, &
     log_pressure, log_radius, log_mass, log_temperature, convective_flag, &
     num_zones, log_teff, ierr)
      use viscosity, only: viscos
      use convection_lib, only: compute_temperature_gradients

      use star_info_lib, only: star, json, i_h1, i_he4, i_metals, i_he3
      use physical_constants
      use eos_lib
      use kap_lib
      use numerics_lib
      use math_lib
      implicit none

      double precision, intent(in) :: fp(json), ft(json)
      double precision, intent(in) :: composition(15,json), log_density(json)
      double precision, intent(in) :: log_luminosity(json), &
           log_pressure(json), log_radius(json), log_mass(json), &
           log_temperature(json)
      logical, intent(out) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: log_teff
      integer, intent(out) :: ierr

      double precision :: atomic_weight(4)
! Indexed by species slot (i_h1, i_he4, i_metals as C12, i_he3).
      data atomic_weight /1.007825d0, 4.002603d0, 12.0d0, 3.01603d0/
! --- locals ---
      logical :: want_derivatives, in_atmosphere
      integer :: im
      double precision :: log10_mass, log10_temperature, log10_pressure, &
           log10_radius, luminosity_lsun, hydrogen_fraction, &
           metal_fraction, log10_density, pressure_rotation_factor, &
           temperature_rotation_factor
! eos_res and kap_res hold the EOS and opacity results (index constants
! in eos_lib and kap_lib). Some slots are inputs too: the density, set
! from the model before each eos_get, and the gradient and c_P guesses
! and ionization fractions, which carry over from zone to zone.
! kap_get takes eos_res%log10_density, the density eos_get just
! returned, not the model's log_density; passing log_density instead
! changes rotating runs (the byte gate catches it).
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results)
      integer :: saha_state
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
           convective_velocity
      logical :: is_convective
      double precision :: dfx1, dfx2, dfx3, dfx4, temp_scratch, amu2, emu2

      ierr = 0

! The temperature gradients, opacity and mean molecular weight of every
! zone, without derivatives, with the interior (not atmosphere) EOS.
      want_derivatives = .false.
      in_atmosphere = .false.
! The Saha solver state starts at 0 and carries over from zone to zone.
      saha_state = 0
      do im = 1,num_zones
         log10_mass = log_mass(im)
         log10_temperature = log_temperature(im)
         log10_pressure = log_pressure(im)
         log10_radius = log_radius(im)
         luminosity_lsun = log_luminosity(im)
         hydrogen_fraction = composition(i_h1,im)
         metal_fraction = composition(i_metals,im)
         log10_density = log_density(im)
         pressure_rotation_factor = fp(im)
         temperature_rotation_factor = ft(im)

         eos_res%log10_density = log10_density
         call eos_get(log10_temperature, log10_pressure, &
              hydrogen_fraction, metal_fraction, eos_res, &
              want_derivatives, in_atmosphere, saha_state, &
              composition_at_zone=composition(:,im), ierr=ierr)
         if (ierr /= 0) return
         call kap_get(eos_res%log10_density, log10_temperature, &
              hydrogen_fraction, metal_fraction, kap_res, &
              eos_res%ion_fraction, ierr=ierr)
         if (ierr /= 0) return
         call compute_temperature_gradients(log10_temperature, log10_pressure, &
              eos_res, kap_res, log10_radius, log10_mass, &
              luminosity_lsun, actual_gradient, radiative_gradient, &
              dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
              convective_velocity, want_derivatives, is_convective, &
              pressure_rotation_factor, temperature_rotation_factor, &
              log_teff, &
              im >= star%iov1 .and. im <= star%iov2, ierr)
         if (ierr /= 0) return
         convective_flag(im) = is_convective
         star%gradr(im) = radiative_gradient
         star%gradT(im) = actual_gradient
         star%grada(im) = eos_res%adiabatic_gradient
! The mean molecular weight of a fully ionized gas,
! mu = mu_ion mu_e / (mu_ion + mu_e) (amu2 = mu_ion, emu2 = mu_e):
!   1/mu_ion = amuenv + sum_k (X_k - X_k,env)/A_k over H, He4, Z (as C12)
!              and He3, where star%amuenv is 1/mu_ion of the surface
!              mixture and X_k,env its abundances;
!   1/mu_e = X/A_H + 2 (X_He3/A_He3 + Y/A_He4) + Z/2.
         dfx1 = composition(i_h1,im) - star%envelope_hydrogen_fraction
         dfx2 = composition(i_he4,im) - star%envelope_helium_fraction
         dfx3 = composition(i_metals,im) - star%envelope_metal_fraction
         dfx4 = composition(i_he3,im) - star%envelope_he3_fraction
         temp_scratch = star%amuenv + dfx1/atomic_weight(i_h1) + &
              dfx2/atomic_weight(i_he4) + dfx3/atomic_weight(i_metals) + &
              dfx4/atomic_weight(i_he3)
         amu2 = 1.0d0/temp_scratch
         temp_scratch = composition(i_h1,im)/atomic_weight(i_h1) + &
              2.0d0*(composition(i_he3,im)/atomic_weight(i_he3) + &
              composition(i_he4,im)/atomic_weight(i_he4)) + 0.5d0*composition(i_metals,im)
         emu2 = 1.0d0/temp_scratch
         star%mu(im) = amu2*emu2/(amu2+emu2)
         star%opacity_zone(im) = kap_res(i_kap)
         star%cp(im) = eos_res%specific_heat_cp
         star%qdt(im) = eos_res%dlnrho_dlnt
         star%conv_vel(im) = convective_velocity
      end do
! The thermal diffusivity and kinematic viscosity.
      call viscos(composition, log_density, log_temperature, num_zones)
end subroutine shell_physics
