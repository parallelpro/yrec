!----------------------------------------------------------------------
! semiconvection
!----------------------------------------------------------------------
! Semiconvection at the edge of a mixed region (Castellani et al. 1971;
! formerly sconvec): extend_by_semiconvection widens the
! first mixed region where the radiative gradient, with the composition
! of the region, exceeds the adiabatic one.
module semiconvection
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use convection_configuration, only: convection_settings
      use composition_layout, only: i_h1, i_metals
      use temperature_gradients, only: compute_temperature_gradients
      use convection_zones, only: merge_overlapping_zones
      use convection_messages, only: log_semiconvective_edge, merged_by_semiconvection
      use eos_lib, only: eos_result, eos_get
      use kap_lib, only: kap_get, num_kap_results
      use physical_constants, only: ln10, c4pi
      use math_lib
      implicit none
      private
      public :: extend_by_semiconvection

! The two edges of a mixed region: the columns of a zone list.
      integer, parameter :: lower_edge = 1, upper_edge = 2
! The 10 of the instability's propagation speed (Castellani et al. 1971,
! eq. 5', p. 347): v_p = (1 - del_ad/del_rad) L/(10 4 pi P r^2).
      real(dp), parameter :: propagation_speed_divisor = 10.0d0

contains

!----------------------------------------------------------------------
! extend_by_semiconvection (formerly sconvec)
!----------------------------------------------------------------------
! Semiconvection at the edges of the first mixed region of the zone list
! mixed_zone_bounds (num_zones_mixed regions, from find_convection_zones):
! its upper edge when it is the convective core, its lower edge when it
! is the surface convection zone, both otherwise; nothing for a fully
! convective model (extend_edge says how an edge moves). Regions that
! then overlap are merged. Only the first region is examined, the YREC
! team's core-only switch (a8e60674, 2025), whatever num_zones_mixed is.
!
! In: the timestep (s); the composition, log10 rho, P, r, m and T and
! L/Lsun (not a log) by shell; log10 Teff (the starspot contrast). The
! shells semiconvection evaluates get the density the EOS returns for
! their own composition in log_density. Out: mixed_zone_bounds and
! num_zones_mixed; ierr, an EOS or opacity failure's code, or the
! mixing-length square root's (compute_temperature_gradients).
subroutine extend_by_semiconvection(timestep, composition, log_density, luminosity_lsun, &
     log_pressure, log_radius, log_mass, log_temperature, num_zones, &
     mixed_zone_bounds, num_zones_mixed, log_teff, ierr)
      real(dp), intent(in) :: timestep, composition(:,:)
      real(dp), intent(inout) :: log_density(:)
      real(dp), intent(in) :: luminosity_lsun(:), log_pressure(:), log_radius(:), &
           log_mass(:), log_temperature(:)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: mixed_zone_bounds(:,:), num_zones_mixed
      real(dp), intent(in) :: log_teff
      integer, intent(out) :: ierr
! The EOS's result, kept from one evaluation to the next (each seeds
! only the density).
      type(eos_result) :: eos_res
      logical :: extend(lower_edge:upper_edge)
      integer :: edge

      ierr = 0
      if (mixed_zone_bounds(1,1) == 1) then
         if (mixed_zone_bounds(1,2) == num_zones) return
         extend = [.false., .true.]
      else if (mixed_zone_bounds(1,2) == num_zones) then
         extend = [.true., .false.]
      else
         extend = [.true., .true.]
      end if
      do edge = lower_edge, upper_edge
         if (.not. extend(edge)) cycle
         call extend_edge(edge, timestep, composition, log_density, luminosity_lsun, &
              log_pressure, log_radius, log_mass, log_temperature, num_zones, log_teff, &
              mixed_zone_bounds(1,edge), eos_res, ierr)
         if (ierr /= 0) return
      end do
      call merge_overlapping_zones(mixed_zone_bounds, num_zones_mixed, &
           merge_shared_shell=.false., cause=merged_by_semiconvection)
end subroutine extend_by_semiconvection

!----------------------------------------------------------------------
! extend_edge: one edge of a mixed region, edge_shell, moved outward
! (down for lower_edge, up for upper_edge) and logged.
! 1. The first shell outside the region, with the region's composition:
!    if its radiative gradient del_rad,mix is below the adiabatic one,
!    the edge stays. Otherwise the same shell with its own composition
!    gives del_rad,local, and f = del_rad,mix/del_rad,local.
! 2. The edge moves out shell by shell while each shell, with its own
!    composition, has f del_rad > del_ad, and while the distance covered,
!    weighted by the jump in the inverse mean molecular weight,
!      sum (1 - mu_mix^-1/mu^-1) dr,
!    stays within v_p dt, the distance the instability propagates in the
!    timestep, v_p = (1 - del_ad/del_rad,mix) L/(10 4 pi P r^2) with L,
!    P and r at the edge (mu^-1 is the last evaluated shell's). Reaching
!    that distance or the model's end is the limit.
! Every shell of the search is evaluated at the first outside shell's
! radius, as in the F77 (the radius enters only the mixing-length cubic,
! not del_rad).
subroutine extend_edge(edge, timestep, composition, log_density, luminosity_lsun, &
     log_pressure, log_radius, log_mass, log_temperature, num_zones, log_teff, &
     edge_shell, eos_res, ierr)
      integer, intent(in) :: edge, num_zones
      real(dp), intent(in) :: timestep, composition(:,:), luminosity_lsun(:), &
           log_pressure(:), log_radius(:), log_mass(:), log_temperature(:), log_teff
      real(dp), intent(inout) :: log_density(:)
      integer, intent(inout) :: edge_shell
      type(eos_result), intent(inout) :: eos_res
      integer, intent(out) :: ierr
      integer :: step, outside, last_shell, shell, new_edge, saha_iterations
      real(dp) :: log_density_outside, log_density_shell, gradient_zone_mix, &
           gradient_local, gradient, gradient_ratio, inverse_mu_zone_mix, &
           superadiabatic_fraction, max_distance, distance, radius_prev, &
           radius_curr, radius_step
      logical :: reached_limit

      if (edge == lower_edge) then
         step = -1
         last_shell = 1
      else
         step = 1
         last_shell = num_zones
      end if
      outside = edge_shell + step
      saha_iterations = 0
! 1. The first shell outside, with the region's composition, then its own.
      log_density_outside = log_density(outside)
      call radiative_gradient_at(log_temperature(outside), log_pressure(outside), &
           log_radius(outside), log_mass(outside), luminosity_lsun(outside), &
           composition(:,edge_shell), log_teff, log_density_outside, eos_res, &
           saha_iterations, gradient_zone_mix, ierr)
      if (ierr /= 0) return
      if (gradient_zone_mix < eos_res%adiabatic_gradient) return
      inverse_mu_zone_mix = eos_res%ion_mean_weight_inverse + &
           eos_res%electron_mean_weight_inverse
      superadiabatic_fraction = 1.0d0 - eos_res%adiabatic_gradient/gradient_zone_mix
      call radiative_gradient_at(log_temperature(outside), log_pressure(outside), &
           log_radius(outside), log_mass(outside), luminosity_lsun(outside), &
           composition(:,outside), log_teff, log_density_outside, eos_res, &
           saha_iterations, gradient_local, ierr)
      if (ierr /= 0) return
      log_density(outside) = log_density_outside
      gradient_ratio = gradient_zone_mix/gradient_local
      max_distance = superadiabatic_fraction*(luminosity_lsun(edge_shell)* &
           convection_settings%solar_luminosity_cgs/(propagation_speed_divisor*c4pi* &
           exp(ln10*log_pressure(edge_shell))))* &
           (timestep/exp(ln10*(log_radius(edge_shell) + log_radius(edge_shell))))
! 2. Out shell by shell, to the distance limit or a stable shell.
      distance = 0.0d0
      reached_limit = .false.
      radius_prev = exp(ln10*log_radius(outside - step))
      do shell = outside + step, last_shell, step
         radius_curr = exp(ln10*log_radius(shell - step))
         radius_step = radius_curr - radius_prev
         distance = distance + (1.0d0 - inverse_mu_zone_mix/ &
              (eos_res%ion_mean_weight_inverse + eos_res%electron_mean_weight_inverse))* &
              radius_step
         radius_prev = radius_curr
         if (distance > max_distance) then
            reached_limit = .true.
            exit
         end if
         log_density_shell = log_density(shell)
         call radiative_gradient_at(log_temperature(shell), log_pressure(shell), &
              log_radius(outside), log_mass(shell), luminosity_lsun(shell), &
              composition(:,shell), log_teff, log_density_shell, eos_res, &
              saha_iterations, gradient, ierr)
         if (ierr /= 0) return
         log_density(shell) = log_density_shell
         if (gradient_ratio*gradient < eos_res%adiabatic_gradient) exit
      end do
! A search that ran off the model's end reached the limit too.
      if (shell < 1 .or. shell > num_zones) reached_limit = .true.
      new_edge = shell - step
      call log_semiconvective_edge(edge_shell, new_edge, reached_limit, &
           gradient_zone_mix, gradient_local)
      edge_shell = new_edge
end subroutine extend_edge

!----------------------------------------------------------------------
! radiative_gradient_at: del_rad at a shell's log10 T, P, r, m and L/Lsun
! with the composition shell_composition (its X, Z and species): the EOS
! (eos_get, its density seeded from log_density and returned in it), the
! opacity at that density (kap_get) and compute_temperature_gradients,
! without derivatives or rotation factors. eos_res holds the EOS's result
! afterwards (del_ad, the inverse mean molecular weights).
subroutine radiative_gradient_at(log_temperature, log_pressure, log_radius, log_mass, &
     luminosity_lsun, shell_composition, log_teff, log_density, eos_res, &
     saha_iterations, radiative_gradient, ierr)
      real(dp), intent(in) :: log_temperature, log_pressure, log_radius, log_mass, &
           luminosity_lsun, shell_composition(:), log_teff
      real(dp), intent(inout) :: log_density
      type(eos_result), intent(inout) :: eos_res
      integer, intent(inout) :: saha_iterations
      real(dp), intent(out) :: radiative_gradient
      integer, intent(out) :: ierr
      real(dp) :: kap_res(num_kap_results), actual_gradient, dgrad_dt, dgrad_dp, &
           dgrad_dr, convective_velocity
      logical :: is_convective

      eos_res%log10_density = log_density
      call eos_get(log_temperature, log_pressure, shell_composition(i_h1), &
           shell_composition(i_metals), eos_res, .false., .true., saha_iterations, &
           composition_at_zone=shell_composition, ierr=ierr)
      if (ierr /= 0) return
      log_density = eos_res%log10_density
      call kap_get(eos_res%log10_density, log_temperature, shell_composition(i_h1), &
           shell_composition(i_metals), kap_res, eos_res%ion_fraction, ierr=ierr)
      if (ierr /= 0) return
      call compute_temperature_gradients(log_temperature, log_pressure, eos_res, kap_res, &
           log_radius, log_mass, luminosity_lsun, actual_gradient, radiative_gradient, &
           dgrad_dt, dgrad_dp, dgrad_dr, convective_velocity, .false., is_convective, &
           1.0d0, 1.0d0, log_teff, .false., ierr)
end subroutine radiative_gradient_at

! Sources: SCONVEC of the Yale code, after Castellani et al. (1971,
! Ap&SS 10, 340). The opacity call is DBG's (12/95) and the
! radius's MHP's (10/02); the core-only switch is the YREC team's
! (a8e60674, 2025-12-17), which removed the loop over every zone.

end module semiconvection
