!----------------------------------------------------------------------
! overshoot
!----------------------------------------------------------------------
! Step overshoot of the convection zones, and the pressure scale height
! it is measured in (formerly oversh and hsubp, and core's
! pressure_scale_height and envelope_overshoot_depth):
!   pressure_scale_height     H_P = P r^2/(G m rho), cm
!   envelope_overshoot_depth  alpha_env H_P below the surface convection
!                             zone, cm (core's Li/Be burn and rezone)
!   extend_by_overshoot       moves the edges of the mixed regions by
!                             alpha H_P (find_convection_zones)
! The switches and the alphas are convection_settings'. Two densities
! enter H_P today: extend_by_overshoot takes the EOS's for the edge
! shell's P, T and composition (edge_scale_height), and
! envelope_overshoot_depth the model's; making them one changes results
! (the B list). log is log10 throughout.
module overshoot
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use convection_configuration, only: convection_settings
      use composition_layout, only: i_h1, i_metals
      use convection_messages, only: log_mixed_regions
      use eos_lib, only: eos_result, eos_get
      use physical_constants, only: ln10
      use math_lib
      implicit none
      private
      public :: pressure_scale_height, envelope_overshoot_depth, extend_by_overshoot

contains

!----------------------------------------------------------------------
! pressure_scale_height: H_P = P r^2/(G m rho) in cm, from log10 P, r,
! rho and m (cgs), with G from convection_settings.
pure real(dp) function pressure_scale_height(log_pressure, log_radius, log_density, &
     log_mass)
      real(dp), intent(in) :: log_pressure, log_radius, log_density, log_mass

      pressure_scale_height = exp(ln10*(log_pressure + 2.0d0*log_radius - log_density - &
           convection_settings%log10_gravitational_constant - log_mass))
end function pressure_scale_height

!----------------------------------------------------------------------
! envelope_overshoot_depth: how far, in cm, envelope overshoot reaches
! below the base of the surface convection zone, shell base_zone of the
! log10 P, r, rho and m profiles: alpha_env H_P there, with the model's
! density; 0 when envelope overshoot is off.
pure real(dp) function envelope_overshoot_depth(log_pressure, log_radius, log_density, &
     log_mass, base_zone)
      real(dp), intent(in) :: log_pressure(:), log_radius(:), log_density(:), log_mass(:)
      integer, intent(in) :: base_zone

      if (convection_settings%envelope_overshoot) then
         envelope_overshoot_depth = convection_settings%overshoot_alpha_envelope* &
              pressure_scale_height(log_pressure(base_zone), log_radius(base_zone), &
              log_density(base_zone), log_mass(base_zone))
      else
         envelope_overshoot_depth = 0.0d0
      end if
end function envelope_overshoot_depth

!----------------------------------------------------------------------
! extend_by_overshoot (formerly oversh)
!----------------------------------------------------------------------
! Moves the edges of the mixed regions of the zone list mixed_zone_bounds
! (num_mixed_zones regions; mixed_zone_bounds_no_overshoot holds the same
! regions unextended) out by a distance alpha H_P, H_P at the edge shell
! (edge_scale_height), for each kind of region whose switch is on:
!   the convective core (core_overshoot): up by alpha_core H_P, at most
!     beta r with limit_core_overshoot (Woo & Demarque 2001, for small
!     cores; beta is core_overshoot_limit);
!   the surface convection zone (envelope_overshoot): down by
!     alpha_env H_P; the region it adds is the overshoot region below the
!     envelope, overshoot_bottom..overshoot_top (the edge), which takes
!     the adiabatic gradient with adiabatic_overshoot;
!   an intermediate zone (intermediate_overshoot): both ways by
!     alpha_intermediate H_P (its downward region is stored the same way).
! An edge moves to the last shell within the distance (extend_downward,
! extend_upward). overshoot_bottom and overshoot_top start at -1 (none).
! The regions with and without overshoot are logged; find_convection_zones
! merges the ones that now meet. ierr: an EOS failure's code.
subroutine extend_by_overshoot(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, num_zones, mixed_zone_bounds, &
     mixed_zone_bounds_no_overshoot, num_mixed_zones, overshoot_bottom, overshoot_top, ierr)
      real(dp), intent(in) :: composition(:,:), log_density(:), log_pressure(:), &
           log_radius(:), log_mass(:), log_temperature(:)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: mixed_zone_bounds(:,:)
      integer, intent(in) :: mixed_zone_bounds_no_overshoot(:,:)
      integer, intent(in) :: num_mixed_zones
      integer, intent(out) :: overshoot_bottom, overshoot_top
      integer, intent(out) :: ierr
      integer :: k, edge
      real(dp) :: scale_height_down, scale_height_up, distance

      overshoot_bottom = -1
      overshoot_top = -1
      ierr = 0
      do k = 1, num_mixed_zones
         if (mixed_zone_bounds(k,1) == 1) then
! The convective core (find_convection_zones returns before calling
! this for a fully convective model, so it ends below the surface).
            if (.not. convection_settings%core_overshoot) cycle
            edge = mixed_zone_bounds(k,2)
            call edge_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge, scale_height_up, ierr)
            if (ierr /= 0) return
            if (convection_settings%limit_core_overshoot) then
               distance = min(scale_height_up*convection_settings%overshoot_alpha_core, &
                    convection_settings%core_overshoot_limit*exp(ln10*log_radius(edge)))
            else
               distance = scale_height_up*convection_settings%overshoot_alpha_core
            end if
            call extend_upward(log_radius, num_zones, distance, &
                 mixed_zone_bounds_no_overshoot(k,2), mixed_zone_bounds(k,2))
         else if (mixed_zone_bounds(k,2) == num_zones) then
! The surface convection zone.
            if (.not. convection_settings%envelope_overshoot) cycle
            edge = mixed_zone_bounds(k,1)
            call edge_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge, scale_height_down, ierr)
            if (ierr /= 0) return
            distance = scale_height_down*convection_settings%overshoot_alpha_envelope
            call extend_downward(log_radius, distance, mixed_zone_bounds_no_overshoot(k,1), &
                 mixed_zone_bounds(k,1))
            overshoot_top = edge
            overshoot_bottom = mixed_zone_bounds(k,1)
         else
! An intermediate zone, neither at the centre nor at the surface.
            if (.not. convection_settings%intermediate_overshoot) cycle
            call edge_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, mixed_zone_bounds(k,1), scale_height_down, ierr)
            if (ierr /= 0) return
            call edge_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, mixed_zone_bounds(k,2), scale_height_up, ierr)
            if (ierr /= 0) return
            edge = mixed_zone_bounds(k,1)
            distance = scale_height_down*convection_settings%overshoot_alpha_intermediate
            call extend_downward(log_radius, distance, mixed_zone_bounds_no_overshoot(k,1), &
                 mixed_zone_bounds(k,1))
            overshoot_top = edge
            overshoot_bottom = mixed_zone_bounds(k,1)
            distance = scale_height_up*convection_settings%overshoot_alpha_intermediate
            call extend_upward(log_radius, num_zones, distance, &
                 mixed_zone_bounds_no_overshoot(k,2), mixed_zone_bounds(k,2))
         end if
      end do
      call log_mixed_regions(mixed_zone_bounds_no_overshoot, mixed_zone_bounds, &
           num_mixed_zones)
end subroutine extend_by_overshoot

!----------------------------------------------------------------------
! extend_downward: a region's lower edge, shell `edge`, moves down to the
! deepest shell less than `distance` (cm) below it, or to the centre.
! In a rotating model with instability transport the edge moves at least
! one shell below edge_no_overshoot, so that the overshoot region holds
! a radiative shell (MHP 11/91).
pure subroutine extend_downward(log_radius, distance, edge_no_overshoot, edge)
      real(dp), intent(in) :: log_radius(:), distance
      integer, intent(in) :: edge_no_overshoot
      integer, intent(inout) :: edge
      real(dp) :: target_radius
      integer :: k

      target_radius = exp(ln10*log_radius(edge)) - distance
      do k = edge - 1, 1, -1
         if (exp(ln10*log_radius(k)) < target_radius) exit
      end do
      edge = k + 1
      if (convection_settings%rotation_active .and. &
           convection_settings%instability_transport_active .and. &
           edge == edge_no_overshoot) edge = edge - 1
end subroutine extend_downward

!----------------------------------------------------------------------
! extend_upward: a region's upper edge, shell `edge`, moves up to the
! highest shell less than `distance` (cm) above it, or to the surface
! (shell num_zones); in a rotating model with instability transport, at
! least one shell above edge_no_overshoot (as extend_downward).
pure subroutine extend_upward(log_radius, num_zones, distance, edge_no_overshoot, edge)
      real(dp), intent(in) :: log_radius(:), distance
      integer, intent(in) :: num_zones, edge_no_overshoot
      integer, intent(inout) :: edge
      real(dp) :: target_radius
      integer :: k

      target_radius = exp(ln10*log_radius(edge)) + distance
      do k = edge + 1, num_zones
         if (exp(ln10*log_radius(k)) > target_radius) exit
      end do
      edge = k - 1
      if (convection_settings%rotation_active .and. &
           convection_settings%instability_transport_active .and. &
           edge == edge_no_overshoot) edge = edge + 1
end subroutine extend_upward

!----------------------------------------------------------------------
! edge_scale_height (formerly hsubp; PD 1984, MHP 1991): H_P (cm) at
! shell `edge`, with the density the EOS gives for its log10 P and T and
! its composition (seeded from the model's). ierr: the EOS's.
subroutine edge_scale_height(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, edge, scale_height, ierr)
      real(dp), intent(in) :: composition(:,:), log_density(:), log_pressure(:), &
           log_radius(:), log_mass(:), log_temperature(:)
      integer, intent(in) :: edge
      real(dp), intent(out) :: scale_height
      integer, intent(out) :: ierr
      type(eos_result) :: eos_res
      integer :: saha_iterations

      saha_iterations = 0
      eos_res = eos_result()
      eos_res%log10_density = log_density(edge)
      call eos_get(log_temperature(edge), log_pressure(edge), composition(i_h1,edge), &
           composition(i_metals,edge), eos_res, .false., .true., saha_iterations, &
           composition_at_zone=composition(:,edge), ierr=ierr)
      if (ierr /= 0) return
      scale_height = pressure_scale_height(log_pressure(edge), log_radius(edge), &
           eos_res%log10_density, log_mass(edge))
end subroutine edge_scale_height

end module overshoot
