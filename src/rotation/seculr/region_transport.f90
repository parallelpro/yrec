!----------------------------------------------------------------------
! region_transport
!----------------------------------------------------------------------
! The angular momentum diffusion across one transport region over a
! rotation sub-step, on the region's equally spaced grid (chi_grid, from
! am_transport_grid): the implicit diffusion equation for omega
! (am_diffusion_coeffs), its solve (solve_angular_momentum_transport),
! and the map of the change back to the model's j/m
! (equal_grid_to_model). A star-blind kernel (tools/check_boundaries.py
! STAR_BLIND_FILES), checked by test_rotation_kernels.
module region_transport
      use composition_layout, only: json
      use rotation_scratch_lib, only: transport_grid_state
      use am_transport_coefficients, only: am_diffusion_coeffs, &
           solve_angular_momentum_transport
      use am_grid, only: equal_grid_to_model
      implicit none
      private
      public :: diffuse_region_angular_momentum

contains

! ---------------------------------------------------------------
! Over sub_timestep, on the n = chi_grid%num_eq_points points of the
! region zone_begin..zone_end, omega at the end of the sub-step (primed)
! solves
!   I_k (omega'_k - omega_k) = (dt/dchi) [D_{k+1} (omega'_{k+1} - omega'_k)
!                                        - D_k (omega'_k - omega'_{k-1})],
! with D_k = chi_grid%eq_am_diffusion_coeff(k) at the edge between
! points k-1 and k and no flux through the region's ends. The outermost
! point also loses I_n (wind_loss_explicit + wind_loss_implicit)/2, or,
! with fix_omega_at_surface, keeps its omega (outer_point_wind_terms in
! magnetic_braking sets these). Each point's change of angular momentum,
! I (omega' - omega), is interpolated back as a fractional change of j/m
! to the region's zones (and an adjoining convection zone), rescaled so
! that the total change is the grid's (equal_grid_to_model), into
! specific_angular_momentum. A zero pivot returns the solve's ierr
! (numerics_singular), a spline failure osplin's, each after a report.
subroutine diffuse_region_angular_momentum(chi_grid, sub_timestep, &
     wind_loss_explicit, wind_loss_implicit, fix_omega_at_surface, &
     shell_mass, zone_begin, zone_end, am_transport_convective_flag, &
     num_zones, specific_angular_momentum, ierr)
      type(transport_grid_state), intent(in) :: chi_grid
      double precision, intent(in) :: sub_timestep
      double precision, intent(in) :: wind_loss_explicit, wind_loss_implicit
      logical, intent(in) :: fix_omega_at_surface
      double precision, intent(in) :: shell_mass(json)
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: specific_angular_momentum(json)
      integer, intent(out) :: ierr
! The tridiagonal system and its solution (omega', which only the
! changes below use); the wind's loss from the outermost point, from
! the coefficients to the solve; each point's change of angular
! momentum and their sum.
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json), surface_wind_loss_term
      double precision :: delta_angular_momentum(json), sum_delta_angular_momentum

      call am_diffusion_coeffs(chi_grid%eq_am_diffusion_coeff, chi_grid%eq_chi_spacing, &
           sub_timestep, chi_grid%eq_moment_of_inertia, chi_grid%eq_angular_momentum, &
           chi_grid%eq_omega, chi_grid%num_eq_points, &
           wind_loss_explicit, wind_loss_implicit, fix_omega_at_surface, &
           sub_diag, diag, super_diag, rhs, surface_wind_loss_term)
      call solve_angular_momentum_transport(chi_grid%num_eq_points, chi_grid%eq_moment_of_inertia, &
           delta_angular_momentum, sum_delta_angular_momentum, &
           sub_diag, diag, super_diag, rhs, solution, &
           surface_wind_loss_term, ierr)
      if (ierr /= 0) return
      call equal_grid_to_model(chi_grid, delta_angular_momentum, chi_grid%eq_angular_momentum, &
           shell_mass, zone_begin, zone_end, am_transport_convective_flag, &
           num_zones, sum_delta_angular_momentum, specific_angular_momentum, ierr)
end subroutine diffuse_region_angular_momentum

end module region_transport
