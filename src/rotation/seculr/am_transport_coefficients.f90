!----------------------------------------------------------------------
! am_transport_coefficients
!----------------------------------------------------------------------
! The implicit difference equations of angular momentum transport on
! the equally spaced grid of a transport region (am_transport_grid):
!   - am_diffusion_coeffs: the tridiagonal system of pure diffusion,
!     solved by solve_angular_momentum_transport (secular_transport).
!   - am_advection_diffusion_coeffs, commented out (see its banner):
!     advection plus diffusion (Zahn 1992; Maeder & Zahn 1998), four
!     unknowns per grid point in one band system solved by banded_solver.
module am_transport_coefficients
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: am_diffusion_coeffs

contains

! ---------------------------------------------------------------
! am_diffusion_coeffs (legacy dcoeft)
! The tridiagonal system for omega at the end of a step of angular
! momentum diffusion, on the equally spaced grid of a transport region
! (coordinate chi, spacing grid_spacing, points i = 1..N with
! N = num_eq_points, from the bottom of the region up).
!
! Shell i has moment of inertia I_i (eq_moment_of_inertia), angular
! momentum J_i (eq_angular_momentum) and angular velocity omega_i
! (eq_omega) at the start of the step. Edge i lies between shells i-1
! and i, and D_i = diffusion_coeff(i) is the transport coefficient there
! (eq_am_diffusion_coeff, which am_transport_grid has already multiplied
! by the geometric factor and the Jacobian d chi/dr). Each shell's
! angular momentum changes by the net flux through its two edges,
! differenced fully implicitly (n the start of the step, n+1 the end):
!   omega_i^{n+1} - omega_i^n = timestep/(grid_spacing I_i)
!       [D_{i+1} (omega_{i+1} - omega_i) - D_i (omega_i - omega_{i-1})]^{n+1},
! that is
!   sub_diag(i) omega_{i-1} + diag(i) omega_i + super_diag(i) omega_{i+1}
!       = rhs(i) = omega_i^n.
! No angular momentum crosses the bottom edge of shell 1 or the top edge
! of shell N. With fix_omega_at_surface (disk locking) the top shell's
! omega is held at rhs(N) instead.
!
! The wind: wind_loss_explicit and wind_loss_implicit are the drops in
! surface omega over the step from the torque at its start and at its
! end (compute_solid_body_braking; zero when the surface convection zone
! is not part of the region). Their mean, as angular momentum,
!   surface_wind_loss_term = -(1/2) (wind_loss_explicit
!       + wind_loss_implicit) I_N,
! is taken from the top shell, rhs(N) = omega_N (1 + surface_wind_loss_term
! /J_N), and returned for solve_angular_momentum_transport to start
! dj(N) from.
!
! The commented-out advection-diffusion transport had a branch here
! adding the advective flux A_i (omega_{i-1} + omega_i)/2 through each
! edge, with A = am_advective_coeff and D = am_diffusive_coeff (from
! am_transport_grid) in place of diffusion_coeff:
!   omega_i^{n+1} - omega_i^n = (timestep/I_i) {[D_{i+1} (omega_{i+1}
!       - omega_i) - D_i (omega_i - omega_{i-1})]/grid_spacing
!       + A_{i+1} (omega_i + omega_{i+1})/2
!       - A_i (omega_{i-1} + omega_i)/2}^{n+1}.
!
! composition_diffusion_coeffs (rotational_mixing) sets up the same
! tridiagonal system for composition.
subroutine am_diffusion_coeffs(diffusion_coeff, grid_spacing, timestep, &
     eq_moment_of_inertia, eq_angular_momentum, eq_omega, num_eq_points, &
     wind_loss_explicit, wind_loss_implicit, fix_omega_at_surface, &
     sub_diag, diag, super_diag, rhs, surface_wind_loss_term)
      use star_info_lib, only: json
      implicit none

      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: grid_spacing, timestep
      double precision, intent(in) :: eq_moment_of_inertia(json), &
           eq_angular_momentum(json), eq_omega(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: wind_loss_explicit, wind_loss_implicit
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: surface_wind_loss_term
      logical, intent(in) :: fix_omega_at_surface
      integer :: i
      double precision :: fact0, fact

      fact0 = timestep/grid_spacing
      do i = 1,num_eq_points
         rhs(i) = eq_omega(i)
      end do
      surface_wind_loss_term = -0.5d0*(wind_loss_explicit+ &
           wind_loss_implicit)*eq_moment_of_inertia(num_eq_points)
      rhs(num_eq_points) = rhs(num_eq_points)*(1.0d0+ &
           surface_wind_loss_term/eq_angular_momentum(num_eq_points))
! Shell 1: no flux through its bottom edge.
      fact = fact0/eq_moment_of_inertia(1)
      sub_diag(1) = 0.0d0
      diag(1) = 1.0d0 + fact*diffusion_coeff(2)
      super_diag(1) = -fact*diffusion_coeff(2)
      do i = 2,num_eq_points-1
         fact = fact0/eq_moment_of_inertia(i)
         sub_diag(i) = -fact*diffusion_coeff(i)
         diag(i) = 1.0d0 + fact*(diffusion_coeff(i)+diffusion_coeff(i+1))
         super_diag(i) = -fact*diffusion_coeff(i+1)
      end do
! Top shell: no flux through its top edge, or omega held fixed.
      fact = fact0/eq_moment_of_inertia(num_eq_points)
      if (.not.fix_omega_at_surface) then
         sub_diag(num_eq_points) = -fact*diffusion_coeff(num_eq_points)
         diag(num_eq_points) = 1.0d0 + fact*diffusion_coeff(num_eq_points)
      else
         sub_diag(num_eq_points) = 0.0d0
         diag(num_eq_points) = 1.0d0
      end if
      super_diag(num_eq_points) = 0.0d0
! (commented out: the advection-diffusion branch, see the header)
! Advection plus diffusion (see the header): facta = timestep/(2 I_i)
! multiplies the advective terms. Shell 1: no flux through its bottom
! edge.
!         fact = fact0/eq_moment_of_inertia(1)
!         facta = 0.5d0*timestep/eq_moment_of_inertia(1)
!         sub_diag(1) = 0.0d0
!         diag(1) = 1.0d0 + fact*rot_scr%am_diffusive_coeff(2) - &
!              facta*rot_scr%am_advective_coeff(2)
!         super_diag(1) = -fact*rot_scr%am_diffusive_coeff(2) - &
!              facta*rot_scr%am_advective_coeff(2)
!         do i = 2,num_eq_points-1
!            fact = fact0/eq_moment_of_inertia(i)
!            facta = 0.5d0*timestep/eq_moment_of_inertia(i)
!            sub_diag(i) = -fact*rot_scr%am_diffusive_coeff(i) + &
!                 facta*rot_scr%am_advective_coeff(i)
!            diag(i) = 1.0d0 + fact*(rot_scr%am_diffusive_coeff(i)+ &
!                 rot_scr%am_diffusive_coeff(i+1)) + &
!                 facta*(rot_scr%am_advective_coeff(i)-rot_scr%am_advective_coeff(i+1))
!            super_diag(i) = -fact*rot_scr%am_diffusive_coeff(i+1) - &
!                 facta*rot_scr%am_advective_coeff(i+1)
!         end do
! Top shell: no flux through its top edge, or omega held fixed.
!         fact = fact0/eq_moment_of_inertia(num_eq_points)
!         facta = 0.5d0*timestep/eq_moment_of_inertia(num_eq_points)
!         if (.not.fix_omega_at_surface) then
!            sub_diag(num_eq_points) = -fact*rot_scr%am_diffusive_coeff(num_eq_points) &
!                 + facta*rot_scr%am_advective_coeff(num_eq_points)
!            diag(num_eq_points) = 1.0d0 + fact*rot_scr%am_diffusive_coeff(num_eq_points) &
!                 + facta*rot_scr%am_advective_coeff(num_eq_points)
!         else
!            sub_diag(num_eq_points) = 0.0d0
!            diag(num_eq_points) = 1.0d0
!         end if
!         super_diag(num_eq_points) = 0.0d0
! (end of the commented-out branch)

      return
end subroutine am_diffusion_coeffs

! ---------------------------------------------------------------
! COMMENTED OUT (2026-10-12, the author's decision): the advection-
! diffusion transport (use_diffusion_advection_transport, legacy LDIFAD):
! am_advection_diffusion_coeffs and banded_solver below, and the
! advection branches of am_diffusion_coeffs, am_transport_grid,
! circulation_velocities and rotation_stability_setup, the d theta/dt
! term (rezone, evolve_angular_momentum) and their scratch arrays
! (rotation_scratch_lib) and henyey_coefficients' derivative copies.
! It never worked in the modernized code: most of the coefficient
! arrays it reads were never filled, and apply_secular_transport stops
! a run that sets the flag. The code is kept below for reference. To
! revive it, uncomment those pieces (each says where it was cut) and
! fill the arrays the header below lists as unset; it has not been
! compiled since 2026-10.
! ---------------------------------------------------------------
! am_advection_diffusion_coeffs (legacy dadcoeft; MHP 06/02)
! Angular momentum transport by the meridional circulation and by
! diffusion (Zahn 1992; Maeder & Zahn 1998),
!   rho d(r^2 omega)/dt = (1/(5 r^2)) d/dr(rho r^4 omega U)
!                         + (1/r^2) d/dr(rho nu r^4 d omega/dr),
! over one step, on the equally spaced grid of a transport region
! (coordinate chi, spacing grid_spacing, shells i = 1..N with
! N = num_eq_points, edge i between shells i-1 and i). The circulation
! velocity U depends on omega and its first three derivatives, so the
! equation is written as four first-order equations per shell and
! solved, fully implicitly, as one band system (banded_solver). The
! unknowns of shell i, in rows 4i-3..4i:
!   omega_i;
!   s_i = [S_{i+1} g_i - S_i g_{i-1}]/grid_spacing (legacy d^2 omega/dr^2),
!     S_k = second_deriv_geom_factor_eqgrid(k) dchi_dr_edge(k) fx_k;
!     s_1 = s_N = 0;
!   g_i = (omega_{i+1} - omega_i)/grid_spacing, d omega/d chi at edge i+1;
!   t_i = dchi_dr_edge(i) [T_{i+1} s_{i+1} - T_i s_i]/grid_spacing
!     (legacy d^3 omega/dr^3), T = third_deriv_geom_factor_eqgrid;
!     t_1 = t_2/3.
! g_N and t_N are not unknowns (zero above the region), which leaves
! 4N - 2 equations; that must not exceed band_nmax (not checked). The
! omega equation of shell i is
!   omega_i^{n+1} - omega_i^n = dt/(I_i grid_spacing) (F_{i+1} - F_i)
! (dt the sub-step, I_i = eq_moment_of_inertia(i)), with the flux
! through edge k
!   F_k = A_k (omega_{k-1} + omega_k)/2 + D_k dchi_dr_edge(k) g_{k-1}
!         + B_k (s_{k-1} + s_k)/2 + C_k t_{k-1},
! A = am_advective_coeff, D = am_diffusive_coeff,
! B_k = am_2nd_deriv_coeff(k) dchi_dr_edge(k) fx_k and
! C_k = am_3rd_deriv_coeff(k) dchi_dr_edge(k) fx_k, and no flux through
! the bottom of shell 1 or the top of shell N. fx_k
! (third_order_ratio_factor) is the current omega at edge k over its
! start-of-step value, which updates the higher-order coefficients for
! the change in omega. The derivative rows have zero right-hand sides,
! so they need no start-of-step values.
!
! The iterations, from the inside out:
!   - coefficient iteration (legacy NN): recompute A and D from the
!     current omega and solve, until the largest relative change of
!     omega is below convergence_tolerance;
!   - d theta/dt iteration (legacy NNN): rescale the implicit wind loss
!     to the current surface omega and rerun the coefficient iteration,
!     until a pass changes omega by at most convergence_tolerance (from
!     the second pass on);
!   - num_substeps sub-steps of timestep/num_substeps;
!   - if the last sub-step has not converged, the step is redone with
!     twice the sub-steps; if it still fails after 5 cuts, ierr = 1.
! The two iterations each run at most max_diffusion_iters times. After each
! solve omega is rescaled so that the total angular momentum is
! conserved apart from the wind loss.
!
! Inputs: grid_spacing; timestep, the full step; eq_moment_of_inertia
! and eq_omega at the start of the step; num_eq_points;
! wind_loss_explicit and wind_loss_implicit, the drops in surface omega
! over the step from the torque at its start and at its end (zero when
! the surface convection zone is not part of the region);
! fix_omega_at_surface, hold the top shell's omega (disk locking).
! Outputs: eq_delta_angular_momentum(i) = I_i (omega_i^new - omega_i^n)
! and their sum, sum_delta_angular_momentum; eq_mixing_diffusion_coeff,
! the composition mixing coefficient from the final circulation
! velocity; diffusion_converged; ierr.
!
! Never reached: apply_secular_transport stops the run when
! use_diffusion_advection_transport would take effect, because the
! rot_scr arrays read here (geometric_factor, dchi_dr_edge,
! eq_velocity_coeff*, am_2nd/3rd_deriv_coeff, ...) are never filled.
!subroutine am_advection_diffusion_coeffs(grid_spacing, timestep, eq_moment_of_inertia, eq_omega, &
!     num_eq_points, wind_loss_explicit, wind_loss_implicit, &
!     eq_delta_angular_momentum, eq_mixing_diffusion_coeff, &
!     sum_delta_angular_momentum, fix_omega_at_surface, diffusion_converged, ierr)
!      use rotation_scratch_lib
!      use star_info_lib, only: json
!      use math_lib
!      implicit none
! Length of the convergence history, indexed by the coefficient
! iteration: it must be at least max_diffusion_iters (not checked).
!      integer, parameter :: iter_history_max = 50
!      double precision, parameter :: tiny = 1.0d-30

!      double precision, intent(in) :: grid_spacing
!      double precision, intent(in) :: timestep
!      double precision, intent(in) :: eq_moment_of_inertia(json), &
!           eq_omega(json)
!      integer, intent(in) :: num_eq_points
!      double precision, intent(in) :: wind_loss_explicit
!      double precision, intent(in) :: wind_loss_implicit
!      double precision, intent(out) :: eq_delta_angular_momentum(json)
!      double precision, intent(inout) :: eq_mixing_diffusion_coeff(json)
!      double precision, intent(out) :: sum_delta_angular_momentum
!      logical, intent(in) :: fix_omega_at_surface
!      logical, intent(out) :: diffusion_converged
!      integer, intent(out) :: ierr

!      double precision :: coeff_matrix(band_nmax,10), rhs(band_nmax), &
!           omega_working(json), max_omega_change_history(iter_history_max)
!      integer :: max_omega_change_zone_history(iter_history_max)
!      double precision :: third_order_ratio_factor(json), domega_dr(json), &
!           omega_mid(json)
!      double precision :: omega_prev_medium_iter(json), &
!           omega_substep_start(json)
!      double precision :: omega_curvature(json)
! Counters: step cuts, sub-steps, the iteration levels, rows and columns.
!      integer :: timestep_cut_count, num_substeps, substep_idx, &
!           theta_iter_idx, coeff_iter_idx, num_equations, i, j, ii, k
! substep_timestep = timestep/num_substeps is the step actually
! differenced; wind_loss_implicit_iter is wind_loss_implicit rescaled to
! the current surface omega.
!      double precision :: substep_timestep, wind_loss_implicit_iter, &
!           substep_time_sum, total_angular_momentum_start, substep_frac, &
!           wind_saturation_threshold
!      double precision :: omega_capped, omega_prev_capped, &
!           domega_dr_prev_substep, advective_term1, &
!           theta_term_n, theta_term_p, diffusive_term2, dt_over_dr, &
!           fact_over_ei, facta_half_dt_over_ei_dr, wind_loss_delta, &
!           max_omega_change, &
!           max_omega_change_total, max_omega_change_medium_iter, &
!           omega_change, omega_change_total, omega_change_medium_iter, &
!           total_angular_momentum_updated, &
!           angular_momentum_correction, omega_mid_new, domega_dr_mid_new, &
!           advective_velocity_term, domega_dr_velocity_term, &
!           domega_dr_velocity_term_alt, curvature_velocity_term, &
!           third_deriv_velocity_term, total_velocity, total_velocity_alt, &
!           mixing_diffusion_raw, wind_loss_delta_full_step
!      integer :: max_omega_change_zone, max_omega_change_total_zone

!      ierr = 0

!      timestep_cut_count = 0
!      num_substeps = 1
! d omega/dr at each edge at the start of the step,
! dchi_dr_edge(i) (omega_i - omega_{i-1})/grid_spacing.
!      do i = 2,num_eq_points
!         domega_dr(i) = rot_scr%dchi_dr_edge(i)*(eq_omega(i)-eq_omega(i-1))/ &
!              grid_spacing
!      end do
!      num_equations = 4*num_eq_points-2
! The step-cutting loop: the whole step, redone with twice the sub-steps
! until its last sub-step converges.
!      timestep_cut: do   ! redone after each cut
!      substep_time_sum = 0.0d0
! Start from the start-of-step omega and its total angular momentum,
! sum of I_i omega_i.
!      wind_loss_implicit_iter = wind_loss_implicit
!      total_angular_momentum_start = 0.0d0
!      do i = 1,num_eq_points
!         omega_working(i) = eq_omega(i)
!         omega_prev_medium_iter(i) = eq_omega(i)
!         omega_substep_start(i) = eq_omega(i)
!         total_angular_momentum_start = total_angular_momentum_start + &
!              eq_omega(i)*eq_moment_of_inertia(i)
!      end do
!      do substep_idx = 1,num_substeps
!         substep_timestep = timestep/dfloat(num_substeps)
!         substep_time_sum = substep_time_sum + substep_timestep
!         substep_frac = substep_timestep/timestep
! The d theta/dt iteration (legacy NNN). The d theta/dt terms are
! switched off (below), so a pass differs from the last only through
! the rescaled implicit wind loss.
!      do theta_iter_idx = 1,rotation_settings%max_diffusion_iters
!         diffusion_converged = .false.
! From the second pass on, if there is a wind, rescale the implicit
! (end-of-step) wind loss to the current surface omega with the law of
! compute_solid_body_braking, delta omega ~ min(omega, omega_sat)^(p-1)
! omega with p = wind_law_omega_exponent:
!   wind_loss_implicit_iter = wind_loss_implicit
!       [min(omega_N, omega_sat)/min(omega_N^n, omega_sat)]^(p-1)
!       omega_N/omega_N^n,
! omega_N the current and omega_N^n the start-of-step surface omega.
! omega_sat is the plain wind_saturation_omega, without the Rossby-number
! or field scaling compute_solid_body_braking may apply (the original's
! Rossby-scaled branch here, LROSSBY, could never execute).
!         if (theta_iter_idx.gt.1 .and. wind_loss_implicit_iter.gt.0.0d0) then
!            wind_saturation_threshold = rotation_settings%wind_saturation_omega
!            omega_capped = min(eq_omega(num_eq_points), &
!                 wind_saturation_threshold)
!            omega_prev_capped = min(omega_working(num_eq_points), &
!                 wind_saturation_threshold)
!            wind_loss_implicit_iter = wind_loss_implicit* &
!                 pow(omega_prev_capped/omega_capped, &
!                 rotation_settings%wind_law_omega_exponent-1.0d0)* &
!                 (omega_working(num_eq_points)/eq_omega(num_eq_points))
!         end if
! The coefficient iteration (legacy NN): recompute the omega-dependent
! coefficients from the current iterate omega_working, then solve.
!      do coeff_iter_idx = 1,rotation_settings%max_diffusion_iters
! omega at each edge,
!   omega_mid(i) = (omega_{i-1} + omega_i)/2 - (c_i - c_{i-1})/8,
! with c_i = omega_curvature(i) = omega_{i+1} - 2 omega_i + omega_{i-1}
! (zero at both ends). This is not the symmetric cubic midpoint value,
! which would subtract (c_i + c_{i-1})/16.
!      omega_curvature(1) = 0.0d0
!      do i = 2,num_eq_points-1
!         omega_curvature(i) = (omega_working(i+1)-2.0d0*omega_working(i)+ &
!              omega_working(i-1))
!      end do
!      omega_curvature(num_eq_points) = 0.0d0
!      do i = 2,num_eq_points
!         omega_mid(i) = 0.5d0*(omega_working(i)+omega_working(i-1)) &
!              -0.125d0*(omega_curvature(i)-omega_curvature(i-1))
! Rescale the shear and GSF coefficients as (d omega/dr)^2, by
! (current/start-of-step gradient)^2. The product is stored back, so the
! factor compounds from one iteration to the next.
!         domega_dr_prev_substep = rot_scr%dchi_dr_edge(i)*(omega_working(i)- &
!              omega_working(i-1))/grid_spacing
!         if (abs(domega_dr(i)).gt.tiny) then
!            rot_scr%shear_diffusion_coeff_eqgrid(i) = rot_scr%shear_diffusion_coeff_eqgrid(i)* &
!                 (domega_dr_prev_substep/domega_dr(i))**2
!            rot_scr%gsf_diffusion_coeff_eqgrid(i) = rot_scr%gsf_diffusion_coeff_eqgrid(i)* &
!                 (domega_dr_prev_substep/domega_dr(i))**2
!         end if
! The circulation velocity at the edge, as a part independent of the
! omega gradient and a part that multiplies it (omega = omega_mid,
! f = difad_velocity_scale, ev = eq_velocity_coeff):
!   U_0 = f ev0 omega^2 (ev1a + omega^2 ev1b)   (advective_term1),
!   U_1 = f ev0 omega^2 (ev2a + ev2b)           (diffusive_term2).
! Then, with the 1/5 of Zahn's advective term,
!   A = am_advective_coeff = geometric_factor U_0/5,
!   D = am_diffusive_coeff = geometric_factor (U_1/5 + D_shear + D_GSF),
! D floored at zero. These overwrite the A and D that am_transport_grid
! set, its vertical mixing term in D included.
!         advective_term1 = rotation_settings%difad_velocity_scale*rot_scr%eq_velocity_coeff0(i)* &
!              omega_mid(i)**2*(rot_scr%eq_velocity_coeff1a(i)+ &
!              omega_mid(i)**2*rot_scr%eq_velocity_coeff1b(i))
! The d theta/dt contributions (theta_term_n to U_1, theta_term_p to
! U_0; legacy VTHN and VTHP) are switched off: both are zero.
!         theta_term_n = 0.0d0
!         theta_term_p = 0.0d0
!         diffusive_term2 = rotation_settings%difad_velocity_scale*rot_scr%eq_velocity_coeff0(i)* &
!              omega_mid(i)**2*(rot_scr%eq_velocity_coeff2a(i)+ &
!              rot_scr%eq_velocity_coeff2b(i))+theta_term_n
!         advective_term1 = advective_term1 + theta_term_p
!         rot_scr%am_advective_coeff(i) = 0.2d0*advective_term1*rot_scr%geometric_factor(i)
!         rot_scr%am_diffusive_coeff(i) = (0.2d0*diffusive_term2+ &
!              rot_scr%shear_diffusion_coeff_eqgrid(i)+rot_scr%gsf_diffusion_coeff_eqgrid(i))* &
!              rot_scr%geometric_factor(i)
!         rot_scr%am_diffusive_coeff(i) = max(0.0d0,rot_scr%am_diffusive_coeff(i))
!      end do
! fx_i = third_order_ratio_factor(i): omega_mid over the start-of-step
! edge omega (the plain mean), for the higher-order terms.
!      do i = 2,num_eq_points
!         third_order_ratio_factor(i) = 2.0d0*omega_mid(i)/(eq_omega(i)+ &
!              eq_omega(i-1))
!      end do
! The right-hand side: omega at the start of the sub-step, and zero for
! the derivative rows.
!      do i = 1,num_eq_points
!         ii = 1 + 4*(i-1)
!         rhs(ii) = omega_substep_start(i)
!         rhs(ii+1) = 0.d0
!         rhs(ii+2) = 0.d0
!         rhs(ii+3) = 0.d0
!      end do
! The wind: over the sub-step the top shell loses
!   wind_loss_delta = -(1/2) (wind_loss_explicit + wind_loss_implicit_iter)
!       I_N substep_frac
! of angular momentum, so its right-hand side becomes
! omega_N + wind_loss_delta/I_N.
!      wind_loss_delta = -0.5d0*(wind_loss_explicit+wind_loss_implicit_iter)* &
!           eq_moment_of_inertia(num_eq_points)*substep_frac
!      rhs(4*num_eq_points-3) = rhs(4*num_eq_points-3)* &
!           (1.0d0+wind_loss_delta/eq_moment_of_inertia(num_eq_points)/ &
!           omega_substep_start(num_eq_points))
! dt/d chi, common to the flux terms.
!      dt_over_dr = substep_timestep/grid_spacing
! Shell 1, rows 1-4: no flux through its bottom edge.
! fact_over_ei = dt/(I d chi) multiplies the D and C terms,
! facta_half_dt_over_ei_dr = dt/(2 I d chi) the edge means of the A and
! B terms.
!      fact_over_ei = dt_over_dr/eq_moment_of_inertia(1)
!      facta_half_dt_over_ei_dr = 0.5d0*substep_timestep/eq_moment_of_inertia(1)/ &
!           grid_spacing
! Clear the band matrix.
!      do i = 1,num_equations
!         do j = 1,10
!            coeff_matrix(i,j) = 0.0d0
!         end do
!      end do
! The omega equation: the advective term at the top edge (edge 2).
!      coeff_matrix(1,5) = 1.0d0  - facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(2)
!      coeff_matrix(1,9) = -facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(2)
! The diffusive term, D_2 dchi_dr_edge(2) g_1.
!      coeff_matrix(1,7) = -fact_over_ei*rot_scr%dchi_dr_edge(2)*rot_scr%am_diffusive_coeff(2)
! The second-derivative term, B_2 s_2/2 (s_1 = 0).
!      coeff_matrix(1,10) = -facta_half_dt_over_ei_dr*rot_scr%dchi_dr_edge(2)* &
!           rot_scr%am_2nd_deriv_coeff(2)*third_order_ratio_factor(2)
! The third-derivative term, C_2 t_1. Its sign is opposite to that of
! the same term (column 8) in the interior rows.
!      coeff_matrix(1,8) = fact_over_ei*rot_scr%dchi_dr_edge(2)*rot_scr%am_3rd_deriv_coeff(2)* &
!           third_order_ratio_factor(2)
! s_1 = 0 at the lower boundary.
!      coeff_matrix(2,5) = 1.0d0
! g_1 = (omega_2 - omega_1)/grid_spacing.
!      coeff_matrix(3,3) = 1.0d0/grid_spacing
!      coeff_matrix(3,5) = 1.0d0
!      coeff_matrix(3,7) = -1.0d0/grid_spacing
! The boundary condition on the highest derivative: t_1 = t_2/3.
!      coeff_matrix(4,5) = 1.0d0
!      coeff_matrix(4,9) = -1.0d0/3.0d0
! Shells ii = 2..N-1. Row i is omega_ii's equation; in it, columns 1-4
! multiply omega, s, g, t of shell ii-1, column 5 omega_ii, columns 6-8
! s, g, t of shell ii, and columns 9-10 omega and s of shell ii+1.
!      do ii = 2,num_eq_points-1
!         fact_over_ei = dt_over_dr/eq_moment_of_inertia(ii)
!         facta_half_dt_over_ei_dr = 0.5d0*substep_timestep/ &
!              eq_moment_of_inertia(ii)/grid_spacing
!         i = 1+(ii-1)*4
! The omega equation: the advective terms at edges ii and ii+1.
!         coeff_matrix(i,1) = facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(ii)
!         coeff_matrix(i,5) = 1.0d0 + facta_half_dt_over_ei_dr* &
!              (rot_scr%am_advective_coeff(ii)-rot_scr%am_advective_coeff(ii+1))
!         coeff_matrix(i,9) = -facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(ii+1)
! The diffusive terms, D dchi_dr_edge g.
!         coeff_matrix(i,3) = fact_over_ei*rot_scr%dchi_dr_edge(ii)* &
!              rot_scr%am_diffusive_coeff(ii)
!         coeff_matrix(i,7) = -fact_over_ei*rot_scr%dchi_dr_edge(ii+1)* &
!              rot_scr%am_diffusive_coeff(ii+1)
! The second-derivative terms, B times the edge means of s.
!         coeff_matrix(i,2) = facta_half_dt_over_ei_dr*rot_scr%dchi_dr_edge(ii)* &
!              rot_scr%am_2nd_deriv_coeff(ii)*third_order_ratio_factor(ii)
!         coeff_matrix(i,6) = facta_half_dt_over_ei_dr* &
!              (rot_scr%am_2nd_deriv_coeff(ii)*rot_scr%dchi_dr_edge(ii)* &
!              third_order_ratio_factor(ii) &
!              - rot_scr%am_2nd_deriv_coeff(ii+1)*rot_scr%dchi_dr_edge(ii+1)* &
!              third_order_ratio_factor(ii+1))
!         coeff_matrix(i,10) = -facta_half_dt_over_ei_dr* &
!              rot_scr%am_2nd_deriv_coeff(ii+1)*rot_scr%dchi_dr_edge(ii+1)* &
!              third_order_ratio_factor(ii+1)
! The third-derivative terms, C t.
!         coeff_matrix(i,4) = fact_over_ei*rot_scr%am_3rd_deriv_coeff(ii)* &
!              rot_scr%dchi_dr_edge(ii)*third_order_ratio_factor(ii)
!         coeff_matrix(i,8) = -fact_over_ei*rot_scr%am_3rd_deriv_coeff(ii+1)* &
!              rot_scr%dchi_dr_edge(ii+1)*third_order_ratio_factor(ii+1)
! s_ii = [S_{ii+1} g_ii - S_ii g_{ii-1}]/grid_spacing.
!         coeff_matrix(i+1,2) = rot_scr%second_deriv_geom_factor_eqgrid(ii)* &
!              rot_scr%dchi_dr_edge(ii)/grid_spacing*third_order_ratio_factor(ii)
!         coeff_matrix(i+1,5) = 1.0d0
!         coeff_matrix(i+1,6) = -rot_scr%second_deriv_geom_factor_eqgrid(ii+1)* &
!              rot_scr%dchi_dr_edge(ii+1)/grid_spacing*third_order_ratio_factor(ii+1)
! g_ii = (omega_{ii+1} - omega_ii)/grid_spacing.
!         coeff_matrix(i+2,3) = 1.0d0/grid_spacing
!         coeff_matrix(i+2,5) = 1.0d0
!         coeff_matrix(i+2,7) = -1.0d0/grid_spacing
! t_ii = dchi_dr_edge(ii) [T_{ii+1} s_{ii+1} - T_ii s_ii]/grid_spacing
! (dchi_dr_edge of the lower edge for both terms).
!         coeff_matrix(i+3,3) = rot_scr%dchi_dr_edge(ii)* &
!              rot_scr%third_deriv_geom_factor_eqgrid(ii)/grid_spacing
!         coeff_matrix(i+3,5) = 1.0d0
!         coeff_matrix(i+3,7) = -rot_scr%dchi_dr_edge(ii)* &
!              rot_scr%third_deriv_geom_factor_eqgrid(ii+1)/grid_spacing
!      end do
! The top shell N, rows 4N-3 and 4N-2: no flux through its top edge.
!      fact_over_ei = dt_over_dr/eq_moment_of_inertia(num_eq_points)
!      facta_half_dt_over_ei_dr = 0.5d0*substep_timestep/ &
!           eq_moment_of_inertia(num_eq_points)/grid_spacing
!      i = 4*num_eq_points - 3
! s_N and s_1 are zero: drop them from the omega equations of shells
! N-1 and 2.
!      coeff_matrix(i-4,10) = 0.0d0
!      coeff_matrix(5,2) = 0.0d0
!      if (.not.fix_omega_at_surface) then
! The omega equation: the advective term at edge N.
!         coeff_matrix(i,1) = facta_half_dt_over_ei_dr* &
!              rot_scr%am_advective_coeff(num_eq_points)
!         coeff_matrix(i,5) = 1.0d0 + facta_half_dt_over_ei_dr* &
!              rot_scr%am_advective_coeff(num_eq_points)
! The diffusive term.
!         coeff_matrix(i,3) = fact_over_ei*rot_scr%dchi_dr_edge(num_eq_points)* &
!              rot_scr%am_diffusive_coeff(num_eq_points)
! The second-derivative term (s_N = 0).
!         coeff_matrix(i,2) = facta_half_dt_over_ei_dr* &
!              rot_scr%dchi_dr_edge(num_eq_points)*rot_scr%am_2nd_deriv_coeff(num_eq_points)* &
!              third_order_ratio_factor(num_eq_points)
! The third-derivative term, without the fx_N the interior rows apply.
!         coeff_matrix(i,4) = fact_over_ei*rot_scr%dchi_dr_edge(num_eq_points)* &
!              rot_scr%am_3rd_deriv_coeff(num_eq_points)
!      else
! omega_N held at its right-hand side.
!         coeff_matrix(i,5) = 1.0d0
!      end if
! s_N = 0: d omega/dr is taken as zero above the region.
!      coeff_matrix(i+1,5) = 1.0d0
! g_N and t_N are not unknowns: they are zero above the region.
! Solve the num_equations = 4N - 2 equations; the solution overwrites rhs.
!      call banded_solver(coeff_matrix,num_equations,rhs, ierr)
!      if (ierr /= 0) return
! The change of angular momentum of each shell since the start of the
! step, and the largest relative changes of omega: since the previous
! coefficient iteration (max_omega_change, legacy DWMAX), since the start
! of the step (max_omega_change_total, DWW) and since the previous
! d theta/dt pass (max_omega_change_medium_iter, DWI).
!      sum_delta_angular_momentum = 0.0d0
!      max_omega_change = (rhs(1)-omega_working(1))/omega_working(1)
!      max_omega_change_zone = 1
!      max_omega_change_total = (rhs(1)-eq_omega(1))/eq_omega(1)
!      max_omega_change_total_zone = 1
!      max_omega_change_medium_iter = (rhs(1)-omega_prev_medium_iter(1))/ &
!           omega_prev_medium_iter(1)
!      do i = 1,num_eq_points
!         ii = 1+4*(i-1)
!         eq_delta_angular_momentum(i) = (rhs(ii)-eq_omega(i))* &
!              eq_moment_of_inertia(i)
!         sum_delta_angular_momentum = sum_delta_angular_momentum + &
!              eq_delta_angular_momentum(i)
!         omega_change = (rhs(ii)-omega_working(i))/omega_working(i)
!         omega_change_total = (rhs(ii)-eq_omega(i))/eq_omega(i)
!         omega_change_medium_iter = (rhs(ii)-omega_prev_medium_iter(i))/ &
!              omega_prev_medium_iter(i)
!         if (abs(omega_change).gt.abs(max_omega_change)) then
!            max_omega_change_zone = i
!            max_omega_change = omega_change
!         end if
!         if (abs(omega_change_total).gt.abs(max_omega_change_total)) then
!            max_omega_change_total_zone = i
!            max_omega_change_total = omega_change_total
!         end if
!         if (abs(omega_change_medium_iter).gt.abs(max_omega_change_medium_iter)) &
!              then
!            max_omega_change_medium_iter = omega_change_medium_iter
!         end if
!      end do
! A relative change above 10% is taken as spurious and not applied: the
! coefficient iteration ends unconverged. Only that loop is left; the
! d theta/dt iteration and any later sub-steps carry on.
!      if (abs(max_omega_change).gt.0.1d0) then
!         write(*,*) substep_idx,theta_iter_idx,coeff_iter_idx
!         write(*,*) 'CORRECTIONS TOO LARGE ',max_omega_change, &
!              max_omega_change_zone
!         diffusion_converged = .false.
!         exit
!      end if
! Take the new omega in full (no damping).
!      total_angular_momentum_updated = 0.0d0
!      do i = 1,num_eq_points
!         ii = 1+4*(i-1)
!         omega_working(i) = omega_working(i)+ &
!              (rhs(ii)-omega_working(i))
!         total_angular_momentum_updated = total_angular_momentum_updated+ &
!              omega_working(i)*eq_moment_of_inertia(i)
!      end do
! Conserve angular momentum: scale omega so that the total, sum of
! I_i omega_i, is the start-of-sub-step total plus the wind loss.
!      angular_momentum_correction = (total_angular_momentum_start+ &
!           wind_loss_delta)/total_angular_momentum_updated
!      write(*,*) total_angular_momentum_start,wind_loss_delta, &
!           total_angular_momentum_updated,angular_momentum_correction
!      sum_delta_angular_momentum = 0.0d0
!      do i = 1,num_eq_points
!         omega_working(i) = omega_working(i)*angular_momentum_correction
!         eq_delta_angular_momentum(i) = (omega_working(i)-eq_omega(i))* &
!              eq_moment_of_inertia(i)
!         sum_delta_angular_momentum = sum_delta_angular_momentum+ &
!              eq_delta_angular_momentum(i)
!      end do
!      write(*,*) sum_delta_angular_momentum,wind_loss_delta
!      max_omega_change_history(coeff_iter_idx) = max_omega_change
!      max_omega_change_zone_history(coeff_iter_idx) = max_omega_change_zone
! Converged when the largest relative change is below
! convergence_tolerance.
!      if (abs(max_omega_change).lt.rotation_settings%convergence_tolerance) then
!         diffusion_converged = .true.
!         exit
!      end if
!      end do
! The convergence history of the coefficient iteration. When the loop
! runs out, coeff_iter_idx is max_diffusion_iters + 1 and the last entry
! printed was never set.
!      write(*,*) substep_idx,theta_iter_idx,coeff_iter_idx
!      write(*,100) max_omega_change_total,max_omega_change_total_zone, &
!           (max_omega_change_history(j),max_omega_change_zone_history(j), &
!           j=1,coeff_iter_idx)
! 100  format(' MAX D(OMEGA)/OMEGA',1pe12.3,' AT PT.',i5, &
!         ' BY ITERATION'/5(1x,e11.3,i4))
! The d theta/dt iteration ends when a pass, from the second on, changed
! omega by at most convergence_tolerance.
!      do k = 1,num_eq_points
!         omega_prev_medium_iter(k) = omega_working(k)
!      end do
!      if (abs(max_omega_change_medium_iter).le.rotation_settings%convergence_tolerance .and. &
!           theta_iter_idx.ge.2) then
!         exit
!      else
!         diffusion_converged = .false.
!      end if
!      end do
! The next sub-step starts from this one's omega and angular momentum.
!      do k = 1,num_eq_points
!         omega_substep_start(k) = omega_working(k)
!      end do
!      total_angular_momentum_start = total_angular_momentum_start+ &
!           wind_loss_delta
!      end do
! Only the last sub-step's convergence is tested: if it failed, redo the
! step with twice the sub-steps; stop if it still fails after 5 cuts.
!      if (.not.diffusion_converged) then
!         timestep_cut_count = timestep_cut_count + 1
!         num_substeps = 2*num_substeps
!         if (timestep_cut_count.gt.5) then
!            write(*,*) 'TIMESTEP CUT MORE THAN 5 TIMES - RUN STOPPED'
            ! A nonzero ierr stops the run.
!            ierr = 1
!            return
!         else
!            write(*,*) 'TIMESTEP CUT #',timestep_cut_count,' IN DADCOEFT'
!         end if
!         cycle timestep_cut
!      end if
!      exit timestep_cut
!      end do timestep_cut
! Print shell 1's solution: omega_1 (and its start-of-step value), s_1,
! g_1, t_1.
!      write(*,914) rhs(1),eq_omega(1),rhs(2),rhs(3),rhs(4)
! 914  format(1p5e12.3)
! The effective circulation velocity at each edge i from the last
! solution (rhs, before the angular momentum rescaling), term by term:
! each term of the flux F_i (header, without fx) over geometric_factor(i)
! times the edge omega (the mean of the two shells'), so that the terms
! add up to U/5. domega_dr_mid_new is the edge mean of s, not d omega/dr.
! total_velocity_alt keeps only the circulation's own gradient term,
! without the shear and GSF diffusion. Edges within 10 of either end
! are printed. The mixing coefficient is then
!   eq_mixing_diffusion_coeff = mixing_geometric_factor
!       [5 |total_velocity_alt| r mixing_velocity_scale es_mixing_scale
!        + mixing_velocity_scale (secular_shear_mixing_scale D_shear
!          + gsf_mixing_scale D_GSF)],
! r = equatorial_radius.
!      do i = 2,num_eq_points
!         ii = 1 + 4*(i-1)
!         omega_mid_new = 0.5d0*(rhs(ii)+rhs(ii-4))
!         domega_dr_mid_new = 0.5d0*(rhs(ii+1)+rhs(ii-3))
!         advective_velocity_term = rot_scr%am_advective_coeff(i)/rot_scr%geometric_factor(i)
!         domega_dr_velocity_term = rot_scr%am_diffusive_coeff(i)/rot_scr%geometric_factor(i)* &
!              rhs(ii-2)*rot_scr%dchi_dr_edge(i)/omega_mid_new
!         domega_dr_velocity_term_alt = 0.2d0*rotation_settings%difad_velocity_scale* &
!              rot_scr%eq_velocity_coeff0(i)*omega_mid_new*rot_scr%dchi_dr_edge(i)* &
!              (rot_scr%eq_velocity_coeff2a(i)+rot_scr%eq_velocity_coeff2b(i))*rhs(ii-2)
!         if (rot_scr%am_diffusive_coeff(i).le.0.0d0) domega_dr_velocity_term_alt = 0.0d0
!         curvature_velocity_term = rot_scr%am_2nd_deriv_coeff(i)/rot_scr%geometric_factor(i)* &
!              rot_scr%dchi_dr_edge(i)*domega_dr_mid_new/omega_mid_new
!         third_deriv_velocity_term = rot_scr%am_3rd_deriv_coeff(i)/rot_scr%geometric_factor(i)* &
!              rhs(ii-1)*rot_scr%dchi_dr_edge(i)/omega_mid_new
!         total_velocity = advective_velocity_term+domega_dr_velocity_term+ &
!              curvature_velocity_term+third_deriv_velocity_term
!         total_velocity_alt = advective_velocity_term+ &
!              domega_dr_velocity_term_alt+curvature_velocity_term+ &
!              third_deriv_velocity_term
!         if (i.le.10 .or. (num_eq_points-i).le.10) then
!         write(*,911) rhs(ii),eq_omega(i),advective_velocity_term, &
!              domega_dr_velocity_term,curvature_velocity_term, &
!              third_deriv_velocity_term,total_velocity, &
!              rhs(ii+1),rhs(ii+2),rhs(ii+3)
!         end if
!         rot_scr%mixing_velocity_estimate(i) = 5.0d0*abs(rotation_settings%mixing_velocity_scale* &
!              rotation_settings%es_mixing_scale*total_velocity_alt*rot_scr%equatorial_radius(i))
!         mixing_diffusion_raw = (rot_scr%mixing_velocity_estimate(i)+ &
!              rotation_settings%mixing_velocity_scale*rotation_settings%secular_shear_mixing_scale* &
!              rot_scr%shear_diffusion_coeff_eqgrid(i)+rotation_settings%mixing_velocity_scale* &
!              rotation_settings%gsf_mixing_scale*rot_scr%gsf_diffusion_coeff_eqgrid(i))
!         eq_mixing_diffusion_coeff(i) = mixing_diffusion_raw* &
!              rot_scr%mixing_geometric_factor(i)
! 911  format(1p10e12.3)
!      end do
! Print the total change of angular momentum and the wind loss over the
! full step.
!      wind_loss_delta_full_step = wind_loss_delta/substep_frac
!      write(*,*) sum_delta_angular_momentum,wind_loss_delta_full_step
!      return
!end subroutine am_advection_diffusion_coeffs

! ---------------------------------------------------------------
! banded_solver (legacy bandw)
! Solves the nm equations coeff_matrix x = rhs assembled by
! am_advection_diffusion_coeffs, by Gaussian elimination without
! pivoting. The matrix is banded, with 4 sub-diagonals and 5
! super-diagonals, and stored by rows: coeff_matrix(i, c) multiplies
! unknown i + c - 5, so column 5 is the diagonal. The solution
! overwrites rhs and coeff_matrix is overwritten. A pivot smaller than
! 1e-20 in magnitude sets ierr = 1.
!subroutine banded_solver(coeff_matrix, nm, rhs, ierr)
!      use rotation_scratch_lib, only: band_nmax
!      implicit none

!      integer, intent(in) :: nm
!      double precision, intent(inout) :: coeff_matrix(band_nmax,10), rhs(band_nmax)
!      integer, intent(out) :: ierr

!      double precision, parameter :: tiny = 1.d-20

!      integer :: i, j, k, imj, imj5, ii, jj

!      ierr = 0
! Forward elimination, row by row from the top: divide row i by its
! pivot and eliminate unknown i from the four rows below it (row j holds
! it in column 5 + i - j). Rows 1..nm-5 have all five super-diagonals.
!      do i = 1,nm-5
!         if (abs(coeff_matrix(i,5)).lt.tiny) then
!            write(*,911) i,(coeff_matrix(i,j),j=1,10)
            ! A nonzero ierr stops the run.
!            ierr = 1
!            return
!         end if
! 911     format(i5,1p10e12.3)
!         rhs(i) = rhs(i)/coeff_matrix(i,5)
!         do j = 6,10
!            coeff_matrix(i,j) = coeff_matrix(i,j)/coeff_matrix(i,5)
!         end do
!         do j = i+1,i+4
!            imj = i - j
!            imj5 = 5 + imj
!            rhs(j) = rhs(j) - coeff_matrix(j,imj5)*rhs(i)
!            do k = imj5+1,10+imj
!               coeff_matrix(j,k) = coeff_matrix(j,k)- &
!                    coeff_matrix(j,imj5)*coeff_matrix(i,k+j-i)
!            end do
!         end do
!      end do
! The last rows, where the end of the matrix cuts the band: row i has
! only ii = nm - i unknowns to its right.
!      do i = nm-4,nm-1
!         ii = nm-i
!         if (abs(coeff_matrix(i,5)).lt.tiny) then
!            write(*,911) i,(coeff_matrix(i,j),j=1,10)
            ! A nonzero ierr stops the run.
!            ierr = 1
!            return
!         end if
!         rhs(i) = rhs(i)/coeff_matrix(i,5)
!         do j = 6,5+ii
!            coeff_matrix(i,j) = coeff_matrix(i,j)/coeff_matrix(i,5)
!         end do
!         do j = i+1,nm
!            imj = i - j
!            imj5 = 5 + imj
!            rhs(j) = rhs(j) - coeff_matrix(j,imj5)*rhs(i)
!            do k = imj5+1,imj5+ii
!               coeff_matrix(j,k) = coeff_matrix(j,k)- &
!                    coeff_matrix(j,imj5)*coeff_matrix(i,k+j-i)
!            end do
!         end do
!      end do
!      if (abs(coeff_matrix(nm,5)).lt.tiny) then
!         write(*,*) 'banded_solver: singular final pivot'
!         ierr = 1
!         return
!      end if
!      rhs(nm) = rhs(nm)/coeff_matrix(nm,5)

! Back substitution, from the bottom: each solved unknown i is removed
! from the (up to) five rows above it.
!      do i = nm,6,-1
!         do j = -1,-5,-1
!            rhs(i+j) = rhs(i+j) - rhs(i)*coeff_matrix(i+j,5-j)
!         end do
!      end do
!      do i = 5,2,-1
!         jj = -(i-1)
!         do j = -1,jj,-1
!            rhs(i+j) = rhs(i+j) - rhs(i)*coeff_matrix(i+j,5-j)
!         end do
!      end do
!      return
!end subroutine banded_solver

end module am_transport_coefficients
