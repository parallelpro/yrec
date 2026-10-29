!----------------------------------------------------------------------
! am_transport_coefficients
!----------------------------------------------------------------------
! The implicit difference equations of angular momentum transport on
! the equally spaced grid of a transport region (am_transport_grid):
!   - am_diffusion_coeffs: the tridiagonal system of pure diffusion,
!     solved by solve_angular_momentum_transport (secular_transport).
!   - am_advection_diffusion_coeffs, retired (commented out 2026-10-12;
!     rotation/retired/advection_diffusion.f90.txt):
!     advection plus diffusion (Zahn 1992; Maeder & Zahn 1998), four
!     unknowns per grid point in one band system solved by banded_solver.
module am_transport_coefficients
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: am_diffusion_coeffs, solve_angular_momentum_transport

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
      use composition_layout, only: json
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
! The advection-diffusion transport (am_advection_diffusion_coeffs and its
! band solver, commented out 2026-10-12 by the author's decision) is in
! rotation/retired/advection_diffusion.f90.txt.

!----------------------------------------------------------------------
! solve_angular_momentum_transport (legacy tridia)
!----------------------------------------------------------------------
! Solves the tridiagonal system of the angular momentum diffusion
! equation on n grid points by the Thomas algorithm. sub_diag, diag,
! super_diag and rhs come from am_diffusion_coeffs: rhs is omega at the
! start of the sub-step, the outermost point's scaled by
! (1 + dj_n_seed/J_n) for the wind. Returns omega at the end in
! solution, and the change of angular momentum of each point and its
! sum,
!   dj(k) = (solution(k) - rhs(k)) ei(k),            k < n,
!   dj(n) = dj_n_seed + (solution(n) - rhs(n)) ei(n),
!   sumdj = sum of dj(k),
! with ei the points' moments of inertia and dj_n_seed the wind's loss
! from the outermost point, -I_n (wind_loss_explicit +
! wind_loss_implicit)/2; dj(n) is thus I_n (omega_end - omega_start)
! when J_n = I_n omega_start. Before the back substitution solution(n)
! is rewritten as rhs(n) (1 + (solution(n) - rhs(n))/rhs(n)), the same
! value up to round-off. A zero pivot returns ierr = 1.
subroutine solve_angular_momentum_transport(n, ei, dj, sumdj, sub_diag, diag, super_diag, rhs, &
     solution, dj_n_seed, ierr)
      use composition_layout, only: json

      integer, intent(in) :: n
      double precision, intent(in) :: ei(json)
      double precision, intent(out) :: dj(json)
      double precision, intent(out) :: sumdj
! The matrix's three diagonals and the right-hand side, from the caller;
! the solution; dj(n)'s starting value.
      double precision, intent(in) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: solution(json)
      double precision, intent(in) :: dj_n_seed
      integer, intent(out) :: ierr

! The elimination's factors, gamma_elim(j) = super_diag(j-1)/bet.
      double precision :: gamma_elim(json)

      double precision :: rhs_orig(json)
      integer :: i, j
      double precision :: bet, fj

      ierr = 0
      dj(n) = dj_n_seed
      do i = 1, n
         rhs_orig(i) = rhs(i)
      end do
      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2, n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) then
            write(*,*) 'solve_angular_momentum_transport: singular matrix'
            ierr = 1
            return
         end if
         solution(j) = (rhs(j) - sub_diag(j)*solution(j-1))/bet
      end do
      dj(n) = dj(n)+(solution(n)-rhs_orig(n))*ei(n)
      sumdj = dj(n)
      fj = 1.0d0+(solution(n)-rhs_orig(n))/rhs_orig(n)
      solution(n) = rhs_orig(n)*fj
      do j = n-1, 1, -1
         solution(j) = solution(j) - gamma_elim(j+1)*solution(j+1)
         dj(j) = (solution(j)-rhs_orig(j))*ei(j)
         sumdj = sumdj + dj(j)
      end do

      return
end subroutine solve_angular_momentum_transport

end module am_transport_coefficients
