!----------------------------------------------------------------------
! microdiff_run
!----------------------------------------------------------------------
! microdiff's solve for one species (hydrogen, the metals, or one light
! element) over one step, on the two equally spaced grids (Bahcall &
! Loeb units). The species' abundance obeys
!     dX/dt = (4 pi/dm) d/dr [ D1 + D2 dX/dr ]   (dm the cell's mass),
! solved in two parts:
!   1. the drift term D1 explicitly, by the two-step Lax-Wendroff scheme
!      (diffusion_scheme): half a step at the midpoints, the
!      coefficients there from the abundance it gives, then the full
!      step at the centres;
!   2. the diffusion term D2 implicitly (solve_diffusion_term): one
!      tridiagonal solve (diffusion_scheme's implicit_diffusion_coeffs,
!      numerics_lib's solve_tridiagonal).
! The coefficients come from microdiff_coefficients. The species' change,
! not its new abundance, goes back to the model, which keeps the
! interpolation's error out of the unchanged part.
module microdiff_run
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use element_diffusion_configuration, only: diffusion_settings
      use diffusion_scheme, only: lax_wendroff_step1, lax_wendroff_step2, &
           implicit_diffusion_coeffs
      implicit none
      private
      public :: species_solve_report, run_microdiff

! What one species' solve did, for its driver's messages
! (element_diffusion_messages): the passes of the implicit solve the
! settings ask for (settling_num_iterations), counted until one's change
! is below settling_tolerance (pass 1 solves; a later pass would repeat
! it and change nothing), whether one was, pass 1's largest signed change
! of the abundance and its equally spaced shell; and, with ierr nonzero,
! the numerics routine that failed.
      type :: species_solve_report
         integer :: num_passes = 0
         logical :: converged = .false.
         real(dp) :: change = 0.0d0
         integer :: change_shell = 0
         character(len=24) :: failed_routine = ''
      end type species_solve_report

contains

! ---------------------------------------------------------------
! One species' change over timestep (Bahcall & Loeb units) on the
! num_eq_points centres (centre_grid) and the midpoints (mid_grid),
! grid_spacing apart, in a star of total_mass. The species is Thoul
! column species_col of species_fraction (centres) and
! species_fraction_mid (midpoints), with the atomic weight and charge
! given (column 3); on return species_fraction(species_col,:) holds its
! change at the centres. hydrogen_dlnc_dr(_mid): d ln C_H/dr, which the
! hydrogen solve computes and the others read. ierr: a numerics
! routine's status (report%failed_routine names it).
subroutine run_microdiff(grid_spacing, timestep, total_mass, num_eq_points, &
     centre_grid, species_fraction, hydrogen_dlnc_dr, &
     mid_grid, species_fraction_mid, hydrogen_dlnc_dr_mid, &
     atomic_weight_diffused, atomic_charge_diffused, species_col, report, ierr)
      use microdiff_mapping, only: microdiff_grid
      use microdiff_coefficients, only: compute_microdiff_coefficients
      real(dp), intent(in) :: grid_spacing, timestep, total_mass
      integer, intent(in) :: num_eq_points
      type(microdiff_grid), intent(in) :: centre_grid, mid_grid
      real(dp), intent(inout) :: species_fraction(:,:), species_fraction_mid(:,:)
      real(dp), intent(inout) :: hydrogen_dlnc_dr(:), hydrogen_dlnc_dr_mid(:)
      real(dp), intent(in) :: atomic_weight_diffused, atomic_charge_diffused
      integer, intent(in) :: species_col
      type(species_solve_report), intent(out) :: report
      integer, intent(out) :: ierr
! the abundance at the centres, and at the start of the step on both grids
      real(dp) :: abundance(json), abundance_start(json), abundance_mid_start(json)
! the drift's change at the midpoints over half the step
      real(dp) :: drift_change_mid(json)
! the drift (D1) and diffusion (D2) coefficients on the two grids
      real(dp) :: drift_coeff(json), drift_coeff_mid(json), &
           diffusion_coeff(json), diffusion_coeff_mid(json)
      integer :: n

      ierr = 0
      n = num_eq_points
      abundance(1:n) = species_fraction(species_col,1:n)
      abundance_start(1:n) = abundance(1:n)
      abundance_mid_start(1:n-1) = species_fraction_mid(species_col,1:n-1)

! The coefficients at the centres.
      call compute_microdiff_coefficients(n, species_fraction, centre_grid, &
           drift_coeff, diffusion_coeff, hydrogen_dlnc_dr, &
           atomic_weight_diffused, atomic_charge_diffused, species_col, ierr)
      if (ierr /= 0) then
         report%failed_routine = 'ludcmp'
         return
      end if
! Lax-Wendroff step 1: the drift over half the step at the midpoints, and
! the coefficients there from the abundance it gives (the midpoints are
! one fewer than the centres).
      call lax_wendroff_step1(timestep, drift_coeff, centre_grid%mass, n, &
           total_mass, drift_change_mid)
      species_fraction_mid(species_col,1:n-1) = abundance_mid_start(1:n-1) &
           + drift_change_mid(1:n-1)
      call compute_microdiff_coefficients(n-1, species_fraction_mid, mid_grid, &
           drift_coeff_mid, diffusion_coeff_mid, hydrogen_dlnc_dr_mid, &
           atomic_weight_diffused, atomic_charge_diffused, species_col, ierr)
      if (ierr /= 0) then
         report%failed_routine = 'ludcmp'
         return
      end if
! Lax-Wendroff step 2: the drift over the step at the centres.
      call lax_wendroff_step2(timestep, drift_coeff_mid, mid_grid%mass, &
           abundance, n, total_mass)
! The diffusion term.
      call solve_diffusion_term(grid_spacing, timestep, total_mass, n, mid_grid%mass, &
           diffusion_coeff_mid, abundance_start, drift_change_mid, abundance, report, ierr)
      if (ierr /= 0) return

      species_fraction(species_col,1:n) = abundance(1:n) - abundance_start(1:n)
end subroutine run_microdiff

! ---------------------------------------------------------------
! The diffusion term over timestep, implicitly: with cell masses
! dm_1 = m(1), dm_i = m(i) - m(i-1), dm_n = M - m(n-1) (m(i) = mid_mass(i),
! the mass inside midpoint i, so the first cell holds any convective core
! and the last the envelope) and alpha_i = 4 pi timestep/(grid_spacing
! dm_i), the system
!     X_i - alpha_i [D_i (X_(i+1) - X_i) - D_(i-1) (X_i - X_(i-1))] = X'_i
! (D the midpoints' diffusion_coeff_mid, X' the abundance after the
! drift) gives the abundance at the end of the step. The midpoint
! coefficients carry no derivative in the abundance here
! (implicit_diffusion_coeffs adds zero times the change in
! change_mid), so the system is solved once (docs/bug_sweep.md M253:
! the F77 iterated, every pass solving the same system). The passes the
! settings ask for are counted for the messages (species_solve_report):
! pass 1's largest signed change, then zero, until one is below
! settling_tolerance; with settling_num_iterations < 1 there is no solve.
! abundance: the drift's result in, the step's result out.
subroutine solve_diffusion_term(grid_spacing, timestep, total_mass, n, mid_mass, &
     diffusion_coeff_mid, abundance_start, change_mid, abundance, report, ierr)
      use physical_constants, only: c4pi
      use numerics_lib, only: solve_tridiagonal
      real(dp), intent(in) :: grid_spacing, timestep, total_mass
      integer, intent(in) :: n
      real(dp), intent(in) :: mid_mass(:), abundance_start(:)
      real(dp), intent(inout) :: diffusion_coeff_mid(:), change_mid(:), abundance(:)
      type(species_solve_report), intent(inout) :: report
      integer, intent(out) :: ierr
! alpha_i, the system's three diagonals, the abundance after the drift
! (the right-hand side) and before the solve, and the coefficients'
! derivative (zero)
      real(dp) :: alpha(json), sub_diag(json), diag(json), super_diag(json), &
           abundance_drift(json), abundance_before(json), zero_derivative(json)
      real(dp) :: alpha_scale, change, dx
      integer :: i, pass, change_shell

      ierr = 0
      alpha_scale = c4pi*timestep/grid_spacing
      alpha(1) = alpha_scale/mid_mass(1)
      do i = 2, n-1
         alpha(i) = alpha_scale/(mid_mass(i) - mid_mass(i-1))
      end do
      alpha(n) = alpha_scale/(total_mass - mid_mass(n-1))
      abundance_drift(1:n) = abundance(1:n)
      change = 0.0d0
      change_shell = 0
      if (diffusion_settings%settling_num_iterations >= 1) then
! the drift's change at the midpoints over the step
         do i = 2, n
            change_mid(i) = 0.5d0*(abundance(i) + abundance(i-1) - abundance_start(i) - &
                 abundance_start(i-1))
         end do
         abundance_before(1:n) = abundance(1:n)
         zero_derivative(1:n-1) = 0.0d0
         call implicit_diffusion_coeffs(alpha, diffusion_coeff_mid, change_mid, &
              zero_derivative, sub_diag, diag, super_diag, n)
         call solve_tridiagonal(sub_diag, diag, super_diag, abundance_drift, n, &
              abundance, ierr)
         if (ierr /= 0) then
            report%failed_routine = 'solve_tridiagonal'
            return
         end if
         do i = 1, n
            dx = abundance(i) - abundance_before(i)
            if (dx > change) then
               change = dx
               change_shell = i
            end if
         end do
      end if
      report%change = change
      report%change_shell = change_shell
      do pass = 1, diffusion_settings%settling_num_iterations
         report%num_passes = pass
         if (pass > 1) change = 0.0d0
         if (change < diffusion_settings%settling_tolerance) then
            report%converged = .true.
            exit
         end if
      end do
end subroutine solve_diffusion_term

end module microdiff_run

! Sources
! - G. Somers (05/2015), after the original routines' grsett (MHP, BC).
! - 2026: one implicit solve where the F77 iterated (docs/bug_sweep.md
!   M253), with its messages kept.
