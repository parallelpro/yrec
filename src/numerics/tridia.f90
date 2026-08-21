!----------------------------------------------------------------------
! tridia
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original tridia.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Thomas-algorithm tridiagonal solve (see also tridiag_gs.f90, the
! same algorithm operating on explicit dummy arguments instead of a
! shared common block, and ctridi.f90, the same algorithm for
! compositional diffusion). The caller fills sub_diag/diag/super_diag/
! rhs before calling; tridia solves for solution(:) in place and also
! returns dj(:), the per-zone contribution of the change in solution
! weighted by ei(:), and sumdj, their sum.
!
! Was originally shared with ctridi.f90 via common/tridi/ (positional
! storage); converted (2026, GUIDELINES.md) to explicit arguments
! since this is real per-call data flow, not global configuration.
! dj_n_seed replaces what was previously smuggled in via
! gamma_elim(n): the caller (dcoeft.f90, via seculr.f90) computes a
! surface wind-angular-momentum-loss term and used to stash it in the
! shared common block's gamma_elim(num_eq_points) slot specifically so
! this routine's very first statement (before gamma_elim is
! overwritten as pure solver scratch below) could pick it up as dj(n)'s
! initial value. That's now an explicit input instead of a COMMON
! side-channel; see dcoeft.f90's matching surface_wind_loss_term
! output argument.
! KC 2025-05-31 removed the unused ej dummy argument (see the
! commented-out original signature below).
!       SUBROUTINE TRIDIA(N,EI,EJ,DJ,SUMDJ)  ! KC 2025-05-31
subroutine tridia(n, ei, dj, sumdj, sub_diag, diag, super_diag, rhs, &
     solution, dj_n_seed)
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: ei(json)
      double precision, intent(out) :: dj(json)
      double precision, intent(out) :: sumdj
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution is solved for here. dj_n_seed is dj(n)'s initial
! value (see header note above).
      double precision, intent(in) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: solution(json)
      double precision, intent(in) :: dj_n_seed

! gamma_elim is solver-internal scratch (SAVE preserved from the
! original common-block version though nothing here actually depends
! on it persisting across calls, now that dj_n_seed carries the one
! value that used to be read from it across calls).
      double precision :: gamma_elim(json)

      double precision :: rhs_orig(json)
      integer :: i, j
      double precision :: bet, fj
      save

      dj(n) = dj_n_seed
      do i = 1, n
         rhs_orig(i) = rhs(i)
      end do
      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2, n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) stop '#TRIDIA:SINGULAR MATRIX'
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
end subroutine tridia
