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
! shared common block). The caller fills common/tridi/'s sub_diag/
! diag/super_diag/rhs before calling; tridia solves for solution(:)
! in place and also returns dj(:), the per-zone contribution of the
! change in solution weighted by ei(:), and sumdj, their sum.
! KC 2025-05-31 removed the unused ej dummy argument (see the
! commented-out original signature below).
!       SUBROUTINE TRIDIA(N,EI,EJ,DJ,SUMDJ)  ! KC 2025-05-31
subroutine tridia(n, ei, dj, sumdj)
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: ei(json)
      double precision, intent(out) :: dj(json)
      double precision, intent(out) :: sumdj

! common/tridi/: tridiagonal-solve work arrays (Thomas algorithm).
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution and gamma_elim are solver-internal work arrays
! that persist across calls via SAVE. Not referenced elsewhere in
! the already-converted reference files, so names are free to choose
! here (COMMON storage is positional, not name-matched).
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json), gamma_elim(json)
      common/tridi/ sub_diag, diag, super_diag, rhs, solution, gamma_elim

      double precision :: rhs_orig(json)
      integer :: i, j
      double precision :: bet, fj
      save

      dj(n) = gamma_elim(n)
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
