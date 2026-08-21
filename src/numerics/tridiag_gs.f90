!----------------------------------------------------------------------
! tridiag_gs
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original tridiag_gs.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0
! regression suite (examples/run_standard_solar_model).
!
! Thomas-algorithm tridiagonal solve (see also tridia.f90, the same
! algorithm operating through a shared common block instead of
! explicit dummy arguments): solves the tridiagonal system with
! sub-diagonal a, diagonal b, super-diagonal c, and right-hand side
! ex_prime, returning the solution in ex.
module tridiag_gs_mod
contains
subroutine tridiag_gs(a, b, c, ex_prime, npt, ex)
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: a(json), b(json), c(json), &
           ex_prime(json)
      integer, intent(in) :: npt
      double precision, intent(out) :: ex(json)

      double precision :: gama(json)
      integer :: j
      double precision :: bet
      save

      bet = b(1)
      ex(1) = ex_prime(1)/bet
      do j = 2, npt
         gama(j) = c(j-1)/bet
         bet = b(j) - a(j)*gama(j)
         if (bet.eq.0) stop '#TRIDIA:SINGULAR MATRIX'
         ex(j) = (ex_prime(j) - a(j)*ex(j-1))/bet
      end do
      do j = npt-1, 1, -1
         ex(j) = ex(j) - gama(j+1)*ex(j+1)
      end do
!  911  format(1p6e10.2)

      return
end subroutine tridiag_gs
end module tridiag_gs_mod
