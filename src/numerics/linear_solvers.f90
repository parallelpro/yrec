!----------------------------------------------------------------------
! linear_solvers
!----------------------------------------------------------------------
! Linear systems (YREC's numerics_lib, split 2026-10):
!   ludcmp, lubksb     LU decomposition with partial pivoting and the
!                      back substitution (Numerical Recipes)
!   simeqc             Gaussian elimination with partial pivoting and
!                      back substitution (IBM SSP's SIMQ), for one or
!                      more right-hand sides
!   solve_tridiagonal  the Thomas algorithm (formerly tridiag_gs and
!                      ctridi, the same algorithm twice)
! A singular system returns numerics_singular (numerics_status).
module linear_solvers
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      use numerics_status, only: numerics_ok, numerics_singular, numerics_bad_input
      implicit none
      private
      public :: ludcmp, lubksb, simeqc, solve_tridiagonal
      public :: tridiagonal_forward_sweep, tridiagonal_back_substitution

contains

!----------------------------------------------------------------------
! ludcmp
!----------------------------------------------------------------------
! Numerical Recipes' ludcmp: replaces a(1:n,1:n) (leading dimension
! np) with the LU decomposition of a row permutation of itself; indx
! records the permutation and d is +1 or -1 for an even or odd number
! of interchanges. The pivot of each column is the largest candidate
! scaled by its row's largest original element (implicit pivoting).
! ierr (numerics_status): numerics_singular for a zero row (nothing is
! done) or for a pivot whose scaled size is at most
! singular_pivot_tolerance, a matrix singular to working precision; the
! factorization then still completes, with a zero pivot replaced by
! 1e-20 as in Numerical Recipes, so that lubksb stays finite.
! numerics_bad_input when a column's candidate pivots are all NaN (the
! factorization stops there).
subroutine ludcmp(a, n, np, indx, d, ierr)
      real(dp), parameter :: tiny = 1.0d-20
! the largest scaled pivot, |pivot|/max|row|, of a matrix singular to
! working precision
      real(dp), parameter :: singular_pivot_tolerance = 1.0d-14

      integer, intent(in) :: n, np
      real(dp), intent(inout) :: a(np,np)
      integer, intent(out) :: indx(n)
      real(dp), intent(out) :: d
      integer, intent(out) :: ierr

      integer :: i, imax, j, k
      real(dp) :: vv(n)
      real(dp) :: aamax, dum, sum

      ierr = numerics_ok

      d = 1.d0
      do i = 1, n
          aamax = 0.d0
          do j = 1, n
              if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))
          end do
          if (aamax == 0d0) then
             ierr = numerics_singular
             return
          end if
          vv(i) = 1.d0/aamax
      end do
      do j = 1, n
          if (j > 1) then
              do i = 1, j - 1
                  sum = a(i,j)
                  if (i > 1) then
                      do k = 1, i - 1
                          sum = sum - a(i,k)*a(k,j)
                      end do
                      a(i,j) = sum
                  end if

              end do
          end if

          aamax = 0.d0
          imax = j
          do i = j, n
              sum = a(i,j)
              if (j > 1) then
                  do k = 1, j - 1
                      sum = sum - a(i,k)*a(k,j)
                  end do
                  a(i,j) = sum
              end if

              dum = vv(i)*abs(sum)
              if (dum >= aamax) then
                  imax = i
                  aamax = dum
              end if

          end do
          if (ieee_is_nan(a(imax,j))) then
             ierr = numerics_bad_input
             return
          end if
          if (aamax <= singular_pivot_tolerance) ierr = numerics_singular
          if (j /= imax) then
              do k = 1, n
                  dum = a(imax,k)
                  a(imax,k) = a(j,k)
                  a(j,k) = dum
              end do
              d = -d
              vv(imax) = vv(j)
          end if

          indx(j) = imax
          if (j /= n) then
              if (a(j,j) == 0d0) a(j,j) = tiny
              dum = 1.d0/a(j,j)
              do i = j + 1, n
                  a(i,j) = a(i,j)*dum
              end do
          end if

      end do
      if (a(n,n) == 0d0) a(n,n) = tiny
end subroutine ludcmp

!----------------------------------------------------------------------
! lubksb
!----------------------------------------------------------------------
! Numerical Recipes' lubksb: solves a x = b with ludcmp's LU form of a
! and its indx; b is overwritten with x.
pure subroutine lubksb(a, n, np, indx, b)
      integer, intent(in) :: n, np
      real(dp), intent(in) :: a(np,np)
      integer, intent(in) :: indx(n)
      real(dp), intent(inout) :: b(n)
      integer :: i, ii, j, ll
      real(dp) :: sum

      ii = 0
      do i = 1, n
          ll = indx(i)
          sum = b(ll)
          b(ll) = b(i)
          if (ii /= 0) then
              do j = ii, i - 1
                  sum = sum - a(i,j)*b(j)
              end do

          else if (sum /= 0d0) then
              ii = i
          end if

          b(i) = sum
      end do
      do i = n, 1, -1
          sum = b(i)
          if (i < n) then
              do j = i + 1, n
                  sum = sum - a(i,j)*b(j)
              end do
          end if

          b(i) = sum/a(i,i)
      end do
end subroutine lubksb

!----------------------------------------------------------------------
! simeqc
!----------------------------------------------------------------------
! Gaussian elimination with partial pivoting and back substitution
! (IBM SSP's SIMQ, extended to several right-hand sides) on the
! num_unknowns x num_cols system a = [A | B]: columns 1..num_unknowns
! hold A, the rest (num_cols > num_unknowns) the right-hand sides,
! which are overwritten with the solutions. A is destroyed. The pivot
! of column j is the first largest |a(i,j)|, i >= j. ierr
! (numerics_status): numerics_singular for a zero pivot. (2026-10: on a
! 2-D array, with the flat-index arithmetic of SIMQ written as
! subscripts in the same order of operations; sized from the arguments,
! it was a fixed system_matrix(56), the 7x8 kemcom case; the singular
! case no longer writes "STOPPED AT 1010" to unit 5, standard input.)
pure subroutine simeqc(a, num_cols, num_unknowns, ierr)
      integer, intent(in) :: num_cols, num_unknowns
      real(dp), intent(inout) :: a(num_unknowns, num_cols)
      integer, intent(out) :: ierr

      integer :: n, j, i, k, pivot_row, row, col
      real(dp) :: pivot, swap

      ierr = numerics_ok
      n = num_unknowns
! forward elimination: row j scaled by its pivot, then subtracted from
! the rows below
      do j = 1, n
         pivot = 0.0d0
         do i = j, n
            if (abs(pivot) >= abs(a(i,j))) cycle
            pivot = a(i,j)
            pivot_row = i
         end do
         if (abs(pivot) == 0.0d0) then
            ierr = numerics_singular
            return
         end if
         do col = j, num_cols
            swap = a(j,col)
            a(j,col) = a(pivot_row,col)
            a(pivot_row,col) = swap
            a(j,col) = a(j,col)/pivot
         end do
         if (j == n) exit
         do row = j+1, n
            do col = j+1, num_cols
               a(row,col) = a(row,col) - a(row,j)*a(j,col)
            end do
         end do
      end do
! back substitution, from row n-1 up: for each column k of the upper
! triangle, the last right-hand side first, then the others from the
! last but one down
      do j = 1, n-1
         row = n - j
         do i = 1, j
            k = n - i + 1
            a(row,num_cols) = a(row,num_cols) - a(row,k)*a(k,num_cols)
            do col = num_cols-1, n+1, -1
               a(row,col) = a(row,col) - a(row,k)*a(k,col)
            end do
         end do
      end do
end subroutine simeqc

!----------------------------------------------------------------------
! solve_tridiagonal
!----------------------------------------------------------------------
! The Thomas algorithm for the tridiagonal system with sub-diagonal a,
! diagonal b, super-diagonal c and right-hand side ex_prime, into ex:
! tridiagonal_forward_sweep, then tridiagonal_back_substitution
! (formerly tridiag_gs, and ctridi with a blanket SAVE of scratch that
! every call writes before it reads). ierr (numerics_status):
! numerics_singular for a zero pivot, the first included.
subroutine solve_tridiagonal(a, b, c, ex_prime, npt, ex, ierr)
      integer, intent(in) :: npt
      real(dp), intent(in) :: a(npt), b(npt), c(npt), ex_prime(npt)
      real(dp), intent(out) :: ex(npt)
      integer, intent(out) :: ierr

      real(dp) :: gama(npt)

      call tridiagonal_forward_sweep(a, b, c, ex_prime, npt, ex, gama, ierr)
      if (ierr /= numerics_ok) return
      call tridiagonal_back_substitution(gama, npt, ex)
end subroutine solve_tridiagonal

!----------------------------------------------------------------------
! tridiagonal_forward_sweep
!----------------------------------------------------------------------
! The Thomas algorithm's elimination of the sub-diagonal a (diagonal b,
! super-diagonal c, right-hand side r): with bet_1 = b(1) and, for
! j = 2..n, gama(j) = c(j-1)/bet_(j-1) and bet_j = b(j) - a(j) gama(j),
! x(j) = (r(j) - a(j) x(j-1))/bet_j. x(n) is then the solution's last
! component; tridiagonal_back_substitution gives the others. gama(1) is
! not set. ierr: numerics_singular for a zero pivot bet_j.
pure subroutine tridiagonal_forward_sweep(a, b, c, r, n, x, gama, ierr)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n), b(n), c(n), r(n)
      real(dp), intent(out) :: x(n), gama(n)
      integer, intent(out) :: ierr

      integer :: j
      real(dp) :: bet

      ierr = numerics_ok
      bet = b(1)
      if (bet == 0) then
         ierr = numerics_singular
         return
      end if
      x(1) = r(1)/bet
      do j = 2, n
         gama(j) = c(j-1)/bet
         bet = b(j) - a(j)*gama(j)
         if (bet == 0) then
            ierr = numerics_singular
            return
         end if
         x(j) = (r(j) - a(j)*x(j-1))/bet
      end do
end subroutine tridiagonal_forward_sweep

!----------------------------------------------------------------------
! tridiagonal_back_substitution
!----------------------------------------------------------------------
! The Thomas algorithm's back substitution, x(j) = x(j) - gama(j+1)
! x(j+1) for j = n-1 down to 1, on tridiagonal_forward_sweep's x and
! gama.
pure subroutine tridiagonal_back_substitution(gama, n, x)
      integer, intent(in) :: n
      real(dp), intent(in) :: gama(n)
      real(dp), intent(inout) :: x(n)

      integer :: j

      do j = n-1, 1, -1
         x(j) = x(j) - gama(j+1)*x(j+1)
      end do
end subroutine tridiagonal_back_substitution

end module linear_solvers
