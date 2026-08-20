!----------------------------------------------------------------------
! lubksb
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lubksb.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! LU back-substitution (Numerical Recipes LUBKSB companion to
! LUDCMP): solves a*x=b for x, given the LU decomposition of a (in
! place, as produced by ludcmp) and its pivot record indx. b is
! overwritten with the solution.
subroutine lubksb(a, n, np, indx, b)
      implicit none
      integer, intent(in) :: n, np
      double precision, intent(in) :: a(np,np)
      integer, intent(in) :: indx(n)
      double precision, intent(inout) :: b(n)

      integer :: i, ii, j, ll
      double precision :: sum
      save

      ii = 0
      do i = 1, n
          ll = indx(i)
          sum = b(ll)
          b(ll) = b(i)
          if (ii.ne.0) then
              do j = ii, i - 1
                  sum = sum - a(i,j)*b(j)
              end do

          else if (sum.ne.0d0) then
              ii = i
          end if

          b(i) = sum
      end do
      do i = n, 1, -1
          sum = b(i)
          if (i.lt.n) then
              do j = i + 1, n
                  sum = sum - a(i,j)*b(j)
              end do
          end if

          b(i) = sum/a(i,i)
      end do
      return

end subroutine lubksb
