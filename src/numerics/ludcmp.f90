!----------------------------------------------------------------------
! ludcmp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ludcmp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! LU decomposition with partial pivoting (Numerical Recipes LUDCMP):
! replaces a by the LU decomposition of a row-wise permutation of
! itself; indx records the row permutation, d returns +/-1 depending
! on whether the number of row interchanges was even or odd (used by
! the caller to get the sign of the determinant).
module ludcmp_mod
contains
subroutine ludcmp(a, n, np, indx, d)
      implicit none
      integer, parameter :: nmax = 100
      double precision, parameter :: tiny = 1.0d-20

      integer, intent(in) :: n, np
      double precision, intent(inout) :: a(np,np)
      integer, intent(out) :: indx(n)
      double precision, intent(out) :: d

      integer :: i, imax, j, k
      double precision :: vv(nmax)
      double precision :: aamax, dum, sum
      save

      d = 1.d0
      do i = 1, n
          aamax = 0.d0
          do j = 1, n
              if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
          end do
          if (aamax.eq.0d0) stop 'Singular matrix.'
          vv(i) = 1.d0/aamax
      end do
      do j = 1, n
          if (j.gt.1) then
              do i = 1, j - 1
                  sum = a(i,j)
                  if (i.gt.1) then
                      do k = 1, i - 1
                          sum = sum - a(i,k)*a(k,j)
                      end do
                      a(i,j) = sum
                  end if

              end do
          end if

          aamax = 0.d0
          do i = j, n
              sum = a(i,j)
              if (j.gt.1) then
                  do k = 1, j - 1
                      sum = sum - a(i,k)*a(k,j)
                  end do
                  a(i,j) = sum
              end if

              dum = vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax = i
                  aamax = dum
              end if

          end do
          if (j.ne.imax) then
              do k = 1, n
                  dum = a(imax,k)
                  a(imax,k) = a(j,k)
                  a(j,k) = dum
              end do
              d = -d
              vv(imax) = vv(j)
          end if

          indx(j) = imax
          if (j.ne.n) then
              if (a(j,j).eq.0d0) a(j,j) = tiny
              dum = 1.d0/a(j,j)
              do i = j + 1, n
                  a(i,j) = a(i,j)*dum
              end do
          end if

      end do
      if (a(n,n).eq.0d0) a(n,n) = tiny
      return

end subroutine ludcmp
end module ludcmp_mod
