!----------------------------------------------------------------------
! linear_solvers
!----------------------------------------------------------------------
! Linear systems (YREC's numerics_lib, split 2026-10):
!   ludcmp, lubksb     LU decomposition with partial pivoting and the
!                      back substitution (Numerical Recipes)
!   simeqc             Gauss-Jordan elimination with partial pivoting,
!                      for one or more right-hand sides
!   solve_tridiagonal  the Thomas algorithm (formerly tridiag_gs and
!                      ctridi, the same algorithm twice)
! A singular system returns ierr = 1.
module linear_solvers
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: ludcmp, lubksb, simeqc, solve_tridiagonal

contains

!----------------------------------------------------------------------
! ludcmp
!----------------------------------------------------------------------
! Numerical Recipes' ludcmp: replaces a(1:n,1:n) (leading dimension
! np) with the LU decomposition of a row permutation of itself; indx
! records the permutation and d is +1 or -1 for an even or odd number
! of interchanges. A zero row returns ierr = 1.
subroutine ludcmp(a, n, np, indx, d, ierr)
      real(dp), parameter :: tiny = 1.0d-20

      integer, intent(in) :: n, np
      real(dp), intent(inout) :: a(np,np)
      integer, intent(out) :: indx(n)
      real(dp), intent(out) :: d
      integer, intent(out) :: ierr

      integer :: i, imax, j, k
      real(dp) :: vv(n)
      real(dp) :: aamax, dum, sum

      ierr = 0

      d = 1.d0
      do i = 1, n
          aamax = 0.d0
          do j = 1, n
              if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
          end do
          if (aamax.eq.0d0) then
             write(*,*) 'ludcmp: singular matrix'
             ierr = 1
             return
          end if
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

!----------------------------------------------------------------------
! lubksb
!----------------------------------------------------------------------
! Numerical Recipes' lubksb: solves a x = b with ludcmp's LU form of a
! and its indx; b is overwritten with x.
subroutine lubksb(a, n, np, indx, b)
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

!----------------------------------------------------------------------
! simeqc
!----------------------------------------------------------------------
! Gauss-Jordan elimination with partial pivoting on a system matrix
! stored as a flat array, column-major: num_unknowns equations (rows)
! and num_cols columns (num_cols > num_unknowns: the augmented
! right-hand sides, overwritten with the solutions). A singular matrix
! returns ierr = 1. (2026-10: sized from the arguments, it was a fixed
! system_matrix(56), the 7x8 kemcom case; the singular case no longer
! writes "STOPPED AT 1010" to unit 5, which is standard input.)
pure subroutine simeqc(system_matrix, num_cols, num_unknowns, ierr)
      integer, intent(in) :: num_cols, num_unknowns
      real(dp), intent(inout) :: system_matrix(num_unknowns*num_cols)
      integer :: jj, j, jy, it, i, ij, imax, ia, ib, ic, id, ix, jx, ny, &
           n1, ig, ih
      real(dp) :: biga, swap_val

      integer, intent(out) :: ierr

      ierr = 0

      jj=-num_unknowns
      do j=1,num_unknowns
      jy=j+1
      jj=jj+num_unknowns+1
      biga=0.0d0
      it=jj-j
      do i=j,num_unknowns
      ij=it+i
      if(dabs(biga).ge.dabs(system_matrix(ij))) cycle
      biga=system_matrix(ij)
      imax=i
      end do
      if (dabs(biga).eq.0.0d0) then
      ierr = 1
      return
      end if
      ia=j+num_unknowns*(j-2)
      it=imax-j
      do i=j,num_cols
      ia=ia+num_unknowns
      ib=ia+it
      swap_val=system_matrix(ia)
      system_matrix(ia)=system_matrix(ib)
      system_matrix(ib)=swap_val
      system_matrix(ia)=system_matrix(ia)/biga
      end do
      if(j.eq.num_unknowns) exit
      ia=num_unknowns*(j-1)
      do ix=jy,num_unknowns
      ib=ia+ix
      it=j-ix
      do jx=jy,num_cols
      ic=num_unknowns*(jx-1)+ix
      id=ic+it
      system_matrix(ic)=system_matrix(ic)-system_matrix(ib)*system_matrix(id)
      end do
      end do
      end do
      ny=num_unknowns-1
      it=num_unknowns*num_unknowns
      do j=1,ny
      ia=it-j
      ic=num_unknowns*num_cols
      ib=ic-j
      do i=1,j
      system_matrix(ib)=system_matrix(ib)-system_matrix(ia)*system_matrix(ic)
      n1=num_cols-1
      ig=ib
      ih=ic
      do
      if(n1.le.num_unknowns) exit
      ig=ig-num_unknowns
      ih=ih-num_unknowns
      system_matrix(ig)=system_matrix(ig)-system_matrix(ia)*system_matrix(ih)
      n1=n1-1
      end do
      ia=ia-num_unknowns
      ic=ic-1
      end do
      end do
      return
end subroutine simeqc

!----------------------------------------------------------------------
! solve_tridiagonal
!----------------------------------------------------------------------
! The Thomas algorithm for the tridiagonal system with sub-diagonal a,
! diagonal b, super-diagonal c and right-hand side ex_prime, into ex
! (formerly tridiag_gs, and ctridi with a blanket SAVE of scratch that
! every call writes before it reads). A zero pivot returns ierr = 1.
subroutine solve_tridiagonal(a, b, c, ex_prime, npt, ex, ierr)
      integer, intent(in) :: npt
      real(dp), intent(in) :: a(npt), b(npt), c(npt), ex_prime(npt)
      real(dp), intent(out) :: ex(npt)
      integer, intent(out) :: ierr

      real(dp) :: gama(npt)
      integer :: j
      real(dp) :: bet

      ierr = 0
      bet = b(1)
      ex(1) = ex_prime(1)/bet
      do j = 2, npt
         gama(j) = c(j-1)/bet
         bet = b(j) - a(j)*gama(j)
         if (bet.eq.0) then
            write(*,*) 'solve_tridiagonal: singular matrix'
            ierr = 1
            return
         end if
         ex(j) = (ex_prime(j) - a(j)*ex(j-1))/bet
      end do
      do j = npt-1, 1, -1
         ex(j) = ex(j) - gama(j+1)*ex(j+1)
      end do
      return
end subroutine solve_tridiagonal

end module linear_solvers
