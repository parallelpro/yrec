!----------------------------------------------------------------------
! bandw
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original bandw.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Custom fixed-10-column banded-matrix eliminator for the angular-
! momentum-transport diffusion system built by am_advection_diffusion_coeffs.f90: a system
! of NM = 4*NTOT-2 equations (4 per shell: omega, d2(omega)/dr2,
! d(omega)/dr, d3(omega)/dr3) with a fixed bandwidth stored 10 columns
! wide per row in coeff_matrix, column 5 being the diagonal (columns
! 1-4 sub-diagonal, 6-10 super-diagonal). This is NOT a standard
! textbook band-diagonal solve (c.f. tridia.f90/tridiag_gs.f90's plain
! tridiagonal Thomas algorithm) -- it is a hand-built Gaussian
! elimination exploiting the specific sparsity pattern that the 4-
! variables-per-shell formulation produces, with special-cased
! treatment of the first/last few rows where the band is truncated by
! the domain boundary. Transliterated as directly as possible from the
! original; the elimination structure itself is not independently
! re-derived here. (The original's commented-out older version of the
! elimination was removed in the 2026 readability sweep.)
subroutine banded_solver(coeff_matrix, nm, rhs, ierr)
      implicit none
      integer, parameter :: nmax = 8000

      integer, intent(in) :: nm
      double precision, intent(inout) :: coeff_matrix(nmax,10), rhs(nmax)
      integer, intent(out) :: ierr

      double precision, parameter :: tiny = 1.d-20

      integer :: i, j, k, imj, imj5, ii, jj

      ierr = 0
! FORWARD ELIMINATION: DIAGONALIZE THE TERMS INVOLVING OMEGA(1)
! AND THEN STEP THROUGH ITS FIRST THROUGH THIRD DERIVATIVES.  THEN
! PROCEED TO THE SECOND SHELL AND REPEAT, GOING TO THE BOTTOM.  THE
! REMAINDER OF THE TERMS ARE THEN REMOVED THROUGH BACKSUBSTITUTION.
      do i = 1,nm-5
         if (abs(coeff_matrix(i,5)).lt.tiny) then
            write(*,911) i,(coeff_matrix(i,j),j=1,10)
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
            ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
            ! preserve the historical stop on a nonzero return.
            ierr = 1
            return
         end if
 911     format(i5,1p10e12.3)
         rhs(i) = rhs(i)/coeff_matrix(i,5)
         do j = 6,10
            coeff_matrix(i,j) = coeff_matrix(i,j)/coeff_matrix(i,5)
         end do
         do j = i+1,i+4
            imj = i - j
            imj5 = 5 + imj
            rhs(j) = rhs(j) - coeff_matrix(j,imj5)*rhs(i)
            do k = imj5+1,10+imj
               coeff_matrix(j,k) = coeff_matrix(j,k)- &
                    coeff_matrix(j,imj5)*coeff_matrix(i,k+j-i)
            end do
         end do
      end do
      do i = nm-4,nm-1
         ii = nm-i
         if (abs(coeff_matrix(i,5)).lt.tiny) then
            write(*,911) i,(coeff_matrix(i,j),j=1,10)
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
            ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
            ! preserve the historical stop on a nonzero return.
            ierr = 1
            return
         end if
         rhs(i) = rhs(i)/coeff_matrix(i,5)
         do j = 6,5+ii
            coeff_matrix(i,j) = coeff_matrix(i,j)/coeff_matrix(i,5)
         end do
         do j = i+1,nm
            imj = i - j
            imj5 = 5 + imj
            rhs(j) = rhs(j) - coeff_matrix(j,imj5)*rhs(i)
            do k = imj5+1,imj5+ii
               coeff_matrix(j,k) = coeff_matrix(j,k)- &
                    coeff_matrix(j,imj5)*coeff_matrix(i,k+j-i)
            end do
         end do
      end do
      if (abs(coeff_matrix(nm,5)).lt.tiny) then
         write(*,*) 'banded_solver: singular final pivot'
         ierr = 1
         return
      end if
      rhs(nm) = rhs(nm)/coeff_matrix(nm,5)

! BACKSUBSTITUTION: EACH SOLVED ROW I REMOVES ITS CONTRIBUTION FROM
! THE (UP TO) FIVE ROWS ABOVE IT.
      do i = nm,6,-1
         do j = -1,-5,-1
            rhs(i+j) = rhs(i+j) - rhs(i)*coeff_matrix(i+j,5-j)
         end do
      end do
      do i = 5,2,-1
         jj = -(i-1)
         do j = -1,jj,-1
            rhs(i+j) = rhs(i+j) - rhs(i)*coeff_matrix(i+j,5-j)
         end do
      end do
      return
end subroutine banded_solver
