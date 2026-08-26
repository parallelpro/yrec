!----------------------------------------------------------------------
! bandw
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original bandw.f; only variable names, source form, and comment
! style were updated (all-caps names in the preserved dead/commented-
! out code below are left as in the original, since that code does
! not execute). Validated against the Stage 0 regression suite
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
! re-derived here. The original carried an extensive alternate/older
! version of the elimination and backsubstitution as commented-out
! code; that dead code is preserved verbatim (as comments, unrenamed)
! below for historical reference -- none of it executes.
subroutine banded_solver(coeff_matrix, nm, rhs, ierr)
      implicit none
      integer, parameter :: json = 5000, nmax = 8000

      integer, intent(in) :: nm
      double precision, intent(inout) :: coeff_matrix(nmax,10), rhs(nmax)

      double precision, parameter :: tiny = 1.d-20

      integer :: i, j, k, imj, imj5, ii, jj
! INITIAL STEP: DIAGONALIZE TERMS INVOLVING OMEGA(1)
! AND THEN STEP THROUGH ITS FIRST THROUGH THIRD
! DERIVATIVES.  THEN PROCEED TO THE SECOND SHELL AND
! REPEAT, GOING TO THE BOTTOM.  THE REMAINDER OF THE
! TERMS CAN THEN BE REMOVED THROUGH BACKSUBSTITUTION.
!      DO I = 1,NM-5,4
! INITIAL STEP: DIAGONALIZE OMEGA(II) TERM.
      integer, intent(out) :: ierr

      ierr = 0

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
!         DO J = 5,10
            coeff_matrix(i,j) = coeff_matrix(i,j)/coeff_matrix(i,5)
         end do
         do j = i+1,i+4
            imj = i - j
            imj5 = 5 + imj
!            WRITE(*,*)I,J,IMJ5
            rhs(j) = rhs(j) - coeff_matrix(j,imj5)*rhs(i)
            do k = imj5+1,10+imj
!               WRITE(*,*)I,J,K,IMJ5,K+J-I
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
!         WRITE(*,*)I,B(I),A(I,5)
!         DO J = 5,5+II
         do j = 6,5+ii
            coeff_matrix(i,j) = coeff_matrix(i,j)/coeff_matrix(i,5)
         end do
!         WRITE(*,*)I,(J,A(I,J),J=6,5+II)
         do j = i+1,nm
            imj = i - j
            imj5 = 5 + imj
            rhs(j) = rhs(j) - coeff_matrix(j,imj5)*rhs(i)
!            WRITE(*,*)I,J,B(I),B(J),IMJ5,A(J,IMJ5)
            do k = imj5+1,imj5+ii
               coeff_matrix(j,k) = coeff_matrix(j,k)- &
                    coeff_matrix(j,imj5)*coeff_matrix(i,k+j-i)
!               WRITE(*,*)I,J,K,A(J,K),IMJ5,A(J,IMJ5),
!     *                   K+J-I,A(I,K+J-I)
            end do
         end do
      end do
      if (abs(coeff_matrix(nm,5)).lt.tiny) stop 999
      rhs(nm) = rhs(nm)/coeff_matrix(nm,5)

! D2W/DR2 TERM
! DW/DR (II) TERM
! OMEGA(II+1) TERM
! 3 NONZERO ENTRIES ABOVE I FOR THE NEXT TWO
! COLUMNS
!         DO J = I1,I2
! DIAGONALIZE ROW
! CORRECT RHS
!            DO K = J+1,J+3
!               K2 = 5 - K + J
!               B(K) = B(K) - A(K,K2)*B(J)
! CORRECT LHS
!               DO JJ = K3-1,K4
! D3W/DR3 TERM
!         K1 = 7
! CORRECT RHS
!         DO J = I3+1,I3+2
!            K2 = 5 - J + I3
!            B(J) = B(J) - A(J,K2)*B(I3)
! CORRECT LHS
!            K3 = K2 + 1
!            DO JJ = K3,K3+1
!            DO JJ = K3-1,K3+1
! FINAL 2 ENTRIES ARE A BIT SIMPLER.
      do i = nm,6,-1
! REMOVE ENTRIES IN LAST COLUMN
!         DO J = I-1,I-5,-1
!            B(J) = B(J) - B(I)*A(J,5+I-J)
         do j = -1,-5,-1
            rhs(i+j) = rhs(i+j) - rhs(i)*coeff_matrix(i+j,5-j)
         end do
!         DO J = I-2,I-5,-1
!            B(J) = B(J) - B(I-1)*A(J,4+I-J)
!         END DO
!         DO J = I-3,I-5,-1
!            B(J) = B(J) - B(I-2)*A(J,3+I-J)
!         END DO
!         DO J = I-4,I-5,-1
!            DO K = J-1,I-5,-1
      end do
      do i = 5,2,-1
         jj = -(i-1)
         do j = -1,jj,-1
            rhs(i+j) = rhs(i+j) - rhs(i)*coeff_matrix(i+j,5-j)
         end do
      end do
!      B(1) = B(1) - B(2)*A(1,6)
      return
end subroutine banded_solver
