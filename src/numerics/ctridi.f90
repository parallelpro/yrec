!----------------------------------------------------------------------
! ctridi
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ctridi.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Solves a tridiagonal matrix system for the fractional abundance of a
! species. This routine is from Numerical Recipes, p.40.
!
! Shares common/tridi/ with tridia.f90 (the Henyey structure-equation
! tridiagonal solver): same six work arrays, same Thomas-algorithm
! layout, reused here for compositional (species-abundance) diffusion
! solves instead. The solution is left in the common block's solution
! array (originally named U here) for the caller to read back.
subroutine ctridi(n)

      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n

! common/tridi/: tridiagonal-solve work arrays (Thomas algorithm).
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution (originally U, the new species abundance) and
! gamma_elim are solver-internal work arrays that persist across
! calls via SAVE. Naming matches tridia.f90, which shares this exact
! common block for the structure-equation solve.
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json), gamma_elim(json)
      common/tridi/ sub_diag, diag, super_diag, rhs, solution, gamma_elim

      save

      double precision :: bet
      integer :: j

      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2,n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) stop '#TRIDIA:SINGULAR MATRIX'
         solution(j) = (rhs(j) - sub_diag(j)*solution(j-1))/bet
      end do
      do j = n-1,1,-1
         solution(j) = solution(j) - gamma_elim(j+1)*solution(j+1)
      end do

      return
end subroutine ctridi
