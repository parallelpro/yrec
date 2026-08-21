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
! Analogue of tridia.f90 (the Henyey structure-equation tridiagonal
! solver): same Thomas-algorithm layout, used here for compositional
! (species-abundance) diffusion solves instead. Was originally shared
! with tridia.f90 via common/tridi/ (positional storage); converted
! (2026, GUIDELINES.md) to explicit arguments since this is real
! per-call data flow (matrix in, solution out), not global
! configuration -- see GUIDELINES.md's module-vs-argument distinction.
subroutine ctridi(n, sub_diag, diag, super_diag, rhs, solution)

      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution (originally U, the new species abundance) is
! solved for here.
      double precision, intent(in) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: solution(json)

! gamma_elim is solver-internal scratch, never read by any caller
! (unlike tridia.f90's gamma_elim(n), which doubled as a genuine
! cross-call input -- see tridia.f90's header note). SAVE preserved
! from the original common-block version though nothing here actually
! depends on it persisting across calls.
      double precision :: gamma_elim(json)
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
