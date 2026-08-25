!----------------------------------------------------------------------
! hsolve
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original hsolve.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This is YREC's core Henyey relaxation solve: backward-sweep
! elimination then forward back-substitution of the block-tridiagonal
! correction system for the two dependent structure variables (P, T)
! at each of the M mesh shells. elim_coeff/elim_rhs are filled in by
! the caller (crrect.f) from the linearized structure-equation
! coefficients at each shell (henyey_coefficients.f90's elim_coeff/elim_rhs are the
! same arrays); surface_bc holds the outer boundary condition. The
! exact physical role of each of the 4 rows of elim_coeff/elim_rhs is
! not confidently rederived here (this is dense, unlabeled linear
! algebra in the original with no explanatory comments) -- names
! below describe structural role (pivot/coupling/correction terms),
! not verified physical meaning; treat with the same care as the
! original when modifying.
subroutine henyey_solve(num_shells, elim_coeff, elim_rhs, surface_bc)
      use star_info_lib, only: json

      implicit none

      integer, intent(in) :: num_shells
      double precision, intent(in) :: elim_coeff(4,2,json)
      double precision, intent(inout) :: elim_rhs(4,json)
      double precision, intent(in) :: surface_bc(6)
      double precision :: pivot_divisor, pivot_ratio, t_coupling_term
      double precision :: correction_p, correction_t
      double precision :: p_store, t_store
      integer :: shell_idx, im

      pivot_divisor = 1.0d0/(surface_bc(1) - elim_coeff(3,1,num_shells))
      pivot_ratio = (surface_bc(4) - elim_coeff(4,1,num_shells))*pivot_divisor
      elim_rhs(1,1) = surface_bc(3) - elim_rhs(3,num_shells)
      elim_rhs(2,1) = (surface_bc(6) - elim_rhs(4,num_shells)) - &
           pivot_ratio*elim_rhs(1,1)
      t_coupling_term = surface_bc(2) - elim_coeff(3,2,num_shells)
      correction_t = elim_rhs(2,1)/(surface_bc(5) - &
           elim_coeff(4,2,num_shells) - pivot_ratio*t_coupling_term)
      correction_p = pivot_divisor*(elim_rhs(1,1) - correction_t*t_coupling_term)
      do shell_idx = 2, num_shells
         im = num_shells + 2 - shell_idx
         p_store = elim_rhs(1,im) - correction_p*elim_coeff(1,1,im) - &
              correction_t*elim_coeff(1,2,im)
         t_store = elim_rhs(2,im) - correction_p*elim_coeff(2,1,im) - &
              correction_t*elim_coeff(2,2,im)
         elim_rhs(3,im) = elim_rhs(3,im) - correction_p*elim_coeff(3,1,im) - &
              correction_t*elim_coeff(3,2,im)
         elim_rhs(4,im) = elim_rhs(4,im) - correction_p*elim_coeff(4,1,im) - &
              correction_t*elim_coeff(4,2,im)
         elim_rhs(1,im) = correction_p
         elim_rhs(2,im) = correction_t
         correction_p = p_store
         correction_t = t_store
      end do
      elim_rhs(1,1) = correction_p
      elim_rhs(2,1) = correction_t
      elim_rhs(3,1) = elim_rhs(3,1) - correction_p*elim_coeff(3,1,1) - &
           correction_t*elim_coeff(3,2,1)
      elim_rhs(4,1) = elim_rhs(4,1) - correction_p*elim_coeff(4,1,1) - &
           correction_t*elim_coeff(4,2,1)

      return
end subroutine henyey_solve
