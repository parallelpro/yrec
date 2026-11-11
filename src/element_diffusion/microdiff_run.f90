!----------------------------------------------------------------------
! microdiff_run
!----------------------------------------------------------------------
! THIS PROGRAM EXECUTES THE DIFFUSING OF ELEMENTS. ALL OF THE RELEVANT
! PARAMETERS ARE DEFINED AT BOTH THE GRID POINTS AND MIDPOINTS.
!
! microdiff's solve for one species (hydrogen, the metals, or one light
! element): the coefficients from microdiff_coefficients, the explicit
! first-derivative term by the two Lax-Wendroff steps, and the implicit
! second-derivative term by implicit_diffusion_coeffs and
! solve_tridiagonal (diffusion_scheme, numerics_lib), on the two equally
! spaced grids (eq, eq_mid).
module microdiff_run
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use element_diffusion_configuration, only: diffusion_settings
      use diffusion_scheme, only: lax_wendroff_step1, lax_wendroff_step2, implicit_diffusion_coeffs
      implicit none

! What one species' solve did, for its driver's messages
! (element_diffusion_messages): the passes of the implicit solve the
! settings ask for, counted until one's change is below
! settling_tolerance (pass 1 solves, each later pass would repeat it
! and changes nothing), whether one was, pass 1's largest signed change
! of the abundance and its equally spaced shell; and, with ierr nonzero,
! the numerics routine that failed.
      type :: species_solve_report
         integer :: num_passes = 0
         logical :: converged = .false.
         real(dp) :: change = 0.0d0
         integer :: change_shell = 0
         character(len=24) :: failed_routine = ''
      end type species_solve_report
contains

subroutine run_microdiff(grid_spacing, timestep, total_mass, num_eq_points, &
     eq, species_fraction, hydrogen_dlnc_dr, &
     eq_mid, species_fraction_mid, hydrogen_dlnc_dr_mid, &
     atomic_weight_diffused, atomic_charge_diffused, species_col, report, ierr)
      use microdiff_mapping, only: microdiff_grid
      use microdiff_coefficients
      use composition_layout, only: json
      use physical_constants
      use numerics_lib, only: solve_tridiagonal
      implicit none
      double precision, intent(in) :: grid_spacing
      double precision, intent(in) :: timestep
      double precision, intent(in) :: total_mass
      integer, intent(in) :: num_eq_points
      type(microdiff_grid), intent(in) :: eq, eq_mid
      double precision, intent(inout) :: species_fraction(3,json)
      double precision, intent(inout) :: hydrogen_dlnc_dr(json)
      double precision, intent(inout) :: species_fraction_mid(3,json)
      double precision, intent(inout) :: hydrogen_dlnc_dr_mid(json)
      double precision, intent(in) :: atomic_weight_diffused, &
           atomic_charge_diffused
      integer, intent(in) :: species_col
      type(species_solve_report), intent(out) :: report
      integer, intent(out) :: ierr
! --- locals ---
      double precision :: diffused_abundance(json), diffused_abundance_mid(json), &
           diffused_abundance_orig(json), diffused_abundance_orig_mid(json), &
           diffused_abundance_prev(json), diffused_abundance_prime(json), &
           diffusion_coeff1(json), diffusion_coeff1_mid(json), &
           diffusion_coeff2(json), diffusion_coeff2_mid(json), &
           diffusion_coeff2_deriv_mid(json)
      double precision :: alpha(json), sub_diag(json), diag(json), &
           super_diag(json)
      integer :: i, iter, num_mid_points, max_change_zone
      double precision :: fac, max_abundance_change, dx

!
! CREATE DUMMY VECTORS FOR THE DIFFUSING ELEMENT. THIS ALLOWS A SINGLE
! ELEMENT TO BE EASILY PASSED INTO SUB-ROUTINES.
!
      do i=1,num_eq_points-1
         diffused_abundance(i) = species_fraction(species_col,i)
         diffused_abundance_mid(i) = species_fraction_mid(species_col,i)
      enddo
      diffused_abundance(num_eq_points) = species_fraction(species_col,num_eq_points)
!
! STORE ORIGINAL RUN OF LIGHT ELEMENT MASS FRACTIONS IN VECTOR ED_ORIG.
! THE CHANGE IN ELEMENTS (ED - ED_ORIG) IS INTERPOLATED BACK TO THE
! ORIGINAL GRID AT THE END OF THE ROUTINE; THE CHANGE IN ABUNDANCE
! (RATHER THAN THE NEW RUN OF ABUNDANCE) IS USED TO MINIMIZE ERRORS
! ARISING FROM THE INTERPOLATION.
      do i = 1,num_eq_points-1
         diffused_abundance_orig(i) = diffused_abundance(i)
         diffused_abundance_orig_mid(i) = diffused_abundance_mid(i)
      enddo
      diffused_abundance_orig(num_eq_points) = diffused_abundance(num_eq_points)
!
      ierr = 0
!
! NOW USE THE EQUALLY SPACED GRID TO CALUCLATE DIFFUSION COEFFICIENTS IN THE
! METHOD OF THE THOUL ROUTINE. MUCH OF THE FOLLOWING CODE WAS TAKEN FROM
! SETUP_LISETT.F.
!
      call compute_microdiff_coefficients(num_eq_points, species_fraction, eq, &
           diffusion_coeff1, diffusion_coeff2, hydrogen_dlnc_dr, &
           atomic_weight_diffused, atomic_charge_diffused, species_col, ierr)
      if (ierr /= 0) then
         report%failed_routine = 'ludcmp'
         return
      end if
!
! FIRST STEP OF TWO STEP LAX-WENDROFF METHOD. COMPUTE NEW ABUNDANCES AT ZONE
! MIDPOINTS USING THE LAX SCHEME :
! X(N+1/2,J+1/2)=1/2(X(N,J+1/2)-1/2(DT/DR)(COD1(N,J+1)-COD1(N,J))
! WHERE N IS THE TIME VARIABLE, J IS THE SPATIAL ONE, AND COD1 IS THE
! DIFFUSION COEFFICIENT.
!
      call lax_wendroff_step1(timestep, diffusion_coeff1, eq%mass, num_eq_points, &
           total_mass, diffused_abundance_mid)
!
! UPDATE ESPEC WITH THE NEW RUN OF ED_H, FOR THE PURPOSE OF
! CALUCLATING NEW DIFFUSION COEFFICIENTS.
!
      do i=1,num_eq_points-1
         species_fraction_mid(species_col,i) = diffused_abundance_orig_mid(i) &
              + diffused_abundance_mid(i)
      enddo
!
! GET NEW DIFFUSION COEFFICIENTS AT THE ZONE MIDPOINTS, USING THE
! PROVISIONAL SOLUTION. NOTE THAT NPT-1 MUST BE PASSED IN, AS THE
! MIDPOINTS HAVE ONE LESS GRID ELEMENT.
!
      num_mid_points = num_eq_points-1
      call compute_microdiff_coefficients(num_mid_points, species_fraction_mid, eq_mid, &
           diffusion_coeff1_mid, diffusion_coeff2_mid, &
           hydrogen_dlnc_dr_mid, atomic_weight_diffused, &
           atomic_charge_diffused, species_col, ierr)
      if (ierr /= 0) then
         report%failed_routine = 'ludcmp'
         return
      end if
!
! USING THE NEW COEFFICIENTS, SOLVE FOR THE NEW RUN OF ABUNDANCES.
! NOTE : THE SR ACTUALLY RETURNS THE *CHANGE* IN THE ABUNDANCE AS A FUNCTION
! OF RADIUS, WHICH IS APPLIED TO THE ORIGINAL RUN OF ELEMENTS, RATHER THAN
! COMPUTING A NEW RUN OF ABUNDANCE AND TRANSFORMING IT BACK.  THIS IS DONE
! TO MINIMIZE ERRORS FROM THE INTERPOLATION.
!
      call lax_wendroff_step2(timestep, diffusion_coeff1_mid, eq_mid%mass, &
           diffused_abundance, num_eq_points, total_mass)
!
!  NOW IMPLICITLY SOLVE FOR THE SECOND TERM (INVOLVING THE SECOND
!  DERIVATIVE OF THE COMPOSITION GRADIENT). FIRST DEFINE VECTOR ED_P,
!  WHICH IS THE ABUNDANCE AT THE START OF THE TIMESTEP, AND ELI_PRIME,
!  THE ABUNDANCE ONE WOULD HAVE WITHOUT THE SECOND TERM IN THE DIFFUSION
!  EQUATION.
      do i = 1,num_eq_points
         diffused_abundance_prev(i) = diffused_abundance_orig(i)
         diffused_abundance_prime(i) = diffused_abundance(i)
      enddo
!  ALPHA IS THE NUMERICAL FACTOR IN FRONT OF THE DIFFUSION COEFFICIENTS.
      fac = c4pi*timestep/grid_spacing
      alpha(1) = fac/eq_mid%mass(1)
      do i = 2,num_eq_points-1
         alpha(i) = fac/(eq_mid%mass(i)-eq_mid%mass(i-1))
      enddo
      alpha(num_eq_points) = fac/(total_mass-eq_mid%mass(num_eq_points-1))
! The implicit solve, once. The midpoint coefficients carry no
! derivative in the abundance here (implicit_diffusion_coeffs adds zero
! times the change), so each pass of the original iteration solved the
! same system for the same right-hand side: the first pass's answer is
! the answer, and every later pass changed nothing (docs/bug_sweep.md
! M253). settling_num_iterations and settling_tolerance keep their
! effect on the messages: the first pass's change is the solve's, every
! later pass's change is zero (at shell 0), and the solve "fails to
! converge" when no pass's change is below the tolerance. With no pass
! (settling_num_iterations < 1) there is no implicit solve.
      max_abundance_change = 0.0d0
      max_change_zone = 0
      if (diffusion_settings%settling_num_iterations >= 1) then
!  FIND CHANGE IN D AT THE ZONE MIDPOINTS, GIVEN CHANGE IN D AT
!  THE ZONE CENTERS.
         do i = 2,num_eq_points
            diffused_abundance_mid(i) = 0.5d0*(diffused_abundance(i)+ &
                 diffused_abundance(i-1)-diffused_abundance_prev(i)- &
                 diffused_abundance_prev(i-1))
         enddo
!  STORE THE RUN OF ABUNDANCES BEFORE THE SOLVE IN diffused_abundance_prev,
!  AND ZERO THE COEFFICIENTS' DERIVATIVE (ONLY THE ORIGINAL ROUTINES HAVE ONE).
         do i = 1,num_eq_points-1
            diffused_abundance_prev(i) = diffused_abundance(i)
            diffusion_coeff2_deriv_mid(i) = 0.0
         enddo
         diffused_abundance_prev(num_eq_points) = diffused_abundance(num_eq_points)
         call implicit_diffusion_coeffs(alpha, diffusion_coeff2_mid, &
              diffused_abundance_mid, diffusion_coeff2_deriv_mid, sub_diag, &
              diag, super_diag, num_eq_points)
!
!  SOLVE THE TRIDIAGONAL MATRIX SYSTEM FOR THE NEW RUN OF LI.
!
         call solve_tridiagonal(sub_diag, diag, super_diag, diffused_abundance_prime, &
              num_eq_points, diffused_abundance, ierr)
         if (ierr /= 0) then
            report%failed_routine = 'solve_tridiagonal'
            return
         end if
!  THE LARGEST (SIGNED) CHANGE OF THE SOLVE.
         do i = 1,num_eq_points
            dx = diffused_abundance(i)-diffused_abundance_prev(i)
            if(dx.gt.max_abundance_change)then
               max_abundance_change = dx
               max_change_zone = i
            endif
         enddo
      end if
! The passes: pass 1's change, then zero, until one is below the
! tolerance.
      report%change = max_abundance_change
      report%change_shell = max_change_zone
      do iter=1,diffusion_settings%settling_num_iterations
         report%num_passes = iter
         if (iter > 1) then
            max_abundance_change = 0.0d0
         end if
         if(max_abundance_change.lt.diffusion_settings%settling_tolerance)then
            report%converged = .true.
            exit
         end if
      end do
!
!  STORE THE RUN OF CHANGES IN THE DIFFUSED ELEMENT.
      do i = 1,num_eq_points
         species_fraction(species_col,i) = diffused_abundance(i)- &
              diffused_abundance_orig(i)
      enddo
      return
end subroutine run_microdiff

end module microdiff_run
