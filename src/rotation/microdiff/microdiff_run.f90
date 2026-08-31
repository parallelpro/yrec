!----------------------------------------------------------------------
! microdiff_run
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original microdiff_run.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! THIS PROGRAM EXECUTES THE DIFFUSING OF ELEMENTS. ALL OF THE RELEVANT
! PARAMETERS ARE DEFINED AT BOTH THE GRID POINTS AND MIDPOINTS.
!
! Part of the microdiff.f90 pipeline (see also microdiff_setup.f90,
! microdiff_mte.f90, microdiff_coefficients.f90, microdiff_etm.f90): performs
! the actual two-step Lax-Wendroff + implicit-second-derivative
! diffusion solve for one species (hydrogen, metals, or a light
! element in turn), calling microdiff_coefficients.f90 for the diffusion
! coefficients, lax_wendrof1.f/lax_wendrof2.f (not part of this batch)
! for the explicit first-derivative term, and implicit_diffusion_coeffs.f90 +
! tridiag_gs.f (not part of this batch) for the implicit second-
! derivative term.
! 2026 de-tramp (ROADMAP item 3): 27 arguments -> 13; the two
! equally-spaced grids are microdiff_grid records (eq, eq_mid).
module microdiff_run_lib
      implicit none
contains

subroutine microdiff_run(grid_spacing, timestep, total_mass, num_eq_points, &
     eq, species_fraction, hydrogen_dlnc_dr, &
     eq_mid, species_fraction_mid, hydrogen_dlnc_dr_mid, &
     atomic_weight_diffused, atomic_charge_diffused, species_col)
      use microdiff_mte_lib, only: microdiff_grid
      use microdiff_coefficients_lib
      use star_info_lib, only: star

      use star_info_lib
      use luout_lib
      use run_log_lib, only: solver_diagnostics
      use phys_const_lib
      use numerics_lib
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





      double precision :: diffused_abundance(json), diffused_abundance_mid(json), &
           diffused_abundance_orig(json), diffused_abundance_orig_mid(json), &
           diffused_abundance_prev(json), diffused_abundance_prime(json), &
           diffusion_coeff1(json), diffusion_coeff1_mid(json), &
           diffusion_coeff2(json), diffusion_coeff2_mid(json), &
           diffusion_coeff2_deriv_mid(json)
      double precision :: alpha(json), sub_diag(json), diag(json), &
           super_diag(json)
      integer :: i, iter, num_mid_points, max_change_zone
! use_generic_diffusion_vectors (originally LDOLI): set true so that
! the Lax-Wendroff routines (lax_wendrof1.f/lax_wendrof2.f, not part
! of this batch) use the current single diffused-species vector rather
! than any legacy metal-diffusion-specific vectors -- exact downstream
! effect not confirmed here, only the value set.
      logical :: use_generic_diffusion_vectors
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
!  SET THE FLAG LDOLI = TRUE, SO THAT THE LAX-WENDROF ROUTINES KNOW TO NOT
!  USE THE OLD METAL DIFFUSION VECTORS.
      use_generic_diffusion_vectors = .true.
!
! NOW USE THE EQUALLY SPACED GRID TO CALUCLATE DIFFUSION COEFFICIENTS IN THE
! METHOD OF THE THOUL ROUTINE. MUCH OF THE FOLLOWING CODE WAS TAKEN FROM
! SETUP_LISETT.F.
!
      call microdiff_coefficients(num_eq_points, species_fraction, eq, &
           diffusion_coeff1, diffusion_coeff2, hydrogen_dlnc_dr, &
           atomic_weight_diffused, atomic_charge_diffused, species_col)
!
! FIRST STEP OF TWO STEP LAX-WENDROFF METHOD. COMPUTE NEW ABUNDANCES AT ZONE
! MIDPOINTS USING THE LAX SCHEME :
! X(N+1/2,J+1/2)=1/2(X(N,J+1/2)-1/2(DT/DR)(COD1(N,J+1)-COD1(N,J))
! WHERE N IS THE TIME VARIABLE, J IS THE SPATIAL ONE, AND COD1 IS THE
! DIFFUSION COEFFICIENT.
!
      call lax_wendroff_step1(timestep, diffusion_coeff1, eq%mass, num_eq_points, &
           total_mass, diffused_abundance_mid, use_generic_diffusion_vectors)
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
      call microdiff_coefficients(num_mid_points, species_fraction_mid, eq_mid, &
           diffusion_coeff1_mid, diffusion_coeff2_mid, &
           hydrogen_dlnc_dr_mid, atomic_weight_diffused, &
           atomic_charge_diffused, species_col)
!
! USING THE NEW COEFFICIENTS, SOLVE FOR THE NEW RUN OF ABUNDANCES.
! NOTE : THE SR ACTUALLY RETURNS THE *CHANGE* IN THE ABUNDANCE AS A FUNCTION
! OF RADIUS, WHICH IS APPLIED TO THE ORIGINAL RUN OF ELEMENTS, RATHER THAN
! COMPUTING A NEW RUN OF ABUNDANCE AND TRANSFORMING IT BACK.  THIS IS DONE
! TO MINIMIZE ERRORS FROM THE INTERPOLATION.
!
      call lax_wendroff_step2(timestep, diffusion_coeff1_mid, eq_mid%mass, &
           diffused_abundance, num_eq_points, total_mass, &
           use_generic_diffusion_vectors)
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
!  START ITERATION LOOP FOR THE NEW RUN OF ABUNDANCES.
      do iter=1,star%ctrl%settling_num_iterations
!  FIND CHANGE IN D AT THE ZONE MIDPOINTS, GIVEN CHANGE IN D AT
!  THE ZONE CENTERS.
         do i = 2,num_eq_points
            diffused_abundance_mid(i) = 0.5d0*(diffused_abundance(i)+ &
                 diffused_abundance(i-1)-diffused_abundance_prev(i)- &
                 diffused_abundance_prev(i-1))
         enddo
!  STORE CURRENT RUN OF ABUNDANCES IN VECTOR ED_P. THE ITERATION
!  LOOP IS COMPLETED ONCE ED - ED_P IS EVERYWHERE LESS THAN THE TOLERANCE
!  GRTOL. ALSO ZERO OUT EQCOD2X_H, WHICH IS ONLY USED IN GRSETT (OLD SETTLING
!  ROUTINE), BUT RETAINED FOR LEGACY.
         do i = 1,num_eq_points-1
            diffused_abundance_prev(i) = diffused_abundance(i)
            diffusion_coeff2_deriv_mid(i) = 0.0
         enddo
         diffused_abundance_prev(num_eq_points) = diffused_abundance(num_eq_points)
!  GET NEW DIFFUSION COEFFICIENTS, TAKING INTO ACCOUNT THE CHANGE IN D
!  FROM THE PREVIOUS ITERATION.
!
         call implicit_diffusion_coeffs(alpha, diffusion_coeff2_mid, &
              diffused_abundance_mid, diffusion_coeff2_deriv_mid, sub_diag, &
              diag, super_diag, num_eq_points)
!
!  SOLVE THE TRIDIAGONAL MATRIX SYSTEM FOR THE NEW RUN OF LI.
!
         call tridiag_gs(sub_diag, diag, super_diag, diffused_abundance_prime, &
              num_eq_points, diffused_abundance)
!
!  CHECK TO SEE IF THE CORRECTIONS TO THE HYDROGEN ABUNDANCE ARE SMALL
!  ENOUGH TO EXIT.
         max_abundance_change = 0.0d0
         max_change_zone = 0
         do i = 1,num_eq_points
            dx = diffused_abundance(i)-diffused_abundance_prev(i)
            if(dx.gt.max_abundance_change)then
               max_abundance_change = dx
               max_change_zone = i
            endif
         enddo
! solver forensics (2026 run-log verbosity sweep)
         if (solver_diagnostics()) &
              write(run_log_unit,90)iter,max_abundance_change,max_change_zone
   90    format(1x,'ITERATION ',i3,' DXMAX ',1pe10.2,' IMAX ',i4)
!  EXIT ITERATION LOOP IF SYSTEM HAS CONVERGED.
         if(max_abundance_change.lt.star%ctrl%settling_tolerance)exit
      end do
      if (iter > star%ctrl%settling_num_iterations) then
      write(terminal_unit,110)star%ctrl%settling_tolerance,star%ctrl%settling_num_iterations, &
           max_abundance_change,max_change_zone
      write(run_log_unit,110)star%ctrl%settling_tolerance,star%ctrl%settling_num_iterations, &
           max_abundance_change,max_change_zone
  110 format(1x,'MICRODIFF FAILED TO CONVERGE TO WITHIN ',1pe9.3,' IN ',i3, &
           'ITERATIONS'/1x,'LAST ITERATION CHANGE IN D ',1pe9.3, &
           ' IN EQUALLY SPACED SHELL ',i5)
      end if
!
!  STORE THE RUN OF CHANGES IN THE DIFFUSED ELEMENT.
      do i = 1,num_eq_points
         species_fraction(species_col,i) = diffused_abundance(i)- &
              diffused_abundance_orig(i)
      enddo
      return
end subroutine microdiff_run

end module microdiff_run_lib
