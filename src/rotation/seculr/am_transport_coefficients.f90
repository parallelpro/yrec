!----------------------------------------------------------------------
! am_transport_coefficients
!----------------------------------------------------------------------
! The coefficients of angular momentum transport: diffusion
! (am_diffusion_coeffs) and advection-diffusion, assembled into a banded
! system and eliminated (am_advection_diffusion_coeffs, banded_solver).
module am_transport_coefficients
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: am_diffusion_coeffs, am_advection_diffusion_coeffs

contains

! ---------------------------------------------------------------
!  DCOEFT SETS UP THE COEFFICIENTS FOR THE DIFFUSION DIFFERENCE EQUATION.
!
!  INPUT VARIABLES:
!  diffusion_coeff (ECOD) : DIFFUSION COEFFICIENTS AT THE EQUALLY SPACED
!     GRID POINTS.
!  grid_spacing (DR) : GRID SPACING.
!  timestep (DT) : TIMESTEP.
!  eq_moment_of_inertia (EI) : RUN OF MOMENTS OF INERTIA OF EQUALLY
!     SPACED GRID POINTS.
!  eq_angular_momentum (EJ) : RUN OF ANGULAR MOMENTUM OF EQUALLY SPACED
!     GRID POINTS (see am_transport_grid.f90's eq_angular_momentum).
!  eq_omega (EW) : RUN OF ANGULAR VELOCITY OF EQUALLY SPACED GRID POINTS.
!  num_eq_points (NTOT) : NUMBER OF EQUALLY SPACED GRID POINTS.
!  wind_loss_explicit (WIND1) : THE ANGULAR MOMENTUM LOSS FROM A
!     SURFACE C.Z. COMPUTED EXPLICITILY.
!  wind_loss_implicit (WIND2) : AS WIND1, BUT COMPUTED IMPLICITLY.
!  *NOTE: IF NOT APPLICABLE, WIND1 AND WIND2 ARE ZEROED OUT.
!
!  OUTPUT VARIABLES (sub_diag/diag/super_diag/rhs dummy arguments) :
!  THE ANGULAR VELOCITY OF SHELL I AT TIME N+1 (W(I,N+1)) IS A FUNCTION OF
!  W(I-1),W(I),AND W(I+1) AS DISCUSSED BELOW.  THIS CAN BE EXPRESSED
!  AS A TRIDIAGONAL MATRIX.
!  sub_diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I-1).
!  diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I).
!  super_diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I+1)
!  rhs(I) : THE ANGULAR VELOCITY TERMS AT THE BEGINNING OF THE TIMESTEP.
!  *NOTE: ANGULAR MOMENTUM LOSS FROM A SURFACE C.Z. IS SUBTRACTED FROM
!     THE LAST ELEMENT OF ARRAY rhs IF APPLICABLE, AND ALSO RETURNED
!     VIA surface_wind_loss_term FOR solve_angular_momentum_transport
!     (secular_transport) TO SEED dj(n) WITH (THIS WAS PREVIOUSLY
!     SMUGGLED THROUGH common/tridi/'s gamma_elim(num_eq_points) SLOT).
!
!  THE DIFFUSION EQUATION WE ARE SOLVING IS
!  dW/dT = 1/(4pi*RHO*R**2) 1/(I/M) d/dR(D*dW/dR)
!  WHERE D IS OUR DIFFUSION COEFFICIENT,W IS THE ANGULAR VELOCITY,
!  I IS THE MOMENT OF INERTIA,M IS MASS,T IS TIME,R IS RADIUS,
!  AND RHO IS DENSITY.
!  USING M = 4pi*RHO*R**2*DR THIS IS DIFFERENCED FOR SHELL I AS
!  W(I,N+1)-W(I,N)=(DT/DR)*(1/I)*(D(I+1/2)*(W(I+1,N+1)-W(I,N+1)) -
!  D(I-1/2)*(W(I,N+1)-W(I-1,N+1)))
!  WHERE INDEX N+1 REFERS TO VALUES AT THE END OF THE TIMESTEP,AND
!  INDEX N REFERS TO VALUES AT THE BEGINNING OF THE TIMESTEP.
!  THIS SYSTEM IS SUPPLEMENTED WITH APPROPRIATE BOUNDARY CONDITIONS.
!
! Simpler (single diffusion coefficient, no third/fourth-order terms)
! companion to am_advection_diffusion_coeffs.f90's 4-band angular-momentum-transport solve,
! and analogue of composition_diffusion_coeffs.f90 (the composition-transport version of the
! same tridiagonal setup).
! sub_diag/diag/super_diag/rhs/surface_wind_loss_term were originally
! shared with the solver (formerly tridia; secular_transport calls both)
! through common/tridi/ (positional storage; surface_wind_loss_term was
! smuggled through the gamma_elim(num_eq_points) slot specifically so
! the solver's first statement could read it before overwriting
! gamma_elim as pure solver scratch).
subroutine am_diffusion_coeffs(diffusion_coeff, grid_spacing, timestep, &
     eq_moment_of_inertia, eq_angular_momentum, eq_omega, num_eq_points, &
     wind_loss_explicit, wind_loss_implicit, fix_omega_at_surface, &
     sub_diag, diag, super_diag, rhs, surface_wind_loss_term)
      use rotation_scratch_lib
      use star_info_lib, only: json
      implicit none

      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: grid_spacing, timestep
      double precision, intent(in) :: eq_moment_of_inertia(json), &
           eq_angular_momentum(json), eq_omega(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: wind_loss_explicit, wind_loss_implicit
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: surface_wind_loss_term
      logical, intent(in) :: fix_omega_at_surface
      integer :: i
      double precision :: fact0, fact, facta

      fact0 = timestep/grid_spacing
      do i = 1,num_eq_points
         rhs(i) = eq_omega(i)
      end do
      surface_wind_loss_term = -0.5d0*(wind_loss_explicit+ &
           wind_loss_implicit)*eq_moment_of_inertia(num_eq_points)
      rhs(num_eq_points) = rhs(num_eq_points)*(1.0d0+ &
           surface_wind_loss_term/eq_angular_momentum(num_eq_points))
      if (.not.rotation_settings%use_diffusion_advection_transport) then
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.
         fact = fact0/eq_moment_of_inertia(1)
         sub_diag(1) = 0.0d0
         diag(1) = 1.0d0 + fact*diffusion_coeff(2)
         super_diag(1) = -fact*diffusion_coeff(2)
         do i = 2,num_eq_points-1
            fact = fact0/eq_moment_of_inertia(i)
            sub_diag(i) = -fact*diffusion_coeff(i)
            diag(i) = 1.0d0 + fact*(diffusion_coeff(i)+diffusion_coeff(i+1))
            super_diag(i) = -fact*diffusion_coeff(i+1)
         end do
!  LAST SHELL B.C. : SAME AS FIRST SHELL B.C.
         fact = fact0/eq_moment_of_inertia(num_eq_points)
         if (.not.fix_omega_at_surface) then
            sub_diag(num_eq_points) = -fact*diffusion_coeff(num_eq_points)
            diag(num_eq_points) = 1.0d0 + fact*diffusion_coeff(num_eq_points)
         else
            sub_diag(num_eq_points) = 0.0d0
            diag(num_eq_points) = 1.0d0
         end if
         super_diag(num_eq_points) = 0.0d0
      else
! IF LDIFAD=T, WE ARE SOLVING A COMBINED DIFFUSION/ADVECTION EQUATION.
! THIS ADDS A TERM D/DR(RHO*R**4*V*W) TO THE ORIGINAL D/DR(RHO*R**4*
! V*R*DW/DR) EQUATION.
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.
         fact = fact0/eq_moment_of_inertia(1)
         facta = 0.5d0*timestep/eq_moment_of_inertia(1)
         sub_diag(1) = 0.0d0
         diag(1) = 1.0d0 + fact*rot_scr%am_diffusive_coeff(2) - &
              facta*rot_scr%am_advective_coeff(2)
         super_diag(1) = -fact*rot_scr%am_diffusive_coeff(2) - &
              facta*rot_scr%am_advective_coeff(2)
         do i = 2,num_eq_points-1
            fact = fact0/eq_moment_of_inertia(i)
            facta = 0.5d0*timestep/eq_moment_of_inertia(i)
            sub_diag(i) = -fact*rot_scr%am_diffusive_coeff(i) + &
                 facta*rot_scr%am_advective_coeff(i)
            diag(i) = 1.0d0 + fact*(rot_scr%am_diffusive_coeff(i)+ &
                 rot_scr%am_diffusive_coeff(i+1)) + &
                 facta*(rot_scr%am_advective_coeff(i)-rot_scr%am_advective_coeff(i+1))
            super_diag(i) = -fact*rot_scr%am_diffusive_coeff(i+1) - &
                 facta*rot_scr%am_advective_coeff(i+1)
         end do
!  LAST SHELL B.C. : SAME AS FIRST SHELL B.C.
         fact = fact0/eq_moment_of_inertia(num_eq_points)
         facta = 0.5d0*timestep/eq_moment_of_inertia(num_eq_points)
         if (.not.fix_omega_at_surface) then
            sub_diag(num_eq_points) = -fact*rot_scr%am_diffusive_coeff(num_eq_points) &
                 + facta*rot_scr%am_advective_coeff(num_eq_points)
            diag(num_eq_points) = 1.0d0 + fact*rot_scr%am_diffusive_coeff(num_eq_points) &
                 + facta*rot_scr%am_advective_coeff(num_eq_points)
         else
            sub_diag(num_eq_points) = 0.0d0
            diag(num_eq_points) = 1.0d0
         end if
         super_diag(num_eq_points) = 0.0d0
      end if

      return
end subroutine am_diffusion_coeffs

! ---------------------------------------------------------------
! MHP 06/02
! DCOEFT SETS UP THE COEFFICIENTS FOR THE DIFFUSION DIFFERENCE EQUATION.
!
! INPUT VARIABLES:
! ECOD : DIFFUSION COEFFICIENTS AT THE EQUALLY SPACED GRID POINTS.
! grid_spacing (DR) : GRID SPACING.
! timestep (DT) : TIMESTEP (the full step; the routine works with
!    substep_timestep = timestep/num_substeps internally).
! eq_moment_of_inertia (EI) : RUN OF MOMENTS OF INERTIA OF EQUALLY
!    SPACED GRID POINTS.
! eq_omega (EW) : RUN OF ANGULAR VELOCITY OF EQUALLY SPACED GRID POINTS.
! num_eq_points (NTOT) : NUMBER OF EQUALLY SPACED GRID POINTS.
! wind_loss_explicit (WIND1) : THE ANGULAR MOMENTUM LOSS FROM A
!    SURFACE C.Z. COMPUTED EXPLICITILY.
! wind_loss_implicit (WIND2) : AS WIND1, BUT COMPUTED IMPLICITLY.
!    Rescaled per theta iteration into the local wind_loss_implicit_iter.
! *NOTE: IF NOT APPLICABLE, WIND1 AND WIND2 ARE ZEROED OUT.
!
! OUTPUT VARIABLES :
! THE ANGULAR VELOCITY OF SHELL I AT TIME N+1 (W(I,N+1)) IS A FUNCTION OF
! W(I-1),W(I),AND W(I+1) AS DISCUSSED BELOW.  THIS CAN BE EXPRESSED
! AS A TRIDIAGONAL MATRIX.
! A(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I-1).
! B(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I).
! C(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I+1)
! D(I) : THE ANGULAR VELOCITY TERMS AT THE BEGINNING OF THE TIMESTEP.
! *NOTE: ANGULAR MOMENTUM LOSS FROM A SURFACE C.Z. IS SUBTRACTED FROM
!    THE LAST ELEMENT OF ARRAY D IF APPLICABLE.
!
! THE DIFFUSION EQUATION WE ARE SOLVING IS
! dW/dT = 1/(4pi*RHO*R**2) 1/(I/M) d/dR(D*dW/dR)
! WHERE D IS OUR DIFFUSION COEFFICIENT,W IS THE ANGULAR VELOCITY,
! I IS THE MOMENT OF INERTIA,M IS MASS,T IS TIME,R IS RADIUS,
! AND RHO IS DENSITY.
! USING M = 4pi*RHO*R**2*DR THIS IS DIFFERENCED FOR SHELL I AS
! W(I,N+1)-W(I,N)=(DT/DR)*(1/I)*(D(I+1/2)*(W(I+1,N+1)-W(I,N+1)) -
! D(I-1/2)*(W(I,N+1)-W(I-1,N+1)))
! WHERE INDEX N+1 REFERS TO VALUES AT THE END OF THE TIMESTEP,AND
! INDEX N REFERS TO VALUES AT THE BEGINNING OF THE TIMESTEP.
! THIS SYSTEM IS SUPPLEMENTED WITH APPROPRIATE BOUNDARY CONDITIONS.
!
! WE ARE SOLVING A 4XN SYSTEM OF EQUATIONS.
! INITIAL DEFINITIONS:
! OMEGA IS STORED IN THE FIRST ENTRY
! D^2 OMEGA/DR^2 TERM IS STORED IN THE SECOND
! D OMEGA/DR TERM IS STORED IN THE THIRD
! D^3 OMEGA/DR^3 TERM IS STORED IN THE FOURTH
! NOTE THAT SINCE THE DERIVATIVES ARE DEFINED IN
! TERMS OF LOWER ORDER QUANTITIES WE DON'T NEED THE
! START OF TIMESTEP VALUES FOR THEM.
subroutine am_advection_diffusion_coeffs(grid_spacing, timestep, eq_moment_of_inertia, eq_omega, &
     num_eq_points, wind_loss_explicit, wind_loss_implicit, &
     eq_delta_angular_momentum, eq_mixing_diffusion_coeff, &
     sum_delta_angular_momentum, fix_omega_at_surface, diffusion_converged, ierr)
      use rotation_scratch_lib
      use star_info_lib, only: json
      use math_lib
      implicit none
! length of the per-coefficient-iteration convergence history printed
! on failure (bounded by rotation_settings%max_diffusion_iters)
      integer, parameter :: iter_history_max = 50
      double precision, parameter :: tiny = 1.0d-30

      double precision, intent(in) :: grid_spacing
      double precision, intent(in) :: timestep
      double precision, intent(in) :: eq_moment_of_inertia(json), &
           eq_omega(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: wind_loss_explicit
      double precision, intent(in) :: wind_loss_implicit
      double precision, intent(out) :: eq_delta_angular_momentum(json)
      double precision, intent(inout) :: eq_mixing_diffusion_coeff(json)
      double precision, intent(out) :: sum_delta_angular_momentum
      logical, intent(in) :: fix_omega_at_surface
      logical, intent(out) :: diffusion_converged
      integer, intent(out) :: ierr

      double precision :: coeff_matrix(band_nmax,10), rhs(band_nmax), &
           omega_working(json), max_omega_change_history(iter_history_max)
      integer :: max_omega_change_zone_history(iter_history_max)
      double precision :: third_order_ratio_factor(json), domega_dr(json), &
           omega_mid(json)
      double precision :: omega_prev_medium_iter(json), &
           omega_substep_start(json)
      double precision :: omega_curvature(json)
! locals
      integer :: timestep_cut_count, num_substeps, substep_idx, &
           theta_iter_idx, coeff_iter_idx, num_equations, i, j, ii, k
! substep_timestep: timestep/num_substeps, the step actually differenced.
! wind_loss_implicit_iter: wind_loss_implicit rescaled to the current
! theta iteration's surface omega.
      double precision :: substep_timestep, wind_loss_implicit_iter, &
           substep_time_sum, total_angular_momentum_start, substep_frac, &
           wind_saturation_threshold
      double precision :: omega_capped, omega_prev_capped, &
           domega_dr_prev_substep, advective_term1, &
           theta_term_n, theta_term_p, diffusive_term2, dt_over_dr, &
           fact_over_ei, facta_half_dt_over_ei_dr, wind_loss_delta, &
           max_omega_change, &
           max_omega_change_total, max_omega_change_medium_iter, &
           omega_change, omega_change_total, omega_change_medium_iter, &
           total_angular_momentum_updated, &
           angular_momentum_correction, omega_mid_new, domega_dr_mid_new, &
           advective_velocity_term, domega_dr_velocity_term, &
           domega_dr_velocity_term_alt, curvature_velocity_term, &
           third_deriv_velocity_term, total_velocity, total_velocity_alt, &
           mixing_diffusion_raw, wind_loss_delta_full_step
      integer :: max_omega_change_zone, max_omega_change_total_zone

      ierr = 0

      timestep_cut_count = 0
      num_substeps = 1
! STORE START OF TIMESTEP GRADIENTS AND
! AVERAGED OMEGAS.
      do i = 2,num_eq_points
         domega_dr(i) = rot_scr%dchi_dr_edge(i)*(eq_omega(i)-eq_omega(i-1))/ &
              grid_spacing
      end do
      num_equations = 4*num_eq_points-2
! LOOP FOR TIMESTEP CUTTING
      timestep_cut: do   ! (was label 5)
      substep_time_sum = 0.0d0
! STORE START OF TIMESTEP OMEGA VALUES
      wind_loss_implicit_iter = wind_loss_implicit
      total_angular_momentum_start = 0.0d0
      do i = 1,num_eq_points
         omega_working(i) = eq_omega(i)
         omega_prev_medium_iter(i) = eq_omega(i)
         omega_substep_start(i) = eq_omega(i)
         total_angular_momentum_start = total_angular_momentum_start + &
              eq_omega(i)*eq_moment_of_inertia(i)
      end do
      do substep_idx = 1,num_substeps
         substep_timestep = timestep/dfloat(num_substeps)
         substep_time_sum = substep_time_sum + substep_timestep
         substep_frac = substep_timestep/timestep
! LOOP FOR ITERATION ON THE D THETA/DT TERM;
! COEFFICIENTS UPDATED ONCE PER NNN
      do theta_iter_idx = 1,rotation_settings%max_diffusion_iters
! LOOP FOR ITERATION ON THE OTHER COEFFICIENTS
! THAT ARE FUNCTIONS OF OMEGA; UPDATED ONCE PER
! NN.
         diffusion_converged = .false.
! UPDATE JDOT TERM
         if (theta_iter_idx.gt.1 .and. wind_loss_implicit_iter.gt.0.0d0) then
! The original dadcoeft.f had an "if (LROSSBY)" Rossby-scaled
! saturation threshold here, but LROSSBY/PMMSOLTAU were plain
! never-assigned locals (not wired to COMMON/PMMWIND/), so the branch
! never executed; the plain threshold is the only behaviour.
            wind_saturation_threshold = rotation_settings%wind_saturation_omega
            omega_capped = min(eq_omega(num_eq_points), &
                 wind_saturation_threshold)
            omega_prev_capped = min(omega_working(num_eq_points), &
                 wind_saturation_threshold)
            wind_loss_implicit_iter = wind_loss_implicit* &
                 pow(omega_prev_capped/omega_capped, &
                 rotation_settings%wind_law_omega_exponent-1.0d0)* &
                 (omega_working(num_eq_points)/eq_omega(num_eq_points))
         end if
      do coeff_iter_idx = 1,rotation_settings%max_diffusion_iters
! COMPUTE THE DIFFUSION COEFFICIENTS FOR
! THE FIRST AND SECOND ORDER TERMS.
      omega_curvature(1) = 0.0d0
      do i = 2,num_eq_points-1
         omega_curvature(i) = (omega_working(i+1)-2.0d0*omega_working(i)+ &
              omega_working(i-1))
      end do
      omega_curvature(num_eq_points) = 0.0d0
      do i = 2,num_eq_points
         omega_mid(i) = 0.5d0*(omega_working(i)+omega_working(i-1)) &
              -0.125d0*(omega_curvature(i)-omega_curvature(i-1))
         domega_dr_prev_substep = rot_scr%dchi_dr_edge(i)*(omega_working(i)- &
              omega_working(i-1))/grid_spacing
         if (abs(domega_dr(i)).gt.tiny) then
            rot_scr%shear_diffusion_coeff_eqgrid(i) = rot_scr%shear_diffusion_coeff_eqgrid(i)* &
                 (domega_dr_prev_substep/domega_dr(i))**2
            rot_scr%gsf_diffusion_coeff_eqgrid(i) = rot_scr%gsf_diffusion_coeff_eqgrid(i)* &
                 (domega_dr_prev_substep/domega_dr(i))**2
         end if
         advective_term1 = rotation_settings%difad_velocity_scale*rot_scr%eq_velocity_coeff0(i)* &
              omega_mid(i)**2*(rot_scr%eq_velocity_coeff1a(i)+ &
              omega_mid(i)**2*rot_scr%eq_velocity_coeff1b(i))
! THE D THETA/DT TERMS (VTHN/VTHP IN THE ORIGINAL) ARE DISABLED: BOTH
! ARE IDENTICALLY ZERO ON EVERY ITERATION.  THE ADDITIONS BELOW ARE
! KEPT AS WRITTEN.
         theta_term_n = 0.0d0
         theta_term_p = 0.0d0
         diffusive_term2 = rotation_settings%difad_velocity_scale*rot_scr%eq_velocity_coeff0(i)* &
              omega_mid(i)**2*(rot_scr%eq_velocity_coeff2a(i)+ &
              rot_scr%eq_velocity_coeff2b(i))+theta_term_n
         advective_term1 = advective_term1 + theta_term_p
         rot_scr%am_advective_coeff(i) = 0.2d0*advective_term1*rot_scr%geometric_factor(i)
         rot_scr%am_diffusive_coeff(i) = (0.2d0*diffusive_term2+ &
              rot_scr%shear_diffusion_coeff_eqgrid(i)+rot_scr%gsf_diffusion_coeff_eqgrid(i))* &
              rot_scr%geometric_factor(i)
         rot_scr%am_diffusive_coeff(i) = max(0.0d0,rot_scr%am_diffusive_coeff(i))
      end do
! USE THE PRIOR RUN TO CORRECT THE
! DIFFUSION COEFFICIENTS FOR CHANGES IN OMEGA;
! FX1 IS USED FOR THE THIRD AND FOURTH ORDER TERMS.
      do i = 2,num_eq_points
         third_order_ratio_factor(i) = 2.0d0*omega_mid(i)/(eq_omega(i)+ &
              eq_omega(i-1))
      end do
! RIGHT HAND SIDE - STARTING OMEGA RUN
      do i = 1,num_eq_points
         ii = 1 + 4*(i-1)
         rhs(ii) = omega_substep_start(i)
         rhs(ii+1) = 0.d0
         rhs(ii+2) = 0.d0
         rhs(ii+3) = 0.d0
      end do
! INCLUDE ANGULAR MOMENTUM LOSS
      wind_loss_delta = -0.5d0*(wind_loss_explicit+wind_loss_implicit_iter)* &
           eq_moment_of_inertia(num_eq_points)*substep_frac
      rhs(4*num_eq_points-3) = rhs(4*num_eq_points-3)* &
           (1.0d0+wind_loss_delta/eq_moment_of_inertia(num_eq_points)/ &
           omega_substep_start(num_eq_points))
! GLOBAL FACTOR FOR THE DIFFUSION COEFFICIENTS
      dt_over_dr = substep_timestep/grid_spacing
! IF LDIFAD=T, WE ARE SOLVING A COMBINED DIFFUSION/ADVECTION EQUATION.
! THIS ADDS A TERM D/DR(RHO*R**4*V*W) TO THE ORIGINAL D/DR(RHO*R**4*
! V*R*DW/DR) EQUATION.  SINCE V DEPENDS ON OMEGA AND ITS FIRST THROUGH
! THIRD DERIVATIVES, WE CAN EFFECTIVELY RECAST THIS AS A SET OF 4 FIRST
! ORDER EQUATIONS.
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.
      fact_over_ei = dt_over_dr/eq_moment_of_inertia(1)
      facta_half_dt_over_ei_dr = 0.5d0*substep_timestep/eq_moment_of_inertia(1)/ &
           grid_spacing
! ZERO OUT INITIAL COEFFICIENT MATRIX
      do i = 1,num_equations
         do j = 1,10
            coeff_matrix(i,j) = 0.0d0
         end do
      end do
! OMEGA TERMS - ADVECTIVE
      coeff_matrix(1,5) = 1.0d0  - facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(2)
      coeff_matrix(1,9) = -facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(2)
! OMEGA TERMS - DW/DR
      coeff_matrix(1,7) = -fact_over_ei*rot_scr%dchi_dr_edge(2)*rot_scr%am_diffusive_coeff(2)
! OMEGA TERMS - D2W/DR2
      coeff_matrix(1,10) = -facta_half_dt_over_ei_dr*rot_scr%dchi_dr_edge(2)* &
           rot_scr%am_2nd_deriv_coeff(2)*third_order_ratio_factor(2)
! OMEGA TERMS - D3W/DR3
      coeff_matrix(1,8) = fact_over_ei*rot_scr%dchi_dr_edge(2)*rot_scr%am_3rd_deriv_coeff(2)* &
           third_order_ratio_factor(2)
! D^2W/DR^2 - SET TO ZERO AT THE LOWER BOUNDARY.
      coeff_matrix(2,5) = 1.0d0
! D OMEGA/DR TERMS
      coeff_matrix(3,3) = 1.0d0/grid_spacing
      coeff_matrix(3,5) = 1.0d0
      coeff_matrix(3,7) = -1.0d0/grid_spacing
! D^3W/DR^3
! APPLY B.C. TO THE HIGHEST ORDER TERM
      coeff_matrix(4,5) = 1.0d0
      coeff_matrix(4,9) = -1.0d0/3.0d0
      do ii = 2,num_eq_points-1
         fact_over_ei = dt_over_dr/eq_moment_of_inertia(ii)
         facta_half_dt_over_ei_dr = 0.5d0*substep_timestep/ &
              eq_moment_of_inertia(ii)/grid_spacing
         i = 1+(ii-1)*4
! OMEGA TERMS - ADVECTIVE
         coeff_matrix(i,1) = facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(ii)
         coeff_matrix(i,5) = 1.0d0 + facta_half_dt_over_ei_dr* &
              (rot_scr%am_advective_coeff(ii)-rot_scr%am_advective_coeff(ii+1))
         coeff_matrix(i,9) = -facta_half_dt_over_ei_dr*rot_scr%am_advective_coeff(ii+1)
! OMEGA TERMS - DW/DR
         coeff_matrix(i,3) = fact_over_ei*rot_scr%dchi_dr_edge(ii)* &
              rot_scr%am_diffusive_coeff(ii)
         coeff_matrix(i,7) = -fact_over_ei*rot_scr%dchi_dr_edge(ii+1)* &
              rot_scr%am_diffusive_coeff(ii+1)
! OMEGA TERMS - D2W/DR2
         coeff_matrix(i,2) = facta_half_dt_over_ei_dr*rot_scr%dchi_dr_edge(ii)* &
              rot_scr%am_2nd_deriv_coeff(ii)*third_order_ratio_factor(ii)
         coeff_matrix(i,6) = facta_half_dt_over_ei_dr* &
              (rot_scr%am_2nd_deriv_coeff(ii)*rot_scr%dchi_dr_edge(ii)* &
              third_order_ratio_factor(ii) &
              - rot_scr%am_2nd_deriv_coeff(ii+1)*rot_scr%dchi_dr_edge(ii+1)* &
              third_order_ratio_factor(ii+1))
         coeff_matrix(i,10) = -facta_half_dt_over_ei_dr* &
              rot_scr%am_2nd_deriv_coeff(ii+1)*rot_scr%dchi_dr_edge(ii+1)* &
              third_order_ratio_factor(ii+1)
! OMEGA TERMS - D3W/DR3
         coeff_matrix(i,4) = fact_over_ei*rot_scr%am_3rd_deriv_coeff(ii)* &
              rot_scr%dchi_dr_edge(ii)*third_order_ratio_factor(ii)
         coeff_matrix(i,8) = -fact_over_ei*rot_scr%am_3rd_deriv_coeff(ii+1)* &
              rot_scr%dchi_dr_edge(ii+1)*third_order_ratio_factor(ii+1)
! D^2W/DR^2 COEFFICIENTS
         coeff_matrix(i+1,2) = rot_scr%second_deriv_geom_factor_eqgrid(ii)* &
              rot_scr%dchi_dr_edge(ii)/grid_spacing*third_order_ratio_factor(ii)
         coeff_matrix(i+1,5) = 1.0d0
         coeff_matrix(i+1,6) = -rot_scr%second_deriv_geom_factor_eqgrid(ii+1)* &
              rot_scr%dchi_dr_edge(ii+1)/grid_spacing*third_order_ratio_factor(ii+1)
! D W/DR COEFFICIENTS
         coeff_matrix(i+2,3) = 1.0d0/grid_spacing
         coeff_matrix(i+2,5) = 1.0d0
         coeff_matrix(i+2,7) = -1.0d0/grid_spacing
! D^3 W/DR^3 COEFFICIENTS
         coeff_matrix(i+3,3) = rot_scr%dchi_dr_edge(ii)* &
              rot_scr%third_deriv_geom_factor_eqgrid(ii)/grid_spacing
         coeff_matrix(i+3,5) = 1.0d0
         coeff_matrix(i+3,7) = -rot_scr%dchi_dr_edge(ii)* &
              rot_scr%third_deriv_geom_factor_eqgrid(ii+1)/grid_spacing
      end do
!  LAST SHELL B.C. : SAME AS FIRST SHELL B.C.
      fact_over_ei = dt_over_dr/eq_moment_of_inertia(num_eq_points)
      facta_half_dt_over_ei_dr = 0.5d0*substep_timestep/ &
           eq_moment_of_inertia(num_eq_points)/grid_spacing
      i = 4*num_eq_points - 3
! ZERO OUT TERMS RELATED TO D2W/DR2 AT THE EDGES
      coeff_matrix(i-4,10) = 0.0d0
      coeff_matrix(5,2) = 0.0d0
      if (.not.fix_omega_at_surface) then
! OMEGA TERMS - ADVECTIVE
         coeff_matrix(i,1) = facta_half_dt_over_ei_dr* &
              rot_scr%am_advective_coeff(num_eq_points)
         coeff_matrix(i,5) = 1.0d0 + facta_half_dt_over_ei_dr* &
              rot_scr%am_advective_coeff(num_eq_points)
! OMEGA TERMS - DW/DR
         coeff_matrix(i,3) = fact_over_ei*rot_scr%dchi_dr_edge(num_eq_points)* &
              rot_scr%am_diffusive_coeff(num_eq_points)
! OMEGA TERMS - D2W/DR2
         coeff_matrix(i,2) = facta_half_dt_over_ei_dr* &
              rot_scr%dchi_dr_edge(num_eq_points)*rot_scr%am_2nd_deriv_coeff(num_eq_points)* &
              third_order_ratio_factor(num_eq_points)
! OMEGA TERMS - D3W/DR3
         coeff_matrix(i,4) = fact_over_ei*rot_scr%dchi_dr_edge(num_eq_points)* &
              rot_scr%am_3rd_deriv_coeff(num_eq_points)
      else
! HOLD OMEGA FIXED
         coeff_matrix(i,5) = 1.0d0
      end if
! D^2W/DR^2 COEFFICIENTS: ASSUME DW/DR = 0 ABOVE BOUNDARY
      coeff_matrix(i+1,5) = 1.0d0
! D W/DR COEFFICIENTS AND D^3 W/DR^3 COEFFICIENTS NOT USED,
! SINCE THEY ARE TREATED AS ZERO ABOVE THE UNSTABLE REGION
! NOW SOLVE THE BAND SYSTEM (NM = TOTAL NUMBER OF EQUATIONS, 4 PER
! SHELL).  THE SOLUTION OVERWRITES RHS.
      call banded_solver(coeff_matrix,num_equations,rhs, ierr)
      if (ierr /= 0) return
! CONVERT RUN OF OMEGA TO DELTA OMEGA
      sum_delta_angular_momentum = 0.0d0
! DWMAX - CHANGE FROM THE PREVIOUS LOW-LEVEL ITERATION
! (NN).
      max_omega_change = (rhs(1)-omega_working(1))/omega_working(1)
      max_omega_change_zone = 1
! DWW - ABSOLUTE CHANGE FROM THE INITIAL VECTOR OF OMEGA
      max_omega_change_total = (rhs(1)-eq_omega(1))/eq_omega(1)
      max_omega_change_total_zone = 1
! DWI - CHANGE FROM THE RESULT OF THE PRIOR MEDIUM
! LEVEL ITERATION (NNN).
      max_omega_change_medium_iter = (rhs(1)-omega_prev_medium_iter(1))/ &
           omega_prev_medium_iter(1)
      do i = 1,num_eq_points
         ii = 1+4*(i-1)
         eq_delta_angular_momentum(i) = (rhs(ii)-eq_omega(i))* &
              eq_moment_of_inertia(i)
         sum_delta_angular_momentum = sum_delta_angular_momentum + &
              eq_delta_angular_momentum(i)
         omega_change = (rhs(ii)-omega_working(i))/omega_working(i)
         omega_change_total = (rhs(ii)-eq_omega(i))/eq_omega(i)
         omega_change_medium_iter = (rhs(ii)-omega_prev_medium_iter(i))/ &
              omega_prev_medium_iter(i)
         if (abs(omega_change).gt.abs(max_omega_change)) then
            max_omega_change_zone = i
            max_omega_change = omega_change
         end if
         if (abs(omega_change_total).gt.abs(max_omega_change_total)) then
            max_omega_change_total_zone = i
            max_omega_change_total = omega_change_total
         end if
         if (abs(omega_change_medium_iter).gt.abs(max_omega_change_medium_iter)) &
              then
            max_omega_change_medium_iter = omega_change_medium_iter
         end if
      end do
! DAMP OUT LARGE (SPURIOUS) CHANGES IN OMEGA
      if (abs(max_omega_change).gt.0.1d0) then
         write(*,*) substep_idx,theta_iter_idx,coeff_iter_idx
         write(*,*) 'CORRECTIONS TOO LARGE ',max_omega_change, &
              max_omega_change_zone
         diffusion_converged = .false.
         exit
      end if
! (THE ORIGINAL'S DAMPING FACTOR ON THE UPDATE WAS HARDWIRED TO 1.)
      total_angular_momentum_updated = 0.0d0
      do i = 1,num_eq_points
         ii = 1+4*(i-1)
         omega_working(i) = omega_working(i)+ &
              (rhs(ii)-omega_working(i))
         total_angular_momentum_updated = total_angular_momentum_updated+ &
              omega_working(i)*eq_moment_of_inertia(i)
      end do
! ANGULAR MOMENTUM CONSERVATION ENFORCED
      angular_momentum_correction = (total_angular_momentum_start+ &
           wind_loss_delta)/total_angular_momentum_updated
      write(*,*) total_angular_momentum_start,wind_loss_delta, &
           total_angular_momentum_updated,angular_momentum_correction
      sum_delta_angular_momentum = 0.0d0
      do i = 1,num_eq_points
         omega_working(i) = omega_working(i)*angular_momentum_correction
         eq_delta_angular_momentum(i) = (omega_working(i)-eq_omega(i))* &
              eq_moment_of_inertia(i)
         sum_delta_angular_momentum = sum_delta_angular_momentum+ &
              eq_delta_angular_momentum(i)
      end do
      write(*,*) sum_delta_angular_momentum,wind_loss_delta
      max_omega_change_history(coeff_iter_idx) = max_omega_change
      max_omega_change_zone_history(coeff_iter_idx) = max_omega_change_zone
! DETERMINE IF RUN HAS CONVERGED
      if (abs(max_omega_change).lt.rotation_settings%convergence_tolerance) then
         diffusion_converged = .true.
         exit
      end if
      end do
      write(*,*) substep_idx,theta_iter_idx,coeff_iter_idx
      write(*,100) max_omega_change_total,max_omega_change_total_zone, &
           (max_omega_change_history(j),max_omega_change_zone_history(j), &
           j=1,coeff_iter_idx)
 100  format(' MAX D(OMEGA)/OMEGA',1pe12.3,' AT PT.',i5, &
         ' BY ITERATION'/5(1x,e11.3,i4))
      do k = 1,num_eq_points
         omega_prev_medium_iter(k) = omega_working(k)
      end do
      if (abs(max_omega_change_medium_iter).le.rotation_settings%convergence_tolerance .and. &
           theta_iter_idx.ge.2) then
         exit
      else
         diffusion_converged = .false.
      end if
      end do
      do k = 1,num_eq_points
         omega_substep_start(k) = omega_working(k)
      end do
      total_angular_momentum_start = total_angular_momentum_start+ &
           wind_loss_delta
      end do
      if (.not.diffusion_converged) then
         timestep_cut_count = timestep_cut_count + 1
         num_substeps = 2*num_substeps
         if (timestep_cut_count.gt.5) then
            write(*,*) 'TIMESTEP CUT MORE THAN 5 TIMES - RUN STOPPED'
            ! a nonzero ierr stops the run
            ierr = 1
            return
         else
            write(*,*) 'TIMESTEP CUT #',timestep_cut_count,' IN DADCOEFT'
         end if
         cycle timestep_cut
      end if
      exit timestep_cut
      end do timestep_cut
      write(*,914) rhs(1),eq_omega(1),rhs(2),rhs(3),rhs(4)
 914  format(1p5e12.3)
      do i = 2,num_eq_points
         ii = 1 + 4*(i-1)
         omega_mid_new = 0.5d0*(rhs(ii)+rhs(ii-4))
         domega_dr_mid_new = 0.5d0*(rhs(ii+1)+rhs(ii-3))
         advective_velocity_term = rot_scr%am_advective_coeff(i)/rot_scr%geometric_factor(i)
         domega_dr_velocity_term = rot_scr%am_diffusive_coeff(i)/rot_scr%geometric_factor(i)* &
              rhs(ii-2)*rot_scr%dchi_dr_edge(i)/omega_mid_new
         domega_dr_velocity_term_alt = 0.2d0*rotation_settings%difad_velocity_scale* &
              rot_scr%eq_velocity_coeff0(i)*omega_mid_new*rot_scr%dchi_dr_edge(i)* &
              (rot_scr%eq_velocity_coeff2a(i)+rot_scr%eq_velocity_coeff2b(i))*rhs(ii-2)
         if (rot_scr%am_diffusive_coeff(i).le.0.0d0) domega_dr_velocity_term_alt = 0.0d0
         curvature_velocity_term = rot_scr%am_2nd_deriv_coeff(i)/rot_scr%geometric_factor(i)* &
              rot_scr%dchi_dr_edge(i)*domega_dr_mid_new/omega_mid_new
         third_deriv_velocity_term = rot_scr%am_3rd_deriv_coeff(i)/rot_scr%geometric_factor(i)* &
              rhs(ii-1)*rot_scr%dchi_dr_edge(i)/omega_mid_new
         total_velocity = advective_velocity_term+domega_dr_velocity_term+ &
              curvature_velocity_term+third_deriv_velocity_term
         total_velocity_alt = advective_velocity_term+ &
              domega_dr_velocity_term_alt+curvature_velocity_term+ &
              third_deriv_velocity_term
         if (i.le.10 .or. (num_eq_points-i).le.10) then
         write(*,911) rhs(ii),eq_omega(i),advective_velocity_term, &
              domega_dr_velocity_term,curvature_velocity_term, &
              third_deriv_velocity_term,total_velocity, &
              rhs(ii+1),rhs(ii+2),rhs(ii+3)
         end if
         rot_scr%mixing_velocity_estimate(i) = 5.0d0*abs(rotation_settings%mixing_velocity_scale* &
              rotation_settings%es_mixing_scale*total_velocity_alt*rot_scr%equatorial_radius(i))
         mixing_diffusion_raw = (rot_scr%mixing_velocity_estimate(i)+ &
              rotation_settings%mixing_velocity_scale*rotation_settings%secular_shear_mixing_scale* &
              rot_scr%shear_diffusion_coeff_eqgrid(i)+rotation_settings%mixing_velocity_scale* &
              rotation_settings%gsf_mixing_scale*rot_scr%gsf_diffusion_coeff_eqgrid(i))
         eq_mixing_diffusion_coeff(i) = mixing_diffusion_raw* &
              rot_scr%mixing_geometric_factor(i)
 911  format(1p10e12.3)
      end do
      wind_loss_delta_full_step = wind_loss_delta/substep_frac
      write(*,*) sum_delta_angular_momentum,wind_loss_delta_full_step
      return
end subroutine am_advection_diffusion_coeffs

! ---------------------------------------------------------------
! Custom fixed-10-column banded-matrix eliminator for the angular-
! momentum-transport diffusion system built by am_advection_diffusion_coeffs.f90: a system
! of NM = 4*NTOT-2 equations (4 per shell: omega, d2(omega)/dr2,
! d(omega)/dr, d3(omega)/dr3) with a fixed bandwidth stored 10 columns
! wide per row in coeff_matrix, column 5 being the diagonal (columns
! 1-4 sub-diagonal, 6-10 super-diagonal). This is NOT a standard
! textbook band-diagonal solve (c.f. numerics' solve_tridiagonal, the plain
! tridiagonal Thomas algorithm) -- it is a hand-built Gaussian
! elimination exploiting the specific sparsity pattern that the 4-
! variables-per-shell formulation produces, with special-cased
! treatment of the first/last few rows where the band is truncated by
! the domain boundary. Transliterated as directly as possible from the
! original; the elimination structure itself is not independently
! re-derived here.
subroutine banded_solver(coeff_matrix, nm, rhs, ierr)
      use rotation_scratch_lib, only: band_nmax
      implicit none

      integer, intent(in) :: nm
      double precision, intent(inout) :: coeff_matrix(band_nmax,10), rhs(band_nmax)
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
            ! a nonzero ierr stops the run
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
            ! a nonzero ierr stops the run
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

end module am_transport_coefficients
