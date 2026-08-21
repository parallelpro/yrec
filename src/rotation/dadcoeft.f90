!----------------------------------------------------------------------
! dadcoeft
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original dadcoeft.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 06/02
! DCOEFT SETS UP THE COEFFICIENTS FOR THE DIFFUSION DIFFERENCE EQUATION.
!
! INPUT VARIABLES:
! ECOD : DIFFUSION COEFFICIENTS AT THE EQUALLY SPACED GRID POINTS.
! grid_spacing (DR) : GRID SPACING.
! timestep (DT) : TIMESTEP.
! eq_moment_of_inertia (EI) : RUN OF MOMENTS OF INERTIA OF EQUALLY
!    SPACED GRID POINTS.
! eq_omega (EW) : RUN OF ANGULAR VELOCITY OF EQUALLY SPACED GRID POINTS.
! num_eq_points (NTOT) : NUMBER OF EQUALLY SPACED GRID POINTS.
! wind_loss_explicit (WIND1) : THE ANGULAR MOMENTUM LOSS FROM A
!    SURFACE C.Z. COMPUTED EXPLICITILY.
! wind_loss_implicit (WIND2) : AS WIND1, BUT COMPUTED IMPLICITLY.
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
!
!       SUBROUTINE DADCOEFT(DR,DT,EI,EJ,EW,NTOT,WIND1,WIND2,DJ,  ! KC 2025-05-31
!                           ECOD2,SUMDJ,LFIX,LOKAD)
subroutine dadcoeft(grid_spacing, timestep, eq_moment_of_inertia, eq_omega, &
     num_eq_points, wind_loss_explicit, wind_loss_implicit, &
     eq_delta_angular_momentum, eq_mixing_diffusion_coeff, &
     sum_delta_angular_momentum, fix_omega_at_surface, diffusion_converged)
      use const_lib
      use light_burn_lib
      use turnover_lib
      implicit none
      integer, parameter :: json = 5000, nmax = 8000

      double precision, intent(in) :: grid_spacing
      double precision, intent(inout) :: timestep
      double precision, intent(in) :: eq_moment_of_inertia(json), &
           eq_omega(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: wind_loss_explicit
      double precision, intent(inout) :: wind_loss_implicit
      double precision, intent(out) :: eq_delta_angular_momentum(json)
      double precision, intent(inout) :: eq_mixing_diffusion_coeff(json)
      double precision, intent(out) :: sum_delta_angular_momentum
      logical, intent(in) :: fix_omega_at_surface
      logical, intent(out) :: diffusion_converged


! DIFFUSION COEFFICIENTS FOR THE TERMS INVOLVING D/DR (OMEGA)
! (ECOD3) AND D/DR (D OMEGA/DR) (ECOD4).
! common/difad/: am_advective_coeff/am_diffusive_coeff (originally
! ECOD3/ECOD4). Both used here. Naming matches rotgrid.f90.
      double precision :: am_advective_coeff(json), am_diffusive_coeff(json)
      common/difad/ am_advective_coeff, am_diffusive_coeff

! common/difad3/: higher-order terms in the meridional-circulation
! velocity (MHP 06/02). facd2/facd3 (second/third-order coefficients,
! used to build vesd2/vesd3), am_2nd_deriv_coeff/am_3rd_deriv_coeff
! (originally ECOD5/ECOD6, the diffusion coefficients for the
! d2(omega)/dr2 and d3(omega)/dr3 terms), geometric_factor (originally
! FGEOM, the rho*r^4*dchi/dr-type Jacobian factor multiplying the
! velocity-based terms to form diffusion coefficients), and the
! eq_velocity_coeff*/velocity_coeff* pairs (originally EV*/FV*, with
! FV1A/FV1B/FV2A/FV2B unused in this batch) which parameterize the
! circulation velocity as a polynomial in omega**2 -- exact physical
! role of each term beyond the algebra below is not confidently
! identified. shear_diffusion_coeff/gsf_diffusion_coeff (DSS/DGSF) are
! written but not read here; shear_diffusion_coeff_eqgrid/
! gsf_diffusion_coeff_eqgrid (ESS/EGSF) are the equally-spaced-grid
! shear/GSF diffusion coefficients, updated and read here. Naming is
! local to this batch.
      double precision :: facd2(json), facd3(json), vesd2(json), &
           vesd3(json), am_2nd_deriv_coeff(json), am_3rd_deriv_coeff(json), &
           geometric_factor(json), velocity_coeff0(json), &
           velocity_coeff1a(json), velocity_coeff1b(json), &
           velocity_coeff2a(json), velocity_coeff2b(json), &
           eq_velocity_coeff0(json), eq_velocity_coeff1a(json), &
           eq_velocity_coeff1b(json), eq_velocity_coeff2a(json), &
           eq_velocity_coeff2b(json), shear_diffusion_coeff(json), &
           gsf_diffusion_coeff(json), shear_diffusion_coeff_eqgrid(json), &
           gsf_diffusion_coeff_eqgrid(json)
      common/difad3/ facd2, facd3, vesd2, vesd3, am_2nd_deriv_coeff, &
           am_3rd_deriv_coeff, geometric_factor, velocity_coeff0, &
           velocity_coeff1a, velocity_coeff1b, velocity_coeff2a, &
           velocity_coeff2b, eq_velocity_coeff0, eq_velocity_coeff1a, &
           eq_velocity_coeff1b, eq_velocity_coeff2a, eq_velocity_coeff2b, &
           shear_diffusion_coeff, gsf_diffusion_coeff, &
           shear_diffusion_coeff_eqgrid, gsf_diffusion_coeff_eqgrid

! time change of theta
! common/difaddt/: ethvn/ethvp (originally ETHVN/ETHVP) are unused
! placeholders in this file -- the code that would read them is
! entirely commented out below, so the d(theta)/dt terms (theta_term_n/
! theta_term_p) are always zero in the live code. omega_avg_start/
! domega_dr_start (EWST/EQWRST) are read (to seed omega_mid_init/
! domega_dr_init on the first sub-step) but the values stored into
! those are not subsequently used either; kept as their own names since
! the intent is clear even though the result is unused. Naming is
! local to this batch.
      double precision :: ethvn(json), ethvp(json), omega_avg_start(json), &
           domega_dr_start(json)
      common/difaddt/ ethvn, ethvp, omega_avg_start, domega_dr_start


! VALUES OF D CHI/DR AT THE ZONE EDGES AND ZONE CENTERS
! common/egridchi/: dchi_dr_edge (QCHIRE) is used here; dchi_dr_center
! (QCHIRC) is not (its one reference is commented out). Naming is
! local to this batch.
      double precision :: dchi_dr_edge(json), dchi_dr_center(json)
      common/egridchi/ dchi_dr_edge, dchi_dr_center

! DEFINITION TERMS FOR THE SECOND AND THIRD DERIVATIVE
! TERMS
! common/egridder/: second/third-derivative geometric factors, raw
! grid (second_deriv_geom_factor/third_deriv_geom_factor, originally
! QQCOD/QQQCOD -- not used in this file) and their equally-spaced-grid
! counterparts (second_deriv_geom_factor_eqgrid/
! third_deriv_geom_factor_eqgrid, originally EQQCOD/EQQQCOD -- used
! here). Naming is local to this batch.
      double precision :: second_deriv_geom_factor_eqgrid(json), &
           third_deriv_geom_factor_eqgrid(json), second_deriv_geom_factor(json), &
           third_deriv_geom_factor(json)
      common/egridder/ second_deriv_geom_factor_eqgrid, &
           third_deriv_geom_factor_eqgrid, second_deriv_geom_factor, &
           third_deriv_geom_factor



! common/difad4/: mixing_geometric_factor/mixing_velocity_estimate/
! equatorial_radius (FGEOMIX/VESN/REQ), all used here. Naming is local
! to this batch.
      double precision :: mixing_geometric_factor(json), &
           mixing_velocity_estimate(json), equatorial_radius(json)
      common/difad4/ mixing_geometric_factor, mixing_velocity_estimate, &
           equatorial_radius

! MHP 8/17 ADDED EXCEN, C_2 TO COMMON BLOCK FOR MATT ET AL. 2012 CENT. TERM
! former common/cwind/: wind_saturation_omega/wind_law_omega_exponent
! (WMAX/EXW) are used here; the rest are unused placeholders (magnetic-
! braking law parameters -- exponents on mass-loss rate/convective
! turnover timescale/radius/mass/luminosity/Rossby number, and overall
! normalization constants, per the Matt et al. 2012-style wind-torque
! law this block evidently parameterizes).
!      COMMON/CWIND/WMAX,EXMD,EXW,EXTAU,EXR,EXM,CONSTFACTOR,STRUCTFACTOR,LJDOT0


!       DIMENSION EI(JSON),EW(JSON),EJ(JSON),DJ(JSON),  ! KC 2025-05-31
      double precision :: coeff_matrix(nmax,10), rhs(nmax), &
           omega_working(json), max_omega_change_history(50)
      integer :: max_omega_change_zone_history(50)
      double precision :: third_order_ratio_factor(json), domega_dr(json), &
           omega_mid(json)
      double precision :: omega_mid_start(json), omega_prev_medium_iter(json), &
           omega_prev_medium_iter_avg(json), omega_mid_prev(json), &
           domega_dr_prev(json), omega_substep_start(json), &
           omega_mid_init(json), domega_dr_init(json)
      double precision :: residual_check(nmax), rhs_orig(nmax), &
           coeff_matrix_orig(nmax,10), omega_curvature(json)
      save

! locals
      integer :: timestep_cut_count, num_substeps, substep_idx, &
           theta_iter_idx, coeff_iter_idx, num_equations, i, j, ii, k
      double precision :: full_timestep, wind_loss_implicit_initial, tiny, &
           substep_time_sum, total_angular_momentum_start, substep_frac, &
           wind_saturation_threshold
      logical :: lrossby
      double precision :: pmmsoltau
      double precision :: omega_capped, omega_prev_capped, omega_mid_it, &
           domega_dr_prev_substep, domega_dr_it, advective_term1, &
           theta_term_n, theta_term_p, diffusive_term2, dt_over_dr, &
           fact_over_ei, facta_half_dt_over_ei_dr, wind_loss_delta, &
           total_delta_angular_momentum_alt, max_omega_change, &
           max_omega_change_total, max_omega_change_medium_iter, &
           omega_change, omega_change_total, omega_change_medium_iter, &
           damping_factor, total_angular_momentum_updated, &
           angular_momentum_correction, omega_mid_new, domega_dr_mid_new, &
           advective_velocity_term, domega_dr_velocity_term, &
           domega_dr_velocity_term_alt, curvature_velocity_term, &
           third_deriv_velocity_term, total_velocity, total_velocity_alt, &
           mixing_diffusion_raw, wind_loss_delta_full_step
      integer :: max_omega_change_zone, max_omega_change_total_zone, &
           max_omega_change_medium_iter_zone

! DCOEFT SETS UP THE COEFFICIENTS FOR THE DIFFUSION DIFFERENCE EQUATION.
      timestep_cut_count = 0
      num_substeps = 1
      full_timestep = timestep
      wind_loss_implicit_initial = wind_loss_implicit
      tiny = 1.0d-30
! STORE START OF TIMESTEP GRADIENTS AND
! AVERAGED OMEGAS.
      do i = 2,num_eq_points
         domega_dr(i) = dchi_dr_edge(i)*(eq_omega(i)-eq_omega(i-1))/ &
              grid_spacing
         omega_mid_start(i) = 0.5d0*(eq_omega(i)+eq_omega(i-1))
      end do
      num_equations = 4*num_eq_points-2
! LOOP FOR TIMESTEP CUTTING
 5    continue
      substep_time_sum = 0.0d0
! STORE START OF TIMESTEP OMEGA VALUES
      wind_loss_implicit = wind_loss_implicit_initial
      total_angular_momentum_start = 0.0d0
      do i = 1,num_eq_points
         omega_working(i) = eq_omega(i)
         omega_prev_medium_iter(i) = eq_omega(i)
         omega_prev_medium_iter_avg(i) = eq_omega(i)
         omega_substep_start(i) = eq_omega(i)
         total_angular_momentum_start = total_angular_momentum_start + &
              eq_omega(i)*eq_moment_of_inertia(i)
      end do
      do i = 2,num_eq_points
         omega_mid_prev(i) = omega_mid_start(i)
         domega_dr_prev(i) = domega_dr(i)
      end do
      do substep_idx = 1,num_substeps
         timestep = full_timestep/dfloat(num_substeps)
         substep_time_sum = substep_time_sum + timestep
         substep_frac = timestep/full_timestep
! LOOP FOR ITERATION ON THE D THETA/DT TERM;
! COEFFICIENTS UPDATED ONCE PER NNN
      do theta_iter_idx = 1,itdif2
! LOOP FOR ITERATION ON THE OTHER COEFFICIENTS
! THAT ARE FUNCTIONS OF OMEGA; UPDATED ONCE PER
! NN.
         diffusion_converged = .false.
! UPDATE JDOT TERM
         if (theta_iter_idx.gt.1 .and. wind_loss_implicit.gt.0.0d0) then
! G Somers 8/17
! ADD ROSSBY SCALING IF DESIRED.
! lrossby/pmmsoltau (originally LROSSBY/PMMSOLTAU, from
! COMMON/PMMWIND/) are NOT declared via any COMMON block in the
! original dadcoeft.f -- an apparent pre-existing bug, since
! COMMON/PMMWIND/ is declared (with these same names) in amcalc.f,
! mwind.f, mcowind.f, and parmin.f. Left exactly as in the original:
! plain (uninitialized, SAVEd) locals, NOT wired to COMMON/PMMWIND/,
! so this IF branch reads whatever value a prior call happened to
! leave in them rather than the intended global Rossby-scaling flag.
            if (lrossby) then
               wind_saturation_threshold = wind_saturation_omega* &
                    pmmsoltau/turnover%convective_turnover_timescale
            else
               wind_saturation_threshold = wind_saturation_omega
            end if
! COMMENT OUT OLD WMAX STUFF
!C MHP 3/09 IF WMAX > 1 THEN ASSUME THAT THE PARAMETER WMAX IS DEFINED BY
!C WMAX = WMAX(SUN)*TAUCZ(SUN) AND THE SATURATION THRESHOLD WSAT = WMAX/TAUCZ(STAR)
!            IF(WMAX.GT.1.0D0)THEN
!               IF(TAUCZ.GT.1.0D0)THEN
!                  WSAT = WMAX/TAUCZ
!               ELSE
!                  WRITE(*,912)WMAX,TAUCZ
! 912      FORMAT('ERROR IN WIND - TAUCZ NOT DEFINED ',1P2E12.3,'STOPPED')
!                  STOP
!               ENDIF
!            ELSE
!               WSAT = WMAX
!            ENDIF
            omega_capped = min(eq_omega(num_eq_points), &
                 wind_saturation_threshold)
            omega_prev_capped = min(omega_working(num_eq_points), &
                 wind_saturation_threshold)
            wind_loss_implicit = wind_loss_implicit_initial* &
                 (omega_prev_capped/omega_capped)** &
                 (wind_law_omega_exponent-1.0d0)* &
                 (omega_working(num_eq_points)/eq_omega(num_eq_points))
         end if
      do coeff_iter_idx = 1,itdif2
! COMPUTE THE DIFFUSION COEFFICIENTS FOR
! THE FIRST AND SECOND ORDER TERMS.
      if (substep_idx.eq.1) then
         do i = 2,num_eq_points
            omega_mid_init(i) = 0.5d0*(omega_avg_start(i)+omega_avg_start(i-1))
            domega_dr_init(i) = domega_dr_start(i)
         end do
      else
         do i = 2,num_eq_points
            omega_mid_init(i) = omega_mid_prev(i)
            domega_dr_init(i) = domega_dr_prev(i)
         end do
      end if
      omega_curvature(1) = 0.0d0
      do i = 2,num_eq_points-1
         omega_curvature(i) = (omega_working(i+1)-2.0d0*omega_working(i)+ &
              omega_working(i-1))
      end do
      omega_curvature(num_eq_points) = 0.0d0
      do i = 2,num_eq_points
         omega_mid(i) = 0.5d0*(omega_working(i)+omega_working(i-1)) &
              -0.125d0*(omega_curvature(i)-omega_curvature(i-1))
!         WM(I) = 0.5D0*(WM(I)+WMINIT(I))
         omega_mid_it = 0.5d0*(omega_prev_medium_iter_avg(i)+ &
              omega_prev_medium_iter_avg(i-1))
         domega_dr_prev_substep = dchi_dr_edge(i)*(omega_working(i)- &
              omega_working(i-1))/grid_spacing
!         QWR2 = 0.5D0*(QWR2+QWRINIT(I))
         domega_dr_it = dchi_dr_edge(i)*(omega_prev_medium_iter_avg(i)- &
              omega_prev_medium_iter_avg(i-1))/grid_spacing
         if (abs(domega_dr(i)).gt.tiny) then
            shear_diffusion_coeff_eqgrid(i) = shear_diffusion_coeff_eqgrid(i)* &
                 (domega_dr_prev_substep/domega_dr(i))**2
            gsf_diffusion_coeff_eqgrid(i) = gsf_diffusion_coeff_eqgrid(i)* &
                 (domega_dr_prev_substep/domega_dr(i))**2
         end if
         advective_term1 = difad_velocity_scale*eq_velocity_coeff0(i)* &
              omega_mid(i)**2*(eq_velocity_coeff1a(i)+ &
              omega_mid(i)**2*eq_velocity_coeff1b(i))
!         VTH = FW*(ETHVN(I)*WM(I)*QWR2-ETHVP(I))/DT
         theta_term_n = 0.0d0
!         VTHN = FW*ETHVN(I)*WM(I)**2/DT
!         VTHP = -FW*ETHVP(I)/DT
!          VTHP = FW/DT0*(ETHVN(I)*WMI*QWRI-ETHVP(I))
!C          VTHP0 = FW/DT0*(ETHVN(I)*WM0(I)*QWR(I)-ETHVP(I))
!C          VTHP1 = FW/DT*ETHVN(I)*(WMI*QWRI-WMP(I)*QWRP(I))
!C          IF(NSTEP.LE.2)THEN
!C             VTHP = VTHP0+VTHP1
!C          ELSE
             theta_term_p = 0.0d0
!C          ENDIF
!         VTHN = FW*WM(I)*(ETHVN(I)*WM(I)-ETHVP(I)/QWR2)/DT
!         VTHN = MAX(0.0D0,VTHN)
! SET D THETA/DT TERM TO ZERO ON THE FIRST ITERATION
         if (theta_iter_idx.eq.1) then
!         IF(NTOT.GT.1)THEN
            theta_term_n=0.0d0
            theta_term_p=0.0d0
         end if
         diffusive_term2 = difad_velocity_scale*eq_velocity_coeff0(i)* &
              omega_mid(i)**2*(eq_velocity_coeff2a(i)+ &
              eq_velocity_coeff2b(i))+theta_term_n
         advective_term1 = advective_term1 + theta_term_p
!         ECOD3(I) = 0.2D0*(V1+VTH)*FGEOM(I)
         am_advective_coeff(i) = 0.2d0*advective_term1*geometric_factor(i)
         am_diffusive_coeff(i) = (0.2d0*diffusive_term2+ &
              shear_diffusion_coeff_eqgrid(i)+gsf_diffusion_coeff_eqgrid(i))* &
              geometric_factor(i)
         am_diffusive_coeff(i) = max(0.0d0,am_diffusive_coeff(i))
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
         rhs_orig(ii) = omega_substep_start(i)
         rhs_orig(ii+1) = 0.d0
         rhs_orig(ii+2) = 0.d0
         rhs_orig(ii+3) = 0.d0
      end do
! INCLUDE ANGULAR MOMENTUM LOSS
      wind_loss_delta = -0.5d0*(wind_loss_explicit+wind_loss_implicit)* &
           eq_moment_of_inertia(num_eq_points)*substep_frac
      rhs(4*num_eq_points-3) = rhs(4*num_eq_points-3)* &
           (1.0d0+wind_loss_delta/eq_moment_of_inertia(num_eq_points)/ &
           omega_substep_start(num_eq_points))
      rhs_orig(4*num_eq_points-3) = rhs(4*num_eq_points-3)
! GLOBAL FACTOR FOR THE DIFFUSION COEFFICIENTS
      dt_over_dr = timestep/grid_spacing
! IF LDIFAD=T, WE ARE SOLVING A COMBINED DIFFUSION/ADVECTION EQUATION.
! THIS ADDS A TERM D/DR(RHO*R**4*V*W) TO THE ORIGINAL D/DR(RHO*R**4*
! V*R*DW/DR) EQUATION.  SINCE V DEPENDS ON OMEGA AND ITS FIRST THROUGH
! THIRD DERIVATIVES, WE CAN EFFECTIVELY RECAST THIS AS A SET OF 4 FIRST
! ORDER EQUATIONS.
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.
      fact_over_ei = dt_over_dr/eq_moment_of_inertia(1)
      facta_half_dt_over_ei_dr = 0.5d0*timestep/eq_moment_of_inertia(1)/ &
           grid_spacing
! ZERO OUT INITIAL COEFFICIENT MATRIX
      do i = 1,num_equations
         do j = 1,10
            coeff_matrix(i,j) = 0.0d0
         end do
      end do
! OMEGA TERMS - ADVECTIVE
      coeff_matrix(1,5) = 1.0d0  - facta_half_dt_over_ei_dr*am_advective_coeff(2)
      coeff_matrix(1,9) = -facta_half_dt_over_ei_dr*am_advective_coeff(2)
! OMEGA TERMS - DW/DR
      coeff_matrix(1,7) = -fact_over_ei*dchi_dr_edge(2)*am_diffusive_coeff(2)
! OMEGA TERMS - D2W/DR2
!      A(1,6) = -FACTA*ECOD5(2)*FX1(2)
!      FPL = EQQCOD(2)*QCHIRE(2)*FX1(2)
      coeff_matrix(1,10) = -facta_half_dt_over_ei_dr*dchi_dr_edge(2)* &
           am_2nd_deriv_coeff(2)*third_order_ratio_factor(2)
!     *          +0.125D0*ECOD3(2)*DR**2/FPL
! OMEGA TERMS - D3W/DR3
      coeff_matrix(1,8) = fact_over_ei*dchi_dr_edge(2)*am_3rd_deriv_coeff(2)* &
           third_order_ratio_factor(2)
! D^2W/DR^2 - SET TO ZERO AT THE LOWER BOUNDARY.
      coeff_matrix(2,5) = 1.0d0
!      A(2,6) = -EQQCOD(2)*QCHIRE(2)/DR*FX1(2)
! D OMEGA/DR TERMS
      coeff_matrix(3,3) = 1.0d0/grid_spacing
      coeff_matrix(3,5) = 1.0d0
      coeff_matrix(3,7) = -1.0d0/grid_spacing
! D^3W/DR^3
! APPLY B.C. TO THE HIGHEST ORDER TERM
!      A(4,3) = EQQQCOD(1)*QCHIRE(2)/DR
      coeff_matrix(4,5) = 1.0d0
!      A(4,7) = -EQQQCOD(2)*QCHIRE(2)/DR
      coeff_matrix(4,9) = -1.0d0/3.0d0
      do ii = 2,num_eq_points-1
         fact_over_ei = dt_over_dr/eq_moment_of_inertia(ii)
         facta_half_dt_over_ei_dr = 0.5d0*timestep/ &
              eq_moment_of_inertia(ii)/grid_spacing
         i = 1+(ii-1)*4
! OMEGA TERMS - ADVECTIVE
         coeff_matrix(i,1) = facta_half_dt_over_ei_dr*am_advective_coeff(ii)
         coeff_matrix(i,5) = 1.0d0 + facta_half_dt_over_ei_dr* &
              (am_advective_coeff(ii)-am_advective_coeff(ii+1))
         coeff_matrix(i,9) = -facta_half_dt_over_ei_dr*am_advective_coeff(ii+1)
! OMEGA TERMS - DW/DR
         coeff_matrix(i,3) = fact_over_ei*dchi_dr_edge(ii)* &
              am_diffusive_coeff(ii)
         coeff_matrix(i,7) = -fact_over_ei*dchi_dr_edge(ii+1)* &
              am_diffusive_coeff(ii+1)
! OMEGA TERMS - D2W/DR2
!         FMI = EQQCOD(II)*QCHIRE(II)*FX1(II)
!         FPL = EQQCOD(II+1)*QCHIRE(II+1)*FX1(II+1)
         coeff_matrix(i,2) = facta_half_dt_over_ei_dr*dchi_dr_edge(ii)* &
              am_2nd_deriv_coeff(ii)*third_order_ratio_factor(ii)
!     *           -0.125D0*ECOD3(II)*DR**2/FMI
         coeff_matrix(i,6) = facta_half_dt_over_ei_dr* &
              (am_2nd_deriv_coeff(ii)*dchi_dr_edge(ii)* &
              third_order_ratio_factor(ii) &
              - am_2nd_deriv_coeff(ii+1)*dchi_dr_edge(ii+1)* &
              third_order_ratio_factor(ii+1))
!     *          +0.125D0*DR**2*(ECOD3(II)/FMI - ECOD3(II+1)/FPL)
         coeff_matrix(i,10) = -facta_half_dt_over_ei_dr* &
              am_2nd_deriv_coeff(ii+1)*dchi_dr_edge(ii+1)* &
              third_order_ratio_factor(ii+1)
!     *           +0.125D0*ECOD3(II+1)*DR**2/FPL
! OMEGA TERMS - D3W/DR3
         coeff_matrix(i,4) = fact_over_ei*am_3rd_deriv_coeff(ii)* &
              dchi_dr_edge(ii)*third_order_ratio_factor(ii)
         coeff_matrix(i,8) = -fact_over_ei*am_3rd_deriv_coeff(ii+1)* &
              dchi_dr_edge(ii+1)*third_order_ratio_factor(ii+1)
! D^2W/DR^2 COEFFICIENTS
         coeff_matrix(i+1,2) = second_deriv_geom_factor_eqgrid(ii)* &
              dchi_dr_edge(ii)/grid_spacing*third_order_ratio_factor(ii)
         coeff_matrix(i+1,5) = 1.0d0
         coeff_matrix(i+1,6) = -second_deriv_geom_factor_eqgrid(ii+1)* &
              dchi_dr_edge(ii+1)/grid_spacing*third_order_ratio_factor(ii+1)
! D W/DR COEFFICIENTS
         coeff_matrix(i+2,3) = 1.0d0/grid_spacing
         coeff_matrix(i+2,5) = 1.0d0
         coeff_matrix(i+2,7) = -1.0d0/grid_spacing
! D^3 W/DR^3 COEFFICIENTS
         coeff_matrix(i+3,3) = dchi_dr_edge(ii)* &
              third_deriv_geom_factor_eqgrid(ii)/grid_spacing
         coeff_matrix(i+3,5) = 1.0d0
         coeff_matrix(i+3,7) = -dchi_dr_edge(ii)* &
              third_deriv_geom_factor_eqgrid(ii+1)/grid_spacing
      end do
!  LAST SHELL B.C. : SAME AS FIRST SHELL B.C.
      fact_over_ei = dt_over_dr/eq_moment_of_inertia(num_eq_points)
      facta_half_dt_over_ei_dr = 0.5d0*timestep/ &
           eq_moment_of_inertia(num_eq_points)/grid_spacing
      i = 4*num_eq_points - 3
!      A(I-1,3) = 0.0D0
!      A(I-1,7) = 0.0D0
!      A(I-1,1) = -1.0D0/3.0D0
! ZERO OUT TERMS RELATED TO D2W/DR2 AT THE EDGES
      coeff_matrix(i-4,10) = 0.0d0
      coeff_matrix(5,2) = 0.0d0
      if (.not.fix_omega_at_surface) then
! OMEGA TERMS - ADVECTIVE
         coeff_matrix(i,1) = facta_half_dt_over_ei_dr* &
              am_advective_coeff(num_eq_points)
         coeff_matrix(i,5) = 1.0d0 + facta_half_dt_over_ei_dr* &
              am_advective_coeff(num_eq_points)
! OMEGA TERMS - DW/DR
         coeff_matrix(i,3) = fact_over_ei*dchi_dr_edge(num_eq_points)* &
              am_diffusive_coeff(num_eq_points)
! OMEGA TERMS - D2W/DR2
!         FMI = EQQCOD(NTOT)*QCHIRE(NTOT)*FX1(NTOT)
         coeff_matrix(i,2) = facta_half_dt_over_ei_dr* &
              dchi_dr_edge(num_eq_points)*am_2nd_deriv_coeff(num_eq_points)* &
              third_order_ratio_factor(num_eq_points)
!     *            -0.125D0*ECOD3(NTOT)*DR**2/FMI
!         A(I,6) = FACTA*QCHIRE(NTOT)*ECOD5(NTOT)*FX1(NTOT)
! OMEGA TERMS - D3W/DR3
         coeff_matrix(i,4) = fact_over_ei*dchi_dr_edge(num_eq_points)* &
              am_3rd_deriv_coeff(num_eq_points)
      else
! HOLD OMEGA FIXED
         coeff_matrix(i,5) = 1.0d0
      end if
! D^2W/DR^2 COEFFICIENTS: ASSUME DW/DR = 0 ABOVE BOUNDARY
!      A(I+1,2) = EQQCOD(NTOT)*QCHIRC(NTOT)/DR*FX1(NTOT)
      coeff_matrix(i+1,5) = 1.0d0
! D W/DR COEFFICIENTS AND D^3 W/DR^3 COEFFICIENTS NOT USED,
! SINCE THEY ARE TREATED AS ZERO ABOVE THE UNSTABLE REGION
! STORE COEFFICIENT MATRIX TO CHECK ON THE MATRIX INVERSION
      do i = 1,10
         do j = 1,num_equations
            coeff_matrix_orig(j,i) = coeff_matrix(j,i)
         end do
      end do
! NOW DECOMPOSE THE BAND MATRIX.  SRS. ARE FROM NUMERICAL
! RECIPES.
! SUBDIAGONAL ROWS
!       M1 = 4  ! KC 2025-05-31
! SUPERDIAGONAL ROWS
!       M2 = 5  ! KC 2025-05-31
! NM = TOTAL NUMBER OF ELEMENTS, 4 PER SHELL
!      WRITE(*,909)(ECOD3(I),ECOD4(I),ECOD5(I),ECOD6(I),
!     *             EV0(I),EV1A(I),EV1B(I),EV2A(I),
!     *             EV2B(I),I=1,NTOT)
! 909  FORMAT(1P9E12.3)
!      WRITE(*,910)((A(J,I),I=1,10),J=1,NM)
!  910  FORMAT(1P10E12.3)
!       CALL BANDW(A,NM,M1,M2,B)  ! KC 2025-05-31
      call bandw(coeff_matrix,num_equations,rhs)
! CHECK ON MATRIX INVERSION
      do i =1,num_equations
         residual_check(i) = 0.0d0
      end do
      do i = 5,10
         residual_check(1) = residual_check(1) + &
              coeff_matrix_orig(1,i)*rhs(i-4)
      end do
      do i = 4,10
         residual_check(2) = residual_check(2) + &
              coeff_matrix_orig(2,i)*rhs(i-3)
      end do
      do i = 3,10
         residual_check(3) = residual_check(3) + &
              coeff_matrix_orig(3,i)*rhs(i-2)
      end do
      do i = 2,10
         residual_check(4) = residual_check(4) + &
              coeff_matrix_orig(4,i)*rhs(i-1)
      end do
      do j = 5,num_equations-5
         do i = 1,10
            residual_check(j) = residual_check(j) + &
                 coeff_matrix_orig(j,i)*rhs(i-5+j)
         end do
      end do
      do i = 1,9
         residual_check(num_equations-4) = residual_check(num_equations-4)+ &
              coeff_matrix_orig(num_equations-4,i)*rhs(i+num_equations-9)
      end do
      do i = 1,8
         residual_check(num_equations-3) = residual_check(num_equations-3)+ &
              coeff_matrix_orig(num_equations-3,i)*rhs(i+num_equations-8)
      end do
      do i = 1,7
         residual_check(num_equations-2) = residual_check(num_equations-2)+ &
              coeff_matrix_orig(num_equations-2,i)*rhs(i+num_equations-7)
      end do
      do i = 1,6
         residual_check(num_equations-1) = residual_check(num_equations-1)+ &
              coeff_matrix_orig(num_equations-1,i)*rhs(i+num_equations-6)
      end do
      do i = 1,5
         residual_check(num_equations) = residual_check(num_equations)+ &
              coeff_matrix_orig(num_equations,i)*rhs(i+num_equations-5)
      end do
! CONVERT RUN OF OMEGA TO DELTA OMEGA
      sum_delta_angular_momentum = 0.0d0
      total_delta_angular_momentum_alt = 0.0d0
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
      max_omega_change_medium_iter_zone = 1
      do i = 1,num_eq_points
         ii = 1+4*(i-1)
         eq_delta_angular_momentum(i) = (rhs(ii)-eq_omega(i))* &
              eq_moment_of_inertia(i)
         total_delta_angular_momentum_alt = total_delta_angular_momentum_alt &
              + (rhs(ii)-omega_substep_start(i))*eq_moment_of_inertia(i)
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
            max_omega_change_medium_iter_zone = i
            max_omega_change_medium_iter = omega_change_medium_iter
         end if
      end do
! DAMP OUT LARGE (SPURIOUS) CHANGES IN OMEGA
      if (abs(max_omega_change).gt.0.1d0) then
         write(*,*) substep_idx,theta_iter_idx,coeff_iter_idx
         write(*,*) 'CORRECTIONS TOO LARGE ',max_omega_change, &
              max_omega_change_zone
         diffusion_converged = .false.
!      WRITE(*,920)((AA(I,J),J=1,10),B(I),C(I),D(I),I=1,NM)
!  920  FORMAT(1P13E10.3)
         goto 950
      else if (abs(max_omega_change).lt.1.0d-2) then
         damping_factor = 1.0d0
      else
!         FX = 0.01D0/ABS(DWMAX)
         damping_factor = 1.0d0
      end if
      total_angular_momentum_updated = 0.0d0
      do i = 1,num_eq_points
         ii = 1+4*(i-1)
         omega_working(i) = omega_working(i)+damping_factor* &
              (rhs(ii)-omega_working(i))
         eq_delta_angular_momentum(i) = eq_delta_angular_momentum(i)* &
              damping_factor
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
!      WRITE(*,912)FX,DWMAX,SUMDJ
! 912  FORMAT(1P2E12.3)
!      WRITE(*,913)(EW(I),EWPREV(I),DJ(I),I=1,NTOT)
! 913  FORMAT(1P6E12.3)
      max_omega_change_history(coeff_iter_idx) = max_omega_change
      max_omega_change_zone_history(coeff_iter_idx) = max_omega_change_zone
! DETERMINE IF RUN HAS CONVERGED
      if (abs(max_omega_change).lt.convergence_tolerance) then
         diffusion_converged = .true.
         goto 900
      end if
      end do
 900  continue
      write(*,*) substep_idx,theta_iter_idx,coeff_iter_idx
      write(*,100) max_omega_change_total,max_omega_change_total_zone, &
           (max_omega_change_history(j),max_omega_change_zone_history(j), &
           j=1,coeff_iter_idx)
 100  format(' MAX D(OMEGA)/OMEGA',1pe12.3,' AT PT.',i5, &
         ' BY ITERATION'/5(1x,e11.3,i4))
      if (theta_iter_idx.le.2) then
         do k = 1,num_eq_points
            omega_prev_medium_iter(k) = omega_working(k)
            omega_prev_medium_iter_avg(k) = omega_working(k)
         end do
      else
         do k = 1,num_eq_points
            omega_prev_medium_iter(k) = omega_working(k)
            omega_prev_medium_iter_avg(k) = 0.5d0*omega_prev_medium_iter(k)+ &
                 0.5d0*omega_working(k)
         end do
      end if
      if (abs(max_omega_change_medium_iter).le.convergence_tolerance .and. &
           theta_iter_idx.ge.2) then
         goto 9999
      else
         diffusion_converged = .false.
      end if
      end do
 9999 continue
      do k = 1,num_eq_points
         omega_substep_start(k) = omega_working(k)
      end do
      total_angular_momentum_start = total_angular_momentum_start+ &
           wind_loss_delta
      do k = 2,num_eq_points
         omega_mid_prev(k) = 0.5d0*(omega_working(k)+omega_working(k-1))
         domega_dr_prev(k) = dchi_dr_edge(k)*(omega_working(k)- &
              omega_working(k-1))/grid_spacing
      end do
      end do
 950  continue
      if (.not.diffusion_converged) then
         timestep_cut_count = timestep_cut_count + 1
         num_substeps = 2*num_substeps
         if (timestep_cut_count.gt.5) then
            write(*,*) 'TIMESTEP CUT MORE THAN 5 TIMES - RUN STOPPED'
            stop 111
         else
            write(*,*) 'TIMESTEP CUT #',timestep_cut_count,' IN DADCOEFT'
         end if
         goto 5
      end if
      write(*,914) rhs(1),eq_omega(1),rhs(2),rhs(3),rhs(4)
 914  format(1p5e12.3)
      do i = 2,num_eq_points
         ii = 1 + 4*(i-1)
         omega_mid_new = 0.5d0*(rhs(ii)+rhs(ii-4))
         domega_dr_mid_new = 0.5d0*(rhs(ii+1)+rhs(ii-3))
         advective_velocity_term = am_advective_coeff(i)/geometric_factor(i)
         domega_dr_velocity_term = am_diffusive_coeff(i)/geometric_factor(i)* &
              rhs(ii-2)*dchi_dr_edge(i)/omega_mid_new
         domega_dr_velocity_term_alt = 0.2d0*difad_velocity_scale* &
              eq_velocity_coeff0(i)*omega_mid_new*dchi_dr_edge(i)* &
              (eq_velocity_coeff2a(i)+eq_velocity_coeff2b(i))*rhs(ii-2)
         if (am_diffusive_coeff(i).le.0.0d0) domega_dr_velocity_term_alt = 0.0d0
         curvature_velocity_term = am_2nd_deriv_coeff(i)/geometric_factor(i)* &
              dchi_dr_edge(i)*domega_dr_mid_new/omega_mid_new
         third_deriv_velocity_term = am_3rd_deriv_coeff(i)/geometric_factor(i)* &
              rhs(ii-1)*dchi_dr_edge(i)/omega_mid_new
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
         mixing_velocity_estimate(i) = 5.0d0*abs(mixing_velocity_scale* &
              es_mixing_scale*total_velocity_alt*equatorial_radius(i))
         mixing_diffusion_raw = (mixing_velocity_estimate(i)+ &
              mixing_velocity_scale*secular_shear_mixing_scale* &
              shear_diffusion_coeff_eqgrid(i)+mixing_velocity_scale* &
              gsf_mixing_scale*gsf_diffusion_coeff_eqgrid(i))
         eq_mixing_diffusion_coeff(i) = mixing_diffusion_raw* &
              mixing_geometric_factor(i)
!         WRITE(*,1111)I,VESN(I),DCMIX,ECOD2(I),REQ(I)
 911  format(1p10e12.3)
!  1111 FORMAT(I5,1P4E12.3)
      end do
      timestep = full_timestep
      wind_loss_delta_full_step = wind_loss_delta/substep_frac
      write(*,*) sum_delta_angular_momentum,wind_loss_delta_full_step
!      IF(.NOT.LOKAD)THEN
!         SUMDJ = 0.0D0
!         DO I = 1,NTOT
!            DJ(I) = (EWPREV(I)-EW(I))*EI(I)
!            SUMDJ = SUMDJ + DJ(I)
!         END DO
!      ENDIF
!      LOKAD = .TRUE.
      return
end subroutine dadcoeft
