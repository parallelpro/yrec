!----------------------------------------------------------------------
! grsett
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original grsett.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! GRSETT SOLVES FOR THE GRAVITATIONAL SETTLING OF HELIUM USING THE
! FORMULISM OF BAHCALL AND LOEB 1989.  SETTLING OF METALS WILL BE ADDED
! LATER.
! GRSETT PERFORMS THE FOLLOWING SERIES OF OPERATIONS :
! 1)  TRANSFORMS TO AN EQUALLY SPACED GRID IN RADIUS (AND MANIPULATES
! THE MODEL VARIABLES INTO BAHCALL AND LOEB UNITS).
!
! THE DIFFUSION EQUATION HAS TWO TERMS : ONE THAT DEPENDS ON DLNP/DR*DX/DR
! AND ONE THAT DEPENDS ON D^2X/DR^2. THE DIFFUSION COEFFICIENTS ARE
! THEMSELVES FUNCTIONS OF BOTH THE THERMAL STRUCTURE AND X.
! GRSETT SOLVES THIS DIFFUSION EQUATION IN TWO STEPS :
! 2) THE FIRST TERM IS SOLVED EXPLICITLY USING THE TWO-STEP
! LAX-WENDROFF SCHEME (NUMERICAL RECIPES, PRESS ET AL. 1986,
! CAMBRIDGE UNIVERSITY PRESS,P.633)
! 3) THE SECOND TERM IS SOLVED IMPLICITLY, USING THE MODIFIED RUN OF
! ABUNDANCE FOUND IN STEP 2) AS THE 'INITIAL' RUN OF ABUNDANCE.
! THE IMPLICIT SCHEME ITERATES ON THE DIFFUSION COEFFICIENTS UNTIL
! THE SOLUTION CONVERGES TO WITHIN THE USER-SPECIFIED TOLERANCE
! settling_tolerance.
!
! 4) GRSETT THEN UPDATES THE ABUNDANCE ARRAY AND EXITS.
!
! INPUT VARIABLES :
!
! timestep - TIMESTEP (SEC)
! composition - RUN OF MASS FRACTIONS OF SPECIES. composition(1,..)=X,
!               composition(2,...)=Y, composition(3,...)=Z
! dlnp_dr - DLNP/DR
! log_radius - LOG RADIUS (CM)
! mass_grams - MASS (GM), UNLOGGED
! log_temperature - LOG TEMPERATURE (K)
! convective_flag - FLAG T/F FOR CONVECTION
! num_zones - NUMBER OF MODEL POINTS
! [COMMON] del_grad(2,...) - DEL (=DLNT/DLNP)
!
! OUTPUT VARIABLES :
!
! NEW RUN OF composition(1,...) AND composition(2,...).
!
! CONVERT MODEL QUANTITIES TO BAHCALL AND LOEB UNITS; LOCATE
! BOUNDARIES OF CONVECTIVE CORE AND ENVELOPE IF APPLICABLE.
! COMPUTE THE DIFFUSION COEFFICIENTS (diffusion_coeff1 AND
! diffusion_coeff2) AND THEIR DERIVATIVES WITH RESPECT TO X
! (diffusion_coeff1_dx AND diffusion_coeff2_dx).
!
subroutine grsett(timestep, composition, dlnp_dr, log_radius, log_density, &
     mass_grams, log_temperature, convective_flag, num_zones, total_mass)

      use star_info_lib, only: star, json
      use luout_lib
      use const_lib
      use numerics_lib
      implicit none

      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json), &
           log_temperature(json)
      double precision, intent(inout) :: mass_grams(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass






! --- locals ---
! Names below are chosen to match the dummy-argument names of the
! callees (setup_grsett.f90/model_to_equal.f90/lax_wendrof1.f90/
! lax_wendrof2.f90/get_imp_diffco.f90/tridiag_gs.f90/
! equal_to_model.f90) at each corresponding call site.
      double precision :: equal_mass(json), equal_radius(json), &
           sub_diag(json), diag(json), super_diag(json), &
           equal_hydrogen_fraction(json)
      double precision :: radius_bl(json), temperature_bl(json)
      double precision :: equal_mass_mid(json), hydrogen_x_prev_iter(json), &
           hydrogen_x_orig(json), alpha(json), &
           equal_hydrogen_fraction_mid(json)
      double precision :: diffusion_coeff1(json), diffusion_coeff2(json), &
           diffusion_coeff1_dx(json), diffusion_coeff2_dx(json)
      double precision :: equal_diffusion_coeff1(json), &
           equal_diffusion_coeff1_mid(json), hydrogen_x_prime(json), &
           equal_diffusion_coeff2_mid(json), &
           equal_diffusion_coeff1_dx_mid(json), &
           equal_diffusion_coeff2_dx_mid(json)
      double precision :: metal_x_orig(json), metal_x_prev_iter(json)
      integer :: zone_begin, zone_end
      logical :: fully_convective_flag
      integer :: eq_idx
      logical :: use_generic_diffusion_vectors
      double precision :: grid_spacing
      integer :: num_equal_points
      double precision :: alpha_prefactor
      integer :: iter_count
      double precision :: max_delta_x
      integer :: max_delta_x_zone
      double precision :: delta_x_local
      double precision :: z_change_first, z_change_last

      call setup_grsett(timestep,dlnp_dr,log_radius,log_density,mass_grams, &
           log_temperature,convective_flag,num_zones,total_mass, &
           diffusion_coeff1,diffusion_coeff2,composition,radius_bl, &
           temperature_bl,zone_begin,zone_end,fully_convective_flag, &
           diffusion_coeff1_dx,diffusion_coeff2_dx)
!
! SKIP SETTLING FOR FULLY CONVECTIVE MODELS.
      if(fully_convective_flag) return
!  TRANSFORM TO AN EQUALLY SPACED GRID IN RADIUS.
!  NOTE : PREFIX E ALONE=VARIABLE AT EQUALLY SPACED GRID POINTS.
!         PREFIX E + SUFFIX _H= VARIABLE AT MIDPOINT BETWEEN EQUALLY
!         SPACED GRID POINTS.(BOTH ARE NEEDED FOR THE DIFFUSION TECHNIQUE).
!
      call model_to_equal(diffusion_coeff1,diffusion_coeff2,composition, &
           radius_bl,mass_grams, &
           zone_begin,zone_end,num_zones,diffusion_coeff1_dx, &
           diffusion_coeff2_dx,grid_spacing, &
           equal_diffusion_coeff1,equal_diffusion_coeff1_mid, &
           equal_diffusion_coeff2_mid,equal_mass,equal_mass_mid, &
           equal_diffusion_coeff1_dx_mid,equal_diffusion_coeff2_dx_mid, &
           equal_radius,equal_hydrogen_fraction, &
           equal_hydrogen_fraction_mid,num_equal_points)
!
! STORE ORIGINAL RUN OF HYDROGEN MASS FRACTION IN VECTOR hydrogen_x_orig.
! THE CHANGE IN HYDROGEN (equal_hydrogen_fraction - hydrogen_x_orig) IS
! INTERPOLATED BACK TO THE ORIGINAL GRID AT THE END OF THE ROUTINE; THE
! CHANGE IN ABUNDANCE (RATHER THAN THE NEW RUN OF ABUNDANCE) IS USED TO
! MINIMIZE ERRORS ARISING FROM THE INTERPOLATION.
      do eq_idx = 1,num_equal_points
         hydrogen_x_orig(eq_idx) = equal_hydrogen_fraction(eq_idx)
      end do
! MHP 3/94 METAL DIFFUSION
      if(star%job%use_diffusion_z)then
         do eq_idx = 1,num_equal_points
            metal_x_orig(eq_idx) = star%rot%metal_abundance_change(eq_idx)
         end do
      endif
!  FIRST STEP OF TWO STEP LAX-WENDROFF METHOD.   COMPUTE NEW X'S AT ZONE
!  MIDPOINTS USING THE LAX SCHEME :
! X(N+1/2,J+1/2)=1/2(X(N,J+1/2)-1/2(DT/DR)(diffusion_coeff1(N,J+1)-diffusion_coeff1(N,J))
! WHERE N IS THE TIME VARIABLE, J IS THE SPATIAL ONE, AND diffusion_coeff1
! IS THE DIFFUSION COEFFICIENT.
!
! SOLVE FOR NEW ABUNDANCES AT THE ZONE MIDPOINTS USING THE ABOVE D.C.
! WITH use_generic_diffusion_vectors = FALSE, SIGNIFYING THAT LAX_WENDROF
! SHOULD DIFFUSE Z.
!
      use_generic_diffusion_vectors = .false.
      call lax_wendrof1(timestep,equal_diffusion_coeff1,equal_mass, &
           num_equal_points,total_mass,equal_hydrogen_fraction_mid, &
           use_generic_diffusion_vectors)
!
! NOW COMPUTE NEW DIFFUSION COEFFICIENTS AT THE ZONE MIDPOINTS USING THE
! PROVISIONAL SOLUTION FOR X AT THE ZONE MIDPOINTS.
!
! GET NEW DIFFUSION COEFFICIENTS.
! ALSO DEFINE VECTOR hydrogen_x_prev_iter, WHICH IS THE HYDROGEN
! ABUNDANCE AT THE END OF THE PREVIOUS ITERATION.
!
      do eq_idx=1,num_equal_points-1
         equal_diffusion_coeff1_mid(eq_idx)=equal_diffusion_coeff1_mid(eq_idx)+ &
              equal_hydrogen_fraction_mid(eq_idx)*equal_diffusion_coeff1_dx_mid(eq_idx)
         hydrogen_x_prev_iter(eq_idx) = hydrogen_x_orig(eq_idx)
      end do
      hydrogen_x_prev_iter(num_equal_points) = hydrogen_x_orig(num_equal_points)
! MHP 3/94 METAL DIFFUSION
      if(star%job%use_diffusion_z)then
         do eq_idx=1,num_equal_points-1
            star%rot%metal_diffusion_coeff1_mid(eq_idx)=star%rot%metal_diffusion_coeff1_mid(eq_idx)+ &
                 star%rot%metal_abundance_change_mid(eq_idx)*star%rot%eq_metal_diffusion_coeff1_mid(eq_idx)
            metal_x_prev_iter(eq_idx) = metal_x_orig(eq_idx)
         end do
         metal_x_prev_iter(num_equal_points) = metal_x_orig(num_equal_points)
      endif
!
! USING THE NEW COEFFICIENTS, SOLVE FOR THE NEW RUN OF HYDROGEN ABUNDANCES.
! NOTE : THE SR ACTUALLY RETURNS THE *CHANGE* IN THE ABUNDANCE AS A FUNCTION
! OF RADIUS, WHICH IS APPLIED TO THE ORIGINAL RUN OF X, RATHER THAN
! COMPUTING A NEW RUN OF ABUNDANCE AND TRANSFORMING IT BACK.  THIS IS DONE
! TO MINIMIZE ERRORS FROM THE INTERPOLATION.
!
      call lax_wendrof2(timestep,equal_diffusion_coeff1_mid,equal_mass_mid, &
           equal_hydrogen_fraction,num_equal_points,total_mass, &
           use_generic_diffusion_vectors)
!
!  NOW IMPLICITLY SOLVE FOR THE SECOND TERM (INVOLVING THE SECOND
!  DERIVATIVE OF THE COMPOSITION GRADIENT).
!  hydrogen_x_prime IS THE HYDROGEN ABUNDANCE ONE WOULD HAVE IN THE
!  ABSENCE OF THE SECOND TERM IN THE DIFFUSION EQUATION.
      do eq_idx = 1,num_equal_points
         hydrogen_x_prime(eq_idx) = equal_hydrogen_fraction(eq_idx)
      end do
!  alpha IS THE NUMERICAL FACTOR IN FRONT OF THE DIFFUSION COEFFICIENTS.
      alpha_prefactor = c4pi*timestep/grid_spacing
      alpha(1) = alpha_prefactor/equal_mass_mid(1)
      do eq_idx = 2,num_equal_points-1
         alpha(eq_idx) = alpha_prefactor/(equal_mass_mid(eq_idx)-equal_mass_mid(eq_idx-1))
      end do
      alpha(num_equal_points) = alpha_prefactor/ &
           (total_mass-equal_mass_mid(num_equal_points-1))
!  START ITERATION LOOP FOR THE NEW RUN OF HYDROGEN.
      do iter_count=1,star%ctrl%settling_num_iterations
!  FIND CHANGE IN X AT THE ZONE MIDPOINTS, GIVEN CHANGE IN X AT
!  THE ZONE CENTERS.
         do eq_idx = 2,num_equal_points
            equal_hydrogen_fraction_mid(eq_idx) = 0.5d0*(equal_hydrogen_fraction(eq_idx)+ &
                 equal_hydrogen_fraction(eq_idx-1)-hydrogen_x_prev_iter(eq_idx)- &
                 hydrogen_x_prev_iter(eq_idx-1))
         end do
!  STORE CURRENT RUN OF HYDROGEN ABUNDANCES IN VECTOR hydrogen_x_prev_iter.
!  THE ITERATION LOOP IS COMPLETED ONCE equal_hydrogen_fraction -
!  hydrogen_x_prev_iter IS EVERYWHERE LESS THAN THE TOLERANCE
!  settling_tolerance.
         do eq_idx = 1,num_equal_points
            hydrogen_x_prev_iter(eq_idx) = equal_hydrogen_fraction(eq_idx)
         end do
!  GET NEW DIFFUSION COEFFICIENTS, TAKING INTO ACCOUNT THE CHANGE IN X
!  FROM THE PREVIOUS ITERATION.
!
         call get_imp_diffco(alpha,equal_diffusion_coeff2_mid, &
              equal_hydrogen_fraction_mid,equal_diffusion_coeff2_dx_mid, &
              sub_diag,diag,super_diag,num_equal_points)
!
!  SOLVE THE TRIDIAGONAL MATRIX SYSTEM FOR THE NEW RUN OF X.
!
         call tridiag_gs(sub_diag,diag,super_diag,hydrogen_x_prime, &
              num_equal_points,equal_hydrogen_fraction)
!
!  CHECK TO SEE IF THE CORRECTIONS TO THE HYDROGEN ABUNDANCE ARE SMALL
!  ENOUGH TO EXIT.
         max_delta_x = 0.0d0
         max_delta_x_zone = 0
         do eq_idx = 1,num_equal_points
            delta_x_local = equal_hydrogen_fraction(eq_idx)-hydrogen_x_prev_iter(eq_idx)
            if(delta_x_local.gt.max_delta_x)then
               max_delta_x = delta_x_local
               max_delta_x_zone = eq_idx
            endif
         end do
!         WRITE(IOWR,90)ITER,DXMAX,IMAX
         write(short_file_unit,90)iter_count,max_delta_x,max_delta_x_zone
 90      format(1x,'ITERATION ',i3,' DXMAX ',1pe10.2,' IMAX ',i4)
!  EXIT ITERATION LOOP IF SYSTEM HAS CONVERGED.
         if(max_delta_x.lt.star%ctrl%settling_tolerance)exit
      end do
      if (iter_count > star%ctrl%settling_num_iterations) then
      write(iowr,110)star%ctrl%settling_tolerance,star%ctrl%settling_num_iterations,max_delta_x, &
           max_delta_x_zone
      write(short_file_unit,110)star%ctrl%settling_tolerance,star%ctrl%settling_num_iterations, &
           max_delta_x,max_delta_x_zone
 110  format(1x,'GRSETT FAILED TO CONVERGE TO WITHIN ',1pe9.3,' IN ',i3, &
           'ITERATIONS'/1x,'LAST ITERATION CHANGE IN X ',1pe9.3, &
           ' IN EQUALLY SPACED SHELL ',i5)
      end if
!  FIND RUN OF CHANGES IN X.
      do eq_idx = 1,num_equal_points
         equal_hydrogen_fraction(eq_idx) = equal_hydrogen_fraction(eq_idx)-hydrogen_x_orig(eq_idx)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
      if(star%job%use_diffusion_z)then
         z_change_first=star%rot%metal_abundance_change(1)-metal_x_orig(1)
         z_change_last=star%rot%metal_abundance_change(num_equal_points)- &
              metal_x_orig(num_equal_points)
         do eq_idx = 1,num_equal_points
            star%rot%metal_abundance_change(eq_idx) = star%rot%metal_abundance_change(eq_idx) - &
                 metal_x_orig(eq_idx)
         end do
      endif
! TRANSFORM BACK TO ORIGINAL MODEL GRID; UPDATE HELIUM ARRAY USING
! X+Y+Z=1.  PRINT DIAGNOSTIC OUTPUT.

      call equal_to_model(timestep,equal_radius,equal_hydrogen_fraction, &
           zone_begin,zone_end,num_equal_points,composition,dlnp_dr, &
           radius_bl,mass_grams,temperature_bl,num_zones,total_mass)

      return
end subroutine grsett
