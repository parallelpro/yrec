!----------------------------------------------------------------------
! original_settling
!----------------------------------------------------------------------
! The original settling routines (use_new_diffusion_routines false):
! gravitational_settling, and its setup of the diffusion coefficients.
module original_settling
      use settling_grid, only: model_to_equal, equal_to_model
      use diffusion_scheme, only: lax_wendroff_step1, lax_wendroff_step2, implicit_diffusion_coeffs
      use thoul_solver, only: thoul_diffusion
      implicit none
      private
      public :: gravitational_settling

contains

! ---------------------------------------------------------------
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
! FORMULISM OF BAHCALL AND LOEB 1989.  METAL SETTLING (MHP 3/94) IS
! CARRIED ALONGSIDE IN star%metal_abundance_change WHEN
! star%job%use_diffusion_z IS SET.
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
! del_grad - DEL (=DLNT/DLNP) per zone (2026 W2: an explicit argument;
!     mix passes star%gradT, rotmix passes mix_scr%delm, which it used
!     to copy into star%gradT around this call)
! convective_flag - FLAG T/F FOR CONVECTION
! num_zones - NUMBER OF MODEL POINTS
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
subroutine gravitational_settling(timestep, composition, dlnp_dr, log_radius, log_density, &
     mass_grams, log_temperature, del_grad, convective_flag, num_zones, total_mass, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use numerics_lib
      implicit none

      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json), &
           log_temperature(json), del_grad(json)
      double precision, intent(inout) :: mass_grams(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      integer, intent(out) :: ierr
! --- locals ---
! Names below are chosen to match the dummy-argument names of the
! callees (gravitational_settling_setup/model_to_equal/lax_wendroff_step1/
! lax_wendroff_step2/implicit_diffusion_coeffs/solve_tridiagonal/
! equal_to_model) at each corresponding call site.
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
      double precision :: metal_x_orig(json)
      integer :: zone_begin, zone_end
      logical :: settling_skipped_flag
      integer :: eq_idx
      logical :: use_generic_diffusion_vectors
      double precision :: grid_spacing
      integer :: num_equal_points
      double precision :: alpha_prefactor
      integer :: iter_count
      double precision :: max_delta_x
      integer :: max_delta_x_zone
      double precision :: delta_x_local

      ierr = 0

      call gravitational_settling_setup(timestep,dlnp_dr,log_radius,log_density,mass_grams, &
           log_temperature,del_grad,convective_flag,num_zones,total_mass, &
           diffusion_coeff1,diffusion_coeff2,composition,radius_bl, &
           temperature_bl,zone_begin,zone_end,settling_skipped_flag, &
           diffusion_coeff1_dx,diffusion_coeff2_dx, ierr)
      if (ierr /= 0) return
!
! SKIP SETTLING FOR FULLY CONVECTIVE (OR H/HE-EXHAUSTED) MODELS.
      if(settling_skipped_flag) return
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
            metal_x_orig(eq_idx) = star%metal_abundance_change(eq_idx)
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
      call lax_wendroff_step1(timestep,equal_diffusion_coeff1,equal_mass, &
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
            rot_scr%metal_diffusion_coeff1_mid(eq_idx)=rot_scr%metal_diffusion_coeff1_mid(eq_idx)+ &
                 rot_scr%metal_abundance_change_mid(eq_idx)*rot_scr%eq_metal_diffusion_coeff1_mid(eq_idx)
         end do
      endif
!
! USING THE NEW COEFFICIENTS, SOLVE FOR THE NEW RUN OF HYDROGEN ABUNDANCES.
! NOTE : THE SR ACTUALLY RETURNS THE *CHANGE* IN THE ABUNDANCE AS A FUNCTION
! OF RADIUS, WHICH IS APPLIED TO THE ORIGINAL RUN OF X, RATHER THAN
! COMPUTING A NEW RUN OF ABUNDANCE AND TRANSFORMING IT BACK.  THIS IS DONE
! TO MINIMIZE ERRORS FROM THE INTERPOLATION.
!
      call lax_wendroff_step2(timestep,equal_diffusion_coeff1_mid,equal_mass_mid, &
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
         call implicit_diffusion_coeffs(alpha,equal_diffusion_coeff2_mid, &
              equal_hydrogen_fraction_mid,equal_diffusion_coeff2_dx_mid, &
              sub_diag,diag,super_diag,num_equal_points)
!
!  SOLVE THE TRIDIAGONAL MATRIX SYSTEM FOR THE NEW RUN OF X.
!
         call solve_tridiagonal(sub_diag,diag,super_diag,hydrogen_x_prime, &
              num_equal_points,equal_hydrogen_fraction, ierr)
         if (ierr /= 0) return
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
! solver forensics (2026 run-log verbosity sweep)
         if (solver_diagnostics()) &
              write(run_log_unit,90)iter_count,max_delta_x,max_delta_x_zone
 90      format(1x,'ITERATION ',i3,' DXMAX ',1pe10.2,' IMAX ',i4)
!  EXIT ITERATION LOOP IF SYSTEM HAS CONVERGED.
         if(max_delta_x.lt.star%ctrl%settling_tolerance)exit
      end do
      if (iter_count > star%ctrl%settling_num_iterations) then
      write(terminal_unit,110)star%ctrl%settling_tolerance,star%ctrl%settling_num_iterations,max_delta_x, &
           max_delta_x_zone
      write(run_log_unit,110)star%ctrl%settling_tolerance,star%ctrl%settling_num_iterations, &
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
         do eq_idx = 1,num_equal_points
            star%metal_abundance_change(eq_idx) = star%metal_abundance_change(eq_idx) - &
                 metal_x_orig(eq_idx)
         end do
      endif
! TRANSFORM BACK TO ORIGINAL MODEL GRID; UPDATE HELIUM ARRAY USING
! X+Y+Z=1.  PRINT DIAGNOSTIC OUTPUT.

      call equal_to_model(timestep,equal_radius,equal_hydrogen_fraction, &
           zone_begin,zone_end,num_equal_points,composition,dlnp_dr, &
           radius_bl,mass_grams,num_zones,total_mass)

      return
end subroutine gravitational_settling

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! setup_grsett
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setup_grsett.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! Prepares a model for the gravitational-settling solve done by
! gravitational_settling.f90: locates the diffusion zone boundaries (excluding any
! convective core/envelope and any hydrogen- or helium-exhausted
! region), converts the model to Bahcall & Loeb (1990) natural units,
! and computes the helium (and, if enabled, iron) settling diffusion
! coefficients diffusion_coeff1/diffusion_coeff2 and their partial
! derivatives with respect to X (diffusion_coeff1_dx/diffusion_
! coeff2_dx) at every zone, using either the simple analytic fit
! (LTHOUL false) or the full Thoul et al. (1994) coefficients (via
! thoul_diffusion, LTHOUL true) with a choice of Coulomb-logarithm prescription
! selected by coulomb_log_choice.
subroutine gravitational_settling_setup(timestep_seconds, dlnp_dr, log_radius, &
     log_density, mass_grams, log_temperature, del_grad, convective_flag, &
     num_zones, total_mass, diffusion_coeff1, diffusion_coeff2, &
     composition, radius_bl, temperature_bl, zone_begin, zone_end, &
     settling_skipped_flag, diffusion_coeff1_dx, diffusion_coeff2_dx, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star, json, i_h1, i_he4, i_metals, i_he3
      use species_table, only: thoul_h1, thoul_he4, thoul_fe, thoul_electron, &
           thoul_col_h, thoul_col_he, thoul_col_metal
      use bahcall_loeb_units, only: set_bahcall_loeb_scales
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(inout) :: timestep_seconds
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json), &
           log_temperature(json)
! del_grad: DEL (=DLNT/DLNP) per zone, from the caller (was star%gradT)
      double precision, intent(in) :: del_grad(json)
      double precision, intent(inout) :: mass_grams(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      double precision, intent(out) :: diffusion_coeff1(json), &
           diffusion_coeff2(json)
      double precision, intent(in) :: composition(15,json)
      double precision, intent(out) :: radius_bl(json), temperature_bl(json)
      integer, intent(out) :: zone_begin, zone_end
      logical, intent(out) :: settling_skipped_flag
      integer, intent(out) :: ierr
      double precision, intent(out) :: diffusion_coeff1_dx(json), &
           diffusion_coeff2_dx(json)
! MHP 8/94 ADDED ATOMIC WEIGHTS AND CHARGES FOR H,HE,FE,ELECTRONS -
! NEEDED FOR FULL THOUL COEFFICIENTS (VALUES FROM species_table)
      double precision :: atomic_weight(4), atomic_charge(4), &
           settling_ap(4), settling_at(4), settling_ac(4,4), &
           coulomb_log(4,4), species_mass_fraction(4)
!     additional variables to be declared for variable ln lambda
!     intermediate variables:
! mass_weighted_conc (originally AC): concentration-weighted atomic
! weight sum used to get the electron density (coulomb_log_choice==4).
      double precision :: charge_to_mass_sum, mass_weighted_conc, &
           ion_number_density, mean_charge_sq, coulomb_xij, &
           electron_number_density, interion_distance, debye_length, &
           coulomb_lambda, concentration(4)
      double precision :: ln_lambda
      data atomic_weight/thoul_h1%weight,thoul_he4%weight,thoul_fe%weight, &
           thoul_electron%weight/
      data atomic_charge/thoul_h1%charge,thoul_he4%charge,thoul_fe%charge, &
           thoul_electron%charge/
      integer :: num_species
      data num_species/4/
! --- locals ---
      integer :: zone_idx
      double precision :: hydrogen_fraction, z_plus_he3_fraction, &
           metal_fraction, hydrogen_fraction_sq, &
           hydrogen_metal_product, thoul_denominator
! helium_ah_coeff (also originally AC): empirical Thoul-fit helium-
! settling amplitude coefficient feeding diffusion_coeff2.
      double precision :: settling_prefactor, iron_settling_ah, helium_ah_coeff
      double precision :: settling_coeff_p, settling_coeff_t, dap_dx, dat_dx, dac_dx
      double precision :: rho_local, temp_local
      integer :: species_a_idx, species_b_idx

!     OUTPUT VARIABLES :
!
!     radius_bl - VECTOR OF UNLOGGED RADII IN BAHCALL AND LOEB UNITS.
!     temperature_bl - VECTOR OF UNLOGGED TEMPERATURES IN BAHCALL AND
!          LOEB UNITS.
!     zone_begin - FIRST ZONE FOR DIFFUSION PURPOSES (EITHER THE FIRST
!          MODEL POINT OR THE OUTERMOST POINT OF A CENTRAL CONVECTION
!          ZONE).
!     zone_end - LAST ZONE FOR DIFFUSION PURPOSES (EITHER THE LAST
!          MODEL POINT OR THE INNERMOST POINT OF A SURFACE CONVECTION
!          ZONE).
!     THE VECTORS dlnp_dr AND mass_grams, AND THE SCALARS
!     timestep_seconds AND total_mass, ARE ALSO CONVERTED TO BAHCALL
!     AND LOEB UNITS.
!     CONSTANTS DEFINED :
!     ln10 = CONVERSION FACTOR FROM LN TO LOG10
!     settling_skipped_flag=T IF SETTLING IS SKIPPED THIS MODEL (FULLY
!     CONVECTIVE, HYDROGEN-EXHAUSTED OR HELIUM-EXHAUSTED).
      ierr = 0
      settling_skipped_flag=.false.
!     CHECK FOR CONVECTIVE CORE.
      if(convective_flag(1))then
         do zone_idx=2,num_zones
            if(.not.convective_flag(zone_idx))exit
         end do
         if (zone_idx > num_zones) then
!        DIFFUSION NOT COMPUTED FOR FULLY CONVECTIVE MODELS.
         settling_skipped_flag=.true.
! print once per suspension; every model only under
! report_solver_diagnostics (2026 run-log verbosity sweep)
         if (solver_diagnostics() .or. &
              .not. star%settling_suspended_reported) then
            write(run_log_unit,15)
   15       format(1x,' FULLY CONVECTIVE MODEL - NO SETTLING')
            star%settling_suspended_reported = .true.
         end if
         return
         end if
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_begin = zone_idx-1
      else
         zone_begin = 1
      endif
! MHP 6/90 CHECK FOR HYDROGEN-EXHAUSTED CORE.
      do zone_idx = zone_begin,num_zones
         if(composition(i_h1,zone_idx).gt.star%ctrl%hydrogen_diffusion_floor)exit
      end do
      if (zone_idx > num_zones) then
!     HYDROGEN-FREE MODEL - EXIT.
! print once per suspension; every model only under
! report_solver_diagnostics (2026 run-log verbosity sweep)
      if (solver_diagnostics() .or. &
           .not. star%settling_suspended_reported) then
         write(run_log_unit,16)star%ctrl%hydrogen_diffusion_floor
   16    format(1x,'X BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
         star%settling_suspended_reported = .true.
      end if
      settling_skipped_flag = .true.
      return
      end if
      zone_begin = zone_idx
!     CHECK FOR CONVECTIVE ENVELOPE.
      if(convective_flag(num_zones))then
         do zone_idx=num_zones-1,2,-1
            if(.not.convective_flag(zone_idx))exit
         end do
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_end = zone_idx+1
      else
         zone_end = num_zones
      endif
!     CHECK FOR HELIUM-EXHAUSTED SURFACE.
!     OUTER POINT IS SET WHEREVER Y>YMIN.
      do zone_idx=zone_end,1,-1
         if(composition(i_he4,zone_idx).gt.star%ctrl%helium_diffusion_min) exit
      end do
      if (zone_idx < 1) then
!     HELIUM-EXHAUSTED MODEL - EXIT.
! print once per suspension; every model only under
! report_solver_diagnostics (2026 run-log verbosity sweep)
      if (solver_diagnostics() .or. &
           .not. star%settling_suspended_reported) then
         write(run_log_unit,17)star%ctrl%helium_diffusion_min
   17    format(1x,'Y BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
         star%settling_suspended_reported = .true.
      end if
      settling_skipped_flag = .true.
      return
      end if
      zone_end = zone_idx
! all suspension checks passed: settling proceeds this model
      if (star%settling_suspended_reported) then
         write(run_log_unit,916)
  916    format(1x,' SETTLING RESUMED')
         star%settling_suspended_reported = .false.
      end if
!     SET THE star%bl_*_scale CONVERSION FACTORS.
      call set_bahcall_loeb_scales()
!     CONVERT LOG(RADIUS) AND LOG(TEMPERATURE) TO NATURAL UNITS.
!     ALSO CONVERT NATURAL UNITS TO BAHCALL AND LOEB UNITS.
      do zone_idx=1,num_zones
         radius_bl(zone_idx)=exp(ln10*log_radius(zone_idx))*star%bl_radius_scale
         temperature_bl(zone_idx)=exp(ln10*log_temperature(zone_idx))*star%bl_temp_scale
         mass_grams(zone_idx)=mass_grams(zone_idx)*star%bl_mass_scale
         dlnp_dr(zone_idx)=dlnp_dr(zone_idx)/star%bl_radius_scale
      end do
      timestep_seconds=timestep_seconds/star%bl_time_scale
      total_mass=total_mass*star%bl_mass_scale
!     SET UP DIFFUSION COEFFICIENTS.
!     MODIFIED BY BC MAY/90 -- VALID FOR ALL X WITH VARIABLE LN LAMBDA
!     GENERAL EQUATION IS
!     DX/DT = (D/DR(D1)+D/DR(D2 DX/DR))/(RHO*R**2), WHERE
!     D1 = R**2/LN LAMBDA  * X  * T**5/2 * (DLNP/DR) * (1-X) *
!          [5/4 + DEL*6*(X-0.32)/(5.4+6.3X-4.5X**2)]
!     D2 = R**2/LN LAMBDA * T**5/2 * (X+3)/(X+1)/(3+5X)
      do zone_idx = 1,num_zones
         hydrogen_fraction = composition(i_h1,zone_idx)
! MHP 10/02 INITIALIZED X - WAS NOT DONE PRIOR TO USAGE IN SHELL 1
         if(star%ctrl%coulomb_log_choice.eq.2)then
!           Noerdlinger's formula (1977 A&A 57,407) for LN LAMBDA:
!           Ln Lambda = -19.7 - ln[4/(3*X + 1)]/2 - ln(rho)/2 + 1.5*ln(T)
            ln_lambda=-1.97d1 - 0.5d0*log(1.0d0/(0.75d0*hydrogen_fraction+0.25d0)) &
                      -0.5d0*log_density(zone_idx)*ln10 + 1.5d0*log_temperature(zone_idx)*ln10
         else if(star%ctrl%coulomb_log_choice.eq.3)then
!           Loeb's formula (1989 Phys. Rev. D 39, 1009) for LN LAMBDA
!           ln lambda = -19.105747 - ln(rho)/2 + 1.5*ln(T)
            ln_lambda=-1.9105747d1-0.50d0*log_density(zone_idx)*ln10+1.5d0*log_temperature(zone_idx)*ln10
         else
            ln_lambda = 2.2d0
         end if
!
         settling_prefactor=star%ctrl%fgry*radius_bl(zone_idx)**2*pow(temperature_bl(zone_idx), 2.5d0)/ln_lambda
         z_plus_he3_fraction = composition(i_metals,zone_idx) + composition(i_he3,zone_idx)
         metal_fraction = composition(i_metals,zone_idx)
         hydrogen_fraction_sq = hydrogen_fraction*hydrogen_fraction
         if(.not.star%ctrl%use_thoul_diffusion)then
            hydrogen_metal_product = hydrogen_fraction*z_plus_he3_fraction
            thoul_denominator=5.4d0+6.3d0*hydrogen_fraction-4.5d0*hydrogen_fraction_sq
            diffusion_coeff1(zone_idx)=settling_prefactor*dlnp_dr(zone_idx)* &
                 (hydrogen_fraction - hydrogen_fraction_sq - hydrogen_metal_product)*(1.25d0+ &
                 del_grad(zone_idx)*6.0d0*(hydrogen_fraction+0.32d0)/thoul_denominator)
            diffusion_coeff2(zone_idx)=settling_prefactor*(hydrogen_fraction+3.0d0)/ &
                 (5.0d0*hydrogen_fraction_sq + 8.0d0*hydrogen_fraction + 3.0d0)
            diffusion_coeff1_dx(zone_idx)=settling_prefactor*dlnp_dr(zone_idx)* &
                 ( (1.0d0-2.0d0*hydrogen_fraction-z_plus_he3_fraction)*(1.25d0+ &
                 (6.0d0*del_grad(zone_idx)*(hydrogen_fraction+0.32d0))/thoul_denominator)+ &
                 (hydrogen_fraction-hydrogen_fraction_sq-hydrogen_metal_product)*6.0d0* &
                 del_grad(zone_idx)*(3.384d0+2.88d0*hydrogen_fraction+4.5d0*hydrogen_fraction_sq)/ &
                 thoul_denominator**2 )
            diffusion_coeff2_dx(zone_idx)=-settling_prefactor*(5.0d0*hydrogen_fraction_sq + &
                 3.0d1*hydrogen_fraction + 2.1d1)/ &
                 (5.0d0*hydrogen_fraction_sq + 8.0d0*hydrogen_fraction + 3.0d0)**2
         else
            species_mass_fraction(thoul_col_h) = composition(i_h1,zone_idx)
            species_mass_fraction(thoul_col_he) = composition(i_he4,zone_idx)
            species_mass_fraction(thoul_col_metal) = composition(i_metals,zone_idx)
            if(.not.star%ctrl%use_thoul_fit)then
               if(star%ctrl%coulomb_log_choice.eq.4)then
                  rho_local = exp(ln10*log_density(zone_idx))
                  temp_local = exp(ln10*log_temperature(zone_idx))
!                 calculate concentrations from mass fractions:
                  charge_to_mass_sum=0.d0
                  do species_a_idx=1,num_species-1
                   charge_to_mass_sum=charge_to_mass_sum+atomic_charge(species_a_idx)* &
                        species_mass_fraction(species_a_idx)/atomic_weight(species_a_idx)
                  enddo
                  do species_a_idx=1,num_species-1
                     concentration(species_a_idx)=species_mass_fraction(species_a_idx)/ &
                          (atomic_weight(species_a_idx)*charge_to_mass_sum)
                  enddo
                  concentration(num_species)=1.d0
!                 calculate density of electrons (NE) from mass density (RHO):
                  mass_weighted_conc=0.d0
                  do species_a_idx=1,num_species
                   mass_weighted_conc=mass_weighted_conc+atomic_weight(species_a_idx)*concentration(species_a_idx)
                  enddo
                  electron_number_density=rho_local/(m_proton_cgs*mass_weighted_conc)
!                 calculate interionic distance (AO):
                  ion_number_density=0.d0
                  do species_a_idx=1,num_species-1
                     ion_number_density=ion_number_density+concentration(species_a_idx)*electron_number_density
                  enddo
                  interion_distance=pow((0.23873d0/ion_number_density), cc13)
!                 calculate Debye length (LAMBDAD):
                  mean_charge_sq=0.d0
                  do species_a_idx=1,num_species
                   mean_charge_sq=mean_charge_sq+concentration(species_a_idx)*atomic_charge(species_a_idx)**2
                  enddo
                  debye_length=6.9010d0*sqrt(temp_local/(electron_number_density*mean_charge_sq))
!                 calculate LAMBDA to use in Coulomb logarithm:
                  coulomb_lambda=max(debye_length,interion_distance)
!                 calculate Coulomb logarithms:
                  do species_a_idx=1,num_species
                   do species_b_idx=1,num_species
                      coulomb_xij=2.3939d3*temp_local*coulomb_lambda/ &
                           abs(atomic_charge(species_a_idx)*atomic_charge(species_b_idx))
                      coulomb_log(species_a_idx,species_b_idx)=0.81245d0 &
                             *log(1.d0+0.18769d0*pow(coulomb_xij, 1.2d0))
                   enddo
                  enddo
               else
                  do species_a_idx = 1,num_species
                     do species_b_idx = 1,num_species
                        coulomb_log(species_a_idx,species_b_idx) = ln_lambda
                     end do
                  end do
               endif
               call thoul_diffusion(num_species,atomic_weight,atomic_charge, &
                    species_mass_fraction,coulomb_log,settling_ap,settling_at,settling_ac, &
                    ierr)
               if (ierr /= 0) return
               settling_coeff_p = -settling_ap(thoul_col_h)
               settling_coeff_t = -del_grad(zone_idx)*settling_at(thoul_col_h)
            else
               settling_coeff_p = 1.58d0 - 2.42d0*hydrogen_fraction + 0.844d0*hydrogen_fraction_sq
               settling_coeff_t = del_grad(zone_idx)*(1.90d0 - 2.69d0*hydrogen_fraction + 0.805d0*hydrogen_fraction_sq)
            endif
            helium_ah_coeff = 1.15d0 - 1.42d0*hydrogen_fraction + 0.647d0*hydrogen_fraction_sq
            dap_dx = -2.42d0 + 1.688d0*hydrogen_fraction
            dat_dx = del_grad(zone_idx)*(-2.69d0 + 1.61d0*hydrogen_fraction)
            dac_dx = -1.42d0 + 1.294d0*hydrogen_fraction
!CFD 10/09 Mimic Mixing to reduce settling.
            diffusion_coeff1(zone_idx) = star%ctrl%constant_mixing_coeff*settling_prefactor*dlnp_dr(zone_idx)* &
                 hydrogen_fraction*(settling_coeff_p+settling_coeff_t)
            diffusion_coeff2(zone_idx) = star%ctrl%constant_mixing_coeff*settling_prefactor*helium_ah_coeff
            diffusion_coeff1_dx(zone_idx) = settling_prefactor*dlnp_dr(zone_idx)* &
                 (settling_coeff_p+settling_coeff_t+hydrogen_fraction*(dap_dx+dat_dx))
            diffusion_coeff2_dx(zone_idx) = settling_prefactor*dac_dx
         endif
!        METAL DIFFUSION, USING THE THOUL ET AL. COEFFICIENTS FOR FULLY
!        IONIZED IRON.
         if(star%job%use_diffusion_z)then
            settling_prefactor=star%job%fgrz*radius_bl(zone_idx)**2*pow(temperature_bl(zone_idx), 2.5d0)/ln_lambda
            if(star%ctrl%use_thoul_diffusion)then
               if(star%ctrl%use_thoul_fit)then
                  settling_coeff_p = -0.157d0 -0.511d0*hydrogen_fraction + 0.389d0*hydrogen_fraction_sq
                  settling_coeff_t = del_grad(zone_idx)*(-1.36d0 - 1.42d0*hydrogen_fraction + &
                       0.549d0*hydrogen_fraction_sq)
               else
                  settling_coeff_p = -settling_ap(thoul_col_metal)
                  settling_coeff_t = -del_grad(zone_idx)*settling_at(thoul_col_metal)
               endif
               iron_settling_ah = -0.0375d0 -0.193d0*hydrogen_fraction + 0.107d0*hydrogen_fraction_sq
!CFD 10/09 Mimic Mixing to reduce settling (constant_mixing_coeff)
!         and add the uncertainties of differential mixing (constant_settling_reduction).
               rot_scr%src_grid_metal_diffusion_coeff1(zone_idx) = star%ctrl%constant_settling_reduction*star%ctrl%constant_mixing_coeff* &
                    settling_prefactor*dlnp_dr(zone_idx)*metal_fraction*(settling_coeff_p+settling_coeff_t)
!              POSITIVE DIFFUSION COEFFICIENTS NEEDED!
               rot_scr%src_grid_metal_diffusion_coeff2(zone_idx) = star%ctrl%constant_mixing_coeff*abs(settling_prefactor*iron_settling_ah)
               rot_scr%src_grid_metal_diffusion_coeff1_dz(zone_idx) = star%ctrl%constant_mixing_coeff*settling_prefactor* &
                    dlnp_dr(zone_idx)*(settling_coeff_p+settling_coeff_t)
               rot_scr%src_grid_metal_diffusion_coeff2_dz(zone_idx) = 0.0d0
            endif
         endif
!
! If using Noerdlinger's formula for LN LAMBDA, have a new term
! in D(D1)/DX
         if(star%ctrl%coulomb_log_choice.eq.2)  then
           diffusion_coeff1_dx(zone_idx)=diffusion_coeff1_dx(zone_idx) + diffusion_coeff1(zone_idx)*1.5d0/ &
           (ln_lambda*(3.0d0*hydrogen_fraction+1.0d0))
       end if
      end do
      return
end subroutine gravitational_settling_setup

end module original_settling
