!----------------------------------------------------------------------
! rotational_mixing
!----------------------------------------------------------------------
! The rotational mixing of the composition: the equally spaced grid
! (composition_grid), the diffusion coefficients and the implicit solve
! (diffuse_composition), driven by diffuse_composition_driver.
module rotational_mixing
      implicit none
      private
      public :: diffuse_composition_driver

contains

! ---------------------------------------------------------------
! Dummy-argument names for the quantities passed straight through to
! composition_grid.f90 (log_density/luminosity_lsun/log_pressure/log_radius/
! log_mass/enclosed_mass/log_total_mass) match that already-converted
! callee's own parameter names, and the whole argument list matches the
! already-converted caller burn_settle_mix.f90's call to ndifcom.
!
!  DIFCOM CALCULATES THE DIFFUSION OF COMPOSITION DUE TO ANGULAR MOMENTUM
!  TRANSPORT.  THIS IS DONE BY TRANSFORMING TO AN EQUALLY SPACED GRID IN
!  RADIUS, SOLVING A DIFFUSION EQUATION, AND TRANSFORMING BACK.
!
!  INPUT VARIABLES:
!  timestep (DT) : DIFFUSION TIMESTEP(SEC).
!  diffusion_coeff (COD2) : DIFFUSION COEFFICIENTS FOR COMPOSITION
!     TRANSPORT AT THE ORIGINAL MODEL POINTS. Passed straight through
!     to mixgrid.
!  equally_spaced_diffusion_coeff (ECOD2) : DIFFUSION COEFFICIENTS FOR
!     COMPOSITION TRANSPORT AT THE EQUALLY SPACED GRID POINTS.
!  equally_spaced_mass (EM) : MASSES OF THE EQUALLY SPACED GRID
!     POINTS(GM).
!     NOTE: FOR CONVECTIVE BOUNDARIES THE MASS OF THE LAST GRID POINT IS
!           THE MASS OF THE ENTIRE CONVECTION ZONE.
!  composition (HCOMP) : ARRAY OF MASS FRACTION OF ALL OF THE SPECIES AT
!     THE ORIGINAL MODEL POINTS.
!  enclosed_mass (HS1) : LOCATION IN MASS(UNLOGGED) OF THE ORIGINAL MODEL
!     POINTS.
!  shell_mass (HS2) : MASSES OF THE ORIGINAL MODEL POINTS(UNLOGGED).
!  velocity (HV) : RUN OF DIFFUSION VELOCITIES.
!  zone_begin,zone_end (IBEG,IEND) :THE FIRST/LAST UNSTABLE POINTS IN
!     THE REGION.
!     NOTE: FOR CONVECTIVE BOUNDARIES THESE ARE ONLY THE FIRST CONVECTIVE
!           POINTS ADJACENT TO AN UNSTABLE RADIATIVE REGION.
!  zone_min (IMIN) : THE INNERMOST RADIATIVE ZONE OUTSIDE OF ZONE 1.
!  zone_max (IMAX) : the outermost radiative zone.
!  convective_flag (LCZ) : FLAG WHICH TELLS WHICH OF THE ORIGINAL MODEL
!     POINTS ARE CONVECTIVE FOR ANGULAR MOMENTUM TRANSPORT PURPOSES
!     (I.E. INCLUDES OVERSHOOT REGIONS). convective_flag=T IF
!     CONVECTIVE.
!  final_iteration_flag (LOK) : FLAG SET T WHEN FINAL ITERATION IS
!     BEING PERFORMED.
!  num_zones (M) : NUMBER OF MODEL POINTS.
!
!  OUTPUT VARIABLES:
!  composition IS UPDATED IN DIFCOM TO GIVE THE NEW RUN OF COMPOSITION
!  AFTER ANGULAR MOMENTUM TRANSPORT.
!
!  BEFORE THE LAST ITERATION(LOK=F),ONLY DIFFUSION OF H,HE4,HE3 CALCULATED
!  TO CALCULATE CHANGE IN MU GRADIENTS CAUSED BY DIFFUSION.
!
!  WHEN LOK IS T, ONLY DIFFUSION OF SPECIES HEAVIER THAN HE4 IS PERFORMED
!  (AS DIFFUSION OF LIGHTER SPECIES HAS ALREADY BEEN DONE).
!  THE EQUALLY SPACED GRID POINTS FROM THE LAST UNSTABLE REGION SOLVED
!  ARE USED INITIALLY, AND THEN THE PROGRAM CHECKS FOR ADDITIONAL UNSTABLE
!  REGIONS INTERIOR TO THAT.
! CALL diffuse_composition FOR THE PREVIOUSLY IDENTIFIED UNSTABLE REGION IF THIS IS NOT
! THE CONVERGED CALL
subroutine diffuse_composition_driver(timestep, diffusion_coeff, equally_spaced_diffusion_coeff, &
     equally_spaced_mass, log_density, luminosity_lsun, log_pressure, &
     log_radius, log_mass, enclosed_mass, shell_mass, log_total_mass, &
     velocity, zone_begin, zone_end, zone_max, zone_min, convective_flag, &
     final_iteration_flag, num_zones, composition, species_begin, species_end, ierr)
      use star_info_lib, only: json
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(inout) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      double precision, intent(in) :: log_density(json), luminosity_lsun(json), &
           log_pressure(json), log_radius(json), log_mass(json), &
           enclosed_mass(json), shell_mass(json), log_total_mass
      double precision, intent(in) :: velocity(json)
      integer, intent(inout) :: zone_begin, zone_end
      integer, intent(in) :: zone_max, zone_min
      logical, intent(in) :: convective_flag(json), final_iteration_flag
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: species_begin, species_end
      integer, intent(out) :: ierr
      logical :: unstable_region_active, unstable_zone_found, &
           two_zone_region
      integer :: search_start, zone_idx

      ierr = 0
      if (.not.final_iteration_flag) then
         call diffuse_composition(timestep, equally_spaced_diffusion_coeff, &
              equally_spaced_mass, zone_begin, zone_end, &
              convective_flag, final_iteration_flag, num_zones, composition, &
              species_begin, species_end, ierr)
         if (ierr /= 0) return
      else
! FIND UNSTABLE REGIONS IN ORDER, AND CALL composition_grid TO SET UP THE
! EQUALLY SPACED GRID AND diffuse_composition TO MIX THEM IN ORDER
!  EACH UNSTABLE REGION IS SOLVED SEPARATELY STARTING HERE.
!  unstable_zone_found IS SET T IF A NON-ZERO VELOCITY IS ENCOUNTERED.
!  zone_begin IS THE ZONE BELOW THE FIRST NON-ZERO V;zone_end IS THE
!  ZONE ABOVE THE LAST NON-ZERO V.
         unstable_region_active = .false.
         search_start = zone_min
         region_loop: do
         unstable_zone_found = .false.
         do zone_idx = search_start,zone_max
            if (velocity(zone_idx).gt.0.0d0) then
               unstable_zone_found = .true.
               if (.not.unstable_region_active) then
!  START OF UNSTABLE REGION
                  unstable_region_active = .true.
                  zone_begin = zone_idx - 1
               end if
            else if (unstable_region_active) then
!  END OF UNSTABLE REGION
               zone_end = zone_idx - 1
               unstable_region_active = .false.
               search_start = zone_idx + 1
               exit
            end if
         end do
         if (zone_idx > zone_max) then
!  IF THE LAST INTERFACE IS UNSTABLE (NON-ZERO V) ENSURE THAT zone_end
!  IS SET PROPERLY.
         if (unstable_region_active) zone_end = zone_max
         search_start = zone_max + 1
         end if
!  IF NO NON-ZERO V'S ENCOUNTERED, EXIT.
         if (unstable_zone_found) then
!  TRANSFORM TO EQUAL GRID SPACING IN R FOR THE REGION.
         call composition_grid(diffusion_coeff, log_density, luminosity_lsun, &
              log_pressure, log_radius, log_mass, enclosed_mass, shell_mass, &
              log_total_mass, zone_begin, zone_end, convective_flag, &
              num_zones, equally_spaced_diffusion_coeff, equally_spaced_mass, &
              two_zone_region)
!  two_zone_region=T IF TWO ZONES IN UNSTABLE REGION;
!  SKIP IF THIS OCCURS.
         if (two_zone_region) then
            if (search_start.le.zone_max) then
               cycle region_loop
            else
               exit region_loop
            end if
         end if
!  PERFORM COMPOSITION DIFFUSION.
!  UNTIL THE FINAL ITERATION, ONLY COMPOSITION DIFFUSION OF SPECIES WHICH
!  AFFECT GRADIENTS IN MEAN MOLECULAR WEIGHT IS COMPUTED (H,HE3,HE4).
!  ON THE FINAL ITERATION, DIFFUSION OF ALL SPECIES IS PERFORMED.
         call diffuse_composition(timestep, equally_spaced_diffusion_coeff, &
              equally_spaced_mass, zone_begin, zone_end, &
              convective_flag, final_iteration_flag, num_zones, composition, &
              species_begin, species_end, ierr)
         if (ierr /= 0) return
!  RETURN FOR NEXT REGION IF APPLICABLE
         if (search_start.gt.zone_max) exit region_loop
         else
         exit region_loop
         end if
         end do region_loop
      end if
      return
end subroutine diffuse_composition_driver

! ---------------------------------------------------------------
! DIFCOM (this routine's original header name) calculates the
! diffusion of composition due to angular momentum transport. This is
! done by transforming to an equally spaced grid which was previously
! defined (by composition_grid.f90) and then transforming back.
!
! INPUT VARIABLES:
! timestep : diffusion timestep (sec).
! equally_spaced_diffusion_coeff : diffusion coefficients for
!    composition transport at the equally spaced grid points.
! equally_spaced_mass : masses of the equally spaced grid points (gm).
!    NOTE: for convective boundaries the mass of the last grid point is
!          the mass of the entire convection zone.
! composition : array of mass fraction of all of the species at the
!    original model points.
! zone_begin,zone_end : the first/last unstable points in the region.
!    NOTE: for convective boundaries these are only the first
!          convective points adjacent to an unstable radiative region.
! convective_flag : flag which tells which of the original model
!    points are convective for angular momentum transport purposes
!    (i.e. includes overshoot regions). true if convective.
! num_zones : number of model points.
!
! OUTPUT VARIABLES:
! composition is updated in this routine to give the new run of
! composition after angular momentum transport.
!
! Before the last iteration (final_iteration_flag=false), only
! diffusion of H,He4,He3 calculated to calculate change in mu
! gradients caused by diffusion. This routine will therefore be called
! either to mix the other species or just the first 4 depending on
! final_iteration_flag.
subroutine diffuse_composition(timestep, equally_spaced_diffusion_coeff, &
     equally_spaced_mass, zone_begin, zone_end, &
     convective_flag, final_iteration_flag, num_zones, composition, &
     species_begin, species_end, ierr)
      use rotation_scratch_lib
      use star_info_lib, only: json, i_h1, i_he4, i_metals, i_he3
      use numerics_lib
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(in) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json), final_iteration_flag
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: species_begin, species_end
      integer, intent(out) :: ierr

! Tridiagonal-solve work arrays (Thomas algorithm): filled in by
! composition_diffusion_coeffs, then consumed by solve_tridiagonal (solution is
! read back afterward).
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json)

      double precision :: equally_spaced_composition(json)
      integer :: varying_species_id(15)
      integer :: num_varying_species, j_idx, zone_idx, ntab, species_num, &
           orig_zone_idx, i0, i1
      double precision :: test_value

      ierr = 0

! BEFORE THE LAST ITERATION(LOK=F),ONLY DIFFUSION OF H,HE4,HE3 CALCULATED
! TO CALCULATE CHANGE IN MU GRADIENTS CAUSED BY DIFFUSION.
! THIS ROUTINE WILL THEREFORE BE CALLED EITHER TO MIX THE OTHER SPECIES
! OR JUST THE FIRST 4 DEPENDING ON LOK.
! DETERMINE WHICH SPECIES VARY OVER THE UNSTABLE REGION.
      num_varying_species = 0
      do j_idx = species_begin, species_end
         test_value = composition(j_idx,zone_begin)
! ABUNDANCE 2(HE4) IS DEFINED AS 1-X-Z-HE3 SO IT IS NOT SOLVED FOR.
         if (j_idx.ne.2) then
            do zone_idx = zone_begin+1, zone_end
               if (dabs(composition(j_idx,zone_idx)-test_value).gt.1.0d-14) &
                    then
                  num_varying_species = num_varying_species + 1
                  varying_species_id(num_varying_species) = j_idx
                  exit
               end if
            end do
         end if
      end do
      if (num_varying_species.ne.0) then
! NOW SOLVE FOR DIFFUSION OF ALL SPECIES THAT VARY OVER THE
! UNSTABLE REGION USING THE SAME DIFFUSION COEFFICIENTS.
      ntab = zone_end - zone_begin + 1
! FIND RUN OF COMPOSITION AT THE EQUALLY SPACED GRID POINTS AT THE
! START OF THE DIFFUSION TIMESTEP.  THIS IS DONE USING AN OSCULATORY
! SPLINE.
      do species_num = 1, num_varying_species
         do zone_idx = 1, ntab
            orig_zone_idx = zone_idx + zone_begin - 1
            rot_scr%xtab(zone_idx) = rot_scr%chi(zone_idx)
            rot_scr%ytab(zone_idx) = composition(varying_species_id(species_num), &
                 orig_zone_idx)
         end do
         call osplin(rot_scr%echi, equally_spaced_composition, rot_scr%xtab, rot_scr%ytab, ntab, &
              rot_scr%ntot)
! SET UP DIFFUSION EQUATION ARRAYS TO SOLVE FOR COMP AT END OF TSTEP
         call composition_diffusion_coeffs(equally_spaced_diffusion_coeff, rot_scr%dchi, timestep, &
              equally_spaced_composition, equally_spaced_mass, rot_scr%ntot, &
              sub_diag, diag, super_diag, rhs)
! SOLVE MATRIX FOR THE RUN OF COMP AT TIME N+1 AT THE NEW GRID.
         call solve_tridiagonal(sub_diag, diag, super_diag, rhs, rot_scr%ntot, solution, ierr)
         if (ierr /= 0) return
! TRANSFORM BACK TO THE ORIGINAL GRID AND UPDATE HCOMP IN THE
! DIFFUSED REGION. U IS THE NEW RUN OF COMPOSITION IN THE REGION AT THE
! EQUALLY SPACED GRID POINTS.
! FOR THE BOUNDARY POINTS IN THE MODEL,EXTRAPOLATING U PAST THE
! LAST EQUALLY SPACED GRID POINTS CAN LEAD TO SERIOUS ERRORS.
! THEREFORE ADD THE *CHANGE* IN COMPOSITION AT THE EQUALLY SPACED GRID
! POINTS TO HCOMP AND DO NOT REPLACE HCOMP WITH U.
         do zone_idx = 1, rot_scr%ntot
            rot_scr%xtab(zone_idx) = rot_scr%echi(zone_idx)
            rot_scr%ytab(zone_idx) = solution(zone_idx) - &
                 equally_spaced_composition(zone_idx)
         end do
         do zone_idx = 1, ntab
            rot_scr%xval(zone_idx) = rot_scr%chi(zone_idx)
         end do
! USE YVAL AS DUMMY ARRAY FOR THE NEW RUN OF COMP AT THE ORIGINAL
! POINT GRID.
         call osplin(rot_scr%xval, rot_scr%yval, rot_scr%xtab, rot_scr%ytab, rot_scr%ntot, ntab)
! CHECK FOR LOWER CONVECTION ZONE
         if (convective_flag(zone_begin).and.zone_begin.gt.1) then
            do zone_idx = zone_begin-1, 1, -1
               if (.not.convective_flag(zone_idx)) then
                  i0 = zone_idx + 1
                  exit
               end if
            end do
            if (zone_idx .lt. 1) i0 = 1
         else
            i0 = zone_begin
         end if
! CHECK FOR UPPER CONVECTION ZONE.
         if (convective_flag(zone_end).and.zone_end.lt.num_zones) then
            do zone_idx = zone_end+1, num_zones
               if (.not.convective_flag(zone_idx)) then
                  i1 = zone_idx - 1
                  exit
               end if
            end do
            if (zone_idx .gt. num_zones) i1 = num_zones
         else
! The original (mixcom.f) never set I1 on
! this branch, so when the diffused region's top is radiative or is
! the surface, the species-mass sum and the He4 renormalisation loop
! below ran over 1..0 (i1 was zero via -finit-local-zero; whatever
! the previous call left in F77). Mirror the i0 branch.
            i1 = zone_end
         end if
! UPDATE COMPOSITION ARRAY.
         do zone_idx = 1, ntab
            j_idx = zone_begin + zone_idx - 1
            composition(varying_species_id(species_num),j_idx) = &
                 composition(varying_species_id(species_num),j_idx) + &
                 rot_scr%yval(zone_idx)
         end do
! UPDATE INNER CZ COMPOSITION IF APPLICABLE.
         if (i0.lt.zone_begin) then
            do zone_idx = zone_begin-1, i0, -1
               composition(varying_species_id(species_num),zone_idx) = &
                    composition(varying_species_id(species_num),zone_begin)
            end do
         end if
! UPDATE OUTER CZ COMPOSITION IF APPLICABLE.
         if (zone_end.lt.i1) then
            do zone_idx = zone_end+1, i1
               composition(varying_species_id(species_num),zone_idx) = &
                    composition(varying_species_id(species_num),zone_end)
            end do
         end if
      end do
! ADJUST HE4 FOR CHANGES IN X, Z, AND HE3.
      if (.not.final_iteration_flag) then
         do zone_idx = i0, i1
            composition(i_he4,zone_idx) = 1.0d0 - composition(i_h1,zone_idx) - &
                 composition(i_metals,zone_idx) - composition(i_he3,zone_idx)
         end do
      end if
      end if
      return
end subroutine diffuse_composition

! ---------------------------------------------------------------
!  SET UP COEFFICIENTS OF DIFFUSION DIFFERENCE EQUATIONS FOR COMPOSITION
!  TRANSPORT DUE TO ROTATIONALLY INDUCED MIXING.
!
!  INPUT VARIABLES:
!  diffusion_coeff (COD) : RUN OF DIFFUSION COEFFICIENTS (COMPUTED IN
!     SR COEFFT).
!  grid_spacing (DR) : GRID SPACING IN RADIUS(CM).
!  timestep (DT) : DIFFUSION TIMESTEP (SEC).
!  eq_composition (ECOMP) : STARTING ARRAY OF COMPOSITION MASS FRACTION
!     AT THE EQUALLY SPACED GRID.
!  eq_mass (EM) : RUN OF MASSES OF THE EQUALLY SPACED GRID POINTS (GM).
!  num_eq_points (NTOT) : NUMBER OF POINTS IN THE EQUALLY SPACED GRID
!     USED TO SOLVE THE DIFFUSION EQUATION.
!  OUTPUT VARIABLES (sub_diag/diag/super_diag/rhs dummy arguments) :
!  THE ABUNDANCE OF SHELL I AT TIME N+1 (COMP(I,N+1)) IS A FUNCTION OF
!  COMP(I-1),COMP(I),AND COMP(I+1) AS DISCUSSED BELOW.  THIS CAN BE EXPRESSED
!  AS A TRIDIAGONAL MATRIX.  TERMS WHICH DEPEND ON COMP(I-1) ARE STORED IN
!  ARRAY sub_diag; THOSE WHICH DEPEND ON COMP(I) ARE STORED IN ARRAY diag;
!  THOSE WHICH DEPEND ON COMP(I+1) ARE STORED IN ARRAY super_diag; AND THE
!  STARTING RUN OF COMPOSITION IS STORED IN ARRAY rhs.
!
!  THE DIFFUSION EQUATION WE ARE SOLVING IS
!  dX/dT = 1/(4pi*RHO*R**2)d/dR(D*dX/dR)
!  WHERE D IS OUR DIFFUSION COEFFICIENT,X IS THE MASS FRACTION OF THE
!  SPECIES IN QUESTION,T IS TIME,R IS RADIUS, AND RHO IS DENSITY.
!  USING M = 4pi*RHO*R**2*DR THIS IS DIFFERENCED FOR SHELL I AS
!  X(I,N+1)-X(I,N)=(DT/DR)*(1/M)*(D(I+1/2)*(X(I+1,N+1)-X(I,N+1)) -
!  D(I-1/2)*(X(I,N+1)-X(I-1,N+1)))
!  WHERE M IS THE MASS OF SHELL , INDEX N+1 REFERS TO VALUES AT THE
!  END OF THE TIMESTEP, INDEX N REFERS TO VALUES AT THE BEGINNING OF
!  THE TIMESTEP.
!  THIS SYSTEM IS SUPPLEMENTED WITH APPROPRIATE BOUNDARY CONDITIONS.
!
! Composition-transport analogue of am_diffusion_coeffs.f90's angular-momentum-
! transport tridiagonal setup; called by diffuse_composition.f90, which passes the
! output matrix straight on to numerics' solve_tridiagonal.
!
! sub_diag/diag/super_diag/rhs: the tridiagonal system diffuse_composition
! solves (solve_tridiagonal).
subroutine composition_diffusion_coeffs(diffusion_coeff, grid_spacing, timestep, eq_composition, &
     eq_mass, num_eq_points, sub_diag, diag, super_diag, rhs)
      use star_info_lib, only: json
      implicit none

      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: grid_spacing, timestep
      double precision, intent(in) :: eq_composition(json), eq_mass(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      integer :: i
      double precision :: fact0, fact

!  STORE THE RUN OF COMPOSITION AT THE BEGINNING OF THE TIMESTEP IN ARRAY rhs.
      do i = 1,num_eq_points
         rhs(i) = eq_composition(i)
      end do
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.(SAME USED FOR LAST SHELL B.C.)
      fact0 = timestep/grid_spacing
      fact = fact0/eq_mass(1)
      sub_diag(1) = 0.0d0
      diag(1) = 1.0d0 + fact*diffusion_coeff(2)
      super_diag(1) = -fact*diffusion_coeff(2)
      do i = 2,num_eq_points-1
         fact = fact0/eq_mass(i)
         sub_diag(i) = -fact*diffusion_coeff(i)
         diag(i) = 1.0d0 + fact*(diffusion_coeff(i)+diffusion_coeff(i+1))
         super_diag(i) = -fact*diffusion_coeff(i+1)
      end do
!  LAST SHELL BOUNDARY CONDITIONS.
      fact = fact0/eq_mass(num_eq_points)
      sub_diag(num_eq_points) = -fact*diffusion_coeff(num_eq_points)
      diag(num_eq_points) = 1.0d0 + fact*diffusion_coeff(num_eq_points)
      super_diag(num_eq_points) = 0.0d0

      return
end subroutine composition_diffusion_coeffs

! ---------------------------------------------------------------
! Builds the equally-spaced-grid quantities needed by diffuse_composition.f90's
! diffusive composition solve: given the first/last unstable zones of
! a region (zone_begin/zone_end), calls equal_spaced_grid to lay down the
! equally spaced coordinate rot_scr%chi, then computes the equally spaced grid
! masses (equally_spaced_mass) and the geometrically weighted
! diffusion coefficients (equally_spaced_diffusion_coeff) at the zone
! edges, including the Jacobian factor for the transformation from
! radius to the rot_scr%chi coordinate.
subroutine composition_grid(diffusion_coeff, log_density, luminosity_lsun, &
     log_pressure, log_radius, log_mass, enclosed_mass, shell_mass, &
     log_total_mass, zone_begin, zone_end, convective_flag, num_zones, &
     equally_spaced_diffusion_coeff, equally_spaced_mass, &
     single_interface_flag)
      use rotation_scratch_lib
      use equal_grid
      use star_info_lib, only: json
      use physical_constants, only: ln10
      use math_lib
      implicit none

      double precision, intent(in) :: diffusion_coeff(json), &
           log_density(json), luminosity_lsun(json), log_pressure(json), &
           log_radius(json), log_mass(json), enclosed_mass(json), &
           shell_mass(json)
      double precision, intent(in) :: log_total_mass
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      logical, intent(out) :: single_interface_flag
      integer :: idx, search_idx, ntab, ntabb
      double precision :: em_top, em_bot

! FLAG THE SPECIAL CASE OF A SINGLE UNSTABLE INTERFACE AND EXIT
      if (zone_end - zone_begin.le.1) then
         single_interface_flag = .true.
         return
      else
         single_interface_flag = .false.
      end if
! DEFINE A GRID OF EQUALLY SPACED POINTS.
      call equal_spaced_grid(luminosity_lsun, log_pressure, log_mass, zone_begin, &
           zone_end, num_zones)
! EQUAL_SPACED_GRID HAS DEFINED A SET OF CO-ORDINATES (CHI) AND EQUALLY
! SPACED MASS POINTS.  NOW FIND THE OTHER QUANTITIES OF INTEREST AT ZONE
! CENTERS:
! TOTAL ZONE MASSES
! INTERMEDIATE POINTS
      do idx = 2, rot_scr%ntot-1
         equally_spaced_mass(idx) = 0.5d0*(rot_scr%eq_enclosed_mass(idx+1) - rot_scr%eq_enclosed_mass(idx-1))
      end do
! SPECIAL TREATMENT OF THE BOUNDARIES; CAN BE CONVECTIVE.
! IF CONVECTIVE SUM OVER ALL SHELLS.  CARE IS NEEDED TO DO BOOK-KEEPING
! PROPERLY AT THE EDGES - TOP IS HALFWAY TO EQUALLY SPACED POINT, NOT
! HALFWAY TO EDGE OF UNEQUALLY SPACED ORIGINAL SET OF POINTS.
!
! CENTER
      em_top = 0.5d0*(rot_scr%eq_enclosed_mass(2) + rot_scr%eq_enclosed_mass(1))
      if (zone_begin.gt.1) then
         em_bot = 0.5d0*(enclosed_mass(zone_begin) + &
              enclosed_mass(zone_begin-1))
      else
         em_bot = 0.0d0
      end if
      equally_spaced_mass(1) = em_top - em_bot
! (INCLUDE EVERY CONVECTIVE SHELL BELOW ZONE_BEGIN.)
      if (zone_begin.gt.1) then
         do search_idx = zone_begin-1, 1, -1
            if (.not.convective_flag(search_idx)) exit
            equally_spaced_mass(1) = equally_spaced_mass(1) + &
                 shell_mass(search_idx)
         end do
      end if
! SURFACE
      em_bot = 0.5d0*(rot_scr%eq_enclosed_mass(rot_scr%ntot) + rot_scr%eq_enclosed_mass(rot_scr%ntot-1))
      if (zone_end.lt.num_zones) then
         em_top = 0.5d0*(enclosed_mass(zone_end) + &
              enclosed_mass(zone_end+1))
      else
         em_top = exp(ln10*log_total_mass)
      end if
      equally_spaced_mass(rot_scr%ntot) = em_top - em_bot
! (INCLUDE EVERY CONVECTIVE SHELL ABOVE ZONE_END.)
      if (zone_end.lt.num_zones) then
         do search_idx = zone_end+1, num_zones
            if (.not.convective_flag(search_idx)) exit
            equally_spaced_mass(rot_scr%ntot) = equally_spaced_mass(rot_scr%ntot) + &
                 shell_mass(search_idx)
         end do
      end if
! NOW SOLVE FOR QUANTITIES NEEDED AT THE ZONE EDGES.  THESE ARE
! RELATED TO THE DIFFUSION COEFFICIENTS.  UNLIKE THE EQUALLY SPACED
! GRID IN R, WE NEED TO INCLUDE A JACOBIAN TERM FOR THE TRANSFORMATION
! OF VARIABLES.
      ntab = zone_end - zone_begin + 1
      call edge_grid_abscissae(rot_scr%chi, rot_scr%echi, rot_scr%dchi, ntab, &
           rot_scr%ntot, rot_scr%xtab, rot_scr%xval, ntabb)
! DIFFUSION COEFFICIENT FOR MIXING - ASSUME CONSTANT BELOW
! BOTTOM INTERFACE OR ABOVE TOP INTERFACE
      call interp_edge_coeff(diffusion_coeff, zone_begin, zone_end, ntab, ntabb, &
           rot_scr%ntot, rot_scr%xval, rot_scr%xtab, rot_scr%ytab, &
           equally_spaced_diffusion_coeff)
! PRODUCT OF RHO R^2 BY D CHI/DR (LOG10, ON THE EQUAL GRID IN YVAL)
      call dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
           enclosed_mass, mix_scr%epsm, luminosity_lsun(num_zones), zone_begin, &
           ntab, rot_scr%ntot, rot_scr%chi, rot_scr%xval, rot_scr%xtab, &
           rot_scr%ytab, rot_scr%yval)
! NOW ADD MULTIPLICATIVE FACTORS TO DIFFUSION COEFFICIENTS
      call multiply_by_exp10(equally_spaced_diffusion_coeff, rot_scr%yval, rot_scr%ntot)
      return
end subroutine composition_grid

end module rotational_mixing
