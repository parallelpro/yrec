!----------------------------------------------------------------------
! mixcom
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mixcom.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
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
! composition_diffusion_coeffs, then consumed by ctridi (solution is
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
         call ctridi(rot_scr%ntot, sub_diag, diag, super_diag, rhs, solution, ierr)
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
! 2026 (bugsweep sec-11): the original (mixcom.f) never set I1 on
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
