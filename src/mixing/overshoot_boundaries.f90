!----------------------------------------------------------------------
! oversh
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original oversh.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine computes the local pressure scale height at both edges
! of a given convective region and locates the boundaries of overshoot
! regions based on the user-specified extent.
subroutine overshoot_boundaries(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, num_zones, mixed_zone_bounds, &
     mixed_zone_bounds_no_overshoot, num_mixed_zones, ierr)
      use star_info_lib, only: star
      use star_info_lib, only: star, json
      use luout_lib
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: composition(15,json), &
           log_density(json), log_pressure(json), log_radius(json), &
           log_mass(json), log_temperature(json)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: mixed_zone_bounds(12,2)
      integer, intent(in) :: mixed_zone_bounds_no_overshoot(12,2)
      integer, intent(in) :: num_mixed_zones
      integer, intent(out) :: ierr
      logical :: up_overshoot_flag, down_overshoot_flag
      integer :: zone_idx, edge_idx, j_idx
      double precision :: pscale_up, pscale_down
      double precision :: cz_radius, overshoot_radius, radius

! IOV1/IOV2 (from common/dpmix/) store the position of overshoot for
! adiabatic extension.
      star%iov1 = -1
      star%iov2 = -1
      ierr = 0
      do zone_idx = 1, num_mixed_zones
! DETERMINE IF THIS REGION IS A CORE CONVECTION ZONE, SURFACE CZ,
! OR INTERMEDIATE CZ. THERE ARE SEPARATE FLAGS GOVERNING WHETHER
! OVERSHOOT WILL BE PERFORMED IN EACH CASE, AND SEPARATE USER
! PARAMETERS GOVERNING THE DEGREE OF OVERSHOOT.
         if (mixed_zone_bounds(zone_idx,1).eq.1) then
! CONVECTIVE CORE
! CHECK FOR A FULLY CONVECTIVE STAR; SKIP THIS SR IF THERE IS ONE.
            if (mixed_zone_bounds(zone_idx,2).eq.num_zones) then
               write(run_log_unit,5)
    5          format(1x,'FULLY CONVECTIVE MODEL - NO OVERSHOOT')
               write(run_log_unit,200) (mixed_zone_bounds_no_overshoot( &
                    1,j_idx), j_idx=1,2)
               return
            end if
! SKIP IF NO CORE OVERSHOOT IS DESIRED.
            if (.not.star%job%core_overshoot_active) cycle
            up_overshoot_flag = .true.
            down_overshoot_flag = .false.
            edge_idx = mixed_zone_bounds(zone_idx,2)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_up, ierr)
            if (ierr /= 0) return
! PSCALU IS THE PRESSURE SCALE HEIGHT ABOVE THE CONVECTIVE REGION;
! ALPHAC IS THE DESIRED OVERSHOOT (IN SCALE HEIGHTS).
!            PSCALU = PSCALU*ALPHAC
! JVS 07/13 ALLOW FOR THE LIMITING OF OVERSHOOTING ABOVE THE CONVECTIVE
! CORES OF LOVE MASS STARS AS PER WOO & DEMARQUE 2001
            if (star%ctrl%lovmax) then
               pscale_up = min(pscale_up*star%ctrl%overshoot_alpha_core, &
                    star%ctrl%betac*exp(ln10*log_radius(edge_idx)))
            else
               pscale_up = pscale_up*star%ctrl%overshoot_alpha_core
            end if

         else if (mixed_zone_bounds(zone_idx,2).eq.num_zones) then
! CONVECTIVE ENVELOPE
! SKIP IF NO ENVELOPE OVERSHOOT IS DESIRED.
            if (.not.star%job%envelope_overshoot_active) cycle
            up_overshoot_flag = .false.
            down_overshoot_flag = .true.
            edge_idx = mixed_zone_bounds(zone_idx,1)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_down, ierr)
            if (ierr /= 0) return
! PSCALD IS THE PRESSURE SCALE HEIGHT BELOW THE CONVECTIVE REGION;
! ALPHAE IS THE DESIRED OVERSHOOT (IN SCALE HEIGHTS).
            pscale_down = pscale_down*star%ctrl%overshoot_alpha_envelope
         else
! INTERMEDIATE CONVECTION ZONE (NOT INCLUDING CENTRAL OR SURFACE POINT).
! SKIP IF NO INTERMEDIATE CONVECTION.
            if (.not.star%job%lovstm) cycle
            up_overshoot_flag = .true.
            down_overshoot_flag = .true.
! PSCALU AND PSCALD HAVE THE SAME MEANING AS ABOVE; OVERSHOOT BOTH BELOW
! AND ABOVE IS PERFORMED BY AN AMOUNT ALPHAM.
            edge_idx = mixed_zone_bounds(zone_idx,1)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_down, ierr)
            if (ierr /= 0) return
            pscale_down = pscale_down*star%ctrl%alpham
            edge_idx = mixed_zone_bounds(zone_idx,2)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_up, ierr)
            if (ierr /= 0) return
            pscale_up = pscale_up*star%ctrl%alpham
         end if
! COMPUTE EXTENSION OF CONVECTION ZONE BELOW SCHWARTZSCHILD BOUNDARY.
         if (down_overshoot_flag) then
            edge_idx = mixed_zone_bounds(zone_idx,1)
            cz_radius = exp(ln10*log_radius(edge_idx))
            overshoot_radius = cz_radius - pscale_down
! THE OVERSHOOT REGION IS EXTENDED THE RADIAL DISTANCE PSCALD DOWN; THE
! LAST POINT LESS THAN PSCALD FROM THE FORMAL EDGE OF THE CZ IS DEFINED
! AS THE NEW EDGE OF THE MIXED REGION.
            do j_idx = edge_idx-1, 1, -1
               radius = exp(ln10*log_radius(j_idx))
               if (radius.lt.overshoot_radius) exit
            end do
! IF THE LOOP COMPLETES, THE OVERSHOOT REGION EXTENDS BELOW THE FIRST
! POINT AND THE CZ WILL EXTEND TO THE CENTER: natural completion of the
! downward loop already leaves j_idx = 0 (the old explicit assignment
! before label 20 was redundant).
! FOR ROTATING MODELS, ENSURE THAT THERE IS AT LEAST ONE RADIATIVE POINT
! IN THE OVERSHOOT REGION.
            mixed_zone_bounds(zone_idx,1) = j_idx + 1
! 11/91 MHP CHANGED TO REQUIRE AN OVERSHOOT ZONE ONLY IF LINSTB=T.
            if (star%job%rotation_active .and. star%job%instability_transport_active .and. &
                 mixed_zone_bounds(zone_idx,1).eq. &
                 mixed_zone_bounds_no_overshoot(zone_idx,1)) &
                 mixed_zone_bounds(zone_idx,1) = &
                 mixed_zone_bounds(zone_idx,1) - 1
! DBG 8/94 STORE POSITION OF OVERSHOOT FOR ADIABATIC EXTENSION
            star%iov2 = edge_idx
            star%iov1 = mixed_zone_bounds(zone_idx,1)
         end if
! COMPUTE EXTENSION OF CONVECTION ZONE ABOVE SCHWARTZSCHILD BOUNDARY.
         if (up_overshoot_flag) then
            edge_idx = mixed_zone_bounds(zone_idx,2)
            cz_radius = exp(ln10*log_radius(edge_idx))
            overshoot_radius = cz_radius + pscale_up
! THE OVERSHOOT REGION IS EXTENDED THE RADIAL DISTANCE PSCALU UP; THE
! LAST POINT LESS THAN PSCALU FROM THE FORMAL EDGE OF THE CZ IS DEFINED
! AS THE NEW EDGE OF THE MIXED REGION.
            do j_idx = edge_idx+1, num_zones
               radius = exp(ln10*log_radius(j_idx))
               if (radius.gt.overshoot_radius) exit
            end do
! IF THE LOOP COMPLETES, THE OVERSHOOT REGION EXTENDS ABOVE THE LAST
! POINT AND THE CZ WILL EXTEND TO THE SURFACE: natural completion of
! the upward loop already leaves j_idx = num_zones + 1 (the old
! explicit assignment before label 40 was redundant).
            mixed_zone_bounds(zone_idx,2) = j_idx - 1
! 11/91 MHP CHANGED TO REQUIRE AN OVERSHOOT ZONE ONLY IF LINSTB=T.
            if (star%job%rotation_active .and. star%job%instability_transport_active .and. &
                 mixed_zone_bounds(zone_idx,2).eq. &
                 mixed_zone_bounds_no_overshoot(zone_idx,2)) &
                 mixed_zone_bounds(zone_idx,2) = &
                 mixed_zone_bounds(zone_idx,2) + 1
         end if
      end do
! OUTPUT : THE OLD AND NEW MIXED REGIONS ARE PRINTED OUT IN ISHORT.
      write(run_log_unit,200) ((mixed_zone_bounds_no_overshoot( &
           zone_idx,j_idx), j_idx=1,2), zone_idx=1,num_mixed_zones)
  200 format(1x,'MIXED REGIONS WITHOUT OVERSHOOT', &
           4('[',i4,'-',i4,' ]'))
      write(run_log_unit,210) ((mixed_zone_bounds(zone_idx,j_idx), &
           j_idx=1,2), zone_idx=1,num_mixed_zones)
  210 format(1x,'MIXED REGIONS WITH OVERSHOOT   ', &
           4('[',i4,'-',i4,' ]'))

      return
end subroutine overshoot_boundaries
