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
subroutine oversh(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, num_zones, mixed_zone_bounds, &
     mixed_zone_bounds_no_overshoot, num_mixed_zones)

      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: composition(15,json), &
           log_density(json), log_pressure(json), log_radius(json), &
           log_mass(json), log_temperature(json)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: mixed_zone_bounds(12,2)
      integer, intent(in) :: mixed_zone_bounds_no_overshoot(12,2)
      integer, intent(in) :: num_mixed_zones



! common/dpmix/: dpenv, alphac, alphae, alpham, betac, lovstc,
! envelope_overshoot_active, lovstm, and lovmax are used here. Naming
! matches liburn.f90.
      double precision :: dpenv, alphac, alphae, alpham, betac
      integer :: iov1, iov2, iovim
      logical :: lovstc, envelope_overshoot_active, lovstm, lsemic, ladov, &
           lovmax
      common/dpmix/ dpenv, alphac, alphae, alpham, betac, iov1, iov2, &
           iovim, lovstc, envelope_overshoot_active, lovstm, lsemic, ladov, &
           lovmax

! common/rot/: only rotation_active and instability_transport_active
! are used here. Naming matches liburn.f90.
      double precision :: wnew, walpcz, acfpft
      integer :: itfp1, itfp2
      logical :: rotation_active, instability_transport_active, lwnew
      common/rot/ wnew, walpcz, acfpft, itfp1, itfp2, rotation_active, &
           instability_transport_active, lwnew

      save

      logical :: up_overshoot_flag, down_overshoot_flag
      integer :: zone_idx, edge_idx, j_idx
      double precision :: pscale_up, pscale_down
      double precision :: cz_radius, overshoot_radius, radius

! IOV1/IOV2 (from common/dpmix/) store the position of overshoot for
! adiabatic extension.
      iov1 = -1
      iov2 = -1
      do 100 zone_idx = 1, num_mixed_zones
! DETERMINE IF THIS REGION IS A CORE CONVECTION ZONE, SURFACE CZ,
! OR INTERMEDIATE CZ. THERE ARE SEPARATE FLAGS GOVERNING WHETHER
! OVERSHOOT WILL BE PERFORMED IN EACH CASE, AND SEPARATE USER
! PARAMETERS GOVERNING THE DEGREE OF OVERSHOOT.
         if (mixed_zone_bounds(zone_idx,1).eq.1) then
! CONVECTIVE CORE
! CHECK FOR A FULLY CONVECTIVE STAR; SKIP THIS SR IF THERE IS ONE.
            if (mixed_zone_bounds(zone_idx,2).eq.num_zones) then
               write(short_file_unit,5)
    5          format(1x,'FULLY CONVECTIVE MODEL - NO OVERSHOOT')
               write(short_file_unit,200) (mixed_zone_bounds_no_overshoot( &
                    1,j_idx), j_idx=1,2)
               return
            end if
! SKIP IF NO CORE OVERSHOOT IS DESIRED.
            if (.not.lovstc) goto 100
            up_overshoot_flag = .true.
            down_overshoot_flag = .false.
            edge_idx = mixed_zone_bounds(zone_idx,2)
            call hsubp(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_up)
! PSCALU IS THE PRESSURE SCALE HEIGHT ABOVE THE CONVECTIVE REGION;
! ALPHAC IS THE DESIRED OVERSHOOT (IN SCALE HEIGHTS).
!            PSCALU = PSCALU*ALPHAC
! JVS 07/13 ALLOW FOR THE LIMITING OF OVERSHOOTING ABOVE THE CONVECTIVE
! CORES OF LOVE MASS STARS AS PER WOO & DEMARQUE 2001
            if (lovmax) then
               pscale_up = min(pscale_up*alphac, &
                    betac*exp(ln10*log_radius(edge_idx)))
            else
               pscale_up = pscale_up*alphac
            end if

         else if (mixed_zone_bounds(zone_idx,2).eq.num_zones) then
! CONVECTIVE ENVELOPE
! SKIP IF NO ENVELOPE OVERSHOOT IS DESIRED.
            if (.not.envelope_overshoot_active) goto 100
            up_overshoot_flag = .false.
            down_overshoot_flag = .true.
            edge_idx = mixed_zone_bounds(zone_idx,1)
            call hsubp(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_down)
! PSCALD IS THE PRESSURE SCALE HEIGHT BELOW THE CONVECTIVE REGION;
! ALPHAE IS THE DESIRED OVERSHOOT (IN SCALE HEIGHTS).
            pscale_down = pscale_down*alphae
         else
! INTERMEDIATE CONVECTION ZONE (NOT INCLUDING CENTRAL OR SURFACE POINT).
! SKIP IF NO INTERMEDIATE CONVECTION.
            if (.not.lovstm) goto 100
            up_overshoot_flag = .true.
            down_overshoot_flag = .true.
! PSCALU AND PSCALD HAVE THE SAME MEANING AS ABOVE; OVERSHOOT BOTH BELOW
! AND ABOVE IS PERFORMED BY AN AMOUNT ALPHAM.
            edge_idx = mixed_zone_bounds(zone_idx,1)
            call hsubp(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_down)
            pscale_down = pscale_down*alpham
            edge_idx = mixed_zone_bounds(zone_idx,2)
            call hsubp(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_up)
            pscale_up = pscale_up*alpham
         end if
! COMPUTE EXTENSION OF CONVECTION ZONE BELOW SCHWARTZSCHILD BOUNDARY.
         if (down_overshoot_flag) then
            edge_idx = mixed_zone_bounds(zone_idx,1)
            cz_radius = exp(ln10*log_radius(edge_idx))
            overshoot_radius = cz_radius - pscale_down
! THE OVERSHOOT REGION IS EXTENDED THE RADIAL DISTANCE PSCALD DOWN; THE
! LAST POINT LESS THAN PSCALD FROM THE FORMAL EDGE OF THE CZ IS DEFINED
! AS THE NEW EDGE OF THE MIXED REGION.
            do 10 j_idx = edge_idx-1, 1, -1
               radius = exp(ln10*log_radius(j_idx))
               if (radius.lt.overshoot_radius) goto 20
   10       continue
! IF THE CODE GETS HERE, THE OVERSHOOT REGION EXTENDS BELOW THE FIRST POINT.
! THE CODE WILL ASSIGN THE FIRST POINT AS THE LOWER EDGE(I.E. THE CZ WILL
! EXTEND TO THE CENTER).
            j_idx = 0
   20       continue
! FOR ROTATING MODELS, ENSURE THAT THERE IS AT LEAST ONE RADIATIVE POINT
! IN THE OVERSHOOT REGION.
            mixed_zone_bounds(zone_idx,1) = j_idx + 1
! 11/91 MHP CHANGED TO REQUIRE AN OVERSHOOT ZONE ONLY IF LINSTB=T.
            if (rotation_active .and. instability_transport_active .and. &
                 mixed_zone_bounds(zone_idx,1).eq. &
                 mixed_zone_bounds_no_overshoot(zone_idx,1)) &
                 mixed_zone_bounds(zone_idx,1) = &
                 mixed_zone_bounds(zone_idx,1) - 1
! DBG 8/94 STORE POSITION OF OVERSHOOT FOR ADIABATIC EXTENSION
            iov2 = edge_idx
            iov1 = mixed_zone_bounds(zone_idx,1)
         end if
! COMPUTE EXTENSION OF CONVECTION ZONE ABOVE SCHWARTZSCHILD BOUNDARY.
         if (up_overshoot_flag) then
            edge_idx = mixed_zone_bounds(zone_idx,2)
            cz_radius = exp(ln10*log_radius(edge_idx))
            overshoot_radius = cz_radius + pscale_up
! THE OVERSHOOT REGION IS EXTENDED THE RADIAL DISTANCE PSCALU UP; THE
! LAST POINT LESS THAN PSCALU FROM THE FORMAL EDGE OF THE CZ IS DEFINED
! AS THE NEW EDGE OF THE MIXED REGION.
            do 30 j_idx = edge_idx+1, num_zones
               radius = exp(ln10*log_radius(j_idx))
               if (radius.gt.overshoot_radius) goto 40
   30       continue
! IF THE CODE GETS HERE, THE OVERSHOOT REGION EXTENDS ABOVE THE LAST POINT.
! THE CODE WILL ASSIGN THE LAST POINT AS THE UPPER EDGE(I.E. THE CZ WILL
! EXTEND TO THE SURFACE).
            j_idx = num_zones + 1
   40       continue
            mixed_zone_bounds(zone_idx,2) = j_idx - 1
! 11/91 MHP CHANGED TO REQUIRE AN OVERSHOOT ZONE ONLY IF LINSTB=T.
            if (rotation_active .and. instability_transport_active .and. &
                 mixed_zone_bounds(zone_idx,2).eq. &
                 mixed_zone_bounds_no_overshoot(zone_idx,2)) &
                 mixed_zone_bounds(zone_idx,2) = &
                 mixed_zone_bounds(zone_idx,2) + 1
         end if
  100 continue
! OUTPUT : THE OLD AND NEW MIXED REGIONS ARE PRINTED OUT IN ISHORT.
      write(short_file_unit,200) ((mixed_zone_bounds_no_overshoot( &
           zone_idx,j_idx), j_idx=1,2), zone_idx=1,num_mixed_zones)
  200 format(1x,'MIXED REGIONS WITHOUT OVERSHOOT', &
           4('[',i4,'-',i4,' ]'))
      write(short_file_unit,210) ((mixed_zone_bounds(zone_idx,j_idx), &
           j_idx=1,2), zone_idx=1,num_mixed_zones)
  210 format(1x,'MIXED REGIONS WITH OVERSHOOT   ', &
           4('[',i4,'-',i4,' ]'))

      return
end subroutine oversh
