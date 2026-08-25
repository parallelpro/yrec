!----------------------------------------------------------------------
! tauintnew
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original tauintnew.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! G Somers 3/17, IMPLEMENTING NEW TAUCZ METHOD.
! THE METHOD WORKS AS FOLLOWS. SEARCH FROM THE CENTER OF THE STAR UP UNTIL YOU
! FIND THE BASE OF THE SURFACE CONVECTION ZONE. THEN SEARCH FROM THE SURFACE
! DOWN UNTIL YOU FIND THE SAME. IF BOTH DIRECTIONS GIVE SAME ANSWER, THEN
! SURFACE CONVECTION ZONE IS UNIQUE, AND THE OLD TOP-DOWN SINKLINE METHOD CAN
! BE EMPLOYED. IF THEY ARE DIFFERENT, THEN MULTIPLE SURFACE CONVECTION CELLS
! EXIST. IN THIS CASE, DETERMINE THE CROSSING TIME OF THE LOWEST CONVECTION
! CELL IN THE ENVELOPE.
!
! IMAX: BASE OF THE MAIN INTERNAL CONVECTION ZONE.
! JMAX: BASE OF THE UPPER-MOST CONVECTION ZONE.
! RBCZ: RADIUS AT WHICH THE BCZ EXISTS
! TAUCZ: OVERTURN TIMESCALE
subroutine tauintnew(shell_mass, convective_flag, log10_radius, &
     log10_pressure, log10_density, local_gravity, num_points, &
     num_interior_points, convective_velocity, radiative_gradient, &
     adiabatic_gradient, radius_at_bcz)
      use atm_table_lib
      use envstruct_lib
      use const_lib
      use star_info_lib, only: star, json
      use numerics_lib
      implicit none

      double precision, intent(in) :: shell_mass(json)
      logical, intent(in) :: convective_flag(json)
      double precision, intent(inout) :: log10_radius(json), &
           log10_pressure(json), log10_density(json), local_gravity(json)
      integer, intent(in) :: num_points, num_interior_points
      double precision, intent(in) :: convective_velocity(json)
      double precision, intent(in) :: radiative_gradient(json), &
           adiabatic_gradient(json)
      double precision, intent(out) :: radius_at_bcz
! --- locals ---
! G Somers; Adding vectors for cubic spline int.
      double precision :: spline_x_delta(4), spline_x_radius(4), &
           spline_y_radius(4), spline_y_velocity(4), spline_y_pscale(4), &
           spline_deriv(4)
      logical :: fully_convective_flag, surface_cz_deep_enough
      integer :: i, j, k, kk, search_start_index
      logical :: in_radiative_zone
      integer :: core_cz_top_index, cz_base_index, upper_cz_base_index, &
           index_gap, cz_top_index
      double precision :: mass_at_cz_top, mass_at_cz_bottom
      double precision :: dd2, dd1, interp_fraction
      double precision :: pressure_scale_height2, radius_test2, val_test
      double precision :: pressure_scale_height1, radius_test1
      double precision :: pressure_scale_height, convective_velocity_bcz
      double precision :: log10_radius_interp, radius_top_cz, cz_width, &
           radius_test
      logical :: spline_taucz_done

! DETERMINE EXTENT OF SURFACE CONVECTION ZONE.
      fully_convective_flag = .false.
! FIND THE TOP OF THE MAIN INTERNAL CONVECTION ZONE.
      search_start_index = 1
! (Restructured 2026: the label-50 retry became the named search loop;
! the label 70/90/100 jumps became the if/else chain below; the label
! 140 jump became the spline_taucz_done flag.)
      search: do
      in_radiative_zone = .true.
      if (convective_flag(search_start_index)) in_radiative_zone = .false.
      do i = search_start_index,num_points,1
         if (.not.in_radiative_zone.and..not.convective_flag(i)) then
            core_cz_top_index = i-1
            in_radiative_zone = .true.
         endif
         if (in_radiative_zone.and.convective_flag(i)) exit
      end do
! IF THE SCAN COMPLETES, EITHER THE STAR IS FULLY CONVECTIVE, OR
! HAS A COMPLETELY RADIATIVE ENVLOPE, MAYBE WITH A CONVECTIVE CORE.
      if (i .gt. num_points .and. in_radiative_zone) then ! RADIATIVE ENVELOPE
         fully_convective_flag = .false.
         surface_cz_deep_enough = .false.
         cz_base_index = num_points
         upper_cz_base_index = num_points
         core_cz_top_index = num_points
      else
      if (i .gt. num_points) then ! FULLY CONVECTIVE
         fully_convective_flag = .true.
         surface_cz_deep_enough = .true.
         i = 1
         j = 0
         core_cz_top_index = 1
      else
! FIND THE BASE OF THE UPPER-MOST CONVECTION ZONE
! JvS 02/2026: ADDED LCZ(J+1) TO TRAP OUT CASES WHERE THE LAST (FEW) POINT(S) IN
! THE MODEL IS(ARE) ALREADY RADIATIVE
      do j = num_points-1,1,-1
         if (convective_flag(j+1) .and. .not.convective_flag(j)) exit
      end do
      endif
      cz_base_index = i
      upper_cz_base_index = j + 1
      index_gap = upper_cz_base_index-cz_base_index
! LOCATING THE CORRECT CZ INCURES VARIOUS PATHOLOGIES. TRAP THEM OUT.
!
! ENSURE THIS CZ IS SUFFICIENTLY FAR FROM THE CONVECTIVE CORE. IN SOME CASES
! THERE ARE A FEW SHELLS RIGHT ABOVE THE CORE THAT ARE CONVECTIVE. SKIP THESE.
! KC 2025-05-30 ADDED IF CHECK TO AVOID RUNTIME ERROR.
      if (cz_base_index .gt. 0 .and. core_cz_top_index .gt. 0) then
         if (num_points.ne.num_interior_points .and. &
              log10_pressure(core_cz_top_index)-log10_pressure(cz_base_index).lt.1) then
            search_start_index = cz_base_index+1
            cycle search
         endif
      endif
! ENSURE THIS CZ IS AT LEAST 3 SHELLS THICK. IF NOT, WE ARE LIKELY SNAGGED
! ON A SINGLE SHELL CONV ZONE DEEP IN THE INTERIOR. RECOMMENSE SEARCH FROM
! ABOVE THIS SHELL.
      if (.not.(convective_flag(cz_base_index).and.convective_flag(cz_base_index+1) &
           .and.convective_flag(cz_base_index+2))) then
         if (.not.convective_flag(cz_base_index+1)) then
            search_start_index = cz_base_index+1
         else
            search_start_index = cz_base_index+2
         endif
         cycle search
      endif
! IN THE CASE WHERE ONLY THE INTERIOR MODEL IS CONSIDERED, USE JMAX.
      if (index_gap.ne.0.and.num_points.eq.num_interior_points) then
         cz_base_index = upper_cz_base_index
         index_gap = 0
      endif
!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
      mass_at_cz_top = shell_mass(num_points)
      if (upper_cz_base_index.gt.1) then
         mass_at_cz_bottom = 0.5d0*(shell_mass(upper_cz_base_index)+shell_mass(upper_cz_base_index-1))
      else
         mass_at_cz_bottom = 0.0d0
      endif
!  LCZSUR=T IF A SURFACE C.Z.DEEP ENOUGH FOR ANGULAR MOMENTUM LOSS EXISTS
      if ((mass_at_cz_top-mass_at_cz_bottom)/star%solar_mass_cgs.gt.0.0d0) then
         surface_cz_deep_enough = .true.
      else
         surface_cz_deep_enough= .false.
      endif
      end if

! CALCULATE TAUCZ TOP DOWN FOR FULLY AND PARTIALLY CONVECTIVE. FOR A FULLY
! RADIATIVE ENVELOPE, SET TAUCZ = 0.0.
      if (surface_cz_deep_enough) then
! INFER RBCZ, THE RADIUS AT THE BASE OF THE SURFACE CZ.
         if (fully_convective_flag) then
!           SET TO CENTER OF STAR, I.E. A VERY SMALL NUMBER
!           IN CGS UNITS, RELATIVE TO THE RADII OF STARS.
            radius_at_bcz = 1.
         else
!           PINPOINT RCZ LINEARLY
            dd2 = radiative_gradient(cz_base_index-1)-adiabatic_gradient(cz_base_index-1)
            dd1 = radiative_gradient(cz_base_index)-adiabatic_gradient(cz_base_index)
            interp_fraction = dd2/(dd2-dd1)
            radius_at_bcz = log10_radius(cz_base_index-1)+interp_fraction* &
                 (log10_radius(cz_base_index)-log10_radius(cz_base_index-1))
         endif
         radius_at_bcz = exp(ln10*radius_at_bcz)
!        IF THERE ARE MULTIPLE CONVECTION CELLS IN THE PREDOMINATELY
!        RADIATIVE ENVELOPE, USE AVERAGING METHOD.
         spline_taucz_done = .false.
         if (index_gap.eq.0) then
! CALCULATE HP
         pressure_scale_height2 = exp(ln10*(log10_pressure(cz_base_index)- &
              log10_density(cz_base_index)))/local_gravity(cz_base_index)
         radius_test2 = exp(ln10*log10_radius(cz_base_index))-radius_at_bcz
         val_test = radius_test2 - pressure_scale_height2
         if (val_test.gt.0.0) then
! INNER MOST POINT PSCA IS BELOW BCZ
            k = cz_base_index
            convective_velocity_bcz = convective_velocity(cz_base_index)
            pressure_scale_height = pressure_scale_height2
            star%convective_turnover_timescale = pressure_scale_height/convective_velocity_bcz
         else
            do k = cz_base_index+1,num_points-1,1
               pressure_scale_height1 = pressure_scale_height2
               radius_test1 = radius_test2
               pressure_scale_height2 = exp(ln10*(log10_pressure(k)-log10_density(k)))/local_gravity(k)
               radius_test2 = exp(ln10*log10_radius(k))-radius_at_bcz
               val_test = radius_test2 - pressure_scale_height2
!              FIND LOCATION WHERE HP = R
               if (val_test.gt.0.0) then
!                 IF K = IMAX+1, MOVE AWAY FROM EDGE CASE.
                  kk = k
                  if (k.eq.cz_base_index+1) kk = k+1
!                 FIND V WITH CUBIC SPLINE
                  spline_x_delta(1) = exp(ln10*log10_radius(kk-2))-radius_at_bcz- &
                            exp(ln10*(log10_pressure(kk-2)-log10_density(kk-2)))/local_gravity(kk-2)
                  spline_x_delta(2) = exp(ln10*log10_radius(kk-1))-radius_at_bcz- &
                            exp(ln10*(log10_pressure(kk-1)-log10_density(kk-1)))/local_gravity(kk-1)
                  spline_x_delta(3) = exp(ln10*log10_radius(kk))-radius_at_bcz- &
                            exp(ln10*(log10_pressure(kk)-log10_density(kk)))/local_gravity(kk)
                  spline_x_delta(4) = exp(ln10*log10_radius(kk+1))-radius_at_bcz- &
                            exp(ln10*(log10_pressure(kk+1)-log10_density(kk+1)))/local_gravity(kk+1)
                  spline_y_velocity(1) = convective_velocity(kk-2)
                  spline_y_velocity(2) = convective_velocity(kk-1)
                  spline_y_velocity(3) = convective_velocity(kk)
                  spline_y_velocity(4) = convective_velocity(kk+1)
                  call kspline(spline_x_delta,spline_y_velocity,spline_deriv)
                  call ksplint(spline_x_delta,spline_y_velocity,spline_deriv,0.0d0,convective_velocity_bcz)
                  spline_y_pscale(1) = exp(ln10*(log10_pressure(kk-2)-log10_density(kk-2)))/local_gravity(kk-2)
                  spline_y_pscale(2) = exp(ln10*(log10_pressure(kk-1)-log10_density(kk-1)))/local_gravity(kk-1)
                  spline_y_pscale(3) = exp(ln10*(log10_pressure(kk)-log10_density(kk)))/local_gravity(kk)
                  spline_y_pscale(4) = exp(ln10*(log10_pressure(kk+1)-log10_density(kk+1)))/local_gravity(kk+1)
                  call kspline(spline_x_delta,spline_y_pscale,spline_deriv)
                  call ksplint(spline_x_delta,spline_y_pscale,spline_deriv,0.0d0,pressure_scale_height)
                  spline_y_radius(1) = log10_radius(kk-2)
                  spline_y_radius(2) = log10_radius(kk-1)
                  spline_y_radius(3) = log10_radius(kk)
                  spline_y_radius(4) = log10_radius(kk+1)
                  call kspline(spline_x_delta,spline_y_radius,spline_deriv)
                  call ksplint(spline_x_delta,spline_y_radius,spline_deriv,0.0d0,log10_radius_interp)
! DEFINE TAUCZ
                  star%convective_turnover_timescale = pressure_scale_height/convective_velocity_bcz
                  spline_taucz_done = .true.
                  exit
               endif
            end do
! KC 2025-05-31 MOVED ENDIF HERE TO AVOID BLOCK MISMATCH.
         endif
         end if
         if (.not. spline_taucz_done) then
! IF CODE GETS HERE, THEN 1 PSCA ABOVE THE BCZ IS OUTSIDE OF THE STAR,
! OR MULTIPLE CONVECTION CELLS RESIDE WITHIN THE RADIATIVE ENVELOPE.
! IN THIS CASE, EVALUATE CVEL AT CENTER OF BOTTOM CZ. FIRST, FIND UPPER EDGE.
            do k = cz_base_index,num_points-1,1
               if (.not.convective_flag(k))exit
            enddo
            cz_top_index = k-1
! PINPOINT TOP OF CZ LINEARLY
            dd2 = radiative_gradient(cz_top_index)-adiabatic_gradient(cz_top_index)
            dd1 = radiative_gradient(cz_top_index+1)-adiabatic_gradient(cz_top_index+1)
            interp_fraction = dd2/(dd2-dd1)
            radius_top_cz = log10_radius(cz_top_index)+interp_fraction*(log10_radius(cz_top_index+1)-log10_radius(cz_top_index))
            radius_top_cz = exp(ln10*radius_top_cz)
! CALCULATE WIDTH OF THE CONVECTION CELL
            cz_width = radius_top_cz - radius_at_bcz
! FIND INDICES SURROUNDING LINEAR CENTER OF CZ
            radius_test = log10(radius_at_bcz + 0.5*cz_width)
            do k = cz_base_index,cz_top_index+1,1
               if (log10_radius(k).gt.radius_test) exit
            enddo
            radius_test = exp(ln10*radius_test)
!           IF K = IMAX+1, MOVE AWAY FROM EDGE CASE.
            kk = k
            if (k.eq.cz_base_index+1) kk = k+1
! CALCULATE CVEL AT THAT CENTRAL LOCATION
            spline_x_radius(1) = exp(ln10*log10_radius(kk-2))
            spline_x_radius(2) = exp(ln10*log10_radius(kk-1))
            spline_x_radius(3) = exp(ln10*log10_radius(kk))
            spline_x_radius(4) = exp(ln10*log10_radius(kk+1))
            spline_y_velocity(1) = convective_velocity(kk-2)
            spline_y_velocity(2) = convective_velocity(kk-1)
            spline_y_velocity(3) = convective_velocity(kk)
            spline_y_velocity(4) = convective_velocity(kk+1)
            call kspline(spline_x_radius,spline_y_velocity,spline_deriv)
            call ksplint(spline_x_radius,spline_y_velocity,spline_deriv,radius_test,convective_velocity_bcz)
! DEFINE TAUCZ
            star%convective_turnover_timescale = cz_width/convective_velocity_bcz
         end if
! KC 2025-05-31 MOVED ENDIF ABOVE TO AVOID BLOCK MISMATCH.
!          ENDIF
!        CONVERT CORE RADIUS INTO SOLAR UNITS
         radius_at_bcz = radius_at_bcz/star%solar_radius_cgs
      else
         star%convective_turnover_timescale = 0.0d0
         radius_at_bcz = 0.0
      endif
!     ENSURE THAT TAUCZ WAS NOT ACCIDENTALLY CALCULATED
!     DEEP IN THE STELLAR INTERIOR. IF YES, REDO CALCULATION.
      if (star%convective_turnover_timescale.gt.1.0e20) then
         search_start_index = cz_base_index + 1
         cycle search
      endif
      exit search
      end do search
      return
end subroutine tauintnew
