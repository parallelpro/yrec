!----------------------------------------------------------------------
! taucal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original taucal.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! JVS 02/12 CALCULATE THE LOCAL CONVECTIVE OVERTURN TIMESCALE AT THE
! BASE OF THE SURFACE CONVECTION ZONE, WITHIN THE ENVELOPE. THIS CODE
! SNAGGED FROM MIDMOD; unlike TAUINT (the analogous interior-model
! routine), this version deals with the possibility of two or more
! surface convection zones separated by a radiative gap within the
! envelope.
subroutine taucal(delta_mass, shell_mass, convective_flag, log10_radius, &
     log10_pressure, log10_density, local_gravity, num_points, &
     convective_velocity, radiative_gradient, adiabatic_gradient)
      use const_lib
      use star_info_lib, only: star, json
      implicit none

      double precision, intent(in) :: delta_mass(json), shell_mass(json)
      logical, intent(in) :: convective_flag(json)
      double precision, intent(inout) :: log10_radius(json), &
           log10_pressure(json), log10_density(json), local_gravity(json)
      integer, intent(in) :: num_points
      double precision, intent(in) :: convective_velocity(json)
      double precision, intent(in) :: radiative_gradient(json), &
           adiabatic_gradient(json)
! --- locals ---
      logical :: fully_convective_flag, exiting_cz_flag, two_cz_flag, &
           cz_exists_flag, surface_cz_deep_enough
      integer :: i, cz_base_index, upper_cz_top_index, j, k
      double precision :: mass_at_cz_top, mass_at_cz_bottom
      double precision :: dd2, dd1, interp_fraction
      double precision :: log10_radius_bcz, radius_bcz
      double precision :: pressure_scale_height_upper, &
           pressure_scale_height_lower, pressure_scale_height_bcz
      double precision :: log10_radius_test, match_tolerance
      double precision :: avg_convective_velocity, mass_sum, &
           velocity_mass_sum, cz_width, convective_velocity_bcz



! JVS 02/12 CALCULATE THE LOCAL CONVECTIVE OVERTURN TIMESCALE AT THE BASE
! OF THE CZ. IN OLDER VERSIONS THIS WAS ONLY DONE FOR ROTATING MODELS;
! THIS MAKES IT SO TAUCZ IS CALCULATED FOR ALL MODELS. THIS PARTICULAR SUBROUTINE
! DOES CALCULATES THE OVERTURN TIMESCALES IN MODELS WHERE THE CZ IS WITHIN
! THE ENVELOPE

! THIS CODE SNAGGED FROM MIDMOD

!  DETERMINE EXTENT OF SURFACE CONVECTION ZONE.
      fully_convective_flag = .false.
! JVS Allows the last point to be stable.
! JVS 05/14 removed requirement for last point to be convective, since this routine
! is only called from ENVINT, the near surface points in the envelope are often
! radiative, despite the fact that a convective envelope exists. Code now checks
! to make sure there is a CZ somewhere in the envelope.
!      IF(LCZ(M) .OR. LCZ(M-1) .OR. LCZ(M-2))THEN



!  SURFACE C.Z. EXISTS.  FIND LOWEST SHELL (IMAX), WHICH IS ALSO THE
!  UPPERMOST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONALLY INDUCED MIXING.
!         DO 71 I = M-1,1,-1
!   81    IMAX = I + 1

! DEAL THE POSSIBILITY OF TWO OR MORE SURFACE CONVECTION ZONES IN THE ENVELOPE
        exiting_cz_flag = .false.
        two_cz_flag = .false.
        cz_exists_flag = .false. ! Flag for the existence of a CZ
        do i = num_points-1,1,-1
            if (.not.convective_flag(i) .and. .not. exiting_cz_flag .and. cz_exists_flag) then ! EXITING CZ
                  cz_base_index = i+1
                  exiting_cz_flag = .true.
            else if (exiting_cz_flag .and. convective_flag(i)) then ! TWO CZS WITH RADIATIVE ZONE BETWEEN
                  exiting_cz_flag = .false.
                  two_cz_flag = .true.
                  upper_cz_top_index = i
            else if( convective_flag(i) .and. .not. cz_exists_flag) then !First convective point from surface
                  cz_exists_flag = .true.
            endif
        end do
         if (cz_base_index .lt. 1) then
              fully_convective_flag = .true.
         endif

      if (cz_exists_flag) then ! REPLACES "IF(LCZ(M)" THEN JVS 05/14


!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
         mass_at_cz_top = shell_mass(num_points)
         if (cz_base_index.gt.1) then
            mass_at_cz_bottom = 0.5d0*(shell_mass(cz_base_index)+shell_mass(cz_base_index-1))
         else
            mass_at_cz_bottom = shell_mass(1)
         endif
!  LCZSUR=T IF A SURFACE C.Z.DEEP ENOUGH FOR ANGULAR MOMENTUM LOSS EXISTS
         if ((mass_at_cz_top-mass_at_cz_bottom)/star%solar_mass_cgs.gt.0.0d0) then
            surface_cz_deep_enough = .true.
         else
            surface_cz_deep_enough= .false.
         endif
      else
!  NO SURFACE C.Z.
         cz_base_index = num_points
         surface_cz_deep_enough = .false.
      endif

! JVS 10/11/13 Concerned that this is inappropriate for LCZSUR = FALSE.
! Should never be the case, but removed it.
!      IF(LCZSUR)THEN
         if (.not.fully_convective_flag) then
            if (.not.rotation_active) then
               local_gravity(cz_base_index)=shell_mass(cz_base_index)*exp(ln10*(cgl-2.0d0*log10_radius(cz_base_index)))
               local_gravity(cz_base_index-1)=shell_mass(cz_base_index-1)*exp(ln10*(cgl-2.0d0*log10_radius(cz_base_index-1)))
            endif
! PINPOINT RCZ
            dd2 = radiative_gradient(cz_base_index-1)-adiabatic_gradient(cz_base_index-1)
            dd1 = radiative_gradient(cz_base_index)-adiabatic_gradient(cz_base_index)
            interp_fraction = dd2/(dd2-dd1)
! INFER HP
            log10_radius_bcz = log10_radius(cz_base_index-1)+interp_fraction* &
                 (log10_radius(cz_base_index)-log10_radius(cz_base_index-1))
            radius_bcz = exp(ln10*log10_radius_bcz)
            pressure_scale_height_upper = exp(ln10*(log10_pressure(cz_base_index)- &
                 log10_density(cz_base_index)))/local_gravity(cz_base_index)
            pressure_scale_height_lower = exp(ln10*(log10_pressure(cz_base_index-1)- &
                 log10_density(cz_base_index-1)))/local_gravity(cz_base_index-1)
            pressure_scale_height_bcz = pressure_scale_height_lower + &
                 interp_fraction*(pressure_scale_height_upper-pressure_scale_height_lower)
            log10_radius_test = dlog10(radius_bcz+pressure_scale_height_bcz)
! JVS 03/14 WHEN THERE ARE TWO CZS, THE DEEPER ONE CAN BE ~HP THICK. IN
! THIS CASE, WE"LL WANT AND AVERAGE CONVECTIVE VELOCITY ACROSS THE REGION
!--------------------------------------------------------------
            match_tolerance = 0.5
!            IF(ABS(1.0 - ((DEXP(CLN*(RTESTL-HR(IMAX))) - 1.0)/
!     *         (DEXP(CLN*(HR(ICZTOP)-HR(IMAX))) - 1.0))) .LT. TOL
!     *          .AND. LTWOCZ) THEN
            if ((abs(1.0 - ((dexp(ln10*(log10_radius_test-log10_radius(cz_base_index))) - 1.0)/ &
               (dexp(ln10*(log10_radius(upper_cz_top_index)-log10_radius(cz_base_index))) - 1.0))) &
               .lt. match_tolerance) .or. two_cz_flag) then
                  ! TAKE AVERAGE CONV VELOCITY
                  avg_convective_velocity = 0.0
!                  DENOM = ABS(HS1(ICZTOP)-HS1(IMAX))/DEXP(CLN*STOTAL) ! MASS IN CZ
                  mass_sum = 0.0
                  velocity_mass_sum = 0.0
                  cz_width = abs(dexp(ln10*log10_radius(upper_cz_top_index))- &
                       dexp(ln10*log10_radius(cz_base_index)))
                  do j=cz_base_index,upper_cz_top_index,1
                        mass_sum = mass_sum + delta_mass(j)
                        velocity_mass_sum = velocity_mass_sum+ 0.5d0* &
                             (convective_velocity(j)+convective_velocity(j+1))*delta_mass(j)
                  enddo
                  avg_convective_velocity = velocity_mass_sum/mass_sum
                  star%turnover%convective_turnover_timescale = cz_width/avg_convective_velocity
            else
!             (ORIGINAL ROUTINE)
!             FIND V
                       do k = cz_base_index+1,num_points
                           if (log10_radius(k).gt.log10_radius_test) then
                                      interp_fraction = (log10_radius_test-log10_radius(k-1))/ &
                                           (log10_radius(k)-log10_radius(k-1))
                                      convective_velocity_bcz = convective_velocity(k-1)+ &
                                           interp_fraction*(convective_velocity(k)-convective_velocity(k-1))
                                      exit
                           endif
                    end do
                       if (k > num_points) then
                    convective_velocity_bcz = convective_velocity(num_points)
                       end if
!                  DEFINE TAUCZ
                  star%turnover%convective_turnover_timescale = pressure_scale_height_bcz/convective_velocity_bcz
            endif

! JVS 10/11/13
!         ELSE
!C INFER HP
!C HP < R AT THE FIRST POINT.  ASSUME V CONSTANT INSIDE AND HP = K/R FOR
!C SLOWLY VARYING DENSITY AND PRESSURE NEAR THE CENTER.
!C FIND LOCATION WHERE HP = R
!                  IF(PSCA2.LE.RTEST2)THEN
!                     FX = (RTEST1-PSCA1)/((PSCA2-RTEST2)-(PSCA1-RTEST1))
!C FIND V
!                     CVEL = SVEL(K-1)+FX*(SVEL(K)-SVEL(K-1))
!                     PSCA = PSCA1+FX*(PSCA2-PSCA1)
!C DEFINE TAUCZ
!C 95            CONTINUE
!            ENDIF
!         ENDIF
      else
         star%turnover%convective_turnover_timescale = 0.0d0
      endif
      print*, 'TauCal'

! END JVS

      return
end subroutine taucal
