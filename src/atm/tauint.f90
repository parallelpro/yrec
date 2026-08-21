!----------------------------------------------------------------------
! tauint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original tauint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! JVS 02/12 CALCULATE THE LOCAL CONVECTIVE OVERTURN TIMESCALE AT THE
! BASE OF THE SURFACE CONVECTION ZONE, for the interior model (as
! opposed to TAUCAL, which does the analogous calculation within the
! envelope). If the calculation overshoots the last interior point, it
! is stitched onto the envelope structure (common/envstruct/) for more
! room. THIS CODE SNAGGED FROM MIDMOD.
subroutine tauint(shell_mass, convective_flag, log10_radius, &
     log10_pressure, log10_density, local_gravity, num_points, &
     convective_velocity, radiative_gradient, adiabatic_gradient)
      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: shell_mass(json)
      logical, intent(in) :: convective_flag(json)
      double precision, intent(inout) :: log10_radius(json), &
           log10_pressure(json), log10_density(json), local_gravity(json)
      integer, intent(in) :: num_points
      double precision, intent(in) :: convective_velocity(json)
      double precision, intent(in) :: radiative_gradient(json), &
           adiabatic_gradient(json)

! common/comp/: not used in this file. Naming matches getopac.f90.
      double precision :: envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv(12), xnew, znew, stotal, senv
      common/comp/envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv
! common/deuter/: not used in this file. Naming matches mix.f90.
      double precision :: deuterium_burning_rate(json), &
           deuterium_burning_rate_start(json), accreted_mass_fraction
      integer :: jcz
      common/deuter/deuterium_burning_rate, deuterium_burning_rate_start, &
           accreted_mass_fraction, jcz
! common/rot/: only rotation_active is used here. Naming matches
! mixcz.f90.
      double precision :: wnew, walpcz, acfpft
      integer :: itfp1, itfp2
      logical :: rotation_active, instability_transport_active, lwnew
      common/rot/wnew, walpcz, acfpft, itfp1, itfp2, rotation_active, &
           instability_transport_active, lwnew
! common/const/: only solar_mass_cgs is used here. Naming matches
! wrtout.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
! common/jtest/: not used in this file. Naming is local to this batch.
      integer :: imax1_placeholder, imax2_placeholder
      logical :: ljvs_placeholder
      common/jtest/imax1_placeholder, imax2_placeholder, ljvs_placeholder
! Needs access to this common block: grants knowledge of envellope
! KC 2025-05-30 reordered common block elements
! JvS 08/25 Updated with new elements
! common/envstruct/: only env_log10_radius/env_convective_velocity/
! num_env_points are used here (stitching the envelope on for extra
! room when the interior search overshoots). Naming is local to this
! batch.
      double precision :: env_log10_pressure(json), env_log10_temperature(json), &
           env_log10_mass(json), env_log10_density(json), env_log10_radius(json), &
           env_hydrogen_fraction(json), env_metal_fraction(json)
      logical :: env_convective_flag(json)
      double precision :: env_gradients(3,json), env_convective_velocity(json), &
           env_beta(json)
      double precision :: env_gamma1(json), env_specific_heat_cp(json), &
           env_ion_fraction(3,json)
      double precision :: env_opacity(json), env_luminosity(json), &
           env_dlnrho_dlnt(json)
      integer :: num_env_points
      common/envstruct/env_log10_pressure, env_log10_temperature, &
           env_log10_mass, env_log10_density, env_log10_radius, &
           env_hydrogen_fraction, env_metal_fraction, env_convective_flag, &
           env_gradients, env_convective_velocity, env_beta, &
           env_gamma1, env_specific_heat_cp, env_ion_fraction, &
           env_opacity, env_luminosity, env_dlnrho_dlnt, num_env_points
! G Somers 3/17, ADDING NEW TAUCZ COMMON BLOCK
! common/ovrtrn/: only convective_turnover_timescale is set here.
! Naming matches mixcz.f90.
      logical :: use_new_turnover_timescale, calc_envelope_flag
      double precision :: convective_turnover_timescale, &
           convective_turnover_timescale_old, pphot, pphot0, fracstep
      common/ovrtrn/use_new_turnover_timescale, calc_envelope_flag, &
           convective_turnover_timescale, convective_turnover_timescale_old, &
           pphot, pphot0, fracstep

      save

! --- locals ---
      logical :: fully_convective_flag, surface_cz_deep_enough
      integer :: i, cz_base_index, k
      double precision :: mass_at_cz_top, mass_at_cz_bottom
      double precision :: dd2, dd1, interp_fraction
      double precision :: log10_radius_bcz, radius_bcz
      double precision :: pressure_scale_height_upper, &
           pressure_scale_height_lower, pressure_scale_height_bcz
      double precision :: log10_radius_test, convective_velocity_bcz
      double precision :: pressure_scale_height2, pressure_scale_height1, &
           radius_test2, radius_test1


! JVS 02/12 CALCULATE THE LOCAL CONVECTIVE OVERTURN TIMESCALE AT THE BASE
! OF THE CZ. IN OLDER VERSIONS THIS WAS ONLY DONE FOR ROTATING MODELS;
! THIS MAKES IT SO TAUCZ IS CALCULATED FOR ALL MODELS.
! THIS CODE SNAGGED FROM MIDMOD
!  DETERMINE EXTENT OF SURFACE CONVECTION ZONE.
      fully_convective_flag = .false.
! JVS Allows the last point to be stable.
      if (convective_flag(num_points) .or. convective_flag(num_points-1)) then
!  SURFACE C.Z. EXISTS.  FIND LOWEST SHELL (IMAX), WHICH IS ALSO THE
!  UPPERMOST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONALLY INDUCED MIXING.
         do 71 i = num_points-1,1,-1
            if (.not.convective_flag(i)) goto 81
   71    continue
         fully_convective_flag = .true.
         i = 0
   81    cz_base_index = i + 1
!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
         mass_at_cz_top = shell_mass(num_points)
         if (cz_base_index.gt.1) then
            mass_at_cz_bottom = 0.5d0*(shell_mass(cz_base_index)+shell_mass(cz_base_index-1))
         else
            mass_at_cz_bottom = 0.0d0
         endif
!  LCZSUR=T IF A SURFACE C.Z.DEEP ENOUGH FOR ANGULAR MOMENTUM LOSS EXISTS
         if ((mass_at_cz_top-mass_at_cz_bottom)/solar_mass_cgs.gt.0.0d0) then
            surface_cz_deep_enough = .true.
         else
            surface_cz_deep_enough= .false.
         endif
      else
!  NO SURFACE C.Z.
         cz_base_index = num_points
         surface_cz_deep_enough = .false.
      endif

      if (surface_cz_deep_enough) then
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
! FIND V
            do k = cz_base_index+1,num_points
               if (log10_radius(k).gt.log10_radius_test) then
                  interp_fraction = (log10_radius_test-log10_radius(k-1))/ &
                       (log10_radius(k)-log10_radius(k-1))
                  convective_velocity_bcz = convective_velocity(k-1)+ &
                       interp_fraction*(convective_velocity(k)-convective_velocity(k-1))
                  goto 85
               endif
            end do
            ! One pressure scale height overshoots edge of interior
            ! calculation. Stitch on the envelope for more room
            do k = 2,num_env_points
               if (env_log10_radius(k).gt.log10_radius_test .and. env_convective_velocity(k) .gt. 0.0) then
                  interp_fraction = (log10_radius_test-log10_radius(num_points))/ &
                       (env_log10_radius(k)-log10_radius(num_points))
                  convective_velocity_bcz = convective_velocity(num_points)+ &
                       interp_fraction*(env_convective_velocity(k)-convective_velocity(num_points))
                  goto 85
               endif
            end do
            convective_velocity_bcz = convective_velocity(num_points)
 85         continue
! DEFINE TAUCZ
            convective_turnover_timescale = pressure_scale_height_bcz/convective_velocity_bcz

         else
! INFER HP
            if (.not.rotation_active) then
               local_gravity(1)=shell_mass(1)*exp(ln10*(cgl-2.0d0*log10_radius(1)))
            endif
            pressure_scale_height2 = exp(ln10*(log10_pressure(1)-log10_density(1)))/local_gravity(1)
            radius_test2 = exp(ln10*log10_radius(1))
            if (pressure_scale_height2.le.radius_test2) then
! HP < R AT THE FIRST POINT.  ASSUME V CONSTANT INSIDE AND HP = K/R FOR
! SLOWLY VARYING DENSITY AND PRESSURE NEAR THE CENTER.
               convective_velocity_bcz = convective_velocity(1)
               pressure_scale_height_bcz = (pressure_scale_height2*radius_test2)**0.5d0
               convective_turnover_timescale = pressure_scale_height_bcz/convective_velocity_bcz
            else
               do k = 2,num_points
                  pressure_scale_height1 = pressure_scale_height2
                  radius_test1 = radius_test2
                  if (.not.rotation_active) then
                     local_gravity(k)=shell_mass(k)*exp(ln10*(cgl-2.0d0*log10_radius(k)))
                  endif
                  pressure_scale_height2 = exp(ln10*(log10_pressure(k)-log10_density(k)))/local_gravity(k)
                  radius_test2 = exp(ln10*log10_radius(k))
! FIND LOCATION WHERE HP = R
                  if (pressure_scale_height2.le.radius_test2) then
                     interp_fraction = (radius_test1-pressure_scale_height1)/ &
                          ((pressure_scale_height2-radius_test2)-(pressure_scale_height1-radius_test1))
! FIND V
                     convective_velocity_bcz = convective_velocity(k-1)+ &
                          interp_fraction*(convective_velocity(k)-convective_velocity(k-1))
                     pressure_scale_height_bcz = pressure_scale_height1+ &
                          interp_fraction*(pressure_scale_height2-pressure_scale_height1)
! DEFINE TAUCZ
                     convective_turnover_timescale = pressure_scale_height_bcz/convective_velocity_bcz
                     goto 95
                  endif
               end do
               k = num_points
               convective_velocity_bcz = convective_velocity(num_points)
               pressure_scale_height_bcz = pressure_scale_height2
               convective_turnover_timescale = pressure_scale_height_bcz/convective_velocity_bcz
 95            continue
            endif
         endif
      else
         convective_turnover_timescale = 0.0d0
      endif

!--------------------------------------------------------------
!                  OPEN(UNIT=100,FILE='diagnostic.out',STATUS='OLD')
!                  DO 1505 I=1,M
!                  IF (LCZ(I)) HOLD = 1.0
!                  IF (.NOT.LCZ(I)) HOLD = 0.0
!                              WRITE(100,1504), HS2(I), HS1(I), HR(I),HOLD,
!     *                  DEL1(I), DEL2(I), SVEL(I),
!1505                  CONTINUE
!1504                  FORMAT(1X,7E16.8)
!                  CLOSE(100)
!----------------------------------------------------------------

! END JVS
      return
end subroutine tauint
