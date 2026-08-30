!$$$$$$
!----------------------------------------------------------------------
! midmod
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original midmod.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! In general a diffusion timestep is less than a model timestep. This
! subroutine takes the structure at the beginning and end of the
! timestep and assigns an intermediate structure using linear
! interpolation in time. Nuclear burning is accounted for here, and
! the changes in angular velocity due to changes in radius and
! conservation of angular momentum are also computed. MIDMOD also
! locates the outer edge of a central convection zone and the inner
! edge of a surface C.Z., and determines the mass and moment of
! inertia of the surface C.Z.
subroutine mid_timestep_model(full_timestep, sub_timestep, time_fraction, first_call, &
     eta_squared_mid, log_density_mid, hg_mid, moment_of_inertia_cz, &
     moment_of_inertia_mid, log_luminosity_mid, log_pressure_mid, &
     log_radius_mid, cz_mass_bottom, cz_mass_top, log_temperature_mid, &
     core_boundary_zone, envelope_boundary_zone, fully_convective_flag, &
     convective_flag_mid, am_transport_convective_flag_mid, &
     surface_cz_active, omega_mid, mean_radius_mid, qiw_mid, &
     radiative_zone_bounds, convective_zone_bounds, num_radiative_zones, &
     num_convective_zones, ierr)
      use rotation_scratch_lib
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, i_grad_actual, i_grad_ad, i_grad_rad, json
      use net_lib
      use phys_const_lib
      use burn_lib

      use math_lib
      implicit none

      double precision, intent(in) :: full_timestep, sub_timestep
      double precision, intent(inout) :: time_fraction
! LC(JSON),LCZ(JSON),ETA2M(JSON),HDM(JSON),HGM(JSON),  ! KC 2025-05-31
      logical, intent(in) :: first_call
      double precision, intent(out) :: eta_squared_mid(json), &
           log_density_mid(json), hg_mid(json)
      double precision, intent(out) :: moment_of_inertia_cz
      double precision, intent(out) :: moment_of_inertia_mid(json)
      double precision, intent(out) :: log_luminosity_mid(json), &
           log_pressure_mid(json)
! log_radius_mid (HRM) is read here before being (re)computed in this
! call, in the LFIRST=.false. branch below that seeds rot_scr%radius_prev from
! the previous MIDMOD call's log_radius_mid -- hence intent(inout)
! rather than intent(out).
      double precision, intent(inout) :: log_radius_mid(json)
      double precision, intent(out) :: cz_mass_bottom, cz_mass_top
      double precision, intent(out) :: log_temperature_mid(json)
      integer, intent(out) :: core_boundary_zone, envelope_boundary_zone
      logical, intent(out) :: fully_convective_flag
      logical, intent(out) :: convective_flag_mid(json)
      logical, intent(out) :: am_transport_convective_flag_mid(json)
      logical, intent(out) :: surface_cz_active
      double precision, intent(out) :: omega_mid(json)
      double precision, intent(out) :: mean_radius_mid(json), qiw_mid(json)
      integer, intent(out) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(out) :: num_radiative_zones, num_convective_zones




! MHP 8/17 added excen, c_2 to common block for Matt et al. 2012 cent. term
! former common/cwind/: not used in this file (only referenced in the
! commented-out taucz calculation below).
!      COMMON/CWIND/WMAX,EXMD,EXW,EXTAU,EXR,EXM,CONSTFACTOR,STRUCTFACTOR,LJDOT0















! MHP 06/02
      logical :: convective_state_changed(json)
! added vector for deuterium burning
      double precision :: deuterium_rate_mid(json), deuterium_rate_mid_start(json)
      integer :: num_species_tracked
      integer :: i, j, k, ii
      double precision :: step_fraction_ratio
      logical :: new_cz_detected
      double precision :: total_epsilon
      logical :: cz_change_active, redistribute_j_flag
      integer :: change_region_start, change_region_end
      integer :: cz_direction_flag
      integer :: cz_zone_bottom, cz_zone_top
      double precision :: cz_moment_of_inertia, cz_angular_momentum
      double precision :: convective_fraction
      double precision :: radius_interp
      double precision :: radius_shell_factor
      double precision :: specific_angular_momentum_corrected
      double precision :: delta_angular_momentum
      double precision :: angular_momentum_ratio
      double precision :: angular_momentum_check_ratio
      integer :: burn_zone_begin, burn_zone_end
      integer :: solid_zone_start, solid_zone_end

! IN GENERAL A DIFFUSION TIMESTEP IS LESS THAN A MODEL TIMESTEP.
! SR MIDMOD TAKES THE STRUCTURE AT THE BEGINNING AND END OF THE TIMESTEP
! AND ASSIGNS AN INTERMEDIATE STRUCTURE USING LINEAR INTERPOLATION IN TIME.
! NUCLEAR BURNING IS ACCOUNTED FOR HERE, AND THE CHANGES IN ANGULAR
! VELOCITY DUE TO CHANGES IN RADIUS AND CONSERVATION OF ANGULAR MOMENTUM
! ARE ALSO COMPUTED.
! MIDMOD ALSO LOCATES THE OUTER EDGE OF A CENTRAL CONVECTION ZONE AND THE INNER
! EDGE OF A SURFACE C.Z., AND DETERMINES THE MASS AND MOMENT OF INERTIA OF
! THE SURFACE C.Z.
      integer, intent(out) :: ierr

      ierr = 0

      num_species_tracked = 11
      if (star%job%use_extended_composition) num_species_tracked = 15
!  INITIALIZE COMPOSITION ARRAY THE FIRST TIME THROUGH.
!  HCOMPM IS THE ARRAY OF CHANGES IN COMPOSITION DUE TO NUCLEAR BURNING
!  DURING THE COURSE OF THE MODEL TIMESTEP.
! MHP 06/02 SAVE CHANGES IN R, CZ DEPTH, DELAD-DELRAD FROM THE
! PRIOR STEP.
      if (first_call) then
         do j = 1,star%nz
            rot_scr%convective_flag_prev(j) = star%convective_flag_start(j)
            rot_scr%radius_prev(j) = star%logR_start(j)
            rot_scr%del_grad_diff_prev(j) = rot_scr%old_del_adiabatic_mix(j)-rot_scr%old_del_radiative_mix(j)
            mix_scr%amum(j) = rot_scr%old_amu(j)
            do i = 1,num_species_tracked
               star%xa(i,j) = star%xa_start(i,j)
            end do
         end do
      else
         do j = 1,star%nz
            rot_scr%radius_prev(j) = log_radius_mid(j)
            rot_scr%del_grad_diff_prev(j) = rot_scr%del_grad_diff_new(j)
            rot_scr%convective_flag_prev(j) = convective_flag_mid(j)
         end do
      endif
      new_cz_detected = .false.
!  INTERPOLATE IN THE MODEL VARIABLES AND AUXILLARY PHYSICS.
      step_fraction_ratio = sub_timestep/full_timestep
      do j = 1,star%nz
         log_density_mid(j) = star%logRho_start(j) + &
              time_fraction*(star%logRho(j)-star%logRho_start(j))
         log_luminosity_mid(j) = star%luminosity_lsun_start(j) + &
              time_fraction*(star%luminosity_lsun(j)-star%luminosity_lsun_start(j))
         log_pressure_mid(j) = star%logP_start(j) + &
              time_fraction*(star%logP(j)-star%logP_start(j))
         log_radius_mid(j) = star%logR_start(j) + &
              time_fraction*(star%logR(j)-star%logR_start(j))
         log_temperature_mid(j) = star%logT_start(j) + &
              time_fraction*(star%logT(j)-star%logT_start(j))
!        DO 30 I = 1,IEND
         hg_mid(j) = star%old_hg(j) + time_fraction*(star%mean_gravity(j) - star%old_hg(j))
         mix_scr%del_adiabatic_mix(j) = rot_scr%old_del_adiabatic_mix(j) + &
              time_fraction*(star%grada(j)-rot_scr%old_del_adiabatic_mix(j))
         mix_scr%delm(j) = rot_scr%old_delm(j) + time_fraction*(star%gradT(j) - rot_scr%old_delm(j))
         mix_scr%del_radiative_mix(j) = rot_scr%old_del_radiative_mix(j) + &
              time_fraction*(star%gradr(j) - rot_scr%old_del_radiative_mix(j))
         mix_scr%esumm(j) = rot_scr%old_esum(j) + time_fraction*(star%eps_total(j) - rot_scr%old_esum(j))
         mix_scr%viscm(j) = rot_scr%old_visc(j) + time_fraction*(star%visc(j) - rot_scr%old_visc(j))
         mix_scr%thdifm(j) = rot_scr%old_thdif(j) + time_fraction*(star%thdif(j) - rot_scr%old_thdif(j))
         mix_scr%cpm(j) = rot_scr%old_cp(j) + time_fraction*(star%cp(j) - rot_scr%old_cp(j))
         mix_scr%qdtm(j) = rot_scr%old_qdt(j) + time_fraction*(star%qdt(j) - rot_scr%old_qdt(j))
         mix_scr%om(j) = rot_scr%old_om(j) + time_fraction*(star%opacity_zone(j) - rot_scr%old_om(j))
         mix_scr%amum(j) = mix_scr%amum(j) + step_fraction_ratio*(star%mu(j) - rot_scr%old_amu(j))
! MHP 6/00 ADDED TOTAL ENERGY GENERATION
         total_epsilon = star%eps_total(j)+star%eps_channels(i_eps_neu,j)+star%eps_channels(i_eps_grav,j)
         mix_scr%epsm(j) = rot_scr%old_eps(j)+time_fraction*(total_epsilon-rot_scr%old_eps(j))
      end do
!  CHECK FOR ADVANCING OR RECEDING CONVECTIVE REGIONS.USE INTERPOLATED
!  RADIATIVE AND ADIABATIC TEMPERATURE GRADIENTS TO DETERMINE WHETHER
!  OR NOT A ZONE IS CONVECTIVE IF IT CHANGES STATE OVER THE COURSE OF A
!  TIMESTEP.
      do i = 1,star%nz
         rot_scr%del_grad_diff_new(i) = mix_scr%del_adiabatic_mix(i)-mix_scr%del_radiative_mix(i)
         if (star%convective_flag(i).eqv.rot_scr%convective_flag_prev(i)) then
            convective_flag_mid(i) = star%convective_flag(i)
            mix_scr%velm(i) = rot_scr%old_vel(i) + time_fraction*(star%conv_vel(i)-rot_scr%old_vel(i))
            convective_state_changed(i) = .false.
         else
            convective_state_changed(i) = .true.
            new_cz_detected = .true.
            if (mix_scr%del_adiabatic_mix(i).lt.mix_scr%del_radiative_mix(i)) then
               convective_flag_mid(i) = .true.
               mix_scr%velm(i) = max(rot_scr%old_vel(i),star%conv_vel(i))
            else
               convective_flag_mid(i) = .false.
               mix_scr%velm(i) = 0.0D0
            endif
         endif
      end do
      new_cz_detected = .false.
! MHP 06/02 IF THE CZ DEPTHS HAVE CHANGED, RESOLVE
! OUT WHEN A GIVEN ZONE "DROPPED OUT" OF THE CONVECTION
! ZONE OR ENTERED IT.  MODIFY THE SPECIFIC ANGULAR MOMENTUM
! ACCORDINGLY.
      if (new_cz_detected) then
         cz_change_active = .false.
         redistribute_j_flag = .false.
         do i = 1,star%nz-1
            if (convective_state_changed(i)) then
               if (.not.cz_change_active) then
! START OF CHANGED REGION
                  cz_change_active = .true.
                  change_region_start = i
               endif
            endif
            if (.not.convective_state_changed(i)) then
               if (cz_change_active) then
! END OF CHANGED REGION
               change_region_end = i - 1
               cz_change_active = .false.
               redistribute_j_flag = .true.
               write(*,*)change_region_start,change_region_end, &
                    rot_scr%del_grad_diff_new(change_region_start), &
                    rot_scr%del_grad_diff_new(change_region_start-1), &
                    rot_scr%del_grad_diff_new(change_region_end+1)
! ONLY CHANGE IF THERE IS A RECEDING CONVECTION ZONE;
! WE CAN ASSUME LOCAL CONSERVATION OF J PRIOR TO HAVING
! A SHELL INCORPORATED INTO A CZ.
               if (rot_scr%del_grad_diff_new(change_region_start).gt.0.0D0) then
                  if (rot_scr%del_grad_diff_new(change_region_start-1).lt.0.0D0) then
! SUBTRACTED FROM A CORE CZ
                     cz_direction_flag = -1
                  else if (rot_scr%del_grad_diff_new(change_region_end+1).lt.0.0D0) then
! SUBTRACTED FROM AN ENVELOPE CZ
                     cz_direction_flag = 1
                  else
                     redistribute_j_flag = .false.
                  endif
               else
                  redistribute_j_flag = .false.
               endif
               endif
            endif
            if (redistribute_j_flag) then
               redistribute_j_flag = .false.
! CORE CZ
               if (cz_direction_flag.lt.0) then
! LOCATE EDGES OF CURRENT LOWER CZ
                  cz_zone_top = change_region_start - 1
                  do ii = change_region_start-1,1,-1
                     if (.not.convective_flag_mid(ii)) then
                        cz_zone_bottom = ii + 1
                        exit
                     endif
                  end do
                  if (ii .lt. 1) cz_zone_bottom = 1
                  do j = change_region_end,change_region_start,-1
                     if (rot_scr%del_grad_diff_new(j).lt.0.0D0) exit
                     cz_moment_of_inertia = 0.0D0
                     cz_angular_momentum = 0.0D0
! FRACTION OF THE TIMESTEP SHELL WAS RADIATIVE
                     time_fraction = rot_scr%del_grad_diff_new(j)/ &
                          (rot_scr%del_grad_diff_new(j)-rot_scr%del_grad_diff_prev(j))
! FRACTION OF THE TIMESTEP THAT THE SHELL WAS CONVECTIVE
                     convective_fraction = 1.0D0 - time_fraction
                     do ii = cz_zone_bottom,j-1
                        radius_interp = rot_scr%radius_prev(ii)+ &
                             convective_fraction*(star%logR(ii)-star%logR_start(ii))
                        cz_moment_of_inertia = cz_moment_of_inertia+ &
                             cc23*star%dm(ii)*exp10((2.0D0*radius_interp))
                        cz_angular_momentum = cz_angular_momentum + &
                             star%j_rot(ii)*star%dm(ii)
                     end do
! GET THE CZ MOMENT OF INERTIA AND TOTAL J AT THE TIME THE
! SHELL WAS RELEASED; THIS GIVES A CORRECTED VALUE FOR
! J/M.  REDISTRIBUTE THE DIFFERENCE IN J BACK INTO THE CZ.
                     radius_shell_factor = 2.0D0*(rot_scr%radius_prev(j)+ &
                          convective_fraction*(star%logR(j)-rot_scr%radius_prev(j)))
                     specific_angular_momentum_corrected = &
                          cc23*(cz_angular_momentum/cz_moment_of_inertia)* &
                          exp10(radius_shell_factor)
                     delta_angular_momentum = &
                          (star%j_rot(j)- &
                          specific_angular_momentum_corrected)*star%dm(j)
                     angular_momentum_ratio = &
                          (cz_angular_momentum-delta_angular_momentum)/ &
                          cz_angular_momentum
                     star%j_rot(j)=specific_angular_momentum_corrected
                     angular_momentum_check_ratio = &
                          specific_angular_momentum_corrected/ &
                          star%j_rot(j)
                     write(*,*)j,time_fraction,angular_momentum_check_ratio, &
                          delta_angular_momentum,angular_momentum_ratio
                     do ii = cz_zone_bottom,j-1
                        star%j_rot(ii) = &
                             star%j_rot(ii)*angular_momentum_ratio
                     end do
                  end do
               else
! ENVELOPE CZ
! LOCATE EDGES OF CURRENT LOWER CZ
                  cz_zone_bottom = change_region_end + 1
                  do ii = change_region_end+2,star%nz
                     if (.not.convective_flag_mid(ii)) then
                        cz_zone_top = ii - 1
                        exit
                     endif
                  end do
                  if (ii .gt. star%nz) cz_zone_top = star%nz
                  do j = change_region_start,change_region_end
                     if (rot_scr%del_grad_diff_new(j).lt.0.0D0) exit
                     cz_moment_of_inertia = 0.0D0
                     cz_angular_momentum = 0.0D0
! FRACTION OF THE TIMESTEP SHELL WAS RADIATIVE
                     time_fraction = rot_scr%del_grad_diff_new(j)/ &
                          (rot_scr%del_grad_diff_new(j)-rot_scr%del_grad_diff_prev(j))
! FRACTION OF THE TIMESTEP THAT THE SHELL WAS CONVECTIVE
                     convective_fraction = 1.0D0 - time_fraction
                     do ii = j+1,cz_zone_top
                        radius_interp = rot_scr%radius_prev(ii)+ &
                             convective_fraction*(star%logR(ii)-star%logR_start(ii))
                        cz_moment_of_inertia = cz_moment_of_inertia+ &
                             cc23*star%dm(ii)*exp10((2.0D0*radius_interp))
                        cz_angular_momentum = cz_angular_momentum + &
                             star%j_rot(ii)*star%dm(ii)
                     end do
! GET THE CZ MOMENT OF INERTIA AND TOTAL J AT THE TIME THE
! SHELL WAS RELEASED; THIS GIVES A CORRECTED VALUE FOR
! J/M.  REDISTRIBUTE THE DIFFERENCE IN J BACK INTO THE CZ.
                     radius_shell_factor = 2.0D0*(rot_scr%radius_prev(j)+ &
                          convective_fraction*(star%logR(j)-rot_scr%radius_prev(j)))
                     specific_angular_momentum_corrected = &
                          cc23*(cz_angular_momentum/cz_moment_of_inertia)* &
                          exp10(radius_shell_factor)
                     delta_angular_momentum = &
                          (star%j_rot(j)- &
                          specific_angular_momentum_corrected)*star%dm(j)
                     angular_momentum_ratio = &
                          (cz_angular_momentum-delta_angular_momentum)/ &
                          cz_angular_momentum
                     star%j_rot(j)=specific_angular_momentum_corrected
                     angular_momentum_check_ratio = &
                          specific_angular_momentum_corrected/ &
                          star%j_rot(j)
                     write(*,*)j,time_fraction,angular_momentum_check_ratio, &
                          delta_angular_momentum,angular_momentum_ratio
                     do ii = j+1,star%nz
                        star%j_rot(ii) = &
                             star%j_rot(ii)*angular_momentum_ratio
                     end do
                  end do
               endif
            endif
         end do
      endif
!  CONVECTIVE OVERSHOOT APPLIED TO NORMAL CONVECTION ZONES.
      call am_convective_regions(star%xa,log_density_mid,log_pressure_mid,log_radius_mid, &
           star%log_mass,log_temperature_mid,convective_flag_mid,star%nz, &
           am_transport_convective_flag_mid,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones, ierr)
      if (ierr /= 0) return
!  REBURN THE ORIGINAL MIXTURE FOR THE SMALL DIFFUSION TIME STEP.
! MHP 6/00 - FOR THERMODYNAMIC CONSISTENCY USE CURRENT T, NOT END OF
! TIMESTEP T (REPLACED HT IN CALL WITH HTM)
!       CALL ROTMIX(DT,HCOMP,HS2,HTM,5,M,MRZONE,MXZONE,  ! KC 2025-05-31
      call rotmix(sub_timestep,star%xa,star%dm,log_temperature_mid, &
           star%nz,radiative_zone_bounds,convective_zone_bounds, &
           num_radiative_zones,num_convective_zones &
           ,star%log_total_mass,log_density_mid,star%log_mass,log_radius_mid, &
           log_pressure_mid,am_transport_convective_flag_mid,star%m, ierr)
      if (ierr /= 0) return
      do i = 1,star%nz
         do j = 1,num_species_tracked
            rot_scr%composition_snapshot(j,i) = star%xa(j,i)
         end do
      end do
! MHP 05/02 ADDED DEUTERIUM BURNING
      if (star%xa_start(12,star%nz).gt.1.0D-14) then
! INCREMENT THE TIMESTEP
         if (first_call) then
            do i = 1,star%nz
               deuterium_rate_mid(i) = star%deuterium_burning_rate_start(i)+ &
                    step_fraction_ratio*(star%deuterium_burning_rate(i)- &
                    star%deuterium_burning_rate_start(i))
               deuterium_rate_mid_start(i) = star%deuterium_burning_rate_start(i)
            end do
         else
            do i = 1,star%nz
               deuterium_rate_mid_start(i) = deuterium_rate_mid(i)
               deuterium_rate_mid(i) = deuterium_rate_mid(i)+ &
                    step_fraction_ratio*(star%deuterium_burning_rate(i)- &
                    star%deuterium_burning_rate_start(i))
            end do
         endif
! RADIATIVE ZONES.
!
         do k = 1,num_radiative_zones
            do j = radiative_zone_bounds(k,1),radiative_zone_bounds(k,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (log_temperature_mid(j).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
               burn_zone_begin = j
               burn_zone_end = j
               call dburnm(burn_zone_begin,burn_zone_end,star%nz,star%dm, &
                    star%xa,sub_timestep,deuterium_rate_mid, &
                    deuterium_rate_mid_start,step_fraction_ratio)
           end do
        end do
         if (k > num_radiative_zones) then
         end if
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
         do k = 1,num_convective_zones
            burn_zone_begin = convective_zone_bounds(k,1)
            burn_zone_end = convective_zone_bounds(k,2)
            call dburnm(burn_zone_begin,burn_zone_end,star%nz,star%dm, &
                 star%xa,sub_timestep,deuterium_rate_mid, &
                 deuterium_rate_mid_start,step_fraction_ratio)
         end do
      endif
!  DETERMINE EXTENT OF CENTRAL CONVECTION ZONE.
!  IMIN IS THE FIRST ZONE ABOVE A CENTRAL CONVECTION ZONE, AND THUS THE
!  FIRST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONAL INSTABILITIES.
      if (am_transport_convective_flag_mid(1)) then
         do i = 2,star%nz
            if (.not.am_transport_convective_flag_mid(i)) exit
         end do
         if (i > (star%nz)) then
         i = star%nz + 1
         end if
         core_boundary_zone = max(2,i-1)
      else
         core_boundary_zone = 2
      endif
!  DETERMINE EXTENT OF SURFACE CONVECTION ZONE.
      fully_convective_flag = .false.
      if (am_transport_convective_flag_mid(star%nz)) then
!  SURFACE C.Z. EXISTS.  FIND LOWEST SHELL (IMAX), WHICH IS ALSO THE
!  UPPERMOST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONALLY INDUCED MIXING.
         do i = star%nz-1,1,-1
            if (.not.am_transport_convective_flag_mid(i)) exit
         end do
         if (i < (1)) then
         fully_convective_flag = .true.
         i = 0
         end if
         envelope_boundary_zone = i + 1
!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
         cz_mass_top = exp(ln10*star%log_total_mass)
         if (envelope_boundary_zone.gt.1) then
            cz_mass_bottom = 0.5D0*(star%m(envelope_boundary_zone)+ &
                 star%m(envelope_boundary_zone-1))
         else
            cz_mass_bottom = 0.0D0
         endif
!  LCZSUR=T IF A SURFACE C.Z.DEEP ENOUGH FOR ANGULAR MOMENTUM LOSS EXISTS
         if ((cz_mass_top-cz_mass_bottom)/star%solar_mass_cgs.gt.0.0D0) then
            surface_cz_active = .true.
         else
            surface_cz_active= .false.
         endif
      else
!  NO SURFACE C.Z.
         envelope_boundary_zone = star%nz
         surface_cz_active = .false.
      endif
! MHP 3/09 ADD ABILITY TO COMPUTE THE CONVECTIVE OVERTURN TIMESCALE 'LOCALLY'.
! THIS IS DEFINED AS HP/V, WHERE HP IS THE PRESSURE SCALE HEIGHT AT THE BASE
! AND V IS THE CONVECTIVE VELOCITY ONE PRESSURE SCALE HEIGHT ABOVE THE BASE.
! FOR FULLY CONVECTIVE STARS THE CENTRAL PRESSURE SCALE HEIGHT DIVERGES, SO
! INSTEAD DEFINE V AT THE POINT WHERE HP = P/RHO*G = P*R^2/RHO*G*M = R.
!      IF(LCZSUR)THEN
!         IF(.NOT.LALLCZ)THEN
!C PINPOINT RCZ
!C            WRITE(*,911)IMAX,FX,DD2,DD1
!C 911        FORMAT(I5,1P3E12.3)
!C INFER HP
!C FIND V
!            DO K = IMAX+1,M
!C 85         CONTINUE
!C DEFINE TAUCZ
!            TAUCZ = PSCA/CVEL
!C            WRITE(*,911)K,PSCA,CVEL,TAUCZ
!         ELSE
!C INFER HP
!C HP < R AT THE FIRST POINT.  ASSUME V CONSTANT INSIDE AND HP = K/R FOR
!C SLOWLY VARYING DENSITY AND PRESSURE NEAR THE CENTER.
!C               WRITE(*,912)PSCA2,RTEST2,PSCA,CVEL,TAUCZ
!C 912           FORMAT(1P5E12.3)
!C FIND LOCATION WHERE HP = R
!                  IF(PSCA2.LE.RTEST2)THEN
!                     FX = (RTEST1-PSCA1)/((PSCA2-RTEST2)-(PSCA1-RTEST1))
!C FIND V
!                     CVEL = VELM(K-1)+FX*(VELM(K)-VELM(K-1))
!                     PSCA = PSCA1+FX*(PSCA2-PSCA1)
!C DEFINE TAUCZ
!C 95            CONTINUE
!C               WRITE(*,911)K,PSCA,CVEL,TAUCZ
! JNT 09/25 FOR 05/15 IMPJMOD=1 IS THE SAME AS LSOLID
      if (.not.star%ctrl%force_solid_body_rotation .and. (star%ctrl%solid_body_mode_flag.ne.1)) then
!  NOW FIND THE RUN OF ROTATION VARIABLES THAT ARE CONSISTENT WITH THE
!  INTERMEDIATE STRUCTURE AND THE RUN OF SPECIFIC ANGULAR MOMENTUM J/M.
!  J/M = I/M * OMEGA.
!  FIRST GUESS AT MOMENT OF INERTIA OF DISTORTED SPERICAL SHELLS:
!  I/M = 2/3 R**2.
         do i = 1,star%nz
            moment_of_inertia_mid(i) = cc23*exp(ln10*2.0D0*log_radius_mid(i))
            omega_mid(i) = star%j_rot(i)/moment_of_inertia_mid(i)
            moment_of_inertia_mid(i) = star%dm(i)*moment_of_inertia_mid(i)
         end do
!  SOLVE FOR THE ANGULAR VELOCITIES OF THE SHELLS GIVEN THE SPECIFIC
!  ANGULAR MOMENTUM (HJM) AND A FIRST GUESS AT THE ANGULAR VELOCITY(OMEGAM)
!  AND MOMENT OF INERTIA (HIM).
         call omega_from_j(log_density_mid,star%j_rot,log_radius_mid, &
              star%log_mass,star%dm,am_transport_convective_flag_mid,star%nz, &
              eta_squared_mid,moment_of_inertia_mid,omega_mid,qiw_mid, &
              mean_radius_mid)
!  FIND TOTAL MOMENT OF INERTIA OF THE SURFACE C.Z. IF APPLICABLE.
         if (surface_cz_active) then
            moment_of_inertia_cz = moment_of_inertia_mid(star%nz)
            do i = star%nz-1,1,-1
               if (.not.am_transport_convective_flag_mid(i)) exit
               moment_of_inertia_cz = moment_of_inertia_cz + moment_of_inertia_mid(i)
            end do
         endif
      else
         solid_zone_start = 1
         solid_zone_end = star%nz
         call solid_body_omega(log_density_mid,star%j_rot,log_radius_mid, &
              star%log_mass,star%dm,solid_zone_start,solid_zone_end, &
              eta_squared_mid,moment_of_inertia_mid,omega_mid,qiw_mid, &
              mean_radius_mid,star%nz)
         if (surface_cz_active) then
            moment_of_inertia_cz = moment_of_inertia_mid(1)
            do i = 2,star%nz
               moment_of_inertia_cz = moment_of_inertia_cz + moment_of_inertia_mid(i)
            end do
         endif
      endif
      return
end subroutine mid_timestep_model
