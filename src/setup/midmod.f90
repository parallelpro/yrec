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
subroutine midmod(full_timestep,sub_timestep,time_fraction,composition, &
     log_density,hg,specific_angular_momentum,log_luminosity,log_pressure, &
     log_radius,log_mass,enclosed_mass, &
! HS2,HSTOT,HT,LC,LCZ,LFIRST,M,ETA2M,HDM,HGM,HICZ,HIM,  ! KC 2025-05-31
     shell_mass,log_total_mass,log_temperature,convective_flag,first_call, &
     num_zones,eta_squared_mid,log_density_mid,hg_mid,moment_of_inertia_cz, &
     moment_of_inertia_mid, &
     log_luminosity_mid,log_pressure_mid,log_radius_mid,cz_mass_bottom, &
     cz_mass_top,log_temperature_mid,core_boundary_zone, &
     envelope_boundary_zone,fully_convective_flag,convective_flag_mid, &
     am_transport_convective_flag_mid,surface_cz_active,omega_mid, &
     mean_radius_mid,qiw_mid,radiative_zone_bounds,convective_zone_bounds, &
     num_radiative_zones,num_convective_zones, ierr)

      use nuclear_lib
      use rotdiff_lib
      use run_diag_lib
      use temp_lib
      use mdphy_lib
      use light_burn_lib
      use turnover_lib
      use scrtch_lib
      use oldmod_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: full_timestep, sub_timestep
      double precision, intent(inout) :: time_fraction
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: log_density(json), hg(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      double precision, intent(in) :: log_luminosity(json), &
           log_pressure(json), log_radius(json), log_mass(json)
      double precision, intent(in) :: enclosed_mass(json), shell_mass(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: log_temperature(json)
! LC(JSON),LCZ(JSON),ETA2M(JSON),HDM(JSON),HGM(JSON),  ! KC 2025-05-31
      logical, intent(in) :: convective_flag(json)
      logical, intent(in) :: first_call
      integer, intent(in) :: num_zones
      double precision, intent(out) :: eta_squared_mid(json), &
           log_density_mid(json), hg_mid(json)
      double precision, intent(out) :: moment_of_inertia_cz
      double precision, intent(out) :: moment_of_inertia_mid(json)
      double precision, intent(out) :: log_luminosity_mid(json), &
           log_pressure_mid(json)
! log_radius_mid (HRM) is read here before being (re)computed in this
! call, in the LFIRST=.false. branch below that seeds rot_diff%radius_prev from
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
      save

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
      if (use_extended_composition) num_species_tracked = 15
!  INITIALIZE COMPOSITION ARRAY THE FIRST TIME THROUGH.
!  HCOMPM IS THE ARRAY OF CHANGES IN COMPOSITION DUE TO NUCLEAR BURNING
!  DURING THE COURSE OF THE MODEL TIMESTEP.
! MHP 06/02 SAVE CHANGES IN R, CZ DEPTH, DELAD-DELRAD FROM THE
! PRIOR STEP.
      if (first_call) then
         do 20 j = 1,num_zones
            rot_diff%convective_flag_prev(j) = prev_model%old_convective_flag(j)
            rot_diff%radius_prev(j) = prev_model%old_radius(j)
            rot_diff%del_grad_diff_prev(j) = rot_diff%old_del_adiabatic_mix(j)-rot_diff%old_del_radiative_mix(j)
            mix_phys%amum(j) = rot_diff%old_amu(j)
            do 10 i = 1,num_species_tracked
               composition(i,j) = prev_model%old_composition(i,j)
   10       continue
   20    continue
      else
         do j = 1,num_zones
            rot_diff%radius_prev(j) = log_radius_mid(j)
            rot_diff%del_grad_diff_prev(j) = rot_diff%del_grad_diff_new(j)
            rot_diff%convective_flag_prev(j) = convective_flag_mid(j)
         end do
      endif
      new_cz_detected = .false.
!  INTERPOLATE IN THE MODEL VARIABLES AND AUXILLARY PHYSICS.
      step_fraction_ratio = sub_timestep/full_timestep
      do 40 j = 1,num_zones
         log_density_mid(j) = prev_model%old_density(j) + &
              time_fraction*(log_density(j)-prev_model%old_density(j))
         log_luminosity_mid(j) = prev_model%old_luminosity(j) + &
              time_fraction*(log_luminosity(j)-prev_model%old_luminosity(j))
         log_pressure_mid(j) = prev_model%old_pressure(j) + &
              time_fraction*(log_pressure(j)-prev_model%old_pressure(j))
         log_radius_mid(j) = prev_model%old_radius(j) + &
              time_fraction*(log_radius(j)-prev_model%old_radius(j))
         log_temperature_mid(j) = prev_model%old_temperature(j) + &
              time_fraction*(log_temperature(j)-prev_model%old_temperature(j))
!        DO 30 I = 1,IEND
!           HCOMP(I,J)=HCOMP(I,J)+FAC2*HCOMPM(I,J)
!           HCOMPP(I,J) = HCOMP(I,J)
!  30    CONTINUE
         hg_mid(j) = run_diag%old_hg(j) + time_fraction*(hg(j) - run_diag%old_hg(j))
         mix_phys%del_adiabatic_mix(j) = rot_diff%old_del_adiabatic_mix(j) + &
              time_fraction*(shell_diag%del_grad(3,j)-rot_diff%old_del_adiabatic_mix(j))
         mix_phys%delm(j) = rot_diff%old_delm(j) + time_fraction*(shell_diag%del_grad(2,j) - rot_diff%old_delm(j))
         mix_phys%del_radiative_mix(j) = rot_diff%old_del_radiative_mix(j) + &
              time_fraction*(shell_diag%del_grad(1,j) - rot_diff%old_del_radiative_mix(j))
         mix_phys%esumm(j) = rot_diff%old_esum(j) + time_fraction*(shell_diag%sesum(j) - rot_diff%old_esum(j))
         mix_phys%viscm(j) = rot_diff%old_visc(j) + time_fraction*(shell_temp%visc(j) - rot_diff%old_visc(j))
         mix_phys%thdifm(j) = rot_diff%old_thdif(j) + time_fraction*(shell_temp%thdif(j) - rot_diff%old_thdif(j))
         mix_phys%cpm(j) = rot_diff%old_cp(j) + time_fraction*(shell_temp%cp(j) - rot_diff%old_cp(j))
         mix_phys%qdtm(j) = rot_diff%old_qdt(j) + time_fraction*(shell_temp%qdt(j) - rot_diff%old_qdt(j))
         mix_phys%om(j) = rot_diff%old_om(j) + time_fraction*(shell_diag%so(j) - rot_diff%old_om(j))
         mix_phys%amum(j) = mix_phys%amum(j) + step_fraction_ratio*(shell_temp%mean_molecular_weight(j) - rot_diff%old_amu(j))
! MHP 6/00 ADDED TOTAL ENERGY GENERATION
         total_epsilon = shell_diag%sesum(j)+shell_diag%seg(6,j)+shell_diag%seg(7,j)
         mix_phys%epsm(j) = rot_diff%old_eps(j)+time_fraction*(total_epsilon-rot_diff%old_eps(j))
   40 continue
!  CHECK FOR ADVANCING OR RECEDING CONVECTIVE REGIONS.USE INTERPOLATED
!  RADIATIVE AND ADIABATIC TEMPERATURE GRADIENTS TO DETERMINE WHETHER
!  OR NOT A ZONE IS CONVECTIVE IF IT CHANGES STATE OVER THE COURSE OF A
!  TIMESTEP.
      do 50 i = 1,num_zones
         rot_diff%del_grad_diff_new(i) = mix_phys%del_adiabatic_mix(i)-mix_phys%del_radiative_mix(i)
         if (convective_flag(i).eqv.rot_diff%convective_flag_prev(i)) then
            convective_flag_mid(i) = convective_flag(i)
            mix_phys%velm(i) = rot_diff%old_vel(i) + time_fraction*(shell_diag%svel(i)-rot_diff%old_vel(i))
            convective_state_changed(i) = .false.
         else
            convective_state_changed(i) = .true.
            new_cz_detected = .true.
            if (mix_phys%del_adiabatic_mix(i).lt.mix_phys%del_radiative_mix(i)) then
               convective_flag_mid(i) = .true.
               mix_phys%velm(i) = max(rot_diff%old_vel(i),shell_diag%svel(i))
            else
               convective_flag_mid(i) = .false.
               mix_phys%velm(i) = 0.0D0
            endif
         endif
   50 continue
      new_cz_detected = .false.
! MHP 06/02 IF THE CZ DEPTHS HAVE CHANGED, RESOLVE
! OUT WHEN A GIVEN ZONE "DROPPED OUT" OF THE CONVECTION
! ZONE OR ENTERED IT.  MODIFY THE SPECIFIC ANGULAR MOMENTUM
! ACCORDINGLY.
      if (new_cz_detected) then
         cz_change_active = .false.
         redistribute_j_flag = .false.
         do i = 1,num_zones-1
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
                    rot_diff%del_grad_diff_new(change_region_start), &
                    rot_diff%del_grad_diff_new(change_region_start-1), &
                    rot_diff%del_grad_diff_new(change_region_end+1)
! ONLY CHANGE IF THERE IS A RECEDING CONVECTION ZONE;
! WE CAN ASSUME LOCAL CONSERVATION OF J PRIOR TO HAVING
! A SHELL INCORPORATED INTO A CZ.
               if (rot_diff%del_grad_diff_new(change_region_start).gt.0.0D0) then
                  if (rot_diff%del_grad_diff_new(change_region_start-1).lt.0.0D0) then
! SUBTRACTED FROM A CORE CZ
                     cz_direction_flag = -1
                  else if (rot_diff%del_grad_diff_new(change_region_end+1).lt.0.0D0) then
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
                        goto 11
                     endif
                  end do
                  cz_zone_bottom = 1
 11               continue
                  do j = change_region_end,change_region_start,-1
                     if (rot_diff%del_grad_diff_new(j).lt.0.0D0) goto 13
                     cz_moment_of_inertia = 0.0D0
                     cz_angular_momentum = 0.0D0
! FRACTION OF THE TIMESTEP SHELL WAS RADIATIVE
                     time_fraction = rot_diff%del_grad_diff_new(j)/ &
                          (rot_diff%del_grad_diff_new(j)-rot_diff%del_grad_diff_prev(j))
! FRACTION OF THE TIMESTEP THAT THE SHELL WAS CONVECTIVE
                     convective_fraction = 1.0D0 - time_fraction
                     do ii = cz_zone_bottom,j-1
                        radius_interp = rot_diff%radius_prev(ii)+ &
                             convective_fraction*(log_radius(ii)-prev_model%old_radius(ii))
                        cz_moment_of_inertia = cz_moment_of_inertia+ &
                             cc23*shell_mass(ii)*10.0D0**(2.0D0*radius_interp)
                        cz_angular_momentum = cz_angular_momentum + &
                             specific_angular_momentum(ii)*shell_mass(ii)
                     end do
! GET THE CZ MOMENT OF INERTIA AND TOTAL J AT THE TIME THE
! SHELL WAS RELEASED; THIS GIVES A CORRECTED VALUE FOR
! J/M.  REDISTRIBUTE THE DIFFERENCE IN J BACK INTO THE CZ.
                     radius_shell_factor = 2.0D0*(rot_diff%radius_prev(j)+ &
                          convective_fraction*(log_radius(j)-rot_diff%radius_prev(j)))
                     specific_angular_momentum_corrected = &
                          cc23*(cz_angular_momentum/cz_moment_of_inertia)* &
                          10.0D0**radius_shell_factor
                     delta_angular_momentum = &
                          (specific_angular_momentum(j)- &
                          specific_angular_momentum_corrected)*shell_mass(j)
                     angular_momentum_ratio = &
                          (cz_angular_momentum-delta_angular_momentum)/ &
                          cz_angular_momentum
                     specific_angular_momentum(j)=specific_angular_momentum_corrected
                     angular_momentum_check_ratio = &
                          specific_angular_momentum_corrected/ &
                          specific_angular_momentum(j)
                     write(*,*)j,time_fraction,angular_momentum_check_ratio, &
                          delta_angular_momentum,angular_momentum_ratio
                     do ii = cz_zone_bottom,j-1
                        specific_angular_momentum(ii) = &
                             specific_angular_momentum(ii)*angular_momentum_ratio
                     end do
                  end do
               else
! ENVELOPE CZ
! LOCATE EDGES OF CURRENT LOWER CZ
                  cz_zone_bottom = change_region_end + 1
                  do ii = change_region_end+2,num_zones
                     if (.not.convective_flag_mid(ii)) then
                        cz_zone_top = ii - 1
                        goto 12
                     endif
                  end do
                  cz_zone_top = num_zones
 12               continue
                  do j = change_region_start,change_region_end
                     if (rot_diff%del_grad_diff_new(j).lt.0.0D0) goto 13
                     cz_moment_of_inertia = 0.0D0
                     cz_angular_momentum = 0.0D0
! FRACTION OF THE TIMESTEP SHELL WAS RADIATIVE
                     time_fraction = rot_diff%del_grad_diff_new(j)/ &
                          (rot_diff%del_grad_diff_new(j)-rot_diff%del_grad_diff_prev(j))
! FRACTION OF THE TIMESTEP THAT THE SHELL WAS CONVECTIVE
                     convective_fraction = 1.0D0 - time_fraction
                     do ii = j+1,cz_zone_top
                        radius_interp = rot_diff%radius_prev(ii)+ &
                             convective_fraction*(log_radius(ii)-prev_model%old_radius(ii))
                        cz_moment_of_inertia = cz_moment_of_inertia+ &
                             cc23*shell_mass(ii)*10.0D0**(2.0D0*radius_interp)
                        cz_angular_momentum = cz_angular_momentum + &
                             specific_angular_momentum(ii)*shell_mass(ii)
                     end do
! GET THE CZ MOMENT OF INERTIA AND TOTAL J AT THE TIME THE
! SHELL WAS RELEASED; THIS GIVES A CORRECTED VALUE FOR
! J/M.  REDISTRIBUTE THE DIFFERENCE IN J BACK INTO THE CZ.
                     radius_shell_factor = 2.0D0*(rot_diff%radius_prev(j)+ &
                          convective_fraction*(log_radius(j)-rot_diff%radius_prev(j)))
                     specific_angular_momentum_corrected = &
                          cc23*(cz_angular_momentum/cz_moment_of_inertia)* &
                          10.0D0**radius_shell_factor
                     delta_angular_momentum = &
                          (specific_angular_momentum(j)- &
                          specific_angular_momentum_corrected)*shell_mass(j)
                     angular_momentum_ratio = &
                          (cz_angular_momentum-delta_angular_momentum)/ &
                          cz_angular_momentum
                     specific_angular_momentum(j)=specific_angular_momentum_corrected
                     angular_momentum_check_ratio = &
                          specific_angular_momentum_corrected/ &
                          specific_angular_momentum(j)
                     write(*,*)j,time_fraction,angular_momentum_check_ratio, &
                          delta_angular_momentum,angular_momentum_ratio
                     do ii = j+1,num_zones
                        specific_angular_momentum(ii) = &
                             specific_angular_momentum(ii)*angular_momentum_ratio
                     end do
                  end do
               endif
 13            continue
            endif
         end do
      endif
!  CONVECTIVE OVERSHOOT APPLIED TO NORMAL CONVECTION ZONES.
      call ovrot(composition,log_density_mid,log_pressure_mid,log_radius_mid, &
           log_mass,log_temperature_mid,convective_flag_mid,num_zones, &
           am_transport_convective_flag_mid,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones)
!  REBURN THE ORIGINAL MIXTURE FOR THE SMALL DIFFUSION TIME STEP.
! MHP 6/00 - FOR THERMODYNAMIC CONSISTENCY USE CURRENT T, NOT END OF
! TIMESTEP T (REPLACED HT IN CALL WITH HTM)
!       CALL ROTMIX(DT,HCOMP,HS2,HTM,5,M,MRZONE,MXZONE,  ! KC 2025-05-31
      call rotmix(sub_timestep,composition,shell_mass,log_temperature_mid, &
           num_zones,radiative_zone_bounds,convective_zone_bounds, &
           num_radiative_zones,num_convective_zones &
           ,log_total_mass,log_density_mid,log_mass,log_radius_mid, &
           log_pressure_mid,am_transport_convective_flag_mid,enclosed_mass, ierr)
      if (ierr /= 0) return
      do i = 1,num_zones
         do j = 1,num_species_tracked
            rot_diff%composition_snapshot(j,i) = composition(j,i)
         end do
      end do
! MHP 05/02 ADDED DEUTERIUM BURNING
      if (prev_model%old_composition(12,num_zones).gt.1.0D-14) then
! INCREMENT THE TIMESTEP
         if (first_call) then
            do i = 1,num_zones
               deuterium_rate_mid(i) = light_burn%deuterium_burning_rate_start(i)+ &
                    step_fraction_ratio*(light_burn%deuterium_burning_rate(i)- &
                    light_burn%deuterium_burning_rate_start(i))
               deuterium_rate_mid_start(i) = light_burn%deuterium_burning_rate_start(i)
            end do
         else
            do i = 1,num_zones
               deuterium_rate_mid_start(i) = deuterium_rate_mid(i)
               deuterium_rate_mid(i) = deuterium_rate_mid(i)+ &
                    step_fraction_ratio*(light_burn%deuterium_burning_rate(i)- &
                    light_burn%deuterium_burning_rate_start(i))
            end do
         endif
! RADIATIVE ZONES.
!
         do k = 1,num_radiative_zones
            do j = radiative_zone_bounds(k,1),radiative_zone_bounds(k,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (log_temperature_mid(j).le.tcut(1)) goto 190
               burn_zone_begin = j
               burn_zone_end = j
               call dburnm(burn_zone_begin,burn_zone_end,num_zones,shell_mass, &
                    composition,sub_timestep,deuterium_rate_mid, &
                    deuterium_rate_mid_start,step_fraction_ratio)
           end do
        end do
 190    continue
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
         do k = 1,num_convective_zones
            burn_zone_begin = convective_zone_bounds(k,1)
            burn_zone_end = convective_zone_bounds(k,2)
            call dburnm(burn_zone_begin,burn_zone_end,num_zones,shell_mass, &
                 composition,sub_timestep,deuterium_rate_mid, &
                 deuterium_rate_mid_start,step_fraction_ratio)
         end do
      endif
!  DETERMINE EXTENT OF CENTRAL CONVECTION ZONE.
!  IMIN IS THE FIRST ZONE ABOVE A CENTRAL CONVECTION ZONE, AND THUS THE
!  FIRST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONAL INSTABILITIES.
      if (am_transport_convective_flag_mid(1)) then
         do 60 i = 2,num_zones
            if (.not.am_transport_convective_flag_mid(i)) goto 65
   60    continue
         i = num_zones + 1
   65    core_boundary_zone = max(2,i-1)
      else
         core_boundary_zone = 2
      endif
!  DETERMINE EXTENT OF SURFACE CONVECTION ZONE.
      fully_convective_flag = .false.
      if (am_transport_convective_flag_mid(num_zones)) then
!  SURFACE C.Z. EXISTS.  FIND LOWEST SHELL (IMAX), WHICH IS ALSO THE
!  UPPERMOST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONALLY INDUCED MIXING.
         do 70 i = num_zones-1,1,-1
            if (.not.am_transport_convective_flag_mid(i)) goto 80
   70    continue
         fully_convective_flag = .true.
         i = 0
   80    envelope_boundary_zone = i + 1
!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
         cz_mass_top = exp(ln10*log_total_mass)
         if (envelope_boundary_zone.gt.1) then
            cz_mass_bottom = 0.5D0*(enclosed_mass(envelope_boundary_zone)+ &
                 enclosed_mass(envelope_boundary_zone-1))
         else
            cz_mass_bottom = 0.0D0
         endif
!  LCZSUR=T IF A SURFACE C.Z.DEEP ENOUGH FOR ANGULAR MOMENTUM LOSS EXISTS
         if ((cz_mass_top-cz_mass_bottom)/solar_mass_cgs.gt.0.0D0) then
            surface_cz_active = .true.
         else
            surface_cz_active= .false.
         endif
      else
!  NO SURFACE C.Z.
         envelope_boundary_zone = num_zones
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
!            DD2 = DELRM(IMAX-1)-DELAM(IMAX-1)
!            DD1 = DELRM(IMAX)-DELAM(IMAX)
!            FX = DD2/(DD2-DD1)
!C            WRITE(*,911)IMAX,FX,DD2,DD1
!C 911        FORMAT(I5,1P3E12.3)
!C INFER HP
!            ENVRL = HRM(IMAX-1)+FX*(HRM(IMAX)-HRM(IMAX-1))
!            ENVR = EXP(CLN*ENVRL)
!            PS2 = EXP(CLN*(HPM(IMAX)-HDM(IMAX)))/HGM(IMAX)
!            PS1 = EXP(CLN*(HPM(IMAX-1)-HDM(IMAX-1)))/HGM(IMAX-1)
!            PSCA = PS1 + FX*(PS2-PS1)
!            RTESTL = DLOG10(ENVR+PSCA)
!C FIND V
!            DO K = IMAX+1,M
!               IF(HRM(K).GT.RTESTL)THEN
!                  FX = (RTESTL-HRM(K-1))/(HRM(K)-HRM(K-1))
!                  CVEL = VELM(K-1)+FX*(VELM(K)-VELM(K-1))
!                  GOTO 85
!               ENDIF
!            END DO
!            CVEL = VELM(M)
!C 85         CONTINUE
!C DEFINE TAUCZ
!            TAUCZ = PSCA/CVEL
!C            WRITE(*,911)K,PSCA,CVEL,TAUCZ
!         ELSE
!C INFER HP
!            PSCA2 = EXP(CLN*(HPM(1)-HDM(1)))/HGM(1)
!            RTEST2 = EXP(CLN*HR(1))
!            IF(PSCA2.LE.RTEST2)THEN
!C HP < R AT THE FIRST POINT.  ASSUME V CONSTANT INSIDE AND HP = K/R FOR
!C SLOWLY VARYING DENSITY AND PRESSURE NEAR THE CENTER.
!               CVEL = VELM(1)
!               PSCA = (PSCA2*RTEST2)**0.5D0
!               TAUCZ = PSCA/CVEL
!C               WRITE(*,912)PSCA2,RTEST2,PSCA,CVEL,TAUCZ
!C 912           FORMAT(1P5E12.3)
!            ELSE
!               DO K = 2,M
!                  PSCA1 = PSCA2
!                  RTEST1 = RTEST2
!                  PSCA2 = EXP(CLN*(HPM(K)-HDM(K)))/HGM(K)
!                  RTEST2 = EXP(CLN*HR(K))
!C FIND LOCATION WHERE HP = R
!                  IF(PSCA2.LE.RTEST2)THEN
!                     FX = (RTEST1-PSCA1)/((PSCA2-RTEST2)-(PSCA1-RTEST1))
!C FIND V
!                     CVEL = VELM(K-1)+FX*(VELM(K)-VELM(K-1))
!                     PSCA = PSCA1+FX*(PSCA2-PSCA1)
!C DEFINE TAUCZ
!                     TAUCZ = PSCA/CVEL
!                     GOTO 95
!                  ENDIF
!               END DO
!               K = M
!               CVEL = VELM(M)
!               PSCA = PSCA2
!               TAUCZ = PSCA/CVEL
!C 95            CONTINUE
!C               WRITE(*,911)K,PSCA,CVEL,TAUCZ
!            ENDIF
!         ENDIF
!      ENDIF
! JNT 09/25 FOR 05/15 IMPJMOD=1 IS THE SAME AS LSOLID
      if (.not.force_solid_body_rotation .and. (solid_body_mode_flag.ne.1)) then
!  NOW FIND THE RUN OF ROTATION VARIABLES THAT ARE CONSISTENT WITH THE
!  INTERMEDIATE STRUCTURE AND THE RUN OF SPECIFIC ANGULAR MOMENTUM J/M.
!  J/M = I/M * OMEGA.
!  FIRST GUESS AT MOMENT OF INERTIA OF DISTORTED SPERICAL SHELLS:
!  I/M = 2/3 R**2.
         do i = 1,num_zones
            moment_of_inertia_mid(i) = cc23*exp(ln10*2.0D0*log_radius_mid(i))
            omega_mid(i) = specific_angular_momentum(i)/moment_of_inertia_mid(i)
            moment_of_inertia_mid(i) = shell_mass(i)*moment_of_inertia_mid(i)
         end do
!  SOLVE FOR THE ANGULAR VELOCITIES OF THE SHELLS GIVEN THE SPECIFIC
!  ANGULAR MOMENTUM (HJM) AND A FIRST GUESS AT THE ANGULAR VELOCITY(OMEGAM)
!  AND MOMENT OF INERTIA (HIM).
         call getrot(log_density_mid,specific_angular_momentum,log_radius_mid, &
              log_mass,shell_mass,am_transport_convective_flag_mid,num_zones, &
              eta_squared_mid,moment_of_inertia_mid,omega_mid,qiw_mid, &
              mean_radius_mid)
!  FIND TOTAL MOMENT OF INERTIA OF THE SURFACE C.Z. IF APPLICABLE.
         if (surface_cz_active) then
            moment_of_inertia_cz = moment_of_inertia_mid(num_zones)
            do i = num_zones-1,1,-1
               if (.not.am_transport_convective_flag_mid(i)) goto 110
               moment_of_inertia_cz = moment_of_inertia_cz + moment_of_inertia_mid(i)
            end do
  110       continue
         endif
      else
         solid_zone_start = 1
         solid_zone_end = num_zones
         call solid(log_density_mid,specific_angular_momentum,log_radius_mid, &
              log_mass,shell_mass,solid_zone_start,solid_zone_end, &
              eta_squared_mid,moment_of_inertia_mid,omega_mid,qiw_mid, &
              mean_radius_mid,num_zones)
         if (surface_cz_active) then
            moment_of_inertia_cz = moment_of_inertia_mid(1)
            do i = 2,num_zones
               moment_of_inertia_cz = moment_of_inertia_cz + moment_of_inertia_mid(i)
            end do
         endif
      endif
      return
end subroutine midmod
