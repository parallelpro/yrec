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
! 2026 de-tramp (ROADMAP item 3): 29 arguments -> 16. The 13 midpoint
! structure arrays (eta_squared_mid ... am_transport_convective_flag_mid)
! moved into rot_scr (rotation_scratch_lib) -- they are the seculr/
! midmod pipeline's shared per-sub-step workspace, not per-call data.
! The remaining arguments are genuine per-call scalars/bounds. Wrapped
! as a module procedure so the compiler enforces the call.
module mid_timestep_model_lib
      implicit none
contains

subroutine mid_timestep_model(full_timestep, sub_timestep, time_fraction, first_call, &
     moment_of_inertia_cz, cz_mass_bottom, cz_mass_top, &
     core_boundary_zone, envelope_boundary_zone, fully_convective_flag, &
     surface_cz_active, &
     radiative_zone_bounds, convective_zone_bounds, num_radiative_zones, &
     num_convective_zones, ierr)
      use rotation_scratch_lib
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, json
      use net_lib
      use phys_const_lib
      use burn_lib
      use math_lib
      implicit none

      double precision, intent(in) :: full_timestep, sub_timestep
      double precision, intent(inout) :: time_fraction
      logical, intent(in) :: first_call
      double precision, intent(out) :: moment_of_inertia_cz
      double precision, intent(out) :: cz_mass_bottom, cz_mass_top
      integer, intent(out) :: core_boundary_zone, envelope_boundary_zone
      logical, intent(out) :: fully_convective_flag
      logical, intent(out) :: surface_cz_active
      integer, intent(out) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(out) :: num_radiative_zones, num_convective_zones
      integer, intent(out) :: ierr
! --- locals ---
! (deuterium_rate_mid/_start live in rot_scr: carried across sub-steps)
      integer :: num_species_tracked
      integer :: i, j, k
      double precision :: step_fraction_ratio
      double precision :: total_epsilon
      integer :: burn_zone_begin, burn_zone_end
      integer :: solid_zone_start, solid_zone_end

      ierr = 0

      num_species_tracked = 11
      if (star%job%use_extended_composition) num_species_tracked = 15
!  INITIALIZE COMPOSITION ARRAY THE FIRST TIME THROUGH.
!  HCOMPM IS THE ARRAY OF CHANGES IN COMPOSITION DUE TO NUCLEAR BURNING
!  DURING THE COURSE OF THE MODEL TIMESTEP.
! MHP 06/02 SAVE THE CONVECTIVE STATE FROM THE PRIOR STEP.
      if (first_call) then
         do j = 1,star%nz
            rot_scr%convective_flag_prev(j) = star%convective_flag_start(j)
            mix_scr%amum(j) = rot_scr%old_amu(j)
            do i = 1,num_species_tracked
               star%xa(i,j) = star%xa_start(i,j)
            end do
         end do
      else
         do j = 1,star%nz
            rot_scr%convective_flag_prev(j) = rot_scr%convective_flag_mid(j)
         end do
      endif
!  INTERPOLATE IN THE MODEL VARIABLES AND AUXILLARY PHYSICS.
      step_fraction_ratio = sub_timestep/full_timestep
      do j = 1,star%nz
         rot_scr%log_density_mid(j) = star%logRho_start(j) + &
              time_fraction*(star%logRho(j)-star%logRho_start(j))
         rot_scr%log_luminosity_mid(j) = star%luminosity_lsun_start(j) + &
              time_fraction*(star%luminosity_lsun(j)-star%luminosity_lsun_start(j))
         rot_scr%log_pressure_mid(j) = star%logP_start(j) + &
              time_fraction*(star%logP(j)-star%logP_start(j))
         rot_scr%log_radius_mid(j) = star%logR_start(j) + &
              time_fraction*(star%logR(j)-star%logR_start(j))
         rot_scr%log_temperature_mid(j) = star%logT_start(j) + &
              time_fraction*(star%logT(j)-star%logT_start(j))
         rot_scr%hg_mid(j) = star%old_hg(j) + time_fraction*(star%mean_gravity(j) - star%old_hg(j))
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
         if (star%convective_flag(i).eqv.rot_scr%convective_flag_prev(i)) then
            rot_scr%convective_flag_mid(i) = star%convective_flag(i)
            mix_scr%velm(i) = rot_scr%old_vel(i) + time_fraction*(star%conv_vel(i)-rot_scr%old_vel(i))
         else
            if (mix_scr%del_adiabatic_mix(i).lt.mix_scr%del_radiative_mix(i)) then
               rot_scr%convective_flag_mid(i) = .true.
               mix_scr%velm(i) = max(rot_scr%old_vel(i),star%conv_vel(i))
            else
               rot_scr%convective_flag_mid(i) = .false.
               mix_scr%velm(i) = 0.0D0
            endif
         endif
      end do
! (THE ORIGINAL'S MHP 06/02 RE-DISTRIBUTION OF J FOR SHELLS THAT
!  DROPPED OUT OF A RECEDING CONVECTION ZONE WAS HARDWIRED OFF AND IS
!  NOT COMPUTED.)
!  CONVECTIVE OVERSHOOT APPLIED TO NORMAL CONVECTION ZONES.
      call am_convective_regions(star%xa,rot_scr%log_density_mid,rot_scr%log_pressure_mid,rot_scr%log_radius_mid, &
           star%log_mass,rot_scr%log_temperature_mid,rot_scr%convective_flag_mid,star%nz, &
           rot_scr%am_transport_convective_flag_mid,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones, ierr)
      if (ierr /= 0) return
!  REBURN THE ORIGINAL MIXTURE FOR THE SMALL DIFFUSION TIME STEP.
! MHP 6/00 - FOR THERMODYNAMIC CONSISTENCY USE CURRENT T, NOT END OF
! TIMESTEP T (REPLACED HT IN CALL WITH HTM)
      call rotmix(sub_timestep,star%xa,star%dm,rot_scr%log_temperature_mid, &
           star%nz,radiative_zone_bounds,convective_zone_bounds, &
           num_radiative_zones,num_convective_zones &
           ,star%log_total_mass,rot_scr%log_density_mid,star%log_mass,rot_scr%log_radius_mid, &
           rot_scr%log_pressure_mid,rot_scr%am_transport_convective_flag_mid,star%m, ierr)
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
               rot_scr%deuterium_rate_mid(i) = star%deuterium_burning_rate_start(i)+ &
                    step_fraction_ratio*(star%deuterium_burning_rate(i)- &
                    star%deuterium_burning_rate_start(i))
               rot_scr%deuterium_rate_mid_start(i) = star%deuterium_burning_rate_start(i)
            end do
         else
            do i = 1,star%nz
               rot_scr%deuterium_rate_mid_start(i) = rot_scr%deuterium_rate_mid(i)
               rot_scr%deuterium_rate_mid(i) = rot_scr%deuterium_rate_mid(i)+ &
                    step_fraction_ratio*(star%deuterium_burning_rate(i)- &
                    star%deuterium_burning_rate_start(i))
            end do
         endif
! RADIATIVE ZONES.
!
         do k = 1,num_radiative_zones
            do j = radiative_zone_bounds(k,1),radiative_zone_bounds(k,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (rot_scr%log_temperature_mid(j).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
               burn_zone_begin = j
               burn_zone_end = j
               call dburnm(burn_zone_begin,burn_zone_end,star%nz,star%dm, &
                    star%xa,sub_timestep,rot_scr%deuterium_rate_mid, &
                    rot_scr%deuterium_rate_mid_start,step_fraction_ratio)
            end do
         end do
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
         do k = 1,num_convective_zones
            burn_zone_begin = convective_zone_bounds(k,1)
            burn_zone_end = convective_zone_bounds(k,2)
            call dburnm(burn_zone_begin,burn_zone_end,star%nz,star%dm, &
                 star%xa,sub_timestep,rot_scr%deuterium_rate_mid, &
                 rot_scr%deuterium_rate_mid_start,step_fraction_ratio)
         end do
      endif
!  DETERMINE EXTENT OF CENTRAL CONVECTION ZONE.
!  IMIN IS THE FIRST ZONE ABOVE A CENTRAL CONVECTION ZONE, AND THUS THE
!  FIRST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONAL INSTABILITIES.
      if (rot_scr%am_transport_convective_flag_mid(1)) then
         do i = 2,star%nz
            if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
         end do
         core_boundary_zone = max(2,i-1)
      else
         core_boundary_zone = 2
      endif
!  DETERMINE EXTENT OF SURFACE CONVECTION ZONE.
      fully_convective_flag = .false.
      if (rot_scr%am_transport_convective_flag_mid(star%nz)) then
!  SURFACE C.Z. EXISTS.  FIND LOWEST SHELL (IMAX), WHICH IS ALSO THE
!  UPPERMOST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONALLY INDUCED MIXING.
         do i = star%nz-1,1,-1
            if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
         end do
!  (LOOP COMPLETED WITHOUT EXIT, I.E. I = 0: FULLY CONVECTIVE.)
         if (i < 1) fully_convective_flag = .true.
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
! JNT 09/25 FOR 05/15 IMPJMOD=1 IS THE SAME AS LSOLID
      if (.not.star%ctrl%force_solid_body_rotation .and. (star%ctrl%solid_body_mode_flag.ne.1)) then
!  NOW FIND THE RUN OF ROTATION VARIABLES THAT ARE CONSISTENT WITH THE
!  INTERMEDIATE STRUCTURE AND THE RUN OF SPECIFIC ANGULAR MOMENTUM J/M.
!  J/M = I/M * OMEGA.
!  FIRST GUESS AT MOMENT OF INERTIA OF DISTORTED SPERICAL SHELLS:
!  I/M = 2/3 R**2.
         do i = 1,star%nz
            rot_scr%moment_of_inertia_mid(i) = cc23*exp(ln10*2.0D0*rot_scr%log_radius_mid(i))
            rot_scr%omega_mid(i) = star%j_rot(i)/rot_scr%moment_of_inertia_mid(i)
            rot_scr%moment_of_inertia_mid(i) = star%dm(i)*rot_scr%moment_of_inertia_mid(i)
         end do
!  SOLVE FOR THE ANGULAR VELOCITIES OF THE SHELLS GIVEN THE SPECIFIC
!  ANGULAR MOMENTUM (HJM) AND A FIRST GUESS AT THE ANGULAR VELOCITY(OMEGAM)
!  AND MOMENT OF INERTIA (HIM).
         call omega_from_j(rot_scr%log_density_mid,star%j_rot,rot_scr%log_radius_mid, &
              star%log_mass,star%dm,rot_scr%am_transport_convective_flag_mid,star%nz, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,rot_scr%omega_mid,rot_scr%qiw_mid, &
              rot_scr%mean_radius_mid)
!  FIND TOTAL MOMENT OF INERTIA OF THE SURFACE C.Z. IF APPLICABLE.
         if (surface_cz_active) then
            moment_of_inertia_cz = rot_scr%moment_of_inertia_mid(star%nz)
            do i = star%nz-1,1,-1
               if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
               moment_of_inertia_cz = moment_of_inertia_cz + rot_scr%moment_of_inertia_mid(i)
            end do
         endif
      else
         solid_zone_start = 1
         solid_zone_end = star%nz
         call solid_body_omega(rot_scr%log_density_mid,star%j_rot,rot_scr%log_radius_mid, &
              star%log_mass,star%dm,solid_zone_start,solid_zone_end, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,rot_scr%omega_mid,rot_scr%qiw_mid, &
              rot_scr%mean_radius_mid,star%nz)
         if (surface_cz_active) then
            moment_of_inertia_cz = rot_scr%moment_of_inertia_mid(1)
            do i = 2,star%nz
               moment_of_inertia_cz = moment_of_inertia_cz + rot_scr%moment_of_inertia_mid(i)
            end do
         endif
      endif
      return
end subroutine mid_timestep_model

end module mid_timestep_model_lib
