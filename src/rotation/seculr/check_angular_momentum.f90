!----------------------------------------------------------------------
! checkj
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original checkj.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 10/02 ECOD, ECOD2 NO LONGER USED; OMITTED FROM CALL
!
! SR CHECKJ PERFORMS SEVERAL FUNCTIONS.
! FIRST, IT CHECKS FOR NEGATIVE SPECIFIC ANGULAR MOMENTA.
! IF THEY ARE ENCOUNTERED, THE TIMESTEP IS CUT.
! SECOND, IT GUARDS AGAINST REVERSAL OF ANGULAR VELOCITY GRADIENTS.
! THIRD, IT COMPARES THE RUN OF OMEGA TO THAT FOR THE PREVIOUS ITERATION.
! IF THE RUN HAS CHANGED BY LITTLE ENOUGH THEN THE DIFFUSION EQUATIONS
! ARE CONSIDERED CONVERGED.  LOK=T IF THE RUN HAS CONVERGED.
! FOURTH, IT CORRECTS THE MOMENT OF INERTIA FOR CHANGES IN OMEGA.
! INPUT VARIABLES:
! FIFTH, IT WRITES OUT DETAILS OF THE DIFFUSION IF REQUESTED.
!
! log_density : RUN OF LOG DENSITY.
! specific_angular_momentum : RUN OF SPECIFIC ANGULAR MOMENTUM AFTER
!    THE CURRENT ITERATION.
! specific_angular_momentum_prev : RUN OF SPECIFIC ANGULAR MOMENTUM
!    AFTER THE LAST ITERATION FOR THE DIFFUSION COEFFICIENTS.
! specific_angular_momentum_start : RUN OF SPECIFIC ANGULAR MOMENTUM AT
!    THE BEGINNING OF THE TIMESTEP.
! log_radius : RUN OF LOG RADIUS.
! log_mass : RUN OF LOG MASS.
! shell_mass : MASS (UNLOGGED) CONTAINED IN EACH OF THE MODEL POINTS.
! diffusion_velocity : RUN OF CHARACTERISTIC DIFFUSION VELOCITIES. IF
!    diffusion_velocity(I)>0 THEN ZONE I IS UNSTABLE WITH RESPECT TO
!    ZONE I-1.
! zone_min,zone_max : THE FIRST AND LAST ZONES AT WHICH
!    diffusion_velocity IS COMPUTED.
!   *NOTE: zone_min = 2 AND zone_max = NUMBER OF MODEL POINTS UNLESS A
!    SURFACE OR CENTRAL CONVECTION ZONE EXISTS.
! iteration_number : ITERATION NUMBER.
! itdif2 : MAXIMUM NUMBER OF ITERATIONS ALLOWED IN A GIVEN
!    DIFFUSION TIMESTEP.
! am_transport_convective_flag : ARRAY SET T IF A ZONE IS CONVECTIVE
!    FOR ANGULAR MOMENTUM REDISTRIBUTION PURPOSES (I.E. INCLUDES
!    OVERSHOOT REGION.)
! LPRT : FLAG SET T IF MODEL I/O DESIRED.
! num_zones : NUMBER OF MODEL POINTS.
!
! OUTPUT VARIABLES:
!
! dt : DIFFUSION TIMESTEP, WHICH CAN BE CUT IF ERRORS IN THE DIFFUSION
!    ARE DETECTED.
! eta_squared,mean_radius : AUXILLARY QUANTITIES USED TO FIND OMEGA
!    GIVEN J/M.
! moment_of_inertia,qiw : RUN OF MOMENTS OF INERTIA AND THEIR
!    DERIVATIVES WITH RESPECT TO OMEGA.
! specific_angular_momentum_prev : ARRAY IS SET EQUAL TO
!    specific_angular_momentum AFTER CONVERGENCE IS CHECKED.
! cut_count : NUMBER OF TIMES DIFFUSION TIMESTEP HAS BEEN CUT.
! converged_flag : FLAG SET T IF DIFFUSION CEOFFICIENTS HAVE CONVERGED.
! redo_flag : FLAG SET T IF A PROBLEM REQUIRES CUTTING THE DIFFUSION
!    TIMESTEP.
! omega : RUN OF ANGULAR VELOCITY.
subroutine check_angular_momentum(log_density, specific_angular_momentum_prev, &
     specific_angular_momentum_start, log_radius, log_mass, shell_mass, &
     diffusion_velocity, zone_min, zone_max, iteration_number, &
! LPRT,M,DT,ETA2,HI,HJM,IREDO,LOK,LREDO,OMEGA,  ! KC 2025-05-31
     am_transport_convective_flag, num_zones, dt, eta_squared, &
     moment_of_inertia, specific_angular_momentum, cut_count, &
     converged_flag, redo_flag, omega, &
! QIW,R0,WSAV,ID,IDM,ECOD,ECOD2,LOKAD)
     qiw, mean_radius, omega_start, print_zone_id, print_zone_count, &
     already_converged_flag, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use luout_lib
      implicit none

      double precision, intent(inout) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum_prev(json)
      double precision, intent(in) :: specific_angular_momentum_start(json)
      double precision, intent(inout) :: log_radius(json), log_mass(json), &
           shell_mass(json)
      double precision, intent(in) :: diffusion_velocity(json)
      integer, intent(in) :: zone_min, zone_max, iteration_number
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: dt
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: converged_flag
      logical, intent(out) :: redo_flag
      double precision, intent(inout) :: omega(json)
      double precision, intent(inout) :: qiw(json), mean_radius(json)
      double precision, intent(in) :: omega_start(json)
      integer, intent(inout) :: print_zone_id(json)
      integer, intent(inout) :: print_zone_count
!     ECOD(JSON),ECOD2(JSON)
      logical, intent(in) :: already_converged_flag
! locals
      double precision :: max_delta_j_by_iter(16)
      integer :: max_delta_j_zone_by_iter(16)
      integer :: zone_index, scan_index
      double precision :: delta_j_fraction
      double precision :: saved_tolerance, saved_acc_tolerance
      integer :: zone_bottom, zone_top
      integer :: loop_start
      double precision :: max_fractional_dj
      integer :: max_dj_zone
      integer :: print_zone_begin, print_zone_end

!  CHECK FOR NEGATIVE J/M.
      integer, intent(out) :: ierr

      ierr = 0

      converged_flag = .false.
      redo_flag = .false.
      if(already_converged_flag)then
         converged_flag = .true.
         redo_flag = .false.
      endif
      do zone_index = 1,num_zones
         if(specific_angular_momentum(zone_index).le.0.0d0) then
            cut_count = cut_count + 1
!  STOP IF TIMESTEP CUT MORE THAN 3 TIMES.
!            IF(IREDO.GT.3)THEN
            if(cut_count.gt.0)then
               write(6,1000) zone_index
               write(run_log_unit,1000) zone_index
 1000          format(1x,39('>'),39('<')/5x,'ERROR IN SR CHECKJ'/ &
                       5x,'NEGATIVE J/M ENCOUNTERED IN ZONE',i5, &
                       ' AND 3 ATTEMPTS AT CUTTING TIMESTEP FAILED'/ &
                       'RUN STOPPED')
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
               ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
               ! preserve the historical stop on a nonzero return.
               ierr = 1
               return
            else
               redo_flag = .true.
               dt = 0.5d0*dt
               write(6,1005)cut_count,zone_index
               write(run_log_unit,1005)cut_count,zone_index
 1005          format(5x,'ERROR IN SR CHECKJ'/5x,'TIMESTEP CUT,',1x, &
                       'NUMBER',i5,' DUE TO NEGATIVE J/M IN ZONE',i5)
               continue
               return
            endif
         endif
      end do
!  CHECK IF THE FRACTIONAL CHANGE IN OMEGA RELATIVE TO THE PREVIOUS
!  ITERATION IS SMALL ENOUGH TO BE CONSIDERED CONVERGED.
!  ALSO LOCATE THE ZONE WHERE THE MAXIMUM CHANGE OCCURS FOR I/O.
      max_delta_j_by_iter(iteration_number) = &
           (specific_angular_momentum(1)-specific_angular_momentum_prev(1))/ &
           specific_angular_momentum_prev(1)
      specific_angular_momentum_prev(1) = specific_angular_momentum(1)
      max_delta_j_zone_by_iter(iteration_number) = 1
      do zone_index = 2,num_zones
         delta_j_fraction = &
              (specific_angular_momentum(zone_index)- &
              specific_angular_momentum_prev(zone_index))/ &
              specific_angular_momentum_prev(zone_index)
         if(abs(delta_j_fraction).gt. &
              abs(max_delta_j_by_iter(iteration_number))) then
            max_delta_j_by_iter(iteration_number) = delta_j_fraction
            max_delta_j_zone_by_iter(iteration_number) = zone_index
         endif
         specific_angular_momentum_prev(zone_index) = &
              specific_angular_momentum(zone_index)
      end do
      if(abs(max_delta_j_by_iter(iteration_number)).le. &
           star%ctrl%convergence_tolerance) then
!         LOK = .FALSE.
!      ELSE
         converged_flag = .true.
      endif
!  FIND THE RUN OF OMEGA THAT CORRESPONDS TO THE NEW RUN OF J/M.
!  THE MOMENT OF INERTIA IS A FUNCTION OF OMEGA WHICH IS SOLVED
!  ITERATIVELY.  BECAUSE THE ERROR IN THE DIFFUSION IS OF ORDER DJOK,
!  RELAX TOLERANCE FOR MOMENT OF INERTIA ITERATION EXCEPT FOR THE
!  FINAL STEP.
      saved_tolerance = rot_scr%moment_of_inertia_tolerance
      saved_acc_tolerance = star%job%acfpft
      if(iteration_number.lt.star%ctrl%itdif2.and..not.converged_flag)then
         rot_scr%moment_of_inertia_tolerance = &
              max(star%ctrl%convergence_tolerance*1.0d-2,saved_tolerance)
         star%job%acfpft = max(star%ctrl%convergence_tolerance*1.0d-2,saved_acc_tolerance)
      endif
      call omega_from_j(log_density,specific_angular_momentum,log_radius,log_mass, &
           shell_mass,am_transport_convective_flag,num_zones,eta_squared, &
           moment_of_inertia,omega,qiw,mean_radius)
      rot_scr%moment_of_inertia_tolerance = saved_tolerance
      star%job%acfpft = saved_acc_tolerance
!  SEARCH FOR REVERSAL OF OMEGA GRADIENTS.  IF ONE EXISTS, ENFORCE
!  SOLID-BODY ROTATION IN THE OFFENDING PAIR OF ZONES.
      zone_index = num_zones
      zone_bottom = num_zones
! (Restructured 2026: the label 20/130 zone scan became the named
! zone_scan loop; the label-70 redo loop a plain do.)
      zone_scan: do
!  POSITIVE OMEGA GRADIENT ENCOUNTERED.
      if(omega(zone_index)-omega(zone_index-1).gt.1.0d0)then
!  IF PREVIOUS GRADIENT WAS POSITIVE, LEAVE ALONE.
         if(star%old_omega(zone_index)-star%old_omega(zone_index-1).gt.1.0d-15)then
            zone_index = zone_bottom-1
            zone_bottom = zone_index
            if(zone_index.gt.1) cycle zone_scan
         else
!  SIGN OF D OMEGA/DR HAS CHANGED,INDICATING AN ERROR IN DIFFUSION.
!  MIX THE OFFENDING ZONES TO SOLID BODY ROTATION.
!  ITOP IS THE UPPERMOST UNSTABLE SHELL.
         zone_top = zone_index
!  IF ADJACENT TO A CONVECTION ZONE, MIX THE CONVECTION ZONE AS WELL.
         if(am_transport_convective_flag(zone_top) .and. zone_top.lt.num_zones) then
            do scan_index = zone_top + 1,num_zones
               if(.not.am_transport_convective_flag(scan_index)) exit
            end do
            zone_top = scan_index - 1
         endif
!  IBOT IS THE BOTTOM UNSTABLE ZONE. CHECK FOR ADJACENT CZ AS ABOVE.
         zone_bottom = zone_index - 1
         if(am_transport_convective_flag(zone_bottom) .and. zone_bottom.gt.1) then
            do scan_index = zone_bottom - 1,1,-1
               if(.not.am_transport_convective_flag(scan_index)) exit
            end do
            zone_bottom = scan_index + 1
         endif
!  ENFORCE A SOLID BODY ROTATION CURVE FROM IBOT TO ITOP.
         call solid_body_omega(log_density,specific_angular_momentum,log_radius,log_mass, &
                    shell_mass,zone_bottom,zone_top,eta_squared, &
                    moment_of_inertia,omega,qiw,mean_radius,num_zones)
!  NOW CHECK TO SEE IF THE REDISTRIBUTION HAS GENERATED ANY NEW REVERSALS.
         do
         redo_flag = .false.
!  CHECK FOR GRADIENT REVERSALS BELOW ZONE IBOT.
         if(zone_bottom.gt.1) then
            if(omega(zone_bottom)-omega(zone_bottom-1).gt.1.0d-15)then
               if(star%old_omega(zone_bottom)-star%old_omega(zone_bottom-1).lt.1.0d-15)then
                  redo_flag = .true.
                  zone_bottom = zone_bottom - 1
                  if(am_transport_convective_flag(zone_bottom) .and. zone_bottom.gt.1) then
                     do scan_index = zone_bottom - 1,1,-1
                        if(.not.am_transport_convective_flag(scan_index)) exit
                     end do
                     zone_bottom = scan_index + 1
                  endif
               endif
            endif
         endif
!  CHECK FOR GRADIENT REVERSALS ABOVE ZONE ITOP.
         if(zone_top.lt.num_zones) then
            if(omega(zone_top+1)-omega(zone_top).gt.1.0d-15)then
               if(star%old_omega(zone_top+1)-star%old_omega(zone_top).lt.1.0d-15)then
                  redo_flag = .true.
                  zone_top = zone_top+1
                  if(am_transport_convective_flag(zone_top) .and. zone_top.lt.num_zones) then
                     do scan_index = zone_top+1,num_zones
                        if(.not.am_transport_convective_flag(scan_index)) exit
                     end do
                     zone_top = scan_index - 1
                  endif
               endif
            endif
         endif
!  IF LREDO=T THEN THE REDISTRIBUTION OF ANGULAR MOMENTUM IN A REVERSED
!  REGION HAS EFFECTED A CHANGE IN OMEGA AT ONE OF BOTH OF THE BOUNDARIES
!  THAT HAS CAUSED A NEW GRADIENT REVERSAL AT THE BOUNDARY.
         if(.not. redo_flag) exit
            call solid_body_omega(log_density,specific_angular_momentum,log_radius,log_mass, &
                       shell_mass,zone_bottom,zone_top,eta_squared, &
                       moment_of_inertia,omega,qiw,mean_radius,num_zones)
         end do
         if(iteration_number.eq.star%ctrl%itdif2.or.converged_flag) &
              write(*,120)zone_bottom,zone_top,iteration_number
  120    format(5x,'OMEGA GRADIENT REVERSAL BETWEEN ZONES ', &
                 i5,' AND ',i5,' ITERATION ',i5)
         endif
      endif
!  RETURN FOR NEXT ZONE. (was label 130)
      zone_index = zone_bottom - 1
      zone_bottom = zone_index
      if(zone_index.le.1) exit zone_scan
      end do zone_scan
!  I/O FOR END OF DIFFUSION STEP.
      if(iteration_number.eq.star%ctrl%itdif2.or.converged_flag)then
!  FIND MAXIMUM FRACTIONAL CHANGE IN J/M OVER TIMESTEP.
         max_fractional_dj = 0.0d0
         max_dj_zone = 0
! MHP 10/02 ICRIT REMOVED
!         IF(ICRIT.EQ.0)THEN
           loop_start = 1
         do zone_index = loop_start,num_zones
            delta_j_fraction = &
                 (specific_angular_momentum(zone_index)- &
                 specific_angular_momentum_start(zone_index))/ &
                 specific_angular_momentum_start(zone_index)
            if(abs(delta_j_fraction).gt.abs(max_fractional_dj)) then
               max_fractional_dj = delta_j_fraction
               max_dj_zone = zone_index
            endif
         end do
         if(.not.already_converged_flag)then
            write(*,160)max_fractional_dj,max_dj_zone, &
                 (max_delta_j_by_iter(scan_index), &
                 max_delta_j_zone_by_iter(scan_index),scan_index=1, &
                 iteration_number)
  160       format(' MAX D(J/M)/(J/M)',1pe12.3,' AT PT.',i5, &
                    ' BY ITERATION'/5(1x,e11.3,i4))
         endif
!
! 2026 retire-legacy: the .FULL (FMODPT) diagnostics block that
! lived here -- per-shell omega/J tables, circulation-velocity and
! transport-coefficient listings -- was hard-disabled in 11/14
! ('G Somers: I am turning off the output to the .full file') and
! is deleted with the file itself. The model-grid transport
! coefficients are profile columns (D_omega, D_mix) instead.
      endif
      return
end subroutine check_angular_momentum
