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
subroutine checkj(log_density, specific_angular_momentum_prev, &
     specific_angular_momentum_start, log_radius, log_mass, shell_mass, &
     diffusion_velocity, zone_min, zone_max, iteration_number, &
! LPRT,M,DT,ETA2,HI,HJM,IREDO,LOK,LREDO,OMEGA,  ! KC 2025-05-31
     am_transport_convective_flag, num_zones, dt, eta_squared, &
     moment_of_inertia, specific_angular_momentum, cut_count, &
     converged_flag, redo_flag, omega, &
! QIW,R0,WSAV,ID,IDM,ECOD,ECOD2,LOKAD)
     qiw, mean_radius, omega_start, print_zone_id, print_zone_count, &
     already_converged_flag)

      use run_diag_lib
      use temp2_lib
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

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



! common/errmom/: moment_of_inertia_tolerance (originally TOLERI), used
! here to relax the GETROT convergence tolerance except on the final
! diffusion iteration. Not referenced in any already-converted file.
      double precision :: moment_of_inertia_tolerance
      common/errmom/ moment_of_inertia_tolerance






! MHP 3/93
! common/quadd/: only circulation_correction_ratio (RAT) is used here.
! Naming matches vcirc.f90.
      double precision :: phisp(json), phirot(json), phidis(json), &
           circulation_correction_ratio(json)
      common/quadd/ phisp, phirot, phidis, circulation_correction_ratio

! MHP 11/94
! common/egrid/: all members used in the (dead-code) diffusion-velocity
! print block. Naming matches rotgrid.f90.
      double precision :: chi(json), echi(json), es1(json), dchi
      integer :: ntot
      common/egrid/ chi, echi, es1, dchi, ntot

! common/difad/: am_advective_coeff/am_diffusive_coeff (originally
! ECOD3/ECOD4), used in the (dead-code) print block. Naming matches
! dcoeft.f90/dadcoeft.f90/rotgrid.f90.
      double precision :: am_advective_coeff(json), am_diffusive_coeff(json)
      common/difad/ am_advective_coeff, am_diffusive_coeff

! common/difad2/: es_advective_velocity/es_advective_velocity_prev/
! es_diffusive_velocity/es_diffusive_velocity_prev (VESA/VESA0/VESD/
! VESD0), used in the (dead-code) print block. Naming matches
! rotgrid.f90/vcirc.f90.
      double precision :: es_advective_velocity(json), &
           es_advective_velocity_prev(json), es_diffusive_velocity(json), &
           es_diffusive_velocity_prev(json)
      common/difad2/ es_advective_velocity, es_advective_velocity_prev, &
           es_diffusive_velocity, es_diffusive_velocity_prev

      save

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
      converged_flag = .false.
      redo_flag = .false.
      if(already_converged_flag)then
         converged_flag = .true.
         redo_flag = .false.
      endif
      do 10 zone_index = 1,num_zones
         if(specific_angular_momentum(zone_index).le.0.0d0) then
            cut_count = cut_count + 1
!  STOP IF TIMESTEP CUT MORE THAN 3 TIMES.
!            IF(IREDO.GT.3)THEN
            if(cut_count.gt.0)then
               write(6,1000) zone_index
               write(short_file_unit,1000) zone_index
 1000          format(1x,39('>'),39('<')/5x,'ERROR IN SR CHECKJ'/ &
                       5x,'NEGATIVE J/M ENCOUNTERED IN ZONE',i5, &
                       ' AND 3 ATTEMPTS AT CUTTING TIMESTEP FAILED'/ &
                       'RUN STOPPED')
               stop
            else
               redo_flag = .true.
               dt = 0.5d0*dt
               write(6,1005)cut_count,zone_index
               write(short_file_unit,1005)cut_count,zone_index
 1005          format(5x,'ERROR IN SR CHECKJ'/5x,'TIMESTEP CUT,',1x, &
                       'NUMBER',i5,' DUE TO NEGATIVE J/M IN ZONE',i5)
               goto 240
            endif
         endif
   10 continue
!  CHECK IF THE FRACTIONAL CHANGE IN OMEGA RELATIVE TO THE PREVIOUS
!  ITERATION IS SMALL ENOUGH TO BE CONSIDERED CONVERGED.
!  ALSO LOCATE THE ZONE WHERE THE MAXIMUM CHANGE OCCURS FOR I/O.
      max_delta_j_by_iter(iteration_number) = &
           (specific_angular_momentum(1)-specific_angular_momentum_prev(1))/ &
           specific_angular_momentum_prev(1)
      specific_angular_momentum_prev(1) = specific_angular_momentum(1)
      max_delta_j_zone_by_iter(iteration_number) = 1
      do 140 zone_index = 2,num_zones
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
  140 continue
      if(abs(max_delta_j_by_iter(iteration_number)).le. &
           convergence_tolerance) then
!         LOK = .FALSE.
!      ELSE
         converged_flag = .true.
      endif
!  FIND THE RUN OF OMEGA THAT CORRESPONDS TO THE NEW RUN OF J/M.
!  THE MOMENT OF INERTIA IS A FUNCTION OF OMEGA WHICH IS SOLVED
!  ITERATIVELY.  BECAUSE THE ERROR IN THE DIFFUSION IS OF ORDER DJOK,
!  RELAX TOLERANCE FOR MOMENT OF INERTIA ITERATION EXCEPT FOR THE
!  FINAL STEP.
      saved_tolerance = moment_of_inertia_tolerance
      saved_acc_tolerance = acfpft
      if(iteration_number.lt.itdif2.and..not.converged_flag)then
         moment_of_inertia_tolerance = &
              max(convergence_tolerance*1.0d-2,saved_tolerance)
         acfpft = max(convergence_tolerance*1.0d-2,saved_acc_tolerance)
      endif
      call getrot(log_density,specific_angular_momentum,log_radius,log_mass, &
           shell_mass,am_transport_convective_flag,num_zones,eta_squared, &
           moment_of_inertia,omega,qiw,mean_radius)
      moment_of_inertia_tolerance = saved_tolerance
      acfpft = saved_acc_tolerance
!  SEARCH FOR REVERSAL OF OMEGA GRADIENTS.  IF ONE EXISTS, ENFORCE
!  SOLID-BODY ROTATION IN THE OFFENDING PAIR OF ZONES.
      zone_index = num_zones
      zone_bottom = num_zones
   20 continue
!  POSITIVE OMEGA GRADIENT ENCOUNTERED.
      if(omega(zone_index)-omega(zone_index-1).gt.1.0d0)then
!  IF PREVIOUS GRADIENT WAS POSITIVE, LEAVE ALONE.
         if(run_diag%old_omega(zone_index)-run_diag%old_omega(zone_index-1).gt.1.0d-15)then
            zone_index = zone_bottom-1
            zone_bottom = zone_index
            if(zone_index.gt.1)then
               goto 20
            else
               goto 130
            endif
         endif
!  SIGN OF D OMEGA/DR HAS CHANGED,INDICATING AN ERROR IN DIFFUSION.
!  MIX THE OFFENDING ZONES TO SOLID BODY ROTATION.
!  ITOP IS THE UPPERMOST UNSTABLE SHELL.
         zone_top = zone_index
!  IF ADJACENT TO A CONVECTION ZONE, MIX THE CONVECTION ZONE AS WELL.
         if(am_transport_convective_flag(zone_top) .and. zone_top.lt.num_zones) then
            do 30 scan_index = zone_top + 1,num_zones
               if(.not.am_transport_convective_flag(scan_index)) goto 40
   30       continue
   40       zone_top = scan_index - 1
         endif
!  IBOT IS THE BOTTOM UNSTABLE ZONE. CHECK FOR ADJACENT CZ AS ABOVE.
         zone_bottom = zone_index - 1
         if(am_transport_convective_flag(zone_bottom) .and. zone_bottom.gt.1) then
            do 50 scan_index = zone_bottom - 1,1,-1
               if(.not.am_transport_convective_flag(scan_index)) goto 60
   50       continue
   60       zone_bottom = scan_index + 1
         endif
!  ENFORCE A SOLID BODY ROTATION CURVE FROM IBOT TO ITOP.
         call solid(log_density,specific_angular_momentum,log_radius,log_mass, &
                    shell_mass,zone_bottom,zone_top,eta_squared, &
                    moment_of_inertia,omega,qiw,mean_radius,num_zones)
!  NOW CHECK TO SEE IF THE REDISTRIBUTION HAS GENERATED ANY NEW REVERSALS.
   70    continue
         redo_flag = .false.
!  CHECK FOR GRADIENT REVERSALS BELOW ZONE IBOT.
         if(zone_bottom.gt.1) then
            if(omega(zone_bottom)-omega(zone_bottom-1).gt.1.0d-15)then
               if(run_diag%old_omega(zone_bottom)-run_diag%old_omega(zone_bottom-1).lt.1.0d-15)then
                  redo_flag = .true.
                  zone_bottom = zone_bottom - 1
                  if(am_transport_convective_flag(zone_bottom) .and. zone_bottom.gt.1) then
                     do 80 scan_index = zone_bottom - 1,1,-1
                        if(.not.am_transport_convective_flag(scan_index)) goto 90
   80                continue
   90                zone_bottom = scan_index + 1
                  endif
               endif
            endif
         endif
!  CHECK FOR GRADIENT REVERSALS ABOVE ZONE ITOP.
         if(zone_top.lt.num_zones) then
            if(omega(zone_top+1)-omega(zone_top).gt.1.0d-15)then
               if(run_diag%old_omega(zone_top+1)-run_diag%old_omega(zone_top).lt.1.0d-15)then
                  redo_flag = .true.
                  zone_top = zone_top+1
                  if(am_transport_convective_flag(zone_top) .and. zone_top.lt.num_zones) then
                     do 100 scan_index = zone_top+1,num_zones
                        if(.not.am_transport_convective_flag(scan_index)) goto 110
  100                continue
  110                zone_top = scan_index - 1
                  endif
               endif
            endif
         endif
!  IF LREDO=T THEN THE REDISTRIBUTION OF ANGULAR MOMENTUM IN A REVERSED
!  REGION HAS EFFECTED A CHANGE IN OMEGA AT ONE OF BOTH OF THE BOUNDARIES
!  THAT HAS CAUSED A NEW GRADIENT REVERSAL AT THE BOUNDARY.
         if(redo_flag) then
            call solid(log_density,specific_angular_momentum,log_radius,log_mass, &
                       shell_mass,zone_bottom,zone_top,eta_squared, &
                       moment_of_inertia,omega,qiw,mean_radius,num_zones)
            goto 70
         endif
         if(iteration_number.eq.itdif2.or.converged_flag) &
              write(*,120)zone_bottom,zone_top,iteration_number
  120    format(5x,'OMEGA GRADIENT REVERSAL BETWEEN ZONES ', &
                 i5,' AND ',i5,' ITERATION ',i5)
      endif
  130 continue
!  RETURN FOR NEXT ZONE.
      zone_index = zone_bottom - 1
      zone_bottom = zone_index
      if(zone_index.gt.1) goto 20
!  I/O FOR END OF DIFFUSION STEP.
      if(iteration_number.eq.itdif2.or.converged_flag)then
!  FIND MAXIMUM FRACTIONAL CHANGE IN J/M OVER TIMESTEP.
         max_fractional_dj = 0.0d0
         max_dj_zone = 0
! MHP 10/02 ICRIT REMOVED
!         IF(ICRIT.EQ.0)THEN
           loop_start = 1
!         ELSE
!           II = ICRIT
!         ENDIF
         do 150 zone_index = loop_start,num_zones
            delta_j_fraction = &
                 (specific_angular_momentum(zone_index)- &
                 specific_angular_momentum_start(zone_index))/ &
                 specific_angular_momentum_start(zone_index)
            if(abs(delta_j_fraction).gt.abs(max_fractional_dj)) then
               max_fractional_dj = delta_j_fraction
               max_dj_zone = zone_index
            endif
  150    continue
         if(.not.already_converged_flag)then
            write(*,160)max_fractional_dj,max_dj_zone, &
                 (max_delta_j_by_iter(scan_index), &
                 max_delta_j_zone_by_iter(scan_index),scan_index=1, &
                 iteration_number)
  160       format(' MAX D(J/M)/(J/M)',1pe12.3,' AT PT.',i5, &
                    ' BY ITERATION'/5(1x,e11.3,i4))
         endif
!
! G Somers 11/14, I AM TURNING OFF THE OUTPUT TO THE .FULL FILE.
! THE AM CHANGES WILL NOT BE RECORDED, BUT THIS CAN BE TRIVIALLY
! EXTRACTED FROM THE EXTENDED .STORE FILE.
!
! SKIP OUTPUT IF NOT DESIRED.
!         IF(.NOT.LPRT)GOTO 240
         if(.true.)goto 240
! G Somers END
!
!  IF NPRTPT IS SET TO A LARGE NUMBER, SKIP DETAILED OUTPUT.
         if(print_point_interval.gt.num_zones)goto 240
         write(imodpt,170)
  170 format(' SHELL',3x,'OMEGA',5x,'DEL OMEGA',6x,'J/M',7x,'DEL J/M')
!  DETERMINE WHICH SHELLS TO PRINT.
!  FIRST POINT ALWAYS PRINTED OUT.
         print_zone_id(1) = 1
         print_zone_count = 2
         print_zone_begin = max(zone_min,print_point_interval)
         print_zone_end = min(zone_max, &
              int(zone_max/print_point_interval)*print_point_interval)
! PRINT OUT EVERY NPRTPT POINTS. WHEN V=0, SKIP POINTS.
         do 180 scan_index = print_zone_begin,print_zone_end,print_point_interval
!            IF(HV(J).EQ.0.0D0)GOTO 180
            print_zone_id(print_zone_count) = scan_index
            print_zone_count = print_zone_count + 1
  180    continue
! OUTERMOST MODEL POINT (OR POINT AT BASE OF SURFACE C.Z.)ALWAYS PRINTED.
         if(print_zone_id(print_zone_count-1).ne.zone_max)then
            print_zone_id(print_zone_count) = zone_max
         else
            print_zone_count = print_zone_count-1
         endif
!  I/O CONCERNING ANGULAR MOMENTUM TRANSPORT.
         do 200 zone_index=1,print_zone_count
            write(imodpt,190)print_zone_id(zone_index), &
                 omega(print_zone_id(zone_index)), &
                 omega(print_zone_id(zone_index))- &
                 omega_start(print_zone_id(zone_index)), &
                 specific_angular_momentum(print_zone_id(zone_index)), &
                 specific_angular_momentum(print_zone_id(zone_index))- &
                 specific_angular_momentum_start(print_zone_id(zone_index))
  190 format(1x,i5,1p4e12.3)
  200    continue
!  I/O CONCERNING DIFFUSION VELOCITIES AND SCALE LENGTHS.
         write(imodpt,210)
  210 format(1x,'SHELL',4x,'VES0',9x,'VES',7x,'VGSF0',8x,'VGSF',9x, &
              'VSS',9x,'RAT',8x,'VTOT',7x,'LENGTH',8x,'VMU')
         do 230 zone_index = 1,print_zone_count
            write(imodpt,220)print_zone_id(zone_index), &
                 circ_vel%es_circulation_velocity_prev(print_zone_id(zone_index)), &
                 circ_vel%es_circulation_velocity(print_zone_id(zone_index)), &
                 circ_vel%gsf_circulation_velocity_prev(print_zone_id(zone_index)), &
                 circ_vel%gsf_circulation_velocity(print_zone_id(zone_index)), &
                 circ_vel%secular_shear_velocity(print_zone_id(zone_index)), &
                 circulation_correction_ratio(print_zone_id(zone_index)), &
                 diffusion_velocity(print_zone_id(zone_index)), &
                 circ_vel%hle(print_zone_id(zone_index)), &
                 circ_vel%mu_gradient_velocity(print_zone_id(zone_index))
  220 format(1x,i5,1p10e12.3)
  230    continue
         if(use_diffusion_advection_transport)then
!            DO I = 1,IDM
!               WRITE(IMODPT,221)ID(I),VES(ID(I)),VESA(ID(I)),
!     *         VESD(ID(I)),
!     *         ECOD(ID(I)),ECOD2(ID(I)),ECOD3(ID(I)),ECOD4(ID(I))
! 221           FORMAT(1X,I5,1P7E12.3)
!            END DO
            if(print_zone_count.eq.ntot)then
            do zone_index = 1,print_zone_count
               write(imodpt,221)zone_index,chi(zone_index), &
                    circ_vel%es_circulation_velocity(zone_index), &
                    es_advective_velocity(zone_index), &
                    es_diffusive_velocity(zone_index),echi(zone_index), &
                    am_advective_coeff(zone_index),am_diffusive_coeff(zone_index)
 221           format(1x,i5,1p7e12.3)
            end do
            else if(print_zone_count.lt.ntot)then
            do zone_index = 1,print_zone_count
               write(imodpt,221)zone_index,chi(zone_index), &
                    circ_vel%es_circulation_velocity(zone_index), &
                    es_advective_velocity(zone_index), &
                    es_diffusive_velocity(zone_index),echi(zone_index), &
                    am_advective_coeff(zone_index),am_diffusive_coeff(zone_index)
            end do
            do zone_index = print_zone_count+1,ntot
               write(imodpt,222)zone_index,echi(zone_index), &
                    am_advective_coeff(zone_index),am_diffusive_coeff(zone_index)
 222           format(1x,i5,48x,1p3e12.3)
            end do
            else
            do zone_index = 1,ntot
               write(imodpt,221)zone_index,chi(zone_index), &
                    circ_vel%es_circulation_velocity(zone_index), &
                    es_advective_velocity(zone_index), &
                    es_diffusive_velocity(zone_index),echi(zone_index), &
                    am_advective_coeff(zone_index),am_diffusive_coeff(zone_index)
            end do
            do zone_index = ntot+1,print_zone_count
               write(imodpt,223)zone_index,chi(zone_index), &
                    circ_vel%es_circulation_velocity(zone_index), &
                    es_advective_velocity(zone_index), &
                    es_diffusive_velocity(zone_index)
 223           format(1x,i5,1p4e12.3)
            end do
            endif
         endif
      endif
  240 continue
      return
end subroutine checkj
