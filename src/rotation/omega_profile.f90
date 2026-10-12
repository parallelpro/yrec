!----------------------------------------------------------------------
! omega_profile
!----------------------------------------------------------------------
! The rotation curve from the angular momentum: omega_from_j over the
! whole star, convection zones as solid bodies (enforce_rotation_profile)
! and radiative shells one at a time (solid_body_omega).
module omega_profile
      implicit none
      private
      public :: omega_from_j, solid_body_omega, enforce_rotation_profile

contains

! ---------------------------------------------------------------
! GETROT takes the angular momentum distribution and finds the
! rotation curve that corresponds to it. Convective regions have
! solid body rotation enforced on them.
subroutine omega_from_j(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, am_transport_convective_flag, num_zones, &
     eta_squared, moment_of_inertia, omega, qiw, mean_radius)
      use star_info_lib, only: json

      implicit none

      double precision, intent(inout) :: log_density(json), &
           specific_angular_momentum(json), log_radius(json), log_mass(json), &
           shell_mass(json)
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json), omega(json), qiw(json), mean_radius(json)
      logical :: in_convection_zone
      integer :: zone_index, zone_start, zone_end

      in_convection_zone = .false.
      do zone_index = 1,num_zones
         if (.not.in_convection_zone) then
            if (.not.am_transport_convective_flag(zone_index)) then
! RADIATIVE SHELL
               zone_start = zone_index
               zone_end = zone_index
               call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                    log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                    moment_of_inertia,omega,qiw,mean_radius,num_zones)
            else
! CONVECTION ZONE START - SET LCONVZ TO T
               in_convection_zone = .true.
               zone_start = zone_index
            endif
         else if (.not.am_transport_convective_flag(zone_index)) then
! THE PREVIOUS SHELL WAS CONVECTIVE.  IF THIS ONE IS RADIATIVE, ENFORCE
! THE APPROPRIATE ROTATION LAW IN THE ENTIRE CONVECTION ZONE PRECEDING IT.
            zone_end = zone_index - 1
! JNT 09/25 FOR 05/15 CALL WCZIMP INSTEAD OF WCZ
            call enforce_rotation_profile(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
! NOW FIND OMEGA IN THE CURRENT SHELL
            zone_start = zone_index
            zone_end = zone_index
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
            in_convection_zone = .false.
         else if (zone_index.eq.num_zones) then
! SPECIAL CASE OF LAST SHELL CONVECTIVE - ENFORCE SOLID-BODY ROTATION
! IN CZ AS ABOVE.
            zone_end = num_zones
! JNT 09/25 FOR 05/15 CALL WCZIMP INSTEAD OF WCZ
            call enforce_rotation_profile(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         endif
      end do
      return
end subroutine omega_from_j

! ---------------------------------------------------------------
! Given the angular momentum per unit mass of zones zone_start..zone_end,
! returns the single (solid-body) rotation rate omega that carries the
! same total angular momentum, iterating on omega with shape (r0, eta2)
! and zone_moments_of_inertia (I, dI/domega) at each guess. r0 and eta2
! are later needed by rotation_shape_factors.
subroutine solid_body_omega(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, zone_start, zone_end, eta_squared, &
     moment_of_inertia, omega, di_domega, mean_radius, num_zones)
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: shape
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      implicit none

      double precision, intent(in) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      double precision, intent(in) :: log_radius(json), log_mass(json), &
           shell_mass(json)
      integer, intent(in) :: zone_start, zone_end
      double precision, intent(inout) :: eta_squared(json)
      double precision, intent(out) :: moment_of_inertia(json)
      double precision, intent(inout) :: omega(json)
      double precision, intent(out) :: di_domega(json)
      double precision, intent(out) :: mean_radius(json)
      integer, intent(in) :: num_zones
! --- locals ---
      double precision :: total_angular_momentum, omega_sum
      integer :: iteration_count
      logical :: disk_locked
      double precision :: omega_guess
      double precision :: total_moment_of_inertia, total_di_domega
      double precision :: new_angular_momentum, delta_angular_momentum, &
           delta_omega
      integer :: zone_idx

      total_angular_momentum = 0.0d0
      omega_sum = 0.0d0
      iteration_count = 0
!  FIND THE TOTAL ANGULAR MOMENTUM OF SHELLS JSTART TO JEND.
!  ALSO MAKE A FIRST GUESS AT OMEGA BY AVERAGING THE PREVIOUS VALUES.
      do zone_idx = zone_start,zone_end
         total_angular_momentum = total_angular_momentum + &
              shell_mass(zone_idx)*specific_angular_momentum(zone_idx)
         omega_sum = omega_sum + omega(zone_idx)
      end do
!  MHP 9/94 OPTION ADDED TO ENFORCE DISK LOCKING UP TO A GIVEN AGE IN
!  THE SURFACE C.Z. ONLY.
      disk_locked = .false.
      if(zone_end.eq.num_zones .and. zone_start.ne.zone_end .and. &
           disk_locking_engaged())disk_locked = .true.
      if(.not.disk_locked)then
         omega_guess = omega_sum/dfloat(zone_end - zone_start + 1)
      else
         omega_guess = star%job%disk_omega_rad_s
      end if
      omega_iter: do
      do zone_idx = zone_start,zone_end
         omega(zone_idx) = omega_guess
      end do
!  DETERMINE THE MOMENTS OF INERTIA OF SHELLS JSTART TO JEND WITH OMEGA
!  EQUAL TO WGUESS.
      call shape(log_density,log_radius,log_mass,zone_start,zone_end,omega, &
           eta_squared,mean_radius)
      call zone_moments_of_inertia(eta_squared,log_radius,log_mass,shell_mass,zone_start, &
           zone_end,omega,mean_radius,moment_of_inertia,di_domega)
      total_moment_of_inertia = 0.0d0
      total_di_domega = 0.0d0
!  FIND TOTAL MOMENT OF INERTIA(CZI) AND TOTAL DI/DOMEGA (CZQIW)
      do zone_idx = zone_start,zone_end
         total_moment_of_inertia = total_moment_of_inertia + &
              moment_of_inertia(zone_idx)
         total_di_domega = total_di_domega + di_domega(zone_idx)
      end do
      new_angular_momentum = omega_guess*total_moment_of_inertia
!  CHECK IF THE TOTAL ANGULAR MOMENTUM FOUND WITH OMEGA = WGUESS IS CLOS
!  ENOUGH.  IF NOT, CALCULATE DELTA OMEGA AND TRY AGAIN WITH A NEW WGUES
      delta_angular_momentum = total_angular_momentum - new_angular_momentum
      if(.not.disk_locked)then
      if(dabs(delta_angular_momentum/total_angular_momentum).gt. &
           rot_scr%moment_of_inertia_tolerance) then
         if(iteration_count.lt.20) then
            iteration_count = iteration_count + 1
            delta_omega = delta_angular_momentum/ &
                 (total_moment_of_inertia + omega_guess*total_di_domega)
            omega_guess = omega_guess + delta_omega
            cycle omega_iter
         end if
      end if
      end if
      exit omega_iter
      end do omega_iter
      do zone_idx = zone_start,zone_end
         specific_angular_momentum(zone_idx) = &
              moment_of_inertia(zone_idx)*omega(zone_idx)/shell_mass(zone_idx)
      end do

      return
end subroutine solid_body_omega

! ---------------------------------------------------------------
! JNT 2025/09/03 copy of 2025/05/14 to make the code enforce specific
! rotation profiles. Extends the original wcz.f (base solid-body/
! constant-J-per-mass/power-law behaviour, now the mode-0 branches
! below) with three additional,
! explicitly-selected imposed-rotation-profile modes (solid_body_mode_flag,
! originally IMPJMOD = 1, 2, or 3) used by evolve_angular_momentum.f90/omega_from_j.f90 in place
! of wcz.
subroutine enforce_rotation_profile(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, istart, iend, eta_squared, moment_of_inertia, &
     omega, qiw, mean_radius, num_zones)
      use star_info_lib, only: star, json
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      double precision, intent(in) :: log_radius(json), log_mass(json)
      double precision, intent(in) :: shell_mass(json)
      integer, intent(inout) :: istart, iend
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json)
      double precision, intent(inout) :: omega(json), qiw(json), &
           mean_radius(json)
      integer, intent(in) :: num_zones
! --- locals ---
      double precision :: cz_total_am, cz_total_mass, power_law_norm
      integer :: zone_idx, zone_start, zone_end

! JNT 2025/09/03 copy of 2025/05/14 TO MAKE THE CODE ENFORCE
! SPECIFIC ROTATION PROFILES

! JNT 2025/09/03 COPY OF 2015/05/07 ADD IEND.LT.M TO ENFORCE SB
! ROTATION IN CONVECTIVE ZONES THAT ARE NOT THE SURFACE CONVECTIVE
! ZONE EVEN IF WALPCZ IS TURNED ON (I.E. WALPCZ NOW ONLY AFFECTS
! SURFACE CONVECTION ZONES)
! JVS 2026/02/13: REVERTED GT TO GE
      if((star%ctrl%walpcz.ge.0.0d0) .or. (star%ctrl%force_solid_body_rotation) .or. &
           (star%ctrl%solid_body_mode_flag.eq.1) .or. (iend.lt.num_zones)) then
!  SOLID BODY ROTATION IN CONVECTIVE REGIONS.
         call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
              log_mass,shell_mass,istart,iend,eta_squared, &
              moment_of_inertia,omega,qiw,mean_radius,num_zones)
! JNT 2025/09/03 COPY OF 2015/05/14 ADD OPTION TO MAKE EVERYTHING
! BELOW THE SURFACE CONVECTION ZONE ROTATE AS A SOLID BODY,
! DETACHED FROM WHAT IS HAPPENING IN THE CONVECTION ZONE
      else if(star%ctrl%solid_body_mode_flag.eq.2)then
         cz_total_am = specific_angular_momentum(istart)*shell_mass(istart)
         cz_total_mass = exp(ln10*star%ctrl%walpcz*log_radius(istart))* &
              moment_of_inertia(istart)
         do zone_idx = istart+1,iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + &
                 exp(ln10*star%ctrl%walpcz*log_radius(zone_idx))* &
                 moment_of_inertia(zone_idx)
         end do

!  ASSIGN NEW RUN OF J/M IN THE C.Z. AND FIND THE NEW RUN OF OMEGA.
         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            omega(zone_idx) = power_law_norm*exp(ln10*star%ctrl%walpcz* &
                 log_radius(zone_idx))
            specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                 moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         end do
         if (istart.gt.1) then
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,1,istart-1,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         endif
! JNT 2025/09/04 REPLICATING 2015/05/14 ADD OPTION TO MAKE
! EVERYTHING BELOW THE SURFACE CONVECTION ZONE ROTATE AS A
! SOLID BODY AT THE RATE OF THE BASE OF THE CONVECTION ZONE.
      else if(star%ctrl%solid_body_mode_flag.eq.3) then
         cz_total_am = specific_angular_momentum(1)*shell_mass(1)
         cz_total_mass = exp(ln10*star%ctrl%walpcz*log_radius(istart))* &
              moment_of_inertia(1)
         if (istart.gt.1) then
            do zone_idx = 2,istart
               cz_total_am = cz_total_am + &
                    specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
               cz_total_mass = cz_total_mass + &
                    exp(ln10*star%ctrl%walpcz*log_radius(istart))* &
                    moment_of_inertia(zone_idx)
            end do
         endif
         do zone_idx = (istart+1),iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + &
                 exp(ln10*star%ctrl%walpcz*log_radius(zone_idx))* &
                 moment_of_inertia(zone_idx)
         end do
!  ASSIGN NEW RUN OF J/M IN THE C.Z. AND FIND THE NEW RUN OF OMEGA.

         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            omega(zone_idx) = power_law_norm*exp(ln10*star%ctrl%walpcz* &
                 log_radius(zone_idx))
            specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                 moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         end do
         if (istart.gt.1) then
            do zone_idx = 1,(istart-1)
               omega(zone_idx) = omega(istart)
               specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                    moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            end do
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,1,istart-1,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         endif

      else if((star%ctrl%walpcz.le.-2.0d0) .and. (star%ctrl%solid_body_mode_flag.eq.0))then
!  CONSTANT SPECIFIC ANGULAR MOMENTUM PER UNIT MASS IN THE C.Z.
!  FIND TOTAL MASS AND ANGULAR MOMENTUM OF C.Z.
         cz_total_am = specific_angular_momentum(istart)*shell_mass(istart)
         cz_total_mass = shell_mass(istart)
         do zone_idx = istart+1,iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + shell_mass(zone_idx)
         end do
!  ASSIGN NEW RUN OF J/M IN THE C.Z. AND FIND THE NEW RUN OF OMEGA.
         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            specific_angular_momentum(zone_idx) = power_law_norm
            omega(zone_idx) = power_law_norm*shell_mass(zone_idx)/ &
                 moment_of_inertia(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         end do
      else
!  GENERAL LAW FOR OMEGA IN C.Z.: OMEGA = C*R**WALPCZ,WHERE C IS A CONSTANT
!  FOR THE ENTIRE C.Z. IF THIS HOLDS, THE RUN OF J/M CAN BE FOUND BY
!  SOLVING FOR THE CONSTANT BY REQUIRING THAT THE TOTAL J OF THE C.Z. BE
!  CONSERVED AND THEN USING J/M=C*R**WALPCZ*(MOMENT OF INERTIA PER UNIT
!  MASS).
!  FIND TOTAL MASS AND ANGULAR MOMENTUM OF C.Z.
         cz_total_am = specific_angular_momentum(istart)*shell_mass(istart)
         cz_total_mass = exp(ln10*star%ctrl%walpcz*log_radius(istart))* &
              moment_of_inertia(istart)
         do zone_idx = istart+1,iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + &
                 exp(ln10*star%ctrl%walpcz*log_radius(zone_idx))* &
                 moment_of_inertia(zone_idx)
         end do
!  ASSIGN NEW RUN OF J/M IN THE C.Z. AND FIND THE NEW RUN OF OMEGA.
         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            omega(zone_idx) = power_law_norm*exp(ln10*star%ctrl%walpcz* &
                 log_radius(zone_idx))
            specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                 moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         end do
      endif
      return
end subroutine enforce_rotation_profile

end module omega_profile
