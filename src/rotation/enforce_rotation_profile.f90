!----------------------------------------------------------------------
! wczimp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wczimp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! JNT 2025/09/03 copy of 2025/05/14 to make the code enforce specific
! rotation profiles. Extends wcz.f90 (see there for the base solid-
! body/constant-J-per-mass/power-law behavior) with three additional,
! explicitly-selected imposed-rotation-profile modes (solid_body_mode_flag,
! originally IMPJMOD = 1, 2, or 3) used by evolve_angular_momentum.f90/omega_from_j.f90 in place
! of wcz.
subroutine enforce_rotation_profile(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, istart, iend, eta_squared, moment_of_inertia, &
     omega, qiw, mean_radius, num_zones)
      use star_info_lib, only: star, json
      use phys_const_lib
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
