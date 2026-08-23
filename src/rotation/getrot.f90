!----------------------------------------------------------------------
! getrot
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getrot.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! GETROT takes the angular momentum distribution and finds the
! rotation curve that corresponds to it. Convective regions have
! solid body rotation enforced on them.
subroutine getrot(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, am_transport_convective_flag, num_zones, &
     eta_squared, moment_of_inertia, omega, qiw, mean_radius)

      implicit none
      integer, parameter :: json = 5000

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
               call solid(log_density,specific_angular_momentum,log_radius, &
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
            call wczimp(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
! NOW FIND OMEGA IN THE CURRENT SHELL
            zone_start = zone_index
            zone_end = zone_index
            call solid(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
            in_convection_zone = .false.
         else if (zone_index.eq.num_zones) then
! SPECIAL CASE OF LAST SHELL CONVECTIVE - ENFORCE SOLID-BODY ROTATION
! IN CZ AS ABOVE.
            zone_end = num_zones
! JNT 09/25 FOR 05/15 CALL WCZIMP INSTEAD OF WCZ
            call wczimp(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         endif
      end do
      return
end subroutine getrot
