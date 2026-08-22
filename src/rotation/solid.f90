!----------------------------------------------------------------------
! solid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original solid.f; only variable names, source form, and comment
! style were updated.
!
! GIVEN THE ANGULAR MOMENTUM PER UNIT MASS, SOLID RETURNS THE ROTATION
! RATE OMEGA.  IT ALSO CALLS MOMI, WHICH CALCULATES THE MOMENT OF
! INERTIA AND R0 AND ETA2; R0 AND ETA2 ARE NEEDED BY FPFT.
subroutine solid(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, zone_start, zone_end, eta_squared, &
     moment_of_inertia, omega, di_domega, mean_radius, num_zones)

      use star_info_lib, only: star
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

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

!  GIVEN THE ANGULAR MOMENTUM PER UNIT MASS, SOLID RETURNS THE ROTATION
!  RATE OMEGA.  IT ALSO CALLS MOMI, WHICH CALCULATES THE MOMENT OF INERT
!  AND R0 AND ETA2; R0 AND ETA2 ARE NEEDED BY FPFT.
      total_angular_momentum = 0.0d0
      omega_sum = 0.0d0
      iteration_count = 0
!  FIND THE TOTAL ANGULAR MOMENTUM OF SHELLS JSTART TO JEND.
!  ALSO MAKE A FIRST GUESS AT OMEGA BY AVERAGING THE PREVIOUS VALUES.
      do zone_idx = zone_start,zone_end
         total_angular_momentum = total_angular_momentum + &
              shell_mass(zone_idx)*specific_angular_momentum(zone_idx)
         omega_sum = omega_sum + omega(zone_idx)
   10 continue
      end do
!  MHP 9/94 OPTION ADDED TO ENFORCE DISK LOCKING UP TO A GIVEN AGE IN
!  THE SURFACE C.Z. ONLY.
      disk_locked = .false.
      if(zone_end.eq.num_zones .and. zone_start.ne.zone_end .and. &
           disk_locking_active .and. &
           disk_lifetime.le.disk_temperature)disk_locked = .true.
      if(.not.disk_locked)then
         omega_guess = omega_sum/dfloat(zone_end - zone_start + 1)
      else
         omega_guess = disk_pressure
      end if
   20 do 30 zone_idx = zone_start,zone_end
         omega(zone_idx) = omega_guess
   30 continue
!  DETERMINE THE MOMENTS OF INERTIA OF SHELLS JSTART TO JEND WITH OMEGA
!  EQUAL TO WGUESS.
      call shape(log_density,log_radius,log_mass,zone_start,zone_end,omega, &
           eta_squared,mean_radius)
!       CALL MOMI(ETA2,HD,HR,HS,HS2,JSTART,JEND,OMEGA,R0,HI,QIW,M)  ! KC 2025-05-31
      call momi(eta_squared,log_radius,log_mass,shell_mass,zone_start, &
           zone_end,omega,mean_radius,moment_of_inertia,di_domega)
      total_moment_of_inertia = 0.0d0
      total_di_domega = 0.0d0
!  FIND TOTAL MOMENT OF INERTIA(CZI) AND TOTAL DI/DOMEGA (CZQIW)
      do zone_idx = zone_start,zone_end
         total_moment_of_inertia = total_moment_of_inertia + &
              moment_of_inertia(zone_idx)
         total_di_domega = total_di_domega + di_domega(zone_idx)
   40 continue
      end do
      new_angular_momentum = omega_guess*total_moment_of_inertia
!  CHECK IF THE TOTAL ANGULAR MOMENTUM FOUND WITH OMEGA = WGUESS IS CLOS
!  ENOUGH.  IF NOT, CALCULATE DELTA OMEGA AND TRY AGAIN WITH A NEW WGUES
      delta_angular_momentum = total_angular_momentum - new_angular_momentum
      if(disk_locked)goto 45
      if(dabs(delta_angular_momentum/total_angular_momentum).gt. &
           star%rot%moment_of_inertia_tolerance) then
         if(iteration_count.lt.20) then
            iteration_count = iteration_count + 1
            delta_omega = delta_angular_momentum/ &
                 (total_moment_of_inertia + omega_guess*total_di_domega)
            omega_guess = omega_guess + delta_omega
            goto 20
         end if
      end if
 45   continue
      do zone_idx = zone_start,zone_end
         specific_angular_momentum(zone_idx) = &
              moment_of_inertia(zone_idx)*omega(zone_idx)/shell_mass(zone_idx)
   50 continue
      end do

      return
end subroutine solid
