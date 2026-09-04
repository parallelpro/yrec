!----------------------------------------------------------------------
! ptime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ptime.f; only variable names, source form, and comment
! style were updated.
!
! This subroutine limits the timestep based on the requirement that the
! changes from one model to the next in T, P, R, and L be less than the
! user parameters atime(itime_dt)-atime(itime_dl) at any model point. This criterion
! dominates in the pre-main sequence, and this subroutine must be used
! (use_structure_dt via LPTIME=T) to perform pre-main-sequence evolution.
subroutine timestep_limit_structure(previous_timestep, luminosity, log_pressure, log_radius, &
     log_temperature, num_points, struct_dt)

      use star_info_lib, only: star, json
      use controls_lib, only: itime_dl, itime_dp, itime_dr, itime_dt, itime_max_dt_frac
      implicit none

! previous_timestep: previous model timestep.
! luminosity: run of luminosity (as a function of mass) in the current
!             model. Confirmed linear (not logged): main.f updates it
!             via HL(I)=HL(I)+HL(I)*rate*dt, an Euler step only
!             consistent with a linear quantity.
! log_pressure: run of log10(pressure) in the current model. Confirmed
!               logged: ptime.f's own comments describe atime(itime_dt)-(itime_dr)
!               as "MAX DELTA LOG T/P/R".
! log_radius: run of log10(radius) in the current model.
! log_temperature: run of log10(temperature) in the current model.
! num_points: number of points in the model.
! struct_dt: timestep for requested change in P, T, R, L.
      double precision, intent(in) :: previous_timestep
      double precision, intent(in) :: luminosity(json), log_pressure(json), &
           log_radius(json), log_temperature(json)
      integer, intent(in) :: num_points
      double precision, intent(out) :: struct_dt
! max_change: maximum absolute time differences for each quantity
! (1 = P, 2 = T, 3 = R, 4 = L).
      double precision :: max_change(4)
      integer :: i
      double precision :: test_p, test_t, test_r, test_l, dt_factor, &
           dt_factor_limit

! find maximum absolute time differences for each quantity
! pressure
      max_change(1)=abs(star%logP_start(1)-log_pressure(1))
! temperature
      max_change(2)=abs(star%logT_start(1)-log_temperature(1))
! radius
      max_change(3)=abs(star%logR_start(1)-log_radius(1))
! luminosity
      if(luminosity(1)+luminosity(2).gt.0.0d0) then
       max_change(4)=abs((star%luminosity_lsun_start(1)-luminosity(1))*2.d0/(luminosity(2)+luminosity(1)))
      else
       max_change(4) = 0.0d0
      endif
      do i = 2,num_points
       test_p = abs(star%logP_start(i)-log_pressure(i))
       if(max_change(1).le.test_p) then
          max_change(1) = test_p
       endif
       test_t = abs(star%logT_start(i)-log_temperature(i))
       if(max_change(2).le.test_t) then
          max_change(2) = test_t
       endif
       test_r = abs(star%logR_start(i)-log_radius(i))
       if(max_change(3).le.test_r) then
          max_change(3) = test_r
       endif
       if(luminosity(i)+luminosity(i-1).gt.0.0d0) then
          test_l = abs((star%luminosity_lsun_start(i)-luminosity(i))*2.0d0/(luminosity(i)+luminosity(i-1)))
       else
          test_l = 0.0d0
       endif
       if(max_change(4).le.test_l) then
          max_change(4) = test_l
       endif
      end do
! now actually limit the timestep by a factor that reduces the
! time changes in all quantities to the ps values or less
      dt_factor = max_change(1)/star%ctrl%atime(itime_dp)
      if (max_change(2)/star%ctrl%atime(itime_dt).gt.dt_factor) dt_factor=max_change(2)/star%ctrl%atime(itime_dt)
      if (max_change(3)/star%ctrl%atime(itime_dr).gt.dt_factor) dt_factor=max_change(3)/star%ctrl%atime(itime_dr)
      if (max_change(4)/(star%ctrl%atime(itime_dl)*2.3026d0).gt.dt_factor) then
        dt_factor=max_change(4)/(star%ctrl%atime(itime_dl)*2.3026d0)
      endif
! if no change from previous model,set struct_dt to timestep
! stored in the previous model.
      if (dt_factor.eq.0.d0) dt_factor=1.0d0
! restrict change in timestep to no more than a factor of atime(itime_max_dt_frac)
! (the global timestep limiter) up or down.
      dt_factor_limit = star%ctrl%atime(itime_max_dt_frac)
      if (dt_factor.gt.dt_factor_limit) dt_factor=dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor=1.0d0/dt_factor_limit
      struct_dt = previous_timestep/dt_factor

      return
end subroutine timestep_limit_structure
