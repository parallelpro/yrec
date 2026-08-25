!----------------------------------------------------------------------
! entime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original entime.f; only variable names, source form, and comment
! style were updated.
!
! This subroutine limits the timestep based on the requirement that the
! changes from one model to the next in Teff and L be less than the
! size of the envelope triangle (tri_delta_teffl and tri_delta_logl).
! This prevents extrapolation outside of the triangle for small
! envelope triangles.
subroutine entime(previous_timestep, luminosity, log_teff, &
     num_points, envelope_dt)

      use star_info_lib, only: star, json
      use const_lib
      implicit none

! previous_timestep: previous model timestep.
! luminosity: run of luminosity (as a function of mass) in the current
!             model.
! log_teff: log10(effective temperature) of the current model
!           (confirmed logged: main.f computes 10.0d0**TEFFL).
! num_points: number of points in the model.
! envelope_dt: timestep for requested change in Teff, L.
      double precision, intent(in) :: previous_timestep
      double precision, intent(in) :: luminosity(json)
      double precision, intent(in) :: log_teff
      integer, intent(in) :: num_points
      double precision, intent(out) :: envelope_dt
! dt_scale: maximum absolute time differences for each quantity
! (1 = Teff, 2 = L). Only elements 1-2 are ever set, matching the
! original DIMENSION HMAX(4) (element 3-4, and the commented-out
! NHMAX(4), were unused there too).
      double precision :: dt_scale(4)
      double precision :: teffl_change, logl_change, dt_factor, &
           dt_factor_limit

! find maximum absolute time differences for each quantity
! temperature
      dt_scale(1) = star%ctrl%tri_delta_teffl
! luminosity
      dt_scale(2) = star%ctrl%tri_delta_logl

      teffl_change = abs(star%prev%log_Teff_start - log_teff)
      logl_change = abs(dlog10(star%prev%luminosity_lsun_start(star%prev%nz_start)) - &
           dlog10(luminosity(num_points)))

! now actually limit the timestep by a factor that reduces the
! time changes in all quantities to the triangle values or less
      dt_factor = teffl_change/star%ctrl%tri_delta_teffl
      if (logl_change/star%ctrl%tri_delta_logl .gt.dt_factor) dt_factor = logl_change/star%ctrl%tri_delta_logl
! if no change from previous model, set envelope_dt to timestep
! stored in the previous model.
      if (dt_factor .eq. 0.0d0) dt_factor = 1.0d0
! use atime(13) as the global factor for limiting timestep changes
      dt_factor_limit = star%ctrl%atime(13)
      if (dt_factor.gt.dt_factor_limit) dt_factor = dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor = 1.0d0/dt_factor_limit
      envelope_dt = previous_timestep/dt_factor

      return
end subroutine entime
