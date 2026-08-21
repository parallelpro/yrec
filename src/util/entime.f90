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

      use const_lib
      implicit none
      integer, parameter :: json = 5000

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

! common/const/: not used in this file's logic; declared only to
! preserve the storage layout shared with every other file that
! references common/const/. Naming matches wrtout.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude



! common/const3/: not used in this file's logic; layout placeholder.
! Naming matches tpgrad.f90.
      double precision :: cdelrl, cmixl, cmixl2, cmixl3, clndp, &
           seconds_per_year
      common/const3/cdelrl, cmixl, cmixl2, cmixl3, clndp, seconds_per_year

! common/ctlim/: only atime(13) (global timestep-change limiter) is
! used here. Naming matches eqstat2.f90/eqburn.f90.
      double precision :: atime(14), tcut(5), saha_log10t_cutoff, &
           tenv0, tenv1, tenv, tgcut
      common/ctlim/ atime, tcut, saha_log10t_cutoff, tenv0, tenv1, tenv, tgcut

! common/oldmod/: only old_luminosity/old_teff/old_num_zones are used
! here. Naming matches eqburn.f90/dburn.f90.
      double precision :: old_pressure(json), old_temperature(json), &
           old_radius(json), old_luminosity(json), old_density(json), &
           old_composition(15,json), old_shell_mass(json), old_teff
      logical :: old_convective_flag(json), old_cz_flag(json)
      integer :: old_num_zones
      common/oldmod/ old_pressure, old_temperature, old_radius, &
           old_luminosity, old_density, old_composition, old_shell_mass, &
           old_convective_flag, old_cz_flag, old_teff, old_num_zones

! common/cenv/: tri_delta_teffl/tri_delta_logl (the envelope triangle
! half-widths in Teff and log L) are used here. Naming matches
! surfbc.f90.
      double precision :: tri_delta_teffl, tri_delta_logl, senv0_placeholder
      logical :: lsenv0_placeholder, lnew0_placeholder
      common/cenv/tri_delta_teffl, tri_delta_logl, senv0_placeholder, &
           lsenv0_placeholder, lnew0_placeholder

      save

! dt_scale: maximum absolute time differences for each quantity
! (1 = Teff, 2 = L). Only elements 1-2 are ever set, matching the
! original DIMENSION HMAX(4) (element 3-4, and the commented-out
! NHMAX(4), were unused there too).
      double precision :: dt_scale(4)
      double precision :: teffl_change, logl_change, dt_factor, &
           dt_factor_limit

! find maximum absolute time differences for each quantity
! temperature
      dt_scale(1) = tri_delta_teffl
! luminosity
      dt_scale(2) = tri_delta_logl

      teffl_change = abs(old_teff - log_teff)
      logl_change = abs(dlog10(old_luminosity(old_num_zones)) - &
           dlog10(luminosity(num_points)))

! now actually limit the timestep by a factor that reduces the
! time changes in all quantities to the triangle values or less
      dt_factor = teffl_change/tri_delta_teffl
      if (logl_change/tri_delta_logl .gt.dt_factor) dt_factor = logl_change/tri_delta_logl
! if no change from previous model, set envelope_dt to timestep
! stored in the previous model.
      if (dt_factor .eq. 0.0d0) dt_factor = 1.0d0
! use atime(13) as the global factor for limiting timestep changes
      dt_factor_limit = atime(13)
      if (dt_factor.gt.dt_factor_limit) dt_factor = dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor = 1.0d0/dt_factor_limit
      envelope_dt = previous_timestep/dt_factor

      return
end subroutine entime
