!----------------------------------------------------------------------
! envint_lib
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass, atm split -- ROADMAP.md): the
! envelope/atmosphere INTEGRATION, moved out of atm_lib. atm_get is
! the star solver's outer-boundary integrator: it integrates the
! atmosphere and envelope down to the fit point, driving the qatm/
! qenv integrands (which call the eos/kap facades) through numerics'
! bsstep, consuming the pure atm surface lookups (surfp/kcsurfp/
! alsurfp -- still atm-domain physics) and writing model state
! (star%pulse, star%run, env_comp bookkeeping, the envelope/
! atmosphere structure arrays). In MESA terms this is star-solver
! code (create_atm/env integration), not the atm physics module --
! which is why it now lives in core/ with qatm/qenv/surfbc and the
! turnover diagnostics. atm/ keeps what MESA's atm has: tables,
! their init, and pure surface lookups.
module envint_lib
      implicit none
contains

subroutine atm_get(luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, log10_star_mass, &
     vertex_index, print_flag, save_boundary_flag, log10_pressure_limit, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     stored_envelope_state, stored_vertex_index, atm_call_count, &
     env_call_count, saha_state, vtx_logp, vtx_logr, vtx_logt, &
     ierr)

! 2026 envint purity split phase B: atm_get is now the thin star-layer
! DRIVER. It snapshots the configuration out of star%job / star%ctrl
! into an envint_config, calls the star-blind kernel
! (core/envint_kernel.f90), and applies the star-state consequences:
! star%pphot from the integrated photosphere, and the Allard
! gray-fallback switch (atm_choice -> 0, use_ttau_relation) when the
! kernel reports it. The public signature is unchanged.
      use envint_kernel
      use star_info_lib, only: star
      implicit none

      double precision, intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity
      double precision, intent(in) :: log10_star_mass
      integer, intent(in) :: vertex_index
      logical, intent(inout) :: print_flag
      logical, intent(in) :: save_boundary_flag
      double precision, intent(in) :: log10_pressure_limit
      double precision, intent(inout) :: log10_radius
      double precision, intent(inout) :: log10_teff
      double precision, intent(inout) :: hydrogen_fraction, metal_fraction
      double precision, intent(inout) :: stored_envelope_state(4)
      integer, intent(inout) :: stored_vertex_index
      integer, intent(inout) :: atm_call_count, env_call_count, saha_state
      double precision, intent(inout) :: vtx_logp(3), vtx_logr(3), vtx_logt(3)
      integer, intent(out) :: ierr

      type(envint_config) :: cfg
      logical :: switched_to_gray
      double precision :: log10_photo_pressure

      cfg%atm_choice = star%job%atm_choice
      cfg%atm_step_begin = star%job%atm_step_begin
      cfg%atm_step_min = star%job%atm_step_min
      cfg%atm_step_max = star%job%atm_step_max
      cfg%env_step_begin = star%job%env_step_begin
      cfg%env_step_min = star%job%env_step_min
      cfg%env_step_max = star%job%env_step_max
      cfg%atm_error_tol = star%ctrl%atm_error_tol
      cfg%env_error_tol = star%ctrl%env_error_tol
      cfg%atm_step_initial = star%ctrl%atm_step_initial
      cfg%atm_hras = star%atm_hras
      cfg%senv = star%senv
      cfg%tenv = star%tenv
      cfg%log10_solar_luminosity = star%log10_solar_luminosity
      cfg%calc_envelope = star%job%calc_envelope_flag

      call integrate_envelope_atmosphere(cfg, switched_to_gray, &
           log10_photo_pressure, &
           luminosity_linear, pressure_rotation_factor, &
           temperature_rotation_factor, log10_gravity, log10_star_mass, &
           vertex_index, print_flag, save_boundary_flag, log10_pressure_limit, &
           log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
           stored_envelope_state, stored_vertex_index, atm_call_count, &
           env_call_count, saha_state, vtx_logp, vtx_logr, vtx_logt, &
           ierr)

! star-state consequences of the integration, applied HERE and only
! here (the kernel never sees star_info):
      if (log10_photo_pressure > -1.0d98) star%pphot = log10_photo_pressure
      if (switched_to_gray) then
! Allard table lookup out of range: the kernel fell back to a gray
! atmosphere for this and all future calls (historical behavior)
         star%job%atm_choice = 0
         star%use_ttau_relation = .true.
      end if
end subroutine atm_get

end module envint_lib
