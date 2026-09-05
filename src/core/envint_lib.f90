!----------------------------------------------------------------------
! envint_lib
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass, atm split -- ROADMAP.md): the
! envelope/atmosphere INTEGRATION, moved out of atm_lib. atm_get is
! the star solver's outer-boundary integrator: it integrates the
! atmosphere and envelope down to the fit point, driving the
! atmosphere_derivs/envelope_derivs integrands (which call the eos/
! kap facades) through numerics' bsstep in core/envint_kernel.f90,
! consuming the pure atm surface lookups (surfp/kcsurfp/alsurfp --
! still atm-domain physics) and writing model state (star%pphot,
! the pt_scr pulse scratch, the envelope/atmosphere structure
! arrays). In MESA terms this is star-solver code (create_atm/env
! integration), not the atm physics module -- which is why it now
! lives in core/ with the derivative routines, surfbc and the
! turnover diagnostics. atm/ keeps what MESA's atm has: tables,
! their init, and pure surface lookups.
module envint_lib
      implicit none

! 2026 W2: fixed Bulirsch-Stoer step-size triple (begin/min/max) a
! driver can hand atm_get in place of the star%job%{atm,env}_step_*
! values -- replaces the former save/override/restore of star%job
! around the call in rebuild_envelope, read_starting_model and
! build_stitched_model.
      type, public :: envint_step_config
            double precision :: step_begin, step_min, step_max
      end type envint_step_config

contains

! all three step sizes equal to h (the only pattern the drivers use)
pure function fixed_envint_step(h) result(steps)
      double precision, intent(in) :: h
      type(envint_step_config) :: steps
      steps%step_begin = h
      steps%step_min = h
      steps%step_max = h
end function fixed_envint_step

! 2026 W3: the surface geometry atm_get's callers derive from (L, Teff)
! before integrating -- shared by rebuild_envelope, surfbc and
! build_stitched_model.
!
! log10 of the surface radius (cm) from log10 L/Lsun and log10 Teff,
! L = 4 pi R**2 sigma Teff**4. NOTE the operand order inside the
! parenthesis: this is the order the three callers above have always
! used; observables_lib's log_r_surface_cm (and wind_lib's
! log10_radius_from_l_teff) evaluate the same identity with the
! - 4 log Teff term last, which can differ in the last bit. Unifying
! them is a numbers-changing decision (R6), so both forms are kept.
double precision function surface_log10_radius_cm(log_luminosity_lsun, log_teff)
      use star_info_lib, only: star
      use phys_const_lib
      double precision, intent(in) :: log_luminosity_lsun, log_teff
      surface_log10_radius_cm = 0.5d0*(log_luminosity_lsun + star%log10_solar_luminosity &
           - 4.0d0*log_teff - c4pil - csigl)
end function surface_log10_radius_cm

! log10 Teff of the unspotted ("ambient") surface for a spotted run
! (G Somers 10/14): the caller has checked that spots are on
! (spot_filling_factor /= 0, spot_temp_contrast /= 1) and that the
! surface shell is convective. The single-precision literals
! (0.25, 4.0, 1.0) are the original ones; build_stitched_model keeps
! its own copy with d0 literals, which is not token-identical.
double precision function ambient_log10_teff(log_teff)
      use star_info_lib, only: star
      use math_lib
      double precision, intent(in) :: log_teff
      ambient_log10_teff = log_teff - 0.25*log10(star%ctrl%spot_filling_factor * &
           pow(star%ctrl%spot_temp_contrast, 4.0) + 1.0 - star%ctrl%spot_filling_factor)
end function ambient_log10_teff

subroutine atm_get(luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, log10_star_mass, &
     print_flag, save_boundary_flag, log10_pressure_limit, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     saha_state, ierr, &
     vertex_index, stored_envelope_state, stored_vertex_index, &
     vtx_logp, vtx_logr, vtx_logt, atm_steps, env_steps)

! 2026 envint purity split phase B: atm_get is now the thin star-layer
! DRIVER. It snapshots the configuration out of star%job / star%ctrl
! into an envint_config, calls the star-blind kernel
! (core/envint_kernel.f90), and applies the star-state consequences:
! star%pphot from the integrated photosphere, and the Allard
! gray-fallback switch (atm_choice -> 0, use_ttau_relation) when the
! kernel reports it.
!
! 2026 W2 (readability, core 5): the fit-vertex bookkeeping
! (vertex_index, stored_envelope_state, stored_vertex_index,
! vtx_logp/logr/logt) is OPTIONAL and only surfbc passes it -- all
! six together; the other drivers call with save_boundary_flag
! .false., under which the kernel never writes them (it only
! compares stored_vertex_index against vertex_index), so they get
! zero-valued scratch here. atm_steps/env_steps, when present,
! replace the star%job%{atm,env}_step_* triples for this call only
! (the drivers that used to overwrite and restore star%job).
      use envint_kernel
      use star_info_lib, only: star
      implicit none

      double precision, intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity
      double precision, intent(in) :: log10_star_mass
      logical, intent(inout) :: print_flag
      logical, intent(in) :: save_boundary_flag
      double precision, intent(in) :: log10_pressure_limit
      double precision, intent(inout) :: log10_radius
      double precision, intent(inout) :: log10_teff
      double precision, intent(inout) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: saha_state
      integer, intent(out) :: ierr
      integer, intent(in), optional :: vertex_index
      double precision, intent(inout), optional :: stored_envelope_state(4)
      integer, intent(inout), optional :: stored_vertex_index
      double precision, intent(inout), optional :: vtx_logp(3), vtx_logr(3), &
           vtx_logt(3)
      type(envint_step_config), intent(in), optional :: atm_steps, env_steps

      type(envint_config) :: cfg
      logical :: switched_to_gray
      double precision :: log10_photo_pressure
! scratch stand-ins for the absent vertex group (never read back)
      double precision :: scratch_state(4), scratch_logp(3), scratch_logr(3), &
           scratch_logt(3)
      integer :: scratch_vertex_index, scratch_stored_vertex_index

      cfg%atm_choice = star%job%atm_choice
      if (present(atm_steps)) then
         cfg%atm_step_begin = atm_steps%step_begin
         cfg%atm_step_min = atm_steps%step_min
         cfg%atm_step_max = atm_steps%step_max
      else
         cfg%atm_step_begin = star%job%atm_step_begin
         cfg%atm_step_min = star%job%atm_step_min
         cfg%atm_step_max = star%job%atm_step_max
      end if
      if (present(env_steps)) then
         cfg%env_step_begin = env_steps%step_begin
         cfg%env_step_min = env_steps%step_min
         cfg%env_step_max = env_steps%step_max
      else
         cfg%env_step_begin = star%job%env_step_begin
         cfg%env_step_min = star%job%env_step_min
         cfg%env_step_max = star%job%env_step_max
      end if
      cfg%atm_error_tol = star%ctrl%atm_error_tol
      cfg%env_error_tol = star%ctrl%env_error_tol
      cfg%atm_step_initial = star%ctrl%atm_step_initial
      cfg%atm_hras = star%atm_hras
      cfg%senv = star%senv
      cfg%tenv = star%tenv
      cfg%log10_solar_luminosity = star%log10_solar_luminosity
      cfg%calc_envelope = star%job%calc_envelope_flag

      if (present(vertex_index)) then
         call run_kernel(vertex_index, stored_envelope_state, &
              stored_vertex_index, vtx_logp, vtx_logr, vtx_logt)
      else
         scratch_vertex_index = 0
         scratch_stored_vertex_index = 0
         scratch_state = 0.0d0
         scratch_logp = 0.0d0
         scratch_logr = 0.0d0
         scratch_logt = 0.0d0
         call run_kernel(scratch_vertex_index, scratch_state, &
              scratch_stored_vertex_index, scratch_logp, scratch_logr, &
              scratch_logt)
      end if

! star-state consequences of the integration, applied HERE and only
! here (the kernel never sees star_info):
      if (log10_photo_pressure > -1.0d98) star%pphot = log10_photo_pressure
      if (switched_to_gray) then
! Allard table lookup out of range: the kernel fell back to a gray
! atmosphere for this and all future calls (historical behavior)
         star%job%atm_choice = 0
         star%use_ttau_relation = .true.
      end if

contains

      subroutine run_kernel(vidx, state, sidx, logp, logr, logt)
            integer, intent(in) :: vidx
            double precision, intent(inout) :: state(4)
            integer, intent(inout) :: sidx
            double precision, intent(inout) :: logp(3), logr(3), logt(3)
            call integrate_envelope_atmosphere(cfg, switched_to_gray, &
                 log10_photo_pressure, &
                 luminosity_linear, pressure_rotation_factor, &
                 temperature_rotation_factor, log10_gravity, log10_star_mass, &
                 vidx, print_flag, save_boundary_flag, log10_pressure_limit, &
                 log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
                 state, sidx, saha_state, logp, logr, logt, ierr)
      end subroutine run_kernel

end subroutine atm_get

end module envint_lib
