!----------------------------------------------------------------------
! evolve_state_lib
!----------------------------------------------------------------------
! Added 2026 (phase five -- the embeddable engine; ROADMAP.md). The
! evolution driver's own cross-level working state: the 26 scalars
! that the model loop in run_yrec genuinely shares with its
! surrounding run/Monte-Carlo levels (timestep control, divergence
! flag, H-shell tracking, previous-model trackers for pulse-output
! pathlength, punch bookkeeping, per-run angular-momentum totals).
!
! Measured, not guessed -- twice: a first audit that counted
! declarations as "outside uses" swept 28 loop-local temporaries in
! here as well, and the Stage-0 byte-diff caught the resulting
! codegen drift in the per-zone prediction arithmetic; the corrected
! executable-lines-only audit yields exactly these 26. The 28
! region-only scalars stay run_yrec locals and will move into
! evolve_step (with save) at extraction.
!
! One module-level instance (`evo`), the star_info single-instance
! pattern. Distinct from the model (star_info), the job configuration
! (star_job), and the physics controls (const_lib): this is driver
! working state, and the phase-C reset audit will re-initialize
! exactly this (plus the model) between run_yrec calls.
module evolve_state_lib
      implicit none
      private

      type, public :: evolve_state
            logical :: has_h_shell, model_diverged_flag, punch_pending_flag, &
                 recompute_envelope_triangle, reset_triangle, &
                 saved_pulse_output_flag
            integer :: h_shell_end_index, h_shell_midpoint_zone, &
                 h_shell_zone_begin, ikut_flag, istore_flag
            double precision :: convective_velocity, dt, &
                 dt_saved, dlnrho_dlnp, dlnrho_dlnt, hydrogen_dt, &
                 max_domega_frac, path_length_sq, prev_age, prev_log_l, &
                 prev_log_teff, timestep_yr, total_angular_momentum, &
                 total_rotational_ke, trial_sign_flag
      end type evolve_state

      ! 2026 (MESA-convention pass): the single instance now lives
      ! inside star_info as star%evo.

! 2026 (phase five, step C): set by yrec_reset_lib's run prologue at
! every run_yrec entry after the first; evolve_step and
! update_output_diagnostics re-initialize their SAVEd locals to
! process-start state (static zero, plus nao's data value) when they
! see their flag, then clear it. Gives repeated run_yrec calls in one
! process the same starting state as fresh processes.
      logical, public, save :: evolve_step_reset_pending = .false.
      logical, public, save :: output_diag_reset_pending = .false.

end module evolve_state_lib
