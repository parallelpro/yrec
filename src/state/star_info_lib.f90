!----------------------------------------------------------------------
! star_info_lib
!----------------------------------------------------------------------
! Added 2026 (phase four, step 3 -- ROADMAP.md "Phase four: the star
! layer"). The single owned representation of the stellar model,
! adapted from MESA's star_info to YREC's single-star-per-process
! design (phase-one decision: one module-level instance, no handles --
! the same pattern as every state/ module).
!
! Fields are exactly the model arrays that used to be declared in
! `program main` and threaded through 20-60-argument positional
! signatures (crrect took 60, starin 50, midmod/wrtout 43, hpoint 35,
! getw 29, mix 21) -- the measured root cause of the star-layer
! entanglement, and the source of all three intent bugs phases 2-3
! surfaced. Names and shapes are main's own, unchanged.
!
! What is deliberately NOT here (yet):
!  * model-level scalars (num_zones, log_total_mass, total_mass_msun,
!    log_teff, the CZ bookkeeping indices): still passed as arguments;
!    several are intent(inout) at rezoning and their locally-copied
!    uses need per-routine investigation before absorbing them.
!  * run control (timesteps, iteration counters, convergence
!    tolerances): driver policy, stays in arguments by design.
!  * the former-COMMON state/ modules (oldmod_lib, scrtch_lib,
!    run_diag_lib, rotdiff_lib, ...): step 4 folds the model-shaped
!    ones in; oldmod_lib becomes simply a second instance (prev).
!  * physics-domain state: eos/kap/atm/nuclear never see star_info
!    (MESA's own boundary); their facades keep plain arguments.
!
! Mid-timestep working copies (the *_mid arrays getw builds for
! seculr/rotmix) are working state of the rotation pipeline, not the
! model; routines that receive those keep explicit arguments.
module star_info_lib
      use oldmod_lib,   only: prev_model_state
      use scrtch_lib,   only: shell_diagnostics_state
      use run_diag_lib, only: run_diagnostics_state
      use rotdiff_lib,  only: rotation_diffusion_state
      implicit none
      private
      integer, parameter, public :: json = 5000

      type, public :: star_info
! structure / thermodynamics, per zone
            double precision :: log_mass(json), luminosity_lsun(json), &
                 log_radius(json), log_pressure(json), &
                 log_temperature(json), log_density(json)
            logical :: convective_flag(json)
            double precision :: composition(15,json), enclosed_mass(json), &
                 shell_mass(json), gravitational_luminosity(json)
! rotation corrections to the structure equations
            double precision :: pressure_rotation_factor(json), &
                 temperature_rotation_factor(json)
! Henyey solver work arrays and corrections
            double precision :: elim_coeff(4,2,json), elim_rhs(4,json)
            double precision :: log_pressure_delta(json), &
                 log_temperature_delta(json)
            double precision :: max_residual(4)
            integer :: max_correction_index(4)
! surface / envelope fit
            double precision :: surface_bc(6), stored_envelope_state(4), &
                 envelope_fit_coeffs(9), luminosity_breakdown(8)
            double precision :: trial_log_luminosity(3), &
                 trial_log_temperature(3), fit_point_pressure(3), &
                 fit_point_temperature(3), fit_point_radius(3)
! rotation state, per zone
            double precision :: omega(json), moment_of_inertia(json), &
                 specific_angular_momentum(json), kinetic_energy_rot(json), &
                 kinetic_energy_rot_old(json), mean_radius(json), &
                 eta_squared(json), qiw(json), mean_gravity(json)
            logical :: am_transport_convective_flag(json)
! nuclear / neutrino diagnostics, per zone
            double precision :: be7_mass_fraction_zone(json), &
                 neutrino_flux_zone(10,json)
            double precision :: reaction_rate_1(json), reaction_rate_2(json), &
                 reaction_rate_3(json), reaction_rate_4(json), &
                 reaction_rate_5(json), reaction_rate_6(json), &
                 reaction_rate_7(json), reaction_rate_8(json), &
                 reaction_rate_9(json), reaction_rate_10(json), &
                 reaction_rate_11(json), reaction_rate_12(json), &
                 reaction_rate_13(json), n15_alpha_branch_fraction(json), &
                 be7_electron_capture_fraction(json)
! mixed/radiative zone bookkeeping
            integer :: mixed_zone_bounds(12,2), &
                 mixed_zone_bounds_no_overshoot(12,2), &
                 radiative_zone_bounds(13,2)
! 2026 (phase four, step 4): the former-COMMON model-state modules,
! folded in as components -- their types stay defined in their own
! state/ files, the single instances now live here. prev is the
! previous-model store (former oldmod_lib prev_model), diag the
! per-shell diagnostics (former scrtch_lib shell_diag), run the
! run-level diagnostics (former run_diag_lib run_diag), rot the
! rotation/diffusion working state (former rotdiff_lib rot_diff).
! Physics-domain files (atm, nuclear, wind internals) that read these
! now visibly reference star%... -- the remaining star-coupling inside
! physics domains is grep-able as star% under those directories.
            type(prev_model_state) :: prev
            type(shell_diagnostics_state) :: diag
            type(run_diagnostics_state) :: run
            type(rotation_diffusion_state) :: rot
      end type star_info

! the one star this process evolves (no handles -- see header)
      type(star_info), public, save :: star

end module star_info_lib
