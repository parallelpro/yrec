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
      use envelope_comp_lib, only: envelope_composition_state
      use turnover_lib, only: turnover_state
      use light_burn_lib, only: light_element_burn_state
      use fluxes_lib, only: neutrino_flux_state
      use engeb_diag_lib, only: engeb_diagnostics_state
      use mdphy_lib, only: mdphy_state
      use temp_lib, only: shell_temp_state
      use temp2_lib, only: circulation_velocity_state
      use pulse_diag_lib, only: pulsation_diagnostics_state
      use star_job_lib, only: star_job
      use evolve_state_lib, only: evolve_state
      implicit none
      private
      integer, parameter, public :: json = 5000

      type, public :: star_info
! structure / thermodynamics, per zone
            double precision :: log_mass(json), luminosity_lsun(json), &
                 logR(json), logP(json), &
                 logT(json), logRho(json)
            logical :: convective_flag(json)
            double precision :: xa(15,json), m(json), &
                 dm(json), gravitational_luminosity(json)
! rotation corrections to the structure equations
            double precision :: fp_rot(json), &
                 ft_rot(json)
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
            double precision :: omega(json), i_rot(json), &
                 j_rot(json), kinetic_energy_rot(json), &
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
! model-level scalars (2026, phase four follow-on): absorbed after the
! per-call-site audit the star_info header called for -- every
! converted routine receives exactly main's storage for these at every
! call site. NOT absorbed, still arguments, because call sites pass
! separate storage: mix's core_cz_edge/envelope_cz_edge (crrect passes
! its own locals; mix writes them) and its
! mixed_zone_bounds_no_overshoot (same reason, see step 3).
            integer :: nz, model_number
            integer :: core_cz_top_index, envelope_cz_bottom_index
            double precision :: log_total_mass, star_mass
            double precision :: log_Teff, log_L
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
! 2026 (phase six, step 1): the nine remaining model-state modules,
! folded in the same way as prev/diag/run/rot -- types stay in their
! own files, the single instances live here. See ROADMAP.md phase six
! for the per-module classification (envstruct/atmstruct deliberately
! stay atm-domain state and are NOT here).
            type(envelope_composition_state) :: env_comp
            type(turnover_state) :: turnover
            type(light_element_burn_state) :: light_burn
            type(neutrino_flux_state) :: flux
            type(engeb_diagnostics_state) :: engeb
            type(mdphy_state) :: mix_phys
            type(shell_temp_state) :: thermo
            type(circulation_velocity_state) :: circ
            type(pulsation_diagnostics_state) :: pulse
! job configuration and driver-step state (2026 MESA-convention pass:
! folded in as nested sub-structs -- MESA's s% job precedent; one
! root means yrec_reset's star snapshot covers them automatically)
            type(star_job) :: job
            type(evolve_state) :: evo
      end type star_info

! the one star this process evolves (no handles -- see header)
      type(star_info), public, save :: star

end module star_info_lib
