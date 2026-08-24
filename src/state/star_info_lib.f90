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
      implicit none
      private
      integer, parameter, public :: json = 5000

! ---- named indices for the positional array slots (2026) ----
! The composition array xa(15,json): species slot per row. Slot 3 is
! the TOTAL metal mass fraction Z (YREC evolves Z as one species);
! the heavy elements in slots 5..11,13..15 are diffusion/burning
! tracers carried inside Z, not independent contributions to it.
      integer, parameter, public :: &
           i_h1  = 1,  i_he4 = 2,  i_metals = 3,  i_he3 = 4, &
           i_c12 = 5,  i_c13 = 6,  i_n14 = 7,     i_n15 = 8, &
           i_o16 = 9,  i_o17 = 10, i_o18 = 11,    i_h2  = 12, &
           i_li6 = 13, i_li7 = 14, i_be9 = 15

! del_grad(3,json): the three temperature gradients per zone.
      integer, parameter, public :: &
           i_grad_rad = 1, i_grad_actual = 2, i_grad_ad = 3

! diag%seg(7,json): specific energy generation channels [erg/g/s].
      integer, parameter, public :: &
           i_eps_pp1 = 1, i_eps_pp2 = 2, i_eps_pp3 = 3, &
           i_eps_cno = 4, i_eps_he3 = 5, i_eps_neu = 6, &
           i_eps_grav = 7

! luminosity_breakdown(8): integrated luminosity per channel [Lsun].
      integer, parameter, public :: &
           i_lum_pp1 = 1, i_lum_pp2 = 2, i_lum_pp3 = 3, &
           i_lum_cno = 4, i_lum_3alpha = 5, i_lum_neu = 6, &
           i_lum_grav = 7, i_lum_he_c = 8

! flux%neutrino_flux / neutrino_flux_total(10) and
! neutrino_flux_zone(10,json): solar-neutrino source per slot
! (wrtout's SNU tables); slots 9-10 are spares.
      integer, parameter, public :: &
           i_nu_pp = 1, i_nu_pep = 2, i_nu_hep = 3, i_nu_be7 = 4, &
           i_nu_b8 = 5, i_nu_n13 = 6, i_nu_o15 = 7, i_nu_f17 = 8

! 2026 state consolidation: the sub-struct TYPE definitions moved
! here from their former one-type-per-file modules in state/ --
! the instances have lived inside star_info since the phase-6 and
! MESA-convention folds; this completes the merge (one root, one
! file). Physics-domain table state stays with its domain.

! ---- from state/oldmod_lib.f90 ----
      type, public :: prev_model_state
! 2026 MESA-convention pass: start-of-step values carry MESA's _start
! suffix, mirroring the live star_info member they snapshot.
            double precision :: logP_start(json), logT_start(json), &
                 logR_start(json), luminosity_lsun_start(json), logRho_start(json)
            double precision :: xa_start(15,json)
! NOT renamed dm_start: despite the name this slot holds LOG-MASS
! coordinates of the pre-rezoning grid (see hpoint's spline blocks) --
! a rezoning scratch, misnamed since the COMMON era.
            double precision :: old_shell_mass(json)
            logical :: convective_flag_start(json), cz_flag_start(json)
            double precision :: log_Teff_start
            integer :: nz_start
      end type prev_model_state

! ---- from state/scrtch_lib.f90 ----
      type, public :: shell_diagnostics_state
            double precision :: sesum(json), seg(7,json), sbeta(json), &
                 seta(json)
            logical :: locons(json)
            double precision :: so(json), del_grad(3,json), &
                 sfxion(3,json), svel(json), scp(json)
      end type shell_diagnostics_state

! ---- from state/turnover_lib.f90 ----
      type, public :: turnover_state
            double precision :: convective_turnover_timescale, &
                 convective_turnover_timescale_old
            double precision :: pphot, pphot0
            double precision :: fracstep
      end type turnover_state

! ---- from state/light_burn_lib.f90 ----
      type, public :: light_element_burn_state
! former common/newrat/: Li6/Li7/Be9 burning rates at the end of the
! timestep, at the current (possibly overshoot-adjusted) depth.
            double precision :: rate_li6(json), rate_li7(json), &
                 rate_be9(json)
! former common/oldrat/: same, at the start of the timestep.
            double precision :: rate_li6_start(json), rate_li7_start(json), &
                 rate_be9_start(json)
! former common/liov/: pressure scale heights used to search downward
! from the CZ base for the true (overshoot-corrected) base location.
            double precision :: pressure_scale_height_start, &
                 pressure_scale_height_end
! former common/prevcz/: previous end-of-timestep values, used as the
! new beginning-of-timestep values.
            double precision :: cz_base_radius_prev, log_rate_li6_prev, &
                 log_rate_li7_prev, log_rate_be9_prev
            integer :: envelope_cz_base_zone_prev
! former common/deuter/: deuterium burning rate (current/start of
! timestep), accreted mass fraction, and the convection-zone base zone
! index used by the deuterium-burning routines.
            double precision :: deuterium_burning_rate(json), &
                 deuterium_burning_rate_start(json)
            double precision :: accreted_mass_fraction
            integer :: jcz
      end type light_element_burn_state

! ---- from state/engeb_diag_lib.f90 ----
      type, public :: engeb_diagnostics_state
! former common/neweps/
            double precision :: alpha_capture_energy, neutrino_loss_rate
! former common/be7/
            double precision :: be7_mass_fraction
! former common/grab/: He3+He3/He3+He4 luminosity and per-shell rate
! diagnostics (JVS 10/11).
            double precision :: he3_he3_energy_rate, &
                 he3_burning_energy_rate
            double precision :: he3_he3_luminosity_zone(json), &
                 he3_burning_luminosity_zone(json)
      end type engeb_diagnostics_state

! ---- from state/fluxes_lib.f90 ----
      type, public :: neutrino_flux_state
            double precision :: neutrino_flux(10), neutrino_flux_total(10)
            double precision :: cl37_snu_rate, ga71_snu_rate
      end type neutrino_flux_state

! ---- from state/mdphy_lib.f90 ----
      type, public :: mdphy_state
            double precision :: amum(json), cpm(json), delm(json)
            double precision :: del_adiabatic_mix(json), &
                 del_radiative_mix(json)
            double precision :: esumm(json), om(json), qdtm(json)
            double precision :: thdifm(json), velm(json), viscm(json)
            double precision :: epsm(json)
      end type mdphy_state

! ---- from state/envelope_comp_lib.f90 ----
      type, public :: envelope_composition_state
            double precision :: envelope_hydrogen_fraction, &
                 envelope_metal_fraction
            double precision :: zenvm, amuenv, fxenv(12)
            double precision :: xnew, znew, stotal, senv
      end type envelope_composition_state

! ---- from state/temp_lib.f90 ----
      type, public :: shell_temp_state
            double precision :: cp(json), mean_molecular_weight(json)
            double precision :: qdt(json), thdif(json), visc(json)
      end type shell_temp_state

! ---- from state/temp2_lib.f90 ----
      type, public :: circulation_velocity_state
            double precision :: es_circulation_velocity(json), &
                 es_circulation_velocity_prev(json)
            double precision :: secular_shear_velocity(json), &
                 secular_shear_velocity_prev(json)
            double precision :: hle(json)
            double precision :: gsf_circulation_velocity(json), &
                 gsf_circulation_velocity_prev(json)
            double precision :: mu_gradient_velocity(json)
      end type circulation_velocity_state

! ---- from state/pulse_diag_lib.f90 ----
      type, public :: pulsation_diagnostics_state
! former common/pulse1/
            double precision :: pulse_dlnrho_dlnp(json), &
                 pulse_dlneps_dlnrho(json)
            double precision :: pulse_dlneps_dlnt(json), &
                 pulse_dlnkap_dlnrho(json)
            double precision :: pulse_dlnkap_dlnt(json), &
                 pulse_specific_heat(json)
            double precision :: pulse_mean_molecular_weight(json), &
                 pulse_dlnrho_dlnt(json)
            double precision :: pulse_electron_mean_molecular_weight(json)
            logical :: lpumod
! former common/pulse2/
            double precision :: qqdp, qqed, qqet, qqod, qqot, qdel, qdela, &
                 qqcp
            double precision :: qrmu, qtl, qpl, qdl, qo, qol, qt, qp
            double precision :: qqdt, qemu, qd, qfs
      end type pulsation_diagnostics_state

! ---- from state/run_diag_lib.f90 ----
      type, public :: run_diagnostics_state
! 2026 MESA-style output: per-model history sources. Computed by
! update_output_diagnostics after all physics for the step is done;
! writers (write_history) only read them. Zero at run start
! (star0 snapshot / static zero), refreshed every converged model.
           double precision :: log_R_surface, log_g_surface
           double precision :: total_moment_of_inertia, cz_moment_of_inertia
           double precision :: rotation_period_days, surf_velocity_kms
           double precision :: h_shell_bot_mass, h_shell_mid_mass, &
                h_shell_top_mass, h_shell_bot_radius, h_shell_mid_radius, &
                h_shell_top_radius
! former common/entrop/
           double precision :: temperature_entropy_term(json), &
                pressure_entropy_term(json), &
                luminosity_entropy_term(json), &
                radius_entropy_term(json)
! former common/rotprt/
           logical :: print_rotation_diagnostics
! former common/theage/
           double precision :: dage
! former common/stch/
           double precision :: composition_final(15,json), &
                log_radius_final(json), &
                log_pressure_final(json), &
                log_density_final(json), &
                log_mass_final(json), &
                log_temperature_final(json)
! former common/calsun/
           double precision :: dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
                log_l_prev, log_r_prev, delta_x, delta_alpha
           logical :: solar_calibration_active
! former common/sound/
           double precision :: adiabatic_index_gamma1(json)
           logical :: sound_speed_output_active
! former common/monte2/
           double precision :: s11_rate(1000), s33_rate(1000), s34_rate(1000), &
                s17_rate(1000), metal_to_h_ratio(1000), &
                helium_fraction_param(1000), diffusion_factor(1000), &
                luminosity_target(1000), age_target(1000)
! former common/cent/
           double precision :: central_log10_temperature, central_log10_pressure, &
                central_log10_density, envelope_mass, envelope_radius
! 2026 (phase four, step 5): output diagnostics formerly computed as
! locals inside io/wrtout.f90, now filled by the star layer
! (core/update_output_diagnostics.f90) and only READ by the writers.
           double precision :: central_beta, central_degeneracy_eta
           double precision :: core_cz_mass
           double precision :: envelope_cz_temperature, envelope_cz_density, &
                envelope_cz_pressure, envelope_cz_o16, envelope_cz_log_radius
! former common/origstart/
           double precision :: orig_specific_angular_momentum(json), &
                orig_composition(15,json)
! former common/envcz/
           double precision :: envelope_cz_base_radius_rsun, rint_placeholder
! former common/comp2/
           double precision :: envelope_helium_fraction, envelope_he3_fraction
! former common/envprt/
           double precision :: current_log10_pressure, current_log10_temperature, &
                current_log10_radius, current_log10_mass, current_log10_density, &
                current_opacity, current_beta, current_gradients(3), &
                current_ion_fraction(3), current_velocity
! former common/oldrot/
           double precision :: old_omega(json), &
                old_specific_angular_momentum(json), &
                old_moment_of_inertia(json), old_hg(json), &
                old_mean_radius(json), old_eta_squared(json)
! former common/i2o/
           character(len=4) :: initial_composition_code
      end type run_diagnostics_state

! ---- from state/rotdiff_lib.f90 ----
      type, public :: rotation_diffusion_state
! former common/advec/
           double precision :: fadv(json), fadv0(json)
! former common/bsburn/
           double precision :: bs_extrapolation_table(11,15,json), &
                bs_extrapolated_composition(15,json), &
                bs_extrapolation_increment(15,json)
! former common/burn/
           double precision :: reaction_rate_by_zone(15,json)
! former common/confac/
           double precision :: bl_radius_scale, bl_mass_scale, &
                bl_temp_scale, bl_time_scale
! former common/difad/
           double precision :: am_advective_coeff(json), &
                am_diffusive_coeff(json)
! former common/difad2/
           double precision :: es_advective_velocity(json), &
                es_advective_velocity_prev(json), &
                es_diffusive_velocity(json), &
                es_diffusive_velocity_prev(json)
! former common/difad3/
           double precision :: facd2(json), facd3(json), &
                vesd2(json), vesd3(json), &
                am_2nd_deriv_coeff(json), &
                am_3rd_deriv_coeff(json), &
                geometric_factor(json), &
                velocity_coeff0(json), &
                velocity_coeff1a(json), &
                velocity_coeff1b(json), &
                velocity_coeff2a(json), &
                velocity_coeff2b(json), &
                eq_velocity_coeff0(json), &
                eq_velocity_coeff1a(json), &
                eq_velocity_coeff1b(json), &
                eq_velocity_coeff2a(json), &
                eq_velocity_coeff2b(json), &
                shear_diffusion_coeff(json), &
                gsf_diffusion_coeff(json), &
                shear_diffusion_coeff_eqgrid(json), &
                gsf_diffusion_coeff_eqgrid(json)
! former common/difad4/
           double precision :: mixing_geometric_factor(json), &
                mixing_velocity_estimate(json), &
                equatorial_radius(json)
! former common/difaddt/
           double precision :: ethvn(json), ethvp(json), &
                omega_avg_start(json), domega_dr_start(json)
! former common/egrid/
           double precision :: chi(json), echi(json), &
                es1(json)
           double precision :: dchi
           integer :: ntot
! former common/egridchi/
           double precision :: dchi_dr_edge(json), &
                dchi_dr_center(json)
! former common/egridder/
           double precision :: second_deriv_geom_factor_eqgrid(json), &
                third_deriv_geom_factor_eqgrid(json), &
                second_deriv_geom_factor(json), &
                third_deriv_geom_factor(json)
! former common/errmom/
           double precision :: moment_of_inertia_tolerance
! former common/gravez/
           double precision :: metal_diffusion_coeff1(json), &
                metal_diffusion_coeff1_mid(json), &
                metal_diffusion_coeff2_mid(json), &
                eq_metal_diffusion_coeff1_mid(json), &
                eq_metal_diffusion_coeff2_mid(json), &
                metal_abundance_change(json), &
                metal_abundance_change_mid(json)
! former common/gscof/ -- dead-everywhere placeholder in both its
! declaring files, kept lowercased pending a confirmed source.
           double precision :: app(json), atp(json), &
                apzp(json), atzp(json)
! former common/intfac/
           double precision :: lagrange_interp_weights(4,json)
! former common/intvar/
           double precision :: interface_luminosity(json), &
                delami(json), delmi(json), dm(json), &
                epsilm(json), interface_gravity_factor(json), &
                hs3(json), pm(json), qdtmi(json), &
                interface_radius(json), tm(json)
! former common/intvr2/
           double precision :: mean_molecular_weight_interface(json), &
                thermal_diffusivity_interface(json), &
                kinematic_viscosity_interface(json), &
                omega_interface(json)
! former common/oldab/
           double precision :: composition_snapshot(15,json)
! former common/oldphy/
           double precision :: old_delm(json), &
                old_del_adiabatic_mix(json), old_amu(json), &
                old_om(json), old_cp(json), &
                old_qdt(json), old_vel(json), &
                old_visc(json), old_thdif(json), &
                old_esum(json), old_del_radiative_mix(json), &
                old_eps(json)
! former common/oldrot2/
           double precision :: tho(json), theta_new(json), &
                theta_mean(json), &
                del_grad_diff_interface(json), &
                es_relaxation_factor(json), theta_prev(json), &
                qwrst(json), wmst(json), qwrmst(json)
! former common/prevmu/
           double precision :: mu_gradient_velocity_prev(json)
! former common/pualpha/
           double precision :: alfmlt, phmlt, cmxmlt
           double precision :: valfmlt(json), vphmlt(json), &
                vcmxmlt(json)
! former common/quadd/
           double precision :: phisp(json), phirot(json), &
                phidis(json), circulation_correction_ratio(json)
! former common/quadru/
           double precision :: quadrupole_moment(json), &
                local_gravity(json)
! former common/splin/
           double precision :: xval(json), yval(json), &
                xtab(json), ytab(json)
! former common/taukh/
           double precision :: fact6(json), &
                es_velocity_coeff1(json), &
                es_velocity_coeff2(json), fgsfj(json), &
                gsf_kippenhahn_coeff(json), &
                es_shear_coeff(json)
! former common/vfact/
           double precision :: fact1(json), fact2(json), &
                fact3(json), fact4(json), &
                mu_gradient_richardson_coeff(json), &
                difad_shear_coeff1(json), difad_shear_coeff2(json)
! former common/dwmax/
           double precision :: max_domega_dr(json), &
                max_domega_dr_old(json)
! former common/gravsz/
           double precision :: src_grid_metal_diffusion_coeff1(json), &
                src_grid_metal_diffusion_coeff2(json), &
                src_grid_metal_diffusion_coeff1_dz(json), &
                src_grid_metal_diffusion_coeff2_dz(json)
! former common/prevmid/
           double precision :: del_grad_diff_prev(json), &
                del_grad_diff_new(json), radius_prev(json)
           logical :: convective_flag_prev(json)
! former common/rotder/
           double precision :: dlnkappa_dlnrho(json), &
                dlnkappa_dlnt(json), dlnepsilon_dlnrho(json), &
                dlnepsilon_dlnt(json), neutrino_loss_fraction(json)
! former common/roten/
           double precision :: rotational_energy_term(json)
! former common/masschg2/: massloss.f90 declared delta_log_pressure/
! delta_log_temperature in swapped order relative to coefft.f90/
! mdot.f90 (a self-documented, harmless pre-existing bug there since
! both are unused placeholders in massloss.f90) -- uses the majority
! (coefft.f90/mdot.f90) order here.
           double precision :: accretion_specific_entropy, &
                envelope_specific_entropy, updated_mass_msun, &
                delta_log_pressure, delta_log_temperature
! former common/masschg3/
           double precision :: solar_wind_mass_loss_rate_msun_yr, &
                wind_reference_omega, wind_max_omega
           logical :: use_rotation_scaled_solar_wind
      end type rotation_diffusion_state

! ---- from state/star_job_lib.f90 ----
      type, public :: star_job
            character(len=256) :: alex06_table_path, allard_table_path, &
                 atm_table_path, fermi_table_path, kurucz_table_path, &
                 kurucz_table2_path, laol_table_path, laol_table2_path, &
                 opal95_table_path, opal92_table_path
            character(len=256) :: zams_a_table_path, zams_b_table_path, &
                 zams_c_table_path, centre1_table_path, centre2_table_path, &
                 centre3_table_path, centre4_table_path, centre5_table_path
            character(len=256) :: opal92_table2_path, pure_z_table_path, &
                 scv_h_table_path, scv_he_table_path, scv_z_table_path
            character(len=256) :: alex95_table_paths(7)
            character(len=256) :: pulse_atm_path, pulse_env_path, &
                 pulse_mod_path
            double precision :: mixture_weights(12)
            integer :: mc_run_start, mc_run_end
      end type star_job

! ---- from state/evolve_state_lib.f90 ----
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


! Module-level (deliberately OUTSIDE the star snapshot/reset):
! evolve_step/output_diag reset flags (set by yrec_reset's
! prologue) and the libyrec namelist-path overrides.
      character(len=256), public, save :: control_nml_override = ' '
      character(len=256), public, save :: physics_nml_override = ' '
      logical, public, save :: evolve_step_reset_pending = .false.
      logical, public, save :: output_diag_reset_pending = .false.

end module star_info_lib
