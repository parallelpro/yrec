!----------------------------------------------------------------------
! run_diag_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles former common/entrop/, common/rotprt/,
! common/theage/, common/stch/, common/calsun/, common/sound/,
! common/monte2/, common/cent/, common/origstart/, common/envcz/,
! common/comp2/, common/envprt/, common/oldrot/, and common/i2o/ --
! small, unrelated per-model diagnostic/state blocks that are
! recomputed or reassigned over the course of a run, not set once at
! startup -- into one derived type, following the pattern established
! by state/pulse_diag_lib.f90 for similarly-sized unrelated blocks.
module run_diag_lib
      implicit none
      integer, parameter :: run_diag_json = 5000

      type, public :: run_diagnostics_state
! former common/entrop/
           double precision :: temperature_entropy_term(run_diag_json), &
                pressure_entropy_term(run_diag_json), &
                luminosity_entropy_term(run_diag_json), &
                radius_entropy_term(run_diag_json)
! former common/rotprt/
           logical :: lprt0_placeholder
! former common/theage/
           double precision :: dage
! former common/stch/
           double precision :: composition_final(15,run_diag_json), &
                log_radius_final(run_diag_json), &
                log_pressure_final(run_diag_json), &
                log_density_final(run_diag_json), &
                log_mass_final(run_diag_json), &
                log_temperature_final(run_diag_json)
! former common/calsun/
           double precision :: dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
                log_l_prev, log_r_prev, delta_x, delta_alpha
           logical :: solar_calibration_active
! former common/sound/
           double precision :: adiabatic_index_gamma1(run_diag_json)
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
           double precision :: orig_specific_angular_momentum(run_diag_json), &
                orig_composition(15,run_diag_json)
! former common/envcz/
           double precision :: convection_zone_radius_placeholder, rint_placeholder
! former common/comp2/
           double precision :: envelope_helium_fraction, envelope_he3_fraction
! former common/envprt/
           double precision :: current_log10_pressure, current_log10_temperature, &
                current_log10_radius, current_log10_mass, current_log10_density, &
                current_opacity, current_beta, current_gradients(3), &
                current_ion_fraction(3), current_velocity
! former common/oldrot/
           double precision :: old_omega(run_diag_json), &
                old_specific_angular_momentum(run_diag_json), &
                old_moment_of_inertia(run_diag_json), old_hg(run_diag_json), &
                old_mean_radius(run_diag_json), old_eta_squared(run_diag_json)
! former common/i2o/
           character(len=4) :: initial_composition_code
      end type run_diagnostics_state
! 2026 (phase four, step 4 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module run_diag_lib
