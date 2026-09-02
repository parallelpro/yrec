!----------------------------------------------------------------------
! rotation_scratch_lib
!----------------------------------------------------------------------
! 2026 solver-scratch cleanup (ROADMAP): the rotation/mixing solver
! WORKSPACE, moved out of star_info -- these are working arrays of
! the seculr/rotmix/midmod pipeline (coefficients, equal-grid
! copies, iteration state), not properties of the star. Members
! with property evidence (read by the output writers or across
! domains) were flattened onto star% instead; what remains here is
! scratch by construction. MESA's own convention keeps solver
! workspace out of s% for the same reason.
!
! Instances follow the state-module pattern (single module-level
! instances). Re-entrancy: yrec_reset snapshots/restores them
! alongside star0, so repeated in-process runs start pristine.
module rotation_scratch_lib
      use star_info_lib, only: json
      implicit none
      private
      public :: rotation_diffusion_state, mdphy_state, &
           circulation_velocity_state, rot_scr, mix_scr, circ_scr

! ---- from state/rotdiff_lib.f90 ----
      type, public :: rotation_diffusion_state
! former common/advec/
! former common/bsburn/
           double precision :: bs_extrapolation_table(11,15,json), &
                bs_extrapolated_composition(15,json), &
                bs_extrapolation_increment(15,json)
! former common/burn/
           double precision :: reaction_rate_by_zone(15,json)
! former common/confac/ bl_* scales flattened to star% (2026:
! output observables, not scratch)
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
           double precision :: omega_avg_start(json), domega_dr_start(json)
! former common/egrid/
           double precision :: chi(json), echi(json), &
                es1(json)
           double precision :: dchi
           integer :: ntot
! former common/egridchi/
           double precision :: dchi_dr_edge(json)
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
                metal_abundance_change_mid(json)
! former common/gscof/ -- dead-everywhere placeholder in both its
! declaring files, kept lowercased pending a confirmed source.
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
! former common/pualpha/ MLT alpha/phi/cmx scalars and per-zone
! vectors flattened to star% (2026: cross-domain/output properties)
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
! delta_log_temperature in swapped order relative to henyey_coefficients.f90/
! mdot.f90 (a self-documented, harmless pre-existing bug there since
! both are unused placeholders in massloss.f90) -- uses the majority
! (henyey_coefficients.f90/mdot.f90) order here.
           double precision :: accretion_specific_entropy, &
                envelope_specific_entropy, updated_mass_msun, &
                delta_log_pressure, delta_log_temperature
! former common/masschg3/
           double precision :: solar_wind_mass_loss_rate_msun_yr, &
                wind_reference_omega, wind_max_omega
           logical :: use_rotation_scaled_solar_wind
! ---- midpoint-in-time structure (2026 de-tramp, ROADMAP item 3) ----
! The stellar structure interpolated to the middle of the current
! diffusion sub-step: filled by mid_timestep_model, consumed by
! secular_transport and evolve_angular_momentum's other callees.
! Formerly 13 locals of evolve_angular_momentum trampled through both
! signatures. Reads only ever follow a write in the same driver call
! (the .not.first_call seeding reads the previous sub-step's values),
! so module lifetime changes nothing; yrec_reset snapshots it anyway.
           double precision :: eta_squared_mid(json), &
                log_density_mid(json), hg_mid(json), &
                moment_of_inertia_mid(json), log_luminosity_mid(json), &
                log_pressure_mid(json), log_radius_mid(json), &
                log_temperature_mid(json), omega_mid(json), &
                mean_radius_mid(json), qiw_mid(json)
           logical :: convective_flag_mid(json), &
                am_transport_convective_flag_mid(json)
! 2026 (bugsweep sec-11): the sub-step deuterium burning rates. The
! .not.first_call branch of mid_timestep_model reads the PREVIOUS
! sub-step's deuterium_rate_mid before overwriting it; as plain
! locals (midmod.f had a blanket SAVE) they restarted from zero on
! every sub-step after the first.
           double precision :: deuterium_rate_mid(json), &
                deuterium_rate_mid_start(json)
      end type rotation_diffusion_state
! ---- from state/mdphy_lib.f90 ----
      type, public :: mdphy_state
            double precision :: amum(json), cpm(json), delm(json)
            double precision :: del_adiabatic_mix(json), &
                 del_radiative_mix(json)
            double precision :: esumm(json), om(json), qdtm(json)
            double precision :: thdifm(json), velm(json), viscm(json)
            double precision :: epsm(json)
      end type mdphy_state
! ---- from state/temp2_lib.f90 ----
      type, public :: circulation_velocity_state
! the current-step circulation velocities (es/ss/gsf) flattened to
! star% (2026 solver-scratch cleanup: they are output observables);
! the *_prev copies below are solver iteration state and stay here.
            double precision :: es_circulation_velocity_prev(json)
            double precision :: secular_shear_velocity_prev(json)
            double precision :: hle(json)
            double precision :: gsf_circulation_velocity_prev(json)
            double precision :: mu_gradient_velocity(json)
      end type circulation_velocity_state

      type(rotation_diffusion_state), save :: rot_scr
      type(mdphy_state), save :: mix_scr
      type(circulation_velocity_state), save :: circ_scr

end module rotation_scratch_lib
