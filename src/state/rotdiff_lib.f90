!----------------------------------------------------------------------
! rotdiff_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles the rotation/angular-momentum-transport
! and diffusion solver's per-model working-state COMMON blocks --
! former common/advec/, bsburn/, burn/, confac/, difad/, difad2/,
! difad3/, difad4/, difaddt/, egrid/, egridchi/, egridder/, errmom/,
! gravez/, gscof/, intfac/, intvar/, intvr2/, oldab/, oldphy/,
! oldrot2/, prevmu/, pualpha/, quadd/, quadru/, splin/, taukh/,
! vfact/ -- into one derived type, following the pulse_diag_lib/
! run_diag_lib precedent for grouping many small, unrelated,
! per-model-recomputed blocks into a single module rather than dozens
! of one-off modules. None of these touch core/parmin.f90 (confirmed);
! every declaring file across rotation/, setup/, io/, atm/, wind/,
! mixing/, misc/ used byte-identical member names/order for a given
! block, so no cross-file disambiguation was needed for this batch.
module rotdiff_lib
      implicit none
      integer, parameter :: rotdiff_json = 5000

      type, public :: rotation_diffusion_state
! former common/advec/
           double precision :: fadv(rotdiff_json), fadv0(rotdiff_json)
! former common/bsburn/
           double precision :: bs_extrapolation_table(11,15,rotdiff_json), &
                bs_extrapolated_composition(15,rotdiff_json), &
                bs_extrapolation_increment(15,rotdiff_json)
! former common/burn/
           double precision :: reaction_rate_by_zone(15,rotdiff_json)
! former common/confac/
           double precision :: bl_radius_scale, bl_mass_scale, &
                bl_temp_scale, bl_time_scale
! former common/difad/
           double precision :: am_advective_coeff(rotdiff_json), &
                am_diffusive_coeff(rotdiff_json)
! former common/difad2/
           double precision :: es_advective_velocity(rotdiff_json), &
                es_advective_velocity_prev(rotdiff_json), &
                es_diffusive_velocity(rotdiff_json), &
                es_diffusive_velocity_prev(rotdiff_json)
! former common/difad3/
           double precision :: facd2(rotdiff_json), facd3(rotdiff_json), &
                vesd2(rotdiff_json), vesd3(rotdiff_json), &
                am_2nd_deriv_coeff(rotdiff_json), &
                am_3rd_deriv_coeff(rotdiff_json), &
                geometric_factor(rotdiff_json), &
                velocity_coeff0(rotdiff_json), &
                velocity_coeff1a(rotdiff_json), &
                velocity_coeff1b(rotdiff_json), &
                velocity_coeff2a(rotdiff_json), &
                velocity_coeff2b(rotdiff_json), &
                eq_velocity_coeff0(rotdiff_json), &
                eq_velocity_coeff1a(rotdiff_json), &
                eq_velocity_coeff1b(rotdiff_json), &
                eq_velocity_coeff2a(rotdiff_json), &
                eq_velocity_coeff2b(rotdiff_json), &
                shear_diffusion_coeff(rotdiff_json), &
                gsf_diffusion_coeff(rotdiff_json), &
                shear_diffusion_coeff_eqgrid(rotdiff_json), &
                gsf_diffusion_coeff_eqgrid(rotdiff_json)
! former common/difad4/
           double precision :: mixing_geometric_factor(rotdiff_json), &
                mixing_velocity_estimate(rotdiff_json), &
                equatorial_radius(rotdiff_json)
! former common/difaddt/
           double precision :: ethvn(rotdiff_json), ethvp(rotdiff_json), &
                omega_avg_start(rotdiff_json), domega_dr_start(rotdiff_json)
! former common/egrid/
           double precision :: chi(rotdiff_json), echi(rotdiff_json), &
                es1(rotdiff_json)
           double precision :: dchi
           integer :: ntot
! former common/egridchi/
           double precision :: dchi_dr_edge(rotdiff_json), &
                dchi_dr_center(rotdiff_json)
! former common/egridder/
           double precision :: second_deriv_geom_factor_eqgrid(rotdiff_json), &
                third_deriv_geom_factor_eqgrid(rotdiff_json), &
                second_deriv_geom_factor(rotdiff_json), &
                third_deriv_geom_factor(rotdiff_json)
! former common/errmom/
           double precision :: moment_of_inertia_tolerance
! former common/gravez/
           double precision :: metal_diffusion_coeff1(rotdiff_json), &
                metal_diffusion_coeff1_mid(rotdiff_json), &
                metal_diffusion_coeff2_mid(rotdiff_json), &
                eq_metal_diffusion_coeff1_mid(rotdiff_json), &
                eq_metal_diffusion_coeff2_mid(rotdiff_json), &
                metal_abundance_change(rotdiff_json), &
                metal_abundance_change_mid(rotdiff_json)
! former common/gscof/ -- dead-everywhere placeholder in both its
! declaring files, kept lowercased pending a confirmed source.
           double precision :: app(rotdiff_json), atp(rotdiff_json), &
                apzp(rotdiff_json), atzp(rotdiff_json)
! former common/intfac/
           double precision :: lagrange_interp_weights(4,rotdiff_json)
! former common/intvar/
           double precision :: interface_luminosity(rotdiff_json), &
                delami(rotdiff_json), delmi(rotdiff_json), dm(rotdiff_json), &
                epsilm(rotdiff_json), interface_gravity_factor(rotdiff_json), &
                hs3(rotdiff_json), pm(rotdiff_json), qdtmi(rotdiff_json), &
                interface_radius(rotdiff_json), tm(rotdiff_json)
! former common/intvr2/
           double precision :: mean_molecular_weight_interface(rotdiff_json), &
                thermal_diffusivity_interface(rotdiff_json), &
                kinematic_viscosity_interface(rotdiff_json), &
                omega_interface(rotdiff_json)
! former common/oldab/
           double precision :: composition_snapshot(15,rotdiff_json)
! former common/oldphy/
           double precision :: old_delm(rotdiff_json), &
                old_del_adiabatic_mix(rotdiff_json), old_amu(rotdiff_json), &
                old_om(rotdiff_json), old_cp(rotdiff_json), &
                old_qdt(rotdiff_json), old_vel(rotdiff_json), &
                old_visc(rotdiff_json), old_thdif(rotdiff_json), &
                old_esum(rotdiff_json), old_del_radiative_mix(rotdiff_json), &
                old_eps(rotdiff_json)
! former common/oldrot2/
           double precision :: tho(rotdiff_json), theta_new(rotdiff_json), &
                theta_mean(rotdiff_json), &
                del_grad_diff_interface(rotdiff_json), &
                es_relaxation_factor(rotdiff_json), theta_prev(rotdiff_json), &
                qwrst(rotdiff_json), wmst(rotdiff_json), qwrmst(rotdiff_json)
! former common/prevmu/
           double precision :: mu_gradient_velocity_prev(rotdiff_json)
! former common/pualpha/
           double precision :: alfmlt, phmlt, cmxmlt
           double precision :: valfmlt(rotdiff_json), vphmlt(rotdiff_json), &
                vcmxmlt(rotdiff_json)
! former common/quadd/
           double precision :: phisp(rotdiff_json), phirot(rotdiff_json), &
                phidis(rotdiff_json), circulation_correction_ratio(rotdiff_json)
! former common/quadru/
           double precision :: quadrupole_moment(rotdiff_json), &
                local_gravity(rotdiff_json)
! former common/splin/
           double precision :: xval(rotdiff_json), yval(rotdiff_json), &
                xtab(rotdiff_json), ytab(rotdiff_json)
! former common/taukh/
           double precision :: fact6(rotdiff_json), &
                es_velocity_coeff1(rotdiff_json), &
                es_velocity_coeff2(rotdiff_json), fgsfj(rotdiff_json), &
                gsf_kippenhahn_coeff(rotdiff_json), &
                es_shear_coeff(rotdiff_json)
! former common/vfact/
           double precision :: fact1(rotdiff_json), fact2(rotdiff_json), &
                fact3(rotdiff_json), fact4(rotdiff_json), &
                mu_gradient_richardson_coeff(rotdiff_json), &
                difad_shear_coeff1(rotdiff_json), difad_shear_coeff2(rotdiff_json)
! former common/dwmax/
           double precision :: max_domega_dr(rotdiff_json), &
                max_domega_dr_old(rotdiff_json)
! former common/gravsz/
           double precision :: src_grid_metal_diffusion_coeff1(rotdiff_json), &
                src_grid_metal_diffusion_coeff2(rotdiff_json), &
                src_grid_metal_diffusion_coeff1_dz(rotdiff_json), &
                src_grid_metal_diffusion_coeff2_dz(rotdiff_json)
! former common/prevmid/
           double precision :: del_grad_diff_prev(rotdiff_json), &
                del_grad_diff_new(rotdiff_json), radius_prev(rotdiff_json)
           logical :: convective_flag_prev(rotdiff_json)
! former common/rotder/
           double precision :: dlnkappa_dlnrho(rotdiff_json), &
                dlnkappa_dlnt(rotdiff_json), dlnepsilon_dlnrho(rotdiff_json), &
                dlnepsilon_dlnt(rotdiff_json), neutrino_loss_fraction(rotdiff_json)
! former common/roten/
           double precision :: rotational_energy_term(rotdiff_json)
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
! 2026 (phase four, step 4 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module rotdiff_lib
