!----------------------------------------------------------------------
! rotation_scratch_lib
!----------------------------------------------------------------------
! The rotation solver's working state: the coefficient arrays, grids,
! iteration state and the structure at the middle of the current
! sub-step (rot_scr); the physics interpolated to the sub-step
! (mid_scr); the previous iterates of the circulation velocities
! (circ_scr). These are working arrays, not properties of the star (what
! the output or other domains read lives on star%). One module-level
! instance each; yrec_reset snapshots and restores them, so an
! in-process rerun starts from the same state.
module rotation_scratch_lib
      use star_info_lib, only: json, star
      implicit none
      private
      public :: rotation_diffusion_state, midstep_physics_state, &
           circulation_velocity_state, rot_scr, mid_scr, circ_scr, &
           disk_locking_engaged

! Row capacity of the 10-column banded system assembled by
! am_advection_diffusion_coeffs and eliminated by banded_solver
! (4*num_eq_points - 2 rows for num_eq_points equal-grid points, i.e.
! room for 2000 points; num_eq_points itself is only bounded by json
! and is not checked
! against this capacity).
      integer, parameter, public :: band_nmax = 8000


      type, public :: rotation_diffusion_state
! the commented-out Bulirsch-Stoer burn-and-mix path's tables, commented
! out with it:
!          double precision :: bs_extrapolation_table(11,15,json), &
!               bs_extrapolated_composition(15,json), &
!               bs_extrapolation_increment(15,json)
! angular momentum transport: the advective and diffusive coefficients
           double precision :: am_advective_coeff(json), &
                am_diffusive_coeff(json)
! the Eddington-Sweet advective and diffusive velocities, now and at the
! previous iteration
           double precision :: es_advective_velocity(json), &
                es_advective_velocity_prev(json), &
                es_diffusive_velocity(json), &
                es_diffusive_velocity_prev(json)
! the advection-diffusion (Maeder & Zahn) coefficients on the model and
! equal grids, and the shear and GSF diffusion coefficients. facd2/facd3
! are filled by setup/rotation_stability_setup and read by nothing.
           double precision :: facd2(json), facd3(json), &
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
! mixing: the geometric factor, the velocity estimate, the equatorial radius
           double precision :: mixing_geometric_factor(json), &
                mixing_velocity_estimate(json), &
                equatorial_radius(json)
! omega and d omega/dr at the start of the step
           double precision :: omega_avg_start(json), domega_dr_start(json)
! The equally spaced grid, set by equal_spaced_grid. Its coordinate is
!   chi = log M/dm + L/(L_tot dl) - log P/dp
! (dm, dl, dp the ranges' scales). chi holds it at the model zones of the
! transport region, eq_chi at the num_eq_points equally spaced points
! (spacing eq_chi_spacing), and eq_enclosed_mass the enclosed mass (g) at
! those points, which am_transport_grid and the composition grid
! difference into the grid's zone masses.
           double precision :: chi(json), eq_chi(json), &
                eq_enclosed_mass(json)
           double precision :: eq_chi_spacing
           integer :: num_eq_points
! d chi/dr at the grid edges
           double precision :: dchi_dr_edge(json)
! the geometric factors of the second- and third-derivative terms, on the
! equal and model grids
           double precision :: second_deriv_geom_factor_eqgrid(json), &
                third_deriv_geom_factor_eqgrid(json), &
                second_deriv_geom_factor(json), &
                third_deriv_geom_factor(json)
! the angular momentum conservation tolerance
           double precision :: moment_of_inertia_tolerance
! 4-point Lagrange interpolation weights
           double precision :: lagrange_interp_weights(4,json)
! The structure at the interfaces, set by rotation_stability_setup:
! element i is the interface between zones i-1 and i, at radius
! radius_interface(i) (cm). Density, pressure, temperature, enclosed mass
! (g) and luminosity (erg/s) are cgs; grada_interface is the adiabatic
! gradient and gradT_interface min(grad_rad, grad_ad); qdt_interface is
! d ln rho/d ln T at constant P (star%qdt's); eps_total_interface the
! energy generation (erg/g/s); gravity_factor_interface the local
! gravity factor.
           double precision :: luminosity_interface(json), &
                grada_interface(json), gradT_interface(json), &
                density_interface(json), &
                eps_total_interface(json), gravity_factor_interface(json), &
                mass_interface(json), pressure_interface(json), &
                qdt_interface(json), &
                radius_interface(json), temperature_interface(json)
! mean molecular weight, thermal diffusivity, viscosity and omega at the
! interfaces
           double precision :: mean_molecular_weight_interface(json), &
                thermal_diffusivity_interface(json), &
                kinematic_viscosity_interface(json), &
                omega_interface(json)
! the composition at the start of the step
           double precision :: composition_snapshot(15,json)
! the previous model's physics on the new mesh (setup/rezone), for the
! interpolation in time to the sub-step: the old values of mid_scr's
! fields
           double precision :: old_gradT(json), &
                old_grada(json), old_mu(json), &
                old_opacity(json), old_cp(json), &
                old_qdt(json), old_conv_vel(json), &
                old_visc(json), old_thdif(json), &
                old_eps_total(json), old_gradr(json), &
                old_eps_plus_neu_grav(json)
! theta (the Maeder & Zahn term) and omega across sub-steps.
! omega_substep_start (WMST): omega at the start of the current
! sub-step (star%old_omega on the first, star%omega afterwards), read at
! the surface by secular_transport for the wind torque. theta_rezoned
! is the previous model's theta term on the new mesh (setup/rezone).
           double precision :: theta_rezoned(json), theta_new(json), &
                theta_mean(json), &
                del_grad_diff_interface(json), &
                es_relaxation_factor(json), theta_prev(json), &
                omega_substep_start(json)
! the circulation's correction ratio
           double precision :: circulation_correction_ratio(json)
! the quadrupole moment and the local gravity
           double precision :: quadrupole_moment(json), &
                local_gravity(json)
! the Eddington-Sweet and GSF (Kippenhahn) velocity coefficients
           double precision :: es_velocity_coeff1(json), &
                es_velocity_coeff2(json), &
                gsf_kippenhahn_coeff(json), &
                es_shear_coeff(json)
! the mu-gradient (Richardson) term and the advection-diffusion shear terms
           double precision :: mu_gradient_richardson_coeff(json), &
                difad_shear_coeff1(json), difad_shear_coeff2(json)
! the convective flags at the previous sub-step
           logical :: convective_flag_prev(json)
! the opacity and energy-generation derivatives, filled by
! henyey_coefficients: dlnkappa_dlnt is logarithmic; deps_dlnt is the
! absolute d eps/d ln T (erg/g/s).
           double precision :: dlnkappa_dlnt(json), &
                deps_dlnt(json), neutrino_loss_fraction(json)
! the rotational energy term
           double precision :: rotational_energy_term(json)
! The structure interpolated to the middle of the current sub-step:
! filled by mid_timestep_model, read by secular_transport and the rest of
! the sub-step. Reads follow a write in the same step (the later
! sub-steps' seeding reads the previous sub-step's values).
! luminosity_lsun_mid is linear L/Lsun (interpolated from
! star%luminosity_lsun), not a log, unlike its log_* neighbours.
! mean_gravity_mid is star%mean_gravity's; di_domega_mid is dI/d omega.
           double precision :: eta_squared_mid(json), &
                log_density_mid(json), mean_gravity_mid(json), &
                moment_of_inertia_mid(json), luminosity_lsun_mid(json), &
                log_pressure_mid(json), log_radius_mid(json), &
                log_temperature_mid(json), omega_mid(json), &
                mean_radius_mid(json), di_domega_mid(json)
           logical :: convective_flag_mid(json), &
                am_transport_convective_flag_mid(json)
! the sub-step deuterium burning rates: mid_timestep_model reads the
! previous sub-step's before overwriting them (as locals they restarted
! from zero every sub-step; fixed in the 2026 bug sweep)
           double precision :: deuterium_rate_mid(json), &
                deuterium_rate_mid_start(json)
      end type rotation_diffusion_state
! The physics interpolated in time to the current sub-step (by
! mid_timestep_model): each field has the name of the star% field it
! comes from (opacity is star%opacity_zone's). eps_plus_neu_grav is
! eps_total plus the neutrino and gravothermal channels as stored.
      type, public :: midstep_physics_state
            double precision :: mu(json), cp(json), gradT(json)
            double precision :: grada(json), gradr(json)
            double precision :: eps_total(json), opacity(json), qdt(json)
            double precision :: thdif(json), conv_vel(json), visc(json)
            double precision :: eps_plus_neu_grav(json)
      end type midstep_physics_state
      type, public :: circulation_velocity_state
! the circulation velocities' previous iterates (the current ones are
! on star%, which the output reads)
            double precision :: es_circulation_velocity_prev(json)
            double precision :: secular_shear_velocity_prev(json)
            double precision :: gsf_circulation_velocity_prev(json)
      end type circulation_velocity_state

      type(rotation_diffusion_state), save :: rot_scr
      type(midstep_physics_state), save :: mid_scr
      type(circulation_velocity_state), save :: circ_scr

contains

! MHP 9/94 disk locking: true while the disk-locking option is on and
! the model is still younger than the disk-locking age (Gyr).
      logical function disk_locking_engaged()
      disk_locking_engaged = star%job%disk_locking_active .and. &
           star%disk_gate_age_gyr.le.star%job%disk_locking_age_gyr
      end function disk_locking_engaged

end module rotation_scratch_lib
