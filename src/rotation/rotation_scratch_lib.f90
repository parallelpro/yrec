!----------------------------------------------------------------------
! rotation_scratch_lib
!----------------------------------------------------------------------
! The rotation solver's working state: the coefficient arrays, grids,
! iteration state and the structure at the middle of the current
! sub-step (rot_scr); the physics interpolated to the sub-step
! (mix_scr); the previous iterates of the circulation velocities
! (circ_scr). These are working arrays, not properties of the star (what
! the output or other domains read lives on star%). One module-level
! instance each; yrec_reset snapshots and restores them, so an
! in-process rerun starts from the same state.
module rotation_scratch_lib
      use star_info_lib, only: json, star
      implicit none
      private
      public :: rotation_diffusion_state, mdphy_state, &
           circulation_velocity_state, rot_scr, mix_scr, circ_scr, &
           disk_locking_engaged

! Row capacity of the 10-column banded system assembled by
! am_advection_diffusion_coeffs and eliminated by banded_solver
! (4*ntot - 2 rows for ntot equal-grid points, i.e. room for 2000
! points; ntot itself is only bounded by json and is not checked
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
! the equally spaced grid in chi. eq_enclosed_mass (originally ES1): enclosed mass in g at the
! equally spaced chi points, set by equal_spaced_grid (spline of log M
! onto echi, then exp10); am_transport_grid/composition_grid difference
! it into the equal-grid zone masses.
           double precision :: chi(json), echi(json), &
                eq_enclosed_mass(json)
           double precision :: dchi
           integer :: ntot
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
! quantities at the zone interfaces
           double precision :: interface_luminosity(json), &
                delami(json), delmi(json), dm(json), &
                epsilm(json), interface_gravity_factor(json), &
                hs3(json), pm(json), qdtmi(json), &
                interface_radius(json), tm(json)
! mean molecular weight, thermal diffusivity, viscosity and omega at the
! interfaces
           double precision :: mean_molecular_weight_interface(json), &
                thermal_diffusivity_interface(json), &
                kinematic_viscosity_interface(json), &
                omega_interface(json)
! the composition at the start of the step
           double precision :: composition_snapshot(15,json)
! the previous model's physics, for the interpolation in time
           double precision :: old_delm(json), &
                old_del_adiabatic_mix(json), old_amu(json), &
                old_om(json), old_cp(json), &
                old_qdt(json), old_vel(json), &
                old_visc(json), old_thdif(json), &
                old_esum(json), old_del_radiative_mix(json), &
                old_eps(json)
! theta (the Maeder & Zahn term) and omega across sub-steps.
! omega_substep_start (WMST): omega at the start of the current
! sub-step (star%old_omega on the first, star%omega afterwards), read at
! the surface by secular_transport for the wind torque. tho is written
! by setup/rezone (the theta term on the new mesh).
           double precision :: tho(json), theta_new(json), &
                theta_mean(json), &
                del_grad_diff_interface(json), &
                es_relaxation_factor(json), theta_prev(json), &
                omega_substep_start(json)
! the potential terms of the quadrupole
           double precision :: phisp(json), phirot(json), &
                phidis(json), circulation_correction_ratio(json)
! the quadrupole moment and the local gravity
           double precision :: quadrupole_moment(json), &
                local_gravity(json)
! the spline tables of the grid interpolation
           double precision :: xval(json), yval(json), &
                xtab(json), ytab(json)
! the Eddington-Sweet and GSF (Kippenhahn) velocity coefficients
           double precision :: fact6(json), &
                es_velocity_coeff1(json), &
                es_velocity_coeff2(json), fgsfj(json), &
                gsf_kippenhahn_coeff(json), &
                es_shear_coeff(json)
! the circulation-velocity factors and the mu-gradient (Richardson) terms
           double precision :: fact1(json), fact2(json), &
                fact3(json), fact4(json), &
                mu_gradient_richardson_coeff(json), &
                difad_shear_coeff1(json), difad_shear_coeff2(json)
! the convective flags at the previous sub-step
           logical :: convective_flag_prev(json)
! the opacity and energy-generation derivatives, filled by
! henyey_coefficients: dlnkappa_* are logarithmic; deps_dlnrho and
! deps_dlnt are the absolute d eps/d ln rho and d eps/d ln T (erg/g/s).
           double precision :: dlnkappa_dlnrho(json), &
                dlnkappa_dlnt(json), deps_dlnrho(json), &
                deps_dlnt(json), neutrino_loss_fraction(json)
! the rotational energy term
           double precision :: rotational_energy_term(json)
! The structure interpolated to the middle of the current sub-step:
! filled by mid_timestep_model, read by secular_transport and the rest of
! the sub-step. Reads follow a write in the same step (the later
! sub-steps' seeding reads the previous sub-step's values).
! luminosity_lsun_mid is linear L/Lsun (interpolated from
! star%luminosity_lsun), not a log, unlike its log_* neighbours.
           double precision :: eta_squared_mid(json), &
                log_density_mid(json), hg_mid(json), &
                moment_of_inertia_mid(json), luminosity_lsun_mid(json), &
                log_pressure_mid(json), log_radius_mid(json), &
                log_temperature_mid(json), omega_mid(json), &
                mean_radius_mid(json), qiw_mid(json)
           logical :: convective_flag_mid(json), &
                am_transport_convective_flag_mid(json)
! the sub-step deuterium burning rates: mid_timestep_model reads the
! previous sub-step's before overwriting them (as locals they restarted
! from zero every sub-step; fixed in the 2026 bug sweep)
           double precision :: deuterium_rate_mid(json), &
                deuterium_rate_mid_start(json)
      end type rotation_diffusion_state
! the physics interpolated to the sub-step
      type, public :: mdphy_state
            double precision :: amum(json), cpm(json), delm(json)
            double precision :: del_adiabatic_mix(json), &
                 del_radiative_mix(json)
            double precision :: esumm(json), om(json), qdtm(json)
            double precision :: thdifm(json), velm(json), viscm(json)
            double precision :: epsm(json)
      end type mdphy_state
      type, public :: circulation_velocity_state
! the circulation velocities' previous iterates (the current ones are
! on star%, which the output reads)
            double precision :: es_circulation_velocity_prev(json)
            double precision :: secular_shear_velocity_prev(json)
            double precision :: gsf_circulation_velocity_prev(json)
      end type circulation_velocity_state

      type(rotation_diffusion_state), save :: rot_scr
      type(mdphy_state), save :: mix_scr
      type(circulation_velocity_state), save :: circ_scr

contains

! MHP 9/94 disk locking: true while the disk-locking option is on and
! the model is still younger than the disk-locking age (Gyr).
      logical function disk_locking_engaged()
      disk_locking_engaged = star%job%disk_locking_active .and. &
           star%disk_gate_age_gyr.le.star%job%disk_locking_age_gyr
      end function disk_locking_engaged

end module rotation_scratch_lib
