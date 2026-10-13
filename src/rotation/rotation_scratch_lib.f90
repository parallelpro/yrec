!----------------------------------------------------------------------
! rotation_scratch_lib
!----------------------------------------------------------------------
! The rotation solver's working state: the coefficient arrays, grids,
! iteration state and the structure interpolated to the current
! sub-step (rot_scr); the physics interpolated to it (mid_scr); the
! previous iterates of the circulation velocities
! (circ_scr). These are working arrays, not properties of the star (what
! the output or other domains read, and what outlives a step, lives on
! star%). One module-level instance each. Nothing in them outlives a
! step: every field is written within a step before it is read (all
! three poisoned with NaN at every step start, 2026-10-13: every pinned
! deck byte-identical), so an in-process rerun needs no restore.
module rotation_scratch_lib
      use rotation_configuration, only: rotation_settings
      use star_info_lib, only: json
      implicit none
      private
      public :: rotation_diffusion_state, midstep_physics_state, &
           circulation_velocity_state, rot_scr, mid_scr, circ_scr

! (commented out with the advection-diffusion transport: the row
! capacity of the 10-column banded system of
! am_advection_diffusion_coeffs and banded_solver, 4*num_eq_points - 2
! rows, room for 2000 points, not checked against num_eq_points)
!      integer, parameter, public :: band_nmax = 8000


      type, public :: rotation_diffusion_state
! the commented-out Bulirsch-Stoer burn-and-mix path's tables, commented
! out with it:
!          double precision :: bs_extrapolation_table(11,15,json), &
!               bs_extrapolated_composition(15,json), &
!               bs_extrapolation_increment(15,json)
! The commented-out advection-diffusion transport's arrays (2026-10-12,
! the author's decision; see am_transport_coefficients), commented out
! with it. Those marked "unset" were never filled, which is why that
! transport never worked; facd2/facd3 were filled and read by nothing.
! The advective and diffusive coefficients on the equal grid; the
! Eddington-Sweet advective and diffusive velocities and their previous
! iterates; the Maeder & Zahn coefficients on the model grid and (unset)
! on the equal grid; the shear and GSF diffusion coefficients on the
! model grid and (unset) on the equal grid; the (unset) geometric
! factors, mixing velocity estimate and equatorial radius; d chi/dr at
! the grid edges (unset); the second- and third-derivative geometric
! factors; the theta term (Maeder & Zahn 1998) on the new mesh, now,
! mean and previous, and its relaxation factor; the order d omega/dr
! shear terms; henyey_coefficients' copies of d ln kappa/d ln T,
! d eps/d ln T and (eps - eps_nu - eps_grav)/eps.
!           double precision :: am_advective_coeff(json), &
!                am_diffusive_coeff(json)
!           double precision :: es_advective_velocity(json), &
!                es_advective_velocity_prev(json), &
!                es_diffusive_velocity(json), &
!                es_diffusive_velocity_prev(json)
!           double precision :: facd2(json), facd3(json), &
!                am_2nd_deriv_coeff(json), &
!                am_3rd_deriv_coeff(json), &
!                geometric_factor(json), &
!                velocity_coeff0(json), &
!                velocity_coeff1a(json), &
!                velocity_coeff1b(json), &
!                velocity_coeff2a(json), &
!                velocity_coeff2b(json), &
!                eq_velocity_coeff0(json), &
!                eq_velocity_coeff1a(json), &
!                eq_velocity_coeff1b(json), &
!                eq_velocity_coeff2a(json), &
!                eq_velocity_coeff2b(json), &
!                shear_diffusion_coeff(json), &
!                gsf_diffusion_coeff(json), &
!                shear_diffusion_coeff_eqgrid(json), &
!                gsf_diffusion_coeff_eqgrid(json)
!           double precision :: mixing_geometric_factor(json), &
!                mixing_velocity_estimate(json), &
!                equatorial_radius(json)
!           double precision :: dchi_dr_edge(json)
!           double precision :: second_deriv_geom_factor_eqgrid(json), &
!                third_deriv_geom_factor_eqgrid(json), &
!                second_deriv_geom_factor(json), &
!                third_deriv_geom_factor(json)
!           double precision :: theta_rezoned(json), theta_new(json), &
!                theta_mean(json), es_relaxation_factor(json), &
!                theta_prev(json)
!           double precision :: difad_shear_coeff1(json), difad_shear_coeff2(json)
!           double precision :: dlnkappa_dlnt(json), &
!                deps_dlnt(json), neutrino_loss_fraction(json)
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
! the sub-step's burned and settled composition, which
! secular_transport restores (slots 1-4) at each later iteration
           double precision :: composition_snapshot(15,json)
! (The start-of-step physics they interpolate from, on the step's mesh,
! is the star's: star%old_gradT, ..., which the rezone stores.)
! omega_substep_start (WMST): omega at the start of the current
! sub-step (star%old_omega on the first, star%omega afterwards), read at
! the surface by secular_transport for the wind torque.
! del_grad_diff_interface: max(grad_ad - grad + d ln mu/d ln P, floor)
! at the interfaces (rotation_stability_setup).
           double precision :: del_grad_diff_interface(json), &
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
! the mu-gradient (Richardson) term
           double precision :: mu_gradient_richardson_coeff(json)
! the convective flags at the previous sub-step
           logical :: convective_flag_prev(json)
! The structure interpolated in time to the current sub-step, at the
! fraction of the step elapsed when the sub-step ends (the "_mid" names
! are the legacy midpoint model's): filled by mid_timestep_model, read
! by secular_transport and the rest of the sub-step. Reads follow a
! write in the same step (the later sub-steps' seeding reads the
! previous sub-step's values).
! (Its omega is star%omega_mid, which outlives the step.)
! luminosity_lsun_mid is linear L/Lsun (interpolated from
! star%luminosity_lsun), not a log, unlike its log_* neighbours.
! mean_gravity_mid is star%mean_gravity's; di_domega_mid is
! zone_moments_of_inertia's di_domega, (omega/2) dI/d omega.
           double precision :: eta_squared_mid(json), &
                log_density_mid(json), mean_gravity_mid(json), &
                moment_of_inertia_mid(json), luminosity_lsun_mid(json), &
                log_pressure_mid(json), log_radius_mid(json), &
                log_temperature_mid(json), &
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
! mid_timestep_model, to the same time as rot_scr's *_mid): each field has the name of the star% field it
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


end module rotation_scratch_lib
