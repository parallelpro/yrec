!----------------------------------------------------------------------
! rotation_scratch_lib
!----------------------------------------------------------------------
! The rotation solver's working state, one record per stage of a
! sub-step:
!   substep_model_state     the structure interpolated in time to the
!                           sub-step (mid_timestep_model), with its
!                           convection zones, composition and starting
!                           omega (sub_model);
!   substep_physics_state   the physics interpolated to it (sub_physics);
!   stability_state         the structure at the interfaces, the
!                           stability coefficients and the quadrupole
!                           moment (rotation_stability_setup,
!                           compute_quadrupole), which the velocities
!                           read (stab);
!   transport_grid_state    the equally spaced grid of the current
!                           transport region (equal_spaced_grid; chi_grid);
!   circulation_velocity_state  the mechanisms' velocities at the
!                           previous iteration (circ_scr).
! These are working arrays, not properties of the star (what the output
! or other domains read, and what outlives a step, lives on star%).
! evolve_angular_momentum allocates one of each for the step and passes
! each routine the records it uses; this module holds only the types.
! Nothing in them outlives a step: every field is written within a step
! before it is read (all of them poisoned at every step start, 2026-10-13
! and 2026-10-28: every pinned deck byte-identical), so a step starts
! from unset records. Some fields carry from one sub-step to the next
! within the step (sub_model's deuterium rates and previous convective
! flags).
module rotation_scratch_lib
      use rotation_configuration, only: rotation_settings
      use star_info_lib, only: json
      use transport_mechanism_layout, only: n_transport_mechanisms
      implicit none
      private
      public :: substep_model_state, substep_physics_state, stability_state, &
           transport_grid_state, circulation_velocity_state

! (commented out with the advection-diffusion transport: the row
! capacity of the 10-column banded system of
! am_advection_diffusion_coeffs and banded_solver, 4*num_eq_points - 2
! rows, room for 2000 points, not checked against num_eq_points)
!      integer, parameter, public :: band_nmax = 8000

! ---------------------------------------------------------------
! The structure interpolated in time to the current sub-step, at the
! fraction of the step elapsed when the sub-step ends (the legacy
! midpoint model, whose fields were the "_mid" arrays): filled by
! mid_timestep_model, read by secular_transport and the rest of the
! sub-step. Reads follow a write in the same step (the later sub-steps'
! seeding reads the previous sub-step's values). Its omega is
! star%omega_mid, which outlives the step. luminosity_lsun is linear
! L/Lsun (interpolated from star%luminosity_lsun), not a log, unlike its
! log_* neighbours; mean_gravity is star%mean_gravity's; di_domega is
! zone_moments_of_inertia's, (omega/2) dI/d omega.
      type :: substep_model_state
           double precision :: eta_squared(json), &
                log_density(json), mean_gravity(json), &
                moment_of_inertia(json), luminosity_lsun(json), &
                log_pressure(json), log_radius(json), &
                log_temperature(json), &
                mean_radius(json), di_domega(json)
           logical :: convective_flag(json), &
                am_transport_convective_flag(json)
! the convective flags at the previous sub-step
           logical :: convective_flag_prev(json)
! the sub-step deuterium burning rates: mid_timestep_model reads the
! previous sub-step's before overwriting them (as locals they restarted
! from zero every sub-step; fixed in the 2026 bug sweep)
           double precision :: deuterium_rate(json), &
                deuterium_rate_start(json)
! the sub-step's burned and settled composition, which
! secular_transport restores (slots 1 to last_mu_species) at each later
! iteration
           double precision :: composition_snapshot(15,json)
! omega_start (legacy WMST): omega at the start of the current sub-step
! (star%old_omega on the first, star%omega afterwards), read at the
! surface by secular_transport for the wind torque. (The start-of-step
! physics the model is interpolated from, on the step's mesh, is the
! star's: star%old_gradT, ..., which the rezone stores.)
           double precision :: omega_start(json)
      end type substep_model_state

! ---------------------------------------------------------------
! The physics interpolated in time to the current sub-step (by
! mid_timestep_model, to the same time as sub_model): each field has the
! name of the star% field it comes from (opacity is star%opacity_zone's).
! eps_plus_neu_grav is eps_total plus the neutrino and gravothermal
! channels as stored.
      type :: substep_physics_state
            double precision :: mu(json), cp(json), gradT(json)
            double precision :: grada(json), gradr(json)
            double precision :: eps_total(json), opacity(json), qdt(json)
            double precision :: thdif(json), conv_vel(json), visc(json)
            double precision :: eps_plus_neu_grav(json)
      end type substep_physics_state

! ---------------------------------------------------------------
! The stability analysis of a sub-step, set once per sub-step by
! rotation_stability_setup and compute_quadrupole and each iteration by
! circulation_velocities' interpolation: element i is the interface
! between zones i-1 and i, at radius radius_interface(i) (cm). Density
! and luminosity (erg/s) are cgs; gravity_factor_interface is the local
! gravity factor; del_grad_diff_interface is max(grad_ad - grad +
! d ln mu/d ln P, floor). The mean molecular weight, thermal
! diffusivity, kinematic viscosity and omega at the interfaces; the
! 4-point Lagrange weights the interpolation uses; the Eddington-Sweet
! and GSF (Kippenhahn) velocity coefficients and the mu-gradient
! (Richardson) term; the quadrupole moment and the local gravity (at the
! zone centres) and the circulation's correction ratio.
      type :: stability_state
           double precision :: lagrange_interp_weights(4,json)
           double precision :: luminosity_interface(json), &
                density_interface(json), gravity_factor_interface(json), &
                radius_interface(json)
           double precision :: mean_molecular_weight_interface(json), &
                thermal_diffusivity_interface(json), &
                kinematic_viscosity_interface(json), &
                omega_interface(json)
           double precision :: del_grad_diff_interface(json)
           double precision :: es_velocity_coeff1(json), &
                es_velocity_coeff2(json), &
                gsf_kippenhahn_coeff(json), &
                es_shear_coeff(json)
           double precision :: mu_gradient_richardson_coeff(json)
           double precision :: quadrupole_moment(json), &
                local_gravity(json)
           double precision :: circulation_correction_ratio(json)
      end type stability_state

! ---------------------------------------------------------------
! The equally spaced grid of the current transport region, set by
! equal_spaced_grid. Its coordinate is
!   chi = log M/dm + L/(L_tot dl) - log P/dp
! (dm, dl, dp the ranges' scales). chi holds it at the model zones of the
! transport region, eq_chi at the num_eq_points equally spaced points
! (spacing eq_chi_spacing), and eq_enclosed_mass the enclosed mass (g) at
! those points, which am_transport_grid and the composition grid
! difference into the grid's zone masses.
      type :: transport_grid_state
           double precision :: chi(json), eq_chi(json), &
                eq_enclosed_mass(json)
           double precision :: eq_chi_spacing
           integer :: num_eq_points
! The region's quantities on the grid (am_transport_grid, or
! composition_grid for the composition): the diffusion coefficients of
! angular momentum and of the composition at the grid edges (with the
! Jacobian factor of the change to chi; legacy ECOD, ECOD2), each grid
! zone's mass (g; EM), moment of inertia, angular momentum (EJ) and
! omega; and single_interface, a region of a single unstable interface,
! which the transport skips (legacy LDUM2).
           double precision :: eq_am_diffusion_coeff(json), &
                eq_mixing_diffusion_coeff(json), eq_moment_of_inertia(json), &
                eq_angular_momentum(json), eq_mass(json), eq_omega(json)
           logical :: single_interface
      end type transport_grid_state

! ---------------------------------------------------------------
      type :: circulation_velocity_state
! each mechanism's velocity at the previous iteration (the current ones
! are the star's mechanism_velocity, which the output reads)
            double precision :: velocity_prev(n_transport_mechanisms, json)
      end type circulation_velocity_state

! ---------------------------------------------------------------
! Kept commented out by the author's decision, for the paths they
! served:
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

end module rotation_scratch_lib
