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
!                           zones' radii, the stability coefficients,
!                           the quadrupole moment (rotation_stability_setup,
!                           compute_quadrupole) and the iterate's
!                           d omega/d ln r, which the velocities read
!                           (stab);
!   transport_grid_state    the equally spaced grid of the current
!                           transport region (equal_spaced_grid; chi_grid);
!   circulation_velocity_state  the mechanisms' velocities at the
!                           previous iteration and the Zahn coupling
!                           factor (circ_scr).
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
! (now in rotation/retired/advection_diffusion.f90.txt)

! ---------------------------------------------------------------
! The structure interpolated in time to the current sub-step, at the
! fraction of the step elapsed when the sub-step ends (the legacy
! midpoint model, whose fields were the "_mid" arrays): filled by
! mid_timestep_model, read by secular_transport and the rest of the
! sub-step. Reads follow a write in the same step (the later sub-steps'
! seeding reads the previous sub-step's values). Its omega is
! star%omega_substep, which outlives the step. luminosity_lsun is linear
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
! (star%omega_start on the first, star%omega afterwards), read at the
! surface by secular_transport for the wind torque. (The start-of-step
! physics the model is interpolated from, on the step's mesh, is the
! star's: star%gradT_start, ..., which the rezone stores.)
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
            double precision :: eps_total(json), opacity(json), dlnrho_dlnt(json)
            double precision :: thermal_diffusivity(json), conv_vel(json), kinematic_viscosity(json)
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
! zone centres) and the circulation's correction ratio. At the zones:
! radius, the unlogged radius (cm). At the interfaces:
! dynamical_shear_omega_limit, the |d omega/d ln r| at which Ri = 1/4
! (rotation_stability_setup), and domega_dlnr, d omega/d ln r of the
! current iterate's omega (secular_transport, each iteration).
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
           double precision :: radius(json), dynamical_shear_omega_limit(json), &
                domega_dlnr(json)
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
! the Zahn (1992) coupling factor between angular-momentum transport and
! mixing, <= 1, at each interface (zahn_coupling_factor; zero unless
! use_zahn_coupling_factor), which scales D_mix
            double precision :: zahn_coupling(json)
      end type circulation_velocity_state

! ---------------------------------------------------------------
! The fields of the commented-out paths (the Bulirsch-Stoer burn-and-mix
! and the advection-diffusion transport) are in
! rotation/retired/advection_diffusion.f90.txt.
end module rotation_scratch_lib
