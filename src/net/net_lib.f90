!----------------------------------------------------------------------
! net_lib -- the public interface of net/
!----------------------------------------------------------------------
! Everything code outside net/ may use from the nuclear network, in one
! place (docs/GUIDELINES.md, naming rule 2; tools/check_boundaries.py
! enforces that the rest of the code uses net/ only through this
! module). It holds no code of its own: it re-exports the public parts
! of net/'s five modules.
!
!   reaction_rates     the rate of every reaction in a zone, the
!                      composition rate vector, the run's settings
!   network            the reaction table (stoichiometry, Q-values,
!                      chains), the implicit and equilibrium burns
!   nuclear_energy     the energy generation of a zone
!   light_elements     deuterium and Li/Be in the trace approximation,
!                      and the trace-element integrator
!   thermal_neutrinos  thermal neutrino losses (Itoh et al. 1996)
!
! net/ reads no model state and no star%ctrl: its configuration is the
! net_settings record, filled at the start of every run by the star
! layer (configure_network in core/burn_step.f90). How the pieces are
! called during a timestep: net/README.md.
!----------------------------------------------------------------------
module net_lib
! the reactions: slots, rates, settings
      use reaction_rates, only: &
           r_pp, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, r_po16, &
           r_ac13, r_ao16_unused, r_ac12, r_an14, r_3alpha, &
           r_c12c12_unused, num_reactions, &
           num_rate_rows, rv_frac_n15_pa, rv_frac_be7_ecap, &
           xs_pep, xs_be7e, xs_be7p, xs_hep, num_cross_section_scales, &
           iq_be7p, nuclear_rates, compute_reaction_rates, &
           compute_network_rates, deuterium_proton_capture_rate, &
           net_settings_t, net_settings, set_net_settings
! the network: species, the implicit and equilibrium burns
      use network, only: &
           s_h1, s_he3, s_he4, s_c12, s_c13, s_n14, s_o16, s_o18, &
           n_solved, n_network, mass_number, &
           burn_implicitly, compute_equilibrium_burn_rates
! the energy generation
      use nuclear_energy, only: nuclear_energy_result, compute_nuclear_energy
! the light elements
      use light_elements, only: &
           deuterium_rate, burn_deuterium_region, surface_accretion, &
           min_deuterium_step, min_deuterium_substep, compute_li_be_rates, &
           surviving_fraction, accreted_survival, linear_in_time, &
           linear_in_log
! thermal neutrino losses
      use thermal_neutrinos, only: compute_thermal_neutrino_losses, &
           compute_mixture_neutrino_losses, sneut, azbar, ifermi12, &
           zfermim12, nu_ionmax
      implicit none
      public
end module net_lib
