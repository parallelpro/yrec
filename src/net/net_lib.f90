!----------------------------------------------------------------------
! net_lib -- the public interface of net/
!----------------------------------------------------------------------
! Everything code outside net/ may use from the nuclear network, in one
! place (docs/GUIDELINES.md, naming rule 2; tools/check_boundaries.py
! enforces that the rest of the code uses net/ only through this
! module). It holds no code of its own: it re-exports what the rest of
! the code uses from net/'s modules, and nothing more.
!
!   net_configuration  the reaction slots and the scale slots, the
!                      run's settings record
!   network            the reaction table, the rate evaluation, the
!                      implicit burn, set_net_settings
!   nuclear_energy     the energy generation of a zone
!   light_elements     deuterium and Li/Be in the trace approximation,
!                      and the trace-element integrator
!   thermal_neutrinos  thermal neutrino losses (Itoh et al. 1996); used
!                      by nuclear_energy only, so nothing is exported
!
! The individual reaction rates, the table and the neutrino routines
! stay inside; net/test/test_net reaches them through the internal
! modules.
! net/ reads no model state and no star%ctrl: its configuration is the
! net_settings record, built at the start of every run by the star
! layer (net_settings_from_star, setup/net_settings_setup.f90) and
! installed with set_net_settings. How the pieces are
! called during a timestep: net/README.md.
!----------------------------------------------------------------------
module net_lib
! the composition rate vector (num_rate_rows rows per zone) and the
! settings record; the reaction slots (r_*: num_reactions, the size of
! the per-reaction scales) and the extra cross-section scales (xs_*:
! pep, the Be7 captures, hep; num_extra_scales); the legacy sstandard
! map
      use net_configuration, only: num_rate_rows, net_settings_t, &
           num_reactions, r_pp, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, &
           r_po16, num_extra_scales, xs_pep, xs_be7e, xs_be7p, xs_hep, &
           num_sstandard, cross_section_scales_from_sstandard
      use network, only: set_net_settings
! the implicit burn, with the network species a caller maps its
! abundances onto, and the network's dX/dt
      use network, only: compute_network_rates, burn_implicitly, n_solved, n_network, &
           network_species_slot, s_h1, s_he3, s_he4, compute_abundance_rates
! the energy generation, by chain (eps_by_chain(chain_*))
      use nuclear_energy, only: nuclear_energy_result, compute_nuclear_energy
      use network, only: n_chains, chain_name, chain_pp1, chain_pp2, chain_pp3, &
           chain_cno, chain_alpha_capture, chain_triple_alpha
! the light elements
      use light_elements, only: &
           deuterium_rate, burn_deuterium_region, surface_accretion, &
           min_deuterium_step, min_deuterium_substep, compute_li_be_rates, &
           surviving_fraction, accreted_survival, linear_in_time, &
           linear_in_log
      implicit none
      public
end module net_lib
