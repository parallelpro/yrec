!----------------------------------------------------------------------
! reaction_layout
!----------------------------------------------------------------------
! The slots of net/'s per-reaction arrays, constants only, shared like
! composition_layout: by net/'s modules, which use them directly (and
! net_lib re-exports them to the star layer); by the rate-scale controls
! (physics%cross_section_scale and its derivative scales sized by
! num_reactions, physics%extra_cross_section_scale by num_extra_scales;
! io/controls_buffer.f90) and the star's run copies of the two scales
! (star_info_lib's cross_section_scale_in_use and
! extra_cross_section_scale_in_use); and by the code that fills those
! (io's map_user_inputs, core's Monte Carlo, through net_lib).
!   - the reaction slots r_*: the index of every per-reaction array (the
!     rates and their derivatives, the energies, the first num_reactions
!     rows of the composition rate vector, the per-reaction scales). Each
!     slot belongs to exactly one row of network's reaction table that
!     computes its rate (network_checks checks it), so a new reaction
!     takes a new slot: raise num_reactions;
!   - the two branching rows of the composition rate vector after the
!     reactions (rv_*);
!   - the extra cross-section scales (xs_*): the rates outside the
!     reaction table, pep, Be7 electron and proton capture, and hep;
!   - the length of the legacy cross-section vector sstandard
!     (num_sstandard; net_configuration maps it onto the scales).
module reaction_layout
      implicit none
      private
! The reaction slots.
      integer, parameter, public :: &
           r_pp = 1, r_he3he3 = 2, r_he3he4 = 3, r_pc12 = 4, r_pc13 = 5, &
           r_pn14 = 6, r_po16 = 7, r_ac13 = 8, r_ac12 = 9, r_an14 = 10, &
           r_3alpha = 11, num_reactions = 11
! The composition rate vector (network's compute_network_rates): the
! reactions at their r_* slots, per s per amu, then the N15(p,alpha) and
! Be7 electron-capture fractions.
      integer, parameter, public :: rv_frac_n15_pa = num_reactions + 1, &
           rv_frac_be7_ecap = num_reactions + 2, num_rate_rows = num_reactions + 2
! The extra cross-section scales: pep (its ratio to pp), the Be7
! electron and proton captures, and hep.
      integer, parameter, public :: xs_pep = 1, xs_be7e = 2, xs_be7p = 3, &
           xs_hep = 4, num_extra_scales = 4
! The legacy cross-section vector sstandard: 17 ratios to Bahcall &
! Ulrich (1988), used when use_new_nuclear_rates is off.
      integer, parameter, public :: num_sstandard = 17
end module reaction_layout
