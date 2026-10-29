!----------------------------------------------------------------------
! net_configuration
!----------------------------------------------------------------------
! The network's layout and the run's settings, in one place:
!   - the reaction slots r_*: the index of every per-reaction array (the
!     rates and their derivatives, the energies, the first num_reactions
!     rows of the composition rate vector, and the per-reaction scales
!     of the settings). Each slot belongs to exactly one row of network's
!     reaction table that computes its rate (network_table_errors checks
!     it), so a new reaction takes a new slot: raise num_reactions;
!   - the two branching rows of the composition rate vector after the
!     reactions (rv_*);
!   - the extra cross-section scales (xs_*): the rates outside the
!     reaction table, pep, Be7 electron and proton capture, and hep;
!   - the settings record net_settings (net_settings_t), the only
!     configuration net/ reads. The star layer builds it from the
!     controls and installs it at the start of every run (network's
!     set_net_settings, which checks it first).
! It also maps the legacy cross-section vector sstandard onto the
! scales (cross_section_scales_from_sstandard).
module net_configuration
      use, intrinsic :: iso_fortran_env, only: dp => real64
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

! The legacy cross-section vector sstandard (17 ratios to Bahcall &
! Ulrich 1988; used when use_new_nuclear_rates is off): the reaction
! slot or the extra slot of each of its entries. Entries 9 (O16 + alpha)
! and 13 (C12 + C12) belong to reactions YREC never computed.
      integer, parameter, public :: num_sstandard = 17
      integer, parameter :: sstandard_reaction_slot(num_sstandard) = &
           [r_pp, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, r_po16, &
           r_ac13, 0, r_ac12, r_an14, r_3alpha, 0, 0, 0, 0, 0]
      integer, parameter :: sstandard_extra_slot(num_sstandard) = &
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
           xs_pep, xs_be7e, xs_be7p, xs_hep]

! Everything net/ needs from the run's configuration.
      type, public :: net_settings_t
! log10 T at or below which: (1) nothing burns; (4) helium burning is
! off; (5) the thermal neutrino losses are skipped. (2) and (3) gated
! the level-1 equilibrium burn, removed 2026-10; they are read and
! unused.
         real(dp) :: logT_cutoffs(5) = 0.0d0
! the weak-screening exponent below which screening is weak
         real(dp) :: weak_screening_threshold = 0.0d0
! multipliers of each reaction's S-factor (cross_section_scale) and of
! its first and second energy derivatives (s_prime_scale,
! s_double_prime_scale), by reaction slot; only the S-factor rate law
! reads them (set_net_settings warns about a scale set on any other)
         real(dp) :: cross_section_scale(num_reactions) = 1.0d0
         real(dp) :: s_prime_scale(num_reactions) = 1.0d0
         real(dp) :: s_double_prime_scale(num_reactions) = 1.0d0
! the multipliers of the rates outside the table (xs_* slots), and of
! the Be7 + p S-factor's derivatives
         real(dp) :: extra_cross_section_scale(num_extra_scales) = 1.0d0
         real(dp) :: be7p_s_prime_scale = 1.0d0, be7p_s_double_prime_scale = 1.0d0
! the multipliers of the N15(p,alpha)C12 and N15(p,gamma)O16 S-factors,
! which set the N15 + p branching
         real(dp) :: n15_pa_scale = 1.0d0, n15_pg_scale = 1.0d0
! Li6(p,a), Li7(p,a), Be9(p,g), Be9(p,d), Be9(p,a) rate multipliers
         real(dp) :: li6_rate_scale = 1.0d0, li7_rate_scale = 1.0d0, &
              be9_pg_rate_scale = 1.0d0, be9_pd_rate_scale = 1.0d0, &
              be9_palpha_rate_scale = 1.0d0
! compute the solar neutrino fluxes; use the Itoh et al. (1996) thermal
! neutrino losses (else Beaudet, Petrosian & Salpeter 1967)
         logical :: calc_neutrinos = .false., use_itoh_neutrino_loss = .false.
      end type net_settings_t
      type(net_settings_t), public, protected :: net_settings

      public :: store_net_settings, cross_section_scales_from_sstandard

contains

!----------------------------------------------------------------------
! store_net_settings: make settings the run's net_settings, unchecked
! (network's set_net_settings checks them first and is what the star
! layer calls).
subroutine store_net_settings(settings)
      type(net_settings_t), intent(in) :: settings
      net_settings = settings
end subroutine store_net_settings

!----------------------------------------------------------------------
! cross_section_scales_from_sstandard: the cross-section scales of the
! reactions and the extra scales from the legacy vector sstandard (for
! the controls read, when use_new_nuclear_rates is off). An entry that
! belongs to no reaction gets a warning if it is not 1: it has no
! effect.
subroutine cross_section_scales_from_sstandard(sstandard, &
     cross_section_scale, extra_cross_section_scale)
      use net_messages, only: net_report, message_length
      real(dp), intent(in) :: sstandard(num_sstandard)
      real(dp), intent(out) :: cross_section_scale(num_reactions), &
           extra_cross_section_scale(num_extra_scales)
      character(len=message_length) :: line(1)
      integer :: i

      cross_section_scale = 1.0d0
      extra_cross_section_scale = 1.0d0
      do i = 1, num_sstandard
         if (sstandard_reaction_slot(i) /= 0) then
            cross_section_scale(sstandard_reaction_slot(i)) = sstandard(i)
         else if (sstandard_extra_slot(i) /= 0) then
            extra_cross_section_scale(sstandard_extra_slot(i)) = sstandard(i)
         else if (sstandard(i) /= 1.0d0) then
            write(line(1), '(a,i0,a,1pe12.5,a)') ' WARNING: net: sstandard(', i, &
                 ') = ', sstandard(i), ' belongs to no reaction YREC computes; it has no effect'
            call net_report(line)
         end if
      end do
end subroutine cross_section_scales_from_sstandard

end module net_configuration
