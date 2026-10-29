!----------------------------------------------------------------------
! network_checks
!----------------------------------------------------------------------
! The reaction table's self-check and the settings install that runs it:
!   set_net_settings      check the table, warn about scales without
!                         effect, and install the run's net_settings
!   network_table_errors  the number of wrong rows and slots, with
!                         table_row_problem and slot_problem saying what
!                         is wrong
!   inert_scale_count     the rows whose scale is set but unread
!   write_network_table   the table, with each row's Q against its
!                         masses, for test_net and documentation
! The tables themselves are network's.
module network_checks
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use net_configuration, only: num_reactions, net_settings_t, &
           store_net_settings
      use screening, only: max_screening_charge
      use reaction_rates, only: shares_rate, s_factor_fit, num_rate_laws
      use network, only: reactions, n_rows, n_solved, n_network, n_chains, &
           chain_name, mass_number, mass_excess_mev, neutron_mass_excess_mev, &
           q_rounding_mev, unbranched, be7_proton_capture, hydrogen_burning, &
           helium_burning, max_change, reaction_t, branch_fraction
      implicit none
      private
      public :: set_net_settings, inert_scale_count, network_table_errors, &
           write_network_table

contains

!----------------------------------------------------------------------
! set_net_settings: check the reaction table and install the run's
! settings (net_configuration's net_settings_t). The star layer calls it
! at the start of every Monte-Carlo run, with the record it builds from
! the controls (setup/net_settings_setup). A wrong table
! (network_table_errors) is reported row by row and slot by slot and
! returns ierr = 1, without installing. A
! reaction scale set to anything but 1 on a reaction whose rate law
! takes none (only the S-factor law reads the per-reaction scales) gets
! a warning: it has no effect.
subroutine set_net_settings(settings, ierr)
      use net_messages, only: net_report, message_length
      type(net_settings_t), intent(in) :: settings
      integer, intent(out) :: ierr
      character(len=message_length) :: line(1)
      integer :: j, k

      ierr = 0
      do j = 1, n_rows
         if (len_trim(table_row_problem(j)) == 0) cycle
         write(line(1), '(4a)') ' ERROR: net: reaction table row "', &
              trim(reactions(j)%name), '": ', trim(table_row_problem(j))
         call net_report(line)
         ierr = 1
      end do
      do k = 1, num_reactions
         if (len_trim(slot_problem(k)) == 0) cycle
         write(line(1), '(a,i0,2a)') ' ERROR: net: reaction slot ', k, ': ', &
              trim(slot_problem(k))
         call net_report(line)
         ierr = 1
      end do
      if (ierr /= 0) return
      do j = 1, n_rows
         if (.not. scale_has_no_effect(settings, j)) cycle
         write(line(1), '(4a)') ' WARNING: net: a cross-section scale of ', &
              trim(reactions(j)%name), ' is not 1, but its rate law', &
              ' takes no scale; it has no effect'
         call net_report(line)
      end do
      call store_net_settings(settings)
end subroutine set_net_settings

! Row j computes its rate with a law that takes no scale (only the
! S-factor law reads the per-reaction scales), yet one of its scales in
! settings is not 1.
pure logical function scale_has_no_effect(settings, j)
      type(net_settings_t), intent(in) :: settings
      integer, intent(in) :: j
      integer :: k

      scale_has_no_effect = .false.
      if (reactions(j)%rate_law == shares_rate .or. &
           reactions(j)%rate_law == s_factor_fit) return
      k = reactions(j)%slot
      scale_has_no_effect = settings%cross_section_scale(k) /= 1.0d0 .or. &
           settings%s_prime_scale(k) /= 1.0d0 .or. &
           settings%s_double_prime_scale(k) /= 1.0d0
end function scale_has_no_effect

! The number of rows scale_has_no_effect reports.
pure integer function inert_scale_count(settings)
      type(net_settings_t), intent(in) :: settings
      integer :: j

      inert_scale_count = 0
      do j = 1, n_rows
         if (scale_has_no_effect(settings, j)) inert_scale_count = inert_scale_count + 1
      end do
end function inert_scale_count

! ---------------------------------------------------------------
! The table's self-check. A row is wrong when it names a reaction slot,
! a chain, a branch, a rate law, a burning group or a species that does
! not exist; shares the rate of a slot no row computes; has a reactant
! outside the implicit solve (an explicit product cannot feed back);
! is an S-factor reaction without a fit or its reference S0; computes a
! rate with screening charges outside 1..max_screening_charge; does not
! conserve baryon number; or deposits more energy than its atomic masses
! release (beyond the q_rounding_mev the Q values are rounded to). A
! reaction slot is wrong unless exactly one row computes its rate (so no
! slot is dead, and none is written twice) and its rows' branch factors
! add up to 1 whatever the branching fractions (a branched reaction has
! both of its outcomes, and only those). network_table_errors counts the
! wrong rows and slots; table_row_problem and slot_problem say what is
! wrong; set_net_settings reports them.
pure integer function network_table_errors()
      integer :: j, k

      network_table_errors = 0
      do j = 1, n_rows
         if (len_trim(table_row_problem(j)) > 0) &
              network_table_errors = network_table_errors + 1
      end do
      do k = 1, num_reactions
         if (len_trim(slot_problem(k)) > 0) &
              network_table_errors = network_table_errors + 1
      end do
end function network_table_errors

! What is wrong with row j ('' if nothing).
pure function table_row_problem(j) result(problem)
      integer, intent(in) :: j
      character(len=64) :: problem
      type(reaction_t) :: row
      integer :: c

      problem = ''
      row = reactions(j)
      if (row%slot < 1 .or. row%slot > num_reactions) then
         problem = 'its rate slot is not a reaction slot'
      else if (row%chain < 1 .or. row%chain > n_chains) then
         problem = 'its chain does not exist'
      else if (row%branch < unbranched .or. row%branch > be7_proton_capture) then
         problem = 'its branch does not exist'
      else if (row%rate_law < shares_rate .or. row%rate_law > num_rate_laws) then
         problem = 'its rate law does not exist'
      else if (row%burns_in /= hydrogen_burning .and. &
           row%burns_in /= helium_burning) then
         problem = 'its burning group does not exist'
      else if (row%rate_law == shares_rate .and. &
           slot_owner_count(row%slot) == 0) then
         problem = 'it shares the rate of a slot no row computes'
      else if (any(row%reactants < 0) .or. any(row%reactants > n_solved)) then
         problem = 'a reactant is not a solved species'
      else if (row%rate_law == s_factor_fit .and. all(row%fit == 0.0d0)) then
         problem = 'it is an S-factor reaction without a fit'
      else if (row%rate_law == s_factor_fit .and. row%s0_kevb <= 0.0d0) then
         problem = 'it is an S-factor reaction without its reference S0'
      else if (row%rate_law /= shares_rate .and. &
           (min(row%screening%z1, row%screening%z2) < 1 .or. &
           max(row%screening%z1, row%screening%z2) > max_screening_charge)) then
         problem = 'its screening charges are not 1 to max_screening_charge'
      end if
      if (len_trim(problem) > 0) return
      do c = 1, 2*max_change - 1, 2
         if (row%change(c) < 0 .or. row%change(c) > n_network) &
              problem = 'its change list names a species that does not exist'
      end do
      if (len_trim(problem) > 0) return
      if (baryon_imbalance(j) /= 0) then
         problem = 'it does not conserve baryon number'
      else if (neutrino_energy(j) < -q_rounding_mev) then
         problem = 'it deposits more energy than its masses release'
      end if
end function table_row_problem

! What is wrong with reaction slot k ('' if nothing): not exactly one
! row computes its rate, or its rows' branch factors do not add up to 1
! (tried at two pairs of branching fractions that sum exactly).
pure function slot_problem(k) result(problem)
      integer, intent(in) :: k
      character(len=64) :: problem
      real(dp), parameter :: frac_n15_pa(2) = [0.25d0, 0.5d0], &
           frac_be7_ecap(2) = [0.75d0, 0.125d0]
      real(dp) :: total
      integer :: j, trial

      problem = ''
      if (slot_owner_count(k) /= 1) then
         write(problem, '(a,i0,a)') 'it has ', slot_owner_count(k), &
              ' table rows computing its rate (it needs one)'
         return
      end if
      do trial = 1, 2
         total = 0.0d0
         do j = 1, n_rows
            if (reactions(j)%slot == k) total = total + &
                 branch_fraction(j, frac_n15_pa(trial), frac_be7_ecap(trial))
         end do
         if (total /= 1.0d0) problem = 'its rows'' branch factors do not add up to 1'
      end do
end function slot_problem

! The number of rows that compute the rate of slot k.
pure integer function slot_owner_count(k)
      integer, intent(in) :: k
      integer :: j

      slot_owner_count = 0
      do j = 1, n_rows
         if (reactions(j)%slot == k .and. &
              reactions(j)%rate_law /= shares_rate) &
              slot_owner_count = slot_owner_count + 1
      end do
end function slot_owner_count

! The energy row j releases, MeV, from the atomic mass excesses of what
! it consumes and makes (neutrinos included; positron annihilation
! included, through the atomic masses).
pure real(dp) function q_from_masses(j)
      integer, intent(in) :: j
      integer :: c

      q_from_masses = -reactions(j)%emitted_mass*neutron_mass_excess_mev
      do c = 1, 2*max_change - 1, 2
         if (reactions(j)%change(c) /= 0) q_from_masses = q_from_masses - &
              reactions(j)%change(c+1)*mass_excess_mev(reactions(j)%change(c))
      end do
end function q_from_masses

! The mean neutrino energy per reaction the table implies, MeV: zero,
! to the 1 keV the table's Q values are rounded to, for a row with no
! beta+ decay or electron capture.
pure real(dp) function neutrino_energy(j)
      integer, intent(in) :: j
      neutrino_energy = q_from_masses(j) - reactions(j)%q_mev
end function neutrino_energy

pure integer function baryon_imbalance(j)
      integer, intent(in) :: j
      integer :: c

      baryon_imbalance = reactions(j)%emitted_mass
      do c = 1, 2*max_change - 1, 2
         if (reactions(j)%change(c) /= 0) baryon_imbalance = baryon_imbalance + &
              nint(mass_number(reactions(j)%change(c)))*reactions(j)%change(c+1)
      end do
end function baryon_imbalance

! ---------------------------------------------------------------
! The table, one row per reaction, for tests and documentation: each
! row's slot, branch, baryon imbalance, Q, chain, the Q its masses give
! and the neutrino energy that implies.
subroutine write_network_table(unit)
      integer, intent(in) :: unit
      character(len=len(chain_name)) :: chain
      integer :: j

      write(unit, '(a)') '# reaction                     slot  branch  baryon imbalance     Q (MeV)  chain' // &
           '           Q masses (MeV)  neutrinos (MeV)'
      do j = 1, n_rows
         chain = '?'
         if (reactions(j)%chain >= 1 .and. reactions(j)%chain <= n_chains) &
              chain = chain_name(reactions(j)%chain)
         write(unit, '(a24, 2x, i8, 2x, i6, 2x, i16, 2x, f10.6, 2x, a14, 2x, f14.6, 2x, f15.6)') &
              reactions(j)%name, reactions(j)%slot, reactions(j)%branch, &
              baryon_imbalance(j), reactions(j)%q_mev, chain, &
              q_from_masses(j), neutrino_energy(j)
      end do
      write(unit, '(a, i0)') '# table errors: ', network_table_errors()
end subroutine write_network_table

end module network_checks
