!----------------------------------------------------------------------
! run_cards
!----------------------------------------------------------------------
! The run plan: each kind card (evolve, rescale, or rescale and evolve)
! turned into its stop flags and rescale targets, and echoed to the
! terminal and the run log.
module run_cards
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use io_units, only: terminal_unit, log_line
      use controls_buffer, only: num_runs, run_kind, num_models, target_end_age, &
           timestep_override, central_deuterium_stop, central_hydrogen_stop, &
           central_helium_stop, rescale_mass, rescale_x, rescale_z, rescale_core_mass, &
           rescale_kind, rescale_params, first_call_flag, &
           end_age_stop_active, timestep_override_active
      implicit none
      private
      public :: interpret_kind_cards
      public :: card_evolve, card_rescale, card_rescale_and_evolve

! The kinds of run card (run_kind, legacy KINDRN): evolve the model;
! rescale it (relaxed num_models times at zero timestep); or rescale
! and evolve it (pre-main sequence, which has no zero-timestep model).
      integer, parameter :: card_evolve = 1, card_rescale = 2, &
           card_rescale_and_evolve = 3

contains

! ---------------------------------------------------------------
! Every card's stop flags and rescale targets, and its line of the run
! plan. (The envelope Z a rescale card asks for is prepare_run's.)
subroutine interpret_kind_cards
      integer :: n

      call log_line('')
      call log_line(' run plan')
      call log_line('')

! card 1 has no previous card's model: it starts from the starting
! model whatever the deck says (as the F77 PARMIN forced LFIRST(1))
      first_call_flag(1) = .true.
      do n = 1, num_runs
         rescale_kind(n) = run_kind(n)
         if (evolves(n)) call set_stop_flags(n)
         if (rescales(n)) rescale_params(:,n) = [rescale_mass(n), rescale_x(n), &
              rescale_z(n), rescale_core_mass(n)]
         call echo_card(n)
      end do
end subroutine interpret_kind_cards

! Does card n evolve the model (and so stop at its stops)?
logical function evolves(n)
      integer, intent(in) :: n
      evolves = run_kind(n) == card_evolve .or. run_kind(n) == card_rescale_and_evolve
end function evolves

! Does card n rescale the model?
logical function rescales(n)
      integer, intent(in) :: n
      rescales = run_kind(n) == card_rescale .or. run_kind(n) == card_rescale_and_evolve
end function rescales

! ---------------------------------------------------------------
! An evolving card stops at its age and central-abundance stops (any
! positive one arms them), and may fix the timestep.
subroutine set_stop_flags(n)
      integer, intent(in) :: n
      end_age_stop_active(n) = target_end_age(n) > 0.0_dp .or. &
           central_deuterium_stop(n) > 0.0_dp .or. central_hydrogen_stop(n) > 0.0_dp &
           .or. central_helium_stop(n) > 0.0_dp
      timestep_override_active(n) = timestep_override(n) > 0.0_dp
end subroutine set_stop_flags

! ---------------------------------------------------------------
! Card n's lines of the run plan, to the terminal and the run log: what
! it does and from which model, its rescale targets, and its stops. An
! evolve or rescale card with no models prints nothing.
subroutine echo_card(n)
      integer, intent(in) :: n
      character(len=200) :: line
      integer :: i

      select case (run_kind(n))
      case (card_evolve)
         if (num_models(n) <= 0) return
         if (first_call_flag(n)) then
            write(line,350) n, num_models(n)
         else
            write(line,351) n, num_models(n)
         end if
      case (card_rescale)
         if (num_models(n) <= 0) return
         if (first_call_flag(n)) then
            write(line,450) n
         else
            write(line,451) n
         end if
      case (card_rescale_and_evolve)
         if (first_call_flag(n)) then
            write(line,550) n
         else
            write(line,551) n
         end if
      case default
         return
      end select
      call echo(line)
      if (rescales(n)) then
         write(line,452) num_models(n), (rescale_params(i,n), i = 1,4)
         call echo(line)
      end if
      if (evolves(n) .and. (end_age_stop_active(n) .or. timestep_override_active(n))) then
         write(line,370) target_end_age(n), timestep_override(n), central_deuterium_stop(n), central_hydrogen_stop(n), central_helium_stop(n)
         call echo(line)
      end if

  350 format(1x,'card',i3,': evolve up to',i6,' models from the starting model')
  351 format(1x,'card',i3,': evolve up to',i6,' models from the previous card''s model')
  450 format(1x,'card',i3,': rescale the starting model')
  451 format(1x,'card',i3,': rescale the previous card''s model')
  550 format(1x,'card',i3,': rescale and evolve the starting model')
  551 format(1x,'card',i3,': rescale and evolve the previous card''s model')
  452 format(9x,'relax for',i3,' models; targets (0 or negative = keep current):', &
           '  M =',f9.6,'  X =',f9.6,'  Z =',f9.6,'  core mass =',f10.6)
  370 format(9x,'stop conditions (0 = unused): age =',es9.2,' yr   fixed dt =',es9.2, &
           ' yr   central: D =',es10.3,'  X =',es10.3,'  Y =',es10.3)
end subroutine echo_card

! One line of the run plan, to the terminal and the run log.
subroutine echo(line)
      character(len=*), intent(in) :: line
      write(terminal_unit,'(a)') trim(line)
      call log_line(trim(line))
end subroutine echo

end module run_cards
