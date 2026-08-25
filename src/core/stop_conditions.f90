!----------------------------------------------------------------------
! stop_conditions
!----------------------------------------------------------------------
! New (2026, core/ readability phase 2). The single home for the
! driver-level stop machinery that was previously duplicated between
! run_yrec.f90 and evolve_step.f90:
!
!  * the step_status protocol constants (evolve_step's contract with
!    run_yrec's model loop);
!  * the end-age predicates -- the "within 1 year of the target age"
!    convention existed as four hand-written copies in evolve_step;
!  * the central-abundance stop table (D / X / Y): the trigger check
!    (evolve_step, once per converged model) and the start-of-run
!    disarm pass (run_yrec) are two walks over the same three rows.
!    Keeping them as one table closes the door on the copy-paste bug
!    fixed in phase 1, where the disarm branch for hydrogen negated
!    the DEUTERIUM threshold and printed the deuterium values (a
!    defect inherited from the original F77).
!
! The three threshold arrays stay separate namelist-bound variables
! in controls_lib (central_deuterium_stop / central_hydrogen_stop /
! central_helium_stop); stop_value/set_stop_value select by row. A
! stop is configured when its threshold is positive; disarming
! negates it (the historical convention -- every check tests .gt.0).
module stop_conditions
      use star_info_lib, only: star, i_h1, i_h2, i_he4
      use const_lib
      use luout_lib
      implicit none
      private
      public :: step_continue, step_kind_card_done, step_leave_run_loop
      public :: reached_end_age, approaching_end_age
      public :: abundance_stop_triggered, disarm_satisfied_stops

! evolve_step -> run_yrec model-loop protocol
      integer, parameter :: step_continue = 0        ! advance accepted
      integer, parameter :: step_kind_card_done = 1  ! age/abundance stop
      integer, parameter :: step_leave_run_loop = 2  ! target radius crossed

! the abundance-stop table: central xa slot and message letter per row
      integer, parameter :: nstops = 3
      integer, parameter :: stop_species(nstops) = [i_h2, i_h1, i_he4]
      character(len=1), parameter :: stop_letter(nstops) = ['D','X','Y']

contains

! ---------------------------------------------------------------
! The current age is at (within 1 year of) the kind card's target
! end age. Used after a step: the kind card is done.
logical function reached_end_age(nk)
      integer, intent(in) :: nk
      reached_end_age = star%job%end_age_stop_active(nk) .and. &
           star%job%target_end_age(nk).gt.0.0d0 .and. &
           (star%job%target_end_age(nk)-star%dage*1.0d9).le.1.0d0
end function reached_end_age

! ---------------------------------------------------------------
! The NEXT step will land on the target end age (same 1-year
! tolerance). Used before a step to pre-arm the final model's pulse
! and sound-speed output.
logical function approaching_end_age(nk)
      integer, intent(in) :: nk
      approaching_end_age = star%job%end_age_stop_active(nk) .and. &
           star%job%target_end_age(nk).gt.0.0d0 .and. &
           abs(star%job%target_end_age(nk)-star%dage*1.0d9-star%timestep_yr) &
           .le.1.0d0
end function approaching_end_age

! ---------------------------------------------------------------
! Post-step trigger check: has any configured central-abundance stop
! been crossed? Prints the historical CENTRAL <el> message per hit.
logical function abundance_stop_triggered(nk)
      integer, intent(in) :: nk
      integer :: k
      abundance_stop_triggered = .false.
      do k = 1, nstops
         if (star%job%end_age_stop_active(nk) .and. stop_value(k,nk).gt.0.0d0 .and. &
             star%xa(stop_species(k),1).lt.stop_value(k,nk)) then
            write(*,'(A,E12.4,A,E12.4)') 'CENTRAL '//stop_letter(k)//' ', &
                 star%xa(stop_species(k),1), ' BELOW STOP VALUE ', &
                 stop_value(k,nk)
            abundance_stop_triggered = .true.
         end if
      end do
end function abundance_stop_triggered

! ---------------------------------------------------------------
! Start-of-run pass (run_yrec, after starin): any stop whose target
! the starting model already satisfies is disarmed (threshold
! negated), with the historical STARTING <el> ... STOP DISABLED
! message to the terminal and the short/log stream. The message
! reports the disarmed (negative) threshold, as the original did.
subroutine disarm_satisfied_stops(nk)
      integer, intent(in) :: nk
      integer :: k
      if (.not. star%job%end_age_stop_active(nk)) return
      do k = 1, nstops
         if (stop_value(k,nk).gt.0.0d0 .and. &
             star%xa(stop_species(k),1).lt.stop_value(k,nk)) then
            call set_stop_value(k, nk, -stop_value(k,nk))
            write(*,'(A,E12.4,A,E12.4,A)') 'STARTING '//stop_letter(k)//' ', &
                 star%xa(stop_species(k),1), ' BELOW STOP VALUE ', &
                 stop_value(k,nk), ' STOP DISABLED.'
            write(short_file_unit,'(A,E12.4,A,E12.4,A)') &
                 'STARTING '//stop_letter(k)//' ', &
                 star%xa(stop_species(k),1), ' BELOW STOP VALUE ', &
                 stop_value(k,nk), ' STOP DISABLED.'
         end if
      end do
end subroutine disarm_satisfied_stops

! ---------------------------------------------------------------
double precision function stop_value(k, nk)
      integer, intent(in) :: k, nk
      select case (k)
      case (1);      stop_value = star%job%central_deuterium_stop(nk)
      case (2);      stop_value = star%job%central_hydrogen_stop(nk)
      case default;  stop_value = star%job%central_helium_stop(nk)
      end select
end function stop_value

subroutine set_stop_value(k, nk, val)
      integer, intent(in) :: k, nk
      double precision, intent(in) :: val
      select case (k)
      case (1);      star%job%central_deuterium_stop(nk) = val
      case (2);      star%job%central_hydrogen_stop(nk) = val
      case default;  star%job%central_helium_stop(nk) = val
      end select
end subroutine set_stop_value

end module stop_conditions
