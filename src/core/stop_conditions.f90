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
! Protocol (2026 single-entry redesign): run_yrec calls
! init_stop_conditions once per kind card, right after the starting
! model is read; evolve_step calls check_stop_conditions once per
! converged model, after the model is written. New stop options
! (2026: the log_L/Teff/log_g/nu_max structure limits below are the
! first) belong as further checks inside check_stop_conditions plus
! their registry rows -- not as new public entries.
!
! The three threshold arrays stay separate namelist-bound variables
! in controls_lib (central_deuterium_stop / central_hydrogen_stop /
! central_helium_stop); stop_value/set_stop_value select by row. A
! stop is configured when its threshold is positive; disarming
! negates it (the historical convention -- every check tests .gt.0).
module stop_conditions
      use star_info_lib, only: star, i_h1, i_h2, i_he4
      use luout_lib
      use math_lib
      implicit none
      private
      public :: step_continue, step_kind_card_done, step_leave_run_loop
      public :: reached_end_age
      public :: check_stop_conditions, init_stop_conditions

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
! Per-model stop check, called by evolve_step after the converged
! model has been written. In order: the end-age stop, the configured
! central-abundance stops (both end the kind card), then -- when a
! star calibration is running -- the target-radius check, whose hit
! leaves the model loop so run_yrec can rescale and retry (the
! chkscal protocol: iteration 1 only primes previous-model state, so
! the calibration comparison starts at iteration 2).
subroutine check_stop_conditions(model_iteration, step_status)
      integer, intent(in) :: model_iteration
      integer, intent(out) :: step_status

      step_status = step_continue

      if (reached_end_age(star%job%nk)) then
         step_status = step_kind_card_done
         return
      end if

      if (abundance_stop_triggered(star%job%nk)) then
         step_status = step_kind_card_done
         return
      end if

! Structure limits (2026): log_L / Teff / log g / nu_max
! against the current kind card's *_upper_limit(nk) /
! *_lower_limit(nk) values -- per-card arrays like the other
! stopping criteria (a limit at its +-1d99 sentinel is off). Not
! checked on pure-rescale relax cards, where the structure is
! mid-relaxation.
      if (star%job%rescale_kind(star%job%nk) /= 2) then
         if (structure_limit_stop_triggered()) then
            step_status = step_kind_card_done
            return
         end if
      end if

! TEST IF MODEL IS NEAR DESIRED Teff AND L. IF NOT RESCALE AND TRY AGAIN.
      if (star%ctrl%calibrate_star_flag .and. .not. star%star_found_flag) then
         if (mod(star%job%nk,2).eq.0 .and. model_iteration.ne.1) then
            call check_star_calibration(star%log_L, star%log_Teff, &
                 star%dage, star%job%nk)
            if (star%just_passed_target_radius_flag) then
               step_status = step_leave_run_loop
            end if
         end if
      end if
end subroutine check_stop_conditions

! ---------------------------------------------------------------
! Structure-limit stops: each configured limit is checked
! against the freshly computed observables (compute_observables runs
! before this, so log_g_surface and the scaling-relation nu_max are
! current). Prints one STOP line per hit to the terminal and the
! run log.
logical function structure_limit_stop_triggered()
      integer :: nk
      integer, parameter :: nlim = 4
      character(len=7), parameter :: qname(nlim) = &
           ['log_L  ', 'Teff   ', 'log_g  ', 'nu_max ']
      double precision :: qval(nlim), qup(nlim), qlo(nlim)
      integer :: k

      nk = star%job%nk
      qval = [star%log_L, exp10(star%log_Teff), star%log_g_surface, &
              star%nu_max]
      qup = [star%job%log_L_upper_limit(nk), star%job%Teff_upper_limit(nk), &
             star%job%log_g_upper_limit(nk), star%job%nu_max_upper_limit(nk)]
      qlo = [star%job%log_L_lower_limit(nk), star%job%Teff_lower_limit(nk), &
             star%job%log_g_lower_limit(nk), star%job%nu_max_lower_limit(nk)]

      structure_limit_stop_triggered = .false.
      do k = 1, nlim
         if (qup(k) < 0.9d99 .and. qval(k) > qup(k)) then
            write(*,10) trim(qname(k)), qval(k), 'above', &
                 trim(qname(k))//'_upper_limit', qup(k)
            write(run_log_unit,10) trim(qname(k)), qval(k), 'above', &
                 trim(qname(k))//'_upper_limit', qup(k)
            structure_limit_stop_triggered = .true.
         else if (qlo(k) > -0.9d99 .and. qval(k) < qlo(k)) then
            write(*,10) trim(qname(k)), qval(k), 'below', &
                 trim(qname(k))//'_lower_limit', qlo(k)
            write(run_log_unit,10) trim(qname(k)), qval(k), 'below', &
                 trim(qname(k))//'_lower_limit', qlo(k)
            structure_limit_stop_triggered = .true.
         end if
      end do
   10 format(1x,'STOP: ',a,' =',es12.4,1x,a,1x,a,' =',es12.4)
end function structure_limit_stop_triggered

! ---------------------------------------------------------------
! The current age is at (within 1 year of) the kind card's target
! end age. Also used directly by evolve_step BEFORE the model is
! written, to keep the final model's saved timestep from being
! zeroed (that use runs earlier than check_stop_conditions can).
logical function reached_end_age(nk)
      integer, intent(in) :: nk
      reached_end_age = star%job%end_age_stop_active(nk) .and. &
           star%job%target_end_age(nk).gt.0.0d0 .and. &
           (star%job%target_end_age(nk)-star%dage*1.0d9).le.1.0d0
end function reached_end_age

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
! Start-of-card setup (run_yrec, after the starting model is read):
! any stop whose target the starting model already satisfies is
! disarmed (threshold negated) so the run doesn't stop at model 1,
! with the historical STARTING <el> ... STOP DISABLED message to the
! terminal and the run log. The message reports the disarmed
! (negative) threshold, as the original did.
subroutine init_stop_conditions(nk)
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
            write(run_log_unit,'(A,E12.4,A,E12.4,A)') &
                 'STARTING '//stop_letter(k)//' ', &
                 star%xa(stop_species(k),1), ' BELOW STOP VALUE ', &
                 stop_value(k,nk), ' STOP DISABLED.'
         end if
      end do
end subroutine init_stop_conditions

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
