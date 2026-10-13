!----------------------------------------------------------------------
! stop_conditions
!----------------------------------------------------------------------
! New (2026, core/ readability phase 2). The single home for the
! driver-level stop machinery that was previously duplicated between
! run_driver (run_yrec) and step_driver (evolve_step):
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
! in controls_buffer (central_deuterium_stop / central_hydrogen_stop /
! central_helium_stop); stop_value/set_stop_value select by row. A
! stop is configured when its threshold is positive; disarming
! negates it (the historical convention -- every check tests .gt.0).
module stop_conditions
      use calibration, only: check_star_calibration
      use star_info_lib, only: i_h1, i_h2, i_he4, star_info
      use io_units
      use math_lib
      implicit none
      private
      public :: step_continue, step_kind_card_done, step_leave_run_loop
      public :: reached_end_age
      public :: check_stop_conditions, init_stop_conditions
      public :: converged_model_is_nan
      public :: check_rotation_initialised

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
! Configuration check shared by run_yrec (after the starting model is
! read) and evolve_step (after a diverged-model reload): a rotating
! run whose model carries no rotation rate is an error (ierr = 1).
subroutine check_rotation_initialised(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(inout) :: ierr
      if ((s%omega(1) .eq. 0) .and. (s%job%rotation_active)) then
1611      format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                 ' Initialize rotation rates or set LROT to', &
                 ' FALSE.')
          print 1611
          ! 2026 (phase five, step B): configuration error returns to the
          ! CLI wrapper (which stops) instead of stopping here.
          ierr = 1
      endif
end subroutine check_rotation_initialised

! ---------------------------------------------------------------
! 2026 (bugsweep Batch 3): NaN guard. Nothing in the solver tests
! for NaN, and a NaN structure passes every convergence test (all
! comparisons are false), so a run that went non-finite (the
! zahb->tahb dt = 0 case before the heburn fix) kept "converging",
! wrote NaN models to the end of its budget and exited 0. Called by
! evolve_step on the converged structure, before it is written;
! a hit is a real error (positive ierr -> exit 1), not a stop.
logical function converged_model_is_nan(s)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      type(star_info), intent(inout) :: s
      integer :: nz, k
      character(len=24) :: what

      converged_model_is_nan = .true.
      nz = s%nz
      what = ' '
      if (ieee_is_nan(s%log_L)) then
         what = 'log_L'
      else if (ieee_is_nan(s%log_Teff)) then
         what = 'log_Teff'
      else if (ieee_is_nan(s%age_gyr) .or. ieee_is_nan(s%dt)) then
         what = 'age/timestep'
      else
         do k = 1, nz
            if (ieee_is_nan(s%logT(k)) .or. ieee_is_nan(s%logRho(k)) &
                 .or. ieee_is_nan(s%logP(k)) .or. ieee_is_nan(s%logR(k)) &
                 .or. ieee_is_nan(s%luminosity_lsun(k))) then
               write(what,'(a,i0)') 'structure at zone ', k
               exit
            end if
         end do
         if (what == ' ') converged_model_is_nan = .false.
      end if
      if (converged_model_is_nan) then
         write(*,10) trim(what), s%model_number + 1
         write(run_log_unit,10) trim(what), s%model_number + 1
         s%termination_reason = 'NaN in '//trim(what)
      end if
   10 format(1x,'ERROR: NaN in ',a,' after model ',i0,' converged; stopping')
end function converged_model_is_nan

! ---------------------------------------------------------------
! Per-model stop check, called by evolve_step after the converged
! model has been written. In order: the end-age stop, the configured
! central-abundance stops (both end the kind card), then -- when a
! star calibration is running -- the target-radius check, whose hit
! leaves the model loop so run_yrec can rescale and retry (the
! check_star_calibration protocol: iteration 1 only primes
! previous-model state, so the calibration comparison starts at
! iteration 2).
subroutine check_stop_conditions(s, model_iteration, step_status)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: model_iteration
      integer, intent(out) :: step_status

      step_status = step_continue

      if (reached_end_age(s, s%job%nk)) then
         write(s%termination_reason,'(a,es10.3,a)') &
              'target age reached (', &
              s%job%target_end_age(s%job%nk), ' yr)'
         step_status = step_kind_card_done
         return
      end if

      if (abundance_stop_triggered(s, s%job%nk)) then
         step_status = step_kind_card_done
         return
      end if

! Structure limits (2026): log_L / Teff / log g / nu_max
! against the current kind card's *_upper_limit(nk) /
! *_lower_limit(nk) values -- per-card arrays like the other
! stopping criteria (a limit at its +-1d99 sentinel is off). Not
! checked on pure-rescale relax cards, where the structure is
! mid-relaxation.
      if (s%job%rescale_kind(s%job%nk) /= 2) then
         if (structure_limit_stop_triggered(s)) then
            step_status = step_kind_card_done
            return
         end if
      end if

! TEST IF MODEL IS NEAR DESIRED Teff AND L. IF NOT RESCALE AND TRY AGAIN.
      if (s%ctrl%calibrate_star_flag .and. .not. s%star_found_flag) then
         if (mod(s%job%nk,2).eq.0 .and. model_iteration.ne.1) then
            call check_star_calibration(s, s%log_L, s%log_Teff, &
                 s%age_gyr, s%job%nk)
            if (s%just_passed_target_radius_flag) then
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
logical function structure_limit_stop_triggered(s)
      type(star_info), intent(inout) :: s
      integer :: nk
      integer, parameter :: nlim = 4
! a limit at its +-1d99 sentinel is off: the test below is against
! this slightly smaller magnitude
      double precision, parameter :: limit_active_below = 0.9d99
      character(len=7), parameter :: qname(nlim) = &
           ['log_L  ', 'Teff   ', 'log_g  ', 'nu_max ']
      double precision :: qval(nlim), qup(nlim), qlo(nlim)
      integer :: k

      nk = s%job%nk
      qval = [s%log_L, exp10(s%log_Teff), s%log_g_surface, &
              s%nu_max]
      qup = [s%job%log_L_upper_limit(nk), s%job%Teff_upper_limit(nk), &
             s%job%log_g_upper_limit(nk), s%job%nu_max_upper_limit(nk)]
      qlo = [s%job%log_L_lower_limit(nk), s%job%Teff_lower_limit(nk), &
             s%job%log_g_lower_limit(nk), s%job%nu_max_lower_limit(nk)]

      structure_limit_stop_triggered = .false.
      do k = 1, nlim
         if (qup(k) < limit_active_below .and. qval(k) > qup(k)) then
            write(*,10) trim(qname(k)), qval(k), 'above', &
                 trim(qname(k))//'_upper_limit', qup(k)
            write(run_log_unit,10) trim(qname(k)), qval(k), 'above', &
                 trim(qname(k))//'_upper_limit', qup(k)
            s%termination_reason = &
                 trim(qname(k))//' above '//trim(qname(k))//'_upper_limit'
            structure_limit_stop_triggered = .true.
         else if (qlo(k) > -limit_active_below .and. qval(k) < qlo(k)) then
            write(*,10) trim(qname(k)), qval(k), 'below', &
                 trim(qname(k))//'_lower_limit', qlo(k)
            write(run_log_unit,10) trim(qname(k)), qval(k), 'below', &
                 trim(qname(k))//'_lower_limit', qlo(k)
            s%termination_reason = &
                 trim(qname(k))//' below '//trim(qname(k))//'_lower_limit'
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
logical function reached_end_age(s, nk)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: nk
      reached_end_age = s%job%end_age_stop_active(nk) .and. &
           s%job%target_end_age(nk).gt.0.0d0 .and. &
           (s%job%target_end_age(nk)-s%age_gyr*1.0d9).le.1.0d0
end function reached_end_age

! ---------------------------------------------------------------
! Post-step trigger check: has any configured central-abundance stop
! been crossed? Prints the historical CENTRAL <el> message per hit.
logical function abundance_stop_triggered(s, nk)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: nk
      integer :: k
      abundance_stop_triggered = .false.
      do k = 1, nstops
         if (s%job%end_age_stop_active(nk) .and. stop_value(s, k,nk).gt.0.0d0 .and. &
             s%xa(stop_species(k),1).lt.stop_value(s, k,nk)) then
            write(*,'(A,E12.4,A,E12.4)') 'CENTRAL '//stop_letter(k)//' ', &
                 s%xa(stop_species(k),1), ' BELOW STOP VALUE ', &
                 stop_value(s, k,nk)
            s%termination_reason = &
                 'central '//stop_letter(k)//' below stop value'
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
subroutine init_stop_conditions(s, nk)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: nk
      integer :: k
! default end-of-card reason; overridden by whichever stop fires
! (read by run_log's end-of-run summary)
      s%termination_reason = 'model budget exhausted'
      if (.not. s%job%end_age_stop_active(nk)) return
      do k = 1, nstops
         if (stop_value(s, k,nk).gt.0.0d0 .and. &
             s%xa(stop_species(k),1).lt.stop_value(s, k,nk)) then
            call set_stop_value(s, k, nk, -stop_value(s, k,nk))
            write(*,'(A,E12.4,A,E12.4,A)') 'STARTING '//stop_letter(k)//' ', &
                 s%xa(stop_species(k),1), ' BELOW STOP VALUE ', &
                 stop_value(s, k,nk), ' STOP DISABLED.'
            write(run_log_unit,'(A,E12.4,A,E12.4,A)') &
                 'STARTING '//stop_letter(k)//' ', &
                 s%xa(stop_species(k),1), ' BELOW STOP VALUE ', &
                 stop_value(s, k,nk), ' STOP DISABLED.'
         end if
      end do
end subroutine init_stop_conditions

! ---------------------------------------------------------------
double precision function stop_value(s, k, nk)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: k, nk
      select case (k)
      case (1);      stop_value = s%job%central_deuterium_stop(nk)
      case (2);      stop_value = s%job%central_hydrogen_stop(nk)
      case default;  stop_value = s%job%central_helium_stop(nk)
      end select
end function stop_value

subroutine set_stop_value(s, k, nk, val)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: k, nk
      double precision, intent(in) :: val
      select case (k)
      case (1);      s%job%central_deuterium_stop(nk) = val
      case (2);      s%job%central_hydrogen_stop(nk) = val
      case default;  s%job%central_helium_stop(nk) = val
      end select
end subroutine set_stop_value

end module stop_conditions
