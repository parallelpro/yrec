!----------------------------------------------------------------------
! yrec_reset
!----------------------------------------------------------------------
! The re-entrancy
! prologue. Fresh-process semantics for repeated run_yrec calls in one
! process. yrec_run_prologue gives run_yrec a star in its process-start
! state: a copy of star_template, a star that nothing writes (static
! storage, so zero-initialized plus the type's declaration defaults --
! exactly what a module-level star was at process start). At every entry
! after the first it also closes the units the last run left open
! (close_run_units: the run's own files, not a host program's units). The controls need
! no snapshot: read_controls resets s%control and s%physics to their
! defaults and re-seeds
! the namelist buffer before every read. Physics-domain table state is
! deliberately NOT reset: tables reloaded by prepare_run are
! overwritten, and the lazy first-use guards (readco and friends) keep
! already-loaded tables, which is correct while both calls use the same
! table configuration -- the acceptance test (test_reentry.py) verifies
! the end-to-end equivalence.
module yrec_reset
      use star_info_lib, only: star_info
      implicit none
      private
      public :: yrec_run_prologue, close_run_units

      logical, save :: first_entry = .true.
! The process-start star: never written, only copied (see the header).
! (The rotation and element-diffusion workspaces need nothing here: their
! drivers allocate them for each step.)
      type(star_info), save :: star_template

contains

! Gives run_yrec its star: s becomes a copy of star_template, after
! closing the units an earlier run left open (close_run_units: the run's own, not a host's). It writes
! the whole record, which the contract's member lists do not show.
!   star reads: none
!   star writes: none
subroutine yrec_run_prologue(s)
      type(star_info), allocatable, intent(inout) :: s

      if (first_entry) then
         first_entry = .false.
      else
         call close_run_units
      end if
      if (.not. allocated(s)) allocate(s)
      s = star_template

      return
end subroutine yrec_run_prologue

! ---------------------------------------------------------------
! The units the run connects, closed if they are open, so that a host
! that reads the run's files after yrec_run finds them complete (M312):
! the two connected by number, the run log (run_messages' run_log_unit)
! and the final model (run_files' last_model_unit); and the files
! open_run_files opens with newunit, the starting model, the Monte
! Carlo samples, summary and model files and the isochrone file, whose
! unit variables are set back to -1 (not connected) as run_files' own
! close does before a later run opens them: a stale number could be the
! one a later newunit hands out. A unit a host program opened stays
! open, and so do the standard units, which the run's units are before
! the controls read assigns them.
subroutine close_run_units
      use run_messages, only: run_log_unit
      use run_files, only: last_model_unit, starting_model_unit, dynamics_unit, &
           mc_summary_unit, mc_model_unit, isochrone_unit
      call close_run_unit(run_log_unit)
      call close_run_unit(last_model_unit)
      call close_newunit_file(starting_model_unit)
      call close_newunit_file(dynamics_unit)
      call close_newunit_file(mc_summary_unit)
      call close_newunit_file(mc_model_unit)
      call close_newunit_file(isochrone_unit)
contains
subroutine close_newunit_file(unit)
      integer, intent(inout) :: unit
      if (unit == -1) return
      call close_run_unit(unit)
      unit = -1
end subroutine close_newunit_file
subroutine close_run_unit(unit)
      use, intrinsic :: iso_fortran_env, only: error_unit, input_unit, output_unit
      integer, intent(in) :: unit
      logical :: is_open
      integer :: ios
      if (unit == input_unit .or. unit == output_unit .or. unit == error_unit) return
      inquire(unit=unit, opened=is_open, iostat=ios)
      if (ios == 0 .and. is_open) close(unit)
end subroutine close_run_unit
end subroutine close_run_units

! ---------------------------------------------------------------
! Sources
!   The re-entrancy prologue came with the embeddable engine (2026, phase
!   five, step C; docs/ROADMAP.md).
end module yrec_reset
