!----------------------------------------------------------------------
! yrec_reset
!----------------------------------------------------------------------
! Added 2026 (phase five, step C -- docs/ROADMAP.md): the re-entrancy
! prologue. Fresh-process semantics for repeated run_yrec calls in one
! process. yrec_run_prologue gives run_yrec a star in its process-start
! state: a copy of star_template, a star that nothing writes (static
! storage, so zero-initialized plus the type's declaration defaults --
! exactly what a module-level star was at process start). At every entry
! after the first it also tells the step/diagnostics routines to
! re-initialize their SAVEd locals and closes any file units the last run
! left open (closing a not-connected unit is a no-op). The controls need
! no snapshot: read_controls resets s%ctrl to its defaults and re-seeds
! the namelist buffer before every read. Physics-domain table state is
! deliberately NOT reset: tables reloaded by prepare_run are
! overwritten, and the lazy first-use guards (readco and friends) keep
! already-loaded tables, which is correct while both calls use the same
! table configuration -- the acceptance test (test_reentry.py) verifies
! the end-to-end equivalence.
module yrec_reset
      use star_info_lib, only: star_info, evolve_step_reset_pending
      implicit none
      private
      public :: yrec_run_prologue

      logical, save :: first_entry = .true.
! The process-start star: never written, only copied (see the header).
! (The rotation and element-diffusion workspaces need nothing here: their
! drivers allocate them for each step.)
      type(star_info), save :: star_template

contains

subroutine yrec_run_prologue(s)
      type(star_info), allocatable, intent(inout) :: s
      integer :: u

      if (first_entry) then
         first_entry = .false.
      else
         evolve_step_reset_pending = .true.
         do u = 7, 99
            close(u)
         end do
      end if
      if (.not. allocated(s)) allocate(s)
      s = star_template

      return
end subroutine yrec_run_prologue

end module yrec_reset
