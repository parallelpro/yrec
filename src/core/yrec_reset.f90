!----------------------------------------------------------------------
! yrec_reset_lib
!----------------------------------------------------------------------
! Added 2026 (phase five, step C -- ROADMAP.md): the re-entrancy
! prologue. Fresh-process semantics for repeated run_yrec calls in
! one process, by pristine snapshot: at the FIRST entry the four
! state roots' startup state is captured (module/static storage is
! zero-initialized plus declaration defaults -- exactly what a new
! process sees); at every LATER entry it is restored, the controls
! rewind to their declaration defaults (so parmin reads over pristine
! values), the step/diagnostics routines are told to re-initialize
! their SAVEd locals, and any file units run 1 left open are closed
! (closing a not-connected unit is a no-op). Physics-domain table
! state is deliberately NOT reset: tables reloaded by star_setup are
! overwritten, and the lazy first-use guards (readco and friends)
! keep already-loaded tables, which is correct while both calls use
! the same table configuration -- the acceptance test
! (test_reentry.py) verifies the end-to-end equivalence.
module yrec_reset_lib
      use star_info_lib, only: star, star_info
      use star_job_lib, only: job, star_job
      use evolve_state_lib, only: evo, evolve_state, &
           evolve_step_reset_pending, output_diag_reset_pending
      use controls_reset_lib
      implicit none

      logical, save :: first_entry = .true.
      type(star_info), save :: star0
      type(star_job), save :: job0
      type(evolve_state), save :: evo0

contains

subroutine yrec_run_prologue

      integer :: u

      if (first_entry) then
         star0 = star
         job0 = job
         evo0 = evo
         call controls_capture
         first_entry = .false.
      else
         star = star0
         job = job0
         evo = evo0
         call controls_restore
         evolve_step_reset_pending = .true.
         output_diag_reset_pending = .true.
         do u = 7, 99
            close(u)
         end do
      end if

      return
end subroutine yrec_run_prologue

end module yrec_reset_lib
