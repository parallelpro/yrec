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
      use star_info_lib, only: star, star_info, evolve_step_reset_pending, &
           observables_reset_pending
      use controls_reset_lib, only: controls_capture, controls_restore
      implicit none

      logical, save :: first_entry = .true.
      ! star0 also carries star%job and star%evo since the 2026
      ! MESA-convention fold -- one snapshot covers all three roots.
      type(star_info), save :: star0

contains

subroutine yrec_run_prologue

      integer :: u

      if (first_entry) then
         star0 = star
         call controls_capture
         first_entry = .false.
      else
         star = star0
         call controls_restore
         evolve_step_reset_pending = .true.
         observables_reset_pending = .true.
         do u = 7, 99
            close(u)
         end do
      end if

      return
end subroutine yrec_run_prologue

end module yrec_reset_lib
