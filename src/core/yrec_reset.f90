!----------------------------------------------------------------------
! yrec_reset_lib
!----------------------------------------------------------------------
! Added 2026 (phase five, step C -- ROADMAP.md): the re-entrancy
! prologue. Fresh-process semantics for repeated run_yrec calls in
! one process, by pristine snapshot: at the FIRST entry star's
! startup state is captured (module/static storage is
! zero-initialized plus declaration defaults -- exactly what a new
! process sees); at every LATER entry it is restored, the
! step/diagnostics routines are told to re-initialize their SAVEd
! locals, and any file units run 1 left open are closed (closing a
! not-connected unit is a no-op). The controls no longer need a
! snapshot here: read_controls itself resets star%ctrl to defaults
! and re-seeds the namelist buffer before every read (2026 phase B). Physics-domain table
! state is deliberately NOT reset: tables reloaded by star_setup are
! overwritten, and the lazy first-use guards (readco and friends)
! keep already-loaded tables, which is correct while both calls use
! the same table configuration -- the acceptance test
! (test_reentry.py) verifies the end-to-end equivalence.
module yrec_reset_lib
      use star_info_lib, only: star, star_info, evolve_step_reset_pending, &
           observables_reset_pending
      use rotation_scratch_lib
      use point_scratch_lib
      implicit none

      logical, save :: first_entry = .true.
! 2026 solver-scratch cleanup: the rotation workspace lives outside
! star_info now (rotation_scratch_lib), so it needs its own pristine
! snapshots alongside star0.
      type(rotation_diffusion_state), save :: rot_scr0
      type(point_physics_scratch), save :: pt_scr0
      type(mdphy_state), save :: mix_scr0
      type(circulation_velocity_state), save :: circ_scr0
      ! star0 also carries star%job and star%evo since the 2026
      ! MESA-convention fold -- one snapshot covers all three roots.
      type(star_info), save :: star0

contains

subroutine yrec_run_prologue

      integer :: u

! 2026 phase B (controls->star% campaign): the controls_capture/
! controls_restore snapshot pair is gone. The namelist BUFFER
! (controls_lib) no longer needs restoring here because
! read_controls re-seeds it from pristine star%ctrl defaults before
! every read; phys_const_lib needs no restoring because, with cmixl
! evicted to star%, every member is recomputed by setups each run.
      if (first_entry) then
         star0 = star
         rot_scr0 = rot_scr
         pt_scr0 = pt_scr
         mix_scr0 = mix_scr
         circ_scr0 = circ_scr
         first_entry = .false.
      else
         star = star0
         rot_scr = rot_scr0
         pt_scr = pt_scr0
         mix_scr = mix_scr0
         circ_scr = circ_scr0
         evolve_step_reset_pending = .true.
         observables_reset_pending = .true.
         do u = 7, 99
            close(u)
         end do
      end if

      return
end subroutine yrec_run_prologue

end module yrec_reset_lib
