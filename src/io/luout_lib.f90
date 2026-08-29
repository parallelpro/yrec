!----------------------------------------------------------------------
! luout_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/luout/: I/O logical unit numbers,
! assigned once at run startup (io/read_controls.f90) and read
! broadly -- global configuration, not per-call data, so (like
! const_lib) this is a module of plain module-level variables rather
! than subroutine arguments.
!
! This module is the residue of the F77 fixed-unit convention. The
! MESA-style writers (history/profile/pulse, io/yrec_output.f90) use
! open(newunit=...) and never appear here; streams migrate out of
! this module as they modernize (2026: the second LAOL/OPAL92
! opacity-table units became newunit locals of their readers).
!
! 2026 descriptive-rename pass (formerly short_file_unit/iowr/
! ilast; idebug retired with the .debug stream):
!   run_log_unit    -- the run log (log_output_file, default run.log)
!   terminal_unit   -- stdout (always 6; kept as a named unit only
!                      because ~20 files write through it)
!   last_model_unit -- the final .mod model (last_model_file)
module luout_lib
      implicit none

      integer :: last_model_unit, run_log_unit, terminal_unit

end module luout_lib
