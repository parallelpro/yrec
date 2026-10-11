!----------------------------------------------------------------------
! io_units
!----------------------------------------------------------------------
! The three fixed units written from everywhere, assigned at each read
! (run_files' assign_unit_numbers). Newer readers and writers use
! open(newunit=...) instead.
!   run_log_unit    -- the run log (log_output_file, default run.log)
!   terminal_unit   -- stdout (always 6; kept as a named unit only
!                      because ~20 files write through it)
!   last_model_unit -- the final .mod model (last_model_file)
module io_units
      implicit none

      integer :: last_model_unit, run_log_unit, terminal_unit

end module io_units
