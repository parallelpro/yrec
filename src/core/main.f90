!----------------------------------------------------------------------
! main
!----------------------------------------------------------------------
! The command-line entry point. The whole job -- parameter and table
! setup, the Monte-Carlo / kind-card run loops, model evolution,
! output and calibration -- lives in core/run_yrec.f90 (the
! embeddable engine, also reachable through core/yrec_capi.f90);
! this program only runs it and turns its error code into an exit
! status.
!
! The original main.f had no PROGRAM statement (an unnamed F77 main
! program); it is named "main" here, matching the file name.
program main
      implicit none
      integer :: ierr

      ierr = 0
      call run_yrec(ierr)
! negative ierr is numerics_termination (numerics_lib) -- the
! historical "solution diverged" stop, which exited 0; preserve that
! as a clean exit. Positive ierr is a real error (exit 1).
      if (ierr > 0) stop 1
end program main
