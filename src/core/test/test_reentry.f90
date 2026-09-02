!----------------------------------------------------------------------
! test_reentry
!----------------------------------------------------------------------
! Phase-five step C acceptance driver (2026 -- ROADMAP.md): calls the
! run_yrec engine TWICE in one process. The second run's outputs must
! be byte-identical to a fresh single-process run's -- the re-entrancy
! contract ("two calls in one process == two processes"), enforced by
! test_reentry.py, which runs this binary and compares against a
! fresh yrec run.
!
! Usage:
!   test_reentry A.inlist            -- run A twice (original form)
!   test_reentry A.inlist B.inlist   -- run A, then B
!
! The two-inlist form (2026, bugsweep sec-11 batch 0) is the one that
! actually exercises control carry-over: A sets a control that B
! omits, and B's outputs must still match a fresh B. The paths go
! through the libyrec override slots (single-inlist style: the same
! file fills both slots) so the engine's own getarg default is
! bypassed for BOTH calls.
program test_reentry
      use star_info_lib, only: control_nml_override, physics_nml_override
      implicit none
      integer :: ierr, nargs
      character(len=256) :: inlist(2)

      nargs = command_argument_count()
      if (nargs < 1 .or. nargs > 2) then
         write(*,'(a)') 'usage: test_reentry A.inlist [B.inlist]'
         stop 2
      end if
      call get_command_argument(1, inlist(1))
      if (nargs == 2) then
         call get_command_argument(2, inlist(2))
      else
         inlist(2) = inlist(1)
      end if

      ierr = 0
      control_nml_override = inlist(1)
      physics_nml_override = inlist(1)
      call run_yrec(ierr)
      if (ierr /= 0) then
         write(*,'(a,i4)') 'test_reentry: FIRST run failed, ierr = ', ierr
         stop 1
      end if
      write(*,'(a)') 'test_reentry: first run done'
      control_nml_override = inlist(2)
      physics_nml_override = inlist(2)
      call run_yrec(ierr)
      if (ierr /= 0) then
         write(*,'(a,i4)') 'test_reentry: SECOND run failed, ierr = ', ierr
         stop 1
      end if
      write(*,'(a)') 'test_reentry: second run done'
      control_nml_override = ' '
      physics_nml_override = ' '
end program test_reentry
