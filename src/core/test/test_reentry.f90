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
!   test_reentry A.nml1 A.nml2               -- run A twice (original form)
!   test_reentry A.nml1 A.nml2 B.nml1 B.nml2  -- run A, then B
!
! The four-argument form (2026, bugsweep sec-11 batch 0) is the one
! that actually exercises control carry-over: A sets a legacy-spelled
! NAMELIST control that B omits, and B's outputs must still match a
! fresh B. The paths go through the libyrec override slots so the
! engine's own getarg default is bypassed for BOTH calls.
program test_reentry
      use star_info_lib, only: control_nml_override, physics_nml_override
      implicit none
      integer :: ierr, nargs
      character(len=256) :: nml(4)

      nargs = command_argument_count()
      if (nargs /= 2 .and. nargs /= 4) then
         write(*,'(a)') 'usage: test_reentry A.nml1 A.nml2 [B.nml1 B.nml2]'
         stop 2
      end if
      call get_command_argument(1, nml(1))
      call get_command_argument(2, nml(2))
      if (nargs == 4) then
         call get_command_argument(3, nml(3))
         call get_command_argument(4, nml(4))
      else
         nml(3) = nml(1)
         nml(4) = nml(2)
      end if

      ierr = 0
      control_nml_override = nml(1)
      physics_nml_override = nml(2)
      call run_yrec(ierr)
      if (ierr /= 0) then
         write(*,'(a,i4)') 'test_reentry: FIRST run failed, ierr = ', ierr
         stop 1
      end if
      write(*,'(a)') 'test_reentry: first run done'
      control_nml_override = nml(3)
      physics_nml_override = nml(4)
      call run_yrec(ierr)
      if (ierr /= 0) then
         write(*,'(a,i4)') 'test_reentry: SECOND run failed, ierr = ', ierr
         stop 1
      end if
      write(*,'(a)') 'test_reentry: second run done'
      control_nml_override = ' '
      physics_nml_override = ' '
end program test_reentry
