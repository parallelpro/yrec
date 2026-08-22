!----------------------------------------------------------------------
! test_reentry
!----------------------------------------------------------------------
! Phase-five step C acceptance driver (2026 -- ROADMAP.md): calls the
! run_yrec engine TWICE in one process with the same command-line
! namelists. The second run's outputs must be byte-identical to a
! fresh single-process run's -- the re-entrancy contract ("two calls
! in one process == two processes"), enforced by test_reentry.py,
! which runs this binary and compares against a fresh yrec run.
program test_reentry
      implicit none
      integer :: ierr

      ierr = 0
      call run_yrec(ierr)
      if (ierr /= 0) then
         write(*,'(a,i4)') 'test_reentry: FIRST run failed, ierr = ', ierr
         stop 1
      end if
      write(*,'(a)') 'test_reentry: first run done'
      call run_yrec(ierr)
      if (ierr /= 0) then
         write(*,'(a,i4)') 'test_reentry: SECOND run failed, ierr = ', ierr
         stop 1
      end if
      write(*,'(a)') 'test_reentry: second run done'
end program test_reentry
