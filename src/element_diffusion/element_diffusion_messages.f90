!----------------------------------------------------------------------
! element_diffusion_messages
!----------------------------------------------------------------------
! Element diffusion's messages, in one place: every text the domain
! writes, with where it goes. The kernels print nothing; the two
! pipelines' drivers (microdiff, gravitational_settling) and
! diffuse_elements call these.
!   note_settling_region  the run log: a model that cannot settle (fully
!                         convective, or X or Y below its floor
!                         everywhere), once per suspension through the
!                         star's latch, every model under solver
!                         diagnostics; and "SETTLING RESUMED" when a
!                         suspension ends;
!   note_solve_pass       the run log under solver diagnostics: each pass
!                         of the implicit solve, its largest change and
!                         where;
!   report_unconverged    the terminal (run_messages' terminal_unit) and
!                         the run log: an implicit solve that missed its
!                         tolerance in its passes (the run goes on);
!   report_numerics_failure
!                         standard output only: the numerics routine
!                         whose failure diffuse_elements returns as ierr,
!                         and the text of its status.
! The texts and destinations are the F77's (the B phase rewords them:
! audit/module-evaluation-7-2026-11-10/B_LISTS.md, element_diffusion 6).
module element_diffusion_messages
      use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
      use run_messages, only: run_log_unit, terminal_unit, solver_diagnostics
      use numerics_lib, only: numerics_status_text
      use settling_region, only: region_found, region_fully_convective, &
           region_no_hydrogen, region_no_helium
      implicit none
      private
      public :: note_settling_region, note_solve_pass, report_unconverged, &
           report_numerics_failure

      integer, parameter :: message_length = 160
      character(len=*), parameter :: unwritten_line = achar(0)

contains

! The settling region's verdict in the run log. A model that cannot
! settle is noted when the latch suspension_reported is clear (or always,
! under solver diagnostics) and sets it; a model that settles after a
! suspension notes "SETTLING RESUMED" and clears it. hydrogen_floor and
! helium_min are the floors the note quotes.
subroutine note_settling_region(region, hydrogen_floor, helium_min, suspension_reported)
      integer, intent(in) :: region
      real(dp), intent(in) :: hydrogen_floor, helium_min
      logical, intent(inout) :: suspension_reported
      character(len=message_length) :: lines(1)

      lines = unwritten_line
      if (region == region_found) then
         if (suspension_reported) then
            write(lines,'(1x,'' SETTLING RESUMED'')')
            call write_diffusion_lines(run_log_unit, lines)
            suspension_reported = .false.
         end if
         return
      end if
      if (.not. (solver_diagnostics() .or. .not. suspension_reported)) return
      select case (region)
      case (region_fully_convective)
         write(lines,'(1x,'' FULLY CONVECTIVE MODEL - NO SETTLING'')')
      case (region_no_hydrogen)
         write(lines,'(1x,''X BELOW '',f9.6,'' IN WHOLE MODEL-NO SETTLING'')') hydrogen_floor
      case (region_no_helium)
         write(lines,'(1x,''Y BELOW '',f9.6,'' IN WHOLE MODEL-NO SETTLING'')') helium_min
      end select
      call write_diffusion_lines(run_log_unit, lines)
      suspension_reported = .true.
end subroutine note_settling_region

! One pass of an implicit solve (pass number, its largest change of the
! abundance and the equally spaced shell where it is), under solver
! diagnostics.
subroutine note_solve_pass(pass, change, shell)
      integer, intent(in) :: pass, shell
      real(dp), intent(in) :: change
      character(len=message_length) :: lines(1)

      if (.not. solver_diagnostics()) return
      lines = unwritten_line
      write(lines,'(1x,''ITERATION '',i3,'' DXMAX '',1pe10.2,'' IMAX '',i4)') pass, change, shell
      call write_diffusion_lines(run_log_unit, lines)
end subroutine note_solve_pass

! An implicit solve that missed its tolerance in num_passes passes:
! routine is the solve's name (MICRODIFF, GRSETT), quantity the
! abundance it solves for (D, X), change and shell the last pass's.
subroutine report_unconverged(routine, quantity, tolerance, num_passes, change, shell)
      character(len=*), intent(in) :: routine, quantity
      real(dp), intent(in) :: tolerance, change
      integer, intent(in) :: num_passes, shell
      character(len=message_length) :: lines(2)

      lines = unwritten_line
      write(lines,'(1x,a,'' FAILED TO CONVERGE TO WITHIN '',1pe9.3,'' IN '',i3, &
           &''ITERATIONS''/1x,''LAST ITERATION CHANGE IN '',a,'' '',1pe9.3, &
           &'' IN EQUALLY SPACED SHELL '',i5)') routine, tolerance, num_passes, &
           quantity, change, shell
      call write_diffusion_lines(terminal_unit, lines)
      call write_diffusion_lines(run_log_unit, lines)
end subroutine report_unconverged

! ' <routine>: <the text of status>' on standard output, as the
! list-directed write it replaces printed it.
subroutine report_numerics_failure(routine, status)
      character(len=*), intent(in) :: routine
      integer, intent(in) :: status
      character(len=message_length) :: lines(1)

      lines(1) = ' '//routine//': '//numerics_status_text(status)
      call write_diffusion_lines(output_unit, lines)
end subroutine report_numerics_failure

subroutine write_diffusion_lines(unit, lines)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: lines(:)
      integer :: k

      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         write(unit,'(a)') trim(lines(k))
      end do
end subroutine write_diffusion_lines

end module element_diffusion_messages
