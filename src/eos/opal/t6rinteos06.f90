!----------------------------------------------------------------------
! t6rinteos06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original t6rinteos06.f; only variable names, source form, and
! comment style were updated.
!
! OPAL 2006 EOS analogue of t6rinteos01.f90 (see there and
! t6rinterp.f90 for the general description). Declared RECURSIVE in
! the original; preserved verbatim even though nothing here actually
! recurses.
!
! Readability W3 (2026): the interpolation itself is t6rint_core.f90,
! shared with t6rinteos01; this wrapper keeps the 2006 out-of-range
! message.
recursive subroutine t6rinteos06(v, slr, slt, ierr)

      use opal_eos_lib
      use luout_lib
      implicit none

      type(opal_eos_vintage), intent(inout) :: v

      double precision, intent(in) :: slr, slt

      integer, intent(out) :: ierr

      ierr = 0

      call t6rint_core(v, slr, slt)
      if (v%esact.gt.1.0d+15) then
         write(run_log_unit,'("T6RINTEOS06: Interpolation indices out", &
              &" of range;please report conditions.")')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if

      return
end subroutine t6rinteos06
