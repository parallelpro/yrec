!----------------------------------------------------------------------
! t6rinteos01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original t6rinteos01.f; only variable names, source form, and
! comment style were updated.
!
! OPAL 2001 EOS analogue of t6rinterp.f90 (see there for the general
! description). Interpolates the already X-interpolated table slice
! (v%x_interp_result) in T6 and density to produce v%esact.
!
! Readability W3 (2026): the interpolation itself is t6rint_core.f90,
! shared with t6rinteos06; this wrapper keeps the 2001 out-of-range
! message.
subroutine t6rinteos01(v, slr, slt, ierr)

      use opal_eos_lib
      use luout_lib
      implicit none

      type(opal_eos_vintage), intent(inout) :: v

      double precision, intent(in) :: slr, slt

      integer, intent(out) :: ierr

      ierr = 0

      call t6rint_core(v, slr, slt)
      if (v%esact.gt.1.0d+15) then
         write(run_log_unit,'("Interpolation indices out of range", &
              &";please report conditions.")')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if

      return
end subroutine t6rinteos01
