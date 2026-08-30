!----------------------------------------------------------------------
! surfopac
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original surfopac.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Sets up the low-temperature/surface opacity tables (OPAL95, OPAL92,
! ALEX95) for the envelope hydrogen fraction, assuming the underlying
! tables have already been read in. Relocated here from atm/ (2026,
! atm/ phase-two reorg): zero atm-domain content -- it only refreshes
! cached table slices in kap/opal95/, kap/opal92/, kap/alex94/ -- and
! its callers (core/read_starting_model.f90, setup/rezone.f90) aren't in atm/
! either; same "misplaced, meval.f90-style" pattern as alsurfp.f90's
! earlier move the other direction (kap/ -> atm/).
subroutine surfopac(hydrogen_fraction, ierr)
      use star_info_lib, only: star
      implicit none

      double precision, intent(in) :: hydrogen_fraction
      integer, intent(out) :: ierr

      ierr = 0
!     THIS SUBROUTINE SETS UP SURF OPACITY TABLES
!     ASSUMES TABLES HAVE ALREADY BEEN READ IN

!     INTERIOR TABLES

!     SETUP OPAL95 TABLES
      if (star%ctrl%use_opal95_tables) then
      call opal95_surface_table(hydrogen_fraction)
      end if
!     SETUP IN OPAL92 TABLES AT ZOPAL1 AND ZOPAL2
      if (star%ctrl%use_opal92_tables) then
        call opal92_surface_table(hydrogen_fraction, ierr)
        if (ierr /= 0) return
      end if

!     LOW TEMP TABLES

!     INTERPOLATE ALEX95 TABLES
      if (star%ctrl%use_alex95_tables) then
       call alex94_surface_table(hydrogen_fraction)
      end if


      return
end subroutine surfopac
