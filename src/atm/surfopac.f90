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
! tables have already been read in.
subroutine surfopac(hydrogen_fraction)
      implicit none

      double precision, intent(in) :: hydrogen_fraction

! --- opacity common blocks - modified 3/09 ---
! common/newopac/: only use_opal95_tables/use_opal92_tables/
! use_alex95_tables are used here. Naming matches getopac.f90.
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max
      logical :: use_alex06_tables, use_laol89_tables, use_opal92_tables, &
           use_opal95_tables, use_kurucz90_tables, use_alex95_tables, &
           use_two_z_tables
      common /newopac/ laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max, use_alex06_tables, &
           use_laol89_tables, use_opal92_tables, use_opal95_tables, &
           use_kurucz90_tables, use_alex95_tables, use_two_z_tables

      save

!     THIS SUBROUTINE SETS UP SURF OPACITY TABLES
!     ASSUMES TABLES HAVE ALREADY BEEN READ IN

!     INTERIOR TABLES

!     SETUP OPAL95 TABLES
      if (use_opal95_tables) then
      call op95xtab(hydrogen_fraction)
      end if
!     SETUP IN OPAL92 TABLES AT ZOPAL1 AND ZOPAL2
      if (use_opal92_tables) then
        call ll4th(hydrogen_fraction)
      end if

!     LOW TEMP TABLES

!     INTERPOLATE ALEX95 TABLES
      if (use_alex95_tables) then
       call alx8th(hydrogen_fraction)
      end if


      return
end subroutine surfopac
