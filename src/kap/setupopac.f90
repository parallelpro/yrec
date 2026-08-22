!----------------------------------------------------------------------
! setupopac
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setupopac.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Reads in the opacity tables selected by the use_*_tables flags (set
! from the run's namelist) and builds the interpolation splines used
! later by kap_lib.f90's kap_get.
subroutine setupopac(envelope_hydrogen_fraction, laol_work_array, &
     alex06_table_path, kurucz_table_path, kurucz_table2_path, &
     laol_table_path, laol_table2_path, opal95_table_path, &
     opal92_table_path, opal92_table2_path, pure_z_table_path, &
     alex95_table_paths, ierr)

      use const_lib
      implicit none

      double precision, intent(in) :: envelope_hydrogen_fraction
      double precision, intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: alex06_table_path, &
           kurucz_table_path, kurucz_table2_path, laol_table_path, &
           laol_table2_path, opal95_table_path, opal92_table_path, &
           opal92_table2_path, pure_z_table_path
      character(len=256), intent(in) :: alex95_table_paths(7)
      integer, intent(out) :: ierr





      save

!     THIS SUBROUTINE READS IN SPECIFIED OPACITY TABLES AND
!     SET UP SPLINES FOR THE TABLES.
!     WHEN LZRAMP=T OR LDIFZ=T THEN READ IN SECOND SET OF
!     OPACITY TABLES AT DIFFERENT Z (E.G. ZOPAL952).
      use_two_z_tables = use_z_ramp .or. use_diffusion_z

!     INTERIOR TABLES

!     READ IN OPAL95 TABLES
      ierr = 0
      if (use_opal95_tables) then
         call ll95tbl(opal95_table_path, ierr)
         if (ierr /= 0) return
         call op95xtab(envelope_hydrogen_fraction)
      end if

!     READ IN OPAL92 TABLES AT ZOPAL1 AND ZOPAL2
      if (use_opal92_tables) then
         call setllo(opal92_table_path, opal92_table2_path)
         call ll4th(envelope_hydrogen_fraction)
      end if
!     READ IN LAOL89 TABLES AT ZLAOL1 AND ZLAOL2
      if (use_laol89_tables) then
         call rdlaol(laol_work_array, laol_table_path, laol_table2_path, ierr)
         if (ierr /= 0) return
         call sulaol
      end if

!     READ IN LAOL89 PURE Z TABLE

      if (use_pure_z_table) then
         call rdzlaol(pure_z_table_path, ierr)
         if (ierr /= 0) return
         call zsulaol
      end if

!     LOW TEMP TABLES

!     READ IN ALEX 2006 TABLES
      if (use_alex06_tables) then
         call readalex06(alex06_table_path, ierr)
         if (ierr /= 0) return
!     READ IN ALEX 1995 TABLES
      else if (use_alex95_tables) then
         call alxtbl(alex95_table_paths, ierr)
         if (ierr /= 0) return
         call alx8th(envelope_hydrogen_fraction)
!     READ IN KURUCZ TABLE AT ZKUR1 AND ZKUR2
      else if (use_kurucz90_tables) then
         call setkrz(kurucz_table_path, kurucz_table2_path)
      end if
      return
end subroutine setupopac
