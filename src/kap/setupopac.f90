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
! later by getopac.
subroutine setupopac(envelope_hydrogen_fraction, laol_work_array, &
     alex06_table_path, kurucz_table_path, kurucz_table2_path, &
     laol_table_path, laol_table2_path, opal95_table_path, &
     opal92_table_path, opal92_table2_path, pure_z_table_path, &
     alex95_table_paths)

      implicit none

      double precision, intent(in) :: envelope_hydrogen_fraction
      double precision, intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: alex06_table_path, &
           kurucz_table_path, kurucz_table2_path, laol_table_path, &
           laol_table2_path, opal95_table_path, opal92_table_path, &
           opal92_table2_path, pure_z_table_path
      character(len=256), intent(in) :: alex95_table_paths(7)

! --- opacity common blocks - modified 3/09 ---
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

! common/zramp/: only use_z_ramp is used here; remaining members are
! placeholders preserving the shared storage layout (see getopac.f90
! for why renaming here doesn't require touching other files).
      double precision :: rsclzc(50), rsclzm1(50), rsclzm2(50)
      integer :: iolaol2, ioopal2, nk
      logical :: use_z_ramp
      common/zramp/ rsclzc, rsclzm1, rsclzm2, iolaol2, ioopal2, nk, &
           use_z_ramp

! common/gravs3/: only use_diffusion_z is used here.
      double precision :: fgry, fgrz
      logical :: lthoul, use_diffusion_z
      common/gravs3/ fgry, fgrz, lthoul, use_diffusion_z

! common/nwlaol/: not read here, but written by the table-loading
! routines called below (rdzlaol/zsulaol etc.) via this same block.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez

      save

!     THIS SUBROUTINE READS IN SPECIFIED OPACITY TABLES AND
!     SET UP SPLINES FOR THE TABLES.
!     WHEN LZRAMP=T OR LDIFZ=T THEN READ IN SECOND SET OF
!     OPACITY TABLES AT DIFFERENT Z (E.G. ZOPAL952).
      use_two_z_tables = use_z_ramp .or. use_diffusion_z

!     INTERIOR TABLES

!     READ IN OPAL95 TABLES
      if (use_opal95_tables) then
         call ll95tbl(opal95_table_path)
         call op95xtab(envelope_hydrogen_fraction)
      end if

!     READ IN OPAL92 TABLES AT ZOPAL1 AND ZOPAL2
      if (use_opal92_tables) then
         call setllo(opal92_table_path, opal92_table2_path)
         call ll4th(envelope_hydrogen_fraction)
      end if
!     READ IN LAOL89 TABLES AT ZLAOL1 AND ZLAOL2
      if (use_laol89_tables) then
         call rdlaol(laol_work_array, laol_table_path, laol_table2_path)
         call sulaol
      end if

!     READ IN LAOL89 PURE Z TABLE

      if (use_pure_z_table) then
         call rdzlaol(pure_z_table_path)
         call zsulaol
      end if

!     LOW TEMP TABLES

!     READ IN ALEX 2006 TABLES
      if (use_alex06_tables) then
         call readalex06(alex06_table_path)
!     READ IN ALEX 1995 TABLES
      else if (use_alex95_tables) then
         call alxtbl(alex95_table_paths)
         call alx8th(envelope_hydrogen_fraction)
!     READ IN KURUCZ TABLE AT ZKUR1 AND ZKUR2
      else if (use_kurucz90_tables) then
         call setkrz(kurucz_table_path, kurucz_table2_path)
      end if
      return
end subroutine setupopac
