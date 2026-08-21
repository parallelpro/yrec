!----------------------------------------------------------------------
! setllo
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setllo.f; only variable names, source form, and comment
! style were updated. FORMAT strings and the READ/OPEN statements
! that must byte-for-byte match the on-disk table file are preserved
! verbatim.
!
! DBG 5/94 Modified to read in second opacity table at different Z.
! Reads the OPAL92 (Lawrence Livermore) opacity table(s) and builds
! the spline interpolation coefficients used later by ll4th.f90 and
! the OPAL92 lookup routines (yllo3d/yllo3d2, not part of this batch).
subroutine setllo(opal92_table_path, opal92_table2_path)

      use const_lib
      implicit none
      integer, parameter :: num_t = 50
      integer, parameter :: num_d = 17
      integer, parameter :: num_x = 3
      integer, parameter :: num_xt = num_t*num_x
! CONST is declared but not referenced anywhere in the original
! setllo.f body; preserved here as an unused placeholder.
      double precision, parameter :: const_unused = 11604.5d0

! MHP 8/25 removed variables not used in subroutine
      character(len=256), intent(in) :: opal92_table_path, opal92_table2_path


! MHP 8/25 removed common block with file names
!     COMMON/LUFNM/ FLAST, FFIRST, FRUN, FSTAND, FFERMI,
!     1    FDEBUG, FTRACK, FSHORT, FMILNE, FMODPT,
!     2    FSTOR, FPMOD, FPENV, FPATM, FDYN,
!     3    FLLDAT, FSNU, FSCOMP, FKUR,
!     4    FMHD1, FMHD2, FMHD3, FMHD4, FMHD5, FMHD6, FMHD7, FMHD8
! GRID ENTRIES FOR TEMPERATURE, AND ABUNDANCE (X)
      double precision :: opal92_grid_logt(num_t), opal92_grid_x(num_x), &
           opal92_grid_logr(num_d)
      common /gllot/ opal92_grid_logt, opal92_grid_x, opal92_grid_logr
! LL OPACITY
      double precision :: opal92_log10_opacity(num_xt, num_d)
      integer :: opal92_num_x, opal92_num_temps
      common /llot/ opal92_log10_opacity, opal92_num_x, opal92_num_temps
! DBG 5/94 for different Z
      double precision :: opal92_grid_logt_z2(num_t), opal92_grid_x_z2(num_x), &
           opal92_grid_logr_z2(num_d)
      common /gllot2/ opal92_grid_logt_z2, opal92_grid_x_z2, opal92_grid_logr_z2
      double precision :: opal92_log10_opacity_z2(num_xt, num_d)
      integer :: opal92_num_x_z2, opal92_num_temps_z2
      common /llot2/ opal92_log10_opacity_z2, opal92_num_x_z2, opal92_num_temps_z2

      double precision :: fgry, fgrz
      logical :: lthoul, use_diffusion_z
      common/gravs3/fgry, fgrz, lthoul, use_diffusion_z

! MHP 8/25 Removed character file names from common block
! common/zramp/: not used here; declared only to preserve the shared
! storage layout (see getopac.f90/setupopac.f90 for these names).
      double precision :: rsclzc(50), rsclzm1(50), rsclzm2(50)
      integer :: iolaol2, ioopal2, nk
      logical :: use_z_ramp
      common/zramp/ rsclzc, rsclzm1, rsclzm2, iolaol2, ioopal2, nk, &
           use_z_ramp

! common/newopac/: only use_two_z_tables is used here.
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

      double precision :: local_grid_y(num_x), local_grid_z(num_x)
      integer :: i, k, density_index, num_temps_read
      double precision :: grid_temp_k

!     OPEN TABLE
      open(unit=laol_table_unit,file=opal92_table_path)
      do 10 i=1,num_x
!        READ GRID POINT FOR ABUNDANCE
!        READ NUMBER OF GRIDS FOR DENSITY, AND TEMPERATURE
        read(laol_table_unit,190,end=97) opal92_grid_x(i), local_grid_z(i)
        local_grid_y(i)=1.0d0-opal92_grid_x(i)-local_grid_z(i)
  190   format(33x,f7.4,2x,f7.4)
         read(laol_table_unit,'()')
!        READ  LOG(DENSITY/TEMPERATURE**3)
            read(laol_table_unit, 200) (opal92_grid_logr(density_index), density_index=1, num_d)
  200   format (6x, 17f7.1)
!        READ GRID VALUES FOR TEMPERATURE, AND OPACITY TABLE
         do 20 k=1, num_t
         read(laol_table_unit,196,end=93) grid_temp_k, &
              (opal92_log10_opacity(k+(i-1)*num_t,density_index),density_index=1,num_d)
         opal92_grid_logt(k)=dlog10(grid_temp_k)
   20    continue
   93    num_temps_read=k-1
  196    format(18f7.3)
!
   10 continue
!     CLOSE THE TABLE WE HAVE READ
   97 close(laol_table_unit,err=99)
      opal92_num_temps = num_temps_read
      opal92_num_x=i-1

! DBG 5/94 Second Opacity Table read here
      if (use_two_z_tables) then
         open(unit=ioopal2,file=opal92_table2_path)
         do 510 i=1,num_x
            read(ioopal2,190,end=597) opal92_grid_x_z2(i), local_grid_z(i)
            local_grid_y(i)=1.0d0-opal92_grid_x_z2(i)-local_grid_z(i)
            read(ioopal2,'()')
            read(ioopal2, 200) (opal92_grid_logr_z2(density_index), density_index=1, num_d)
            do 520 k=1, num_t
               read(ioopal2,196,end=593) grid_temp_k, &
                    (opal92_log10_opacity_z2(k+(i-1)*num_t,density_index),density_index=1,num_d)
               opal92_grid_logt_z2(k)=dlog10(grid_temp_k)
  520       continue
  593       num_temps_read=k-1
  510    continue
  597    close(ioopal2,err=99)
         opal92_num_temps_z2 = num_temps_read
         opal92_num_x_z2=i-1
      end if
!
      call ylloc

      return
   99 stop 'ERROR IN FILE CLOSING'
end subroutine setllo
