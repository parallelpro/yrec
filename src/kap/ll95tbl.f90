!----------------------------------------------------------------------
! ll95tbl
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ll95tbl.f; only variable names, source form, and comment
! style were updated. FORMAT strings and the READ/OPEN statements
! that must byte-for-byte match the on-disk table file are preserved
! verbatim.
!
! MHP 7/98 MODIFIED TO READ IN ALL METAL ABUNDACES FOR OPACITY TABLES.
! Reads the full set of OPAL95 opacity tables (all tabulated Z and X),
! skipping the header, then builds the fixed-Z table via op95ztab.f
! (not part of this batch).
subroutine ll95tbl(opal95_table_path)

      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      character(len=256), intent(in) :: opal95_table_path


! FULL SET OF TABLES: OPACITY AS A FUNCTION OF Z AND X, T, RHO/T6**3
! TABLES ARE INCREMENTED IN SETS OF NZ*NX.  SO THE TABLES FOR THE
! THIRD METAL ABUNDANCE (3 X 10**-4)BEGIN AT TABLE 21 AND END AT TABLE 30.
! FOR THE HIGH VALUES OF Z, THE NUMBER OF X TABLES IS NOT THE SAME (I.E.
! X<0.9 FOR Z=0.1).
! FOR EACH COMPOSITION A FULL GRID IN (T,RHO/T6**3) IS RETAINED.
      double precision :: opal95_grid_logt(num_t), opal95_grid_x(num_x), &
           opal95_grid_logr(num_d), opal95_grid_z(num_z), &
           opal95_full_opacity(num_xz,num_t,num_d)
      integer :: opal95_num_x_at_z(num_z), opal95_table_start_index(num_z)
      common /llot95a/ opal95_grid_logt, opal95_grid_x, opal95_grid_logr, &
           opal95_grid_z, opal95_full_opacity, opal95_num_x_at_z, &
           opal95_table_start_index
! INDICES FOR INTERPOLATION IN Z,X,T, AND R
      integer :: opal95_index_z, opal95_index_x(4,4), opal95_index_t, &
           opal95_index_rho(4)
      common/op95indx/ opal95_index_z, opal95_index_x, opal95_index_t, &
           opal95_index_rho

! common/newopac/: ll95tbl.f declares an OUT-OF-SYNC 16-member version
! of this block (10 doubles + 6 logicals) vs. the canonical 17-member
! layout (10 doubles + 7 logicals) used by getopac.f90/setupopac.f90
! and every other file in this batch -- a pre-existing bug in the
! original ll95tbl.f (this file is the ONLY one in the opacity family
! with this mismatched declaration; confirmed by comparing against
! setkrz.f, yalo3d.f, alxtbl.f, ll4th.f, setllo.f, sulaol.f, rdlaol.f,
! and rdzlaol.f, which all share the canonical layout). Only the 5th
! double (opal95_single_table_z / ZOPAL951) is actually referenced in
! this file's body, and it sits at the same byte offset in both
! layouts, so the mismatch is harmless in practice here -- but it is
! preserved exactly (not "fixed") per the conversion's COMMON-block
! rules (same order/type/dimensions as originally declared in THIS
! file). All other members below are unused placeholders local to
! this file's (buggy) view of /newopac/.
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, opal95_table_z2_placeholder, &
           alex_table_z1_placeholder, alex_table_z2_placeholder, &
           kurucz_table_z1_placeholder, kurucz_table_z2_placeholder
      logical :: use_laol89_tables_placeholder, use_opal92_tables_placeholder, &
           use_opal95_tables_placeholder, use_kurucz90_tables_placeholder, &
           use_alex95_tables_placeholder, use_two_z_tables_placeholder
      common /newopac/ laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, opal95_table_z2_placeholder, &
           alex_table_z1_placeholder, alex_table_z2_placeholder, &
           kurucz_table_z1_placeholder, kurucz_table_z2_placeholder, &
           use_laol89_tables_placeholder, use_opal92_tables_placeholder, &
           use_opal95_tables_placeholder, use_kurucz90_tables_placeholder, &
           use_alex95_tables_placeholder, use_two_z_tables_placeholder


! NUMBER OF COMPOSITION TABLES AT LOWER Z FOR EACH ABUNDANCE
      data opal95_table_start_index/0,10,20,30,40,50,60,70,80,90,100,109,118/
! NUMBER OF X TABLES AT EACH Z
      data opal95_num_x_at_z/10,10,10,10,10,10,10,10,10,10,9,9,8/
! TABULTED SET OF X
      data opal95_grid_x/0.0d0, 0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0,0.9d0, &
           0.95d0,1.0d0/
! TABULATED SET OF Z
      data opal95_grid_z/0.0d0, 0.0001d0, 0.0003d0, 0.001d0, 0.002d0, &
           0.004d0, 0.01d0, 0.02d0, 0.03d0, &
           0.04d0, 0.06d0, 0.08d0, 0.10d0/
! INDICES FOR INTERPOLATION
      data opal95_index_z,opal95_index_x,opal95_index_t,opal95_index_rho &
           /1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
! LOCAL VECTOR, USED TO SKIP LONG HEADER.
      character(len=132) :: header_line

      integer :: fmt_start, iz, ix, nn, n, i, j, nmax
      double precision :: xx, zz, xxt, target_z

!     READ IN OPACITY TABLES, SKIPPING HEADER.  HEY, IT WORKS.
      open(opal95_table_unit,file=opal95_table_path,status='OLD',access='SEQUENTIAL')
      fmt_start=1

   10 read(opal95_table_unit,1,end=9999) header_line
      if ( header_line(fmt_start:fmt_start+4).ne.'TABLE' ) go to 10
    1 format(a)
      read(header_line,'(36X,F7.4,11X,F7.4)') xx,zz
      iz = 1
      ix = 1
      nn = 1
!     ENTRY POINT FOR GETTING NEW TABLES.
   15 continue
      n = opal95_table_start_index(iz)+ix
      if (ix.lt.num_x) then
         xxt = opal95_grid_x(ix)
      else
         xxt = 1.0d0 - opal95_grid_z(iz)
      endif
      if (opal95_grid_z(iz).ne.zz .or. xxt.ne.xx) then
       write(short_file_unit,*)' OPAL95: Z ERROR INCOMPATIBLE TABLE'
         stop
      endif

!     READ IN HEADER INFO: GRID IN RHO/T6**3
      read(opal95_table_unit,20,end=1000) (opal95_grid_logr(i),i=1,num_d)
   20 format(///,4x,19f7.1,/)
!     READ IN FULL TABLE: LOG CAPPA AS A FUNCTION OF LOG T AND
!     LOG R = RHO/T6**3
      do i = 1,57
         read(opal95_table_unit,30,end=9999) opal95_grid_logt(i), &
              (opal95_full_opacity(n,i,j),j=1,num_d)
      end do
   30 format(f4.2,19f7.3)
! MHP 12/97 NOW TREAT CORNER WITHOUT DATA.
      read(opal95_table_unit,31,end=9999) opal95_grid_logt(58), &
           (opal95_full_opacity(n,58,j),j=1,num_d-1)
   31 format(f4.2,18f7.3)

      opal95_full_opacity(n,58,19) = 9.999d0
      do i = 59,60
         read(opal95_table_unit,32,end=9999) opal95_grid_logt(i), &
              (opal95_full_opacity(n,i,j),j=1,num_d-2)
   32    format(f4.2,17f7.3)
         do j = 18,19
            opal95_full_opacity(n,i,j) = 9.999d0
         end do
      end do
      do i = 61,64
         read(opal95_table_unit,33,end=9999) opal95_grid_logt(i), &
              (opal95_full_opacity(n,i,j),j=1,num_d-3)
   33    format(f4.2,16f7.3)
         do j = 17,19
            opal95_full_opacity(n,i,j) = 9.999d0
         end do
      end do
      do i = 65,69
         read(opal95_table_unit,34,end=9999) opal95_grid_logt(i), &
              (opal95_full_opacity(n,i,j),j=1,num_d-4)
   34    format(f4.2,15f7.3)
         do j = 16,19
            opal95_full_opacity(n,i,j) = 9.999d0
         end do
      end do
      i = 70
      read(opal95_table_unit,35,end=9999) opal95_grid_logt(i), &
           (opal95_full_opacity(n,i,j),j=1,num_d-5)
   35 format(f4.2,14f7.3)
      do j = 15,19
         opal95_full_opacity(n,i,j) = 9.999d0
      end do

!     EXIT IF CORRECT NUMBER OF TABLES READ IN.
      if (nn.ge.num_xz) goto 1000
      nn = nn + 1
!     NEED TO ACCOUNT FOR FEWER X VALUES AT HIGHER Z.
      if (ix.le.8) then
         nmax = num_z
      else if (ix.eq.9) then
         nmax = num_z - 3
      else
         nmax = num_z - 1
      endif
!     READ IN NEXT METAL ABUNDANCE (AT FIXED X) OR READ IN FIRST METAL
!     ABUNDANCE AT NEXT X.
      if (iz.lt.nmax) then
         iz = iz + 1
      else
         iz = 1
         ix = ix + 1
      endif

!     RETURN TO READ IN NEXT TABLE.
      read(opal95_table_unit,900) xx,zz
  900 format(/36x,f7.4,11x,f7.4)
      goto 15
 1000 continue

      close(opal95_table_unit)
!     NOW GENERATE A TABLE AT A FIXED VALUE OF Z.
!     NOTE THAT FOR METAL DIFFUSION A 4-D INTERPOLATION (IN X,Z,T,RHO)
!     IS PERFORMED RATHER THAN A LINEAR INTERPOLATION BETWEEN TWO FIXED
!     Z TABLES.
      target_z = opal95_single_table_z

      call op95ztab(target_z)
 9999 continue
      return
end subroutine ll95tbl
