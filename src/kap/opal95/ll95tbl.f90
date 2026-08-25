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
subroutine ll95tbl(opal95_table_path, ierr)
      use star_info_lib, only: star

      use opacity_table_lib
      use luout_lib
      implicit none
      integer, intent(out) :: ierr
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      character(len=256), intent(in) :: opal95_table_path



! common/newopac/: ll95tbl.f declares an OUT-OF-SYNC 16-member version
! of this block (10 doubles + 6 logicals) vs. the canonical 17-member
! layout (10 doubles + 7 logicals) used by kap_lib.f90's kap_get/setupopac.f90
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
! opal95_single_table_z is now use-associated from const_lib -- its
! out-of-sync 5th-slot position happens to be identical in both this
! file's buggy 16-member layout and the canonical 17-member layout, so
! use-association is correct here despite the historic mismatch.


! opal95_table_start_index/opal95_num_x_at_z/opal95_grid_x/
! opal95_grid_z/opal95_index_z/opal95_index_x/opal95_index_t/
! opal95_index_rho defaults moved to opacity_table_lib.f90: DATA can
! no longer target them here now that they're use-associated.
! LOCAL VECTOR, USED TO SKIP LONG HEADER.
      character(len=132) :: header_line

      integer :: fmt_start, iz, ix, nn, n, i, j, nmax
      integer :: read_status
      double precision :: xx, zz, xxt, target_z

!     READ IN OPACITY TABLES, SKIPPING HEADER.  HEY, IT WORKS.
      ierr = 0
      open(star%ctrl%opal95_table_unit,file=opal95_table_path,status='OLD',access='SEQUENTIAL')
      fmt_start=1

      do
         read(star%ctrl%opal95_table_unit,1,end=9999) header_line
      if (header_line(fmt_start:fmt_start+4).eq.'TABLE') exit
      end do
    1 format(a)
      read(header_line,'(36X,F7.4,11X,F7.4)') xx,zz
      iz = 1
      ix = 1
      nn = 1
!     ENTRY POINT FOR GETTING NEW TABLES.
      table_loop: do
      n = opacity_table%opal95_table_start_index(iz)+ix
      if (ix.lt.num_x) then
         xxt = opacity_table%opal95_grid_x(ix)
      else
         xxt = 1.0d0 - opacity_table%opal95_grid_z(iz)
      endif
      if (opacity_table%opal95_grid_z(iz).ne.zz .or. xxt.ne.xx) then
       write(short_file_unit,*)' OPAL95: Z ERROR INCOMPATIBLE TABLE'
         ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
         ierr = 1
         return
      endif

!     READ IN HEADER INFO: GRID IN RHO/T6**3
      read(star%ctrl%opal95_table_unit,20,iostat=read_status) (opacity_table%opal95_grid_logr(i),i=1,num_d)
      if (read_status .lt. 0) exit table_loop   ! was end=1000
      if (read_status .gt. 0) stop 'OPAL95 TABLE READ ERROR'
   20 format(///,4x,19f7.1,/)
!     READ IN FULL TABLE: LOG CAPPA AS A FUNCTION OF LOG T AND
!     LOG R = RHO/T6**3
      do i = 1,57
         read(star%ctrl%opal95_table_unit,30,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,num_d)
      end do
   30 format(f4.2,19f7.3)
! MHP 12/97 NOW TREAT CORNER WITHOUT DATA.
      read(star%ctrl%opal95_table_unit,31,end=9999) opacity_table%opal95_grid_logt(58), &
           (opacity_table%opal95_full_opacity(n,58,j),j=1,num_d-1)
   31 format(f4.2,18f7.3)

      opacity_table%opal95_full_opacity(n,58,19) = 9.999d0
      do i = 59,60
         read(star%ctrl%opal95_table_unit,32,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,num_d-2)
   32    format(f4.2,17f7.3)
         do j = 18,19
            opacity_table%opal95_full_opacity(n,i,j) = 9.999d0
         end do
      end do
      do i = 61,64
         read(star%ctrl%opal95_table_unit,33,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,num_d-3)
   33    format(f4.2,16f7.3)
         do j = 17,19
            opacity_table%opal95_full_opacity(n,i,j) = 9.999d0
         end do
      end do
      do i = 65,69
         read(star%ctrl%opal95_table_unit,34,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,num_d-4)
   34    format(f4.2,15f7.3)
         do j = 16,19
            opacity_table%opal95_full_opacity(n,i,j) = 9.999d0
         end do
      end do
      i = 70
      read(star%ctrl%opal95_table_unit,35,end=9999) opacity_table%opal95_grid_logt(i), &
           (opacity_table%opal95_full_opacity(n,i,j),j=1,num_d-5)
   35 format(f4.2,14f7.3)
      do j = 15,19
         opacity_table%opal95_full_opacity(n,i,j) = 9.999d0
      end do

!     EXIT IF CORRECT NUMBER OF TABLES READ IN.
      if (nn.lt.num_xz) then
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
      read(star%ctrl%opal95_table_unit,900) xx,zz
  900 format(/36x,f7.4,11x,f7.4)
      cycle table_loop
      end if
      exit table_loop
      end do table_loop

      close(star%ctrl%opal95_table_unit)
!     NOW GENERATE A TABLE AT A FIXED VALUE OF Z.
!     NOTE THAT FOR METAL DIFFUSION A 4-D INTERPOLATION (IN X,Z,T,RHO)
!     IS PERFORMED RATHER THAN A LINEAR INTERPOLATION BETWEEN TWO FIXED
!     Z TABLES.
      target_z = star%ctrl%opal95_single_table_z

      call op95ztab(target_z, ierr)
      if (ierr /= 0) return
 9999 continue
      return
end subroutine ll95tbl
