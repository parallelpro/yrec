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
! skipping the header, then builds the fixed-Z table via
! opal95_fixed_z_table. Any premature end of file (fewer than n_opal95_tables
! tables, or a truncated table body) is reported and returns ierr = 1.
subroutine ll95tbl(opal95_table_path, ierr)
      use star_info_lib, only: star

      use opacity_table_lib
      use luout_lib
      implicit none
      integer, intent(out) :: ierr

      character(len=256), intent(in) :: opal95_table_path

! opal95_table_start_index/opal95_num_x_at_z/opal95_grid_x/
! opal95_grid_z/opal95_index_z/opal95_index_x/opal95_index_t/
! opal95_index_rho defaults live in opacity_table_lib.f90.
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
      if (ix.lt.n_opal95_x) then
         xxt = opacity_table%opal95_grid_x(ix)
      else
         xxt = 1.0d0 - opacity_table%opal95_grid_z(iz)
      endif
      if (opacity_table%opal95_grid_z(iz).ne.zz .or. xxt.ne.xx) then
       write(run_log_unit,*)' OPAL95: Z ERROR INCOMPATIBLE TABLE'
         ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
         ierr = 1
         return
      endif

!     READ IN HEADER INFO: GRID IN RHO/T6**3
      read(star%ctrl%opal95_table_unit,20,iostat=read_status) (opacity_table%opal95_grid_logr(i),i=1,n_opal95_d)
      if (read_status .lt. 0) then
         write(*,*) 'll95tbl: OPAL95 opacity table file ended after ', nn-1, &
              ' of ', n_opal95_tables, ' tables'
         close(star%ctrl%opal95_table_unit)
         ierr = 1
         return
      end if
      if (read_status .gt. 0) then
         write(*,*) 'll95tbl: malformed OPAL95 opacity table (header read failed)'
         ierr = 1
         return
      end if
   20 format(///,4x,19f7.1,/)
!     READ IN FULL TABLE: LOG CAPPA AS A FUNCTION OF LOG T AND
!     LOG R = RHO/T6**3
      do i = 1,57
         read(star%ctrl%opal95_table_unit,30,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d)
      end do
   30 format(f4.2,19f7.3)
! MHP 12/97 NOW TREAT CORNER WITHOUT DATA.
      read(star%ctrl%opal95_table_unit,31,end=9999) opacity_table%opal95_grid_logt(58), &
           (opacity_table%opal95_full_opacity(n,58,j),j=1,n_opal95_d-1)
   31 format(f4.2,18f7.3)

      opacity_table%opal95_full_opacity(n,58,n_opal95_d) = opal95_missing_opacity
      do i = 59,60
         read(star%ctrl%opal95_table_unit,32,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-2)
   32    format(f4.2,17f7.3)
         do j = n_opal95_d-1,n_opal95_d
            opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      do i = 61,64
         read(star%ctrl%opal95_table_unit,33,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-3)
   33    format(f4.2,16f7.3)
         do j = n_opal95_d-2,n_opal95_d
            opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      do i = 65,69
         read(star%ctrl%opal95_table_unit,34,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-4)
   34    format(f4.2,15f7.3)
         do j = n_opal95_d-3,n_opal95_d
            opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      i = n_opal95_t
      read(star%ctrl%opal95_table_unit,35,end=9999) opacity_table%opal95_grid_logt(i), &
           (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-5)
   35 format(f4.2,14f7.3)
      do j = n_opal95_d-4,n_opal95_d
         opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
      end do

!     EXIT IF CORRECT NUMBER OF TABLES READ IN.
      if (nn.lt.n_opal95_tables) then
      nn = nn + 1
!     NEED TO ACCOUNT FOR FEWER X VALUES AT HIGHER Z.
      if (ix.le.8) then
         nmax = n_opal95_z
      else if (ix.eq.9) then
         nmax = n_opal95_z - 3
      else
         nmax = n_opal95_z - 1
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

      call opal95_fixed_z_table(target_z, ierr)
      return
!     PREMATURE END OF FILE INSIDE A TABLE BODY (end=9999 above).
 9999 write(*,*) 'll95tbl: OPAL95 opacity table file ended inside table ', nn
      close(star%ctrl%opal95_table_unit)
      ierr = 1
      return
end subroutine ll95tbl
