!----------------------------------------------------------------------
! alxtbl
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original alxtbl.f; only variable names, source form, and comment
! style were updated.
!
! ALEXANDER OPACITIES   APJ,437,879,1994
! IN EACH 7 FILE FOR X = 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8
! DESCENDING ORDER OF Z = 0.0, 0.00001, 0.00003, 0.0001, 0.0003,
!     (15)               0.001, 0.002, 0.004, 0.01, 0.02, 0.03,
!                        0.04, 0.06, 0.08, 0.10
! DESCENDING ORDER OF LOGT  4.10 4.05 4.00 3.95 3.90 3.85 3.80 3.75
!     (23)                 3.70 3.65 3.60 3.55 3.50 3.45 3.40 3.35
!                          3.30 3.25 3.20 3.15 3.10 3.05 3.00
! ASCENDING ORDER OF  LOGR -7.0 -6.5 -6.0 -5.5 -5.0 -4.5 -4.0 -3.5
!     (17)                -3.0 -2.5 -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0
! WHERE R = RHO/T_6**3
! READ IN ALL Z VALUES AT FIXED X AND GENERATE OPACITY ENTRIES AT
! EITHER ONE OR TWO FIXED VALUES OF Z.  THEN GENERATE SURFACE X TABLE
! IF MODEL HAS NO DIFFUSION.  INTERPOLATION BETWEEN COMPOSITIONS IS
! 4-POINT LAGRANGIAN.
!
! MHP 8/25 Added file name to subroutine call
subroutine alxtbl(alex95_table_paths)

      implicit none
      integer, parameter :: num_x = 7
      integer, parameter :: num_z = 15
      integer, parameter :: num_xz = 105
      integer, parameter :: num_t = 23
      integer, parameter :: num_d = 17
      integer, parameter :: num_xt = 8

      character(len=256), intent(in) :: alex95_table_paths(7)

! MHP 8/25 Removed all character strings from common blocks
      integer :: alex95_table_unit
      common /alexo/ alex95_table_unit

! common/newopac/: only alex_table_z1 is used here.
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

! common/galot/: ALEX95 low-T opacity table grids.
      double precision :: alex95_grid_logt(num_t), alex95_grid_x(num_x), &
           alex95_grid_logr(num_d), alex95_grid_z(num_z)
      common /galot/ alex95_grid_logt, alex95_grid_x, alex95_grid_logr, &
           alex95_grid_z

! common/alot/: ALEX95 low-T opacity table and cached-index state.
! NOTE: the file-read loop below reads directly into alex95_index_t
! (the common's cached T-grid index), matching the original ALXTBL,
! which reads its per-record row-index field from the table file into
! the variable named IT -- the SAME "IT" as common/ALOT/'s T-grid
! index cache (an accidental name collision in the original F77, not
! something introduced by this conversion). This leaves
! alex95_index_t holding whatever value was in the last-read record's
! row-index field, not a valid T-grid index, when this routine
! returns; downstream callers (yalo3d.f90/getalex... path) clamp IT
! into range before using it, so the collision is preserved exactly
! but is understood to be harmless in practice.
      double precision :: alex95_opacity(num_xt, num_t, num_d), &
           alex95_cached_x, alex95_cached_z
      integer :: alex95_index_x, alex95_index_t, alex95_index_r
      common /alot/ alex95_opacity, alex95_cached_x, alex95_cached_z, &
           alex95_index_x, alex95_index_t, alex95_index_r

! common/alotall/: full (X,Z) grid of ALEX95 low-T opacity tables.
      double precision :: alex95_full_opacity(num_xz, num_t, num_d)
      common /alotall/ alex95_full_opacity

! MHP 10/02 ISHORT UNDECLARED
      integer :: ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr

! --- local arrays ---
      double precision :: row_opacity_temp(num_d)
      data alex95_grid_x/0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0/
      data alex95_grid_z/0.0d0, 0.00001d0, 0.00003d0, 0.0001d0, 0.0003d0, &
           0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
           0.04d0, 0.06d0, 0.08d0, 0.10d0/
      data alex95_grid_logt/3.00d0,3.05d0,3.10d0,3.15d0,3.20d0,3.25d0,3.30d0, &
           3.35d0,3.40d0,3.45d0,3.50d0,3.55d0,3.60d0,3.65d0, &
           3.70d0,3.75d0,3.80d0,3.85d0,3.90d0,3.95d0,4.00d0, &
           4.05d0,4.10d0/
      data alex95_grid_logr/-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0,-4.5d0,-4.0d0, &
           -3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0, &
           0.0d0, 0.5d0, 1.0d0/
! INITIALIZE INDEX OF PREVIOUS CLOSEST POINTS AND ENVELOPE ABUNDANCES
      data alex95_cached_x,alex95_cached_z,alex95_index_x,alex95_index_t, &
           alex95_index_r/0.0d0,0.0d0,4,12,9/
      save

      integer :: table_index, i, ii, j, k, row_density_count
      double precision :: header_x, header_z, row_temp, row_logr0
      double precision :: target_z

      do table_index=1,num_x
!        OPEN EACH OF THE TABLES AT FIXED X WITH A RANGE OF Z.
         open(unit=alex95_table_unit,file=alex95_table_paths(table_index), &
              status='OLD', form='FORMATTED')
!        READ IN INITIAL X AND Z; ENSURE THAT THEY HAVE THE EXPECTED VALUES.
         do i = num_z,1,-1
!           INDEX FOR STORING TABLES : 15 SETS OF METAL ABUNDANCE WITH 7 SETS
!           OF HYDROGEN FOR EACH.  THESE ARE STORED IN THE ORDER (X,Z) OF
!           (0,0),(0.1,0),...(0.8,0),(0,0.00001),...
            ii = table_index + 7*(i-1)
!           HEADER INFORMATION: X AND Z
            read(alex95_table_unit,10) header_x, header_z
   10       format(18x,f6.2,2x,f7.2)
            if (header_x.ne.alex95_grid_x(table_index) .or. &
                 header_z.ne.alex95_grid_z(i)) then
               write(short_file_unit,20) header_x, alex95_grid_x(table_index), &
                    header_z, alex95_grid_z(i)
   20          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
               stop
            endif
!           OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T,
!           STARTING R.  STORE IN A NUMZ*NUMT*NUMR ARRAY.
            do j = num_t,1,-1
               read(alex95_table_unit,30) alex95_index_t, row_density_count, &
                    row_temp, row_logr0, (row_opacity_temp(k),k=1,num_d)
               do k = 1, num_d
                  alex95_full_opacity(ii,j,k) = row_opacity_temp(k)
               end do
   30          format(i2,i3,f6.3,f5.1,8f8.3/9f8.3)
               if (row_density_count.ne.17 .or. row_temp.ne.alex95_grid_logt(j) &
                    .or. row_logr0.ne.alex95_grid_logr(1)) then
                  write(short_file_unit,40) row_density_count, row_temp, &
                       alex95_grid_logt(j), row_logr0, alex95_grid_logr(1)
   40             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T,RHO',i3,4f7.2, &
                       ' RUN STOPPED')
                  stop
               endif
            end do
         end do
         close(unit=alex95_table_unit)
      end do
      target_z = alex_table_z1
      call alxztab(target_z)
      return
end subroutine alxtbl
