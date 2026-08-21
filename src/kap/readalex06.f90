!----------------------------------------------------------------------
! readalex06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original readalex06.f; only variable names, source form, and
! comment style were updated. FORMAT strings and the READ/OPEN
! statements that must byte-for-byte match the on-disk table file are
! preserved verbatim.
!
! MHP 3/09  SUBROUTINE FOR READING IN ALEXANDER LOW T OPACITIES.
!
! ALEXANDER OPACITIES FERGUSON ET AL. 2005  APJ,623,585
! FOR EACH  X = 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 1-Z
! INCREASING ORDER OF Z = 0.0, 0.00001, 0.00003, 0.0001, 0.0003,
!     (16)               0.001, 0.002, 0.004, 0.01, 0.02, 0.03,
!                        0.04, 0.05, 0.06, 0.08, 0.10
! DESCENDING ORDER OF LOGT 4.5-3.5 IN 0.05 DEX INCREMENTS
!     (85)                 3.5-2.9 IN 0.01 DEX INCREMENTS
!                          2.9-2.7 IN 0.05 DEX INCREMENTS
! ASCENDING ORDER OF  LOGR -8.0 -7.5 -7.0 -6.5 -6.0 -5.5 -5.0 -4.5 -4.0 -3.5
!     (19)                -3.0 -2.5 -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0
! WHERE R = RHO/T_6**3
! READ IN ALL Z VALUES AT FIXED X.  THEN GENERATE SURFACE X,Z TABLE.
! INTERPOLATION BETWEEN COMPOSITIONS IS 4-POINT LAGRANGIAN.
!
! MHP 8/25 Added file name to subroutine call
subroutine readalex06(alex06_table_path)

      use luout_lib
      implicit none
      integer, parameter :: num_x = 9
      integer, parameter :: num_z = 16
      integer, parameter :: num_xz = 143
      integer, parameter :: num_t = 85
      integer, parameter :: num_d = 19

      character(len=256), intent(in) :: alex06_table_path

!     ALEX LOW T OPACITY COMMON BLOCKS
! MHP 8/25 Removed file names from common block
      integer :: alex06_table_unit
      common /alex06/ alex06_table_unit

      double precision :: alex06_grid_logt(num_t), alex06_grid_x(num_x), &
           alex06_grid_logr(num_d), alex06_grid_z(num_z)
      common /galot06/ alex06_grid_logt, alex06_grid_x, alex06_grid_logr, &
           alex06_grid_z

      double precision :: alex06_opacity(num_t, num_d), alex06_cached_x, &
           alex06_cached_z
      integer :: alex06_index_x, alex06_index_t, alex06_index_r
      common /alot06/ alex06_opacity, alex06_cached_x, alex06_cached_z, &
           alex06_index_x, alex06_index_t, alex06_index_r

      double precision :: alex06_full_opacity(num_xz, num_t, num_d)
      common /alot06all/ alex06_full_opacity


!     LOCAL ARRAYS
      double precision :: row_logr_check(num_d), header_x, header_z, row_temp
      data alex06_grid_x/0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0,0.9d0, &
           1.0d0/
      data alex06_grid_z/0.0d0, 0.00001d0, 0.00003d0, 0.0001d0, 0.0003d0, &
           0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
           0.04d0, 0.05d0,0.06d0, 0.08d0, 0.10d0/
      data alex06_grid_logt/2.70d0,2.75d0,2.80d0,2.85d0,2.90d0,2.91d0,2.92d0, &
           2.93d0,2.94d0,2.95d0,2.96d0,2.97d0,2.98d0,2.99d0,3.00d0, &
           3.01d0,3.02d0,3.03d0,3.04d0,3.05d0,3.06d0,3.07d0,3.08d0, &
           3.09d0,3.10d0,3.11d0,3.12d0,3.13d0,3.14d0,3.15d0,3.16d0, &
           3.17d0,3.18d0,3.19d0,3.20d0,3.21d0,3.22d0,3.23d0,3.24d0, &
           3.25d0,3.26d0,3.27d0,3.28d0,3.29d0,3.30d0,3.31d0,3.32d0, &
           3.33d0,3.34d0,3.35d0,3.36d0,3.37d0,3.38d0,3.39d0,3.40d0, &
           3.41d0,3.42d0,3.43d0,3.44d0,3.45d0,3.46d0,3.47d0,3.48d0, &
           3.49d0,3.50d0,3.55d0,3.60d0,3.65d0,3.70d0,3.75d0,3.80d0, &
           3.85d0,3.90d0,3.95d0,4.00d0,4.05d0,4.10d0,4.15d0,4.20d0, &
           4.25d0,4.30d0,4.35d0,4.40d0,4.45d0,4.50d0/
      data alex06_grid_logr/-8.0d0,-7.5d0,-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0, &
           -4.5d0,-4.0d0,-3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0, &
           -0.5d0, 0.0d0, 0.5d0, 1.0d0/
!     INITIALIZE INDEX OF PREVIOUS CLOSEST POINTS AND ENVELOPE ABUNDANCES
      data alex06_cached_x,alex06_cached_z,alex06_index_x,alex06_index_t, &
           alex06_index_r/0.0d0,0.0d0,4,43,10/
      save

      integer :: i, ii, j, jj, k, kk

      open(unit=alex06_table_unit,file=alex06_table_path)
!     READ IN INITIAL X AND Z; ENSURE THAT THEY HAVE THE EXPECTED VALUES.
      do i = 1,num_x-1
         do ii = 1, num_z
!           INDEX FOR STORING TABLES : 15 SETS OF METAL ABUNDANCE WITH 7 SETS
!           OF HYDROGEN FOR EACH.  THESE ARE STORED IN THE ORDER (X,Z) OF
!           (0,0),(0,0.1),...(0.9,0),(0.9,0.1).  The final tables are defined
!           at X = 1-Z except for the 0.1 Z case (already read in at X=0.9).
!           HEADER INFORMATION: X AND Z; CHECK FOR CONSISTENCY
            read(alex06_table_unit,10) header_x, header_z
   10       format(30x,f9.6,7x,f9.6//)
            if (header_x.ne.alex06_grid_x(i) .or. header_z.ne.alex06_grid_z(ii)) then
               write(*,15) alex06_grid_x(i), header_x, alex06_grid_z(ii), header_z
   15          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
!               STOP
            endif
            read(alex06_table_unit,20) (row_logr_check(k),k=1,num_d)
   20       format(6x,19f7.3)
            do kk=1,16
               if (row_logr_check(kk).ne.alex06_grid_logr(kk)) then
                  write(*,25) kk, alex06_grid_logr(kk), row_logr_check(kk)
   25             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL R',i3,2f7.3,' RUN STOPPED')
                  stop
               endif
            end do
            jj = (i-1)*num_z+ii
!           OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T.
!           STORE IN A NUMXZ*NUMT*NUMR ARRAY.
            do j = num_t,1,-1
               read(alex06_table_unit,30) row_temp, (alex06_full_opacity(jj,j,k),k=1,num_d)
   30          format(f5.3,1x,19f7.3)
               if (row_temp.ne.alex06_grid_logt(j)) then
                  write(*,35) j, alex06_grid_logt(j), row_temp
   35             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T',i3,2f7.3,' RUN STOPPED')
                  stop
               endif
            end do
         end do
      end do
!     FINAL SET OF X TABLES ARE DEFINED AT X = 1 - Z AND THERE IS NOT ONE FOR
!     Z = 0.1 (ALREADY READ IN AT X=0.9).
      do ii = 1, num_z-1
!        INDEX FOR STORING TABLES : 15 SETS OF METAL ABUNDANCE WITH 7 SETS
!        OF HYDROGEN FOR EACH.  THESE ARE STORED IN THE ORDER (X,Z) OF
!        (0,0),(0,0.1),...(0.9,0),(0.9,0.1).  The final tables are defined
!        at X = 1-Z except for the 0.1 Z case (already read in at X=0.9).
!        HEADER INFORMATION: X AND Z; CHECK FOR CONSISTENCY
         read(alex06_table_unit,10) header_x, header_z
         if (header_x.ne.(1.0d0 - alex06_grid_z(ii)) .or. &
              header_z.ne.alex06_grid_z(ii)) then
            write(*,15) 1.0d0 - alex06_grid_z(ii), header_x, alex06_grid_z(ii), header_z
            stop
         endif
         read(alex06_table_unit,20) (row_logr_check(k),k=1,num_d)
         do kk=1,16
            if (row_logr_check(kk).ne.alex06_grid_logr(kk)) then
               write(*,25) kk, alex06_grid_logr(kk), row_logr_check(kk)
               stop
            endif
         end do
         jj = (num_x-1)*num_z+ii
!        OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T.
!        STORE IN A NUMXZ*NUMT*NUMR ARRAY.
         do j = num_t,1,-1
            read(alex06_table_unit,30) row_temp, (alex06_full_opacity(jj,j,k),k=1,num_d)
            if (row_temp.ne.alex06_grid_logt(j)) then
               write(*,35) j, alex06_grid_logt(j), row_temp
               stop
            endif
         end do
      end do
      close(unit=alex06_table_unit)
!     INITIALIZE FIXED Z,X TABLE
      do i = 1,num_t
         do j = 1,num_d
            alex06_opacity(i,j)=alex06_full_opacity(1,i,j)
         end do
      end do
      return
end subroutine readalex06
