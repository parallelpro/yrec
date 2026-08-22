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
subroutine readalex06(alex06_table_path, ierr)

      use opacity_table_lib
      use const_lib
      use luout_lib
      implicit none
      integer, intent(out) :: ierr
      integer, parameter :: num_x = 9
      integer, parameter :: num_z = 16
      integer, parameter :: num_xz = 143
      integer, parameter :: num_t = 85
      integer, parameter :: num_d = 19

      character(len=256), intent(in) :: alex06_table_path






!     LOCAL ARRAYS
      double precision :: row_logr_check(num_d), header_x, header_z, row_temp
! alex06_grid_x/alex06_grid_z/alex06_grid_logt/alex06_grid_logr/
! alex06_cached_x/alex06_cached_z/alex06_index_x/alex06_index_t/
! alex06_index_r defaults moved to opacity_table_lib.f90: DATA can no
! longer target them here now that they're use-associated.
      save

      integer :: i, ii, j, jj, k, kk

      ierr = 0
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
            if (header_x.ne.opacity_table%alex06_grid_x(i) .or. header_z.ne.opacity_table%alex06_grid_z(ii)) then
               write(*,15) opacity_table%alex06_grid_x(i), header_x, opacity_table%alex06_grid_z(ii), header_z
   15          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
!               STOP
            endif
            read(alex06_table_unit,20) (row_logr_check(k),k=1,num_d)
   20       format(6x,19f7.3)
            do kk=1,16
               if (row_logr_check(kk).ne.opacity_table%alex06_grid_logr(kk)) then
                  write(*,25) kk, opacity_table%alex06_grid_logr(kk), row_logr_check(kk)
   25             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL R',i3,2f7.3,' RUN STOPPED')
                  ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
                  ierr = 1
                  return
               endif
            end do
            jj = (i-1)*num_z+ii
!           OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T.
!           STORE IN A NUMXZ*NUMT*NUMR ARRAY.
            do j = num_t,1,-1
               read(alex06_table_unit,30) row_temp, (opacity_table%alex06_full_opacity(jj,j,k),k=1,num_d)
   30          format(f5.3,1x,19f7.3)
               if (row_temp.ne.opacity_table%alex06_grid_logt(j)) then
                  write(*,35) j, opacity_table%alex06_grid_logt(j), row_temp
   35             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T',i3,2f7.3,' RUN STOPPED')
                  ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
                  ierr = 1
                  return
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
         if (header_x.ne.(1.0d0 - opacity_table%alex06_grid_z(ii)) .or. &
              header_z.ne.opacity_table%alex06_grid_z(ii)) then
            write(*,15) 1.0d0 - opacity_table%alex06_grid_z(ii), header_x, opacity_table%alex06_grid_z(ii), header_z
            ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
            ierr = 1
            return
         endif
         read(alex06_table_unit,20) (row_logr_check(k),k=1,num_d)
         do kk=1,16
            if (row_logr_check(kk).ne.opacity_table%alex06_grid_logr(kk)) then
               write(*,25) kk, opacity_table%alex06_grid_logr(kk), row_logr_check(kk)
               ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
               ierr = 1
               return
            endif
         end do
         jj = (num_x-1)*num_z+ii
!        OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T.
!        STORE IN A NUMXZ*NUMT*NUMR ARRAY.
         do j = num_t,1,-1
            read(alex06_table_unit,30) row_temp, (opacity_table%alex06_full_opacity(jj,j,k),k=1,num_d)
            if (row_temp.ne.opacity_table%alex06_grid_logt(j)) then
               write(*,35) j, opacity_table%alex06_grid_logt(j), row_temp
               ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
               ierr = 1
               return
            endif
         end do
      end do
      close(unit=alex06_table_unit)
!     INITIALIZE FIXED Z,X TABLE
      do i = 1,num_t
         do j = 1,num_d
            opacity_table%alex06_opacity(i,j)=opacity_table%alex06_full_opacity(1,i,j)
         end do
      end do
      return
end subroutine readalex06
