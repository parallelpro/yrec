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
subroutine read_alex94_tables(alex95_table_paths, ierr)
      use star_info_lib, only: star

      use opacity_table_lib
      use luout_lib
      implicit none
      integer, intent(out) :: ierr
      integer, parameter :: num_x = 7
      integer, parameter :: num_z = 15
      integer, parameter :: num_xz = 105
      integer, parameter :: num_t = 23
      integer, parameter :: num_d = 17
      integer, parameter :: num_xt = 8

      character(len=256), intent(in) :: alex95_table_paths(7)







! --- local arrays ---
      double precision :: row_opacity_temp(num_d)
! alex95_grid_x/alex95_grid_z/alex95_grid_logt/alex95_grid_logr/
! alex95_cached_x/alex95_cached_z/alex95_index_x/alex95_index_t/
! alex95_index_r defaults moved to opacity_table_lib.f90: DATA can no
! longer target them here now that they're use-associated.
      integer :: table_index, i, ii, j, k, row_density_count
      double precision :: header_x, header_z, row_temp, row_logr0
      double precision :: target_z

      do table_index=1,num_x
!        OPEN EACH OF THE TABLES AT FIXED X WITH A RANGE OF Z.
         open(unit=star%ctrl%alex95_table_unit,file=alex95_table_paths(table_index), &
              status='OLD', form='FORMATTED')
!        READ IN INITIAL X AND Z; ENSURE THAT THEY HAVE THE EXPECTED VALUES.
         do i = num_z,1,-1
!           INDEX FOR STORING TABLES : 15 SETS OF METAL ABUNDANCE WITH 7 SETS
!           OF HYDROGEN FOR EACH.  THESE ARE STORED IN THE ORDER (X,Z) OF
!           (0,0),(0.1,0),...(0.8,0),(0,0.00001),...
            ii = table_index + 7*(i-1)
!           HEADER INFORMATION: X AND Z
            read(star%ctrl%alex95_table_unit,10) header_x, header_z
   10       format(18x,f6.2,2x,f7.2)
            if (header_x.ne.opacity_table%alex95_grid_x(table_index) .or. &
                 header_z.ne.opacity_table%alex95_grid_z(i)) then
               write(short_file_unit,20) header_x, opacity_table%alex95_grid_x(table_index), &
                    header_z, opacity_table%alex95_grid_z(i)
   20          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
               ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
               ierr = 1
               return
            endif
!           OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T,
!           STARTING R.  STORE IN A NUMZ*NUMT*NUMR ARRAY.
            do j = num_t,1,-1
               read(star%ctrl%alex95_table_unit,30) opacity_table%alex95_index_t, row_density_count, &
                    row_temp, row_logr0, (row_opacity_temp(k),k=1,num_d)
               do k = 1, num_d
                  opacity_table%alex95_full_opacity(ii,j,k) = row_opacity_temp(k)
               end do
   30          format(i2,i3,f6.3,f5.1,8f8.3/9f8.3)
               if (row_density_count.ne.17 .or. row_temp.ne.opacity_table%alex95_grid_logt(j) &
                    .or. row_logr0.ne.opacity_table%alex95_grid_logr(1)) then
                  write(short_file_unit,40) row_density_count, row_temp, &
                       opacity_table%alex95_grid_logt(j), row_logr0, opacity_table%alex95_grid_logr(1)
   40             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T,RHO',i3,4f7.2, &
                       ' RUN STOPPED')
                  ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
                  ierr = 1
                  return
               endif
            end do
         end do
         close(unit=star%ctrl%alex95_table_unit)
      end do
      ierr = 0
      target_z = star%ctrl%alex_table_z1
      call alex94_fixed_z_table(target_z)
      return
end subroutine read_alex94_tables
