!----------------------------------------------------------------------
! op95ztab
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original op95ztab.f; only variable names, source form, and comment
! style were updated.
!
! MHP 7/98 COMPUTE OPACITY TABLE AT FIXED Z FROM THE FULL OP95 SET.
! Builds the fixed-Z (X, T, rho) opacity table by cubic interpolation
! across the 4 nearest Z tables in the full OPAL95 table set
! (common/llot95a/, filled by ll95tbl).
subroutine op95ztab(metal_fraction, ierr)

      use opacity_table_lib
      use numerics_lib
      implicit none
      integer, intent(out) :: ierr
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      double precision, intent(in) :: metal_fraction

      double precision :: table_z_nodes(4), z_weight(4), z_weight_hix(4)
      integer :: i, j, k, z_table_index, z_table_index2
      integer :: table1_index, table2_index, table3_index, table4_index

      ierr = 0
      opacity_table%opal95_fixed_z = metal_fraction
!  FIND 4 NEAREST TABLES IN Z TO DESIRED VALUE.
      do i = 3,num_z-1
         if (opacity_table%opal95_grid_z(i).ge.metal_fraction) then
            z_table_index = i - 2
            exit
         endif
      end do
      if (i > (num_z-1)) then
! DESIRED Z > 0.1D0; STOP
      write(*,5)metal_fraction
    5 format(1x,'DESIRED Z',f10.6,'OUTSIDE OP95 TABLE RANGE'/ &
           ' RUN STOPPED')
      ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
      ierr = 1
      return
      end if
   10 continue
!  FIND INTERPOLATION FACTORS IN Z.
      do i = 1,4
         table_z_nodes(i) = opacity_table%opal95_grid_z(z_table_index+i-1)
      end do
      call intrp2(table_z_nodes, z_weight, metal_fraction)
!  FOR X TABLES 1-8, ALL VALUES OF Z PRESENT; PERFORM INTERPOLATION.
      do i = 1,8
! INDICES FOR 4 DESIRED COMPOSITIONS.
         table1_index = opacity_table%opal95_table_start_index(z_table_index)+i
         table2_index = opacity_table%opal95_table_start_index(z_table_index+1)+i
         table3_index = opacity_table%opal95_table_start_index(z_table_index+2)+i
         table4_index = opacity_table%opal95_table_start_index(z_table_index+3)+i
         do j = 1,num_t
            do k = 1,num_d
               opacity_table%opal95_fixed_z_opacity(i,j,k) = z_weight(1)*opacity_table%opal95_full_opacity(table1_index,j,k)+ &
               z_weight(2)*opacity_table%opal95_full_opacity(table2_index,j,k)+z_weight(3)*opacity_table%opal95_full_opacity(table3_index,j,k)+ &
               z_weight(4)*opacity_table%opal95_full_opacity(table4_index,j,k)
            end do
         end do
      end do
!  FOR X TABLES 9 AND 10, HIGH VALUES OF Z ARE NOT PRESENT
!  OMIT TABLE 9 (X=0.95) IF DESIRED Z > 0.04
      if (.not. (metal_fraction.ge.0.04d0)) then
!  OTHERWISE, CHECK TO ENSURE THAT THE 4 Z TABLES USED HAVE Z < 0.04
!  ADJUST INTERPOLATION FACTORS IF NEEDED
      if (opacity_table%opal95_grid_z(z_table_index+3).gt.0.04d0) then
         z_table_index2 = num_z - 6
         do i = 1,4
            table_z_nodes(i) = opacity_table%opal95_grid_z(z_table_index2+i-1)
         end do
         call intrp2(table_z_nodes, z_weight_hix, metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      endif
      table1_index = opacity_table%opal95_table_start_index(z_table_index2)+9
      table2_index = opacity_table%opal95_table_start_index(z_table_index2+1)+9
      table3_index = opacity_table%opal95_table_start_index(z_table_index2+2)+9
      table4_index = opacity_table%opal95_table_start_index(z_table_index2+3)+9
      do j = 1,num_t
         do k = 1,num_d
            opacity_table%opal95_fixed_z_opacity(9,j,k) = z_weight_hix(1)*opacity_table%opal95_full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*opacity_table%opal95_full_opacity(table2_index,j,k)+z_weight_hix(3)*opacity_table%opal95_full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*opacity_table%opal95_full_opacity(table4_index,j,k)
         end do
      end do
      end if
   20 continue
!  OMIT TABLE 10 (X = 1-Z) IF DESIRED Z >= 0.1
      if (.not. (metal_fraction.ge.0.1d0)) then
!  CHECK TO ENSURE THAT Z=0.1 TABLE IS NOT ONE OF THE 4 TABLES; ADJUST
!  INTERPOLATION FACTORS IF NEEDED.
      if (opacity_table%opal95_grid_z(z_table_index+3).ge.0.1d0) then
         z_table_index2 = num_z - 4
         do i = 1,4
            table_z_nodes(i) = opacity_table%opal95_grid_z(z_table_index2+i-1)
         end do
         call intrp2(table_z_nodes, z_weight_hix, metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      endif
      end if
   30 continue
      table1_index = opacity_table%opal95_table_start_index(z_table_index2)+10
      table2_index = opacity_table%opal95_table_start_index(z_table_index2+1)+10
      table3_index = opacity_table%opal95_table_start_index(z_table_index2+2)+10
      table4_index = opacity_table%opal95_table_start_index(z_table_index2+3)+10
      do j = 1,num_t
         do k = 1,num_d
            opacity_table%opal95_fixed_z_opacity(10,j,k) = z_weight_hix(1)*opacity_table%opal95_full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*opacity_table%opal95_full_opacity(table2_index,j,k)+z_weight_hix(3)*opacity_table%opal95_full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*opacity_table%opal95_full_opacity(table4_index,j,k)
         end do
      end do
      return
end subroutine op95ztab
