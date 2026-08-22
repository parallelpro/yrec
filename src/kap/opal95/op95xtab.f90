!----------------------------------------------------------------------
! op95xtab
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original op95xtab.f; only variable names, source form, and comment
! style were updated.
!
! MHP 7/98 COMPUTE OPACITY TABLE AT FIXED X FROM THE OP95 TABLE AT
! THE MODEL Z. Builds the surface (T, rho) opacity table for the
! given hydrogen abundance by cubic interpolation across the 4
! nearest X tables in the fixed-Z table set (common/llot95/, filled
! by op95ztab).
subroutine op95xtab(hydrogen_fraction)

      use opacity_table_lib
      use numerics_lib
      implicit none
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      double precision, intent(in) :: hydrogen_fraction

      double precision :: table_x_nodes(4), x_weight(4)
      integer :: i, j, k, x_table_index
      integer :: table1_index, table2_index, table3_index, table4_index

      opacity_table%opal95_surface_x = hydrogen_fraction
!  FIND 4 NEAREST TABLES IN X TO DESIRED VALUE.
      if (hydrogen_fraction.le.0.8d0) then
! DON'T NEED TO WORRY ABOUT MISSING X TABLES AT HIGH Z.
         do i = 3,num_x-1
            if (opacity_table%opal95_grid_x(i).ge.hydrogen_fraction) then
               x_table_index = i - 2
               exit
            endif
         end do
         if (i > (num_x-1)) then
         x_table_index = num_x - 3
         end if
   10    continue
         do i = 1,4
            table_x_nodes(i) = opacity_table%opal95_grid_x(x_table_index+i-1)
         end do
      else if (opacity_table%opal95_fixed_z.le.0.04d0) then
! HIGH X TABLES PRESENT AT LOW Z.
         do i = 3,num_x-1
            if (opacity_table%opal95_grid_x(i).ge.hydrogen_fraction) then
               x_table_index = i - 2
               exit
            endif
         end do
         if (i > (num_x-1)) then
         x_table_index = num_x - 3
         end if
   20    continue
         do i = 1,4
            table_x_nodes(i) = opacity_table%opal95_grid_x(x_table_index+i-1)
         end do
      else if (opacity_table%opal95_fixed_z.ge.0.1d0) then
! USE TABLES 6-9
         x_table_index = 6
         do i = 1,4
            table_x_nodes(i) = opacity_table%opal95_grid_x(x_table_index+i-1)
         end do
      else
! IF Z IS BETWEEN 0.04 AND 0.1, TABLES 1-8 AND 10 EXIST.
! SINCE WE HAVE ALREADY DETERMINED THAT X > 0.8, USE 6-8 AND 10.
         x_table_index = 11
         table_x_nodes(1) = opacity_table%opal95_grid_x(6)
         table_x_nodes(2) = opacity_table%opal95_grid_x(7)
         table_x_nodes(3) = opacity_table%opal95_grid_x(8)
         table_x_nodes(4) = opacity_table%opal95_grid_x(10)
      endif
!
!  FIND INTERPOLATION FACTORS IN Z.
      call intrp2(table_x_nodes, x_weight, hydrogen_fraction)
! INDICES FOR 4 DESIRED COMPOSITIONS.
      if (x_table_index.lt.10) then
         table1_index = x_table_index
         table2_index = x_table_index + 1
         table3_index = x_table_index + 2
         table4_index = x_table_index + 3
      else
         table1_index = 6
         table2_index = 7
         table3_index = 8
         table4_index = 10
      endif
      do j = 1,num_t
         do k = 1,num_d
            opacity_table%opal95_surface_opacity(j,k) = x_weight(1)*opacity_table%opal95_fixed_z_opacity(table1_index,j,k) + &
            x_weight(2)*opacity_table%opal95_fixed_z_opacity(table2_index,j,k) + x_weight(3)*opacity_table%opal95_fixed_z_opacity(table3_index,j,k) + &
            x_weight(4)*opacity_table%opal95_fixed_z_opacity(table4_index,j,k)
         end do
      end do
      return
end subroutine op95xtab
