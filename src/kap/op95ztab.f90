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
subroutine op95ztab(metal_fraction)

      use intrp2_mod
      implicit none
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      double precision, intent(in) :: metal_fraction

! GRID ENTRIES FOR TEMPERATURE, ABUNDANCE (X), AND RHO/T6**3
! OPACITY AS A FUNCTION OF X, T, AND RHO/T6**3
      double precision :: opal95_fixed_z_opacity(num_x,num_t,num_d), opal95_fixed_z
      common /llot95/opal95_fixed_z_opacity, opal95_fixed_z
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
      double precision :: table_z_nodes(4), z_weight(4), z_weight_hix(4)
      save

      integer :: i, j, k, z_table_index, z_table_index2
      integer :: table1_index, table2_index, table3_index, table4_index

      opal95_fixed_z = metal_fraction
!  FIND 4 NEAREST TABLES IN Z TO DESIRED VALUE.
      do i = 3,num_z-1
         if (opal95_grid_z(i).ge.metal_fraction) then
            z_table_index = i - 2
            goto 10
         endif
      end do
! DESIRED Z > 0.1D0; STOP
      write(*,5)metal_fraction
    5 format(1x,'DESIRED Z',f10.6,'OUTSIDE OP95 TABLE RANGE'/ &
           ' RUN STOPPED')
      stop
   10 continue
!  FIND INTERPOLATION FACTORS IN Z.
      do i = 1,4
         table_z_nodes(i) = opal95_grid_z(z_table_index+i-1)
      end do
      call intrp2(table_z_nodes, z_weight, metal_fraction)
!  FOR X TABLES 1-8, ALL VALUES OF Z PRESENT; PERFORM INTERPOLATION.
      do i = 1,8
! INDICES FOR 4 DESIRED COMPOSITIONS.
         table1_index = opal95_table_start_index(z_table_index)+i
         table2_index = opal95_table_start_index(z_table_index+1)+i
         table3_index = opal95_table_start_index(z_table_index+2)+i
         table4_index = opal95_table_start_index(z_table_index+3)+i
         do j = 1,num_t
            do k = 1,num_d
               opal95_fixed_z_opacity(i,j,k) = z_weight(1)*opal95_full_opacity(table1_index,j,k)+ &
               z_weight(2)*opal95_full_opacity(table2_index,j,k)+z_weight(3)*opal95_full_opacity(table3_index,j,k)+ &
               z_weight(4)*opal95_full_opacity(table4_index,j,k)
            end do
         end do
      end do
!  FOR X TABLES 9 AND 10, HIGH VALUES OF Z ARE NOT PRESENT
!  OMIT TABLE 9 (X=0.95) IF DESIRED Z > 0.04
      if (metal_fraction.ge.0.04d0) goto 20
!  OTHERWISE, CHECK TO ENSURE THAT THE 4 Z TABLES USED HAVE Z < 0.04
!  ADJUST INTERPOLATION FACTORS IF NEEDED
      if (opal95_grid_z(z_table_index+3).gt.0.04d0) then
         z_table_index2 = num_z - 6
         do i = 1,4
            table_z_nodes(i) = opal95_grid_z(z_table_index2+i-1)
         end do
         call intrp2(table_z_nodes, z_weight_hix, metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      endif
      table1_index = opal95_table_start_index(z_table_index2)+9
      table2_index = opal95_table_start_index(z_table_index2+1)+9
      table3_index = opal95_table_start_index(z_table_index2+2)+9
      table4_index = opal95_table_start_index(z_table_index2+3)+9
      do j = 1,num_t
         do k = 1,num_d
            opal95_fixed_z_opacity(9,j,k) = z_weight_hix(1)*opal95_full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*opal95_full_opacity(table2_index,j,k)+z_weight_hix(3)*opal95_full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*opal95_full_opacity(table4_index,j,k)
         end do
      end do
   20 continue
!  OMIT TABLE 10 (X = 1-Z) IF DESIRED Z >= 0.1
      if (metal_fraction.ge.0.1d0) goto 30
!  CHECK TO ENSURE THAT Z=0.1 TABLE IS NOT ONE OF THE 4 TABLES; ADJUST
!  INTERPOLATION FACTORS IF NEEDED.
      if (opal95_grid_z(z_table_index+3).ge.0.1d0) then
         z_table_index2 = num_z - 4
         do i = 1,4
            table_z_nodes(i) = opal95_grid_z(z_table_index2+i-1)
         end do
         call intrp2(table_z_nodes, z_weight_hix, metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      endif
   30 continue
      table1_index = opal95_table_start_index(z_table_index2)+10
      table2_index = opal95_table_start_index(z_table_index2+1)+10
      table3_index = opal95_table_start_index(z_table_index2+2)+10
      table4_index = opal95_table_start_index(z_table_index2+3)+10
      do j = 1,num_t
         do k = 1,num_d
            opal95_fixed_z_opacity(10,j,k) = z_weight_hix(1)*opal95_full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*opal95_full_opacity(table2_index,j,k)+z_weight_hix(3)*opal95_full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*opal95_full_opacity(table4_index,j,k)
         end do
      end do
      return
end subroutine op95ztab
