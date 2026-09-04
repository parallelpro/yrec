!----------------------------------------------------------------------
! alxztab
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original alxztab.f; only variable names, source form, and comment
! style were updated.
!
! Generates the fixed-Z table (one table per X, interpolated in Z)
! from the full set of Alexander 1994 low-temperature opacity tables,
! by 4-point Lagrangian interpolation in Z.
subroutine alex94_fixed_z_table(metal_fraction)

      use opacity_table_lib
      use numerics_lib
      implicit none
      integer, parameter :: num_x = 7
      integer, parameter :: num_z = 15
      integer, parameter :: num_t = 23
      integer, parameter :: num_d = 17

      double precision, intent(in) :: metal_fraction

      double precision :: interp_z(4), weight_z(4)
      integer :: i, jj, kk, idz, ii1, ii2, ii3, ii4

!     LOCATE FOUR NEAREST TABLES IN Z
      do i=3,num_z-1
         if (metal_fraction.lt.opacity_table%alex95_grid_z(i)) then
            idz = i - 2
            exit
         endif
      end do
      if (i > (num_z-1)) then
      idz = num_z - 3
      end if
      do i = 1,4
         interp_z(i) = opacity_table%alex95_grid_z(idz+i-1)
      end do
!     GET INTERPOLATION FACTORS FOR Z.
      call intrp2(interp_z, weight_z, metal_fraction)
!     INTERPOLATE IN Z AT FIXED X.
      do i = 1,num_x
!        INDICES FOR TABLES: SETS OF 7 X AT FIXED Z, STARTING FROM
!        Z= 0 AND GOING TO Z = 0.1.
         ii1 = i + num_x*(idz-1)
         ii2 = ii1 + num_x
         ii3 = ii2 + num_x
         ii4 = ii3 + num_x
         do jj = 1,num_t
            do kk = 1,num_d
               opacity_table%alex95_opacity(i,jj,kk) = weight_z(1)*opacity_table%alex95_full_opacity(ii1,jj,kk)+ &
                    weight_z(2)*opacity_table%alex95_full_opacity(ii2,jj,kk) + &
                    weight_z(3)*opacity_table%alex95_full_opacity(ii3,jj,kk) + &
                    weight_z(4)*opacity_table%alex95_full_opacity(ii4,jj,kk)
            end do
         end do
      end do
      opacity_table%alex95_cached_z = metal_fraction
      return
end subroutine alex94_fixed_z_table
