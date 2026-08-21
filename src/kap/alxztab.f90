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
subroutine alxztab(metal_fraction)

      use numerics_lib
      implicit none
      integer, parameter :: num_x = 7
      integer, parameter :: num_z = 15
      integer, parameter :: num_xz = 105
      integer, parameter :: num_t = 23
      integer, parameter :: num_d = 17
      integer, parameter :: num_xt = 8

      double precision, intent(in) :: metal_fraction

! common/galot/: ALEX95 low-T opacity table grids.
      double precision :: alex95_grid_logt(num_t), alex95_grid_x(num_x), &
           alex95_grid_logr(num_d), alex95_grid_z(num_z)
      common /galot/ alex95_grid_logt, alex95_grid_x, alex95_grid_logr, &
           alex95_grid_z

! common/alot/: ALEX95 low-T opacity table and cached-index state.
      double precision :: alex95_opacity(num_xt, num_t, num_d), &
           alex95_cached_x, alex95_cached_z
      integer :: alex95_index_x, alex95_index_t, alex95_index_r
      common /alot/ alex95_opacity, alex95_cached_x, alex95_cached_z, &
           alex95_index_x, alex95_index_t, alex95_index_r

! common/alotall/: full (X,Z) grid of ALEX95 low-T opacity tables.
      double precision :: alex95_full_opacity(num_xz, num_t, num_d)
      common /alotall/ alex95_full_opacity

      double precision :: interp_z(4), weight_z(4)
      save

      integer :: i, jj, kk, idz, ii1, ii2, ii3, ii4

!     LOCATE FOUR NEAREST TABLES IN Z
      do i=3,num_z-1
         if (metal_fraction.lt.alex95_grid_z(i)) then
            idz = i - 2
            goto 10
         endif
      end do
      idz = num_z - 3
   10 continue
      do i = 1,4
         interp_z(i) = alex95_grid_z(idz+i-1)
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
               alex95_opacity(i,jj,kk) = weight_z(1)*alex95_full_opacity(ii1,jj,kk)+ &
                    weight_z(2)*alex95_full_opacity(ii2,jj,kk) + &
                    weight_z(3)*alex95_full_opacity(ii3,jj,kk) + &
                    weight_z(4)*alex95_full_opacity(ii4,jj,kk)
            end do
         end do
      end do
      alex95_cached_z = metal_fraction
      return
end subroutine alxztab
