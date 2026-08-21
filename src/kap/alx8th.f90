!----------------------------------------------------------------------
! alx8th
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original alx8th.f; only variable names, source form, and comment
! style were updated.
!
! Generates the surface-X (envelope-abundance) table for the
! Alexander 1994 low-temperature opacity tables by 4-point Lagrangian
! interpolation in X, storing it in table slot 8 of common/alot/.
subroutine alx8th(hydrogen_fraction)

      use intrp2_mod
      implicit none
      integer, parameter :: num_x = 7
      integer, parameter :: num_z = 15
      integer, parameter :: num_t = 23
      integer, parameter :: num_d = 17
      integer, parameter :: num_xt = 8

      double precision, intent(in) :: hydrogen_fraction

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

      double precision :: interp_x(4), weight_x(4)
      save

      integer :: i, j

      alex95_cached_x = hydrogen_fraction
!     FIND 4 NEAREST TABLES IN X.
      if (hydrogen_fraction.lt.alex95_grid_x(4)) then
         if (hydrogen_fraction.gt.alex95_grid_x(3)) then
            alex95_index_x = 2
         else
            alex95_index_x = 1
         endif
      else
         if (hydrogen_fraction.gt.alex95_grid_x(5)) then
            alex95_index_x = 4
         else
            alex95_index_x = 3
         endif
      endif
      do i = 1,4
         interp_x(i) = alex95_grid_x(alex95_index_x+i-1)
      end do
      call intrp2(interp_x, weight_x, hydrogen_fraction)
      do i = 1,num_t
         do j = 1,num_d
            alex95_opacity(8,i,j) = weight_x(1)*alex95_opacity(alex95_index_x,i,j)+ &
                 weight_x(2)*alex95_opacity(alex95_index_x+1,i,j) + &
                 weight_x(3)*alex95_opacity(alex95_index_x+2,i,j) + &
                 weight_x(4)*alex95_opacity(alex95_index_x+3,i,j)
         end do
      end do
      return
end subroutine alx8th
