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
! interpolation in X, storing it in table slot 8 of
! opacity_table%alex95_full_opacity.
subroutine alex94_surface_table(hydrogen_fraction)

      use opacity_table_lib
      use numerics_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction

      double precision :: interp_x(4), weight_x(4)
      integer :: i, j

      opacity_table%alex95_cached_x = hydrogen_fraction
!     FIND 4 NEAREST TABLES IN X.
      opacity_table%alex95_index_x = alex_x_stencil_start(hydrogen_fraction)
      do i = 1,4
         interp_x(i) = opacity_table%alex95_grid_x(opacity_table%alex95_index_x+i-1)
      end do
      call intrp2(interp_x, weight_x, hydrogen_fraction)
      do i = 1,n_alex95_t
         do j = 1,n_alex95_d
            opacity_table%alex95_opacity(8,i,j) = weight_x(1)*opacity_table%alex95_opacity(opacity_table%alex95_index_x,i,j)+ &
                 weight_x(2)*opacity_table%alex95_opacity(opacity_table%alex95_index_x+1,i,j) + &
                 weight_x(3)*opacity_table%alex95_opacity(opacity_table%alex95_index_x+2,i,j) + &
                 weight_x(4)*opacity_table%alex95_opacity(opacity_table%alex95_index_x+3,i,j)
         end do
      end do
      return
end subroutine alex94_surface_table
