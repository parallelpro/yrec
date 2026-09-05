!----------------------------------------------------------------------
! t6rint_core
!----------------------------------------------------------------------
! Readability W3 (2026): the T6/density interpolation shared by the
! OPAL 2001 and 2006 EOS (t6rinteos01.f90 / t6rinteos06.f90). Those two
! routines were identical statement for statement apart from the state
! instance they addressed, the RECURSIVE attribute of the 2006 copy,
! and the text of the out-of-range message; the arithmetic lives here
! once, and each vintage's wrapper keeps its own message. Not used by
! the 1995 EOS (t6rinterp.f90 orders its esactq2 block differently).
!
! Interpolates the already X-interpolated table slice
! (v%x_interp_result) in density (slr) and T6 (slt) with the stencil
! indices esac01/esac06 stored in v (density_index_*, t6_index_*,
! *_interp_order), leaving the result in v%esact. The caller tests
! v%esact against the out-of-range sentinel.
subroutine t6rint_core(v, slr, slt)

      use opal_eos_lib
      implicit none

      type(opal_eos_vintage), intent(inout) :: v

      double precision, intent(in) :: slr, slt

! --- locals ---
      integer :: hi_loop_count, recompute_flag, cache_slot, t6_grid_idx
      double precision :: esactq, esact2, esactq2, dix, dix2
      double precision, external :: quad

      hi_loop_count = 0
      recompute_flag = 0

      do t6_grid_idx = v%t6_index_1, v%t6_index_1 + v%t6_interp_order
         cache_slot = 1
         hi_loop_count = hi_loop_count + 1
         v%rho_interp_lo(hi_loop_count) = quad(v%quad, recompute_flag, cache_slot, &
              slr, v%x_interp_result(t6_grid_idx,v%density_index_1), &
              v%x_interp_result(t6_grid_idx,v%density_index_2), &
              v%x_interp_result(t6_grid_idx,v%density_index_3), &
              v%density_grid(v%density_index_1), v%density_grid(v%density_index_2), &
              v%density_grid(v%density_index_3))
         if (v%density_interp_order.eq.3) then
            cache_slot = 2
            v%rho_interp_hi(hi_loop_count) = quad(v%quad, recompute_flag, cache_slot, &
                 slr, v%x_interp_result(t6_grid_idx,v%density_index_2), &
                 v%x_interp_result(t6_grid_idx,v%density_index_3), &
                 v%x_interp_result(t6_grid_idx,v%density_index_4), &
                 v%density_grid(v%density_index_2), v%density_grid(v%density_index_3), &
                 v%density_grid(v%density_index_4))
         end if
         recompute_flag = 1
      end do

      recompute_flag = 0
      cache_slot = 1
! ..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      v%esact = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(1), &
           v%rho_interp_lo(2), v%rho_interp_lo(3), v%t6_grid(v%t6_index_1), &
           v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      if (v%density_interp_order.eq.3) then
! .....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
         esactq = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_hi(1), &
              v%rho_interp_hi(2), v%rho_interp_hi(3), v%t6_grid(v%t6_index_1), &
              v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      end if
      if (v%t6_interp_order.eq.3) then
! .....    eos(i) in lower-left 3x3.
         esact2 = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(2), &
              v%rho_interp_lo(3), v%rho_interp_lo(4), v%t6_grid(v%t6_index_2), &
              v%t6_grid(v%t6_index_3), v%t6_grid(v%t6_index_4))
! .....    eos(i) smoothed in left 3x4
         dix = (v%t6_grid(v%t6_index_3) - slt)*v%t6_grid_spacing_inv(v%t6_index_3)
         v%esact = v%esact*dix + esact2*(1.0d0 - dix)
! endif   ! moved to loc a
         if (v%density_interp_order.eq.3) then

! .....     eos(i) in upper-right 3x3.
            esactq2 = quad(v%quad, recompute_flag, cache_slot, slt, &
                 v%rho_interp_hi(2), v%rho_interp_hi(3), v%rho_interp_hi(4), &
                 v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3), &
                 v%t6_grid(v%t6_index_4))
            esactq = esactq*dix + esactq2*(1.0d0 - dix)
         end if
      end if  ! loc a
!
      if (v%density_interp_order.eq.3) then
         dix2 = (v%density_grid(v%density_index_3) - slr)* &
              v%density_grid_spacing_inv(v%density_index_3)
         if (v%t6_interp_order.eq.3) then
! .....        eos(i) smoothed in both log(T6) and log(R)
            v%esact = v%esact*dix2 + esactq*(1.0d0 - dix2)
         end if
      end if

      return
end subroutine t6rint_core
