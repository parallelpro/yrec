!----------------------------------------------------------------------
! t6rinterp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original t6rinterp.f; only variable names, source form, and comment
! style were updated.
!
! Second stage of the OPAL 1995 EOS interpolation: given the already
! X-interpolated table slice (opal_eos%x_interp_result, computed by esac.f90's
! call site right before this), interpolates in T6 and density on the
! surrounding 3x3 (or 4x4, smoothed by mixing overlapping quadratics)
! grid to produce opal_eos%esact, the interpolated value of the current
! thermodynamic variable at (slt, slr) = (T6, density).
subroutine t6rinterp(slr, slt, ierr)

      use opal_eos_lib
      use luout_lib
      implicit none

      double precision, intent(in) :: slr, slt

      integer, parameter :: mx = 5, mv = 10, nr = 77, nt = 56







      save

! --- locals ---
      integer :: hi_loop_count, recompute_flag, cache_slot, t6_grid_idx
      double precision :: esactq, esact2, esactq2, dix, dix2
      double precision, external :: quad

      integer, intent(out) :: ierr

      ierr = 0

      hi_loop_count = 0
      recompute_flag = 0

      do t6_grid_idx = opal_eos%t6_index_1, opal_eos%t6_index_1 + opal_eos%t6_interp_order
         cache_slot = 1
         hi_loop_count = hi_loop_count + 1
         opal_eos%rho_interp_lo(hi_loop_count) = quad(recompute_flag, cache_slot, &
              slr, opal_eos%x_interp_result(t6_grid_idx,opal_eos%density_index_1), &
              opal_eos%x_interp_result(t6_grid_idx,opal_eos%density_index_2), &
              opal_eos%x_interp_result(t6_grid_idx,opal_eos%density_index_3), &
              opal_eos%density_grid(opal_eos%density_index_1), opal_eos%density_grid(opal_eos%density_index_2), &
              opal_eos%density_grid(opal_eos%density_index_3))
         if (opal_eos%density_interp_order.eq.3) then
            cache_slot = 2
            opal_eos%rho_interp_hi(hi_loop_count) = quad(recompute_flag, cache_slot, &
                 slr, opal_eos%x_interp_result(t6_grid_idx,opal_eos%density_index_2), &
                 opal_eos%x_interp_result(t6_grid_idx,opal_eos%density_index_3), &
                 opal_eos%x_interp_result(t6_grid_idx,opal_eos%density_index_4), &
                 opal_eos%density_grid(opal_eos%density_index_2), opal_eos%density_grid(opal_eos%density_index_3), &
                 opal_eos%density_grid(opal_eos%density_index_4))
         end if
         recompute_flag = 1
      end do

      recompute_flag = 0
      cache_slot = 1
! ..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      opal_eos%esact = quad(recompute_flag, cache_slot, slt, opal_eos%rho_interp_lo(1), &
           opal_eos%rho_interp_lo(2), opal_eos%rho_interp_lo(3), opal_eos%t6_grid(opal_eos%t6_index_1), &
           opal_eos%t6_grid(opal_eos%t6_index_2), opal_eos%t6_grid(opal_eos%t6_index_3))
      if (opal_eos%density_interp_order.eq.3) then
! .....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
         esactq = quad(recompute_flag, cache_slot, slt, opal_eos%rho_interp_hi(1), &
              opal_eos%rho_interp_hi(2), opal_eos%rho_interp_hi(3), opal_eos%t6_grid(opal_eos%t6_index_1), &
              opal_eos%t6_grid(opal_eos%t6_index_2), opal_eos%t6_grid(opal_eos%t6_index_3))
      end if
      if (opal_eos%t6_interp_order.eq.3) then
! .....    eos(i) in lower-left 3x3.
         esact2 = quad(recompute_flag, cache_slot, slt, opal_eos%rho_interp_lo(2), &
              opal_eos%rho_interp_lo(3), opal_eos%rho_interp_lo(4), opal_eos%t6_grid(opal_eos%t6_index_2), &
              opal_eos%t6_grid(opal_eos%t6_index_3), opal_eos%t6_grid(opal_eos%t6_index_4))
! .....    eos(i) smoothed in left 3x4
         dix = (opal_eos%t6_grid(opal_eos%t6_index_3) - slt)*opal_eos%t6_grid_spacing_inv(opal_eos%t6_index_3)
         opal_eos%esact = opal_eos%esact*dix + esact2*(1.0d0 - dix)
      end if
! NOTE: preserved verbatim from the original -- this block is a
! sibling of (not nested inside) the "opal_eos%t6_interp_order.eq.3" block
! above, so if opal_eos%density_interp_order.eq.3 but opal_eos%t6_interp_order.ne.3,
! dix here is whatever it was left as by a previous call (t6rinterp
! is SAVE'd), not freshly computed in this call.
      if (opal_eos%density_interp_order.eq.3) then

! .....     eos(i) in upper-right 3x3.
         esactq2 = quad(recompute_flag, cache_slot, slt, &
              opal_eos%rho_interp_hi(2), opal_eos%rho_interp_hi(3), opal_eos%rho_interp_hi(4), &
              opal_eos%t6_grid(opal_eos%t6_index_2), opal_eos%t6_grid(opal_eos%t6_index_3), &
              opal_eos%t6_grid(opal_eos%t6_index_4))
         esactq = esactq*dix + esactq2*(1.0d0 - dix)
      end if
!
      if (opal_eos%density_interp_order.eq.3) then
         dix2 = (opal_eos%density_grid(opal_eos%density_index_3) - slr)* &
              opal_eos%density_grid_spacing_inv(opal_eos%density_index_3)
         if (opal_eos%t6_interp_order.eq.3) then
! .....        eos(i) smoothed in both log(T6) and log(R)
            opal_eos%esact = opal_eos%esact*dix2 + esactq*(1.0d0 - dix2)
         end if
      end if
      if (opal_eos%esact.gt.1.0d+15) then
         write(short_file_unit,'("INTERPOLATION INDICES OUT OF RANGE", &
              &";PLEASE REPORT CONDITIONS.")')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if

      return
end subroutine t6rinterp
