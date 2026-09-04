!----------------------------------------------------------------------
! t6rinteos01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original t6rinteos01.f; only variable names, source form, and
! comment style were updated.
!
! OPAL 2001 EOS analogue of t6rinterp.f90 (see there for the general
! description). Interpolates the already X-interpolated table slice
! (opal_eos%x_interp_result_01) in T6 and density to produce opal_eos%esact_01.
subroutine t6rinteos01(slr, slt, ierr)

      use opal_eos_lib
      use luout_lib
      implicit none

      double precision, intent(in) :: slr, slt

! --- locals ---
      integer :: hi_loop_count, recompute_flag, cache_slot, t6_grid_idx
      double precision :: esactq, esact2, esactq2, dix, dix2
      double precision, external :: quadeos01

      integer, intent(out) :: ierr

      ierr = 0

      hi_loop_count = 0
      recompute_flag = 0

      do t6_grid_idx = opal_eos%t6_index_1_01, opal_eos%t6_index_1_01 + opal_eos%t6_interp_order_01
         cache_slot = 1
         hi_loop_count = hi_loop_count + 1
         opal_eos%rho_interp_lo_01(hi_loop_count) = quadeos01(recompute_flag, cache_slot, &
              slr, opal_eos%x_interp_result_01(t6_grid_idx,opal_eos%density_index_1_01), &
              opal_eos%x_interp_result_01(t6_grid_idx,opal_eos%density_index_2_01), &
              opal_eos%x_interp_result_01(t6_grid_idx,opal_eos%density_index_3_01), &
              opal_eos%density_grid_01(opal_eos%density_index_1_01), opal_eos%density_grid_01(opal_eos%density_index_2_01), &
              opal_eos%density_grid_01(opal_eos%density_index_3_01))
         if (opal_eos%density_interp_order_01.eq.3) then
            cache_slot = 2
            opal_eos%rho_interp_hi_01(hi_loop_count) = quadeos01(recompute_flag, cache_slot, &
                 slr, opal_eos%x_interp_result_01(t6_grid_idx,opal_eos%density_index_2_01), &
                 opal_eos%x_interp_result_01(t6_grid_idx,opal_eos%density_index_3_01), &
                 opal_eos%x_interp_result_01(t6_grid_idx,opal_eos%density_index_4_01), &
                 opal_eos%density_grid_01(opal_eos%density_index_2_01), opal_eos%density_grid_01(opal_eos%density_index_3_01), &
                 opal_eos%density_grid_01(opal_eos%density_index_4_01))
         end if
         recompute_flag = 1
      end do

      recompute_flag = 0
      cache_slot = 1
! ..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      opal_eos%esact_01 = quadeos01(recompute_flag, cache_slot, slt, opal_eos%rho_interp_lo_01(1), &
           opal_eos%rho_interp_lo_01(2), opal_eos%rho_interp_lo_01(3), opal_eos%t6_grid_01(opal_eos%t6_index_1_01), &
           opal_eos%t6_grid_01(opal_eos%t6_index_2_01), opal_eos%t6_grid_01(opal_eos%t6_index_3_01))
      if (opal_eos%density_interp_order_01.eq.3) then
! .....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
         esactq = quadeos01(recompute_flag, cache_slot, slt, opal_eos%rho_interp_hi_01(1), &
              opal_eos%rho_interp_hi_01(2), opal_eos%rho_interp_hi_01(3), opal_eos%t6_grid_01(opal_eos%t6_index_1_01), &
              opal_eos%t6_grid_01(opal_eos%t6_index_2_01), opal_eos%t6_grid_01(opal_eos%t6_index_3_01))
      end if
      if (opal_eos%t6_interp_order_01.eq.3) then
! .....    eos(i) in lower-left 3x3.
         esact2 = quadeos01(recompute_flag, cache_slot, slt, opal_eos%rho_interp_lo_01(2), &
              opal_eos%rho_interp_lo_01(3), opal_eos%rho_interp_lo_01(4), opal_eos%t6_grid_01(opal_eos%t6_index_2_01), &
              opal_eos%t6_grid_01(opal_eos%t6_index_3_01), opal_eos%t6_grid_01(opal_eos%t6_index_4_01))
! .....    eos(i) smoothed in left 3x4
         dix = (opal_eos%t6_grid_01(opal_eos%t6_index_3_01) - slt)*opal_eos%t6_grid_spacing_inv_01(opal_eos%t6_index_3_01)
         opal_eos%esact_01 = opal_eos%esact_01*dix + esact2*(1.0d0 - dix)
! endif   ! moved to loc a
         if (opal_eos%density_interp_order_01.eq.3) then

! .....     eos(i) in upper-right 3x3.
            esactq2 = quadeos01(recompute_flag, cache_slot, slt, &
                 opal_eos%rho_interp_hi_01(2), opal_eos%rho_interp_hi_01(3), opal_eos%rho_interp_hi_01(4), &
                 opal_eos%t6_grid_01(opal_eos%t6_index_2_01), opal_eos%t6_grid_01(opal_eos%t6_index_3_01), &
                 opal_eos%t6_grid_01(opal_eos%t6_index_4_01))
            esactq = esactq*dix + esactq2*(1.0d0 - dix)
         end if
      end if  ! loc a
!
      if (opal_eos%density_interp_order_01.eq.3) then
         dix2 = (opal_eos%density_grid_01(opal_eos%density_index_3_01) - slr)* &
              opal_eos%density_grid_spacing_inv_01(opal_eos%density_index_3_01)
         if (opal_eos%t6_interp_order_01.eq.3) then
! .....        eos(i) smoothed in both log(T6) and log(R)
            opal_eos%esact_01 = opal_eos%esact_01*dix2 + esactq*(1.0d0 - dix2)
         end if
      end if
      if (opal_eos%esact_01.gt.1.0d+15) then
         write(run_log_unit,'("Interpolation indices out of range", &
              &";please report conditions.")')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if

      return
end subroutine t6rinteos01
