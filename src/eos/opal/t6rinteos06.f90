!----------------------------------------------------------------------
! t6rinteos06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original t6rinteos06.f; only variable names, source form, and
! comment style were updated.
!
! OPAL 2006 EOS analogue of t6rinteos01.f90 (see there and
! t6rinterp.f90 for the general description). Declared RECURSIVE in
! the original; preserved verbatim even though nothing here actually
! recurses.
recursive subroutine t6rinteos06(v, slr, slt, ierr)

      use opal_eos_lib
      use luout_lib
      implicit none

      type(opal_eos_vintage), intent(inout) :: v

      double precision, intent(in) :: slr, slt

! --- locals ---
      integer :: hi_loop_count, recompute_flag, cache_slot, t6_grid_idx
      double precision :: esactq, esact2, esactq2, dix, dix2
      double precision, external :: quadeos06
      integer, intent(out) :: ierr

      ierr = 0

      hi_loop_count = 0
      recompute_flag = 0

      do t6_grid_idx = v%t6_index_1, v%t6_index_1 + v%t6_interp_order
         cache_slot = 1
         hi_loop_count = hi_loop_count + 1
         v%rho_interp_lo(hi_loop_count) = quadeos06(v%quad, recompute_flag, cache_slot, &
              slr, v%x_interp_result(t6_grid_idx,v%density_index_1), &
              v%x_interp_result(t6_grid_idx,v%density_index_2), &
              v%x_interp_result(t6_grid_idx,v%density_index_3), &
              v%density_grid(v%density_index_1), v%density_grid(v%density_index_2), &
              v%density_grid(v%density_index_3))
         if (v%density_interp_order.eq.3) then
            cache_slot = 2
            v%rho_interp_hi(hi_loop_count) = quadeos06(v%quad, recompute_flag, cache_slot, &
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
      v%esact = quadeos06(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(1), &
           v%rho_interp_lo(2), v%rho_interp_lo(3), v%t6_grid(v%t6_index_1), &
           v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      if (v%density_interp_order.eq.3) then
! .....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
         esactq = quadeos06(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_hi(1), &
              v%rho_interp_hi(2), v%rho_interp_hi(3), v%t6_grid(v%t6_index_1), &
              v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      end if
      if (v%t6_interp_order.eq.3) then
! .....    eos(i) in lower-left 3x3.
         esact2 = quadeos06(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(2), &
              v%rho_interp_lo(3), v%rho_interp_lo(4), v%t6_grid(v%t6_index_2), &
              v%t6_grid(v%t6_index_3), v%t6_grid(v%t6_index_4))
! .....    eos(i) smoothed in left 3x4
         dix = (v%t6_grid(v%t6_index_3) - slt)*v%t6_grid_spacing_inv(v%t6_index_3)
         v%esact = v%esact*dix + esact2*(1.0d0 - dix)
! endif   ! moved to loc a
         if (v%density_interp_order.eq.3) then

! .....     eos(i) in upper-right 3x3.
            esactq2 = quadeos06(v%quad, recompute_flag, cache_slot, slt, &
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
      if (v%esact.gt.1.0d+15) then
         write(run_log_unit,'("T6RINTEOS06: Interpolation indices out", &
              &" of range;please report conditions.")')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if

      return
end subroutine t6rinteos06
