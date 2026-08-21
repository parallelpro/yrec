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
! (x_interp_result) in T6 and density to produce esact.
subroutine t6rinteos01(slr, slt)

      use luout_lib
      implicit none

      double precision, intent(in) :: slr, slt

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 191

! common/eeeos/: not used in this file; placeholders (see esac01.f90).
      double precision :: x_interp_workspace(mx,nt,nr), x_grid_copy(mx)
      common/eeeos/ x_interp_workspace, x_grid_copy

! common/aaeos/: rho_interp_hi/rho_interp_lo are used below; xxh is a
! placeholder.
      double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
      common/aaeos/ rho_interp_hi, rho_interp_lo, xxh

! common/aeos/: see esac01.f90 for the full description.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr)
      integer :: x_loop_index, x_index_lo
      double precision :: x_grid(mx)
      common/aeos/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_loop_index, &
           x_index_lo, x_grid

! common/bbeos/: density index window (l1..l4), t6 index window
! (k1..k4), and the interpolation-order flags t6_interp_order
! (original ip) and density_interp_order (original iq).
      integer :: density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order
      common/bbeos/ density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order

! common/eeos/: esact, eos_output(mv).
      double precision :: esact, eos_output(mv)
      common/eeos/ esact, eos_output


      save

! --- locals ---
      integer :: hi_loop_count, recompute_flag, cache_slot, t6_grid_idx
      double precision :: esactq, esact2, esactq2, dix, dix2
      double precision, external :: quadeos01

      hi_loop_count = 0
      recompute_flag = 0

      do t6_grid_idx = t6_index_1, t6_index_1 + t6_interp_order
         cache_slot = 1
         hi_loop_count = hi_loop_count + 1
         rho_interp_lo(hi_loop_count) = quadeos01(recompute_flag, cache_slot, &
              slr, x_interp_result(t6_grid_idx,density_index_1), &
              x_interp_result(t6_grid_idx,density_index_2), &
              x_interp_result(t6_grid_idx,density_index_3), &
              density_grid(density_index_1), density_grid(density_index_2), &
              density_grid(density_index_3))
         if (density_interp_order.eq.3) then
            cache_slot = 2
            rho_interp_hi(hi_loop_count) = quadeos01(recompute_flag, cache_slot, &
                 slr, x_interp_result(t6_grid_idx,density_index_2), &
                 x_interp_result(t6_grid_idx,density_index_3), &
                 x_interp_result(t6_grid_idx,density_index_4), &
                 density_grid(density_index_2), density_grid(density_index_3), &
                 density_grid(density_index_4))
         end if
         recompute_flag = 1
      end do

      recompute_flag = 0
      cache_slot = 1
! ..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      esact = quadeos01(recompute_flag, cache_slot, slt, rho_interp_lo(1), &
           rho_interp_lo(2), rho_interp_lo(3), t6_grid(t6_index_1), &
           t6_grid(t6_index_2), t6_grid(t6_index_3))
      if (density_interp_order.eq.3) then
! .....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
         esactq = quadeos01(recompute_flag, cache_slot, slt, rho_interp_hi(1), &
              rho_interp_hi(2), rho_interp_hi(3), t6_grid(t6_index_1), &
              t6_grid(t6_index_2), t6_grid(t6_index_3))
      end if
      if (t6_interp_order.eq.3) then
! .....    eos(i) in lower-left 3x3.
         esact2 = quadeos01(recompute_flag, cache_slot, slt, rho_interp_lo(2), &
              rho_interp_lo(3), rho_interp_lo(4), t6_grid(t6_index_2), &
              t6_grid(t6_index_3), t6_grid(t6_index_4))
! .....    eos(i) smoothed in left 3x4
         dix = (t6_grid(t6_index_3) - slt)*t6_grid_spacing_inv(t6_index_3)
         esact = esact*dix + esact2*(1.0d0 - dix)
! endif   ! moved to loc a
         if (density_interp_order.eq.3) then

! .....     eos(i) in upper-right 3x3.
            esactq2 = quadeos01(recompute_flag, cache_slot, slt, &
                 rho_interp_hi(2), rho_interp_hi(3), rho_interp_hi(4), &
                 t6_grid(t6_index_2), t6_grid(t6_index_3), &
                 t6_grid(t6_index_4))
            esactq = esactq*dix + esactq2*(1.0d0 - dix)
         end if
      end if  ! loc a
!
      if (density_interp_order.eq.3) then
         dix2 = (density_grid(density_index_3) - slr)* &
              density_grid_spacing_inv(density_index_3)
         if (t6_interp_order.eq.3) then
! .....        eos(i) smoothed in both log(T6) and log(R)
            esact = esact*dix2 + esactq*(1.0d0 - dix2)
         end if
      end if
      if (esact.gt.1.0d+15) then
         write(short_file_unit,'("Interpolation indices out of range", &
              &";please report conditions.")')
         stop
      end if

      return
end subroutine t6rinteos01
