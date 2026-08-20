!----------------------------------------------------------------------
! blkdta000
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original blkdta000.f; only variable names, source form, and comment
! style were updated.
!
! BLOCK DATA unit initializing the OPAL 1995 EOS interpolation grids
! (common/a/, common/b/) at program startup: the fixed X-composition
! grid (x_grid), the T6-index lower bound per density row
! (t6_index_lo), and the column-order index (eos_var_order). Member
! names for common/aa/, common/a/, common/b/ are reused verbatim from
! the already-converted esac.f90/readco.f90/radsub.f90/rhoofp.f90/
! t6rinterp.f90 (COMMON storage is positional, not name-matched, but
! this file shares the same underlying storage and so must match
! layout: order, type, and dimensions).
block data
      implicit none
      integer, parameter :: mx=5, mv=10, nr=77, nt=56

! common/aa/: not used in this file; placeholders (see t6rinterp.f90).
      double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
      common/aa/ rho_interp_hi, rho_interp_lo, xxh

! common/a/: only x_grid is initialized here; the rest are
! placeholders (see esac.f90/readco.f90/t6rinterp.f90 for the full
! description of this block's layout).
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr), &
           x_grid(mx)
      integer :: x_loop_index, x_index_lo
      common/a/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_grid, &
           x_loop_index, x_index_lo

! common/b/: t6_index_lo and eos_var_order are initialized here;
! z_table/eos_index_inverse are placeholders (see readco.f90).
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), t6_index_lo(nr)
      common/b/ z_table, eos_index_inverse, eos_var_order, t6_index_lo

      integer :: i

      data (x_grid(i),i=1,mx)/0.0d0,0.2d0,0.4d0,0.6d0,0.8d0/
      data (t6_index_lo(i),i=1,nr)/37*56,7*54, &
     &39,37,36,34,33,31,31,30,29,28,27,26,26,25,24,23,22,21,21,20, &
     &19,18,17,16,15,15,13,13,11,11,9,9,7/
      data (eos_var_order(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/
      save
end block data
