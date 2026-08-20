!------------------------    GROUP: SR_P     -------------------------------
!----------------------------------------------------------------------
! intpt
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intpt.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). NOT exercised by that suite
! (only called from mhdpx2.f90, the MHD equation-of-state path, which
! the suite's 4 test cases do not select) -- verified by build +
! code review only.
!
! Bicubic-style interpolation in a (log10_pressure, log10_temperature)
! table_data(table_dim_t,table_dim_r,num_vars) grid, via two passes of
! the external 4-point Lagrangian interpolator lir: first at fixed
! t_indices columns to interpolate in pressure (work1 -> work2), then
! across the 4 selected t_indices to interpolate in temperature
! (work2 -> interp_vars). Dummy-argument names match the actual
! arguments used at the intpt call sites in mhdpx2.f90.
subroutine intpt(log10_pressure, log10_temperature, table_data, &
     table_dim_t, table_dim_r, num_vars, table_log10t, num_t, num_r, &
     work1, work2, y_work, interp_vars)
      implicit none
      integer, intent(in) :: table_dim_t, table_dim_r, num_vars
      double precision, intent(in) :: log10_pressure, log10_temperature
      double precision, intent(in) :: table_data(table_dim_t,table_dim_r,num_vars)
      double precision, intent(in) :: table_log10t(table_dim_t)
      integer, intent(in) :: num_t, num_r
      double precision, intent(inout) :: work1(num_vars,4), work2(num_vars,4)
      double precision, intent(inout) :: y_work(num_vars)
      double precision, intent(out) :: interp_vars(num_vars)

      integer :: r_lo_guess(4), r_indices(4,4), t_indices(4)
      double precision :: x_nodes(4)
! common/ccout2/: no member is used anywhere in this batch of files;
! all are unrenamed placeholders preserving the storage layout.
! Naming matches meqos.f90.
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch
! lir_order is INTEGER*4 in the original (overriding this file's
! IMPLICIT LOGICAL*4(L) for the single name L), holding a flag passed
! to the external routine lir; its exact meaning there is not
! established from this file alone.
      integer :: lir_order
      save

      integer :: n, i, m, j, t_col, iv, t_idx, r_idx, t_idx_max, r_idx_max
      double precision :: p_min, p_max
      integer :: lir_num_vars, lir_leading_dim, lir_num_points, lir_interp_mode

      do 100 n=1,num_t
         if(table_log10t(n).ge.log10_temperature) goto 101
         t_indices(1)=n
 100  continue
 101  if(t_indices(1).ge.2) t_indices(1)=t_indices(1)-1
      t_idx_max=num_t-3
      if(t_indices(1).gt.t_idx_max) t_indices(1)=t_idx_max
      do i=2,4
         t_indices(i)=t_indices(1)+i-1
      end do
      do i=1,4
         r_lo_guess(i)=1
         t_idx=t_indices(i)
         p_min=table_data(t_idx, 1,2)
         p_max=table_data(t_idx,num_r,2)
         if(log10_pressure.gt.p_max) then
            return
         end if
         do 200 m=1,num_r
            if(table_data(t_idx,m,2).ge.log10_pressure) goto 201
            r_lo_guess(i)=m
 200     continue
 201     if(r_lo_guess(i).ge.2) r_lo_guess(i)=r_lo_guess(i)-1
         r_idx_max=num_r-3
         if(r_lo_guess(i).gt.r_idx_max) r_lo_guess(i)=r_idx_max
      end do
      do i=1,4
         do j=1,4
            r_indices(j,i)=r_lo_guess(i)+j-1
         end do
      end do
      do t_col=1,4
         t_idx=t_indices(t_col)
         do i=1,4
            r_idx=r_indices(i,t_col)
            x_nodes(i)=table_data(t_idx,r_idx,2)
         end do
         do i=1,4
            r_idx=r_indices(i,t_col)
            do iv=1,num_vars
               work1(iv,i)=table_data(t_idx,r_idx,iv)
            end do
         end do

         lir_num_vars=num_vars
         lir_leading_dim=num_vars
         lir_num_points=4
         lir_order=1
         lir_interp_mode=1
         call lir(log10_pressure, x_nodes, y_work, work1, lir_num_vars, &
              lir_leading_dim, lir_num_points, lir_order, lir_interp_mode)
         do iv=1,num_vars
            work2(iv,t_col) = y_work(iv)
         end do

      end do
      do i=1,4
         t_idx=t_indices(i)
         x_nodes(i) = table_log10t(t_idx)
      end do

      lir_num_vars=num_vars
      lir_leading_dim=num_vars
      lir_num_points=4
      lir_order=1
      lir_interp_mode=1

      call lir(log10_temperature, x_nodes, interp_vars, work2, lir_num_vars, &
           lir_leading_dim, lir_num_points, lir_order, lir_interp_mode)

      return
end subroutine intpt
