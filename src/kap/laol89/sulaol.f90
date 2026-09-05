!----------------------------------------------------------------------
! sulaol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original sulaol.f; only variable names, source form, and comment
! style were updated.
!
! DBG 4/94 Modified to do ZRAMP stuff.
! Builds the natural cubic spline coefficients (in log rho, at each
! tabulated X and T) for the LAOL89 opacity table(s) read by
! rdlaol.f90, for use by gtlaol.f90.
!
! 2026 wave 3 (R5): the spline loop, identical for both tables, is
! the internal spline_laol_table below; the log10 conversion of each
! table's T grid stays here because the second table's loop bound
! was (and is) the first table's num_t.
subroutine sulaol
      use star_info_lib, only: star

      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none

! MHP 8/25 Removed unused variables
      double precision :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho)
      integer :: it, ix, ir, num_valid_rho

      do it=1, opacity_table%laol(1)%num_t
         opacity_table%laol(1)%grid_t(it) = log10(opacity_table%laol(1)%grid_t(it))
      end do
      call spline_laol_table(opacity_table%laol(1))
! DBG 4/94 Do SPLINE on second opacity table if ZRAMP
      if (star%use_two_z_tables) then
! NOTE: the log10 conversion of the second table's T grid runs over
! the FIRST table's num_t, as the original did (an R6 decision
! whether that is a latent bug when the two tables differ in T
! extent); it is therefore kept here rather than in the helper.
       do it=1, opacity_table%laol(1)%num_t
          opacity_table%laol(2)%grid_t(it) = log10(opacity_table%laol(2)%grid_t(it))
       end do
         call spline_laol_table(opacity_table%laol(2))
      end if
      return

contains

! spline_laol_table: the natural cubic spline in log rho of every
! (X, T) row of tbl -- token-identical in both former halves.
subroutine spline_laol_table(tbl)
      type(laol_table_set), intent(inout) :: tbl

      do ix=1, tbl%num_x
         do it=1, tbl%num_t
            num_valid_rho=0
            do ir=1, tbl%num_rho
                tbl%slaol_opacity(ix,ir,it) = 0.0d0
                tbl%slaol_log_rho(ix,ir,it) = 0.0d0
                tbl%slaol_d2opacity(ix,ir,it) = 0.0d0
                if (tbl%opacity(ix,ir,it) .ne. 0.0d0) then
                   num_valid_rho = num_valid_rho+1
                   row_log10_opacity(num_valid_rho) = log10(tbl%opacity(ix,ir,it))
                   row_log_rho(num_valid_rho) = log10(tbl%grid_rho(ir))
                end if
            end do
            if (num_valid_rho .ge. 4) then
               tbl%slaol_num_points(ix,it)=num_valid_rho
               call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    1.0d30, 1.0d30, row_d2opacity)
               do ir=1,num_valid_rho
                   tbl%slaol_opacity(ix,ir,it) = row_log10_opacity(ir)
                   tbl%slaol_log_rho(ix,ir,it) = row_log_rho(ir)
                   tbl%slaol_d2opacity(ix,ir,it) = row_d2opacity(ir)
               end do
            else
               tbl%slaol_num_points(ix,it) = 0
            end if
         end do
      end do
end subroutine spline_laol_table

end subroutine sulaol
