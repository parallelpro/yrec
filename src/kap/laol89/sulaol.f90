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
! rdlaol.f90, for use by gtlaol.f90/gtlaol2.f90.
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

      do it=1, opacity_table%laol_num_t
         opacity_table%laol_grid_t(it) = log10(opacity_table%laol_grid_t(it))
      end do
      do ix=1, opacity_table%laol_num_x
         do it=1, opacity_table%laol_num_t
            num_valid_rho=0
            do ir=1, opacity_table%laol_num_rho
                opacity_table%slaol_opacity(ix,ir,it) = 0.0d0
                opacity_table%slaol_log_rho(ix,ir,it) = 0.0d0
                opacity_table%slaol_d2opacity(ix,ir,it) = 0.0d0
                if (opacity_table%laol_opacity(ix,ir,it) .ne. 0.0d0) then
                   num_valid_rho = num_valid_rho+1
                   row_log10_opacity(num_valid_rho) = log10(opacity_table%laol_opacity(ix,ir,it))
                   row_log_rho(num_valid_rho) = log10(opacity_table%laol_grid_rho(ir))
                end if
            end do
            if (num_valid_rho .ge. 4) then
               opacity_table%slaol_num_points(ix,it)=num_valid_rho
               call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    1.0d30, 1.0d30, row_d2opacity)
               do ir=1,num_valid_rho
                   opacity_table%slaol_opacity(ix,ir,it) = row_log10_opacity(ir)
                   opacity_table%slaol_log_rho(ix,ir,it) = row_log_rho(ir)
                   opacity_table%slaol_d2opacity(ix,ir,it) = row_d2opacity(ir)
               end do
            else
               opacity_table%slaol_num_points(ix,it) = 0
            end if
         end do
      end do
! DBG 4/94 Do SPLINE on second opacity table if ZRAMP
      if (star%use_two_z_tables) then
       do it=1, opacity_table%laol_num_t
          opacity_table%laol2_grid_t(it) = log10(opacity_table%laol2_grid_t(it))
       end do
         do ix=1, opacity_table%laol2_num_x
            do it=1, opacity_table%laol2_num_t
               num_valid_rho=0
               do ir=1, opacity_table%laol2_num_rho
                  opacity_table%slaol2_opacity(ix,ir,it) = 0.0d0
                  opacity_table%slaol2_log_rho(ix,ir,it) = 0.0d0
                  opacity_table%slaol2_d2opacity(ix,ir,it) = 0.0d0
                  if (opacity_table%laol2_opacity(ix,ir,it) .ne. 0.0d0) then
                     num_valid_rho = num_valid_rho+1
                     row_log10_opacity(num_valid_rho) = log10(opacity_table%laol2_opacity(ix,ir,it))
                     row_log_rho(num_valid_rho) = log10(opacity_table%laol2_grid_rho(ir))
                  end if
               end do
               if (num_valid_rho .ge. 4) then
                  opacity_table%slaol2_num_points(ix,it)=num_valid_rho
                  call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                       1.0d30, 1.0d30, row_d2opacity)
                  do ir=1,num_valid_rho
                     opacity_table%slaol2_opacity(ix,ir,it) = row_log10_opacity(ir)
                     opacity_table%slaol2_log_rho(ix,ir,it) = row_log_rho(ir)
                     opacity_table%slaol2_d2opacity(ix,ir,it) = row_d2opacity(ir)
                  end do
               else
                  opacity_table%slaol2_num_points(ix,it) = 0
               end if
            end do
         end do
      end if
      return
end subroutine sulaol
