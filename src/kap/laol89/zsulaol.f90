!----------------------------------------------------------------------
! zsulaol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original zsulaol.f; only variable names, source form, and comment
! style were updated.
!
! DBG 4/94 Modified to do ZRAMP stuff.
! Builds the natural cubic spline coefficients (in log rho, at each
! tabulated T) for the pure-Z LAOL89 opacity table read by
! rdzlaol.f90, for use by gtpurz.f90.
subroutine zsulaol

      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none
      double precision :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho)
      integer :: it, ir, num_valid_rho

      do it=1, opacity_table%zlaol_num_t
         opacity_table%zlaol_logt_grid(it) = log10(opacity_table%zlaol_logt_grid(it))
      end do
        do it=1, opacity_table%zlaol_num_t
            num_valid_rho=0
            do ir=1, opacity_table%zlaol_num_rho
                  opacity_table%zslaol_opacity(ir,it) = 0.0d0
                  opacity_table%zslaol_log_rho(ir,it) = 0.0d0
                  opacity_table%zslaol_d2opacity(ir,it) = 0.0d0
                  if (opacity_table%zlaol_opacity(ir,it) .ne. 0.0d0) then
                     num_valid_rho = num_valid_rho+1
                     row_log10_opacity(num_valid_rho) = log10(opacity_table%zlaol_opacity(ir,it))
                     row_log_rho(num_valid_rho) = log10(opacity_table%zlaol_logrho_grid(ir))
                  end if
            end do
            if (num_valid_rho .ge. 4) then
               opacity_table%zslaol_num_points(it)=num_valid_rho
               call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    1.0d30, 1.0d30, row_d2opacity)
               do ir=1,num_valid_rho
                     opacity_table%zslaol_opacity(ir,it) = row_log10_opacity(ir)
                     opacity_table%zslaol_log_rho(ir,it) = row_log_rho(ir)
                     opacity_table%zslaol_d2opacity(ir,it) = row_d2opacity(ir)
               end do
            else
               opacity_table%zslaol_num_points(it) = 0
            end if
        end do
      return
end subroutine zsulaol
