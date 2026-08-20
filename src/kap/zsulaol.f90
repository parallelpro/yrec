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

      implicit none
      double precision :: row_log10_opacity(104), row_log_rho(104), &
           row_d2opacity(104)

! DBG 12/95 ARRAYS FOR PURE Z TABLE
      double precision :: zlaol_opacity(104,52), zlaol_logt_grid(52), &
           zlaol_logrho_grid(104)
      integer :: zlaol_num_rho, zlaol_num_t
      common/zlaol/ zlaol_opacity, zlaol_logt_grid, zlaol_logrho_grid, &
           zlaol_num_rho, zlaol_num_t

      double precision :: zslaol_opacity(104,52), zslaol_log_rho(104,52), &
           zslaol_d2opacity(104,52)
      integer :: zslaol_num_points(52)
      common/zslaol/ zslaol_opacity, zslaol_log_rho, zslaol_d2opacity, &
           zslaol_num_points

      save

      integer :: it, ir, num_valid_rho

      do it=1, zlaol_num_t
         zlaol_logt_grid(it) = log10(zlaol_logt_grid(it))
      end do
        do it=1, zlaol_num_t
            num_valid_rho=0
            do ir=1, zlaol_num_rho
                  zslaol_opacity(ir,it) = 0.0d0
                  zslaol_log_rho(ir,it) = 0.0d0
                  zslaol_d2opacity(ir,it) = 0.0d0
                  if (zlaol_opacity(ir,it) .ne. 0.0d0) then
                     num_valid_rho = num_valid_rho+1
                     row_log10_opacity(num_valid_rho) = log10(zlaol_opacity(ir,it))
                     row_log_rho(num_valid_rho) = log10(zlaol_logrho_grid(ir))
                  end if
            end do
            if (num_valid_rho .ge. 4) then
               zslaol_num_points(it)=num_valid_rho
               call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    1.0d30, 1.0d30, row_d2opacity)
               do ir=1,num_valid_rho
                     zslaol_opacity(ir,it) = row_log10_opacity(ir)
                     zslaol_log_rho(ir,it) = row_log_rho(ir)
                     zslaol_d2opacity(ir,it) = row_d2opacity(ir)
               end do
            else
               zslaol_num_points(it) = 0
            end if
        end do
      return
end subroutine zsulaol
