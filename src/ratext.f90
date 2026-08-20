!----------------------------------------------------------------------
! ratext
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ratext.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Diagonal rational-function extrapolator used to extrapolate a
! sequence of sub-stepped burning-rate estimates (indexed by
! decreasing step size) to the zero-step-size limit. Same algorithm as
! subroutine RZEXTR in Numerical Recipes, p.566. Called repeatedly by
! liburn/liburn2 with increasing est_index as the lithium/beryllium
! burning sub-stepping is refined.
subroutine ratext(est_index, x_est, y_est, y_extrap, y_err, num_vars, &
     max_use)
      implicit none
      integer, parameter :: imax=11, nmax=15, ncol=7

      integer, intent(in) :: est_index, num_vars, max_use
      double precision, intent(in) :: x_est
      double precision, intent(in) :: y_est(num_vars)
      double precision, intent(out) :: y_extrap(num_vars), y_err(num_vars)

      double precision :: x_hist(imax), tableau(nmax,ncol), fx(ncol)
      save

! SAME AS SR RZEXTR FROM NUMERICAL RECIPES, P.566.
!
! SAVE CURRENT INDEPENDENT VARIABLE.
      integer :: var_idx, k, k2, num_use
      double precision :: yy, v, c, b1, b, ddy

      x_hist(est_index) = x_est
      if (est_index.eq.1) then
         do var_idx = 1, num_vars
            y_extrap(var_idx) = y_est(var_idx)
            tableau(var_idx,1) = y_est(var_idx)
            y_err(var_idx) = y_est(var_idx)
         end do
      else
!        USE AT MOST max_use PREVIOUS MEMBERS.
         num_use = min(est_index,max_use)
         do k = 1, num_use-1
            fx(k+1) = x_hist(est_index-k)/x_est
         end do
!        EVALUATE NEXT DIAGONAL IN TABLEAU.
         do var_idx=1,num_vars
            yy = y_est(var_idx)
            v = tableau(var_idx,1)
            c = yy
            tableau(var_idx,1) = yy
            do k2 = 2, num_use
               b1 = fx(k2)*v
               b = b1 - c
!              CARE NEEDED TO AVOID DIVISION BY ZERO.
               if (b.ne.0d0) then
                  b = (c - v)/b
                  ddy = c*b
                  c = b1*b
               else
                  ddy = v
               end if
               v = tableau(var_idx,k2)
               tableau(var_idx,k2) = ddy
               yy = yy + ddy
            end do
            y_err(var_idx) = ddy
            y_extrap(var_idx) = yy
         end do
      end if
      return
end subroutine ratext
