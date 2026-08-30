!----------------------------------------------------------------------
! gtlaol2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original gtlaol2.f; only variable names, source form, and comment
! style were updated.
!
! DBG 5/94 Identical to gtlaol.f90 except gets opacity from the
! second (Z2) LAOL89 table.
subroutine gtlaol2(log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)

      use opacity_table_lib
      use luout_lib
      use numerics_lib
      use math_lib
      implicit none
      integer :: jerr_gate
      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      double precision :: row_log10_opacity(104), row_log_rho(104), &
           row_d2opacity(104), dlnkap_dlnrho_by_t(52), dlnkap_dlnrho_by_x(4)
      double precision :: logt_interp_opacity(52), logt_values(52), &
           logt_d2opacity(52), dlnkap_dlnt_by_x(4)
      double precision :: opacity_by_x(4), x_values(4)
      double precision :: log_extrap_tolerance, local_x, local_logt, &
           local_logrho
      integer :: t_locate_guess, t_index, x_grid_index, x_grid_index_hi, &
           num_valid_x, x_loop_index, num_valid_t, t_range_lo, t_range_hi, &
           num_valid_rho, rho_loop_index, spline_index_lo, spline_index_hi
      double precision :: log10_opacity_value, slope

!     1. CUBIC SPLINE INTERPOLATE IN DENSITY
!     2. CUBIC SPLINE INTERPOLATE IN TEMPERATURE
!     3. LINEAR INTERPOLATE IN COMPOSITION
!     IF WITHIN TOLLAOL OF EDGE THEN LINEAR EXTRAPOLATE
!
!     TOLLAOL PERMITS SOME EXTRAPLOATION BEYOND TABLE EDGE.
      ierr = 0
      log_extrap_tolerance = log(tollaol)
      local_x = hydrogen_fraction
      local_logt = log10_temperature
      local_logrho = log10_density
      call locate(opacity_table%ot2, opacity_table%nt2, local_logt, t_locate_guess)
      call locate(opacity_table%oxa2, opacity_table%nxyz2, local_x, x_grid_index)
      if (x_grid_index .eq. opacity_table%nxyz2) then
          x_grid_index = opacity_table%nxyz2-1
      end if
      if (x_grid_index .eq. 0) then
          x_grid_index = 1
      end if
      num_valid_x = 0
      x_grid_index_hi = x_grid_index + 1
      do x_loop_index=x_grid_index, x_grid_index_hi
         num_valid_t = 0
!        GET RANGE OF FOUR TT SURROUNDING T
         call xrng4(t_locate_guess, opacity_table%nt2, t_range_lo, t_range_hi)
         do t_index=t_range_lo, t_range_hi
            num_valid_rho = opacity_table%slaol2_num_points(x_loop_index,t_index)
            if (num_valid_rho .ge. 4) then
               do rho_loop_index=1, num_valid_rho
                  row_log10_opacity(rho_loop_index) = &
                       opacity_table%slaol2_opacity(x_loop_index,rho_loop_index,t_index)
                  row_log_rho(rho_loop_index) = &
                       opacity_table%slaol2_log_rho(x_loop_index,rho_loop_index,t_index)
                  row_d2opacity(rho_loop_index) = &
                       opacity_table%slaol2_d2opacity(x_loop_index,rho_loop_index,t_index)
               end do
               if (local_logrho.gt.row_log_rho(1) .and. &
                    local_logrho.lt.row_log_rho(num_valid_rho)) then
                  call splint(row_log_rho, row_log10_opacity, num_valid_rho, &
                       row_d2opacity, local_logrho, log10_opacity_value, &
                       spline_index_lo, spline_index_hi, jerr_gate)
                  ! 2026 numerics-gate opt-in: interpolation failure returns via
                  ! ierr (diagnostic printed at the gate) instead of stopping.
                  if (jerr_gate /= 0) then
                     ierr = jerr_gate
                     return
                  end if
                  num_valid_t = num_valid_t+1
                  logt_interp_opacity(num_valid_t) = log10_opacity_value
                  logt_values(num_valid_t) = opacity_table%ot2(t_index)
                  dlnkap_dlnrho_by_t(num_valid_t) = &
                       (row_log10_opacity(spline_index_hi)-row_log10_opacity(spline_index_lo))/ &
                       (row_log_rho(spline_index_hi)-row_log_rho(spline_index_lo))
               else if (local_logrho.gt.row_log_rho(1)-log_extrap_tolerance .and. &
                    local_logrho.le.row_log_rho(1)) then
!                 PUT LINEAR EXTRAPOLATION ROUTINES HERE
                  slope = (row_log10_opacity(2)-row_log10_opacity(1))/ &
                       (row_log_rho(2)-row_log_rho(1))
                  log10_opacity_value = row_log10_opacity(1)+slope*(local_logrho-row_log_rho(1))
                  num_valid_t = num_valid_t+1
                  logt_interp_opacity(num_valid_t) = log10_opacity_value
                  logt_values(num_valid_t) = opacity_table%ot2(t_index)
                  dlnkap_dlnrho_by_t(num_valid_t) = slope
               else if (local_logrho.ge.row_log_rho(num_valid_rho) .and. &
                    local_logrho.lt.row_log_rho(num_valid_rho)+log_extrap_tolerance) then
!                 PUT LINEAR EXTRAPOLATION ROUTINES HERE
                  slope = (row_log10_opacity(num_valid_rho-1)-row_log10_opacity(num_valid_rho))/ &
                       (row_log_rho(num_valid_rho-1)-row_log_rho(num_valid_rho))
                  log10_opacity_value = row_log10_opacity(num_valid_rho)+ &
                       slope*(local_logrho-row_log_rho(num_valid_rho))
                  num_valid_t = num_valid_t+1
                  logt_interp_opacity(num_valid_t) = log10_opacity_value
                  logt_values(num_valid_t) = opacity_table%ot2(t_index)
                  dlnkap_dlnrho_by_t(num_valid_t) = slope
               end if
            else
               write(run_log_unit,120) log10_density, log10_temperature
  120          format(' OUTSIDE OPACITY TABLE #2, IN DENSITY.  ', &
                    'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
               ierr = 1
               return
            end if
         end do
         if (num_valid_t .ge. 4) then
            call cspline(logt_values, logt_interp_opacity, num_valid_t, &
                 1.0d30, 1.0d30, logt_d2opacity)
            call splint(logt_values, logt_interp_opacity, num_valid_t, &
                 logt_d2opacity, local_logt, log10_opacity_value, &
                 spline_index_lo, spline_index_hi, jerr_gate)
            if (jerr_gate /= 0) then
               ierr = jerr_gate
               return
            end if
            num_valid_x = num_valid_x + 1
            opacity_by_x(num_valid_x) = log10_opacity_value
            slope=(dlnkap_dlnrho_by_t(spline_index_hi)-dlnkap_dlnrho_by_t(spline_index_lo))/ &
                 (logt_values(spline_index_hi)-logt_values(spline_index_lo))
            dlnkap_dlnrho_by_x(num_valid_x)=dlnkap_dlnrho_by_t(spline_index_lo)+ &
                 slope*(local_logt-logt_values(spline_index_lo))
            dlnkap_dlnt_by_x(num_valid_x) = &
                 (logt_interp_opacity(spline_index_hi)-logt_interp_opacity(spline_index_lo))/ &
                 (logt_values(spline_index_hi)-logt_values(spline_index_lo))
            x_values(num_valid_x) = opacity_table%oxa2(x_loop_index)
         else
            write(run_log_unit,121) log10_density, log10_temperature
  121       format(' OUTSIDE OPACITY TABLE #2, IN TEMPERATURE.  ', &
                 'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
            ierr = 1
            return
         end if
      end do
      if (num_valid_x .ge. 2) then
         slope = (opacity_by_x(2)-opacity_by_x(1))/(x_values(2)-x_values(1))
         log10_opacity_value = opacity_by_x(1)+slope*(hydrogen_fraction-x_values(1))
         slope = (dlnkap_dlnt_by_x(2)-dlnkap_dlnt_by_x(1))/(x_values(2)-x_values(1))
         dlnkap_dlnt = dlnkap_dlnt_by_x(1) + slope*(hydrogen_fraction-x_values(1))
         slope = (dlnkap_dlnrho_by_x(2)-dlnkap_dlnrho_by_x(1))/(x_values(2)-x_values(1))
         dlnkap_dlnrho = dlnkap_dlnrho_by_x(1)+slope*(hydrogen_fraction-x_values(1))
         if (log10_opacity_value .gt. 35) then
            opacity = 1.0d35
            log10_opacity = 35.0d0
         else
            opacity = exp10(log10_opacity_value)
            log10_opacity = log10_opacity_value
         end if
      else
         write(run_log_unit,122) log10_density, log10_temperature
  122    format(' OUTSIDE OPACITY TABLE #2.  ', &
              'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
          ierr = 1
          return
      end if
      return
end subroutine gtlaol2
