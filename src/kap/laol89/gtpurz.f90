!----------------------------------------------------------------------
! gtpurz
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original gtpurz.f; only variable names, source form, and comment
! style were updated.
!
! Get LAOL89 opacity for the pure-Z (Z=1) table: cubic-spline
! interpolate in density, then in temperature (no composition
! interpolation needed since the table is pure Z).
subroutine gtpurz(log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)

      use star_info_lib, only: star
      use opacity_table_lib
      use luout_lib
      use numerics_lib
      use math_lib
      implicit none
      integer :: jerr_gate
      double precision, intent(in) :: log10_density, log10_temperature
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      double precision :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho), dlnkap_dlnrho_by_t(n_laol_t)
      double precision :: logt_interp_opacity(n_laol_t), logt_values(n_laol_t), &
           logt_d2opacity(n_laol_t)
      double precision :: log_extrap_tolerance
      integer :: t_locate_guess, t_index, t_range_lo, t_range_hi, &
           num_valid_t, num_valid_rho, rho_loop_index, spline_index_lo, &
           spline_index_hi
      double precision :: log10_opacity_value, slope

!     1. CUBIC SPLINE INTERPOLATE IN DENSITY
!     2. CUBIC SPLINE INTERPOLATE IN TEMPERATURE
!     IF WITHIN TOLLAOL OF EDGE THEN LINEAR EXTRAPOLATE
!
!     TOLLAOL PERMITS SOME EXTRAPLOATION BEYOND TABLE EDGE.
      ierr = 0
      log_extrap_tolerance = log(star%ctrl%tollaol)
      call locate(opacity_table%zlaol_logt_grid, opacity_table%zlaol_num_t, log10_temperature, t_locate_guess)
      num_valid_t = 0
!     GET RANGE OF FOUR TT SURROUNDING T
      call xrng4(t_locate_guess, opacity_table%zlaol_num_t, t_range_lo, t_range_hi)
      do t_index=t_range_lo, t_range_hi
         num_valid_rho = opacity_table%zslaol_num_points(t_index)
         if (num_valid_rho .ge. 4) then
            do rho_loop_index=1, num_valid_rho
               row_log10_opacity(rho_loop_index) = opacity_table%zslaol_opacity(rho_loop_index,t_index)
               row_log_rho(rho_loop_index) = opacity_table%zslaol_log_rho(rho_loop_index,t_index)
               row_d2opacity(rho_loop_index) = opacity_table%zslaol_d2opacity(rho_loop_index,t_index)
            end do
            if (log10_density.gt.row_log_rho(1) .and. &
                 log10_density.lt.row_log_rho(num_valid_rho)) then
               call splint(row_log_rho, row_log10_opacity, num_valid_rho, &
                    row_d2opacity, log10_density, log10_opacity_value, &
                    spline_index_lo, spline_index_hi, jerr_gate)
               ! 2026 numerics-gate opt-in: interpolation failure returns via
               ! ierr (diagnostic printed at the gate) instead of stopping.
               if (jerr_gate /= 0) then
                  ierr = jerr_gate
                  return
               end if
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = opacity_table%zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = &
                    (row_log10_opacity(spline_index_hi)-row_log10_opacity(spline_index_lo))/ &
                    (row_log_rho(spline_index_hi)-row_log_rho(spline_index_lo))
            else if (log10_density.gt.row_log_rho(1)-log_extrap_tolerance .and. &
                 log10_density.le.row_log_rho(1)) then
!              PUT LINEAR EXTRAPOLATION ROUTINES HERE
               slope = (row_log10_opacity(2)-row_log10_opacity(1))/ &
                    (row_log_rho(2)-row_log_rho(1))
               log10_opacity_value = row_log10_opacity(1)+slope*(log10_density-row_log_rho(1))
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = opacity_table%zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = slope
            else if (log10_density.ge.row_log_rho(num_valid_rho) .and. &
                 log10_density.lt.row_log_rho(num_valid_rho)+log_extrap_tolerance) then
!              PUT LINEAR EXTRAPOLATION ROUTINES HERE
               slope = (row_log10_opacity(num_valid_rho-1)-row_log10_opacity(num_valid_rho))/ &
                    (row_log_rho(num_valid_rho-1)-row_log_rho(num_valid_rho))
               log10_opacity_value = row_log10_opacity(num_valid_rho)+ &
                    slope*(log10_density-row_log_rho(num_valid_rho))
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = opacity_table%zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = slope
            end if
         else
            write(run_log_unit,120) log10_density, log10_temperature
  120       format(' OUTSIDE Z OPACITY TABLE, IN DENSITY.  ', &
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
              logt_d2opacity, log10_temperature, log10_opacity_value, &
              spline_index_lo, spline_index_hi, jerr_gate)
         if (jerr_gate /= 0) then
            ierr = jerr_gate
            return
         end if
         slope=(dlnkap_dlnrho_by_t(spline_index_hi)-dlnkap_dlnrho_by_t(spline_index_lo))/ &
              (logt_values(spline_index_hi)-logt_values(spline_index_lo))
         dlnkap_dlnrho = dlnkap_dlnrho_by_t(spline_index_lo)+ &
              slope*(log10_temperature-logt_values(spline_index_lo))
         dlnkap_dlnt = (logt_interp_opacity(spline_index_hi)-logt_interp_opacity(spline_index_lo))/ &
              (logt_values(spline_index_hi)-logt_values(spline_index_lo))
      else
         write(run_log_unit,121) log10_density, log10_temperature
  121    format(' OUTSIDE Z OPACITY TABLE, IN TEMPERATURE.  ', &
              'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      end if
      if (log10_opacity_value .gt. 35) then
         opacity = laol_opacity_cap
         log10_opacity = 35.0d0
      else
         opacity = exp10(log10_opacity_value)
         log10_opacity = log10_opacity_value
      end if
      return
end subroutine gtpurz
