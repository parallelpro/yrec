!----------------------------------------------------------------------
! gtlaol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original gtlaol.f; only variable names, source form, and comment
! style were updated.
!
! Get LAOL89 opacity: cubic-spline interpolate in density, then in
! temperature, then linearly interpolate (or, within tollaol of the
! table edge, linearly extrapolate) in hydrogen abundance.
subroutine gtlaol(log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use luout_lib
      use numerics_lib
      implicit none
      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

! MHP 8/25 Removed character file names from common block
      double precision :: row_log10_opacity(104), row_log_rho(104), &
           row_d2opacity(104), dlnkap_dlnrho_by_t(52), dlnkap_dlnrho_by_x(4)
      double precision :: logt_interp_opacity(52), logt_values(52), &
           logt_d2opacity(52), dlnkap_dlnt_by_x(4)
      double precision :: opacity_by_x(4), x_values(4)


! MHP 8/25 Removed character file names from common block
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez

      double precision :: slaol_opacity(12,104,52), slaol_log_rho(12,104,52), &
           slaol_d2opacity(12,104,52)
      integer :: slaol_num_points(12,52)
      common/slaol/ slaol_opacity, slaol_log_rho, slaol_d2opacity, &
           slaol_num_points

      save

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
      log_extrap_tolerance = log(tollaol)
      local_x = hydrogen_fraction
      local_logt = log10_temperature
      local_logrho = log10_density
      call locate(ot, numt, local_logt, t_locate_guess)
      call locate(oxa, numofxyz, local_x, x_grid_index)
      if (x_grid_index .eq. numofxyz) then
          x_grid_index = numofxyz-1
      end if
      if (x_grid_index .eq. 0) then
          x_grid_index = 1
      end if
      num_valid_x = 0
      x_grid_index_hi = x_grid_index + 1
      do x_loop_index=x_grid_index, x_grid_index_hi
         num_valid_t = 0
!        GET RANGE OF FOUR TT SURROUNDING T
         call xrng4(t_locate_guess, numt, t_range_lo, t_range_hi)
         do t_index=t_range_lo, t_range_hi
            num_valid_rho = slaol_num_points(x_loop_index,t_index)
            if (num_valid_rho .ge. 4) then
               do rho_loop_index=1, num_valid_rho
                  row_log10_opacity(rho_loop_index) = &
                       slaol_opacity(x_loop_index,rho_loop_index,t_index)
                  row_log_rho(rho_loop_index) = &
                       slaol_log_rho(x_loop_index,rho_loop_index,t_index)
                  row_d2opacity(rho_loop_index) = &
                       slaol_d2opacity(x_loop_index,rho_loop_index,t_index)
               end do
               if (local_logrho.gt.row_log_rho(1) .and. &
                    local_logrho.lt.row_log_rho(num_valid_rho)) then
                  call splint(row_log_rho, row_log10_opacity, num_valid_rho, &
                       row_d2opacity, local_logrho, log10_opacity_value, &
                       spline_index_lo, spline_index_hi)
                  num_valid_t = num_valid_t+1
                  logt_interp_opacity(num_valid_t) = log10_opacity_value
                  logt_values(num_valid_t) = ot(t_index)
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
                  logt_values(num_valid_t) = ot(t_index)
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
                  logt_values(num_valid_t) = ot(t_index)
                  dlnkap_dlnrho_by_t(num_valid_t) = slope
               end if
            else
               write(short_file_unit,120) log10_density, log10_temperature
  120          format(' OUTSIDE OPACITY TABLE, IN DENSITY.  ', &
                    'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
               stop
            end if
         end do
         if (num_valid_t .ge. 4) then
            call cspline(logt_values, logt_interp_opacity, num_valid_t, &
                 1.0d30, 1.0d30, logt_d2opacity)
            call splint(logt_values, logt_interp_opacity, num_valid_t, &
                 logt_d2opacity, local_logt, log10_opacity_value, &
                 spline_index_lo, spline_index_hi)
            num_valid_x = num_valid_x + 1
            opacity_by_x(num_valid_x) = log10_opacity_value
            slope=(dlnkap_dlnrho_by_t(spline_index_hi)-dlnkap_dlnrho_by_t(spline_index_lo))/ &
                 (logt_values(spline_index_hi)-logt_values(spline_index_lo))
            dlnkap_dlnrho_by_x(num_valid_x)=dlnkap_dlnrho_by_t(spline_index_lo)+ &
                 slope*(local_logt-logt_values(spline_index_lo))
            dlnkap_dlnt_by_x(num_valid_x) = &
                 (logt_interp_opacity(spline_index_hi)-logt_interp_opacity(spline_index_lo))/ &
                 (logt_values(spline_index_hi)-logt_values(spline_index_lo))
            x_values(num_valid_x) = oxa(x_loop_index)
         else
            write(short_file_unit,121) log10_density, log10_temperature
  121       format(' OUTSIDE OPACITY TABLE, IN TEMPERATURE.  ', &
                 'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
            stop
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
            opacity = 10.0d0**log10_opacity_value
            log10_opacity = log10_opacity_value
         end if
      else
         write(short_file_unit,122) log10_density, log10_temperature
  122    format(' OUTSIDE OPACITY TABLE.  ', &
              'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
          stop
      end if
      return
end subroutine gtlaol
