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
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use numerics_lib
      implicit none
      double precision, intent(in) :: log10_density, log10_temperature
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

! MHP 8/25 Removed unused variables
      double precision :: row_log10_opacity(104), row_log_rho(104), &
           row_d2opacity(104), dlnkap_dlnrho_by_t(52)
      double precision :: logt_interp_opacity(52), logt_values(52), &
           logt_d2opacity(52)

! common/luout/: only short_file_unit is used here.
      integer :: ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr

! MHP 8/25 Removed character file names from common block
! common/nwlaol/: only tollaol is used here.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez

! common/slaol/: not used here; declared only to preserve the shared
! storage layout (see gtlaol.f90 for these names).
      double precision :: slaol_opacity(12,104,52), slaol_log_rho(12,104,52), &
           slaol_d2opacity(12,104,52)
      integer :: slaol_num_points(12,52)
      common/slaol/ slaol_opacity, slaol_log_rho, slaol_d2opacity, &
           slaol_num_points

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

      double precision :: log_extrap_tolerance, local_logt, local_logrho
      integer :: t_locate_guess, t_index, t_range_lo, t_range_hi, &
           num_valid_t, num_valid_rho, rho_loop_index, spline_index_lo, &
           spline_index_hi
      double precision :: log10_opacity_value, slope

!     1. CUBIC SPLINE INTERPOLATE IN DENSITY
!     2. CUBIC SPLINE INTERPOLATE IN TEMPERATURE
!     IF WITHIN TOLLAOL OF EDGE THEN LINEAR EXTRAPOLATE
!
!     TOLLAOL PERMITS SOME EXTRAPLOATION BEYOND TABLE EDGE.
      log_extrap_tolerance = log(tollaol)
      local_logt = log10_temperature
      local_logrho = log10_density
      call locate(zlaol_logt_grid, zlaol_num_t, local_logt, t_locate_guess)
      num_valid_t = 0
!     GET RANGE OF FOUR TT SURROUNDING T
      call xrng4(t_locate_guess, zlaol_num_t, t_range_lo, t_range_hi)
      do t_index=t_range_lo, t_range_hi
         num_valid_rho = zslaol_num_points(t_index)
         if (num_valid_rho .ge. 4) then
            do rho_loop_index=1, num_valid_rho
               row_log10_opacity(rho_loop_index) = zslaol_opacity(rho_loop_index,t_index)
               row_log_rho(rho_loop_index) = zslaol_log_rho(rho_loop_index,t_index)
               row_d2opacity(rho_loop_index) = zslaol_d2opacity(rho_loop_index,t_index)
            end do
            if (local_logrho.gt.row_log_rho(1) .and. &
                 local_logrho.lt.row_log_rho(num_valid_rho)) then
               call splint(row_log_rho, row_log10_opacity, num_valid_rho, &
                    row_d2opacity, local_logrho, log10_opacity_value, &
                    spline_index_lo, spline_index_hi)
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = &
                    (row_log10_opacity(spline_index_hi)-row_log10_opacity(spline_index_lo))/ &
                    (row_log_rho(spline_index_hi)-row_log_rho(spline_index_lo))
            else if (local_logrho.gt.row_log_rho(1)-log_extrap_tolerance .and. &
                 local_logrho.le.row_log_rho(1)) then
!              PUT LINEAR EXTRAPOLATION ROUTINES HERE
               slope = (row_log10_opacity(2)-row_log10_opacity(1))/ &
                    (row_log_rho(2)-row_log_rho(1))
               log10_opacity_value = row_log10_opacity(1)+slope*(local_logrho-row_log_rho(1))
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = slope
            else if (local_logrho.ge.row_log_rho(num_valid_rho) .and. &
                 local_logrho.lt.row_log_rho(num_valid_rho)+log_extrap_tolerance) then
!              PUT LINEAR EXTRAPOLATION ROUTINES HERE
               slope = (row_log10_opacity(num_valid_rho-1)-row_log10_opacity(num_valid_rho))/ &
                    (row_log_rho(num_valid_rho-1)-row_log_rho(num_valid_rho))
               log10_opacity_value = row_log10_opacity(num_valid_rho)+ &
                    slope*(local_logrho-row_log_rho(num_valid_rho))
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = slope
            end if
         else
            write(short_file_unit,120) log10_density, log10_temperature
  120       format(' OUTSIDE Z OPACITY TABLE, IN DENSITY.  ', &
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
         slope=(dlnkap_dlnrho_by_t(spline_index_hi)-dlnkap_dlnrho_by_t(spline_index_lo))/ &
              (logt_values(spline_index_hi)-logt_values(spline_index_lo))
         dlnkap_dlnrho = dlnkap_dlnrho_by_t(spline_index_lo)+ &
              slope*(local_logt-logt_values(spline_index_lo))
         dlnkap_dlnt = (logt_interp_opacity(spline_index_hi)-logt_interp_opacity(spline_index_lo))/ &
              (logt_values(spline_index_hi)-logt_values(spline_index_lo))
      else
         write(short_file_unit,121) log10_density, log10_temperature
  121    format(' OUTSIDE Z OPACITY TABLE, IN TEMPERATURE.  ', &
              'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
         stop
      end if
      if (log10_opacity_value .gt. 35) then
         opacity = 1.0d35
         log10_opacity = 35.0d0
      else
         opacity = 10.0d0**log10_opacity_value
         log10_opacity = log10_opacity_value
      end if
      return
end subroutine gtpurz
