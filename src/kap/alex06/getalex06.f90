!----------------------------------------------------------------------
! getalex06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getalex06.f; only variable names, source form, and comment
! style were updated.
!
! MHP 3/09 Computes Alexander 2006 (Ferguson et al. 2005) low-
! temperature opacities by 4-point Lagrangian interpolation in log T
! and log R = rho/T6**3, reloading the fixed-(X,Z) table via
! alex06tab.f90 whenever the requested composition changes.
subroutine getalex06(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use opacity_table_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: num_x = 9
      integer, parameter :: num_z = 16
      integer, parameter :: num_xz = 143
      integer, parameter :: num_t = 85
      integer, parameter :: num_d = 19

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt




      double precision :: interp_nodes(4), weight_t(4), dweight_t(4), &
           weight_r(4), dweight_r(4), opacity_row(4), dlnkap_dlnr_row(4)
      save

      double precision :: delta_z, delta_x, logr, saved_r
      logical :: extrapolate_linear
      integer :: i, ii

      delta_z = abs(metal_fraction-opacity_table%alex06_cached_z)
      delta_x = abs(hydrogen_fraction-opacity_table%alex06_cached_x)
!     ENSURE THAT OPACITY TABLE HAS THE SAME X,Z VALUE AS THE POINT
      if (delta_z.gt.1.0d-8 .or. delta_x.gt.1.0d-8) then
         opacity_table%alex06_cached_z = metal_fraction
         opacity_table%alex06_cached_x = hydrogen_fraction
         call alex06tab
      endif
!     COMPUTE R FOR GRID POINT
      logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     Insure that index IT is within the required array bounds
      if (opacity_table%alex06_index_t .lt. 1) opacity_table%alex06_index_t=1
      if ((opacity_table%alex06_index_t+2) .gt. num_t) opacity_table%alex06_index_t=num_t-2
!     FIND NEAREST GRID POINTS IN T.
      if (log10_temperature.lt.opacity_table%alex06_grid_logt(opacity_table%alex06_index_t+2)) then
         do i = opacity_table%alex06_index_t+1,2,-1
            if (log10_temperature.gt.opacity_table%alex06_grid_logt(i)) then
               opacity_table%alex06_index_t = i - 1
               goto 10
            endif
         end do
         opacity_table%alex06_index_t = 1
   10    continue
      else
         do i = opacity_table%alex06_index_t+3,num_t
            if (log10_temperature.lt.opacity_table%alex06_grid_logt(i)) then
               opacity_table%alex06_index_t = i - 2
               opacity_table%alex06_index_t = min(num_t-3,opacity_table%alex06_index_t)
               goto 20
            endif
         end do
         opacity_table%alex06_index_t = num_t - 3
   20    continue
      endif
!     INTERPOLATION FACTORS IN LOG T
      do i = 1,4
         interp_nodes(i) = opacity_table%alex06_grid_logt(opacity_table%alex06_index_t+i-1)
      end do
      call interp(interp_nodes, weight_t, dweight_t, log10_temperature)
!     Insure that index ID is within the required array bounds
      if (opacity_table%alex06_index_r .lt. 1) opacity_table%alex06_index_r=1
      if ((opacity_table%alex06_index_r+2) .gt. num_d) opacity_table%alex06_index_r=num_d-2
!     FIND NEAREST GRID POINTS IN R = RHO/T6**3
      if (logr.lt.opacity_table%alex06_grid_logr(opacity_table%alex06_index_r+2)) then
         do i = opacity_table%alex06_index_r+1,2,-1
            if (logr.gt.opacity_table%alex06_grid_logr(i)) then
               opacity_table%alex06_index_r = i - 1
               goto 30
            endif
         end do
         opacity_table%alex06_index_r = 1
   30    continue
      else
         do i = opacity_table%alex06_index_r+3,num_d
            if (logr.lt.opacity_table%alex06_grid_logr(i)) then
               opacity_table%alex06_index_r = i - 2
               opacity_table%alex06_index_r = min(num_d-3,opacity_table%alex06_index_r)
               goto 40
            endif
         end do
         opacity_table%alex06_index_r = num_d - 3
   40    continue
      endif
!     INTERPOLATION FACTORS IN LOG R
      if (logr.gt.opacity_table%alex06_grid_logr(num_d)) then
         extrapolate_linear = .true.
         saved_r = logr
         logr = opacity_table%alex06_grid_logr(num_d)
      else
         extrapolate_linear = .false.
      endif
      do i = 1,4
         interp_nodes(i) = opacity_table%alex06_grid_logr(opacity_table%alex06_index_r+i-1)
      end do
      call interp(interp_nodes, weight_r, dweight_r, logr)
!     INTERPOLATE IN LOG R AT FIXED T
      do i = 1,4
         ii = opacity_table%alex06_index_t + i - 1
         opacity_row(i) = weight_r(1)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r) + &
              weight_r(2)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+1) + &
              weight_r(3)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+2) + &
              weight_r(4)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+3)
         dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r) + &
              dweight_r(2)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+1) + &
              dweight_r(3)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+2) + &
              dweight_r(4)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+3)
      end do
      if (extrapolate_linear) then
         do i = 1,4
            ii = opacity_table%alex06_index_t+i-1
            opacity_row(i) = opacity_row(i)+(saved_r-logr)* &
                 (opacity_table%alex06_opacity(ii,num_d)-opacity_table%alex06_opacity(ii,num_d-1))/ &
                 (opacity_table%alex06_grid_logr(num_d)-opacity_table%alex06_grid_logr(num_d-1))
         end do
         logr = saved_r
      endif
!     INTERPOLATE IN T
      log10_opacity = weight_t(1)*opacity_row(1)+weight_t(2)*opacity_row(2)+ &
           weight_t(3)*opacity_row(3)+weight_t(4)*opacity_row(4)
!     D LN CAPPA/D LN T AT FIXED R
      dlnkap_dlnt = dweight_t(1)*opacity_row(1)+dweight_t(2)*opacity_row(2)+ &
           dweight_t(3)*opacity_row(3)+dweight_t(4)*opacity_row(4)
!     INTERPOLATE IN D LN CAPPA/ D LN R AT FIXED T
      dlnkap_dlnrho = weight_t(1)*dlnkap_dlnr_row(1) + weight_t(2)*dlnkap_dlnr_row(2) + &
           weight_t(3)*dlnkap_dlnr_row(3) + weight_t(4)*dlnkap_dlnr_row(4)
!     CORRECT FROM DERIVATIVE AT FIXED R TO DERIVATIVE AT FIXED RHO
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = exp(ln10*log10_opacity)
      return
end subroutine getalex06
