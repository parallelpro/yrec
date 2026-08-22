!----------------------------------------------------------------------
! yalo3d
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original yalo3d.f; only variable names, source form, and comment
! style were updated.
!
! MHP 4/98 Computes Alexander 1994 low-temperature opacities by
! 4-point Lagrangian interpolation in log T and log R = rho/T6**3,
! and (away from the envelope X, Z) also in X.
subroutine yalo3d(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use opacity_table_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: num_x = 7
      integer, parameter :: num_z = 15
      integer, parameter :: num_t = 23
      integer, parameter :: num_d = 17
      integer, parameter :: num_xt = 8

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt





      double precision :: interp_nodes(4), weight_t(4), dweight_t(4), &
           weight_r(4), dweight_r(4), weight_x(4), opacity_row(4), &
           dlnkap_dlnr_row(4), opacity_x(4), dlnkap_dlnt_x(4), &
           dlnkap_dlnr_x(4)
      double precision :: delta_z, logr, saved_r
      logical :: extrapolate_linear
      integer :: i, j, jj, ii

      delta_z = abs(metal_fraction-opacity_table%alex95_cached_z)
!     ENSURE THAT OPACITY TABLE HAS THE SAME Z VALUE AS THE ENVELOPE.
      if (delta_z.gt.1.0d-8) then
         call alxztab(metal_fraction)
         call alx8th(hydrogen_fraction)
      endif
!     COMPUTE R FOR GRID POINT
      logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     FIND NEAREST GRID POINTS IN T.
!     Insure that index IT is within the required array bounds  llp  8/19/08
      if (opacity_table%alex95_index_t .lt. 1) opacity_table%alex95_index_t=1
      if ((opacity_table%alex95_index_t+2) .gt. num_t) opacity_table%alex95_index_t=num_t-2
      if (log10_temperature.lt.opacity_table%alex95_grid_logt(opacity_table%alex95_index_t+2)) then
         do i = opacity_table%alex95_index_t+1,2,-1
            if (log10_temperature.gt.opacity_table%alex95_grid_logt(i)) then
               opacity_table%alex95_index_t = i - 1
               exit
            endif
         end do
         if (i < (2)) then
         opacity_table%alex95_index_t = 1
         end if
   10    continue
      else
         do i = opacity_table%alex95_index_t+3,num_t
            if (log10_temperature.lt.opacity_table%alex95_grid_logt(i)) then
               opacity_table%alex95_index_t = i - 2
               opacity_table%alex95_index_t = min(num_t-3,opacity_table%alex95_index_t)
               exit
            endif
         end do
         if (i > num_t) then
         opacity_table%alex95_index_t = num_t - 3
         end if
   20    continue
      endif
!     INTERPOLATION FACTORS IN LOG T
      do i = 1,4
         interp_nodes(i) = opacity_table%alex95_grid_logt(opacity_table%alex95_index_t+i-1)
      end do
      call interp(interp_nodes, weight_t, dweight_t, log10_temperature)
!     FIND NEAREST GRID POINTS IN R = RHO/T6**3
      if (logr.lt.opacity_table%alex95_grid_logr(opacity_table%alex95_index_r+2)) then
         do i = opacity_table%alex95_index_r+1,2,-1
            if (logr.gt.opacity_table%alex95_grid_logr(i)) then
               opacity_table%alex95_index_r = i - 1
               exit
            endif
         end do
         if (i < (2)) then
         opacity_table%alex95_index_r = 1
         end if
   30    continue
      else
         do i = opacity_table%alex95_index_r+3,num_d
            if (logr.lt.opacity_table%alex95_grid_logr(i)) then
               opacity_table%alex95_index_r = i - 2
               opacity_table%alex95_index_r = min(num_d-3,opacity_table%alex95_index_r)
               exit
            endif
         end do
         if (i > num_d) then
         opacity_table%alex95_index_r = num_d - 3
         end if
   40    continue
      endif
!     INTERPOLATION FACTORS IN LOG R
      if (logr.gt.opacity_table%alex95_grid_logr(num_d).and. &
           abs(hydrogen_fraction-opacity_table%alex95_cached_x).lt.1.0d-8) then
         extrapolate_linear = .true.
         saved_r = logr
         logr = opacity_table%alex95_grid_logr(num_d)
      else
         extrapolate_linear = .false.
      endif
      do i = 1,4
         interp_nodes(i) = opacity_table%alex95_grid_logr(opacity_table%alex95_index_r+i-1)
      end do
      call interp(interp_nodes, weight_r, dweight_r, logr)
!     NOW EITHER INTERPOLATE IN SURFACE X TABLE OR CALCULATE OPACITY AT
!     4 DIFFERENT VALUES OF X AND INTERPOLATE IN X.
      if (abs(hydrogen_fraction-opacity_table%alex95_cached_x).lt.1.0d-8) then
!        SURFACE ABUNDANCE TABLE
!        INTERPOLATE IN LOG R AT FIXED T
         do i = 1,4
            ii = opacity_table%alex95_index_t+i - 1
            opacity_row(i) = weight_r(1)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r) + &
                 weight_r(2)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r+1) + &
                 weight_r(3)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r+2) + &
                 weight_r(4)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r+3)
            dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r) + &
                 dweight_r(2)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r+1) + &
                 dweight_r(3)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r+2) + &
                 dweight_r(4)*opacity_table%alex95_opacity(8,ii,opacity_table%alex95_index_r+3)
         end do
         if (extrapolate_linear) then
            do i = 1,4
               ii = opacity_table%alex95_index_t+i-1
               opacity_row(i) = opacity_row(i)+(saved_r-logr)* &
                    (opacity_table%alex95_opacity(8,ii,num_d)-opacity_table%alex95_opacity(8,ii,num_d-1))/ &
                    (opacity_table%alex95_grid_logr(num_d)-opacity_table%alex95_grid_logr(num_d-1))
            end do
            logr = saved_r
         endif
!        INTERPOLATE IN T
         log10_opacity = weight_t(1)*opacity_row(1)+weight_t(2)*opacity_row(2)+ &
              weight_t(3)*opacity_row(3)+weight_t(4)*opacity_row(4)
!        D LN CAPPA/D LN T AT FIXED R
         dlnkap_dlnt = dweight_t(1)*opacity_row(1)+dweight_t(2)*opacity_row(2)+ &
              dweight_t(3)*opacity_row(3)+dweight_t(4)*opacity_row(4)
!        INTERPOLATE IN D LN CAPPA/ D LN R AT FIXED T
         dlnkap_dlnrho = weight_t(1)*dlnkap_dlnr_row(1) + weight_t(2)*dlnkap_dlnr_row(2) + &
              weight_t(3)*dlnkap_dlnr_row(3) + weight_t(4)*dlnkap_dlnr_row(4)
!        CORRECT FROM DERIVATIVE AT FIXED R TO DERIVATIVE AT FIXED RHO
         dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
         opacity = exp(ln10*log10_opacity)
      else
!        FIND 4 NEAREST TABLES IN X.
         if (hydrogen_fraction.lt.opacity_table%alex95_grid_x(4)) then
            if (hydrogen_fraction.gt.opacity_table%alex95_grid_x(3)) then
               opacity_table%alex95_index_x = 2
            else
               opacity_table%alex95_index_x = 1
            endif
         else
            if (hydrogen_fraction.gt.opacity_table%alex95_grid_x(5)) then
               opacity_table%alex95_index_x = 4
            else
               opacity_table%alex95_index_x = 3
            endif
         endif
         do i = 1,4
            interp_nodes(i) = opacity_table%alex95_grid_x(opacity_table%alex95_index_x+i-1)
         end do
         call intrp2(interp_nodes, weight_x, hydrogen_fraction)
!        INTERPOLATE IN LOG R AT FIXED T
         do j = 1,4
            jj = opacity_table%alex95_index_x+j-1
            do i = 1,4
               ii = opacity_table%alex95_index_t+i - 1
               opacity_row(i) = weight_r(1)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r) + &
                    weight_r(2)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r+1) + &
                    weight_r(3)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r+2) + &
                    weight_r(4)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r+3)
               dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r) + &
                    dweight_r(2)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r+1) + &
                    dweight_r(3)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r+2) + &
                    dweight_r(4)*opacity_table%alex95_opacity(jj,ii,opacity_table%alex95_index_r+3)
            end do
!           INTERPOLATE IN T
            opacity_x(j) = weight_t(1)*opacity_row(1)+weight_t(2)*opacity_row(2)+ &
                 weight_t(3)*opacity_row(3)+weight_t(4)*opacity_row(4)
!           D LN CAPPA/D LN T AT FIXED R
            dlnkap_dlnt_x(j) = dweight_t(1)*opacity_row(1) + dweight_t(2)*opacity_row(2) + &
                 dweight_t(3)*opacity_row(3) + dweight_t(4)*opacity_row(4)
!           INTERPOLATE IN D LN CAPPA/ D LN R AT FIXED T
            dlnkap_dlnr_x(j) = weight_t(1)*dlnkap_dlnr_row(1) + weight_t(2)*dlnkap_dlnr_row(2) + &
                 weight_t(3)*dlnkap_dlnr_row(3) + weight_t(4)*dlnkap_dlnr_row(4)
         end do
!        INTERPOLATE IN X
         log10_opacity = weight_x(1)*opacity_x(1) + weight_x(2)*opacity_x(2) + &
              weight_x(3)*opacity_x(3) + weight_x(4)*opacity_x(4)
         dlnkap_dlnt = weight_x(1)*dlnkap_dlnt_x(1) + weight_x(2)*dlnkap_dlnt_x(2) + &
              weight_x(3)*dlnkap_dlnt_x(3) + weight_x(4)*dlnkap_dlnt_x(4)
         dlnkap_dlnrho = weight_x(1)*dlnkap_dlnr_x(1) + weight_x(2)*dlnkap_dlnr_x(2) + &
              weight_x(3)*dlnkap_dlnr_x(3) + weight_x(4)*dlnkap_dlnr_x(4)
!        CORRECT FROM DERIVATIVE AT FIXED R TO DERIVATIVE AT FIXED RHO
         dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
         opacity = exp(ln10*log10_opacity)
      endif
      return
end subroutine yalo3d
