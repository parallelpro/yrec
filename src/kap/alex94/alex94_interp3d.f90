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
subroutine alex94_interp3d(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use opacity_table_lib
      use phys_const_lib
      use numerics_lib
      use math_lib
      implicit none

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

      delta_z = abs(metal_fraction-opacity_table%alex94_cached_z)
!     ENSURE THAT OPACITY TABLE HAS THE SAME Z VALUE AS THE ENVELOPE.
      if (delta_z.gt.alex_composition_tol) then
         call alex94_fixed_z_table(metal_fraction)
         call alex94_surface_table(hydrogen_fraction)
      endif
!     COMPUTE R FOR GRID POINT
      logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     FIND NEAREST GRID POINTS IN T.
!     Insure that index IT is within the required array bounds  llp  8/19/08
      if (opacity_table%alex94_index_t .lt. 1) opacity_table%alex94_index_t=1
      if ((opacity_table%alex94_index_t+2) .gt. n_alex94_t) opacity_table%alex94_index_t=n_alex94_t-2
      call stencil4_locate(opacity_table%alex94_grid_logt, n_alex94_t, log10_temperature, opacity_table%alex94_index_t)
!     INTERPOLATION FACTORS IN LOG T
      do i = 1,4
         interp_nodes(i) = opacity_table%alex94_grid_logt(opacity_table%alex94_index_t+i-1)
      end do
      call interp(interp_nodes, weight_t, dweight_t, log10_temperature)
!     FIND NEAREST GRID POINTS IN R = RHO/T6**3
      call stencil4_locate(opacity_table%alex94_grid_logr, n_alex94_d, logr, opacity_table%alex94_index_r)
!     INTERPOLATION FACTORS IN LOG R
      if (logr.gt.opacity_table%alex94_grid_logr(n_alex94_d).and. &
           abs(hydrogen_fraction-opacity_table%alex94_cached_x).lt.alex_composition_tol) then
         extrapolate_linear = .true.
         saved_r = logr
         logr = opacity_table%alex94_grid_logr(n_alex94_d)
      else
         extrapolate_linear = .false.
      endif
      do i = 1,4
         interp_nodes(i) = opacity_table%alex94_grid_logr(opacity_table%alex94_index_r+i-1)
      end do
      call interp(interp_nodes, weight_r, dweight_r, logr)
!     NOW EITHER INTERPOLATE IN SURFACE X TABLE OR CALCULATE OPACITY AT
!     4 DIFFERENT VALUES OF X AND INTERPOLATE IN X.
      if (abs(hydrogen_fraction-opacity_table%alex94_cached_x).lt.alex_composition_tol) then
!        SURFACE ABUNDANCE TABLE
!        INTERPOLATE IN LOG R AT FIXED T
         do i = 1,4
            ii = opacity_table%alex94_index_t+i - 1
            opacity_row(i) = weight_r(1)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r) + &
                 weight_r(2)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+1) + &
                 weight_r(3)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+2) + &
                 weight_r(4)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+3)
            dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r) + &
                 dweight_r(2)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+1) + &
                 dweight_r(3)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+2) + &
                 dweight_r(4)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+3)
         end do
         if (extrapolate_linear) then
            do i = 1,4
               ii = opacity_table%alex94_index_t+i-1
               opacity_row(i) = opacity_row(i)+(saved_r-logr)* &
                    (opacity_table%alex94_opacity(8,ii,n_alex94_d)-opacity_table%alex94_opacity(8,ii,n_alex94_d-1))/ &
                    (opacity_table%alex94_grid_logr(n_alex94_d)-opacity_table%alex94_grid_logr(n_alex94_d-1))
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
         opacity_table%alex94_index_x = alex_x_stencil_start(hydrogen_fraction)
         do i = 1,4
            interp_nodes(i) = opacity_table%alex94_grid_x(opacity_table%alex94_index_x+i-1)
         end do
         call intrp2(interp_nodes, weight_x, hydrogen_fraction)
!        INTERPOLATE IN LOG R AT FIXED T
         do j = 1,4
            jj = opacity_table%alex94_index_x+j-1
            do i = 1,4
               ii = opacity_table%alex94_index_t+i - 1
               opacity_row(i) = weight_r(1)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r) + &
                    weight_r(2)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+1) + &
                    weight_r(3)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+2) + &
                    weight_r(4)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+3)
               dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r) + &
                    dweight_r(2)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+1) + &
                    dweight_r(3)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+2) + &
                    dweight_r(4)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+3)
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
end subroutine alex94_interp3d
