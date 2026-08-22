!----------------------------------------------------------------------
! op953d
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original op953d.f; only variable names, source form, and comment
! style were updated.
!
! MHP 7/98 GET OPAL95 OPACITY FROM TABLE AT FIXED Z
! 3D interpolation (X, T, rho) in the fixed-Z table produced by
! op95ztab. Called from getopal95 when Z matches the model metal
! fraction but X does not match the surface composition.
subroutine op953d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use opacity_table_lib
      implicit none
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      double precision :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      double precision :: logcappa_at_x(4), dlogcappa_dlogt_at_x(4), &
           dlogcappa_dlogr_at_x(4)
      integer :: i, j, x_table_index, temp_index, density_index
      double precision :: delta_logr

! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF X
      do j = 1,4
         x_table_index = opacity_table%opal95_index_x(1,j)
! FOR EACH X, GET CAPPA FOR 4 VALUES OF T
         do i = 1,4
            temp_index = opacity_table%opal95_index_t + i - 1
! FOR EACH T, GET CAPPA FOR 4 VALUES OF R
            density_index = opacity_table%opal95_index_rho(i)
! LOG CAPPA AT DESIRED RHO FOR EACH OF THE 4 DESIRED T.
            logcappa_at_t(i) = opacity_table%opal95_weight_rho(i,1)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index) &
                 + opacity_table%opal95_weight_rho(i,2)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+1) &
                 + opacity_table%opal95_weight_rho(i,3)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+2) &
                 + opacity_table%opal95_weight_rho(i,4)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+3)
! D LOG CAPPA/D LOG R FOR EACH OF THE 4 DESIRED T.
            dlogcappa_dlogr_at_t(i) = opacity_table%opal95_dweight_rho(i,1)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index) &
                 + opacity_table%opal95_dweight_rho(i,2)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+1) &
                 + opacity_table%opal95_dweight_rho(i,3)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+2) &
                 + opacity_table%opal95_dweight_rho(i,4)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+3)
         end do
! CHECK ON WHETHER RHO IS OUTSIDE OF THE TABLE AND NEEDS EXTRAPOLATION
         if (opacity_table%opal95_extrap_lo) then
! FACTOR IN R
            delta_logr = opacity_table%opal95_logr - opacity_table%opal95_logr_lo_edge
! CORRECT CAPPA BY USING THE DERIVATIVE AT THE BOUNDARY.
            do i = 1,4
               logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
            end do
         else if (opacity_table%opal95_extrap_hi) then
            do i = 1,4
               if (opacity_table%opal95_extrap_hi_row(i)) then
                  delta_logr = opacity_table%opal95_logr - opacity_table%opal95_logr_hi_edge(i)
                  logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
               endif
            end do
         endif
! INTERPOLATE FOR LOG CAPPA IN T.
         logcappa_at_x(j) = opacity_table%opal95_weight_t(1)*logcappa_at_t(1) + opacity_table%opal95_weight_t(2)*logcappa_at_t(2) &
              + opacity_table%opal95_weight_t(3)*logcappa_at_t(3) + opacity_table%opal95_weight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG T
         dlogcappa_dlogt_at_x(j) = opacity_table%opal95_dweight_t(1)*logcappa_at_t(1) + opacity_table%opal95_dweight_t(2)*logcappa_at_t(2) &
              + opacity_table%opal95_dweight_t(3)*logcappa_at_t(3) + opacity_table%opal95_dweight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG R
         dlogcappa_dlogr_at_x(j) = opacity_table%opal95_weight_t(1)*dlogcappa_dlogr_at_t(1) + opacity_table%opal95_weight_t(2)*dlogcappa_dlogr_at_t(2) &
              + opacity_table%opal95_weight_t(3)*dlogcappa_dlogr_at_t(3) + opacity_table%opal95_weight_t(4)*dlogcappa_dlogr_at_t(4)
      end do
! INTERPOLATE FOR LOG CAPPA IN X.
      log10_opacity = opacity_table%opal95_weight_x(1,1)*logcappa_at_x(1) + opacity_table%opal95_weight_x(1,2)*logcappa_at_x(2) + &
           opacity_table%opal95_weight_x(1,3)*logcappa_at_x(3) + opacity_table%opal95_weight_x(1,4)*logcappa_at_x(4)
! INTERPOLATE FOR QOT IN X.
      dlnkap_dlnt = opacity_table%opal95_weight_x(1,1)*dlogcappa_dlogt_at_x(1) + opacity_table%opal95_weight_x(1,2)*dlogcappa_dlogt_at_x(2) + &
           opacity_table%opal95_weight_x(1,3)*dlogcappa_dlogt_at_x(3) + opacity_table%opal95_weight_x(1,4)*dlogcappa_dlogt_at_x(4)
! INTERPOLATE FOR QOD IN X.
      dlnkap_dlnrho = opacity_table%opal95_weight_x(1,1)*dlogcappa_dlogr_at_x(1) + opacity_table%opal95_weight_x(1,2)*dlogcappa_dlogr_at_x(2) + &
           opacity_table%opal95_weight_x(1,3)*dlogcappa_dlogr_at_x(3) + opacity_table%opal95_weight_x(1,4)*dlogcappa_dlogr_at_x(4)
! CORRECT FROM DERIVATE AT FIXED R TO DERIVATIVE AT FIXED RHO.
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = 1.0d1**log10_opacity
      return
end subroutine op953d
