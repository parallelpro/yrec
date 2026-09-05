!----------------------------------------------------------------------
! opal95_interp4d
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original op954d.f; only variable names, source form, and comment
! style were updated.
!
! MHP 7/98 GET OPAL95 OPACITY FROM FULL SET OF OPACITY TABLES.
! 4D interpolation (Z, X, T, rho) across the full OPAL95 table set.
! Called from getopal95 when neither X nor Z match a precomputed
! reduced table.
! 2026 wave 3 (R5): the T/rho stencil (indices, weights, log R
! extrapolation state) that getopal95 computes per call used to reach
! this routine as hidden module state (opacity_table%opal95_index_*,
! opal95_weight_*, opal95_logr*, opal95_extrap_*); it now arrives as
! explicit intent(in) arguments, passed straight from those members by
! getopal95. The opacity tables themselves stay in opacity_table.
subroutine opal95_interp4d(index_z, weight_z, index_x, weight_x, &
     index_t, index_rho, weight_t, dweight_t, weight_rho, dweight_rho, &
     logr, logr_lo_edge, logr_hi_edge, extrap_lo, extrap_hi, extrap_hi_row, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use opacity_table_lib
      use numerics_lib, only: lagrange4
      use math_lib
      implicit none

! Z stencil, and the X stencil (index_x(k,:), weight_x(k,:)) at each
! of its four Z rows.
      integer, intent(in) :: index_z, index_x(4,4)
      double precision, intent(in) :: weight_z(4), weight_x(4,4)
      integer, intent(in) :: index_t, index_rho(4)
      double precision, intent(in) :: weight_t(4), dweight_t(4), &
           weight_rho(4,4), dweight_rho(4,4)
      double precision, intent(in) :: logr, logr_lo_edge, logr_hi_edge(4)
      logical, intent(in) :: extrap_lo, extrap_hi, extrap_hi_row(4)
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      double precision :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      double precision :: logcappa_at_x(4), dlogcappa_dlogt_at_x(4), &
           dlogcappa_dlogr_at_x(4)
      double precision :: logcappa_at_z(4), dlogcappa_dlogt_at_z(4), &
           dlogcappa_dlogr_at_z(4)
      integer :: i, j, k, z_table_offset, xz_table_index, temp_index, &
           density_index
      double precision :: delta_logr

! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF Z
      do k = 1,4
         z_table_offset = opacity_table%opal95_table_start_index(index_z+k-1)
! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF X
         do j = 1,4
            xz_table_index = index_x(k,j) + z_table_offset
! FOR EACH X, GET CAPPA FOR 4 VALUES OF T
            do i = 1,4
               temp_index = index_t + i - 1
! FOR EACH T, GET CAPPA FOR 4 VALUES OF R
               density_index = index_rho(i)
! LOG CAPPA AT DESIRED RHO FOR EACH OF THE 4 DESIRED T.
               logcappa_at_t(i) = weight_rho(i,1)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index) + &
                    weight_rho(i,2)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+1) &
                    + weight_rho(i,3)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+2) &
                    + weight_rho(i,4)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+3)
! D LOG CAPPA/D LOG R FOR EACH OF THE 4 DESIRED T.
               dlogcappa_dlogr_at_t(i) = dweight_rho(i,1)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index) &
                    + dweight_rho(i,2)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+1) &
                    + dweight_rho(i,3)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+2) &
                    + dweight_rho(i,4)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+3)
            end do
! CHECK ON WHETHER RHO IS OUTSIDE OF THE TABLE AND NEEDS EXTRAPOLATION
            if (extrap_lo) then
! FACTOR IN R
               delta_logr = logr - logr_lo_edge
! CORRECT CAPPA BY USING THE DERIVATIVE AT THE BOUNDARY.
               do i = 1,4
                  logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
               end do
            else if (extrap_hi) then
               do i = 1,4
                  if (extrap_hi_row(i)) then
                     delta_logr = logr - logr_hi_edge(i)
                     logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
                  endif
               end do
            endif
! INTERPOLATE FOR LOG CAPPA IN T. (Inline rather than lagrange4: the
! original accumulates three terms, stores, then adds the fourth.)
            logcappa_at_x(j) = weight_t(1)*logcappa_at_t(1) + weight_t(2)*logcappa_at_t(2) &
                 + weight_t(3)*logcappa_at_t(3)
            logcappa_at_x(j) = logcappa_at_x(j) + weight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG T
            dlogcappa_dlogt_at_x(j) = dweight_t(1)*logcappa_at_t(1) + dweight_t(2)*logcappa_at_t(2) &
                 + dweight_t(3)*logcappa_at_t(3)
            dlogcappa_dlogt_at_x(j) = dlogcappa_dlogt_at_x(j) + dweight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG R
            dlogcappa_dlogr_at_x(j) = weight_t(1)*dlogcappa_dlogr_at_t(1) + weight_t(2)*dlogcappa_dlogr_at_t(2) &
                 + weight_t(3)*dlogcappa_dlogr_at_t(3)
            dlogcappa_dlogr_at_x(j) = dlogcappa_dlogr_at_x(j) + weight_t(4)*dlogcappa_dlogr_at_t(4)
         end do
! INTERPOLATE FOR LOG CAPPA IN X. (Inline rather than lagrange4: the
! weights are the strided row weight_x(k,1:4).)
         logcappa_at_z(k) = weight_x(k,1)*logcappa_at_x(1) + weight_x(k,2)*logcappa_at_x(2) + &
              weight_x(k,3)*logcappa_at_x(3) + &
              weight_x(k,4)*logcappa_at_x(4)
! INTERPOLATE FOR QOT IN X.
         dlogcappa_dlogt_at_z(k) = weight_x(k,1)*dlogcappa_dlogt_at_x(1) + weight_x(k,2)*dlogcappa_dlogt_at_x(2) + &
              weight_x(k,3)*dlogcappa_dlogt_at_x(3) + &
              weight_x(k,4)*dlogcappa_dlogt_at_x(4)
! INTERPOLATE FOR QOD IN X.
         dlogcappa_dlogr_at_z(k) = weight_x(k,1)*dlogcappa_dlogr_at_x(1) + weight_x(k,2)*dlogcappa_dlogr_at_x(2) + &
              weight_x(k,3)*dlogcappa_dlogr_at_x(3) + &
              weight_x(k,4)*dlogcappa_dlogr_at_x(4)
      end do
! INTERPOLATE FOR LOG CAPPA IN Z.
      log10_opacity = lagrange4(weight_z, logcappa_at_z)
! INTERPOLATE FOR QOT IN Z.
      dlnkap_dlnt = lagrange4(weight_z, dlogcappa_dlogt_at_z)
! INTERPOLATE FOR QOD IN Z.
      dlnkap_dlnrho = lagrange4(weight_z, dlogcappa_dlogr_at_z)
! CORRECT FROM DERIVATE AT FIXED R TO DERIVATIVE AT FIXED RHO.
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = pow(1.0d1, log10_opacity)
      return
end subroutine opal95_interp4d
