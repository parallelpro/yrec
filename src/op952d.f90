!----------------------------------------------------------------------
! op952d
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original op952d.f; only variable names, source form, and comment
! style were updated.
!
! MHP 7/98 GET OPAL95 OPACITY FROM TABLE AT SURFACE X,Z
! 2D interpolation (T, rho) in the surface-X table produced by
! op95xtab. Called from getopal95 when both X and Z match the model
! surface composition exactly.
subroutine op952d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      implicit none
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

! AS ABOVE FOR DESIRED SURFACE VALUE OF X.
      double precision :: opal95_surface_opacity(num_t,num_d), opal95_surface_x
      common /llot95e/opal95_surface_opacity, opal95_surface_x
! INDICES FOR INTERPOLATION IN Z,X,T, AND R
      integer :: opal95_index_z, opal95_index_x(4,4), opal95_index_t, &
           opal95_index_rho(4)
      common/op95indx/ opal95_index_z, opal95_index_x, opal95_index_t, &
           opal95_index_rho
! INTERPOLATION FACTORS FOR Z,X,T, AND R, AS WELL AS DERIVATIVE
! FACTORS FOR T AND RHO.
      double precision :: opal95_weight_z(4), opal95_weight_x(4,4), &
           opal95_weight_t(4), opal95_dweight_t(4), opal95_weight_rho(4,4), &
           opal95_dweight_rho(4,4)
      common/op95fact/ opal95_weight_z, opal95_weight_x, opal95_weight_t, &
           opal95_dweight_t, opal95_weight_rho, opal95_dweight_rho
! DATA FOR LINEAR EXTRAPOLATION WHEN OUTSIDE TABLE IN RHO
      double precision :: opal95_logr, opal95_logr_lo_edge, opal95_logr_hi_edge(4)
      logical :: opal95_extrap_lo, opal95_extrap_hi, opal95_extrap_hi_row(4)
      common/op95ext/ opal95_logr, opal95_logr_lo_edge, opal95_logr_hi_edge, &
           opal95_extrap_lo, opal95_extrap_hi, opal95_extrap_hi_row
      double precision :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      save

      integer :: i, temp_index, density_index
      double precision :: delta_logr

! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF T.
      do i = 1,4
         temp_index = opal95_index_t + i - 1
! FOR EACH T, INTERPOLATE IN 4 VALUES OF R
         density_index = opal95_index_rho(i)
! LOG CAPPA AT DESIRED RHO FOR EACH OF THE 4 DESIRED T.
         logcappa_at_t(i) = opal95_weight_rho(i,1)*opal95_surface_opacity(temp_index,density_index) &
              + opal95_weight_rho(i,2)*opal95_surface_opacity(temp_index,density_index+1) &
              + opal95_weight_rho(i,3)*opal95_surface_opacity(temp_index,density_index+2) &
              + opal95_weight_rho(i,4)*opal95_surface_opacity(temp_index,density_index+3)
! D LOG CAPPA/D LOG R FOR EACH OF THE 4 DESIRED T.
         dlogcappa_dlogr_at_t(i) = opal95_dweight_rho(i,1)*opal95_surface_opacity(temp_index,density_index) &
              + opal95_dweight_rho(i,2)*opal95_surface_opacity(temp_index,density_index+1) &
              + opal95_dweight_rho(i,3)*opal95_surface_opacity(temp_index,density_index+2) &
              + opal95_dweight_rho(i,4)*opal95_surface_opacity(temp_index,density_index+3)
      end do
! CHECK ON WHETHER RHO IS OUTSIDE OF THE TABLE AND NEEDS EXTRAPOLATION
      if (opal95_extrap_lo) then
! FACTOR IN R
         delta_logr = opal95_logr - opal95_logr_lo_edge
! CORRECT CAPPA BY USING THE DERIVATIVE AT THE BOUNDARY.
         do i = 1,4
            logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
         end do
      else if (opal95_extrap_hi) then
         do i = 1,4
            if (opal95_extrap_hi_row(i)) then
               delta_logr = opal95_logr - opal95_logr_hi_edge(i)
               logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
            endif
         end do
      endif
! INTERPOLATE FOR LOG CAPPA IN T.
      log10_opacity = opal95_weight_t(1)*logcappa_at_t(1) + opal95_weight_t(2)*logcappa_at_t(2) &
           + opal95_weight_t(3)*logcappa_at_t(3) + opal95_weight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG T
      dlnkap_dlnt = opal95_dweight_t(1)*logcappa_at_t(1) + opal95_dweight_t(2)*logcappa_at_t(2) &
           + opal95_dweight_t(3)*logcappa_at_t(3) + opal95_dweight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG R
      dlnkap_dlnrho = opal95_weight_t(1)*dlogcappa_dlogr_at_t(1) + opal95_weight_t(2)*dlogcappa_dlogr_at_t(2) &
           + opal95_weight_t(3)*dlogcappa_dlogr_at_t(3) + opal95_weight_t(4)*dlogcappa_dlogr_at_t(4)
! CORRECT FROM DERIVATE AT FIXED R TO DERIVATIVE AT FIXED RHO.
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = 1.0d1**log10_opacity
      return
end subroutine op952d
