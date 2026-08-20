!----------------------------------------------------------------------
! op954d
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
subroutine op954d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      implicit none
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126

      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

! FULL SET OF TABLES: OPACITY AS A FUNCTION OF Z AND X, T, RHO/T6**3
! TABLES ARE INCREMENTED IN SETS OF NZ*NX.  SO THE TABLES FOR THE
! THIRD METAL ABUNDANCE (3 X 10**-4)BEGIN AT TABLE 21 AND END AT TABLE 30.
! FOR THE HIGH VALUES OF Z, THE NUMBER OF X TABLES IS NOT THE SAME (I.E.
! X<0.9 FOR Z=0.1).
! FOR EACH COMPOSITION A FULL GRID IN (T,RHO/T6**3) IS RETAINED.
      double precision :: opal95_grid_logt(num_t), opal95_grid_x(num_x), &
           opal95_grid_logr(num_d), opal95_grid_z(num_z), &
           opal95_full_opacity(num_xz,num_t,num_d)
      integer :: opal95_num_x_at_z(num_z), opal95_table_start_index(num_z)
      common /llot95a/ opal95_grid_logt, opal95_grid_x, opal95_grid_logr, &
           opal95_grid_z, opal95_full_opacity, opal95_num_x_at_z, &
           opal95_table_start_index
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
      double precision :: logcappa_at_x(4), dlogcappa_dlogt_at_x(4), &
           dlogcappa_dlogr_at_x(4)
      double precision :: logcappa_at_z(4), dlogcappa_dlogt_at_z(4), &
           dlogcappa_dlogr_at_z(4)
      save

      integer :: i, j, k, z_table_offset, xz_table_index, temp_index, &
           density_index
      double precision :: delta_logr

! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF Z
      do k = 1,4
         z_table_offset = opal95_table_start_index(opal95_index_z+k-1)
! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF X
         do j = 1,4
            xz_table_index = opal95_index_x(k,j) + z_table_offset
! FOR EACH X, GET CAPPA FOR 4 VALUES OF T
            do i = 1,4
               temp_index = opal95_index_t + i - 1
! FOR EACH T, GET CAPPA FOR 4 VALUES OF R
               density_index = opal95_index_rho(i)
! LOG CAPPA AT DESIRED RHO FOR EACH OF THE 4 DESIRED T.
               logcappa_at_t(i) = opal95_weight_rho(i,1)*opal95_full_opacity(xz_table_index,temp_index,density_index) + &
                    opal95_weight_rho(i,2)*opal95_full_opacity(xz_table_index,temp_index,density_index+1) &
                    + opal95_weight_rho(i,3)*opal95_full_opacity(xz_table_index,temp_index,density_index+2) &
                    + opal95_weight_rho(i,4)*opal95_full_opacity(xz_table_index,temp_index,density_index+3)
! D LOG CAPPA/D LOG R FOR EACH OF THE 4 DESIRED T.
               dlogcappa_dlogr_at_t(i) = opal95_dweight_rho(i,1)*opal95_full_opacity(xz_table_index,temp_index,density_index) &
                    + opal95_dweight_rho(i,2)*opal95_full_opacity(xz_table_index,temp_index,density_index+1) &
                    + opal95_dweight_rho(i,3)*opal95_full_opacity(xz_table_index,temp_index,density_index+2) &
                    + opal95_dweight_rho(i,4)*opal95_full_opacity(xz_table_index,temp_index,density_index+3)
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
            logcappa_at_x(j) = opal95_weight_t(1)*logcappa_at_t(1) + opal95_weight_t(2)*logcappa_at_t(2) &
                 + opal95_weight_t(3)*logcappa_at_t(3)
            logcappa_at_x(j) = logcappa_at_x(j) + opal95_weight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG T
            dlogcappa_dlogt_at_x(j) = opal95_dweight_t(1)*logcappa_at_t(1) + opal95_dweight_t(2)*logcappa_at_t(2) &
                 + opal95_dweight_t(3)*logcappa_at_t(3)
            dlogcappa_dlogt_at_x(j) = dlogcappa_dlogt_at_x(j) + opal95_dweight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG R
            dlogcappa_dlogr_at_x(j) = opal95_weight_t(1)*dlogcappa_dlogr_at_t(1) + opal95_weight_t(2)*dlogcappa_dlogr_at_t(2) &
                 + opal95_weight_t(3)*dlogcappa_dlogr_at_t(3)
            dlogcappa_dlogr_at_x(j) = dlogcappa_dlogr_at_x(j) + opal95_weight_t(4)*dlogcappa_dlogr_at_t(4)
         end do
! INTERPOLATE FOR LOG CAPPA IN X.
         logcappa_at_z(k) = opal95_weight_x(k,1)*logcappa_at_x(1) + opal95_weight_x(k,2)*logcappa_at_x(2) + &
              opal95_weight_x(k,3)*logcappa_at_x(3) + &
              opal95_weight_x(k,4)*logcappa_at_x(4)
! INTERPOLATE FOR QOT IN X.
         dlogcappa_dlogt_at_z(k) = opal95_weight_x(k,1)*dlogcappa_dlogt_at_x(1) + opal95_weight_x(k,2)*dlogcappa_dlogt_at_x(2) + &
              opal95_weight_x(k,3)*dlogcappa_dlogt_at_x(3) + &
              opal95_weight_x(k,4)*dlogcappa_dlogt_at_x(4)
! INTERPOLATE FOR QOD IN X.
         dlogcappa_dlogr_at_z(k) = opal95_weight_x(k,1)*dlogcappa_dlogr_at_x(1) + opal95_weight_x(k,2)*dlogcappa_dlogr_at_x(2) + &
              opal95_weight_x(k,3)*dlogcappa_dlogr_at_x(3) + &
              opal95_weight_x(k,4)*dlogcappa_dlogr_at_x(4)
      end do
! INTERPOLATE FOR LOG CAPPA IN Z.
      log10_opacity = opal95_weight_z(1)*logcappa_at_z(1) + opal95_weight_z(2)*logcappa_at_z(2) + &
           opal95_weight_z(3)*logcappa_at_z(3) + &
           opal95_weight_z(4)*logcappa_at_z(4)
! INTERPOLATE FOR QOT IN Z.
      dlnkap_dlnt = opal95_weight_z(1)*dlogcappa_dlogt_at_z(1) + opal95_weight_z(2)*dlogcappa_dlogt_at_z(2) + &
           opal95_weight_z(3)*dlogcappa_dlogt_at_z(3) + &
           opal95_weight_z(4)*dlogcappa_dlogt_at_z(4)
! INTERPOLATE FOR QOD IN Z.
      dlnkap_dlnrho = opal95_weight_z(1)*dlogcappa_dlogr_at_z(1) + opal95_weight_z(2)*dlogcappa_dlogr_at_z(2) + &
           opal95_weight_z(3)*dlogcappa_dlogr_at_z(3) + &
           opal95_weight_z(4)*dlogcappa_dlogr_at_z(4)
! CORRECT FROM DERIVATE AT FIXED R TO DERIVATIVE AT FIXED RHO.
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = 1.0d1**log10_opacity
      return
end subroutine op954d
