!----------------------------------------------------------------------
! yllo2d2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original yllo2d2.f; only variable names, source form, and comment
! style were updated. common/gllot2/, common/llot2/, common/llot42/,
! and common/lintpl2/ member names match those established in
! ll4th.f90/setllo.f90.
!
! TWO DIMENSIONAL INTERPOLATION FOR OPACITY AND THE DERIVATIVES.
! DBG 5/94 Modified to use second opacity table at different Z
! (OPAL92 second-Z tables). abund_index, temp_index, and dens_index
! are the nearest grid point of abundance, temperature, and density
! found by yllo3d2. opacity is opacity, log10_opacity is
! dlog(opacity), dlnkap_dlnrho is the partial derivative of opacity
! wrt density, dlnkap_dlnt is the partial derivative of opacity wrt
! temperature.
subroutine yllo2d2(temperature, density, abund_index, temp_index, &
     dens_index, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

      use numerics_lib
      implicit none
      integer, parameter :: num_t = 50
      integer, parameter :: num_d = 17
      integer, parameter :: num_x = 3
      integer, parameter :: num_xt = num_t*num_x
      integer, parameter :: num_4d = 4*num_d

      double precision, intent(in) :: temperature, density
      integer, intent(in) :: abund_index, temp_index
      integer, intent(inout) :: dens_index
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      double precision :: opal92_grid_logt_z2(num_t), opal92_grid_x_z2(num_x), &
           opal92_grid_logr_z2(num_d)
      common /gllot2/ opal92_grid_logt_z2, opal92_grid_x_z2, opal92_grid_logr_z2
      double precision :: opal92_log10_opacity_z2(num_xt, num_d)
      integer :: opal92_num_x_z2, opal92_num_temps_z2
      common /llot2/ opal92_log10_opacity_z2, opal92_num_x_z2, opal92_num_temps_z2
      double precision :: opal92_surface_x_z2, opal92_surface_z_z2, &
           opal92_surface_spline_coeffs_z2(num_t,num_4d)
      integer :: opal92_surface_x_index_z2
      common /llot42/opal92_surface_x_z2, opal92_surface_z_z2, &
           opal92_surface_spline_coeffs_z2, opal92_surface_x_index_z2
      double precision :: opal92_spline_coeffs_z2(num_xt, num_4d)
      integer :: opal92_density_start_index_z2(num_xt), opal92_density_count_z2(num_xt)
      common /lintpl2/ opal92_spline_coeffs_z2, opal92_density_start_index_z2, &
           opal92_density_count_z2
      double precision :: xt(num_t), yto(num_t)
      double precision :: aqod(num_t)
      integer :: jt, it, its, itf
      logical :: lmore
      save

      integer :: mm1, index1, ndss, ndf, knot, index2
      double precision :: dx, c1, c2, c3, c4, ol0, qodi
      double precision :: ol00, unused_ddensity_dtemp

      lmore = .true.
!!! FOR SIX GRID POINTS OF TEMPERATURE
      its = temp_index - 2
      if (its.le.0) its = 1
      itf = temp_index + 3
      if (itf.gt.opal92_num_temps_z2) itf = opal92_num_temps_z2
      if (abund_index.lt.4) then
          mm1 = abund_index
       else
          mm1 = opal92_surface_x_index_z2
      endif
      jt = 0
      do 300 it = its,itf
        index1 = it + (mm1-1)*opal92_num_temps_z2
        ndss = opal92_density_start_index_z2(index1)
        if (ndss.ne.1) stop ' OPAL95 2D2 CHECK NDSS '
        ndf = ndss + opal92_density_count_z2(index1) - 1
        call findex(opal92_grid_logr_z2, ndf, density, dens_index)
        if (dens_index.lt.0) then
! OUT SIDE THEN  LINEAR EXTRAPOLATION
          dens_index = -dens_index
          knot = dens_index - ndss + 1
          index2 = 4*(knot-1)
          dx = density - opal92_grid_logr_z2(dens_index)
          if (abund_index.lt.4) then
            c1 = opal92_spline_coeffs_z2(index1,index2+1)
            c2 = opal92_spline_coeffs_z2(index1,index2+2)
          else
            c1 = opal92_surface_spline_coeffs_z2(it,index2+1)
            c2 = opal92_surface_spline_coeffs_z2(it,index2+2)
          endif
          ol0 = c2*dx + c1
          qodi = c2
        else
! IN SIDE THEN  SPLINE INTERPOLATION
          knot = dens_index - ndss + 1
          index2 = 4*(knot-1)
          dx = density - opal92_grid_logr_z2(dens_index)
          if (abund_index.lt.4) then
            c1 = opal92_spline_coeffs_z2(index1,index2+1)
            c2 = opal92_spline_coeffs_z2(index1,index2+2)
            c3 = opal92_spline_coeffs_z2(index1,index2+3)
            c4 = opal92_spline_coeffs_z2(index1,index2+4)
          else
            c1 = opal92_surface_spline_coeffs_z2(it,index2+1)
            c2 = opal92_surface_spline_coeffs_z2(it,index2+2)
            c3 = opal92_surface_spline_coeffs_z2(it,index2+3)
            c4 = opal92_surface_spline_coeffs_z2(it,index2+4)
          endif
! INTERPOLATION FOR OPACITY(OL) IN THE ENTRY D AND THE EACH T-GRID
! ESTIMATES THE PARTIAL DERIVATIVE OF OL WRT D
! EVALUATES THE INTERPOLATION VALUE IN THE SUB-RANGE WE DETERMINED.
          ol0 = ((c4*dx+c3)*dx+c2)*dx + c1
          qodi = (3.0d0*c4*dx+2.0d0*c3)*dx + c2
        endif
!
        jt = jt + 1
        xt(jt) = opal92_grid_logt_z2(it)
        yto(jt) = ol0
        aqod(jt) = qodi
  300 continue
      if (xt(1).gt.temperature.or.xt(jt).lt.temperature) stop ' EXTRAPOLATION FAILS '
!! INTERPOLATION FOR THE OPACITY IN THE ENTRY T AND D.
!! GET THE PARTIAL DERIVATIVE OF OL WRT T.
      call intpol(xt, yto, jt, temperature, ol00, dlnkap_dlnt)
      log10_opacity = ol00
      opacity = 10.0d0**log10_opacity
! QOTF = D LN(O)/D LN(T)
!! FIND THE PARTIAL DERIVATIVE VALUE OF OL WRT D IN THE GIVEN T AND D
      call intpol(xt, aqod, jt, temperature, dlnkap_dlnrho, unused_ddensity_dtemp)
! QODF = D LN(O)/D LN(D)
      return
end subroutine yllo2d2
