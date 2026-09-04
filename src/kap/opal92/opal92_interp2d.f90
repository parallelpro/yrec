!----------------------------------------------------------------------
! opal92_interp2d
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original yllo2d.f; only variable names, source form, and comment
! style were updated.
!
! TWO DIMENSIONAL INTERPOLATION FOR OPACITY AND THE DERIVATIVES
! (OPAL92 tables). abund_index, temp_index, and dens_index are the
! nearest grid point of abundance, temperature, and density found by
! opal92_interp3d. opacity is opacity, log10_opacity is log(opacity),
! dlnkap_dlnrho is the partial derivative of opacity wrt density,
! dlnkap_dlnt is the partial derivative of opacity wrt temperature.
subroutine opal92_interp2d(temperature, density, abund_index, temp_index, &
     dens_index, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)

      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none

      double precision, intent(in) :: temperature, density
! 2026 ierr campaign: interpolation/extrapolation failures return
! via ierr (kap_eval gates); the historical stops are gone.
      integer, intent(out) :: ierr
      integer, intent(in) :: abund_index, temp_index
      integer, intent(inout) :: dens_index
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      double precision :: xt(n_opal92_t), yto(n_opal92_t)
      double precision :: aqod(n_opal92_t)
      integer :: jt, it, its, itf
      integer :: mm1, index1, ndss, ndf, knot, index2
      double precision :: dx, c1, c2, c3, c4, ol0, qodi
      double precision :: ol00, unused_ddensity_dtemp

      ierr = 0
! FOR SIX GRID POINTS OF TEMPERATURE
      its = temp_index - 2
      if (its.le.0) its = 1
      itf = temp_index + 3
      if (itf.gt.opacity_table%opal92_num_temps) itf = opacity_table%opal92_num_temps
      if (abund_index.lt.4) then
         mm1 = abund_index
      else
         mm1 = opacity_table%opal92_surface_x_index
      endif
      jt = 0
      do it = its,itf
         index1 = it + (mm1-1)*opacity_table%opal92_num_temps
         ndss = opacity_table%opal92_density_start_index(index1)
         if (ndss.ne.1) then
            write(*,*) 'opal92_interp2d: CHECK NDSS'
            ierr = 1
            return
         end if
         ndf = ndss + opacity_table%opal92_density_count(index1) - 1
         call findex(opacity_table%opal92_grid_logr, ndf, density, dens_index)
         if (dens_index.lt.0) then
! OUT SIDE THEN  LINEAR EXTRAPOLATION
            dens_index = -dens_index
            knot = dens_index - ndss + 1
            index2 = 4*(knot-1)
            dx = density - opacity_table%opal92_grid_logr(dens_index)
            if (abund_index.lt.4) then
               c1 = opacity_table%opal92_spline_coeffs(index1,index2+1)
               c2 = opacity_table%opal92_spline_coeffs(index1,index2+2)
            else
               c1 = opacity_table%opal92_surface_spline_coeffs(it,index2+1)
               c2 = opacity_table%opal92_surface_spline_coeffs(it,index2+2)
            endif
            ol0 = c2*dx + c1
            qodi = c2
         else
! IN SIDE THEN  SPLINE INTERPOLATION
            knot = dens_index - ndss + 1
            index2 = 4*(knot-1)
            dx = density - opacity_table%opal92_grid_logr(dens_index)
            if (abund_index.lt.4) then
               c1 = opacity_table%opal92_spline_coeffs(index1,index2+1)
               c2 = opacity_table%opal92_spline_coeffs(index1,index2+2)
               c3 = opacity_table%opal92_spline_coeffs(index1,index2+3)
               c4 = opacity_table%opal92_spline_coeffs(index1,index2+4)
            else
               c1 = opacity_table%opal92_surface_spline_coeffs(it,index2+1)
               c2 = opacity_table%opal92_surface_spline_coeffs(it,index2+2)
               c3 = opacity_table%opal92_surface_spline_coeffs(it,index2+3)
               c4 = opacity_table%opal92_surface_spline_coeffs(it,index2+4)
            endif
! INTERPOLATION FOR OPACITY(OL) IN THE ENTRY D AND THE EACH T-GRID
! ESTIMATES THE PARTIAL DERIVATIVE OF OL WRT D
! EVALUATES THE INTERPOLATION VALUE IN THE SUB-RANGE WE DETERMINED.
            ol0 = ((c4*dx+c3)*dx+c2)*dx + c1
            qodi = (3.0d0*c4*dx+2.0d0*c3)*dx + c2
         endif
!
         jt = jt + 1
         xt(jt) = opacity_table%opal92_grid_logt(it)
         yto(jt) = ol0
         aqod(jt) = qodi
      end do
      if (xt(1).gt.temperature.or.xt(jt).lt.temperature) then
         write(*,*) 'opal92_interp2d: extrapolation fails'
         ierr = 1
         return
      end if
!! INTERPOLATION FOR THE OPACITY IN THE ENTRY T AND D.
!! GET THE PARTIAL DERIVATIVE OF OL WRT T.
      call intpol(xt, yto, jt, temperature, ol00, dlnkap_dlnt, ierr)
      if (ierr /= 0) return
      log10_opacity = ol00
      opacity = exp10(log10_opacity)
! QOTF = D LN(O)/D LN(T)
!! FIND THE PARTIAL DERIVATIVE VALUE OF OL WRT D IN THE GIVEN T AND D
      call intpol(xt, aqod, jt, temperature, dlnkap_dlnrho, unused_ddensity_dtemp, ierr)
      if (ierr /= 0) return
! QODF = D LN(O)/D LN(D)

      return
end subroutine opal92_interp2d
