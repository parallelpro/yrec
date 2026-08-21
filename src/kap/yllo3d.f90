! THIS IS THE INTERPOLATION FACILITY FOR THE LIVERMORE OPACITY TABLES
! USING CUBIC SPLINE INTERPLOATION SCHEME.
!*** YLLO3D AND SETLLO ARE THE ROUTINES TO BE CALLED!
!      SUBROUTINE YLLO3D(DL,TL,X,OF,OLF,QODF,QOTF)
!      SUBROUTINE YLLO2D(T,D,M1,M2,M3,O,OL,QODF,QOTF)
!      SUBROUTINE YLLO3D2(DL,TL,X,OF,OLF,QODF,QOTF)
!      SUBROUTINE YLLO2D2(T,D,M1,M2,M3,O,OL,QODF,QOTF)
!      SUBROUTINE SETLLO(IOCTRL)                     READ LLOT
!      SUBROUTINE YLLOC                              SPLINE COEFF.
!      SUBROUTINE LL4TH(X) SURFACE TABLES
!
!----------------------------------------------------------------------
! yllo3d
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original yllo3d.f; only variable names, source form, and comment
! style were updated. common/gllot/ and common/llot/ member names
! match those established in ll4th.f90/setllo.f90.
!
! 3D interpolation (X, T, rho) in the OPAL92 opacity tables. Finds
! the nearest X table(s), converts to the table's internal (T6,
! rho/T6**3) coordinates, calls yllo2d for the 2D interpolation, and
! linearly extrapolates/interpolates in X if needed.
subroutine yllo3d(log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      use findex_mod
      implicit none
      integer, parameter :: num_t = 50
      integer, parameter :: num_d = 17
      integer, parameter :: num_x = 3
      integer, parameter :: num_xt = num_t*num_x
      integer, parameter :: num_4d = 4*num_d

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      double precision :: opal92_grid_logt(num_t), opal92_grid_x(num_x), &
           opal92_grid_logr(num_d)
      common /gllot/ opal92_grid_logt, opal92_grid_x, opal92_grid_logr
      double precision :: opal92_log10_opacity(num_xt, num_d)
      integer :: opal92_num_x, opal92_num_temps
      common /llot/ opal92_log10_opacity, opal92_num_x, opal92_num_temps
      double precision :: opal92_surface_x, opal92_surface_z, &
           opal92_surface_spline_coeffs(num_t,num_4d)
      integer :: opal92_surface_x_index
      common /llot4/opal92_surface_x, opal92_surface_z, &
           opal92_surface_spline_coeffs, opal92_surface_x_index
      integer :: abund_index, temp_index, dens_index
      common /kipmll/abund_index, temp_index, dens_index
      data abund_index, temp_index, dens_index/1,1,1/
      save

      logical :: single_x_table
      double precision :: t6, rhot3
      integer :: im1
      double precision :: o0, ol0, qod0, qot0, o1, ol1, qod1, qot1
      double precision :: qodi, qoti, grdnt

      single_x_table = .false.
! INDEPENDENT PARAMETER IN LIVERMORE OPACITY TABLE;
! T6 = LN(T/10E6)
      t6 = log10_temperature - 6.0d0
! RHOT3 = LN(RHO/(T/10E6)**3)
      rhot3 = log10_density - 3.0d0*t6

      if (dabs(opal92_surface_x-hydrogen_fraction).le.1.0d-5) then
         abund_index = 4
         single_x_table = .true.
         go to 131
      endif
      do 130 im1 = 1,num_x
         if (dabs(opal92_grid_x(im1)-hydrogen_fraction).le.1.0d-5) then
            abund_index = im1
            single_x_table = .true.
            go to 131
         endif
 130  continue
      call findex(opal92_grid_x, num_x, hydrogen_fraction, abund_index)
      if (abund_index.lt.0) abund_index = -abund_index
      if (abund_index.le.1.and.rhot3.gt.-1.0d0) abund_index = 2
      if (abund_index.ge.3) abund_index = 2
      if (abund_index.le.0) stop ' ERROR IN X GRID'
 131  continue
      call findex(opal92_grid_logt, opal92_num_temps, t6, temp_index)
      if (temp_index.lt.0.and.opal92_grid_logt(opal92_num_temps).eq.t6) temp_index = -temp_index
      if (temp_index.lt.0) stop ' T OUT OF TABLE '
      call yllo2d(t6, rhot3, abund_index, temp_index, dens_index, o0, ol0, qod0, qot0)
      if (single_x_table) then
! USE ONLY ONE X TABLE
         log10_opacity = ol0
         opacity = o0
         qodi = qod0
         qoti = qot0
      else
! LINEAR EXTRAPOLATION IN X
         call yllo2d(t6, rhot3, abund_index+1, temp_index, dens_index, o1, ol1, qod1, qot1)
         grdnt = (hydrogen_fraction-opal92_grid_x(abund_index))/(opal92_grid_x(abund_index+1)-opal92_grid_x(abund_index))
         log10_opacity = (ol1-ol0)*grdnt + ol0
         qodi = (qod1-qod0)*grdnt + qod0
         qoti = (qot1-qot0)*grdnt + qot0
         opacity = 10.0d0**log10_opacity
      endif
! CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
      dlnkap_dlnrho = qodi
      dlnkap_dlnt = qoti - 3.0d0*dlnkap_dlnrho

      return
end subroutine yllo3d
