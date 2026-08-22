!----------------------------------------------------------------------
! yllo3d2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original yllo3d2.f; only variable names, source form, and comment
! style were updated. common/gllot2/ and common/llot2/ member names
! match those established in ll4th.f90/setllo.f90.
!
! DBG 5/94 opacity at different Z. 3D interpolation (X, T, rho) in
! the second (different-Z) set of OPAL92 opacity tables. Mirrors
! yllo3d but reads the "2" common blocks and calls yllo2d2.
subroutine yllo3d2(log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      use opacity_table_lib
      use numerics_lib
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

! former common/kipmll2/: abund_index/temp_index/dens_index now
! use-associated from opacity_table_lib as abund_index_z2/
! temp_index_z2/dens_index_z2.
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

      if (dabs(opacity_table%opal92_surface_x_z2-hydrogen_fraction).le.1.0d-5) then
         opacity_table%abund_index_z2 = 4
         single_x_table = .true.
         go to 131
      endif
      do im1 = 1,num_x
         if (dabs(opacity_table%opal92_grid_x_z2(im1)-hydrogen_fraction).le.1.0d-5) then
            opacity_table%abund_index_z2 = im1
            single_x_table = .true.
            go to 131
         endif
 130  continue
      end do
      call findex(opacity_table%opal92_grid_x_z2, num_x, hydrogen_fraction, opacity_table%abund_index_z2)
      if (opacity_table%abund_index_z2.lt.0) opacity_table%abund_index_z2 = -opacity_table%abund_index_z2
      if (opacity_table%abund_index_z2.le.1.and.rhot3.gt.-1.0d0) opacity_table%abund_index_z2 = 2
      if (opacity_table%abund_index_z2.ge.3) opacity_table%abund_index_z2 = 2
      if (opacity_table%abund_index_z2.le.0) stop ' ERROR IN X2 GRID'
 131  continue
      call findex(opacity_table%opal92_grid_logt_z2, opacity_table%opal92_num_temps_z2, t6, opacity_table%temp_index_z2)
      if (opacity_table%temp_index_z2.lt.0.and.opacity_table%opal92_grid_logt_z2(opacity_table%opal92_num_temps_z2).eq.t6) opacity_table%temp_index_z2 = -opacity_table%temp_index_z2
      if (opacity_table%temp_index_z2.lt.0) stop ' T OUT OF TABLE '
      call yllo2d2(t6, rhot3, opacity_table%abund_index_z2, opacity_table%temp_index_z2, opacity_table%dens_index_z2, o0, ol0, qod0, qot0)
      if (single_x_table) then
!>>> USE ONLY ONE X TABLE
         log10_opacity = ol0
         opacity = o0
         qodi = qod0
         qoti = qot0
      else
!>>> LINEAR EXTRAPOLATION IN X
         call yllo2d2(t6, rhot3, opacity_table%abund_index_z2+1, opacity_table%temp_index_z2, opacity_table%dens_index_z2, o1, ol1, qod1, qot1)
         grdnt = (hydrogen_fraction-opacity_table%opal92_grid_x_z2(opacity_table%abund_index_z2))/ &
              (opacity_table%opal92_grid_x_z2(opacity_table%abund_index_z2+1)-opacity_table%opal92_grid_x_z2(opacity_table%abund_index_z2))
         log10_opacity = (ol1-ol0)*grdnt + ol0
         qodi = (qod1-qod0)*grdnt + qod0
         qoti = (qot1-qot0)*grdnt + qot0
         opacity = 10.0d0**log10_opacity
      endif
! CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
      dlnkap_dlnrho = qodi
      dlnkap_dlnt = qoti - 3.0d0*dlnkap_dlnrho

      return
end subroutine yllo3d2
