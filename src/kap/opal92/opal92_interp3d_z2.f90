!----------------------------------------------------------------------
! opal92_interp3d_z2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original yllo3d2.f; only variable names, source form, and comment
! style were updated.
!
! DBG 5/94 opacity at different Z. 3D interpolation (X, T, rho) in
! the second (different-Z) set of OPAL92 opacity tables. Mirrors
! opal92_interp3d but reads the *_z2 table members and calls
! opal92_interp2d_z2.
subroutine opal92_interp3d_z2(log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)
      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
! 2026 ierr campaign: failures return via ierr (kap_eval gates).
      integer, intent(out) :: ierr
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

! The warm-start cursors opal92_index_x_z2/opal92_index_t_z2/opal92_index_rho_z2
! are members of opacity_table (opacity_table_lib).
      logical :: single_x_table
      double precision :: t6, rhot3
      integer :: im1
      double precision :: o0, ol0, qod0, qot0, o1, ol1, qod1, qot1
      double precision :: qodi, qoti, grdnt

      ierr = 0
      single_x_table = .false.
! INDEPENDENT PARAMETER IN LIVERMORE OPACITY TABLE;
! T6 = LOG10(T/10E6)
      t6 = log10_temperature - 6.0d0
! RHOT3 = LOG10(RHO/(T/10E6)**3)
      rhot3 = log10_density - 3.0d0*t6

      if (dabs(opacity_table%opal92_surface_x_z2-hydrogen_fraction).le.opal92_x_match_tol) then
         opacity_table%opal92_index_x_z2 = 4
         single_x_table = .true.
      endif
      if (.not. single_x_table) then
      do im1 = 1,n_opal92_x
         if (dabs(opacity_table%opal92_grid_x_z2(im1)-hydrogen_fraction).le.opal92_x_match_tol) then
            opacity_table%opal92_index_x_z2 = im1
            single_x_table = .true.
            exit
         endif
      end do
      end if
      if (.not. single_x_table) then
      call findex(opacity_table%opal92_grid_x_z2, n_opal92_x, hydrogen_fraction, opacity_table%opal92_index_x_z2)
      if (opacity_table%opal92_index_x_z2.lt.0) opacity_table%opal92_index_x_z2 = -opacity_table%opal92_index_x_z2
      if (opacity_table%opal92_index_x_z2.le.1.and.rhot3.gt.-1.0d0) opacity_table%opal92_index_x_z2 = 2
      if (opacity_table%opal92_index_x_z2.ge.3) opacity_table%opal92_index_x_z2 = 2
      if (opacity_table%opal92_index_x_z2.le.0) then
         write(*,*) 'opal92_interp3d_z2: error in X2 grid'
         ierr = 1
         return
      end if
      end if
      call findex(opacity_table%opal92_grid_logt_z2, opacity_table%opal92_num_temps_z2, t6, opacity_table%opal92_index_t_z2)
      if (opacity_table%opal92_index_t_z2.lt.0.and.opacity_table%opal92_grid_logt_z2(opacity_table%opal92_num_temps_z2).eq.t6) opacity_table%opal92_index_t_z2 = -opacity_table%opal92_index_t_z2
      if (opacity_table%opal92_index_t_z2.lt.0) then
         write(*,*) 'opal92_interp3d_z2: T out of table'
         ierr = 1
         return
      end if
      call opal92_interp2d_z2(t6, rhot3, opacity_table%opal92_index_x_z2, opacity_table%opal92_index_t_z2, opacity_table%opal92_index_rho_z2, o0, ol0, qod0, qot0, ierr)
      if (ierr /= 0) return
      if (single_x_table) then
!>>> USE ONLY ONE X TABLE
         log10_opacity = ol0
         opacity = o0
         qodi = qod0
         qoti = qot0
      else
!>>> LINEAR EXTRAPOLATION IN X
         call opal92_interp2d_z2(t6, rhot3, opacity_table%opal92_index_x_z2+1, opacity_table%opal92_index_t_z2, opacity_table%opal92_index_rho_z2, o1, ol1, qod1, qot1, ierr)
      if (ierr /= 0) return
         grdnt = (hydrogen_fraction-opacity_table%opal92_grid_x_z2(opacity_table%opal92_index_x_z2))/ &
              (opacity_table%opal92_grid_x_z2(opacity_table%opal92_index_x_z2+1)-opacity_table%opal92_grid_x_z2(opacity_table%opal92_index_x_z2))
         log10_opacity = (ol1-ol0)*grdnt + ol0
         qodi = (qod1-qod0)*grdnt + qod0
         qoti = (qot1-qot0)*grdnt + qot0
         opacity = exp10(log10_opacity)
      endif
! CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
      dlnkap_dlnrho = qodi
      dlnkap_dlnt = qoti - 3.0d0*dlnkap_dlnrho

      return
end subroutine opal92_interp3d_z2
