!----------------------------------------------------------------------
! opal92_interp3d
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original yllo3d.f; only variable names, source form, and comment
! style were updated.
!
! THIS IS THE INTERPOLATION FACILITY FOR THE LIVERMORE (OPAL92)
! OPACITY TABLES USING CUBIC SPLINE INTERPOLATION. Companion routines
! in this directory: read_opal92_tables (read), opal92_table_prep
! (spline coefficients), opal92_surface_table (surface X,Z table),
! opal92_interp2d (2D in T,rho).
!
! 3D interpolation (X, T, rho) in the OPAL92 opacity tables. Finds
! the nearest X table(s), converts to the table's internal (T6,
! rho/T6**3) coordinates, calls opal92_interp2d for the 2D
! interpolation, and linearly extrapolates/interpolates in X if
! needed.
!
! 2026 wave 3 (R5): tbl is the table set (opacity_table%opal92(1) or
! (2)); the former opal92_interp3d_z2.f90 was a member-for-member
! rename of this file (its two stdout messages named it; its grdnt
! line was only continued differently) and is deleted.
subroutine opal92_interp3d(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)
      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none

      type(opal92_table_set), intent(inout) :: tbl
      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
! 2026 ierr campaign: failures return via ierr (kap_eval gates).
      integer, intent(out) :: ierr
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

! tbl%index_x/index_t/index_rho (the warm-start cursors) default to
! 1 in opacity_table_lib.f90 (opal92_table_set).
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

      if (dabs(tbl%surface_x-hydrogen_fraction).le.opal92_x_match_tol) then
         tbl%index_x = 4
         single_x_table = .true.
      endif
      if (.not. single_x_table) then
      do im1 = 1,n_opal92_x
         if (dabs(tbl%grid_x(im1)-hydrogen_fraction).le.opal92_x_match_tol) then
            tbl%index_x = im1
            single_x_table = .true.
            exit
         endif
      end do
      end if
      if (.not. single_x_table) then
      call findex(tbl%grid_x, n_opal92_x, hydrogen_fraction, tbl%index_x)
      if (tbl%index_x.lt.0) tbl%index_x = -tbl%index_x
      if (tbl%index_x.le.1.and.rhot3.gt.-1.0d0) tbl%index_x = 2
      if (tbl%index_x.ge.3) tbl%index_x = 2
      if (tbl%index_x.le.0) then
         write(*,*) 'opal92_interp3d: error in X grid'
         ierr = 1
         return
      end if
      end if
      call findex(tbl%grid_logt, tbl%num_temps, t6, tbl%index_t)
      if (tbl%index_t.lt.0.and.tbl%grid_logt(tbl%num_temps).eq.t6) tbl%index_t = -tbl%index_t
      if (tbl%index_t.lt.0) then
         write(*,*) 'opal92_interp3d: T out of table'
         ierr = 1
         return
      end if
      call opal92_interp2d(tbl, t6, rhot3, tbl%index_x, tbl%index_t, tbl%index_rho, o0, ol0, qod0, qot0, ierr)
      if (ierr /= 0) return
      if (single_x_table) then
! USE ONLY ONE X TABLE
         log10_opacity = ol0
         opacity = o0
         qodi = qod0
         qoti = qot0
      else
! LINEAR EXTRAPOLATION IN X
         call opal92_interp2d(tbl, t6, rhot3, tbl%index_x+1, tbl%index_t, tbl%index_rho, o1, ol1, qod1, qot1, ierr)
      if (ierr /= 0) return
         grdnt = (hydrogen_fraction-tbl%grid_x(tbl%index_x))/(tbl%grid_x(tbl%index_x+1)-tbl%grid_x(tbl%index_x))
         log10_opacity = (ol1-ol0)*grdnt + ol0
         qodi = (qod1-qod0)*grdnt + qod0
         qoti = (qot1-qot0)*grdnt + qot0
         opacity = exp10(log10_opacity)
      endif
! CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
      dlnkap_dlnrho = qodi
      dlnkap_dlnt = qoti - 3.0d0*dlnkap_dlnrho

      return
end subroutine opal92_interp3d
