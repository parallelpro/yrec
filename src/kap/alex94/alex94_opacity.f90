!----------------------------------------------------------------------
! alex94_opacity
!----------------------------------------------------------------------
! The low-temperature molecular opacities of Alexander & Ferguson
! (1994/95): seven files, one per X, each with tables in Z, log T and
! log R.
!     read_alex94_tables          read the seven files, then build the
!                                 fixed-Z tables at the run's Z
!     build_alex94_surface_table  the fixed-X table at the envelope X
!     compute_alex94_opacity      the opacity and its log derivatives
! The tables, the compositions they were built at and the lookup's
! cursors are an alex94_table_set (opacity_tables), which kap_lib holds
! and passes in. The controls name this family alex95 (use_alex95_tables,
! alex95_table_paths), the code after the paper's year. The sources are
! at the end of the file.
module alex94_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: alex_composition_tol, alex94_surface_slot, &
           alex94_table_set, n_alex94_d, n_alex94_t, n_alex94_x, n_alex94_z, &
           stencil4_locate, kap_settings_t, kap_table_paths_t, opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use numerics_lib, only: interp, lagrange4
      use physical_constants, only: ln10
      use alex94_reader, only: read_alex94_tables
      use math_lib
      implicit none
      private
      public :: alex94_source
      public :: build_alex94_fixed_z_table, build_alex94_surface_table, compute_alex94_opacity

contains


!----------------------------------------------------------------------
! build_alex94_fixed_z_table (formerly alex94_fixed_z_table)
!----------------------------------------------------------------------
! Builds the fixed-Z tables tbl%opacity(1:7,:,:), one per X, by 4-point
! Lagrange interpolation in Z of the full set, and records their Z.
subroutine build_alex94_fixed_z_table(tbl, metal_fraction)
      type(alex94_table_set), intent(inout) :: tbl


      real(dp), intent(in) :: metal_fraction

      real(dp) :: interp_z(4), weight_z(4)
      integer :: i, jj, kk, idz, ii1, ii2, ii3, ii4

!     Locate four nearest tables in Z
      do i=3,n_alex94_z-1
         if (metal_fraction < tbl%grid_z(i)) then
            idz = i - 2
            exit
         end if
      end do
      if (i > (n_alex94_z-1)) then
      idz = n_alex94_z - 3
      end if
      do i = 1,4
         interp_z(i) = tbl%grid_z(idz+i-1)
      end do
!     Get interpolation factors for Z.
      call interp(interp_z, weight_z, x_eval=metal_fraction)
!     Interpolate in Z at fixed X.
      do i = 1,n_alex94_x
!        the slots of the four Z tables of this X (seven X per Z, from
!        Z = 0 up)
         ii1 = i + n_alex94_x*(idz-1)
         ii2 = ii1 + n_alex94_x
         ii3 = ii2 + n_alex94_x
         ii4 = ii3 + n_alex94_x
         do jj = 1,n_alex94_t
            do kk = 1,n_alex94_d
               tbl%opacity(i,jj,kk) = weight_z(1)*tbl%full_opacity(ii1,jj,kk)+ &
                    weight_z(2)*tbl%full_opacity(ii2,jj,kk) + &
                    weight_z(3)*tbl%full_opacity(ii3,jj,kk) + &
                    weight_z(4)*tbl%full_opacity(ii4,jj,kk)
            end do
         end do
      end do
      tbl%cached_z = metal_fraction
end subroutine build_alex94_fixed_z_table

!----------------------------------------------------------------------
! build_alex94_surface_table (formerly alex94_surface_table)
!----------------------------------------------------------------------
! Builds the surface table, tbl%opacity(alex94_surface_slot,:,:), at
! hydrogen_fraction by 4-point Lagrange interpolation in X of the
! fixed-Z tables, and records its X.
subroutine build_alex94_surface_table(tbl, hydrogen_fraction)
      type(alex94_table_set), intent(inout) :: tbl


      real(dp), intent(in) :: hydrogen_fraction

      real(dp) :: interp_x(4), weight_x(4)
      integer :: i, j

      tbl%cached_x = hydrogen_fraction
!     Find 4 nearest tables in X.
      tbl%index_x = alex_x_stencil_start(tbl%grid_x, hydrogen_fraction)
      do i = 1,4
         interp_x(i) = tbl%grid_x(tbl%index_x+i-1)
      end do
      call interp(interp_x, weight_x, x_eval=hydrogen_fraction)
      do i = 1,n_alex94_t
         do j = 1,n_alex94_d
            tbl%opacity(alex94_surface_slot,i,j) = weight_x(1)*tbl%opacity(tbl%index_x,i,j)+ &
                 weight_x(2)*tbl%opacity(tbl%index_x+1,i,j) + &
                 weight_x(3)*tbl%opacity(tbl%index_x+2,i,j) + &
                 weight_x(4)*tbl%opacity(tbl%index_x+3,i,j)
         end do
      end do
end subroutine build_alex94_surface_table

!----------------------------------------------------------------------
! compute_alex94_opacity (formerly alex94_interp3d)
!----------------------------------------------------------------------
! The opacity at (log rho, log T, X, Z) and its log derivatives, by
! 4-point Lagrange interpolation in log T and log R = log rho - 3 log T6:
! in the surface table at its X (extrapolating linearly above the top
! log R), else also in X among the fixed-Z tables. A Z other than the
! fixed-Z tables' rebuilds them, and the surface table, at the lookup's
! composition.
subroutine compute_alex94_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      type(alex94_table_set), intent(inout) :: tbl


      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      real(dp) :: interp_nodes(4), weight_t(4), dweight_t(4), &
           weight_r(4), dweight_r(4), weight_x(4), opacity_row(4), &
           dlnkap_dlnr_row(4), opacity_x(4), dlnkap_dlnt_x(4), &
           dlnkap_dlnr_x(4)
      real(dp) :: delta_z, logr, saved_r
      logical :: extrapolate_linear
      integer :: i, j, jj, ii

      delta_z = abs(metal_fraction-tbl%cached_z)
!     tables at this Z
      if (delta_z > alex_composition_tol) then
         call build_alex94_fixed_z_table(tbl, metal_fraction)
         call build_alex94_surface_table(tbl, hydrogen_fraction)
      end if
!     Compute R for grid point
      logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     Find nearest grid points in T.
!     (the T cursor clamped into the stencil's range first)
      if (tbl%index_t < 1) tbl%index_t=1
      if ((tbl%index_t+2) > n_alex94_t) tbl%index_t=n_alex94_t-2
      call stencil4_locate(tbl%grid_logt, n_alex94_t, log10_temperature, tbl%index_t)
!     Interpolation factors in log T
      do i = 1,4
         interp_nodes(i) = tbl%grid_logt(tbl%index_t+i-1)
      end do
      call interp(interp_nodes, weight_t, dweight_t, log10_temperature)
!     Find nearest grid points in R = rho/T6**3
      call stencil4_locate(tbl%grid_logr, n_alex94_d, logr, tbl%index_r)
!     Interpolation factors in log R
      if (logr > tbl%grid_logr(n_alex94_d).and. &
           abs(hydrogen_fraction-tbl%cached_x) < alex_composition_tol) then
         extrapolate_linear = .true.
         saved_r = logr
         logr = tbl%grid_logr(n_alex94_d)
      else
         extrapolate_linear = .false.
      end if
      do i = 1,4
         interp_nodes(i) = tbl%grid_logr(tbl%index_r+i-1)
      end do
      call interp(interp_nodes, weight_r, dweight_r, logr)
!     Now either interpolate in surface X table or calculate opacity at
!     4 different values of X and interpolate in X.
      if (abs(hydrogen_fraction-tbl%cached_x) < alex_composition_tol) then
!        Surface abundance table
!        interpolate in log R at fixed T
         do i = 1,4
            ii = tbl%index_t+i - 1
            opacity_row(i) = weight_r(1)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r) + &
                 weight_r(2)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r+1) + &
                 weight_r(3)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r+2) + &
                 weight_r(4)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r+3)
            dlnkap_dlnr_row(i) = dweight_r(1)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r) + &
                 dweight_r(2)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r+1) + &
                 dweight_r(3)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r+2) + &
                 dweight_r(4)*tbl%opacity(alex94_surface_slot,ii,tbl%index_r+3)
         end do
         if (extrapolate_linear) then
            do i = 1,4
               ii = tbl%index_t+i-1
               opacity_row(i) = opacity_row(i)+(saved_r-logr)* &
                    (tbl%opacity(alex94_surface_slot,ii,n_alex94_d)-tbl%opacity(alex94_surface_slot,ii,n_alex94_d-1))/ &
                    (tbl%grid_logr(n_alex94_d)-tbl%grid_logr(n_alex94_d-1))
            end do
            logr = saved_r
         end if
!        Interpolate in T
         log10_opacity = lagrange4(weight_t, opacity_row)
!        D ln cappa/d ln T at fixed R
         dlnkap_dlnt = lagrange4(dweight_t, opacity_row)
!        Interpolate in d ln cappa/ d ln R at fixed T
         dlnkap_dlnrho = lagrange4(weight_t, dlnkap_dlnr_row)
!        Correct from derivative at fixed R to derivative at fixed rho
         dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
         opacity = exp(ln10*log10_opacity)
      else
!        Find 4 nearest tables in X.
         tbl%index_x = alex_x_stencil_start(tbl%grid_x, hydrogen_fraction)
         do i = 1,4
            interp_nodes(i) = tbl%grid_x(tbl%index_x+i-1)
         end do
         call interp(interp_nodes, weight_x, x_eval=hydrogen_fraction)
!        Interpolate in log R at fixed T
         do j = 1,4
            jj = tbl%index_x+j-1
            do i = 1,4
               ii = tbl%index_t+i - 1
               opacity_row(i) = weight_r(1)*tbl%opacity(jj,ii,tbl%index_r) + &
                    weight_r(2)*tbl%opacity(jj,ii,tbl%index_r+1) + &
                    weight_r(3)*tbl%opacity(jj,ii,tbl%index_r+2) + &
                    weight_r(4)*tbl%opacity(jj,ii,tbl%index_r+3)
               dlnkap_dlnr_row(i) = dweight_r(1)*tbl%opacity(jj,ii,tbl%index_r) + &
                    dweight_r(2)*tbl%opacity(jj,ii,tbl%index_r+1) + &
                    dweight_r(3)*tbl%opacity(jj,ii,tbl%index_r+2) + &
                    dweight_r(4)*tbl%opacity(jj,ii,tbl%index_r+3)
            end do
!           Interpolate in T
            opacity_x(j) = lagrange4(weight_t, opacity_row)
!           D ln cappa/d ln T at fixed R
            dlnkap_dlnt_x(j) = lagrange4(dweight_t, opacity_row)
!           Interpolate in d ln cappa/ d ln R at fixed T
            dlnkap_dlnr_x(j) = lagrange4(weight_t, dlnkap_dlnr_row)
         end do
!        Interpolate in X
         log10_opacity = lagrange4(weight_x, opacity_x)
         dlnkap_dlnt = lagrange4(weight_x, dlnkap_dlnt_x)
         dlnkap_dlnrho = lagrange4(weight_x, dlnkap_dlnr_x)
!        Correct from derivative at fixed R to derivative at fixed rho
         dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
         opacity = exp(ln10*log10_opacity)
      end if
end subroutine compute_alex94_opacity

! alex_x_stencil_start
! First of the four Alexander-1994 X tables bracketing hydrogen_fraction
! (the 7-column X grid is walked as a fixed ladder: 1, 2, 3 or 4).
pure integer function alex_x_stencil_start(grid_x, hydrogen_fraction)
      real(dp), intent(in) :: grid_x(:), hydrogen_fraction

      if (hydrogen_fraction < grid_x(4)) then
         if (hydrogen_fraction > grid_x(3)) then
            alex_x_stencil_start = 2
         else
            alex_x_stencil_start = 1
         end if
      else
         if (hydrogen_fraction > grid_x(5)) then
            alex_x_stencil_start = 4
         else
            alex_x_stencil_start = 3
         end if
      end if
end function alex_x_stencil_start

!----------------------------------------------------------------------
! The family's descriptor (opacity_source) and its procedures
!----------------------------------------------------------------------
! A molecular family (the controls call it alex95: use_alex95_tables,
! alex95_table_paths): fixed-Z tables at alex_table_z1, rebuilt when a
! lookup's Z differs, and a surface-X table at the envelope X.
function alex94_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'Alexander 1994'
      source%selected = settings%use_alex95_tables
      source%covers_any_z = .true.
      source%load => load_alex94
      source%refresh => refresh_alex94
      source%lookup => lookup_alex94
end function alex94_source

subroutine load_alex94(settings, paths, envelope_hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! (the reader reports its own errors)
      message = ''
      call read_alex94_tables(tables%alex94, paths%alex94, ierr)
      if (ierr /= 0) return
      call build_alex94_fixed_z_table(tables%alex94, settings%alex_table_z1)
      call build_alex94_surface_table(tables%alex94, envelope_hydrogen_fraction)
end subroutine load_alex94

subroutine refresh_alex94(settings, hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      real(dp), intent(in) :: hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      ierr = 0
      message = ''
      call build_alex94_surface_table(tables%alex94, hydrogen_fraction)
end subroutine refresh_alex94

subroutine lookup_alex94(settings, tables, point, res, in_table, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      ierr = 0
      message = ''
      in_table = .true.
      call compute_alex94_opacity(tables%alex94, point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, point%metal_fraction, &
           res%opacity, res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt)
end subroutine lookup_alex94

!----------------------------------------------------------------------
! Sources
!----------------------------------------------------------------------
! The tables: Alexander & Ferguson (1994, ApJ 437, 879), seven files in
! input/opacity/molecular/legacy/alex/. The code: the F77 ALXTBL,
! ALXZTAB, ALX8TH and YALO3D (MHP 4/98; the path argument 8/25; the
! cursor clamp, llp 8/19/08), converted to Fortran 90 in 2026, the tables
! now a table set passed in and the reads checked (M220).

end module alex94_opacity
