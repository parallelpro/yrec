!----------------------------------------------------------------------
! alex06_opacity
!----------------------------------------------------------------------
! The low-temperature molecular opacities of Ferguson et al. (2005), the
! 'Alexander 2006' tables: one file of tables in X and Z, each in log T
! and log R.
!     read_alex06_tables       read the file
!     compute_alex06_opacity   the opacity and its log derivatives, from
!                              a table at the zone's (X, Z) built by
!                              build_alex06_composition_table and
!                              reused while the composition holds
! The tables, the built table, its composition (the cache) and the
! lookup's cursors are an alex06_table_set (opacity_tables), which
! kap_lib holds and passes in. The sources are at the end of the file.
module alex06_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: alex_composition_tol, alex06_table_set, n_alex06_d, &
           n_alex06_t, n_alex06_x, n_alex06_z, stencil4_locate, kap_settings_t, &
           kap_table_paths_t, opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use numerics_lib, only: interp, lagrange4
      use physical_constants, only: ln10
      use alex06_reader, only: read_alex06_tables
      use math_lib
      implicit none
      private
      public :: alex06_source
      public :: compute_alex06_opacity

contains


!----------------------------------------------------------------------
! build_alex06_composition_table (formerly alex06tab)
!----------------------------------------------------------------------
! Builds tbl%opacity, the table (in log T and log R) at one composition,
! by 4-point Lagrange interpolation of the full set in Z, then in X. A
! composition off the tables (X < 0 or > 1 - Z; Z < 0 or more than one
! node step above the top Z) returns ierr = 1 with the reason in
! message, the table untouched.
! The table is built at (hydrogen_fraction, metal_fraction); the caller
! (compute_alex06_opacity) records that composition as the cache's only
! when the build succeeds.
subroutine build_alex06_composition_table(tbl, hydrogen_fraction, metal_fraction, message, ierr)
      type(alex06_table_set), intent(inout) :: tbl

      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      real(dp) :: interp_nodes(4), weight_z(4), weight_x(4)
      real(dp) :: opacity_by_x(4,n_alex06_t,n_alex06_d)
      real(dp) :: x_max, z_max, interp_target
      integer :: i, j, iz, k, kk, kk2, kk3, kk4

      ierr = 0
      message = ''
      x_max = 1.0d0 - metal_fraction
!     Check that the requested composition is inside table bounds
      if (hydrogen_fraction < 0.0d0 .or. hydrogen_fraction > x_max) then
         write(message,5) hydrogen_fraction, metal_fraction
    5    format('the composition X =',f8.5,', Z =',f8.5,' is outside the Alexander 2006 tables')
         ierr = 1
         return
      end if
!     Permit extrapolation in Z by up to 1 table element
      z_max = tbl%grid_z(n_alex06_z)+(tbl%grid_z(n_alex06_z)-tbl%grid_z(n_alex06_z-1))
      if (metal_fraction < 0.0d0 .or. metal_fraction > z_max) then
         write(message,5) hydrogen_fraction, metal_fraction
         ierr = 1
         return
      end if
!     Find 4 nearest tables in Z.
      do i = 3,n_alex06_z-2
         if (metal_fraction <= tbl%grid_z(i)) then
            iz = i - 2
            exit
         end if
      end do
      if (i > (n_alex06_z-2)) then
      iz = n_alex06_z - 3
      end if
!     Find 4 nearest tables in X.
      do i = 3,n_alex06_x-2
         if (hydrogen_fraction <= tbl%grid_x(i)) then
            tbl%index_x = i - 2
            exit
         end if
      end do
      if (i > (n_alex06_x-2)) then
!     No table for X > 0.9 if Z =0.10 or more
      if (metal_fraction >= 0.1d0) then
         tbl%index_x = n_alex06_x - 4
      else
         tbl%index_x = n_alex06_x - 3
      end if
      end if
!     Interpolation factors for Z
      do i = 1,4
         interp_nodes(i) = tbl%grid_z(iz+i-1)
      end do
      interp_target = metal_fraction
      call interp(interp_nodes, weight_z, x_eval=interp_target)
!     The Z = 0.1 case has one table fewer (its X = 1 - Z table is the
!     X = 0.9 one), which only the fourth X position meets: the first
!     three take the four consecutive Z slots of their X.
      do k = 1,3
         kk = n_alex06_z*(tbl%index_x+k-2)+iz
         do i = 1,n_alex06_t
            do j = 1,n_alex06_d
               opacity_by_x(k,i,j) = lagrange4(weight_z, tbl%full_opacity(kk:kk+3,i,j))
            end do
         end do
      end do
!     At the top X position and the top Z stencil, the Z = 0.1 table of
!     the X = 1 - Z family is the X = 0.9 table (slot n_alex06_z (X index
!     + 2), its first).
      if (tbl%index_x == n_alex06_x-3.and.iz == n_alex06_z-3) then
!        Use different indexing for the last table
         kk = n_alex06_z*(tbl%index_x+2)+iz
         kk2 = kk + 1
         kk3 = kk2 + 1
         kk4 = n_alex06_z*(tbl%index_x+2)
      else
         kk = n_alex06_z*(tbl%index_x+2)+iz
         kk2 = kk + 1
         kk3 = kk2 + 1
         kk4 = kk3 + 1
      end if
!     (Inline rather than lagrange4: kk4 is not kk+3 in the Z=0.10 case.)
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            opacity_by_x(4,i,j) = weight_z(1)*tbl%full_opacity(kk,i,j)+ &
                 weight_z(2)*tbl%full_opacity(kk2,i,j) + &
                 weight_z(3)*tbl%full_opacity(kk3,i,j) + &
                 weight_z(4)*tbl%full_opacity(kk4,i,j)
         end do
      end do
!     Now do X interpolation
!     interpolation factors for X
      do i = 1,3
         interp_nodes(i) = tbl%grid_x(tbl%index_x+i-1)
      end do
! the fourth table is the X = 1 - Z family's (at Z = 0.1, the X = 0.9
! table) when the stencil ends at the last X position
      if (tbl%index_x == n_alex06_x-3) then
         interp_nodes(4) = 1.0d0-metal_fraction
      else
         interp_nodes(4) = tbl%grid_x(tbl%index_x+3)
      end if
      interp_target = hydrogen_fraction
      call interp(interp_nodes, weight_x, x_eval=interp_target)
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            tbl%opacity(i,j) = lagrange4(weight_x, opacity_by_x(1:4,i,j))
         end do
      end do
end subroutine build_alex06_composition_table

!----------------------------------------------------------------------
! compute_alex06_opacity (formerly getalex06)
!----------------------------------------------------------------------
! The opacity at (log rho, log T, X, Z) and its log derivatives, by
! 4-point Lagrange interpolation in log T and log R = log rho - 3 log T6
! of the table at the composition, rebuilt (build_alex06_composition_table)
! when X or Z moves more than alex_composition_tol from the cache's. Above
! the top log R the opacity is extrapolated linearly from the last two
! columns. A composition off the tables returns ierr with the reason in
! message.
subroutine compute_alex06_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, message, ierr)
      type(alex06_table_set), intent(inout) :: tbl


      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      real(dp) :: interp_nodes(4), weight_t(4), dweight_t(4), &
           weight_r(4), dweight_r(4), opacity_row(4), dlnkap_dlnr_row(4)
      real(dp) :: delta_z, delta_x, logr, saved_r
      logical :: extrapolate_linear
      integer :: i, ii

      ierr = 0
      message = ''
      delta_z = abs(metal_fraction-tbl%cached_z)
      delta_x = abs(hydrogen_fraction-tbl%cached_x)
!     Ensure that opacity table has the same X,Z value as the point
      if (delta_z > alex_composition_tol .or. delta_x > alex_composition_tol) then
         call build_alex06_composition_table(tbl, hydrogen_fraction, metal_fraction, message, ierr)
         if (ierr /= 0) return
         tbl%cached_z = metal_fraction
         tbl%cached_x = hydrogen_fraction
      end if
!     Compute R for grid point
      logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     Insure that index IT is within the required array bounds
      if (tbl%index_t < 1) tbl%index_t=1
      if ((tbl%index_t+2) > n_alex06_t) tbl%index_t=n_alex06_t-2
!     Find nearest grid points in T.
      call stencil4_locate(tbl%grid_logt, n_alex06_t, log10_temperature, tbl%index_t)
!     Interpolation factors in log T
      do i = 1,4
         interp_nodes(i) = tbl%grid_logt(tbl%index_t+i-1)
      end do
      call interp(interp_nodes, weight_t, dweight_t, log10_temperature)
!     Insure that index ID is within the required array bounds
      if (tbl%index_r < 1) tbl%index_r=1
      if ((tbl%index_r+2) > n_alex06_d) tbl%index_r=n_alex06_d-2
!     Find nearest grid points in R = rho/T6**3
      call stencil4_locate(tbl%grid_logr, n_alex06_d, logr, tbl%index_r)
!     Interpolation factors in log R
      if (logr > tbl%grid_logr(n_alex06_d)) then
         extrapolate_linear = .true.
         saved_r = logr
         logr = tbl%grid_logr(n_alex06_d)
      else
         extrapolate_linear = .false.
      end if
      do i = 1,4
         interp_nodes(i) = tbl%grid_logr(tbl%index_r+i-1)
      end do
      call interp(interp_nodes, weight_r, dweight_r, logr)
!     Interpolate in log R at fixed T
      do i = 1,4
         ii = tbl%index_t + i - 1
         opacity_row(i) = weight_r(1)*tbl%opacity(ii,tbl%index_r) + &
              weight_r(2)*tbl%opacity(ii,tbl%index_r+1) + &
              weight_r(3)*tbl%opacity(ii,tbl%index_r+2) + &
              weight_r(4)*tbl%opacity(ii,tbl%index_r+3)
         dlnkap_dlnr_row(i) = dweight_r(1)*tbl%opacity(ii,tbl%index_r) + &
              dweight_r(2)*tbl%opacity(ii,tbl%index_r+1) + &
              dweight_r(3)*tbl%opacity(ii,tbl%index_r+2) + &
              dweight_r(4)*tbl%opacity(ii,tbl%index_r+3)
      end do
      if (extrapolate_linear) then
         do i = 1,4
            ii = tbl%index_t+i-1
            opacity_row(i) = opacity_row(i)+(saved_r-logr)* &
                 (tbl%opacity(ii,n_alex06_d)-tbl%opacity(ii,n_alex06_d-1))/ &
                 (tbl%grid_logr(n_alex06_d)-tbl%grid_logr(n_alex06_d-1))
         end do
         logr = saved_r
      end if
!     Interpolate in T
      log10_opacity = lagrange4(weight_t, opacity_row)
!     D ln cappa/d ln T at fixed R
      dlnkap_dlnt = lagrange4(dweight_t, opacity_row)
!     Interpolate in d ln cappa/ d ln R at fixed T
      dlnkap_dlnrho = lagrange4(weight_t, dlnkap_dlnr_row)
!     Correct from derivative at fixed R to derivative at fixed rho
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = exp(ln10*log10_opacity)
end subroutine compute_alex06_opacity

!----------------------------------------------------------------------
! The family's descriptor (opacity_source) and its procedures
!----------------------------------------------------------------------
! A molecular family that interpolates in X and Z itself (a table built
! at each new composition).
function alex06_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'Alexander 2006'
      source%selected = settings%use_alex06_tables
      source%covers_any_z = .true.
      source%load => load_alex06
      source%lookup => lookup_alex06
end function alex06_source

subroutine load_alex06(settings, paths, envelope_hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! (the reader reports its own errors)
      message = ''
      call read_alex06_tables(tables%alex06, paths%alex06, ierr)
end subroutine load_alex06

subroutine lookup_alex06(settings, tables, point, res, in_table, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      call compute_alex06_opacity(tables%alex06, point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, point%metal_fraction, &
           res%opacity, res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
      in_table = ierr == 0
end subroutine lookup_alex06

!----------------------------------------------------------------------
! Sources
!----------------------------------------------------------------------
! The tables: Ferguson et al. (2005, ApJ 623, 585), the low-temperature
! opacities YREC calls Alexander 2006 (input/opacity/molecular/
! alexmol06gs98.tab, at the GS98 mixture). The code: MHP's F77 READALEX06,
! ALEX06TAB and GETALEX06 (3/09; the path argument 8/25), converted to
! Fortran 90 in 2026. Since then: a failed build returns ierr (the F77
! stopped), the cache is set only after a successful build and reset with
! the table (M219), the header checks all return ierr and cover all 19
! columns (M230), the fourth X node is 1 - Z where the table there is
! (M225), the reads are checked (M170), and the tables are a table set
! passed in.

end module alex06_opacity
