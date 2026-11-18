!----------------------------------------------------------------------
! kurucz_opacity
!----------------------------------------------------------------------
! The low-temperature molecular opacities of Kurucz (1990), one table
! per metallicity in log T and log rho (a second at another Z when the
! run blends two).
!     build_kurucz_splines     each T row's spline in log rho, from the
!                              table kurucz_reader read
!     compute_kurucz_opacity   the opacity and its log derivatives from
!                              the splines in log rho, interpolated in
!                              log T; in_table is .false. outside the
!                              table (the first table also has a
!                              log rho <= -3, log T <= 4.1 range check)
! The tables, their splines and the lookup's cursors are kurucz_table_sets
! (opacity_tables), which kap_lib holds and passes in. The sources are at
! the end of the file.
module kurucz_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kurucz_max_num_densities, kurucz_max_num_temps, &
           kurucz_table_set, kap_settings_t, kap_table_paths_t, opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use opacity_blending, only: blend_results_in_z
      use numerics_lib, only: findex, natural_spline, splint, &
           compute_spline_coefficients, numerics_status_text
      use kurucz_reader, only: read_kurucz_tables
      use math_lib
      implicit none
      private
      public :: kurucz90_source
      public :: build_kurucz_splines, compute_kurucz_opacity

! The first table is used only up to this log rho (g/cm^3) and log T (K)
! (kurucz_table_set%check_range)
      real(dp), parameter :: kurucz_logrho_max = -3.0d0, kurucz_logt_max = 4.1d0
! A T row needs at least this many densities with data to be splined
! (build_kurucz_splines), and a lookup at least this many T rows with
! data at its density for the spline in log T
      integer, parameter :: min_row_densities = 25, min_spline_rows = 4
! the spline coefficient work array's interval count, more than any
! row's densities
      integer, parameter :: max_spline_intervals = 100

contains


!----------------------------------------------------------------------
! build_kurucz_splines
!----------------------------------------------------------------------
! Builds, for each T row of tbl, the natural cubic spline of log kappa in
! log rho over the row's densities with data, as coefficients per
! interval (tbl%spline_coeffs), with the row's first density index and
! count. A row that does not start at the first density, or has fewer
! than min_row_densities (25) points, returns ierr = 1 with the reason in
! message.
subroutine build_kurucz_splines(tbl, message, ierr)
      type(kurucz_table_set), intent(inout) :: tbl
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      real(dp) :: spline_work(4,max_spline_intervals), density_nodes(kurucz_max_num_densities), &
           row_log10_opacity(kurucz_max_num_densities), &
           row_second_derivs(kurucz_max_num_densities)
      integer :: row, column, num_filled, interval, i

      ierr = 0
      message = ''
      do row = 1,tbl%num_temps
!        the row's densities with data (kappa > 0)
         num_filled = 0
         do column = 1,kurucz_max_num_densities
            if (tbl%opacity(row,column) <= 0.0d0) cycle
            if (num_filled <= 0) tbl%density_start_index(row) = column
            num_filled = num_filled + 1
            density_nodes(num_filled) = tbl%log10_rho(row,column)
            row_log10_opacity(num_filled) = log10(tbl%opacity(row,column))
         end do
         tbl%density_count(row) = num_filled
         if (tbl%density_start_index(row) /= 1) then
            write(message,'(a,f6.3,a)') 'the Kurucz 1990 table''s row at log T', &
                 tbl%grid_logt(row), ' does not start at its first density'
            ierr = 1
            return
         end if
         if (tbl%density_count(row) < min_row_densities) then
            write(message,'(a,f6.3,a)') 'the Kurucz 1990 table''s row at log T', &
                 tbl%grid_logt(row), ' has fewer than 25 densities'
            ierr = 1
            return
         end if
! the row's natural cubic spline in log rho, as coefficients per interval
         call natural_spline(density_nodes, row_log10_opacity, num_filled, row_second_derivs)
         call compute_spline_coefficients(density_nodes, row_log10_opacity, &
              row_second_derivs, num_filled, spline_work)
         do interval = 1,num_filled
            do i = 1,4
               tbl%spline_coeffs(row,i + (interval-1)*4) = spline_work(i,interval)
            end do
         end do
      end do

end subroutine build_kurucz_splines

!----------------------------------------------------------------------
! compute_kurucz_opacity
!----------------------------------------------------------------------
! The opacity at (log rho, log T) and its log derivatives from table set
! tbl (the first or the second Z's), by cubic splines: in log rho on each
! of up to six T rows around the point (kurucz_rows_at_density: every
! row if fewer than four of the six have data at this density), then
! natural splines in log T through the rows' log kappa and through
! their dln kappa/dln rho. tbl%ix_t and tbl%ix_rho are the warm-start
! cursors (the T row and the density knot) the previous call left.
! in_table is .false. for a point outside the table (and, for the first
! table, above log rho = -3 or log T = 4.1), where the caller uses the
! interior opacity alone. A point with fewer than four T rows of data
! at its density, or a failed spline, returns ierr with the reason in
! message.
subroutine compute_kurucz_opacity(tbl, log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, in_table, message, ierr)
      type(kurucz_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! the rows with data at the density: their log T, log kappa and
! dln kappa/dln rho there, and the T splines' second derivatives
      real(dp) :: row_logt(kurucz_max_num_temps), row_log10_opacity(kurucz_max_num_temps), &
           row_dlnkap_dlnrho(kurucz_max_num_temps), second_derivs(kurucz_max_num_temps)
      integer :: num_rows, first_row, last_row, spline_ierr

      call findex(tbl%grid_logt, tbl%num_temps, log10_temperature, &
           tbl%ix_t)

      ierr = 0
      message = ''
      in_table = .false.
      if (tbl%check_range) then
         if (.not.(log10_density <= kurucz_logrho_max .and. &
              log10_temperature <= kurucz_logt_max)) return
      end if
!     six T rows around the point
      first_row = tbl%ix_t-2
      if (first_row <= 0) first_row = 1
      last_row = tbl%ix_t+3
      if (last_row > tbl%num_temps) last_row = tbl%num_temps
      call kurucz_rows_at_density(tbl, first_row, last_row, log10_density, num_rows, &
           row_logt, row_log10_opacity, row_dlnkap_dlnrho)
      if (num_rows < min_spline_rows) then
!        every row of the table
         call kurucz_rows_at_density(tbl, 1, tbl%num_temps, log10_density, num_rows, &
              row_logt, row_log10_opacity, row_dlnkap_dlnrho)
         if (num_rows < min_spline_rows) then
            write(message,'(a,f8.4,a,f8.4)') 'the Kurucz 1990 table has no value at log rho =', &
                 log10_density, ', log T =', log10_temperature
            ierr = 1
            return
         end if
      end if
      if (row_logt(1) > log10_temperature .or. &
           row_logt(num_rows) < log10_temperature) return
      in_table = .true.
!     log kappa and its T derivative from the spline in log T
      call natural_spline(row_logt, row_log10_opacity, num_rows, second_derivs)
      call splint(row_logt, row_log10_opacity, num_rows, second_derivs, log10_temperature, &
           log10_opacity, spline_ierr, dy=dlnkap_dlnt)
      if (spline_ierr /= 0) then
         message = 'the Kurucz 1990 table''s spline in log T failed: ' // &
              trim(numerics_status_text(spline_ierr))
         ierr = spline_ierr
         return
      end if
      opacity = exp10(log10_opacity)
!     and the rho derivative, from its own spline in log T
      call natural_spline(row_logt, row_dlnkap_dlnrho, num_rows, second_derivs)
      call splint(row_logt, row_dlnkap_dlnrho, num_rows, second_derivs, log10_temperature, &
           dlnkap_dlnrho, spline_ierr)
      if (spline_ierr /= 0) then
         message = 'the Kurucz 1990 table''s spline in log T failed: ' // &
              trim(numerics_status_text(spline_ierr))
         ierr = spline_ierr
         return
      end if
end subroutine compute_kurucz_opacity

! The T rows first_row to last_row of tbl that have data at
! log10_density (num_rows of them, in order): each row's log T, and its
! log kappa and dln kappa/dln rho there from the row's cubic spline in
! log rho.
subroutine kurucz_rows_at_density(tbl, first_row, last_row, log10_density, num_rows, &
     row_logt, row_log10_opacity, row_dlnkap_dlnrho)
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(in) :: first_row, last_row
      real(dp), intent(in) :: log10_density
      integer, intent(out) :: num_rows
      real(dp), intent(inout) :: row_logt(:), row_log10_opacity(:), row_dlnkap_dlnrho(:)
      integer :: row, knot, coeff_base
      real(dp) :: delta_logrho, c1, c2, c3, c4
      logical :: found

      num_rows = 0
      do row = first_row, last_row
         call locate_kurucz_density(tbl, row, log10_density, found)
         if (.not. found) cycle
!        (tbl%ix_rho, tbl%ix_rho+1) is the interval of the row's densities
!        that holds log rho; its spline coefficients
         knot = tbl%ix_rho - tbl%density_start_index(row) + 1
         coeff_base = 4*(knot-1)
         delta_logrho = log10_density - tbl%log10_rho(row, tbl%ix_rho)
         c1 = tbl%spline_coeffs(row, coeff_base+1)
         c2 = tbl%spline_coeffs(row, coeff_base+2)
         c3 = tbl%spline_coeffs(row, coeff_base+3)
         c4 = tbl%spline_coeffs(row, coeff_base+4)
         num_rows = num_rows+1
         row_logt(num_rows) = tbl%grid_logt(row)
         row_log10_opacity(num_rows) = ((c4*delta_logrho+c3)*delta_logrho+c2)* &
              delta_logrho+c1
         row_dlnkap_dlnrho(num_rows) = (3.0d0*c4*delta_logrho+2.0d0*c3)*delta_logrho+c2
      end do
end subroutine kurucz_rows_at_density

! The density knot of T row `row` at or below log10_density, searched
! from the cursor tbl%ix_rho (reset to the row's last density when it is
! outside the row) and left in it; found is .false. when the row has no
! data at that density (below its first density or above its last).
subroutine locate_kurucz_density(tbl, row, log10_density, found)
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(in) :: row
      real(dp), intent(in) :: log10_density
      logical, intent(out) :: found
      integer :: first_density, last_density, knot, scan

      found = .false.
      first_density = tbl%density_start_index(row)
      last_density = first_density + tbl%density_count(row) - 1
      if (tbl%ix_rho < first_density .or. &
           tbl%ix_rho > last_density) tbl%ix_rho = last_density
      knot = tbl%ix_rho
      if (log10_density < tbl%log10_rho(row, knot)) then
         do scan = knot-1, first_density, -1
            if (tbl%log10_rho(row, scan) <= log10_density) then
               knot = scan
               found = .true.
               exit
            end if
         end do
      else
         do scan = knot, last_density-1
            if (tbl%log10_rho(row, scan+1) > log10_density) then
               knot = scan
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            if (tbl%log10_rho(row, last_density) >= log10_density) then
               knot = last_density
               found = .true.
            end if
         end if
      end if
      if (found) tbl%ix_rho = knot
end subroutine locate_kurucz_density

!----------------------------------------------------------------------
! The family's descriptor (opacity_source) and its procedures
!----------------------------------------------------------------------
! A molecular family of one table per Z (kurucz_table_z1, and with
! use_two_z_tables kurucz_table_z2, blended linearly in Z), with no
! value outside its table (in_table .false.: the interior opacity is
! used alone there).
function kurucz90_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'Kurucz 1990'
      source%blends_two_z = .true.
      source%table_z = [settings%kurucz_table_z1, settings%kurucz_table_z2]
      source%table_z_controls = [character(len=16) :: 'kurucz_table_z1', 'kurucz_table_z2']
      source%selected = settings%use_kurucz90_tables
      source%load => load_kurucz90
      source%lookup => lookup_kurucz90
end function kurucz90_source

subroutine load_kurucz90(settings, paths, envelope_hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! (the reader reports its own errors)
      message = ''
      call read_kurucz_tables(tables%kurucz, paths%kurucz, paths%kurucz_z2, &
           settings%use_two_z_tables, ierr)
      if (ierr /= 0) return
      call build_kurucz_splines(tables%kurucz(1), message, ierr)
      if (ierr /= 0) return
      if (settings%use_two_z_tables) call build_kurucz_splines(tables%kurucz(2), message, ierr)
end subroutine load_kurucz90

subroutine lookup_kurucz90(settings, tables, point, res, in_table, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
      type(kap_result) :: res_z2

      call compute_kurucz_opacity(tables%kurucz(1), point%log10_density, &
           point%log10_temperature, res%opacity, res%log10_opacity, &
           res%dlnkap_dlnrho, res%dlnkap_dlnt, in_table, message, ierr)
      if (ierr /= 0) return
      if (in_table .and. settings%use_two_z_tables) then
         call compute_kurucz_opacity(tables%kurucz(2), point%log10_density, &
              point%log10_temperature, res_z2%opacity, res_z2%log10_opacity, &
              res_z2%dlnkap_dlnrho, res_z2%dlnkap_dlnt, in_table, message, ierr)
         if (ierr /= 0) return
      end if
      if (in_table .and. settings%use_two_z_tables) then
         call blend_results_in_z(res, res_z2, point%metal_fraction, &
              settings%kurucz_table_z1, settings%kurucz_table_z2)
      end if
end subroutine lookup_kurucz90

!----------------------------------------------------------------------
! Sources
!----------------------------------------------------------------------
! The tables: Kurucz (1990) Rosseland mean opacities, one file per
! [M/H] in input/opacity/molecular/legacy/kurmol/. The code: YCK's F77
! spline set-up (3/91), the second Z table (DBG 12/95) and the KURUCZ
! lookup, converted to Fortran 90 in 2026: the lookup's alternate return
! became in_table, the two Z tables' copies one routine on a table set
! passed in, the reads checked and bounded (M220, M229).

end module kurucz_opacity
