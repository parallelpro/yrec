!----------------------------------------------------------------------
! kurucz_opacity
!----------------------------------------------------------------------
! The low-temperature molecular opacities of Kurucz (1990), one table
! per metallicity in log T and log rho (a second at another Z when the
! run blends two).
!     read_kurucz_tables       read and spline the table(s)
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
      use kap_messages, only: kap_log, kap_log_values, kap_report, kap_terminal, &
           message_length, unwritten_line
      use kurucz_reader, only: read_kurucz_tables
      use math_lib
      implicit none
      private
      public :: kurucz90_source
      public :: build_kurucz_splines, compute_kurucz_opacity

contains


!----------------------------------------------------------------------
! build_kurucz_splines
!----------------------------------------------------------------------
! Builds, for each T row of tbl, the natural cubic spline of log kappa in
! log rho over the row's densities with data, as coefficients per
! interval (tbl%spline_coeffs), with the row's first density index and
! count. A row that does not start at the first density, or has fewer
! than 25 points, returns ierr = 1.
subroutine build_kurucz_splines(tbl, ierr)
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(out) :: ierr
! (np: the coefficient array's interval count, more than any row's)
      integer, parameter :: np = 100

      real(dp) :: spline_work(4,np), density_nodes(kurucz_max_num_densities), &
           row_log10_opacity(kurucz_max_num_densities), &
           row_second_derivs(kurucz_max_num_densities)
      integer :: it, jd, id, index2, j, i
      real(dp) :: chkd, chko
      character(len=message_length) :: lines(1)

      ierr = 0
      do it = 1,tbl%num_temps
         jd = 0
         do id = 1,kurucz_max_num_densities
            chkd = tbl%log10_rho(it,id)
            chko = tbl%opacity(it,id)
            if (chko <= 0.0d0) cycle
            if (jd <= 0) tbl%density_start_index(it) = id
            jd = jd + 1
            density_nodes(jd) = chkd
            row_log10_opacity(jd) = log10(chko)
         end do
         tbl%density_count(it) = jd
         if (tbl%density_start_index(it) /= 1) then
            lines = unwritten_line
            write(lines,*) 'build_kurucz_splines: kurucz table density grid does not start at index 1'
            call kap_terminal(lines)
            ierr = 1
            return
         end if
         if (tbl%density_count(it) < 25) then
            lines = unwritten_line
            write(lines,*) 'build_kurucz_splines: kurucz table has fewer than 25 density points'
            call kap_terminal(lines)
            ierr = 1
            return
         end if
! the row's natural cubic spline in log rho, as coefficients per interval
         call natural_spline(density_nodes, row_log10_opacity, jd, row_second_derivs)
         call compute_spline_coefficients(density_nodes, row_log10_opacity, &
              row_second_derivs, jd, spline_work)
         do j = 1,jd
            do i = 1,4
               index2 = i + (j-1)*4
               tbl%spline_coeffs(it,index2) = spline_work(i,j)
            end do
         end do
      end do

end subroutine build_kurucz_splines

!----------------------------------------------------------------------
! compute_kurucz_opacity (formerly kurucz)
!----------------------------------------------------------------------
! The opacity at (log rho, log T) and its log derivatives from table set
! tbl (the first or the second Z's), by cubic splines: in log rho on each
! T row, then in log T across the rows. in_table is .false. for a point
! outside the table (and, for the first table, above log rho = -3 or
! log T = 4.1), where the caller uses the interior opacity alone.
subroutine compute_kurucz_opacity(tbl, log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, in_table, ierr)

!     Two dimensional interpolation for opacity.
!     tbl%ix_t / tbl%ix_rho are the warm-start
!     cursors (nearest temperature row and density knot) left by the
!     previous call.  For up to six temperature rows around the point
!     the row's cubic spline in log rho gives log opacity and its rho
!     derivative; those are then interpolated in log T (natural_spline and splint) to give
!     log10_opacity, dlnkap_dlnt and dlnkap_dlnrho.

      integer :: jerr_gate

      type(kurucz_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      logical, intent(out) :: in_table
      integer, intent(out) :: ierr

      real(dp) :: temp_subset_logt(kurucz_max_num_temps), &
           temp_subset_log10_opacity(kurucz_max_num_temps)
      real(dp) :: temp_subset_dlnkap_dlnrho(kurucz_max_num_temps), &
           temp_subset_second_derivs(kurucz_max_num_temps)
      ! The Kurucz tables are only used below this density and temperature.
      real(dp), parameter :: kurucz_logrho_max = -3.0d0, kurucz_logt_max = 4.1d0
      integer :: num_valid_temps, temp_index, temp_index_start, &
           temp_index_end
      logical :: search_full_range
      logical :: density_found
! --- locals ---
      integer :: row_density_start, row_density_end, &
           density_pointer, density_scan_index, knot_index, &
           coeff_base_index
      real(dp) :: delta_logrho, c1, c2, c3, c4
      real(dp) :: log10_opacity_at_rho, dlnkap_dlnrho_at_rho
      real(dp) :: log10_opacity_interp
      character(len=message_length) :: lines(1)

      call findex(tbl%grid_logt, tbl%num_temps, log10_temperature, &
           tbl%ix_t)

      ierr = 0
      in_table = .false.
      if (tbl%check_range) then
      if (.not.(log10_density <= kurucz_logrho_max .and. &
           log10_temperature <= kurucz_logt_max)) return
      end if
      search_full_range = .true.
!     For six grid points of temperature
      temp_index_start = tbl%ix_t-2
      if (temp_index_start <= 0) temp_index_start = 1
      temp_index_end = tbl%ix_t+3
      if (temp_index_end > tbl%num_temps) temp_index_end = tbl%num_temps
      full_search: do
      num_valid_temps = 0
      do temp_index = temp_index_start, temp_index_end
         density_found = .false.
         row_density_start = tbl%density_start_index(temp_index)
         row_density_end = row_density_start + &
              tbl%density_count(temp_index) - 1
!        the density knot below log rho on this row
         if (tbl%ix_rho < row_density_start .or. &
              tbl%ix_rho > row_density_end) tbl%ix_rho = row_density_end
         density_pointer = tbl%ix_rho
         if (log10_density < tbl%log10_rho(temp_index, density_pointer)) then
            do density_scan_index = density_pointer-1, row_density_start, -1
               if (tbl%log10_rho(temp_index, density_scan_index) <= &
                    log10_density) then
                  density_pointer = density_scan_index
                  density_found = .true.
                  exit
               end if
            end do
         else
            do density_scan_index = density_pointer, row_density_end-1
               if (tbl%log10_rho(temp_index, density_scan_index+1) > &
                    log10_density) then
                  density_pointer = density_scan_index
                  density_found = .true.
                  exit
               end if
            end do
            if (.not. density_found) then
            if (tbl%log10_rho(temp_index, row_density_end) >= &
                 log10_density) then
               density_pointer = row_density_end
               density_found = .true.
            end if
            end if
         end if
         if (.not. density_found) cycle
         tbl%ix_rho = density_pointer
         knot_index = tbl%ix_rho - row_density_start + 1
         coeff_base_index = 4*(knot_index-1)
!        Now, (knot,knot+1) is sub-range of rho which contains d.
         delta_logrho = log10_density - tbl%log10_rho(temp_index, tbl%ix_rho)
         c1 = tbl%spline_coeffs(temp_index, coeff_base_index+1)
         c2 = tbl%spline_coeffs(temp_index, coeff_base_index+2)
         c3 = tbl%spline_coeffs(temp_index, coeff_base_index+3)
         c4 = tbl%spline_coeffs(temp_index, coeff_base_index+4)
!        Interpolation for opacity(ol) in the entry d and the each T-grid
!        estimates the partial derivative of ol wrt d
!        evaluates the interpolation value in the sub-range we determined.
         log10_opacity_at_rho = ((c4*delta_logrho+c3)*delta_logrho+c2)* &
              delta_logrho+c1
         dlnkap_dlnrho_at_rho = (3.0d0*c4*delta_logrho+2.0d0*c3)*delta_logrho+c2
         num_valid_temps = num_valid_temps+1
         temp_subset_logt(num_valid_temps) = tbl%grid_logt(temp_index)
         temp_subset_log10_opacity(num_valid_temps) = log10_opacity_at_rho
         temp_subset_dlnkap_dlnrho(num_valid_temps) = dlnkap_dlnrho_at_rho
      end do
      if (num_valid_temps <= 3) then
         if (search_full_range) then
            temp_index_start = 1
            temp_index_end = tbl%num_temps
            search_full_range = .false.
            cycle full_search
         end if
         call kap_log_values([log10_density, log10_temperature], &
              'ERROR KURUCZ OP: NO TABLE VALUE ')
         ierr = 1
         return
      end if
      exit full_search
      end do full_search
      if (temp_subset_logt(1) > log10_temperature .or. &
           temp_subset_logt(num_valid_temps) < log10_temperature) return
      in_table = .true.
!     Interpolation for the opacity in the entry T and d.
!     Get the partial derivative of ol wrt T.
      call natural_spline(temp_subset_logt, temp_subset_log10_opacity, &
           num_valid_temps, temp_subset_second_derivs)
      call splint(temp_subset_logt, temp_subset_log10_opacity, num_valid_temps, &
           temp_subset_second_derivs, log10_temperature, log10_opacity_interp, &
           jerr_gate, dy=dlnkap_dlnt)
      if (jerr_gate /= 0) then
         lines = unwritten_line
         write(lines,*) 'splint: '//numerics_status_text(jerr_gate)
         call kap_log(lines)
         ierr = jerr_gate
         return
      end if
      log10_opacity = log10_opacity_interp
      opacity = exp10(log10_opacity)
!     (dlnkap_dlnt from the same spline's slope)
!     find the partial derivative value of ol wrt d in the given T and d
      call natural_spline(temp_subset_logt, temp_subset_dlnkap_dlnrho, &
           num_valid_temps, temp_subset_second_derivs)
      call splint(temp_subset_logt, temp_subset_dlnkap_dlnrho, num_valid_temps, &
           temp_subset_second_derivs, log10_temperature, dlnkap_dlnrho, jerr_gate)
      if (jerr_gate /= 0) then
         lines = unwritten_line
         write(lines,*) 'splint: '//numerics_status_text(jerr_gate)
         call kap_log(lines)
         ierr = jerr_gate
         return
      end if
end subroutine compute_kurucz_opacity

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

subroutine load_kurucz90(settings, paths, envelope_hydrogen_fraction, tables, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      integer, intent(out) :: ierr

      call read_kurucz_tables(tables%kurucz, paths%kurucz, paths%kurucz_z2, &
           settings%use_two_z_tables, ierr)
      if (ierr /= 0) return
      call build_kurucz_splines(tables%kurucz(1), ierr)
      if (ierr /= 0) return
      if (settings%use_two_z_tables) call build_kurucz_splines(tables%kurucz(2), ierr)
end subroutine load_kurucz90

subroutine lookup_kurucz90(settings, tables, point, res, in_table, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      integer, intent(out) :: ierr
      type(kap_result) :: res_z2

      call compute_kurucz_opacity(tables%kurucz(1), point%log10_density, &
           point%log10_temperature, res%opacity, res%log10_opacity, &
           res%dlnkap_dlnrho, res%dlnkap_dlnt, in_table, ierr)
      if (ierr /= 0) return
      if (in_table .and. settings%use_two_z_tables) then
         call compute_kurucz_opacity(tables%kurucz(2), point%log10_density, &
              point%log10_temperature, res_z2%opacity, res_z2%log10_opacity, &
              res_z2%dlnkap_dlnrho, res_z2%dlnkap_dlnt, in_table, ierr)
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
