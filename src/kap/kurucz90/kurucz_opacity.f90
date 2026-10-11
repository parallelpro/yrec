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
! The tables live in opacity_tables (opacity_table%kurucz(1:2)).
module kurucz_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kap_settings, kurucz_max_num_densities, &
           kurucz_max_num_temps, kurucz_table_set, opacity_table
      use numerics_lib, only: findex, intpol, ysplin
      use io_units, only: run_log_unit
      use math_lib
      implicit none
      private
      public :: read_kurucz_tables, compute_kurucz_opacity

contains

!----------------------------------------------------------------------
! read_kurucz_tables
!----------------------------------------------------------------------
! Reads the Kurucz opacity table (and, if a second Z table is
! requested, a second Kurucz table) into opacity_table%kurucz(1)/(2),
! then calls build_kurucz_splines on each for the spline
! interpolation coefficients used later by kurucz.f90.
!
! 2026 wave 3 (R5): the two verbatim read loops (first table on
! star%ctrl%kurucz_table_unit, second on kurucz_table2_unit) became
! the internal read_kurucz_table below, called once per set; the
! order read 1, read 2, splines 1, splines 2 is the original's. The
! loop's working scalars stay host-associated so that whatever value
! they carried from the first table into the second is carried
! exactly as before.
subroutine read_kurucz_tables(kurucz_table_path, kurucz_table2_path, ierr)

      character(len=256), intent(in) :: kurucz_table_path, kurucz_table2_path
      integer, intent(out) :: ierr
      integer :: num_read, density_index
      integer :: read_status
      real(dp) :: prev_grid_temp, grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col

      ierr = 0
! The early-out range check in kurucz() belongs to the first table
! only (the former kurucz2.f90 never had it).
      opacity_table%kurucz(2)%check_range = .false.

      call read_kurucz_table(kurucz_table_path, &
           opacity_table%kurucz(1), ierr)
      if (ierr /= 0) return

! DBG 12/95 read in second Z table if requested
      if (kap_settings%use_two_z_tables) then
         call read_kurucz_table(kurucz_table2_path, &
              opacity_table%kurucz(2), ierr)
         if (ierr /= 0) return
      end if

      call build_kurucz_splines(opacity_table%kurucz(1), ierr)
      if (ierr /= 0) return
! DBG 12/95 second Z table
      if (kap_settings%use_two_z_tables) then
         call build_kurucz_splines(opacity_table%kurucz(2), ierr)
      end if

      return

contains

! read_kurucz_table: read one Kurucz table file into tbl.
subroutine read_kurucz_table(path, tbl, ierr)
      character(len=256), intent(in) :: path
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      integer :: unit

    1 format(2f5.2,5f7.3,3f9.5,f8.3)
!     OPEN TABLE
      open(newunit=unit, file=path, status='OLD')
      read(unit,'(/)')

      num_read = 0
      prev_grid_temp = 0.0d0
      table_read: do
      read(unit,1,iostat=read_status) grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col
      if (read_status > 0) cycle table_read   ! was err=110: reread next record
      if (read_status < 0) exit table_read    ! was end=120
      if (prev_grid_temp /= grid_temp) then
         num_read = num_read+1
         if (num_read > kurucz_max_num_temps) then
            write(*,*) 'read_kurucz_tables: too many temperature rows in ', &
                 trim(path)
            ierr = 1
            return
         end if
         tbl%grid_logt(num_read) = grid_temp
         prev_grid_temp = grid_temp
         density_index = 1
      endif
      if (num_read < 1) then
         write(*,*) 'read_kurucz_tables: no header row read from ', &
              trim(path)
         ierr = 1
         return
      end if
      tbl%log10_rho(num_read, density_index) = density
      tbl%log10_opacity(num_read, density_index) = exp10(log10_opacity0)
      density_index = density_index+1
      end do table_read

      tbl%num_temps = num_read
!     Close the table we have read
      close(unit,err=99)
      return
   99 continue
      write(*,*) 'read_kurucz_tables: error closing table file'
      ierr = 1
end subroutine read_kurucz_table

end subroutine read_kurucz_tables

!----------------------------------------------------------------------
! build_kurucz_splines
!----------------------------------------------------------------------
! YCK 3/91. Builds the cubic-spline coefficients (in density) for
! one Kurucz90 opacity table set tbl. Reads the tbl%grid_logt/
! log10_rho/log10_opacity grids filled by read_kurucz_tables.f90 and
! writes the per-row spline coefficients, density start index and
! count used by kurucz.f90.
!
! 2026 wave 3 (R5): was one routine that did this for
! opacity_table%kurucz_* and then (if use_two_z_tables) repeated the
! identical loop for kurucz2_*; read_kurucz_tables now calls it once
! per table set in the same order. The two stdout diagnostics that
! said "kurucz table" / "second kurucz table" are now generic.
subroutine build_kurucz_splines(tbl, ierr)
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(out) :: ierr
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      real(dp) :: spline_work(4,np), density_nodes(kurucz_max_num_densities)
      integer :: it, jd, id, index2, j, i
      real(dp) :: chkd, chko

      ierr = 0
      do it = 1,tbl%num_temps
         jd = 0
         do id = 1,kurucz_max_num_densities
            chkd = tbl%log10_rho(it,id)
            chko = tbl%log10_opacity(it,id)
            if (chko <= 0.0d0) cycle
            if (jd <= 0) tbl%density_start_index(it) = id
            jd = jd + 1
            density_nodes(jd) = chkd
            spline_work(1,jd) = log10(chko)
         end do
         tbl%density_count(it) = jd
         if (tbl%density_start_index(it) /= 1) then
            write(*,*) 'build_kurucz_splines: kurucz table density grid does not start at index 1'
            ierr = 1
            return
         end if
         if (tbl%density_count(it) < 25) then
            write(*,*) 'build_kurucz_splines: kurucz table has fewer than 25 density points'
            ierr = 1
            return
         end if
         call ysplin(density_nodes, spline_work, jd)
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
! kurucz
! Interpolation facility for Kurucz's opacity tables using a cubic
! spline interpolation scheme (two-dimensional: temperature and
! density). in_table is .false. for a point outside the table, where the
! caller uses the interior opacity alone. (2026-10: this was an F77
! alternate return, `return 1` to the caller's label 100.)
!
! 2026 wave 3 (R5): tbl is the table set to interpolate in --
! opacity_table%kurucz(1) (the former kurucz.f90) or
! opacity_table%kurucz(2) (the former kurucz2.f90, deleted; it was a
! member-for-member rename of this file except for the early-out
! range check below, which it lacked and which tbl%check_range now
! gates).
subroutine compute_kurucz_opacity(tbl, log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, in_table, ierr)

!     Two dimensional interpolation for opacity.
!     tbl%ix_t / tbl%ix_rho ARE THE WARM-START
!     cursors (nearest temperature row and density knot) left by the
!     previous call.  For up to six temperature rows around the point
!     the row's cubic spline in log rho gives log opacity and its rho
!     derivative; those are then interpolated in log T (intpol) to give
!     log10_opacity, dlnkap_dlnt AND dlnkap_dlnrho.

      integer :: jerr_gate

      type(kurucz_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      logical, intent(out) :: in_table
      integer, intent(out) :: ierr

      real(dp) :: temp_subset_logt(kurucz_max_num_temps), &
           temp_subset_log10_opacity(kurucz_max_num_temps)
      real(dp) :: temp_subset_dlnkap_dlnrho(kurucz_max_num_temps)
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
      real(dp) :: log10_opacity_interp, dlnkap_dlnrho_unused

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
!        FIND THE 'M3'
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
               endif
            end do
         else
            do density_scan_index = density_pointer, row_density_end-1
               if (tbl%log10_rho(temp_index, density_scan_index+1) > &
                    log10_density) then
                  density_pointer = density_scan_index
                  density_found = .true.
                  exit
               endif
            end do
            if (.not. density_found) then
            if (tbl%log10_rho(temp_index, row_density_end) >= &
                 log10_density) then
               density_pointer = row_density_end
               density_found = .true.
            endif
            endif
         endif
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
         endif
         write(run_log_unit,*) 'ERROR KURUCZ OP: NO TABLE VALUE ', &
              log10_density, log10_temperature
! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      endif
      exit full_search
      end do full_search
      if (temp_subset_logt(1) > log10_temperature .or. &
           temp_subset_logt(num_valid_temps) < log10_temperature) return
      in_table = .true.
!     Interpolation for the opacity in the entry T and d.
!     Get the partial derivative of ol wrt T.
      call intpol(temp_subset_logt, temp_subset_log10_opacity, &
           num_valid_temps, log10_temperature, log10_opacity_interp, dlnkap_dlnt, jerr_gate)
      ! 2026 numerics-gate opt-in: interpolation failure returns via
      ! ierr (diagnostic printed at the gate) instead of stopping.
      if (jerr_gate /= 0) then
         ierr = jerr_gate
         return
      end if
      log10_opacity = log10_opacity_interp
      opacity = exp10(log10_opacity)
!     QOTF = D LN(O)/D LN(T)
!     find the partial derivative value of ol wrt d in the given T and d
      call intpol(temp_subset_logt, temp_subset_dlnkap_dlnrho, &
           num_valid_temps, log10_temperature, dlnkap_dlnrho, dlnkap_dlnrho_unused, jerr_gate)
      if (jerr_gate /= 0) then
         ierr = jerr_gate
         return
      end if
!     QODF = D LN(O)/D LN(D)
end subroutine compute_kurucz_opacity

end module kurucz_opacity
