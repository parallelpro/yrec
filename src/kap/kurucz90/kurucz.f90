!----------------------------------------------------------------------
! kurucz
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original kurucz.f; only variable names, source form, and comment
! style were updated.
!
! Interpolation facility for Kurucz's opacity tables using a cubic
! spline interpolation scheme (two-dimensional: temperature and
! density). Uses the classic F77 "alternate return" (the trailing `*`
! dummy argument / `return 1`) to signal "point outside the table" to
! the caller (see kap_lib.f90's kap_get, label 100).
subroutine kurucz(log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr, *)

!     TWO DIMENSIONAL INTERPOLATION FOR OPACITY.
!     opacity_table%kurucz_ix_t / kurucz_ix_rho ARE THE WARM-START
!     CURSORS (NEAREST TEMPERATURE ROW AND DENSITY KNOT) LEFT BY THE
!     PREVIOUS CALL.  FOR UP TO SIX TEMPERATURE ROWS AROUND THE POINT
!     THE ROW'S CUBIC SPLINE IN LOG RHO GIVES LOG OPACITY AND ITS RHO
!     DERIVATIVE; THOSE ARE THEN INTERPOLATED IN LOG T (intpol) TO GIVE
!     log10_opacity, dlnkap_dlnt AND dlnkap_dlnrho.

      use opacity_table_lib
      use luout_lib
      use numerics_lib
      use math_lib
      implicit none
      integer :: jerr_gate

      double precision, intent(in) :: log10_density, log10_temperature
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      double precision :: temp_subset_logt(kurucz_max_num_temps), &
           temp_subset_log10_opacity(kurucz_max_num_temps)
      double precision :: temp_subset_dlnkap_dlnrho(kurucz_max_num_temps)
      ! The Kurucz tables are only used below this density and temperature.
      double precision, parameter :: kurucz_logrho_max = -3.0d0, kurucz_logt_max = 4.1d0
      integer :: num_valid_temps, temp_index, temp_index_start, &
           temp_index_end
      logical :: search_full_range
      logical :: density_found
! --- locals ---
      integer :: row_density_start, row_density_end, &
           density_pointer, density_scan_index, knot_index, &
           coeff_base_index
      double precision :: delta_logrho, c1, c2, c3, c4
      double precision :: log10_opacity_at_rho, dlnkap_dlnrho_at_rho
      double precision :: log10_opacity_interp, dlnkap_dlnrho_unused

      call findex(opacity_table%kurucz_grid_logt, opacity_table%kurucz_num_temps, log10_temperature, &
           opacity_table%kurucz_ix_t)

      ierr = 0
      if (.not.(log10_density.le.kurucz_logrho_max .and. &
           log10_temperature.le.kurucz_logt_max)) return 1
      search_full_range = .true.
!     FOR SIX GRID POINTS OF TEMPERATURE
      temp_index_start = opacity_table%kurucz_ix_t-2
      if (temp_index_start.le.0) temp_index_start = 1
      temp_index_end = opacity_table%kurucz_ix_t+3
      if (temp_index_end.gt.opacity_table%kurucz_num_temps) temp_index_end = opacity_table%kurucz_num_temps
      full_search: do
      num_valid_temps = 0
      do temp_index = temp_index_start, temp_index_end
         density_found = .false.
         row_density_start = opacity_table%kurucz_density_start_index(temp_index)
         row_density_end = row_density_start + &
              opacity_table%kurucz_density_count(temp_index) - 1
!        FIND THE 'M3'
         if (opacity_table%kurucz_ix_rho.lt.row_density_start .or. &
              opacity_table%kurucz_ix_rho.gt.row_density_end) opacity_table%kurucz_ix_rho = row_density_end
         density_pointer = opacity_table%kurucz_ix_rho
         if (log10_density.lt.opacity_table%kurucz_log10_rho(temp_index, density_pointer)) then
            do density_scan_index = density_pointer-1, row_density_start, -1
               if (opacity_table%kurucz_log10_rho(temp_index, density_scan_index).le. &
                    log10_density) then
                  density_pointer = density_scan_index
                  density_found = .true.
                  exit
               endif
            end do
         else
            do density_scan_index = density_pointer, row_density_end-1
               if (opacity_table%kurucz_log10_rho(temp_index, density_scan_index+1).gt. &
                    log10_density) then
                  density_pointer = density_scan_index
                  density_found = .true.
                  exit
               endif
            end do
            if (.not. density_found) then
            if (opacity_table%kurucz_log10_rho(temp_index, row_density_end).ge. &
                 log10_density) then
               density_pointer = row_density_end
               density_found = .true.
            endif
            endif
         endif
         if (.not. density_found) cycle
         opacity_table%kurucz_ix_rho = density_pointer
         knot_index = opacity_table%kurucz_ix_rho - row_density_start + 1
         coeff_base_index = 4*(knot_index-1)
!        NOW, (KNOT,KNOT+1) IS SUB-RANGE OF RHO WHICH CONTAINS D.
         delta_logrho = log10_density - opacity_table%kurucz_log10_rho(temp_index, opacity_table%kurucz_ix_rho)
         c1 = opacity_table%kurucz_spline_coeffs(temp_index, coeff_base_index+1)
         c2 = opacity_table%kurucz_spline_coeffs(temp_index, coeff_base_index+2)
         c3 = opacity_table%kurucz_spline_coeffs(temp_index, coeff_base_index+3)
         c4 = opacity_table%kurucz_spline_coeffs(temp_index, coeff_base_index+4)
!        INTERPOLATION FOR OPACITY(OL) IN THE ENTRY D AND THE EACH T-GRID
!        ESTIMATES THE PARTIAL DERIVATIVE OF OL WRT D
!        EVALUATES THE INTERPOLATION VALUE IN THE SUB-RANGE WE DETERMINED.
         log10_opacity_at_rho = ((c4*delta_logrho+c3)*delta_logrho+c2)* &
              delta_logrho+c1
         dlnkap_dlnrho_at_rho = (3.0d0*c4*delta_logrho+2.0d0*c3)*delta_logrho+c2
         num_valid_temps = num_valid_temps+1
         temp_subset_logt(num_valid_temps) = opacity_table%kurucz_grid_logt(temp_index)
         temp_subset_log10_opacity(num_valid_temps) = log10_opacity_at_rho
         temp_subset_dlnkap_dlnrho(num_valid_temps) = dlnkap_dlnrho_at_rho
      end do
      if (num_valid_temps.le.3) then
         if (search_full_range) then
            temp_index_start = 1
            temp_index_end = opacity_table%kurucz_num_temps
            search_full_range = .false.
            cycle full_search
         endif
         write(run_log_unit,*) 'ERROR KURUCZ OP: NO TABLE VALUE ', &
              log10_density, log10_temperature
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      endif
      exit full_search
      end do full_search
      if (temp_subset_logt(1).gt.log10_temperature .or. &
           temp_subset_logt(num_valid_temps).lt.log10_temperature) return 1
!     INTERPOLATION FOR THE OPACITY IN THE ENTRY T AND D.
!     GET THE PARTIAL DERIVATIVE OF OL WRT T.
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
!     FIND THE PARTIAL DERIVATIVE VALUE OF OL WRT D IN THE GIVEN T AND D
      call intpol(temp_subset_logt, temp_subset_dlnkap_dlnrho, &
           num_valid_temps, log10_temperature, dlnkap_dlnrho, dlnkap_dlnrho_unused, jerr_gate)
      if (jerr_gate /= 0) then
         ierr = jerr_gate
         return
      end if
!     QODF = D LN(O)/D LN(D)
      return
end subroutine kurucz
