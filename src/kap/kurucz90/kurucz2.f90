!----------------------------------------------------------------------
! kurucz2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original kurucz2.f; only variable names, source form, and comment
! style were updated.
!
! Interpolation facility for Kurucz's opacity tables using a cubic
! spline interpolation scheme -- second Z table (see kurucz.f90).
! Uses the classic F77 "alternate return" (the trailing `*` dummy
! argument / `return 1`) to signal "point outside the table" to the
! caller (see kap_lib.f90's kap_get, label 100). Unlike kurucz.f90, this
! routine does NOT pre-check log10_density/log10_temperature against
! the -3.0/4.1 bounds before searching -- that early-out check is
! absent from the original kurucz2.f and is preserved as-is here.
subroutine kurucz2(log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr, *)

!     TWO DIMENSIONAL INTERPOLATION FOR OPACITY
!     M1, M2, AND M3 ARE NEAREST GRID POINT OF ABUNDANCE, TEMPERATURE,
!     AND DENSITY WE GOT AT 'GETD'.
!     O IS OPACITY.
!     OL IS log(O)
!     QODF IS THE PARTIAL DERIVATIVE OF O WRT D.
!     QOTF IS THE PARTIAL DERIVATIVE OF O WRT T.
!     THE VARIABLES WHICH HAVE 'TEM'-HEAD OR 'TEM'-TAIL ARE TEMPORARILLY
!     REQUIRED FOR DERIVATIVE VALUE ESTIMATION.  THESE ARE NOT THE VALUES
!     WE WANT.
!     THE VARIABLES WHICH HAVE 'I'-TAIL ARE OBTAINED FROM INTERPOLATION
!     ROUTINE.  WE WILL NOT USE THESE EVALUATIONS.
!     THE VARIABLES WHICH HAVE 'F'-TAIL ARE WHAT WILL BE RETURNED.
!     THIS VALUES IS ESTIMATED BY NUMERICAL DIFFERENTIATION.

      use opacity_table_lib
      use luout_lib
      use numerics_lib
      use math_lib
      implicit none
      integer :: jerr_gate
      integer, parameter :: max_num_temps = 60
      integer, parameter :: max_num_densities = 50
      integer, parameter :: num_x_tables = 1
      integer, parameter :: num_x_temp_entries = max_num_temps*num_x_tables
      integer, parameter :: num_spline_coeffs = 4*max_num_densities

      double precision, intent(in) :: log10_density, log10_temperature
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      double precision :: temp_subset_logt(max_num_temps), &
           temp_subset_log10_opacity(max_num_temps)
      double precision :: temp_subset_dlnkap_dlnrho(max_num_temps)
      integer :: num_valid_temps, temp_index, temp_index_start, &
           temp_index_end
      logical :: search_full_range
      logical :: density_found






! kurucz2_ix_x/kurucz2_ix_t/kurucz2_ix_rho defaults moved to
! opacity_table_lib.f90: DATA can no longer target them here now that
! they're use-associated.
! --- locals ---
      double precision :: local_logt, local_logrho
      integer :: t_row_index, row_density_start, row_density_end, &
           density_pointer, density_scan_index, knot_index, &
           coeff_base_index
      double precision :: delta_logrho, c1, c2, c3, c4
      double precision :: log10_opacity_at_rho, dlnkap_dlnrho_at_rho
      double precision :: log10_opacity_interp, dlnkap_dlnrho_unused

      call findex(opacity_table%kurucz2_grid_logt, opacity_table%kurucz2_num_temps, log10_temperature, &
           opacity_table%kurucz2_ix_t)

      ierr = 0
      local_logt = log10_temperature
      local_logrho = log10_density
      search_full_range = .true.
!     FOR SIX GRID POINTS OF TEMPERATURE
      temp_index_start = opacity_table%kurucz2_ix_t-2
      if (temp_index_start.le.0) temp_index_start = 1
      temp_index_end = opacity_table%kurucz2_ix_t+3
      if (temp_index_end.gt.opacity_table%kurucz2_num_temps) temp_index_end = opacity_table%kurucz2_num_temps
      full_search: do
      num_valid_temps = 0
      do temp_index = temp_index_start, temp_index_end
         t_row_index = temp_index
         density_found = .false.
         row_density_start = opacity_table%kurucz2_density_start_index(t_row_index)
         row_density_end = row_density_start + &
              opacity_table%kurucz2_density_count(t_row_index) - 1
!        FIND THE 'M3'
         if (opacity_table%kurucz2_ix_rho.lt.row_density_start .or. &
              opacity_table%kurucz2_ix_rho.gt.row_density_end) opacity_table%kurucz2_ix_rho = row_density_end
         density_pointer = opacity_table%kurucz2_ix_rho
         if (local_logrho.lt.opacity_table%kurucz2_log10_rho(temp_index, density_pointer)) then
            do density_scan_index = density_pointer-1, row_density_start, -1
               if (opacity_table%kurucz2_log10_rho(temp_index, density_scan_index).le. &
                    local_logrho) then
                  density_pointer = density_scan_index
                  density_found = .true.
                  exit
               endif
            end do
         else
            do density_scan_index = density_pointer, row_density_end-1
               if (opacity_table%kurucz2_log10_rho(temp_index, density_scan_index+1).gt. &
                    local_logrho) then
                  density_pointer = density_scan_index
                  density_found = .true.
                  exit
               endif
            end do
            if (.not. density_found) then
            if (opacity_table%kurucz2_log10_rho(temp_index, row_density_end).ge. &
                 local_logrho) then
               density_pointer = row_density_end
               density_found = .true.
            endif
            endif
         endif
         if (.not. density_found) cycle
         opacity_table%kurucz2_ix_rho = density_pointer
         knot_index = opacity_table%kurucz2_ix_rho - row_density_start + 1
         coeff_base_index = 4*(knot_index-1)
!        NOW, (KNOT,KNOT+1) IS SUB-RANGE OF RHO WHICH CONTAINS D.
         delta_logrho = local_logrho - opacity_table%kurucz2_log10_rho(temp_index, opacity_table%kurucz2_ix_rho)
         c1 = opacity_table%kurucz2_spline_coeffs(t_row_index, coeff_base_index+1)
         c2 = opacity_table%kurucz2_spline_coeffs(t_row_index, coeff_base_index+2)
         c3 = opacity_table%kurucz2_spline_coeffs(t_row_index, coeff_base_index+3)
         c4 = opacity_table%kurucz2_spline_coeffs(t_row_index, coeff_base_index+4)
!        INTERPOLATION FOR OPACITY(OL) IN THE ENTRY D AND THE EACH T-GRID
!        ESTIMATES THE PARTIAL DERIVATIVE OF OL WRT D
!        EVALUATES THE INTERPOLATION VALUE IN THE SUB-RANGE WE DETERMINED.
         log10_opacity_at_rho = ((c4*delta_logrho+c3)*delta_logrho+c2)* &
              delta_logrho+c1
         dlnkap_dlnrho_at_rho = (3.0d0*c4*delta_logrho+2.0d0*c3)*delta_logrho+c2
         num_valid_temps = num_valid_temps+1
         temp_subset_logt(num_valid_temps) = opacity_table%kurucz2_grid_logt(temp_index)
         temp_subset_log10_opacity(num_valid_temps) = log10_opacity_at_rho
         temp_subset_dlnkap_dlnrho(num_valid_temps) = dlnkap_dlnrho_at_rho
      end do
      if (num_valid_temps.le.3) then
         if (search_full_range) then
            temp_index_start = 1
            temp_index_end = opacity_table%kurucz2_num_temps
            search_full_range = .false.
            cycle full_search
         endif
         write(run_log_unit,*) 'ERROR KURUCZ OP: NO TABLE VALUE ', &
              local_logrho, local_logt
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      endif
      exit full_search
      end do full_search
      if (temp_subset_logt(1).gt.local_logt .or. &
           temp_subset_logt(num_valid_temps).lt.local_logt) return 1
!     INTERPOLATION FOR THE OPACITY IN THE ENTRY T AND D.
!     GET THE PARTIAL DERIVATIVE OF OL WRT T.
      call intpol(temp_subset_logt, temp_subset_log10_opacity, &
           num_valid_temps, local_logt, log10_opacity_interp, dlnkap_dlnt, jerr_gate)
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
           num_valid_temps, local_logt, dlnkap_dlnrho, dlnkap_dlnrho_unused, jerr_gate)
      if (jerr_gate /= 0) then
         ierr = jerr_gate
         return
      end if
!     QODF = D LN(O)/D LN(D)
      return
end subroutine kurucz2
