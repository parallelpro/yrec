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
! the caller (see getopac.f90, label 100).
subroutine kurucz(log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, *)

!     TWO DIMENSIONAL INTERPOLATION FOR OPACITY
!     M1, M2, AND M3 ARE NEAREST GRID POINT OF ABUNDANCE, TEMPERATURE,
!     AND DENSITY WE GOT AT 'GETD'.
!     O IS OPACITY.
!     OL IS DLOG(O)
!     QODF IS THE PARTIAL DERIVATIVE OF O WRT D.
!     QOTF IS THE PARTIAL DERIVATIVE OF O WRT T.
!     THE VARIABLES WHICH HAVE 'TEM'-HEAD OR 'TEM'-TAIL ARE TEMPORARILLY
!     REQUIRED FOR DERIVATIVE VALUE ESTIMATION.  THESE ARE NOT THE VALUES
!     WE WANT.
!     THE VARIABLES WHICH HAVE 'I'-TAIL ARE OBTAINED FROM INTERPOLATION
!     ROUTINE.  WE WILL NOT USE THESE EVALUATIONS.
!     THE VARIABLES WHICH HAVE 'F'-TAIL ARE WHAT WILL BE RETURNED.
!     THIS VALUES IS ESTIMATED BY NUMERICAL DIFFERENTIATION.

      use findex_mod
      implicit none
      integer, parameter :: max_num_temps = 60
      integer, parameter :: max_num_densities = 50
      integer, parameter :: num_x_tables = 1
      integer, parameter :: num_x_temp_entries = max_num_temps*num_x_tables
      integer, parameter :: num_spline_coeffs = 4*max_num_densities

      double precision, intent(in) :: log10_density, log10_temperature
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      double precision :: temp_subset_logt(max_num_temps), &
           temp_subset_log10_opacity(max_num_temps)
      double precision :: temp_subset_dlnkap_dlnrho(max_num_temps)
      integer :: num_valid_temps, temp_index, temp_index_start, &
           temp_index_end
      logical :: search_full_range

      double precision :: kurucz_grid_logt(max_num_temps)
      common /gkrz/ kurucz_grid_logt

      double precision :: kurucz_log10_opacity(num_x_temp_entries, max_num_densities), &
           kurucz_log10_rho(num_x_temp_entries, max_num_densities)
      integer :: kurucz_num_temps
      common /krz/ kurucz_log10_opacity, kurucz_log10_rho, kurucz_num_temps

      integer :: kurucz_ix_x, kurucz_ix_t, kurucz_ix_rho
      common /kipm/ kurucz_ix_x, kurucz_ix_t, kurucz_ix_rho

      double precision :: kurucz_spline_coeffs(num_x_temp_entries, num_spline_coeffs)
      integer :: kurucz_density_start_index(num_x_temp_entries), &
           kurucz_density_count(num_x_temp_entries)
      common /intpl2/ kurucz_spline_coeffs, kurucz_density_start_index, &
           kurucz_density_count

! common/luout/: only short_file_unit is used here.
      integer :: ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr

      data kurucz_ix_x, kurucz_ix_t, kurucz_ix_rho /1, 1, 1/
      save

! --- locals ---
      double precision :: local_logt, local_logrho
      integer :: t_row_index, row_density_start, row_density_end, &
           density_pointer, density_scan_index, knot_index, &
           coeff_base_index
      double precision :: delta_logrho, c1, c2, c3, c4
      double precision :: log10_opacity_at_rho, dlnkap_dlnrho_at_rho
      double precision :: log10_opacity_interp, dlnkap_dlnrho_unused

      call findex(kurucz_grid_logt, kurucz_num_temps, log10_temperature, &
           kurucz_ix_t)

      local_logt = log10_temperature
      local_logrho = log10_density
      if (.not.(log10_density.le.-3.0d0 .and. &
           log10_temperature.le.4.1d0)) return 1
      search_full_range = .true.
!     FOR SIX GRID POINTS OF TEMPERATURE
      temp_index_start = kurucz_ix_t-2
      if (temp_index_start.le.0) temp_index_start = 1
      temp_index_end = kurucz_ix_t+3
      if (temp_index_end.gt.kurucz_num_temps) temp_index_end = kurucz_num_temps
  333 continue
      num_valid_temps = 0
      do 300 temp_index = temp_index_start, temp_index_end
         t_row_index = temp_index
         row_density_start = kurucz_density_start_index(t_row_index)
         row_density_end = row_density_start + &
              kurucz_density_count(t_row_index) - 1
!        FIND THE 'M3'
         if (kurucz_ix_rho.lt.row_density_start .or. &
              kurucz_ix_rho.gt.row_density_end) kurucz_ix_rho = row_density_end
         density_pointer = kurucz_ix_rho
         if (local_logrho.lt.kurucz_log10_rho(temp_index, density_pointer)) then
            do 211 density_scan_index = density_pointer-1, row_density_start, -1
               if (kurucz_log10_rho(temp_index, density_scan_index).le. &
                    local_logrho) then
                  density_pointer = density_scan_index
                  goto 213
               endif
  211       continue
            go to 300
         else
            do 212 density_scan_index = density_pointer, row_density_end-1
               if (kurucz_log10_rho(temp_index, density_scan_index+1).gt. &
                    local_logrho) then
                  density_pointer = density_scan_index
                  goto 213
               endif
  212       continue
            if (kurucz_log10_rho(temp_index, row_density_end).ge. &
                 local_logrho) then
               density_pointer = row_density_end
               goto 213
            endif
            go to 300
         endif
  213    kurucz_ix_rho = density_pointer
         knot_index = kurucz_ix_rho - row_density_start + 1
         coeff_base_index = 4*(knot_index-1)
!        NOW, (KNOT,KNOT+1) IS SUB-RANGE OF RHO WHICH CONTAINS D.
         delta_logrho = local_logrho - kurucz_log10_rho(temp_index, kurucz_ix_rho)
         c1 = kurucz_spline_coeffs(t_row_index, coeff_base_index+1)
         c2 = kurucz_spline_coeffs(t_row_index, coeff_base_index+2)
         c3 = kurucz_spline_coeffs(t_row_index, coeff_base_index+3)
         c4 = kurucz_spline_coeffs(t_row_index, coeff_base_index+4)
!        INTERPOLATION FOR OPACITY(OL) IN THE ENTRY D AND THE EACH T-GRID
!        ESTIMATES THE PARTIAL DERIVATIVE OF OL WRT D
!        EVALUATES THE INTERPOLATION VALUE IN THE SUB-RANGE WE DETERMINED.
         log10_opacity_at_rho = ((c4*delta_logrho+c3)*delta_logrho+c2)* &
              delta_logrho+c1
         dlnkap_dlnrho_at_rho = (3.0d0*c4*delta_logrho+2.0d0*c3)*delta_logrho+c2
         num_valid_temps = num_valid_temps+1
         temp_subset_logt(num_valid_temps) = kurucz_grid_logt(temp_index)
         temp_subset_log10_opacity(num_valid_temps) = log10_opacity_at_rho
         temp_subset_dlnkap_dlnrho(num_valid_temps) = dlnkap_dlnrho_at_rho
  300 continue
      if (num_valid_temps.le.3) then
         if (search_full_range) then
            temp_index_start = 1
            temp_index_end = kurucz_num_temps
            search_full_range = .false.
            go to 333
         endif
         write(short_file_unit,*) 'ERROR KURUCZ OP: NO TABLE VALUE ', &
              local_logrho, local_logt
         stop
      endif
      if (temp_subset_logt(1).gt.local_logt .or. &
           temp_subset_logt(num_valid_temps).lt.local_logt) return 1
!     INTERPOLATION FOR THE OPACITY IN THE ENTRY T AND D.
!     GET THE PARTIAL DERIVATIVE OF OL WRT T.
      call intpol(temp_subset_logt, temp_subset_log10_opacity, &
           num_valid_temps, local_logt, log10_opacity_interp, dlnkap_dlnt)
      log10_opacity = log10_opacity_interp
      opacity = 10.0d0**log10_opacity
!     QOTF = D LN(O)/D LN(T)
!     FIND THE PARTIAL DERIVATIVE VALUE OF OL WRT D IN THE GIVEN T AND D
      call intpol(temp_subset_logt, temp_subset_dlnkap_dlnrho, &
           num_valid_temps, local_logt, dlnkap_dlnrho, dlnkap_dlnrho_unused)
!     QODF = D LN(O)/D LN(D)
      return
end subroutine kurucz
