!----------------------------------------------------------------------
! kurucz_atmosphere
!----------------------------------------------------------------------
! The surface boundary condition from the Kurucz (1993) and the
! Castelli & Kurucz (2003) tables of log P at (log Teff, log g):
!     read_surface_pressure_table   read a table (at atm_init)
!     interpolate_surface_pressure  log P at (log Teff, log g): cubic
!                                   splines in log Teff and in log g
!                                   (3-point or linear near the table's
!                                   high-gravity edge)
! The tables live in atm_tables (atm_table%kurucz_surface_p,
! atm_table%castelli_surface_p).
module kurucz_atmosphere
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use atm_tables, only: atm_table, atm_table_ng, atm_table_nt, &
           surface_p_table
      use io_units, only: run_log_unit, terminal_unit
      use numerics_lib, only: inter3, natural_spline, splint
      implicit none
      private
      public :: read_surface_pressure_table, interpolate_surface_pressure

contains

!----------------------------------------------------------------------
! interpolate_surface_pressure (formerly surface_p_interp)
!----------------------------------------------------------------------
! log P at (log Teff, log g) from a surface-pressure table: cubic
! splines through 4 points in log Teff and in log g (3-point Lagrangian
! or linear in log g near the table's high-gravity edge). The result,
! with log T = log Teff, goes to atm_table%atm_log10_pressure and
! atm_log10_temperature. A point more than one table step outside the
! table returns ierr = 1.
subroutine interpolate_surface_pressure(tbl, log10_teff, log10_gravity, print_flag, ierr)

      integer :: jerr_gate

      type(surface_p_table), intent(inout) :: tbl
      real(dp), intent(in) :: log10_teff, log10_gravity
      logical, intent(in) :: print_flag
! --- locals ---
      real(dp) :: teff_nodes(4), gravity_nodes(4), &
           pressure_at_nodes(4), unused_deriv(3), gravity_weights(3), &
           gravity_nodes_3pt(3), pressure_table_vals(4), &
           teff_spline_deriv(4), gravity_spline_deriv(4)
      real(dp) :: fx, interpolated_value
      integer :: row, row_base, node, k, kk, kkk

! Interpolates in temperature with a 4-point natural cubic spline
! (natural_spline/splint), and in gravity the same way if 4 or more points are
! available (3-point Lagrangian or linear at the table edge).  It will quit if the desired data point
! has Teff or log G more than one table point from the data.
!
! Check to ensure that data is within table.
      integer, intent(out) :: ierr

! Lower table edge, and the log Teff rows above which only 3 or 2 log g
! columns are tabulated (so 3-point Lagrangian or linear interpolation).
      real(dp), parameter :: table_logteff_min = 3.5d0, table_logg_min = -0.5d0
      real(dp), parameter :: logteff_three_point = 4.5d0, logteff_two_point = 4.55d0

      ierr = 0

      if (log10_teff < table_logteff_min .or. log10_gravity < table_logg_min) then
         write(terminal_unit,911) log10_teff, log10_gravity
         write(run_log_unit,911) log10_teff, log10_gravity
  911    format(1X,'DESIRED ATMOSPHERE OUTSIDE TABLE RANGE'/ &
              ' LOG TEFF',F10.6,' LOG G',F10.6/' RUN STOPPED')
         ierr = 1
         return
      endif
! Temperature interpolation factors.
      do row = 1,tbl%num_teff
         if (log10_teff <= tbl%teff(row)) exit
      end do
      if (row > tbl%num_teff) then
      row = tbl%num_teff
      end if
      row_base = max(1,row-2)
      row_base = min(tbl%num_teff-3,row_base)
      do k = 1,4
         teff_nodes(k) = tbl%teff(row_base+k-1)
      end do
! Gravity interpolation factors.
      do row = row_base,row_base+3
         node = row-row_base+1
! Check if 4 log values available - otherwise, use 3 point Lagrangian
! or linear interpolation.
         if (tbl%teff(row) > logteff_three_point) then
            if (tbl%teff(row) > logteff_two_point) then
! Linear interpolation
               fx = (log10_gravity-tbl%logg(atm_table_ng-1))/ &
                    (tbl%logg(atm_table_ng)-tbl%logg(atm_table_ng-1))
               pressure_at_nodes(node) = &
                    tbl%log10_pressure(row,atm_table_ng-1) + &
                    fx*(tbl%log10_pressure(row,atm_table_ng) - &
                    tbl%log10_pressure(row,atm_table_ng-1))
            else
! 3-point Lagrangian interpolation.
               do k = 1,3
                  gravity_nodes_3pt(k) = tbl%logg(atm_table_ng-3+k)
               end do
               call inter3(gravity_nodes_3pt, gravity_weights, &
                    unused_deriv, log10_gravity)
               pressure_at_nodes(node) = &
                    tbl%log10_pressure(row,atm_table_ng-2)*gravity_weights(1) + &
                    tbl%log10_pressure(row,atm_table_ng-1)*gravity_weights(2) + &
                    tbl%log10_pressure(row,atm_table_ng)*gravity_weights(3)
            endif
            cycle
         endif
         if (log10_gravity >= tbl%logg(atm_table_ng-1)) then
! Desired log G above second to top table log G - use top 4 log G values.
            do kk = 1,4
               gravity_nodes(kk) = tbl%logg(atm_table_ng-4+kk)
               pressure_table_vals(kk) = &
                    tbl%log10_pressure(row,atm_table_ng-4+kk)
            end do
            call natural_spline(gravity_nodes, pressure_table_vals, 4, gravity_spline_deriv)
            call splint(gravity_nodes, pressure_table_vals, 4, &
                 gravity_spline_deriv, log10_gravity, interpolated_value, jerr_gate)
            ! 2026 numerics-gate opt-in: interpolation failure returns via
            ! ierr (diagnostic printed at the gate) instead of stopping.
            if (jerr_gate /= 0) then
               ierr = jerr_gate
               return
            end if
            pressure_at_nodes(node) = interpolated_value
            cycle
         endif
! General case - find 4 nearest points in gravity that are in the table.
! G Somers, I changed NG to IMINMAX in the next line. This prevents the
! code from using -999 to interpolate in some instances.
         do k = tbl%gmax_index(row)-3,tbl%gmin_index(row),-1
            if (log10_gravity < tbl%logg(k+2) .and. &
                 log10_gravity >= tbl%logg(k+1)) then
               kk = max(tbl%gmin_index(row),k)
               kk = min(atm_table_ng-3,kk)
               do kkk = 1,4
                  gravity_nodes(kkk) = tbl%logg(kk+kkk-1)
                  pressure_table_vals(kkk) = &
                       tbl%log10_pressure(row,kk+kkk-1)
               end do
               call natural_spline(gravity_nodes, pressure_table_vals, 4, gravity_spline_deriv)
               call splint(gravity_nodes, pressure_table_vals, 4, &
                    gravity_spline_deriv, log10_gravity, interpolated_value, jerr_gate)
               if (jerr_gate /= 0) then
                  ierr = jerr_gate
                  return
               end if
               pressure_at_nodes(node) = interpolated_value
               exit
            endif
         end do
         if (k < tbl%gmin_index(row)) then
! Desired log G below 2nd table entry -use first 4 points.
         do k = 1,4
            gravity_nodes(k) = tbl%logg(k+tbl%gmin_index(row)-1)
            pressure_table_vals(k) = &
                 tbl%log10_pressure(row,k+tbl%gmin_index(row)-1)
         end do
         call natural_spline(gravity_nodes, pressure_table_vals, 4, gravity_spline_deriv)
         call splint(gravity_nodes, pressure_table_vals, 4, &
              gravity_spline_deriv, log10_gravity, interpolated_value, jerr_gate)
         if (jerr_gate /= 0) then
            ierr = jerr_gate
            return
         end if
         pressure_at_nodes(node) = interpolated_value
         end if
      end do
! Interpolate in temperature to find correct log P.
      call natural_spline(teff_nodes, pressure_at_nodes, 4, teff_spline_deriv)
      call splint(teff_nodes, pressure_at_nodes, 4, teff_spline_deriv, &
           log10_teff, interpolated_value, jerr_gate)
      if (jerr_gate /= 0) then
         ierr = jerr_gate
         return
      end if
      atm_table%atm_log10_pressure = interpolated_value
      atm_table%atm_log10_temperature = log10_teff
! Write out information to the model file.
      if (print_flag) then
        write(run_log_unit,70)
70      format('********PRESSURE AT T=TEFF INTERPOLATED FROM TABULATED' &
              ,  ' VALUES********')
        write(run_log_unit,71) log10_teff, atm_table%atm_log10_pressure
71      format(' ',20X,'LOG (Teff) =',F10.5,' LOG P =',F10.5)
      endif
end subroutine interpolate_surface_pressure


!----------------------------------------------------------------------
! read_surface_pressure_table
!----------------------------------------------------------------------
! Reads a surface-pressure table into tbl: the metal fraction of the
! header (read past), num_teff log Teff rows (atm_table_nt for the
! Kurucz table, atm_table_ntc for Castelli/Kurucz), the atm_table_ng
! log g columns and log P; then finds each row's valid log g range
! (G Somers 5/15: a row may end in -999 entries at high gravity).
subroutine read_surface_pressure_table(table_path, tbl, num_teff, ierr)
      character(len=*), intent(in) :: table_path
      type(surface_p_table), intent(inout) :: tbl
      integer, intent(in) :: num_teff
      integer, intent(out) :: ierr
      integer :: table_unit, teff_idx, logg_idx
      logical :: found_valid_pressure
      real(dp) :: kurucz_table_z

      ierr = 0
      call allocate_surface_p(tbl, num_teff)
! Open surface pressure table
      open(newunit=table_unit, file=table_path, status='OLD', iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'atm_init: cannot open the atmosphere table ', trim(table_path)
         ierr = 1
         return
      end if
! GET ABUNDANCE:
      read(table_unit,200) kurucz_table_z
  200   format(1x,10x,1pe16.8,/)
! Get values of log teff:
      if (num_teff == atm_table_nt) then
         read(table_unit,201) (tbl%teff(teff_idx),teff_idx=1,num_teff)
      else
         read(table_unit,206) (tbl%teff(teff_idx),teff_idx=1,num_teff)
      end if
  201   format(1p4e16.8,13(/1p4e16.8),/1pe16.8,/)
  206   format(1p4e16.8,18(/1p4e16.8),/)
! Get values of log G:
      read(table_unit,202) (tbl%logg(logg_idx),logg_idx=1,atm_table_ng)
  202   format(1p4e16.8,/1p4e16.8,/1p3e16.8,/)
! Get grid of log pressure values:
      do teff_idx=1,num_teff
        read(table_unit,203) (tbl%log10_pressure(teff_idx,logg_idx),logg_idx=1,atm_table_ng)
  203     format(1p4e16.8,/1p4e16.8,/1p3e16.8)
      end do
      rewind table_unit
      close(table_unit)
!       G Somers 5/15; added a catch that allows the code to work
!                      if the highest gravity P value for a specific
!                      temperature is -999 (find_gravity_edges).
      call find_gravity_edges(tbl)

contains

! allocate_surface_p: on first use give tbl num_teff rows x
! atm_table_ng columns, zero-filled like the static arrays it
! replaces (an all -999 row leaves gmax_index untouched, so it must
! read 0 as before). Later calls (a re-read of the same table) reuse
! the storage as the static arrays did.
pure subroutine allocate_surface_p(tbl, num_teff)
      type(surface_p_table), intent(inout) :: tbl
      integer, intent(in) :: num_teff
      if (allocated(tbl%teff)) return
      allocate(tbl%teff(num_teff), tbl%logg(atm_table_ng), &
           tbl%log10_pressure(num_teff,atm_table_ng), &
           tbl%gmin_index(num_teff), tbl%gmax_index(num_teff))
      tbl%teff = 0.0d0
      tbl%logg = 0.0d0
      tbl%log10_pressure = 0.0d0
      tbl%gmin_index = 0
      tbl%gmax_index = 0
      tbl%num_teff = num_teff
end subroutine allocate_surface_p

! find_gravity_edges: G Somers 5/15; per log Teff row, the lowest and
! highest log g columns with a valid (positive) tabulated pressure.
! LCATCH (found_valid_pressure) is set to true once the first valid P
! is read, and it can set the minimum P thereafter.
subroutine find_gravity_edges(tbl)
      type(surface_p_table), intent(inout) :: tbl
         do teff_idx = 1,tbl%num_teff
            found_valid_pressure = .false.
            do logg_idx = atm_table_ng,1,-1
               if(tbl%log10_pressure(teff_idx,logg_idx) <= 0.0d0)then
!                 check if the first non-999 value has been reached.
!                 if so, set IMIN.
                  if(found_valid_pressure)then
                     tbl%gmin_index(teff_idx) = logg_idx + 1
                     exit
                  endif
               else
!                 if we have reached a non-negative pressure value,
!                 turn off the catch so IMIN can be set. also record
!                 the highest gravity with a pressure for later int.
                  if(.not.found_valid_pressure) tbl%gmax_index(teff_idx) = logg_idx
                  found_valid_pressure = .true.
               endif
            end do
            if (logg_idx < 1) tbl%gmin_index(teff_idx) = 1
!           if all of the P values at a given T are -999, set IMIN
!           to the number of gravity terms. in responce, the code
!           should break when trying to find surface P.
            if(.not.found_valid_pressure) tbl%gmin_index(teff_idx) = atm_table_ng
         end do
end subroutine find_gravity_edges

end subroutine read_surface_pressure_table

end module kurucz_atmosphere
