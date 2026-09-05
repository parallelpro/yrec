!----------------------------------------------------------------------
! surfp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original surfp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Interpolates a tabulated surface-pressure table (Kurucz ATMPL or
! Kurucz/Castelli ATMPLC) to find log(P) at the given log(Teff),
! log(g): a 4-point Lagrangian interpolation in temperature, and
! cubic-spline (or 3-point/linear fallback near the table's
! high-gravity edge) interpolation in gravity.
!
! 2026 wave 3 (R5): surface_p_interp is the interpolation on a
! surface_p_table; surfp (Kurucz, atm_table%kurucz_surface_p) and
! kcsurfp (Kurucz/Castelli, atm_table%castelli_surface_p; formerly
! its own tables/kcsurfp.f90) are the two entries core/envint_kernel
! calls. `diff -w surfp.f90 kcsurfp.f90` differed only in the table
! names, the row count (atm_table_nt vs atm_table_ntc, now
! tbl%num_teff) and line continuation.
subroutine surface_p_interp(tbl, log10_teff, log10_gravity, print_flag, ierr)

      use atm_table_lib
      use luout_lib
      use numerics_lib
      implicit none
      integer :: jerr_gate
! tbl%num_teff AND atm_table_ng FOR TABULATED SURFACE PRESSURES.

      type(surface_p_table), intent(inout) :: tbl
      double precision, intent(in) :: log10_teff, log10_gravity
      logical, intent(in) :: print_flag
! --- locals ---
      double precision :: teff_nodes(4), gravity_nodes(4), &
           pressure_at_nodes(4), unused_deriv(3), gravity_weights(3), &
           gravity_nodes_3pt(3), pressure_table_vals(4), &
           teff_spline_deriv(4), gravity_spline_deriv(4)
      double precision :: fx, interpolated_value
      integer :: row, row_base, node, k, kk, kkk

! INTERPOLATES IN TEMPERATURE WITH A 4-POINT CUBIC SPLINE (kspline/
! ksplint), AND IN GRAVITY THE SAME WAY IF 4 OR MORE POINTS ARE
! AVAILABLE (3-POINT LAGRANGIAN OR LINEAR AT THE TABLE EDGE).  IT WILL QUIT IF THE DESIRED DATA POINT
! HAS TEFF OR LOG G MORE THAN ONE TABLE POINT FROM THE DATA.
!
! CHECK TO ENSURE THAT DATA IS WITHIN TABLE.
      integer, intent(out) :: ierr

! Lower table edge, and the log Teff rows above which only 3 or 2 log g
! columns are tabulated (so 3-point Lagrangian or linear interpolation).
      double precision, parameter :: table_logteff_min = 3.5d0, table_logg_min = -0.5d0
      double precision, parameter :: logteff_three_point = 4.5d0, logteff_two_point = 4.55d0

      ierr = 0

      if (log10_teff.lt.table_logteff_min .or. log10_gravity.lt.table_logg_min) then
         write(terminal_unit,911) log10_teff, log10_gravity
         write(run_log_unit,911) log10_teff, log10_gravity
  911    format(1X,'DESIRED ATMOSPHERE OUTSIDE TABLE RANGE'/ &
              ' LOG TEFF',F10.6,' LOG G',F10.6/' RUN STOPPED')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the atm_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      endif
! TEMPERATURE INTERPOLATION FACTORS.
      do row = 1,tbl%num_teff
         if (log10_teff.le.tbl%teff(row)) exit
      end do
      if (row > tbl%num_teff) then
      row = tbl%num_teff
      end if
      row_base = max(1,row-2)
      row_base = min(tbl%num_teff-3,row_base)
      do k = 1,4
         teff_nodes(k) = tbl%teff(row_base+k-1)
      end do
! GRAVITY INTERPOLATION FACTORS.
      do row = row_base,row_base+3
         node = row-row_base+1
! CHECK IF 4 LOG VALUES AVAILABLE - OTHERWISE, USE 3 POINT LAGRANGIAN
! OR LINEAR INTERPOLATION.
         if (tbl%teff(row).gt.logteff_three_point) then
            if (tbl%teff(row).gt.logteff_two_point) then
! LINEAR INTERPOLATION
               fx = (log10_gravity-tbl%logg(atm_table_ng-1))/ &
                    (tbl%logg(atm_table_ng)-tbl%logg(atm_table_ng-1))
               pressure_at_nodes(node) = &
                    tbl%log10_pressure(row,atm_table_ng-1) + &
                    fx*(tbl%log10_pressure(row,atm_table_ng) - &
                    tbl%log10_pressure(row,atm_table_ng-1))
            else
! 3-POINT LAGRANGIAN INTERPOLATION.
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
         if (log10_gravity.ge.tbl%logg(atm_table_ng-1)) then
! DESIRED LOG G ABOVE SECOND TO TOP TABLE LOG G - USE TOP 4 LOG G VALUES.
            do kk = 1,4
               gravity_nodes(kk) = tbl%logg(atm_table_ng-4+kk)
               pressure_table_vals(kk) = &
                    tbl%log10_pressure(row,atm_table_ng-4+kk)
            end do
            call kspline(gravity_nodes, pressure_table_vals, &
                 gravity_spline_deriv)
            call ksplint(gravity_nodes, pressure_table_vals, &
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
! GENERAL CASE - FIND 4 NEAREST POINTS IN GRAVITY THAT ARE IN THE TABLE.
! G Somers, I changed NG to IMINMAX in the next line. This prevents the
! code from using -999 to interpolate in some instances.
         do k = tbl%gmax_index(row)-3,tbl%gmin_index(row),-1
            if (log10_gravity.lt.tbl%logg(k+2) .and. &
                 log10_gravity.ge.tbl%logg(k+1)) then
               kk = max(tbl%gmin_index(row),k)
               kk = min(atm_table_ng-3,kk)
               do kkk = 1,4
                  gravity_nodes(kkk) = tbl%logg(kk+kkk-1)
                  pressure_table_vals(kkk) = &
                       tbl%log10_pressure(row,kk+kkk-1)
               end do
               call kspline(gravity_nodes, pressure_table_vals, &
                    gravity_spline_deriv)
               call ksplint(gravity_nodes, pressure_table_vals, &
                    gravity_spline_deriv, log10_gravity, interpolated_value, jerr_gate)
               if (jerr_gate /= 0) then
                  ierr = jerr_gate
                  return
               end if
               pressure_at_nodes(node) = interpolated_value
               exit
            endif
         end do
         if (k .lt. tbl%gmin_index(row)) then
! DESIRED LOG G BELOW 2ND TABLE ENTRY -USE FIRST 4 POINTS.
         do k = 1,4
            gravity_nodes(k) = tbl%logg(k+tbl%gmin_index(row)-1)
            pressure_table_vals(k) = &
                 tbl%log10_pressure(row,k+tbl%gmin_index(row)-1)
         end do
         call kspline(gravity_nodes, pressure_table_vals, gravity_spline_deriv)
         call ksplint(gravity_nodes, pressure_table_vals, &
              gravity_spline_deriv, log10_gravity, interpolated_value, jerr_gate)
         if (jerr_gate /= 0) then
            ierr = jerr_gate
            return
         end if
         pressure_at_nodes(node) = interpolated_value
         end if
      end do
! INTERPOLATE IN TEMPERATURE TO FIND CORRECT LOG P.
      call kspline(teff_nodes, pressure_at_nodes, teff_spline_deriv)
      call ksplint(teff_nodes, pressure_at_nodes, teff_spline_deriv, &
           log10_teff, interpolated_value, jerr_gate)
      if (jerr_gate /= 0) then
         ierr = jerr_gate
         return
      end if
      atm_table%atm_log10_pressure = interpolated_value
      atm_table%atm_log10_temperature = log10_teff
! WRITE OUT INFORMATION TO THE MODEL FILE.
      if (print_flag) then
        write(run_log_unit,70)
70      format('********PRESSURE AT T=TEFF INTERPOLATED FROM TABULATED' &
              ,  ' VALUES********')
        write(run_log_unit,71) log10_teff, atm_table%atm_log10_pressure
71      format(' ',20X,'LOG (Teff) =',F10.5,' LOG P =',F10.5)
      endif
      return
end subroutine surface_p_interp

! surfp: log P from the Kurucz (1993) table, atm_choice 3.
subroutine surfp(log10_teff, log10_gravity, print_flag, ierr)
      use atm_table_lib
      implicit none
      double precision, intent(in) :: log10_teff, log10_gravity
      logical, intent(in) :: print_flag
      integer, intent(out) :: ierr
      call surface_p_interp(atm_table%kurucz_surface_p, log10_teff, &
           log10_gravity, print_flag, ierr)
end subroutine surfp

! kcsurfp: JNT 06/2014, same interpolation on the newer Kurucz/
! Castelli table, atm_choice 5.
subroutine kcsurfp(log10_teff, log10_gravity, print_flag, ierr)
      use atm_table_lib
      implicit none
      double precision, intent(in) :: log10_teff, log10_gravity
      logical, intent(in) :: print_flag
      integer, intent(out) :: ierr
      call surface_p_interp(atm_table%castelli_surface_p, log10_teff, &
           log10_gravity, print_flag, ierr)
end subroutine kcsurfp
