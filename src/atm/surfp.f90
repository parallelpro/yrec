!----------------------------------------------------------------------
! surfp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original surfp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Interpolates the (Kurucz) surface-pressure table ATMPL to find
! log(P) at the given log(Teff), log(g): a 4-point Lagrangian
! interpolation in temperature, and cubic-spline (or 3-point/linear
! fallback near the table's high-gravity edge) interpolation in
! gravity.
subroutine surfp(log10_teff, log10_gravity, print_flag)

      use luout_lib
      use numerics_lib
      implicit none
! PARAMETERS NT AND NG FOR TABULATED SURFACE PRESSURES.
      integer, parameter :: nt = 57, ng = 11
      integer, parameter :: ntc = 76, ngc = 11

      double precision, intent(in) :: log10_teff, log10_gravity
      logical, intent(in) :: print_flag

! common/atmprt/: only atm_log10_pressure/atm_log10_temperature (AP,
! AT) are set here. Naming matches alsurfp.f90.
      double precision :: atm_tau, atm_log10_pressure, &
           atm_log10_temperature, atm_log10_density, atm_opacity, &
           atm_ion_fraction(3)
      common/atmprt/atm_tau, atm_log10_pressure, atm_log10_temperature, &
           atm_log10_density, atm_opacity, atm_ion_fraction
! common/atmos2/: the Kurucz surface-pressure table, all used here.
! Naming is local to this batch.
      double precision :: kurucz_log10_pressure_table(nt,ng), &
           kurucz_teff_table(nt), kurucz_logg_table(ng), kurucz_table_z
      integer :: atm_table_file_unit
      common/atmos2/kurucz_log10_pressure_table, kurucz_teff_table, &
           kurucz_logg_table, kurucz_table_z, atm_table_file_unit
! common/fac/: kurucz_gmin_index/kurucz_gmax_index/teff_interp_start_index/
! gravity_interp_indices are used here; castelli_gmin_index/
! castelli_gmax_index are unused placeholders (KCSURFP's table sizes).
! Naming is local to this batch.
      integer :: kurucz_gmin_index(nt), kurucz_gmax_index(nt), &
           teff_interp_start_index, gravity_interp_indices(4), &
           castelli_gmin_index(ntc), castelli_gmax_index(ntc)
      common/fac/kurucz_gmin_index, kurucz_gmax_index, &
           teff_interp_start_index, gravity_interp_indices, &
           castelli_gmin_index, castelli_gmax_index

      save

! --- locals ---
      double precision :: teff_nodes(4), gravity_nodes(4), &
           pressure_at_nodes(4), unused_deriv(3), gravity_weights(3), &
           gravity_nodes_3pt(3), pressure_table_vals(4), &
           teff_spline_deriv(4), gravity_spline_deriv(4)
      double precision :: fx, interpolated_value
      integer :: row, row_base, node, k, kk, kkk

! SURFPL INTERPOLATES IN TEMPERATURE USING A 4-POINT LAGRANGIAN
! INTERPOLATOR, AND INTERPOLATES IN GRAVITY THE SAME WAY IF 4 OR
! MORE POINTS ARE AVAILABLE.  IT WILL QUIT IF THE DESIRED DATA POINT
! HAS TEFF OR LOG G MORE THAN ONE TABLE POINT FROM THE DATA.
!
! CHECK TO ENSURE THAT DATA IS WITHIN TABLE.
      if (log10_teff.lt.3.5d0 .or. log10_gravity.lt.-0.5d0) then
         write(iowr,911) log10_teff, log10_gravity
         write(short_file_unit,911) log10_teff, log10_gravity
  911    format(1X,'DESIRED ATMOSPHERE OUTSIDE TABLE RANGE'/ &
              ' LOG TEFF',F10.6,' LOG G',F10.6/' RUN STOPPED')
         stop
      endif
! TEMPERATURE INTERPOLATION FACTORS.
      do row = 1,nt
         if (log10_teff.le.kurucz_teff_table(row)) goto 10
      end do
      row = nt
   10 continue
      row_base = max(1,row-2)
      row_base = min(nt-3,row_base)
      do k = 1,4
         teff_nodes(k) = kurucz_teff_table(row_base+k-1)
      end do
      teff_interp_start_index = row_base
! GRAVITY INTERPOLATION FACTORS.
      do 20 row = row_base,row_base+3
         node = row-row_base+1
! CHECK IF 4 LOG VALUES AVAILABLE - OTHERWISE, USE 3 POINT LAGRANGIAN
! OR LINEAR INTERPOLATION.
         if (kurucz_teff_table(row).gt.4.5d0) then
            if (kurucz_teff_table(row).gt.4.55d0) then
! LINEAR INTERPOLATION
               fx = (log10_gravity-kurucz_logg_table(ng-1))/ &
                    (kurucz_logg_table(ng)-kurucz_logg_table(ng-1))
               pressure_at_nodes(node) = &
                    kurucz_log10_pressure_table(row,ng-1) + &
                    fx*(kurucz_log10_pressure_table(row,ng) - &
                    kurucz_log10_pressure_table(row,ng-1))
               gravity_interp_indices(node) = ng-1
            else
! 3-POINT LAGRANGIAN INTERPOLATION.
               do k = 1,3
                  gravity_nodes_3pt(k) = kurucz_logg_table(ng-3+k)
               end do
               call inter3(gravity_nodes_3pt, gravity_weights, &
                    unused_deriv, log10_gravity)
               pressure_at_nodes(node) = &
                    kurucz_log10_pressure_table(row,ng-2)*gravity_weights(1) + &
                    kurucz_log10_pressure_table(row,ng-1)*gravity_weights(2) + &
                    kurucz_log10_pressure_table(row,ng)*gravity_weights(3)
               gravity_interp_indices(node) = ng-2
            endif
            goto 20
         endif
         if (log10_gravity.ge.kurucz_logg_table(ng-1)) then
! DESIRED LOG G ABOVE SECOND TO TOP TABLE LOG G - USE TOP 4 LOG G VALUES.
            do kk = 1,4
               gravity_nodes(kk) = kurucz_logg_table(ng-4+kk)
               pressure_table_vals(kk) = &
                    kurucz_log10_pressure_table(row,ng-4+kk)
            end do
            call kspline(gravity_nodes, pressure_table_vals, &
                 gravity_spline_deriv)
            call ksplint(gravity_nodes, pressure_table_vals, &
                 gravity_spline_deriv, log10_gravity, interpolated_value)
            pressure_at_nodes(node) = interpolated_value
            gravity_interp_indices(node) = ng-3
            goto 20
         endif
! GENERAL CASE - FIND 4 NEAREST POINTS IN GRAVITY THAT ARE IN THE TABLE.
! G Somers, I changed NG to IMINMAX in the next line. This prevents the
! code from using -999 to interpolate in some instances.
         do k = kurucz_gmax_index(row)-3,kurucz_gmin_index(row),-1
            if (log10_gravity.lt.kurucz_logg_table(k+2) .and. &
                 log10_gravity.ge.kurucz_logg_table(k+1)) then
               kk = max(kurucz_gmin_index(row),k)
               kk = min(ng-3,kk)
               do kkk = 1,4
                  gravity_nodes(kkk) = kurucz_logg_table(kk+kkk-1)
                  pressure_table_vals(kkk) = &
                       kurucz_log10_pressure_table(row,kk+kkk-1)
               end do
               call kspline(gravity_nodes, pressure_table_vals, &
                    gravity_spline_deriv)
               call ksplint(gravity_nodes, pressure_table_vals, &
                    gravity_spline_deriv, log10_gravity, interpolated_value)
               pressure_at_nodes(node) = interpolated_value
               gravity_interp_indices(node) = kk
               goto 20
            endif
         end do
! DESIRED LOG G BELOW 2ND TABLE ENTRY -USE FIRST 4 POINTS.
         do k = 1,4
            gravity_nodes(k) = kurucz_logg_table(k+kurucz_gmin_index(row)-1)
            pressure_table_vals(k) = &
                 kurucz_log10_pressure_table(row,k+kurucz_gmin_index(row)-1)
         end do
         call kspline(gravity_nodes, pressure_table_vals, gravity_spline_deriv)
         call ksplint(gravity_nodes, pressure_table_vals, &
              gravity_spline_deriv, log10_gravity, interpolated_value)
         pressure_at_nodes(node) = interpolated_value
         gravity_interp_indices(node) = kurucz_gmin_index(row)
   20 continue
! INTERPOLATE IN TEMPERATURE TO FIND CORRECT LOG P.
      call kspline(teff_nodes, pressure_at_nodes, teff_spline_deriv)
      call ksplint(teff_nodes, pressure_at_nodes, teff_spline_deriv, &
           log10_teff, interpolated_value)
      atm_log10_pressure = interpolated_value
      atm_log10_temperature = log10_teff
! WRITE OUT INFORMATION TO THE MODEL FILE.
      if (print_flag) then
        write(short_file_unit,70)
        write(istor,70)
70      format('********PRESSURE AT T=TEFF INTERPOLATED FROM TABULATED' &
              ,  ' VALUES********')
        write(short_file_unit,71) log10_teff, atm_log10_pressure
        write(istor,71) log10_teff, atm_log10_pressure
71      format(' ',20X,'LOG (Teff) =',F10.5,' LOG P =',F10.5)
      endif
      return
end subroutine surfp
