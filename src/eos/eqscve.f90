!----------------------------------------------------------------------
! eqscve
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqscve.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Interpolates in the surface (envelope) abundance table of the
! Saumon-Chabrier-Van Horn (SCV) equation of state -- i.e. the table
! built specifically for the envelope composition (X, Z). If the
! requested composition differs from the envelope composition, falls
! back to eqscvg (the general-abundance SCV table). Locates the SCV
! table cell in (log10 T, log10 Pgas), interpolates density and its
! logarithmic derivatives, the specific heat, adiabatic gradient, and
! ionization fractions, with an additional smoothing pass across
! adjacent table cells (the 4-pt Lagrange interpolation of the
! original version was replaced with a 4-pt natural spline; see the
! historical comments retained below).
subroutine eqscve(log10_temperature, temperature, pressure, &
     log10_density, density, hydrogen_fraction, metal_fraction, beta, &
     ion_fraction, dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
     adiabatic_gradient, valid_table_point)

      use numerics_lib
      implicit none

      integer, parameter :: nts = 63, nps = 76

! common/const1/: only ln10 is used here; the rest are placeholders
! preserving the shared storage layout (see getopac.f90).
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! common/const2/: no member is used in this file. Names match where
! gas_constant/radiation_constant_over_3 are actually read (eqstat2.f90).
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn
      common/const2/ gas_constant, radiation_constant_over_3, ca3l, csig, &
           csigl, cgl, cmkh, cmkhn

! common/scveos/: the SCV equation-of-state tables and the persistent
! (search-hunt) table indices, all used here. mhp 5/97 added this
! common block for SCV EOS tables.
      double precision :: table_log10_temperature(nts), &
           hydrogen_table(nts,nps,12), helium_table(nts,nps,12), &
           entropy_of_mixing_table(nts,nps), metal_table(nts,nps,13), &
           envelope_table(nts,nps,12)
      integer :: num_pressure_points(nts)
      logical :: use_scv_eos
      integer :: scv_temp_index, scv_pressure_index
      common/scveos/ table_log10_temperature, hydrogen_table, &
           helium_table, entropy_of_mixing_table, metal_table, &
           envelope_table, num_pressure_points, use_scv_eos, &
           scv_temp_index, scv_pressure_index

! common/comp/: only envelope_hydrogen_fraction/envelope_metal_fraction
! are used here; the rest are placeholders preserving the shared
! storage layout. Names are chosen to match their usage in
! eqstat2.f90, where the rest are read.
      double precision :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, zenvm, envelope_amu, &
           envelope_species_fractions(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, envelope_amu, envelope_species_fractions, xnew, znew, &
           stotal, senv

      double precision, intent(in) :: log10_temperature, temperature, &
           pressure
      double precision, intent(out) :: log10_density, density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction, beta
      double precision, intent(out) :: ion_fraction(3)
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      logical, intent(out) :: valid_table_point

! REPLACE 4-PT LAGRANGE INTERPOLATE WITH 4-PT NATURAL SPLINE
!      DIMENSION ATOMWT(4),QR(4),FT(4),FTD(4),
!     *     FP(4),FPD(4),TEMPT(5),TEMP(4,5),FXION(3)
      double precision :: press_interp_weights(4), &
           press_interp_weight_derivs(4), temp_interp_weights(4), &
           temp_interp_weight_derivs(4)
!       DIMENSION ATOMWT(4),QR(4),YTAB(4),Y2(4),  ! KC 2025-05-31
      double precision :: interp_nodes(4), ytab_work(4), &
           spline_second_deriv(4), interp_result(5), temp_grid(4,5), &
           interp_result_temp_shift(5), interp_result_press_shift(5), &
           interp_result_both_shift(5), interp_result_temp_smoothed(5), &
           interp_result_press_smoothed(5)
!       DATA ATOMWT/0.9921D0,0.24975D0,0.08322D0,0.4995D0/
!       DATA CMH,CMHE,CBOLTZ/1.67357D-24,6.646442D-24,1.380658D-16/
      double precision, parameter :: tol_pressure_smooth = 0.08d0, &
           tol_temp_smooth = 0.032d0

! --- locals ---
      double precision :: beta_complement, log10_gas_pressure
      double precision :: temp_dist_above, temp_dist_below, &
           temp_smooth_weight
      double precision :: press_dist_above, press_dist_below, &
           press_smooth_weight
      double precision :: spline_value
      double precision :: helium_fraction, radiation_pressure, gas_pressure
      double precision :: change_t1, change_t2
      double precision :: dlnrho_dlnp_gas, dlnrho_dlnt_gas, dlnp_dlnt_gas, &
           dlnp_dlnt, specific_heat_cp_gas, du_dt, &
           radiation_energy_density_term
      double precision :: xtf_h2, xtf_he, xtf_h1, xtf_hep, xtf_h_e, xtf_hp, &
           xtf_he_e, xtf_hepp, xhp, xhep, xhepp
      logical :: temp_needs_smoothing, press_needs_smoothing
      integer :: temp_smooth_direction, press_smooth_direction
      integer :: i, ii, j, jj, jjj, k

      save

      if (abs(hydrogen_fraction-envelope_hydrogen_fraction).gt.1.0d-5 &
           .or. abs(metal_fraction-envelope_metal_fraction).gt.1.0d-5) then
!          CALL EQSCVG(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,RMU,
!      *               AMU,EMU,ETA,QDT,QDP,QCP,DELA,LCALC)  ! KC 2025-05-31
         call eqscvg(log10_temperature, temperature, pressure, &
              log10_density, density, hydrogen_fraction, metal_fraction, &
              beta, ion_fraction, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, valid_table_point)
         return
      end if
! find gas pressure (which is the quantity which is tabulated).
      beta_complement = 1.0d0 - beta
      log10_gas_pressure = dlog10(beta*pressure)
!  check if the point is within the table
      if (log10_temperature.lt.table_log10_temperature(1) .or. &
           log10_temperature.gt.table_log10_temperature(nts) .or. &
           log10_gas_pressure.lt.4.0d0) then
         valid_table_point = .false.   ! Error exit - no valid table entry
         return
      end if
!  find nearest points in temperature.
      if (log10_temperature.lt.table_log10_temperature(scv_temp_index+1)) then
! search down to find nearest 4 table elements
         do i = scv_temp_index, 1, -1
            if (log10_temperature.gt.table_log10_temperature(i)) then
               ii = i - 1
               goto 10
            end if
         end do
         ii = 1
  10     continue
      else
! search up for nearest 4 table elements
         do i = scv_temp_index+2, nts
            if (log10_temperature.lt.table_log10_temperature(i)) then
               ii = i - 2
               goto 20
            end if
         end do
         ii = nts - 3
  20     continue
      end if
      scv_temp_index = max(1,ii)
      scv_temp_index = min(nts-3,scv_temp_index)
      temp_needs_smoothing = .false.
      if (scv_temp_index.eq.ii) then
         temp_dist_above = log10_temperature - &
              table_log10_temperature(scv_temp_index+1)
         temp_dist_below = table_log10_temperature(scv_temp_index+2) - &
              log10_temperature
         if (temp_dist_above.lt.tol_temp_smooth) then
            if (scv_temp_index.gt.1) then
               temp_needs_smoothing = .true.
               temp_smooth_direction = -1
               temp_smooth_weight = (temp_dist_above+tol_temp_smooth)/ &
                    (2.0d0*tol_temp_smooth)
            end if
         else if (temp_dist_below.lt.tol_temp_smooth) then
            if (scv_temp_index.lt.nts-3) then
               temp_needs_smoothing = .true.
               temp_smooth_direction = 1
               temp_smooth_weight = 0.5d0*temp_dist_below/tol_temp_smooth
            end if
         end if
      end if
!  find nearest points in pressure.
      jjj = min(num_pressure_points(scv_temp_index)-3, scv_pressure_index)
      if (log10_gas_pressure.lt.envelope_table(scv_temp_index,jjj+1,1)) then
! search down to find nearest 4 table elements
         do j = jjj, 1, -1
            if (log10_gas_pressure.gt.envelope_table(scv_temp_index,j,1)) then
               jj = j - 1
               goto 30
            end if
         end do
         jj = 1
  30     continue
         scv_pressure_index = max(1,jj)
         scv_pressure_index = min(num_pressure_points(scv_temp_index)-3, &
              scv_pressure_index)
      else
! search up for nearest 4 table elements.  note search is done at lowest
! temperature point (with the minimum range in p).
         do j = jjj+2, num_pressure_points(scv_temp_index)
            if (log10_gas_pressure.lt.envelope_table(scv_temp_index,j,1)) then
               jj = j - 2
               goto 40
            end if
         end do
! point is outside table; return.
         valid_table_point = .false.   ! Error exit - no valid table entry
         return
  40     continue
         scv_pressure_index = min(num_pressure_points(scv_temp_index)-3, jj)
      end if
      valid_table_point = .true.
      press_needs_smoothing = .false.
      if (scv_pressure_index.eq.jj) then
         press_dist_above = log10_gas_pressure - &
              envelope_table(scv_temp_index,scv_pressure_index+1,1)
         press_dist_below = envelope_table(scv_temp_index, &
              scv_pressure_index+2,1) - log10_gas_pressure
         if (press_dist_above.lt.tol_pressure_smooth) then
            if (scv_pressure_index.gt.1) then
               press_needs_smoothing = .true.
               press_smooth_direction = -1
               press_smooth_weight = (press_dist_above+tol_pressure_smooth)/ &
                    (2.0d0*tol_pressure_smooth)
            end if
         else if (press_dist_below.lt.tol_pressure_smooth) then
            if (scv_pressure_index.lt.num_pressure_points(scv_temp_index)-3) then
               press_needs_smoothing = .true.
               press_smooth_direction = 1
               press_smooth_weight = 0.5d0*press_dist_below/tol_pressure_smooth
            end if
         end if
      end if
!      DO K = 1,4
!         QR(K) = TLOGX(K+IDT-1)
!      END DO
!      CALL INTERP(QR,FT,FTD,TL)
!      DO K = 1,4
!         QR(K) = TABLENV(IDT,IDP+K-1,1)
!      END DO
!      CALL INTERP(QR,FP,FPD,PP)
      helium_fraction = 1.0d0 - hydrogen_fraction - metal_fraction
! include radiation pressure in the equation of state.
      radiation_pressure = beta_complement*pressure
      gas_pressure = beta*pressure
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         do k = 1,4
            interp_nodes(k) = envelope_table(scv_temp_index, &
                 scv_pressure_index+k-1,1)
         end do
         do j = 1,5
            do k = 1,4
               ytab_work(k) = envelope_table(ii,scv_pressure_index+k-1,j+1)
            end do
            call kspline(interp_nodes, ytab_work, spline_second_deriv)
            call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
                 log10_gas_pressure, spline_value)
            temp_grid(i,j) = spline_value
!            TEMP(I,J)=FP(1)*TABLENV(II,IDP,J+1) +
!     *      FP(2)*TABLENV(II,IDP+1,J+1) + FP(3)*TABLENV(II,IDP+2,J+1)
!     *      + FP(4)*TABLENV(II,IDP+3,J+1)
         end do
      end do
! interpolate in temperature
      do k = 1,4
         interp_nodes(k) = table_log10_temperature(k+scv_temp_index-1)
      end do
      do j = 1,5
         do k = 1,4
            ytab_work(k) = temp_grid(k,j)
         end do
         call kspline(interp_nodes, ytab_work, spline_second_deriv)
         call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
              log10_temperature, spline_value)
         interp_result(j) = spline_value
      end do
! perform interpolation at adjacent temperature table points
      if (temp_needs_smoothing) then
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1+temp_smooth_direction
         do k = 1,4
            interp_nodes(k) = envelope_table(scv_temp_index, &
                 scv_pressure_index+k-1,1)
         end do
         do j = 1,5
            do k = 1,4
               ytab_work(k) = envelope_table(ii,scv_pressure_index+k-1,j+1)
            end do
            call kspline(interp_nodes, ytab_work, spline_second_deriv)
            call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
                 log10_gas_pressure, spline_value)
            temp_grid(i,j) = spline_value
!            TEMP(I,J)=FP(1)*TABLENV(II,IDP,J+1) +
!     *      FP(2)*TABLENV(II,IDP+1,J+1) + FP(3)*TABLENV(II,IDP+2,J+1)
!     *      + FP(4)*TABLENV(II,IDP+3,J+1)
         end do
      end do
! interpolate in temperature
      do k = 1,4
         interp_nodes(k) = table_log10_temperature(k+scv_temp_index+ &
              temp_smooth_direction-1)
      end do
      do j = 1,5
         do k = 1,4
            ytab_work(k) = temp_grid(k,j)
         end do
         call kspline(interp_nodes, ytab_work, spline_second_deriv)
         call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
              log10_temperature, spline_value)
         interp_result_temp_shift(j) = spline_value
      end do
      if (temp_smooth_direction.eq.-1) then
         do j = 1,5
            interp_result_temp_smoothed(j) = interp_result_temp_shift(j) + &
                 temp_smooth_weight*(interp_result(j) - &
                 interp_result_temp_shift(j))
         end do
      else
         do j = 1,5
            interp_result_temp_smoothed(j) = interp_result(j) + &
                 temp_smooth_weight*(interp_result_temp_shift(j) - &
                 interp_result(j))
         end do
      end if
      end if
! perform interpolation at adjacent pressure table points
      if (press_needs_smoothing) then
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         do k = 1,4
            interp_nodes(k) = envelope_table(scv_temp_index, &
                 scv_pressure_index+k-1+press_smooth_direction,1)
         end do
         do j = 1,5
            do k = 1,4
               ytab_work(k) = envelope_table(ii,scv_pressure_index+k-1+ &
                    press_smooth_direction,j+1)
            end do
            call kspline(interp_nodes, ytab_work, spline_second_deriv)
            call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
                 log10_gas_pressure, spline_value)
            temp_grid(i,j) = spline_value
!            TEMP(I,J)=FP(1)*TABLENV(II,IDP,J+1) +
!     *      FP(2)*TABLENV(II,IDP+1,J+1) + FP(3)*TABLENV(II,IDP+2,J+1)
!     *      + FP(4)*TABLENV(II,IDP+3,J+1)
         end do
      end do
! interpolate in temperature
      do k = 1,4
         interp_nodes(k) = table_log10_temperature(k+scv_temp_index-1)
      end do
      do j = 1,5
         do k = 1,4
            ytab_work(k) = temp_grid(k,j)
         end do
         call kspline(interp_nodes, ytab_work, spline_second_deriv)
         call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
              log10_temperature, spline_value)
         interp_result_press_shift(j) = spline_value
      end do
      if (press_smooth_direction.eq.-1) then
         do j = 1,5
            interp_result_press_smoothed(j) = interp_result_press_shift(j) + &
                 press_smooth_weight*(interp_result(j) - &
                 interp_result_press_shift(j))
         end do
      else
         do j = 1,5
            interp_result_press_smoothed(j) = interp_result(j) + &
                 press_smooth_weight*(interp_result_press_shift(j) - &
                 interp_result(j))
         end do
      end if
      end if
! perform interpolation at adjacent t+p table points
      if (press_needs_smoothing .and. temp_needs_smoothing) then
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1 + temp_smooth_direction
         do k = 1,4
            interp_nodes(k) = envelope_table(scv_temp_index, &
                 scv_pressure_index+k-1+press_smooth_direction,1)
         end do
         do j = 1,5
            do k = 1,4
               ytab_work(k) = envelope_table(ii,scv_pressure_index+k-1+ &
                    press_smooth_direction,j+1)
            end do
            call kspline(interp_nodes, ytab_work, spline_second_deriv)
            call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
                 log10_gas_pressure, spline_value)
            temp_grid(i,j) = spline_value
!            TEMP(I,J)=FP(1)*TABLENV(II,IDP,J+1) +
!     *      FP(2)*TABLENV(II,IDP+1,J+1) + FP(3)*TABLENV(II,IDP+2,J+1)
!     *      + FP(4)*TABLENV(II,IDP+3,J+1)
         end do
      end do
! interpolate in temperature
      do k = 1,4
         interp_nodes(k) = table_log10_temperature(k+scv_temp_index-1 + &
              temp_smooth_direction)
      end do
      do j = 1,5
         do k = 1,4
            ytab_work(k) = temp_grid(k,j)
         end do
         call kspline(interp_nodes, ytab_work, spline_second_deriv)
         call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
              log10_temperature, spline_value)
         interp_result_both_shift(j) = spline_value
      end do
      end if
      if (temp_needs_smoothing) then
! add changes for both t and p interpolation
         if (press_needs_smoothing) then
!            WRITE(*,911)(TEMPT(J),J=1,5)
!  911        FORMAT(1X,'ORIG ',1P5E16.7)
!            WRITE(*,912)(TEMPT1(J)-TEMPT(J),J=1,5)
!  912        FORMAT(1X,'INT T',1P5E16.7)
!            WRITE(*,913)(TEMPT2(J)-TEMPT(J),J=1,5)
!  913        FORMAT(1X,'INT P',1P5E16.7)
!            WRITE(*,914)(TEMPT3(J)-TEMPT(J),J=1,5)
!  914        FORMAT(1X,'INT PT',1P5E16.7)
            do j = 1,5
               if (temp_smooth_direction.eq.-1) then
! interpolate in t at fixed p
                  change_t1 = interp_result_temp_shift(j) + &
                       temp_smooth_weight*(interp_result(j) - &
                       interp_result_temp_shift(j))
! interpolate in t at different p
                  change_t2 = interp_result_both_shift(j) + &
                       temp_smooth_weight*(interp_result_press_shift(j) - &
                       interp_result_both_shift(j))
               else
                  change_t1 = interp_result(j) + &
                       temp_smooth_weight*(interp_result_temp_shift(j) - &
                       interp_result(j))
                  change_t2 = interp_result_press_shift(j) + &
                       temp_smooth_weight*(interp_result_both_shift(j) - &
                       interp_result_press_shift(j))
               end if
               if (press_smooth_direction.eq.-1) then
                  interp_result(j) = change_t2+press_smooth_weight* &
                       (change_t1-change_t2)
! change in p at fixed t
!                  CHGP1 = TEMPT2+FSP*(TEMPT(J)-TEMPT2(J))
! change in p at different t
!                  CHGP2 = TEMPT3(J)+FSP*(TEMPT1(J)-TEMPT3(J))
               else
                  interp_result(j) = change_t1+press_smooth_weight* &
                       (change_t2-change_t1)
!                  CHGP1 = TEMPT(J)+FSP*(TEMPT2(J)-TEMPT(J))
!                  CHGP2 = TEMPT1(J)+FSP*(TEMPT3(J)-TEMPT1(J))
               end if
!               WRITE(*,915)ISMP,ISMT,TEMPT(J),CHGT1,CHGT2,FSP,FST
!  915           FORMAT(2I2,1P5E16.7)
            end do
          else
! add t interpolation changes only
!            WRITE(*,911)(TEMPT(J),J=1,5)
!            WRITE(*,912)(TEMPT1(J)-TEMPT(J),J=1,5)
             do j = 1,5
                interp_result(j) = interp_result_temp_smoothed(j)
             end do
          end if
       else if (press_needs_smoothing) then
! add p interpolation changes only
!            WRITE(*,911)(TEMPT(J),J=1,5)
!            WRITE(*,912)(TEMPT2(J)-TEMPT(J),J=1,5)
          do j = 1,5
             interp_result(j) = interp_result_press_smoothed(j)
          end do
       end if
!      DO J = 1,5
!         TEMPT(J)=FT(1)*TEMP(1,J) + FT(2)*TEMP(2,J) + FT(3)*TEMP(3,J)
!     *   + FT(4)*TEMP(4,J)
!      END DO
!  density
      log10_density = interp_result(1)
      density = exp(ln10*log10_density)
! d ln rho/ d ln p = qdp
! for gas pressure qdp(tot) = qdp(x)*x*rho/rho(x)+qdp(y)*y*rho/rho(y)
! for radiation pressure qdp = 0, so qdp(tot) = qdp(gas)*p/pgas
      dlnrho_dlnp_gas = interp_result(3)
      dlnrho_dlnp = dlnrho_dlnp_gas/beta
! d ln rho/ d ln t = qdt (note : d ln p/ d ln t = qpt)
! for gas pressure, correct as per qdp
! for radiation pressure, use qdt = -qdp*qpt.  correct qpt for
! radiation pressure and use the corrected qdp, qpt to get qdt.
      dlnrho_dlnt_gas = interp_result(2)
      dlnp_dlnt_gas = -dlnrho_dlnt_gas/dlnrho_dlnp_gas
      dlnp_dlnt = dlnp_dlnt_gas*beta + 4.0d0*beta_complement
      dlnrho_dlnt = -dlnp_dlnt*dlnrho_dlnp
!  cp = s*(d ln s/ d ln t)|p is tabulated. use
!  cp = du/dt + p*(d ln rho/d ln t)**2/rho/t/(d ln rho/ d ln p)
!  to find du/dt.  then include the effects of radiation pressure
!  on du/dt and the otehr terms and get a corrected cp.
! cp (gas pressure only).
      specific_heat_cp_gas = interp_result(4)
! now find du/dt from the original table.
      du_dt = interp_result(5)
!      QUT = QCP0 - PGAS*QDT0**2/QDP0/D/T
! correct du/dt for radiation
      radiation_energy_density_term = 3.0d0*radiation_pressure/density
      du_dt = du_dt + 4.0d0*radiation_energy_density_term/temperature
! correct cp for radiation pressure
      specific_heat_cp = du_dt - pressure*dlnrho_dlnt*dlnp_dlnt/density/ &
           temperature
! adiabatic temperature gradient
      adiabatic_gradient = -pressure*dlnrho_dlnt/density/temperature/ &
           specific_heat_cp

! Get fractions of total particles (including electrons), as follows:
!   XTF_H2  the fraction that is neutral hydrogen molecules, and
!   XTF_He  the fraction thst is neutral helium).
!   These are in column 2 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.

      do k = 1,4
         interp_nodes(k) = table_log10_temperature(k+scv_temp_index-1)
      end do
      call interp(interp_nodes, temp_interp_weights, &
           temp_interp_weight_derivs, log10_temperature)
      do k = 1,4
         interp_nodes(k) = envelope_table(scv_temp_index, &
              scv_pressure_index+k-1,1)
      end do
      call interp(interp_nodes, press_interp_weights, &
           press_interp_weight_derivs, log10_gas_pressure)

      do i = 1,4
         ii = scv_temp_index+i-1
         temp_grid(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,2) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,2) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,2) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,2)
         temp_grid(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,2) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,2) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,2) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,2)
      end do
! interpolate in temperature
      xtf_h2 = temp_interp_weights(1)*temp_grid(1,1) + &
           temp_interp_weights(2)*temp_grid(2,1) + &
           temp_interp_weights(3)*temp_grid(3,1) &
           + temp_interp_weights(4)*temp_grid(4,1)
      xtf_he = temp_interp_weights(1)*temp_grid(1,2) + &
           temp_interp_weights(2)*temp_grid(2,2) + &
           temp_interp_weights(3)*temp_grid(3,2) &
           + temp_interp_weights(4)*temp_grid(4,2)

! Get more fractions of total particles (including electrons), as follows:
!   XTF_H1  the fraction that is neutral hydrogen atoms, and
!   XTF_HeP  the fraction thst is singly ionized helium).
!   These are in column 3 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_grid(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,3) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,3) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,3) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,3)
         temp_grid(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,3) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,3) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,3) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,3)
      end do
! interpolate in temperature
      xtf_h1 = temp_interp_weights(1)*temp_grid(1,1) + &
           temp_interp_weights(2)*temp_grid(2,1) + &
           temp_interp_weights(3)*temp_grid(3,1) &
           + temp_interp_weights(4)*temp_grid(4,1)
      xtf_hep = temp_interp_weights(1)*temp_grid(1,2) + &
           temp_interp_weights(2)*temp_grid(2,2) + &
           temp_interp_weights(3)*temp_grid(3,2) &
           + temp_interp_weights(4)*temp_grid(4,2)

! At tis time we can calculate the ramaining particles fractions using
! conservation of particles and coservation of charge.  The variable
! names are:
!  XTF_H_e  the hycrogen related fraction that is electrons
!  XTF_HP   the hydrogen ralated fraction that is H+ ions
!  XTF_He_e the helium related fraction that is electrons
!  XTF_HePP the helium related fraction that is doubly ionized helium (He++)
! Particle and charge conservation yields:
      xtf_h_e = .5d0*(1d0 - xtf_h2 - xtf_h1)
      xtf_hp = xtf_h_e
      xtf_he_e = (1d0/3d0)*(2d0 - 2d0*xtf_he - xtf_hep)
      xtf_hepp = 1d0 - xtf_he - xtf_hep - xtf_he_e

! We are seeking the fraction XHP of hydrogen nuclei that are singly
! ionized, the fraction XHeP of helium huclei that are singly ionized
! and XHePP, the helium fraction that is doubly ionized.
      xhp = xtf_hp / (xtf_h2 + xtf_h1 + xtf_hp)
      xhep = xtf_hep / (xtf_he + xtf_hep + xtf_hepp)
      xhepp = xtf_hepp / (xtf_he + xtf_hep + xtf_hepp)

! In the same order, these are the desired elements of array ion_fraction(3)
      ion_fraction(1) = xhp
      ion_fraction(2) = xhep
      ion_fraction(3) = xhepp

      return  ! Normal exit. Valid table entry. valid_table_point is true.

end subroutine eqscve
