!----------------------------------------------------------------------
! eqscvg
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqscvg.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Interpolates in the general-abundance tables of the Saumon-
! Chabrier-Van Horn (SCV) equation of state. Given log10(T), the
! total pressure P, and the composition (X, Z), locates the SCV
! table cell in (log10 T, log10 Pgas) and interpolates density and
! its logarithmic derivatives, plus the specific heat, adiabatic
! gradient, and ionization fractions.
subroutine eqscvg(log10_temperature, temperature, pressure, &
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

! common/comp/: none of these members are used in this file; declared
! only to preserve the storage layout shared with every other file
! that references common/comp/. Names are chosen to match their usage
! in eqstat2.f90, where they are read.
      double precision :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, zenvm, envelope_amu, &
           envelope_species_fractions(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, envelope_amu, envelope_species_fractions, xnew, znew, &
           stotal, senv

      double precision, intent(in) :: log10_temperature, temperature
      double precision, intent(in) :: pressure
      double precision, intent(out) :: log10_density, density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction, beta
      double precision, intent(out) :: ion_fraction(3)
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      logical, intent(out) :: valid_table_point

      double precision :: interp_nodes(4), temp_interp_weights(4), &
           temp_interp_weight_derivs(4), press_interp_weights(4), &
           press_interp_weight_derivs(4), temp_work(4,3)

! --- locals ---
      double precision :: beta_complement, log10_gas_pressure, &
           helium_fraction, radiation_pressure, gas_pressure
      double precision :: log10_density_pure_h, density_pure_h, &
           log10_density_pure_he, density_pure_he, log10_density_pure_z, &
           density_pure_z
      double precision :: volume_frac_h, volume_frac_he, volume_frac_z
      double precision :: dlnrho_dlnp_pure_h, dlnrho_dlnp_pure_he, &
           dlnrho_dlnp_pure_z, dlnrho_dlnp_gas
      double precision :: dlnrho_dlnt_pure_h, dlnrho_dlnt_pure_he, &
           dlnrho_dlnt_pure_z, dlnrho_dlnt_gas, dlnp_dlnt_gas, dlnp_dlnt
      double precision :: dlns_dlnt_p, specific_heat_cp_gas, du_dt, &
           radiation_energy_density_term
      double precision :: log10_entropy_pure_h, entropy_pure_h, &
           log10_entropy_pure_he, entropy_pure_he, entropy_of_mixing, &
           dlnsmix_dlnt, total_entropy
      double precision :: dlns_dlnt_pure_h, dlns_dlnt_pure_he, &
           log10_du_dt_pure_z, du_dt_pure_z
      double precision :: xtf_h2, xtf_he, xtf_h1, xtf_hep, xtf_h_e, xtf_hp, &
           xtf_he_e, xtf_hepp, xhp, xhep, xhepp
      integer :: i, ii, j, jj, jjj, k

      save

      beta_complement = 1.0d0 - beta
      log10_gas_pressure = dlog10(beta*pressure)
! check if the point is within the table
      if (log10_temperature.lt.table_log10_temperature(1) .or. &
           log10_temperature.gt.table_log10_temperature(nts) .or. &
           log10_gas_pressure.lt.4.0d0) then
         valid_table_point = .false.
         return
      end if
! find nearest points in temperature.
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
         scv_temp_index = max(1,ii)
         scv_temp_index = min(nts-3,scv_temp_index)
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
         scv_temp_index = max(1,ii)
         scv_temp_index = min(nts-3,scv_temp_index)
      end if
! find nearest points in pressure.
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
!         WRITE(*,5)TL,PL
         valid_table_point = .false.
         return
  40     continue
         scv_pressure_index = min(num_pressure_points(scv_temp_index)-3, jj)
      end if
      valid_table_point = .true.
      do k = 1,4
         interp_nodes(k) = table_log10_temperature(k+scv_temp_index-1)
      end do
      call interp(interp_nodes, temp_interp_weights, &
           temp_interp_weight_derivs, log10_temperature)
      do k = 1,4
         interp_nodes(k) = envelope_table(scv_temp_index, &
              scv_pressure_index+k-1, 1)
      end do
      call interp(interp_nodes, press_interp_weights, &
           press_interp_weight_derivs, log10_gas_pressure)
      helium_fraction = 1.0d0 - hydrogen_fraction - metal_fraction
! include radiation pressure in the equation of state.
      radiation_pressure = beta_complement*pressure
      gas_pressure = beta*pressure
! density for x=1, y=1, z=1
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_work(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,4) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,4) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,4) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,4)
         temp_work(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,4) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,4) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,4) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,4)
         temp_work(i,3) = press_interp_weights(1)*metal_table(ii,scv_pressure_index,4) + &
         press_interp_weights(2)*metal_table(ii,scv_pressure_index+1,4) + &
         press_interp_weights(3)*metal_table(ii,scv_pressure_index+2,4) &
         + press_interp_weights(4)*metal_table(ii,scv_pressure_index+3,4)
      end do
! interpolate in temperature
      log10_density_pure_h = temp_interp_weights(1)*temp_work(1,1) + &
           temp_interp_weights(2)*temp_work(2,1) + &
           temp_interp_weights(3)*temp_work(3,1) &
           + temp_interp_weights(4)*temp_work(4,1)
      density_pure_h = exp(ln10*log10_density_pure_h)
      log10_density_pure_he = temp_interp_weights(1)*temp_work(1,2) + &
           temp_interp_weights(2)*temp_work(2,2) + &
           temp_interp_weights(3)*temp_work(3,2) &
           + temp_interp_weights(4)*temp_work(4,2)
      density_pure_he = exp(ln10*log10_density_pure_he)
      log10_density_pure_z = temp_interp_weights(1)*temp_work(1,3) + &
           temp_interp_weights(2)*temp_work(2,3) + &
           temp_interp_weights(3)*temp_work(3,3) &
           + temp_interp_weights(4)*temp_work(4,3)
      density_pure_z = exp(ln10*log10_density_pure_z)
! density
      density = 1.0d0/(hydrogen_fraction/density_pure_h + &
           helium_fraction/density_pure_he + metal_fraction/density_pure_z)
      log10_density = log10(density)
      volume_frac_h = hydrogen_fraction*density/density_pure_h
      volume_frac_he = helium_fraction*density/density_pure_he
      volume_frac_z = metal_fraction*density/density_pure_z
! d ln rho/ d ln p = qdp
! for gas pressure qdp(tot) = qdp(x)*x*rho/rho(x)+qdp(y)*y*rho/rho(y)
! for radiation pressure qdp = 0, so qdp(tot) = qdp(gas)*p/pgas
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_work(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,8) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,8) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,8) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,8)
         temp_work(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,8) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,8) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,8) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,8)
         temp_work(i,3) = press_interp_weights(1)*metal_table(ii,scv_pressure_index,13) + &
         press_interp_weights(2)*metal_table(ii,scv_pressure_index+1,13) + &
         press_interp_weights(3)*metal_table(ii,scv_pressure_index+2,13) &
         + press_interp_weights(4)*metal_table(ii,scv_pressure_index+3,13)
      end do
! interpolate in temperature
      dlnrho_dlnp_pure_h = temp_interp_weights(1)*temp_work(1,1) + &
           temp_interp_weights(2)*temp_work(2,1) + &
           temp_interp_weights(3)*temp_work(3,1) &
           + temp_interp_weights(4)*temp_work(4,1)
      dlnrho_dlnp_pure_he = temp_interp_weights(1)*temp_work(1,2) + &
           temp_interp_weights(2)*temp_work(2,2) + &
           temp_interp_weights(3)*temp_work(3,2) &
           + temp_interp_weights(4)*temp_work(4,2)
      dlnrho_dlnp_pure_z = temp_interp_weights(1)*temp_work(1,3) + &
           temp_interp_weights(2)*temp_work(2,3) + &
           temp_interp_weights(3)*temp_work(3,3) &
           + temp_interp_weights(4)*temp_work(4,3)
      dlnrho_dlnp_gas = dlnrho_dlnp_pure_h*volume_frac_h + &
           dlnrho_dlnp_pure_he*volume_frac_he + dlnrho_dlnp_pure_z*volume_frac_z
      dlnrho_dlnp = dlnrho_dlnp_gas/beta
! d ln rho/ d ln t = qdt (note : d ln p/ d ln t = qpt)
! for gas pressure, correct as per qdp
! for radiation pressure, use qdt = -qdp*qpt.  correct qpt for
! radiation pressure and use the corrected qdp, qpt to get qdt.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_work(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,7) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,7) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,7) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,7)
         temp_work(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,7) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,7) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,7) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,7)
         temp_work(i,3) = press_interp_weights(1)*metal_table(ii,scv_pressure_index,10) + &
         press_interp_weights(2)*metal_table(ii,scv_pressure_index+1,10) + &
         press_interp_weights(3)*metal_table(ii,scv_pressure_index+2,10) &
         + press_interp_weights(4)*metal_table(ii,scv_pressure_index+3,10)
      end do
! interpolate in temperature
      dlnrho_dlnt_pure_h = temp_interp_weights(1)*temp_work(1,1) + &
           temp_interp_weights(2)*temp_work(2,1) + &
           temp_interp_weights(3)*temp_work(3,1) &
           + temp_interp_weights(4)*temp_work(4,1)
      dlnrho_dlnt_pure_he = temp_interp_weights(1)*temp_work(1,2) + &
           temp_interp_weights(2)*temp_work(2,2) + &
           temp_interp_weights(3)*temp_work(3,2) &
           + temp_interp_weights(4)*temp_work(4,2)
      dlnrho_dlnt_pure_z = temp_interp_weights(1)*temp_work(1,3) + &
           temp_interp_weights(2)*temp_work(2,3) + &
           temp_interp_weights(3)*temp_work(3,3) &
           + temp_interp_weights(4)*temp_work(4,3)
      dlnrho_dlnt_gas = dlnrho_dlnt_pure_h*volume_frac_h + &
           dlnrho_dlnt_pure_he*volume_frac_he + dlnrho_dlnt_pure_z*volume_frac_z
      dlnp_dlnt_gas = -dlnrho_dlnt_gas/dlnrho_dlnp_gas
      dlnp_dlnt = dlnp_dlnt_gas*beta + 4.0d0*beta_complement
      dlnrho_dlnt = -dlnp_dlnt*dlnrho_dlnp
! cp = s*(d ln s/ d ln t)|p is tabulated. use
! cp = du/dt + p*(d ln rho/d ln t)**2/rho/t/(d ln rho/ d ln p)
! to find du/dt.  then include the effects of radiation pressure
! on du/dt and the other terms and get a corrected cp.
! cp (gas pressure only).
! note that unlike the envelope table, we need to interpolate in
! entropy s, d ln s/d lnt, and entropy of mixing, then calculate
! du/dt
! entropy and entropy of mixing
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_work(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,5) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,5) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,5) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,5)
         temp_work(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,5) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,5) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,5) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,5)
         temp_work(i,3) = press_interp_weights(1)*entropy_of_mixing_table(ii,scv_pressure_index) + &
         press_interp_weights(2)*entropy_of_mixing_table(ii,scv_pressure_index+1) + &
         press_interp_weights(3)*entropy_of_mixing_table(ii,scv_pressure_index+2) &
         + press_interp_weights(4)*entropy_of_mixing_table(ii,scv_pressure_index+3)
      end do
! interpolate in temperature
      log10_entropy_pure_h = temp_interp_weights(1)*temp_work(1,1) + &
           temp_interp_weights(2)*temp_work(2,1) + &
           temp_interp_weights(3)*temp_work(3,1) &
           + temp_interp_weights(4)*temp_work(4,1)
      entropy_pure_h = exp(ln10*log10_entropy_pure_h)
      log10_entropy_pure_he = temp_interp_weights(1)*temp_work(1,2) + &
           temp_interp_weights(2)*temp_work(2,2) + &
           temp_interp_weights(3)*temp_work(3,2) &
           + temp_interp_weights(4)*temp_work(4,2)
      entropy_pure_he = exp(ln10*log10_entropy_pure_he)
      entropy_of_mixing = temp_interp_weights(1)*temp_work(1,3) + &
           temp_interp_weights(2)*temp_work(2,3) + &
           temp_interp_weights(3)*temp_work(3,3) &
           + temp_interp_weights(4)*temp_work(4,3)
      dlnsmix_dlnt = temp_interp_weight_derivs(1)*log10(temp_work(1,3)) + &
           temp_interp_weight_derivs(2)*log10(temp_work(2,3)) &
           + temp_interp_weight_derivs(3)*log10(temp_work(3,3)) &
           + temp_interp_weight_derivs(4)*log10(temp_work(4,3))
! total_entropy is computed but not used further below, matching the
! original (dead) computation of S0 in eqscvg.f.
      total_entropy = hydrogen_fraction*entropy_pure_h + &
           helium_fraction*entropy_pure_he + entropy_of_mixing
! d ln s/ d ln t (x and y) and du/dt (z)
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_work(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,9) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,9) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,9) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,9)
         temp_work(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,9) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,9) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,9) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,9)
         temp_work(i,3) = press_interp_weights(1)*metal_table(ii,scv_pressure_index,7) + &
         press_interp_weights(2)*metal_table(ii,scv_pressure_index+1,7) + &
         press_interp_weights(3)*metal_table(ii,scv_pressure_index+2,7) &
         + press_interp_weights(4)*metal_table(ii,scv_pressure_index+3,7)
      end do
! interpolate in temperature
      dlns_dlnt_pure_h = temp_interp_weights(1)*temp_work(1,1) + &
           temp_interp_weights(2)*temp_work(2,1) + &
           temp_interp_weights(3)*temp_work(3,1) &
           + temp_interp_weights(4)*temp_work(4,1)
      dlns_dlnt_pure_he = temp_interp_weights(1)*temp_work(1,2) + &
           temp_interp_weights(2)*temp_work(2,2) + &
           temp_interp_weights(3)*temp_work(3,2) &
           + temp_interp_weights(4)*temp_work(4,2)
      log10_du_dt_pure_z = temp_interp_weights(1)*temp_work(1,3) + &
           temp_interp_weights(2)*temp_work(2,3) + &
           temp_interp_weights(3)*temp_work(3,3) &
           + temp_interp_weights(4)*temp_work(4,3)
      du_dt_pure_z = exp(ln10*log10_du_dt_pure_z)
      dlns_dlnt_p = hydrogen_fraction*entropy_pure_h*dlns_dlnt_pure_h + &
           helium_fraction*entropy_pure_he*dlns_dlnt_pure_he + &
           entropy_of_mixing*dlnsmix_dlnt
      specific_heat_cp_gas = dlns_dlnt_p + metal_fraction*(du_dt_pure_z + &
           gas_pressure*dlnrho_dlnt_pure_z**2/dlnrho_dlnp_pure_z/ &
           density_pure_z/temperature)
! now find du/dt from the original table.
      du_dt = specific_heat_cp_gas - gas_pressure*dlnrho_dlnt_gas**2/ &
           dlnrho_dlnp_gas/density/temperature
! correct du/dt for radiation
      radiation_energy_density_term = 3.0d0*radiation_pressure/density
      du_dt = du_dt + 4.0d0*radiation_energy_density_term/temperature
! correct cp for radiation pressure
      specific_heat_cp = du_dt + pressure*dlnrho_dlnt**2/dlnrho_dlnp/ &
           density/temperature
! adiabatic temperature gradient
      adiabatic_gradient = -pressure*dlnrho_dlnt/density/temperature/ &
           specific_heat_cp

! Get fractions of total particles (including electrons), as follows:
!   XTF_H2  the fraction that is neutral hydrogen molecules, and
!   XTF_He  the fraction thst is neutral helium).
!   These are in column 2 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_work(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,2) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,2) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,2) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,2)
         temp_work(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,2) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,2) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,2) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,2)
      end do
! interpolate in temperature
      xtf_h2 = temp_interp_weights(1)*temp_work(1,1) + &
           temp_interp_weights(2)*temp_work(2,1) + &
           temp_interp_weights(3)*temp_work(3,1) &
           + temp_interp_weights(4)*temp_work(4,1)
      xtf_he = temp_interp_weights(1)*temp_work(1,2) + &
           temp_interp_weights(2)*temp_work(2,2) + &
           temp_interp_weights(3)*temp_work(3,2) &
           + temp_interp_weights(4)*temp_work(4,2)

! Get more fractions of total particles (including electrons), as follows:
!   XTF_H1  the fraction that is neutral hydrogen atoms, and
!   XTF_HeP  the fraction thst is singly ionized helium).
!   These are in column 3 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = scv_temp_index+i-1
         temp_work(i,1) = press_interp_weights(1)*hydrogen_table(ii,scv_pressure_index,3) + &
         press_interp_weights(2)*hydrogen_table(ii,scv_pressure_index+1,3) + &
         press_interp_weights(3)*hydrogen_table(ii,scv_pressure_index+2,3) &
         + press_interp_weights(4)*hydrogen_table(ii,scv_pressure_index+3,3)
         temp_work(i,2) = press_interp_weights(1)*helium_table(ii,scv_pressure_index,3) + &
         press_interp_weights(2)*helium_table(ii,scv_pressure_index+1,3) + &
         press_interp_weights(3)*helium_table(ii,scv_pressure_index+2,3) &
         + press_interp_weights(4)*helium_table(ii,scv_pressure_index+3,3)
      end do
! interpolate in temperature
      xtf_h1 = temp_interp_weights(1)*temp_work(1,1) + &
           temp_interp_weights(2)*temp_work(2,1) + &
           temp_interp_weights(3)*temp_work(3,1) &
           + temp_interp_weights(4)*temp_work(4,1)
      xtf_hep = temp_interp_weights(1)*temp_work(1,2) + &
           temp_interp_weights(2)*temp_work(2,2) + &
           temp_interp_weights(3)*temp_work(3,2) &
           + temp_interp_weights(4)*temp_work(4,2)

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

      return    ! Primay exit. Success. valid_table_point is true.

end subroutine eqscvg
