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

      use phys_const_lib
      use star_info_lib
      use numerics_lib
      use scv_eos_lib
      implicit none

      integer, parameter :: nts = 63, nps = 76





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
      beta_complement = 1.0d0 - beta
      log10_gas_pressure = log10(beta*pressure)
! check if the point is within the table
      if (log10_temperature.lt.tlogx(1) .or. &
           log10_temperature.gt.tlogx(nts) .or. &
           log10_gas_pressure.lt.4.0d0) then
         valid_table_point = .false.
         return
      end if
! find nearest points in temperature.
      if (log10_temperature.lt.tlogx(idtt+1)) then
! search down to find nearest 4 table elements
         do i = idtt, 1, -1
            if (log10_temperature.gt.tlogx(i)) then
               ii = i - 1
               exit
            end if
         end do
         if (i < (1)) then
         ii = 1
         end if
         idtt = max(1,ii)
         idtt = min(nts-3,idtt)
      else
! search up for nearest 4 table elements
         do i = idtt+2, nts
            if (log10_temperature.lt.tlogx(i)) then
               ii = i - 2
               exit
            end if
         end do
         if (i > nts) then
         ii = nts - 3
         end if
         idtt = max(1,ii)
         idtt = min(nts-3,idtt)
      end if
! find nearest points in pressure.
      jjj = min(nptsx(idtt)-3, idp)
      if (log10_gas_pressure.lt.tablenv(idtt,jjj+1,1)) then
! search down to find nearest 4 table elements
         do j = jjj, 1, -1
            if (log10_gas_pressure.gt.tablenv(idtt,j,1)) then
               jj = j - 1
               exit
            end if
         end do
         if (j < (1)) then
         jj = 1
         end if
         idp = max(1,jj)
         idp = min(nptsx(idtt)-3, &
              idp)
      else
! search up for nearest 4 table elements.  note search is done at lowest
! temperature point (with the minimum range in p).
         do j = jjj+2, nptsx(idtt)
            if (log10_gas_pressure.lt.tablenv(idtt,j,1)) then
               jj = j - 2
               exit
            end if
         end do
         if (j > (nptsx(idtt))) then
! point is outside table; return.
!         WRITE(*,5)TL,PL
         valid_table_point = .false.
         return
         end if
         idp = min(nptsx(idtt)-3, jj)
      end if
      valid_table_point = .true.
      do k = 1,4
         interp_nodes(k) = tlogx(k+idtt-1)
      end do
      call interp(interp_nodes, temp_interp_weights, &
           temp_interp_weight_derivs, log10_temperature)
      do k = 1,4
         interp_nodes(k) = tablenv(idtt, &
              idp+k-1, 1)
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
         ii = idtt+i-1
         temp_work(i,1) = press_interp_weights(1)*tablex(ii,idp,4) + &
         press_interp_weights(2)*tablex(ii,idp+1,4) + &
         press_interp_weights(3)*tablex(ii,idp+2,4) &
         + press_interp_weights(4)*tablex(ii,idp+3,4)
         temp_work(i,2) = press_interp_weights(1)*tabley(ii,idp,4) + &
         press_interp_weights(2)*tabley(ii,idp+1,4) + &
         press_interp_weights(3)*tabley(ii,idp+2,4) &
         + press_interp_weights(4)*tabley(ii,idp+3,4)
         temp_work(i,3) = press_interp_weights(1)*tablez(ii,idp,4) + &
         press_interp_weights(2)*tablez(ii,idp+1,4) + &
         press_interp_weights(3)*tablez(ii,idp+2,4) &
         + press_interp_weights(4)*tablez(ii,idp+3,4)
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
         ii = idtt+i-1
         temp_work(i,1) = press_interp_weights(1)*tablex(ii,idp,8) + &
         press_interp_weights(2)*tablex(ii,idp+1,8) + &
         press_interp_weights(3)*tablex(ii,idp+2,8) &
         + press_interp_weights(4)*tablex(ii,idp+3,8)
         temp_work(i,2) = press_interp_weights(1)*tabley(ii,idp,8) + &
         press_interp_weights(2)*tabley(ii,idp+1,8) + &
         press_interp_weights(3)*tabley(ii,idp+2,8) &
         + press_interp_weights(4)*tabley(ii,idp+3,8)
         temp_work(i,3) = press_interp_weights(1)*tablez(ii,idp,13) + &
         press_interp_weights(2)*tablez(ii,idp+1,13) + &
         press_interp_weights(3)*tablez(ii,idp+2,13) &
         + press_interp_weights(4)*tablez(ii,idp+3,13)
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
         ii = idtt+i-1
         temp_work(i,1) = press_interp_weights(1)*tablex(ii,idp,7) + &
         press_interp_weights(2)*tablex(ii,idp+1,7) + &
         press_interp_weights(3)*tablex(ii,idp+2,7) &
         + press_interp_weights(4)*tablex(ii,idp+3,7)
         temp_work(i,2) = press_interp_weights(1)*tabley(ii,idp,7) + &
         press_interp_weights(2)*tabley(ii,idp+1,7) + &
         press_interp_weights(3)*tabley(ii,idp+2,7) &
         + press_interp_weights(4)*tabley(ii,idp+3,7)
         temp_work(i,3) = press_interp_weights(1)*tablez(ii,idp,10) + &
         press_interp_weights(2)*tablez(ii,idp+1,10) + &
         press_interp_weights(3)*tablez(ii,idp+2,10) &
         + press_interp_weights(4)*tablez(ii,idp+3,10)
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
         ii = idtt+i-1
         temp_work(i,1) = press_interp_weights(1)*tablex(ii,idp,5) + &
         press_interp_weights(2)*tablex(ii,idp+1,5) + &
         press_interp_weights(3)*tablex(ii,idp+2,5) &
         + press_interp_weights(4)*tablex(ii,idp+3,5)
         temp_work(i,2) = press_interp_weights(1)*tabley(ii,idp,5) + &
         press_interp_weights(2)*tabley(ii,idp+1,5) + &
         press_interp_weights(3)*tabley(ii,idp+2,5) &
         + press_interp_weights(4)*tabley(ii,idp+3,5)
         temp_work(i,3) = press_interp_weights(1)*smix(ii,idp) + &
         press_interp_weights(2)*smix(ii,idp+1) + &
         press_interp_weights(3)*smix(ii,idp+2) &
         + press_interp_weights(4)*smix(ii,idp+3)
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
         ii = idtt+i-1
         temp_work(i,1) = press_interp_weights(1)*tablex(ii,idp,9) + &
         press_interp_weights(2)*tablex(ii,idp+1,9) + &
         press_interp_weights(3)*tablex(ii,idp+2,9) &
         + press_interp_weights(4)*tablex(ii,idp+3,9)
         temp_work(i,2) = press_interp_weights(1)*tabley(ii,idp,9) + &
         press_interp_weights(2)*tabley(ii,idp+1,9) + &
         press_interp_weights(3)*tabley(ii,idp+2,9) &
         + press_interp_weights(4)*tabley(ii,idp+3,9)
         temp_work(i,3) = press_interp_weights(1)*tablez(ii,idp,7) + &
         press_interp_weights(2)*tablez(ii,idp+1,7) + &
         press_interp_weights(3)*tablez(ii,idp+2,7) &
         + press_interp_weights(4)*tablez(ii,idp+3,7)
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
         ii = idtt+i-1
         temp_work(i,1) = press_interp_weights(1)*tablex(ii,idp,2) + &
         press_interp_weights(2)*tablex(ii,idp+1,2) + &
         press_interp_weights(3)*tablex(ii,idp+2,2) &
         + press_interp_weights(4)*tablex(ii,idp+3,2)
         temp_work(i,2) = press_interp_weights(1)*tabley(ii,idp,2) + &
         press_interp_weights(2)*tabley(ii,idp+1,2) + &
         press_interp_weights(3)*tabley(ii,idp+2,2) &
         + press_interp_weights(4)*tabley(ii,idp+3,2)
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
         ii = idtt+i-1
         temp_work(i,1) = press_interp_weights(1)*tablex(ii,idp,3) + &
         press_interp_weights(2)*tablex(ii,idp+1,3) + &
         press_interp_weights(3)*tablex(ii,idp+2,3) &
         + press_interp_weights(4)*tablex(ii,idp+3,3)
         temp_work(i,2) = press_interp_weights(1)*tabley(ii,idp,3) + &
         press_interp_weights(2)*tabley(ii,idp+1,3) + &
         press_interp_weights(3)*tabley(ii,idp+2,3) &
         + press_interp_weights(4)*tabley(ii,idp+3,3)
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
