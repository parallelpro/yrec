!----------------------------------------------------------------------
! setscv
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setscv.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This is to evaluate the EOS with the new tables and call eqstat to
! compute the new equation of state and compare to the old one...
!
! Builds the envelope-mixture (X/Y/Z) SCVH-style equation-of-state
! table tablenv from the separately tabulated pure-hydrogen (tablex),
! pure-helium (tabley), and pure-metal (tablez) tables: first
! completes du/dt and the entropy-of-mixing correction (smix) on the
! native H/He table grid, then combines the three species by the
! additive-volume rule at the envelope composition (X=envelope_hydrogen_
! fraction, Z=star%env_comp%envelope_metal_fraction) to fill tablenv columns 1-6,
! then numerically differentiates those columns (in log T and log P)
! to fill tablenv columns 7-12.
subroutine setscv

      use star_info_lib, only: star
      use const_lib
      use numerics_lib
      use scv_eos_lib
      implicit none
      integer, parameter :: nts = 63, nps = 76





      double precision :: interp_x(3), t_interp_weight(3), &
           t_interp_dweight(3), p_interp_weight(3), p_interp_dweight(3)
      double precision :: hydrogen_atom_mass, helium_atom_mass, &
           boltzmann_constant
      data hydrogen_atom_mass, helium_atom_mass, boltzmann_constant &
           /1.67357d-24, 6.646442d-24, 1.380658d-16/
! --- locals ---
      integer :: t_idx, p_idx, k_idx, ii, jj
! log_t_work/log_p_work are reused as scratch across the file: log10
! scale in the first two loops below, natural-log scale (CLN*...) in
! the final derivative loop -- this mirrors the original TL/PL reuse.
      double precision :: log_t_work, log_p_work, ln_t_work, temp_value
      double precision :: helium_fraction_local, one_minus_y_local, &
           hydrogen_fraction_local, metal_fraction_local
      double precision :: density_h, density_he, density_z, density_mix
      double precision :: hydrogen_number_density, helium_number_density, &
           hydrogen_electron_density, helium_electron_density, &
           total_electron_density
      double precision :: mixing_beta, mixing_gamma, hydrogen_xd, &
           helium_yd, mix_de, mix_dei, total_number_density
      double precision :: pressure_value, radiation_pressure, gas_pressure, &
           rho_mix, log_rho_mix
      double precision :: dlnrho_dlnp_gas, dlnrho_dlnt_gas, dlnp_dlnt_gas
      double precision :: entropy_h, entropy_he, entropy_total, &
           dsmix_dlnt, dlns_dlnt, cp_gas, du_dt
      double precision :: dqdt_dlnp, dqdt_dlnt, dlncp_dlnp, dlncp_dlnt, &
           dqut_dlnt, dqut_dlnp, dlnp_dlnrho, ln_rho_value

!  READ IN EQUATION OF STATE TABLES FOR HYDROGEN AND HELIUM
      do t_idx = 1, nts
         do p_idx = 1, nptsx(t_idx)
!  COMPUTE DU/DT - NEEDED FOR SPECIFIC HEAT DERIVATIVES
!  CP = T (DS/DT)|P = S(D LN S/D LN T)|P (EQ 1)
! AND CP = DU/DT + P*(D LN RHO/D LN T)**2/(RHO*T*D LN RHO/D LN P) (EQ 2)
! DU/DT CAN BE INFERRED FROM EQ 1 AND THE SECOND TERM IN EQ 2, ALL OF
! WHICH ARE CONTAINED IN THE TABLES.
            tablex(t_idx,p_idx,12) = &
                 exp(ln10*tablex(t_idx,p_idx,5))*tablex(t_idx,p_idx,9) - &
                 exp(ln10*(tablex(t_idx,p_idx,1)-tablex(t_idx,p_idx,4)- &
                 tlogx(t_idx)))*tablex(t_idx,p_idx,7)**2/tablex(t_idx,p_idx,8)
            tabley(t_idx,p_idx,12) = &
                 exp(ln10*tabley(t_idx,p_idx,5))*tabley(t_idx,p_idx,9) - &
                 exp(ln10*(tabley(t_idx,p_idx,1)-tabley(t_idx,p_idx,4)- &
                 tlogx(t_idx)))*tabley(t_idx,p_idx,7)**2/tabley(t_idx,p_idx,8)
! COMPUTE THE ENTROPY OF MIXING.
! NUMBER DENSITY OF HYDROGEN AND HELIUM
            helium_fraction_local = 1.0d0-star%env_comp%envelope_hydrogen_fraction
            one_minus_y_local = 1.0d0 - helium_fraction_local
            density_h = exp(ln10*tablex(t_idx,p_idx,4))
            density_he = exp(ln10*tabley(t_idx,p_idx,4))
            density_mix = 1.0d0/(one_minus_y_local/density_h + &
                 helium_fraction_local/density_he)
            hydrogen_number_density = 2.0d0*one_minus_y_local*density_mix/ &
                 hydrogen_atom_mass/(1.0d0+3.0d0*tablex(t_idx,p_idx,2)+ &
                 tablex(t_idx,p_idx,3))
            helium_number_density = 3.0d0*helium_fraction_local*density_mix/ &
                 helium_atom_mass/(1.0d0 + 2.0d0*tabley(t_idx,p_idx,2) + &
                 tabley(t_idx,p_idx,3))
! HYDROGEN AND HELIUM ELECTRON NUMBER DENSITIES
            hydrogen_electron_density = 0.5d0*hydrogen_number_density* &
                 (1.0d0 - tablex(t_idx,p_idx,2) - tablex(t_idx,p_idx,3))
            helium_electron_density = cc13*helium_number_density* &
                 (2.0d0 - 2.0d0*tabley(t_idx,p_idx,2) - tabley(t_idx,p_idx,3))
            total_electron_density = hydrogen_electron_density + &
                 helium_electron_density
            hydrogen_electron_density = max(hydrogen_electron_density,1.0d0)
            helium_electron_density = max(helium_electron_density,1.0d0)
            total_electron_density = max(total_electron_density,1.0d0)
            if (hydrogen_number_density.eq.0.0d0 .or. &
                 helium_number_density.eq.0.0d0) then
               smix(t_idx,p_idx) = 0.0d0
            else
               mixing_beta = (hydrogen_atom_mass/helium_atom_mass)* &
                    (helium_fraction_local/one_minus_y_local)
               mixing_gamma = 1.5d0*(1.0d0 + tablex(t_idx,p_idx,3)+ 3.0d0* &
                    tablex(t_idx,p_idx,2))/(1.0d0+2.0d0*tabley(t_idx,p_idx,2)+ &
                    tabley(t_idx,p_idx,3))
               hydrogen_xd = 1.0d0 - tablex(t_idx,p_idx,2) - tablex(t_idx,p_idx,3)
               helium_yd = 2.0d0 - 2.0d0*tabley(t_idx,p_idx,2) - &
                    tabley(t_idx,p_idx,3)
               if (hydrogen_xd.gt.0.0d0) then
                 mix_de = 1.5d0*helium_yd*mixing_beta*mixing_gamma/hydrogen_xd
               else
                 mix_de = 0.0d0
               end if
               if (helium_yd.gt.0.0d0) then
                  mix_dei = hydrogen_xd/(1.5d0*helium_yd*mixing_beta*mixing_gamma)
               else
                  mix_dei = 0.0d0
               end if
               total_number_density = hydrogen_number_density + &
                    helium_number_density
               smix(t_idx,p_idx) = one_minus_y_local/hydrogen_atom_mass*2.0d0/ &
                    (1.0d0+tablex(t_idx,p_idx,3)+3.0d0*tablex(t_idx,p_idx,2))* &
                    (dlog(1.0d0+mixing_beta*mixing_gamma)- &
                    hydrogen_electron_density/total_number_density* &
                    dlog(1.0d0 + mix_de)+mixing_beta*mixing_gamma* &
                    dlog(1.0d0+1.0d0/mixing_beta/mixing_gamma)- &
                    helium_electron_density/total_number_density* &
                    dlog(1.0d0+mix_dei))
               smix(t_idx,p_idx) = boltzmann_constant*smix(t_idx,p_idx)
            end if
         end do
      end do
!  NOW COMPUTE EQUATION OF STATE VARIABLES FOR THE SURFACE MIXTURE.
      do t_idx=1,nts
         log_t_work = tlogx(t_idx)
         temp_value = (10.0d0**log_t_work)
! TEMPERATURE INTERPOLATION FACTORS
        if (t_idx.eq.1) then
           idtt = 1
        else if (t_idx.eq.nts) then
           idtt = t_idx - 2
        else
           idtt = t_idx - 1
        end if
        do k_idx = 1,3
           interp_x(k_idx) = ln10*tlogx(k_idx+idtt-1)
        end do
        ln_t_work = ln10*log_t_work
        call inter3(interp_x,t_interp_weight,t_interp_dweight,ln_t_work)
         do p_idx=1,nptsx(t_idx)
            log_p_work = tablex(t_idx,p_idx,1)
            log_rho_mix = tablex(t_idx,p_idx,4)
            hydrogen_fraction_local = star%env_comp%envelope_hydrogen_fraction
            metal_fraction_local = star%env_comp%envelope_metal_fraction
            helium_fraction_local = 1.0d0 - hydrogen_fraction_local - &
                 metal_fraction_local
            pressure_value = (10.0d0**log_p_work)
! INCLUDE RADIATION PRESSURE IN THE EQUATION OF STATE.
            radiation_pressure = radiation_constant_over_3*(temp_value**2)**2
            gas_pressure = pressure_value
            pressure_value = (pressure_value + radiation_pressure)
            log_p_work = log10(pressure_value)
!  DENSITY : ADD INVERSELY : 1/RHO = X/RHO(X) + Y/RHO(Y) + Z/RHO(Z)
            density_h = exp(ln10*tablex(t_idx,p_idx,4))
            density_he = exp(ln10*tabley(t_idx,p_idx,4))
            density_z = exp(ln10*tablez(t_idx,p_idx,4))
            rho_mix = 1.0d0/(hydrogen_fraction_local/density_h + &
                 helium_fraction_local/density_he + &
                 metal_fraction_local/density_z)
            log_rho_mix = log10(rho_mix)
! D LN RHO/ D LN P = QDP
! FOR GAS PRESSURE QDP(TOT) = QDP(X)*X*RHO/RHO(X)+QDP(Y)*Y*RHO/RHO(Y)
! FOR RADIATION PRESSURE QDP = 0, SO QDP(TOT) = QDP(GAS)*P/PGAS
            dlnrho_dlnp_gas = tablex(t_idx,p_idx,8)*hydrogen_fraction_local* &
                 (rho_mix/density_h)+tabley(t_idx,p_idx,8)* &
                 helium_fraction_local*(rho_mix/density_he) &
                 + tablez(t_idx,p_idx,13)*metal_fraction_local* &
                 (rho_mix/density_z)
! D LN RHO/ D LN T = QDT (NOTE : D LN P/ D LN T = QPT)
! FOR GAS PRESSURE, CORRECT AS PER QDP
! FOR RADIATION PRESSURE, USE QDT = QDP*QPT.  CORRECT QPT FOR
! RADIATION PRESSURE AND USE THE CORRECTED QDP, QPT TO GET QDT.
            dlnrho_dlnt_gas = hydrogen_fraction_local*(rho_mix/density_h)* &
                 tablex(t_idx,p_idx,7)+helium_fraction_local* &
                 (rho_mix/density_he)*tabley(t_idx,p_idx,7) &
                 +metal_fraction_local*(rho_mix/density_z)* &
                 tablez(t_idx,p_idx,10)
            dlnp_dlnt_gas = -dlnrho_dlnt_gas/dlnrho_dlnp_gas
            entropy_h = exp(ln10*tablex(t_idx,p_idx,5))
            entropy_he = exp(ln10*tabley(t_idx,p_idx,5))
! ENTROPY.  OBEYS ADDITIVE VOLUME RULE, BUT ALSO NEED TO INCLUDE THE
! ENTROPY OF MIXING.
            entropy_total = hydrogen_fraction_local*entropy_h + &
                 helium_fraction_local*entropy_he + smix(t_idx,p_idx)
            dsmix_dlnt = t_interp_dweight(1)*smix(idtt,p_idx)+ &
                 t_interp_dweight(2)*smix(idtt+1,p_idx)+ &
                 t_interp_dweight(3)*smix(idtt+2,p_idx)
! D LN S/ D LN T
            dlns_dlnt = (hydrogen_fraction_local*entropy_h* &
                 tablex(t_idx,p_idx,9) + helium_fraction_local*entropy_he* &
                 tabley(t_idx,p_idx,9) + dsmix_dlnt)/entropy_total
!  CP = S*(D LN S/ D LN T)|P IS TABULATED. USE
!  CP = DU/DT + P*(D LN RHO/D LN T)**2/RHO/T/(D LN RHO/ D LN P)
!  TO INCLUDE THE EFFECTS OF RADIATION PRESSURE.
! CP (GAS PRESSURE ONLY).
            cp_gas = entropy_total*dlns_dlnt + metal_fraction_local* &
                 (exp(ln10*tablez(t_idx,p_idx,7)) + &
                 gas_pressure*tablez(t_idx,p_idx,10)**2/ &
                 tablez(t_idx,p_idx,13)/density_z/temp_value)
! NOW FIND DU/DT FROM THE ORIGINAL TABLE.
            du_dt = cp_gas - gas_pressure*dlnrho_dlnt_gas**2/ &
                 dlnrho_dlnp_gas/rho_mix/temp_value
! NOW STORE THE RELEVANT VARIABLES IN A TABLE FOR THE ENVELOPE
! MIXTURE.
            tablenv(t_idx,p_idx,1) = tablex(t_idx,p_idx,1)
            tablenv(t_idx,p_idx,2) = log_rho_mix
            tablenv(t_idx,p_idx,3) = dlnrho_dlnt_gas
            tablenv(t_idx,p_idx,4) = dlnrho_dlnp_gas
            tablenv(t_idx,p_idx,5) = cp_gas
            tablenv(t_idx,p_idx,6) = du_dt
         end do
      end do
! NOW FIND THE FOLLOWING DERIVATIVES NUMERICALLY :
! D QDT/ D LN T
! D QDT/ D LN P
! D LN CP / D LN P
! D LN CP / D LN T
! D LN DEL(AD) / D LN P
! D LN DEL(AD) / D LN T
      do t_idx = 1, nts
! TEMPERATURE INTERPOLATION FACTORS
        if (t_idx.eq.1) then
           idtt = 1
        else if (t_idx.eq.nts) then
           idtt = t_idx - 2
        else
           idtt = t_idx - 1
        end if
        do k_idx = 1,3
           interp_x(k_idx) = ln10*tlogx(k_idx+idtt-1)
        end do
        log_t_work = ln10*tlogx(t_idx)
        call inter3(interp_x,t_interp_weight,t_interp_dweight,log_t_work)
        temp_value = exp(log_t_work)
        do p_idx = 1, nptsx(t_idx)
! PRESSURE INTERPOLATION FACTORS
           if (p_idx.eq.1) then
              idp = 1
           else if (p_idx.eq.nptsx(t_idx)) then
              idp = p_idx - 2
           else
              idp = p_idx - 1
           end if
           do k_idx = 1,3
              interp_x(k_idx) = ln10*tablenv(t_idx,idp+k_idx-1,1)
           end do
           log_p_work = ln10*tablenv(t_idx,p_idx,1)
           call inter3(interp_x,p_interp_weight,p_interp_dweight,log_p_work)
! DERIVATIVES OF D LN RHO/ D LN T
           dqdt_dlnp = tablenv(t_idx,idp,3)*p_interp_dweight(1)+ &
                tablenv(t_idx,idp+1,3)*p_interp_dweight(2) &
                +tablenv(t_idx,idp+2,3)*p_interp_dweight(3)
           dqdt_dlnt = tablenv(idtt,p_idx,3)*t_interp_dweight(1)+ &
                tablenv(idtt+1,p_idx,3)*t_interp_dweight(2) &
                +tablenv(idtt+2,p_idx,3)*t_interp_dweight(3)
! DERIVATIVES OF LN CP
           do k_idx = 1,3
              interp_x(k_idx) = dlog(tablenv(t_idx,idp+k_idx-1,5))
           end do
           dlncp_dlnp = interp_x(1)*p_interp_dweight(1)+ &
                interp_x(2)*p_interp_dweight(2)+interp_x(3)*p_interp_dweight(3)
           do k_idx = 1,3
              ii = idtt+k_idx-1
            jj = min(nptsx(ii),p_idx)
            interp_x(k_idx) = dlog(tablenv(ii,jj,5))
           end do
           dlncp_dlnt = interp_x(1)*t_interp_dweight(1)+ &
                interp_x(2)*t_interp_dweight(2)+interp_x(3)*t_interp_dweight(3)
! DERIVATIVES OF DU/DT
           do k_idx = 1,3
              ii = idtt+k_idx-1
            jj = min(nptsx(ii),p_idx)
            interp_x(k_idx) = tablenv(ii,jj,6)
              interp_x(k_idx) = tablenv(idtt+k_idx-1,p_idx,6)
           end do
           dqut_dlnt = interp_x(1)*t_interp_dweight(1)+ &
                interp_x(2)*t_interp_dweight(2)+interp_x(3)*t_interp_dweight(3)
           do k_idx = 1,3
              interp_x(k_idx) = tablenv(t_idx,idp+k_idx-1,6)
           end do
           dqut_dlnp = interp_x(1)*p_interp_dweight(1)+ &
                interp_x(2)*p_interp_dweight(2)+interp_x(3)*p_interp_dweight(3)
! D/D LN RHO (D LN P/D LN RHO)
           do k_idx = 1,3
              interp_x(k_idx) = ln10*tablenv(t_idx,idp+k_idx-1,2)
           end do
           ln_rho_value = ln10*tablenv(t_idx,p_idx,2)
           call inter3(interp_x,p_interp_weight,p_interp_dweight,ln_rho_value)
           dlnp_dlnrho = p_interp_dweight(1)/tablenv(t_idx,idp,4)+ &
                p_interp_dweight(2)/tablenv(t_idx,idp+1,4) &
                +p_interp_dweight(3)/tablenv(t_idx,idp+2,4)
           tablenv(t_idx,p_idx,7) = dqdt_dlnp
           tablenv(t_idx,p_idx,8) = dqdt_dlnt
           tablenv(t_idx,p_idx,9) = dlnp_dlnrho
           tablenv(t_idx,p_idx,10) = dlncp_dlnt
           tablenv(t_idx,p_idx,11) = dqut_dlnp
           tablenv(t_idx,p_idx,12) = dqut_dlnt
        end do
      end do

      return
end subroutine setscv
