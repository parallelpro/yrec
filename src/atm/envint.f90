!----------------------------------------------------------------------
! envint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original envint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Computes one envelope solution (P, T, R at the model's fitting
! point) for a given (Teff, L) vertex: first integrates the gray (or
! Eddington/Krishna-Swamy) atmosphere in optical depth from a starting
! guess down to tau=2/3 via QATM (bypassed entirely for tabulated
! Kurucz/Allard atmospheres, which look up the boundary pressure
! directly via SURFP/KCSURFP/ALSURFP), then integrates the envelope
! structure equations in pressure from tau=2/3 to the fitting mass
! point via QENV. Both integrations use the Bulirsch-Stoer stepper
! BSSTEP. Also saves the full atmosphere/envelope structure (for
! profile output and pulsation) and, unless the caller has switched to
! the newer TAUINTNEW-based method, locates the surface convection
! zone within the envelope via TAUCAL.
subroutine envint(luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, log10_star_mass, &
     vertex_index, print_flag, save_boundary_flag, log10_pressure_limit, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     stored_envelope_state, stored_vertex_index, atm_call_count, &
     env_call_count, saha_state, vtx_logp, vtx_logr, vtx_logt, &
     pulse_print_flag)

      use eos_lib
      use atm_table_lib
      use rotdiff_lib
      use run_diag_lib
      use atmstruct_lib
      use pulse_diag_lib
      use envstruct_lib
      use envelope_comp_lib
      use light_burn_lib
      use turnover_lib
      use luout_lib
      use const_lib
      use intpar_lib
      use numerics_lib
      implicit none
      integer, parameter :: json=5000
! PARAMETERS NT AND NG FOR TABULATED SURFACE PRESSURES OF KURUCZ.
      integer, parameter :: nt=57, ng=11
! JNT 06/14 ADDED NTC/NGC FOR KTTAU=5
      integer, parameter :: ntc=76, ngc=11
! MHP 8/97 ADDED NTA AND NGA FOR ALLARD ATMOSPHERE TABLES
      integer, parameter :: nta=54, nga=5

! luminosity_linear/pressure_rotation_factor/temperature_rotation_factor/
! log10_gravity/print_flag/hydrogen_fraction/metal_fraction are
! intent(inout), not intent(in): they're relayed unchanged through
! this file's own calls to bsstep (numerics_lib) down to mmid and then
! to the deriv callback (qatm/qenv/etc, whose actual identity varies
! per call site) -- both bsstep and mmid declare their corresponding
! dummy arguments intent(inout), since the callback may modify them.
! Declaring them intent(in) here was a pre-existing bug (silently
! tolerated under the old implicit-interface calling convention,
! surfaced once bsstep moved into numerics_lib and gained an explicit
! interface -- see GUIDELINES.md); every call site already passes real
! variables for these positions, so widening the intent here changes
! nothing about how any of them are called.
      double precision, intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity
      double precision, intent(in) :: log10_star_mass
      integer, intent(in) :: vertex_index
      logical, intent(inout) :: print_flag
      logical, intent(in) :: save_boundary_flag
      double precision, intent(in) :: log10_pressure_limit
      double precision, intent(inout) :: log10_radius
      double precision, intent(inout) :: log10_teff
      double precision, intent(inout) :: hydrogen_fraction, metal_fraction
      double precision, intent(inout) :: stored_envelope_state(4)
      integer, intent(inout) :: stored_vertex_index
      integer, intent(inout) :: atm_call_count, env_call_count, saha_state
      double precision, intent(inout) :: vtx_logp(3), vtx_logr(3), vtx_logt(3)
      logical, intent(in) :: pulse_print_flag

! DBG CHANGED MAXSTEP FROM 200 TO 2000 TO GIVE ATMOSPHERE INTEGRATER A CHANCE.
      integer, parameter :: maxstp = 2000
      double precision, parameter :: tiny = 1.0d-30
      external qatm, qenv
      double precision :: hra
      external hra

      double precision :: ion_fraction(3)
      double precision :: taucal_delta_mass(json), taucal_shell_mass(json), &
           taucal_local_gravity(json), taucal_radiative_gradient(json), &
           taucal_adiabatic_gradient(json)


! G Somers END


! MHP 1/01 CHANGED END OF FILE INDICATOR IN ATMOSPHERE/ENVELOPE FILES TO
! VECTOR FROM SCALAR
      double precision :: xyz(22)
      data xyz/22*99.99d0/
      save

! --- locals ---
      logical :: want_derivatives, in_atmosphere, conductive_opacity_flag
      logical :: allard_lookup_failed
      integer :: jj
      double precision :: log10_temperature, temperature, log10_pressure, &
           pressure, log10_density, density
      double precision :: atm_density_guess
      double precision :: beta, beta_inverse, beta14, specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision :: dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
           adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      double precision :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      double precision :: indep_var
      double precision :: chi_rho, chi_t, specific_heat_cv, gamma1
      double precision :: pulse_energy_sum, prev_tau, prev_opacity, &
           prev_density, delta_tau_step, electron_pressure, &
           pulse_radiative_gradient, pulse_gradient, opacity_now, &
           density_now, tau_now
      integer :: num_eqs, num_ok, num_bad, step_index, i
      double precision :: tolerance, x_limit, h_max, h_min, h_step, &
           step_tolerance, h_did, h_next
      double precision :: y(3), dydx(3), y_scale(3), y_start(3), &
           err_sum(3), step_err(3)
      double precision :: initial_mass_coord, surface_radius_linear
      logical :: store_flag_set, pressure_limit_test_flag
      logical :: cz_in_envelope
      double precision :: mass_diff_remaining, interp_weight
      integer :: inversion_index1, inversion_index2
      double precision :: swap_temp
      logical :: swap_temp_logical
      integer :: cz_start_index
      double precision :: log10_cz_radius
      double precision :: dd2, dd1, interp_fraction
      double precision :: unused_chdelj, unused_chdeld
      double precision :: radius_linear, depth_fraction, local_log10_gravity, &
           local_gravity_linear
      double precision :: x_start, taucz_env_accum, delta_radius_cz

! DBG PULSE TURN ON DERIVATIVE CALCULATOR
      if (pulse_print_flag.and.print_flag) then
          pulse_diag%lpumod = .true.
          want_derivatives = .true.
          in_atmosphere = .false.
      else
          pulse_diag%lpumod = .false.
          want_derivatives = .false.
          in_atmosphere = .true.
      end if
      conductive_opacity_flag = .false.

! JVS 10/07/13 Always calculate derivatives
      want_derivatives = .true.

! G Somers 11/14 WRITE ATMOSHPHERE HEADER TO .STORE FILE, AND ADDED
! I/O FLAGS TO THE ATMOSPHERE CALLS
      if(print_flag.and.lstatm)then
         if(.not.lstch)write(istor,60)
      endif
 60   format(/,'******** ATMOSPHERE BEGIN ********')

! GET PRESSURE AT T=Teff BY INTERPOLATION IN TABLE ATMPL.
      if (atm_choice .eq. 3) then
! KURUCZ ATMOSPHERES
         if(lstch) lstatm=.false.
         call surfp(log10_teff,log10_gravity,print_flag.and.lstatm)
         goto 200
! JNT 06/14
! GET PRESSURE AT T=Teff BY INTERPOLATION IN TABLE ATMPLC.
      else if (atm_choice .eq. 5) then
! KURUCZ ATMOSPHERES
         if(lstch) lstatm=.false.
         call kcsurfp(log10_teff,log10_gravity,print_flag.and.lstatm)
         goto 200
! We have Kurucz atmosphere boundary conditions
      else if (atm_choice .eq. 4) then
! ALLARD & HAUSCHILDT ATMOSPHERES
         if(lstch) lstatm=.false.
         call alsurfp(log10_teff,log10_gravity,print_flag.and.lstatm,allard_lookup_failed)
! Changed to Allard atmosphere code
         if(allard_lookup_failed) then
            atm_choice=0
            use_ttau_relation = .true.
! Set to gray atmosphere (KTTAU=0), as
! TeffL is above Allard max, or GL is out of range.
            write(*,*) 'ENVINT: Change to gray atmosphere (KTTAU=0)'
            write(short_file_unit,*)'ENVINT: Change to gray atmosphere (KTTAU=0)'
            goto 2
         endif
         goto 200
! We have Allard atmosphere boundary conditions
      endif
    2 continue
! Start gray atmosphere bounary conditions
! GUESS THE TEMPERATURE FOR AN OPTICAL DEPTH NEAR ZERO.
      idt = 15
      do 5 jj = 1,4
       idd(jj) = 5
    5 continue
      err_sum(1) = 0.0d0
      if(atm_choice .eq. 0) then
            log10_temperature = log10_teff - 0.031235d0 + 0.25d0*dlog10(cc23)
      else if (atm_choice .eq. 1) then
            log10_temperature = log10_teff - 0.031235d0 + 0.25d0*dlog10(0.550d0)
      else if (atm_choice .eq. 2) then
            log10_temperature = log10_teff + hra(cc23) - atm_hras
      end if
!                 For kttau = 0,1,or 2, very occasionally the integration
!                 fails because the starting point (X0) is past the end
!                 point (XLIM). When this happens, we divide the effective
!                 starting density (atmd0) by 10 and retry.
      atm_density_guess = atm_step_initial
 1998 continue
! Return point if X0 > XLIM
      temperature = dexp(ln10*log10_temperature)
! FIND THE PRESSURE CORRESPONDING TO THIS T AND THE DENSITY CHOSEN
! FOR THE START OF THE ATMOSPHERE INTEGRATION.

!      PL = DLOG10((CGAS*ATMD0 + CA3*T**3)*T)
      log10_pressure = dlog10((gas_constant*atm_density_guess + &
           radiation_constant_over_3*temperature**3)*temperature)

! NOW FIND THE OPTICAL DEPTH(X0) WHERE THE ATMOSPHERE INTEGRATION BEGINS.
      call eos_get(log10_temperature,temperature,log10_pressure,pressure, &
           log10_density,density,hydrogen_fraction,metal_fraction,beta, &
           beta_inverse,beta14,ion_fraction,specific_gas_constant, &
           ion_mean_weight_inverse,electron_mean_weight_inverse, &
           electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
           specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt,adiabatic_gradient_dt,adiabatic_gradient_dp, &
           specific_heat_cp_dt,specific_heat_cp_dp,want_derivatives, &
           in_atmosphere,saha_state)
! DBG 12/95 GET OPACITY
      call getopac(log10_density, log10_temperature, hydrogen_fraction, &
           metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, &
           dlnkap_dlnt, ion_fraction)
      indep_var = log10_pressure - log10_gravity + dlog10(opacity)
      y(1) = log10_pressure
      dydx(1) = dexp(ln10*(log10_gravity+indep_var-log10_opacity-log10_pressure))
! G Somers 11/14 ADDED I/O FLAG AND CHANGED WRITE OUTS TO .STORE.
      if(print_flag.and.lstatm) then
       if(.not.lstch) write(istor,10)
         beta = 1.0d0 - radiation_constant_over_3*(temperature**2)**2/pressure
         chi_rho = 1.0d0/dlnrho_dlnp
         chi_t = -chi_rho*dlnrho_dlnt
         specific_heat_cv = specific_heat_cp - exp(ln10*(log10_pressure-log10_density-log10_temperature))*chi_t**2/chi_rho
         gamma1 = chi_rho*specific_heat_cp/specific_heat_cv
           if(.not.lstch)then
             write(istor,20)indep_var,log10_pressure,log10_temperature,log10_density,opacity, &
                         (ion_fraction(i),i=1,3),saha_state,atm_call_count, &
                         gamma1,dlnrho_dlnp,dlnrho_dlnt,beta,specific_heat_cp,specific_heat_cv
           endif
   10    format(6X,'TAU',9X,'P',10X,'T',10X,'D',11X,'O',9X, &
              'HII  HEII HEIII   SAHA   KATM')
   20    format(1X,4F11.7,1PE14.7,0P3F6.3,2I7,4F12.8,1P2E12.5)
      endif
! DBG PULSE INITIAL POINT FOR PULSATION
      if (pulse_print_flag.and.print_flag) then
        pulse_diag%qqed = 0.0d0
        pulse_energy_sum = 0.0d0
        pulse_diag%qqet = 0.0d0
        pulse_diag%qfs = 1.0d0
        prev_tau = dexp(ln10*indep_var)
        prev_opacity = opacity
        prev_density = dexp(ln10*log10_density)
        delta_tau_step = 0.0d0
        electron_pressure = gas_constant * temperature * density * electron_mean_weight_inverse
! FROM FIRST LINES OF TPGRAD
        pulse_radiative_gradient = opacity*luminosity_linear*dexp(ln10*(log10_pressure-log10_star_mass-4.0d0*log10_temperature+ &
                 log10_solar_luminosity-cgl+cdelrl))*temperature_rotation_factor/pressure_rotation_factor
        if (pulse_radiative_gradient-adiabatic_gradient .le. 1.0d-6) then
            pulse_gradient = pulse_radiative_gradient
        else
            pulse_gradient = adiabatic_gradient
        end if
        if(pulsation_file_version.eq.1) then
! MHP 10/02 TYPO - QED SHOULD HAVE BEEN QQED
!     *         PL, QESUM,O,QDP,QED,
           write(opal_atm_unit,5001)delta_tau_step,pulse_diag%qfs,luminosity_linear,log10_temperature,log10_density, &
               log10_pressure, pulse_energy_sum,opacity,dlnrho_dlnp,pulse_diag%qqed, &
               pulse_diag%qqet,dlnkap_dlnrho,dlnkap_dlnt,pulse_gradient,adiabatic_gradient, &
               specific_heat_cp,specific_gas_constant,dlnrho_dlnt,electron_pressure
!     *         PL, TAUP ,O,QDP,QED,
          else if (pulsation_file_version.eq.2 .or. pulsation_file_version.eq.3) then
           write(opal_atm_unit,6001)delta_tau_step,pulse_diag%qfs,luminosity_linear,log10_temperature,log10_density, &
               log10_pressure, prev_tau ,opacity,dlnrho_dlnp,pulse_diag%qqed, &
               pulse_diag%qqet,dlnkap_dlnrho,dlnkap_dlnt,pulse_gradient,adiabatic_gradient, &
               specific_heat_cp,specific_gas_constant,dlnrho_dlnt,electron_pressure
          end if

 5001 format(5E16.9,/,5E16.9,/,5E16.9,/,4E16.9)
 6001 format(5E23.16,/,5E23.16,/,5E23.16,/,4E23.16)
 6003 format(6E23.16,/,6E23.16,/,6E23.16,/,4E23.16)
      end if



! INTEGRATE DP/DTAU FROM THIS STARTING TAU TO TAU = 2/3.
! SET NUMERICAL PARAMETERS UP.
      num_eqs = 1
      tolerance = atm_error_tol
      if (atm_choice .eq. 1) then
! KRISHNA-SWAMY T TAU HAS DIFFERENT ZERO THAN EDDINGTON T TAU
! TAU = 0.312156330 AT TEFF.
            x_limit = -0.505627854d0
      else
! TAU = 2/3 AT TEFF.
            x_limit = -0.176091259d0
      end if

      if  (indep_var .gt. x_limit) then ! Check that starting point is before endpoint
         write(short_file_unit,*)"ENVINT: X0>XLIM, X0,XLIM:",indep_var,x_limit
       write(short_file_unit,*)"ENVINT: get new X0 by dividing ATMD0 by 10"
       atm_density_guess = atm_density_guess / 10d0   ! If not before, divide starting density by 10
       goto 1998              ! and retry.
      endif

      h_max = atm_step_max
      h_min = atm_step_min
      h_step = atm_step_begin
      step_tolerance = tolerance_fraction
      num_ok = 0
      num_bad = 0

! INTEGRATION LOOP.
      do 40 step_index = 1,maxstp
! YSCAL IS THE ARRAY THAT THE NUMERICAL ERRORS ARE SCALED AGAINST.
       do 30 i = 1,num_eqs
          y_scale(i) = dabs(y(i)) + dabs(h_step*dydx(i))+tiny
   30    continue
! ENSURE THAT STEP DOESN'T EXCEED MAXIMUM STEP SIZE OR GO PAST
! THE LIMIT OF THE INTEGRATION.
       if((indep_var-x_limit)*(indep_var+h_step-x_limit).lt.0.0d0) h_step = x_limit - indep_var
       if(h_step.gt.h_max) h_step = h_max
! INTEGRATE THE ATMOSPHERE FROM TAU TO TAU + H
! H IS THE ATTEMPTED STEP,HDID IS THE ONE PERFORMED, AND HNEXT IS THE
! PREDICTED NEXT STEP.
       call bsstep(y,dydx,num_eqs,indep_var,h_step,tolerance,y_scale,h_did, &
            h_next,qatm, luminosity_linear,pressure_rotation_factor, &
            temperature_rotation_factor,log10_gravity,in_atmosphere, &
            want_derivatives,conductive_opacity_flag,print_flag,log10_radius, &
            log10_teff,hydrogen_fraction,metal_fraction,atm_call_count, &
            saha_state,step_err)
! FIND DP/DTAU AT THE START OF THE NEXT STEP.
       err_sum(1) = err_sum(1) + step_err(1)
       call qatm(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
            temperature_rotation_factor,log10_gravity,in_atmosphere, &
            want_derivatives,conductive_opacity_flag,print_flag,log10_radius, &
            log10_teff,hydrogen_fraction,metal_fraction,atm_call_count,saha_state)
! G Somers 11/14 ADDED I/O FLAG AND CHANGED WRITE OUTS TO .STORE.
       if(print_flag.and.lstatm) then
            beta = 1.0d0 - radiation_constant_over_3*exp(ln10*(4.0d0*atm_table%atm_log10_temperature-atm_table%atm_log10_pressure))
            chi_rho = 1.0d0/pulse_diag%qqdp
            chi_t = -chi_rho*pulse_diag%qqdt
            specific_heat_cv = pulse_diag%qqcp - exp(ln10*(atm_table%atm_log10_pressure-atm_table%atm_log10_density-atm_table%atm_log10_temperature))*chi_t**2/chi_rho
            gamma1 = chi_rho*pulse_diag%qqcp/specific_heat_cv
          if(.not.lstch)then
            write(istor,20)atm_table%atm_tau,atm_table%atm_log10_pressure,atm_table%atm_log10_temperature,atm_table%atm_log10_density,atm_table%atm_opacity, &
                     (atm_table%atm_ion_fraction(i),i=1,3), &
                     saha_state,atm_call_count,gamma1,pulse_diag%qqdp,pulse_diag%qqdt,beta,pulse_diag%qqcp,specific_heat_cv
          endif
! JvS: SAVE STRUCTURE TO COMMON BLOCK
          atmo_struct%atmo_log10_pressure(step_index) = atm_table%atm_log10_pressure
          atmo_struct%atmo_log10_temperature(step_index) = atm_table%atm_log10_temperature
          atmo_struct%atmo_log10_density(step_index) = atm_table%atm_log10_density
          atmo_struct%atmo_beta(step_index) = beta
          atmo_struct%atmo_gamma1(step_index) = gamma1
          atmo_struct%atmo_dlnrho_dlnt(step_index) = pulse_diag%qqdt
          atmo_struct%atmo_ion_fraction(1,step_index) = atm_table%atm_ion_fraction(1)
          atmo_struct%atmo_ion_fraction(2,step_index) = atm_table%atm_ion_fraction(2)
          atmo_struct%atmo_ion_fraction(3,step_index) = atm_table%atm_ion_fraction(3)
          atmo_struct%atmo_opacity(step_index) = atm_table%atm_opacity
          atmo_struct%atmo_specific_heat_cp(step_index) = pulse_diag%qqcp
       endif
       if(h_did.eq.h_step) then
          num_ok = num_ok + 1
       else
          num_bad = num_bad + 1
       endif
! DBG PULSE ATMOSPHERE VALUES FOR PULSATION
! JVS 02/11 - Added LCLCD option to IF statement
       if ((pulse_print_flag.and.print_flag) .or. lclcd_placeholder .or. lstch) then
          pulse_diag%qqed = 0.0d0
          pulse_energy_sum = 0.0d0
          pulse_diag%qqet = 0.0d0
          pulse_diag%qfs = 1.0d0
          opacity_now = pulse_diag%qo
          density_now = dexp(ln10*pulse_diag%qdl)
          tau_now = dexp(ln10*atm_table%atm_tau)
! SEE J.P. COX PRINC. OF STELL. STRUC. P590
          delta_tau_step =  (tau_now - prev_tau)/(((density_now*opacity_now)+(prev_density*prev_opacity))/2)
          electron_pressure = gas_constant * pulse_diag%qt * pulse_diag%qd * pulse_diag%qemu
          prev_opacity = opacity_now
!FROM FIRST LINES OF TPGRAD
          pulse_radiative_gradient = pulse_diag%qo*luminosity_linear*dexp(ln10*(pulse_diag%qpl-log10_star_mass-4.0d0*pulse_diag%qtl+log10_solar_luminosity-cgl+ &
                 cdelrl))*temperature_rotation_factor/pressure_rotation_factor
          if (pulse_radiative_gradient-pulse_diag%qdela .le. 1.0d-6) then
            pulse_gradient = pulse_radiative_gradient
          else
            pulse_gradient = pulse_diag%qdela
          end if
          if (pulse_print_flag.and.print_flag) then
                if (pulsation_file_version.eq.1) then
                         write(opal_atm_unit,5001)delta_tau_step,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl, &
                           pulse_diag%qpl, pulse_energy_sum ,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed, &
                           pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_gradient,pulse_diag%qdela, &
                           pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
                        else if (pulsation_file_version.eq.2.or.pulsation_file_version.eq.3) then
                         write(opal_atm_unit,6001)delta_tau_step,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl, &
                           pulse_diag%qpl, tau_now ,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed, &
                           pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_gradient,pulse_diag%qdela, &
                           pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
                       end if
            end if
! JvS SAVE TO COMMON ATMSTRUCT COMMON BLOCK
          atmo_struct%atmo_delta_depth(step_index) = delta_tau_step
          atmo_struct%atmo_gradients(1,step_index) = pulse_radiative_gradient
          atmo_struct%atmo_gradients(2,step_index) = pulse_gradient
          atmo_struct%atmo_gradients(3,step_index) = pulse_diag%qdela
          atmo_struct%num_atm_points = step_index

          prev_tau = tau_now
          prev_density = density_now
       end if
! DBG END


! CHECK IF INTEGRATION COMPLETE
       if(dabs(indep_var - x_limit).le.step_tolerance) then
          write(short_file_unit,35)num_ok,num_bad,err_sum(1)
   35       format(1X,'ATMOSPHERE INTEGRATION COMPLETE',1X, &
                 'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
                 I5/5X,'MAXIMUM RELATIVE ERROR IN P ',1PE22.13)
          goto 200
       endif
       if(h_next.lt.h_min) h_next = h_min
       h_step = h_next
40    continue
! INTEGRATION HAS FAILED TO FINISH IN MAXSTP STEPS;
! PRINT NASTY MESSAGE AND QUIT.
      write(iowr,50)
      write(short_file_unit,50)
   50 format(5X,'ATMOSPHERE INTEGRATION FAILED AFTER MAXSTP',1X, &
             'INTEGRATIONS.I QUIT.')
      stop
! ENVELOPE INTEGRATION
! HERE P IS THE INDEPENDENT VARIABLE AND M,R,AND T ARE
! DEPENDENT VARIABLES.  INTEGRATE FROM TAU = 2/3 TO THE LAST
! MASS POINT IN THE MODEL.


  200 continue   ! Kurucz and Allard (KTTAU=3 and 4) bypass atmosphere
                 !  integration and come here
! G Somers 3/17, IF INTERESTED ONLY IN PPHOT, BREAK HERE.
      turnover%pphot = atm_table%atm_log10_pressure
      if(.not.calc_envelope_flag) goto 555

! G Somers 11/14 WRITE ENVELOPE HEADER
      if(print_flag.and.lstenv)then
         if(.not.lstch) write(istor,61)
      endif
 61   format(/,'******** ENVELOPE BEGIN ********')

! DBG PULSE WRITE END OF DATA INDICATOR
      if (pulse_print_flag.and.print_flag) then
!         XYZ = 99.99D0
         if(pulsation_file_version.eq.1) then
            write(opal_atm_unit, 5001) (xyz(i),i=1,19)
         else if (pulsation_file_version.eq.2.or.pulsation_file_version.eq.3) then
            write(opal_atm_unit, 6001) (xyz(i),i=1,19)
         end if
! 5002    FORMAT(E16.9)
! 6002    FORMAT(E23.16)
      endif
! DBG
!  IF ENVELOPE MASS(SENV) SMALL ENOUGH,SKIP ENVELOPE INTEGRATION.
! DBG 2/92 CHANGED FROM 1.0D-10 to 1.0D-12
      if(env_comp%senv.gt.-1.0d-12) then
       if(save_boundary_flag) then
          vtx_logp(vertex_index) = atm_table%atm_log10_pressure
          vtx_logr(vertex_index) = log10_radius
          vtx_logt(vertex_index) = atm_table%atm_log10_temperature
          if(print_flag)then
            if(.not.lstch)write(istor,230)vtx_logp(vertex_index),vtx_logt(vertex_index),vtx_logr(vertex_index),env_comp%senv
          endif
       endif
 230     format(4X,3F16.12,8X,F16.12)
       goto 300
      endif
!  INITIALIZE VARIABLES AND SET NUMERICAL PARAMETERS.
      in_atmosphere = .false.
      do 235 i = 1,3
       err_sum(i) = 0.0d0
 235  continue
      num_ok = 0
      num_bad = 0
      indep_var = atm_table%atm_log10_pressure
      initial_mass_coord = 0.0d0
      y(1) = initial_mass_coord
      y(2) = atm_table%atm_log10_temperature
      y(3) = log10_radius
      num_eqs = 3
      tolerance = env_error_tol
      h_max = env_step_max
      h_min = env_step_min
      h_step = env_step_begin
      step_tolerance = dabs(tolerance_fraction*env_comp%senv)
      surface_radius_linear = dexp(ln10*log10_radius)
      if(stored_vertex_index.eq.vertex_index) stored_vertex_index = 0
      store_flag_set = .false.
      if (save_boundary_flag) then
         pressure_limit_test_flag = .false.
      else
         pressure_limit_test_flag = .true.
      end if
!  FIND DY/DX AT THE START OF THE STEP.
      call qenv(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
           temperature_rotation_factor,log10_gravity,in_atmosphere, &
           want_derivatives,conductive_opacity_flag,print_flag,log10_radius, &
           log10_teff,hydrogen_fraction,metal_fraction,env_call_count,saha_state)
! DBG PULSE WRITE FIRST POINT OF ENVELOPE
      if (pulse_print_flag.and.print_flag) then
!         REWIND IOPENV
         pulse_diag%qqed = 0.0d0
         pulse_energy_sum = 0.0d0
         pulse_diag%qqet = 0.0d0
         electron_pressure = gas_constant * pulse_diag%qt * pulse_diag%qd * pulse_diag%qemu
         if(pulsation_file_version.eq.1) then
          write(opal_envelope_unit,5001)log10_radius,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl, &
                 pulse_diag%qpl, pulse_energy_sum,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed, &
                 pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_diag%qdel,pulse_diag%qdela, &
                 pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
         else if (pulsation_file_version.eq.2) then
          write(opal_envelope_unit,6001)log10_radius,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl, &
                 pulse_diag%qpl, pulse_energy_sum,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed, &
                 pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_diag%qdel,pulse_diag%qdela, &
                 pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
         else if (pulsation_file_version.eq.3) then
! DBG 7/95 Appended mixing length info at end of first three lines
          write(opal_envelope_unit,6003)log10_radius,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl,rot_diff%alfmlt, &
                 pulse_diag%qpl, pulse_energy_sum,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed,rot_diff%phmlt, &
                 pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_diag%qdel,pulse_diag%qdela,rot_diff%cmxmlt, &
                 pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
         end if
      end if   !uncertain!
! DBG
! G Somers 11/14 ADDED I/O FLAG AND CHANGED WRITE OUTS TO .STORE.
      if(print_flag.and.lstenv) then
       if(.not.lstch) write(istor,240) 'GRAV  ','P   ','T   ','DEPTH    ','M      ', &
             'D    ','O   ','BETA','DELR  ','DELA','DEL ','HII ', &
             'HEII','HEIII','V   ','GAM1   ','QQDP   '
 240     format(1X,3A10,2A14,2A10,A7,A9,5A6,A9,2A12)
!     X' ',4X,'GRAV',7X,'P',8X,'T',8X,'DEPTH',6X,'M',
!     *        9X,'D',6X,'O',7X,'BETA',3X,'DELR',5X,'DELA',3X,
!     *        'DEL',3X,'HII',3X,'HEII',1X,'HEIII',3X,'V')

!c 260        FORMAT(1X,1PE10.3,0P2F10.7,1P2E14.7,0P1F10.7,
!     *           1PE10.2,0PF7.3,1PE9.2,0P3F6.3,2F6.3,1PE9.2,
!     *           0P2F12.8)


       radius_linear = dexp(ln10*log10_radius)
       depth_fraction = (surface_radius_linear-radius_linear)/surface_radius_linear
       local_log10_gravity = cgl+log10_star_mass-2.0d0*log10_radius
       local_gravity_linear = dexp(ln10*local_log10_gravity)
         chi_rho = 1.0d0/pulse_diag%qqdp
         chi_t = -chi_rho*pulse_diag%qqdt
         specific_heat_cv = pulse_diag%qqcp - exp(ln10*(run_diag%current_log10_pressure-run_diag%current_log10_density- &
              run_diag%current_log10_temperature))*chi_t**2/chi_rho
         gamma1 = chi_rho*pulse_diag%qqcp/specific_heat_cv
       if(.not.lstch) write(istor,260)local_gravity_linear,run_diag%current_log10_pressure, &
            run_diag%current_log10_temperature,depth_fraction,run_diag%current_log10_mass, &
              run_diag%current_log10_density,run_diag%current_opacity,run_diag%current_beta,(run_diag%current_gradients(i),i=1,3), &
                   (run_diag%current_ion_fraction(i),i=1,3),run_diag%current_velocity, &
              gamma1,pulse_diag%qqdp
      endif
! STORE STARTING VALUES OF THE INTEGRATION
! 07/02 INITIALIZE NUMBER OF STORED ENVELOPE POINTS TO 1
      cz_in_envelope = .false.
      taucz_env_accum = 0.0d0
      env_struct%env_log10_density(1) = run_diag%current_log10_density
      env_struct%env_log10_pressure(1) = run_diag%current_log10_pressure
      env_struct%env_log10_radius(1) = run_diag%current_log10_radius
      env_struct%env_log10_mass(1) = run_diag%current_log10_mass
      env_struct%env_log10_temperature(1) = run_diag%current_log10_temperature
      env_struct%env_hydrogen_fraction(1) = hydrogen_fraction
      env_struct%env_metal_fraction(1) = metal_fraction
      env_struct%env_convective_flag(1) = run_diag%current_velocity.gt.0.0d0
! JVS 03/28
      env_struct%env_gradients(1,1) = run_diag%current_gradients(1)
      env_struct%env_gradients(2,1) = run_diag%current_gradients(2)
      env_struct%env_gradients(3,1) = run_diag%current_gradients(3)
      env_struct%env_convective_velocity(1) = run_diag%current_velocity
      env_struct%env_beta(1) = run_diag%current_beta
! JVS end
! JVS 08/25
! Always save these
      chi_rho = 1.0d0/pulse_diag%qqdp
      chi_t = -chi_rho*pulse_diag%qqdt
      specific_heat_cv = pulse_diag%qqcp - exp(ln10*(run_diag%current_log10_pressure-run_diag%current_log10_density- &
           run_diag%current_log10_temperature))*chi_t**2/chi_rho
      env_struct%env_gamma1(1) = chi_rho*pulse_diag%qqcp/specific_heat_cv
      env_struct%env_specific_heat_cp(1) = pulse_diag%qqcp
      env_struct%env_ion_fraction(1,1) = run_diag%current_ion_fraction(1)
      env_struct%env_ion_fraction(2,1) = run_diag%current_ion_fraction(2)
      env_struct%env_ion_fraction(3,1) = run_diag%current_ion_fraction(3)
      env_struct%env_opacity(1) = run_diag%current_opacity
      env_struct%env_luminosity(1) = luminosity_linear
      env_struct%env_dlnrho_dlnt(1) = pulse_diag%qqdt
! JVS 10/10
      unused_chdelj = run_diag%current_gradients(2)
      unused_chdeld = pulse_diag%qdela
      if(env_struct%env_convective_flag(1))cz_in_envelope = .true.
      env_struct%num_env_points = 1
      do 220 step_index = 1,maxstp
       x_start = indep_var
       do 210 i = 1,num_eqs
          y_start(i) = y(i)
          y_scale(i) = dabs(y(i)) + dabs(h_step*dydx(i))+tiny
 210     continue
       swap_temp = y(1) + h_step*dydx(1)
       if(env_comp%senv - y(1).gt.0.0d0 .or. env_comp%senv - swap_temp.gt.0.0d0) then
!  IF THE INTEGRATION HAS OVERSHOT THE FITTING POINT, OR THE NEXT
!  STEP WILL DO SO,LIMIT STEP SIZE OR INTEGRATE BACKWARDS TO THE
!  CORRECT MASS.
          h_step = (env_comp%senv - y(1))/dydx(1)
       endif
!  ENSURE THAT STEP DOESN'T EXCEED MAXIMUM STEP SIZE
       if(h_step.lt.0.0d0) then
          if(h_step.lt.-h_max) h_step = -h_max
       else
          if(h_step.gt.h_max)  h_step = h_max
       endif
!  PLIM IS AN ESTIMATE OF THE ENDING PRESSURE FOR THE INTEGRATION,
!  BASED ON THE PRESSURE OF THE LAST MODEL POINT.
!  THE FIRST TIME THE INTEGRATOR TRIES TO PASS IT,LIMIT THE STEP
!  IN PRESSURE.
       if(pressure_limit_test_flag) then
          if(indep_var + h_step.gt.log10_pressure_limit) then
             h_step = log10_pressure_limit - indep_var
             pressure_limit_test_flag = .false.
          endif
       endif
!  INTEGRATE THE EQUATIONS FROM X0 TO X0 + H
!  H IS THE ATTEMPTED STEP,HDID IS THE ONE PERFORMED, AND HNEXT IS THE
!  PREDICTED NEXT STEP.
       call bsstep(y,dydx,num_eqs,indep_var,h_step,tolerance,y_scale,h_did,h_next,qenv, &
              luminosity_linear,pressure_rotation_factor,temperature_rotation_factor, &
              log10_gravity,in_atmosphere,want_derivatives,conductive_opacity_flag, &
              print_flag,log10_radius,log10_teff,hydrogen_fraction,metal_fraction, &
              env_call_count,saha_state,step_err)
       do 255 i = 1,3
          err_sum(i) = err_sum(i) + step_err(i)
 255     continue
!  FIND DY/DX AT THE START OF THE NEXT STEP.
       call qenv(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
              temperature_rotation_factor,log10_gravity,in_atmosphere, &
              want_derivatives,conductive_opacity_flag,print_flag,log10_radius, &
              log10_teff,hydrogen_fraction,metal_fraction,env_call_count,saha_state)
! DBG PULSE
         if (pulse_print_flag.and.print_flag) then
          pulse_diag%qqed = 0.0d0
          pulse_energy_sum = 0.0d0
          pulse_diag%qqet = 0.0d0
          electron_pressure = gas_constant * pulse_diag%qt * pulse_diag%qd * pulse_diag%qemu
          if(pulsation_file_version.eq.1) then
               write(opal_envelope_unit,5001)log10_radius,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl, &
                    pulse_diag%qpl, pulse_energy_sum,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed, &
                    pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_diag%qdel,pulse_diag%qdela, &
                    pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
            else if (pulsation_file_version.eq.2) then
               write(opal_envelope_unit,6001)log10_radius,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl, &
                    pulse_diag%qpl, pulse_energy_sum,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed, &
                    pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_diag%qdel,pulse_diag%qdela, &
                    pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
            else if (pulsation_file_version.eq.3) then
! DBG 7/95 Appended mixing length info at end of first three lines
               write(opal_envelope_unit,6003)log10_radius,pulse_diag%qfs,luminosity_linear,pulse_diag%qtl,pulse_diag%qdl,rot_diff%alfmlt, &
                    pulse_diag%qpl, pulse_energy_sum,pulse_diag%qo,pulse_diag%qqdp,pulse_diag%qqed,rot_diff%phmlt, &
                    pulse_diag%qqet,pulse_diag%qqod,pulse_diag%qqot,pulse_diag%qdel,pulse_diag%qdela,rot_diff%cmxmlt, &
                    pulse_diag%qqcp,pulse_diag%qrmu,pulse_diag%qqdt,electron_pressure
            end if
         end if
! DBG END
! G Somers 11/14 ADDED I/O FLAG AND CHANGED WRITE OUTS TO .STORE.
       if(print_flag.and.lstenv) then
          radius_linear = dexp(ln10*log10_radius)
          depth_fraction = (surface_radius_linear-radius_linear)/surface_radius_linear
          local_log10_gravity = cgl+log10_star_mass-2.0d0*log10_radius
          local_gravity_linear = dexp(ln10*local_log10_gravity)
            chi_rho = 1.0d0/pulse_diag%qqdp
            chi_t = -chi_rho*pulse_diag%qqdt
            specific_heat_cv = pulse_diag%qqcp - exp(ln10*(run_diag%current_log10_pressure-run_diag%current_log10_density- &
                 run_diag%current_log10_temperature))*chi_t**2/chi_rho
            gamma1 = chi_rho*pulse_diag%qqcp/specific_heat_cv
          if(.not.lstch) write(istor,260)local_gravity_linear,run_diag%current_log10_pressure, &
               run_diag%current_log10_temperature,depth_fraction,run_diag%current_log10_mass, &
                 run_diag%current_log10_density,run_diag%current_opacity,run_diag%current_beta,(run_diag%current_gradients(i),i=1,3), &
                      (run_diag%current_ion_fraction(i),i=1,3),run_diag%current_velocity, &
                 gamma1,pulse_diag%qqdp
 260        format(1X,1PE10.3,0P2F10.7,1P2E14.7,0P1F10.7, &
                 1PE10.2,0PF7.3,1PE9.2,0P3F6.3,2F6.3,1PE9.2, &
                 0P2F12.8)
       endif
       if(h_did.eq.h_step) then
          num_ok = num_ok + 1
       else
          num_bad = num_bad + 1
       endif
!  CHECK IF INTEGRATION COMPLETE
       mass_diff_remaining = env_comp%senv - y(1)
! 07/02 STORE ENVELOPE TERMS IF THE INTEGRATION
! HAS NOT OVERSHOT THE FITTING POINT.
         if(mass_diff_remaining.le.step_tolerance)then
            env_struct%num_env_points = env_struct%num_env_points + 1
            env_struct%env_log10_density(env_struct%num_env_points) = run_diag%current_log10_density
            env_struct%env_log10_pressure(env_struct%num_env_points) = run_diag%current_log10_pressure
            env_struct%env_log10_radius(env_struct%num_env_points) = run_diag%current_log10_radius
            env_struct%env_log10_mass(env_struct%num_env_points) = run_diag%current_log10_mass
            env_struct%env_log10_temperature(env_struct%num_env_points) = run_diag%current_log10_temperature
            env_struct%env_hydrogen_fraction(env_struct%num_env_points) = hydrogen_fraction
            env_struct%env_metal_fraction(env_struct%num_env_points) = metal_fraction
            env_struct%env_convective_flag(env_struct%num_env_points) = run_diag%current_velocity.gt.0.0d0
! JVS 08/13 ADD RUN FOR CZ CALCULATION
            env_struct%env_gradients(1,env_struct%num_env_points) = run_diag%current_gradients(1)
            env_struct%env_gradients(2,env_struct%num_env_points) = run_diag%current_gradients(2)
            env_struct%env_gradients(3,env_struct%num_env_points) = run_diag%current_gradients(3)
            env_struct%env_convective_velocity(env_struct%num_env_points) = run_diag%current_velocity
            env_struct%env_beta(env_struct%num_env_points) = run_diag%current_beta
! JVS 08/25
! Always save these
            chi_rho = 1.0d0/pulse_diag%qqdp
            chi_t = -chi_rho*pulse_diag%qqdt
            specific_heat_cv = pulse_diag%qqcp - exp(ln10*(run_diag%current_log10_pressure-run_diag%current_log10_density- &
                 run_diag%current_log10_temperature))*chi_t**2/chi_rho
            env_struct%env_gamma1(env_struct%num_env_points) = chi_rho*pulse_diag%qqcp/specific_heat_cv
            env_struct%env_specific_heat_cp(env_struct%num_env_points) = pulse_diag%qqcp
            env_struct%env_ion_fraction(1,env_struct%num_env_points) = run_diag%current_ion_fraction(1)
            env_struct%env_ion_fraction(2,env_struct%num_env_points) = run_diag%current_ion_fraction(2)
            env_struct%env_ion_fraction(3,env_struct%num_env_points) = run_diag%current_ion_fraction(3)
            env_struct%env_opacity(env_struct%num_env_points) = run_diag%current_opacity
            env_struct%env_luminosity(env_struct%num_env_points) = luminosity_linear
            env_struct%env_dlnrho_dlnt(env_struct%num_env_points) = pulse_diag%qqdt

            if(.not.cz_in_envelope)then
               if(env_struct%env_convective_flag(env_struct%num_env_points))then
                  cz_in_envelope = .true.
               endif
            else if(run_diag%current_velocity.gt.0.0d0)then
               delta_radius_cz = 10.0d0**env_struct%env_log10_radius(env_struct%num_env_points-1) - &
                    10.0d0**env_struct%env_log10_radius(env_struct%num_env_points)
               taucz_env_accum = taucz_env_accum + delta_radius_cz/run_diag%current_velocity
            endif
         endif
       if(dabs(mass_diff_remaining).le.step_tolerance)then
!            WRITE(*,*)TAUCZENV/CSECYR
          if(save_boundary_flag) then
             interp_weight = mass_diff_remaining/(y_start(1)-y(1))
             vtx_logp(vertex_index) = indep_var + interp_weight*(x_start - indep_var)
             vtx_logr(vertex_index) = y(3) + interp_weight*(y_start(3) - y(3))
             vtx_logt(vertex_index) = y(2) + interp_weight*(y_start(2) - y(2))
             if(print_flag)write(short_file_unit,230)vtx_logp(vertex_index),vtx_logt(vertex_index),vtx_logr(vertex_index),env_comp%senv
          endif
          goto 300
       else if(.not.store_flag_set) then
          if(y(2).ge.tenv .and. save_boundary_flag) then
             store_flag_set = .true.
             stored_vertex_index = vertex_index
             stored_envelope_state(1) = indep_var
             stored_envelope_state(2) = y(2)
             stored_envelope_state(3) = y(3)
             stored_envelope_state(4) = y(1)
          endif
       endif
       if(h_next.lt.h_min) h_next = h_min
       h_step = h_next
 220  continue
! INTEGRATION HAS FAILED TO FINISH IN MAXSTP STEPS;
! PRINT NASTY MESSAGE AND QUIT.
      write(iowr,911)
      write(short_file_unit,911)
 911  format(5X,'ENVELOPE INTEGRATION FAILED AFTER MAXSTP TRIES.',1X, &
           'I QUIT')
      stop
 300  continue
! 07/02 NOW INVERT THE ENVELOPE VECTOR.
      if(env_comp%senv.lt.-1.0d-12)then
         do i = 1,env_struct%num_env_points
            inversion_index1 = i
            inversion_index2 = env_struct%num_env_points - i + 1
            if(inversion_index1.ge.inversion_index2)goto 310
            swap_temp = env_struct%env_log10_density(inversion_index1)
            env_struct%env_log10_density(inversion_index1) = env_struct%env_log10_density(inversion_index2)
            env_struct%env_log10_density(inversion_index2) = swap_temp
            swap_temp = env_struct%env_log10_pressure(inversion_index1)
            env_struct%env_log10_pressure(inversion_index1) = env_struct%env_log10_pressure(inversion_index2)
            env_struct%env_log10_pressure(inversion_index2) = swap_temp
            swap_temp = env_struct%env_log10_radius(inversion_index1)
            env_struct%env_log10_radius(inversion_index1) = env_struct%env_log10_radius(inversion_index2)
            env_struct%env_log10_radius(inversion_index2) = swap_temp
            swap_temp = env_struct%env_log10_mass(inversion_index1)
            env_struct%env_log10_mass(inversion_index1) = env_struct%env_log10_mass(inversion_index2)
            env_struct%env_log10_mass(inversion_index2) = swap_temp
            swap_temp = env_struct%env_log10_temperature(inversion_index1)
            env_struct%env_log10_temperature(inversion_index1) = env_struct%env_log10_temperature(inversion_index2)
            env_struct%env_log10_temperature(inversion_index2) = swap_temp
            swap_temp = env_struct%env_hydrogen_fraction(inversion_index1)
            env_struct%env_hydrogen_fraction(inversion_index1) = env_struct%env_hydrogen_fraction(inversion_index2)
            env_struct%env_hydrogen_fraction(inversion_index2) = swap_temp
            swap_temp = env_struct%env_metal_fraction(inversion_index1)
            env_struct%env_metal_fraction(inversion_index1) = env_struct%env_metal_fraction(inversion_index2)
            env_struct%env_metal_fraction(inversion_index2) = swap_temp
            swap_temp_logical = env_struct%env_convective_flag(inversion_index1)
            env_struct%env_convective_flag(inversion_index1) = env_struct%env_convective_flag(inversion_index2)
            env_struct%env_convective_flag(inversion_index2) = swap_temp_logical
!  08/25 JVS
            swap_temp = env_struct%env_opacity(inversion_index1)
            env_struct%env_opacity(inversion_index1) = env_struct%env_opacity(inversion_index2)
            env_struct%env_opacity(inversion_index2) = swap_temp
            swap_temp = env_struct%env_luminosity(inversion_index1)
            env_struct%env_luminosity(inversion_index1) = env_struct%env_luminosity(inversion_index2)
            swap_temp = env_struct%env_dlnrho_dlnt(inversion_index1)
            env_struct%env_dlnrho_dlnt(inversion_index1) = env_struct%env_dlnrho_dlnt(inversion_index2)
            env_struct%env_dlnrho_dlnt(inversion_index2) = swap_temp

! 08/13 JVS ADDED DEL VECTORS
            swap_temp = env_struct%env_gradients(1,inversion_index1)
            env_struct%env_gradients(1,inversion_index1) = env_struct%env_gradients(1,inversion_index2)
            env_struct%env_gradients(1,inversion_index2) = swap_temp
            swap_temp = env_struct%env_gradients(2,inversion_index1)
            env_struct%env_gradients(2,inversion_index1) = env_struct%env_gradients(2,inversion_index2)
            env_struct%env_gradients(2,inversion_index2) = swap_temp
            swap_temp = env_struct%env_gradients(3,inversion_index1)
            env_struct%env_gradients(3,inversion_index1) = env_struct%env_gradients(3,inversion_index2)
            env_struct%env_gradients(3,inversion_index2) = swap_temp
            swap_temp = env_struct%env_convective_velocity(inversion_index1)
            env_struct%env_convective_velocity(inversion_index1) = env_struct%env_convective_velocity(inversion_index2)
            env_struct%env_convective_velocity(inversion_index2) = swap_temp
            swap_temp = env_struct%env_beta(inversion_index1)
            env_struct%env_beta(inversion_index1) = env_struct%env_beta(inversion_index2)
            env_struct%env_beta(inversion_index2) = swap_temp
!  08/25 JVS
            swap_temp = env_struct%env_gamma1(inversion_index1)
            env_struct%env_gamma1(inversion_index1) = env_struct%env_gamma1(inversion_index2)
            env_struct%env_gamma1(inversion_index2) = swap_temp
            swap_temp = env_struct%env_specific_heat_cp(inversion_index1)
            env_struct%env_specific_heat_cp(inversion_index1) = env_struct%env_specific_heat_cp(inversion_index2)
            env_struct%env_specific_heat_cp(inversion_index2) = swap_temp
            swap_temp = env_struct%env_ion_fraction(1,inversion_index1)
            env_struct%env_ion_fraction(1,inversion_index1) = env_struct%env_ion_fraction(1,inversion_index2)
            env_struct%env_ion_fraction(1,inversion_index2) = swap_temp
            swap_temp = env_struct%env_ion_fraction(2,inversion_index1)
            env_struct%env_ion_fraction(2,inversion_index1) = env_struct%env_ion_fraction(2,inversion_index2)
            env_struct%env_ion_fraction(2,inversion_index2) = swap_temp
            swap_temp = env_struct%env_ion_fraction(3,inversion_index1)
            env_struct%env_ion_fraction(3,inversion_index1) = env_struct%env_ion_fraction(3,inversion_index2)
            env_struct%env_ion_fraction(3,inversion_index2) = swap_temp

         end do
 310     continue
! JVS 07/12 Save the last envelope point pressure
!      PPHOT = ENVP(NUMENV) ! G Somers 3/17, MOVED PPHOT DEF HIGHER UP
! END JVS
      endif
! DBG PULSE WRITE END OF DATA INDICATOR
      if (pulse_print_flag.and.print_flag) then
!         XYZ = 99.99D0
         if(pulsation_file_version.eq.1) then
            write(opal_envelope_unit, 5001) (xyz(i),i=1,19)
         else if (pulsation_file_version.eq.2 .or. pulsation_file_version.eq.3)then
            write(opal_envelope_unit, 6003) (xyz(i),i=1,22)
         end if
      endif

      write(short_file_unit,215)num_ok,num_bad,mass_diff_remaining,y(1),(err_sum(jj),jj=1,3)
 215  format(1X,'ENVELOPE INTEGRATION COMPLETE',1X, &
           'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
           I5/5X,'SENV-LAST M=',1PE22.13,'  LAST M=',E22.13/ &
           5X,'MAX RELATIVE ERRORS:M ',1PE14.5,'  T ',E14.5, &
           '  R ',E14.5)

! JVS 08/13
! IF THE CZ IS IN THE ENVELOPE (.I.E. BEYOND THE FITTING POINT) TRACK ITS
! LOCATION IN MASS AND RADIUS FOR USE WITH AM LOSS ROUTINES
!
! G Somers 3/17, SKIP TAUCAL CALL IF USING NEW TAUCZ ROUTINES
        if(use_new_turnover_timescale)goto 555
! G Somers END
      cz_start_index = 0
      do i=1,env_struct%num_env_points
            if (cz_start_index .eq. 0 .and. env_struct%env_convective_flag(i) ) cz_start_index = i
      end do

      if (cz_start_index .gt. 1) then
            ! Calculate the the location of the base sof the surface CZ in radius
            dd2 = env_struct%env_gradients(1,cz_start_index-1)-env_struct%env_gradients(2,cz_start_index-1)
            dd1 = env_struct%env_gradients(1,cz_start_index)-env_struct%env_gradients(2,cz_start_index)
            interp_fraction = dd2/(dd2-dd1)

            log10_cz_radius = env_struct%env_log10_radius(cz_start_index-1)+interp_fraction* &
                 (env_struct%env_log10_radius(cz_start_index)-env_struct%env_log10_radius(cz_start_index-1))-log10_solar_radius
            run_diag%convection_zone_radius_placeholder = exp(ln10*log10_cz_radius)

            do i=1,env_struct%num_env_points
                  taucal_shell_mass(i) = dexp(ln10*(env_struct%env_log10_mass(i)+log10_star_mass))

                  if (i.eq.env_struct%num_env_points) then
                        taucal_delta_mass(i) = 0.0
                  else
                        taucal_delta_mass(i)=dexp(ln10*env_struct%env_log10_mass(i+1))-dexp(ln10*env_struct%env_log10_mass(i))
                  endif

                  taucal_radiative_gradient(i) = env_struct%env_gradients(1,i)
                  taucal_adiabatic_gradient(i) = env_struct%env_gradients(2,i)
            end do

!             CALL TAUCAL(ENVX,ENVS2,ENVS1,LCENV,ENVR,ENVP,ENVD,ENVG,NUMENV,  ! KC 2025-05-31
            call taucal(taucal_delta_mass,taucal_shell_mass,env_struct%env_convective_flag,env_struct%env_log10_radius, &
                  env_struct%env_log10_pressure,env_struct%env_log10_density,taucal_local_gravity,env_struct%num_env_points, &
                  env_struct%env_convective_velocity, taucal_radiative_gradient,taucal_adiabatic_gradient)
      endif

555      continue
      return
end subroutine envint
