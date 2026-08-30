!----------------------------------------------------------------------
! envint_lib
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass, atm split -- ROADMAP.md): the
! envelope/atmosphere INTEGRATION, moved out of atm_lib. atm_get is
! the star solver's outer-boundary integrator: it integrates the
! atmosphere and envelope down to the fit point, driving the qatm/
! qenv integrands (which call the eos/kap facades) through numerics'
! bsstep, consuming the pure atm surface lookups (surfp/kcsurfp/
! alsurfp -- still atm-domain physics) and writing model state
! (star%pulse, star%run, env_comp bookkeeping, the envelope/
! atmosphere structure arrays). In MESA terms this is star-solver
! code (create_atm/env integration), not the atm physics module --
! which is why it now lives in core/ with qatm/qenv/surfbc and the
! turnover diagnostics. atm/ keeps what MESA's atm has: tables,
! their init, and pure surface lookups.
module envint_lib
      implicit none
contains

subroutine atm_get(luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, log10_star_mass, &
     vertex_index, print_flag, save_boundary_flag, log10_pressure_limit, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     stored_envelope_state, stored_vertex_index, atm_call_count, &
     env_call_count, saha_state, vtx_logp, vtx_logr, vtx_logt, &
     ierr)

      use eos_lib
      use kap_lib
      use atm_table_lib
      use star_info_lib, only: star, json
      use point_scratch_lib
      use atmstruct_lib
      use envstruct_lib
      use luout_lib
      use phys_const_lib
      use intpar_lib
      use numerics_lib
      use run_log_lib, only: solver_diagnostics
      implicit none
! PARAMETERS NT AND NG FOR TABULATED SURFACE PRESSURES OF KURUCZ.
      integer, parameter :: nt=57, ng=11
! JNT 06/14 ADDED NTC/NGC FOR KTTAU=5
      integer, parameter :: ntc=76, ngc=11
! MHP 8/97 ADDED NTA AND NGA FOR ALLARD ATMOSPHERE TABLES
      integer, parameter :: nta=54, nga=5

! luminosity_linear/pressure_rotation_factor/temperature_rotation_factor/
! log10_gravity/hydrogen_fraction/metal_fraction are
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
! 2026 (.store convergence + retirement): print_flag no longer
! reaches the integrand callbacks (they save unconditionally) and
! envint writes no output files anymore. What remains is init-echo
! verbosity: the fitted envelope-vertex line to the .short log.
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
! 2026 (ROADMAP.md stage 3): OPTIONAL ierr, the transitional form of
! MESA's ierr-not-stop discipline (same contract as kap_lib's and
! eos_lib's facades). Passed: any error that used to stop inside the
! atm internals is returned here instead, ierr /= 0, no stop. Omitted:
! exactly the historical diagnostics, then a stop -- now in the funnel
! at the end of this subroutine rather than at the point of failure.
      integer, intent(out) :: ierr

! DBG CHANGED MAXSTEP FROM 200 TO 2000 TO GIVE ATMOSPHERE INTEGRATER A CHANCE.
      integer, parameter :: maxstp = 2000
      logical :: tabulated_bc
      double precision, parameter :: tiny = 1.0d-30
      external atmosphere_derivs, envelope_derivs
      double precision :: harvard_t_tau
      external harvard_t_tau

      double precision :: ion_fraction(3)


! G Somers END


! MHP 1/01 CHANGED END OF FILE INDICATOR IN ATMOSPHERE/ENVELOPE FILES TO
! VECTOR FROM SCALAR
      double precision :: xyz(22)
      data xyz/22*99.99d0/
! --- locals ---
      logical :: want_derivatives, in_atmosphere, conductive_opacity_flag
      logical :: allard_lookup_failed
      integer :: jj, jerr
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
! 2026 named-index call buffers for eos_get/kap_get. Unlike the
! other migrated sites, the relay scalars above STAY: they are host-
! associated across the contained routines (and the atm_retry loop),
! so the call is wrapped with a symmetric prepack/unpack instead.
      double precision :: eos_res(num_eos_results), kap_res(num_kap_results)
      double precision :: indep_var
      double precision :: chi_rho, chi_t, specific_heat_cv, gamma1
      double precision :: prev_tau, prev_opacity, &
           prev_density, delta_tau_step, &
           pulse_radiative_gradient, pulse_gradient, opacity_now, &
           density_now, tau_now
      integer :: num_eqs, num_ok, num_bad, step_index, i
      double precision :: tolerance, x_limit, h_max, h_min, h_step, &
           step_tolerance, h_did, h_next
      double precision :: y(3), dydx(3), y_scale(3), y_start(3), &
           err_sum(3), step_err(3)
      double precision :: initial_mass_coord
      logical :: store_flag_set, pressure_limit_test_flag
      logical :: cz_in_envelope
      double precision :: mass_diff_remaining, interp_weight
      integer :: inversion_index1, inversion_index2
      double precision :: swap_temp
      logical :: swap_temp_logical
      double precision :: unused_chdelj, unused_chdeld
      double precision :: x_start, taucz_env_accum, delta_radius_cz

      ierr = 0
      jerr = 0

      call prepare_surface_boundary
      if (ierr /= 0) return
      call integrate_atmosphere
      if (ierr /= 0) return
! ENVELOPE INTEGRATION
! HERE P IS THE INDEPENDENT VARIABLE AND M,R,AND T ARE
! DEPENDENT VARIABLES.  INTEGRATE FROM TAU = 2/3 TO THE LAST
! MASS POINT IN THE MODEL.


                 !  integration and come here
! G Somers 3/17, IF INTERESTED ONLY IN PPHOT, BREAK HERE.
      star%pphot = atm_table%atm_log10_pressure
      if (.not.star%job%calc_envelope_flag) then
         continue
         return
      end if

      call integrate_envelope
      if (ierr /= 0) return
! 2026: the legacy taucal path (track_envelope_cz, gated by
! LNEWTCZ=.false.) is retired -- the turnover timescale is computed
! by core/turnover/turnover_timescale.f90 from the assembled
! interior+envelope structure, never as a side effect of an
! envelope integration.
      return

! error funnel: reached only when a callee (or one of the two
! integration-failure checks above) reported jerr /= 0. With ierr
! present the caller takes responsibility; without it, preserve the
! historical stop (the diagnostic already printed at the point of
! failure).

contains

! ---------------------------------------------------------------
! Surface boundary condition: pulse-derivative flags, then the
! tabulated surface pressure lookup (Kurucz / Kurucz-Castelli /
! Allard, falling back to gray when the Allard lookup fails).
! Sets tabulated_bc; ierr on a failed table interpolation.
subroutine prepare_surface_boundary
! DBG PULSE TURN ON DERIVATIVE CALCULATOR
! 2026 retire-legacy: with the .pmod/.penv/.patm pulse trio gone,
! the pulse derivative mode (want_derivatives + in-envelope eos
! flags) is never enabled here; the always-taken branch remains.
! (.store convergence: lpumod itself is deleted -- the derivs
! routines now save their output scratch unconditionally.)
      want_derivatives = .false.
      in_atmosphere = .true.
      conductive_opacity_flag = .false.

! JVS 10/07/13 Always calculate derivatives
      want_derivatives = .true.

! 2026 (.store convergence): the "ATMOSPHERE BEGIN" header and the
! mid-integration atmosphere/envelope text tables that used to
! stream into the legacy .store from here are retired -- the same
! structures are materialized in env_struct/atmo_struct and written
! by the profile/pulse/.store writers. envint no longer writes to
! istor at all.

! GET PRESSURE AT T=Teff BY INTERPOLATION IN TABLE ATMPL.
      tabulated_bc = .false.
      if (star%job%atm_choice .eq. 3) then
! KURUCZ ATMOSPHERES
         call surfp(log10_teff,log10_gravity,.false.,jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         tabulated_bc = .true.
! JNT 06/14
! GET PRESSURE AT T=Teff BY INTERPOLATION IN TABLE ATMPLC.
      else if (star%job%atm_choice .eq. 5) then
! KURUCZ ATMOSPHERES
         call kcsurfp(log10_teff,log10_gravity,.false.,jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         tabulated_bc = .true.
! We have Kurucz atmosphere boundary conditions
      else if (star%job%atm_choice .eq. 4) then
! ALLARD & HAUSCHILDT ATMOSPHERES
         call alsurfp(log10_teff,log10_gravity,.false.,allard_lookup_failed,jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
! Changed to Allard atmosphere code
         if(allard_lookup_failed) then
            star%job%atm_choice=0
            star%use_ttau_relation = .true.
! Set to gray atmosphere (KTTAU=0), as
! TeffL is above Allard max, or GL is out of range.
            write(*,'(1x,a)') 'note: Allard table lookup out of range;' &
                 // ' switching to a gray atmosphere'
            write(run_log_unit,'(1x,a)') 'note: Allard table lookup' &
                 // ' out of range; switching to a gray atmosphere'
         else
            tabulated_bc = .true.
         endif
! We have Allard atmosphere boundary conditions
      endif
end subroutine prepare_surface_boundary

! ---------------------------------------------------------------
! Gray-atmosphere start (guess T at tau ~ 0, matching P, rho and
! starting tau) and the dP/dtau integration down to the
! photosphere (tau = 2/3 Eddington, 0.312 Krishna-Swamy), storing
! atmo_struct on the way; the atm_retry loop restarts with a
! reduced step on failure, ierr after maxstp steps. Skipped
! entirely when a tabulated boundary supplied the pressure.
subroutine integrate_atmosphere
      use math_lib
! 2026 (.store convergence): the count is zeroed on every entry so a
! skipped integration (tabulated surface boundary) leaves an
! accurate "no atmosphere points" state rather than the previous
! integration's stale count. atmo_struct is filled unconditionally
! below -- it is the atmosphere the stitched model materializes.
      atmo_struct%num_atm_points = 0
      if (.not. tabulated_bc) then
! Start gray atmosphere bounary conditions
! GUESS THE TEMPERATURE FOR AN OPTICAL DEPTH NEAR ZERO.
      err_sum(1) = 0.0d0
      if(star%job%atm_choice .eq. 0) then
            log10_temperature = log10_teff - 0.031235d0 + 0.25d0*log10(cc23)
      else if (star%job%atm_choice .eq. 1) then
            log10_temperature = log10_teff - 0.031235d0 + 0.25d0*log10(0.550d0)
      else if (star%job%atm_choice .eq. 2) then
            log10_temperature = log10_teff + harvard_t_tau(cc23) - star%atm_hras
      end if
!                 For kttau = 0,1,or 2, very occasionally the integration
!                 fails because the starting point (X0) is past the end
!                 point (XLIM). When this happens, we divide the effective
!                 starting density (atmd0) by 10 and retry.
      atm_density_guess = star%ctrl%atm_step_initial
      atm_retry: do
! Return point if X0 > XLIM (was label 1998)
      temperature = exp(ln10*log10_temperature)
! FIND THE PRESSURE CORRESPONDING TO THIS T AND THE DENSITY CHOSEN
! FOR THE START OF THE ATMOSPHERE INTEGRATION.

!      PL = log10((CGAS*ATMD0 + CA3*T**3)*T)
      log10_pressure = log10((gas_constant*atm_density_guess + &
           radiation_constant_over_3*temperature**3)*temperature)

! NOW FIND THE OPTICAL DEPTH(X0) WHERE THE ATMOSPHERE INTEGRATION BEGINS.
! prepack every result slot from its host scalar, call, unpack every
! one back: identity for slots eqstat leaves alone, and the host
! scalars keep carrying the historical inout guesses across calls.
      eos_res(i_temperature) = temperature
      eos_res(i_pressure) = pressure
      eos_res(i_log10_density) = log10_density
      eos_res(i_density) = density
      eos_res(i_beta) = beta
      eos_res(i_beta_inverse) = beta_inverse
      eos_res(i_beta14) = beta14
      eos_res(i_fxion:i_fxion+2) = ion_fraction
      eos_res(i_gas_constant) = specific_gas_constant
      eos_res(i_mu_ion_inv) = ion_mean_weight_inverse
      eos_res(i_mu_e_inv) = electron_mean_weight_inverse
      eos_res(i_eta) = electron_degeneracy_parameter
      eos_res(i_dlnrho_dlnt) = dlnrho_dlnt
      eos_res(i_dlnrho_dlnp) = dlnrho_dlnp
      eos_res(i_cp) = specific_heat_cp
      eos_res(i_grada) = adiabatic_gradient
      eos_res(i_dlnrho_dlnt_dt) = dlnrho_dlnt_dt
      eos_res(i_dlnrho_dlnp_dt) = dlnrho_dlnp_dt
      eos_res(i_grada_dt) = adiabatic_gradient_dt
      eos_res(i_grada_dp) = adiabatic_gradient_dp
      eos_res(i_cp_dt) = specific_heat_cp_dt
      eos_res(i_cp_dp) = specific_heat_cp_dp
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_state, ierr=jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
      temperature = eos_res(i_temperature)
      pressure = eos_res(i_pressure)
      log10_density = eos_res(i_log10_density)
      density = eos_res(i_density)
      beta = eos_res(i_beta)
      beta_inverse = eos_res(i_beta_inverse)
      beta14 = eos_res(i_beta14)
      ion_fraction = eos_res(i_fxion:i_fxion+2)
      specific_gas_constant = eos_res(i_gas_constant)
      ion_mean_weight_inverse = eos_res(i_mu_ion_inv)
      electron_mean_weight_inverse = eos_res(i_mu_e_inv)
      electron_degeneracy_parameter = eos_res(i_eta)
      dlnrho_dlnt = eos_res(i_dlnrho_dlnt)
      dlnrho_dlnp = eos_res(i_dlnrho_dlnp)
      specific_heat_cp = eos_res(i_cp)
      adiabatic_gradient = eos_res(i_grada)
      dlnrho_dlnt_dt = eos_res(i_dlnrho_dlnt_dt)
      dlnrho_dlnp_dt = eos_res(i_dlnrho_dlnp_dt)
      adiabatic_gradient_dt = eos_res(i_grada_dt)
      adiabatic_gradient_dp = eos_res(i_grada_dp)
      specific_heat_cp_dt = eos_res(i_cp_dt)
      specific_heat_cp_dp = eos_res(i_cp_dp)
! DBG 12/95 GET OPACITY (log10_density just unpacked = eqstat's
! returned density, the historical inout dataflow)
      call kap_get(log10_density, log10_temperature, hydrogen_fraction, &
           metal_fraction, kap_res, ion_fraction, ierr=jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
      opacity = kap_res(i_kap)
      log10_opacity = kap_res(i_log10_kap)
      dlnkap_dlnrho = kap_res(i_dlnkap_dlnrho)
      dlnkap_dlnt = kap_res(i_dlnkap_dlnt)
      indep_var = log10_pressure - log10_gravity + log10(opacity)
      y(1) = log10_pressure
      dydx(1) = exp(ln10*(log10_gravity+indep_var-log10_opacity-log10_pressure))
! DBG PULSE INITIAL POINT FOR PULSATION



! INTEGRATE DP/DTAU FROM THIS STARTING TAU TO TAU = 2/3.
! SET NUMERICAL PARAMETERS UP.
      num_eqs = 1
      tolerance = star%ctrl%atm_error_tol
      if (star%job%atm_choice .eq. 1) then
! KRISHNA-SWAMY T TAU HAS DIFFERENT ZERO THAN EDDINGTON T TAU
! TAU = 0.312156330 AT TEFF.
            x_limit = -0.505627854d0
      else
! TAU = 2/3 AT TEFF.
            x_limit = -0.176091259d0
      end if

      if  (indep_var .gt. x_limit) then ! Check that starting point is before endpoint
         write(run_log_unit,*)"ENVINT: X0>XLIM, X0,XLIM:",indep_var,x_limit
       write(run_log_unit,*)"ENVINT: get new X0 by dividing ATMD0 by 10"
       atm_density_guess = atm_density_guess / 10d0   ! If not before, divide starting density by 10
       cycle atm_retry        ! and retry.
      endif

      h_max = star%job%atm_step_max
      h_min = star%job%atm_step_min
      h_step = star%job%atm_step_begin
      step_tolerance = tolerance_fraction
      num_ok = 0
      num_bad = 0

! INTEGRATION LOOP.
      do step_index = 1,maxstp
! YSCAL IS THE ARRAY THAT THE NUMERICAL ERRORS ARE SCALED AGAINST.
       do i = 1,num_eqs
          y_scale(i) = dabs(y(i)) + dabs(h_step*dydx(i))+tiny
       end do
! ENSURE THAT STEP DOESN'T EXCEED MAXIMUM STEP SIZE OR GO PAST
! THE LIMIT OF THE INTEGRATION.
       if((indep_var-x_limit)*(indep_var+h_step-x_limit).lt.0.0d0) h_step = x_limit - indep_var
       if(h_step.gt.h_max) h_step = h_max
! INTEGRATE THE ATMOSPHERE FROM TAU TO TAU + H
! H IS THE ATTEMPTED STEP,HDID IS THE ONE PERFORMED, AND HNEXT IS THE
! PREDICTED NEXT STEP.
       call bsstep(y,dydx,num_eqs,indep_var,h_step,tolerance,y_scale,h_did, &
            h_next,atmosphere_derivs, luminosity_linear,pressure_rotation_factor, &
            temperature_rotation_factor,log10_gravity,in_atmosphere, &
            want_derivatives,conductive_opacity_flag,log10_radius, &
            log10_teff,hydrogen_fraction,metal_fraction,atm_call_count, &
            saha_state,step_err,ierr=jerr)
! 2026 numerics-gate opt-in: a bsstep failure (the historical
! "solution diverged" stop, its diagnostic already printed by the
! gate) becomes the graceful numerics-termination code when the
! caller opted in via ierr; without ierr the historical stop stands.
       if (jerr /= 0) then
          ierr = numerics_termination
          return
       end if
! FIND DP/DTAU AT THE START OF THE NEXT STEP.
       err_sum(1) = err_sum(1) + step_err(1)
       call atmosphere_derivs(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
            temperature_rotation_factor,log10_gravity,in_atmosphere, &
            want_derivatives,conductive_opacity_flag,log10_radius, &
            log10_teff,hydrogen_fraction,metal_fraction,atm_call_count,saha_state, jerr)
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
! 2026 (.store convergence): the structure save below was gated on
! the print flag (so only the legacy .store stitch, which called
! atm_get with printing on, ever populated atmo_struct), and the
! depth/gradient save further down was gated on lstch. Both now run
! unconditionally: every atmosphere integration materializes the
! full atmo_struct, so the stitched model (and any writer reading
! it) always sees the atmosphere. Values are identical to what the
! print-on path stored.
       beta = 1.0d0 - radiation_constant_over_3*exp(ln10*(4.0d0*atm_table%atm_log10_temperature-atm_table%atm_log10_pressure))
       chi_rho = 1.0d0/pt_scr%qqdp
       chi_t = -chi_rho*pt_scr%qqdt
       specific_heat_cv = pt_scr%qqcp - exp(ln10*(atm_table%atm_log10_pressure-atm_table%atm_log10_density-atm_table%atm_log10_temperature))*chi_t**2/chi_rho
       gamma1 = chi_rho*pt_scr%qqcp/specific_heat_cv
! JvS: SAVE STRUCTURE TO COMMON BLOCK
       atmo_struct%atmo_log10_pressure(step_index) = atm_table%atm_log10_pressure
       atmo_struct%atmo_log10_temperature(step_index) = atm_table%atm_log10_temperature
       atmo_struct%atmo_log10_density(step_index) = atm_table%atm_log10_density
       atmo_struct%atmo_beta(step_index) = beta
       atmo_struct%atmo_gamma1(step_index) = gamma1
       atmo_struct%atmo_dlnrho_dlnt(step_index) = pt_scr%qqdt
       atmo_struct%atmo_ion_fraction(1,step_index) = atm_table%atm_ion_fraction(1)
       atmo_struct%atmo_ion_fraction(2,step_index) = atm_table%atm_ion_fraction(2)
       atmo_struct%atmo_ion_fraction(3,step_index) = atm_table%atm_ion_fraction(3)
       atmo_struct%atmo_opacity(step_index) = atm_table%atm_opacity
       atmo_struct%atmo_specific_heat_cp(step_index) = pt_scr%qqcp
       if(h_did.eq.h_step) then
          num_ok = num_ok + 1
       else
          num_bad = num_bad + 1
       endif
! Geometric step length between successive tau points (J.P. Cox,
! Princ. of Stell. Struc. p590). prev_tau/prev_density/prev_opacity
! carry across steps; their values entering the FIRST step are
! whatever the previous integration left (a historical quirk of the
! .store stitch, preserved -- it only affects the outermost
! atmosphere point's depth).
       opacity_now = pt_scr%qo
       density_now = exp(ln10*pt_scr%qdl)
       tau_now = exp(ln10*atm_table%atm_tau)
       delta_tau_step =  (tau_now - prev_tau)/(((density_now*opacity_now)+(prev_density*prev_opacity))/2)
       prev_opacity = opacity_now
!FROM FIRST LINES OF TPGRAD
       pulse_radiative_gradient = pt_scr%qo*luminosity_linear*exp(ln10*(pt_scr%qpl-log10_star_mass-4.0d0*pt_scr%qtl+star%log10_solar_luminosity-cgl+ &
              cdelrl))*temperature_rotation_factor/pressure_rotation_factor
       if (pulse_radiative_gradient-pt_scr%qdela .le. 1.0d-6) then
         pulse_gradient = pulse_radiative_gradient
       else
         pulse_gradient = pt_scr%qdela
       end if
! JvS SAVE TO COMMON ATMSTRUCT COMMON BLOCK
       atmo_struct%atmo_delta_depth(step_index) = delta_tau_step
       atmo_struct%atmo_gradients(1,step_index) = pulse_radiative_gradient
       atmo_struct%atmo_gradients(2,step_index) = pulse_gradient
       atmo_struct%atmo_gradients(3,step_index) = pt_scr%qdela
       atmo_struct%num_atm_points = step_index
       prev_tau = tau_now
       prev_density = density_now


! CHECK IF INTEGRATION COMPLETE
       if(dabs(indep_var - x_limit).le.step_tolerance) then
          if (solver_diagnostics()) write(run_log_unit,35)num_ok,num_bad,err_sum(1)
   35       format(1X,'ATMOSPHERE INTEGRATION COMPLETE',1X, &
                 'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
                 I5/5X,'MAXIMUM RELATIVE ERROR IN P ',1PE22.13)
          exit
       endif
       if(h_next.lt.h_min) h_next = h_min
       h_step = h_next
      end do
      if (step_index .gt. maxstp) then
! INTEGRATION HAS FAILED TO FINISH IN MAXSTP STEPS;
! PRINT NASTY MESSAGE AND QUIT.
      write(terminal_unit,50)
      write(run_log_unit,50)
   50 format(5X,'ATMOSPHERE INTEGRATION FAILED AFTER MAXSTP',1X, &
             'INTEGRATIONS.I QUIT.')
! 2026 (ROADMAP.md stage 3): stop converted to the ierr funnel below.
      jerr = 1
      ierr = jerr
      return
      end if
      exit atm_retry
      end do atm_retry
      end if
end subroutine integrate_atmosphere

! ---------------------------------------------------------------
! Integrate the envelope in pressure from the photosphere down to
! the fitting mass (M, R, T dependent), storing env_struct
! (inverted to inward-out at the end), with step limiting at the
! fitting point; ierr after maxstp steps.
subroutine integrate_envelope
      use math_lib
! DBG PULSE WRITE END OF DATA INDICATOR
! DBG
!  IF ENVELOPE MASS(SENV) SMALL ENOUGH,SKIP ENVELOPE INTEGRATION.
! DBG 2/92 CHANGED FROM 1.0D-10 to 1.0D-12
      if(star%senv.gt.-1.0d-12) then
       if(save_boundary_flag) then
          vtx_logp(vertex_index) = atm_table%atm_log10_pressure
          vtx_logr(vertex_index) = log10_radius
          vtx_logt(vertex_index) = atm_table%atm_log10_temperature
       endif
 230     format(4X,3F16.12,8X,F16.12)
      else
!  INITIALIZE VARIABLES AND SET NUMERICAL PARAMETERS.
      in_atmosphere = .false.
      do i = 1,3
       err_sum(i) = 0.0d0
      end do
      num_ok = 0
      num_bad = 0
      indep_var = atm_table%atm_log10_pressure
      initial_mass_coord = 0.0d0
      y(1) = initial_mass_coord
      y(2) = atm_table%atm_log10_temperature
      y(3) = log10_radius
      num_eqs = 3
      tolerance = star%ctrl%env_error_tol
      h_max = star%job%env_step_max
      h_min = star%job%env_step_min
      h_step = star%job%env_step_begin
      step_tolerance = dabs(tolerance_fraction*star%senv)
      if(stored_vertex_index.eq.vertex_index) stored_vertex_index = 0
      store_flag_set = .false.
      if (save_boundary_flag) then
         pressure_limit_test_flag = .false.
      else
         pressure_limit_test_flag = .true.
      end if
!  FIND DY/DX AT THE START OF THE STEP.
      call envelope_derivs(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
           temperature_rotation_factor,log10_gravity,in_atmosphere, &
           want_derivatives,conductive_opacity_flag,log10_radius, &
           log10_teff,hydrogen_fraction,metal_fraction,env_call_count,saha_state, jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
! DBG PULSE WRITE FIRST POINT OF ENVELOPE
! DBG
! STORE STARTING VALUES OF THE INTEGRATION
! 07/02 INITIALIZE NUMBER OF STORED ENVELOPE POINTS TO 1
      cz_in_envelope = .false.
      taucz_env_accum = 0.0d0
      env_struct%env_log10_density(1) = pt_scr%current_log10_density
      env_struct%env_log10_pressure(1) = pt_scr%current_log10_pressure
      env_struct%env_log10_radius(1) = pt_scr%current_log10_radius
      env_struct%env_log10_mass(1) = pt_scr%current_log10_mass
      env_struct%env_log10_temperature(1) = pt_scr%current_log10_temperature
      env_struct%env_hydrogen_fraction(1) = hydrogen_fraction
      env_struct%env_metal_fraction(1) = metal_fraction
      env_struct%env_convective_flag(1) = pt_scr%current_velocity.gt.0.0d0
! JVS 03/28
      env_struct%env_gradients(1,1) = pt_scr%current_gradients(1)
      env_struct%env_gradients(2,1) = pt_scr%current_gradients(2)
      env_struct%env_gradients(3,1) = pt_scr%current_gradients(3)
      env_struct%env_convective_velocity(1) = pt_scr%current_velocity
      env_struct%env_beta(1) = pt_scr%current_beta
! JVS end
! JVS 08/25
! Always save these
      chi_rho = 1.0d0/pt_scr%qqdp
      chi_t = -chi_rho*pt_scr%qqdt
      specific_heat_cv = pt_scr%qqcp - exp(ln10*(pt_scr%current_log10_pressure-pt_scr%current_log10_density- &
           pt_scr%current_log10_temperature))*chi_t**2/chi_rho
      env_struct%env_gamma1(1) = chi_rho*pt_scr%qqcp/specific_heat_cv
      env_struct%env_specific_heat_cp(1) = pt_scr%qqcp
      env_struct%env_ion_fraction(1,1) = pt_scr%current_ion_fraction(1)
      env_struct%env_ion_fraction(2,1) = pt_scr%current_ion_fraction(2)
      env_struct%env_ion_fraction(3,1) = pt_scr%current_ion_fraction(3)
      env_struct%env_opacity(1) = pt_scr%current_opacity
      env_struct%env_luminosity(1) = luminosity_linear
      env_struct%env_dlnrho_dlnt(1) = pt_scr%qqdt
! JVS 10/10
      unused_chdelj = pt_scr%current_gradients(2)
      unused_chdeld = pt_scr%qdela
      if(env_struct%env_convective_flag(1))cz_in_envelope = .true.
      env_struct%num_env_points = 1
      do step_index = 1,maxstp
       x_start = indep_var
       do i = 1,num_eqs
          y_start(i) = y(i)
          y_scale(i) = dabs(y(i)) + dabs(h_step*dydx(i))+tiny
       end do
       swap_temp = y(1) + h_step*dydx(1)
       if(star%senv - y(1).gt.0.0d0 .or. star%senv - swap_temp.gt.0.0d0) then
!  IF THE INTEGRATION HAS OVERSHOT THE FITTING POINT, OR THE NEXT
!  STEP WILL DO SO,LIMIT STEP SIZE OR INTEGRATE BACKWARDS TO THE
!  CORRECT MASS.
          h_step = (star%senv - y(1))/dydx(1)
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
       call bsstep(y,dydx,num_eqs,indep_var,h_step,tolerance,y_scale,h_did,h_next,envelope_derivs, &
              luminosity_linear,pressure_rotation_factor,temperature_rotation_factor, &
              log10_gravity,in_atmosphere,want_derivatives,conductive_opacity_flag, &
              log10_radius,log10_teff,hydrogen_fraction,metal_fraction, &
              env_call_count,saha_state,step_err,ierr=jerr)
! 2026 numerics-gate opt-in: same contract as the atmosphere-side
! bsstep call above.
       if (jerr /= 0) then
          ierr = numerics_termination
          return
       end if
       do i = 1,3
          err_sum(i) = err_sum(i) + step_err(i)
       end do
!  FIND DY/DX AT THE START OF THE NEXT STEP.
       call envelope_derivs(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
              temperature_rotation_factor,log10_gravity,in_atmosphere, &
              want_derivatives,conductive_opacity_flag,log10_radius, &
              log10_teff,hydrogen_fraction,metal_fraction,env_call_count,saha_state, jerr)
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
! DBG PULSE
! DBG END
       if(h_did.eq.h_step) then
          num_ok = num_ok + 1
       else
          num_bad = num_bad + 1
       endif
!  CHECK IF INTEGRATION COMPLETE
       mass_diff_remaining = star%senv - y(1)
! 07/02 STORE ENVELOPE TERMS IF THE INTEGRATION
! HAS NOT OVERSHOT THE FITTING POINT.
         if(mass_diff_remaining.le.step_tolerance)then
            env_struct%num_env_points = env_struct%num_env_points + 1
            env_struct%env_log10_density(env_struct%num_env_points) = pt_scr%current_log10_density
            env_struct%env_log10_pressure(env_struct%num_env_points) = pt_scr%current_log10_pressure
            env_struct%env_log10_radius(env_struct%num_env_points) = pt_scr%current_log10_radius
            env_struct%env_log10_mass(env_struct%num_env_points) = pt_scr%current_log10_mass
            env_struct%env_log10_temperature(env_struct%num_env_points) = pt_scr%current_log10_temperature
            env_struct%env_hydrogen_fraction(env_struct%num_env_points) = hydrogen_fraction
            env_struct%env_metal_fraction(env_struct%num_env_points) = metal_fraction
            env_struct%env_convective_flag(env_struct%num_env_points) = pt_scr%current_velocity.gt.0.0d0
! JVS 08/13 ADD RUN FOR CZ CALCULATION
            env_struct%env_gradients(1,env_struct%num_env_points) = pt_scr%current_gradients(1)
            env_struct%env_gradients(2,env_struct%num_env_points) = pt_scr%current_gradients(2)
            env_struct%env_gradients(3,env_struct%num_env_points) = pt_scr%current_gradients(3)
            env_struct%env_convective_velocity(env_struct%num_env_points) = pt_scr%current_velocity
            env_struct%env_beta(env_struct%num_env_points) = pt_scr%current_beta
! JVS 08/25
! Always save these
            chi_rho = 1.0d0/pt_scr%qqdp
            chi_t = -chi_rho*pt_scr%qqdt
            specific_heat_cv = pt_scr%qqcp - exp(ln10*(pt_scr%current_log10_pressure-pt_scr%current_log10_density- &
                 pt_scr%current_log10_temperature))*chi_t**2/chi_rho
            env_struct%env_gamma1(env_struct%num_env_points) = chi_rho*pt_scr%qqcp/specific_heat_cv
            env_struct%env_specific_heat_cp(env_struct%num_env_points) = pt_scr%qqcp
            env_struct%env_ion_fraction(1,env_struct%num_env_points) = pt_scr%current_ion_fraction(1)
            env_struct%env_ion_fraction(2,env_struct%num_env_points) = pt_scr%current_ion_fraction(2)
            env_struct%env_ion_fraction(3,env_struct%num_env_points) = pt_scr%current_ion_fraction(3)
            env_struct%env_opacity(env_struct%num_env_points) = pt_scr%current_opacity
            env_struct%env_luminosity(env_struct%num_env_points) = luminosity_linear
            env_struct%env_dlnrho_dlnt(env_struct%num_env_points) = pt_scr%qqdt

            if(.not.cz_in_envelope)then
               if(env_struct%env_convective_flag(env_struct%num_env_points))then
                  cz_in_envelope = .true.
               endif
            else if(pt_scr%current_velocity.gt.0.0d0)then
               delta_radius_cz = exp10(env_struct%env_log10_radius(env_struct%num_env_points-1)) - &
                    exp10(env_struct%env_log10_radius(env_struct%num_env_points))
               taucz_env_accum = taucz_env_accum + delta_radius_cz/pt_scr%current_velocity
            endif
         endif
       if(dabs(mass_diff_remaining).le.step_tolerance)then
!            WRITE(*,*)TAUCZENV/CSECYR
          if(save_boundary_flag) then
             interp_weight = mass_diff_remaining/(y_start(1)-y(1))
             vtx_logp(vertex_index) = indep_var + interp_weight*(x_start - indep_var)
             vtx_logr(vertex_index) = y(3) + interp_weight*(y_start(3) - y(3))
             vtx_logt(vertex_index) = y(2) + interp_weight*(y_start(2) - y(2))
             if(print_flag .and. solver_diagnostics())write(run_log_unit,230)vtx_logp(vertex_index),vtx_logt(vertex_index),vtx_logr(vertex_index),star%senv
          endif
          exit
       else if(.not.store_flag_set) then
          if(y(2).ge.star%tenv .and. save_boundary_flag) then
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
      end do
      if (step_index .gt. maxstp) then
! INTEGRATION HAS FAILED TO FINISH IN MAXSTP STEPS;
! PRINT NASTY MESSAGE AND QUIT.
      write(terminal_unit,911)
      write(run_log_unit,911)
 911  format(5X,'ENVELOPE INTEGRATION FAILED AFTER MAXSTP TRIES.',1X, &
           'I QUIT')
! 2026 (ROADMAP.md stage 3): stop converted to the ierr funnel below.
      jerr = 1
      ierr = jerr
      return
      end if
      end if
! 07/02 NOW INVERT THE ENVELOPE VECTOR.
      if(star%senv.lt.-1.0d-12)then
         do i = 1,env_struct%num_env_points
            inversion_index1 = i
            inversion_index2 = env_struct%num_env_points - i + 1
            if(inversion_index1.ge.inversion_index2)exit
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
! JVS 07/12 Save the last envelope point pressure
!      PPHOT = ENVP(NUMENV) ! G Somers 3/17, MOVED PPHOT DEF HIGHER UP
! END JVS
      endif
! DBG PULSE WRITE END OF DATA INDICATOR

      if (solver_diagnostics()) write(run_log_unit,215)num_ok,num_bad,mass_diff_remaining,y(1),(err_sum(jj),jj=1,3)
 215  format(1X,'ENVELOPE INTEGRATION COMPLETE',1X, &
           'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
           I5/5X,'SENV-LAST M=',1PE22.13,'  LAST M=',E22.13/ &
           5X,'MAX RELATIVE ERRORS:M ',1PE14.5,'  T ',E14.5, &
           '  R ',E14.5)

end subroutine integrate_envelope


end subroutine atm_get

end module envint_lib
