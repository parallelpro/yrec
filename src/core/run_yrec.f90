!----------------------------------------------------------------------
! run_yrec
!----------------------------------------------------------------------
! Added 2026 (phase five, step A -- the embeddable engine; see
! ROADMAP.md "Next milestone"). This is the entire body of the former
! `program main`, moved verbatim: namelist/controls read (parmin),
! table setup (setups), and the Monte-Carlo/run loop containing the
! full evolution driver. `program main` is now a thin CLI wrapper.
!
! The blanket `save` below is load-bearing: program-unit variables are
! implicitly static, and this body's locals (including its `data`
! initializations and rescale bookkeeping) rely on that; as subroutine
! locals they would otherwise be automatic. With save they are exactly
! as static as before -- which also means run_yrec is NOT yet
! re-entrant (one call per process); phase C of this milestone (the
! SAVE/reset audit) is what will make repeated calls legal.
!
! ierr: reserved for phase B, when the driver-side stops in this body
! (and the `if (jerr /= 0) stop` seams from phase-three stage 3)
! become error returns. In phase A it is set to 0 and the historical
! stops remain exactly where they were.
!
subroutine run_yrec(ierr)

! the array size, i.e. max # of shells is specified in the parameter
! statement. it defines JSON. to change the array size do a global
! change on "JSON=2000" or whatever.
      use net_lib
      use star_info_lib, only: star, i_c12, i_c13, i_h1, i_h2, i_he3, i_he4, i_lum_grav, i_metals, i_n14, i_n15, i_nu_b8, i_nu_be7, i_nu_f17, i_nu_hep, i_nu_n13, i_nu_o15, i_nu_pep, i_nu_pp, i_o16, i_o17, i_o18
      use yrec_output, only: output_run_header
      use luout_lib
      use const_lib
      use yrec_reset_lib, only: yrec_run_prologue
      use burn_lib
      implicit none
      integer :: step_status
      integer, parameter :: json = 5000
      integer, parameter :: numtt = 70
      integer, parameter :: numd = 19
      integer, parameter :: numx = 10
      integer, parameter :: numz = 13
      integer, parameter :: numxz = 126

! DBGLAOL - to save space make tables single precision
! MHP 8/25 Removed unused variables and added pass-through variables
! See the CROSS-CALLEE NAMING NOTE above for the file-path locals
! (declared further below, in the --- locals --- section, using
! setups.f90's/pdist.f90's own descriptive dummy-argument spelling).


































!     MHP 10/24 FLAG FOR END OF RUN



























! MHP 10/24 NEW VARIABLES FOR STOP CRITERIA ON CENTRAL ABUNDANCE are
! carried in common/sett/ above.

! latest values (Bahcall and Pinsonneault 1996)-actual values set in
! subroutine PARMIN
      double precision :: bp96_scale_factor(17)
      data bp96_scale_factor/0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0, &
           1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.92088,0.1625/
! MHP 3/96 added data for base solar age, L
!       DATA SUNAGE,SUNL/4.57D09,3.844D33/  ! KC 2025-05-31
      double precision :: reference_solar_luminosity
      data reference_solar_luminosity/3.844D33/

! --- locals ---
      integer :: monte_carlo_run_number
      double precision :: age_scale_factor
      integer :: convergence_iterations
      double precision :: initial_x_guess, initial_alpha_guess
      logical :: saved_use_structure_dt_limits
      integer :: saved_atm_choice
      logical :: compute_neutrino_fluxes
      double precision :: prev_mass_bound, curr_mass_bound, next_mass_bound
      double precision :: dlnrho_dlnt_unused, dlnrho_dlnp_unused
      integer:: i, j, k
      double precision :: shell_log_density, shell_log_temperature, &
           hydrogen_fraction, helium_fraction, metal_fraction, he3_fraction, &
           c12_fraction, c13_fraction, n14_fraction, n15_fraction, &
           o16_fraction, o17_fraction, o18_fraction, deuterium_fraction
      double precision :: pp_chain_energy_gen, he3he4_be7_electron_energy_gen, &
           he3he4_be7_proton_energy_gen, cno_cycle_energy_gen, &
           triple_alpha_energy_gen, dlnepsilon_dlnrho, dlnepsilon_dlnt, &
           total_energy_gen_rate
      integer :: shell_index
      double precision :: t6_million_k, log_electron_density, &
           zone_mass_fraction, zone_radius_fraction
      integer :: wrtlst_unit, model_iteration
      double precision :: log_r_rsun, current_zx, surface_z_over_x
      double precision :: initial_helium_fraction, initial_metal_fraction
      double precision :: central_temperature_mk, central_pressure_scaled, &
           central_density_linear
      character(len=256) :: alex06_table_path, allard_table_path, &
           atm_table_path, fermi_table_path, kurucz_table_path, &
           kurucz_table2_path, laol_table_path, laol_table2_path, &
           opal95_table_path, opal92_table_path
      character(len=256) :: zams_a_table_path, zams_b_table_path, &
           zams_c_table_path, centre1_table_path, centre2_table_path, &
           centre3_table_path, centre4_table_path, centre5_table_path
      character(len=256) :: opal92_table2_path, pure_z_table_path, &
           scv_h_table_path, scv_he_table_path, scv_z_table_path
      character(len=256) :: alex95_table_paths(7)
      double precision :: monte_helium_diffusion_fraction

      save
   ! INTENTIONAL: run-level state across MC runs; reset via yrec_reset
!*******
! START
!*******
      ! 2026 (ROADMAP.md stage 3): library errors return here via ierr;
      ! this driver-side call site preserves the historical stop.

      integer, intent(out) :: ierr
! load-bearing: see header
      save
   ! INTENTIONAL: run-level state across MC runs; reset via yrec_reset
      ierr = 0

! 2026 (phase five, step C): fresh-process semantics for repeated
! calls -- see core/yrec_reset.f90.
      call yrec_run_prologue

      call setversion()

      iowr = 9
! LPUNCH is TRUE once first model is calculated
      star%evo%punch_pending_flag = .false.
! 2026 (phase five): controls read and setup are now star-layer
! routines operating on the star_job structure (state/star_job_lib).
      call read_controls(ierr)
      if (ierr /= 0) return
      call star_setup(ierr)
      if (ierr /= 0) return

      do monte_carlo_run_number = star%job%mc_run_start,star%job%mc_run_end
! for monte carlo run, input values of parameters being changed.
      if (lmonte) then
         cross_section_scale(1) = star%run%s11_rate(monte_carlo_run_number)*bp96_scale_factor(1)
         cross_section_scale(2) = star%run%s33_rate(monte_carlo_run_number)*bp96_scale_factor(2)
         cross_section_scale(3) = star%run%s34_rate(monte_carlo_run_number)*bp96_scale_factor(3)
         cross_section_scale(16) = star%run%s17_rate(monte_carlo_run_number)*bp96_scale_factor(16)
         monte_helium_diffusion_fraction = star%run%helium_fraction_param(monte_carlo_run_number)
         fgrz = star%run%diffusion_factor(monte_carlo_run_number)
         solar_luminosity_cgs = reference_solar_luminosity*star%run%luminosity_target(monte_carlo_run_number)
         log10_solar_luminosity = dlog10(solar_luminosity_cgs)
         ln_solar_luminosity = ln10/solar_luminosity_cgs
         age_scale_factor = star%run%age_target(monte_carlo_run_number)
! timestep and final age are altered in SR SETCAL; input #s should be
! scaled for a solar age of 4.57 Gyr
         target_end_age(2)=1.0D8
         target_end_age(3)=4.57D9
      else
         age_scale_factor = 1.0D0
      endif
! DBG PULSE: save LPULSE flag, set LPULSE to F except on last model of
! last run, then set LPULSE to saved value of LPULSE.
      star%evo%saved_pulse_output_flag = pulsation_output_active
! 02/11 JVS uncommented LPULSE=.FALSE.
!      LPULSE = .FALSE.
! MHP 1/93 add option to automatically calibrate solar model.
! MHP 3/96 added counter for # of iterations per converged model and
! starting estimate of ALPHA and X
      if (calibrate_solar_model) then
         call setcal(age_scale_factor)
         convergence_iterations = 1
         initial_x_guess = rescale_params(2,1)
         initial_alpha_guess = mixing_length_array(1)
         saved_use_structure_dt_limits = use_structure_dt_limits   ! save LPTIME for reuse during calibration
         saved_atm_choice  = atm_choice    ! save KTTAU for reuse during calibration
      else
         convergence_iterations = 0
      endif
! DBG 12/94 add option to automatically calculate a stellar model
! of specified Teff and L
      if (calibrate_star_flag) then
         call setscal
      endif

!**********
!     RUN THROUGH THE KIND CARDS IN ORDER
!**********
      run_loop: do nk = 1, num_runs
         star%run%sound_speed_output_active = .false.
!         LPULSE=.FALSE.
         initial_envelope_x = initial_x_array(nk)
         initial_envelope_z = initial_z_array(nk)
         cmixl = mixing_length_array(nk)
       change_envelope_mass_flag = has_senv0_array(nk)
       requested_envelope_mass = senv0_array(nk)
       star%evo%reset_triangle = .false.
       star%evo%model_diverged_flag = .false.
! MHP 10/02 ZERO OUT INITIAL ANGULAR MOMENTUM
         star%evo%total_angular_momentum = 0.0D0
         star%evo%total_rotational_ke = 0.0D0
! read in the initial model here
! STARIN also calls RSCALE to perform rescaling if requested
!        CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
       call starin(star%evo%timestep_yr, star%evo%dt, star%evo%hydrogen_dt, star%evo%trial_sign_flag, &
            star%evo%ikut_flag, star%evo%istore_flag, star%evo%model_diverged_flag, &
            star%evo%recompute_envelope_triangle, nk, star%evo%dlnrho_dlnp, star%evo%dlnrho_dlnt, &
            star%evo%total_angular_momentum, star%evo%total_rotational_ke, &
            star%evo%convective_velocity, star%job%mixture_weights, ierr)
       if (ierr /= 0) return

      if ((star%omega(1) .eq. 0) .and. (rotation_active)) then

1611      format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                 ' Initialize rotation rates or set LROT to', &
                 ' FALSE.')
          print 1611
          ! 2026 (phase five, step B): configuration error returns to the
          ! CLI wrapper (which stops) instead of stopping here.
          ierr = 1
          return
      endif
!     MHP 10/24 CHECK STOP CONDITIONS AND DISABLE THEM IF THE STARTING VALUES ARE BELOW THE TARGET THRESHOLD
         if (end_age_stop_active(nk)) then
            if (central_deuterium_stop(nk).gt.0.0D0 .and. &
                 star%xa(i_h2,1).lt.central_deuterium_stop(nk)) then
               central_deuterium_stop(nk)=-central_deuterium_stop(nk)
               write(*,101)star%xa(i_h2,1),central_deuterium_stop(nk)
               write(short_file_unit,101)star%xa(i_h2,1),central_deuterium_stop(nk)
 101           format('STARTING D ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_hydrogen_stop(nk).gt.0.0D0 .and. &
                 star%xa(i_h1,1).lt.central_hydrogen_stop(nk)) then
               central_deuterium_stop(nk)=-central_hydrogen_stop(nk)
               write(*,102)star%xa(i_h2,1),central_deuterium_stop(nk)
               write(short_file_unit,102)star%xa(i_h2,1),central_deuterium_stop(nk)
 102           format('STARTING X ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_helium_stop(nk).gt.0.0D0 .and. &
                 star%xa(i_he4,1).lt.central_helium_stop(nk)) then
               central_helium_stop(nk)=-central_helium_stop(nk)
               write(*,103)star%xa(i_h2,1),central_deuterium_stop(nk)
               write(short_file_unit,103)star%xa(i_h2,1),central_deuterium_stop(nk)
 103           format('STARTING Y ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
         endif
!     MHP 2/04 NEUTRINO TABLE
!      LNUTAB = .TRUE.
      compute_neutrino_fluxes = .false.
      if (compute_neutrino_fluxes) then
! SET UP WEIGHTS AND MASSES.
! HS1 = LOCATION IN GM (UNLOGGED) OF SHELL CENTERS.
! HS2 = MASS IN GM OF EACH SHELL.
      curr_mass_bound = exp(ln10*star%log_mass(1))
      prev_mass_bound = - curr_mass_bound
      do i = 2,star%nz
         next_mass_bound = prev_mass_bound
         prev_mass_bound = curr_mass_bound
         curr_mass_bound = exp(ln10*star%log_mass(i))
         star%m(i-1) = prev_mass_bound
         star%dm(i-1) = 0.5D0*(curr_mass_bound-next_mass_bound)
      end do
      star%m(star%nz) = curr_mass_bound
      star%dm(star%nz) = exp(ln10*star%log_total_mass) - 0.5D0*(prev_mass_bound+curr_mass_bound)
      dlnrho_dlnt_unused = -1.0D0
      dlnrho_dlnp_unused = 1.0D0
      do j = 1,10
         star%flux%neutrino_flux_total(j) = 0.0D0
         do k = 1,star%nz
            star%neutrino_flux_zone(j,k) = 0.0D0
         end do
      end do
! ASSIGN LOCAL VARIABLES FOR SR CALL FROM GLOBAL VECTORS.
      do i = 1,star%nz
         shell_log_density = star%logRho(i)
         shell_log_temperature = star%logT(i)
! SKIP CALCULATIONS FOR LOW TEMPERATURES.
         if (shell_log_temperature.lt.6.0D0) exit
         hydrogen_fraction = star%xa(i_h1,i)
         helium_fraction = star%xa(i_he4,i)
         metal_fraction = star%xa(i_metals,i)
         he3_fraction = star%xa(i_he3,i)
         c12_fraction = star%xa(i_c12,i)
         c13_fraction = star%xa(i_c13,i)
         n14_fraction = star%xa(i_n14,i)
         n15_fraction = star%xa(i_n15,i)
         o16_fraction = star%xa(i_o16,i)
         o17_fraction = star%xa(i_o17,i)
         o18_fraction = star%xa(i_o18,i)
         call engeb(pp_chain_energy_gen,he3he4_be7_electron_energy_gen, &
              he3he4_be7_proton_energy_gen,cno_cycle_energy_gen, &
              triple_alpha_energy_gen,dlnepsilon_dlnrho,dlnepsilon_dlnt, &
              total_energy_gen_rate,shell_log_density, &
!      *TL,PDT,PDP,X,Y,Z,XHE3,XC12,XC13,XN14,XN15,XO16,XO17,
!      *XO18,XH2,XLI6,XLI7,XBE9,I,HR1,HR2,HR3,HR4,HR5,HR6,HR7,  ! KC 2025-05-31
              shell_log_temperature,hydrogen_fraction,helium_fraction, &
              he3_fraction,c12_fraction,c13_fraction,n14_fraction,o16_fraction, &
              o18_fraction,deuterium_fraction,shell_index,star%reaction_rate_1, &
              star%reaction_rate_2,star%reaction_rate_3,star%reaction_rate_4,star%reaction_rate_5, &
              star%reaction_rate_6,star%reaction_rate_7,star%reaction_rate_8,star%reaction_rate_9, &
              star%reaction_rate_10,star%reaction_rate_11,star%reaction_rate_12, &
              star%reaction_rate_13,star%n15_alpha_branch_fraction, &
              star%be7_electron_capture_fraction)
! BE7 MASS FRACTION.
         star%be7_mass_fraction_zone(i) = star%engeb%be7_mass_fraction
! CONVERT FROM ERG/GM/S TO ERG/S FOR EACH SHELL BY MULTIPLYING
! BY THE MASS OF EACH SHELL IN GM (HS2).
         do j = 1,10
            star%neutrino_flux_zone(j,i) = star%flux%neutrino_flux(j)*star%dm(i)
            star%flux%neutrino_flux_total(j) = star%flux%neutrino_flux_total(j) + star%neutrino_flux_zone(j,i)
         end do
         write(*,911)i,star%dm(i),(star%neutrino_flux_zone(j,i),j=1,10)
 911     format(I5,1P11E10.3)
      end do
! WRITE OUT TOTAL NEUTRINO FLUXES.
! ***NOTE THAT THESE ARE IN UNITS OF 10**10. ***
      write(76,222)(star%flux%neutrino_flux_total(i),i=1,10)
! NORMALIZE FLUXES.
      do j = 1,10
         do i = 1,star%nz
            star%neutrino_flux_zone(j,i) = star%neutrino_flux_zone(j,i)/star%flux%neutrino_flux_total(j)
         end do
      end do
      do i = 1,star%nz
! TEMPERATURE IN UNITS OF 10**6 K.
         t6_million_k = exp(ln10*(star%logT(i)-6.0D0))
         if (t6_million_k.lt.5.0D0) exit
! ELECTRON DENSITY.
         log_electron_density = star%logRho(i)+log10((1.0D0+star%xa(i_h1,i))/2.0D0)
! MASS FRACTION.
         zone_mass_fraction = star%dm(i)/1.9891D33
! RADIUS FRACTION.
         zone_radius_fraction = exp(ln10*star%logR(i))/solar_radius_cgs
! FLUXES ARE PRINTED IN THE SAME ORDER AS BAHCALL AND PINSONNEAULT.
         write(76,145)zone_radius_fraction,t6_million_k,log_electron_density, &
         zone_mass_fraction,star%be7_mass_fraction_zone(i),star%neutrino_flux_zone(i_nu_pp,i), &
         star%neutrino_flux_zone(i_nu_b8,i), &
         star%neutrino_flux_zone(i_nu_n13,i), &
         star%neutrino_flux_zone(i_nu_o15,i),star%neutrino_flux_zone(i_nu_f17,i),star%neutrino_flux_zone(i_nu_be7,i), &
         star%neutrino_flux_zone(i_nu_pep,i),star%neutrino_flux_zone(i_nu_hep,i)
  145    format(F9.5,F7.3,F6.3,1P10E10.3)
      end do
  222    format(1P10E10.3)
!         IF(M.GT.1)STOP999
      endif
! save mass in solar units
         pulsation_mass_msun=star%star_mass
! MHP 08/02 STORE STARTING CZ PROPERTIES
         star%light_burn%jcz = star%envelope_cz_bottom_index
         star%turnover%convective_turnover_timescale = 0.0D0
! write out headers of the appropriate output files
      call output_run_header(star%star_mass)
! DBG PULSE OUT 7/92
! initialize variables for calculating when to dump pulse output
         star%evo%prev_log_l = star%log_L
         star%evo%prev_log_teff = star%log_Teff
         star%evo%prev_age = star%run%dage
         star%evo%path_length_sq = 0.0D0

       if (helium_flash_active) then
! timestep cutting requires a model stored in logical unit ILAST
! or it will crash - so copy initial model to unit ILAST
          if (star%evo%punch_pending_flag) then
             wrtlst_unit = ilast
             call wrtlst(wrtlst_unit,star%xa,star%logRho,star%luminosity_lsun, &
                  star%logP,star%logR,star%log_mass,star%logT,star%convective_flag, &
                  star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
                  star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs, &
                  star%evo%trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index, &
                  star%envelope_cz_bottom_index,star%model_number,star%nz, &
                  star%star_mass,star%log_Teff,star%log_L,star%log_total_mass,star%run%dage, &
                  star%evo%timestep_yr,star%omega)
          endif
       endif

! locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if applicable).
         call findsh(star%xa,star%luminosity_lsun,star%convective_flag,star%nz, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%evo%h_shell_zone_begin, &
              star%evo%h_shell_end_index,star%evo%h_shell_midpoint_zone,star%evo%has_h_shell)
! determine timestep for model
! JVS 04/14 Added Teffl to passed variables
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *               JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(star%evo%dt,star%evo%hydrogen_dt,star%nz,star%logRho,star%luminosity_lsun, &
            star%m,star%dm,star%logT,star%xa,star%core_cz_top_index, &
            star%evo%h_shell_midpoint_zone,star%luminosity_breakdown,star%run%dage,star%evo%timestep_yr,nk, &
            star%logP,star%logR,star%omega,star%evo%max_domega_frac,star%evo%h_shell_zone_begin, &
            star%log_Teff)

       star%evo%dt_saved = star%evo%dt
! zero out entropy terms.
         do i = 1,star%nz
            star%run%temperature_entropy_term(i) = 0.0D0
            star%run%pressure_entropy_term(i) = 0.0D0
            star%run%luminosity_entropy_term(i) = 0.0D0
            star%run%radius_entropy_term(i) = 0.0D0
         end do

! zero out light element burning rates in the surface CZ.
         if (use_extended_composition) then
            star%light_burn%log_rate_li6_prev = 0.0D0
            star%light_burn%log_rate_li7_prev = 0.0D0
            star%light_burn%log_rate_be9_prev = 0.0D0
         endif

! for a given kind card, evolve NMODLS(NK) times
! if rescaling is being performed, NMODLS(NK) is the number of times
! the new model is being relaxed
       do model_iteration = 1,num_models(nk)
! 2026 (phase five): one model advance per iteration, extracted to
! core/evolve_step.f90 (see its header for the step_status contract).
       call evolve_step(model_iteration, step_status, ierr)
       if (ierr /= 0) return
       if (step_status == 1) exit
       if (step_status == 2) cycle run_loop   ! (was goto 200, the run-loop terminator)
       end do

! G Somers 11/14, CHANGE CALL TO PUTSTORE INSTEAD OF WRTLST.
! STORE LAST MODEL IN ISTOR IF LSTORE, LSTPCH, AND LPUNCH ARE .TRUE.
         if (lstore.and.lstpch.and.star%evo%punch_pending_flag) then
          call putstore(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,star%evo%trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index,star%envelope_cz_bottom_index,star%model_number,star%nz,star%star_mass,star%log_Teff,star%log_L,star%log_total_mass, &
                 star%run%dage,star%evo%timestep_yr,star%omega,star%m,star%eta_squared,star%mean_radius,star%fp_rot,star%ft_rot,star%j_rot,star%i_rot)
            star%evo%punch_pending_flag = .false.
       endif
! 110  CONTINUE
! G Somers END

! MHP 1/93 CHECK AUTOMATIC CALIBRATATION OF SOLAR MODEL.
!c MHP 5/96 changed solar calibration to perform solar models in 3 kind cards
         if (calibrate_solar_model) then
! JVS Turn off calcad - speeds things up
            compute_acoustic_depth=.false.
            if (mod(nk,3).eq.0) then
               log_r_rsun = 0.5D0*(star%log_L+log10_solar_luminosity-c4pil-csigl-4.0D0*star%log_Teff)-log10_solar_radius
! MHP 06/13 Add solar Z/X to observables
               current_zx = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
               call chkcal(star%log_L,log_r_rsun,nk,current_zx)
!               CALL CHKCAL(BL,RLL,NK)
               use_structure_dt_limits = saved_use_structure_dt_limits  ! Restore LPTIME to original value for next cycle
               atm_choice  = saved_atm_choice    ! Restore KTTAU to original value for next cycle
               if (star%run%solar_calibration_active) then
                  exit
               else
!c MHP 8/96 added counter for # of runs needed for calibration
                  convergence_iterations = convergence_iterations + 1
! MHP 6/97 STOP AFTER 10 ATTEMPTS AT CALIBRATION
!                  IF(ICONV.GE.11) GOTO 250
                  if (convergence_iterations.ge.15) exit
                  if (pulsation_output_active) then
! DBG 6/93 Need to delete pulse output because have not got ultimate
! model yet.
! MHP 8/25 Replaced delete file with rewind file. This is functionally the same and avoids the need to pass the character string for the file name from parmin.
                     rewind(opal_model_unit)
                     rewind(opal_envelope_unit)
                     rewind(opal_atm_unit)
!                     CLOSE(IOPMOD, STATUS='DELETE')
!                     CLOSE(IOPENV, STATUS='DELETE')
!                     CLOSE(IOPATM, STATUS='DELETE')
!                     OPEN(IOPMOD, FILE=FPMOD,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
!                     OPEN(IOPENV, FILE=FPENV,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
!                     OPEN(IOPATM, FILE=FPATM,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
                  end if
               end if
            endif
         endif

! DBG 12/94 NO MORE RUNS NEEDED. HAVE CALIBRATED STELLAR MODEL
         if (calibrate_star_flag .and. star_found_flag.and.(mod(nk,2).eq.0)) exit

! END RUN LOOP
      end do run_loop
! EXIT RUN LOOP

! FOR MONTE CARLO, REWIND OUTPUT FILES AND WRITE OUT SNU FLUXES AND
! MODEL PARAMTERS TO AN OUTPUT FILE.
! RUN FAILED TO CONVERGE.  WRITE FINAL INFO WITH WARNING NOTE.
! 2026 MESA-style output: this whole chain is legacy-file machinery
! (it rewinds the legacy units between calibration/MC cycles and
! writes the .snu summaries) -- legacy mode only. In MESA mode the
! history file simply accumulates every calibration cycle instead.
      if (use_legacy_output) then
      if (lmonte .and. convergence_iterations.ge.11 .and. .not.star%run%solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)
         write(neutrino_unit,1525)star%log_L,log_r_rsun
 1525    format(5X,'DID NOT CONVERGE WITHIN 10 ATTEMPTS L,R',2F10.6)
! MONTE CARLO #, CONVERGED MIXING LENGTH AND INITIAL H, SURFACE X,
! SURFACE Z, Z/X, CENTRAL X, CENTRAL Z
         write(neutrino_unit,1519) monte_carlo_run_number,mixing_length_array(nk),rescale_params(2,nk-2),star%xa(i_h1,star%nz), &
              star%xa(i_metals,star%nz),surface_z_over_x,star%xa(i_h1,1),star%xa(i_metals,1)
 1519    format(1X,I5,3F10.6,4E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,star%run%dlum_dx,star%run%drad_dx,star%run%dlum_dalpha,star%run%drad_dalpha
! SUMMARY OF STRUCTURE : TC, RHOC, PC
         write(neutrino_unit, 1517)star%run%central_log10_temperature,star%run%central_log10_pressure,star%run%central_log10_density, &
              star%xa(i_h1,1),star%xa(i_metals,1)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) star%flux%cl37_snu_rate,star%flux%ga71_snu_rate,(star%flux%neutrino_flux_total(i),i=1,8)
!          CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *        DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *        JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *        CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
         call wrtmonte(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%run%dage, &
              star%evo%timestep_yr,star%star_mass,star%log_Teff,star%log_L, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
              star%envelope_fit_coeffs,star%evo%trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
      else if (calibrate_solar_model .and. lsnu .and. star%run%solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)

         surface_z_over_x = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
! HEADER FILE:  MONTE CARLO PARAMETERS
         if (lmonte) then
            write(neutrino_unit,1520)monte_carlo_run_number,star%run%s11_rate(monte_carlo_run_number),star%run%s33_rate(monte_carlo_run_number),star%run%s34_rate(monte_carlo_run_number),star%run%s17_rate(monte_carlo_run_number), &
                 star%run%metal_to_h_ratio(monte_carlo_run_number),star%run%helium_fraction_param(monte_carlo_run_number),star%run%diffusion_factor(monte_carlo_run_number),star%run%luminosity_target(monte_carlo_run_number),star%run%age_target(monte_carlo_run_number)
         endif
 1520    format(I7,1P9E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,star%run%dlum_dx,star%run%drad_dx,star%run%dlum_dalpha,star%run%drad_dalpha
 1518    format(1X,I2,2F10.6,1P4E11.4)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) star%flux%cl37_snu_rate,star%flux%ga71_snu_rate,(star%flux%neutrino_flux_total(i),i=1,10)
 1516    format(1X,2F8.3,1P10E10.3)
! SUMMARY OF STRUCTURE : TC, RHOC, PC, XC, ZC (ADD MU C)
         central_temperature_mk = 10.0D0**(star%run%central_log10_temperature-6.0D0)
         central_pressure_scaled = 10.0D0**(star%run%central_log10_pressure-17.0D0)
         central_density_linear = 10.0D0**star%run%central_log10_density
         write(neutrino_unit, 1517)central_temperature_mk,central_density_linear,central_pressure_scaled,star%xa(i_h1,1),star%xa(i_metals,1)
 1517    format(1X,F7.3,F7.2,F6.3,2F8.5)
! INITIAL ALPHA,Y,Z,ALPHA; FINAL R, L
         initial_helium_fraction = 1.0D0 - rescale_params(2,nk-2) - rescale_params(3,nk-2)
         initial_metal_fraction = rescale_params(3,nk-2)
         write(neutrino_unit,1521)mixing_length_array(nk),initial_helium_fraction,initial_metal_fraction,star%log_L,log_r_rsun
 1521    format(F7.4,2F8.5,1P2E10.3)
! CZ DEPTH (R,M), SURFACE Y, Z, Z/X (ADD T CZ BASE, RHO CZ BASE)
         write(neutrino_unit,1522)star%run%envelope_radius,star%run%envelope_mass,star%xa(i_he4,star%nz),star%xa(i_metals,star%nz),surface_z_over_x
 1522    format(F8.5,F9.6,2F8.5,F9.6)
! ENERGY GENERATION FRACTIONS PP I,II,III,CNO,EGRAV
         write(neutrino_unit,1523)(star%luminosity_breakdown(j),j=1,4),star%luminosity_breakdown(i_lum_grav)
 1523    format(1P5E10.3)
         if (lmonte) then
!             CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *           DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *           JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *           CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
            call wrtmonte(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%run%dage, &
                 star%evo%timestep_yr,star%star_mass,star%log_Teff,star%log_L, &
                 star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,star%evo%trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
         endif
      endif
      end if
      end do

! 2026 (phase five, step B): the normal end-of-job stop became this
! clean return (ierr stays 0); the CLI wrapper simply ends.
      return
end subroutine run_yrec
