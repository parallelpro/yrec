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
      use nuclear_lib
      use star_info_lib, only: star
      use fluxes_lib
      use engeb_diag_lib
      use light_burn_lib
      use turnover_lib
      use star_info_lib, only: star
      use luout_lib
      use const_lib
      use star_info_lib, only: star
      use star_job_lib, only: job
      implicit none
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


































      double precision :: trial_sign_flag
!     MHP 10/24 FLAG FOR END OF RUN
      logical :: end_kind_flag























      integer :: nao
      data nao/1/




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
      logical :: saved_pulse_output_flag
      integer :: convergence_iterations
      double precision :: initial_x_guess, initial_alpha_guess
      logical :: saved_use_structure_dt_limits
      integer :: saved_atm_choice
      logical :: reset_triangle, model_diverged_flag
      double precision :: total_angular_momentum, total_rotational_ke
      integer :: ikut_flag, istore_flag
      double precision :: delta_time, hydrogen_dt, timestep_yr
      double precision :: convective_velocity
      logical :: compute_neutrino_fluxes
      double precision :: prev_mass_bound, curr_mass_bound, next_mass_bound
      double precision :: dlnrho_dlnt_unused, dlnrho_dlnp_unused
      integer :: i, j, k, ii
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
      double precision :: prev_log_l, prev_log_teff, prev_age, path_length_sq
      integer :: wrtlst_unit, model_iteration
      integer :: h_shell_zone_begin, h_shell_midpoint_zone, &
           h_shell_end_index
      logical :: has_h_shell
      double precision :: delta_time_saved
      logical :: evolve_model_flag
      integer :: num_species
      double precision :: target_envelope_mass
      logical :: new_atmosphere_fit_needed
      integer :: envelope_cz_zone_prev, envelope_cz_zone_end
      logical :: converged, in_atmosphere, want_derivatives, mixing_active, &
           conductive_opacity_flag
      integer :: iteration_level, max_iterations, iterations_done
      double precision :: dlnrho_dlnt, dlnrho_dlnp
      logical :: recompute_surface_bc, recompute_envelope_triangle
      logical :: use_correct_gradients
      integer :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      double precision :: log_gravity
      double precision :: teff_kelvin_unused
      double precision :: log_r_rsun, current_zx, surface_z_over_x
      double precision :: initial_helium_fraction, initial_metal_fraction
      double precision :: central_temperature_mk, central_pressure_scaled, &
           central_density_linear
      logical :: wind_loss_active
      double precision :: max_domega_frac
      integer :: itrot
      double precision :: delta_temp_step, delta_pressure_step, &
           delta_lum_step, delta_radius_step
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
      logical :: punch_pending_flag

      save

!*******
! START
!*******
      ! 2026 (ROADMAP.md stage 3): library errors return here via ierr;
      ! this driver-side call site preserves the historical stop.
      integer :: jerr

      integer, intent(out) :: ierr
! load-bearing: see header
      save

      ierr = 0

      call setversion()

      iowr = 9
! LPUNCH is TRUE once first model is calculated
      punch_pending_flag = .false.
! 2026 (phase five): controls read and setup are now star-layer
! routines operating on the star_job structure (state/star_job_lib).
      call read_controls(ierr)
      if (ierr /= 0) return
      call star_setup(ierr)
      if (ierr /= 0) return

      do 500 monte_carlo_run_number = job%mc_run_start,job%mc_run_end
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
      saved_pulse_output_flag = pulsation_output_active
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
      do 200 nk = 1, num_runs
         star%run%sound_speed_output_active = .false.
!         LPULSE=.FALSE.
         initial_envelope_x = initial_x_array(nk)
         initial_envelope_z = initial_z_array(nk)
         cmixl = mixing_length_array(nk)
       change_envelope_mass_flag = has_senv0_array(nk)
       requested_envelope_mass = senv0_array(nk)
       reset_triangle = .false.
       model_diverged_flag = .false.
! MHP 10/02 ZERO OUT INITIAL ANGULAR MOMENTUM
         total_angular_momentum = 0.0D0
         total_rotational_ke = 0.0D0
! read in the initial model here
! STARIN also calls RSCALE to perform rescaling if requested
!        CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
       call starin(timestep_yr, delta_time, hydrogen_dt, trial_sign_flag, &
            ikut_flag, istore_flag, model_diverged_flag, &
            recompute_envelope_triangle, nk, dlnrho_dlnp, dlnrho_dlnt, &
            total_angular_momentum, total_rotational_ke, &
            convective_velocity, job%mixture_weights, ierr)
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
                 star%composition(12,1).lt.central_deuterium_stop(nk)) then
               central_deuterium_stop(nk)=-central_deuterium_stop(nk)
               write(*,101)star%composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,101)star%composition(12,1),central_deuterium_stop(nk)
 101           format('STARTING D ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_hydrogen_stop(nk).gt.0.0D0 .and. &
                 star%composition(1,1).lt.central_hydrogen_stop(nk)) then
               central_deuterium_stop(nk)=-central_hydrogen_stop(nk)
               write(*,102)star%composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,102)star%composition(12,1),central_deuterium_stop(nk)
 102           format('STARTING X ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_helium_stop(nk).gt.0.0D0 .and. &
                 star%composition(2,1).lt.central_helium_stop(nk)) then
               central_helium_stop(nk)=-central_helium_stop(nk)
               write(*,103)star%composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,103)star%composition(12,1),central_deuterium_stop(nk)
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
      do i = 2,star%num_zones
         next_mass_bound = prev_mass_bound
         prev_mass_bound = curr_mass_bound
         curr_mass_bound = exp(ln10*star%log_mass(i))
         star%enclosed_mass(i-1) = prev_mass_bound
         star%shell_mass(i-1) = 0.5D0*(curr_mass_bound-next_mass_bound)
      end do
      star%enclosed_mass(star%num_zones) = curr_mass_bound
      star%shell_mass(star%num_zones) = exp(ln10*star%log_total_mass) - 0.5D0*(prev_mass_bound+curr_mass_bound)
      dlnrho_dlnt_unused = -1.0D0
      dlnrho_dlnp_unused = 1.0D0
      do j = 1,10
         flux_diag%neutrino_flux_total(j) = 0.0D0
         do k = 1,star%num_zones
            star%neutrino_flux_zone(j,k) = 0.0D0
         end do
      end do
! ASSIGN LOCAL VARIABLES FOR SR CALL FROM GLOBAL VECTORS.
      do i = 1,star%num_zones
         shell_log_density = star%log_density(i)
         shell_log_temperature = star%log_temperature(i)
! SKIP CALCULATIONS FOR LOW TEMPERATURES.
         if (shell_log_temperature.lt.6.0D0) goto 666
         hydrogen_fraction = star%composition(1,i)
         helium_fraction = star%composition(2,i)
         metal_fraction = star%composition(3,i)
         he3_fraction = star%composition(4,i)
         c12_fraction = star%composition(5,i)
         c13_fraction = star%composition(6,i)
         n14_fraction = star%composition(7,i)
         n15_fraction = star%composition(8,i)
         o16_fraction = star%composition(9,i)
         o17_fraction = star%composition(10,i)
         o18_fraction = star%composition(11,i)
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
         star%be7_mass_fraction_zone(i) = engeb_diag%be7_mass_fraction
! CONVERT FROM ERG/GM/S TO ERG/S FOR EACH SHELL BY MULTIPLYING
! BY THE MASS OF EACH SHELL IN GM (HS2).
         do j = 1,10
            star%neutrino_flux_zone(j,i) = flux_diag%neutrino_flux(j)*star%shell_mass(i)
            flux_diag%neutrino_flux_total(j) = flux_diag%neutrino_flux_total(j) + star%neutrino_flux_zone(j,i)
         end do
         write(*,911)i,star%shell_mass(i),(star%neutrino_flux_zone(j,i),j=1,10)
 911     format(I5,1P11E10.3)
      end do
  666 continue
! WRITE OUT TOTAL NEUTRINO FLUXES.
! ***NOTE THAT THESE ARE IN UNITS OF 10**10. ***
      write(76,222)(flux_diag%neutrino_flux_total(i),i=1,10)
! NORMALIZE FLUXES.
      do j = 1,10
         do i = 1,star%num_zones
            star%neutrino_flux_zone(j,i) = star%neutrino_flux_zone(j,i)/flux_diag%neutrino_flux_total(j)
         end do
      end do
      do i = 1,star%num_zones
! TEMPERATURE IN UNITS OF 10**6 K.
         t6_million_k = exp(ln10*(star%log_temperature(i)-6.0D0))
         if (t6_million_k.lt.5.0D0) goto 141
! ELECTRON DENSITY.
         log_electron_density = star%log_density(i)+log10((1.0D0+star%composition(1,i))/2.0D0)
! MASS FRACTION.
         zone_mass_fraction = star%shell_mass(i)/1.9891D33
! RADIUS FRACTION.
         zone_radius_fraction = exp(ln10*star%log_radius(i))/solar_radius_cgs
! FLUXES ARE PRINTED IN THE SAME ORDER AS BAHCALL AND PINSONNEAULT.
         write(76,145)zone_radius_fraction,t6_million_k,log_electron_density, &
         zone_mass_fraction,star%be7_mass_fraction_zone(i),star%neutrino_flux_zone(1,i), &
         star%neutrino_flux_zone(5,i), &
         star%neutrino_flux_zone(6,i), &
         star%neutrino_flux_zone(7,i),star%neutrino_flux_zone(8,i),star%neutrino_flux_zone(4,i), &
         star%neutrino_flux_zone(2,i),star%neutrino_flux_zone(3,i)
  145    format(F9.5,F7.3,F6.3,1P10E10.3)
      end do
  141 continue
  222    format(1P10E10.3)
!         IF(M.GT.1)STOP999
      endif
! save mass in solar units
         pulsation_mass_msun=star%total_mass_msun
! MHP 08/02 STORE STARTING CZ PROPERTIES
         light_burn%jcz = star%envelope_cz_bottom_index
         turnover%convective_turnover_timescale = 0.0D0
! write out headers of the appropriate output files
      call wrthead(star%total_mass_msun)
! DBG PULSE OUT 7/92
! initialize variables for calculating when to dump pulse output
         prev_log_l = star%log10_luminosity
         prev_log_teff = star%log_teff
         prev_age = star%run%dage
         path_length_sq = 0.0D0

       if (helium_flash_active) then
! timestep cutting requires a model stored in logical unit ILAST
! or it will crash - so copy initial model to unit ILAST
          if (punch_pending_flag) then
             wrtlst_unit = ilast
             call wrtlst(wrtlst_unit,star%composition,star%log_density,star%luminosity_lsun, &
                  star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag, &
                  star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
                  star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs, &
                  trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index, &
                  star%envelope_cz_bottom_index,star%model_number,star%num_zones, &
                  star%total_mass_msun,star%log_teff,star%log10_luminosity,star%log_total_mass,star%run%dage, &
                  timestep_yr,star%omega)
          endif
       endif

! locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if applicable).
         call findsh(star%composition,star%luminosity_lsun,star%convective_flag,star%num_zones, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,h_shell_zone_begin, &
              h_shell_end_index,h_shell_midpoint_zone,has_h_shell)
! determine timestep for model
! JVS 04/14 Added Teffl to passed variables
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *               JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(delta_time,hydrogen_dt,star%num_zones,star%log_density,star%luminosity_lsun, &
            star%enclosed_mass,star%shell_mass,star%log_temperature,star%composition,star%core_cz_top_index, &
            h_shell_midpoint_zone,star%luminosity_breakdown,star%run%dage,timestep_yr,nk, &
            star%log_pressure,star%log_radius,star%omega,max_domega_frac,h_shell_zone_begin, &
            star%log_teff)

       delta_time_saved = delta_time
! zero out entropy terms.
         do 99 i = 1,star%num_zones
            star%run%temperature_entropy_term(i) = 0.0D0
            star%run%pressure_entropy_term(i) = 0.0D0
            star%run%luminosity_entropy_term(i) = 0.0D0
            star%run%radius_entropy_term(i) = 0.0D0
   99    continue

! zero out light element burning rates in the surface CZ.
         if (use_extended_composition) then
            light_burn%log_rate_li6_prev = 0.0D0
            light_burn%log_rate_li7_prev = 0.0D0
            light_burn%log_rate_be9_prev = 0.0D0
         endif

! for a given kind card, evolve NMODLS(NK) times
! if rescaling is being performed, NMODLS(NK) is the number of times
! the new model is being relaxed
       do 100 model_iteration = 1,num_models(nk)
! rewind ISHORT if LRWSH is true (keeps ISHORT small)
          if (lrwsh_placeholder) then
             rewind(short_file_unit)
          endif
! DBG PULSE:  if last model of last run then set LPULSE to LSAVPU
            if (model_iteration.eq.num_models(nk) .and. nk .eq. num_runs) then
                 pulsation_output_active = saved_pulse_output_flag
            end if

! JVS 02/11: Also allow pulse output at particular ages along the way
!
!  If the step is bracketing an age of interest, turn on output. This will
! for the step before and step after the age in AGEOUT. Once the info has
! been printed out, AGEOUT is set to the next age.
!
! Turn on calcad:
      if (acoustic_depth_output) then
            compute_acoustic_depth=.true.
      else
            compute_acoustic_depth = .false.
      endif
! If output has been turned on for a previous step, keep it on for the next
! step, but then turn it off.
      if (acoustic_depth_output) then
            if (ljwrt_placeholder) then
                  print*, 'LJWRT on'
                  pulsation_output_active = saved_pulse_output_flag
                  nao=nao+1
                  lclcd_placeholder =.true.
                  ljlast_placeholder =.false.
                  ljwrt_placeholder=.false.
            else if (.not.ljwrt_placeholder) then
                  lclcd_placeholder=.false.
            endif
! If this is the step before one of the ages of interest, print everything out.
! Also, save model structure.
            if (nao.lt.6) then
                  if (star%run%dage+timestep_yr/1.0D9-ageout_placeholder(nao) .le. 0.0D0 .and. &
                  star%run%dage+2.0D0*timestep_yr/1.0D9-ageout_placeholder(nao) .ge. 0.0D0 .and. .not. ljwrt_placeholder) then
                        print*, 'AGEOUT reached'
                        pulsation_output_active = saved_pulse_output_flag
                        lclcd_placeholder = .true.
                        ljlast_placeholder = .true.
                        ljwrt_placeholder=.true.
                  endif
            endif
       endif


! JVS end
!
! DBG PULSE:  if endage reached then set LPULSE to LSAVPU
! MHP 10/24 GENERALIZE CHECK
         if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (abs(target_end_age(nk)-star%run%dage*1.0D9-timestep_yr) .le. 1.0D0)) then
                 pulsation_output_active = saved_pulse_output_flag
! MHP 7/96 compute sound speed for solar model
                 star%run%sound_speed_output_active = .true.
            end if


!FD echo LSOUND
!        print*,'MAIN LSOUND = ',LSOUND
!FD end
            if (po_output_enabled) then
! MHP 8/25 changed to add file names as declared variables
             call pdist(prev_log_l,prev_log_teff,prev_age,path_length_sq,star%log10_luminosity,star%log_teff,model_iteration,job%pulse_atm_path, &
             job%pulse_env_path,job%pulse_mod_path)
          endif

! STARIN called here for timestep cutting
   15       if (model_diverged_flag) then
!              CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
             call starin(timestep_yr, delta_time, hydrogen_dt, &
                  trial_sign_flag, ikut_flag, istore_flag, &
                  model_diverged_flag, recompute_envelope_triangle, nk, &
                  dlnrho_dlnp, dlnrho_dlnt, total_angular_momentum, &
                  total_rotational_ke, convective_velocity, &
                  job%mixture_weights, ierr)
             if (ierr /= 0) return
             if ((star%omega(1) .eq. 0) .and. (rotation_active)) then
18               format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                        ' Initialize rotation rates or set LROT to', &
                        ' FALSE.')
                 print 18
                 ! 2026 (phase five, step B): configuration error returns to the
                 ! CLI wrapper (which stops) instead of stopping here.
                 ierr = 1
                 return
             endif
          endif
          punch_pending_flag = .true.

! skip this section if model not to be aged
! MHP 7/98
! need to add logic to permit resacling + time evolution for
! pre-main sequence models
            if (rescale_kind(nk).ne.2 .and. star%model_number.ge.0) then
               evolve_model_flag = .true.
            else if (star%model_number.ge.0 .and. star%log_temperature(1).lt.6.6D0) then
               evolve_model_flag = .true.
            else
               evolve_model_flag = .false.
            endif
            new_atmosphere_fit_needed = .false.
            if (evolve_model_flag) then
! ADD MASS LOSS CALCULATION
               call massloss(star%log10_luminosity,star%run%dage,delta_time,star%composition,star%log_density,star%specific_angular_momentum,star%log_pressure,star%log_radius, &
                             star%log_mass,star%enclosed_mass,star%shell_mass,star%log_total_mass,star%log_temperature,star%envelope_cz_bottom_index,recompute_envelope_triangle, &
                             star%num_zones,star%omega,star%total_mass_msun,star%log_teff,target_envelope_mass,new_atmosphere_fit_needed)
! STORE COMPOSITION MATRIX AT THE BEGINNING OF THE TIMESTEP.
               num_species = 11
               if (use_extended_composition) num_species=15
               do 33 i = 1,star%num_zones
                  do 32 j = 1,num_species
                     star%prev%old_composition(j,i) = star%composition(j,i)
   32             continue
   33          continue
               iteration_level=1
! mixed_zone_bounds_no_overshoot stays an ARGUMENT of mix (not read as
! star% inside it) because crrect passes its own local array there --
! storage deliberately separate from main's. main passes the star copy
! explicitly.
               call mix(delta_time, iteration_level, timestep_yr, &
                    star%core_cz_top_index, star%envelope_cz_bottom_index, &
                    star%mixed_zone_bounds_no_overshoot, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
             timestep_yr = delta_time/seconds_per_year
             star%run%dage = star%run%dage + 1.0D-9*timestep_yr
            endif
!***MHP 1/04 OPACITY TEST
!      IDT = 15
!      DO JJJ = 1,4
!         IDD(JJJ) = 5
!      END DO
!      XXX = 0.7D0
!      ZZZ = 0.02D0
!      Do JJJ = 1,2000
!         READ(75,*)TL,DL,XX,ZZ,OO
!         IF(TT.GT.1.0D9)STOP911
!         TL = LOG10(TT)
!         DL = LOG10(DD)
! DBG 12/95 GET OPACITY
!         CALL GETOPAC(DL, TL,XXX,ZZZ, O, OL, QOD, QOT)
!         DIFF = (O-OO)/O
!         RL = DL - 3*TL +18.0D0
!         WRITE(76,1554)TL,RL,XXX,ZZZ,O,OO,DIFF
! 1554    format(4f11.6,3e20.10)
!      END DO
!*** END TEST
! rezone new model, except rezoning not performed for He flash calculations
          if (.not.helium_flash_active) then
             call hpoint(istore_flag, reset_triangle, h_shell_zone_begin, &
                  has_h_shell, total_angular_momentum, &
                  total_rotational_ke, ierr)
             if (ierr /= 0) return
! STORE NEW CZ BASE
               light_burn%jcz = star%envelope_cz_bottom_index
            else
! save old model for PTIME
               do i=1, star%num_zones
                  star%prev%old_pressure(i) = star%log_pressure(i)
                  star%prev%old_temperature(i) = star%log_temperature(i)
                  star%prev%old_radius(i) = star%log_radius(i)
                  star%prev%old_luminosity(i) = star%luminosity_lsun(i)
                  star%prev%old_density(i) = star%log_density(i)
               end do
! JVS 04/14 Save Teff as well
               star%prev%old_teff = star%log_teff
!  JVS 05/25 Added model number to list of saved values
           star%prev%old_num_zones = star%num_zones

          endif
! store starting distribution of rotational kinetic energy.
            if (rotation_active) then
               do i = 1,star%num_zones
                  star%kinetic_energy_rot_old(i) = star%kinetic_energy_rot(i)
               end do
            endif
! changed for lithium burning with overshoot.
! store starting depth of C.Z. for light element burning.
            if (use_extended_composition) then
               light_burn%cz_base_radius_prev = 0.0D0
               envelope_cz_zone_prev = star%envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  light_burn%pressure_scale_height_start = alphae*exp(clndp*(star%log_pressure(star%envelope_cz_bottom_index)+2.0D0*star%log_radius(star%envelope_cz_bottom_index) &
                           -star%log_density(star%envelope_cz_bottom_index)-cgl-star%log_mass(star%envelope_cz_bottom_index)))
               else
                  light_burn%pressure_scale_height_start = 0.0D0
               endif
! find burning rates at the beginning of the time step.
               call lirate88(star%composition,star%log_density,star%log_temperature,star%num_zones,1)
            endif
! begin correction routines
! set flags for CRRECT
! CRRECT checks surface boundary conditions in every iteration
! if LNEW0 = T, new envelope triangle calculated the 1st iteration
! (i.e. old triangle ignored)
! LFINI = T if model has converged
! LARGE = T if model has diverged
          if (lnew0) recompute_envelope_triangle = .true.
            if (.not.evolve_model_flag) delta_time = -dabs(delta_time)
            fcorr = dabs(fcorr0) - fcorri
            iterations_done = 0
            model_diverged_flag = .false.
            converged = .false.
            if (.not.lnews .or. delta_time.le.0.0D0) then
               do 20 i = 1,star%num_zones
! zero entropy terms
                  star%log_temperature_delta(i) = 0.0D0
                  star%log_pressure_delta(i) = 0.0D0
                  star%run%temperature_entropy_term(i) = 0.0D0
                  star%run%pressure_entropy_term(i) = 0.0D0
                  star%run%luminosity_entropy_term(i) = 0.0D0
                  star%run%radius_entropy_term(i) = 0.0D0
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
 20            continue
            else
! use the rate of change in the previous model to estimate the new
! run of structure variables.
               do 30 i = 1,star%num_zones
                  delta_temp_step = star%run%temperature_entropy_term(i)*delta_time
                  delta_pressure_step = star%run%pressure_entropy_term(i)*delta_time
                  delta_lum_step = star%luminosity_lsun(i)*star%run%luminosity_entropy_term(i)*delta_time
                  delta_radius_step = star%run%radius_entropy_term(i)*delta_time
                  star%log_temperature_delta(i) = delta_temp_step
                  star%log_pressure_delta(i) = delta_pressure_step
                  star%log_temperature(i) = star%log_temperature(i) + delta_temp_step
                  star%log_pressure(i) = star%log_pressure(i) + delta_pressure_step
                  star%luminosity_lsun(i) = star%luminosity_lsun(i) + delta_lum_step
                  star%log_radius(i) = star%log_radius(i) + delta_radius_step
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
 30            continue
            endif

! FIRST LEVEL OF ITERATIONS
! USE ENVELOPE TRIANGLE OF THE PREVIOUS MODEL;
! FOR THE FIRST MODEL OF A RUN,THE TRIANGLE IS GENERATED HERE.
            max_iterations = niter1
            recompute_surface_bc = .false.
! CALL TO CRRECT - ADDED ITERATION LEVEL
            iteration_level = 1
            call crrect(delta_time, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, dlnrho_dlnt, &
                 dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
! SECOND LEVEL OF ITERATIONS
! CHECK ENVELOPE TRIANGLE BEFORE ITERATING FOR SOLUTION
            if (model_diverged_flag) goto 15
            recompute_surface_bc = .true.
            max_iterations = niter2
            iteration_level = 2
            call crrect(delta_time, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, dlnrho_dlnt, &
                 dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
            if (model_diverged_flag) goto 15
! 7/91 STORE CHANGES IN THE STRUCTURE. THESE CHANGES ARE USED TO GET AN
! IMPROVED FIRST GUESS AT THE STRUCTURE FOR THE NEXT MODEL IF LNEWS=T.
            if (delta_time.gt.0.0D0) then
               do 27 ii = 1,star%num_zones
                  star%run%temperature_entropy_term(ii)=star%log_temperature_delta(ii)/delta_time
                  star%run%pressure_entropy_term(ii)=star%log_pressure_delta(ii)/delta_time
                  star%run%luminosity_entropy_term(ii)=2.0D0*(star%luminosity_lsun(ii)-star%prev%old_luminosity(ii))/(star%luminosity_lsun(ii)+star%prev%old_luminosity(ii))/delta_time
                  star%run%radius_entropy_term(ii)=(star%log_radius(ii)-star%prev%old_radius(ii))/delta_time
 27            continue
            endif
! THIRD LEVEL OF ITERATIONS
            recompute_surface_bc = .false.
            max_iterations = niter3
            iteration_level = 3
            call crrect(delta_time, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, dlnrho_dlnt, &
                 dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
            if (model_diverged_flag) goto 15
            if (.not.rotation_active) then
               itdif1 = 1
            endif
! MHP 05/02
! IF THE CODE IS ITERATING BETWEEN THE STRUCTURE AND ROTATION
! SOLUTIONS, ENSURE THAT THE START-OF-TIMESTEP QUANTITIES
! HCOMPP (COMPOSITION) AND HJMSAV (ANGULAR MOMENTUM) ARE ONLY
! OVERWRITTEN ON THE LAST RUN THROUGH.
            if (itdif1.gt.1) then
               do i = 1,star%num_zones
                  star%run%orig_specific_angular_momentum(i) = star%specific_angular_momentum(i)
                  do j = 1,15
                     star%run%orig_composition(j,i) = star%prev%old_composition(j,i)
                  end do
               end do
            endif
            do itrot = 1, itdif1
! MHP 05/02 RESTORE ORIGINAL "START OF TIMESTEP"
! VALUES FOR THE COMPOSITION MATRIX
               if (itrot.gt.1) then
                  do i = 1,star%num_zones
                     do j = 1,15
                        star%prev%old_composition(j,i) = star%run%orig_composition(j,i)
                     end do
                  end do
               endif
! 7/91 THE FOURTH LEVEL OF ITERATION REPEATS THE ITERATION BETWEEN THE
! MIXING AND THE STRUCTURE VARIABLES.  IT SHOULD NOT BE USED FOR MODELS
! WHERE SEMI-CONVECTION IS IMPORTANT (ITERATING BETWEEN THE BURNING AND
! STRUCTURE GENERATES OSCILLATIONS). IT SHOULD BE USED FOR HIGH-PRECISION
! WORK (E.G. SOLAR MODELS).
! Surface boundary conditions checked again since we've changed the
! star%composition (and hence the structure) of the model in ITLVL=3
! (to be implemented when I know the rest of it works!)
            max_iterations = niter4
            recompute_surface_bc=.false.
            iteration_level = 4
            call crrect(delta_time, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, dlnrho_dlnt, &
                 dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
!  25         CONTINUE
            if (.not.converged) then
! MODEL FAILED TO CONVERGE WITHIN(NITER1+NITER2+NITER3+NITER4)ITERATIONS
               model_diverged_flag = .true.
               goto 15
            endif

! MODEL HAS CONVERGED
! ENSURE CONVECTION ZONES ARE FULLY MIXED.
! MHP 02/12 INFER CONVECTIVE OVERTURN TIMESCALE (USED IN MDOT)
! JVS 02/12 CALL MIXCZ(HCOMP,HS2,LC,M)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
! C G Somers 6/14, SET IMIX = .TRUE. SO THE CORRECT GRADS ARE USED.
!       IMIX = .TRUE.
!       CALL MIXCZ(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,IMIX)
! G Somers 6/14, SET LIMIX = .TRUE. SO THE CORRECT GRADS ARE USED.
      use_correct_gradients = .true.
!       CALL MIXCZ(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,LIMIX)  ! KC 2025-05-31
      call mixcz(star%composition,star%shell_mass,star%convective_flag,star%num_zones)
! G Somers END

! MHP 9/94 STORE TOTAL AGE IN SAGE
            disk_lifetime = star%run%dage
            if (rotation_active) then
! RESTORE ORIGINAL START OF TIMESTEP VALUES
! TO THE ANGULAR MOMENTUM DISTRIBUTION
               if (itrot.gt.1) then
                  do i = 1,star%num_zones
                     star%specific_angular_momentum(i) = star%run%orig_specific_angular_momentum(i)
                  end do
               endif
! MHP 9/94 ADDED FLAG TO TURN ON ROTATION OUTPUT WHEN END OF KIND
! CARD REACHED.
! MHP 10/24 GENERALIZE CHECK
         if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (abs(target_end_age(nk)-star%run%dage*1.0D9-timestep_yr) .le. 1.0D0)) then
!               IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)THEN
                  star%run%lprt0_placeholder = .true.
               else
                  star%run%lprt0_placeholder = .false.
               endif
! FIND THE NEW RUN OF OMEGA
! JENV0 ADDED TO SR CALL.
               wind_loss_active = ljdot0
               call getw(delta_time, max_domega_frac, wind_loss_active, &
                    envelope_cz_zone_prev, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
               call fpft(star%log_density,star%log_radius,star%log_mass,star%num_zones,star%omega,star%eta_squared,star%pressure_rotation_factor,star%temperature_rotation_factor,star%mean_gravity,star%mean_radius)
            endif
            end do
! LOCATE THE HYDROGEN-BURNING SHELL AND THE BOUNDARIES OF THE CENTRAL
! AND SURFACE CONVECTION ZONES (IF APPLICABLE).
       call findsh(star%composition,star%luminosity_lsun,star%convective_flag,star%num_zones, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,h_shell_zone_begin,h_shell_end_index,h_shell_midpoint_zone, &
              has_h_shell)
! PERFORM LIGHT ELEMENT BURNING
         if (use_extended_composition .and. star%model_number.ge.0 .and. delta_time.gt.0.0D0) then
! ONLY FOR MODELS WITHOUT ROTATION, OR WITHOUT ROTATIONAL MIXING.
            if (.not.rotation_active .or. .not.instability_transport_active) then
! FIND CONVECTION ZONE DEPTH AT THE END OF THE TIME STEP.
               call convec(star%composition,star%log_density,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,star%num_zones,star%radiative_zone_bounds,star%mixed_zone_bounds, &
                            star%mixed_zone_bounds_no_overshoot,star%core_cz_top_index,star%envelope_cz_bottom_index,num_radiative_zones,num_mixed_zones,num_mixed_zones_no_overshoot)
! CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
               envelope_cz_zone_end = star%envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  light_burn%pressure_scale_height_end = alphae*exp(clndp*(star%log_pressure(star%envelope_cz_bottom_index)+2.0D0*star%log_radius(star%envelope_cz_bottom_index) &
                           -star%log_density(star%envelope_cz_bottom_index)-cgl-star%log_mass(star%envelope_cz_bottom_index)))
               else
                  light_burn%pressure_scale_height_end = 0.0D0
               endif
! FIND BURNING RATES AT THE END OF THE TIME STEP.
               call lirate88(star%composition,star%log_density,star%log_temperature,star%num_zones,2)
!                CALL LIBURN(DELTS,HCOMP,HD,HR,HS1,HS2,HT,JENV1,JENV0,M)  ! KC 2025-05-31
               call liburn(delta_time,star%composition,star%log_radius,star%enclosed_mass,star%shell_mass,star%log_temperature,envelope_cz_zone_end,envelope_cz_zone_prev,star%num_zones)
            endif
         endif
! MHP 07/02 RESTORE PRIOR FITTING POINT IF MASS ACCRETION IS BEING
! INCLUDED
         if (new_atmosphere_fit_needed) then
            call getnewenv(target_envelope_mass,star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%enclosed_mass,star%shell_mass, &
!     *                     HSTOT,HT,LC,ETA2,HG,HI,HJM,QIW,R0,  ! KC 2025-05-31
                            star%log_total_mass,star%log_temperature,star%convective_flag,star%eta_squared,star%moment_of_inertia,star%specific_angular_momentum,star%qiw,star%mean_radius, &
                            star%kinetic_energy_rot,star%log10_luminosity,total_angular_momentum,total_rotational_ke,star%log_teff,star%num_zones,recompute_envelope_triangle)
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
            call fpft(star%log_density,star%log_radius,star%log_mass,star%num_zones,star%omega,star%eta_squared,star%pressure_rotation_factor,star%temperature_rotation_factor,star%mean_gravity,star%mean_radius)
            new_atmosphere_fit_needed = .false.
         endif
! DETERMINE TIMESTEP FOR NEXT MODEL
! HTIMER ALSO LOCATES THE H-BURNING SHELL
! JVS 04/14 added teffl to passed htimer variables
       delta_time = dabs(delta_time)
       delta_time_saved = delta_time
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *        JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(delta_time,hydrogen_dt,star%num_zones,star%log_density,star%luminosity_lsun,star%enclosed_mass,star%shell_mass,star%log_temperature,star%composition,star%core_cz_top_index, &
              h_shell_midpoint_zone,star%luminosity_breakdown,star%run%dage,timestep_yr,nk,star%log_pressure,star%log_radius,star%omega, &
              max_domega_frac,h_shell_zone_begin,star%log_teff)
! IF EVOLVING TO A GIVEN AGE AND KIND CARD IS DONE, AVOID ZEROING OUT
! TIMESTEP WRITTEN TO MODEL (AS THIS MAKES CONTINUING A SEQUENCE AWKWARD.)
!     INSTEAD WRITE THE PREVIOUS MODEL TIMESTEP TO MODEL.
! ONLY IF A FIXED END AGE IS USED, NOT FOR OTHER STOPS
       if (end_age_stop_active(nk) .and. target_end_age(nk).gt.0.0D0) then
          if (target_end_age(nk)-star%run%dage*1.0D9.le.1.0D0) then
             delta_time = max(delta_time_saved,1.0D-3*star%run%dage*seconds_per_year)
             timestep_yr = delta_time/seconds_per_year
          else
             delta_time_saved = delta_time
          endif
       else
          delta_time_saved = delta_time
       endif
       if (rescale_kind(nk).ne.2) star%model_number = star%model_number+1
! 2026 (phase four, step 5): compute the output diagnostics in the
! star layer (fills star%run%*, star%luminosity_breakdown
! renormalization, turnover% via gettau); wrtout below only reads.
       call update_output_diagnostics(ierr)
       if (ierr /= 0) return
! WRTOUT IS THE OUTPUT DRIVER ROUTINE
       call wrtout(timestep_yr, log_gravity, has_h_shell, &
            h_shell_zone_begin, h_shell_midpoint_zone, h_shell_end_index, &
            trial_sign_flag, punch_pending_flag, total_angular_momentum, &
            total_rotational_ke)

! MHP 10/24 GENERALIZED STOP CONDITIONS
!     IF EVOLVING TO A GIVEN AGE AND AGE IS REACHED, KIND CARD IS DONE
!       IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)GOTO 110
       if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (target_end_age(nk)-star%run%dage*1.0D9).le.1.0D0) goto 110
! MHP 10/24 CHECK ALL STOP CONDITIONS, EXIT IF ANY SATISFIED
         end_kind_flag = .false.
         if (end_age_stop_active(nk).and.central_deuterium_stop(nk).gt.0.0D0 .and. &
              star%composition(12,1).lt.central_deuterium_stop(nk)) then
            write(*,104)star%composition(12,1),central_deuterium_stop(nk)
 104        format('CENTRAL D ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_hydrogen_stop(nk).gt.0.0D0 .and. &
              star%composition(1,1).lt.central_hydrogen_stop(nk)) then
            write(*,105)star%composition(1,1),central_hydrogen_stop(nk)
 105        format('CENTRAL X ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_helium_stop(nk).gt.0.0D0 .and. &
              star%composition(2,1).lt.central_helium_stop(nk)) then
            write(*,106)star%composition(2,1),central_helium_stop(nk)
 106        format('CENTRAL Y ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
! IF EXITING, SET I/O FLAGS PROPERLY AND EXIT LOOP
         if (end_kind_flag) then
            pulsation_output_active = saved_pulse_output_flag
            star%run%sound_speed_output_active = .true.
            star%run%lprt0_placeholder = .true.
            goto 110
         endif
! TEST IF MODEL IS NEAR DESIRED Teff AND L. IF NOT RESCALE AND TRY AGAIN.
         if (calibrate_star_flag .and. .not. star_found_flag) then
            if (mod(nk,2).eq.0) then
             if (model_iteration.eq.1) then
                teff_kelvin_unused = 10.0D0**star%log_teff
             else
                call chkscal(star%log10_luminosity, star%log_teff, star%run%dage, nk)
                if (just_passed_target_radius_flag) goto 200
             end if
          endif
       endif

! END OF RUN
  100    continue

! G Somers 11/14, CHANGE CALL TO PUTSTORE INSTEAD OF WRTLST.
! STORE LAST MODEL IN ISTOR IF LSTORE, LSTPCH, AND LPUNCH ARE .TRUE.
  110    if (lstore.and.lstpch.and.punch_pending_flag) then
          call putstore(star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index,star%envelope_cz_bottom_index,star%model_number,star%num_zones,star%total_mass_msun,star%log_teff,star%log10_luminosity,star%log_total_mass, &
                 star%run%dage,timestep_yr,star%omega,star%enclosed_mass,star%eta_squared,star%mean_radius,star%pressure_rotation_factor,star%temperature_rotation_factor,star%specific_angular_momentum,star%moment_of_inertia)
            punch_pending_flag = .false.
       endif
! 110  CONTINUE
! G Somers END

! MHP 1/93 CHECK AUTOMATIC CALIBRATATION OF SOLAR MODEL.
!c MHP 5/96 changed solar calibration to perform solar models in 3 kind cards
         if (calibrate_solar_model) then
! JVS Turn off calcad - speeds things up
            compute_acoustic_depth=.false.
            if (mod(nk,3).eq.0) then
               log_r_rsun = 0.5D0*(star%log10_luminosity+log10_solar_luminosity-c4pil-csigl-4.0D0*star%log_teff)-log10_solar_radius
! MHP 06/13 Add solar Z/X to observables
               current_zx = star%composition(3,star%num_zones)/star%composition(1,star%num_zones)
               call chkcal(star%log10_luminosity,log_r_rsun,nk,current_zx)
!               CALL CHKCAL(BL,RLL,NK)
               use_structure_dt_limits = saved_use_structure_dt_limits  ! Restore LPTIME to original value for next cycle
               atm_choice  = saved_atm_choice    ! Restore KTTAU to original value for next cycle
               if (star%run%solar_calibration_active) then
                  go to 250
               else
!c MHP 8/96 added counter for # of runs needed for calibration
                  convergence_iterations = convergence_iterations + 1
! MHP 6/97 STOP AFTER 10 ATTEMPTS AT CALIBRATION
!                  IF(ICONV.GE.11) GOTO 250
                  if (convergence_iterations.ge.15) goto 250
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
         if (calibrate_star_flag .and. star_found_flag.and.(mod(nk,2).eq.0)) goto 250

! END RUN LOOP
 200  continue
! EXIT RUN LOOP
 250  continue

! FOR MONTE CARLO, REWIND OUTPUT FILES AND WRITE OUT SNU FLUXES AND
! MODEL PARAMTERS TO AN OUTPUT FILE.
! RUN FAILED TO CONVERGE.  WRITE FINAL INFO WITH WARNING NOTE.
      if (lmonte .and. convergence_iterations.ge.11 .and. .not.star%run%solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)
         write(neutrino_unit,1525)star%log10_luminosity,log_r_rsun
 1525    format(5X,'DID NOT CONVERGE WITHIN 10 ATTEMPTS L,R',2F10.6)
! MONTE CARLO #, CONVERGED MIXING LENGTH AND INITIAL H, SURFACE X,
! SURFACE Z, Z/X, CENTRAL X, CENTRAL Z
         write(neutrino_unit,1519) monte_carlo_run_number,mixing_length_array(nk),rescale_params(2,nk-2),star%composition(1,star%num_zones), &
              star%composition(3,star%num_zones),surface_z_over_x,star%composition(1,1),star%composition(3,1)
 1519    format(1X,I5,3F10.6,4E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,star%run%dlum_dx,star%run%drad_dx,star%run%dlum_dalpha,star%run%drad_dalpha
! SUMMARY OF STRUCTURE : TC, RHOC, PC
         write(neutrino_unit, 1517)star%run%central_log10_temperature,star%run%central_log10_pressure,star%run%central_log10_density, &
              star%composition(1,1),star%composition(3,1)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) flux_diag%cl37_snu_rate,flux_diag%ga71_snu_rate,(flux_diag%neutrino_flux_total(i),i=1,8)
!          CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *        DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *        JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *        CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
         call wrtmonte(star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,star%num_zones,star%run%dage, &
              timestep_yr,star%total_mass_msun,star%log_teff,star%log10_luminosity, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
              star%envelope_fit_coeffs,trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
      else if (calibrate_solar_model .and. lsnu .and. star%run%solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)

         surface_z_over_x = star%composition(3,star%num_zones)/star%composition(1,star%num_zones)
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
         write(neutrino_unit, 1516) flux_diag%cl37_snu_rate,flux_diag%ga71_snu_rate,(flux_diag%neutrino_flux_total(i),i=1,10)
 1516    format(1X,2F8.3,1P10E10.3)
! SUMMARY OF STRUCTURE : TC, RHOC, PC, XC, ZC (ADD MU C)
         central_temperature_mk = 10.0D0**(star%run%central_log10_temperature-6.0D0)
         central_pressure_scaled = 10.0D0**(star%run%central_log10_pressure-17.0D0)
         central_density_linear = 10.0D0**star%run%central_log10_density
         write(neutrino_unit, 1517)central_temperature_mk,central_density_linear,central_pressure_scaled,star%composition(1,1),star%composition(3,1)
 1517    format(1X,F7.3,F7.2,F6.3,2F8.5)
! INITIAL ALPHA,Y,Z,ALPHA; FINAL R, L
         initial_helium_fraction = 1.0D0 - rescale_params(2,nk-2) - rescale_params(3,nk-2)
         initial_metal_fraction = rescale_params(3,nk-2)
         write(neutrino_unit,1521)mixing_length_array(nk),initial_helium_fraction,initial_metal_fraction,star%log10_luminosity,log_r_rsun
 1521    format(F7.4,2F8.5,1P2E10.3)
! CZ DEPTH (R,M), SURFACE Y, Z, Z/X (ADD T CZ BASE, RHO CZ BASE)
         write(neutrino_unit,1522)star%run%envelope_radius,star%run%envelope_mass,star%composition(2,star%num_zones),star%composition(3,star%num_zones),surface_z_over_x
 1522    format(F8.5,F9.6,2F8.5,F9.6)
! ENERGY GENERATION FRACTIONS PP I,II,III,CNO,EGRAV
         write(neutrino_unit,1523)(star%luminosity_breakdown(j),j=1,4),star%luminosity_breakdown(7)
 1523    format(1P5E10.3)
         if (lmonte) then
!             CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *           DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *           JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *           CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
            call wrtmonte(star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,star%num_zones,star%run%dage, &
                 timestep_yr,star%total_mass_msun,star%log_teff,star%log10_luminosity, &
                 star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
         endif
      endif
 500  end do

! 2026 (phase five, step B): the normal end-of-job stop became this
! clean return (ierr stays 0); the CLI wrapper simply ends.
      return
end subroutine run_yrec
