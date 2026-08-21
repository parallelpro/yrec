!----------------------------------------------------------------------
! starin
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original starin.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!       SUBROUTINE STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,FP,  ! KC 2025-05-31
!      * FT,FTRI,HCOMP,HD,HI,HJM,HKEROT,HL,HP,HR,HS,HS1,HS2,HSTOT,HT,IKUT,
!      * ISTORE,JCORE,JENV,LARGE,LC,LNEW,M,MODEL,NK,OMEGA,PS,QDP,QDT,QIW,  ! KC 2025-05-31
!      * R0,RS,SJTOT,SKEROT,SMASS,TEFFL,TLUMX,TRIL,TRIT,TS,VEL,HG,V)
!
! This is YREC's initial-model setup/read driver, called once at the
! start of every model. If a fresh starting model is requested
! (first_call_flag(run_index) true), it reads the model file (YREC7 or
! MODEL2 format, detected from the file's 4-character keyword) via
! getyrec7/getmodel2, optionally extends the innermost shell inward
! (common/core/), rescales it (rscale), optionally changes the
! envelope fitting point mass (calling envint/meqos/eqstat/getopac/
! tpgrad to keep the new last interior shell's radiative/convective
! flag and density consistent), sets up the rotation curve (fpft/
! momi) if rotation is active, sets up the surface mixture and
! opacity tables (surfopac/setscv), and finally calls physic/ovrot/
! gettau to get the initial temperature-gradient/opacity/convective-
! zone structure. If instead the in-memory model is reused
! (first_call_flag(run_index) false), only the envelope-mass-fitting-
! point rescale (rscale) and the final physic/ovrot/gettau calls are
! done.
!
! CROSS-CALLEE NAMING NOTE: several dummy arguments here are threaded
! into more than one already-converted callee, and those callees do
! not always agree on a name for the same physical data (this file's
! own dummy names were free to choose, per the project's incremental
! conversion order, but the callees' names were fixed by earlier
! batches). Judgment calls made below, all verified against the
! actual physics/usage in this file:
!   - HL is the linear luminosity (L/Lsun) -- STOTAL/SENV bookkeeping
!     and DLOG10(HL(M))-style usage elsewhere in the codebase confirm
!     this. Named luminosity_lsun (matches crrect.f90's slot name for
!     the same physical quantity). getyrec7.f90/getmodel2.f90/
!     physic.f90's slot name for the same array is "log_luminosity",
!     which is a misnomer inherited from those files' own earlier
!     conversions; out of scope to fix here. The scalar per-point
!     value of this same array (originally B) is named
!     shell_luminosity_lsun to avoid clashing with the array name.
!   - HS1 is named enclosed_mass (matches crrect.f90's slot name for
!     the same array, and gettau.f90's own parameter name for it).
!   - HS2 is named shell_mass (matches crrect.f90/rscale.f90/
!     momi.f90).
!   - R0 is named mean_radius (matches momi.f90, which computes it,
!     and crrect.f90); fpft.f90's own slot name for the same array is
!     "r0".
!   - QIW is named di_domega (matches momi.f90, which computes it);
!     crrect.f90 instead keeps this slot named "qiw" verbatim.
!   - HG is named mean_gravity (matches fpft.f90, which computes it);
!     physic.f90's own slot name for the same array is "hg".
!   - M (number of mesh points) is named num_shells (matches
!     getyrec7.f90/getmodel2.f90); physic.f90/ovrot.f90/gettau.f90/
!     rscale.f90 call the same count "num_zones", momi.f90 calls it
!     "zone_end", crrect.f90 calls it "num_points".
!   - SMASS is named total_mass_msun (matches getyrec7.f90/
!     getmodel2.f90); rscale.f90's own slot name for the same value is
!     "star_mass".
!   - FP/FT (rotational P/T correction factors) are named
!     pressure_rotation_factor/temperature_rotation_factor (matches
!     crrect.f90); the per-point scalars used in the single-shell
!     EQSTAT/TPGRAD physics block below (originally FPL/FTL) are named
!     point_pressure_rotation_factor/point_temperature_rotation_factor
!     to avoid clashing with the array names (physic.f90 instead
!     avoids this clash by keeping its own array dummies short, "fp"/
!     "ft", and giving the fuller names to its per-point locals).
!   - Per-point scalar physics locals below otherwise follow
!     physic.f90's own naming exactly (log10_temperature/temperature/
!     log10_pressure/pressure/log10_density/density/log10_radius/
!     log10_mass/opacity/beta/.../is_convective etc.), since this
!     block is explicitly "stolen from PHYSIC" per the original
!     source's own comment. Note physic.f90 uses a "log10_" prefix for
!     these per-point scalars specifically so they do not collide with
!     its "log_"-prefixed array dummy names; the same convention is
!     used here for the same reason (this file's own array dummies are
!     named log_pressure/log_temperature/log_density/log_radius/
!     log_mass).
!   - common/core/ (LCORE/MCORE/FCORE) and common/newmx/ (the CNO-
!     mixture/isotope-ratio controls) are referenced only in
!     parmin.f90 among already-converted files; parmin.f90 keeps their
!     members at cryptic lowercased spelling because many are
!     NAMELIST-exposed. Per this batch's own conversion instructions,
!     this file is not bound by that naming and gives them descriptive
!     names instead (see the common block declarations below for the
!     mapping back to parmin.f90's spelling).
!   - common/alexmix/ (XALEX/ZALEX) and common/vnewcb/ (VNEW) are
!     likewise only otherwise established in parmin.f90.  vnew is kept
!     as-is (already a plain, adequately descriptive name, and it is
!     the name used directly in this file's own pre-existing
!     comments); alexmix's members are given descriptive names since
!     they are unused placeholders in this file.
subroutine starin(log10_luminosity, envelope_fit_coeffs, age_gyr, &
     timestep_yr, delta_time, delta_time_abs, eta_squared, &
     pressure_rotation_factor, temperature_rotation_factor, trial_sign_flag, &
     composition, log_density, moment_of_inertia, specific_angular_momentum, &
     kinetic_energy_rot, luminosity_lsun, log_pressure, log_radius, &
     log_mass, enclosed_mass, shell_mass, log_total_mass, log_temperature, &
     ikut_flag, istore_flag, envelope_zone_index, model_failed_flag, &
     convective_flag, envelope_recomputed_flag, num_shells, model_number, &
     run_index, omega, fit_point_pressure, dlnrho_dlnp, dlnrho_dlnt, &
     di_domega, mean_radius, fit_point_radius, total_angular_momentum, &
     total_rotational_ke, total_mass_msun, log_teff, luminosity_breakdown, &
     trial_log_luminosity, trial_log_temperature, fit_point_temperature, &
     convective_velocity, mean_gravity, species_mix_weights)

      use envstruct_lib
      use envelope_comp_lib
      use turnover_lib
      use scrtch_lib
      use oldmod_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000
      integer, parameter :: nts = 63, nps = 76

      double precision, intent(inout) :: log10_luminosity
      double precision, intent(inout) :: envelope_fit_coeffs(9)
      double precision, intent(inout) :: age_gyr, timestep_yr
      double precision, intent(inout) :: delta_time, delta_time_abs
      double precision, intent(inout) :: eta_squared(json)
      double precision, intent(inout) :: pressure_rotation_factor(json), &
           temperature_rotation_factor(json)
      double precision, intent(inout) :: trial_sign_flag
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: log_density(json)
      double precision, intent(inout) :: moment_of_inertia(json), &
           specific_angular_momentum(json)
      double precision, intent(inout) :: kinetic_energy_rot(json)
      double precision, intent(inout) :: luminosity_lsun(json)
      double precision, intent(inout) :: log_pressure(json), &
           log_radius(json), log_mass(json)
      double precision, intent(out) :: enclosed_mass(json), shell_mass(json)
      double precision, intent(inout) :: log_total_mass
      double precision, intent(inout) :: log_temperature(json)
      integer, intent(out) :: ikut_flag
      integer, intent(inout) :: istore_flag
      integer, intent(in) :: envelope_zone_index
      logical, intent(in) :: model_failed_flag
      logical, intent(inout) :: convective_flag(json)
      logical, intent(inout) :: envelope_recomputed_flag
      integer, intent(inout) :: num_shells, model_number
      integer, intent(in) :: run_index
      double precision, intent(inout) :: omega(json)
      double precision, intent(inout) :: fit_point_pressure(3)
      double precision, intent(out) :: dlnrho_dlnp, dlnrho_dlnt
      double precision, intent(inout) :: di_domega(json)
      double precision, intent(inout) :: mean_radius(json)
      double precision, intent(inout) :: fit_point_radius(3)
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke
      double precision, intent(inout) :: total_mass_msun
      double precision, intent(inout) :: log_teff
      double precision, intent(inout) :: luminosity_breakdown(8)
      double precision, intent(inout) :: trial_log_luminosity(3), &
           trial_log_temperature(3)
      double precision, intent(inout) :: fit_point_temperature(3)
      double precision, intent(out) :: convective_velocity
      double precision, intent(inout) :: mean_gravity(json)
      double precision, intent(inout) :: species_mix_weights(12)

! DBGLAOL
      integer*4 :: katm, kenv, saha_state
! common/nwlaol/: olaol/oxa/ot/orho/tollaol (the LAOL pure-Z table
! data itself) are placeholders for layout; llaol is used here. Naming
! matches getopac.f90/envint.f90.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez
!      CHARACTER*256 OPECALEX(7)
      character(len=4) :: format_tag
      character(len=6) :: eos_code
      character(len=4) :: atm_code, alok_code, hik_code, compmix_code
! common/comp2/: both used here. Naming matches checkc.f90/physic.f90.
      double precision :: envelope_helium_fraction, envelope_he3_fraction
      common/comp2/ envelope_helium_fraction, envelope_he3_fraction
! common/envprt/: not used in this file's logic; layout placeholder.
! Naming matches envint.f90/getnewenv.f90/qenv.f90.
      double precision :: current_log10_pressure, current_log10_temperature, &
           current_log10_radius, current_log10_mass, current_log10_density, &
           current_opacity, current_beta, current_gradients(3), &
           current_ion_fraction(3), current_velocity
      common/envprt/ current_log10_pressure, current_log10_temperature, &
           current_log10_radius, current_log10_mass, current_log10_density, &
           current_opacity, current_beta, current_gradients, &
           current_ion_fraction, current_velocity
! common/oldrot/: only old_omega is used here (the rest are unused
! placeholders). Naming matches hpoint.f90/midmod.f90.
      double precision :: old_omega(json), old_specific_angular_momentum(json), &
           old_moment_of_inertia(json), old_hg(json), old_mean_radius(json), &
           old_eta_squared(json)
      common/oldrot/ old_omega, old_specific_angular_momentum, &
           old_moment_of_inertia, old_hg, old_mean_radius, old_eta_squared
! OPACITY COMMON BLOCKS - modified 3/09
! common/newopac/: not used in this file's logic; layout placeholder.
! Naming matches getopac.f90/surfopac.f90.
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max
      logical :: use_alex06_tables, use_laol89_tables, use_opal92_tables, &
           use_opal95_tables, use_kurucz90_tables, use_alex95_tables, &
           use_two_z_tables
      common /newopac/ laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max, use_alex06_tables, use_laol89_tables, &
           use_opal92_tables, use_opal95_tables, use_kurucz90_tables, &
           use_alex95_tables, use_two_z_tables
! MHP 8/25 Removed all character strings from common blocks
! common/alexo/: not used in this file's logic; layout placeholder.
! Naming matches alxtbl.f90 (parmin.f90 keeps this slot's cryptic
! original name "ialxo").
      integer :: alex95_table_unit
      common /alexo/ alex95_table_unit
! common/alexmix/: not used in this file's logic; layout placeholder.
! Only established elsewhere in parmin.f90 (as xalex/zalex); this file
! is not bound by that and uses descriptive names instead.
      double precision :: alex_mixture_x, alex_mixture_z
      common /alexmix/ alex_mixture_x, alex_mixture_z
! MHP  5/97 ADDED COMMON BLOCK FOR SCV EOS TABLES
! common/scveos/: only use_scv_eos is used here (the rest are unused
! placeholders). Naming matches setscv.f90/eqstat2.f90, except idtt:
! setscv.f90 names this slot "idt", which would collide here with
! common/optab/'s idt (also declared in this file); kept as the
! original cryptic spelling (lowercased) to avoid that collision.
      double precision :: tlogx(nts), tablex(nts,nps,12), &
           tabley(nts,nps,12), smix(nts,nps), tablez(nts,nps,13), &
           tablenv(nts,nps,12)
      integer :: nptsx(nts), idtt, idp
      logical :: use_scv_eos
      common/scveos/ tlogx, tablex, tabley, smix, tablez, tablenv, nptsx, &
           use_scv_eos, idtt, idp
! LLP  3/19/03 Add COMMON block /I2O/ for info directly transferred from
!      input to output model - starting with a code for th initial model
!      compostion (COMPMIX)
! common/i2o/: used here (compmix_code is passed to getyrec7/
! getmodel2). Naming matches this file's own atm_code/eos_code/
! hik_code/alok_code convention.
      common /i2o/ compmix_code
      double precision :: atomic_weight(12)
      data atomic_weight /23.0d0,26.99d0,24.32d0,55.86d0,28.1d0,12.015d0, &
           1.008d0,16.0d0,14.01d0,39.96d0,20.19d0,4.004d0/

! MHP 10/24 ENSURE THAT ONLY HOMOGENEOUS MODELS HAVE THE MIXTURE ALTERED
      double precision :: reference_composition(15)

! --- locals ---
      integer :: iread
      double precision :: mixing_length0
      logical :: mixing_length_matches
! NOTE: LEXCP0 in the original source is never assigned before use
! (line "WRITE(ISHORT,1040) CMIXL,CMIXL0,LEXCOM,LEXCP0" below); it
! appears to be a pre-existing bug (likely meant to read
! use_extended_composition0). Preserved exactly as in the original.
      logical :: lexcp0
      logical :: use_extended_composition0
      double precision :: fraction_diff
      integer :: i, j, k, kk
      double precision :: total_carbon_cno_fraction, &
           total_nitrogen_cno_fraction, total_oxygen_cno_fraction
      double precision :: carbon_scale_ratio, nitrogen_scale_ratio, &
           oxygen_scale_ratio
      double precision :: sum_c12_c13, sum_o16_o18
      integer :: num_shells_extended
      double precision :: core_shell_spacing, central_log_density, &
           actual_gradient, central_shell_luminosity, density_estimate_offset
      integer :: first_original_shell
      double precision :: trial_log_pressure, temp_scratch
      logical :: want_derivatives, local_conductive_opacity_flag, &
           in_atmosphere
      double precision :: hydrogen_fraction, metal_fraction
      double precision :: log10_pressure, pressure, log10_temperature, &
           temperature, log10_density, density
      double precision :: log10_radius, log10_mass, shell_luminosity_lsun
      double precision :: point_pressure_rotation_factor, &
           point_temperature_rotation_factor
      double precision :: beta, beta_inverse, beta14, ion_fraction(3), &
           specific_gas_constant, ion_mean_weight_inverse, &
           electron_mean_weight_inverse, electron_degeneracy_parameter, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      double precision :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      double precision :: radiative_gradient, dgrad_dt_component, &
           dgrad_dp_component, dgrad_dr_component
      logical :: is_convective
      integer :: old_last_shell
      double precision :: saved_env_step_max, saved_env_step_min, &
           saved_env_step_begin
      logical :: save_boundary_flag, print_flag, pulse_print_flag
      double precision :: log10_gravity
      integer :: vertex_index
      double precision :: log10_pressure_limit
      integer :: envint_unused_flag
      double precision :: spot_adjusted_log_teff
      double precision :: envint_dummy1(4), envint_dummy2(3), &
           envint_dummy3(3), envint_dummy4(3)
      double precision :: pressure_offset, density_offset, &
           temperature_offset, radius_offset
      integer :: env_point_index
      double precision :: lower_mass_coord, target_mass_coord, &
           upper_mass_coord, envelope_interp_fraction
      double precision :: old_senv, target_log_mass_at_fit
      double precision :: interior_interp_fraction
      integer :: num_species
      double precision :: prev_mass, curr_mass, next_mass
      double precision :: angular_momentum_sum, rotational_ke_sum, &
           shell_angular_momentum
      double precision :: mixture_weight_sum, mixture_scale_factor
      integer :: radiative_zone_bounds(13,2), convective_zone_bounds(12,2)
      logical :: am_transport_convective_flag(json)
      integer :: num_radiative_zones, num_convective_zones
      integer :: core_cz_top_index0, envelope_cz_bottom_index0
      logical :: rotation_active0
      logical :: use_diffusion_y0, use_diffusion_z0, disk_locking_active0, &
           instability_transport_active0, ljdot00
      logical :: lovstc0, envelope_overshoot_active0, lovstm0, &
           use_pure_z_table0, lsemic0
      double precision :: disk_pressure0, disk_temperature0, &
           wind_saturation_omega0

      save

! If flag LARGE is set, model has failed to converge.  Terminate the run.
      if (model_failed_flag) then
          write(short_file_unit,1000)
 1000       format(1x,39('>'),40('<')/, &
                  "STARIN:        ***** RUN STOPPED *****")
          write(short_file_unit,1010)
 1010       format("STARIN: ***** MODEL FAILED TO CONVERGE *****")
            stop
      endif

! THIS SUBROUTINE READS IN THE INITIAL STELLAR MODEL
! INITIAL MODEL IS STORED IN LOGICAL UNIT IFIRST

! Set flags for reading an input model
      ikut_flag = 0
      iread = first_unit

! INITIALIZE VARIABLES
      convective_velocity = 0.0d0
      dlnrho_dlnt = -1.0d0
      dlnrho_dlnp = 1.0d0

! Flag LFIRST(NK) tells where to get the starting stellar model for the current
! step (step NK).  If LFIRST(NK) is true, read in the starting stellarmodel from
! the file specified by LU IFIRST.  If LFIRST(NK) is false, as starting model use
! the stellar model currently stored in memory.

      if (.not.first_call_flag(run_index)) then
! Use the model currently in memory as the starting model.
! DBG 2/92 CHANGED SO WILL RESCALE ENVELOPE MASS ON EACH NEW RUN
         if (rescale_kind(run_index).ne.1) call rscale(luminosity_lsun, &
              composition,log_mass,log_total_mass,num_shells,run_index, &
              total_mass_msun,convective_flag)
! Now skip over the reading and processing of an input model file
         goto 3000
      endif


!     Read in the starting model from LU IFIRST and process it.

!  Get file format of input model

      rewind iread
      read(iread,10) format_tag
 10   format(a4)

! ATEMP now contains a keyword describing the format of the input stellar
! model.  We decide what kind of model format it has and process accordingly.

      if (format_tag .eq. 'NMOD') then
         write(short_file_unit,12)
 12      format('STARIN:  Input model has YREC7 format')
         call getyrec7(log10_luminosity,envelope_fit_coeffs,mixing_length0, &
              age_gyr,timestep_yr,trial_sign_flag,composition,log_density, &
              luminosity_lsun,log_pressure,log_radius,log_mass, &
              log_total_mass,log_temperature,iread,short_file_unit, &
              core_cz_top_index0,envelope_cz_bottom_index0,convective_flag, &
              use_extended_composition0,rotation_active0,num_shells, &
              model_number,omega,fit_point_pressure,fit_point_radius, &
              total_mass_msun,log_teff,luminosity_breakdown, &
              trial_log_luminosity,trial_log_temperature, &
              fit_point_temperature, &
              atm_code,eos_code,hik_code,use_diffusion_y0,use_diffusion_z0, &
              disk_locking_active0,instability_transport_active0,ljdot00, &
              alok_code, &
              lovstc0,envelope_overshoot_active0,lovstm0,use_pure_z_table0, &
              lsemic0,compmix_code,disk_pressure0, &
              disk_temperature0,wind_saturation_omega0)
! First three lines above are YREC7 inputs
! Last three lines are MODEL2 add-ons

      else if (format_tag .eq. 'MOD2 ') then
         write(short_file_unit,16)
 16      format('STARIN:  Input model has MODEL2 format')
         call getmodel2(log10_luminosity,envelope_fit_coeffs,mixing_length0, &
              age_gyr,timestep_yr,trial_sign_flag,composition,log_density, &
              luminosity_lsun,log_pressure,log_radius,log_mass, &
              log_total_mass,log_temperature,iread, &
              core_cz_top_index0,envelope_cz_bottom_index0,convective_flag, &
              use_extended_composition0,rotation_active0,num_shells, &
              model_number,omega,fit_point_pressure,fit_point_radius, &
              total_mass_msun,log_teff,luminosity_breakdown, &
              trial_log_luminosity,trial_log_temperature, &
              fit_point_temperature, &
              atm_code,eos_code,hik_code,use_diffusion_y0,use_diffusion_z0, &
              disk_locking_active0,instability_transport_active0,ljdot00, &
              alok_code, &
              lovstc0,envelope_overshoot_active0,lovstm0,use_pure_z_table0, &
              lsemic0,compmix_code,disk_pressure0, &
              disk_temperature0,wind_saturation_omega0)
! First three lines above are YREC7 inputs
! Last three lines are MODEL2 add-ons

      else
         write(short_file_unit,20)
 20      format('STARIN: ***** RUN TERMINATED, INVALID INPUT', &
                ' MODEL FILE.  *****')
         stop
      endif

! Model has now been read in. Some post-processing is required.

      delta_time = seconds_per_year*timestep_yr
      delta_time_abs = dabs(delta_time)
      env_comp%stotal = log_total_mass

! CHECK TO ENSURE THAT MIX LENGTH, SURFACE B.C. AND CONVECTION ZONE
! THEORY OF THE MODEL ARE THE SAME AS USER PARAMETERS
! CHECKED ONLY FOR EVOLVED MODELS(MODEL NUMBER > 0)
! 1/92 Changed to not stop just give warning.
      if (model_number.gt.0) then
       mixing_length_matches = .true.
! Jan 12, 1989 : in STARIN changed test to see if model has correct mixing
! length from 1.0e-6 to 2.0e-3 because models only store
! mixing length to four sig digits.
       if (mixing_length0.gt.0.0d0) mixing_length_matches = &
            (dabs(cmixl-mixing_length0).lt.2.0d-3)
! MHP 9/03 FIXED TYPO
       if (.not.mixing_length_matches .or. use_extended_composition0.neqv. &
            use_extended_composition) then
          write(short_file_unit,1040) cmixl,mixing_length0, &
               use_extended_composition,lexcp0
          write(iowr,1040) cmixl,mixing_length0,use_extended_composition, &
               lexcp0
 1040       format(1x,'ERROR IN SUBROUTINE STARIN'/1x,'USER PARAMETERS', &
             ' OF WRONG TYPE FOR INITIAL MODEL'/1x,'MIXING LENGTH - USER' &
             ,' DESIRES',f7.3,' MODEL MIX LENGTH',f7.3/1x, &
             'EXTENDED COMP-USER DESIRES ',l1,' MODEL USED ',l1)
! =>RUN STOPPED DUE TO INCONSISTENCY BETWEEN MODEL AND RUN PARMS
!            STOP
       endif
      endif

! ENVELOPE DATA (Now bypassed)
! LNEW0 HAS BEEN READ IN, IF TRUE THEN RECOMPUTE ENVELOPE EVERY MODEL
! STORED ENVELOPE RECORDS ONLY USED FOR HE FLASH CALCS
!      IF(LNEW0) THEN
!       LNEW = .TRUE.
!      ELSE
!       LNEW = .NOT.LKUTHE
!      ENDIF
!      DO 80 I = 1,3
!       IF((.NOT.LNEW).AND.IABS(IO).NE.I) LNEW = .TRUE.
! 80   CONTINUE

      istore_flag = 0
      trial_sign_flag = 1.0d0
! Require recompute of envelope. No assurance of reliable
! input triangle
       envelope_recomputed_flag = .true.


! GET XNEW AND ZNEW FROM HENYEY POINTS

      env_comp%xnew = composition(1,num_shells)
      env_comp%znew = composition(3,num_shells)

! FOURTH PART:  - LOG J/M STORED

      if (rotation_active) then
       if (lwnew) then
! GENERATE A SOLID BODY ROTATION CURVE WITH OMEGA = WNEW;
! THIS IS DONE TO CONVERT A NON-ROTATING MODEL TO A ROTATING ONE.
          do 540 i = 1,num_shells
             omega(i) = wnew
 540        continue
       endif
      else
         do 570 i = 1,json
            pressure_rotation_factor(i) = 1.0d0
            temperature_rotation_factor(i) = 1.0d0
 570     continue
      endif
! KEEP IREAD OPEN
      rewind iread
! End of the reading and processing of an input model file.
      if (.not.first_call_flag(run_index)) goto 3000
!      IF(.NOT.LFIRST(NK).OR.NK.GT.1)GOTO 3000
!     MHP 10/24 MACHINERY TO ALTER THE HEAVY ELEMENT MIXTURE
!     THIS IS ONLY DONE if the first MODEL IS being READ IN, AND ONLY FOR A
! CHEMICALLY HOMOGENEOUS MODEL. IT CAN OVER-WRITE MASS FRACTIONS 4-15 WITH USER-SPECIFIED VALUES
! ISETMIX=1 -> CAN ADJUST CNO FRACTIONS ISETISO=1-> CHANGE ISOTOPE RATIOS
      if (change_cno_mixture_active .or. change_isotope_ratios_active) then
! ENSURE STARTING MODEL IS HOMOGENEOUS BEFORE EITHER IS CHANGED
         do i = 1,15
            reference_composition(i)=composition(i,1)
         end do
         do j = 2,num_shells
            do i = 1,15
               fraction_diff = abs(composition(i,j)-reference_composition(i))
               if (fraction_diff.gt.1.0d-6) then
                  write(*,592)i,j,fraction_diff
                  write(short_file_unit,592)i,j,fraction_diff
 592              format('SPECIES ',i3,' IN SHELL ',i5, &
                    ' DIFFERS FROM CENTER BY ',e12.4, &
                    ' MIX NOT MODIFIED IN EVOLVED MODEL')
                  goto 602
               endif
            end do
         end do
      endif
! LOOP FOR CHANGING CNO MIX
      if (change_cno_mixture_active) then
!     INFER CURRENT TOTAL CNO FRACTIONS AND SCALE ALL ISOTOPES BY THE RATIO BETWEEN
!     DESIRED AND CURRENT FRACTIONS. RELATIVE ISOTOPES ARE ADJUSTED IN THE ISOTOPE SECTION BELOW.
         total_carbon_cno_fraction = (reference_composition(5)+ &
              reference_composition(6))/reference_composition(3)
         total_nitrogen_cno_fraction = (reference_composition(7)+ &
              reference_composition(8))/reference_composition(3)
         total_oxygen_cno_fraction = (reference_composition(9)+ &
              reference_composition(10)+reference_composition(11))/ &
              reference_composition(3)
         carbon_scale_ratio = target_carbon_cno_fraction/ &
              total_carbon_cno_fraction
         nitrogen_scale_ratio = target_nitrogen_cno_fraction/ &
              total_nitrogen_cno_fraction
         oxygen_scale_ratio = target_oxygen_cno_fraction/ &
              total_oxygen_cno_fraction
         write(*,*)target_carbon_cno_fraction,target_nitrogen_cno_fraction, &
              target_oxygen_cno_fraction
         do i = 5,6
            do j = 1,num_shells
               composition(i,j)=carbon_scale_ratio*composition(i,j)
            end do
         end do
         do i = 7,8
            do j = 1,num_shells
               composition(i,j)=nitrogen_scale_ratio*composition(i,j)
            end do
         end do
         do i = 9,11
            do j = 1,num_shells
               composition(i,j)=oxygen_scale_ratio*composition(i,j)
            end do
         end do
         write(*,594)(reference_composition(k),k=5,11), &
              (composition(k,1),k=5,11)
         write(short_file_unit,594)(reference_composition(k),k=5,11), &
              (composition(k,1),k=5,11)
 594     format('CNO MIX CHANGED IN STARIN. OLD C12 C13 N14' &
           ' N15 O16 O17 O18 ',7e12.4,' NEW ',7e12.4)
      endif
! DESIRED ISOTOPE RATIOS AND LIGHT ELEMENT ABUNDANCES ASSIGNED.
!     AT PRESENT B10,B11,N15,O17 ARE NOT USED AND THUS NOT ALTERED.
!     START WITH LIGHT ELEMENTS
      if (change_isotope_ratios_active) then
         sum_c12_c13 = composition(5,1)+composition(6,1)
         sum_o16_o18 = composition(9,1)+composition(11,1)
         do j = 1,num_shells
            composition(4,j)=initial_he3_fraction
            composition(5,j)= c12_to_c13_ratio*sum_c12_c13/ &
                 (1.0d0+c12_to_c13_ratio)
            composition(6,j)= sum_c12_c13/(1.0d0+c12_to_c13_ratio)
            composition(9,j)= o16_to_o18_ratio*sum_o16_o18/ &
                 (1.0d0+o16_to_o18_ratio)
            composition(11,j)= sum_o16_o18/(1.0d0+o16_to_o18_ratio)
            composition(12,j)=initial_h2_fraction
            composition(13,j)=initial_li6_fraction
            composition(14,j)=initial_li7_fraction
            composition(15,j)=initial_be9_fraction
         end do
         write(*,593)(reference_composition(k),k=4,15), &
              (composition(k,1),k=4,15)
         write(short_file_unit,593)(reference_composition(k),k=4,15), &
              (composition(k,1),k=4,15)
 593     format('CNO ISOTOPES AND LIGHT ELEMENTS CHANGED IN ', &
              'STARIN. OLD HE3 C12 C13 N14 N15 O16 O17 O18 H2 LI6 ', &
               'LI7 BE9',12e12.4,' NEW ',12e12.4)
      endif
 602  continue
 3000 continue

!     The following code enables us to extend the model from the current
!     inner most shell to a point ncloser to center, if flag LCORE is set.
! 9/98 MHP REPLACE AS FOLLOWS:
! ADD INTEGER NUMBER OF POINTS USING THE DESIRED SPACING IN MASS,HPTTOL.
! USE STELLAR STRUCTURE EQUATIONS CONSISTENT WITH THE ASSUMPTIONS
! IN THE MODEL, NAMELY CONSTANT ENERGY GENERATION RATE AND DENSITY
! INTERIOR TO THE LOCATION OF THE CENTRAL FITTING POINT.
! MHP 9/14 CHANGED SO THAT MOVING THE CORE FITTING IS ATTACHED TO ANY RUN
! WHICH READS IN THE STARTING MODEL; THIS AVOIDS OVER-WRITING THE CHANGE
! IN AUTO-CALIBRATED SOLAR MODELS
      if (extend_core_inward .and. first_call_flag(run_index)) then
!      IF(LCORE .AND. NK .EQ. 1) THEN
! AVOID SHUFFLING POINTS BY ASSIGNING NEW CENTRAL POINTS IN INTEGER
! MULTIPLES OF THE CENTRAL POINT SPACING.
!     MCORE is number of shells to extrapolate to new core.
!     FCORE is factor to reduce inner mass shell.
          num_core_shells_added = int(dlog10(core_mass_reduction_factor)/ &
               chi_grid_scale(2))+1
          core_mass_reduction_factor = dble(num_core_shells_added)* &
               chi_grid_scale(2)
          num_shells_extended = num_shells + num_core_shells_added
          if (num_shells_extended .gt. json) then
             write(short_file_unit,476)"STARIN: Unable to extend core inward ", &
                   "- JSON too small"
             write(short_file_unit,477) "STARIN: Required size =", &
                    num_shells_extended, ", JSON = ", json
             write(short_file_unit,478) "STARIN:  ***** RUN TERMINATED *****"
  476        format(2a)
  477        format(a, i8, a, i8)
  478        format(a)
          endif
          core_shell_spacing = chi_grid_scale(2)
! shift data for remaining points by the required number
          do i=num_shells,1, -1
             log_mass(i+num_core_shells_added) = log_mass(i)
             log_radius(i+num_core_shells_added) = log_radius(i)
             luminosity_lsun(i+num_core_shells_added) = luminosity_lsun(i)
             log_pressure(i+num_core_shells_added) = log_pressure(i)
             log_temperature(i+num_core_shells_added) = log_temperature(i)
             convective_flag(i+num_core_shells_added) = convective_flag(i)
             do j=1, 15
                composition(j,i+num_core_shells_added) = composition(j,i)
             end do
             omega(i+num_core_shells_added) = omega(i)
          end do
          first_original_shell = num_core_shells_added+1
! MARCH INWARD IN MASS FROM THE INNERMOST MODEL POINT.
! ASSUME EPSILON=CONSTANT AND DEL=CONSTANT
          central_log_density = log_density(1)
          actual_gradient = (log_temperature(2)-log_temperature(1))/ &
               (log_pressure(2)-log_pressure(1))
          central_shell_luminosity = luminosity_lsun(1)
! MHP 4/12 FACTOR FOR ESTIMATING RHO FROM P AND T
          density_estimate_offset = log_pressure(first_original_shell)- &
               log_density(first_original_shell)- &
               log_temperature(first_original_shell)
          do i = num_core_shells_added,1,-1
             log_mass(i) = log_mass(i+1)-core_shell_spacing
! USE M  = 4PI/3*RHOC*R**3 TO GET R AS A FUNCTION OF M
             log_radius(i) = cc13*(log_mass(i)-c4pi3l-central_log_density)
! USE EPSILON=CONSTANT TO GET L
             luminosity_lsun(i) = exp(ln10*(log_mass(i)- &
                  log_mass(first_original_shell)))*central_shell_luminosity
! USE HYDROSTATIC EQUILIBRIUM TO GET P
             trial_log_pressure = log_pressure(i+1)
             temp_scratch = exp(ln10*(cgl+2.0d0*log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*log_radius(i)))
             trial_log_pressure = trial_log_pressure+0.5d0*temp_scratch* &
                  core_shell_spacing
             temp_scratch = exp(ln10*(cgl+2.0d0*log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*log_radius(i)))
             log_pressure(i) = log_pressure(i+1)+temp_scratch* &
                  core_shell_spacing
! ASSUME R/C FLAG IS THE SAME AS FOR THE FIRST POINT
             convective_flag(i) = convective_flag(first_original_shell)
! ASSUME DEL= CONSTANT IN THE CORE
             log_temperature(i) = log_temperature(i+1)+temp_scratch* &
                  actual_gradient*core_shell_spacing
! ASSUME OMEGA = CONSTANT
             omega(i) = omega(first_original_shell)
! ASSUME COMPOSITION IS UNIFORM
             do j=1, 15
                composition(j,i) = composition(j,first_original_shell)
             end do
! CALL EQUATION OF STATE TO GET CONSISTENT DENSITY
!           LDERIV = .FALSE.
!           LOCOND = .FALSE.
!           LATMO = .FALSE.
!             KSAHA = 0
!           X = HCOMP(1,I)
!           Z = HCOMP(3,I)
!           PL = HP(I)
!             P = EXP(CLN*PL)
!             TL = HT(I)
!             T = EXP(CLN*TL)
!           DL = HD(I+1)
!             D = EXP(CLN*DL)
!           FPL = 1.0D0
!           FTL = 1.0D0
!           CALL EQSTAT(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,
!     *                   RMU,AMU,EMU,ETA,QDT,QDP,QCP,DELA,QDTT,QDTP,
!     *                   QAT,QAP,QCPT,QCPP,LDERIV,LATMO,KSAHA)
!             HD(I) = DL
! MHP 4/12 REPLACED (BROKEN) CALL TO EQSTAT WITH LOCAL ESTIMATE FOR RHO
             log_density(i) = log_pressure(i) - log_temperature(i) - &
                  density_estimate_offset
          end do
          num_shells = num_shells_extended
      end if
! End of code to extend core inward

! PERFORM RESCALING OF FIRST MODEL IF REQUIRED
      if (rescale_kind(run_index).ne.1) call rscale(luminosity_lsun, &
           composition,log_mass,log_total_mass,num_shells,run_index, &
           total_mass_msun,convective_flag)

! What is the metallicity of your model? (surface X and Z)
! put in SETUPOPAC here, and take out of setups.
! DBG 2/92 CHANGED SO THAT WORKS ON FIRST MODEL OF EACH NEW RUN
! CHANGE FITTING POINT OF THE ENVELOPE INTEGRATION IF REQUESTED;
! PROCEDURE IS:
! IF THE REQUESTED ENVELOPE MASS(SENV0) IS GREATER THAN THE
! CURRENT ONE, DELETE POINTS UNTIL THE MASS IS EQUAL.
! LINEAR INTERPOLATION BETWEEN THE LAST POINT ABOVE THE FITTING
! MASS AND THE FIRST POINT BELOW IT IS USED;COMPOSITION IS
! ASSUMED EQUAL TO THE LAST OLD POINT BELOW THE NEW FITTING POINT.
! IF THE REQUESTED ENVELOPE MASS IS SMALLER THAN THE CURRENT ONE,
! AN ENVELOPE WITH TEFF AND L EQUAL TO THE STORED MODEL VALUE
! IS INTEGRATED FROM THE SURFACE TO THE DESIRED FITTING POINT.
! 1 NEW POINT IS ADDED, AND THE COMPOSITION OF THE NEW POINT
! IS ASSUMED EQUAL TO THAT OF THE LAST OLD POINT.
      if (change_envelope_mass_flag) then
       if (requested_envelope_mass.gt.0.0d0) requested_envelope_mass = &
            -requested_envelope_mass
! DBG 2/92 CHANGED MINIMUM FROM 1.0D-10 TO 1.0D-12
! RESTRICT MIMIMUM ENVELOPE MASS;1.0D-12 CORRESPONDS TO TAU=2/3
! FOR THIS PURPOSE(BASE OF ATMOSPHERE).
       if (requested_envelope_mass.gt.-1.0d-12) requested_envelope_mass = &
            -1.d-12
       env_comp%senv = log_mass(num_shells) - log_total_mass
       old_senv = env_comp%senv
       if (env_comp%senv.eq.requested_envelope_mass) goto 599
       num_species = 11
       if (use_extended_composition) num_species = 15
       if (requested_envelope_mass.lt.env_comp%senv) then
! NEW ENVELOPE DEEPER THAN THE OLD ONE
          target_log_mass_at_fit = log_total_mass+requested_envelope_mass
          do 575 i = num_shells-1,1,-1
             if (log_mass(i).lt.target_log_mass_at_fit) goto 580
 575        continue
! ENVELOPE MASS DESIRED WITHIN FIRST POINT;PRINT NASTY MESSAGE
! AND ABORT.
          write(short_file_unit,576)requested_envelope_mass
 576        format(5x,'ERROR IN SUBROUTINE STARIN'/5x,'DESIRED', &
              ' ENVELOPE MASS',1pe22.13,' TOO LARGE'/5x,'ENVELOPE', &
              ' MASS NOT CHANGED')
          goto 599
 580        num_shells = i + 1
          env_comp%senv = requested_envelope_mass
          interior_interp_fraction = (target_log_mass_at_fit-log_mass(i))/ &
               (log_mass(i+1) - log_mass(i))
          log_mass(num_shells) = target_log_mass_at_fit
          log_density(num_shells) = log_density(i) + &
               interior_interp_fraction*(log_density(i+1) - log_density(i))
          luminosity_lsun(num_shells) = luminosity_lsun(i) + &
               interior_interp_fraction*(luminosity_lsun(i+1) - &
               luminosity_lsun(i))
          log_pressure(num_shells) = log_pressure(i) + &
               interior_interp_fraction*(log_pressure(i+1) - log_pressure(i))
          log_radius(num_shells) = log_radius(i) + &
               interior_interp_fraction*(log_radius(i+1) - log_radius(i))
          log_temperature(num_shells) = log_temperature(i) + &
               interior_interp_fraction*(log_temperature(i+1) - &
               log_temperature(i))
          do 585 j = 1,num_species
             composition(j,num_shells) = composition(j,i)
 585        continue
          env_comp%xnew = composition(1,num_shells)
          env_comp%znew = composition(3,num_shells)
          if (rotation_active) omega(num_shells) = omega(i) + &
               interior_interp_fraction*(omega(i+1)-omega(i))
          if (convective_flag(i).and.convective_flag(i+1)) then
             convective_flag(num_shells) = .true.
          else if (.not.convective_flag(i).and. .not.convective_flag(i+1)) &
               then
             convective_flag(num_shells) = .false.
          else
! CALL BASIC PHYSICS ROUTINES TO DETERMINE IF THE NEW LAST SHELL IS
! CONVECTIVE OR RADIATIVE.
             want_derivatives = .false.
             local_conductive_opacity_flag = .false.
             in_atmosphere = .true.
             saha_state = 0
             hydrogen_fraction = composition(1,num_shells)
             metal_fraction = composition(3,num_shells)
             log10_pressure = log_pressure(num_shells)
             log10_temperature = log_temperature(num_shells)
             log10_density = log_density(num_shells)
             shell_luminosity_lsun = luminosity_lsun(num_shells)
             log10_radius = log_radius(num_shells)
             log10_mass = log_mass(num_shells)
             point_pressure_rotation_factor = 1.0d0
             point_temperature_rotation_factor = 1.0d0
             idt = 15
             do 588 kk = 1,4
              idd(kk) = 5
 588           continue
               if (use_mhd_eos) then
                  call meqos(log10_temperature,temperature,log10_pressure, &
                       pressure,log10_density,density,hydrogen_fraction, &
                       metal_fraction,beta,beta_inverse,beta14,ion_fraction, &
                       specific_gas_constant,ion_mean_weight_inverse, &
                       electron_mean_weight_inverse, &
                       electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
                       specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
!      *                 QDTP,QAT,QAP,QCPT,QCPP,LDERIV,LATMO,KSAHA)  ! KC 2025-05-31
                       dlnrho_dlnp_dt,adiabatic_gradient_dt, &
                       adiabatic_gradient_dp,specific_heat_cp_dt, &
                       specific_heat_cp_dp)
                  if (use_debye_huckel_correction) then
                     debye_huckel_x = composition(1,num_shells)
                     debye_huckel_y = composition(2,num_shells)+composition(4,num_shells)
                     debye_huckel_z_total = composition(3,num_shells)
                     debye_huckel_z(1) = composition(5,num_shells)+composition(6,num_shells)
                     debye_huckel_z(2) = composition(7,num_shells)+composition(8,num_shells)
                     debye_huckel_z(3) = composition(9,num_shells)+ &
                          composition(10,num_shells)+composition(11,num_shells)
                  end if
                  call eqstat(log10_temperature,temperature,log10_pressure, &
                       pressure,log10_density,density,hydrogen_fraction, &
                       metal_fraction,beta,beta_inverse,beta14,ion_fraction, &
                       specific_gas_constant,ion_mean_weight_inverse, &
                       electron_mean_weight_inverse, &
                       electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
                       specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
                       dlnrho_dlnp_dt,adiabatic_gradient_dt, &
                       adiabatic_gradient_dp,specific_heat_cp_dt, &
                       specific_heat_cp_dp,want_derivatives,in_atmosphere, &
                       saha_state)
               end if
               call getopac(log10_density, log10_temperature, &
                    hydrogen_fraction, metal_fraction, opacity, &
                    log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ion_fraction)
               iovim = -1
               call tpgrad(log10_temperature,temperature,log10_pressure, &
                    pressure,density,log10_radius,log10_mass, &
                    shell_luminosity_lsun,opacity,dlnrho_dlnt,dlnrho_dlnp, &
                    dlnkap_dlnt,dlnkap_dlnrho, &
                    specific_heat_cp,actual_gradient,radiative_gradient, &
                    adiabatic_gradient,dlnrho_dlnt_dt,dlnrho_dlnp_dt, &
                    adiabatic_gradient_dt,adiabatic_gradient_dp, &
                    dgrad_dt_component,dgrad_dp_component,dgrad_dr_component, &
                    specific_heat_cp_dt,specific_heat_cp_dp, &
                    convective_velocity,want_derivatives,is_convective, &
                    point_pressure_rotation_factor, &
                    point_temperature_rotation_factor,log_teff)
               log_density(num_shells) = log10_density
               convective_flag(num_shells) = is_convective
          endif
       else
! DESIRED ENVELOPE MASS LESS THAN CURRENT VALUE.
            old_last_shell = num_shells
            saved_env_step_max = env_step_max
            saved_env_step_min = env_step_min
            saved_env_step_begin = env_step_begin
            env_step_max = chi_grid_scale(8)
            env_step_min = chi_grid_scale(8)
            env_step_begin = chi_grid_scale(8)
!          SENV = SENV0
          save_boundary_flag = .false.
          print_flag = .true.
          katm = 0
          kenv = 0
          saha_state = 0
          shell_luminosity_lsun = dexp(ln10*log10_luminosity)
          log10_radius = 0.5d0*(log10_luminosity + solar_luminosity_cgs - &
               4.0d0*log_teff - c4pil - csigl)
          log10_gravity = cgl + env_comp%stotal - log10_radius - log10_radius
          hydrogen_fraction = composition(1,num_shells)
          metal_fraction = composition(3,num_shells)
          point_pressure_rotation_factor = 1.0d0
          point_temperature_rotation_factor = 1.0d0
          vertex_index=0
          log10_pressure_limit = log_pressure(num_shells)
! DBG PULSE: DO NOT DO PULSE OUTPUT
            pulse_print_flag = .false.
            if (use_debye_huckel_correction) then
               debye_huckel_x = composition(1,num_shells)
               debye_huckel_y = composition(2,num_shells)+composition(4,num_shells)
               debye_huckel_z_total = composition(3,num_shells)
               debye_huckel_z(1) = composition(5,num_shells)+composition(6,num_shells)
               debye_huckel_z(2) = composition(7,num_shells)+composition(8,num_shells)
               debye_huckel_z(3) = composition(9,num_shells)+composition(10,num_shells)+ &
                    composition(11,num_shells)
            end if
! MHP 10/02  define ISTORE - used in ENVINT
            envint_unused_flag = 0
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
          if (envelope_zone_index.eq.num_shells.and.spot_filling_factor.ne. &
               0.0.and.spot_temp_contrast.ne.1.0) then
               spot_adjusted_log_teff = log_teff - 0.25*log10(&
                    spot_filling_factor * spot_temp_contrast**4.0 + 1.0 - &
                    spot_filling_factor)
          else
             spot_adjusted_log_teff = log_teff
          endif
          call envint(shell_luminosity_lsun,point_pressure_rotation_factor, &
                 point_temperature_rotation_factor,log10_gravity,env_comp%stotal, &
                 vertex_index,print_flag,save_boundary_flag, &
                 log10_pressure_limit,log10_radius,spot_adjusted_log_teff, &
                 hydrogen_fraction,metal_fraction,envint_dummy1, &
                 envint_unused_flag,katm,kenv,saha_state,envint_dummy2, &
                 envint_dummy3,envint_dummy4,pulse_print_flag)
! G Somers END
            env_step_max = saved_env_step_max
            env_step_min = saved_env_step_min
            env_step_begin = saved_env_step_begin
          env_comp%senv = requested_envelope_mass
            if (num_shells+env_struct%num_env_points.ge.json) stop 9999
! ENFORCE CONSISTENCY WITH THE INTERIOR SOLUTION;
! ADJUST THE (P, RHO, T, R) POINTS TO BE CONSISTENT
! WITH THE LAST MODEL POINT.
            pressure_offset = log_pressure(num_shells) - env_struct%env_log10_pressure(1)
            density_offset = log_density(num_shells) - env_struct%env_log10_density(1)
            temperature_offset = log_temperature(num_shells) - &
                 env_struct%env_log10_temperature(1)
            radius_offset = log_radius(num_shells) - env_struct%env_log10_radius(1)
            do j = 1,env_struct%num_env_points - 1
               env_struct%env_log10_density(j) = env_struct%env_log10_density(j+1)+density_offset
               env_struct%env_log10_pressure(j) = env_struct%env_log10_pressure(j+1)+pressure_offset
               env_struct%env_log10_radius(j) = env_struct%env_log10_radius(j+1)+radius_offset
               env_struct%env_log10_mass(j) = env_struct%env_log10_mass(j+1)
               env_struct%env_log10_temperature(j) = env_struct%env_log10_temperature(j+1)+ &
                    temperature_offset
               env_struct%env_hydrogen_fraction(j) = env_struct%env_hydrogen_fraction(j+1)
               env_struct%env_metal_fraction(j) = env_struct%env_metal_fraction(j+1)
            end do
            env_struct%num_env_points = env_struct%num_env_points - 1
            do j = num_shells+1,num_shells+env_struct%num_env_points
               env_point_index = j-num_shells
! LUMINOSITY ASSUMED CONSTANT
               luminosity_lsun(j) = luminosity_lsun(num_shells)
! INCLUDE NEW POINTS UP TO THE DIFFERENT DESIRED FITTING POINT
               if (env_struct%env_log10_mass(env_point_index).le.env_comp%senv) then
                  log_density(j) = env_struct%env_log10_density(env_point_index)
                  log_pressure(j) = env_struct%env_log10_pressure(env_point_index)
                  log_radius(j) = env_struct%env_log10_radius(env_point_index)
                  log_mass(j) = env_struct%env_log10_mass(env_point_index) + env_comp%stotal
                  log_temperature(j) = env_struct%env_log10_temperature(env_point_index)
                  composition(1,j) = env_struct%env_hydrogen_fraction(env_point_index)
                  composition(3,j) = env_struct%env_metal_fraction(env_point_index)
                  do k = 4,num_species
                     composition(k,j) = composition(k,num_shells)
                  end do
                  composition(2,j)=1.0d0-composition(1,j)-composition(3,j)- &
                       composition(4,j)
                  convective_flag(j) = env_struct%env_convective_flag(env_point_index)
               else
! POINTS BEYOND THIS ARE ABOVE THE NEW DESIRED FITTING POINT;
! INTERPOLATE LINEARLY, SET NEW NUMBER OF TOTAL POINTS, AND EXIT
                  if (env_point_index.eq.1) then
! INTERPOLATE BETWEEN THE LAST INTERIOR POINT AND THE FIRST ENVELOPE POINT
                     lower_mass_coord = log_mass(num_shells)
                     target_mass_coord = env_comp%stotal + env_comp%senv
                     upper_mass_coord = env_struct%env_log10_mass(env_point_index) + &
                          env_comp%stotal
                     if (upper_mass_coord-lower_mass_coord.lt.1.0d-14) &
                          stop 9998
                     envelope_interp_fraction = (target_mass_coord- &
                          lower_mass_coord)/(upper_mass_coord- &
                          lower_mass_coord)
                     log_density(j) = log_density(num_shells)+ &
                          envelope_interp_fraction*(env_struct%env_log10_density( &
                          env_point_index)-log_density(num_shells))
                     log_pressure(j) = log_pressure(num_shells)+ &
                          envelope_interp_fraction*(env_struct%env_log10_pressure( &
                          env_point_index)-log_pressure(num_shells))
                     log_radius(j) = log_radius(num_shells)+ &
                          envelope_interp_fraction*(env_struct%env_log10_radius( &
                          env_point_index)-log_radius(num_shells))
                     log_mass(j) = target_mass_coord
                     log_temperature(j) = log_temperature(num_shells)+ &
                          envelope_interp_fraction*(env_struct%env_log10_temperature( &
                          env_point_index)-log_temperature(num_shells))
                     composition(1,j) = composition(1,num_shells)+ &
                          envelope_interp_fraction*(composition(1,num_shells) &
                          -env_struct%env_hydrogen_fraction(env_point_index))
                     composition(3,j) = composition(3,num_shells)+ &
                          envelope_interp_fraction*(composition(3,num_shells) &
                          -env_struct%env_metal_fraction(env_point_index))
                     do k = 4,num_species
                        composition(k,j) = composition(k,num_shells)
                     end do
                     composition(2,j)=1.0d0-composition(1,j)- &
                          composition(3,j)-composition(4,j)
                     if (env_struct%env_convective_flag(env_point_index).or. &
                          convective_flag(num_shells)) then
                        convective_flag(j) = .true.
                     else
                        convective_flag(j) = .false.
                     endif
                  else
! INTERPOLATE BETWEEN THE LAST 2 ENVELOPE POINTS
                     lower_mass_coord = env_struct%env_log10_mass(env_point_index-1) + &
                          env_comp%stotal
                     target_mass_coord = env_comp%stotal + env_comp%senv
                     upper_mass_coord = env_struct%env_log10_mass(env_point_index) + &
                          env_comp%stotal
                     if (upper_mass_coord-lower_mass_coord.lt.1.0d-14) &
                          stop 9998
                     envelope_interp_fraction = (target_mass_coord- &
                          lower_mass_coord)/(upper_mass_coord- &
                          lower_mass_coord)
                     log_density(j) = env_struct%env_log10_density(env_point_index-1)+ &
                          envelope_interp_fraction*(env_struct%env_log10_density( &
                          env_point_index)-env_struct%env_log10_density(env_point_index-1))
                     log_pressure(j) = env_struct%env_log10_pressure(env_point_index-1)+ &
                          envelope_interp_fraction*(env_struct%env_log10_pressure( &
                          env_point_index)-env_struct%env_log10_pressure(env_point_index-1))
                     log_radius(j) = env_struct%env_log10_radius(env_point_index-1)+ &
                          envelope_interp_fraction*(env_struct%env_log10_radius( &
                          env_point_index)-env_struct%env_log10_radius(env_point_index-1))
                     log_mass(j) = target_mass_coord
                     log_temperature(j) = env_struct%env_log10_temperature( &
                          env_point_index-1)+envelope_interp_fraction*( &
                          env_struct%env_log10_temperature(env_point_index)- &
                          env_struct%env_log10_temperature(env_point_index-1))
                     composition(1,j) = env_struct%env_hydrogen_fraction( &
                          env_point_index-1)+envelope_interp_fraction*( &
                          env_struct%env_hydrogen_fraction(env_point_index)- &
                          env_struct%env_hydrogen_fraction(env_point_index-1))
                     composition(3,j) = env_struct%env_metal_fraction( &
                          env_point_index-1)+envelope_interp_fraction*( &
                          env_struct%env_metal_fraction(env_point_index)- &
                          env_struct%env_metal_fraction(env_point_index-1))
                     do k = 4,num_species
                        composition(k,j) = composition(k,num_shells)
                     end do
                     composition(2,j)=1.0d0-composition(1,j)- &
                          composition(3,j)-composition(4,j)
                     if (env_struct%env_convective_flag(env_point_index).or. &
                          env_struct%env_convective_flag(env_point_index-1)) then
                        convective_flag(j) = .true.
                     else
                        convective_flag(j) = .false.
                     endif
                  endif
                  num_shells = j
                  goto 587
               endif
            end do
! ASSIGN THE BOUNDARY AT THE PHOTOSPHERE FOR ENVELOPE MASS BELOW 1.0D-12.
            num_shells = num_shells + env_struct%num_env_points
 587        continue
            if (rotation_active) then
               do j = old_last_shell+1,num_shells
                  omega(j) = omega(old_last_shell)
                  specific_angular_momentum(j) = cc23*omega(old_last_shell)* &
                       10.0d0**(2.0d0*log_radius(j))
               end do
            endif
            write(*,910)
 910  format(1x,'NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS'/ &
            ' J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F')
      write(*,911)(j,log_density(j),luminosity_lsun(j),log_pressure(j), &
                   log_radius(j),log_mass(j)-env_comp%stotal, &
                   log_temperature(j),convective_flag(j), j = old_last_shell,num_shells)
 911  format(i5,1p6e16.8,l2)
!          M = M + 1
!          HS(M) = HSTOT + SENV
!          HD(M) = ED
!          HL(M) = HL(M-1)
!          HP(M) = EP
!          HR(M) = ER
!          HT(M) = ET
!          LC(M) = EVEL.GT.0.0D0
!          DO 590 J = 1,JEND
!             HCOMP(J,M) = HCOMP(J,M-1)
! 590        CONTINUE
!          XNEW = HCOMP(1,M)
!          ZNEW = HCOMP(3,M)
!          IF(LROT) OMEGA(M) = OMEGA(I) + FS*(OMEGA(I+1)-OMEGA(I))
       endif
       envelope_recomputed_flag = .true.
       write(short_file_unit,597)old_senv,env_comp%senv
 597     format(5x,'***** NEW ENVELOPE MASS CALCULATED *****'/8x, &
              'OLD SENV ',1pe22.13,'  NEW SENV',e22.13)
 599     continue
      endif

! SET UP WEIGHTS AND MASSES
! HS1 IS THE UNLOGGED HS; HS2 IS THE MASS OF THE SHELL(ALSO NOT LOG).
      next_mass = dexp(ln10*log_mass(1))
      curr_mass = - next_mass
      do 120 i = 2,num_shells
       prev_mass = curr_mass
       curr_mass = next_mass
       next_mass = dexp(ln10*log_mass(i))
       enclosed_mass(i-1) = curr_mass
       shell_mass(i-1) = 0.5d0*(next_mass-prev_mass)
 120  continue
      enclosed_mass(num_shells) = next_mass
      shell_mass(num_shells) = dexp(ln10*log_total_mass) - 0.5d0*(curr_mass+ &
           next_mass)

      if (rotation_active) then
! CALCULATE FP,FT,R0 AND ETA2 GIVEN OMEGA
       call fpft(log_density,log_radius,log_mass,num_shells,omega, &
            eta_squared,pressure_rotation_factor,temperature_rotation_factor, &
            mean_gravity,mean_radius)
! FIND MOMENT OF INERTIA(HI)
!        CALL MOMI(ETA2,HD,HR,HS,HS2,1,M,OMEGA,R0,HI,QIW,M)  ! KC 2025-05-31
       call momi(eta_squared,log_radius,log_mass,shell_mass,1,num_shells, &
            omega,mean_radius,moment_of_inertia,di_domega)
! GIVEN OMEGA AND I, FIND ANGULAR MOMENTUM AND ROTATIONAL K.E.
       angular_momentum_sum = 0.0d0
       rotational_ke_sum = 0.0d0
       do 550 i = 1,num_shells
          shell_angular_momentum = omega(i)*moment_of_inertia(i)
          specific_angular_momentum(i) = shell_angular_momentum/shell_mass(i)
          kinetic_energy_rot(i) = 0.5d0*omega(i)*shell_angular_momentum
          angular_momentum_sum = angular_momentum_sum+shell_angular_momentum
          rotational_ke_sum = rotational_ke_sum + kinetic_energy_rot(i)
 550     continue
       write(short_file_unit,560)total_angular_momentum, &
            angular_momentum_sum,total_rotational_ke,rotational_ke_sum
 560     format(1x,'TOTAL J OF STAR - PREVIOUS ',1pe21.13,' NEW ', &
              1pe21.13/1x,'TOTAL ROTATIONAL K.E. OF STAR - PREVIOUS ', &
              1pe21.13,' NEW ',1pe21.13)
       total_angular_momentum = angular_momentum_sum
       total_rotational_ke = rotational_ke_sum
      endif

      if (run_index.gt.1) goto 630
! SET UP MASS FRACTIONS AND NUMBER FRACTIONS OF ELEMENTS IN
! ENVELOPE.
! DBG 1/96 V (ENVELOPE MASS FRACTIONS WAS NORMALLY READ IN VIA
! RDLAOL. NOW THAT THESE OPACITIES ARE OBSOLETE VNEW IS INTRODUCED
! HOLD THE RELATIVE MASS FRACTIONS OF THE ELEMENTS (SEE PARMIN).
! TO MAINTAIN BACKWARD COMPATIBILITY IF LLAOL=T THEN USE V READ IN
! VIA RDLAOL OTHERWISE USE VNEW.
      if (.not.llaol) then
         do i=1, 12
            species_mix_weights(i)=vnew(i)
         end do
      end if

! COMPUTE SURFACE MIX VALUES.
! ZENVM = Z OTHER THAN CNO CYCLE ELEMENTS;
! AMUENV = MEAN ATOMIC WEIGHT OF SURFACE Z;
! FXENV = NUMBER DENSITY OF SPECIES .
      env_comp%envelope_hydrogen_fraction = env_comp%xnew
      env_comp%envelope_metal_fraction = env_comp%znew
      envelope_helium_fraction = 1.0d0 - env_comp%envelope_hydrogen_fraction - &
           env_comp%envelope_metal_fraction - composition(4,num_shells)
      envelope_he3_fraction = composition(4,num_shells)
! EVERYTHING BUT V(7)=H, AND V(12)=HE
      mixture_weight_sum = species_mix_weights(1)+species_mix_weights(2)+ &
           species_mix_weights(3)+species_mix_weights(4)+ &
           species_mix_weights(5)+species_mix_weights(6)+ &
           species_mix_weights(8)+species_mix_weights(9)+ &
           species_mix_weights(10)+species_mix_weights(11)
      env_comp%zenvm = env_comp%envelope_metal_fraction*(mixture_weight_sum - &
           species_mix_weights(6)-species_mix_weights(8)- &
           species_mix_weights(9))/mixture_weight_sum
      mixture_scale_factor = env_comp%envelope_metal_fraction/mixture_weight_sum
      species_mix_weights(7) = env_comp%envelope_hydrogen_fraction/ &
           mixture_scale_factor
      species_mix_weights(12) = (1.0d0-env_comp%envelope_hydrogen_fraction- &
           env_comp%envelope_metal_fraction)/mixture_scale_factor
      mixture_weight_sum = 0.0d0
      do 610 i = 1,12
       species_mix_weights(i) = mixture_scale_factor*species_mix_weights(i)/ &
            atomic_weight(i)
       mixture_weight_sum = mixture_weight_sum + species_mix_weights(i)
 610  continue
      env_comp%amuenv = mixture_weight_sum
      mixture_scale_factor = 1.0d0/env_comp%amuenv
! DBG 1/96 FXENV ARE NUMBER FRACTIONS OF ELEMENTS REQURIED
! BY EOS ROUTINES (SEE EQSTAT AND EQSAHA)
      do 620 i = 1,12
       env_comp%fxenv(i) = species_mix_weights(i)*mixture_scale_factor
 620  continue
!     FIND SURFACE COMPOSITION OPACITY TABLE
!     FIRST FIND INTERPOLATING FACTOR FOR COMPOSITION
 630  continue
! DBG 11/95 GENERATE NEW SURFACE OPACITY TABLES
      call surfopac(env_comp%envelope_hydrogen_fraction)
      if (use_scv_eos) then
         call setscv
      endif

! CLONE P,T,R,L ARRAY TO DUMMY ARRAY HPOLD.
! HPOLD IS USED TO LIMIT THE TIMESTEP BASED ON CHANGES FROM
! MODEL TO MODEL IN P,T,R,L.
      do 710 i = 1,num_shells
         prev_model%old_pressure(i) = log_pressure(i)
         prev_model%old_temperature(i) = log_temperature(i)
         prev_model%old_radius(i) = log_radius(i)
         prev_model%old_luminosity(i) = luminosity_lsun(i)
!  JVS 04/14 Added Teff to the list of saved values
         prev_model%old_teff = log_teff
!  JVS 05/25 Added model number to list of saved values
       prev_model%old_num_zones = num_shells
 710  continue
      if (rotation_active) then
         do 720 i = 1,num_shells
          old_omega(i) = omega(i)
 720     continue
      endif

! 8/17 G Somers
!  FIND BASIC PHYSICAL QUANTITIES. THIS CODE STOLEN FROM PHYSIC
!  FIND ACTUAL AND ADIABATIC TEMPERATURE GRADIENTS,OPACITY,AND
!  MEAN MOLECULAR WEIGHT FOR ALL RADIATIVE SHELLS.
!      if(.False.)then
!      LDERIV = .FALSE.
!      LOCOND = .FALSE.
!      LATMO = .FALSE.
!      IDT = 15
!      DO 725 I = 1,4
!         IDD(I) = 5
! 725  CONTINUE
!      DO 730 IM = 1,M
!         SL = HS(IM)
!         TL = HT(IM)
!         PL = HP(IM)
!         RL = HR(IM)
!         B = HL(IM)
!         X = HCOMP(1,IM)
!         Z = HCOMP(3,IM)
!         DL = HD(IM)
!         FPL = FP(IM)
!         FTL = FT(IM)
!
!         IF(LMHD) THEN
!            CALL MEQOS(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,RMU,
!     *           AMU,EMU,ETA,QDT,QDP,QCP,DELA,QDTT,QDTP,QAT,QAP,QCPT,
!     *           QCPP,LDERIV,LATMO,KSAHA)
!         ELSE
!            IF (LDH) THEN
!               XXDH = HCOMP(1,IM)
!               YYDH = HCOMP(2,IM)+HCOMP(4,IM)
!               ZZDH = HCOMP(3,IM)
!               ZDH(1) = HCOMP(5,IM)+HCOMP(6,IM)
!               ZDH(2) = HCOMP(7,IM)+HCOMP(8,IM)
!               ZDH(3) = HCOMP(9,IM)+HCOMP(10,IM)+HCOMP(11,IM)
!            END IF
!            CALL EQSTAT(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,RMU,
!     *           AMU,EMU,ETA,QDT,QDP,QCP,DELA,QDTT,QDTP,QAT,QAP,QCPT,
!     *           QCPP,LDERIV,LATMO,KSAHA)
!         ENDIF
!         CALL GETOPAC (DL,TL,X,Z,O,OL,QOD,QOT,FXION)
!         IOVIM=IM
!         CALL TPGRAD(TL,T,PL,P,D,RL,SL,B,O,QDT,QDP,QOT,QOD,QCP,DEL,
!     *        DELR,DELA,QDTT,QDTP,QAT,QAP,QACT,QACP,QACR,QCPT,QCPP,
!     *        VEL,LDERIV,LCONV,FPL,FTL,TEFFL)
!         SDEL(1,IM) = DELR
!         SDEL(2,IM) = DEL
!         SDEL(3,IM) = DELA
!C JVS 10/13 Always want SVEL
!       SVEL(IM) = VEL
! 730  CONTINUE
!      endif

!       CALL PHYSIC(FP,FT,HCOMP,HD,HG,HL,HP,HR,HS,HT,LC,LCZ,M,TEFFL)  ! KC 2025-05-31
      call physic(pressure_rotation_factor,temperature_rotation_factor, &
           composition,log_density,mean_gravity,luminosity_lsun,log_pressure, &
           log_radius,log_mass,log_temperature,convective_flag,num_shells, &
           log_teff)
      call ovrot(composition,log_density,log_pressure,log_radius,log_mass, &
                 log_temperature,convective_flag,num_shells, &
                 am_transport_convective_flag,radiative_zone_bounds, &
                 convective_zone_bounds,num_radiative_zones, &
                 num_convective_zones)
! INITIALIZE TAUCZ, PPHOT, AND FRACSTEP
!       CALL GETTAU(HCOMP,HR,HP,HD,HG,HS1,HT,FP,FT,TEFFL,  ! KC 2025-05-31
      call gettau(composition,log_radius,log_pressure,log_density, &
                  enclosed_mass,log_temperature,pressure_rotation_factor, &
                  temperature_rotation_factor,log_teff, &
                  log_total_mass,log10_luminosity,num_shells,convective_flag, &
                  env_struct%env_log10_radius)
      turnover%convective_turnover_timescale_old = turnover%convective_turnover_timescale
      turnover%pphot0 = turnover%pphot
      turnover%fracstep = 0.5

      return
end subroutine starin
