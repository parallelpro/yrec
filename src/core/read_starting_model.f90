!----------------------------------------------------------------------
! read_starting_model (formerly starin)
!----------------------------------------------------------------------
! YREC's initial-model setup/read driver, called at the start of every
! kind card (and again after a diverged step, with the cut timestep).
! If a fresh starting model is requested (first_call_flag(run_index)
! true), it reads the model file (.mod, YREC7 or MODEL2 layout,
! detected from the file's keyword) via read_mod_model / read_yrec7 /
! read_model2, optionally extends the innermost shell inward
! (extend_core_inward), rescales it (rescale_model), optionally
! changes the envelope fitting point mass (calling atm_get / eos_get /
! kap_get / temperature_gradients to keep the new last interior
! shell's radiative/convective flag and density consistent), sets up
! the rotation curve (rotation_shape_factors / zone_moments_of_inertia)
! if rotation is active, sets up the surface mixture and opacity
! tables, and finally calls shell_physics / am_convective_regions /
! compute_turnover_timescale for the initial temperature-gradient /
! opacity / convective-zone structure. If instead the in-memory model
! is reused (first_call_flag(run_index) false), only the envelope-mass
! rescale (rescale_model) and the final physics calls are done.
!
! Naming: the star% arrays follow henyey_iterate's names for the same
! data (luminosity_lsun is the LINEAR L/Lsun, m/dm the unlogged shell
! masses, mean_radius/qiw/mean_gravity the rotation geometry); the
! per-point scalars of the single-shell physics block below follow
! shell_physics' log10_* naming.
subroutine read_starting_model(timestep_yr, delta_time, delta_time_abs, &
     trial_sign_flag, ikut_flag, istore_flag, model_failed_flag, &
     envelope_recomputed_flag, run_index, dlnrho_dlnp, dlnrho_dlnt, &
     total_angular_momentum, total_rotational_ke, convective_velocity, &
     species_mix_weights, ierr)
      use temperature_gradients_lib
      use star_info_lib, only: star, i_be9, i_c12, i_c13, i_h1, i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18, json, &
           n_species_basic, n_species_extended, n_mix_species, ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, ix_h, ix_o, ix_n, ix_ar, ix_ne, ix_he, &
           max_convective_zones, max_radiative_zones
      use controls_lib, only: ichi_dm_max, ichi_dp_env_max
      use envint_lib, only: atm_get
      use envstruct_lib
      use luout_lib
      use phys_const_lib
      use eos_lib
      use kap_lib
      use opacity_table_lib
      use yale_eos_lib
      use scv_eos_lib
      use stitched_model_lib, only: build_stitched_model
      implicit none
      integer :: jerr_atm

      double precision, intent(inout) :: timestep_yr
      double precision, intent(inout) :: delta_time, delta_time_abs
      double precision, intent(inout) :: trial_sign_flag
      integer, intent(out) :: ikut_flag
      integer, intent(inout) :: istore_flag
      logical, intent(in) :: model_failed_flag
      logical, intent(inout) :: envelope_recomputed_flag
      integer, intent(in) :: run_index
      double precision, intent(out) :: dlnrho_dlnp, dlnrho_dlnt
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke
      double precision, intent(out) :: convective_velocity
      double precision, intent(inout) :: species_mix_weights(n_mix_species)

      integer :: katm, kenv, saha_state
      character(len=4) :: format_tag
      character(len=6) :: eos_code
      character(len=4) :: atm_code, alok_code, hik_code
! (the initial-composition code read with the model is
! star%initial_composition_code)
! atomic weights of the surface mixture in ix_* slot order
! (Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He)
      double precision :: atomic_weight(n_mix_species)
      data atomic_weight /23.0d0,26.99d0,24.32d0,55.86d0,28.1d0,12.015d0, &
           1.008d0,16.0d0,14.01d0,39.96d0,20.19d0,4.004d0/

! ENSURE THAT ONLY HOMOGENEOUS MODELS HAVE THE MIXTURE ALTERED
      double precision :: reference_composition(n_species_extended)

! --- locals ---
      integer :: iread
      double precision :: mixing_length0
      logical :: mixing_length_matches
      logical :: mixture_ok
! NOTE: LEXCP0 in the original source is never assigned before use
! (line "WRITE(ISHORT,1040) CMIXL,CMIXL0,LEXCOM,LEXCP0" below); it
! appears to be a pre-existing bug (likely meant to read
! use_extended_composition0). Preserved exactly as in the original.
      logical :: lexcp0
      logical :: use_extended_composition0
      double precision :: fraction_diff
      integer :: i, j, k
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
      logical :: want_derivatives, in_atmosphere
      double precision :: hydrogen_fraction, metal_fraction
      double precision :: log10_pressure, log10_temperature, log10_density
      double precision :: log10_radius, log10_mass, shell_luminosity_lsun
      double precision :: point_pressure_rotation_factor, &
           point_temperature_rotation_factor
! 2026 named-index results for the convective-flag test of the new
! last shell (the envint/atm_get sites below keep their scalars).
      double precision :: eos_res(num_eos_results), kap_res(num_kap_results)
      double precision :: radiative_gradient, dgrad_dt_component, &
           dgrad_dp_component, dgrad_dr_component
      logical :: is_convective
      integer :: old_last_shell
      double precision :: saved_env_step_max, saved_env_step_min, &
           saved_env_step_begin
      logical :: save_boundary_flag, print_flag
      double precision :: log10_gravity
      integer :: vertex_index
      double precision :: log10_pressure_limit
      integer :: atm_get_unused_flag
      double precision :: spot_adjusted_log_teff
      double precision :: atm_get_dummy1(4), atm_get_dummy2(3), &
           atm_get_dummy3(3), atm_get_dummy4(3)
      double precision :: pressure_offset, density_offset, &
           temperature_offset, radius_offset
      integer :: env_point_index
      double precision :: lower_mass_coord, target_mass_coord, &
           upper_mass_coord, envelope_interp_fraction
      double precision :: old_senv, target_log_mass_at_fit
      double precision :: interior_interp_fraction
      integer :: num_species
      double precision :: angular_momentum_sum, rotational_ke_sum, &
           shell_angular_momentum
      double precision :: mixture_weight_sum, mixture_scale_factor
      integer :: radiative_zone_bounds(max_radiative_zones,2), &
           convective_zone_bounds(max_convective_zones,2)
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
! If flag LARGE is set, model has failed to converge.  Terminate the run.
      ! 2026 (ROADMAP.md stage 3): library errors return here via ierr;
      ! this driver-side call site preserves the historical stop.
      integer :: jerr

      integer, intent(out) :: ierr

      ierr = 0

      if (model_failed_flag) then
          write(run_log_unit,1000)
 1000       format(1x,39('>'),40('<')/, &
                  "STARIN:        ***** RUN STOPPED *****")
          write(run_log_unit,1010)
 1010       format("STARIN: ***** MODEL FAILED TO CONVERGE *****")
            ! 2026 (phase five, step B): stop converted to ierr; run_yrec
            ! returns the error and the CLI wrapper (main) stops.
            ierr = 1
            return
      endif

! THIS SUBROUTINE READS IN THE INITIAL STELLAR MODEL
! INITIAL MODEL IS STORED IN LOGICAL UNIT IFIRST

! Set flags for reading an input model
      ikut_flag = 0
      iread = star%ctrl%first_unit

! INITIALIZE VARIABLES
      convective_velocity = 0.0d0
      dlnrho_dlnt = -1.0d0
      dlnrho_dlnp = 1.0d0

      call acquire_starting_model
      if (ierr /= 0) return
      call extend_core_toward_center
      if (ierr /= 0) return
      call rescale_and_refit_envelope
      if (ierr /= 0) return
      call build_shell_masses
      call initialize_rotation_state
      if (ierr /= 0) return
      call update_surface_mixture
      call snapshot_step_start
      if (ierr /= 0) return
      call shell_physics(star%fp_rot,star%ft_rot, &
           star%xa,star%logRho,star%mean_gravity,star%luminosity_lsun,star%logP, &
           star%logR,star%log_mass,star%logT,star%convective_flag,star%nz, &
           star%log_Teff, jerr)
      if (jerr /= 0) then
      ! 2026 (phase five, step B): propagate instead of stopping
         ierr = jerr
         return
      end if
      call am_convective_regions(star%xa,star%logRho,star%logP,star%logR,star%log_mass, &
                 star%logT,star%convective_flag,star%nz, &
                 am_transport_convective_flag,radiative_zone_bounds, &
                 convective_zone_bounds,num_radiative_zones, &
                 num_convective_zones, ierr)
      if (ierr /= 0) return
! INITIALIZE TAUCZ, PPHOT, AND FRACSTEP (2026: build the stitched
! model first -- the turnover walker reads it; this also seeds
! star%pphot. The historical call passed env_log10_radius as the
! scalar BCZ-radius out-arg, clobbering env_log10_radius(1) by
! sequence association; the sane out-arg is used now.)
      call build_stitched_model
      call compute_turnover_timescale(star%envelope_radius, ierr)
      if (ierr /= 0) return
      star%convective_turnover_timescale_old = star%convective_turnover_timescale
      star%pphot0 = star%pphot
      star%fracstep = 0.5

      return

contains

! ---------------------------------------------------------------
! Obtain the starting model: reuse the model in memory (rescaling
! the envelope mass when requested), or read + process the input
! model file (YREC7 or MODEL2 layout, detected from its keyword),
! check its mixing length / surface BC / CZ settings against the
! user parameters, set up the rotation curve, and -- first model
! of a run only -- apply the heavy-element mixture / isotope-ratio
! alterations (isetmix / isetiso). Sets ierr on failure.
subroutine acquire_starting_model
! Flag LFIRST(NK) tells where to get the starting stellar model for the current
! step (step NK).  If LFIRST(NK) is true, read in the starting stellarmodel from
! the file specified by LU IFIRST.  If LFIRST(NK) is false, as starting model use
! the stellar model currently stored in memory.

      if (.not.star%job%first_call_flag(run_index)) then
! Use the model currently in memory as the starting model.
! DBG 2/92 CHANGED SO WILL RESCALE ENVELOPE MASS ON EACH NEW RUN
         if (star%job%rescale_kind(run_index).ne.1) call rescale_model(star%luminosity_lsun, &
              star%xa,star%log_mass,star%log_total_mass,star%nz,run_index, &
              star%star_mass,star%convective_flag, ierr)
         if (ierr /= 0) return
! The reading and processing of an input model file is skipped
! (the else branch below).
      else


!     Read in the starting model from LU IFIRST and process it.

!  Get file format of input model

      rewind iread
      read(iread,10) format_tag
 10   format(a4)

! ATEMP now contains a keyword describing the format of the input stellar
! model.  We decide what kind of model format it has and process accordingly.

      if (format_tag .eq. 'YMOD') then
         write(run_log_unit,11)
 11      format(1x,'loading starting model (YREC .mod format)')
         call read_mod_model(iread, timestep_yr, mixing_length0, &
              use_extended_composition0, rotation_active0, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      else if (format_tag .eq. 'NMOD') then
         write(run_log_unit,12)
 12      format(1x,'loading starting model (YREC7 format)')
! (2026 de-tramp: the star%-shaped slots are read/written by
! read_yrec7 directly; only genuine locals -- the model file's own
! header values, latched below -- remain arguments.)
         call read_yrec7(mixing_length0,timestep_yr,trial_sign_flag,iread, &
              core_cz_top_index0,envelope_cz_bottom_index0, &
              use_extended_composition0,rotation_active0, &
              atm_code,eos_code,hik_code,use_diffusion_y0,use_diffusion_z0, &
              disk_locking_active0,instability_transport_active0,ljdot00, &
              alok_code,lovstc0,envelope_overshoot_active0,lovstm0, &
              use_pure_z_table0,lsemic0,disk_pressure0, &
              disk_temperature0,wind_saturation_omega0, ierr)
         if (ierr /= 0) return
! First three lines above are YREC7 inputs
! Last three lines are MODEL2 add-ons

      else if (format_tag .eq. 'MOD2 ') then
         write(run_log_unit,16)
 16      format(1x,'loading starting model (MODEL2 format)')
! (2026 de-tramp: same treatment as the read_yrec7 call above.)
         call read_model2(mixing_length0,timestep_yr,trial_sign_flag,iread, &
              core_cz_top_index0,envelope_cz_bottom_index0, &
              use_extended_composition0,rotation_active0, &
              atm_code,eos_code,hik_code,use_diffusion_y0,use_diffusion_z0, &
              disk_locking_active0,instability_transport_active0,ljdot00, &
              alok_code,lovstc0,envelope_overshoot_active0,lovstm0, &
              use_pure_z_table0,lsemic0,disk_pressure0, &
              disk_temperature0,wind_saturation_omega0)
! First three lines above are YREC7 inputs
! Last three lines are MODEL2 add-ons

      else
         write(run_log_unit,20)
 20      format('STARIN: ***** RUN TERMINATED, INVALID INPUT', &
                ' MODEL FILE.  *****')
         ! 2026 (phase five, step B): stop converted to ierr; run_yrec
         ! returns the error and the CLI wrapper (main) stops.
         ierr = 1
         return
      endif

! Model has now been read in. Some post-processing is required.

      delta_time = seconds_per_year*timestep_yr
      delta_time_abs = dabs(delta_time)
      star%stotal = star%log_total_mass

! CHECK TO ENSURE THAT MIX LENGTH, SURFACE B.C. AND CONVECTION ZONE
! THEORY OF THE MODEL ARE THE SAME AS USER PARAMETERS
! CHECKED ONLY FOR EVOLVED MODELS(MODEL NUMBER > 0)
! 1/92 Changed to not stop just give warning.
      if (star%model_number.gt.0) then
       mixing_length_matches = .true.
! Jan 12, 1989 : in STARIN changed test to see if model has correct mixing
! length from 1.0e-6 to 2.0e-3 because models only store
! mixing length to four sig digits.
       if (mixing_length0.gt.0.0d0) mixing_length_matches = &
            (dabs(star%mixing_length_alpha-mixing_length0).lt.2.0d-3)
! MHP 9/03 FIXED TYPO
       if (.not.mixing_length_matches .or. use_extended_composition0.neqv. &
            star%job%use_extended_composition) then
          write(run_log_unit,1040) star%mixing_length_alpha,mixing_length0, &
               star%job%use_extended_composition,lexcp0
          write(terminal_unit,1040) star%mixing_length_alpha,mixing_length0,star%job%use_extended_composition, &
               lexcp0
 1040       format(1x,'ERROR IN SUBROUTINE STARIN'/1x,'USER PARAMETERS', &
             ' OF WRONG TYPE FOR INITIAL MODEL'/1x,'MIXING LENGTH - USER' &
             ,' DESIRES',f7.3,' MODEL MIX LENGTH',f7.3/1x, &
             'EXTENDED COMP-USER DESIRES ',l1,' MODEL USED ',l1)
! (the run is not stopped: only this warning is issued)
       endif
      endif

      istore_flag = 0
      trial_sign_flag = 1.0d0
! Require recompute of envelope. No assurance of reliable
! input triangle
       envelope_recomputed_flag = .true.


! GET XNEW AND ZNEW FROM HENYEY POINTS

      star%xnew = star%xa(i_h1,star%nz)
      star%znew = star%xa(i_metals,star%nz)

! FOURTH PART:  - LOG J/M STORED

      if (star%job%rotation_active) then
       if (star%ctrl%set_initial_omega) then
! GENERATE A SOLID BODY ROTATION CURVE WITH OMEGA = WNEW;
! THIS IS DONE TO CONVERT A NON-ROTATING MODEL TO A ROTATING ONE.
          do i = 1,star%nz
             star%omega(i) = star%ctrl%initial_omega
          end do
       endif
      else
         do i = 1,json
            star%fp_rot(i) = 1.0d0
            star%ft_rot(i) = 1.0d0
         end do
      endif
! KEEP IREAD OPEN
      rewind iread
! End of the reading and processing of an input model file.
      if (star%job%first_call_flag(run_index)) then
!      IF(.NOT.LFIRST(NK).OR.NK.GT.1)GOTO 3000
!     MHP 10/24 MACHINERY TO ALTER THE HEAVY ELEMENT MIXTURE
!     THIS IS ONLY DONE if the first MODEL IS being READ IN, AND ONLY FOR A
! CHEMICALLY HOMOGENEOUS MODEL. IT CAN OVER-WRITE MASS FRACTIONS 4-15 WITH USER-SPECIFIED VALUES
! ISETMIX=1 -> CAN ADJUST CNO FRACTIONS ISETISO=1-> CHANGE ISOTOPE RATIOS
      mixture_ok = .true.
      if (star%ctrl%change_cno_mixture_active .or. star%ctrl%change_isotope_ratios_active) then
! ENSURE STARTING MODEL IS HOMOGENEOUS BEFORE EITHER IS CHANGED
         do i = 1,15
            reference_composition(i)=star%xa(i,1)
         end do
         homogeneity: do j = 2,star%nz
            do i = 1,15
               fraction_diff = abs(star%xa(i,j)-reference_composition(i))
               if (fraction_diff.gt.1.0d-6) then
                  write(*,592)i,j,fraction_diff
                  write(run_log_unit,592)i,j,fraction_diff
 592              format('SPECIES ',i3,' IN SHELL ',i5, &
                    ' DIFFERS FROM CENTER BY ',e12.4, &
                    ' MIX NOT MODIFIED IN EVOLVED MODEL')
                  mixture_ok = .false.
                  exit homogeneity
               endif
            end do
         end do homogeneity
      endif
      if (mixture_ok) then
! LOOP FOR CHANGING CNO MIX
      if (star%ctrl%change_cno_mixture_active) then
!     INFER CURRENT TOTAL CNO FRACTIONS AND SCALE ALL ISOTOPES BY THE RATIO BETWEEN
!     DESIRED AND CURRENT FRACTIONS. RELATIVE ISOTOPES ARE ADJUSTED IN THE ISOTOPE SECTION BELOW.
         total_carbon_cno_fraction = (reference_composition(i_c12)+ &
              reference_composition(i_c13))/reference_composition(i_metals)
         total_nitrogen_cno_fraction = (reference_composition(i_n14)+ &
              reference_composition(i_n15))/reference_composition(i_metals)
         total_oxygen_cno_fraction = (reference_composition(i_o16)+ &
              reference_composition(i_o17)+reference_composition(i_o18))/ &
              reference_composition(i_metals)
         carbon_scale_ratio = star%ctrl%target_carbon_cno_fraction/ &
              total_carbon_cno_fraction
         nitrogen_scale_ratio = star%ctrl%target_nitrogen_cno_fraction/ &
              total_nitrogen_cno_fraction
         oxygen_scale_ratio = star%ctrl%target_oxygen_cno_fraction/ &
              total_oxygen_cno_fraction
         do i = i_c12,i_c13
            do j = 1,star%nz
               star%xa(i,j)=carbon_scale_ratio*star%xa(i,j)
            end do
         end do
         do i = i_n14,i_n15
            do j = 1,star%nz
               star%xa(i,j)=nitrogen_scale_ratio*star%xa(i,j)
            end do
         end do
         do i = i_o16,i_o18
            do j = 1,star%nz
               star%xa(i,j)=oxygen_scale_ratio*star%xa(i,j)
            end do
         end do
         write(*,594)(reference_composition(k),k=i_c12,i_o18), &
              (star%xa(k,1),k=i_c12,i_o18)
         write(run_log_unit,594)(reference_composition(k),k=i_c12,i_o18), &
              (star%xa(k,1),k=i_c12,i_o18)
 594     format(1x,'CNO mixture applied to the starting model:',/, &
             4x,'old (C12 C13 N14 N15 O16 O17 O18):',7es12.4,/, &
             4x,'new (C12 C13 N14 N15 O16 O17 O18):',7es12.4)
      endif
! DESIRED ISOTOPE RATIOS AND LIGHT ELEMENT ABUNDANCES ASSIGNED.
!     AT PRESENT B10,B11,N15,O17 ARE NOT USED AND THUS NOT ALTERED.
!     START WITH LIGHT ELEMENTS
      if (star%ctrl%change_isotope_ratios_active) then
         sum_c12_c13 = star%xa(i_c12,1)+star%xa(i_c13,1)
         sum_o16_o18 = star%xa(i_o16,1)+star%xa(i_o18,1)
         do j = 1,star%nz
            star%xa(i_he3,j)=star%ctrl%initial_he3_fraction
            star%xa(i_c12,j)= star%ctrl%c12_to_c13_ratio*sum_c12_c13/ &
                 (1.0d0+star%ctrl%c12_to_c13_ratio)
            star%xa(i_c13,j)= sum_c12_c13/(1.0d0+star%ctrl%c12_to_c13_ratio)
            star%xa(i_o16,j)= star%ctrl%o16_to_o18_ratio*sum_o16_o18/ &
                 (1.0d0+star%ctrl%o16_to_o18_ratio)
            star%xa(i_o18,j)= sum_o16_o18/(1.0d0+star%ctrl%o16_to_o18_ratio)
            star%xa(i_h2,j)=star%ctrl%initial_h2_fraction
            star%xa(i_li6,j)=star%ctrl%initial_li6_fraction
            star%xa(i_li7,j)=star%ctrl%initial_li7_fraction
            star%xa(i_be9,j)=star%ctrl%initial_be9_fraction
         end do
         write(*,593)(reference_composition(k),k=i_he3,i_be9), &
              (star%xa(k,1),k=i_he3,i_be9)
         write(run_log_unit,593)(reference_composition(k),k=i_he3,i_be9), &
              (star%xa(k,1),k=i_he3,i_be9)
 593     format(1x,'isotope/light-element mixture applied to the', &
             ' starting model:',/, &
             4x,'old (He3 C12 C13 N14 N15 O16 O17 O18 H2 Li6 Li7 Be9):', &
             12es12.4,/, &
             4x,'new (He3 C12 C13 N14 N15 O16 O17 O18 H2 Li6 Li7 Be9):', &
             12es12.4)
      endif
      end if
      end if
      endif
end subroutine acquire_starting_model

! ---------------------------------------------------------------
! If extend_core_inward is set, add central points inward of the
! innermost shell using constant-epsilon, constant-density
! stellar-structure estimates (spacing per hpttol).
subroutine extend_core_toward_center
      use math_lib

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
      if (star%ctrl%extend_core_inward .and. star%job%first_call_flag(run_index)) then
!      IF(LCORE .AND. NK .EQ. 1) THEN
! AVOID SHUFFLING POINTS BY ASSIGNING NEW CENTRAL POINTS IN INTEGER
! MULTIPLES OF THE CENTRAL POINT SPACING.
!     MCORE is number of shells to extrapolate to new core.
!     FCORE is factor to reduce inner mass shell.
          star%job%num_core_shells_added = int(log10(star%job%core_mass_reduction_factor)/ &
               star%ctrl%chi_grid_scale(ichi_dm_max))+1
          star%job%core_mass_reduction_factor = dble(star%job%num_core_shells_added)* &
               star%ctrl%chi_grid_scale(ichi_dm_max)
          num_shells_extended = star%nz + star%job%num_core_shells_added
          if (num_shells_extended .gt. json) then
             write(run_log_unit,476)"STARIN: Unable to extend core inward ", &
                   "- JSON too small"
             write(run_log_unit,477) "STARIN: Required size =", &
                    num_shells_extended, ", JSON = ", json
             write(run_log_unit,478) "STARIN:  ***** RUN TERMINATED *****"
  476        format(2a)
  477        format(a, i8, a, i8)
  478        format(a)
          endif
          core_shell_spacing = star%ctrl%chi_grid_scale(ichi_dm_max)
! shift data for remaining points by the required number
          do i=star%nz,1, -1
             star%log_mass(i+star%job%num_core_shells_added) = star%log_mass(i)
             star%logR(i+star%job%num_core_shells_added) = star%logR(i)
             star%luminosity_lsun(i+star%job%num_core_shells_added) = star%luminosity_lsun(i)
             star%logP(i+star%job%num_core_shells_added) = star%logP(i)
             star%logT(i+star%job%num_core_shells_added) = star%logT(i)
             star%convective_flag(i+star%job%num_core_shells_added) = star%convective_flag(i)
             do j=1, 15
                star%xa(j,i+star%job%num_core_shells_added) = star%xa(j,i)
             end do
             star%omega(i+star%job%num_core_shells_added) = star%omega(i)
          end do
          first_original_shell = star%job%num_core_shells_added+1
! MARCH INWARD IN MASS FROM THE INNERMOST MODEL POINT.
! ASSUME EPSILON=CONSTANT AND DEL=CONSTANT
          central_log_density = star%logRho(1)
          actual_gradient = (star%logT(2)-star%logT(1))/ &
               (star%logP(2)-star%logP(1))
          central_shell_luminosity = star%luminosity_lsun(1)
! MHP 4/12 FACTOR FOR ESTIMATING RHO FROM P AND T
          density_estimate_offset = star%logP(first_original_shell)- &
               star%logRho(first_original_shell)- &
               star%logT(first_original_shell)
          do i = star%job%num_core_shells_added,1,-1
             star%log_mass(i) = star%log_mass(i+1)-core_shell_spacing
! USE M  = 4PI/3*RHOC*R**3 TO GET R AS A FUNCTION OF M
             star%logR(i) = cc13*(star%log_mass(i)-c4pi3l-central_log_density)
! USE EPSILON=CONSTANT TO GET L
             star%luminosity_lsun(i) = exp(ln10*(star%log_mass(i)- &
                  star%log_mass(first_original_shell)))*central_shell_luminosity
! USE HYDROSTATIC EQUILIBRIUM TO GET P
             trial_log_pressure = star%logP(i+1)
             temp_scratch = exp(ln10*(cgl+2.0d0*star%log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*star%logR(i)))
             trial_log_pressure = trial_log_pressure+0.5d0*temp_scratch* &
                  core_shell_spacing
             temp_scratch = exp(ln10*(cgl+2.0d0*star%log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*star%logR(i)))
             star%logP(i) = star%logP(i+1)+temp_scratch* &
                  core_shell_spacing
! ASSUME R/C FLAG IS THE SAME AS FOR THE FIRST POINT
             star%convective_flag(i) = star%convective_flag(first_original_shell)
! ASSUME DEL= CONSTANT IN THE CORE
             star%logT(i) = star%logT(i+1)+temp_scratch* &
                  actual_gradient*core_shell_spacing
! ASSUME OMEGA = CONSTANT
             star%omega(i) = star%omega(first_original_shell)
! ASSUME COMPOSITION IS UNIFORM
             do j=1, 15
                star%xa(j,i) = star%xa(j,first_original_shell)
             end do
! LOCAL ESTIMATE FOR RHO (replaced a broken equation-of-state call)
             star%logRho(i) = star%logP(i) - star%logT(i) - &
                  density_estimate_offset
          end do
          star%nz = num_shells_extended
      end if
! End of code to extend core inward
end subroutine extend_core_toward_center

! ---------------------------------------------------------------
! Rescale the first model when the kind card asks for it, then move
! the envelope fitting point to the requested envelope mass --
! deleting points (deeper new envelope) or integrating a fresh
! envelope down to the new fitting mass (shallower), with EOS
! re-evaluation at the surface point. Sets ierr on failure.
subroutine rescale_and_refit_envelope
      use math_lib

! PERFORM RESCALING OF FIRST MODEL IF REQUIRED
      if (star%job%rescale_kind(run_index).ne.1) call rescale_model(star%luminosity_lsun, &
           star%xa,star%log_mass,star%log_total_mass,star%nz,run_index, &
           star%star_mass,star%convective_flag, ierr)
      if (ierr /= 0) return

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
      if (star%job%change_envelope_mass_flag) then
      envelope_rescale: do
       if (star%job%requested_envelope_mass.gt.0.0d0) star%job%requested_envelope_mass = &
            -star%job%requested_envelope_mass
! DBG 2/92 CHANGED MINIMUM FROM 1.0D-10 TO 1.0D-12
! RESTRICT MIMIMUM ENVELOPE MASS;1.0D-12 CORRESPONDS TO TAU=2/3
! FOR THIS PURPOSE(BASE OF ATMOSPHERE).
       if (star%job%requested_envelope_mass.gt.-1.0d-12) star%job%requested_envelope_mass = &
            -1.d-12
       star%senv = star%log_mass(star%nz) - star%log_total_mass
       old_senv = star%senv
       if (star%senv.eq.star%job%requested_envelope_mass) exit envelope_rescale
       num_species = n_species_basic
       if (star%job%use_extended_composition) num_species = n_species_extended
       if (star%job%requested_envelope_mass.lt.star%senv) then
! NEW ENVELOPE DEEPER THAN THE OLD ONE
          target_log_mass_at_fit = star%log_total_mass+star%job%requested_envelope_mass
          do i = star%nz-1,1,-1
             if (star%log_mass(i).lt.target_log_mass_at_fit) exit
          end do
          if (i < (1)) then
! ENVELOPE MASS DESIRED WITHIN FIRST POINT;PRINT NASTY MESSAGE
! AND ABORT.
          write(run_log_unit,576)star%job%requested_envelope_mass
 576        format(5x,'ERROR IN SUBROUTINE STARIN'/5x,'DESIRED', &
              ' ENVELOPE MASS',1pe22.13,' TOO LARGE'/5x,'ENVELOPE', &
              ' MASS NOT CHANGED')
          exit envelope_rescale
          end if
            star%nz = i + 1
          star%senv = star%job%requested_envelope_mass
          interior_interp_fraction = (target_log_mass_at_fit-star%log_mass(i))/ &
               (star%log_mass(i+1) - star%log_mass(i))
          star%log_mass(star%nz) = target_log_mass_at_fit
          star%logRho(star%nz) = star%logRho(i) + &
               interior_interp_fraction*(star%logRho(i+1) - star%logRho(i))
          star%luminosity_lsun(star%nz) = star%luminosity_lsun(i) + &
               interior_interp_fraction*(star%luminosity_lsun(i+1) - &
               star%luminosity_lsun(i))
          star%logP(star%nz) = star%logP(i) + &
               interior_interp_fraction*(star%logP(i+1) - star%logP(i))
          star%logR(star%nz) = star%logR(i) + &
               interior_interp_fraction*(star%logR(i+1) - star%logR(i))
          star%logT(star%nz) = star%logT(i) + &
               interior_interp_fraction*(star%logT(i+1) - &
               star%logT(i))
          do j = 1,num_species
             star%xa(j,star%nz) = star%xa(j,i)
          end do
          star%xnew = star%xa(i_h1,star%nz)
          star%znew = star%xa(i_metals,star%nz)
          if (star%job%rotation_active) star%omega(star%nz) = star%omega(i) + &
               interior_interp_fraction*(star%omega(i+1)-star%omega(i))
          if (star%convective_flag(i).and.star%convective_flag(i+1)) then
             star%convective_flag(star%nz) = .true.
          else if (.not.star%convective_flag(i).and. .not.star%convective_flag(i+1)) &
               then
             star%convective_flag(star%nz) = .false.
          else
! CALL BASIC PHYSICS ROUTINES TO DETERMINE IF THE NEW LAST SHELL IS
! CONVECTIVE OR RADIATIVE.
             want_derivatives = .false.
             in_atmosphere = .true.
             saha_state = 0
             hydrogen_fraction = star%xa(i_h1,star%nz)
             metal_fraction = star%xa(i_metals,star%nz)
             log10_pressure = star%logP(star%nz)
             log10_temperature = star%logT(star%nz)
             log10_density = star%logRho(star%nz)
             shell_luminosity_lsun = star%luminosity_lsun(star%nz)
             log10_radius = star%logR(star%nz)
             log10_mass = star%log_mass(star%nz)
             point_pressure_rotation_factor = 1.0d0
             point_temperature_rotation_factor = 1.0d0
               eos_res(i_log10_density) = log10_density
               call eos_get(log10_temperature, log10_pressure, &
                    hydrogen_fraction, metal_fraction, eos_res, &
                    want_derivatives, in_atmosphere, saha_state, &
                    composition_at_zone=star%xa(:,star%nz), ierr=ierr)
               if (ierr /= 0) return
! kap at eqstat's returned density -- the historical inout dataflow
               call kap_get(eos_res(i_log10_density), log10_temperature, &
                    hydrogen_fraction, metal_fraction, kap_res, &
                    eos_res(i_fxion:i_fxion+2), ierr=ierr)
               if (ierr /= 0) return
               star%iovim = -1
               call temperature_gradients(log10_temperature, log10_pressure, &
                    eos_res, kap_res, log10_radius, log10_mass, &
                    shell_luminosity_lsun, actual_gradient, radiative_gradient, &
                    dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
                    convective_velocity, want_derivatives, is_convective, &
                    point_pressure_rotation_factor, &
                    point_temperature_rotation_factor, star%log_Teff, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
               log10_density = eos_res(i_log10_density)
               star%logRho(star%nz) = log10_density
               star%convective_flag(star%nz) = is_convective
          endif
       else
! DESIRED ENVELOPE MASS LESS THAN CURRENT VALUE.
            old_last_shell = star%nz
            saved_env_step_max = star%job%env_step_max
            saved_env_step_min = star%job%env_step_min
            saved_env_step_begin = star%job%env_step_begin
            star%job%env_step_max = star%ctrl%chi_grid_scale(ichi_dp_env_max)
            star%job%env_step_min = star%ctrl%chi_grid_scale(ichi_dp_env_max)
            star%job%env_step_begin = star%ctrl%chi_grid_scale(ichi_dp_env_max)
!          SENV = SENV0
          save_boundary_flag = .false.
          print_flag = .true.
          katm = 0
          kenv = 0
          saha_state = 0
          shell_luminosity_lsun = exp(ln10*star%log_L)
          log10_radius = 0.5d0*(star%log_L + star%solar_luminosity_cgs - &
               4.0d0*star%log_Teff - c4pil - csigl)
          log10_gravity = cgl + star%stotal - log10_radius - log10_radius
          hydrogen_fraction = star%xa(i_h1,star%nz)
          metal_fraction = star%xa(i_metals,star%nz)
          point_pressure_rotation_factor = 1.0d0
          point_temperature_rotation_factor = 1.0d0
          vertex_index=0
          log10_pressure_limit = star%logP(star%nz)
! DBG PULSE: DO NOT DO PULSE OUTPUT
            if (use_debye_huckel_correction) then
               debye_huckel_x = star%xa(i_h1,star%nz)
               debye_huckel_y = star%xa(i_he4,star%nz)+star%xa(i_he3,star%nz)
               debye_huckel_z_total = star%xa(i_metals,star%nz)
               debye_huckel_z(1) = star%xa(i_c12,star%nz)+star%xa(i_c13,star%nz)
               debye_huckel_z(2) = star%xa(i_n14,star%nz)+star%xa(i_n15,star%nz)
               debye_huckel_z(3) = star%xa(i_o16,star%nz)+star%xa(i_o17,star%nz)+ &
                    star%xa(i_o18,star%nz)
            end if
! MHP 10/02  define ISTORE - used in ENVINT
            atm_get_unused_flag = 0
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
          if (star%envelope_cz_bottom_index.eq.star%nz.and.star%ctrl%spot_filling_factor.ne. &
               0.0.and.star%ctrl%spot_temp_contrast.ne.1.0) then
               spot_adjusted_log_teff = star%log_Teff - 0.25*log10(&
                    star%ctrl%spot_filling_factor * pow(star%ctrl%spot_temp_contrast, 4.0) + 1.0 - &
                    star%ctrl%spot_filling_factor)
          else
             spot_adjusted_log_teff = star%log_Teff
          endif
          call atm_get(shell_luminosity_lsun,point_pressure_rotation_factor, &
                 point_temperature_rotation_factor,log10_gravity,star%stotal, &
                 vertex_index,print_flag,save_boundary_flag, &
                 log10_pressure_limit,log10_radius,spot_adjusted_log_teff, &
                 hydrogen_fraction,metal_fraction,atm_get_dummy1, &
                 atm_get_unused_flag,katm,kenv,saha_state,atm_get_dummy2, &
                 atm_get_dummy3,atm_get_dummy4,ierr=jerr_atm)
! 2026 numerics-gate opt-in: envelope/atmosphere integration failures
! surface here (incl. numerics_termination) instead of stopping in
! atm_get; the host's ierr is already threaded to run_yrec.
          if (jerr_atm /= 0) then
             ierr = jerr_atm
             return
          end if
! G Somers END
            star%job%env_step_max = saved_env_step_max
            star%job%env_step_min = saved_env_step_min
            star%job%env_step_begin = saved_env_step_begin
          star%senv = star%job%requested_envelope_mass
! 2026 io/driver stops -> ierr (was: stop 9999)
            if (star%nz+env_struct%num_env_points.ge.json) then
               write(*,*) 'read_starting_model: grid overflow (json too small)'
               ierr = 1
               return
            end if
! ENFORCE CONSISTENCY WITH THE INTERIOR SOLUTION;
! ADJUST THE (P, RHO, T, R) POINTS TO BE CONSISTENT
! WITH THE LAST MODEL POINT.
            pressure_offset = star%logP(star%nz) - env_struct%env_log10_pressure(1)
            density_offset = star%logRho(star%nz) - env_struct%env_log10_density(1)
            temperature_offset = star%logT(star%nz) - &
                 env_struct%env_log10_temperature(1)
            radius_offset = star%logR(star%nz) - env_struct%env_log10_radius(1)
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
            do j = star%nz+1,star%nz+env_struct%num_env_points
               env_point_index = j-star%nz
! LUMINOSITY ASSUMED CONSTANT
               star%luminosity_lsun(j) = star%luminosity_lsun(star%nz)
! INCLUDE NEW POINTS UP TO THE DIFFERENT DESIRED FITTING POINT
               if (env_struct%env_log10_mass(env_point_index).le.star%senv) then
                  star%logRho(j) = env_struct%env_log10_density(env_point_index)
                  star%logP(j) = env_struct%env_log10_pressure(env_point_index)
                  star%logR(j) = env_struct%env_log10_radius(env_point_index)
                  star%log_mass(j) = env_struct%env_log10_mass(env_point_index) + star%stotal
                  star%logT(j) = env_struct%env_log10_temperature(env_point_index)
                  star%xa(i_h1,j) = env_struct%env_hydrogen_fraction(env_point_index)
                  star%xa(i_metals,j) = env_struct%env_metal_fraction(env_point_index)
                  do k = 4,num_species
                     star%xa(k,j) = star%xa(k,star%nz)
                  end do
                  star%xa(i_he4,j)=1.0d0-star%xa(i_h1,j)-star%xa(i_metals,j)- &
                       star%xa(i_he3,j)
                  star%convective_flag(j) = env_struct%env_convective_flag(env_point_index)
               else
! POINTS BEYOND THIS ARE ABOVE THE NEW DESIRED FITTING POINT;
! INTERPOLATE LINEARLY, SET NEW NUMBER OF TOTAL POINTS, AND EXIT
                  if (env_point_index.eq.1) then
! INTERPOLATE BETWEEN THE LAST INTERIOR POINT AND THE FIRST ENVELOPE POINT
                     lower_mass_coord = star%log_mass(star%nz)
                     target_mass_coord = star%stotal + star%senv
                     upper_mass_coord = env_struct%env_log10_mass(env_point_index) + &
                          star%stotal
                     if (upper_mass_coord-lower_mass_coord.lt.1.0d-14) then
! 2026 (phase five, step B): was a bare `stop 9998`; now prints and
! returns the error (run_yrec propagates, the CLI wrapper stops).
                        write(run_log_unit,*) 'STARIN: degenerate', &
                             ' envelope interpolation interval (was STOP 9998)'
                        ierr = 1
                        return
                     end if
                     envelope_interp_fraction = (target_mass_coord- &
                          lower_mass_coord)/(upper_mass_coord- &
                          lower_mass_coord)
                     star%logRho(j) = star%logRho(star%nz)+ &
                          envelope_interp_fraction*(env_struct%env_log10_density( &
                          env_point_index)-star%logRho(star%nz))
                     star%logP(j) = star%logP(star%nz)+ &
                          envelope_interp_fraction*(env_struct%env_log10_pressure( &
                          env_point_index)-star%logP(star%nz))
                     star%logR(j) = star%logR(star%nz)+ &
                          envelope_interp_fraction*(env_struct%env_log10_radius( &
                          env_point_index)-star%logR(star%nz))
                     star%log_mass(j) = target_mass_coord
                     star%logT(j) = star%logT(star%nz)+ &
                          envelope_interp_fraction*(env_struct%env_log10_temperature( &
                          env_point_index)-star%logT(star%nz))
                     star%xa(i_h1,j) = star%xa(i_h1,star%nz)+ &
                          envelope_interp_fraction*(star%xa(i_h1,star%nz) &
                          -env_struct%env_hydrogen_fraction(env_point_index))
                     star%xa(i_metals,j) = star%xa(i_metals,star%nz)+ &
                          envelope_interp_fraction*(star%xa(i_metals,star%nz) &
                          -env_struct%env_metal_fraction(env_point_index))
                     do k = 4,num_species
                        star%xa(k,j) = star%xa(k,star%nz)
                     end do
                     star%xa(i_he4,j)=1.0d0-star%xa(i_h1,j)- &
                          star%xa(i_metals,j)-star%xa(i_he3,j)
                     if (env_struct%env_convective_flag(env_point_index).or. &
                          star%convective_flag(star%nz)) then
                        star%convective_flag(j) = .true.
                     else
                        star%convective_flag(j) = .false.
                     endif
                  else
! INTERPOLATE BETWEEN THE LAST 2 ENVELOPE POINTS
                     lower_mass_coord = env_struct%env_log10_mass(env_point_index-1) + &
                          star%stotal
                     target_mass_coord = star%stotal + star%senv
                     upper_mass_coord = env_struct%env_log10_mass(env_point_index) + &
                          star%stotal
                     if (upper_mass_coord-lower_mass_coord.lt.1.0d-14) then
! 2026 (phase five, step B): was a bare `stop 9998`; now prints and
! returns the error (run_yrec propagates, the CLI wrapper stops).
                        write(run_log_unit,*) 'STARIN: degenerate', &
                             ' envelope interpolation interval (was STOP 9998)'
                        ierr = 1
                        return
                     end if
                     envelope_interp_fraction = (target_mass_coord- &
                          lower_mass_coord)/(upper_mass_coord- &
                          lower_mass_coord)
                     star%logRho(j) = env_struct%env_log10_density(env_point_index-1)+ &
                          envelope_interp_fraction*(env_struct%env_log10_density( &
                          env_point_index)-env_struct%env_log10_density(env_point_index-1))
                     star%logP(j) = env_struct%env_log10_pressure(env_point_index-1)+ &
                          envelope_interp_fraction*(env_struct%env_log10_pressure( &
                          env_point_index)-env_struct%env_log10_pressure(env_point_index-1))
                     star%logR(j) = env_struct%env_log10_radius(env_point_index-1)+ &
                          envelope_interp_fraction*(env_struct%env_log10_radius( &
                          env_point_index)-env_struct%env_log10_radius(env_point_index-1))
                     star%log_mass(j) = target_mass_coord
                     star%logT(j) = env_struct%env_log10_temperature( &
                          env_point_index-1)+envelope_interp_fraction*( &
                          env_struct%env_log10_temperature(env_point_index)- &
                          env_struct%env_log10_temperature(env_point_index-1))
                     star%xa(i_h1,j) = env_struct%env_hydrogen_fraction( &
                          env_point_index-1)+envelope_interp_fraction*( &
                          env_struct%env_hydrogen_fraction(env_point_index)- &
                          env_struct%env_hydrogen_fraction(env_point_index-1))
                     star%xa(i_metals,j) = env_struct%env_metal_fraction( &
                          env_point_index-1)+envelope_interp_fraction*( &
                          env_struct%env_metal_fraction(env_point_index)- &
                          env_struct%env_metal_fraction(env_point_index-1))
                     do k = 4,num_species
                        star%xa(k,j) = star%xa(k,star%nz)
                     end do
                     star%xa(i_he4,j)=1.0d0-star%xa(i_h1,j)- &
                          star%xa(i_metals,j)-star%xa(i_he3,j)
                     if (env_struct%env_convective_flag(env_point_index).or. &
                          env_struct%env_convective_flag(env_point_index-1)) then
                        star%convective_flag(j) = .true.
                     else
                        star%convective_flag(j) = .false.
                     endif
                  endif
                  star%nz = j
                  exit
               endif
            end do
! ASSIGN THE BOUNDARY AT THE PHOTOSPHERE FOR ENVELOPE MASS BELOW 1.0D-12.
! (On the exit path above num_zones was just set to j, so this guard is
! false there; on fall-through num_zones is unchanged and it is true.)
            if (j .gt. star%nz + env_struct%num_env_points) then
            star%nz = star%nz + env_struct%num_env_points
            end if
            if (star%job%rotation_active) then
               do j = old_last_shell+1,star%nz
                  star%omega(j) = star%omega(old_last_shell)
                  star%j_rot(j) = cc23*star%omega(old_last_shell)* &
                       exp10((2.0d0*star%logR(j)))
               end do
            endif
            write(*,910)
 910  format(1x,'NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS'/ &
            ' J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F')
      write(*,911)(j,star%logRho(j),star%luminosity_lsun(j),star%logP(j), &
                   star%logR(j),star%log_mass(j)-star%stotal, &
                   star%logT(j),star%convective_flag(j), j = old_last_shell,star%nz)
 911  format(i5,1p6e16.8,l2)
       endif
       envelope_recomputed_flag = .true.
       write(run_log_unit,597)old_senv,star%senv
 597     format(5x,'***** NEW ENVELOPE MASS CALCULATED *****'/8x, &
              'OLD SENV ',1pe22.13,'  NEW SENV',e22.13)
      exit envelope_rescale
      end do envelope_rescale
      endif

end subroutine rescale_and_refit_envelope

! ---------------------------------------------------------------
! Shell masses (unlogged) and per-shell dm from the (possibly
! refitted) log-mass grid.
subroutine build_shell_masses
      use observables_lib, only: shell_masses_from_log_mass
      call shell_masses_from_log_mass(star%log_mass, star%log_total_mass, &
           star%nz, star%m, star%dm)
end subroutine build_shell_masses

! ---------------------------------------------------------------
! Rotation state for the starting model: shape factors, moments of
! inertia, and the J / rotational-KE totals (with the historical
! previous-vs-new echo).
subroutine initialize_rotation_state
      if (star%job%rotation_active) then
! CALCULATE FP,FT,R0 AND ETA2 GIVEN OMEGA
       call rotation_shape_factors(star%logRho,star%logR,star%log_mass,star%nz,star%omega, &
            star%eta_squared,star%fp_rot,star%ft_rot, &
            star%mean_gravity,star%mean_radius,ierr)
       if (ierr /= 0) return
! FIND MOMENT OF INERTIA(HI)
       call zone_moments_of_inertia(star%eta_squared,star%logR,star%log_mass,star%dm,1,star%nz, &
            star%omega,star%mean_radius,star%i_rot,star%qiw)
! GIVEN OMEGA AND I, FIND ANGULAR MOMENTUM AND ROTATIONAL K.E.
       angular_momentum_sum = 0.0d0
       rotational_ke_sum = 0.0d0
       do i = 1,star%nz
          shell_angular_momentum = star%omega(i)*star%i_rot(i)
          star%j_rot(i) = shell_angular_momentum/star%dm(i)
          star%kinetic_energy_rot(i) = 0.5d0*star%omega(i)*shell_angular_momentum
          angular_momentum_sum = angular_momentum_sum+shell_angular_momentum
          rotational_ke_sum = rotational_ke_sum + star%kinetic_energy_rot(i)
       end do
       write(run_log_unit,560)total_angular_momentum, &
            angular_momentum_sum,total_rotational_ke,rotational_ke_sum
 560     format(1x,'TOTAL J OF STAR - PREVIOUS ',1pe21.13,' NEW ', &
              1pe21.13/1x,'TOTAL ROTATIONAL K.E. OF STAR - PREVIOUS ', &
              1pe21.13,' NEW ',1pe21.13)
       total_angular_momentum = angular_momentum_sum
       total_rotational_ke = rotational_ke_sum
      endif
end subroutine initialize_rotation_state

! ---------------------------------------------------------------
! Surface-mixture bookkeeping: on the first kind card, rebuild the
! envelope mass/number fractions (zenvm/amuenv/fxenv) and push them
! to the eos domain; regenerate the surface opacity tables (and the
! SCV envelope EOS table when in use) for the current surface X.
subroutine update_surface_mixture
      if (run_index.le.1) then
! SET UP MASS FRACTIONS AND NUMBER FRACTIONS OF ELEMENTS IN
! ENVELOPE.
! DBG 1/96 V (ENVELOPE MASS FRACTIONS WAS NORMALLY READ IN VIA
! RDLAOL. NOW THAT THESE OPACITIES ARE OBSOLETE VNEW IS INTRODUCED
! HOLD THE RELATIVE MASS FRACTIONS OF THE ELEMENTS (SEE PARMIN).
! TO MAINTAIN BACKWARD COMPATIBILITY IF LLAOL=T THEN USE V READ IN
! VIA RDLAOL OTHERWISE USE VNEW.
      if (.not.llaol) then
         do i=1, n_mix_species
            species_mix_weights(i)=star%ctrl%mixture_weights_seed(i)
         end do
      end if

! COMPUTE SURFACE MIX VALUES.
! ZENVM = Z OTHER THAN CNO CYCLE ELEMENTS;
! AMUENV = MEAN ATOMIC WEIGHT OF SURFACE Z;
! FXENV = NUMBER DENSITY OF SPECIES .
      star%envelope_hydrogen_fraction = star%xnew
      star%envelope_metal_fraction = star%znew
      star%envelope_helium_fraction = 1.0d0 - star%envelope_hydrogen_fraction - &
           star%envelope_metal_fraction - star%xa(i_he3,star%nz)
      star%envelope_he3_fraction = star%xa(i_he3,star%nz)
! SUM OF THE TEN METALS (EVERYTHING BUT H AND HE).
      mixture_weight_sum = species_mix_weights(ix_na)+species_mix_weights(ix_al)+ &
           species_mix_weights(ix_mg)+species_mix_weights(ix_fe)+ &
           species_mix_weights(ix_si)+species_mix_weights(ix_c)+ &
           species_mix_weights(ix_o)+species_mix_weights(ix_n)+ &
           species_mix_weights(ix_ar)+species_mix_weights(ix_ne)
      star%zenvm = star%envelope_metal_fraction*(mixture_weight_sum - &
           species_mix_weights(ix_c)-species_mix_weights(ix_o)- &
           species_mix_weights(ix_n))/mixture_weight_sum
      mixture_scale_factor = star%envelope_metal_fraction/mixture_weight_sum
      species_mix_weights(ix_h) = star%envelope_hydrogen_fraction/ &
           mixture_scale_factor
      species_mix_weights(ix_he) = (1.0d0-star%envelope_hydrogen_fraction- &
           star%envelope_metal_fraction)/mixture_scale_factor
      mixture_weight_sum = 0.0d0
      do i = 1,n_mix_species
       species_mix_weights(i) = mixture_scale_factor*species_mix_weights(i)/ &
            atomic_weight(i)
       mixture_weight_sum = mixture_weight_sum + species_mix_weights(i)
      end do
      star%amuenv = mixture_weight_sum
      mixture_scale_factor = 1.0d0/star%amuenv
! DBG 1/96 FXENV ARE NUMBER FRACTIONS OF ELEMENTS REQURIED
! BY EOS ROUTINES (SEE EQSTAT AND EQSAHA)
      do i = 1,n_mix_species
       star%fxenv(i) = species_mix_weights(i)*mixture_scale_factor
      end do
! push the recomputed mixture to the eos domain (physics-purity pass)
      call eos_set_mixture(star%envelope_hydrogen_fraction, &
           star%envelope_metal_fraction, star%amuenv, &
           star%fxenv)
!     FIND SURFACE COMPOSITION OPACITY TABLE
!     FIRST FIND INTERPOLATING FACTOR FOR COMPOSITION
      end if
! DBG 11/95 GENERATE NEW SURFACE OPACITY TABLES
      call kap_update_surface_tables(star%envelope_hydrogen_fraction, ierr=ierr)
      if (ierr /= 0) return
      if (use_scv_eos) then
         call build_scv_envelope_table
      endif

end subroutine update_surface_mixture

! ---------------------------------------------------------------
! Snapshot the step-start structure (the *_start arrays that the
! timestep limiters diff against) and, when rotating, old_omega.
subroutine snapshot_step_start
! CLONE P,T,R,L ARRAY TO DUMMY ARRAY HPOLD.
! HPOLD IS USED TO LIMIT THE TIMESTEP BASED ON CHANGES FROM
! MODEL TO MODEL IN P,T,R,L.
      do i = 1,star%nz
         star%logP_start(i) = star%logP(i)
         star%logT_start(i) = star%logT(i)
         star%logR_start(i) = star%logR(i)
         star%luminosity_lsun_start(i) = star%luminosity_lsun(i)
      end do
! Teff and the zone count are saved as well
      star%log_Teff_start = star%log_Teff
      star%nz_start = star%nz
      if (star%job%rotation_active) then
         do i = 1,star%nz
          star%old_omega(i) = star%omega(i)
         end do
      endif
end subroutine snapshot_step_start

end subroutine read_starting_model
