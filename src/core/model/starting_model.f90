!----------------------------------------------------------------------
! starting_model
!----------------------------------------------------------------------
! The starting model of each run card (read_starting_model): read from
! file or reused from memory, rescaled and refitted as the card asks,
! and set up with the state the first step needs.
module starting_model
      use turnover_timescale, only: compute_turnover_timescale
      use shell_physics, only: compute_shell_physics
      use model_rescaling, only: rescale_model
      implicit none
      private
      public :: read_starting_model

contains

!----------------------------------------------------------------------
! read_starting_model (legacy starin)
!----------------------------------------------------------------------
! Sets up the starting model of run card run_index, at the start of
! every card, in named steps:
!   acquire_starting_model      read the model file (.mod, YREC7 or
!                               MODEL2 layout) when first_call_flag is
!                               set, else reuse the model in memory;
!   extend_core_toward_center   optionally add points inward of the
!                               innermost shell (file models only);
!   rescale_and_refit_envelope  rescale X, Z or the mass on a rescale
!                               card, and optionally move the envelope
!                               fitting point;
!   build_shell_masses, initialize_rotation_state,
!   update_surface_mixture, snapshot_step_start
!                               shell masses, rotation totals, surface
!                               mixture and opacity tables, step-start
!                               copies.
! Then shell_physics, am_convective_regions and
! compute_turnover_timescale give the gradients, opacities, convection
! zones and turnover time of the starting model.
! Called with model_failed_flag set (a step failed to converge), it
! only stops the run with ierr = 1.
!
! Arguments: timestep_yr, delta_time (s), delta_time_abs = |delta_time|
! and triangle_orientation come from the model file; envelope_recomputed_flag
! is set when the envelope must be recomputed; total_angular_momentum
! and total_rotational_ke return the rotating model's totals;
! species_mix_weights returns the surface mixture as X_i/A_i (first
! card only).
!
! Naming: the star% arrays use henyey_iterate's names (luminosity_lsun
! is the linear L/Lsun; m is the unlogged mass coordinate and dm the
! mass of each shell, in g; mean_radius, qiw and mean_gravity are the
! rotation geometry); the
! single-point physics in set_new_last_shell_convection uses
! shell_physics' log10_* names.
subroutine read_starting_model(timestep_yr, delta_time, delta_time_abs, &
     triangle_orientation, ikut_flag, stored_vertex_index, model_failed_flag, &
     envelope_recomputed_flag, run_index, dlnrho_dlnp, dlnrho_dlnt, &
     total_angular_momentum, total_rotational_ke, convective_velocity, &
     species_mix_weights, ierr)
      use am_regions, only: am_convective_regions
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use convection_lib, only: compute_temperature_gradients, set_convection_envelope_metals
      use star_info_lib, only: star, i_be9, i_c12, i_c13, i_h1, i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18, json, &
           n_species_basic, n_species_extended, n_mix_species, ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, ix_h, ix_o, ix_n, ix_ar, ix_ne, ix_he, &
           max_convective_zones, max_radiative_zones
      use controls_buffer, only: ichi_dm_max, ichi_dp_env_max
      use envelope_integration, only: atm_get, envint_step_config, fixed_envint_step, &
           ambient_log10_teff
      use envelope_refit, only: append_envelope_points
      use envelope_structure
      use io_units
      use model_files, only: read_mod_model, read_yrec7, read_model2
      use physical_constants
      use eos_lib
      use kap_lib
      use stitched_model, only: build_stitched_model
      implicit none
      integer :: jerr_atm

      double precision, intent(inout) :: timestep_yr
      double precision, intent(inout) :: delta_time, delta_time_abs
      double precision, intent(inout) :: triangle_orientation
      integer, intent(out) :: ikut_flag
      integer, intent(inout) :: stored_vertex_index
      logical, intent(in) :: model_failed_flag
      logical, intent(inout) :: envelope_recomputed_flag
      integer, intent(in) :: run_index
      double precision, intent(out) :: dlnrho_dlnp, dlnrho_dlnt
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke
      double precision, intent(out) :: convective_velocity
      double precision, intent(inout) :: species_mix_weights(n_mix_species)

      integer :: saha_state
      character(len=4) :: format_tag
      character(len=6) :: eos_code
      character(len=4) :: atm_code, alok_code, hik_code
! The initial-composition code read with the model is kept in
! star%initial_composition_code.
! Atomic weights of the surface-mixture elements, in ix_* slot order
! (Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He).
      double precision :: atomic_weight(n_mix_species)
      data atomic_weight /23.0d0,26.99d0,24.32d0,55.86d0,28.1d0,12.015d0, &
           1.008d0,16.0d0,14.01d0,39.96d0,20.19d0,4.004d0/

! Central composition of the starting model: the reference for the
! homogeneity test that must pass before the mixture is altered.
      double precision :: reference_composition(n_species_extended)

! --- locals ---
      integer :: iread
      double precision :: mixing_length0
      logical :: mixing_length_matches
      logical :: mixture_ok
! lexcp0 is never assigned: the mismatch warning in
! check_model_against_controls prints it where it presumably means
! use_extended_composition0 (a bug kept from starin).
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
! EOS and opacity results for the convective test of a new last shell.
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results)
      double precision :: radiative_gradient, dgrad_dt_component, &
           dgrad_dp_component, dgrad_dr_component
      logical :: is_convective
      integer :: old_last_shell
      type(envint_step_config) :: env_steps
      logical :: save_boundary_flag, print_flag
      double precision :: log10_gravity
      double precision :: log10_pressure_limit
      double precision :: spot_adjusted_log_teff
      double precision :: pressure_offset, density_offset, &
           temperature_offset, radius_offset
      double precision :: old_senv, target_log_mass_at_fit
      double precision :: interior_interp_fraction
      integer :: num_species
      double precision :: angular_momentum_sum, rotational_ke_sum, &
           shell_angular_momentum
      double precision :: mixture_weight_sum, mixture_scale_factor
      integer :: radiative_zone_bounds(max_radiative_zones,2), &
           convective_zone_bounds(max_convective_zones,2)
      integer :: num_radiative_zones, num_convective_zones
      integer :: core_cz_top_index0, envelope_cz_bottom_index0
      logical :: rotation_active0
      logical :: use_diffusion_y0, use_diffusion_z0, disk_locking_active0, &
           instability_transport_active0, ljdot00
      logical :: lovstc0, envelope_overshoot_active0, lovstm0, &
           use_pure_z_table0, lsemic0
      double precision :: disk_pressure0, disk_temperature0, &
           wind_saturation_omega0
! Error code from the library calls, passed on through ierr.
      integer :: jerr

      integer, intent(out) :: ierr

      ierr = 0

! model_failed_flag (legacy LARGE) means the last model failed to
! converge: stop the run.
      if (model_failed_flag) then
          write(run_log_unit,1000)
 1000       format(1x,39('>'),40('<')/, &
                  "STARIN:        ***** RUN STOPPED *****")
          write(run_log_unit,1010)
 1010       format("STARIN: ***** MODEL FAILED TO CONVERGE *****")
            ! run_yrec returns the error; the command-line wrapper
            ! (main) stops.
            ierr = 1
            return
      endif

! The model file is on unit star%ctrl%first_unit (legacy IFIRST).
      ikut_flag = 0
      iread = star%ctrl%first_unit

! Initial values: no convective velocity, and the ideal-gas
! dlnrho/dlnT = -1 and dlnrho/dlnP = 1.
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
      call compute_shell_physics(star%fp_rot,star%ft_rot, &
           star%xa,star%logRho,star%luminosity_lsun,star%logP, &
           star%logR,star%log_mass,star%logT,star%convective_flag,star%nz, &
           star%log_Teff, jerr)
      if (jerr /= 0) then
      ! Pass the error on.
         ierr = jerr
         return
      end if
      call am_convective_regions(star%xa,star%logRho,star%logP,star%logR,star%log_mass, &
                 star%logT,star%convective_flag,star%nz, &
                 star%am_transport_convective_flag,radiative_zone_bounds, &
                 convective_zone_bounds,num_radiative_zones, &
                 num_convective_zones, star%iov1, star%iov2, ierr)
      if (ierr /= 0) return
! Initialize the convective turnover time (legacy TAUCZ), the
! photospheric pressure pphot and fracstep. The stitched model
! (interior plus a fresh envelope) is built first: it sets star%pphot,
! and compute_turnover_timescale walks it, returning the radius of the
! convection-zone base in star%envelope_radius. (starin passed the
! envelope radius array there, overwriting its first element.) The
! wind takes pphot and the turnover time at the fraction fracstep of
! the step, between the step-start copies (pphot0, *_old) and the
! current values.
      call build_stitched_model
      call compute_turnover_timescale(star%envelope_radius, ierr)
      if (ierr /= 0) return
      star%convective_turnover_timescale_old = star%convective_turnover_timescale
      star%pphot0 = star%pphot
      star%fracstep = 0.5

      return

contains

! ---------------------------------------------------------------
! Obtain the starting model. With first_call_flag(run_index) false,
! reuse the model in memory (legacy LFIRST). With it true, read the
! model file on unit iread, warn if its mixing length or composition
! set differs from the controls, reset the step state, set the
! rotation curve and, in a chemically homogeneous model, alter the
! heavy-element mixture and isotope ratios if asked. Sets ierr on
! failure.
subroutine acquire_starting_model

! A reused model (first_call_flag false) stays as it is in memory. On a
! rescale card it is rescaled once, in rescale_and_refit_envelope, as a
! model read from file is. (starin rescaled a reused model here as well,
! so it was rescaled twice.)
      if (star%job%first_call_flag(run_index)) then
         call read_starting_model_file
         if (ierr /= 0) return
         call check_model_against_controls
         call reset_state_after_read
         call set_starting_rotation_curve
! Leave unit iread open, rewound.
         rewind iread
! Alter the heavy-element mixture (MHP 10/24), only in a chemically
! homogeneous model read from file. change_cno_mixture_active (legacy
! ISETMIX) sets the C, N and O fractions of Z; change_isotope_ratios_active
! (legacy ISETISO) sets the isotope ratios and light-element
! abundances. Between them they can overwrite mass fractions 4-15.
! (The test below is always true inside this branch; starin once also
! required the first card.)
         if (star%job%first_call_flag(run_index)) then
            call check_mixture_homogeneity
            if (mixture_ok) then
               call rescale_cno_mixture
               call set_isotope_ratios
            end if
         end if
      endif
end subroutine acquire_starting_model

! ---------------------------------------------------------------
! Read the model file on unit iread, in the layout its first keyword
! names (YMOD: YREC .mod; NMOD: YREC7; MOD2: MODEL2), and set the
! timestep and star%stotal from it. Sets ierr on an unknown layout or
! a failed read.
subroutine read_starting_model_file

      rewind iread
      read(iread,10) format_tag
 10   format(a4)

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
! read_yrec7 fills the star% model directly; its arguments are the
! file's own header values (mixing length, timestep, flags, physics
! codes), which check_model_against_controls compares with the
! controls.
         call read_yrec7(mixing_length0,timestep_yr,triangle_orientation,iread, &
              core_cz_top_index0,envelope_cz_bottom_index0, &
              use_extended_composition0,rotation_active0, &
              atm_code,eos_code,hik_code,use_diffusion_y0,use_diffusion_z0, &
              disk_locking_active0,instability_transport_active0,ljdot00, &
              alok_code,lovstc0,envelope_overshoot_active0,lovstm0, &
              use_pure_z_table0,lsemic0,disk_pressure0, &
              disk_temperature0,wind_saturation_omega0, ierr)
         if (ierr /= 0) return
! The header values up to rotation_active0 are YREC7's; those from
! atm_code on were added by the MODEL2 layout.

      else if (format_tag .eq. 'MOD2 ') then
         write(run_log_unit,16)
 16      format(1x,'loading starting model (MODEL2 format)')
! Same header values as read_yrec7 (no ierr).
         call read_model2(mixing_length0,timestep_yr,triangle_orientation,iread, &
              core_cz_top_index0,envelope_cz_bottom_index0, &
              use_extended_composition0,rotation_active0, &
              atm_code,eos_code,hik_code,use_diffusion_y0,use_diffusion_z0, &
              disk_locking_active0,instability_transport_active0,ljdot00, &
              alok_code,lovstc0,envelope_overshoot_active0,lovstm0, &
              use_pure_z_table0,lsemic0,disk_pressure0, &
              disk_temperature0,wind_saturation_omega0)

      else
         write(run_log_unit,20)
 20      format('STARIN: ***** RUN TERMINATED, INVALID INPUT', &
                ' MODEL FILE.  *****')
         ! run_yrec returns the error; the command-line wrapper (main)
         ! stops.
         ierr = 1
         return
      endif

! The timestep in s, and the log total mass.
      delta_time = seconds_per_year*timestep_yr
      delta_time_abs = dabs(delta_time)
      star%stotal = star%log_total_mass
end subroutine read_starting_model_file

! ---------------------------------------------------------------
! Warn when an evolved model (model_number > 0) has a mixing length or
! extended-composition flag that differs from the controls. A warning
! only: the run goes on.
subroutine check_model_against_controls
      if (star%model_number.gt.0) then
       mixing_length_matches = .true.
! The tolerance 2e-3 allows for model files storing the mixing length
! to four significant digits; a model mixing length <= 0 is not
! checked.
       if (mixing_length0.gt.0.0d0) mixing_length_matches = &
            (dabs(star%mixing_length_alpha-mixing_length0).lt.2.0d-3)
! The warning prints lexcp0, which is never set (see its declaration),
! for the model's flag.
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
       endif
      endif
end subroutine check_model_against_controls

! ---------------------------------------------------------------
! Reset the store and trial-sign flags, require a fresh envelope, and
! take the surface X and Z from the outermost interior point.
subroutine reset_state_after_read
      stored_vertex_index = 0
      triangle_orientation = 1.0d0
! Recompute the envelope: the envelope triangle stored with the model
! may not be reliable.
       envelope_recomputed_flag = .true.


! Surface X and Z (xnew, znew) from the outermost interior point.
      star%xnew = star%xa(i_h1,star%nz)
      star%znew = star%xa(i_metals,star%nz)
end subroutine reset_state_after_read

! ---------------------------------------------------------------
! Set the solid-body rotation curve when asked (rotating runs), or
! unit rotation factors fp and ft (non-rotating runs).
subroutine set_starting_rotation_curve

      if (star%job%rotation_active) then
       if (star%ctrl%set_initial_omega) then
! Solid-body rotation, omega = initial_omega (legacy WNEW) in every
! shell: this turns a non-rotating model into a rotating one.
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
end subroutine set_starting_rotation_curve

! ---------------------------------------------------------------
! Check that the starting model is chemically homogeneous (every shell
! within 1e-6 of the center in each of species 1-15) before its mixture
! is altered; mixture_ok is false, with a message, if not.
subroutine check_mixture_homogeneity
      mixture_ok = .true.
      if (star%ctrl%change_cno_mixture_active .or. star%ctrl%change_isotope_ratios_active) then
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
end subroutine check_mixture_homogeneity

! ---------------------------------------------------------------
! Scale the C, N and O isotopes in every shell so that C, N and O make
! up the requested fractions of Z. Z itself is not changed.
subroutine rescale_cno_mixture
      if (star%ctrl%change_cno_mixture_active) then
! The current fractions, from the (homogeneous) central composition:
!   f_C = (X(C12) + X(C13))/Z,  f_N = (X(N14) + X(N15))/Z,
!   f_O = (X(O16) + X(O17) + X(O18))/Z.
! Every isotope of an element is scaled by target/f, which keeps the
! isotope ratios; set_isotope_ratios may then change them.
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
end subroutine rescale_cno_mixture

! ---------------------------------------------------------------
! Set the requested He3, H2, Li6, Li7 and Be9 abundances and the
! C12/C13 and O16/O18 ratios in every shell. The central sums
! S_C = X(C12) + X(C13) and S_O = X(O16) + X(O18) are kept and split
! by the requested ratio r:
!   X(C12) = r_C S_C/(1 + r_C),  X(C13) = S_C/(1 + r_C),
! and likewise for O16 and O18. N15 and O17 are not altered, and B10
! and B11 are not species of the model.
subroutine set_isotope_ratios
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
end subroutine set_isotope_ratios

! ---------------------------------------------------------------
! With extend_core_inward (legacy LCORE), add N points inward of the
! innermost shell, in whole steps of the maximum log-mass spacing
! dlogm = chi_grid_scale(ichi_dm_max) (legacy HPTTOL(2)) so that the
! grid needs no reshuffling (MHP 9/98). The new points follow what the
! model assumes inside its central fitting point: constant density,
! constant energy generation rate and constant d log T/d log P.
! It is done on every card that reads its model from file, not only the
! first, so that auto-calibrated solar models keep the change (MHP 9/14).
subroutine extend_core_toward_center
      use math_lib

      if (star%ctrl%extend_core_inward .and. star%job%first_call_flag(run_index)) then
! The innermost mass is to drop by the factor core_mass_reduction_factor
! (legacy FCORE): that takes N = int(log10(factor)/dlogm) + 1 steps
! (num_core_shells_added, legacy MCORE). The control is then
! overwritten with the log10 shift N dlogm actually made.
          star%job%num_core_shells_added = int(log10(star%job%core_mass_reduction_factor)/ &
               star%ctrl%chi_grid_scale(ichi_dm_max))+1
          star%job%core_mass_reduction_factor = dble(star%job%num_core_shells_added)* &
               star%ctrl%chi_grid_scale(ichi_dm_max)
          num_shells_extended = star%nz + star%job%num_core_shells_added
! Too many points for the json-point arrays: the run stops with ierr
! (starin printed RUN TERMINATED but went on writing past json; fixed
! 2026-10).
          if (num_shells_extended .gt. json) then
             write(run_log_unit,476)"STARIN: Unable to extend core inward ", &
                   "- JSON too small"
             write(run_log_unit,477) "STARIN: Required size =", &
                    num_shells_extended, ", JSON = ", json
             write(run_log_unit,478) "STARIN:  ***** RUN TERMINATED *****"
  476        format(2a)
  477        format(a, i8, a, i8)
  478        format(a)
             ierr = 1
             return
          endif
          core_shell_spacing = star%ctrl%chi_grid_scale(ichi_dm_max)
! Move every old point up by N slots (logRho is not moved).
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
! March inward in mass from the innermost old point (now slot N+1).
! Slots 1 and 2 still hold the old innermost values here, so rho_c is
! the old central density and del = d log T/d log P is taken between
! the two innermost old points (for N = 1, slot 2 has already been
! overwritten and del is 0/0).
          central_log_density = star%logRho(1)
          actual_gradient = (star%logT(2)-star%logT(1))/ &
               (star%logP(2)-star%logP(1))
          central_shell_luminosity = star%luminosity_lsun(1)
! Offset for estimating rho from P and T at fixed mean molecular weight
! (MHP 4/12): log rho = log P - log T - offset, with the offset taken
! at slot N+1.
          density_estimate_offset = star%logP(first_original_shell)- &
               star%logRho(first_original_shell)- &
               star%logT(first_original_shell)
          do i = star%job%num_core_shells_added,1,-1
             star%log_mass(i) = star%log_mass(i+1)-core_shell_spacing
! Constant density, m = (4 pi/3) rho_c r^3:
!   log r = (1/3)(log m - log(4 pi/3) - log rho_c).
             star%logR(i) = cc13*(star%log_mass(i)-c4pi3l-central_log_density)
! Constant epsilon: L = L_1 m/m_1, L_1 and m_1 at the innermost old
! point.
             star%luminosity_lsun(i) = exp(ln10*(star%log_mass(i)- &
                  star%log_mass(first_original_shell)))*central_shell_luminosity
! Hydrostatic equilibrium, d log P/d log m = -G m^2/(4 pi r^4 P), taken
! inward in two stages: log P half a step in from point i+1, then the
! full step with that pressure:
!   log P_i = log P_(i+1) + (G m_i^2/(4 pi r_i^4 P_half)) dlogm.
             trial_log_pressure = star%logP(i+1)
             temp_scratch = exp(ln10*(cgl+2.0d0*star%log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*star%logR(i)))
             trial_log_pressure = trial_log_pressure+0.5d0*temp_scratch* &
                  core_shell_spacing
             temp_scratch = exp(ln10*(cgl+2.0d0*star%log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*star%logR(i)))
             star%logP(i) = star%logP(i+1)+temp_scratch* &
                  core_shell_spacing
! Convective or radiative as the innermost old point.
             star%convective_flag(i) = star%convective_flag(first_original_shell)
! Constant del: log T_i = log T_(i+1) + del (log P_i - log P_(i+1)).
             star%logT(i) = star%logT(i+1)+temp_scratch* &
                  actual_gradient*core_shell_spacing
! Omega and composition as at the innermost old point.
             star%omega(i) = star%omega(first_original_shell)
             do j=1, 15
                star%xa(j,i) = star%xa(j,first_original_shell)
             end do
! The local density estimate, in place of an EOS call that did not
! work here (MHP 4/12).
             star%logRho(i) = star%logP(i) - star%logT(i) - &
                  density_estimate_offset
          end do
          star%nz = num_shells_extended
      end if
end subroutine extend_core_toward_center

! ---------------------------------------------------------------
! On a rescale card (rescale_kind /= 1), rescale the model's X, Z or
! mass (rescale_model). Then, with change_envelope_mass_flag, move the
! envelope fitting point to requested_envelope_mass = log10(m_fit/M)
! (legacy SENV0), on the first model of every card (DBG 2/92):
! - deeper than now: drop the points outside the new fitting mass and
!   put the new last point at it, linear in log m between the old
!   points on either side, with the composition of the one inside;
! - shallower: integrate an envelope with the model's Teff and L from
!   the surface down to the old last point, and append its points
!   inside the new fitting mass as shells, with the old last point's
!   composition apart from X and Z.
! Sets ierr on failure.
subroutine rescale_and_refit_envelope
      use math_lib

      if (star%job%rescale_kind(run_index).ne.1) call rescale_model(star%luminosity_lsun, &
           star%xa,star%log_mass,star%log_total_mass,star%nz,run_index, &
           star%star_mass,star%convective_flag, ierr)
      if (ierr /= 0) return

! A once-through loop: exit leaves the envelope change early.
      if (star%job%change_envelope_mass_flag) then
      envelope_rescale: do
       call clamp_requested_envelope_mass
       if (star%senv.eq.star%job%requested_envelope_mass) exit envelope_rescale
       num_species = n_species_basic
       if (star%job%use_extended_composition) num_species = n_species_extended
       if (star%job%requested_envelope_mass.lt.star%senv) then
! The new envelope is deeper than the old one.
          call find_shell_below_new_fit
          if (i < (1)) then
! The requested fitting mass lies inside the innermost shell: report
! it and leave the envelope mass unchanged.
          write(run_log_unit,576)star%job%requested_envelope_mass
 576        format(5x,'ERROR IN SUBROUTINE STARIN'/5x,'DESIRED', &
              ' ENVELOPE MASS',1pe22.13,' TOO LARGE'/5x,'ENVELOPE', &
              ' MASS NOT CHANGED')
          exit envelope_rescale
          end if
          call interpolate_new_fitting_point
          call set_new_last_shell_convection
          if (ierr /= 0) return
       else
! The new envelope is shallower than the old one.
          call integrate_envelope_to_new_fit
          if (ierr /= 0) return
          call append_envelope_as_shells
          if (ierr /= 0) return
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
! Make the requested envelope mass, log10(m_fit/M), negative and no
! shallower than -1e-12, and record the current one,
! senv = log m(nz) - log M, in senv and old_senv.
subroutine clamp_requested_envelope_mass
       if (star%job%requested_envelope_mass.gt.0.0d0) star%job%requested_envelope_mass = &
            -star%job%requested_envelope_mass
! The limit -1e-12 stands for the base of the atmosphere, tau = 2/3
! (DBG 2/92; it was -1e-10).
       if (star%job%requested_envelope_mass.gt.-1.0d-12) star%job%requested_envelope_mass = &
            -1.d-12
       star%senv = star%log_mass(star%nz) - star%log_total_mass
       old_senv = star%senv
end subroutine clamp_requested_envelope_mass

! ---------------------------------------------------------------
! Find the outermost shell i inside the new (deeper) fitting mass,
! log m_fit = log M + requested_envelope_mass; i = 0 when the fitting
! mass lies inside the innermost shell.
subroutine find_shell_below_new_fit
          target_log_mass_at_fit = star%log_total_mass+star%job%requested_envelope_mass
          do i = star%nz-1,1,-1
             if (star%log_mass(i).lt.target_log_mass_at_fit) exit
          end do
end subroutine find_shell_below_new_fit

! ---------------------------------------------------------------
! Make shell i+1 the new last interior shell, at the fitting mass,
! with the composition of shell i and the structure (and omega) linear
! in log m between shells i and i+1:
!   f = (log m_fit - log m_i)/(log m_(i+1) - log m_i),
!   y = y_i + f (y_(i+1) - y_i) for log rho, L, log P, log R, log T.
! The new surface X and Z (xnew, znew) are its X and Z.
subroutine interpolate_new_fitting_point
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
end subroutine interpolate_new_fitting_point

! ---------------------------------------------------------------
! Decide whether the new last shell is convective: as its neighbours
! when they agree; otherwise from the temperature gradients at the
! point (eos_get with the atmosphere EOS, kap_get,
! compute_temperature_gradients, no rotation), which also give its
! density.
subroutine set_new_last_shell_convection
          if (star%convective_flag(i).and.star%convective_flag(i+1)) then
             star%convective_flag(star%nz) = .true.
          else if (.not.star%convective_flag(i).and. .not.star%convective_flag(i+1)) &
               then
             star%convective_flag(star%nz) = .false.
          else
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
               eos_res%log10_density = log10_density
               call eos_get(log10_temperature, log10_pressure, &
                    hydrogen_fraction, metal_fraction, eos_res, &
                    want_derivatives, in_atmosphere, saha_state, &
                    composition_at_zone=star%xa(:,star%nz), ierr=ierr)
               if (ierr /= 0) return
! The opacity is taken at the density eos_get returns.
               call kap_get(eos_res%log10_density, log10_temperature, &
                    hydrogen_fraction, metal_fraction, kap_res, &
                    eos_res%ion_fraction, ierr=ierr)
               if (ierr /= 0) return
               call compute_temperature_gradients(log10_temperature, log10_pressure, &
                    eos_res, kap_res, log10_radius, log10_mass, &
                    shell_luminosity_lsun, actual_gradient, radiative_gradient, &
                    dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
                    convective_velocity, want_derivatives, is_convective, &
                    point_pressure_rotation_factor, &
                    point_temperature_rotation_factor, star%log_Teff, .false., jerr)
               if (jerr /= 0) then
               ! Pass the error on.
                  ierr = jerr
                  return
               end if
               log10_density = eos_res%log10_density
               star%logRho(star%nz) = log10_density
               star%convective_flag(star%nz) = is_convective
          endif
end subroutine set_new_last_shell_convection

! ---------------------------------------------------------------
! Integrate an envelope with the model's Teff and L from the surface
! down to the pressure of the old last interior shell (atm_get), for
! the new, shallower fitting mass; sets senv to it. Sets ierr on
! failure.
subroutine integrate_envelope_to_new_fit
      use math_lib
            old_last_shell = star%nz
! atm_get integrates with one fixed step for start, minimum and
! maximum: chi_grid_scale(ichi_dp_env_max) (legacy HPTTOL(8)).
            env_steps = fixed_envint_step(star%ctrl%chi_grid_scale(ichi_dp_env_max))
          save_boundary_flag = .false.
          print_flag = .true.
          saha_state = 0
! Surface radius and gravity from L = 4 pi R^2 sigma Teff^4 and
! g = G M/R^2:
!   log R = (log L/Lsun + log Lsun - 4 log Teff - log(4 pi) - log sigma)/2,
! with log Lsun = star%log10_solar_luminosity (starin's CLSUNL; the
! 2026-08 modernization had put Lsun in erg/s there, fixed 2026-10-12).
          shell_luminosity_lsun = exp(ln10*star%log_L)
          log10_radius = 0.5d0*(star%log_L + star%log10_solar_luminosity - &
               4.0d0*star%log_Teff - c4pil - csigl)
          log10_gravity = cgl + star%stotal - log10_radius - log10_radius
          hydrogen_fraction = star%xa(i_h1,star%nz)
          metal_fraction = star%xa(i_metals,star%nz)
          point_pressure_rotation_factor = 1.0d0
          point_temperature_rotation_factor = 1.0d0
          log10_pressure_limit = star%logP(star%nz)
! The Debye-Hueckel correction takes the last shell's composition.
            call eos_set_debye_huckel_composition(star%xa(:,star%nz))
! Spotted stars (G. Somers 10/14): with spot_filling_factor f /= 0,
! spot_temp_contrast x /= 1 and envelope_cz_bottom_index = nz, the
! envelope is integrated at the Teff of the unspotted surface
! (ambient_log10_teff, legacy ATEFFL):
!   Teff_amb^4 = Teff^4/(f x^4 + 1 - f).
          if (star%envelope_cz_bottom_index.eq.star%nz.and.star%ctrl%spot_filling_factor.ne. &
               0.0.and.star%ctrl%spot_temp_contrast.ne.1.0) then
               spot_adjusted_log_teff = ambient_log10_teff(star%log_Teff)
          else
             spot_adjusted_log_teff = star%log_Teff
          endif
          call atm_get(shell_luminosity_lsun,point_pressure_rotation_factor, &
                 point_temperature_rotation_factor,log10_gravity,star%stotal, &
                 print_flag,save_boundary_flag, &
                 log10_pressure_limit,log10_radius,spot_adjusted_log_teff, &
                 hydrogen_fraction,metal_fraction,saha_state, &
                 ierr=jerr_atm,env_steps=env_steps)
! atm_get returns integration failures (numerics_termination among
! them) in jerr_atm; pass them on.
          if (jerr_atm /= 0) then
             ierr = jerr_atm
             return
          end if
          star%senv = star%job%requested_envelope_mass
end subroutine integrate_envelope_to_new_fit

! ---------------------------------------------------------------
! Make the fresh envelope continuous with the interior and append its
! points inside the new fitting mass as shells; in rotating models they
! take the old last shell's omega. Echoes the new points. Sets ierr on
! failure.
subroutine append_envelope_as_shells
      use math_lib
! Stop if the envelope points would overflow the json-point arrays.
            if (star%nz+env_struct%num_env_points.ge.json) then
               write(*,*) 'read_starting_model: grid overflow (json too small)'
               ierr = 1
               return
            end if
! Envelope point 1 is the deepest, at the old last shell's pressure.
! Shift log P, log rho, log T and log R of the envelope by the
! interior-minus-envelope difference there, so that the envelope joins
! the interior continuously, and drop point 1. The convective flag moves
! with the rest (the F77 starin.f left it unshifted, one point off;
! fixed 2026-10).
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
               env_struct%env_convective_flag(j) = env_struct%env_convective_flag(j+1)
            end do
            env_struct%num_env_points = env_struct%num_env_points - 1
! Append the envelope points inside the new fitting mass as shells
! and interpolate the last one onto it (append_envelope_points, shared
! with rebuild_envelope; here a point exactly at the fitting mass is
! kept as a shell).
            call append_envelope_points(.true., num_species, star%xa, &
                 star%logRho, star%luminosity_lsun, star%logP, star%logR, &
                 star%log_mass, star%logT, star%convective_flag, star%nz, ierr)
            if (ierr /= 0) then
! A degenerate interpolation interval in mass is fatal; run_yrec
! returns the error.
               write(run_log_unit,*) 'STARIN: degenerate', &
                    ' envelope interpolation interval (was STOP 9998)'
               return
            end if
! The new shells rotate with the old last shell's omega, with the
! specific angular momentum of a thin shell, j = (2/3) omega r^2.
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
end subroutine append_envelope_as_shells

! ---------------------------------------------------------------
! The unlogged mass coordinate m and the mass dm of each shell (g),
! from the (possibly extended or refitted) log-mass grid.
subroutine build_shell_masses
      use observables, only: shell_masses_from_log_mass
      call shell_masses_from_log_mass(star%log_mass, star%log_total_mass, &
           star%nz, star%m, star%dm)
end subroutine build_shell_masses

! ---------------------------------------------------------------
! Rotation state of the starting model: shape factors, moments of
! inertia, and the totals of angular momentum and rotational kinetic
! energy, echoed beside the totals held on entry.
subroutine initialize_rotation_state
      if (star%job%rotation_active) then
! Shape factors fp_rot and ft_rot, mean gravity and radius, and
! eta_squared, from omega.
       call rotation_shape_factors(star%logRho,star%logR,star%log_mass,star%nz,star%omega, &
            star%eta_squared,star%fp_rot,star%ft_rot, &
            star%mean_gravity,star%mean_radius,ierr)
       if (ierr /= 0) return
! Moment of inertia of each shell, i_rot (legacy HI).
       call zone_moments_of_inertia(star%eta_squared,star%logR,star%log_mass,star%dm,1,star%nz, &
            star%omega,star%mean_radius,star%i_rot,star%qiw)
! Per shell: angular momentum J_i = omega_i I_i, specific angular
! momentum j_rot = J_i/dm_i and kinetic energy (1/2) omega_i J_i;
! summed over the star.
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
! Surface-mixture bookkeeping. On the first card (run_index <= 1), set
! the surface X, Y, Z and He3, build from them and the controls' metal
! mixture the number fractions fxenv of the 12 mixture elements, Z
! without C, N and O (zenvm) and amuenv, and pass them to the EOS. On
! every card, regenerate the surface opacity tables, and the SCV
! envelope EOS table when SCV is on, for the current surface
! composition.
subroutine update_surface_mixture
      if (run_index.le.1) then
! The relative mass fractions w_i of the ten metals are the controls'
! mixture_weights_seed.
      do i=1, n_mix_species
         species_mix_weights(i)=star%ctrl%mixture_weights_seed(i)
      end do

! Surface composition: X_env = xnew, Z_env = znew, X(He3) of the last
! shell, Y_env = 1 - X_env - Z_env - X(He3).
      star%envelope_hydrogen_fraction = star%xnew
      star%envelope_metal_fraction = star%znew
      star%envelope_helium_fraction = 1.0d0 - star%envelope_hydrogen_fraction - &
           star%envelope_metal_fraction - star%xa(i_he3,star%nz)
      star%envelope_he3_fraction = star%xa(i_he3,star%nz)
! S = sum of w_i over the ten metals (all but H and He).
      mixture_weight_sum = species_mix_weights(ix_na)+species_mix_weights(ix_al)+ &
           species_mix_weights(ix_mg)+species_mix_weights(ix_fe)+ &
           species_mix_weights(ix_si)+species_mix_weights(ix_c)+ &
           species_mix_weights(ix_o)+species_mix_weights(ix_n)+ &
           species_mix_weights(ix_ar)+species_mix_weights(ix_ne)
! The part of Z not in C, N and O: zenvm = Z_env (S - w_C - w_O - w_N)/S.
      star%zenvm = star%envelope_metal_fraction*(mixture_weight_sum - &
           species_mix_weights(ix_c)-species_mix_weights(ix_o)- &
           species_mix_weights(ix_n))/mixture_weight_sum
      call set_convection_envelope_metals(star%zenvm)
! With f = Z_env/S, the metals' mass fractions are X_i = f w_i; H and He
! get w = X_env/f and (1 - X_env - Z_env)/f, so that f w is their mass
! fraction too. Then w_i becomes X_i/A_i, and
!   amuenv = sum_i X_i/A_i,
! the inverse mean atomic weight 1/mu_ion of the whole surface mixture.
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
! Number fractions of the elements, fxenv(i) = (X_i/A_i)/amuenv, which
! the EOS needs (DBG 1/96): compute_blended_eos and compute_saha_eos
! (legacy EQSAHA).
      do i = 1,n_mix_species
       star%fxenv(i) = species_mix_weights(i)*mixture_scale_factor
      end do
! Pass the new mixture to the EOS.
      call eos_set_mixture(star%envelope_hydrogen_fraction, &
           star%envelope_metal_fraction, star%amuenv, &
           star%fxenv)
      end if
! Regenerate the surface opacity tables for the surface X (DBG 11/95),
! and the SCV envelope EOS table when SCV is on.
      call kap_update_surface_tables(star%envelope_hydrogen_fraction, ierr=ierr)
      if (ierr /= 0) return
      call eos_update_envelope_tables(star%envelope_hydrogen_fraction, &
           star%envelope_metal_fraction)

end subroutine update_surface_mixture

! ---------------------------------------------------------------
! Copy the step-start structure: log P, log T, log R and L into the
! *_start arrays (legacy HPOLD), against which the timestep control
! limits the change from model to model; with them Teff and the surface
! L, and, when rotating, omega into old_omega.
subroutine snapshot_step_start
      do i = 1,star%nz
         star%logP_start(i) = star%logP(i)
         star%logT_start(i) = star%logT(i)
         star%logR_start(i) = star%logR(i)
         star%luminosity_lsun_start(i) = star%luminosity_lsun(i)
      end do
      star%log_Teff_start = star%log_Teff
      star%surface_luminosity_lsun_start = star%luminosity_lsun(star%nz)
      if (star%job%rotation_active) then
         do i = 1,star%nz
          star%old_omega(i) = star%omega(i)
         end do
      endif
end subroutine snapshot_step_start

end subroutine read_starting_model

end module starting_model
