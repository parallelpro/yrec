!----------------------------------------------------------------------
! starting_model
!----------------------------------------------------------------------
! The starting model of each run card (read_starting_model): read from
! file or reused from memory, rescaled and refitted as the card asks,
! and set up with the state the first step needs.
module starting_model
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use turnover_timescale, only: compute_turnover_timescale
      use shell_physics, only: compute_shell_physics
      use model_rescaling, only: rescale_model
      use star_info_lib, only: star_info
      use run_cards, only: card_evolve
      use am_regions, only: am_convective_regions
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use convection_lib, only: compute_temperature_gradients
      use settings_install, only: install_surface_mixture_settings, install_surface_x_settings
      use star_info_lib, only: i_be9, i_c12, i_c13, i_h1, i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18, json, &
           n_species_basic, n_species_extended, n_mix_species, ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, ix_h, ix_o, ix_n, ix_ar, ix_ne, ix_he, &
           max_convective_zones, max_radiative_zones
      use control_slots, only: ichi_dm_max, ichi_dp_env_max
      use envelope_integration, only: atm_get, envint_step_config, fixed_envint_step, &
           ambient_log10_teff
      use envelope_refit, only: append_envelope_points
      use core_messages, only: core_log, core_report, core_stdout, core_stdout_and_log, &
           message_length, unwritten_line
      use run_files, only: starting_model_unit
      use model_files, only: model_file_header, read_model_format_tag, read_mod_model, &
           read_yrec7, read_model2
      use physical_constants, only: c4pi3l, c4pil, cc13, cc23, cgl, csigl, ln10, seconds_per_year
      use eos_lib, only: eos_get, eos_result
      use kap_lib, only: kap_get, num_kap_results
      use stitched_model, only: build_stitched_model
      implicit none
      private
      public :: read_starting_model

! Atomic weights of the surface-mixture elements, in ix_* slot order
! (Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He).
      real(dp), parameter :: atomic_weight(n_mix_species) = [23.0d0, 26.99d0, &
           24.32d0, 55.86d0, 28.1d0, 12.015d0, 1.008d0, 16.0d0, 14.01d0, 39.96d0, 20.19d0, 4.004d0]

! An evolved model's mixing length matches the card's within this
! [-]; a species differs from the centre's by more than this in a model
! that is not homogeneous [mass fraction].
      real(dp), parameter :: mixing_length_match_tol = 2.0d-3, homogeneity_tol = 1.0d-6

! What read_starting_model's steps pass on to each other:
!   run_index              the card (the star's current_card);
!   iread                  the model file's unit (run_files'
!                          starting_model_unit, legacy IFIRST);
!   header                 the model file's header values that the star
!                          does not hold;
!   reference_composition  the central composition of the starting model:
!                          the reference for the homogeneity test that
!                          must pass before the mixture is altered;
!   mixture_ok             that test's verdict;
!   saha_iterations        the Saha iteration count the EOS calls share;
!   num_species            the species carried (11, or 15 with the
!                          extended composition);
!   old_last_shell, old_senv, target_log_mass_at_fit
!                          the envelope refit's old last shell, old
!                          fitting mass and target log m;
!   convective_velocity    the convective velocity of the new last shell
!                          (set_new_last_shell_convection).
      type :: starting_model_state
         integer :: run_index = 0
         integer :: iread = 0
         type(model_file_header) :: header
         real(dp) :: reference_composition(n_species_extended) = 0.0d0
         logical :: mixture_ok = .false.
         integer :: saha_iterations = 0
         integer :: num_species = 0
         integer :: old_last_shell = 0
         real(dp) :: old_senv = 0.0d0
         real(dp) :: target_log_mass_at_fit = 0.0d0
         real(dp) :: convective_velocity = 0.0d0
      end type starting_model_state

contains

!----------------------------------------------------------------------
! read_starting_model (legacy starin)
!----------------------------------------------------------------------
! Sets up the starting model of run card run_index (the star's
! current_card), at the start of every card, in named steps:
!   acquire_starting_model      read the model file (.mod, YREC7 or
!                               MODEL2 layout) when card_from_starting_model is
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
! Called with the star's model_diverged_flag set (a step failed to
! converge), it only stops the run with ierr = 1.
!
! What it sets on the star besides the model: timestep_yr and dt (s)
! from the model file; triangle_orientation, set to 1 after every read
! (reset_state_after_read), as the F77 STARIN set FTRI = 1, whatever
! orientation the file's envelope records give; recompute_envelope_triangle,
! set when the envelope must be recomputed; and total_angular_momentum and
! total_rotational_ke, the rotating model's totals.
!
! Naming: the star% arrays use henyey_iterate's names (luminosity_lsun
! is the linear L/Lsun; m is the unlogged mass coordinate and dm the
! mass of each shell, in g; mean_radius, qiw and mean_gravity are the
! rotation geometry); the
! single-point physics in set_new_last_shell_convection uses
! shell_physics' log10_* names.
!   star reads: convective_turnover_timescale, current_card, fp_rot,
!     ft_rot, log_mass, log_Teff, logP, logR, logRho, logT,
!     luminosity_lsun, model_diverged_flag, nz, photosphere_logP, xa
!   star writes: am_transport_convective_flag, convective_flag,
!     convective_turnover_timescale_prev,
!     envelope_overshoot_bottom_index, envelope_overshoot_top_index,
!     photosphere_logP_prev, substep_time_fraction
!   star reads: none
!   star writes: none
subroutine read_starting_model(s, ierr)
      implicit none
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      type(starting_model_state) :: sm
      character(len=message_length) :: lines(4)
      integer :: radiative_zone_bounds(max_radiative_zones,2)
      integer :: convective_zone_bounds(max_convective_zones,2)
      integer :: num_radiative_zones
      integer :: num_convective_zones
      integer :: jerr

      ierr = 0
      sm%run_index = s%current_card

! model_failed_flag (legacy LARGE) means the last model failed to
! converge: stop the run.
      if (s%model_diverged_flag) then
          lines = unwritten_line
          write(lines,1000)
          call core_log(lines)
 1000       format(1x,39('>'),40('<')/, &
                  "STARIN:        ***** RUN STOPPED *****")
          lines = unwritten_line
          write(lines,1010)
          call core_log(lines)
 1010       format("STARIN: ***** MODEL FAILED TO CONVERGE *****")
            ! run_yrec returns the error; the command-line wrapper
            ! (main) stops.
            ierr = 1
            return
      end if

! The model file is on run_files' starting_model_unit (legacy IFIRST).
      sm%iread = starting_model_unit

! Initial value: no convective velocity.
      sm%convective_velocity = 0.0d0

      call acquire_starting_model(s, sm, ierr)
      if (ierr /= 0) return
      call extend_core_toward_center(s, sm, ierr)
      if (ierr /= 0) return
      call rescale_and_refit_envelope(s, sm, ierr)
      if (ierr /= 0) return
      call build_shell_masses(s)
      call initialize_rotation_state(s, ierr)
      if (ierr /= 0) return
      call update_surface_mixture(s, sm, ierr)
      call snapshot_step_start(s)
      if (ierr /= 0) return
      call compute_shell_physics(s, s%fp_rot,s%ft_rot, &
           s%xa,s%logRho,s%luminosity_lsun,s%logP, &
           s%logR,s%log_mass,s%logT,s%convective_flag,s%nz, &
           s%log_Teff, jerr)
      if (jerr /= 0) then
      ! Pass the error on.
         ierr = jerr
         return
      end if
      call am_convective_regions(s%xa,s%logRho,s%logP,s%logR,s%log_mass, &
                 s%logT,s%convective_flag,s%nz, &
                 s%am_transport_convective_flag,radiative_zone_bounds, &
                 convective_zone_bounds,num_radiative_zones, &
                 num_convective_zones, s%envelope_overshoot_bottom_index, s%envelope_overshoot_top_index, ierr)
      if (ierr /= 0) return
! Initialize the convective turnover time (legacy TAUCZ), the
! photospheric pressure photosphere_logP and substep_time_fraction. The
! stitched model (interior plus a fresh envelope) is built first: it sets
! star%photosphere_logP, and compute_turnover_timescale walks it,
! returning the radius of the convection-zone base in
! star%envelope_radius_rsun. (starin passed the envelope radius array
! there, overwriting its first element.) The wind takes photosphere_logP
! and the turnover time at the fraction substep_time_fraction of the step,
! between the step-start copies (*_prev) and the current values.
      call build_stitched_model(s, ierr)
      if (ierr /= 0) return
      call compute_turnover_timescale(s, ierr)
      if (ierr /= 0) return
      s%convective_turnover_timescale_prev = s%convective_turnover_timescale
      s%photosphere_logP_prev = s%photosphere_logP
      s%substep_time_fraction = 0.5

      return

end subroutine read_starting_model



! ---------------------------------------------------------------
! Obtain the starting model. With card_from_starting_model(run_index) false,
! reuse the model in memory (legacy LFIRST). With it true, read the
! model file on unit iread, warn if its mixing length or composition
! set differs from the controls, reset the step state, set the
! rotation curve and, in a chemically homogeneous model, alter the
! heavy-element mixture and isotope ratios if asked. Sets ierr on
! failure.
!   star reads: card_from_starting_model
!   star writes: none
subroutine acquire_starting_model(s, sm, ierr)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr

! A reused model (card_from_starting_model false) stays as it is in memory. On a
! rescale card it is rescaled once, in rescale_and_refit_envelope, as a
! model read from file is. (starin rescaled a reused model here as well,
! so it was rescaled twice.)
      if (s%card_from_starting_model(sm%run_index)) then
         call read_starting_model_file(s, sm, ierr)
         if (ierr /= 0) return
         call check_model_against_controls(s, sm)
         call reset_state_after_read(s)
         call set_starting_rotation_curve(s)
! Leave unit iread open, rewound.
         rewind sm%iread
! Alter the heavy-element mixture (MHP 10/24), only in a chemically
! homogeneous model read from file. change_cno_mixture_active (legacy
! ISETMIX) sets the C, N and O fractions of Z; change_isotope_ratios_active
! (legacy ISETISO) sets the isotope ratios and light-element
! abundances. Between them they can overwrite mass fractions 4-15.
! (The test below is always true inside this branch; starin once also
! required the first card.)
         if (s%card_from_starting_model(sm%run_index)) then
            call check_mixture_homogeneity(s, sm)
            if (sm%mixture_ok) then
               call rescale_cno_mixture(s, sm)
               call set_isotope_ratios(s, sm)
            end if
         end if
      end if
end subroutine acquire_starting_model


! ---------------------------------------------------------------
! Read the model file on unit iread, in the layout its first keyword
! names (YMOD: YREC .mod; NMOD: YREC7; MOD2: MODEL2), and set the
! timestep and star%log_total_mass from it. The readers fill the star%
! model directly and return the file's other header values in header
! (mixing length, timestep, flags, physics codes), which
! check_model_against_controls compares with the controls. (The file's
! triangle orientation is not used: reset_state_after_read sets it to 1
! and asks for a fresh envelope.) Sets ierr on an unknown layout or a
! failed read (an empty file among them).
!   star reads: none
!   star writes: dt, timestep_yr
subroutine read_starting_model_file(s, sm, ierr)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr
      character(len=4) :: format_tag
      character(len=message_length) :: lines(4)

      call read_model_format_tag(sm%iread, format_tag, ierr)
      if (ierr /= 0) return

      if (format_tag == 'YMOD') then
         lines = unwritten_line
         write(lines,11)
         call core_log(lines)
 11      format(1x,'loading starting model (YREC .mod format)')
         call read_mod_model(s, sm%iread, sm%header, ierr)
         if (ierr /= 0) return
      else if (format_tag == 'NMOD') then
         lines = unwritten_line
         write(lines,12)
         call core_log(lines)
 12      format(1x,'loading starting model (YREC7 format)')
         call read_yrec7(s, sm%iread, sm%header, ierr)
         if (ierr /= 0) return
      else if (format_tag == 'MOD2 ') then
         lines = unwritten_line
         write(lines,16)
         call core_log(lines)
 16      format(1x,'loading starting model (MODEL2 format)')
         call read_model2(s, sm%iread, sm%header, ierr)
         if (ierr /= 0) return
      else
         lines = unwritten_line
         write(lines,20)
         call core_log(lines)
 20      format('STARIN: ***** RUN TERMINATED, INVALID INPUT', &
                ' MODEL FILE.  *****')
         ! run_yrec returns the error; the command-line wrapper (main)
         ! stops.
         ierr = 1
         return
      end if

! The timestep in s, and the log total mass.
      s%timestep_yr = sm%header%timestep_yr
      s%dt = seconds_per_year*s%timestep_yr
end subroutine read_starting_model_file


! ---------------------------------------------------------------
! Warn when an evolved model (model_number > 0) has a mixing length or
! extended-composition flag that differs from the controls. A warning
! only: the run goes on.
!   star reads: mixing_length_alpha, model_number,
!     physics%use_extended_composition
!   star writes: none
subroutine check_model_against_controls(s, sm)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      logical :: mixing_length_matches
! lexcp0 is never assigned: the mismatch warning below prints it where it
! presumably means the model's flag, header%use_extended_composition (a
! bug kept from starin; core.md finding 3, on the B list).
      logical :: lexcp0
      character(len=message_length) :: lines(4)
      if (s%model_number > 0) then
       mixing_length_matches = .true.
! The tolerance 2e-3 allows for model files storing the mixing length
! to four significant digits; a model mixing length <= 0 is not
! checked.
       if (sm%header%mixing_length > 0.0d0) mixing_length_matches = &
            (abs(s%mixing_length_alpha-sm%header%mixing_length) < mixing_length_match_tol)
! The warning prints lexcp0, which is never set (see its declaration),
! for the model's flag.
       if (.not.mixing_length_matches .or. sm%header%use_extended_composition.neqv. &
            s%physics%use_extended_composition) then
          lines = unwritten_line
          write(lines,1040) s%mixing_length_alpha,sm%header%mixing_length, &
               s%physics%use_extended_composition,lexcp0
          call core_report(lines)
 1040       format(1x,'ERROR IN SUBROUTINE STARIN'/1x,'USER PARAMETERS', &
             ' OF WRONG TYPE FOR INITIAL MODEL'/1x,'MIXING LENGTH - USER' &
             ,' DESIRES',f7.3,' MODEL MIX LENGTH',f7.3/1x, &
             'EXTENDED COMP-USER DESIRES ',l1,' MODEL USED ',l1)
       end if
      end if
end subroutine check_model_against_controls


! ---------------------------------------------------------------
! Reset the store and trial-sign flags, require a fresh envelope, and
! take the surface X and Z from the outermost interior point.
!   star reads: nz, xa
!   star writes: recompute_envelope_triangle, stored_vertex_index,
!     surface_hydrogen_fraction, surface_metal_fraction,
!     triangle_orientation
subroutine reset_state_after_read(s)
      type(star_info), intent(inout) :: s
      s%stored_vertex_index = 0
      s%triangle_orientation = 1.0d0
! Recompute the envelope: the envelope triangle stored with the model
! may not be reliable.
       s%recompute_envelope_triangle = .true.


! surface_hydrogen_fraction, surface_metal_fraction: the outermost point's X, Z.
      s%surface_hydrogen_fraction = s%xa(i_h1,s%nz)
      s%surface_metal_fraction = s%xa(i_metals,s%nz)
end subroutine reset_state_after_read


! ---------------------------------------------------------------
! Set the solid-body rotation curve when asked (rotating runs), or
! unit rotation factors fp and ft (non-rotating runs).
!   star reads: nz, physics%initial_omega, physics%rotation_active,
!     physics%set_initial_omega
!   star writes: fp_rot, ft_rot, omega
subroutine set_starting_rotation_curve(s)
      type(star_info), intent(inout) :: s
      integer :: i

      if (s%physics%rotation_active) then
       if (s%physics%set_initial_omega) then
! Solid-body rotation, omega = initial_omega (legacy WNEW) in every
! shell: this turns a non-rotating model into a rotating one.
          do i = 1,s%nz
             s%omega(i) = s%physics%initial_omega
          end do
       end if
      else
         do i = 1,json
            s%fp_rot(i) = 1.0d0
            s%ft_rot(i) = 1.0d0
         end do
      end if
end subroutine set_starting_rotation_curve


! ---------------------------------------------------------------
! Check that the starting model is chemically homogeneous (every shell
! within 1e-6 of the center in each of species 1-15) before its mixture
! is altered; mixture_ok is false, with a message, if not.
!   star reads: control%change_cno_mixture_active,
!     control%change_isotope_ratios_active, nz, xa
!   star writes: none
subroutine check_mixture_homogeneity(s, sm)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      real(dp) :: fraction_diff
      integer :: i
      integer :: j
      character(len=message_length) :: lines(4)
      sm%mixture_ok = .true.
      if (s%control%change_cno_mixture_active .or. s%control%change_isotope_ratios_active) then
         do i = 1, n_species_extended
            sm%reference_composition(i)=s%xa(i,1)
         end do
         homogeneity: do j = 2,s%nz
            do i = 1, n_species_extended
               fraction_diff = abs(s%xa(i,j)-sm%reference_composition(i))
               if (fraction_diff > homogeneity_tol) then
                  lines = unwritten_line
                  write(lines,592)i,j,fraction_diff
                  call core_stdout_and_log(lines)
 592              format('SPECIES ',i3,' IN SHELL ',i5, &
                    ' DIFFERS FROM CENTER BY ',e12.4, &
                    ' MIX NOT MODIFIED IN EVOLVED MODEL')
                  sm%mixture_ok = .false.
                  exit homogeneity
               end if
            end do
         end do homogeneity
      end if
end subroutine check_mixture_homogeneity


! ---------------------------------------------------------------
! Scale the C, N and O isotopes in every shell so that C, N and O make
! up the requested fractions of Z. Z itself is not changed.
!   star reads: control%change_cno_mixture_active,
!     control%target_carbon_cno_fraction,
!     control%target_nitrogen_cno_fraction,
!     control%target_oxygen_cno_fraction, nz
!   star writes: xa
subroutine rescale_cno_mixture(s, sm)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer :: i
      integer :: j
      integer :: k
      real(dp) :: total_carbon_cno_fraction
      real(dp) :: total_nitrogen_cno_fraction
      real(dp) :: total_oxygen_cno_fraction
      real(dp) :: carbon_scale_ratio
      real(dp) :: nitrogen_scale_ratio
      real(dp) :: oxygen_scale_ratio
      character(len=message_length) :: lines(4)
      if (s%control%change_cno_mixture_active) then
! The current fractions, from the (homogeneous) central composition:
!   f_C = (X(C12) + X(C13))/Z,  f_N = (X(N14) + X(N15))/Z,
!   f_O = (X(O16) + X(O17) + X(O18))/Z.
! Every isotope of an element is scaled by target/f, which keeps the
! isotope ratios; set_isotope_ratios may then change them.
         total_carbon_cno_fraction = (sm%reference_composition(i_c12)+ &
              sm%reference_composition(i_c13))/sm%reference_composition(i_metals)
         total_nitrogen_cno_fraction = (sm%reference_composition(i_n14)+ &
              sm%reference_composition(i_n15))/sm%reference_composition(i_metals)
         total_oxygen_cno_fraction = (sm%reference_composition(i_o16)+ &
              sm%reference_composition(i_o17)+sm%reference_composition(i_o18))/ &
              sm%reference_composition(i_metals)
         carbon_scale_ratio = s%control%target_carbon_cno_fraction/ &
              total_carbon_cno_fraction
         nitrogen_scale_ratio = s%control%target_nitrogen_cno_fraction/ &
              total_nitrogen_cno_fraction
         oxygen_scale_ratio = s%control%target_oxygen_cno_fraction/ &
              total_oxygen_cno_fraction
         do i = i_c12,i_c13
            do j = 1,s%nz
               s%xa(i,j)=carbon_scale_ratio*s%xa(i,j)
            end do
         end do
         do i = i_n14,i_n15
            do j = 1,s%nz
               s%xa(i,j)=nitrogen_scale_ratio*s%xa(i,j)
            end do
         end do
         do i = i_o16,i_o18
            do j = 1,s%nz
               s%xa(i,j)=oxygen_scale_ratio*s%xa(i,j)
            end do
         end do
         lines = unwritten_line
         write(lines,594)(sm%reference_composition(k),k=i_c12,i_o18), &
              (s%xa(k,1),k=i_c12,i_o18)
         call core_stdout_and_log(lines)
 594     format(1x,'CNO mixture applied to the starting model:',/, &
             4x,'old (C12 C13 N14 N15 O16 O17 O18):',7es12.4,/, &
             4x,'new (C12 C13 N14 N15 O16 O17 O18):',7es12.4)
      end if
end subroutine rescale_cno_mixture


! ---------------------------------------------------------------
! Set the requested He3, H2, Li6, Li7 and Be9 abundances and the
! C12/C13 and O16/O18 ratios in every shell. The central sums
! S_C = X(C12) + X(C13) and S_O = X(O16) + X(O18) are kept and split
! by the requested ratio r:
!   X(C12) = r_C S_C/(1 + r_C),  X(C13) = S_C/(1 + r_C),
! and likewise for O16 and O18. N15 and O17 are not altered, and B10
! and B11 are not species of the model.
!   star reads: control%c12_to_c13_ratio,
!     control%change_isotope_ratios_active,
!     control%initial_be9_fraction, control%initial_h2_fraction,
!     control%initial_he3_fraction, control%initial_li6_fraction,
!     control%initial_li7_fraction, control%o16_to_o18_ratio, nz
!   star writes: xa
subroutine set_isotope_ratios(s, sm)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer :: j
      integer :: k
      real(dp) :: sum_c12_c13
      real(dp) :: sum_o16_o18
      character(len=message_length) :: lines(4)
      if (s%control%change_isotope_ratios_active) then
         sum_c12_c13 = s%xa(i_c12,1)+s%xa(i_c13,1)
         sum_o16_o18 = s%xa(i_o16,1)+s%xa(i_o18,1)
         do j = 1,s%nz
            s%xa(i_he3,j)=s%control%initial_he3_fraction
            s%xa(i_c12,j)= s%control%c12_to_c13_ratio*sum_c12_c13/ &
                 (1.0d0+s%control%c12_to_c13_ratio)
            s%xa(i_c13,j)= sum_c12_c13/(1.0d0+s%control%c12_to_c13_ratio)
            s%xa(i_o16,j)= s%control%o16_to_o18_ratio*sum_o16_o18/ &
                 (1.0d0+s%control%o16_to_o18_ratio)
            s%xa(i_o18,j)= sum_o16_o18/(1.0d0+s%control%o16_to_o18_ratio)
            s%xa(i_h2,j)=s%control%initial_h2_fraction
            s%xa(i_li6,j)=s%control%initial_li6_fraction
            s%xa(i_li7,j)=s%control%initial_li7_fraction
            s%xa(i_be9,j)=s%control%initial_be9_fraction
         end do
         lines = unwritten_line
         write(lines,593)(sm%reference_composition(k),k=i_he3,i_be9), &
              (s%xa(k,1),k=i_he3,i_be9)
         call core_stdout_and_log(lines)
 593     format(1x,'isotope/light-element mixture applied to the', &
             ' starting model:',/, &
             4x,'old (He3 C12 C13 N14 N15 O16 O17 O18 H2 Li6 Li7 Be9):', &
             12es12.4,/, &
             4x,'new (He3 C12 C13 N14 N15 O16 O17 O18 H2 Li6 Li7 Be9):', &
             12es12.4)
      end if
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
!   star reads: card_from_starting_model, physics%chi_grid_scale,
!     physics%extend_core_inward
!   star writes: convective_flag, core_mass_reduction, log_mass, logP,
!     logR, logRho, logT, luminosity_lsun, num_core_points_added, nz,
!     omega, xa
subroutine extend_core_toward_center(s, sm, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr
      integer :: i
      integer :: j
      integer :: num_shells_extended
      real(dp) :: core_shell_spacing
      real(dp) :: central_log_density
      real(dp) :: actual_gradient
      real(dp) :: central_shell_luminosity
      real(dp) :: density_estimate_offset
      integer :: first_original_shell
      real(dp) :: trial_log_pressure
      real(dp) :: temp_scratch
      character(len=message_length) :: lines(4)

      if (s%physics%extend_core_inward .and. s%card_from_starting_model(sm%run_index)) then
! The innermost mass is to drop by the factor core_mass_reduction_factor
! (legacy FCORE): that takes N = int(log10(factor)/dlogm) + 1 steps
! (num_core_points_added, legacy MCORE). core_mass_reduction is then
! overwritten with the log10 shift N dlogm actually made.
          s%num_core_points_added = int(log10(s%core_mass_reduction)/ &
               s%physics%chi_grid_scale(ichi_dm_max))+1
          s%core_mass_reduction = dble(s%num_core_points_added)* &
               s%physics%chi_grid_scale(ichi_dm_max)
          num_shells_extended = s%nz + s%num_core_points_added
! Too many points for the json-point arrays: the run stops with ierr
! (starin printed RUN TERMINATED but went on writing past json; fixed
! 2026-10).
          if (num_shells_extended > json) then
             lines = unwritten_line
             write(lines,476)"STARIN: Unable to extend core inward ", &
                   "- JSON too small"
             call core_log(lines)
             lines = unwritten_line
             write(lines,477) "STARIN: Required size =", &
                    num_shells_extended, ", JSON = ", json
             call core_log(lines)
             lines = unwritten_line
             write(lines,478) "STARIN:  ***** RUN TERMINATED *****"
             call core_log(lines)
  476        format(2a)
  477        format(a, i8, a, i8)
  478        format(a)
             ierr = 1
             return
          end if
          core_shell_spacing = s%physics%chi_grid_scale(ichi_dm_max)
! Move every old point up by N slots (logRho is not moved).
          do i=s%nz,1, -1
             s%log_mass(i+s%num_core_points_added) = s%log_mass(i)
             s%logR(i+s%num_core_points_added) = s%logR(i)
             s%luminosity_lsun(i+s%num_core_points_added) = s%luminosity_lsun(i)
             s%logP(i+s%num_core_points_added) = s%logP(i)
             s%logT(i+s%num_core_points_added) = s%logT(i)
             s%convective_flag(i+s%num_core_points_added) = s%convective_flag(i)
             do j = 1, n_species_extended
                s%xa(j,i+s%num_core_points_added) = s%xa(j,i)
             end do
             s%omega(i+s%num_core_points_added) = s%omega(i)
          end do
          first_original_shell = s%num_core_points_added+1
! March inward in mass from the innermost old point (now slot N+1).
! Slots 1 and 2 still hold the old innermost values here, so rho_c is
! the old central density and del = d log T/d log P is taken between
! the two innermost old points (for N = 1, slot 2 has already been
! overwritten and del is 0/0).
          central_log_density = s%logRho(1)
          actual_gradient = (s%logT(2)-s%logT(1))/ &
               (s%logP(2)-s%logP(1))
          central_shell_luminosity = s%luminosity_lsun(1)
! Offset for estimating rho from P and T at fixed mean molecular weight
! (MHP 4/12): log rho = log P - log T - offset, with the offset taken
! at slot N+1.
          density_estimate_offset = s%logP(first_original_shell)- &
               s%logRho(first_original_shell)- &
               s%logT(first_original_shell)
          do i = s%num_core_points_added,1,-1
             s%log_mass(i) = s%log_mass(i+1)-core_shell_spacing
! Constant density, m = (4 pi/3) rho_c r^3:
!   log r = (1/3)(log m - log(4 pi/3) - log rho_c).
             s%logR(i) = cc13*(s%log_mass(i)-c4pi3l-central_log_density)
! Constant epsilon: L = L_1 m/m_1, L_1 and m_1 at the innermost old
! point.
             s%luminosity_lsun(i) = exp(ln10*(s%log_mass(i)- &
                  s%log_mass(first_original_shell)))*central_shell_luminosity
! Hydrostatic equilibrium, d log P/d log m = -G m^2/(4 pi r^4 P), taken
! inward in two stages: log P half a step in from point i+1, then the
! full step with that pressure:
!   log P_i = log P_(i+1) + (G m_i^2/(4 pi r_i^4 P_half)) dlogm.
             trial_log_pressure = s%logP(i+1)
             temp_scratch = exp(ln10*(cgl+2.0d0*s%log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*s%logR(i)))
             trial_log_pressure = trial_log_pressure+0.5d0*temp_scratch* &
                  core_shell_spacing
             temp_scratch = exp(ln10*(cgl+2.0d0*s%log_mass(i)-c4pil- &
                  trial_log_pressure-4.0d0*s%logR(i)))
             s%logP(i) = s%logP(i+1)+temp_scratch* &
                  core_shell_spacing
! Convective or radiative as the innermost old point.
             s%convective_flag(i) = s%convective_flag(first_original_shell)
! Constant del: log T_i = log T_(i+1) + del (log P_i - log P_(i+1)).
             s%logT(i) = s%logT(i+1)+temp_scratch* &
                  actual_gradient*core_shell_spacing
! Omega and composition as at the innermost old point.
             s%omega(i) = s%omega(first_original_shell)
             do j = 1, n_species_extended
                s%xa(j,i) = s%xa(j,first_original_shell)
             end do
! The local density estimate, in place of an EOS call that did not
! work here (MHP 4/12).
             s%logRho(i) = s%logP(i) - s%logT(i) - &
                  density_estimate_offset
          end do
          s%nz = num_shells_extended
      end if
end subroutine extend_core_toward_center


! ---------------------------------------------------------------
! On a rescale card (card_kind /= card_evolve), rescale the model's X, Z or
! mass (rescale_model). Then, with fit_point_move_requested, move the
! envelope fitting point to requested_log_fit_point_mass_fraction = log10(m_fit/M)
! (legacy SENV0), on the first model of every card (DBG 2/92):
! - deeper than now: drop the points outside the new fitting mass and
!   put the new last point at it, linear in log m between the old
!   points on either side, with the composition of the one inside;
! - shallower: integrate an envelope with the model's Teff and L from
!   the surface down to the old last point, and append its points
!   inside the new fitting mass as shells, with the old last point's
!   composition apart from X and Z.
! Sets ierr on failure.
!   star reads: card_kind, fit_point_move_requested,
!     log_fit_point_mass_fraction, physics%use_extended_composition,
!     requested_log_fit_point_mass_fraction
!   star writes: recompute_envelope_triangle
subroutine rescale_and_refit_envelope(s, sm, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr
      integer :: i
      character(len=message_length) :: lines(4)

      if (s%card_kind(sm%run_index) /= card_evolve) call rescale_model(s, ierr)
      if (ierr /= 0) return

! A once-through loop: exit leaves the envelope change early.
      if (s%fit_point_move_requested) then
      envelope_rescale: do
       call clamp_requested_envelope_mass(s, sm)
       if (s%log_fit_point_mass_fraction == s%requested_log_fit_point_mass_fraction) exit envelope_rescale
       sm%num_species = n_species_basic
       if (s%physics%use_extended_composition) sm%num_species = n_species_extended
       if (s%requested_log_fit_point_mass_fraction < s%log_fit_point_mass_fraction) then
! The new envelope is deeper than the old one.
          call find_shell_below_new_fit(s, sm)
          if (i < (1)) then
! The requested fitting mass lies inside the innermost shell: report
! it and leave the envelope mass unchanged.
          lines = unwritten_line
          write(lines,576)s%requested_log_fit_point_mass_fraction
          call core_log(lines)
 576        format(5x,'ERROR IN SUBROUTINE STARIN'/5x,'DESIRED', &
              ' ENVELOPE MASS',1pe22.13,' TOO LARGE'/5x,'ENVELOPE', &
              ' MASS NOT CHANGED')
          exit envelope_rescale
          end if
          call interpolate_new_fitting_point(s, sm)
          call set_new_last_shell_convection(s, sm, ierr)
          if (ierr /= 0) return
       else
! The new envelope is shallower than the old one.
          call integrate_envelope_to_new_fit(s, sm, ierr)
          if (ierr /= 0) return
          call append_envelope_as_shells(s, sm, ierr)
          if (ierr /= 0) return
       end if
       s%recompute_envelope_triangle = .true.
       lines = unwritten_line
       write(lines,597)sm%old_senv,s%log_fit_point_mass_fraction
       call core_log(lines)
 597     format(5x,'***** NEW ENVELOPE MASS CALCULATED *****'/8x, &
              'OLD SENV ',1pe22.13,'  NEW SENV',e22.13)
      exit envelope_rescale
      end do envelope_rescale
      end if

end subroutine rescale_and_refit_envelope


! ---------------------------------------------------------------
! Make the requested envelope mass, log10(m_fit/M), negative and no
! shallower than -1e-12, and record the current one, log m(nz) - log M,
! in log_fit_point_mass_fraction and the local old_senv.
!   star reads: log_mass, log_total_mass, nz
!   star writes: log_fit_point_mass_fraction,
!     requested_log_fit_point_mass_fraction
subroutine clamp_requested_envelope_mass(s, sm)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
       if (s%requested_log_fit_point_mass_fraction > 0.0d0) s%requested_log_fit_point_mass_fraction = &
            -s%requested_log_fit_point_mass_fraction
! The limit -1e-12 stands for the base of the atmosphere, tau = 2/3
! (DBG 2/92; it was -1e-10).
       if (s%requested_log_fit_point_mass_fraction > -1.0d-12) s%requested_log_fit_point_mass_fraction = &
            -1.d-12
       s%log_fit_point_mass_fraction = s%log_mass(s%nz) - s%log_total_mass
       sm%old_senv = s%log_fit_point_mass_fraction
end subroutine clamp_requested_envelope_mass


! ---------------------------------------------------------------
! Find the outermost shell i inside the new (deeper) fitting mass,
! log m_fit = log M + requested_log_fit_point_mass_fraction; i = 0 when the fitting
! mass lies inside the innermost shell.
!   star reads: log_mass, log_total_mass, nz,
!     requested_log_fit_point_mass_fraction
!   star writes: none
subroutine find_shell_below_new_fit(s, sm)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer :: i
          sm%target_log_mass_at_fit = s%log_total_mass+s%requested_log_fit_point_mass_fraction
          do i = s%nz-1,1,-1
             if (s%log_mass(i) < sm%target_log_mass_at_fit) exit
          end do
end subroutine find_shell_below_new_fit


! ---------------------------------------------------------------
! Make shell i+1 the new last interior shell, at the fitting mass,
! with the composition of shell i and the structure (and omega) linear
! in log m between shells i and i+1:
!   f = (log m_fit - log m_i)/(log m_(i+1) - log m_i),
!   y = y_i + f (y_(i+1) - y_i) for log rho, L, log P, log R, log T.
! Its X and Z become surface_hydrogen_fraction and surface_metal_fraction.
!   star reads: physics%rotation_active,
!     requested_log_fit_point_mass_fraction
!   star writes: log_fit_point_mass_fraction, log_mass, logP, logR,
!     logRho, logT, luminosity_lsun, nz, omega,
!     surface_hydrogen_fraction, surface_metal_fraction, xa
subroutine interpolate_new_fitting_point(s, sm)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer :: i
      integer :: j
      real(dp) :: interior_interp_fraction
            s%nz = i + 1
          s%log_fit_point_mass_fraction = s%requested_log_fit_point_mass_fraction
          interior_interp_fraction = (sm%target_log_mass_at_fit-s%log_mass(i))/ &
               (s%log_mass(i+1) - s%log_mass(i))
          s%log_mass(s%nz) = sm%target_log_mass_at_fit
          s%logRho(s%nz) = s%logRho(i) + &
               interior_interp_fraction*(s%logRho(i+1) - s%logRho(i))
          s%luminosity_lsun(s%nz) = s%luminosity_lsun(i) + &
               interior_interp_fraction*(s%luminosity_lsun(i+1) - &
               s%luminosity_lsun(i))
          s%logP(s%nz) = s%logP(i) + &
               interior_interp_fraction*(s%logP(i+1) - s%logP(i))
          s%logR(s%nz) = s%logR(i) + &
               interior_interp_fraction*(s%logR(i+1) - s%logR(i))
          s%logT(s%nz) = s%logT(i) + &
               interior_interp_fraction*(s%logT(i+1) - &
               s%logT(i))
          do j = 1,sm%num_species
             s%xa(j,s%nz) = s%xa(j,i)
          end do
          s%surface_hydrogen_fraction = s%xa(i_h1,s%nz)
          s%surface_metal_fraction = s%xa(i_metals,s%nz)
          if (s%physics%rotation_active) s%omega(s%nz) = s%omega(i) + &
               interior_interp_fraction*(s%omega(i+1)-s%omega(i))
end subroutine interpolate_new_fitting_point


! ---------------------------------------------------------------
! Decide whether the new last shell is convective: as its neighbours
! when they agree; otherwise from the temperature gradients at the
! point (eos_get with the atmosphere EOS, kap_get,
! compute_temperature_gradients, no rotation), which also give its
! density.
!   star reads: log_mass, log_Teff, logP, logR, logT, luminosity_lsun,
!     nz, xa
!   star writes: convective_flag, logRho
subroutine set_new_last_shell_convection(s, sm, ierr)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr
      integer :: i
      real(dp) :: actual_gradient
      logical :: want_derivatives
      logical :: in_atmosphere
      real(dp) :: hydrogen_fraction
      real(dp) :: metal_fraction
      real(dp) :: log10_pressure
      real(dp) :: log10_temperature
      real(dp) :: log10_density
      real(dp) :: log10_radius
      real(dp) :: log10_mass
      real(dp) :: shell_luminosity_lsun
      real(dp) :: point_pressure_rotation_factor
      real(dp) :: point_temperature_rotation_factor
      type(eos_result) :: eos_res
      real(dp) :: kap_res(num_kap_results)
      real(dp) :: radiative_gradient
      real(dp) :: dgrad_dt_component
      real(dp) :: dgrad_dp_component
      real(dp) :: dgrad_dr_component
      logical :: is_convective
      integer :: jerr
          if (s%convective_flag(i).and.s%convective_flag(i+1)) then
             s%convective_flag(s%nz) = .true.
          else if (.not.s%convective_flag(i).and. .not.s%convective_flag(i+1)) &
               then
             s%convective_flag(s%nz) = .false.
          else
             want_derivatives = .false.
             in_atmosphere = .true.
             sm%saha_iterations = 0
             hydrogen_fraction = s%xa(i_h1,s%nz)
             metal_fraction = s%xa(i_metals,s%nz)
             log10_pressure = s%logP(s%nz)
             log10_temperature = s%logT(s%nz)
             log10_density = s%logRho(s%nz)
             shell_luminosity_lsun = s%luminosity_lsun(s%nz)
             log10_radius = s%logR(s%nz)
             log10_mass = s%log_mass(s%nz)
             point_pressure_rotation_factor = 1.0d0
             point_temperature_rotation_factor = 1.0d0
               eos_res%log10_density = log10_density
               call eos_get(log10_temperature, log10_pressure, &
                    hydrogen_fraction, metal_fraction, eos_res, &
                    want_derivatives, in_atmosphere, sm%saha_iterations, &
                    composition_at_zone=s%xa(:,s%nz), ierr=ierr)
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
                    sm%convective_velocity, want_derivatives, is_convective, &
                    point_pressure_rotation_factor, &
                    point_temperature_rotation_factor, s%log_Teff, .false., jerr)
               if (jerr /= 0) then
               ! Pass the error on.
                  ierr = jerr
                  return
               end if
               log10_density = eos_res%log10_density
               s%logRho(s%nz) = log10_density
               s%convective_flag(s%nz) = is_convective
          end if
end subroutine set_new_last_shell_convection


! ---------------------------------------------------------------
! Integrate an envelope with the model's Teff and L from the surface
! down to the pressure of the old last interior shell (atm_get), for
! the new, shallower fitting mass; sets log_fit_point_mass_fraction to
! it. Sets ierr on failure.
!   star reads: envelope_cz_bottom_index, log10_solar_luminosity,
!     log_L_lsun, log_Teff, log_total_mass, logP, nz,
!     physics%chi_grid_scale, physics%spot_filling_factor,
!     physics%spot_temp_contrast, requested_log_fit_point_mass_fraction,
!     xa
!   star writes: log_fit_point_mass_fraction
subroutine integrate_envelope_to_new_fit(s, sm, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr
      integer :: jerr_atm
      real(dp) :: hydrogen_fraction
      real(dp) :: metal_fraction
      real(dp) :: log10_radius
      real(dp) :: shell_luminosity_lsun
      real(dp) :: point_pressure_rotation_factor
      real(dp) :: point_temperature_rotation_factor
      type(envint_step_config) :: env_steps
      logical :: save_boundary_flag
      logical :: print_flag
      real(dp) :: log10_gravity
      real(dp) :: log10_pressure_limit
      real(dp) :: spot_adjusted_log_teff
            sm%old_last_shell = s%nz
! atm_get integrates with one fixed step for start, minimum and
! maximum: chi_grid_scale(ichi_dp_env_max) (legacy HPTTOL(8)).
            env_steps = fixed_envint_step(s%physics%chi_grid_scale(ichi_dp_env_max))
          save_boundary_flag = .false.
          print_flag = .true.
          sm%saha_iterations = 0
! Surface radius and gravity from L = 4 pi R^2 sigma Teff^4 and
! g = G M/R^2:
!   log R = (log L/Lsun + log Lsun - 4 log Teff - log(4 pi) - log sigma)/2,
! with log Lsun = star%log10_solar_luminosity (starin's CLSUNL; the
! 2026-08 modernization had put Lsun in erg/s there, fixed 2026-10-12).
          shell_luminosity_lsun = exp(ln10*s%log_L_lsun)
          log10_radius = 0.5d0*(s%log_L_lsun + s%log10_solar_luminosity - &
               4.0d0*s%log_Teff - c4pil - csigl)
          log10_gravity = cgl + s%log_total_mass - log10_radius - log10_radius
          hydrogen_fraction = s%xa(i_h1,s%nz)
          metal_fraction = s%xa(i_metals,s%nz)
          point_pressure_rotation_factor = 1.0d0
          point_temperature_rotation_factor = 1.0d0
          log10_pressure_limit = s%logP(s%nz)
! Spotted stars (G. Somers 10/14): with spot_filling_factor f /= 0,
! spot_temp_contrast x /= 1 and envelope_cz_bottom_index = nz, the
! envelope is integrated at the Teff of the unspotted surface
! (ambient_log10_teff, legacy ATEFFL):
!   Teff_amb^4 = Teff^4/(f x^4 + 1 - f).
          if (s%envelope_cz_bottom_index == s%nz.and.s%physics%spot_filling_factor /= &
               0.0.and.s%physics%spot_temp_contrast /= 1.0) then
               spot_adjusted_log_teff = ambient_log10_teff(s, s%log_Teff)
          else
             spot_adjusted_log_teff = s%log_Teff
          end if
          call atm_get(s, shell_luminosity_lsun,point_pressure_rotation_factor, &
                 point_temperature_rotation_factor,log10_gravity,s%log_total_mass, &
                 print_flag,save_boundary_flag, &
                 log10_pressure_limit,log10_radius,spot_adjusted_log_teff, &
                 hydrogen_fraction,metal_fraction,s%xa(:,s%nz),sm%saha_iterations, &
                 ierr=jerr_atm,env_steps=env_steps)
! atm_get returns integration failures (envelope_termination among
! them) in jerr_atm; pass them on.
          if (jerr_atm /= 0) then
             ierr = jerr_atm
             return
          end if
          s%log_fit_point_mass_fraction = s%requested_log_fit_point_mass_fraction
end subroutine integrate_envelope_to_new_fit


! ---------------------------------------------------------------
! Make the fresh envelope continuous with the interior and append its
! points inside the new fitting mass as shells; in rotating models they
! take the old last shell's omega. Echoes the new points. Sets ierr on
! failure.
!   star reads: log_total_mass, physics%rotation_active
!   star writes: convective_flag, env_points, j_rot, log_mass, logP,
!     logR, logRho, logT, luminosity_lsun, nz, omega, xa
subroutine append_envelope_as_shells(s, sm, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr
      integer :: j
      real(dp) :: pressure_offset
      real(dp) :: density_offset
      real(dp) :: temperature_offset
      real(dp) :: radius_offset
      character(len=message_length), allocatable :: point_lines(:)
      character(len=message_length) :: lines(4)
! Stop if the envelope points would overflow the json-point arrays.
            if (s%nz+s%env_points%num_env_points >= json) then
               lines = unwritten_line
               write(lines,*) 'read_starting_model: grid overflow (json too small)'
               call core_stdout(lines)
               ierr = 1
               return
            end if
! Envelope point 1 is the deepest, at the old last shell's pressure.
! Shift log P, log rho, log T and log R of the envelope by the
! interior-minus-envelope difference there, so that the envelope joins
! the interior continuously, and drop point 1. The convective flag moves
! with the rest (the F77 starin.f left it unshifted, one point off;
! fixed 2026-10).
            pressure_offset = s%logP(s%nz) - s%env_points%env_log10_pressure(1)
            density_offset = s%logRho(s%nz) - s%env_points%env_log10_density(1)
            temperature_offset = s%logT(s%nz) - &
                 s%env_points%env_log10_temperature(1)
            radius_offset = s%logR(s%nz) - s%env_points%env_log10_radius(1)
            do j = 1,s%env_points%num_env_points - 1
               s%env_points%env_log10_density(j) = s%env_points%env_log10_density(j+1)+density_offset
               s%env_points%env_log10_pressure(j) = s%env_points%env_log10_pressure(j+1)+pressure_offset
               s%env_points%env_log10_radius(j) = s%env_points%env_log10_radius(j+1)+radius_offset
               s%env_points%env_log10_mass(j) = s%env_points%env_log10_mass(j+1)
               s%env_points%env_log10_temperature(j) = s%env_points%env_log10_temperature(j+1)+ &
                    temperature_offset
               s%env_points%env_hydrogen_fraction(j) = s%env_points%env_hydrogen_fraction(j+1)
               s%env_points%env_metal_fraction(j) = s%env_points%env_metal_fraction(j+1)
               s%env_points%env_convective_flag(j) = s%env_points%env_convective_flag(j+1)
            end do
            s%env_points%num_env_points = s%env_points%num_env_points - 1
! Append the envelope points inside the new fitting mass as shells
! and interpolate the last one onto it (append_envelope_points, shared
! with rebuild_envelope; here a point exactly at the fitting mass is
! kept as a shell).
            call append_envelope_points(s, .true., sm%num_species, s%xa, &
                 s%logRho, s%luminosity_lsun, s%logP, s%logR, &
                 s%log_mass, s%logT, s%convective_flag, s%nz, ierr)
            if (ierr /= 0) then
! A degenerate interpolation interval in mass is fatal; run_yrec
! returns the error.
               lines = unwritten_line
               write(lines,*) 'STARIN: degenerate', &
                    ' envelope interpolation interval (was STOP 9998)'
               call core_log(lines)
               return
            end if
! The new shells rotate with the old last shell's omega, with the
! specific angular momentum of a thin shell, j = (2/3) omega r^2.
            if (s%physics%rotation_active) then
               do j = sm%old_last_shell+1,s%nz
                  s%omega(j) = s%omega(sm%old_last_shell)
                  s%j_rot(j) = cc23*s%omega(sm%old_last_shell)* &
                       exp10((2.0d0*s%logR(j)))
               end do
            end if
            lines = unwritten_line
            write(lines,910)
            call core_stdout(lines)
 910  format(1x,'NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS'/ &
            ' J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F')
      allocate(point_lines(s%nz-sm%old_last_shell+1))
      point_lines = unwritten_line
      write(point_lines,911)(j,s%logRho(j),s%luminosity_lsun(j),s%logP(j), &
                   s%logR(j),s%log_mass(j)-s%log_total_mass, &
                   s%logT(j),s%convective_flag(j), j = sm%old_last_shell,s%nz)
      call core_stdout(point_lines)
 911  format(i5,1p6e16.8,l2)
end subroutine append_envelope_as_shells


! ---------------------------------------------------------------
! The unlogged mass coordinate m and the mass dm of each shell (g),
! from the (possibly extended or refitted) log-mass grid.
!   star reads: log_mass, log_total_mass, nz
!   star writes: dm, m
subroutine build_shell_masses(s)
      use observables, only: shell_masses_from_log_mass
      type(star_info), intent(inout) :: s
      call shell_masses_from_log_mass(s%log_mass, s%log_total_mass, &
           s%nz, s%m, s%dm)
end subroutine build_shell_masses


! ---------------------------------------------------------------
! Rotation state of the starting model: shape factors, moments of
! inertia, and the totals of angular momentum and rotational kinetic
! energy, echoed beside the totals held on entry.
!   star reads: dm, log_mass, logR, logRho, nz, omega,
!     physics%rotation_active
!   star writes: eta_squared, fp_rot, ft_rot, i_rot, j_rot,
!     kinetic_energy_rot, mean_gravity, mean_radius,
!     total_angular_momentum, total_rotational_ke
subroutine initialize_rotation_state(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(inout) :: ierr
      integer :: i
      real(dp) :: angular_momentum_sum
      real(dp) :: rotational_ke_sum
      real(dp) :: shell_angular_momentum
      character(len=message_length) :: lines(4)
! di_domega of omega_from_j, which this caller does not use
      real(dp) :: qiw(json)
      if (s%physics%rotation_active) then
! Shape factors fp_rot and ft_rot, mean gravity and radius, and
! eta_squared, from omega.
       call rotation_shape_factors(s%logRho,s%logR,s%log_mass,s%nz,s%omega, &
            s%eta_squared,s%fp_rot,s%ft_rot, &
            s%mean_gravity,s%mean_radius,ierr)
       if (ierr /= 0) return
! Moment of inertia of each shell, i_rot (legacy HI).
       call zone_moments_of_inertia(s%eta_squared,s%logR,s%log_mass,s%dm,1,s%nz, &
            s%omega,s%mean_radius,s%i_rot,qiw)
! Per shell: angular momentum J_i = omega_i I_i, specific angular
! momentum j_rot = J_i/dm_i and kinetic energy (1/2) omega_i J_i;
! summed over the star.
       angular_momentum_sum = 0.0d0
       rotational_ke_sum = 0.0d0
       do i = 1,s%nz
          shell_angular_momentum = s%omega(i)*s%i_rot(i)
          s%j_rot(i) = shell_angular_momentum/s%dm(i)
          s%kinetic_energy_rot(i) = 0.5d0*s%omega(i)*shell_angular_momentum
          angular_momentum_sum = angular_momentum_sum+shell_angular_momentum
          rotational_ke_sum = rotational_ke_sum + s%kinetic_energy_rot(i)
       end do
       lines = unwritten_line
       write(lines,560)s%total_angular_momentum, &
            angular_momentum_sum,s%total_rotational_ke,rotational_ke_sum
       call core_log(lines)
 560     format(1x,'TOTAL J OF STAR - PREVIOUS ',1pe21.13,' NEW ', &
              1pe21.13/1x,'TOTAL ROTATIONAL K.E. OF STAR - PREVIOUS ', &
              1pe21.13,' NEW ',1pe21.13)
       s%total_angular_momentum = angular_momentum_sum
       s%total_rotational_ke = rotational_ke_sum
      end if
end subroutine initialize_rotation_state


! ---------------------------------------------------------------
! Surface-mixture bookkeeping. On the first card (run_index <= 1), set
! the surface X, Y, Z and He3, build from them and the controls' metal
! mixture the number fractions envelope_number_fractions of the 12
! mixture elements, Z without C, N and O (envelope_metals_no_cno) and
! envelope_inverse_atomic_weight, and install them (settings_install's
! install_surface_mixture_settings: convection's metals and the EOS
! mixture, from which the EOS builds its SCV envelope table when SCV is
! on). On every card, regenerate the surface opacity tables for the
! envelope X (install_surface_x_settings).
!   star reads: nz, physics%mixture_weights_seed,
!     surface_hydrogen_fraction, surface_metal_fraction, xa
!   star writes: envelope_he3_fraction, envelope_helium_fraction,
!     envelope_hydrogen_fraction, envelope_inverse_atomic_weight,
!     envelope_metal_fraction, envelope_metals_no_cno,
!     envelope_number_fractions
subroutine update_surface_mixture(s, sm, ierr)
      type(star_info), intent(inout) :: s
      type(starting_model_state), intent(inout) :: sm
      integer, intent(inout) :: ierr
      real(dp) :: species_mix_weights(n_mix_species)
      integer :: i
      real(dp) :: mixture_weight_sum
      real(dp) :: mixture_scale_factor
      if (sm%run_index <= 1) then
! The relative mass fractions w_i of the ten metals are the controls'
! mixture_weights_seed.
      do i=1, n_mix_species
         species_mix_weights(i)=s%physics%mixture_weights_seed(i)
      end do

! Surface: X_env = surface_hydrogen_fraction, Z_env = surface_metal_fraction,
! X(He3) of the last shell, Y_env = 1 - X_env - Z_env - X(He3).
      s%envelope_hydrogen_fraction = s%surface_hydrogen_fraction
      s%envelope_metal_fraction = s%surface_metal_fraction
      s%envelope_helium_fraction = 1.0d0 - s%envelope_hydrogen_fraction - &
           s%envelope_metal_fraction - s%xa(i_he3,s%nz)
      s%envelope_he3_fraction = s%xa(i_he3,s%nz)
! S = sum of w_i over the ten metals (all but H and He).
      mixture_weight_sum = species_mix_weights(ix_na)+species_mix_weights(ix_al)+ &
           species_mix_weights(ix_mg)+species_mix_weights(ix_fe)+ &
           species_mix_weights(ix_si)+species_mix_weights(ix_c)+ &
           species_mix_weights(ix_o)+species_mix_weights(ix_n)+ &
           species_mix_weights(ix_ar)+species_mix_weights(ix_ne)
! Z outside C, N and O: envelope_metals_no_cno = Z_env (S - w_C - w_O - w_N)/S.
      s%envelope_metals_no_cno = s%envelope_metal_fraction*(mixture_weight_sum - &
           species_mix_weights(ix_c)-species_mix_weights(ix_o)- &
           species_mix_weights(ix_n))/mixture_weight_sum
! With f = Z_env/S, the metals' mass fractions are X_i = f w_i; H and He
! get w = X_env/f and (1 - X_env - Z_env)/f, so that f w is their mass
! fraction too. Then w_i becomes X_i/A_i, and
!   envelope_inverse_atomic_weight = sum_i X_i/A_i,
! the inverse mean atomic weight 1/mu_ion of the whole surface mixture.
      mixture_scale_factor = s%envelope_metal_fraction/mixture_weight_sum
      species_mix_weights(ix_h) = s%envelope_hydrogen_fraction/ &
           mixture_scale_factor
      species_mix_weights(ix_he) = (1.0d0-s%envelope_hydrogen_fraction- &
           s%envelope_metal_fraction)/mixture_scale_factor
      mixture_weight_sum = 0.0d0
      do i = 1,n_mix_species
       species_mix_weights(i) = mixture_scale_factor*species_mix_weights(i)/ &
            atomic_weight(i)
       mixture_weight_sum = mixture_weight_sum + species_mix_weights(i)
      end do
      s%envelope_inverse_atomic_weight = mixture_weight_sum
      mixture_scale_factor = 1.0d0/s%envelope_inverse_atomic_weight
! Number fractions of the elements, envelope_number_fractions(i) =
! (X_i/A_i)/envelope_inverse_atomic_weight, which the EOS needs (DBG
! 1/96): compute_blended_eos and compute_saha_eos (legacy EQSAHA).
      do i = 1,n_mix_species
       s%envelope_number_fractions(i) = species_mix_weights(i)*mixture_scale_factor
      end do
! Install the new mixture: convection's metals and the EOS mixture.
      call install_surface_mixture_settings(s)
      end if
! Regenerate the surface opacity tables for the surface X (DBG 11/95).
      call install_surface_x_settings(s%envelope_hydrogen_fraction, ierr)
      if (ierr /= 0) return

end subroutine update_surface_mixture


! ---------------------------------------------------------------
! Copy the step-start structure: log P, log T, log R and L into the
! *_start arrays (legacy HPOLD), against which the timestep control
! limits the change from model to model; with them Teff and the surface
! L, and, when rotating, omega into omega_start.
!   star reads: log_Teff, logP, logR, logT, luminosity_lsun, nz, omega,
!     physics%rotation_active
!   star writes: log_Teff_start, logP_start, logR_start, logT_start,
!     luminosity_lsun_start, omega_start, surface_luminosity_lsun_start
subroutine snapshot_step_start(s)
      type(star_info), intent(inout) :: s
      integer :: i
      do i = 1,s%nz
         s%logP_start(i) = s%logP(i)
         s%logT_start(i) = s%logT(i)
         s%logR_start(i) = s%logR(i)
         s%luminosity_lsun_start(i) = s%luminosity_lsun(i)
      end do
      s%log_Teff_start = s%log_Teff
      s%surface_luminosity_lsun_start = s%luminosity_lsun(s%nz)
      if (s%physics%rotation_active) then
         do i = 1,s%nz
          s%omega_start(i) = s%omega(i)
         end do
      end if
end subroutine snapshot_step_start


end module starting_model
