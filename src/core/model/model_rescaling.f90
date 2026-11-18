!----------------------------------------------------------------------
! model_rescaling
!----------------------------------------------------------------------
! Rescaling a model at the start of a rescale run card (rescale_model):
! its X and Z (every model), its total mass and one species
! (main-sequence models) or its core and total mass (horizontal-branch
! models), and optionally a Z ramp toward the center. On a later card
! the envelope mixture follows a new Z. rescale_model is a driver over
! steps that take the star and a rescale_state record; the arithmetic on
! the model's arrays is the star-blind kernel rescaling_kernels.
module model_rescaling
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      use physical_constants, only: ln10
      use shell_boundaries, only: find_shell_boundaries
      use star_info_lib, only: star_info, i_h1
      use control_slots, only: i_target_mass, i_target_x, i_target_z, &
           i_target_core_mass
      use core_messages, only: core_log, core_report, message_length, unwritten_line
      use rescaling_kernels, only: hb_masses, scale_hydrogen, shift_metals, ramp_metals, &
           set_helium_by_difference, hb_region_masses, mass_stretch_factor, &
           stretch_masses_above, move_core_mass, compute_envelope_mixture
      implicit none
      private
      public :: rescale_model

! A model with central X above main_sequence_x_min is rescaled as a
! main-sequence model, else as a horizontal-branch one [mass fraction];
! a species value relative to H is on the scale log N(H) = 12 [dex]; a
! horizontal-branch mass rescale stretches the mass outside the H shell
! when at least hb_old_method_envelope_msun lies outside the core, else
! outside the core [Msun].
      real(dp), parameter :: main_sequence_x_min = 1.0d-10, log_abundance_scale_h = 12.0d0, &
           hb_old_method_envelope_msun = 0.01d0

! What rescale_model's steps pass on to each other:
!   run_index             the card (the star's current_card);
!   shell_begin, shell_end
!                         a horizontal-branch model's hydrogen shell, the
!                         points shell_begin to shell_end
!                         (measure_hb_masses);
!   old                   its masses before the rescale [Msun]
!                         (measure_hb_masses); rescale_hb_total_mass
!                         brings old%envelope to the new total mass, for
!                         rescale_hb_core_mass's message.
      type :: rescale_state
         integer :: run_index = 0
         integer :: shell_begin = 0
         integer :: shell_end = 0
         type(hb_masses) :: old
      end type rescale_state

contains

!----------------------------------------------------------------------
! rescale_model (legacy rscale)
!----------------------------------------------------------------------
! Rescales the model at the start of run card run_index to the targets
! in star%card_rescale_targets(:,run_index), by control_slots' rows: the
! total mass in Msun (i_target_mass), the surface X (i_target_x), the
! surface Z (i_target_z) and the core mass in Msun (i_target_core_mass).
! A target is applied when positive. In order:
! - rescale_surface_x: every shell's X scaled by X_new/X_surf (every
!   model);
! - rescale_interior_z: every shell's Z shifted by Z_new - Z_surf, the
!   CNO isotopes scaled with it (every model);
! - main-sequence models (central X > main_sequence_x_min):
!   rescale_one_species, one species set to a given abundance, and
!   rescale_main_sequence_mass, the mass scaled by M_new/M_old throughout;
! - horizontal-branch models (no central hydrogen, an H-burning shell
!   around a helium core): measure_hb_masses; rescale_hb_total_mass, the
!   total mass changed by stretching the mass outside the H shell;
!   rescale_hb_core_mass, the core mass changed at fixed H-shell and
!   total mass, the envelope taking up the difference; report_hb_masses;
! - set_helium_from_x_and_z: Y = 1 - X - Z - X(He3) in every shell;
! - ramp_interior_z: with use_z_ramp, Z set to card_z_ramp_core_z in the
!   core and ramped to the surface Z (DBG 5/94).
! The model is the star's: xa, log_mass and log_total_mass (log10 m and
! M in g), star_mass (M in Msun), and luminosity_lsun and
! convective_flag, which locate the H shell; run_index is current_card.
! surface_hydrogen_fraction and surface_metal_fraction hold the old
! surface X and Z on entry, and the targets on return when X or Z is
! rescaled. ierr = 1 when a horizontal-branch rescale would leave no
! envelope.
!   star reads: current_card, xa
!   star writes: none
subroutine rescale_model(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      type(rescale_state) :: rs

      ierr = 0
      rs%run_index = s%current_card

      call rescale_surface_x(s, rs)
      call rescale_interior_z(s, rs)
      if (s%xa(i_h1,1) > main_sequence_x_min) then
         call rescale_one_species(s, rs)
         call rescale_main_sequence_mass(s, rs)
      else
         call measure_hb_masses(s, rs)
         call rescale_hb_total_mass(s, rs, ierr)
         if (ierr /= 0) return
         call rescale_hb_core_mass(s, rs, ierr)
         if (ierr /= 0) return
         call report_hb_masses(s, rs)
      end if
      call set_helium_from_x_and_z(s)
      call ramp_interior_z(s, rs)
end subroutine rescale_model

! ---------------------------------------------------------------
! With a positive X target, scale every shell's X by X_new/X_surf,
! X_surf being the old surface X (surface_hydrogen_fraction), capped at
! X_new (scale_hydrogen). Valid for main-sequence and horizontal-branch
! models alike. A requested X above 1 is refused with a warning. The
! envelope X, star%envelope_hydrogen_fraction (legacy XENV), is left
! alone: it sets the initial mean molecular weight for the Yale EOS and
! must not change after the first run (MHP 7/99).
!   star reads: card_rescale_targets, nz
!   star writes: current_card_envelope_x, surface_hydrogen_fraction, xa
subroutine rescale_surface_x(s, rs)
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(in) :: rs
      real(dp) :: new_x
      character(len=message_length) :: lines(4)

      new_x = s%card_rescale_targets(i_target_x,rs%run_index)
      if (new_x <= 0.0d0) return
      if (new_x <= 1.0d0) then
         call scale_hydrogen(s%xa(i_h1,1:s%nz), new_x, s%surface_hydrogen_fraction)
         s%surface_hydrogen_fraction = new_x
         s%current_card_envelope_x = new_x
      else
         lines = unwritten_line
         write(lines,1000) rs%run_index, new_x
         call core_report(lines)
 1000    format(1x,'WARNING: rescale (card ',i0,'): the requested X =',f9.6, &
              ' is above 1; X is not rescaled')
      end if
end subroutine rescale_surface_x

! ---------------------------------------------------------------
! With a positive Z target, shift every shell's Z by dZ = Z_new - Z_surf,
! Z_surf being the old surface Z (surface_metal_fraction), and scale its
! CNO isotopes with it (shift_metals). The target is applied when
! positive, as the run plan's "0 or negative = keep current" says
! (Sources, 1). A requested Z above 1 is refused with a warning.
! Unlike X, the envelope Z (star%envelope_metal_fraction, legacy ZENV)
! follows the new Z, since many routines read it (DBG 4/95), and the EOS
! mixture with it. On a card after the first, the whole envelope mixture
! is computed again for the new Z and installed: the EOS's inverse mean
! atomic weight and number fractions, Y_env, and convection's Z without
! C, N and O (its O16 floor), with the first card's X_env and He3 (X_env
! sets the Yale EOS's reference mixture and stays the first card's, MHP
! 7/99). On the first card, update_surface_mixture computes the whole
! mixture from the rescaled model next (settings_install's
! install_surface_z_settings meanwhile).
!   star reads: card_rescale_targets, envelope_he3_fraction,
!     envelope_hydrogen_fraction, nz, physics%mixture_weights_seed
!   star writes: current_card_envelope_z, envelope_helium_fraction,
!     envelope_inverse_atomic_weight, envelope_metal_fraction,
!     envelope_metals_no_cno, envelope_number_fractions,
!     surface_metal_fraction, xa
subroutine rescale_interior_z(s, rs)
      use settings_install, only: install_surface_z_settings, &
           install_surface_mixture_settings
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(in) :: rs
      real(dp) :: new_z, delta_z
      character(len=message_length) :: lines(4)

      new_z = s%card_rescale_targets(i_target_z,rs%run_index)
      if (new_z <= 0.0d0) return
      if (new_z <= 1.0d0) then
         delta_z = new_z - s%surface_metal_fraction
         call shift_metals(s%xa(:,1:s%nz), delta_z)
         s%surface_metal_fraction = new_z
         s%current_card_envelope_z = new_z
         s%envelope_metal_fraction = s%current_card_envelope_z
         if (rs%run_index > 1) then
            call compute_envelope_mixture(s%envelope_hydrogen_fraction, &
                 s%envelope_metal_fraction, s%envelope_he3_fraction, &
                 s%physics%mixture_weights_seed, s%envelope_helium_fraction, &
                 s%envelope_metals_no_cno, s%envelope_inverse_atomic_weight, &
                 s%envelope_number_fractions)
            call install_surface_mixture_settings(s)
         else
            call install_surface_z_settings(s)
         end if
      else
         lines = unwritten_line
         write(lines,1002) rs%run_index, new_z
         call core_report(lines)
 1002    format(1x,'WARNING: rescale (card ',i0,'): the requested Z =',f9.6, &
              ' is above 1; Z is not rescaled')
      end if
end subroutine rescale_interior_z

! ---------------------------------------------------------------
! With rescale_species_active (legacy LNEWCP), set one species other
! than X, Y and Z, slot new_species_index, to new_species_value in
! every shell (main-sequence models). new_species_value is a mass
! fraction or, with value_relative_to_h (legacy LREL), a log abundance
! relative to the surface X on a scale where X is 12, by mass:
!   X_k = 10^(value - 12) X_surf,
! so 3.0 means 1e-9 times the surface X. The control itself is left as
! given, so every rescale card converts it afresh (Sources, 2). An
! abundance of 1 or more is refused with a warning.
!   star reads: nz, physics%new_species_index,
!     physics%new_species_value, physics%rescale_species_active,
!     physics%value_relative_to_h
!   star writes: xa
subroutine rescale_one_species(s, rs)
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(in) :: rs
      real(dp) :: species_target
      integer :: zone_idx
      character(len=message_length) :: lines(4)

      if (.not. (s%physics%rescale_species_active .and. &
           s%physics%new_species_value >= 0.0d0)) return
      species_target = s%physics%new_species_value
      if (s%physics%value_relative_to_h) species_target = &
           exp(ln10*(s%physics%new_species_value-log_abundance_scale_h))*s%xa(i_h1,s%nz)
      if (species_target < 1.0d0) then
         do zone_idx = 1,s%nz
            s%xa(s%physics%new_species_index,zone_idx) = species_target
         end do
      else
         lines = unwritten_line
         write(lines,1004) rs%run_index, species_target, s%physics%new_species_index, &
              s%physics%new_species_index
         call core_report(lines)
 1004    format(1x,'WARNING: rescale (card ',i0,'): the requested abundance',1pe11.3, &
              ' of species ',i0,' is 1 or more; species ',i0,' is not rescaled')
      end if
end subroutine rescale_one_species

! ---------------------------------------------------------------
! With a positive mass target, scale the total mass and every mass point
! by M_new/M_old (main-sequence models): log M and every log m shift by
! log10(M_new/M_old); nothing else in the model changes.
!   star reads: card_rescale_targets, nz
!   star writes: log_mass, log_total_mass, star_mass
subroutine rescale_main_sequence_mass(s, rs)
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(in) :: rs
      real(dp) :: new_mass_msun, log_mass_shift
      integer :: zone_idx

      new_mass_msun = s%card_rescale_targets(i_target_mass,rs%run_index)
      if (new_mass_msun <= 0.0d0) return
      log_mass_shift = log10(new_mass_msun/s%star_mass)
      s%log_total_mass = s%log_total_mass + log_mass_shift
      s%star_mass = new_mass_msun
      do zone_idx = 1,s%nz
         s%log_mass(zone_idx) = s%log_mass(zone_idx) + log_mass_shift
      end do
end subroutine rescale_main_sequence_mass

! ---------------------------------------------------------------
! When a horizontal-branch model's mass or core mass is to be rescaled,
! find its H-burning shell, points shell_begin..shell_end
! (find_shell_boundaries, gn 10/05), and its old core, shell and
! envelope masses (hb_region_masses), echoed in the run log.
!   star reads: convective_flag, log_mass, log_total_mass,
!     luminosity_lsun, nz, physics%timestep_factors, solar_mass_cgs,
!     star_mass, xa
!   star writes: none
subroutine measure_hb_masses(s, rs)
      use control_slots, only: itime_core_min
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(inout) :: rs
      integer :: core_edge, envelope_edge, shell_mid
      logical :: has_h_shell
      character(len=message_length) :: lines(4)

      if (.not. hb_mass_rescale_requested(s, rs)) return
      call find_shell_boundaries(s%xa(i_h1,:), s%luminosity_lsun, s%convective_flag, s%nz, &
           s%physics%timestep_factors(itime_core_min), core_edge, envelope_edge, &
           rs%shell_begin, rs%shell_end, shell_mid, has_h_shell)
      rs%old = hb_region_masses(s%log_mass(1:s%nz), s%log_total_mass, rs%shell_begin, &
           rs%shell_end, s%solar_mass_cgs)
      lines = unwritten_line
      write(lines,63) rs%old%core, rs%old%shell, rs%old%envelope, s%star_mass
      call core_log(lines)
   63 format(1x,'HB-OLD MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL',f9.6)
end subroutine measure_hb_masses

! ---------------------------------------------------------------
! Whether the card rescales a horizontal-branch model's total mass or its
! core mass.
!   star reads: card_rescale_targets
!   star writes: none
logical function hb_mass_rescale_requested(s, rs)
      type(star_info), intent(in) :: s
      type(rescale_state), intent(in) :: rs

      hb_mass_rescale_requested = &
           s%card_rescale_targets(i_target_mass,rs%run_index) > 0.0d0 .or. &
           s%card_rescale_targets(i_target_core_mass,rs%run_index) > 0.0d0
end function hb_mass_rescale_requested

! ---------------------------------------------------------------
! With a positive mass target, rescale a horizontal-branch model's total
! mass to M_new by stretching the mass above a fixed base m_b
! (stretch_masses_above):
!   m_i -> m_b + s (m_i - m_b),  s = (M_new - m_b)/(M_old - m_b),
! and likewise M. The base is the top of the H shell, m(shell_end), when
! at least hb_old_method_envelope_msun lies outside the core ("old
! method"), else the top of the core, m(shell_begin-1) ("new method").
! ierr = 1 if the envelope mass would not be positive: the new mass
! outside the H shell in the old method; in the new method the old mass
! outside the core, M_old - M_core, not the new one. rs%old%envelope
! becomes the envelope's mass at the new total mass.
!   star reads: card_rescale_targets, nz, solar_mass_cgs
!   star writes: log_mass, log_total_mass, star_mass
subroutine rescale_hb_total_mass(s, rs, ierr)
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(inout) :: rs
      integer, intent(inout) :: ierr
      real(dp) :: new_mass_msun, delta_env_mass, env_mass_new, env_mass_check, &
           mass_scale_factor
      integer :: base

      new_mass_msun = s%card_rescale_targets(i_target_mass,rs%run_index)
      if (new_mass_msun <= 0.0d0) return
      delta_env_mass = new_mass_msun - s%star_mass
      env_mass_new = rs%old%envelope + delta_env_mass

      if ((((pow(10d0, s%log_total_mass))/s%solar_mass_cgs)-rs%old%core) >= &
           hb_old_method_envelope_msun) then
         if (env_mass_new <= 0.0d0) then
            call report_negative_hb_envelope(rs%run_index, rs%old%envelope, env_mass_new, &
                 new_mass_msun, s%star_mass)
            ierr = 1
            return
         end if
         base = rs%shell_end
         mass_scale_factor = mass_stretch_factor(s%log_mass(base), new_mass_msun, &
              s%star_mass, s%solar_mass_cgs)
         call stretch_masses_above(s%log_mass(1:s%nz), s%log_total_mass, base, &
              mass_scale_factor)
         s%star_mass = new_mass_msun
         rs%old%envelope = env_mass_new
      else
         env_mass_check = ((exp10(s%log_total_mass))/s%solar_mass_cgs)-rs%old%core
         if (env_mass_check <= 0.0d0) then
            call report_negative_hb_envelope(rs%run_index, rs%old%envelope, env_mass_new, &
                 new_mass_msun, s%star_mass)
            ierr = 1
            return
         end if
         base = rs%shell_begin - 1
         mass_scale_factor = mass_stretch_factor(s%log_mass(base), new_mass_msun, &
              s%star_mass, s%solar_mass_cgs)
         call stretch_masses_above(s%log_mass(1:s%nz), s%log_total_mass, base, &
              mass_scale_factor)
         s%star_mass = new_mass_msun
         rs%old%envelope = (exp(ln10*s%log_total_mass)- &
              exp(ln10*s%log_mass(rs%shell_end)))/s%solar_mass_cgs
      end if
end subroutine rescale_hb_total_mass

! ---------------------------------------------------------------
! The error of a total-mass rescale that would leave no envelope, with
! the envelope's old and new masses and the new and old total masses
! [Msun]; the caller returns ierr = 1, and the run stops.
subroutine report_negative_hb_envelope(run_index, old_envelope_msun, new_envelope_msun, &
     new_mass_msun, old_mass_msun)
      integer, intent(in) :: run_index
      real(dp), intent(in) :: old_envelope_msun, new_envelope_msun, new_mass_msun, &
           old_mass_msun
      character(len=message_length) :: lines(4)

      lines = unwritten_line
      write(lines,69) run_index, old_envelope_msun, new_envelope_msun, new_mass_msun, &
           old_mass_msun
      call core_report(lines)
   69 format(1x,'ERROR: rescale (card ',i0,'): the new envelope mass outside', &
           ' the H shell would be negative (old',1pe10.2,', new',e10.2, &
           ' Msun; total mass new',e10.2,', old',e10.2,' Msun); run stopped')
end subroutine report_negative_hb_envelope

! ---------------------------------------------------------------
! With a positive core-mass target, rescale a horizontal-branch model's
! core mass to M_core,new at fixed H-shell mass and total mass, moving
! mass between the core and the envelope (move_core_mass). ierr = 1 if
! the envelope inside the fitting point, rs%old%inner_envelope, would not
! keep a positive mass.
!   star reads: card_rescale_targets, nz, solar_mass_cgs
!   star writes: log_mass
subroutine rescale_hb_core_mass(s, rs, ierr)
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(in) :: rs
      integer, intent(inout) :: ierr
      real(dp) :: new_core_msun, delta_core_mass, delta_env_mass, env_mass_new
      character(len=message_length) :: lines(4)

      new_core_msun = s%card_rescale_targets(i_target_core_mass,rs%run_index)
      if (new_core_msun <= 0.0d0) return
      delta_core_mass = new_core_msun - rs%old%core
      delta_env_mass = - delta_core_mass
      env_mass_new = rs%old%inner_envelope + delta_env_mass
      if (env_mass_new <= 0.0d0) then
         lines = unwritten_line
         write(lines,71) rs%run_index, rs%old%envelope, env_mass_new, new_core_msun, &
              rs%old%core
         call core_report(lines)
   71    format(1x,'ERROR: rescale (card ',i0,'): the core rescaling would leave a', &
              ' negative envelope mass (old',1pe10.2,', new',e10.2, &
              ' Msun; core mass new',e10.2,', old',e10.2,' Msun); run stopped')
!        run_yrec returns the error; the command-line wrapper (main) stops.
         ierr = 1
         return
      end if
      call move_core_mass(s%log_mass(1:s%nz), rs%shell_begin, rs%shell_end, &
           new_core_msun, rs%old%core, s%solar_mass_cgs)
end subroutine rescale_hb_core_mass

! ---------------------------------------------------------------
! Echo the rescaled horizontal-branch core, shell, envelope and total
! masses in the run log.
!   star reads: log_mass, log_total_mass, nz, solar_mass_cgs, star_mass
!   star writes: none
subroutine report_hb_masses(s, rs)
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(in) :: rs
      type(hb_masses) :: new
      character(len=message_length) :: lines(4)

      if (.not. hb_mass_rescale_requested(s, rs)) return
      new = hb_region_masses(s%log_mass(1:s%nz), s%log_total_mass, rs%shell_begin, &
           rs%shell_end, s%solar_mass_cgs)
      lines = unwritten_line
      write(lines,101) new%core, new%shell, new%envelope, s%star_mass
      call core_log(lines)
  101 format(1x,'HB-RESCALED MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL ',f9.6)
end subroutine report_hb_masses

! ---------------------------------------------------------------
! Set every shell's Y from its new X, Z and He3: Y = 1 - X - Z - X(He3).
!   star reads: nz
!   star writes: xa
subroutine set_helium_from_x_and_z(s)
      type(star_info), intent(inout) :: s

      call set_helium_by_difference(s%xa(:,1:s%nz))
end subroutine set_helium_from_x_and_z

! ---------------------------------------------------------------
! With use_z_ramp (legacy LZRAMP; DBG 5/94) and the card's
! card_z_ramp_core_z, card_z_ramp_inner_mass_fraction and
! card_z_ramp_outer_mass_fraction all positive, rescale the interior Z as
! a function of the mass fraction q = m/M (ramp_metals): the core value
! for q below the inner fraction, then linear up to the surface Z,
! surface_metal_fraction, at the outer fraction; Z beyond it is left
! alone.
!   star reads: card_z_ramp_core_z, card_z_ramp_inner_mass_fraction,
!     card_z_ramp_outer_mass_fraction, control%use_z_ramp, log_mass, nz,
!     solar_mass_cgs, star_mass, surface_metal_fraction
!   star writes: xa
subroutine ramp_interior_z(s, rs)
      type(star_info), intent(inout) :: s
      type(rescale_state), intent(in) :: rs

      if (s%control%use_z_ramp .and. (s%card_z_ramp_core_z(rs%run_index) > 0d0) .and. &
           (s%card_z_ramp_inner_mass_fraction(rs%run_index) > 0d0) .and. &
           (s%card_z_ramp_outer_mass_fraction(rs%run_index) > 0d0)) then
         call ramp_metals(s%xa(:,1:s%nz), s%log_mass(1:s%nz), &
              s%star_mass*s%solar_mass_cgs, s%card_z_ramp_core_z(rs%run_index), &
              s%card_z_ramp_inner_mass_fraction(rs%run_index), &
              s%card_z_ramp_outer_mass_fraction(rs%run_index), s%surface_metal_fraction)
      end if
end subroutine ramp_interior_z

! ---------------------------------------------------------------
! Sources
! 1. The F77 RSCALE (git show 6cd56739:src/rscale.f) applied the Z target
!    at >= 0, so an unset rescale_z, whose default is 0, set Z = 0;
!    fixed 2026-10. On a later card the envelope mixture kept the first
!    card's metals until 2026-11 (bug_sweep M85, M96).
! 2. The converted species value used to overwrite the control, and a
!    second rescale converted it again, to about 1e-12 X.
end module model_rescaling
