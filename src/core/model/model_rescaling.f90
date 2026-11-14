!----------------------------------------------------------------------
! model_rescaling
!----------------------------------------------------------------------
! Rescaling a model at the start of a rescale run card (rescale_model):
! its X and Z (every model), its total mass and one species
! (main-sequence models) or its core and total mass (horizontal-branch
! models), and optionally a Z ramp toward the center. On a later card
! the envelope mixture follows a new Z (compute_envelope_mixture, the
! star-free kernel the first card's update_surface_mixture shares).
module model_rescaling
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use shell_boundaries, only: find_shell_boundaries
      use star_info_lib, only: star_info
      use composition_layout, only: n_mix_species, ix_na, ix_al, ix_mg, ix_fe, &
           ix_si, ix_c, ix_h, ix_o, ix_n, ix_ar, ix_ne, ix_he
      use control_slots, only: i_target_mass, i_target_x, i_target_z, &
           i_target_core_mass
      implicit none
      private
      public :: rescale_model, compute_envelope_mixture

! A model with central X above main_sequence_x_min is rescaled as a
! main-sequence model, else as a horizontal-branch one [mass fraction];
! a species value relative to H is on the scale log N(H) = 12 [dex]; a
! horizontal-branch mass rescale stretches the mass outside the H shell
! when at least hb_old_method_envelope_msun lies outside the core, else
! outside the core [Msun]; z_floor keeps a Z ratio finite [mass fraction].
      real(dp), parameter :: main_sequence_x_min = 1.0d-10, log_abundance_scale_h = 12.0d0, &
           hb_old_method_envelope_msun = 0.01d0, z_floor = 1.d-30

! Atomic weights of the surface-mixture elements, in ix_* slot order
! (Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He) [amu].
      real(dp), parameter :: mixture_atomic_weight(n_mix_species) = [23.0d0, 26.99d0, &
           24.32d0, 55.86d0, 28.1d0, 12.015d0, 1.008d0, 16.0d0, 14.01d0, 39.96d0, 20.19d0, 4.004d0]

contains

!----------------------------------------------------------------------
! rescale_model (legacy rscale)
!----------------------------------------------------------------------
! Rescales the model at the start of run card run_index to the targets
! in star%card_rescale_targets(:,run_index), by control_slots' rows: the
! total mass in Msun (i_target_mass), the surface X (i_target_x), the
! surface Z (i_target_z) and the core mass in Msun (i_target_core_mass).
! A target is applied when positive. In order:
! - X: every shell's X scaled by X_new/X_surf (every model);
! - Z: every shell's Z shifted by Z_new - Z_surf, the CNO isotopes
!   scaled with it (every model);
! - main-sequence models (central X > 1e-10): one species set to a
!   given abundance, and the mass scaled by M_new/M_old throughout;
! - horizontal-branch models (no central hydrogen, an H-burning shell
!   around a helium core): the total mass changed by stretching the
!   mass outside the H shell, and the core mass changed at fixed
!   H-shell and total mass, the envelope taking up the difference;
! - Y = 1 - X - Z - X(He3) in every shell;
! - with use_z_ramp, Z set to rsclzc in the core and ramped to the
!   surface Z (DBG 5/94).
! The model is the star's: xa, log_mass and log_total_mass (log10 m and
! M in g), star_mass (M in Msun), and luminosity_lsun and
! convective_flag, which locate the H shell; run_index is current_card.
! surface_hydrogen_fraction and surface_metal_fraction hold the old
! surface X and Z on entry, and the targets on return when X or Z is rescaled.
!   star reads: card_rescale_targets, card_z_ramp_core_z,
!     card_z_ramp_inner_mass_fraction, card_z_ramp_outer_mass_fraction,
!     control%use_z_ramp, convective_flag, current_card,
!     envelope_he3_fraction, envelope_hydrogen_fraction,
!     luminosity_lsun, nz, physics%mixture_weights_seed,
!     physics%new_species_index, physics%new_species_value,
!     physics%rescale_species_active, physics%timestep_factors,
!     physics%value_relative_to_h, solar_mass_cgs
!   star writes: current_card_envelope_x, current_card_envelope_z,
!     envelope_helium_fraction, envelope_inverse_atomic_weight,
!     envelope_metal_fraction, envelope_metals_no_cno,
!     envelope_number_fractions, log_mass, log_total_mass, star_mass,
!     surface_hydrogen_fraction, surface_metal_fraction, xa
subroutine rescale_model(s, ierr)
      use physical_constants, only: ln10
      use star_info_lib, only: i_h1, i_he3, i_he4, i_metals
      use control_slots, only: itime_core_min
      use core_messages, only: core_log, core_report, message_length, unwritten_line
      use settings_install, only: install_surface_z_settings, install_surface_mixture_settings
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      integer :: run_index
! --- locals ---
      real(dp) :: x_rescale_factor
      integer :: zone_idx, species_idx
      real(dp) :: delta_z, z_rescale_factor
      real(dp) :: log_mass_shift
      real(dp) :: core_mass_old, env_mass_old, env_mass_old2, &
           shell_mass_old, delta_env_mass, delta_core_mass, env_mass_new, &
           mass_scale_factor, env_mass_check, shell_mass_prev, &
           mass_shift_grams, unlogged_mass_temp, core_mass_new, &
           shell_mass_new
      real(dp) :: z_ramp_slope, mass_fraction_local, z_ramp_value
      real(dp) :: species_target
      integer :: core_edge, envelope_edge, shell_begin, shell_end, shell_mid
      logical :: has_h_shell

      integer, intent(out) :: ierr

      ierr = 0
      run_index = s%current_card

      call rescale_surface_x
      call rescale_interior_z
      if(s%xa(i_h1,1) > main_sequence_x_min) then
! Main-sequence model: one species, and the mass.
         call rescale_one_species
         call rescale_main_sequence_mass
      else
! Horizontal-branch model: the total mass changes only the points
! outside the H-burning shell; the core mass changes the points inside
! it, and the mass outside the shell shrinks or grows to match.
         call measure_hb_masses
         call rescale_hb_total_mass
         if (ierr /= 0) return
         call rescale_hb_core_mass
         if (ierr /= 0) return
         call report_hb_masses
      end if
      call set_helium_from_x_and_z
      call ramp_interior_z
      return

contains

! ---------------------------------------------------------------
! Scale every shell's X by X_new/X_surf, X_surf being the old surface X
! (surface_hydrogen_fraction), capping each at X_new:
!   X_i -> min(X_new, X_i X_new/max(X_surf, 1e-20)).
! Valid for main-sequence and horizontal-branch models alike. A
! requested X above 1 is refused with a message.
subroutine rescale_surface_x
      character(len=message_length) :: lines(4)
      if(s%card_rescale_targets(i_target_x,run_index) > 0.0d0) then
         if(s%card_rescale_targets(i_target_x,run_index) <= 1.0d0) then
            x_rescale_factor = s%card_rescale_targets(i_target_x,run_index)/max(s%surface_hydrogen_fraction,1.0d-20)
            do zone_idx = 1,s%nz
               s%xa(i_h1,zone_idx) = min(s%card_rescale_targets(i_target_x,run_index), &
                    s%xa(i_h1,zone_idx)*x_rescale_factor)
            end do
            s%surface_hydrogen_fraction = s%card_rescale_targets(i_target_x,run_index)
            s%current_card_envelope_x = s%card_rescale_targets(i_target_x,run_index)
! The envelope X, star%envelope_hydrogen_fraction (legacy XENV), is
! deliberately left alone: it sets the initial mean molecular weight
! for the Yale EOS and must not change after the first run (MHP 7/99).
         else
! Requested X above 1: X is not changed.
            lines = unwritten_line
            write(lines,1000)run_index,s%card_rescale_targets(i_target_x,run_index)
            call core_report(lines)
 1000       format(1x,'WARNING: rescale (card ',i0,'): the requested X =',f9.6, &
           ' is above 1; X is not rescaled')
         end if
      end if
end subroutine rescale_surface_x

! ---------------------------------------------------------------
! Shift every shell's Z by dZ = Z_new - Z_surf, Z_surf being the old
! surface Z (surface_metal_fraction), and scale the CNO isotopes (slots
! C12..O18) by each shell's new over old Z:
!   Z_i -> max(0, Z_i + dZ),
!   X_k,i -> X_k,i max(0, Z_i + dZ)/(Z_i + 1e-30),  k = C12..O18.
! He3 and the light elements are left alone. The envelope Z follows, and
! on a card after the first the whole envelope mixture with it (below).
! The target is applied when
! positive, as the run plan's "0 or negative = keep current" says (the
! F77 rscale.f applied it at >= 0, so an unset rescale_z, whose default
! is 0, set Z = 0; fixed 2026-10). A requested Z above 1 is refused
! with a message.
subroutine rescale_interior_z
      character(len=message_length) :: lines(4)
      if(s%card_rescale_targets(i_target_z,run_index) > 0.0d0) then
         if(s%card_rescale_targets(i_target_z,run_index) <= 1.0d0) then
            delta_z = s%card_rescale_targets(i_target_z,run_index) - s%surface_metal_fraction
            do zone_idx = 1,s%nz
               z_rescale_factor=max(0.d0,s%xa(i_metals,zone_idx)+delta_z)/ &
                    (s%xa(i_metals,zone_idx)+z_floor)
               s%xa(i_metals,zone_idx) = max(0.0d0,s%xa(i_metals,zone_idx)+delta_z)
               do species_idx = 5,11
                  s%xa(species_idx,zone_idx) = s%xa(species_idx,zone_idx)*z_rescale_factor
               end do
            end do
            s%surface_metal_fraction = s%card_rescale_targets(i_target_z,run_index)
            s%current_card_envelope_z = s%card_rescale_targets(i_target_z,run_index)
! Unlike X, the envelope Z (star%envelope_metal_fraction, legacy ZENV)
! follows the new Z, since many routines read it (DBG 4/95), and the
! EOS mixture with it. On a card after the first, the whole envelope
! mixture is computed again for the new Z and installed: the EOS's
! inverse mean atomic weight and number fractions, Y_env, and
! convection's Z without C, N and O (its O16 floor), with the first
! card's X_env and He3 (X_env sets the Yale EOS's reference mixture and
! stays the first card's, MHP 7/99); they kept the first card's metals
! (M85, M96). On the first card, update_surface_mixture computes the
! whole mixture from the rescaled model next (settings_install's
! install_surface_z_settings meanwhile).
            s%envelope_metal_fraction = s%current_card_envelope_z
            if (run_index > 1) then
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
! Requested Z above 1: Z is not changed.
            lines = unwritten_line
            write(lines,1002)run_index,s%card_rescale_targets(i_target_z,run_index)
            call core_report(lines)
 1002       format(1x,'WARNING: rescale (card ',i0,'): the requested Z =',f9.6, &
           ' is above 1; Z is not rescaled')
         end if
      end if
end subroutine rescale_interior_z

! ---------------------------------------------------------------
! With rescale_species_active (legacy LNEWCP), set one species other
! than X, Y and Z, slot new_species_index, to new_species_value in
! every shell (main-sequence models).
subroutine rescale_one_species
      character(len=message_length) :: lines(4)
         if(s%physics%rescale_species_active.and.s%physics%new_species_value >= 0.0d0) then
! new_species_value is a mass fraction or, with value_relative_to_h
! (legacy LREL), a log abundance relative to the surface X on a scale
! where X is 12, by mass:
!   X_k = 10^(value - 12) X_surf,
! so 3.0 means 1e-9 times the surface X. The control itself is left
! as given, so every rescale card converts it afresh. (The converted
! value used to overwrite it, and a second rescale converted it again,
! to about 1e-12 X.)
            species_target = s%physics%new_species_value
            if(s%physics%value_relative_to_h) species_target = &
                 exp(ln10*(s%physics%new_species_value-log_abundance_scale_h))*s%xa(i_h1,s%nz)
            if(species_target < 1.0d0) then
               do zone_idx = 1,s%nz
                  s%xa(s%physics%new_species_index,zone_idx) = species_target
               end do
            else
! An abundance of 1 or more is refused; the species is not changed.
               lines = unwritten_line
               write(lines,1004)run_index,species_target,s%physics%new_species_index, &
                    s%physics%new_species_index
               call core_report(lines)
 1004       format(1x,'WARNING: rescale (card ',i0,'): the requested abundance',1pe11.3, &
           ' of species ',i0,' is 1 or more; species ',i0,' is not rescaled')
            end if
         end if
end subroutine rescale_one_species

! ---------------------------------------------------------------
! Scale the total mass and every mass point by M_new/M_old
! (main-sequence models): log M and every log m shift by
! log10(M_new/M_old); nothing else in the model changes.
subroutine rescale_main_sequence_mass
         if(s%card_rescale_targets(i_target_mass,run_index) > 0.0d0) then
            log_mass_shift = log10(s%card_rescale_targets(i_target_mass,run_index)/s%star_mass)
            s%log_total_mass = s%log_total_mass + log_mass_shift
            s%star_mass = s%card_rescale_targets(i_target_mass,run_index)
            do zone_idx = 1,s%nz
               s%log_mass(zone_idx) = s%log_mass(zone_idx) + log_mass_shift
            end do
         end if
end subroutine rescale_main_sequence_mass

! ---------------------------------------------------------------
! When a horizontal-branch model's mass or core mass is to be rescaled,
! find its H-burning shell, points shell_begin..shell_end
! (locate_shell_boundaries, gn 10/05), and its old masses in Msun:
!   core           m(shell_begin-1),
!   shell          m(shell_end) - m(shell_begin-1),
!   envelope       M - m(shell_end),
!   env_mass_old2  m(nz) - m(shell_end), the envelope inside the
!                  fitting point.
subroutine measure_hb_masses
      character(len=message_length) :: lines(4)
         if(s%card_rescale_targets(i_target_mass,run_index) > 0.0d0.or.s%card_rescale_targets(i_target_core_mass,run_index) > 0.0d0) then
            call find_shell_boundaries(s%xa(i_h1,:), s%luminosity_lsun, s%convective_flag, s%nz, &
                 s%physics%timestep_factors(itime_core_min), core_edge, envelope_edge, &
                 shell_begin, shell_end, shell_mid, has_h_shell)
            core_mass_old = exp(ln10*s%log_mass(shell_begin-1))/s%solar_mass_cgs
            env_mass_old = (exp(ln10*s%log_total_mass)-exp(ln10*s%log_mass(shell_end)))/s%solar_mass_cgs
            env_mass_old2 = (exp(ln10*s%log_mass(s%nz))-exp(ln10*s%log_mass(shell_end)))/s%solar_mass_cgs
            shell_mass_old = (exp(ln10*s%log_mass(shell_end))-exp(ln10*s%log_mass(shell_begin-1)))/s%solar_mass_cgs
            delta_env_mass = 0.0d0
            delta_core_mass = 0.0d0
            lines = unwritten_line
            write(lines,63)core_mass_old,shell_mass_old,env_mass_old,s%star_mass
            call core_log(lines)
   63       format(1x,'HB-OLD MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL',f9.6)
         end if
end subroutine measure_hb_masses

! ---------------------------------------------------------------
! Rescale a horizontal-branch model's total mass to M_new by stretching
! the mass above a fixed base m_b:
!   m_i -> m_b + s (m_i - m_b),  s = (M_new - m_b)/(M_old - m_b),
! and likewise M. The base is the top of the H shell, m(shell_end),
! when at least 0.01 Msun lies outside the core ("old method"), else
! the top of the core, m(shell_begin-1) ("new method"). Sets ierr if
! the envelope mass would not be positive.
subroutine rescale_hb_total_mass
      character(len=message_length) :: lines(4)
         if(s%card_rescale_targets(i_target_mass,run_index) > 0.0d0) then

            delta_env_mass = s%card_rescale_targets(i_target_mass,run_index) - s%star_mass
            env_mass_new = env_mass_old + delta_env_mass

            if((((pow(10d0, s%log_total_mass))/s%solar_mass_cgs)-core_mass_old) >= hb_old_method_envelope_msun) then
! The new envelope mass outside the H shell must be positive.
            if(env_mass_new <= 0.0d0)then
               lines = unwritten_line
               write(lines,69)run_index,env_mass_old,env_mass_old+delta_env_mass, &
                    s%card_rescale_targets(i_target_mass,run_index),s%star_mass
               call core_report(lines)
   69          format(1x,'ERROR: rescale (card ',i0,'): the new envelope mass outside', &
           ' the H shell would be negative (old',1pe10.2,', new',e10.2, &
           ' Msun; total mass new',e10.2,', old',e10.2,' Msun); run stopped')
               ! run_yrec returns the error; the command-line wrapper
               ! (main) stops.
               ierr = 1
               return
            end if

! Base m_b = m(shell_end); s from the masses in Msun.
        mass_scale_factor =(s%card_rescale_targets(i_target_mass,run_index)-exp(ln10*s%log_mass(shell_end))/s%solar_mass_cgs)/ &
                          (s%star_mass-exp(ln10*s%log_mass(shell_end))/s%solar_mass_cgs)

        s%log_total_mass=log10(exp10(s%log_mass(shell_end))+mass_scale_factor*(exp10(s%log_total_mass)-exp10(s%log_mass(shell_end))))

            s%star_mass = s%card_rescale_targets(i_target_mass,run_index)
            do zone_idx = shell_end+1,s%nz
        s%log_mass(zone_idx)=log10(exp10(s%log_mass(shell_end))+ &
             mass_scale_factor*(exp10(s%log_mass(zone_idx))-exp10(s%log_mass(shell_end))))
            end do
            env_mass_old = env_mass_new

         else

! The mass outside the core must be positive. (This tests the old
! mass, M_old - M_core, not the new one.)
            env_mass_check = ((exp10(s%log_total_mass))/s%solar_mass_cgs)-core_mass_old
            if(env_mass_check <= 0.0d0)then
               lines = unwritten_line
               write(lines,69)run_index,env_mass_old,env_mass_old+delta_env_mass, &
                    s%card_rescale_targets(i_target_mass,run_index),s%star_mass
               call core_report(lines)
               ! run_yrec returns the error; the command-line wrapper
               ! (main) stops.
               ierr = 1
               return
            end if

! Base m_b = m(shell_begin-1); s from the masses in Msun.
        mass_scale_factor =(s%card_rescale_targets(i_target_mass,run_index)-exp(ln10*s%log_mass(shell_begin-1))/s%solar_mass_cgs)/ &
                          (s%star_mass-exp(ln10*s%log_mass(shell_begin-1))/s%solar_mass_cgs)
       s%log_total_mass=log10(exp10(s%log_mass(shell_begin-1))+ &
            mass_scale_factor*(exp10(s%log_total_mass)-exp10(s%log_mass(shell_begin-1))))

            s%star_mass = s%card_rescale_targets(i_target_mass,run_index)
            do zone_idx = shell_begin,s%nz
      s%log_mass(zone_idx)=log10(exp10(s%log_mass(shell_begin-1))+ &
           mass_scale_factor*(exp10(s%log_mass(zone_idx))-exp10(s%log_mass(shell_begin-1))))
            end do
            env_mass_old = (exp(ln10*s%log_total_mass)-exp(ln10*s%log_mass(shell_end)))/s%solar_mass_cgs
         end if
         end if
end subroutine rescale_hb_total_mass

! ---------------------------------------------------------------
! Rescale a horizontal-branch model's core mass to M_core,new at fixed
! H-shell mass and total mass, moving mass between the core and the
! envelope:
! - the core's log m shift by log10(M_core,new/M_core,old);
! - the H-shell points shift by dM_core = M_core,new - M_core,old in m;
! - the envelope points above the shell move toward m(nz), which stays
!   fixed (as does M):
!     m_i -> m_nz - s (m_nz - m_i),
!     s = (m_nz - m_new(shell_end))/(m_nz - m_old(shell_end)).
! Sets ierr if the envelope inside the fitting point, env_mass_old2,
! would not keep a positive mass.
subroutine rescale_hb_core_mass
      character(len=message_length) :: lines(4)
         if(s%card_rescale_targets(i_target_core_mass,run_index) > 0.0d0) then

            shell_mass_prev= (exp(ln10*s%log_mass(shell_end)))
            delta_core_mass = s%card_rescale_targets(i_target_core_mass,run_index) - core_mass_old
            delta_env_mass = - delta_core_mass
            env_mass_new = env_mass_old2 + delta_env_mass
            if(env_mass_new <= 0.0d0)then
               lines = unwritten_line
               write(lines,71)run_index,env_mass_old,env_mass_new, &
                    s%card_rescale_targets(i_target_core_mass,run_index),core_mass_old
               call core_report(lines)
   71          format(1x,'ERROR: rescale (card ',i0,'): the core rescaling would leave a', &
           ' negative envelope mass (old',1pe10.2,', new',e10.2, &
           ' Msun; core mass new',e10.2,', old',e10.2,' Msun); run stopped')
               ! run_yrec returns the error; the command-line wrapper
               ! (main) stops.
               ierr = 1
               return
            end if
            log_mass_shift = log10(s%card_rescale_targets(i_target_core_mass,run_index)/core_mass_old)
            do zone_idx = 1,shell_begin-1
               s%log_mass(zone_idx) = s%log_mass(zone_idx) + log_mass_shift
            end do

! Every H-shell point moves by dM_core (g), which keeps the shell's
! mass.
            mass_shift_grams = delta_core_mass*s%solar_mass_cgs
            do zone_idx = shell_begin,shell_end
               unlogged_mass_temp = exp(ln10*s%log_mass(zone_idx))
               s%log_mass(zone_idx) = log10(unlogged_mass_temp + mass_shift_grams)
            end do
! Shrink or stretch the envelope toward m(nz) to keep the total mass.
        mass_scale_factor =(exp(ln10*s%log_mass(s%nz))-exp(ln10*s%log_mass(shell_end)))/ &
                          (exp(ln10*s%log_mass(s%nz))-shell_mass_prev)

            do zone_idx = shell_end+1,s%nz
        s%log_mass(zone_idx)=log10(exp10(s%log_mass(s%nz))- &
             mass_scale_factor*(exp10(s%log_mass(s%nz))-exp10(s%log_mass(zone_idx))))
            end do

         end if
end subroutine rescale_hb_core_mass

! ---------------------------------------------------------------
! Echo the rescaled horizontal-branch core, shell, envelope and total
! masses.
subroutine report_hb_masses
      character(len=message_length) :: lines(4)
         if(s%card_rescale_targets(i_target_mass,run_index) > 0.0d0.or.s%card_rescale_targets(i_target_core_mass,run_index) > 0.0d0)then
            core_mass_new = exp(ln10*s%log_mass(shell_begin-1))/s%solar_mass_cgs
            env_mass_new = (exp(ln10*s%log_total_mass)-exp(ln10*s%log_mass(shell_end)))/s%solar_mass_cgs
            shell_mass_new = (exp(ln10*s%log_mass(shell_end))-exp(ln10*s%log_mass(shell_begin-1)))/s%solar_mass_cgs
            lines = unwritten_line
            write(lines,101)core_mass_new,shell_mass_new,env_mass_new,s%star_mass
            call core_log(lines)
  101       format(1x,'HB-RESCALED MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL ',f9.6)

         end if
end subroutine report_hb_masses

! ---------------------------------------------------------------
! Set every shell's Y from its new X, Z and He3: Y = 1 - X - Z - X(He3).
subroutine set_helium_from_x_and_z
      do zone_idx = 1,s%nz
         s%xa(i_he4,zone_idx) = 1.0d0-s%xa(i_h1,zone_idx)-s%xa(i_metals,zone_idx)-s%xa(i_he3,zone_idx)
      end do
end subroutine set_helium_from_x_and_z

! ---------------------------------------------------------------
! With use_z_ramp (legacy LZRAMP; DBG 5/94) and rsclzc, rsclzm1 and
! rsclzm2 all positive, rescale the interior Z as a function of the
! mass fraction q = m/M: Z = rsclzc for q < rsclzm1, then linear up to
! the surface Z, Z_surf = surface_metal_fraction, at q = rsclzm2,
!   Z(q) = rsclzc + (q - rsclzm1)(Z_surf - rsclzc)/(rsclzm2 - rsclzm1);
! Z beyond rsclzm2 is left alone. X takes up the change, Y being kept:
! X = 1 - Z - X(He3) - Y. The CNO isotopes (C12..O18) are multiplied by
! Z_new/(Z_old + 1e-30), as in rescale_interior_z. (The factor was
! (Z_old - Z_new)/Z_old, the fraction of Z removed: right only for a
! halved Z, zero for an unchanged one, negative for a raised one.)
subroutine ramp_interior_z
      if(s%control%use_z_ramp.and.(s%card_z_ramp_core_z(run_index) > 0d0).and.(s%card_z_ramp_inner_mass_fraction(run_index) > 0d0) &
           .and.(s%card_z_ramp_outer_mass_fraction(run_index) > 0d0))then
            z_ramp_slope = (s%surface_metal_fraction - s%card_z_ramp_core_z(run_index))/(s%card_z_ramp_outer_mass_fraction(run_index)-s%card_z_ramp_inner_mass_fraction(run_index))
            do zone_idx = 1,s%nz
               mass_fraction_local = exp10(s%log_mass(zone_idx))/(s%star_mass*s%solar_mass_cgs)
               if (mass_fraction_local < s%card_z_ramp_inner_mass_fraction(run_index)) then
                  z_rescale_factor = s%card_z_ramp_core_z(run_index)/(s%xa(i_metals,zone_idx)+z_floor)
                  s%xa(i_metals,zone_idx) = s%card_z_ramp_core_z(run_index)
                  s%xa(i_h1,zone_idx) = 1.0d0-s%xa(i_metals,zone_idx)-s%xa(i_he3,zone_idx)-s%xa(i_he4,zone_idx)
                  do species_idx = 5,11
                      s%xa(species_idx,zone_idx) = s%xa(species_idx,zone_idx)*z_rescale_factor
                  end do
               else if (mass_fraction_local < s%card_z_ramp_outer_mass_fraction(run_index)) then
                  z_ramp_value = (mass_fraction_local-s%card_z_ramp_inner_mass_fraction(run_index))*z_ramp_slope+s%card_z_ramp_core_z(run_index)
                  z_rescale_factor = z_ramp_value/(s%xa(i_metals,zone_idx)+z_floor)
                  s%xa(i_metals,zone_idx) = z_ramp_value
                  s%xa(i_h1,zone_idx) = 1.0d0-s%xa(i_metals,zone_idx)-s%xa(i_he3,zone_idx)-s%xa(i_he4,zone_idx)
                  do species_idx = 5,11
                      s%xa(species_idx,zone_idx) = s%xa(species_idx,zone_idx)*z_rescale_factor
                  end do
               end if
            end do
      end if
end subroutine ramp_interior_z

end subroutine rescale_model

! ---------------------------------------------------------------
! The envelope mixture from the envelope's X, Z and He3 mass fractions
! (x_env, z_env, he3_env) and the controls' relative mass fractions w_i
! of the ten metals (seed_weights, mixture_weights_seed, n_mix_species
! of them in ix_* slot order; H's and He's slots are not read), into
! number_fractions (n_mix_species):
!   y_env = 1 - x_env - z_env - he3_env;
!   with S the sum of w_i over the ten metals, Z outside C, N and O,
!     metals_no_cno = z_env (S - w_C - w_O - w_N)/S;
!   with f = z_env/S, the metals' mass fractions X_i = f w_i, H's x_env
!   and He's 1 - x_env - z_env (He3 with He4); inverse_atomic_weight is
!   sum_i X_i/A_i, 1/mu_ion of the mixture, and number fractions are
!   the (X_i/A_i)/inverse_atomic_weight (the EOS's, compute_blended_eos
!   and compute_saha_eos; DBG 1/96).
! The arithmetic is update_surface_mixture's (AMUENV), in its order.
pure subroutine compute_envelope_mixture(x_env, z_env, he3_env, seed_weights, &
     y_env, metals_no_cno, inverse_atomic_weight, number_fractions)
      real(dp), intent(in) :: x_env, z_env, he3_env, seed_weights(:)
      real(dp), intent(out) :: y_env, metals_no_cno, inverse_atomic_weight, &
           number_fractions(:)
      real(dp) :: weights(n_mix_species), weight_sum, scale_factor
      integer :: i

      do i = 1, n_mix_species
         weights(i) = seed_weights(i)
      end do
      y_env = 1.0d0 - x_env - z_env - he3_env
      weight_sum = weights(ix_na)+weights(ix_al)+ &
           weights(ix_mg)+weights(ix_fe)+ &
           weights(ix_si)+weights(ix_c)+ &
           weights(ix_o)+weights(ix_n)+ &
           weights(ix_ar)+weights(ix_ne)
      metals_no_cno = z_env*(weight_sum - &
           weights(ix_c)-weights(ix_o)- &
           weights(ix_n))/weight_sum
      scale_factor = z_env/weight_sum
      weights(ix_h) = x_env/scale_factor
      weights(ix_he) = (1.0d0-x_env-z_env)/scale_factor
      weight_sum = 0.0d0
      do i = 1, n_mix_species
         weights(i) = scale_factor*weights(i)/mixture_atomic_weight(i)
         weight_sum = weight_sum + weights(i)
      end do
      inverse_atomic_weight = weight_sum
      scale_factor = 1.0d0/inverse_atomic_weight
      do i = 1, n_mix_species
         number_fractions(i) = weights(i)*scale_factor
      end do
end subroutine compute_envelope_mixture

end module model_rescaling
