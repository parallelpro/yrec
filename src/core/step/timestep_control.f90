!----------------------------------------------------------------------
! timestep_control
!----------------------------------------------------------------------
! The next model's timestep: compute_timestep takes the most stringent
! of the limits from the structure change, the change of omega,
! hydrogen and helium burning and the HR-diagram motion, within the
! growth limit, the fixed-timestep override and the end age. The limits
! are the star-blind kernels of timestep_limits.f90; this module gathers
! their inputs from the star and stores the result.
module timestep_control
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: compute_timestep

contains

!----------------------------------------------------------------------
! compute_timestep (legacy htimer.f)
!----------------------------------------------------------------------
! The timestep for the next model, in steps:
!   1. a fixed timestep (card_fixed_timestep_active for this run card)
!      replaces all the limits; caveat user;
!   2. otherwise switch_off_structure_limits, then the limits of
!      timestep_limits, each set by the change over the step just taken
!      (dt_unlimited for a criterion that is off or does not apply):
!        structure_dt  the change of log P, log T, log R and L at any
!                      point, while the structure-based limits are on;
!        rotation_dt   the change of omega, rotating models only (it
!                      also sets star%max_domega_frac);
!        hydrogen_dt   hydrogen burning;
!        helium_dt     helium burning (with the triple-alpha rate at the
!                      hottest zone in a flash, flash_triple_alpha_rate);
!        envelope_dt   the motion in the HR diagram against the envelope
!                      triangle, when on and the structure-based limits
!                      are off;
!      and the tightest of them, within the growth limit
!      atime(itime_max_dt_frac) times a previous timestep above 1 s
!      (take_tightest_limit);
!   3. when evolving to a fixed age (card_stops_active with
!      card_max_age_yr > 0), the cut that keeps the model from passing
!      it (clip_to_end_age).
! It reads the model from the star and sets star%dt to the new timestep
! (s) and star%timestep_yr to it in years. A negative star%dt (a stored
! model's) skips steps 1-2, and only its magnitude is kept.
!   star reads: card_fixed_timestep_active, card_fixed_timestep_yr,
!     core_cz_top_index, current_card, h_shell_bot_index,
!     h_shell_mid_index, log_Teff, log_Teff_start, logP, logP_start,
!     logR, logR_start, logRho, logT, logT_start,
!     luminosity_breakdown_lsun, luminosity_lsun, luminosity_lsun_start,
!     m, nz, omega, omega_start, surface_luminosity_lsun_start, xa,
!     xa_start
!   star writes: dt, max_domega_frac, timestep_yr
subroutine compute_timestep(s)
      use timestep_limits, only: timestep_settings, structure_change_limit, &
           omega_change_limit, hydrogen_burning_limit, helium_burning_limit, &
           hr_motion_limit
      use star_info_lib, only: i_lum_pp1, i_lum_pp2, i_lum_pp3, i_lum_cno, &
           i_lum_3alpha
      use physical_constants, only: seconds_per_year
      type(star_info), intent(inout) :: s
      type(timestep_settings) :: set
      real(dp) :: structure_dt, rotation_dt, hydrogen_dt, helium_dt, envelope_dt

      if (s%dt >= 0.0d0) then
         if (s%card_fixed_timestep_active(s%current_card)) then
            s%dt = s%card_fixed_timestep_yr(s%current_card)*seconds_per_year
         else
            call switch_off_structure_limits(s)
            set = timestep_settings_of(s)
            structure_dt = structure_change_limit(set, s%dt, s%logP_start, &
                 s%logT_start, s%logR_start, s%luminosity_lsun_start, s%logP, &
                 s%logT, s%logR, s%luminosity_lsun, s%nz)
            call omega_change_limit(set, s%dt, s%omega_start, s%omega, s%nz, &
                 rotation_dt, s%max_domega_frac)
! The hydrogen-burning luminosity, L_H = L_pp1 + L_pp2 + L_pp3 + L_CNO
! (Lsun).
            hydrogen_dt = hydrogen_burning_limit(set, &
                 s%luminosity_breakdown_lsun(i_lum_pp1) + s%luminosity_breakdown_lsun(i_lum_pp2) + &
                 s%luminosity_breakdown_lsun(i_lum_pp3) + s%luminosity_breakdown_lsun(i_lum_cno), &
                 s%xa, s%xa_start, s%logRho, s%logT, s%m, s%luminosity_lsun, &
                 s%core_cz_top_index, s%h_shell_mid_index, s%nz)
            helium_dt = helium_burning_limit(set, s%luminosity_breakdown_lsun(i_lum_3alpha), &
                 flash_triple_alpha_rate(s), s%xa, s%logRho, s%m, s%luminosity_lsun, &
                 s%core_cz_top_index, s%h_shell_bot_index)
            envelope_dt = hr_motion_limit(set, s%dt, s%log_Teff_start, &
                 s%surface_luminosity_lsun_start, s%log_Teff, s%luminosity_lsun(s%nz))
            call take_tightest_limit()
         end if
      end if
      s%dt = abs(s%dt)
      s%timestep_yr = s%dt/seconds_per_year
      call clip_to_end_age(s)

contains

! ---------------------------------------------------------------
! The timestep: the tightest of the limits, in this order, within the
! growth limit (tightest_timestep). A new limit is one more entry here.
subroutine take_tightest_limit()
      use timestep_limits, only: tightest_timestep
      s%dt = tightest_timestep(set, s%dt, [ &
           hydrogen_dt, &
           helium_dt, &
           rotation_dt, &
           structure_dt, &
           envelope_dt])
end subroutine take_tightest_limit

end subroutine compute_timestep

!----------------------------------------------------------------------
! The structure-based limits are switched off for the rest of the run
! (MHP 9/01) once the centre is hotter than log Tc = 7.1 while the
! gravothermal luminosity L_grav is negative; the switch goes to
! standard output.
!   star reads: logT, luminosity_breakdown_lsun
!   star writes: structure_dt_limits_active
subroutine switch_off_structure_limits(s)
      use star_info_lib, only: i_lum_grav
      use core_messages, only: core_stdout, message_length, unwritten_line
      type(star_info), intent(inout) :: s
! The central log T above which, with L_grav < 0, they go off.
      real(dp), parameter :: structure_limits_log_tc_max = 7.1d0
      character(len=message_length) :: lines(1)
      if (s%structure_dt_limits_active) then
         if (s%logT(1) > structure_limits_log_tc_max .and. &
              s%luminosity_breakdown_lsun(i_lum_grav) < 0.0d0) then
            lines = unwritten_line
            write(lines, 100) s%logT(1), s%luminosity_breakdown_lsun(i_lum_grav)
 100        format(1x,'timestep control: structure-based limits', &
                 ' disabled (log Tc =',f7.4,', L_grav =',es10.2,' Lsun)')
            call core_stdout(lines)
            s%structure_dt_limits_active = .false.
         end if
      end if
end subroutine switch_off_structure_limits

!----------------------------------------------------------------------
! The settings of the limits, from the star's controls and state
! (timestep_limits' timestep_settings).
!   star reads: physics%max_domega_global, physics%rotation_active,
!     physics%timestep_factors, physics%tri_delta_logl,
!     physics%tri_delta_teffl, physics%use_envelope_triangle_dt,
!     solar_luminosity_cgs, solar_mass_cgs, structure_dt_limits_active
!   star writes: none
function timestep_settings_of(s) result(set)
      use timestep_limits, only: timestep_settings
      type(star_info), intent(inout) :: s
      type(timestep_settings) :: set
      set%timestep_factors = s%physics%timestep_factors
      set%max_domega_global = s%physics%max_domega_global
      set%tri_delta_teffl = s%physics%tri_delta_teffl
      set%tri_delta_logl = s%physics%tri_delta_logl
      set%solar_luminosity_cgs = s%solar_luminosity_cgs
      set%solar_mass_cgs = s%solar_mass_cgs
      set%structure_limits_active = s%structure_dt_limits_active
      set%rotation_active = s%physics%rotation_active
      set%use_envelope_triangle_dt = s%physics%use_envelope_triangle_dt
end function timestep_settings_of

!----------------------------------------------------------------------
! The triple-alpha rate (erg/g/s) at the hottest zone, which sets the
! helium flash's limit (compute_zone_energy's triple-alpha term, with
! deuterium zeroed unless the extended composition is used); zero when
! helium does not burn (timestep_limits' min_burning_luminosity) or the
! burning is not a flash.
!   star reads: logRho, logT, luminosity_breakdown_lsun, nz,
!     physics%use_extended_composition, xa
!   star writes: none
real(dp) function flash_triple_alpha_rate(s) result(rate)
      use timestep_limits, only: helium_flash_applies, hottest_zone, &
           min_burning_luminosity
      use burn_step, only: compute_zone_energy
      use net_lib, only: nuclear_energy_result, chain_triple_alpha
      use star_info_lib, only: i_lum_3alpha, i_h2, n_species_extended
      type(star_info), intent(inout) :: s
      real(dp) :: zone_xa(n_species_extended)
      type(nuclear_energy_result) :: burn
      integer :: zone
      rate = 0.0d0
      if (s%luminosity_breakdown_lsun(i_lum_3alpha) <= min_burning_luminosity) return
      if (.not. helium_flash_applies(s%logRho)) return
      zone = hottest_zone(s%logT, s%nz)
      zone_xa = s%xa(:,zone)
      if (.not. s%physics%use_extended_composition) zone_xa(i_h2) = 0.0d0
      call compute_zone_energy(s, s%logRho(zone), s%logT(zone), zone_xa, zone, burn)
      rate = burn%eps_by_chain(chain_triple_alpha)
end function flash_triple_alpha_rate

!----------------------------------------------------------------------
! When evolving to a given age, shorten the timestep so that the model
! does not pass it: dt = min(dt, card_max_age_yr - age), in years. Only
! for a true stop age (card_max_age_yr > 0), since card_stops_active also
! covers other stop conditions (MHP 10/24).
!   star reads: age_gyr, card_max_age_yr, card_stops_active,
!     current_card
!   star writes: dt, timestep_yr
subroutine clip_to_end_age(s)
      use physical_constants, only: seconds_per_year
      type(star_info), intent(inout) :: s
      real(dp) :: time_left_years
      if (s%card_stops_active(s%current_card) .and. s%card_max_age_yr(s%current_card) > 0.0d0) then
         time_left_years = s%card_max_age_yr(s%current_card) - s%age_gyr*1.0d9
         if (time_left_years < s%timestep_yr) then
            s%timestep_yr = time_left_years
            s%dt = s%timestep_yr*seconds_per_year
         end if
      end if
end subroutine clip_to_end_age

end module timestep_control
