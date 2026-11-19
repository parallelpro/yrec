!----------------------------------------------------------------------
! mass_change
!----------------------------------------------------------------------
! The mass-change step of a timestep: mass loss or accretion at the
! surface, when
! mass_change_active is set. The rate (Msun/yr, positive for accretion)
! is mass_accretion_rate, or the Reimers law with use_reimers_wind. (A
! solar wind scaled by the surface rotation, MHP 8/10, is retired:
! step/retired/rotation_scaled_solar_wind.f90.txt.) The step
!   - limits the timestep so that |dM| is at most the fraction
!     accretion_total_fraction of the star, the fraction
!     accretion_cz_fraction of its surface convection zone and half the
!     mass outside the outer zone, and with accretion ends the step at
!     the disk-locking age, where the mass change stops;
!   - with accretion, mixes the accreted composition into the
!     convection zone;
!   - with mass loss from a rotating star, removes the angular momentum
!     the lost matter carries (a thermal wind);
!   - with mass loss, compresses the zone masses of the convection
!     zone; then sets the new total mass and the shell masses dm.
! (The thermal-energy and entropy bookkeeping for accretion, llp
! 06/15/2009, whose shift of the convection zone was zero, is retired:
! step/retired/accretion_entropy.f90.txt; docs/bug_sweep.md M323.)
! The model is read from and written to star; the timestep is star%dt.
! The step's timestep limits, the convection zone's new masses, the
! wind's angular momentum and the accreted composition are computed by
! the star-blind kernels of mass_change_kernel.f90. change_mass and
! apply_mass_change call their steps in order; the steps are module
! procedures with the star and what they need as arguments.
module mass_change
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info, i_h1, i_metals
      use mass_change_kernel, only: wind_angular_momentum, &
           convection_zone_mass, change_shell_masses, &
           wind_specific_angular_momentum, remove_wind_angular_momentum, &
           mix_accreted_composition, mass_change_timestep_limit, &
           clip_step_to_disk_locking_age, disk_has_dispersed
      use physical_constants, only: log10_radius_from_l_teff, seconds_per_year
      use math_lib
      use core_messages, only: core_log, core_report, message_length, unwritten_line
      implicit none
      private
      public :: change_mass

contains

!----------------------------------------------------------------------
! change_mass
!----------------------------------------------------------------------
! Sets the rate and measures the convective envelope; then calls
! apply_mass_change. Outputs:
! old_log_envelope_mass_fraction = log10(m(nz)/M) before the change
! (set only when the step runs), and new_atmosphere_fit_needed, true
! when matter was accreted. With accretion and disk locking, a model at
! or past the disk-locking age switches the mass change off and returns
! without a step.
!   star reads: disk_gate_age_gyr, physics%disk_locking_active,
!     physics%disk_locking_age_gyr, physics%mass_accretion_rate
!   star writes: mass_change_active
subroutine change_mass(s, old_log_envelope_mass_fraction, new_atmosphere_fit_needed, ierr)
      type(star_info), intent(inout) :: s
      real(dp), intent(out) :: old_log_envelope_mass_fraction
      logical, intent(out) :: new_atmosphere_fit_needed
      integer, intent(out) :: ierr
      logical :: disk_exhausted
      real(dp) :: mass_loss_rate_msun_yr, total_radius_cm
      real(dp) :: cz_total_mass_below_fitting

! Start from the constant rate mass_accretion_rate (Msun/yr); the step
! runs only with mass_change_active.
      ierr = 0
      new_atmosphere_fit_needed = .false.
      mass_loss_rate_msun_yr = s%physics%mass_accretion_rate
      if(.not.s%mass_change_active) return
      call set_mass_change_rate(s, mass_loss_rate_msun_yr, total_radius_cm)
! Accretion with disk locking, from a model at or past the disk-locking
! age (a run started from a model saved after it, or the age left at its
! default 0): the disk is gone, so the mass change stops here, before
! anything is changed, and the timestep is left alone (M319).
      if (s%physics%disk_locking_active .and. mass_loss_rate_msun_yr > 0.0d0) then
         if (disk_has_dispersed(s%disk_gate_age_gyr, s%physics%disk_locking_age_gyr)) then
            s%mass_change_active = .false.
            return
         end if
      end if
      call measure_convective_envelope(s, total_radius_cm, &
           cz_total_mass_below_fitting)
      if(mass_loss_rate_msun_yr > 0.0d0) new_atmosphere_fit_needed = .true.
      call apply_mass_change(s, total_radius_cm, mass_loss_rate_msun_yr, &
           cz_total_mass_below_fitting, old_log_envelope_mass_fraction, &
           disk_exhausted, ierr)
! The disk-locking age was reached this step: the mass change stops for
! the rest of the run.
      if(disk_exhausted) s%mass_change_active = .false.

end subroutine change_mass

! ---------------------------------------------------------------
! The star's radius (total_radius_cm, cm), and the rate (rate_msun_yr,
! Msun/yr, positive for accretion): the constant mass_accretion_rate it
! holds on entry, or the Reimers law with use_reimers_wind.
!   star reads: log10_solar_luminosity, log_L_lsun, log_Teff,
!     physics%reimers_scaling_factor, physics%use_reimers_wind,
!     solar_radius_cgs, star_mass
!   star writes: none
subroutine set_mass_change_rate(s, rate_msun_yr, total_radius_cm)
      type(star_info), intent(inout) :: s
      real(dp), intent(inout) :: rate_msun_yr
      real(dp), intent(out) :: total_radius_cm
! R from L = 4 pi R^2 sigma Teff^4.
      total_radius_cm = exp10(log10_radius_from_l_teff(s%log_L_lsun, s%log_Teff, &
           s%log10_solar_luminosity))
! The Reimers law replaces the constant rate:
!   Mdot [Msun/yr] = reimers_scaling_factor (L/Lsun) (R/Rsun)/(M/Msun),
! with reimers_scaling_factor = -4e-13 by default (eta = 1, mass loss).
      if(s%physics%use_reimers_wind)then
         rate_msun_yr = s%physics%reimers_scaling_factor* &
              exp10(s%log_L_lsun)* &
              (total_radius_cm/s%solar_radius_cgs)/s%star_mass
      end if
end subroutine set_mass_change_rate

! ---------------------------------------------------------------
! With a surface convection zone, sets envelope_cz_bottom_index_start to
! the zone's base (and, with the solver diagnostics on, prints the
! turnover timescale in yr and R/Rsun, from total_radius_cm); without
! one, zeroes the turnover timescale and sets it to nz. Then sums the
! convective envelope's mass, cz_mass (g).
!   star reads: dm, envelope_cz_bottom_index, nz, solar_radius_cgs
!   star writes: convective_turnover_timescale,
!     envelope_cz_bottom_index_start
subroutine measure_convective_envelope(s, total_radius_cm, cz_mass)
      use run_messages, only: solver_diagnostics
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: total_radius_cm
      real(dp), intent(out) :: cz_mass
      integer :: nz, jenv
      character(len=message_length) :: lines(1)
      nz = s%nz
      jenv = s%envelope_cz_bottom_index
! The turnover timescale was computed beforehand, for the previous model
! (compute_turnover_timescale; MHP 2/12). (The F77 code printed it
! unlabelled every step; it is now a diagnostic.)
      if(jenv < nz)then
         if (solver_diagnostics()) then
            lines = unwritten_line
            write(lines,'(a,1p,2e12.4)') &
                 ' mass change: turnover time (yr), R/Rsun', &
                 s%convective_turnover_timescale/seconds_per_year, &
                 total_radius_cm/s%solar_radius_cgs
            call core_log(lines)
         end if
         s%envelope_cz_bottom_index_start = jenv
      else
         s%convective_turnover_timescale = 0.0d0
         s%envelope_cz_bottom_index_start = nz
      end if
! The mass of the convective envelope, M_cz = sum of dm(jenv..nz)
! (convection_zone_mass; MHP 8/10), the step's one M_cz. dm(nz) reaches
! the total mass, so this includes the mass outside the fitting point.
      cz_mass = convection_zone_mass(s%dm(jenv:nz))
end subroutine measure_convective_envelope

!----------------------------------------------------------------------
! apply_mass_change
!----------------------------------------------------------------------
! Limits the timestep, then removes or adds the step's mass,
! dM = Mdot dt, at the surface: the mixing for accretion, the angular
! momentum carried off by mass loss, and the
! convection zone's masses. The arithmetic is in the star-blind kernels
! of mass_change_kernel; this driver calls its steps in order (the
! module procedures below, each with the star and what it needs as
! arguments), and they gather the zone from the star, call the kernels
! and store the results. Inputs: total_radius_cm (R), the rate
! mass_loss_rate_msun_yr (Msun/yr, positive for accretion) and
! cz_total_mass_below_fitting (the convective envelope's mass M_cz, sum
! of its shell masses: the step's timestep limit, the wind's angular
! momentum and the compression take it).
! Outputs: old_log_envelope_mass_fraction = log10(m(nz)/M) before the
! change; disk_exhausted, true when the disk-locking age was reached
! (the timestep is clipped to land on it); ierr = 1 when the envelope
! boundary is the surface zone.
!   star reads: dt, log_mass, log_total_mass, nz, solar_mass_cgs
!   star writes: none
subroutine apply_mass_change(s, total_radius_cm, mass_loss_rate_msun_yr, &
     cz_total_mass_below_fitting, old_log_envelope_mass_fraction, &
     disk_exhausted, ierr)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: total_radius_cm
      real(dp), intent(in) :: mass_loss_rate_msun_yr
      real(dp), intent(in) :: cz_total_mass_below_fitting
      real(dp), intent(out) :: old_log_envelope_mass_fraction
      logical, intent(out) :: disk_exhausted
      integer, intent(out) :: ierr
! the magnitude of the rate (g/s) and the step's mass change dM (g,
! negative for loss)
      real(dp) :: mass_loss_rate_cgs, delta_mass_cgs
! the convection zone's shell masses before the change (g), which the
! wind's angular momentum and the accreted composition are measured on
      real(dp), allocatable :: cz_shell_masses_before(:)

      ierr = 0
      disk_exhausted = .false.
      old_log_envelope_mass_fraction = s%log_mass(s%nz) - s%log_total_mass
! The solar wind scaled by the surface rotation (MHP 8/10) is retired:
! step/retired/rotation_scaled_solar_wind.f90.txt.
      mass_loss_rate_cgs = abs(mass_loss_rate_msun_yr)*s%solar_mass_cgs/ &
           seconds_per_year
      call require_zone_below_surface(s, ierr)
      if (ierr /= 0) return
      call limit_mass_change_timestep(s, mass_loss_rate_msun_yr, &
           mass_loss_rate_cgs, cz_total_mass_below_fitting, disk_exhausted)
! The mass change of the step, dM = Mdot dt, in g (negative for loss).
      delta_mass_cgs = mass_loss_rate_msun_yr*s%solar_mass_cgs*s%dt/ &
           seconds_per_year
      call update_masses(s, delta_mass_cgs, cz_shell_masses_before)
      call spin_down_cz_by_wind(s, total_radius_cm, delta_mass_cgs, &
           cz_shell_masses_before, cz_total_mass_below_fitting)
      call report_new_mass(s, delta_mass_cgs)
      call mix_in_accreted_composition(s, delta_mass_cgs, &
           cz_shell_masses_before)
end subroutine apply_mass_change

! ---------------------------------------------------------------
! ierr = 1 when the envelope boundary is the surface zone.
!   star reads: envelope_cz_bottom_index, nz
!   star writes: none
subroutine require_zone_below_surface(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      character(len=message_length) :: lines(1)
      ierr = 0
      if (s%envelope_cz_bottom_index == s%nz) then
         lines = unwritten_line
         lines(1) = ' ERROR: mass change: the surface convection zone''s base is the' &
              //' outermost shell, so there is no zone to take the mass change'
         call core_report(lines)
         ierr = 1
      end if
end subroutine require_zone_below_surface

! ---------------------------------------------------------------
! Limits the timestep for the mass change at the rate rate_cgs = |Mdot|
! (g/s), by the fractions accretion_total_fraction of the star and
! accretion_cz_fraction of the convective envelope (cz_mass, g) and by
! half the mass outside the fitting point (mass_change_timestep_limit),
! with a line in the run log when it shortens the step; a zero rate
! sets no limit. With accretion only (rate_msun_yr > 0) and disk
! locking, a step that would pass the disk-locking age is clipped to end
! there (clip_step_to_disk_locking_age), and disk_exhausted says so: the
! mass change stops after this step.
!   star reads: disk_gate_age_gyr, log_mass, log_total_mass, nz,
!     physics%accretion_cz_fraction, physics%accretion_total_fraction,
!     physics%disk_locking_active, physics%disk_locking_age_gyr,
!     solar_mass_cgs, star_mass
!   star writes: dt
subroutine limit_mass_change_timestep(s, rate_msun_yr, rate_cgs, cz_mass, &
     disk_exhausted)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: rate_msun_yr, rate_cgs, cz_mass
      logical, intent(out) :: disk_exhausted
      real(dp) :: timestep_limit
      character(len=message_length) :: lines(1)

      disk_exhausted = .false.
      if (rate_cgs > 0.0d0) then
         timestep_limit = mass_change_timestep_limit(rate_cgs, &
              s%physics%accretion_total_fraction, s%physics%accretion_cz_fraction, &
              s%star_mass, s%solar_mass_cgs, cz_mass, &
              exp10(s%log_total_mass) - exp10(s%log_mass(s%nz)))
         if (timestep_limit < s%dt) then
            lines = unwritten_line
            write(lines,10) s%dt, timestep_limit
 10         format(' mass change: timestep reduced from',1pe11.3,' s to',e11.3,' s')
            call core_log(lines)
            s%dt = timestep_limit
         end if
      end if
      if (s%physics%disk_locking_active .and. rate_msun_yr > 0.0d0) &
           call clip_step_to_disk_locking_age(s%disk_gate_age_gyr, &
           s%physics%disk_locking_age_gyr, s%dt, disk_exhausted)
end subroutine limit_mass_change_timestep

! ---------------------------------------------------------------
! The new total mass after the step's change delta_mass_cgs = dM (g),
! the masses of the convection zone (change_shell_masses: with mass loss
! the zone-center masses compress about the zone's base; the shell
! masses dm, the last reaching the new total) and
! log_fit_point_mass_fraction; the surface boundary condition is to be
! recomputed. cz_shell_masses_before returns the zone's shell masses
! before the change, for spin_down_cz_by_wind and
! mix_in_accreted_composition.
!   star reads: envelope_cz_bottom_index, nz, solar_mass_cgs
!   star writes: dm, log_fit_point_mass_fraction, log_mass,
!     log_total_mass, m, recompute_envelope_triangle, star_mass
subroutine update_masses(s, delta_mass_cgs, cz_shell_masses_before)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: delta_mass_cgs
      real(dp), allocatable, intent(out) :: cz_shell_masses_before(:)
      integer :: nz, jenv
      real(dp) :: total_mass_grams_old, total_mass_grams_new

      nz = s%nz
      jenv = s%envelope_cz_bottom_index
      cz_shell_masses_before = s%dm(jenv:nz)
! The new total mass, in Msun (star_mass) and as log10 M in g
! (log_total_mass).
      s%star_mass = s%star_mass + delta_mass_cgs/s%solar_mass_cgs
      total_mass_grams_old = exp10(s%log_total_mass)
      total_mass_grams_new = total_mass_grams_old + delta_mass_cgs
      s%log_total_mass = log10(total_mass_grams_new)
      call change_shell_masses(jenv, delta_mass_cgs, total_mass_grams_new, &
           s%m(1:nz), s%log_mass(1:nz), s%dm(1:nz))
! log_fit_point_mass_fraction = log10(m(nz)/M), m(nz) the fitting point's mass.
      s%log_fit_point_mass_fraction = s%log_mass(nz) - s%log_total_mass
! The envelope triangle (the surface boundary condition) is to be
! recomputed.
      s%recompute_envelope_triangle = .true.
end subroutine update_masses

! ---------------------------------------------------------------
! With mass loss (delta_mass_cgs = dM < 0) from a rotating star: the
! wind carries off the surface's specific angular momentum under the
! convection zone's rotation law (wind_specific_angular_momentum, at the
! radius total_radius_cm), and the zone keeps J_cz + dJ on its new
! masses (remove_wind_angular_momentum, from the shell masses before
! the change). A line in the run log: dM, dJ, omega(nz), the new and old
! J_cz, and the convective envelope's mass cz_mass.
!   star reads: dm, envelope_cz_bottom_index, logR, nz, omega,
!     physics%cz_rotation_law_exponent, physics%rotation_active
!   star writes: j_rot
subroutine spin_down_cz_by_wind(s, total_radius_cm, delta_mass_cgs, &
     cz_shell_masses_before, cz_mass)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: total_radius_cm, delta_mass_cgs, &
           cz_shell_masses_before(:), cz_mass
      integer :: nz, jenv
      real(dp) :: wind_j_surface
      type(wind_angular_momentum) :: wind_am
      character(len=message_length) :: lines(1)

      if (.not. (s%physics%rotation_active .and. delta_mass_cgs < 0.0d0)) return
      nz = s%nz
      jenv = s%envelope_cz_bottom_index
      wind_j_surface = wind_specific_angular_momentum(s%physics%cz_rotation_law_exponent, &
           total_radius_cm, s%omega(nz), s%logR(nz), s%j_rot(nz))
      call remove_wind_angular_momentum(wind_j_surface, delta_mass_cgs, &
           cz_shell_masses_before, s%dm(jenv:nz), s%j_rot(jenv:nz), &
           wind_am)
      lines = unwritten_line
      write(lines,22)delta_mass_cgs,wind_am%lost,s%omega(nz), &
           wind_am%zone_after,wind_am%zone_before,cz_mass
 22   format(' mass change, wind: dM, dJ, omega(nz), J_cz after and before, M_cz', &
           1P6E13.4)
      call core_log(lines)
end subroutine spin_down_cz_by_wind

! ---------------------------------------------------------------
! The run log: the new mass and the step's change delta_mass_cgs, in
! Msun.
!   star reads: solar_mass_cgs, star_mass
!   star writes: none
subroutine report_new_mass(s, delta_mass_cgs)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: delta_mass_cgs
      real(dp) :: delta_mass_msun
      character(len=message_length) :: lines(1)
      delta_mass_msun = delta_mass_cgs/s%solar_mass_cgs
      lines = unwritten_line
      write(lines,20)s%star_mass,delta_mass_msun
 20   format(' mass change: new mass',1pe19.10,' Msun, change',e19.10,' Msun')
      call core_log(lines)
end subroutine report_new_mass

! ---------------------------------------------------------------
! With accretion (delta_mass_cgs = dM > 0): the accreted composition
! mixed into the convection zone (mix_accreted_composition, on the
! zone's shell masses before the change, with X_acc =
! accreted_composition), and the accreted mass recorded.
!   star reads: envelope_cz_bottom_index, nz,
!     physics%accreted_composition
!   star writes: accreted_mass, xa
subroutine mix_in_accreted_composition(s, delta_mass_cgs, &
     cz_shell_masses_before)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: delta_mass_cgs, cz_shell_masses_before(:)
      if (delta_mass_cgs > 0.0d0) then
         call mix_accreted_composition(cz_shell_masses_before, &
              delta_mass_cgs, s%physics%accreted_composition, &
              s%xa(:,s%envelope_cz_bottom_index:s%nz))
         s%accreted_mass = delta_mass_cgs
      end if
end subroutine mix_in_accreted_composition

! ---------------------------------------------------------------
! Sources
!   The F77's MASSLOSS and MDOT (wind/massloss, wind/mdot; MHP 10/01, the
!   accretion efficiency MHP 5/02, the disk locking and M_cz MHP 8/10, the
!   accretion entropy llp 06/15/2009, retired in 2026-11: its shift was
!   zero, and its surface lookup and two EOS calls could only stop the
!   step, bug_sweep M323). Fixed in 2026: the Reimers law (the
!   legacy code evaluated reimers_scaling_factor L/(g R) in CGS, a rate in
!   g/s, but stored it as Msun/yr, and mdot converted it to CGS again, a
!   rate about 1e25 too large); a zero rate, which divided by zero in the
!   timestep limit; an accreting model past the disk-locking age, which
!   took a step of negative length (docs/bug_sweep.md M319). The steps
!   became module procedures with arguments in round 3 (2026-11).
end module mass_change
