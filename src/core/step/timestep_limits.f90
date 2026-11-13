!----------------------------------------------------------------------
! timestep_limits
!----------------------------------------------------------------------
! The limits on the next model's timestep, each a star-blind kernel: the
! model's arrays and the settings record (timestep_settings) in, one
! limit in seconds out, dt_unlimited when its criterion is off or does
! not apply. compute_timestep (timestep_control.f90) gathers their
! inputs from the star, calls each, and takes the tightest
! (tightest_timestep). The limits, each set by the change over the step
! just taken:
!   structure_change_limit  the change of log P, log T, log R and L at
!                           any point (legacy ptime)
!   omega_change_limit      the change of omega (legacy wtime)
!   hydrogen_burning_limit  hydrogen burning (legacy xtime)
!   helium_burning_limit    helium burning (legacy ytime)
!   hr_motion_limit         the motion in the HR diagram against the
!                           envelope triangle (legacy entime; JVS 04/14)
! A new limit is one more kernel here, its call and one entry in the
! list compute_timestep's take_tightest_limit hands tightest_timestep,
! and a check in core/test/test_timestep_limits.f90.
module timestep_limits
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use control_slots, only: itime_core_min, itime_dx_core_tot, &
           itime_dx_core_frac, itime_dy_core_tot, itime_dy_core_frac, &
           itime_dx_total, itime_dx_shell, itime_dt, itime_dp, itime_dr, &
           itime_dl, itime_dy_shell, itime_max_dt_frac, itime_dy_total
      use composition_layout, only: i_h1, i_he4, i_metals
      implicit none
      private
      public :: timestep_settings, structure_change_limit, omega_change_limit, &
           hydrogen_burning_limit, helium_burning_limit, hr_motion_limit, &
           tightest_timestep, hottest_zone, helium_flash_applies

! "No limit from this criterion" (s).
      real(dp), parameter, public :: dt_unlimited = 1.0d20
! A limit below this (s) counts as no limit: as from a zero timestep
! stored in the model, it would stall the run (MHP 6/90).
      real(dp), parameter, public :: dt_floor = 1.0d0
! The number of timestep factors, physics%timestep_factors (legacy
! atime): control_slots' itime_* slots, of which itime_dy_total is last.
      integer, parameter, public :: num_timestep_factors = itime_dy_total
! A burning luminosity (Lsun) at or below this sets no limit.
      real(dp), parameter, public :: min_burning_luminosity = 1.0d-34
! The energy released per gram of hydrogen and of helium burned (erg/g).
      real(dp), parameter :: h_burn_energy_per_gram = 6.00d18
      real(dp), parameter :: he_burn_energy_per_gram = 5.85d17
! The helium flash: a central log rho at or above this is taken as a
! giant in the flash rather than a horizontal-branch star; its limit is
! the time for the triple-alpha rate at the temperature maximum to
! release flash_energy_per_gram (about 1.7e-3 in Y), with the rate
! floored at min_triple_alpha_rate (erg/g/s).
      real(dp), parameter :: flash_log_central_density = 5.0d0
      real(dp), parameter :: flash_energy_per_gram = 1.0d15
      real(dp), parameter :: min_triple_alpha_rate = 1.0d-22
! ln 10 as the legacy ptime wrote it (four digits; physical_constants'
! ln10 differs by 7e-6 relative, which would move the structure limit:
! a B item).
      real(dp), parameter :: ln10_ptime = 2.3026d0

! The settings of the limits: the timestep factors (legacy atime,
! physics%timestep_factors); the largest fractional change of omega per
! step (physics%max_domega_global, legacy DTWIND); the envelope
! triangle's size (physics%tri_delta_teffl, tri_delta_logl); the solar
! L and M (cgs); and which limits apply: the structure-based ones
! (star%structure_dt_limits_active, after compute_timestep's switch-off),
! the rotation one (physics%rotation_active) and the HR-motion one
! (physics%use_envelope_triangle_dt, applied when the structure-based
! limits are off).
      type, public :: timestep_settings
            real(dp) :: timestep_factors(num_timestep_factors) = 0.0_dp
            real(dp) :: max_domega_global = 0.0_dp
            real(dp) :: tri_delta_teffl = 0.0_dp, tri_delta_logl = 0.0_dp
            real(dp) :: solar_luminosity_cgs = 0.0_dp, solar_mass_cgs = 0.0_dp
            logical :: structure_limits_active = .false.
            logical :: rotation_active = .false.
            logical :: use_envelope_triangle_dt = .false.
      end type timestep_settings

contains

!----------------------------------------------------------------------
! The factor f by which the timestep shrinks, clipped so that the
! timestep changes by at most a = timestep_factors(itime_max_dt_frac),
! up or down: f in [1/a, a].
pure real(dp) function clipped_factor(set, dt_factor)
      type(timestep_settings), intent(in) :: set
      real(dp), intent(in) :: dt_factor
      real(dp) :: dt_factor_limit
      clipped_factor = dt_factor
      dt_factor_limit = set%timestep_factors(itime_max_dt_frac)
      if (clipped_factor > dt_factor_limit) clipped_factor = dt_factor_limit
      if (clipped_factor < 1.0d0/dt_factor_limit) clipped_factor = 1.0d0/dt_factor_limit
end function clipped_factor

!----------------------------------------------------------------------
! structure_change_limit (legacy ptime.f)
!----------------------------------------------------------------------
! The limit from the change of the structure since the start of the
! step. The largest changes over the n points,
!   dP = max_i |log P_start,i - log P_i|, and likewise dT and dR,
!   dL = max_i |L_start,i - L_i| / ((L_i + L_(i-1))/2),
! with L_2 in place of L_0 at the centre (and 0 where L_i + L_(i-1) <= 0),
! give the factor
!   f = max(dP/atime(itime_dp), dT/atime(itime_dt), dR/atime(itime_dr),
!           dL/(ln 10 atime(itime_dl))),
! where atime(itime_dl) limits Delta log L and Delta L/L ~ ln 10 Delta
! log L. f = 1 when nothing changed, f is clipped (clipped_factor), and
! the limit is previous_dt/f. It dominates on the pre-main sequence. Off
! (dt_unlimited) unless set%structure_limits_active; below dt_floor it
! counts as no limit. log P, log T, log R in log10 cgs, L in Lsun.
pure real(dp) function structure_change_limit(set, previous_dt, &
     log_pressure_start, log_temperature_start, log_radius_start, &
     luminosity_start, log_pressure, log_temperature, log_radius, &
     luminosity, n) result(dt)
      type(timestep_settings), intent(in) :: set
      real(dp), intent(in) :: previous_dt
      real(dp), intent(in) :: log_pressure_start(:), log_temperature_start(:), &
           log_radius_start(:), luminosity_start(:), log_pressure(:), &
           log_temperature(:), log_radius(:), luminosity(:)
      integer, intent(in) :: n
! max_change: the largest change of each quantity over all points
! (1 = log P, 2 = log T, 3 = log R, 4 = L, relative).
      real(dp) :: max_change(4), test_p, test_t, test_r, test_l, dt_factor
      integer :: i

      if (.not. set%structure_limits_active) then
         dt = dt_unlimited
         return
      end if
! The changes at the centre, then the largest over the other points.
      max_change(1) = abs(log_pressure_start(1) - log_pressure(1))
      max_change(2) = abs(log_temperature_start(1) - log_temperature(1))
      max_change(3) = abs(log_radius_start(1) - log_radius(1))
      if (luminosity(1) + luminosity(2) > 0.0d0) then
         max_change(4) = abs((luminosity_start(1) - luminosity(1))*2.d0/ &
              (luminosity(2) + luminosity(1)))
      else
         max_change(4) = 0.0d0
      end if
      do i = 2, n
         test_p = abs(log_pressure_start(i) - log_pressure(i))
         if (max_change(1) <= test_p) max_change(1) = test_p
         test_t = abs(log_temperature_start(i) - log_temperature(i))
         if (max_change(2) <= test_t) max_change(2) = test_t
         test_r = abs(log_radius_start(i) - log_radius(i))
         if (max_change(3) <= test_r) max_change(3) = test_r
         if (luminosity(i) + luminosity(i-1) > 0.0d0) then
            test_l = abs((luminosity_start(i) - luminosity(i))*2.0d0/ &
                 (luminosity(i) + luminosity(i-1)))
         else
            test_l = 0.0d0
         end if
         if (max_change(4) <= test_l) max_change(4) = test_l
      end do
! The factor that brings the largest change of every quantity down to
! its limit, the changes taken to scale with the timestep.
      dt_factor = max_change(1)/set%timestep_factors(itime_dp)
      if (max_change(2)/set%timestep_factors(itime_dt) > dt_factor) &
           dt_factor = max_change(2)/set%timestep_factors(itime_dt)
      if (max_change(3)/set%timestep_factors(itime_dr) > dt_factor) &
           dt_factor = max_change(3)/set%timestep_factors(itime_dr)
      if (max_change(4)/(set%timestep_factors(itime_dl)*ln10_ptime) > dt_factor) &
           dt_factor = max_change(4)/(set%timestep_factors(itime_dl)*ln10_ptime)
! No change since the previous model: keep the previous timestep.
      if (dt_factor == 0.d0) dt_factor = 1.0d0
      dt = previous_dt/clipped_factor(set, dt_factor)
      if (dt < dt_floor) dt = dt_unlimited
end function structure_change_limit

!----------------------------------------------------------------------
! omega_change_limit (legacy wtime.f)
!----------------------------------------------------------------------
! For rotating models, the limit from the change of omega since the
! start of the step:
!   max_domega_frac = max_i 2 |omega_i - omega_start,i|/(omega_i + omega_start,i),
!   f = max_domega_frac/max_domega_global   (legacy DTWIND),
! clipped (clipped_factor), and rotation_dt = previous_dt/f;
! dt_unlimited when omega did not change, and below dt_floor. Without
! rotation the limit is off and max_domega_frac keeps its value.
! max_domega_frac sets the number of angular momentum transport
! sub-steps (evolve_angular_momentum).
pure subroutine omega_change_limit(set, previous_dt, omega_start, omega, n, &
     rotation_dt, max_domega_frac)
      type(timestep_settings), intent(in) :: set
      real(dp), intent(in) :: previous_dt, omega_start(:), omega(:)
      integer, intent(in) :: n
      real(dp), intent(out) :: rotation_dt
      real(dp), intent(inout) :: max_domega_frac
      real(dp) :: test_domega, dt_factor
      integer :: i

      rotation_dt = dt_unlimited
      if (.not. set%rotation_active) return
      max_domega_frac = 2.0d0*abs(omega(1) - omega_start(1))/(omega(1) + omega_start(1))
      do i = 2, n
         test_domega = 2.0d0*abs(omega(i) - omega_start(i))/(omega(i) + omega_start(i))
         if (test_domega > max_domega_frac) max_domega_frac = test_domega
      end do
      dt_factor = max_domega_frac/set%max_domega_global
! No change of omega: no limit.
      if (dt_factor == 0.d0) return
      rotation_dt = previous_dt/clipped_factor(set, dt_factor)
      if (rotation_dt < dt_floor) rotation_dt = dt_unlimited
end subroutine omega_change_limit

!----------------------------------------------------------------------
! hydrogen_burning_limit (legacy xtime.f)
!----------------------------------------------------------------------
! The limit from hydrogen burning, with E_H = 6.0e18 erg/g the energy
! released per gram of hydrogen burned; off when the hydrogen-burning
! luminosity L_H (Lsun) is at most min_burning_luminosity.
! Core burning (central X >= atime(itime_core_min)): the time to burn
!   Delta X = min(atime(itime_dx_core_tot), atime(itime_dx_core_frac) X_cc)
! in the mass m_cc inside the convective-core edge (the centre without a
! convective core) at the luminosity L_cc there,
!   dt = E_H Delta X m_cc/(L_cc Lsun),
! X_cc being X at the core edge.
! Shell burning (central X below that): the smaller of the time to burn
! atime(itime_dx_total) Msun of material with the surface X at L_H,
!   E_H atime(itime_dx_total) X_surf Msun/(L_H Lsun),
! and the time for X at the shell midpoint to fall by
! atime(itime_dx_shell) (when that is positive),
!   atime(itime_dx_shell)/|dX/dt|,
! with dX/dt from the network's reactions (compute_abundance_rates) on
! the zone's rates and its start-of-step abundances (composition_start).
! (The legacy code took dX/dt from the equilibrium burn, with He3 and
! the CN cycle in equilibrium.) Masses in g, L in Lsun, composition in
! the star's slots.
real(dp) function hydrogen_burning_limit(set, hydrogen_luminosity, &
     composition, composition_start, log_density, log_temperature, &
     enclosed_mass, luminosity, core_edge, shell_mid, n) result(dt)
      use net_lib, only: compute_network_rates, compute_abundance_rates, &
           num_rate_rows, n_network, network_species_slot, s_h1
      type(timestep_settings), intent(in) :: set
      real(dp), intent(in) :: hydrogen_luminosity
      real(dp), intent(in) :: composition(:,:), composition_start(:,:), &
           log_density(:), log_temperature(:), enclosed_mass(:), luminosity(:)
      integer, intent(in) :: core_edge, shell_mid, n
! rate_vec, dxdt_per_a: the network's rate vector and dX_i/dt / A_i of
! the network species at the shell midpoint.
      real(dp) :: rate_vec(num_rate_rows), dxdt_per_a(n_network)
! delta_x: the hydrogen budget of the branch taken, Delta X in the core
! branch and atime(itime_dx_total) X_surf Msun/Lsun in the shell branch.
      real(dp) :: delta_x, dx_dt

      dt = dt_unlimited
      if (hydrogen_luminosity <= min_burning_luminosity) return
! Core burning. With a convective core, the limit is the time to burn
! the given amount of hydrogen in the whole core, not only at the centre.
      if (composition(i_h1,1) >= set%timestep_factors(itime_core_min)) then
         delta_x = min(set%timestep_factors(itime_dx_core_tot), &
              set%timestep_factors(itime_dx_core_frac)*composition(i_h1,core_edge))
         dt = (h_burn_energy_per_gram/set%solar_luminosity_cgs)* &
              (enclosed_mass(core_edge)/luminosity(core_edge))*delta_x
         return
      end if
! Shell burning: limit the total mass of hydrogen burned,
      delta_x = set%timestep_factors(itime_dx_total)*composition(i_h1,n)* &
           (set%solar_mass_cgs/set%solar_luminosity_cgs)
      dt = h_burn_energy_per_gram*delta_x/hydrogen_luminosity
! and the depletion of X at the shell midpoint, with dX/dt from the
! network (A = 1 for hydrogen).
      call compute_network_rates(log_density(shell_mid), log_temperature(shell_mid), &
           composition(:,shell_mid), rate_vec)
      call compute_abundance_rates(rate_vec, &
           composition_start(network_species_slot, shell_mid), dxdt_per_a)
      dx_dt = dxdt_per_a(s_h1)
      if (dx_dt < 0.0d0 .and. set%timestep_factors(itime_dx_shell) > 0.0d0) &
           dt = min(dt, abs(set%timestep_factors(itime_dx_shell)/dx_dt))
end function hydrogen_burning_limit

!----------------------------------------------------------------------
! helium_burning_limit (legacy ytime.f; the helium-shell branch is new)
!----------------------------------------------------------------------
! The limit from helium burning, with E_He = 5.85e17 erg/g the energy
! released per gram of helium burned; off when the triple-alpha
! luminosity (Lsun) is at most min_burning_luminosity.
! Helium flash (helium_flash_applies: central log rho >= 5): the time
! for the triple-alpha rate at the temperature maximum,
! triple_alpha_rate (erg/g/s, the caller's, from the hottest zone,
! hottest_zone), to release 1e15 erg/g (about 1.7e-3 in Y),
!   dt = 1e15/max(triple_alpha_rate, 1e-22).
! Core burning (Y_core = 1 - Z at the convective-core edge >=
! atime(itime_core_min); X is negligible in a helium-burning core): the
! time to burn
!   Delta Y = min(atime(itime_dy_core_tot), atime(itime_dy_core_frac) Y_core)
! in the mass m_cc inside the core edge at the luminosity L_cc there,
!   dt = E_He Delta Y m_cc/(L_cc Lsun).
! Shell burning (core helium exhausted), mirroring the shell branch of
! hydrogen_burning_limit: at the zone just below the H-burning shell's
! bottom (the top of the helium core), with its m, L and Y, the smaller
! of
!   E_He atime(itime_dy_total) Msun/(L Lsun)   (Msun of helium burned per step),
!   E_He atime(itime_dy_shell) Y m/(L Lsun)    (fraction of Y burned),
! each used only when its factor is positive, and only where L > 0;
! dt_unlimited otherwise. (2026 deviation, bugsweep sec-10/11: ytime.f
! reused the previous call's timestep here, times m/L again, which gave
! about 0 s and, once the SAVE that kept it was gone, 0 and then NaN:
! the run_from_zahb_to_tahb failure.)
real(dp) function helium_burning_limit(set, triple_alpha_luminosity, &
     triple_alpha_rate, composition, log_density, enclosed_mass, &
     luminosity, core_edge, shell_bottom) result(dt)
      type(timestep_settings), intent(in) :: set
      real(dp), intent(in) :: triple_alpha_luminosity, triple_alpha_rate
      real(dp), intent(in) :: composition(:,:), log_density(:), &
           enclosed_mass(:), luminosity(:)
      integer, intent(in) :: core_edge, shell_bottom
      real(dp) :: helium_energy_gen, core_helium_fraction
      integer :: shell_zone

      dt = dt_unlimited
      if (triple_alpha_luminosity <= min_burning_luminosity) return
      if (helium_flash_applies(log_density)) then
         helium_energy_gen = triple_alpha_rate
         if (helium_energy_gen < min_triple_alpha_rate) helium_energy_gen = min_triple_alpha_rate
         dt = flash_energy_per_gram/helium_energy_gen
         return
      end if
      core_helium_fraction = 1d0 - composition(i_metals,core_edge)
      if (core_helium_fraction >= set%timestep_factors(itime_core_min)) then
         dt = min(set%timestep_factors(itime_dy_core_tot), &
              set%timestep_factors(itime_dy_core_frac)*core_helium_fraction)
         dt = (he_burn_energy_per_gram/set%solar_luminosity_cgs)*dt* &
              (enclosed_mass(core_edge)/luminosity(core_edge))
         return
      end if
      shell_zone = max(shell_bottom-1, 1)
      if (luminosity(shell_zone) > 0.0d0) then
         if (set%timestep_factors(itime_dy_total) > 0.0d0) dt = min(dt, &
              (he_burn_energy_per_gram/set%solar_luminosity_cgs)*set%timestep_factors(itime_dy_total)* &
              (set%solar_mass_cgs/luminosity(shell_zone)))
         if (set%timestep_factors(itime_dy_shell) > 0.0d0) dt = min(dt, &
              (he_burn_energy_per_gram/set%solar_luminosity_cgs)*set%timestep_factors(itime_dy_shell)* &
              composition(i_he4,shell_zone)* &
              (enclosed_mass(shell_zone)/luminosity(shell_zone)))
      end if
end function helium_burning_limit

! Whether helium burning is a flash: central log rho >= 5.
pure logical function helium_flash_applies(log_density)
      real(dp), intent(in) :: log_density(:)
      helium_flash_applies = log_density(1) >= flash_log_central_density
end function helium_flash_applies

! The zone of the largest log T among the first n (the first of equal
! ones; zone 1 when none exceeds 0).
pure integer function hottest_zone(log_temperature, n)
      real(dp), intent(in) :: log_temperature(:)
      integer, intent(in) :: n
      real(dp) :: max_temp
      integer :: zone_idx
      max_temp = 0.0d0
      hottest_zone = 1
      do zone_idx = 1, n
         if (log_temperature(zone_idx) > max_temp) then
            hottest_zone = zone_idx
            max_temp = log_temperature(zone_idx)
         end if
      end do
end function hottest_zone

!----------------------------------------------------------------------
! hr_motion_limit (legacy entime.f)
!----------------------------------------------------------------------
! The limit from the motion in the HR diagram since the start of the
! step, which keeps a small envelope triangle from being used for
! extrapolation:
!   f = max(|Delta log Teff|/tri_delta_teffl, |Delta log L|/tri_delta_logl),
! f = 1 when nothing changed, clipped (clipped_factor), and
! dt = previous_dt/f. Applied (else dt_unlimited) when
! set%use_envelope_triangle_dt is on and the structure-based limits are
! off. The start-of-step surface L (Lsun) is the one saved before the
! rezone (seed_new_point_run in rezoning).
real(dp) function hr_motion_limit(set, previous_dt, log_teff_start, &
     surface_luminosity_start, log_teff, surface_luminosity) result(dt)
      use math_lib
      type(timestep_settings), intent(in) :: set
      real(dp), intent(in) :: previous_dt, log_teff_start, &
           surface_luminosity_start, log_teff, surface_luminosity
      real(dp) :: teffl_change, logl_change, dt_factor

      dt = dt_unlimited
      if (.not. set%use_envelope_triangle_dt .or. set%structure_limits_active) return
      teffl_change = abs(log_teff_start - log_teff)
      logl_change = abs(log10(surface_luminosity_start) - log10(surface_luminosity))
      dt_factor = teffl_change/set%tri_delta_teffl
      if (logl_change/set%tri_delta_logl > dt_factor) dt_factor = logl_change/set%tri_delta_logl
! No change since the previous model: keep the previous timestep.
      if (dt_factor == 0.0d0) dt_factor = 1.0d0
      dt = previous_dt/clipped_factor(set, dt_factor)
end function hr_motion_limit

!----------------------------------------------------------------------
! The timestep: the smallest of the limits and, for a previous timestep
! above dt_floor, of its allowed growth,
!   dt = min(atime(itime_max_dt_frac) previous_dt, limits),
! the smallest limit alone otherwise.
pure real(dp) function tightest_timestep(set, previous_dt, limits) result(dt)
      type(timestep_settings), intent(in) :: set
      real(dp), intent(in) :: previous_dt, limits(:)
      if (previous_dt > dt_floor) then
         dt = set%timestep_factors(itime_max_dt_frac)*previous_dt
         dt = min(dt, minval(limits))
      else
         dt = minval(limits)
      end if
end function tightest_timestep

end module timestep_limits
