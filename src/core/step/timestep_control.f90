!----------------------------------------------------------------------
! timestep_control
!----------------------------------------------------------------------
! The next model's timestep: compute_timestep takes the most stringent of
! the limits from the structure change, the change of omega, hydrogen and
! helium burning and the HR-diagram motion, within the growth limit, the
! fixed-timestep override and the end age. The limits are this module's
! own procedures.
module timestep_control
      implicit none
      private
      public :: compute_timestep

contains

!----------------------------------------------------------------------
! compute_timestep (legacy htimer.f)
!----------------------------------------------------------------------
! The timestep for the next model: the smallest of these limits, each
! set by the change over the step just taken (a criterion that is off or
! does not apply gives dt_unlimited):
!   structure_dt  the change of log P, log T, log R and L at any point
!                 (timestep_limit_structure, legacy ptime), while the
!                 structure-based limits are on;
!   rotation_dt   the change of omega, rotating models only
!                 (timestep_limit_omega, legacy wtime);
!   hydrogen_dt   hydrogen burning (timestep_limit_hburn, legacy xtime);
!   helium_dt     helium burning (timestep_limit_heburn, legacy ytime);
!   envelope_dt   the motion in the HR diagram against the envelope
!                 triangle (timestep_limit_hr, legacy entime), when on
!                 and the structure-based limits are off;
! and, for a previous timestep above 1 s, atime(itime_max_dt_frac) times
! the previous timestep.
! A fixed timestep (timestep_override_active for this kind card)
! replaces all of these; caveat user. When evolving to a fixed age
! (end_age_stop_active with target_end_age > 0), the timestep is cut so
! that the model does not pass that age.
! A negative previous_timestep marks a model being relaxed rather than
! evolved: the limits are skipped and only its magnitude is kept.
! On return previous_timestep holds the new timestep (s), also in
! chosen_dt and, in years, timestep_years.
subroutine compute_timestep(previous_timestep, chosen_dt, num_points, log_density, &
     luminosity, enclosed_mass, log_temperature, composition, &
     convective_core_edge_zone, h_shell_midpoint_zone, &
     luminosity_components, age_gyr, timestep_years, kind_card_index, &
     log_pressure, log_radius, omega, max_domega_frac, h_shell_zone_begin, &
     log_teff, disable_structure_dt_limits)
      use star_info_lib, only: star, json
      use controls_buffer, only: itime_max_dt_frac
      use physical_constants
      implicit none

      double precision, intent(inout) :: previous_timestep
! chosen_dt: the timestep actually adopted (seconds); the caller stores
! it in star%chosen_dt.
      double precision, intent(out) :: chosen_dt
      integer, intent(in) :: num_points
      double precision, intent(in) :: log_density(json), luminosity(json), &
           enclosed_mass(json), log_temperature(json)
      double precision, intent(in) :: composition(15,json)
      integer, intent(in) :: convective_core_edge_zone, &
           h_shell_midpoint_zone
      double precision, intent(in) :: luminosity_components(8)
      double precision, intent(in) :: age_gyr
      double precision, intent(out) :: timestep_years
      integer, intent(in) :: kind_card_index
      double precision, intent(in) :: log_pressure(json), log_radius(json), &
           omega(json)
      double precision, intent(out) :: max_domega_frac
      integer, intent(in) :: h_shell_zone_begin
      double precision, intent(in) :: log_teff
! disable_structure_dt_limits: set when this call switches the
! structure-based limits off; the driver then clears
! star%job%use_structure_dt_limits.
      logical, intent(out) :: disable_structure_dt_limits
      double precision :: structure_dt, rotation_dt, hydrogen_dt, helium_dt, &
           hydrogen_luminosity, envelope_dt, time_left_years
! structure_limits_active: whether the structure-based limits apply to
! this step: star%job%use_structure_dt_limits, unless
! decide_structure_limits switches them off (log Tc > 7.1 with
! L_grav < 0). The job flag itself is left to the driver.
      logical :: structure_limits_active
! dt_unlimited: "no limit from this criterion" sentinel (seconds).
      double precision, parameter :: dt_unlimited = 1.0d20

      disable_structure_dt_limits = .false.
      if (previous_timestep.ge.0.0d0) then
! A fixed timestep for this kind card overrides the limits.
      if(star%job%timestep_override_active(kind_card_index)) then
       chosen_dt = star%job%timestep_override(kind_card_index)*seconds_per_year
       previous_timestep = chosen_dt
      else
       call decide_structure_limits()
       call limit_by_structure_change()
       call limit_by_omega_change()
       call limit_by_hydrogen_burning()
       call limit_by_helium_burning()
       call limit_by_hr_motion()
       call take_tightest_limit()
      end if
      end if
      previous_timestep = abs(previous_timestep)
      timestep_years = previous_timestep/seconds_per_year
      call clip_to_end_age()
      return

contains

! ---------------------------------------------------------------
! Whether the structure-based limits apply to this step: they are switched
! off (MHP 9/01) once the centre is hotter than log Tc = 7.1 while the
! gravothermal luminosity L_grav (luminosity_components(7)) is negative.
subroutine decide_structure_limits()
      structure_limits_active = star%job%use_structure_dt_limits
      if(structure_limits_active) then
         if(log_temperature(1).gt.7.1d0 .and. luminosity_components(7).lt. 0.0d0) then
            write(*, 100) log_temperature(1), luminosity_components(7)
 100        format(1x,'timestep control: structure-based limits', &
                 ' disabled (log Tc =',f7.4,', L_grav =',es10.2,' Lsun)')
            structure_limits_active = .false.
            disable_structure_dt_limits = .true.
         endif
      endif
end subroutine decide_structure_limits

! ---------------------------------------------------------------
! The limit from the change of the structure variables
! (timestep_limit_structure, legacy ptime), unlimited when the
! structure-based limits are off. With no change since the previous
! model (as on the first call of a kind card) it is the timestep stored
! in the model.
subroutine limit_by_structure_change()
      if(structure_limits_active) then
       call timestep_limit_structure(previous_timestep,luminosity,log_pressure,log_radius,log_temperature,num_points,structure_dt)
      else
       structure_dt = dt_unlimited
      endif
! A limit below 1 s (as from a zero timestep stored in the model) would
! stall the run, so it counts as no limit (MHP 6/90).
      if(structure_dt.lt.1.0d0) structure_dt=dt_unlimited
end subroutine limit_by_structure_change

! ---------------------------------------------------------------
! The limit from the change of omega (timestep_limit_omega, legacy
! wtime), rotating models only; unlimited when omega did not change. It
! also returns max_domega_frac, the largest fractional change of omega
! since the previous model, which sets the number of angular momentum
! transport sub-steps (evolve_angular_momentum).
subroutine limit_by_omega_change()
      if(star%job%rotation_active) then
       call timestep_limit_omega(previous_timestep,num_points,omega,rotation_dt,max_domega_frac)
      else
       rotation_dt = dt_unlimited
      endif
! A limit below 1 s counts as no limit, as for the structure (MHP 6/90).
      if(rotation_dt.lt.1.0d0) rotation_dt=dt_unlimited
end subroutine limit_by_omega_change

! ---------------------------------------------------------------
! The limit from hydrogen burning (timestep_limit_hburn, legacy xtime),
! for stars that burn hydrogen.
subroutine limit_by_hydrogen_burning()
! The hydrogen-burning luminosity, L_H = L_pp1 + L_pp2 + L_pp3 + L_CNO
! (Lsun); no limit when it is below 1e-34.
      hydrogen_luminosity = luminosity_components(1) + luminosity_components(2) + &
           luminosity_components(3) + luminosity_components(4)
      if(hydrogen_luminosity.gt.1.0d-34) then
       call timestep_limit_hburn(log_density,composition,luminosity,enclosed_mass, &
            log_temperature,hydrogen_luminosity,convective_core_edge_zone, &
            h_shell_midpoint_zone,num_points,hydrogen_dt)
      else
       hydrogen_dt = dt_unlimited
      endif
end subroutine limit_by_hydrogen_burning

! ---------------------------------------------------------------
! The limit from helium burning (timestep_limit_heburn, legacy ytime),
! for stars that burn helium (triple-alpha luminosity
! luminosity_components(5) above 1e-34 Lsun).
subroutine limit_by_helium_burning()
      if(luminosity_components(5).gt.1.0d-34) then
       call timestep_limit_heburn(composition,log_density,luminosity, &
            enclosed_mass,log_temperature,convective_core_edge_zone, &
            num_points,helium_dt,h_shell_zone_begin)
      else
       helium_dt = dt_unlimited
      endif
end subroutine limit_by_helium_burning

! ---------------------------------------------------------------
! The limit from the motion in the HR diagram (timestep_limit_hr, legacy
! entime; JVS 04/14): the model may move no more than tri_delta_logl in
! log L or tri_delta_teffl in log Teff per step, the size of the envelope
! triangle. Applied when use_envelope_triangle_dt is on and the
! structure-based limits are off.
subroutine limit_by_hr_motion()
      if(star%ctrl%use_envelope_triangle_dt .and. .not. structure_limits_active) then
       call timestep_limit_hr(previous_timestep,luminosity,log_teff,num_points,envelope_dt)
      else
       envelope_dt = dt_unlimited
      endif
end subroutine limit_by_hr_motion

! ---------------------------------------------------------------
! The timestep: the smallest of the limits and, for a previous timestep
! above 1 s, of its allowed growth,
!   dt = min(atime(itime_max_dt_frac) dt_prev, hydrogen_dt, helium_dt,
!            rotation_dt, structure_dt, envelope_dt).
subroutine take_tightest_limit()
      if(previous_timestep.gt.1.0d0) then
         previous_timestep = star%ctrl%atime(itime_max_dt_frac)*previous_timestep
         previous_timestep = min(previous_timestep,hydrogen_dt,helium_dt,rotation_dt,structure_dt,envelope_dt)
      else
! A zero (or sub-second) previous timestep: the smallest limit alone.
         previous_timestep = min(hydrogen_dt,helium_dt,rotation_dt,structure_dt,envelope_dt)
      endif
      chosen_dt = previous_timestep
end subroutine take_tightest_limit

! ---------------------------------------------------------------
! When evolving to a given age, shorten the timestep so the model does not
! pass it: dt = min(dt, target_end_age - age), in years. Only for a true
! stop age (target_end_age > 0), since end_age_stop_active also covers
! other stop conditions (MHP 10/24).
subroutine clip_to_end_age()
      if(star%job%end_age_stop_active(kind_card_index) .and. star%job%target_end_age(kind_card_index).gt.0.0d0) then
       time_left_years = star%job%target_end_age(kind_card_index) - age_gyr*1.0d9
       if(time_left_years.lt. timestep_years) then
         timestep_years = time_left_years
         chosen_dt = timestep_years*seconds_per_year
         previous_timestep = chosen_dt
       endif
      endif
end subroutine clip_to_end_age

end subroutine compute_timestep

!----------------------------------------------------------------------
! timestep_limit_structure (legacy ptime.f)
!----------------------------------------------------------------------
! The timestep limit from the change of the structure since the start of
! the step (star%*_start). The largest changes over all points,
!   dP = max_i |log P_start,i - log P_i|, and likewise dT and dR,
!   dL = max_i |L_start,i - L_i| / ((L_i + L_(i-1))/2),
! with L_2 in place of L_0 at the centre, give the factor
!   f = max(dP/atime(itime_dp), dT/atime(itime_dt), dR/atime(itime_dr),
!           dL/(ln 10 atime(itime_dl))),
! where atime(itime_dl) limits Delta log L and Delta L/L ~ ln 10 Delta
! log L. f = 1 when nothing changed; f is clipped to [1/a, a] with
! a = atime(itime_max_dt_frac); and struct_dt = previous_timestep/f.
! This limit dominates on the pre-main sequence, where it must be on
! (star%job%use_structure_dt_limits, legacy LPTIME).
subroutine timestep_limit_structure(previous_timestep, luminosity, log_pressure, log_radius, &
     log_temperature, num_points, struct_dt)

      use star_info_lib, only: star, json
      use controls_buffer, only: itime_dl, itime_dp, itime_dr, itime_dt, itime_max_dt_frac
      implicit none

! previous_timestep: the previous model's timestep (s).
! luminosity: L/Lsun at each point of the current model (linear, not
!   log).
! log_pressure, log_radius, log_temperature: log P, log R and log T at
!   each point of the current model.
! num_points: the number of points.
! struct_dt: the structure-based timestep limit (s).
      double precision, intent(in) :: previous_timestep
      double precision, intent(in) :: luminosity(json), log_pressure(json), &
           log_radius(json), log_temperature(json)
      integer, intent(in) :: num_points
      double precision, intent(out) :: struct_dt
! max_change: the largest change of each quantity over all points
! (1 = log P, 2 = log T, 3 = log R, 4 = L, relative).
      double precision :: max_change(4)
      integer :: i
      double precision :: test_p, test_t, test_r, test_l, dt_factor, &
           dt_factor_limit

! The changes at the centre, then the largest over the other points:
! log P,
      max_change(1)=abs(star%logP_start(1)-log_pressure(1))
! log T,
      max_change(2)=abs(star%logT_start(1)-log_temperature(1))
! log R,
      max_change(3)=abs(star%logR_start(1)-log_radius(1))
! and L, relative to the mean L of the point and its neighbour.
      if(luminosity(1)+luminosity(2).gt.0.0d0) then
       max_change(4)=abs((star%luminosity_lsun_start(1)-luminosity(1))*2.d0/(luminosity(2)+luminosity(1)))
      else
       max_change(4) = 0.0d0
      endif
      do i = 2,num_points
       test_p = abs(star%logP_start(i)-log_pressure(i))
       if(max_change(1).le.test_p) then
          max_change(1) = test_p
       endif
       test_t = abs(star%logT_start(i)-log_temperature(i))
       if(max_change(2).le.test_t) then
          max_change(2) = test_t
       endif
       test_r = abs(star%logR_start(i)-log_radius(i))
       if(max_change(3).le.test_r) then
          max_change(3) = test_r
       endif
       if(luminosity(i)+luminosity(i-1).gt.0.0d0) then
          test_l = abs((star%luminosity_lsun_start(i)-luminosity(i))*2.0d0/(luminosity(i)+luminosity(i-1)))
       else
          test_l = 0.0d0
       endif
       if(max_change(4).le.test_l) then
          max_change(4) = test_l
       endif
      end do
! The factor f that brings the largest change of every quantity down to
! its limit, the changes taken to scale with the timestep (2.3026 is
! ln 10).
      dt_factor = max_change(1)/star%ctrl%atime(itime_dp)
      if (max_change(2)/star%ctrl%atime(itime_dt).gt.dt_factor) dt_factor=max_change(2)/star%ctrl%atime(itime_dt)
      if (max_change(3)/star%ctrl%atime(itime_dr).gt.dt_factor) dt_factor=max_change(3)/star%ctrl%atime(itime_dr)
      if (max_change(4)/(star%ctrl%atime(itime_dl)*2.3026d0).gt.dt_factor) then
        dt_factor=max_change(4)/(star%ctrl%atime(itime_dl)*2.3026d0)
      endif
! No change since the previous model: keep the previous timestep.
      if (dt_factor.eq.0.d0) dt_factor=1.0d0
! Let the timestep change by at most a factor atime(itime_max_dt_frac)
! (the global timestep limiter) up or down.
      dt_factor_limit = star%ctrl%atime(itime_max_dt_frac)
      if (dt_factor.gt.dt_factor_limit) dt_factor=dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor=1.0d0/dt_factor_limit
      struct_dt = previous_timestep/dt_factor

      return
end subroutine timestep_limit_structure

!----------------------------------------------------------------------
! timestep_limit_omega (legacy wtime.f)
!----------------------------------------------------------------------
! For rotating models, the timestep limit from the change of omega since
! the start of the step (star%old_omega):
!   max_domega_frac = max_i 2 |omega_i - omega_old,i|/(omega_i + omega_old,i),
!   f = max_domega_frac/max_domega_global   (legacy DTWIND),
! clipped to [1/a, a] with a = atime(itime_max_dt_frac), and
! rotation_dt = previous_timestep/f; dt_unlimited when omega did not
! change.
subroutine timestep_limit_omega(previous_timestep, num_points, omega, rotation_dt, &
     max_domega_frac)

      use star_info_lib, only: star, json
      use controls_buffer, only: itime_max_dt_frac
      implicit none

! previous_timestep: the previous model's timestep (s).
! num_points: the number of points.
! omega: the angular velocity at each point of the current model.
! rotation_dt: the rotation-based timestep limit (s).
! max_domega_frac: the largest fractional change of omega.
      double precision, intent(in) :: previous_timestep
      integer, intent(in) :: num_points
      double precision, intent(in) :: omega(json)
      double precision, intent(out) :: rotation_dt
      double precision, intent(out) :: max_domega_frac
      integer :: start_index, i
      double precision :: test_domega, dt_factor, dt_factor_limit
! dt_unlimited: "no limit from this criterion" sentinel (seconds).
      double precision, parameter :: dt_unlimited = 1.0d20

      start_index = 1
      max_domega_frac = 2.0d0*abs(omega(start_index)-star%old_omega(start_index))/ &
           (omega(start_index)+star%old_omega(start_index))
      do i = start_index+1,num_points
         test_domega=2.0d0*abs(omega(i)-star%old_omega(i))/(omega(i)+star%old_omega(i))
         if(test_domega.gt.max_domega_frac) then
            max_domega_frac = test_domega
         endif
      end do
      dt_factor = max_domega_frac/star%job%max_domega_global
! No change of omega: no limit.
      if (dt_factor.eq.0.d0)then
          rotation_dt = dt_unlimited
          return
      endif
! Let the timestep change by at most a factor atime(itime_max_dt_frac)
! (the global timestep limiter) up or down.
      dt_factor_limit = star%ctrl%atime(itime_max_dt_frac)
      if (dt_factor.gt.dt_factor_limit) dt_factor=dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor=1.0d0/dt_factor_limit
      rotation_dt = previous_timestep/dt_factor
      return
end subroutine timestep_limit_omega

!----------------------------------------------------------------------
! timestep_limit_hburn (legacy xtime.f)
!----------------------------------------------------------------------
! The timestep limit from hydrogen burning, with E_H = 6.0e18 erg/g the
! energy released per gram of hydrogen burned.
! Core burning (central X >= atime(itime_core_min)): the time to burn
!   Delta X = min(atime(itime_dx_core_tot), atime(itime_dx_core_frac) X_cc)
! in the mass m_cc inside the convective-core edge (the centre without a
! convective core) at the luminosity L_cc there,
!   hydrogen_dt = E_H Delta X m_cc/(L_cc Lsun),
! X_cc being X at the core edge.
! Shell burning (central X below that): the smaller of the time to burn
! atime(itime_dx_total) Msun of material with the surface X at the
! hydrogen-burning luminosity L_H (Lsun),
!   E_H atime(itime_dx_total) X_surf Msun/(L_H Lsun),
! and the time for X at the shell midpoint to fall by
! atime(itime_dx_shell),
!   atime(itime_dx_shell)/|dX/dt|,
! with dX/dt from the network's reactions (compute_abundance_rates) on
! the zone's rates and its start-of-step abundances. (The legacy code
! took dX/dt from the equilibrium burn, with He3 and the CN cycle in
! equilibrium.)
subroutine timestep_limit_hburn(log_density, composition, luminosity, enclosed_mass, &
     log_temperature, hydrogen_luminosity, &
     convective_core_edge_zone, h_shell_midpoint_zone, num_points, &
     hydrogen_dt)
      use star_info_lib, only: star, json, i_h1
      use net_lib, only: compute_network_rates, compute_abundance_rates, &
           num_rate_rows, n_network, network_species_slot, s_h1
      use controls_buffer, only: itime_core_min, itime_dx_core_frac, itime_dx_core_tot, &
           itime_dx_shell, itime_dx_total
      use physical_constants
      implicit none

      double precision, intent(in) :: log_density(json)
      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: luminosity(json)
      double precision, intent(in) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json)
      double precision, intent(in) :: hydrogen_luminosity
      integer, intent(in) :: convective_core_edge_zone, &
           h_shell_midpoint_zone, num_points
      double precision, intent(out) :: hydrogen_dt
! rate_vec, dxdt_per_a: the network's rate vector and dX_i/dt / A_i of
! the network species at the shell midpoint.
      double precision :: rate_vec(num_rate_rows), dxdt_per_a(n_network)
! delta_x: the hydrogen budget of the branch taken, Delta X in the core
! branch and atime(itime_dx_total) X_surf Msun/Lsun in the shell branch.
      double precision :: delta_x
! local_log_density, local_log_temperature: log rho and log T at the
! shell midpoint, for compute_network_rates.
      double precision :: local_log_density, local_log_temperature
      integer :: zone_end
! The energy released per gram of hydrogen burned (erg/g).
      double precision, parameter :: h_burn_energy_per_gram = 6.00d18
      double precision :: dx_dt, shell_dt_x_depletion

! Core burning. With a convective core, the limit is the time to burn
! the given amount of hydrogen in the whole core, not only at the centre.
      if(composition(i_h1,1).ge.star%ctrl%atime(itime_core_min)) then
       delta_x = min(star%ctrl%atime(itime_dx_core_tot),star%ctrl%atime(itime_dx_core_frac)*composition(i_h1,convective_core_edge_zone))
       hydrogen_dt =(h_burn_energy_per_gram/star%solar_luminosity_cgs)* &
            (enclosed_mass(convective_core_edge_zone)/luminosity(convective_core_edge_zone))*delta_x
       return
      endif
! Shell burning: limit the total mass of hydrogen burned,
      delta_x = star%ctrl%atime(itime_dx_total)*composition(i_h1,num_points)*(star%solar_mass_cgs/star%solar_luminosity_cgs)
      hydrogen_dt = h_burn_energy_per_gram*delta_x/hydrogen_luminosity
! and the depletion of X at the shell midpoint, with dX/dt from the
! network (A = 1 for hydrogen).
      zone_end=h_shell_midpoint_zone
      local_log_density = log_density(zone_end)
      local_log_temperature = log_temperature(zone_end)
      call compute_network_rates(local_log_density,local_log_temperature,composition(:,zone_end), &
           rate_vec)
      call compute_abundance_rates(rate_vec, &
           star%xa_start(network_species_slot, zone_end), dxdt_per_a)
      dx_dt = dxdt_per_a(s_h1)
      if(dx_dt.lt.0.0d0 .and. star%ctrl%atime(itime_dx_shell).gt.0.0d0) then
         shell_dt_x_depletion = abs(star%ctrl%atime(itime_dx_shell)/dx_dt)
         hydrogen_dt = min(hydrogen_dt,shell_dt_x_depletion)
      endif
      return
end subroutine timestep_limit_hburn

!----------------------------------------------------------------------
! timestep_limit_heburn (legacy ytime.f; the helium-shell branch is new)
!----------------------------------------------------------------------
! The timestep limit from helium burning, with E_He = 5.85e17 erg/g the
! energy released per gram of helium burned.
! Helium flash (log rho_c >= 5, taken as a giant in the flash rather
! than a horizontal-branch star): the time for the triple-alpha rate at
! the temperature maximum to release 1e15 erg/g (about 1.7e-3 in Y),
!   helium_dt = 1e15/eps_3alpha(T_max),   eps_3alpha >= 1e-22.
! Core burning (Y_core = 1 - Z at the convective-core edge
! >= atime(itime_core_min)): the time to burn
!   Delta Y = min(atime(itime_dy_core_tot), atime(itime_dy_core_frac) Y_core)
! in the mass m_cc inside the convective-core edge at the luminosity
! L_cc there,
!   helium_dt = E_He Delta Y m_cc/(L_cc Lsun).
! Shell burning (core helium exhausted): see that branch below.
subroutine timestep_limit_heburn(composition, log_density, luminosity, &
     enclosed_mass, log_temperature, convective_core_edge_zone, &
     num_points, helium_dt, h_shell_zone_begin)
      use net_lib, only: nuclear_energy_result
      use star_info_lib, only: star, json, i_he4, i_h2, i_metals
      use controls_buffer, only: itime_core_min, itime_dy_core_frac, itime_dy_core_tot, &
           itime_dy_shell, itime_dy_total
      use burn_step, only: compute_zone_energy
      implicit none

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_density(json), luminosity(json), &
           enclosed_mass(json), log_temperature(json)
      integer, intent(in) :: convective_core_edge_zone, num_points
! helium_dt: the helium-burning timestep limit (s), assigned on every
! path.
      double precision, intent(out) :: helium_dt
      integer, intent(in) :: h_shell_zone_begin
      double precision :: max_temp, local_log_density, local_log_temperature
      integer :: max_temp_zone, zone_idx, shell_zone
! zone_xa: the hottest zone's composition, with deuterium zeroed unless
! the extended composition is used.
      double precision :: zone_xa(15)
! burn: the zone's energy record, of which only the triple-alpha term
! (burn%triple_alpha_energy_gen) is used here.
      type(nuclear_energy_result) :: burn
      double precision :: helium_energy_gen
      double precision :: core_helium_fraction
! The energy released per gram of helium burned (erg/g).
      double precision, parameter :: he_burn_energy_per_gram = 5.85d17
! dt_unlimited: "no limit from this criterion" sentinel (seconds).
      double precision, parameter :: dt_unlimited = 1.0d20

      if (log_density(1).ge.5.0d0) then
! Helium flash: find the temperature maximum,
       max_temp = 0.0d0
       max_temp_zone = 1
       do zone_idx = 1,num_points
          if (log_temperature(zone_idx).gt.max_temp) then
             max_temp_zone = zone_idx
             max_temp = log_temperature(zone_idx)
          endif
       end do
! and the triple-alpha energy generation there.
       local_log_density = log_density(max_temp_zone)
       local_log_temperature = log_temperature(max_temp_zone)
       zone_xa = composition(:,max_temp_zone)
       if (.not.star%job%use_extended_composition) zone_xa(i_h2) = 0.0d0
       call compute_zone_energy(local_log_density, local_log_temperature, &
            zone_xa, max_temp_zone, burn)
         helium_energy_gen = burn%triple_alpha_energy_gen
       if(helium_energy_gen.lt.1.d-22) helium_energy_gen = 1.d-22
       helium_dt=1.0d15/helium_energy_gen
      else
! No flash. core_helium_fraction = 1 - Z at the convective-core edge,
! which approximates the core Y since X is negligible in a
! helium-burning core.
       core_helium_fraction = 1d0 - composition(i_metals,convective_core_edge_zone)
       if(core_helium_fraction.ge.star%ctrl%atime(itime_core_min)) then
            helium_dt = min(star%ctrl%atime(itime_dy_core_tot),star%ctrl%atime(itime_dy_core_frac)*core_helium_fraction)
            helium_dt = (he_burn_energy_per_gram/star%solar_luminosity_cgs)*helium_dt* &
               (enclosed_mass(convective_core_edge_zone)/luminosity(convective_core_edge_zone))
         else
! Core helium exhausted: the helium-shell limit, mirroring the shell
! branch of timestep_limit_hburn. At the zone just below the H-burning
! shell (the top of the helium core), with its m, L and Y, helium_dt is
! the smaller of
!   E_He atime(itime_dy_total) Msun/(L Lsun)   (time_dy_total: Msun of
!                                              helium burned per step),
!   E_He atime(itime_dy_shell) Y m/(L Lsun)    (time_dy_shell: fraction
!                                              of Y burned),
! each used only when its control is positive; dt_unlimited otherwise.
! (2026 deviation, bugsweep sec-10/11: ytime.f reused the previous
! call's timestep here, times m/L again, which gave about 0 s and, once
! the SAVE that kept it was gone, 0 and then NaN: the
! run_from_zahb_to_tahb failure.)
            shell_zone = max(h_shell_zone_begin-1, 1)
            helium_dt = dt_unlimited
            if (luminosity(shell_zone) .gt. 0.0d0) then
               if (star%ctrl%atime(itime_dy_total) .gt. 0.0d0) helium_dt = min(helium_dt, &
                    (he_burn_energy_per_gram/star%solar_luminosity_cgs)*star%ctrl%atime(itime_dy_total)* &
                    (star%solar_mass_cgs/luminosity(shell_zone)))
               if (star%ctrl%atime(itime_dy_shell) .gt. 0.0d0) helium_dt = min(helium_dt, &
                    (he_burn_energy_per_gram/star%solar_luminosity_cgs)*star%ctrl%atime(itime_dy_shell)* &
                    composition(i_he4,shell_zone)* &
                    (enclosed_mass(shell_zone)/luminosity(shell_zone)))
            end if
         endif
      endif
      return
end subroutine timestep_limit_heburn

!----------------------------------------------------------------------
! timestep_limit_hr (legacy entime.f)
!----------------------------------------------------------------------
! The timestep limit from the motion in the HR diagram since the start
! of the step, which keeps a small envelope triangle from being used
! for extrapolation:
!   f = max(|Delta log Teff|/tri_delta_teffl, |Delta log L|/tri_delta_logl),
! f = 1 when nothing changed, clipped to [1/a, a] with
! a = atime(itime_max_dt_frac), and envelope_dt = previous_timestep/f.
subroutine timestep_limit_hr(previous_timestep, luminosity, log_teff, &
     num_points, envelope_dt)

      use star_info_lib, only: star, json
      use controls_buffer, only: itime_max_dt_frac
      use math_lib
      implicit none

! previous_timestep: the previous model's timestep (s).
! luminosity: L/Lsun at each point of the current model.
! log_teff: log Teff of the current model.
! num_points: the number of points.
! envelope_dt: the HR-motion timestep limit (s).
      double precision, intent(in) :: previous_timestep
      double precision, intent(in) :: luminosity(json)
      double precision, intent(in) :: log_teff
      integer, intent(in) :: num_points
      double precision, intent(out) :: envelope_dt
      double precision :: teffl_change, logl_change, dt_factor, &
           dt_factor_limit

! The changes of log Teff and of the surface log L since the start of the
! step; the start-of-step surface L is the one saved before the rezone
! (seed_new_point_run in rezoning). (It used to be read as
! luminosity_lsun_start(nz_start), nz_start the point count before the
! rezone, but the rezone re-splines luminosity_lsun_start onto the new
! points, so that read an interior point, or a stale one beyond the new
! count, whenever the count changed.)
      teffl_change = abs(star%log_Teff_start - log_teff)
      logl_change = abs(log10(star%surface_luminosity_lsun_start) - &
           log10(luminosity(num_points)))

! The factor f that brings both changes down to the triangle's size.
      dt_factor = teffl_change/star%ctrl%tri_delta_teffl
      if (logl_change/star%ctrl%tri_delta_logl .gt.dt_factor) dt_factor = logl_change/star%ctrl%tri_delta_logl
! No change since the previous model: keep the previous timestep.
      if (dt_factor .eq. 0.0d0) dt_factor = 1.0d0
! Let the timestep change by at most a factor atime(itime_max_dt_frac)
! up or down.
      dt_factor_limit = star%ctrl%atime(itime_max_dt_frac)
      if (dt_factor.gt.dt_factor_limit) dt_factor = dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor = 1.0d0/dt_factor_limit
      envelope_dt = previous_timestep/dt_factor

      return
end subroutine timestep_limit_hr

end module timestep_control
