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

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! htimer
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original htimer.f; only variable names, source form, and comment
! style were updated.
!
! htimer determines the timestep for the next model.
! There are several possible means of determining the timestep; it
! chooses the most stringent one.
! First of all, the timestep may not be more than atime(13) times the
!      previous one.
! In sr ptime, the timestep based on structure changes from one
!      model to the next is computed if desired (structure_dt).
! In sr wtime, the timestep based on changes in the angular velocity
!      from one model to the next is computed (rotating models only)
!      (rotation_dt).
! In sr xtime, the timestep based on hydrogen burning is computed
!      (hydrogen_dt).
! In sr ytime, the timestep based on helium burning is computed
!      (helium_dt).
! Details of the computation of these timesteps can be found in the
!      appropriate subroutines.
! The final timestep (previous_timestep, also returned as chosen_dt)
!      is the minimum of the above timesteps.
! **note: if a criterion is inapplicable or turned off by user
!      parameters then the 'timestep' for that criterion is
!      arbitrarily set to a large number (dt_unlimited).
! If a fixed timestep is being enforced (timestep_override_active=T),
! the above criteria are short-circuited and a timestep
! (timestep_override) is enforced. Caveat user.
! If the model is being evolved to a fixed age (end_age_stop_active=T),
! the timestep will be reduced if the one used would put the model
! past the desired age (target_end_age).
!
! A negative timestep indicates that a model is being relaxed
! and not evolved. Check for negative dt and exit if found.
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
! disable_structure_dt_limits: set when this call decided to switch the
! structure-based limits off; the driver clears
! star%job%use_structure_dt_limits on seeing it (W3).
      logical, intent(out) :: disable_structure_dt_limits
      double precision :: structure_dt, rotation_dt, hydrogen_dt, helium_dt, &
           hydrogen_luminosity, envelope_dt, time_left_years
! structure_limits_active: whether the structure-based limits apply to
! this step -- star%job%use_structure_dt_limits, switched off once the
! centre passes log Tc = 7.1 with L_grav < 0.  The routine reads only
! this local after that decision and reports it through
! disable_structure_dt_limits; the driver flips the star%job% flag.
      logical :: structure_limits_active
! dt_unlimited: "no limit from this criterion" sentinel (seconds).
      double precision, parameter :: dt_unlimited = 1.0d20

      disable_structure_dt_limits = .false.
      if (previous_timestep.ge.0.0d0) then
! if user is fixing tstep, set dt to given value and exit
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
! off once log Tc passes 7.1 with L_grav < 0.
subroutine decide_structure_limits()
! mhp 9/01  turn off structure-based timestep setting above a critical
!           temperature; this is done when
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
! The limit from the change of the structure variables (ptime), unlimited
! when the structure-based limits are off.
subroutine limit_by_structure_change()
!  find timestep based on changes in structure variables from
!       one model to the next.
!  note that this returns the timestep stored in the model on the
!       first call to htimer for each kind card.
      if(structure_limits_active) then
       call timestep_limit_structure(previous_timestep,luminosity,log_pressure,log_radius,log_temperature,num_points,structure_dt)
      else
       structure_dt = dt_unlimited
      endif
!  mhp 6/90 a zero timestep in the model will fritz out the code, so
!  avoid this by setting the structure timestep arbitrarily large.
      if(structure_dt.lt.1.0d0) structure_dt=dt_unlimited
end subroutine limit_by_structure_change

! ---------------------------------------------------------------
! The limit from the change of omega (wtime), rotating models only; it also
! returns the largest fractional change of omega.
subroutine limit_by_omega_change()
!  find timestep for rotating models based on changes in angular velocity
!       from one model to the next.
!  this sr also returns the maximum change in omega from the previous
!       model to the current one, which is used in the diffusion routine
!       to determine the number of diffusion timesteps required.
!  note that this returns the timestep stored in the model on the
!       first call to htimer for each kind card.
      if(star%job%rotation_active) then
       call timestep_limit_omega(previous_timestep,num_points,omega,rotation_dt,max_domega_frac)
      else
       rotation_dt = dt_unlimited
      endif
!  mhp 6/90 a zero timestep in the model will fritz out the code, so
!  avoid this by setting the rotation timestep arbitrarily large.
      if(rotation_dt.lt.1.0d0) rotation_dt=dt_unlimited
end subroutine limit_by_omega_change

! ---------------------------------------------------------------
! The limit from hydrogen burning (xtime), for stars that burn hydrogen.
subroutine limit_by_hydrogen_burning()
!  find the timestep based on hydrogen burning.
!  determine the total luminosity due to hydrogen burning (in solar units).
      hydrogen_luminosity = luminosity_components(1) + luminosity_components(2) + &
           luminosity_components(3) + luminosity_components(4)
!  skip h-burning timestep criteria for stars without h burning.
      if(hydrogen_luminosity.gt.1.0d-34) then
       call timestep_limit_hburn(log_density,composition,luminosity,enclosed_mass, &
            log_temperature,hydrogen_luminosity,convective_core_edge_zone, &
            h_shell_midpoint_zone,num_points,hydrogen_dt)
      else
       hydrogen_dt = dt_unlimited
      endif
end subroutine limit_by_hydrogen_burning

! ---------------------------------------------------------------
! The limit from helium burning (ytime), for stars that burn helium.
subroutine limit_by_helium_burning()
!  find the timestep based on helium burning.
!  skip for stars without he burning.
      if(luminosity_components(5).gt.1.0d-34) then
       call timestep_limit_heburn(composition,log_density,luminosity, &
            enclosed_mass,log_temperature,convective_core_edge_zone, &
            num_points,helium_dt,h_shell_zone_begin)
      else
       helium_dt = dt_unlimited
      endif
end subroutine limit_by_helium_burning

! ---------------------------------------------------------------
! The limit from the motion in the HR diagram (the envelope triangle), when
! it is on and the structure-based limits are off.
subroutine limit_by_hr_motion()
!  04/14 jvs added timestep governor based on the size of the envelope
!  triangle. timesteps are restricted such that the model does not move
!  more than tri_delta_logl or tri_delta_teffl
      if(star%ctrl%use_envelope_triangle_dt .and. .not. structure_limits_active) then
       call timestep_limit_hr(previous_timestep,luminosity,log_teff,num_points,envelope_dt)
      else
       envelope_dt = dt_unlimited
      endif
end subroutine limit_by_hr_motion

! ---------------------------------------------------------------
! The timestep: the smallest of the limits and, for a non-zero previous
! timestep, of its allowed growth.
subroutine take_tightest_limit()
!     limit increase in time step from one model to the next for models
!     with a non-zero timestep.
      if(previous_timestep.gt.1.0d0) then
         previous_timestep = star%ctrl%atime(itime_max_dt_frac)*previous_timestep
!  now set the timestep to be the minimum of the entropy based timestep,
!  the nuclear burning timesteps, and the previous timestep*atime(13).
!  04/14 jvs added envelope_dt
         previous_timestep = min(previous_timestep,hydrogen_dt,helium_dt,rotation_dt,structure_dt,envelope_dt)
      else
!  set the timestep to be the minimum of the the above timesteps,
!  but neglect the previous timestep.
!  04/14 jvs added envelope_dt
         previous_timestep = min(hydrogen_dt,helium_dt,rotation_dt,structure_dt,envelope_dt)
      endif
      chosen_dt = previous_timestep
end subroutine take_tightest_limit

! ---------------------------------------------------------------
! When evolving to a given age, shorten the timestep so the model does not
! pass it.
subroutine clip_to_end_age()
!     mhp 10/24 flag includes other stop conditions,
!  so only do this if there is a true stop age.
!     if evolving to a given age, ensure age not exceeded
!      IF(LENDAG(NK)) THEN
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

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! ptime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ptime.f; only variable names, source form, and comment
! style were updated.
!
! This subroutine limits the timestep based on the requirement that the
! changes from one model to the next in T, P, R, and L be less than the
! user parameters atime(itime_dt)-atime(itime_dl) at any model point. This criterion
! dominates in the pre-main sequence, and this subroutine must be used
! (use_structure_dt via LPTIME=T) to perform pre-main-sequence evolution.
subroutine timestep_limit_structure(previous_timestep, luminosity, log_pressure, log_radius, &
     log_temperature, num_points, struct_dt)

      use star_info_lib, only: star, json
      use controls_buffer, only: itime_dl, itime_dp, itime_dr, itime_dt, itime_max_dt_frac
      implicit none

! previous_timestep: previous model timestep.
! luminosity: run of luminosity (as a function of mass) in the current
!             model. Confirmed linear (not logged): main.f updates it
!             via HL(I)=HL(I)+HL(I)*rate*dt, an Euler step only
!             consistent with a linear quantity.
! log_pressure: run of log10(pressure) in the current model. Confirmed
!               logged: ptime.f's own comments describe atime(itime_dt)-(itime_dr)
!               as "MAX DELTA LOG T/P/R".
! log_radius: run of log10(radius) in the current model.
! log_temperature: run of log10(temperature) in the current model.
! num_points: number of points in the model.
! struct_dt: timestep for requested change in P, T, R, L.
      double precision, intent(in) :: previous_timestep
      double precision, intent(in) :: luminosity(json), log_pressure(json), &
           log_radius(json), log_temperature(json)
      integer, intent(in) :: num_points
      double precision, intent(out) :: struct_dt
! max_change: maximum absolute time differences for each quantity
! (1 = P, 2 = T, 3 = R, 4 = L).
      double precision :: max_change(4)
      integer :: i
      double precision :: test_p, test_t, test_r, test_l, dt_factor, &
           dt_factor_limit

! find maximum absolute time differences for each quantity
! pressure
      max_change(1)=abs(star%logP_start(1)-log_pressure(1))
! temperature
      max_change(2)=abs(star%logT_start(1)-log_temperature(1))
! radius
      max_change(3)=abs(star%logR_start(1)-log_radius(1))
! luminosity
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
! now actually limit the timestep by a factor that reduces the
! time changes in all quantities to the ps values or less
      dt_factor = max_change(1)/star%ctrl%atime(itime_dp)
      if (max_change(2)/star%ctrl%atime(itime_dt).gt.dt_factor) dt_factor=max_change(2)/star%ctrl%atime(itime_dt)
      if (max_change(3)/star%ctrl%atime(itime_dr).gt.dt_factor) dt_factor=max_change(3)/star%ctrl%atime(itime_dr)
      if (max_change(4)/(star%ctrl%atime(itime_dl)*2.3026d0).gt.dt_factor) then
        dt_factor=max_change(4)/(star%ctrl%atime(itime_dl)*2.3026d0)
      endif
! if no change from previous model,set struct_dt to timestep
! stored in the previous model.
      if (dt_factor.eq.0.d0) dt_factor=1.0d0
! restrict change in timestep to no more than a factor of atime(itime_max_dt_frac)
! (the global timestep limiter) up or down.
      dt_factor_limit = star%ctrl%atime(itime_max_dt_frac)
      if (dt_factor.gt.dt_factor_limit) dt_factor=dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor=1.0d0/dt_factor_limit
      struct_dt = previous_timestep/dt_factor

      return
end subroutine timestep_limit_structure

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! wtime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wtime.f; only variable names, source form, and comment
! style were updated.
!
! For rotating models, limit the timestep based on max change in
! omega. The user parameter max_domega_global (originally DTWIND)
! governs the maximum change allowed between models.
subroutine timestep_limit_omega(previous_timestep, num_points, omega, rotation_dt, &
     max_domega_frac)

      use star_info_lib, only: star, json
      use controls_buffer, only: itime_max_dt_frac
      implicit none

! previous_timestep: previous model timestep.
! num_points: number of model points.
! omega: run of angular velocity in the current model.
! rotation_dt: rotation-based timestep.
! max_domega_frac: maximum fractional change in omega.
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
! if no change from previous model,set rotation_dt to timestep
! stored in the previous model.
      if (dt_factor.eq.0.d0)then
          rotation_dt = dt_unlimited
          return
      endif
! restrict change in timestep to no more than a factor of atime(itime_max_dt_frac)
! (the global timestep limiter) up or down.
      dt_factor_limit = star%ctrl%atime(itime_max_dt_frac)
      if (dt_factor.gt.dt_factor_limit) dt_factor=dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor=1.0d0/dt_factor_limit
      rotation_dt = previous_timestep/dt_factor
      return
end subroutine timestep_limit_omega

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! xtime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original xtime.f; only variable names, source form, and comment
! style were updated.
!
! This subroutine finds the timestep based on hydrogen burning.
! For stars with central X > atime(itime_core_min) the timestep is the time needed
! to burn the minimum of atime(itime_dx_core_tot) of X at the center or the fraction
! atime(itime_dx_core_frac) of the central X.
! Stars with central X < atime(itime_core_min) are considered by the program to
! have a hydrogen shell burning source. The timestep is the minimum of
! the time required to burn atime(itime_dx_shell) of X at the shell midpoint or to
! burn the mass fraction atime(itime_dx_total) of X in the entire star.
! dX/dt at the shell midpoint comes from the network's reactions
! (compute_abundance_rates) on the zone's rates and its start-of-step
! abundances. (2026-10: it came from the equilibrium burn, with He3 and
! the CN cycle in equilibrium, removed with the level-1 burn.)
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
! the composition rate vector and dX_i/dt / A_i of the network species
! at the shell midpoint
      double precision :: rate_vec(num_rate_rows), dxdt_per_a(n_network)
! delta_x: hydrogen mass-fraction change budget for the current
! branch (core-burning limit in the first branch, whole-star shell-
! burning limit in the second -- a single reused scratch variable in
! the original, exactly as here).
      double precision :: delta_x
! local_log_density/local_log_temperature: log10(rho)/log10(T) at the
! shell midpoint, fed to rates().
      double precision :: local_log_density, local_log_temperature
      integer :: zone_end
! energy released per gram of hydrogen burned (erg/g)
      double precision, parameter :: h_burn_energy_per_gram = 6.00d18
      double precision :: dx_dt, shell_dt_x_depletion

!  h-core burning time criterion
! **note that for a convective core, the timestep is based on the time
!   needed to burn the given fraction of hydrogen in the core and not
!   just in the central shell.
      if(composition(i_h1,1).ge.star%ctrl%atime(itime_core_min)) then
       delta_x = min(star%ctrl%atime(itime_dx_core_tot),star%ctrl%atime(itime_dx_core_frac)*composition(i_h1,convective_core_edge_zone))
       hydrogen_dt =(h_burn_energy_per_gram/star%solar_luminosity_cgs)* &
            (enclosed_mass(convective_core_edge_zone)/luminosity(convective_core_edge_zone))*delta_x
       return
      endif
!  h-shell burning criterion
!  limit total mass of hydrogen burned.
      delta_x = star%ctrl%atime(itime_dx_total)*composition(i_h1,num_points)*(star%solar_mass_cgs/star%solar_luminosity_cgs)
      hydrogen_dt = h_burn_energy_per_gram*delta_x/hydrogen_luminosity
!  limit x-depletion at shell mid-point: dX/dt there from the network
!  (A = 1 for hydrogen)
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

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! ytime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ytime.f; only variable names, source form, and comment
! style were updated -- except the He-shell branch (2026, see below).
!
! Stars with a helium luminosity have a timestep limit based on the
! time required to burn atime(itime_dy_core_tot) of Y at the maximum in temperature or
! the fraction atime(itime_dy_core_frac) of Y at this point.
! Test for helium flash burning; high central density is taken as a
! sign that the model is a giant undergoing a He flash rather than a
! horizontal branch star.
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
! helium_dt: assigned on every path (2026 -- the He-shell branch no
! longer reads the incoming value; see the note there).
      double precision, intent(out) :: helium_dt
      integer, intent(in) :: h_shell_zone_begin
      double precision :: max_temp, local_log_density, local_log_temperature
      integer :: max_temp_zone, zone_idx, shell_zone
! the hottest zone's composition (deuterium counted only with the
! extended composition, as it always has been here)
      double precision :: zone_xa(15)
! only the helium-burning term (burn%triple_alpha_energy_gen) of the
! shell's energy record is used here.
      type(nuclear_energy_result) :: burn
      double precision :: helium_energy_gen
      double precision :: core_helium_fraction
! energy released per gram of helium burned (erg/g)
      double precision, parameter :: he_burn_energy_per_gram = 5.85d17
! dt_unlimited: "no limit from this criterion" sentinel (seconds).
      double precision, parameter :: dt_unlimited = 1.0d20

      if (log_density(1).ge.5.0d0) then
!     search for temperature maximum
       max_temp = 0.0d0
       max_temp_zone = 1
       do zone_idx = 1,num_points
          if (log_temperature(zone_idx).gt.max_temp) then
             max_temp_zone = zone_idx
             max_temp = log_temperature(zone_idx)
          endif
       end do
!     calculate helium burning rate at tmax
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
!  set limits on he burning for non-helium flash stars
!  core_helium_fraction: computed as 1 - Z(core); in a He-burning core
!  X is normally negligible there, so this approximates the core Y.
       core_helium_fraction = 1d0 - composition(i_metals,convective_core_edge_zone)
       if(core_helium_fraction.ge.star%ctrl%atime(itime_core_min)) then
            helium_dt = min(star%ctrl%atime(itime_dy_core_tot),star%ctrl%atime(itime_dy_core_frac)*core_helium_fraction)
            helium_dt = (he_burn_energy_per_gram/star%solar_luminosity_cgs)*helium_dt* &
               (enclosed_mass(convective_core_edge_zone)/luminosity(convective_core_edge_zone))
         else
! 2026 (bugsweep sec-10/11): core helium exhausted -- the He-shell
! branch. The original (ytime.f) read DELTSY on the right-hand side
! here, i.e. the PREVIOUS call's timestep (kept alive only by the
! caller's blanket SAVE), and multiplied it by M/L again: a
! dimensionless-nonsense recurrence that gave ~0 s on the first
! model after exhaustion and, once the SAVE was gone, exactly 0 ->
! dt = 0 -> NaN (the run_from_zahb_to_tahb failure at Y_c < atime(itime_core_min)).
! Compute a fresh limit from the He-shell controls instead, mirroring
! the H-shell branch of timestep_limit_hburn: atime(itime_dy_total) (time_dy_total,
! Msun of He burned per step) and atime(itime_dy_shell) (time_dy_shell, fraction
! of Y burned at the shell) -- both were mapped by map_user_inputs but
! never consumed until now. The reference point is the shell just
! below the H-burning shell (the He-rich core), as before.
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

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! entime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original entime.f; only variable names, source form, and comment
! style were updated.
!
! This subroutine limits the timestep based on the requirement that the
! changes from one model to the next in Teff and L be less than the
! size of the envelope triangle (tri_delta_teffl and tri_delta_logl).
! This prevents extrapolation outside of the triangle for small
! envelope triangles.
subroutine timestep_limit_hr(previous_timestep, luminosity, log_teff, &
     num_points, envelope_dt)

      use star_info_lib, only: star, json
      use controls_buffer, only: itime_max_dt_frac
      use math_lib
      implicit none

! previous_timestep: previous model timestep.
! luminosity: run of luminosity (as a function of mass) in the current
!             model.
! log_teff: log10(effective temperature) of the current model
!           (confirmed logged: main.f computes 10.0d0**TEFFL).
! num_points: number of points in the model.
! envelope_dt: timestep for requested change in Teff, L.
      double precision, intent(in) :: previous_timestep
      double precision, intent(in) :: luminosity(json)
      double precision, intent(in) :: log_teff
      integer, intent(in) :: num_points
      double precision, intent(out) :: envelope_dt
      double precision :: teffl_change, logl_change, dt_factor, &
           dt_factor_limit

! find absolute changes in log Teff and log L since the previous model
      teffl_change = abs(star%log_Teff_start - log_teff)
      logl_change = abs(log10(star%luminosity_lsun_start(star%nz_start)) - &
           log10(luminosity(num_points)))

! now actually limit the timestep by a factor that reduces the
! time changes in all quantities to the triangle values or less
      dt_factor = teffl_change/star%ctrl%tri_delta_teffl
      if (logl_change/star%ctrl%tri_delta_logl .gt.dt_factor) dt_factor = logl_change/star%ctrl%tri_delta_logl
! if no change from previous model, set envelope_dt to timestep
! stored in the previous model.
      if (dt_factor .eq. 0.0d0) dt_factor = 1.0d0
! use atime(itime_max_dt_frac) as the global factor for limiting timestep changes
      dt_factor_limit = star%ctrl%atime(itime_max_dt_frac)
      if (dt_factor.gt.dt_factor_limit) dt_factor = dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor = 1.0d0/dt_factor_limit
      envelope_dt = previous_timestep/dt_factor

      return
end subroutine timestep_limit_hr

end module timestep_control
