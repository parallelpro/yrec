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
