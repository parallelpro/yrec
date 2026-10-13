!----------------------------------------------------------------------
! burn_step (formerly burn_lib)
!----------------------------------------------------------------------
! The star-layer side of nuclear burning: the code that reads the model,
! calls the net/ physics and stores the results in star. The physics
! itself lives in net/: the rates (net_lib's compute_reaction_rates), the
! energy generation (nuclear_energy), the network and its implicit and
! equilibrium burns (network), and the light elements (light_elements).
! Here:
!   configure_network     net/'s settings, from the run's controls
!   surface_accretion_now the step's surface accretion, for the
!                         deuterium burn
!   compute_zone_energy   one shell's energy generation
!   burn_regions          the regions a burn visits, in order
!   compute_burn_rates, network_rates_of_model   the rate vectors
!   burn_network, burn_deuterium, burn_network_region
!                         the composition step's burns
!   store_start/end_li_be_rates, burn_li_be, average_li_be_over_cz
!                         Li6, Li7 and Be9 over a step
module burn_step
      use net_lib, only: num_rate_rows, nuclear_energy_result, &
           surface_accretion
      use convection_geometry, only: locate_surface_mixing_base
! exp and log: crmath's under USE_CRMATH
      use math_lib
      use star_info_lib, only: star_info
      implicit none
      private
      public :: configure_network, surface_accretion_now, &
           compute_zone_energy, burn_regions, network_rates_of_model, &
           compute_burn_rates, burn_network, burn_deuterium, burn_network_region, store_start_li_be_rates, &
           store_end_li_be_rates, burn_li_be, average_li_be_over_cz

contains



!----------------------------------------------------------------------
! configure_network
!----------------------------------------------------------------------
! Fills net/'s settings (net_lib's net_settings_t) from the run's
! controls and rate factors, and installs them (set_net_settings). run_yrec
! calls it at the start of every run, after the Monte Carlo rate
! factors are applied; net/ reads no other configuration.
!   star reads: c12_alpha_scale, cross_section_scale,
!     ctrl%be9_palpha_rate_scale, ctrl%be9_pd_rate_scale,
!     ctrl%be9_pg_rate_scale, ctrl%calc_neutrinos, ctrl%li6_rate_scale,
!     ctrl%li7_rate_scale, ctrl%nuclear_logT_cutoffs,
!     ctrl%use_itoh_neutrino_loss, ctrl%weak_screening_threshold,
!     o16_gamma_scale, qqs0ee_scale, qs0e_scale
!   star writes: none
subroutine configure_network(s)
      use net_lib, only: net_settings_t, set_net_settings
      type(star_info), intent(inout) :: s
      type(net_settings_t) :: settings

      settings%logT_cutoffs = s%ctrl%nuclear_logT_cutoffs
      settings%weak_screening_threshold = s%ctrl%weak_screening_threshold
      settings%cross_section_scale = s%cross_section_scale
      settings%qs0e_scale = s%qs0e_scale
      settings%qqs0ee_scale = s%qqs0ee_scale
      settings%o16_gamma_scale = s%o16_gamma_scale
      settings%c12_alpha_scale = s%c12_alpha_scale
      settings%li6_rate_scale = s%ctrl%li6_rate_scale
      settings%li7_rate_scale = s%ctrl%li7_rate_scale
      settings%be9_pg_rate_scale = s%ctrl%be9_pg_rate_scale
      settings%be9_pd_rate_scale = s%ctrl%be9_pd_rate_scale
      settings%be9_palpha_rate_scale = s%ctrl%be9_palpha_rate_scale
      settings%calc_neutrinos = s%ctrl%calc_neutrinos
      settings%use_itoh_neutrino_loss = s%ctrl%use_itoh_neutrino_loss
      call set_net_settings(settings)
end subroutine configure_network

!----------------------------------------------------------------------
! surface_accretion_now
!----------------------------------------------------------------------
! The surface accretion of a step, as burn_deuterium_region takes it:
! whether accretion is on, whether matter flows in (mass_accretion_rate
! > 0), the deuterium fraction of the accreted gas, and the mass accreted
! over the step (mass, g).
!   star reads: ctrl%accreted_composition, ctrl%mass_accretion_rate,
!     job%use_mass_accretion
!   star writes: none
function surface_accretion_now(s, mass) result(accretion)
      use star_info_lib, only: i_h2
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: mass
      type(surface_accretion) :: accretion

      accretion%on = s%job%use_mass_accretion
      accretion%inflow = s%ctrl%mass_accretion_rate.gt.0.0d0
      accretion%deuterium_fraction = s%ctrl%accreted_composition(i_h2)
      accretion%mass = mass
end function surface_accretion_now

!----------------------------------------------------------------------
! compute_zone_energy (formerly engeb)
!----------------------------------------------------------------------
! The energy generation of one shell: net's compute_nuclear_energy, with
! the two inputs that depend on where the shell sits in the model filled
! in: whether it lies in the surface convection zone (shell_index >=
! star%jcz), and the convective turnover time; together they cap
! deuterium burning there. Everything comes back in e; nothing is
! stored.
!   star reads: convective_turnover_timescale, jcz
!   star writes: none
subroutine compute_zone_energy(s, log_density, log_temperature, xa, &
     shell_index, e)
      use net_lib, only: compute_nuclear_energy
      type(star_info), intent(inout) :: s
! xa: the shell's composition, star%xa layout
      double precision, intent(in) :: log_density, log_temperature, xa(:)
      integer, intent(in) :: shell_index
      type(nuclear_energy_result), intent(out) :: e

      call compute_nuclear_energy(log_density, log_temperature, xa, &
           shell_index.ge.s%jcz, s%convective_turnover_timescale, e)
end subroutine compute_zone_energy


!----------------------------------------------------------------------
! average_li_be_over_cz
!----------------------------------------------------------------------
! The Li6, Li7 and Be9 abundances (x_start, x_end) and the natural log of
! their burning rates (log_rate_start, log_rate_end), mass-averaged over
! the surface convection zone, for burn_li_be (and the commented-out
! burn_li_be_midpoint, rotation/burn_li_be_midpoint.f90):
!   start of the step: from cz_base_zone_old up, with the start-of-step
!     rates;
!   end of the step: from cz_base_zone up, with the end-of-step rates and
!     the base zone k = cz_base_zone adjusted to the true edge of the
!     region (when 1 < k < num_zones): its mass changed by
!     cz_base_frac (m(k) - m(k-1)), and its rate moved toward the zone
!     below's,
!       rate(k) + (cz_base_frac/2) (rate(k-1) - rate(k)).
!     cz_mass_end is the region's mass then.
! (Deviations from the legacy code, 2026-10: the base-zone adjustment
! stays local, where the legacy code kept the adjusted rate in
! star%rate_*(cz_base_zone) and so carried it into the next rotation
! sub-step's start rates; and a legacy branch that reused the previous
! step's rate when star%log_rate_li6_prev > 0, which the natural log of
! these rates per second never is, is removed.)
!   star reads: rate_be9, rate_be9_start, rate_li6, rate_li6_start,
!     rate_li7, rate_li7_start
!   star writes: none
subroutine average_li_be_over_cz(s, composition, mass_coordinate, shell_mass, &
     cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, &
     x_start, log_rate_start, x_end, cz_mass_end, log_rate_end)
      use star_info_lib, only: json, i_li6, i_li7, i_be9
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(in) :: shell_mass(json)
      integer, intent(in) :: cz_base_zone, cz_base_zone_old, num_zones
      double precision, intent(in) :: cz_base_frac
      double precision, intent(out) :: x_start(3), log_rate_start(3), &
           x_end(3), cz_mass_end, log_rate_end(3)

      integer, parameter :: species(3) = [i_li6, i_li7, i_be9]
      double precision :: rate_sum(3), rate(3), cz_mass_start, mass
      integer :: zone_idx

! The start of the step, at the start-of-step depth.
      x_start = 0.0d0
      rate_sum = 0.0d0
      cz_mass_start = 0.0d0
      do zone_idx = cz_base_zone_old,num_zones
         x_start = x_start + composition(species,zone_idx)*shell_mass(zone_idx)
         rate_sum = rate_sum + [s%rate_li6_start(zone_idx), &
              s%rate_li7_start(zone_idx), s%rate_be9_start(zone_idx)]* &
              shell_mass(zone_idx)
         cz_mass_start = cz_mass_start + shell_mass(zone_idx)
      end do
      x_start = x_start/cz_mass_start
      log_rate_start = log(rate_sum/cz_mass_start)
! The end of the step, at the end-of-step depth, with the base zone
! adjusted to the true edge of the region (cz_base_frac, from
! locate_surface_mixing_base).
      x_end = 0.0d0
      rate_sum = 0.0d0
      cz_mass_end = 0.0d0
      do zone_idx = cz_base_zone,num_zones
         mass = shell_mass(zone_idx)
         rate = [s%rate_li6(zone_idx), s%rate_li7(zone_idx), &
              s%rate_be9(zone_idx)]
         if (zone_idx.eq.cz_base_zone .and. cz_base_zone.gt.1 .and. &
              cz_base_zone.lt.num_zones) then
            mass = mass + cz_base_frac* &
                 (mass_coordinate(cz_base_zone)-mass_coordinate(cz_base_zone-1))
            rate = rate + 0.5d0*cz_base_frac*([s%rate_li6(cz_base_zone-1), &
                 s%rate_li7(cz_base_zone-1), s%rate_be9(cz_base_zone-1)] - rate)
         end if
         rate_sum = rate_sum + rate*mass
         x_end = x_end + composition(species,zone_idx)*mass
         cz_mass_end = cz_mass_end + mass
      end do
      x_end = x_end/cz_mass_end
      log_rate_end = log(rate_sum/cz_mass_end)
end subroutine average_li_be_over_cz



!----------------------------------------------------------------------
! burn_li_be
!----------------------------------------------------------------------
! Li6, Li7 and Be9 burning over a step (formerly liburn). Each species
! decays at its proton-capture rate lambda, interpolated in log between
! the start- and end-of-step rates (store_start_li_be_rates,
! store_end_li_be_rates), lambda(t) = lambda_0 (lambda_1/lambda_0)^(t/dt),
! and survives as X_end = X_start exp(-integral of lambda dt), in closed
! form (light_elements' surviving_fraction; it replaced the legacy
! sub-stepping with rational extrapolation to zero step, which converged
! on the same integral to 1e-6). In order:
!   1. the base of the surface mixed region (convection zone plus
!      envelope overshoot) at the start and end of the step
!      (locate_surface_mixing_base);
!   2. radiative shells below both bases, each on its own rates (shells
!      hotter than log T = 7 are cleared);
!   3. the convection zone as one unit, on its mass-averaged rates
!      (average_li_be_over_cz), starting from the abundance mixed over
!      the deeper zone when it deepened, and taking in accreted gas
!      (accreted_survival) when accretion is on;
!   4. shells the convection zone left during the step, on a rate
!      weighted between their radiative rate and the zone's by the time
!      spent in each.
! Abundances below 1e-24 are set to zero.
! The model arrays are arguments because the two callers burn different
! models: burn_light_elements the converged end-of-step model,
! evolve_angular_momentum each rotation sub-step's midpoint (rot_scr).
!   star reads: accreted_mass, ctrl%accreted_composition,
!     ctrl%mass_accretion_rate, job%use_mass_accretion, rate_be9,
!     rate_be9_start, rate_li6, rate_li6_start, rate_li7, rate_li7_start
!   star writes: none
subroutine burn_li_be(s, timestep, composition, radius, mass_coordinate, &
     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
      use star_info_lib, only: json, i_li6, i_li7, i_be9
      use net_lib, only: surviving_fraction, accreted_survival, &
           linear_in_time, linear_in_log
      type(star_info), intent(inout) :: s

! timestep: s
      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: radius(json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones

      integer, parameter :: species(3) = [i_li6, i_li7, i_be9]
      integer :: zone_idx, isp, cz_base_zone, cz_base_zone_old, min_zone, max_zone
      double precision :: cz_base_frac
      double precision :: rate_start(3), rate_end(3), rate
      double precision :: x_cz_start(3), x_cz_end(3), log_rate_cz_start(3), &
           log_rate_cz_end(3), cz_mass_end
      logical :: accreting
      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
           radiative_frac

! 1. The base of the surface mixed region at the start and end of the
! step. The lithium burned in a surface convection zone depends
! sensitively on the temperature at its base, so locating it accurately
! matters. min_zone is the deeper of the two bases, max_zone the
! shallower.
      call locate_surface_mixing_base(s, radius, env_cz_zone, env_cz_zone_old, &
           num_zones, cz_base_zone, cz_base_zone_old, cz_base_frac)
      min_zone = min(cz_base_zone,cz_base_zone_old)
      max_zone = max(cz_base_zone,cz_base_zone_old)

! 2. The radiative interior, up to the deeper of the two bases. The loop
! stops at the first shell too cool to burn (rate_be9 <= 1e-32 at the
! start or end) and skips shells already depleted.
      do zone_idx = 1,min_zone-1
         if(s%rate_be9(zone_idx).le.1.0d-32 .or. s%rate_be9_start(zone_idx).le.1.0d-32)exit
         if(all(composition(species,zone_idx).lt.1.0d-24))cycle
         if(log_temperature(zone_idx).gt.7.0d0)then
            composition(species,zone_idx) = 0.0d0
            cycle
         endif
         rate_start = [s%rate_li6_start(zone_idx), s%rate_li7_start(zone_idx), &
              s%rate_be9_start(zone_idx)]
         rate_end = [s%rate_li6(zone_idx), s%rate_li7(zone_idx), &
              s%rate_be9(zone_idx)]
         do isp = 1,3
            composition(species(isp),zone_idx) = composition(species(isp),zone_idx)* &
                 surviving_fraction(rate_start(isp), rate_end(isp), timestep, linear_in_log)
            if(composition(species(isp),zone_idx).lt.1.0d-24) &
                 composition(species(isp),zone_idx) = 0.0d0
         end do
      end do

! 3. The convection zone as one unit. When its base is too cool to burn
! (rate_be9 <= 1e-32 at the start or end), nothing more burns: return.
      if (s%rate_be9_start(cz_base_zone_old).le.1.0d-32.or.s%rate_be9(cz_base_zone).le.1.0d-32) then
         return
      end if
      call average_li_be_over_cz(s, composition, mass_coordinate, shell_mass, &
           cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, &
           x_cz_start, log_rate_cz_start, x_cz_end, cz_mass_end, log_rate_cz_end)
! When the zone deepened, start from the abundance mixed over the deeper
! zone.
      if(cz_base_zone.lt.cz_base_zone_old) x_cz_start = x_cz_end
! The zone's abundance survives the step on its mass-averaged rates.
! With mass accretion (MHP 05/02), the accreted gas burns from its
! arrival on (accreted_survival, S_acc) and is mixed in by mass:
!   x_cz = (x_cz M_cz + x_acc S_acc M_acc) / (M_cz + M_acc).
! This assumes that all the accreted matter lands in a convection zone
! of the end-of-step depth: not strictly true, but not a bad
! approximation either.
      accreting = s%job%use_mass_accretion .and. &
           s%ctrl%mass_accretion_rate.gt.0.0d0
      do isp = 1,3
         x_cz_end(isp) = x_cz_start(isp)*surviving_fraction(exp(log_rate_cz_start(isp)), &
              exp(log_rate_cz_end(isp)), timestep, linear_in_log)
         if (accreting) x_cz_end(isp) = (x_cz_end(isp)*cz_mass_end + &
              s%ctrl%accreted_composition(species(isp))* &
              accreted_survival(exp(log_rate_cz_start(isp)), exp(log_rate_cz_end(isp)), &
              timestep, linear_in_log)*s%accreted_mass)/ &
              (cz_mass_end + s%accreted_mass)
         if(x_cz_end(isp).lt.1.0d-24) x_cz_end(isp) = 0.0d0
      end do
! The new abundance is set from the shallower base (max_zone) out; when
! the zone deepened, the shells between the two bases are left as they
! are in this call.
      do zone_idx = max_zone,num_zones
         composition(species,zone_idx) = x_cz_end
      end do

! 4. Shells that began convective and ended radiative (a retreating
! convection zone).
      if (cz_base_zone.le.cz_base_zone_old) then
         return
      end if
! The base moved out in mass from mass_coord_beg (the midpoint of the old
! base zone and the one below) to mass_coord_end (likewise for the new).
! Assume that the time a shell spent in the convection zone is
! proportional to its place in mass between the two: a shell 1/3 of the
! way from the old base to the new spent 1/3 of the step in the zone and
! 2/3 out of it. The radiative fraction of the step is then
!   radiative_frac = (mass_coord_end - m)/(mass_coord_end - mass_coord_beg),
! 1 at the old base and 0 at the new one. (The sign was fixed in the 2026
! bugsweep.)
      if(cz_base_zone_old.gt.1)then
         mass_coord_beg = 0.5d0*(mass_coordinate(cz_base_zone_old)+mass_coordinate(cz_base_zone_old-1))
      else
         mass_coord_beg = 0.0d0
      endif
      mass_coord_end = 0.5d0*(mass_coordinate(cz_base_zone)+mass_coordinate(cz_base_zone-1))
      delta_mass = mass_coord_end - mass_coord_beg
      do zone_idx = cz_base_zone_old,cz_base_zone-1
! Stop at the first shell too cool to burn.
         if(s%rate_be9(zone_idx).le.1.0d-32)exit
         radiative_frac = (mass_coord_end - mass_coordinate(zone_idx))/delta_mass
         rate_end = [s%rate_li6(zone_idx), s%rate_li7(zone_idx), &
              s%rate_be9(zone_idx)]
         do isp = 1,3
! The rate, weighted in the log between the shell's own end-of-step rate
! and the convection zone's start-of-step rate, and held over the step:
!   ln rate = radiative_frac ln rate_end + (1 - radiative_frac) ln rate_cz.
            rate = exp(radiative_frac*log(rate_end(isp)) + &
                 (1.0d0-radiative_frac)*log_rate_cz_start(isp))
            composition(species(isp),zone_idx) = composition(species(isp),zone_idx)* &
                 surviving_fraction(rate, rate, timestep, linear_in_time)
            if(composition(species(isp),zone_idx).lt.1.0d-24) &
                 composition(species(isp),zone_idx) = 0.0d0
         end do
      end do
end subroutine burn_li_be




!----------------------------------------------------------------------
! store_li_be_rates (formerly lirate88)
!----------------------------------------------------------------------
! The Li6, Li7 and Be9 proton-capture rates per s of every shell, for
! burn_li_be: rate = rho X lambda(T9), with lambda from net's
! compute_li_be_rates (Caughlan & Fowler 1988). With start_of_step, of
! the saved start-of-step model (star%logRho_start, logT_start, xa_start)
! into star%rate_*_start; otherwise of the model in the arguments into
! star%rate_*. The loop stops at the first shell cooler than log T = 6 in
! both the argument model and the start-of-step one; from there out both
! sets are zeroed. (Fixed 2026-10: the legacy callers passed the 1/2 flag
! the wrong way round, so both sets came from the start-of-step model.)
!   star reads: logRho_start, logT_start, xa_start
!   star writes: rate_be9, rate_be9_start, rate_li6, rate_li6_start,
!     rate_li7, rate_li7_start
subroutine store_li_be_rates(s, composition, log_density, log_temperature, &
     num_zones, start_of_step)
      use star_info_lib, only: json, i_h1
      use physical_constants, only: ln10
      use net_lib, only: compute_li_be_rates
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_density(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: num_zones
      logical, intent(in) :: start_of_step

      double precision, parameter :: min_log_temperature = 6.0d0
      integer :: zone_idx
      double precision :: rhox, t9, li6_rate, li7_rate, be9_rate

      do zone_idx = 1,num_zones
         if (log_temperature(zone_idx) < min_log_temperature .and. &
              s%logT_start(zone_idx) < min_log_temperature) exit
         if (start_of_step) then
            rhox = exp(ln10*s%logRho_start(zone_idx))*s%xa_start(i_h1,zone_idx)
            t9 = exp(ln10*(s%logT_start(zone_idx)-9.0d0))
         else
            rhox = exp(ln10*log_density(zone_idx))*composition(i_h1,zone_idx)
            t9 = exp(ln10*(log_temperature(zone_idx)-9.0d0))
         end if
         call compute_li_be_rates(t9, li6_rate, li7_rate, be9_rate)
         if (start_of_step) then
            s%rate_li6_start(zone_idx) = li6_rate*rhox
            s%rate_li7_start(zone_idx) = li7_rate*rhox
            s%rate_be9_start(zone_idx) = be9_rate*rhox
         else
            s%rate_li6(zone_idx) = li6_rate*rhox
            s%rate_li7(zone_idx) = li7_rate*rhox
            s%rate_be9(zone_idx) = be9_rate*rhox
         end if
      end do
      s%rate_li6(zone_idx:num_zones) = 0.0d0
      s%rate_li7(zone_idx:num_zones) = 0.0d0
      s%rate_be9(zone_idx:num_zones) = 0.0d0
      s%rate_li6_start(zone_idx:num_zones) = 0.0d0
      s%rate_li7_start(zone_idx:num_zones) = 0.0d0
      s%rate_be9_start(zone_idx:num_zones) = 0.0d0
end subroutine store_li_be_rates

!----------------------------------------------------------------------
! store_start_li_be_rates / store_end_li_be_rates
!----------------------------------------------------------------------
! The two uses of store_li_be_rates. Start of step: the rates of the
! saved start-of-step model into star%rate_*_start. End of step: the
! rates of the model passed in (the converged model, or a rotation
! sub-step's midpoint model) into star%rate_*.
!   star reads: logRho, logT, nz, xa
!   star writes: none
subroutine store_start_li_be_rates(s)
      type(star_info), intent(inout) :: s
      call store_li_be_rates(s, s%xa, s%logRho, s%logT, s%nz, .true.)
end subroutine store_start_li_be_rates

!   star reads: none
!   star writes: none
subroutine store_end_li_be_rates(s, composition, log_density, log_temperature, num_zones)
      use star_info_lib, only: json
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_temperature(json)
      integer, intent(in) :: num_zones
      call store_li_be_rates(s, composition, log_density, log_temperature, &
           num_zones, .false.)
end subroutine store_end_li_be_rates


!----------------------------------------------------------------------
! burn_network_region (formerly kemcom)
!----------------------------------------------------------------------
! New abundances at the end of the timestep for one radiative shell
! (zone_begin = zone_end) or one convection zone (zone_begin < zone_end,
! homogenized), by an implicit burning step of the network in
! net/network.f90. This routine does only the model bookkeeping around
! it; the reactions, the rate law and the Newton iteration live in
! network.
!
!  log_temperature     run of log10 T; shells below
!                      ctrl%nuclear_logT_cutoffs(1) are only mixed
!  rate_by_zone(r_*/rv_*,k) net_lib's rate vectors, per s per amu
!  shell_mass          the mass weights for a convection zone
!  composition         out: the new abundances (star%xa layout); the
!                      solved shell's previous neighbour is the first
!                      guess, so shells are solved in increasing order
!  timestep            the step, s
! The start-of-step abundances are star%xa_start. Z changes only once
! hydrogen is exhausted (X < 5e-7 at the top of the region), as
! 1 - X - Y - He3.
! log_temperature, rate_by_zone and composition are arguments because
! the two callers differ: burn_network passes the model's, burn_and_settle the
! rotation sub-step's midpoint temperatures (rot_scr). If the Newton
! iteration fails, ierr = 1 and the run stops.
!   star reads: ctrl%absolute_tolerance, ctrl%max_burn_iterations,
!     ctrl%min_abundance, ctrl%nuclear_logT_cutoffs,
!     ctrl%relative_tolerance, xa_start
!   star writes: none
subroutine burn_network_region(s, log_temperature, zone_begin, zone_end, &
     rate_by_zone, shell_mass, composition, timestep, ierr)

      use star_info_lib, only: json, i_h1, i_metals, n_species_basic
      use net_lib, only: burn_implicitly, n_solved, n_network, &
           network_species_slot, s_h1, s_he3, s_he4
      use io_units, only: run_log_unit
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: zone_begin, zone_end
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
! timestep: s
      double precision, intent(in) :: timestep
      integer, intent(out) :: ierr

      double precision :: region_start(n_species_basic), rate_vec(num_rate_rows)
      double precision :: x_start(n_network), x(n_network), region_mass, metals
      integer :: k, guess_zone

      ierr = 0

! The start-of-step abundances and rates of the region (mass-weighted
! averages over a convection zone).
      if (zone_begin /= zone_end) then
         region_mass = 0.0d0
         region_start = 0.0d0
         rate_vec = 0.0d0
         do k = zone_begin, zone_end
            region_mass = region_mass + shell_mass(k)
            region_start = region_start + s%xa_start(1:n_species_basic,k)*shell_mass(k)
            rate_vec = rate_vec + shell_mass(k)*rate_by_zone(:,k)
         end do
         region_start = region_start/region_mass
         rate_vec = rate_vec/region_mass
      else
         region_start = s%xa_start(1:n_species_basic,zone_begin)
         rate_vec = rate_by_zone(:,zone_begin)
      end if

! Too cool to burn (below ctrl%nuclear_logT_cutoffs(1)): only homogenize.
      if (log_temperature(zone_begin) < s%ctrl%nuclear_logT_cutoffs(1)) then
         do k = zone_begin, zone_end
            composition(1:n_species_basic,k) = region_start
         end do
         return
      end if

! The first guess: the shell just solved below (or this one, at the
! centre).
      guess_zone = max(zone_begin - 1, 1)
      x_start = region_start(network_species_slot)
      x(1:n_solved) = composition(network_species_slot(1:n_solved), guess_zone)
      x(n_solved+1:) = x_start(n_solved+1:)

      call burn_implicitly(x_start, rate_vec, timestep, &
           s%ctrl%min_abundance, s%ctrl%absolute_tolerance, &
           s%ctrl%relative_tolerance, s%ctrl%max_burn_iterations, x, ierr)
      if (ierr /= 0) then
         write (run_log_unit,'(a,i0,a,i0,a)') ' burn_network_region: no '// &
              'convergence for the new abundances of shell ', zone_begin, &
              ' in ', s%ctrl%max_burn_iterations, &
              ' iterations (ctrl%max_burn_iterations); run stopped'
         return
      end if

! Z: once hydrogen is exhausted at the top of the region (X < 5e-7),
! Z = 1 - X - Y - He3; otherwise the region's start-of-step Z.
      if (s%xa_start(i_h1,zone_end) < 5.0d-7) then
         metals = 1.0d0 - x(s_h1) - x(s_he3) - x(s_he4)
      else
         metals = region_start(i_metals)
      end if
      do k = zone_begin, zone_end
         composition(network_species_slot,k) = x
         composition(i_metals,k) = metals
      end do
end subroutine burn_network_region

! =================================================================
! The burns of the composition step (core/step/composition_step.f90).
! =================================================================

!----------------------------------------------------------------------
! burn_regions
!----------------------------------------------------------------------
! The regions a burn visits, in the order it visits them: each shell of
! each radiative region in turn, up to the first shell of that region at
! or below the nuclear cutoff ctrl%nuclear_logT_cutoffs(1), then each
! mixed region as a unit (burned whatever its temperature). Column r is
! region r's first and last zone. The temperatures and region tables
! are arguments: the composition step passes the model's, the rotation
! sub-steps their midpoint model's.
!   star reads: ctrl%nuclear_logT_cutoffs
!   star writes: none
function burn_regions(s, log_temperature, radiative_zone_bounds, &
     num_radiative_zones, mixed_zone_bounds, num_mixed_zones) result(regions)
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: log_temperature(:)
      integer, intent(in) :: radiative_zone_bounds(:,:), mixed_zone_bounds(:,:)
      integer, intent(in) :: num_radiative_zones, num_mixed_zones
      integer, allocatable :: regions(:,:)
      integer :: list(2, size(log_temperature) + num_mixed_zones)
      integer :: n, r, k

      n = 0
      do r = 1, num_radiative_zones
         do k = radiative_zone_bounds(r,1), radiative_zone_bounds(r,2)
            if (log_temperature(k) <= s%ctrl%nuclear_logT_cutoffs(1)) exit
            n = n + 1
            list(:,n) = [k, k]
         end do
      end do
      do r = 1, num_mixed_zones
         n = n + 1
         list(:,n) = mixed_zone_bounds(r,:)
      end do
      regions = list(:,1:n)
end function burn_regions

!----------------------------------------------------------------------
! network_rates_of_model
!----------------------------------------------------------------------
! The composition rate vector (net_lib's compute_network_rates) of every
! shell of a model hotter than the nuclear cutoff; zero in the cooler
! shells. The model is in the arguments: compute_burn_rates passes the
! star's, mid_timestep_model a rotation sub-step's midpoint (rot_scr).
!   star reads: ctrl%nuclear_logT_cutoffs
!   star writes: none
subroutine network_rates_of_model(s, log_density, log_temperature, &
     composition, num_zones, rate_by_zone)
      use star_info_lib, only: json
      use net_lib, only: compute_network_rates
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: log_density(json), &
           log_temperature(json), composition(15,json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: rate_by_zone(num_rate_rows, json)
      integer :: zone_idx, clear_idx

      do zone_idx = 1, num_zones
! Stop at the first shell at or below the nuclear cutoff.
         if (log_temperature(zone_idx).le.s%ctrl%nuclear_logT_cutoffs(1)) exit
         call compute_network_rates(log_density(zone_idx), log_temperature(zone_idx), &
              composition(:,zone_idx), rate_by_zone(:, zone_idx))
      end do
! Zero the rates from there out (zone_idx = num_zones + 1 if no shell was
! that cool).
      do clear_idx = zone_idx, num_zones
         rate_by_zone(:, clear_idx) = 0.0d0
      end do
end subroutine network_rates_of_model

!----------------------------------------------------------------------
! compute_burn_rates
!----------------------------------------------------------------------
! The composition rate vector of every shell of the star above the
! nuclear cutoff (network_rates_of_model) and its deuterium burning rate
! (light_elements' deuterium_rate, stored in star%deuterium_burning_rate
! and, with store_start_rates, in star%deuterium_burning_rate_start).
! The rotation sub-steps compute their own rate vectors on the
! sub-step's midpoint model (mid_timestep_model, for burn_and_settle).
!   star reads: convective_turnover_timescale,
!     ctrl%accreted_composition, ctrl%mass_accretion_rate,
!     ctrl%nuclear_logT_cutoffs, jcz, job%use_mass_accretion, logRho,
!     logT, nz, xa
!   star writes: deuterium_burning_rate, deuterium_burning_rate_start
subroutine compute_burn_rates(s, store_start_rates, rate_by_zone)
      use star_info_lib, only: json, i_h1, i_h2
      use net_lib, only: deuterium_rate
      type(star_info), intent(inout) :: s
! store_start_rates: also keep the deuterium rates as the start-of-step
! rates (star%deuterium_burning_rate_start; the composition step's
! predictor pass).
      logical, intent(in) :: store_start_rates
      double precision, intent(out) :: rate_by_zone(num_rate_rows, json)
      double precision :: deuterium_test
      integer :: zone_idx, clear_idx

      call network_rates_of_model(s, s%logRho, s%logT, s%xa, s%nz, &
           rate_by_zone)
! The deuterium burning rates (MHP 5/02), computed only while the surface
! deuterium, or with inflowing accretion the larger of it and the
! accreted gas's, exceeds 1e-11; zero from the first shell at or below
! the nuclear cutoff out.
      if (s%job%use_mass_accretion .and. s%ctrl%mass_accretion_rate.gt.0.0d0) then
         deuterium_test = max(s%xa(i_h2,s%nz), &
              s%ctrl%accreted_composition(i_h2))
      else
         deuterium_test = s%xa(i_h2,s%nz)
      end if
      do zone_idx = 1, s%nz
         if (s%logT(zone_idx).le.s%ctrl%nuclear_logT_cutoffs(1)) exit
         if (deuterium_test.gt.1.0d-11) then
            s%deuterium_burning_rate(zone_idx) = deuterium_rate( &
                 s%logRho(zone_idx), s%logT(zone_idx), s%xa(i_h1,zone_idx), &
                 zone_idx.ge.s%jcz, s%convective_turnover_timescale)
            if (store_start_rates) s%deuterium_burning_rate_start(zone_idx) = &
                 s%deuterium_burning_rate(zone_idx)
         else
            s%deuterium_burning_rate(zone_idx) = 0.0d0
         end if
      end do
      do clear_idx = zone_idx, s%nz
         s%deuterium_burning_rate(clear_idx) = 0.0d0
      end do
end subroutine compute_burn_rates

!----------------------------------------------------------------------
! burn_network
!----------------------------------------------------------------------
! The implicit network burn over the step (burn_network_region, i.e.
! network's burn_implicitly) of each burn region, which also
! homogenizes the mixed regions.
!   star reads: dm, logT, mixed_zone_bounds
!   star writes: xa
subroutine burn_network(s, radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, rate_by_zone, timestep, ierr)
      use star_info_lib, only: json, max_radiative_zones
      type(star_info), intent(inout) :: s
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
! timestep: s
      double precision, intent(in) :: timestep
      integer, intent(out) :: ierr
      integer, allocatable :: regions(:,:)
      integer :: r

      ierr = 0
      regions = burn_regions(s, s%logT, radiative_zone_bounds, &
           num_radiative_zones, s%mixed_zone_bounds, num_mixed_zones)
      do r = 1, size(regions,2)
         call burn_network_region(s, s%logT, regions(1,r), regions(2,r), &
              rate_by_zone, s%dm, s%xa, timestep, ierr)
         if (ierr /= 0) return
      end do
end subroutine burn_network

!----------------------------------------------------------------------
! burn_deuterium
!----------------------------------------------------------------------
! Deuterium over the step (light_elements's burn_deuterium_region) in
! each burn region.
!   star reads: accreted_mass, deuterium_burning_rate,
!     deuterium_burning_rate_start, dm, logT, mixed_zone_bounds, nz,
!     xa_start
!   star writes: xa
subroutine burn_deuterium(s, radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, timestep)
      use star_info_lib, only: max_radiative_zones
      use net_lib, only: burn_deuterium_region, min_deuterium_step
      type(star_info), intent(inout) :: s
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
! timestep: s
      double precision, intent(in) :: timestep
      integer, allocatable :: regions(:,:)
      integer :: r

      regions = burn_regions(s, s%logT, radiative_zone_bounds, &
           num_radiative_zones, s%mixed_zone_bounds, num_mixed_zones)
      do r = 1, size(regions,2)
         call burn_deuterium_region(regions(1,r), regions(2,r), s%nz, &
              s%dm, s%xa, timestep, s%xa_start, &
              s%deuterium_burning_rate_start, s%deuterium_burning_rate, &
              surface_accretion_now(s, s%accreted_mass), min_deuterium_step)
      end do
end subroutine burn_deuterium

end module burn_step
