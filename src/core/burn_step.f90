!----------------------------------------------------------------------
! burn_step (formerly burn_lib)
!----------------------------------------------------------------------
! The star-layer side of nuclear burning (2026-10 reorganization): the
! code that reads the model, calls the net/ physics, and stores the
! results in star. The physics itself lives in net/: the rates
! (net_lib's compute_reaction_rates), the energy generation
! (nuclear_energy), the network and its implicit and equilibrium
! burns (network), and the light elements (light_elements).
! Here:
!   configure_network     net/'s settings, from the run's controls
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
subroutine configure_network()
      use star_info_lib, only: star
      use net_lib, only: net_settings_t, set_net_settings
      type(net_settings_t) :: settings

      settings%logT_cutoffs = star%ctrl%nuclear_logT_cutoffs
      settings%weak_screening_threshold = star%ctrl%weak_screening_threshold
      settings%cross_section_scale = star%cross_section_scale
      settings%qs0e_scale = star%qs0e_scale
      settings%qqs0ee_scale = star%qqs0ee_scale
      settings%o16_gamma_scale = star%o16_gamma_scale
      settings%c12_alpha_scale = star%c12_alpha_scale
      settings%li6_rate_scale = star%ctrl%li6_rate_scale
      settings%li7_rate_scale = star%ctrl%li7_rate_scale
      settings%be9_pg_rate_scale = star%ctrl%be9_pg_rate_scale
      settings%be9_pd_rate_scale = star%ctrl%be9_pd_rate_scale
      settings%be9_palpha_rate_scale = star%ctrl%be9_palpha_rate_scale
      settings%calc_neutrinos = star%ctrl%calc_neutrinos
      settings%use_itoh_neutrino_loss = star%ctrl%use_itoh_neutrino_loss
      call set_net_settings(settings)
end subroutine configure_network

!----------------------------------------------------------------------
! surface_accretion_now
!----------------------------------------------------------------------
! The gas accreted at the surface over a step in which `mass` grams
! were accreted, as burn_deuterium_region takes it.
!   star reads: ctrl%accreted_composition, ctrl%mass_accretion_rate,
!     job%use_mass_accretion
!   star writes: none
function surface_accretion_now(mass) result(accretion)
      use star_info_lib, only: star, i_h2
      double precision, intent(in) :: mass
      type(surface_accretion) :: accretion

      accretion%on = star%job%use_mass_accretion
      accretion%inflow = star%ctrl%mass_accretion_rate.gt.0.0d0
      accretion%deuterium_fraction = star%ctrl%accreted_composition(i_h2)
      accretion%mass = mass
end function surface_accretion_now

!----------------------------------------------------------------------
! compute_zone_energy
!----------------------------------------------------------------------
! compute_zone_energy (formerly engeb)
! The energy generation of one shell: nuclear_energy, with the two
! inputs that depend on where the shell sits in the model filled in --
! whether it lies in the surface convection zone (shell_index >=
! star%jcz) and the convective turnover time, which together cap
! deuterium burning there. Everything comes back in e; nothing is
! stored. (2026-10: it used to return the neutrino loss, alpha-capture
! and He3 energies, the solar neutrino fluxes and the Be7 abundance
! through star% scalars, and to store per-shell reaction rates and
! branchings that nothing read.)
!   star reads: convective_turnover_timescale, jcz
!   star writes: none
subroutine compute_zone_energy(log_density, log_temperature, xa, &
     shell_index, e)
      use star_info_lib, only: star
      use net_lib, only: compute_nuclear_energy
! xa: the shell's composition, star%xa layout
      double precision, intent(in) :: log_density, log_temperature, xa(:)
      integer, intent(in) :: shell_index
      type(nuclear_energy_result), intent(out) :: e

      call compute_nuclear_energy(log_density, log_temperature, xa, &
           shell_index.ge.star%jcz, star%convective_turnover_timescale, e)
end subroutine compute_zone_energy


!----------------------------------------------------------------------
! average_li_be_over_cz: shared by burn_li_be and the compiled-out
! burn_li_be_midpoint (rotation/burn_li_be_midpoint.f90) (2026
! readability R3). The Li6, Li7, Be9 abundances (x_start, x_end) and the
! log of the burning rates (log_rate_start, log_rate_end), mass-averaged
! over the surface convection zone:
!   start of the step: from cz_base_zone_old up, start-of-step rates;
!   end of the step: from cz_base_zone up, end-of-step rates, with the
!     base zone covering the whole CZ -- its mass extended by
!     cz_base_frac of its width and its rate interpolated halfway as far
!     toward the zone below. cz_mass_end is the CZ's mass then.
! (2026-10: the base-zone adjustment is local; it used to be written into
! shell_mass, restored, and into star%rate_*(cz_base_zone), kept -- the
! rotation sub-steps then carried the adjusted rate into the next
! sub-step's start rates. The start rates also had a branch that reused
! the previous step's CZ rate when star%log_rate_li6_prev > 0, which a
! natural log of a rate per second never is: removed.)
!   star reads: rate_be9, rate_be9_start, rate_li6, rate_li6_start,
!     rate_li7, rate_li7_start
!   star writes: none
subroutine average_li_be_over_cz(composition, mass_coordinate, shell_mass, &
     cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, &
     x_start, log_rate_start, x_end, cz_mass_end, log_rate_end)
      use star_info_lib, only: star, json, i_li6, i_li7, i_be9
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

! START OF THE STEP, AT THE START-OF-STEP DEPTH.
      x_start = 0.0d0
      rate_sum = 0.0d0
      cz_mass_start = 0.0d0
      do zone_idx = cz_base_zone_old,num_zones
         x_start = x_start + composition(species,zone_idx)*shell_mass(zone_idx)
         rate_sum = rate_sum + [star%rate_li6_start(zone_idx), &
              star%rate_li7_start(zone_idx), star%rate_be9_start(zone_idx)]* &
              shell_mass(zone_idx)
         cz_mass_start = cz_mass_start + shell_mass(zone_idx)
      end do
      x_start = x_start/cz_mass_start
      log_rate_start = log(rate_sum/cz_mass_start)
! END OF THE STEP, AT THE END-OF-STEP DEPTH, WITH THE BASE ZONE ADJUSTED
! TO THE TRUE EDGE OF THE CZ (cz_base_frac, from locate_surface_mixing_base).
      x_end = 0.0d0
      rate_sum = 0.0d0
      cz_mass_end = 0.0d0
      do zone_idx = cz_base_zone,num_zones
         mass = shell_mass(zone_idx)
         rate = [star%rate_li6(zone_idx), star%rate_li7(zone_idx), &
              star%rate_be9(zone_idx)]
         if (zone_idx.eq.cz_base_zone .and. cz_base_zone.gt.1 .and. &
              cz_base_zone.lt.num_zones) then
            mass = mass + cz_base_frac* &
                 (mass_coordinate(cz_base_zone)-mass_coordinate(cz_base_zone-1))
            rate = rate + 0.5d0*cz_base_frac*([star%rate_li6(cz_base_zone-1), &
                 star%rate_li7(cz_base_zone-1), star%rate_be9(cz_base_zone-1)] - rate)
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
! decays at its proton-capture rate, interpolated in log between the
! start- and end-of-step rates (store_start_li_be_rates,
! store_end_li_be_rates), with the closed-form survival of
! light_elements (surviving_fraction; 2026-10 -- it replaced
! sub-stepping plus rational extrapolation to zero step, which converged
! on the same integral to 1e-6). In order:
!   1. the base of the surface mixed region (CZ plus envelope overshoot)
!      at the start and end of the step (locate_surface_mixing_base);
!   2. radiative shells below the CZ, each on its own rates (shells
!      hotter than log T = 7 are cleared);
!   3. the CZ as one unit, on its mass-averaged rates
!      (average_li_be_over_cz), starting from the deeper CZ's mixed
!      abundance when the CZ deepened, and taking in accreted gas
!      (accreted_survival) when accretion is on;
!   4. shells the CZ left during the step, on a rate weighted between
!      their radiative rate and the CZ's by the time spent in each.
! Abundances below 1e-24 are set to zero.
! The model arrays are arguments because the two callers burn different
! models: burn_light_elements the converged end-of-step model,
! evolve_angular_momentum each rotation sub-step's midpoint (rot_scr).
!   star reads: accreted_mass, ctrl%accreted_composition,
!     ctrl%mass_accretion_rate, job%use_mass_accretion, rate_be9,
!     rate_be9_start, rate_li6, rate_li6_start, rate_li7, rate_li7_start
!   star writes: none
subroutine burn_li_be(timestep, composition, radius, mass_coordinate, &
     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
      use star_info_lib, only: star, json, i_li6, i_li7, i_be9
      use net_lib, only: surviving_fraction, accreted_survival, &
           linear_in_time, linear_in_log

! timestep: s
      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: radius(json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones

      integer, parameter :: species(3) = [i_li6, i_li7, i_be9]
      integer :: zone_idx, s, cz_base_zone, cz_base_zone_old, min_zone, max_zone
      double precision :: cz_base_frac
      double precision :: rate_start(3), rate_end(3), rate
      double precision :: x_cz_start(3), x_cz_end(3), log_rate_cz_start(3), &
           log_rate_cz_end(3), cz_mass_end
      logical :: accreting
      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
           radiative_frac

! 1. THE BASE OF THE SURFACE MIXED REGION AT THE START AND END OF THE STEP.
! THE DEGREE OF LITHIUM BURNING IN A SURFACE CZ DEPENDS SENSITIVELY ON THE
! TEMPERATURE AT ITS BASE - SO ACCURATELY LOCATING IT IS IMPORTANT.
      call locate_surface_mixing_base(radius, env_cz_zone, env_cz_zone_old, &
           num_zones, cz_base_zone, cz_base_zone_old, cz_base_frac)
      min_zone = min(cz_base_zone,cz_base_zone_old)
      max_zone = max(cz_base_zone,cz_base_zone_old)

! 2. RADIATIVE INTERIOR, UP TO THE SHALLOWER OF THE TWO CZ BASES.
      do zone_idx = 1,min_zone-1
         if(star%rate_be9(zone_idx).le.1.0d-32 .or. star%rate_be9_start(zone_idx).le.1.0d-32)exit
         if(all(composition(species,zone_idx).lt.1.0d-24))cycle
         if(log_temperature(zone_idx).gt.7.0d0)then
            composition(species,zone_idx) = 0.0d0
            cycle
         endif
         rate_start = [star%rate_li6_start(zone_idx), star%rate_li7_start(zone_idx), &
              star%rate_be9_start(zone_idx)]
         rate_end = [star%rate_li6(zone_idx), star%rate_li7(zone_idx), &
              star%rate_be9(zone_idx)]
         do s = 1,3
            composition(species(s),zone_idx) = composition(species(s),zone_idx)* &
                 surviving_fraction(rate_start(s), rate_end(s), timestep, linear_in_log)
            if(composition(species(s),zone_idx).lt.1.0d-24) &
                 composition(species(s),zone_idx) = 0.0d0
         end do
      end do

! 3. THE CONVECTION ZONE (SKIPPED IF IT IS BELOW THE BURNING THRESHOLD).
      if (star%rate_be9_start(cz_base_zone_old).le.1.0d-32.or.star%rate_be9(cz_base_zone).le.1.0d-32) then
         return
      end if
      call average_li_be_over_cz(composition, mass_coordinate, shell_mass, &
           cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, &
           x_cz_start, log_rate_cz_start, x_cz_end, cz_mass_end, log_rate_cz_end)
! FOR THE STARTING ABUNDANCE, USE THE MIXED ABUNDANCE FOR THE DEEPER CZ.
      if(cz_base_zone.lt.cz_base_zone_old) x_cz_start = x_cz_end
! MHP 05/02 MASS ACCRETION: THE ACCRETED GAS BURNS FROM ITS ARRIVAL ON
! (accreted_survival) AND IS MIXED INTO THE CZ BY MASS. THIS ASSUMES THAT
! ALL ACCRETED MATTER LANDS IN A CZ WITH THE END-OF-TIMESTEP DEPTH -- NOT
! STRICTLY TRUE, BUT NOT A BAD APPROXIMATION EITHER.
      accreting = star%job%use_mass_accretion .and. &
           star%ctrl%mass_accretion_rate.gt.0.0d0
      do s = 1,3
         x_cz_end(s) = x_cz_start(s)*surviving_fraction(exp(log_rate_cz_start(s)), &
              exp(log_rate_cz_end(s)), timestep, linear_in_log)
         if (accreting) x_cz_end(s) = (x_cz_end(s)*cz_mass_end + &
              star%ctrl%accreted_composition(species(s))* &
              accreted_survival(exp(log_rate_cz_start(s)), exp(log_rate_cz_end(s)), &
              timestep, linear_in_log)*star%accreted_mass)/ &
              (cz_mass_end + star%accreted_mass)
         if(x_cz_end(s).lt.1.0d-24) x_cz_end(s) = 0.0d0
      end do
      do zone_idx = max_zone,num_zones
         composition(species,zone_idx) = x_cz_end
      end do

! 4. SHELLS THAT BEGAN CONVECTIVE AND ENDED RADIATIVE (A RETREATING CZ).
      if (cz_base_zone.le.cz_base_zone_old) then
         return
      end if
! FIND STARTING AND ENDING LOCATION IN MASS OF THE CZ BASE, AND
! ASSUME THAT THE FRACTION OF TIME SPENT IN THE CZ IS PROPORTIONAL TO
! THE LOCATION IN MASS (I.E. A POINT 1/3 FROM THE OLD TO THE NEW CZ
! BASE IN MASS SPENT 1/3 OF THE TIME IN THE CZ AND 2/3 OUT OF IT):
! radiative_frac is 1 at the new base and 0 at the old one (2026
! bugsweep Batch 2 fixed its sign).
      if(cz_base_zone_old.gt.1)then
         mass_coord_beg = 0.5d0*(mass_coordinate(cz_base_zone_old)+mass_coordinate(cz_base_zone_old-1))
      else
         mass_coord_beg = 0.0d0
      endif
      mass_coord_end = 0.5d0*(mass_coordinate(cz_base_zone)+mass_coordinate(cz_base_zone-1))
      delta_mass = mass_coord_end - mass_coord_beg
      do zone_idx = cz_base_zone_old,cz_base_zone-1
! SKIP IF SHELL TEMPERATURE DROPS BELOW BURNING THRESHOLD.
         if(star%rate_be9(zone_idx).le.1.0d-32)exit
         radiative_frac = (mass_coord_end - mass_coordinate(zone_idx))/delta_mass
         rate_end = [star%rate_li6(zone_idx), star%rate_li7(zone_idx), &
              star%rate_be9(zone_idx)]
         do s = 1,3
! USE FRAD*RADIATIVE RATE AND (1-FRAD)*CONVECTIVE RATE, IN THE LOG.
            rate = exp(radiative_frac*log(rate_end(s)) + &
                 (1.0d0-radiative_frac)*log_rate_cz_start(s))
            composition(species(s),zone_idx) = composition(species(s),zone_idx)* &
                 surviving_fraction(rate, rate, timestep, linear_in_time)
            if(composition(species(s),zone_idx).lt.1.0d-24) &
                 composition(species(s),zone_idx) = 0.0d0
         end do
      end do
end subroutine burn_li_be




!----------------------------------------------------------------------
! store_li_be_rates (formerly lirate88)
!----------------------------------------------------------------------
! The Li6, Li7 and Be9 proton-capture rates per s (net's
! compute_li_be_rates, Caughlan & Fowler 1988, times rho X) of every
! shell, for burn_li_be: with start_of_step, of the saved start-of-step
! model (star%logRho_start, logT_start, xa_start) into star%rate_*_start;
! otherwise of the model in the arguments into star%rate_*. The loop
! stops at the first shell cooler than log T = 6 in both the argument
! model and the start-of-step one; from there out both sets are zeroed.
! (2026-10: the two uses are named, store_start_li_be_rates and
! store_end_li_be_rates; callers passed a 1/2 flag, swapped, so both sets
! came from the start-of-step model.)
!   star reads: logRho_start, logT_start, xa_start
!   star writes: rate_be9, rate_be9_start, rate_li6, rate_li6_start,
!     rate_li7, rate_li7_start
subroutine store_li_be_rates(composition, log_density, log_temperature, &
     num_zones, start_of_step)
      use star_info_lib, only: star, json, i_h1
      use physical_constants, only: ln10
      use net_lib, only: compute_li_be_rates

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
              star%logT_start(zone_idx) < min_log_temperature) exit
         if (start_of_step) then
            rhox = exp(ln10*star%logRho_start(zone_idx))*star%xa_start(i_h1,zone_idx)
            t9 = exp(ln10*(star%logT_start(zone_idx)-9.0d0))
         else
            rhox = exp(ln10*log_density(zone_idx))*composition(i_h1,zone_idx)
            t9 = exp(ln10*(log_temperature(zone_idx)-9.0d0))
         end if
         call compute_li_be_rates(t9, li6_rate, li7_rate, be9_rate)
         if (start_of_step) then
            star%rate_li6_start(zone_idx) = li6_rate*rhox
            star%rate_li7_start(zone_idx) = li7_rate*rhox
            star%rate_be9_start(zone_idx) = be9_rate*rhox
         else
            star%rate_li6(zone_idx) = li6_rate*rhox
            star%rate_li7(zone_idx) = li7_rate*rhox
            star%rate_be9(zone_idx) = be9_rate*rhox
         end if
      end do
      star%rate_li6(zone_idx:num_zones) = 0.0d0
      star%rate_li7(zone_idx:num_zones) = 0.0d0
      star%rate_be9(zone_idx:num_zones) = 0.0d0
      star%rate_li6_start(zone_idx:num_zones) = 0.0d0
      star%rate_li7_start(zone_idx:num_zones) = 0.0d0
      star%rate_be9_start(zone_idx:num_zones) = 0.0d0
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
subroutine store_start_li_be_rates()
      use star_info_lib, only: star
      call store_li_be_rates(star%xa, star%logRho, star%logT, star%nz, .true.)
end subroutine store_start_li_be_rates

!   star reads: none
!   star writes: none
subroutine store_end_li_be_rates(composition, log_density, log_temperature, num_zones)
      use star_info_lib, only: json
      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_temperature(json)
      integer, intent(in) :: num_zones
      call store_li_be_rates(composition, log_density, log_temperature, &
           num_zones, .false.)
end subroutine store_end_li_be_rates


!----------------------------------------------------------------------
! burn_network_region (formerly kemcom; in core/burn_step since 2026-10)
!----------------------------------------------------------------------
! New abundances at the end of the timestep for one radiative shell
! (zone_begin = zone_end) or one convection zone (zone_begin < zone_end,
! homogenized), by an implicit burning step of the network in
! net/network.f90. This routine does only the model bookkeeping
! around it; the reactions, the rate law and the Newton iteration live
! in network (2026-10: they replaced the hand-written 7-species
! equations and Jacobian that stood here).
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
! the two callers differ: burn_network passes the model's, rotmix the
! rotation sub-step's midpoint temperatures (rot_scr). If the Newton
! iteration fails, ierr = 1 and the run stops.
!   star reads: ctrl%absolute_tolerance, ctrl%max_burn_iterations,
!     ctrl%min_abundance, ctrl%nuclear_logT_cutoffs,
!     ctrl%relative_tolerance, xa_start
!   star writes: none
subroutine burn_network_region(log_temperature, zone_begin, zone_end, &
     rate_by_zone, shell_mass, composition, timestep, ierr)

      use star_info_lib, only: star, json, i_h1, i_metals, n_species_basic
      use net_lib, only: burn_implicitly, n_solved, n_network, &
           network_species_slot, s_h1, s_he3, s_he4
      use io_units, only: run_log_unit

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

! start-of-step abundances and rates of the region (mass-weighted
! averages over a convection zone)
      if (zone_begin /= zone_end) then
         region_mass = 0.0d0
         region_start = 0.0d0
         rate_vec = 0.0d0
         do k = zone_begin, zone_end
            region_mass = region_mass + shell_mass(k)
            region_start = region_start + star%xa_start(1:n_species_basic,k)*shell_mass(k)
            rate_vec = rate_vec + shell_mass(k)*rate_by_zone(:,k)
         end do
         region_start = region_start/region_mass
         rate_vec = rate_vec/region_mass
      else
         region_start = star%xa_start(1:n_species_basic,zone_begin)
         rate_vec = rate_by_zone(:,zone_begin)
      end if

! too cool to burn: only homogenize
      if (log_temperature(zone_begin) < star%ctrl%nuclear_logT_cutoffs(1)) then
         do k = zone_begin, zone_end
            composition(1:n_species_basic,k) = region_start
         end do
         return
      end if

! first guess: the shell just solved below (or this one, at the centre)
      guess_zone = max(zone_begin - 1, 1)
      x_start = region_start(network_species_slot)
      x(1:n_solved) = composition(network_species_slot(1:n_solved), guess_zone)
      x(n_solved+1:) = x_start(n_solved+1:)

      call burn_implicitly(x_start, rate_vec, timestep, &
           star%ctrl%min_abundance, star%ctrl%absolute_tolerance, &
           star%ctrl%relative_tolerance, star%ctrl%max_burn_iterations, x, ierr)
      if (ierr /= 0) then
         write (run_log_unit,'(a,i0,a,i0,a)') ' burn_network_region: no '// &
              'convergence for the new abundances of shell ', zone_begin, &
              ' in ', star%ctrl%max_burn_iterations, &
              ' iterations (ctrl%max_burn_iterations); run stopped'
         return
      end if

      if (star%xa_start(i_h1,zone_end) < 5.0d-7) then
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
! The burns of the composition step (core/composition_step.f90).
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
function burn_regions(log_temperature, radiative_zone_bounds, &
     num_radiative_zones, mixed_zone_bounds, num_mixed_zones) result(regions)
      use star_info_lib, only: star
      double precision, intent(in) :: log_temperature(:)
      integer, intent(in) :: radiative_zone_bounds(:,:), mixed_zone_bounds(:,:)
      integer, intent(in) :: num_radiative_zones, num_mixed_zones
      integer, allocatable :: regions(:,:)
      integer :: list(2, size(log_temperature) + num_mixed_zones)
      integer :: n, r, k

      n = 0
      do r = 1, num_radiative_zones
         do k = radiative_zone_bounds(r,1), radiative_zone_bounds(r,2)
            if (log_temperature(k) <= star%ctrl%nuclear_logT_cutoffs(1)) exit
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
subroutine network_rates_of_model(log_density, log_temperature, &
     composition, num_zones, rate_by_zone)
      use star_info_lib, only: star, json
      use net_lib, only: compute_network_rates
      double precision, intent(in) :: log_density(json), &
           log_temperature(json), composition(15,json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: rate_by_zone(num_rate_rows, json)
      integer :: zone_idx, clear_idx

      do zone_idx = 1, num_zones
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
         if (log_temperature(zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
         call compute_network_rates(log_density(zone_idx), log_temperature(zone_idx), &
              composition(:,zone_idx), rate_by_zone(:, zone_idx))
      end do
! ZERO THE RATES BELOW THE CUTOFF (zone_idx = num_zones+1 IF NO SHELL
! WAS COLD).
      do clear_idx = zone_idx, num_zones
         rate_by_zone(:, clear_idx) = 0.0d0
      end do
end subroutine network_rates_of_model

!----------------------------------------------------------------------
! compute_burn_rates
!----------------------------------------------------------------------
! The composition rate vector of every shell of the star above the
! nuclear cutoff (network_rates_of_model) and its deuterium burning rate
! (light_elements's deuterium_rate, stored in
! star%deuterium_burning_rate and, with store_start_rates, _start).
! (2026-10: the rate vectors were also copied to
! rot_scr%reaction_rate_by_zone for rotmix, which now computes its own
! on the sub-step's midpoint model.)
!   star reads: convective_turnover_timescale,
!     ctrl%accreted_composition, ctrl%mass_accretion_rate,
!     ctrl%nuclear_logT_cutoffs, jcz, job%use_mass_accretion, logRho,
!     logT, nz, xa
!   star writes: deuterium_burning_rate, deuterium_burning_rate_start
subroutine compute_burn_rates(store_start_rates, rate_by_zone)
      use star_info_lib, only: star, json, i_h1, i_h2
      use net_lib, only: deuterium_rate
! store_start_rates: also keep the deuterium rates as the start-of-step
! rates (star%deuterium_burning_rate_start; the composition step's
! predictor pass)
      logical, intent(in) :: store_start_rates
      double precision, intent(out) :: rate_by_zone(num_rate_rows, json)
      double precision :: deuterium_test
      integer :: zone_idx, clear_idx

      call network_rates_of_model(star%logRho, star%logT, star%xa, star%nz, &
           rate_by_zone)
! MHP 5/02 THE DEUTERIUM BURNING RATES
      if (star%job%use_mass_accretion .and. star%ctrl%mass_accretion_rate.gt.0.0d0) then
         deuterium_test = max(star%xa(i_h2,star%nz), &
              star%ctrl%accreted_composition(i_h2))
      else
         deuterium_test = star%xa(i_h2,star%nz)
      end if
      do zone_idx = 1, star%nz
         if (star%logT(zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
         if (deuterium_test.gt.1.0d-11) then
            star%deuterium_burning_rate(zone_idx) = deuterium_rate( &
                 star%logRho(zone_idx), star%logT(zone_idx), star%xa(i_h1,zone_idx), &
                 zone_idx.ge.star%jcz, star%convective_turnover_timescale)
            if (store_start_rates) star%deuterium_burning_rate_start(zone_idx) = &
                 star%deuterium_burning_rate(zone_idx)
         else
            star%deuterium_burning_rate(zone_idx) = 0.0d0
         end if
      end do
      do clear_idx = zone_idx, star%nz
         star%deuterium_burning_rate(clear_idx) = 0.0d0
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
subroutine burn_network(radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, rate_by_zone, timestep, ierr)
      use star_info_lib, only: star, json, max_radiative_zones
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
! timestep: s
      double precision, intent(in) :: timestep
      integer, intent(out) :: ierr
      integer, allocatable :: regions(:,:)
      integer :: r

      ierr = 0
      regions = burn_regions(star%logT, radiative_zone_bounds, &
           num_radiative_zones, star%mixed_zone_bounds, num_mixed_zones)
      do r = 1, size(regions,2)
         call burn_network_region(star%logT, regions(1,r), regions(2,r), &
              rate_by_zone, star%dm, star%xa, timestep, ierr)
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
subroutine burn_deuterium(radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, timestep)
      use star_info_lib, only: star, max_radiative_zones
      use net_lib, only: burn_deuterium_region, min_deuterium_step
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
! timestep: s
      double precision, intent(in) :: timestep
      integer, allocatable :: regions(:,:)
      integer :: r

      regions = burn_regions(star%logT, radiative_zone_bounds, &
           num_radiative_zones, star%mixed_zone_bounds, num_mixed_zones)
      do r = 1, size(regions,2)
         call burn_deuterium_region(regions(1,r), regions(2,r), star%nz, &
              star%dm, star%xa, timestep, star%xa_start, &
              star%deuterium_burning_rate_start, star%deuterium_burning_rate, &
              surface_accretion_now(star%accreted_mass), min_deuterium_step)
      end do
end subroutine burn_deuterium

end module burn_step
