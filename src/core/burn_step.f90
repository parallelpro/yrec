!----------------------------------------------------------------------
! burn_step (formerly burn_lib)
!----------------------------------------------------------------------
! The star-layer side of nuclear burning (2026-10 reorganization): the
! code that reads the model, calls the net/ physics, and stores the
! results in star. The physics itself lives in net/: the rates
! (net_lib's nuclear_reaction_rates), the energy generation
! (nuclear_energy_lib), the network and its implicit and equilibrium
! burns (network_lib), and the light elements (light_element_lib).
! Here: compute_zone_energy, which fills in the two energy inputs that
! depend on the shell's place in the model; store_li_be_rates, which
! stores the Li/Be rates per shell; burn_li_be/burn_li_be_midpoint,
! the Li/Be depletion over a step, with the convection-zone base
! location (overshoot included) and the CZ averaging they need.
module burn_step
      use net_lib
! nuclear_energy_result (net/nuclear_energy.f90) is re-exported: it is
! what compute_zone_energy returns.
      use nuclear_energy_lib, only: nuclear_energy, nuclear_energy_result
      use convection_geometry, only: locate_surface_mixing_base
      implicit none
      private :: average_li_be_over_cz, store_li_be_rates

contains












!----------------------------------------------------------------------
! configure_network
!----------------------------------------------------------------------
! Fills net/'s settings (net_lib's net_settings_t) from the run's
! controls and rate factors, and installs them (net_setup). run_yrec
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
      call net_setup(settings)
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
      use light_element_lib, only: surface_accretion
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
subroutine compute_zone_energy(log_density, log_temperature, &
     hydrogen_fraction, helium_fraction, he3_fraction, c12_fraction, &
     c13_fraction, n14_fraction, o16_fraction, o18_fraction, &
     deuterium_fraction, shell_index, e)
      use star_info_lib, only: star
      implicit none
      double precision, intent(in) :: log_density, log_temperature, &
           hydrogen_fraction, helium_fraction, he3_fraction, c12_fraction, &
           c13_fraction, n14_fraction, o16_fraction, o18_fraction, &
           deuterium_fraction
      integer, intent(in) :: shell_index
      type(nuclear_energy_result), intent(out) :: e

      call nuclear_energy(log_density, log_temperature, hydrogen_fraction, &
           helium_fraction, he3_fraction, c12_fraction, c13_fraction, &
           n14_fraction, o16_fraction, o18_fraction, deuterium_fraction, &
           shell_index.ge.star%jcz, star%convective_turnover_timescale, e)
end subroutine compute_zone_energy


!----------------------------------------------------------------------
! average_li_be_over_cz: shared by burn_li_be and burn_li_be_midpoint (2026
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
      use math_lib
      implicit none
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
! light_element_lib (surviving_fraction; 2026-10 -- it replaced
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
      use light_element_lib, only: surviving_fraction, accreted_survival, &
           linear_in_time, linear_in_log
      use math_lib
      implicit none

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


! ---------------------------------------------------------------
! Compiled out (2026-10): the Bulirsch-Stoer mixing-plus-burning path
! (burn_settle_mix -> burn_mix_extrapolated, burn_li_be_midpoint) is switched off
! by a hard-coded burs_extrapolation_active = .false. in
! evolve_angular_momentum, with no control to enable it, so it never
! ran. Kept, not deleted: build with -DYREC_BULIRSCH_STOER_BURN to
! compile it again (and set burs_extrapolation_active).
! ---------------------------------------------------------------
#ifdef YREC_BULIRSCH_STOER_BURN
!----------------------------------------------------------------------
! burn_li_be_midpoint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original liburn2.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Determine lithium-6, lithium-7, and beryllium-9 burning. The burning
! rates depend on the local T and rho, the abundance, and (in
! convection zones) the mixing.
!
! In radiative regions burning is done implicitly, using a single
! trapezoidal (midpoint-in-the-log) estimate of the burning rate over
! the timestep -- unlike burn_li_be, this variant does not iterate
! sub-step refinement with a rational-function extrapolation.
!
! In the convective region the program likewise uses a single
! time-averaged (midpoint-in-the-log) rate over the timestep, again
! without the iterative sub-step refinement used by burn_li_be.
!
! The degree of lithium burning in a surface CZ depends sensitively
! on the temperature at its base -- so accurately locating it is
! important. Determine the true location (cz_base_frac) of the base
! of the CZ at the end of the timestep, and the location of the edge
! of overshoot regions if applicable.
!
! 11/91 HR added to call.
!
! burn_li_be_midpoint (formerly liburn2)
!   star reads: rate_be9, rate_be9_start, rate_li6, rate_li6_start,
!     rate_li7, rate_li7_start
!   star writes: none
subroutine burn_li_be_midpoint(timestep, composition, radius, mass_coordinate, &
     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
      use rotation_scratch_lib
      use star_info_lib, only: star, json, i_li6, i_li7, i_be9
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: radius(json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones

      double precision :: li6_substep_depletion, li7_substep_depletion, &
           be9_substep_depletion
      integer :: zone_idx, min_zone, max_zone
      double precision :: cz_base_frac
      integer :: cz_base_zone, cz_base_zone_old
      double precision :: log_rate_li6, log_rate_li7, log_rate_be9
      double precision :: li6_cz_start, li7_cz_start, be9_cz_start
      double precision :: log_rate_li6_cz_start, log_rate_li7_cz_start, &
           log_rate_be9_cz_start
      double precision :: li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end
      double precision :: x_cz_start(3), x_cz_end(3), log_rate_cz_start(3), &
           log_rate_cz_end(3)
      double precision :: log_rate_li6_cz_end, log_rate_li7_cz_end, &
           log_rate_be9_cz_end
      double precision :: cz_log_rate_li6, cz_log_rate_li7, cz_log_rate_be9
      double precision :: li6_depletion, li7_depletion, be9_depletion
      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
           radiative_frac

! THE DEGREE OF LITHIUM BURNING IN A SURFACE CZ DEPENDS SENSITIVELY
! ON THE TEMPERATURE AT ITS BASE - SO ACCURATELY LOCATING IS IMPORTANT.
! DETERMINE THE TRUE LOCATION (FX) OF THE BASE OF THE CZ AT THE END
! OF THE TIMESTEP, AND THE LOCATION OF THE EDGE OF OVERSHOOT REGIONS
! IF APPLICABLE.
      call locate_surface_mixing_base(radius, env_cz_zone, env_cz_zone_old, num_zones, &
           cz_base_zone, cz_base_zone_old, cz_base_frac)
! RADIATIVE INTERIOR.
      min_zone = min(cz_base_zone,cz_base_zone_old)
      max_zone = max(cz_base_zone,cz_base_zone_old)
      do zone_idx = 1,min_zone-1
         if(star%rate_be9(zone_idx).le.1.0d-32 .or. star%rate_be9_start(zone_idx).le.1.0d-32)cycle
         if(composition(i_li6,zone_idx).lt.1.0d-24.and.composition(i_li7,zone_idx).lt.1.0d-24 &
         .and.composition(i_be9,zone_idx).lt.1.0d-24)cycle
         if(log_temperature(zone_idx).gt.7.0d0)then
            composition(i_li6,zone_idx) = 0.0d0
            composition(i_li7,zone_idx) = 0.0d0
            composition(i_be9,zone_idx) = 0.0d0
            cycle
         endif
         log_rate_li6 = 0.5d0*(log(star%rate_li6(zone_idx)) + log(star%rate_li6_start(zone_idx)))
         log_rate_li7 = 0.5d0*(log(star%rate_li7(zone_idx)) + log(star%rate_li7_start(zone_idx)))
         log_rate_be9 = 0.5d0*(log(star%rate_be9(zone_idx)) + log(star%rate_be9_start(zone_idx)))
         li6_substep_depletion = timestep*exp(log_rate_li6)
         li7_substep_depletion = timestep*exp(log_rate_li7)
         be9_substep_depletion = timestep*exp(log_rate_be9)
! THE REACTION RATES ARE OF THE FORM
!    DLI/DT = F(RHO,T,X)*LI =>DLN LI = DT*F(RHO,T,X)
! SOLVE FOR D LN (SPECIES)/DT AND ZERO OUT IF THE DEPLETION IS TOO HIGH.
         if(li6_substep_depletion.lt.3.0d1)then
            composition(i_li6,zone_idx) = composition(i_li6,zone_idx)/ &
                 exp(li6_substep_depletion)
         else
            composition(i_li6,zone_idx) = 0.0d0
         endif
         if(li7_substep_depletion.lt.3.0d1)then
            composition(i_li7,zone_idx) = composition(i_li7,zone_idx)/ &
                 exp(li7_substep_depletion)
         else
            composition(i_li7,zone_idx) = 0.0d0
         endif
         if(be9_substep_depletion.lt.3.0d1)then
            composition(i_be9,zone_idx) = composition(i_be9,zone_idx)/ &
                 exp(be9_substep_depletion)
         else
            composition(i_be9,zone_idx) = 0.0d0
         endif
! WRITE NEW ABUNDANCES AND EXIT.
         if(composition(i_li6,zone_idx).lt.1.0d-24)composition(i_li6,zone_idx)=0.0d0
         if(composition(i_li7,zone_idx).lt.1.0d-24)composition(i_li7,zone_idx)=0.0d0
         if(composition(i_be9,zone_idx).lt.1.0d-24)composition(i_be9,zone_idx)=0.0d0
      end do
! CONVECTION ZONE.
!
! SKIP IF WHOLE CZ IS BELOW THE BURNING THRESHOLD.
      if (star%rate_be9_start(cz_base_zone_old).le.1.0d-32.or.star%rate_be9(cz_base_zone).le.1.0d-32) then
         return
      end if
      call average_li_be_over_cz(composition, mass_coordinate, shell_mass, &
           cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, &
           x_cz_start, log_rate_cz_start, x_cz_end, cz_mass_end, log_rate_cz_end)
      li6_cz_start = x_cz_start(1)
      li7_cz_start = x_cz_start(2)
      be9_cz_start = x_cz_start(3)
      log_rate_li6_cz_start = log_rate_cz_start(1)
      log_rate_li7_cz_start = log_rate_cz_start(2)
      log_rate_be9_cz_start = log_rate_cz_start(3)
      li6_cz_end = x_cz_end(1)
      li7_cz_end = x_cz_end(2)
      be9_cz_end = x_cz_end(3)
      log_rate_li6_cz_end = log_rate_cz_end(1)
      log_rate_li7_cz_end = log_rate_cz_end(2)
      log_rate_be9_cz_end = log_rate_cz_end(3)
! FOR THE STARTING ABUNDANCE, USE THE MIXED ABUNDANCE FOR THE
! DEEPER CZ.
      if(cz_base_zone.lt.cz_base_zone_old)then
         li6_cz_start = li6_cz_end
         li7_cz_start = li7_cz_end
         be9_cz_start = be9_cz_end
      endif
! INITIALIZE ABUNDANCES.
      li6_cz_end = li6_cz_start
      li7_cz_end = li7_cz_start
      be9_cz_end = be9_cz_start
! TIME-AVERAGED RATES USING LINEAR INTERPOLATION IN THE LOG.
      cz_log_rate_li6 = 0.5d0*(log_rate_li6_cz_end+log_rate_li6_cz_start)
      cz_log_rate_li7 = 0.5d0*(log_rate_li7_cz_end+log_rate_li7_cz_start)
      cz_log_rate_be9 = 0.5d0*(log_rate_be9_cz_end+log_rate_be9_cz_start)
      li6_depletion = timestep*exp(cz_log_rate_li6)
      li7_depletion = timestep*exp(cz_log_rate_li7)
      be9_depletion = timestep*exp(cz_log_rate_be9)
      if(li6_depletion.lt.3.0d1 .and. li6_cz_end.gt.1.0d-24)then
         li6_cz_end = li6_cz_end/exp(li6_depletion)
      else
         li6_cz_end = 0.0d0
      endif
      if(li6_cz_end.lt.1.0d-24)li6_cz_end = 0.0d0
      if(li7_depletion.lt.3.0d1 .and. li7_cz_end.gt.1.0d-24)then
         li7_cz_end = li7_cz_end/exp(li7_depletion)
      else
         li7_cz_end = 0.0d0
      endif
      if(li7_cz_end.lt.1.0d-24)li7_cz_end = 0.0d0
      if(be9_depletion.lt.3.0d1 .and. be9_cz_end.gt.1.0d-24)then
         be9_cz_end = be9_cz_end/exp(be9_depletion)
      else
         be9_cz_end = 0.0d0
      endif
      if(be9_cz_end.lt.1.0d-24)be9_cz_end = 0.0d0
      do zone_idx = max_zone,num_zones
         composition(i_li6,zone_idx) = li6_cz_end
         composition(i_li7,zone_idx) = li7_cz_end
         composition(i_be9,zone_idx) = be9_cz_end
      end do
! NOW SOLVE FOR ABUNDANCES IN THE REGION WHICH BEGAN CONVECTIVE AND
! ENDED RADIATIVE.
      if (cz_base_zone.le.cz_base_zone_old) then
         return
      end if
! FIND STARTING AND ENDING LOCATION IN MASS OF THE CZ BASE, AND
! ASSUME THAT THE FRACTION OF TIME SPENT IN THE CZ IS PROPORTIONAL TO
! THE LOCATION IN MASS (I.E. A POINT 1/3 FROM THE OLD TO THE NEW CZ
! BASE IN MASS SPENT 1/3 OF THE TIME IN THE CZ AND 2/3 OUT OF IT).
      if(cz_base_zone_old.gt.1)then
         mass_coord_beg = 0.5d0*(mass_coordinate(cz_base_zone_old)+mass_coordinate(cz_base_zone_old-1))
      else
         mass_coord_beg = 0.0d0
      endif
      mass_coord_end = 0.5d0*(mass_coordinate(cz_base_zone)+mass_coordinate(cz_base_zone-1))
! 2026 (bugsweep Batch 2): this block runs only for a retreating CZ
! (cz_base_zone > cz_base_zone_old), so the original
! delta_mass = beg - end was negative and FRAD = (m - beg)/delta_mass
! ran from 0 at the old base to -1 at the new one. Per the comment
! above FRAD must be 1 at the old base and 0 at the new base:
! FRAD = (end - m)/(end - beg).
      delta_mass = mass_coord_end - mass_coord_beg
      do zone_idx = cz_base_zone_old,cz_base_zone-1
! MHP 9/91 CHANGE TO AVOID DIVISION BY ZERO.
! SKIP IF SHELL TEMPERATURE DROPS BELOW BURNING THRESHOLD.
         if(star%rate_be9(zone_idx).le.1.0d-32)exit
         radiative_frac = (mass_coord_end - mass_coordinate(zone_idx))/delta_mass
! USE FRAD*RADIATIVE RATE AND (1-FRAD)*CONVECTIVE RATE.
         li6_depletion = timestep*exp(radiative_frac*log(star%rate_li6(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li6_cz_start)
         li7_depletion = timestep*exp(radiative_frac*log(star%rate_li7(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li7_cz_start)
         be9_depletion = timestep*exp(radiative_frac*log(star%rate_be9(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_be9_cz_start)
!***REMEMBER TO ADD FAILSAFES FOR LARGE DEPLETION***
         composition(i_li6,zone_idx) = composition(i_li6,zone_idx)/exp(li6_depletion)
         if(composition(i_li6,zone_idx).lt.1.0d-24)composition(i_li6,zone_idx)=0.0d0
         composition(i_li7,zone_idx) = composition(i_li7,zone_idx)/exp(li7_depletion)
         if(composition(i_li7,zone_idx).lt.1.0d-24)composition(i_li7,zone_idx)=0.0d0
         composition(i_be9,zone_idx) = composition(i_be9,zone_idx)/exp(be9_depletion)
         if(composition(i_be9,zone_idx).lt.1.0d-24)composition(i_be9,zone_idx)=0.0d0
      end do
      return
end subroutine burn_li_be_midpoint
#endif


!----------------------------------------------------------------------
! store_li_be_rates
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lirate88.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine computes the burning rates for Li6, Li7, Be9 at the
! beginning (use_current_model = 2: star%rate_*_start, from the
! start-of-step model) or end (use_current_model = 1: star%rate_*,
! from the arguments) of a timestep, for burn_li_be/burn_li_be_midpoint.
!
! Burning rates from Caughlin and Fowler (1988): net's li_be_rates
! (2026-10: the formulas moved there; this loop stores them).
! Callers: flag 2 at the start of the step (start-of-step model ->
! star%rate_*_start), flag 1 at the end (the model passed in ->
! star%rate_*). 2026-10: the calls had the flags the other way round
! (inherited from F77), so both sets came from the start-of-step model.
! store_li_be_rates (formerly lirate88)
!   star reads: logRho_start, logT_start, xa_start
!   star writes: rate_be9, rate_be9_start, rate_li6, rate_li6_start,
!     rate_li7, rate_li7_start
subroutine store_li_be_rates(composition, log_density, log_temperature, num_zones, &
     use_current_model)
      use star_info_lib, only: star, json, i_be9
      use phys_const_lib
      use math_lib
      use light_element_lib, only: li_be_rates
      implicit none

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_density(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: num_zones, use_current_model

      double precision :: tlim
      data tlim/6.0d0/
      integer :: zone_idx, tail_idx
      double precision :: rhox, t9, fli6, fli7, fsbe9

      do zone_idx = 1,num_zones
         if(log_temperature(zone_idx).lt.tlim.and.star%logT_start(zone_idx).lt.tlim)exit
         if(use_current_model.eq.1)then
            rhox = exp(ln10*log_density(zone_idx))*composition(1,zone_idx)
            t9=exp(ln10*(log_temperature(zone_idx)-9.0d0))
         else
            rhox = exp(ln10*star%logRho_start(zone_idx))*star%xa_start(1,zone_idx)
            t9=exp(ln10*(star%logT_start(zone_idx)-9.0d0))
         endif
         call li_be_rates(t9, fli6, fli7, fsbe9)
         if(use_current_model.eq.1)then
            star%rate_li6(zone_idx) = fli6*rhox
            star%rate_li7(zone_idx) = fli7*rhox
            star%rate_be9(zone_idx) = fsbe9*rhox
         else
            star%rate_li6_start(zone_idx) = fli6*rhox
            star%rate_li7_start(zone_idx) = fli7*rhox
            star%rate_be9_start(zone_idx) = fsbe9*rhox
         endif
      end do
      do tail_idx = zone_idx,num_zones
         star%rate_li6(tail_idx)=0.0d0
         star%rate_li7(tail_idx)=0.0d0
         star%rate_be9(tail_idx)=0.0d0
         star%rate_li6_start(tail_idx)=0.0d0
         star%rate_li7_start(tail_idx)=0.0d0
         star%rate_be9_start(tail_idx)=0.0d0
      end do
      return
end subroutine store_li_be_rates

!----------------------------------------------------------------------
! store_start_li_be_rates / store_end_li_be_rates
!----------------------------------------------------------------------
! The two uses of store_li_be_rates, named (2026-10: callers passed its 1/2 flag
! directly, and had them swapped). Start of step: the rates of the
! saved start-of-step model (star%logRho_start, logT_start, xa_start)
! into star%rate_*_start. End of step: the rates of the model passed in
! (the converged model, or a rotation sub-step's midpoint model) into
! star%rate_*.
!   star reads: logRho, logT, nz, xa
!   star writes: none
subroutine store_start_li_be_rates()
      use star_info_lib, only: star
      call store_li_be_rates(star%xa, star%logRho, star%logT, star%nz, 2)
end subroutine store_start_li_be_rates

!   star reads: none
!   star writes: none
subroutine store_end_li_be_rates(composition, log_density, log_temperature, num_zones)
      use star_info_lib, only: json
      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_temperature(json)
      integer, intent(in) :: num_zones
      call store_li_be_rates(composition, log_density, log_temperature, num_zones, 1)
end subroutine store_end_li_be_rates


!----------------------------------------------------------------------
! burn_network_region (formerly kemcom; in core/burn_step since 2026-10)
!----------------------------------------------------------------------
! New abundances at the end of the timestep for one radiative shell
! (zone_begin = zone_end) or one convection zone (zone_begin < zone_end,
! homogenized), by an implicit burning step of the network in
! net/network_lib.f90. This routine does only the model bookkeeping
! around it; the reactions, the rate law and the Newton iteration live
! in network_lib (2026-10: they replaced the hand-written 7-species
! equations and Jacobian that stood here).
!
!  log_temperature     run of log10 T; shells below
!                      ctrl%nuclear_logT_cutoffs(1) are only mixed
!  rate_by_zone(rr_*,k) net_lib's rate vectors, per s per amu
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
! rotation sub-step's midpoint temperatures (rot_scr).
!   star reads: ctrl%absolute_tolerance, ctrl%max_burn_iterations,
!     ctrl%min_abundance, ctrl%nuclear_logT_cutoffs,
!     ctrl%relative_tolerance, xa_start
!   star writes: none
subroutine burn_network_region(log_temperature, zone_begin, zone_end, &
     rate_by_zone, shell_mass, composition, timestep, ierr)

      use star_info_lib, only: star, json, i_h1, i_he3, i_he4, i_metals, &
           i_c12, i_c13, i_n14, i_o16, i_o18, n_species_basic
      use network_lib, only: implicit_burn_step, n_solved, n_network, &
           s_h1, s_he3, s_he4
      use luout_lib
      implicit none

      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: zone_begin, zone_end
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
! timestep: s
      double precision, intent(in) :: timestep
      integer, intent(out) :: ierr

! star%xa slot of each network species (network_lib's s_* order)
      integer, parameter :: xa_slot(n_network) = &
           [i_h1, i_he3, i_he4, i_c12, i_c13, i_n14, i_o16, i_o18]
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
      x_start = region_start(xa_slot)
      x(1:n_solved) = composition(xa_slot(1:n_solved), guess_zone)
      x(n_solved+1:) = x_start(n_solved+1:)

      call implicit_burn_step(x_start, rate_vec, timestep, &
           star%ctrl%min_abundance, star%ctrl%absolute_tolerance, &
           star%ctrl%relative_tolerance, star%ctrl%max_burn_iterations, x, ierr)
      if (ierr /= 0) then
         write (run_log_unit,1000) zone_begin
 1000    format(1X,39('>'),40('<')/1X,'ERROR IN SUBROUTINE KEMCOM'/ &
              1X,'UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL',I4/1X, &
              'RUN STOPPED AFTER 50 ATTEMPTS')
         return
      end if

      if (star%xa_start(i_h1,zone_end) < 5.0d-7) then
         metals = 1.0d0 - x(s_h1) - x(s_he3) - x(s_he4)
      else
         metals = region_start(i_metals)
      end if
      do k = zone_begin, zone_end
         composition(xa_slot,k) = x
         composition(i_metals,k) = metals
      end do
end subroutine burn_network_region

! =================================================================
! The burning steps of the composition step (core/composition_step.f90,
! formerly mix; moved here 2026-10 with bodies unchanged).
! =================================================================

!----------------------------------------------------------------------
! network_rates_of_model
!----------------------------------------------------------------------
! The composition rate vector (net_lib's network_rate_vector) of every
! shell of a model hotter than the nuclear cutoff; zero in the cooler
! shells. The model is in the arguments: compute_burn_rates passes the
! star's, mid_timestep_model a rotation sub-step's midpoint (rot_scr).
!   star reads: ctrl%nuclear_logT_cutoffs
!   star writes: none
subroutine network_rates_of_model(log_density, log_temperature, &
     composition, num_zones, rate_by_zone)
      use star_info_lib, only: star, json
      double precision, intent(in) :: log_density(json), &
           log_temperature(json), composition(15,json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: rate_by_zone(num_rate_rows, json)
      integer :: zone_idx, clear_idx

      do zone_idx = 1, num_zones
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
         if (log_temperature(zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
         call network_rate_vector(log_density(zone_idx), log_temperature(zone_idx), &
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
! (light_element_lib's deuterium_rate, stored in
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
      use light_element_lib, only: deuterium_rate
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
              star%ctrl%accreted_composition(12))
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
! network_lib's implicit_burn_step): shell by shell in the radiative
! regions above the nuclear cutoff, and once per mixed region, which it
! also homogenizes.
!   star reads: ctrl%nuclear_logT_cutoffs, dm, logT, mixed_zone_bounds
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
      integer :: radiative_region_idx, inner_zone_idx, mixed_zone_idx, &
           zone_begin, zone_end

      ierr = 0
!
!  NOW IMPLICITLY SOLVE FOR THE NEW ABUNDANCES AT THE END OF THE
!  TIMESTEP.  THIS IS DONE SHELL BY SHELL FOR RADIATIVE REGIONS,
!  AND FOR EACH CONVECTION ZONE AS A UNIT.
!
! RADIATIVE ZONES.
!
      do radiative_region_idx = 1, num_radiative_zones
         do inner_zone_idx = radiative_zone_bounds(radiative_region_idx,1), &
              radiative_zone_bounds(radiative_region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
            if (star%logT(inner_zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
            zone_begin = inner_zone_idx
            zone_end = inner_zone_idx
            call burn_network_region(star%logT, zone_begin, zone_end, &
                 rate_by_zone, star%dm, star%xa, timestep, ierr)
            if (ierr /= 0) return
         end do
      end do
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
      do mixed_zone_idx = 1, num_mixed_zones
         zone_begin = star%mixed_zone_bounds(mixed_zone_idx,1)
         zone_end = star%mixed_zone_bounds(mixed_zone_idx,2)
         call burn_network_region(star%logT, zone_begin, zone_end, &
              rate_by_zone, star%dm, star%xa, timestep, ierr)
         if (ierr /= 0) return
      end do
end subroutine burn_network

!----------------------------------------------------------------------
! burn_equilibrium
!----------------------------------------------------------------------
! The explicit burn of iteration level 1 (network_lib's equilibrium_burn_rates, He3 and
! the CN cycle in equilibrium): overwrites X, and under helium burning
! Z, C12 and O16, shell by shell in radiative regions and per mixed
! region.
!   star reads: ctrl%nuclear_logT_cutoffs, dm, logT, mixed_zone_bounds,
!     xa_start
!   star writes: xa
subroutine burn_equilibrium(radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, rate_by_zone, timestep)
      use star_info_lib, only: star, json, max_radiative_zones, i_h1, &
           i_metals, i_c12, i_o16
      use network_lib, only: equilibrium_burn_rates
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
! timestep: s
      double precision, intent(in) :: timestep
      integer :: radiative_region_idx, inner_zone_idx, mixed_zone_idx, &
           zone_begin, zone_end, zone_idx
      double precision :: dc_dt, do_dt, dx_dt, dy_dt
! region abundances returned by equilibrium_burn_rates
      double precision :: hydrogen_fraction, metal_fraction, c12_fraction, &
           o16_fraction

!
! IF ITLVL=1 THEN THE RATES OF HYDROGEN AND HELIUM BURNING ARE
!   COMPUTED EXPLICITLY ASSUMING EQUILIBRIUM HE3 AND CN.
!   IF THIS IS THE CASE, CALL OLD ENERGY GENERATION ROUTINE AND
!   OVERWRITE THE HYDROGEN ABUNDANCES (H BURNING) AND Z,C12,O16
!   ABUNDANCES (HE BURNING).
!
!  EXPLICITLY SOLVE FOR THE NEW HYDROGEN AND HELIUM BURNING RATES
!  AT THE BEGINNING OF THE TIMESTEP, ASSUMING EQUILIBRIUM HE3 AND CN
!  CYCLE ABUNDANCES.  THIS IS DONE SHELL BY SHELL FOR RADIATIVE REGIONS,
!  AND FOR EACH CONVECTION ZONE AS A UNIT.
!
! RADIATIVE ZONES.
!
         do radiative_region_idx = 1, num_radiative_zones
            do inner_zone_idx = &
                 radiative_zone_bounds(radiative_region_idx,1), &
                 radiative_zone_bounds(radiative_region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (star%logT(inner_zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
               zone_begin = inner_zone_idx
               zone_end = inner_zone_idx
               call equilibrium_burn_rates(rate_by_zone, star%xa_start, star%dm, &
                    star%logT, zone_begin, zone_end, dc_dt, do_dt, &
                    dx_dt, dy_dt, c12_fraction, o16_fraction, &
                    hydrogen_fraction, metal_fraction)
!  USE THE EXPLICIT HYDROGEN BURNING RATE.
               if (dx_dt.ne.0.0d0) then
                  star%xa(i_h1,inner_zone_idx) = hydrogen_fraction + &
                       dx_dt*timestep
               end if
!  USE THE HELIUM BURNING RATE FROM EQBURN AND THE CARBON,ALPHA
!  BURNING RATE
               if (dy_dt.ne.0.0d0) then
                  star%xa(i_metals,inner_zone_idx) = metal_fraction - &
                       dy_dt*timestep
                  star%xa(i_c12,inner_zone_idx) = c12_fraction + &
                       dc_dt*timestep
                  star%xa(i_o16,inner_zone_idx) = o16_fraction + &
                       do_dt*timestep
               end if
            end do
         end do
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
         do mixed_zone_idx = 1, num_mixed_zones
            zone_begin = star%mixed_zone_bounds(mixed_zone_idx,1)
            zone_end = star%mixed_zone_bounds(mixed_zone_idx,2)
            call equilibrium_burn_rates(rate_by_zone, star%xa_start, star%dm, star%logT, &
                 zone_begin, zone_end, dc_dt, do_dt, dx_dt, dy_dt, &
                 c12_fraction, o16_fraction, hydrogen_fraction, &
                 metal_fraction)
!  USE THE EXPLICIT HYDROGEN BURNING RATE.
            if (dx_dt.ne.0.0d0) then
               hydrogen_fraction = hydrogen_fraction + dx_dt*timestep
               do zone_idx = zone_begin, zone_end
                  star%xa(i_h1,zone_idx) = hydrogen_fraction
               end do
            end if
!  USE THE HELIUM BURNING RATE FROM EQBURN AND THE CARBON,ALPHA
!  BURNING RATE
            if (dy_dt.ne.0.0d0) then
               metal_fraction = metal_fraction - dy_dt*timestep
               c12_fraction = c12_fraction + dc_dt*timestep
               o16_fraction = o16_fraction + do_dt*timestep
               do zone_idx = zone_begin, zone_end
                  star%xa(i_metals,zone_idx) = metal_fraction
                  star%xa(i_c12,zone_idx) = c12_fraction
                  star%xa(i_o16,zone_idx) = o16_fraction
               end do
            end if
         end do
end subroutine burn_equilibrium

!----------------------------------------------------------------------
! burn_deuterium
!----------------------------------------------------------------------
! Deuterium over the step (light_element_lib's burn_deuterium_region), shell by shell in
! radiative regions above the nuclear cutoff and per mixed region.
!   star reads: accreted_mass, ctrl%nuclear_logT_cutoffs,
!     deuterium_burning_rate, deuterium_burning_rate_start, dm, logT,
!     mixed_zone_bounds, nz, xa_start
!   star writes: xa
subroutine burn_deuterium(radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, timestep)
      use star_info_lib, only: star, max_radiative_zones
      use light_element_lib, only: burn_deuterium_region, min_deuterium_step
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
! timestep: s
      double precision, intent(in) :: timestep
      integer :: radiative_region_idx, inner_zone_idx, mixed_zone_idx, &
           zone_begin, zone_end

! MHP 5/02 ADDED DEUTERIUM BURNING
! RADIATIVE ZONES.
!
         do radiative_region_idx = 1, num_radiative_zones
            do inner_zone_idx = &
                 radiative_zone_bounds(radiative_region_idx,1), &
                 radiative_zone_bounds(radiative_region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (star%logT(inner_zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
               zone_begin = inner_zone_idx
               zone_end = inner_zone_idx
               call burn_deuterium_region(zone_begin, zone_end, star%nz, star%dm, &
                    star%xa, timestep, star%xa_start, star%deuterium_burning_rate_start, &
                    star%deuterium_burning_rate, surface_accretion_now(star%accreted_mass), &
                    min_deuterium_step)
            end do
         end do
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
         do mixed_zone_idx = 1, num_mixed_zones
            zone_begin = star%mixed_zone_bounds(mixed_zone_idx,1)
            zone_end = star%mixed_zone_bounds(mixed_zone_idx,2)
            call burn_deuterium_region(zone_begin, zone_end, star%nz, star%dm, &
                 star%xa, timestep, star%xa_start, star%deuterium_burning_rate_start, &
                 star%deuterium_burning_rate, surface_accretion_now(star%accreted_mass), &
                    min_deuterium_step)
         end do
end subroutine burn_deuterium

end module burn_step
