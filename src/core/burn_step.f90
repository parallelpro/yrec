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
      implicit none
      private :: locate_cz_base, average_li_be_over_cz, store_li_be_rates

contains












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
!   star reads: jcz, convective_turnover_timescale
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
! locate_cz_base: shared by burn_li_be and burn_li_be_midpoint (2026 readability R3;
! the block was byte-identical in both). Locates the true base of the
! surface convection zone at the end of the step (cz_base_zone, with
! cz_base_frac the offset of the true base from the zone midpoint in
! units of the zone width: negative = higher, positive = lower) and
! the base at the start of the step (cz_base_zone_old), including
! envelope overshoot when active, and records cz_base_zone in
! star%envelope_cz_base_zone_prev.
! cz_base_frac is intent(inout): it is assigned only when the CZ base
! is an interior zone (1 < env_cz_zone < num_zones), exactly as in the
! original inline block, and is never read by the callers otherwise.
subroutine locate_cz_base(radius, env_cz_zone, env_cz_zone_old, num_zones, &
     cz_base_zone, cz_base_zone_old, cz_base_frac)
      use rotation_scratch_lib
      use star_info_lib, only: star, json
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: radius(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones
      integer, intent(out) :: cz_base_zone, cz_base_zone_old
      double precision, intent(inout) :: cz_base_frac

      integer :: zone_idx
      double precision :: del_diff, del_diff_below
      double precision :: search_radius, shell_radius, delta_radius, &
           cz_base_radius

      if(env_cz_zone.gt.1.and.env_cz_zone.lt.num_zones)then
         if(star%job%rotation_active.and.star%job%instability_transport_active)then
            del_diff = mix_scr%del_adiabatic_mix(env_cz_zone) - &
                 mix_scr%del_radiative_mix(env_cz_zone)
            del_diff_below = mix_scr%del_adiabatic_mix(env_cz_zone-1) - &
                 mix_scr%del_radiative_mix(env_cz_zone-1)
         else
! EVALUATE DEL(AD) - DEL(RAD) AT THE LAST CONVECTIVE POINT AND THE ONE
! BELOW IT.
            del_diff = star%grada(env_cz_zone)-star%gradr(env_cz_zone)
            del_diff_below = star%grada(env_cz_zone-1)-star%gradr(env_cz_zone-1)
         endif
! USE LINEAR INTERPOLATION TO FIND THE DISTANCE OF THE TRUE LOCATION
! OF THE BASE FROM THE ZONE MIDPOINT. IF FX IS NEGATIVE,THEN THE TRUE
! BASE IS HIGHER; IF IT IS POSITIVE, THE TRUE BASE IS LOWER.
         cz_base_frac = max(-0.5d0,0.5d0-del_diff_below/ &
              (del_diff_below-del_diff))
         cz_base_frac = min(0.5d0,cz_base_frac)
         if(.not.star%job%envelope_overshoot_active)then
            cz_base_zone = env_cz_zone
            cz_base_zone_old = env_cz_zone_old
         else
! STARTING CZ DEPTH (cz_base_radius_prev = 0 means no previous step
! has stored one: locate it from the start-of-step model).
            if(star%cz_base_radius_prev.eq.0.0d0)then
               star%cz_base_radius_prev = 0.5d0*(exp(ln10*star%logR_start(env_cz_zone_old)) &
                        +exp(ln10*star%logR_start(env_cz_zone_old-1)))
               search_radius = star%cz_base_radius_prev - star%pressure_scale_height_start
               do zone_idx = env_cz_zone_old-1,1,-1
                  shell_radius = exp(ln10*star%logR_start(zone_idx))
                  if(shell_radius.lt.search_radius)then
                     cz_base_zone_old = zone_idx + 1
                     exit
                  endif
               end do
               if (zone_idx .lt. 1) cz_base_zone_old = 1
            else
               cz_base_zone_old = star%envelope_cz_base_zone_prev
            endif
! ENDING CZ DEPTH : DETERMINE OVERSHOOT FROM TRUE CZ BASE.
            delta_radius = exp(ln10*radius(env_cz_zone))-exp(ln10*radius(env_cz_zone-1))
            cz_base_radius = 0.5d0*(exp(ln10*radius(env_cz_zone))+exp(ln10*radius(env_cz_zone-1))) &
                    -cz_base_frac*delta_radius
            star%cz_base_radius_prev = cz_base_radius
            search_radius = cz_base_radius - star%pressure_scale_height_end
            do zone_idx = env_cz_zone-1,1,-1
               shell_radius = exp(ln10*radius(zone_idx))
               if(shell_radius.lt.search_radius)then
                  cz_base_zone = zone_idx + 1
                  delta_radius = exp(ln10*star%logR_start(zone_idx+1))-shell_radius
                  cz_base_frac = 0.5d0-((search_radius-shell_radius)/delta_radius)
                  cz_base_frac = max(-0.5d0,cz_base_frac)
                  cz_base_frac = min(0.5d0,cz_base_frac)
                  exit
               endif
            end do
            if (zone_idx .lt. 1) cz_base_zone = 1
         endif
      else
         cz_base_zone = env_cz_zone
         cz_base_zone_old = env_cz_zone_old
      endif
      star%envelope_cz_base_zone_prev = cz_base_zone
      return
end subroutine locate_cz_base

!----------------------------------------------------------------------
! average_li_be_over_cz: shared by burn_li_be and burn_li_be_midpoint (2026 readability
! R3; the block was byte-identical in both). Mass-averages the Li6,
! Li7, Be9 abundances and the log of the mass-averaged burning rates
! over the surface convection zone at the start of the step (from
! cz_base_zone_old) and at the end (from cz_base_zone, with the base
! zone's mass and rates adjusted by cz_base_frac so that it covers the
! whole CZ; shell_mass is restored before return, star%rate_li6/li7/be9
! at cz_base_zone keep the adjustment). Start-of-step log rates come
! from star%log_rate_*_prev when a previous step stored them.
subroutine average_li_be_over_cz(composition, mass_coordinate, shell_mass, &
     cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, &
     li6_cz_start, li7_cz_start, be9_cz_start, &
     log_rate_li6_cz_start, log_rate_li7_cz_start, log_rate_be9_cz_start, &
     li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end, &
     log_rate_li6_cz_end, log_rate_li7_cz_end, log_rate_be9_cz_end)
      use star_info_lib, only: star, json, i_li6, i_li7, i_be9
      use math_lib
      implicit none

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(inout) :: shell_mass(json)
      integer, intent(in) :: cz_base_zone, cz_base_zone_old, num_zones
      double precision, intent(in) :: cz_base_frac
      double precision, intent(out) :: li6_cz_start, li7_cz_start, be9_cz_start
      double precision, intent(out) :: log_rate_li6_cz_start, &
           log_rate_li7_cz_start, log_rate_be9_cz_start
      double precision, intent(out) :: li6_cz_end, li7_cz_end, be9_cz_end, &
           cz_mass_end
      double precision, intent(out) :: log_rate_li6_cz_end, &
           log_rate_li7_cz_end, log_rate_be9_cz_end

      integer :: zone_idx
      double precision :: cz_mass_start, shell_mass_save

! FIND RATES AT THE BEGINNING OF THE TIMESTEP (USING THE DEPTH AT THE START).
      li6_cz_start = 0.0d0
      li7_cz_start = 0.0d0
      be9_cz_start = 0.0d0
      cz_mass_start = 0.0d0
      do zone_idx = cz_base_zone_old,num_zones
         li6_cz_start = li6_cz_start+composition(i_li6,zone_idx)*shell_mass(zone_idx)
         li7_cz_start = li7_cz_start+composition(i_li7,zone_idx)*shell_mass(zone_idx)
         be9_cz_start = be9_cz_start+composition(i_be9,zone_idx)*shell_mass(zone_idx)
         cz_mass_start = cz_mass_start + shell_mass(zone_idx)
      end do
      li6_cz_start = li6_cz_start/cz_mass_start
      li7_cz_start = li7_cz_start/cz_mass_start
      be9_cz_start = be9_cz_start/cz_mass_start
! log_rate_li6_prev <= 0 means no previous step has stored CZ rates.
      if(star%log_rate_li6_prev.le.0.0d0)then
! COMPUTE MASS-WEIGHTED AVERAGE RATES AT THE START OF THE STEP.
         log_rate_li6_cz_start = 0.0d0
         log_rate_li7_cz_start = 0.0d0
         log_rate_be9_cz_start = 0.0d0
         do zone_idx = cz_base_zone_old,num_zones
            log_rate_li6_cz_start = log_rate_li6_cz_start + star%rate_li6_start(zone_idx)*shell_mass(zone_idx)
            log_rate_li7_cz_start = log_rate_li7_cz_start + star%rate_li7_start(zone_idx)*shell_mass(zone_idx)
            log_rate_be9_cz_start = log_rate_be9_cz_start + star%rate_be9_start(zone_idx)*shell_mass(zone_idx)
         end do
         log_rate_li6_cz_start = log(log_rate_li6_cz_start/cz_mass_start)
         log_rate_li7_cz_start = log(log_rate_li7_cz_start/cz_mass_start)
         log_rate_be9_cz_start = log(log_rate_be9_cz_start/cz_mass_start)
      else
! USE THE RATE FROM THE END OF THE PREVIOUS TIMESTEP.
         log_rate_li6_cz_start = star%log_rate_li6_prev
         log_rate_li7_cz_start = star%log_rate_li7_prev
         log_rate_be9_cz_start = star%log_rate_be9_prev
      endif
! USE THE LOCATION OF THE TRUE EDGE OF THE CONVECTION ZONE (FX, FOUND AT
! BEGINNING OF SR) TO ADJUST THE BURNING RATE AND MASS OF THE BOTTOM POINT
! SUCH THAT IT INCLUDES THE ENTIRE C.Z.
      shell_mass_save = shell_mass(cz_base_zone)
      if(cz_base_zone.gt.1.and.cz_base_zone.lt.num_zones)then
         shell_mass(cz_base_zone) = shell_mass(cz_base_zone)+cz_base_frac* &
              (mass_coordinate(cz_base_zone)-mass_coordinate(cz_base_zone-1))
         star%rate_li6(cz_base_zone) = star%rate_li6(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%rate_li6(cz_base_zone-1)-star%rate_li6(cz_base_zone))
         star%rate_li7(cz_base_zone) = star%rate_li7(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%rate_li7(cz_base_zone-1)-star%rate_li7(cz_base_zone))
         star%rate_be9(cz_base_zone) = star%rate_be9(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%rate_be9(cz_base_zone-1)-star%rate_be9(cz_base_zone))
      endif
! FIND RATES AT THE END OF THE TIMESTEP (USING THE DEPTH AT THE END).
! ALSO STORE INITIAL ABUNDANCES(FLI60,FLI70,FBE90).
! FM IS THE TOTAL MASS IN THE CZ.
      li6_cz_end = 0.0d0
      li7_cz_end = 0.0d0
      be9_cz_end = 0.0d0
      cz_mass_end = 0.0d0
      log_rate_li6_cz_end = 0.0d0
      log_rate_li7_cz_end = 0.0d0
      log_rate_be9_cz_end = 0.0d0
      do zone_idx = cz_base_zone,num_zones
         log_rate_li6_cz_end = log_rate_li6_cz_end + star%rate_li6(zone_idx)*shell_mass(zone_idx)
         log_rate_li7_cz_end = log_rate_li7_cz_end + star%rate_li7(zone_idx)*shell_mass(zone_idx)
         log_rate_be9_cz_end = log_rate_be9_cz_end + star%rate_be9(zone_idx)*shell_mass(zone_idx)
         li6_cz_end = li6_cz_end+composition(i_li6,zone_idx)*shell_mass(zone_idx)
         li7_cz_end = li7_cz_end+composition(i_li7,zone_idx)*shell_mass(zone_idx)
         be9_cz_end = be9_cz_end+composition(i_be9,zone_idx)*shell_mass(zone_idx)
         cz_mass_end = cz_mass_end + shell_mass(zone_idx)
      end do
      li6_cz_end = li6_cz_end/cz_mass_end
      li7_cz_end = li7_cz_end/cz_mass_end
      be9_cz_end = be9_cz_end/cz_mass_end
! MASS-WEIGHTED AVERAGE RATES AT THE END OF THE STEP.
      log_rate_li6_cz_end = log(log_rate_li6_cz_end/cz_mass_end)
      log_rate_li7_cz_end = log(log_rate_li7_cz_end/cz_mass_end)
      log_rate_be9_cz_end = log(log_rate_be9_cz_end/cz_mass_end)
! RESTORE ORIGINAL MASS TO BASE OF C.Z.
      shell_mass(cz_base_zone) = shell_mass_save
      return
end subroutine average_li_be_over_cz



!----------------------------------------------------------------------
! burn_li_be
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original liburn.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Determine lithium-6, lithium-7, and beryllium-9 burning in a
! non-rotating model (see burn_li_be_midpoint below for the simpler single-pass
! rotating-model variant). The burning rates depend on the local T and
! rho, the abundance, and (in convection zones) the mixing.
!
! In radiative regions burning is done implicitly, correcting the
! rates for the change in abundance as a function of time, by
! repeatedly halving the sub-step size and using a rational-function
! extrapolation (ratext) to converge on the zero-step-size limit.
!
! In the convective region the program iterates between burning and
! mixing in the same way, to correct the rates for the change in
! abundance as a function of time.
!
! 11/91 HR added to call.
!
! shell_mass genuinely needs intent(inout), not intent(in): it's
! temporarily perturbed for the CZ-base mass adjustment
! (shell_mass_save = shell_mass(cz_base_zone) ... shell_mass(cz_base_zone)
! = shell_mass_save) and the perturbed value is read by the
! immediately-following rate/abundance sums before being restored, so
! the write is real, not dead code. Externally the net effect on the
! caller's array is still zero (every path that writes also restores
! before return), which is why rotation/evolve_angular_momentum.f90 --
! one of this routine's two callers -- declares its own shell_mass
! intent(inout) to match this explicit interface.
! burn_li_be (formerly liburn)
subroutine burn_li_be(timestep, composition, radius, mass_coordinate, &
     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
      use rotation_scratch_lib
      use star_info_lib, only: star, json, i_li6, i_li7, i_be9
      use luout_lib
      use phys_const_lib
      use numerics_lib
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
      double precision :: light_element_save(3,json)
      integer :: substep_counts(11)
      double precision :: extrap_tol(3), extrap_y(3), extrap_err(3), &
           extrap_result(3)
      data substep_counts/2,4,6,8,12,16,24,32,48,64,96/
      data extrap_tol/1.0d-6,1.0d-6,1.0d-6/
      integer :: zone_idx, refine_idx, substep_idx, species_idx
      integer :: cz_base_zone, cz_base_zone_old, min_zone, max_zone
      double precision :: cz_base_frac
      double precision :: substep_frac, substep_dt
      logical :: converged
      double precision :: log_rate_li6, log_rate_li7, log_rate_be9
      double precision :: extrap_x
      double precision :: li6_cz_start, li7_cz_start, be9_cz_start
      double precision :: log_rate_li6_cz_start, log_rate_li7_cz_start, &
           log_rate_be9_cz_start
      double precision :: li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end
      double precision :: log_rate_li6_cz_end, log_rate_li7_cz_end, &
           log_rate_be9_cz_end
      logical :: accretion_active
      double precision :: li6_accreted, li7_accreted, be9_accreted
      double precision :: cz_log_rate_li6, cz_log_rate_li7, cz_log_rate_be9
      double precision :: li6_depletion, li7_depletion, be9_depletion
      double precision :: li6_added, li7_added, be9_added
      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
           radiative_frac

! SAVE ORIGINAL ABUNDANCES.
      do zone_idx = 1,num_zones
         light_element_save(1,zone_idx) = composition(i_li6,zone_idx)
         light_element_save(2,zone_idx) = composition(i_li7,zone_idx)
         light_element_save(3,zone_idx) = composition(i_be9,zone_idx)
      end do
! THE DEGREE OF LITHIUM BURNING IN A SURFACE CZ DEPENDS SENSITIVELY
! ON THE TEMPERATURE AT ITS BASE - SO ACCURATELY LOCATING IS IMPORTANT.
! DETERMINE THE TRUE LOCATION (FX) OF THE BASE OF THE CZ AT THE END
! OF THE TIMESTEP, AND THE LOCATION OF THE EDGE OF OVERSHOOT REGIONS
! IF APPLICABLE.
      call locate_cz_base(radius, env_cz_zone, env_cz_zone_old, num_zones, &
           cz_base_zone, cz_base_zone_old, cz_base_frac)
! RADIATIVE INTERIOR.
      min_zone = min(cz_base_zone,cz_base_zone_old)
      max_zone = max(cz_base_zone,cz_base_zone_old)
      do zone_idx = 1,min_zone-1
         if(star%rate_be9(zone_idx).le.1.0d-32 .or. star%rate_be9_start(zone_idx).le.1.0d-32)exit
         if(composition(i_li6,zone_idx).lt.1.0d-24.and.composition(i_li7,zone_idx).lt.1.0d-24 .and.composition(i_be9,zone_idx).lt.1.0d-24)cycle
         if(log_temperature(zone_idx).gt.7.0d0)then
            composition(i_li6,zone_idx) = 0.0d0
            composition(i_li7,zone_idx) = 0.0d0
            composition(i_be9,zone_idx) = 0.0d0
            cycle
         endif
         do refine_idx = 1,11
! PERFORM BURNING IN substep_counts SMALLER TIMESTEPS.
!  DEFINE FRAC - THE FRACTION OF THE TOTAL TIMESTEP IN EACH SMALL ONE -
!  AND DELT - THE FRACTIONAL TIMESTEP IN YEARS.
            substep_frac = 1.0d0/dfloat(substep_counts(refine_idx))
            substep_dt = substep_frac*timestep
            converged = .false.
! RESTORE INITIAL ABUNDANCES.
            composition(i_li6,zone_idx) = light_element_save(1,zone_idx)
            composition(i_li7,zone_idx) = light_element_save(2,zone_idx)
            composition(i_be9,zone_idx) = light_element_save(3,zone_idx)
! STORE STARTING REACTION RATES.
            log_rate_li6 = log(star%rate_li6_start(zone_idx))-0.5d0*substep_frac* &
                   (log(star%rate_li6(zone_idx))-log(star%rate_li6_start(zone_idx)))
            log_rate_li7 = log(star%rate_li7_start(zone_idx))-0.5d0*substep_frac* &
                   (log(star%rate_li7(zone_idx))-log(star%rate_li7_start(zone_idx)))
            log_rate_be9 = log(star%rate_be9_start(zone_idx))-0.5d0*substep_frac* &
                   (log(star%rate_be9(zone_idx))-log(star%rate_be9_start(zone_idx)))
            do substep_idx = 1,substep_counts(refine_idx)
! INCREMENT THE REACTION RATES.
               log_rate_li6 = log_rate_li6+substep_frac* &
                    (log(star%rate_li6(zone_idx))-log(star%rate_li6_start(zone_idx)))
               log_rate_li7 = log_rate_li7+substep_frac* &
                    (log(star%rate_li7(zone_idx))-log(star%rate_li7_start(zone_idx)))
               log_rate_be9 = log_rate_be9+substep_frac* &
                    (log(star%rate_be9(zone_idx))-log(star%rate_be9_start(zone_idx)))
               li6_substep_depletion = substep_dt*exp(log_rate_li6)
               li7_substep_depletion = substep_dt*exp(log_rate_li7)
               be9_substep_depletion = substep_dt*exp(log_rate_be9)
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
            end do
! STORE ABUNDANCES AS A FUNCTION OF TIMESTEP IN VECTOR YEST.
            extrap_y(1) = composition(i_li6,zone_idx)
            extrap_y(2) = composition(i_li7,zone_idx)
            extrap_y(3) = composition(i_be9,zone_idx)
            extrap_x = substep_frac**2
! USE A RATIONAL POLYNOMIAL EXTRAPOLATOR TO EXTRAPOLATE TO THE SOLUTION FOR
! ZERO TIMESTEP.
            call ratext(refine_idx,extrap_x,extrap_y,extrap_result,extrap_err,3,7)
! CHECK IF THE SOLUTION HAS CONVERGED.
            if(refine_idx.gt.1)then
               converged = .true.
               do species_idx=1,3
                  if(extrap_y(species_idx).lt.1.0d-24)cycle
                  extrap_err(species_idx) = abs(extrap_err(species_idx)/extrap_y(species_idx))
                  if(extrap_err(species_idx).gt.extrap_tol(species_idx))converged=.false.
               end do
               if(converged)then
!                 WRITE(ISHORT,912)I,J,(YEXT(K2)K2=1,3)
!   912             FORMAT(1X,'CONVERGED',I5,' LEVEL ',I2,' LI6 ',1P,
!      *            E12.5,' LI7 ',E12.5,' BE9 ',E12.5)
                  exit
               endif
            endif
         end do
         if (.not. converged) then
! THE REFINEMENT LOOP RAN OUT WITHOUT CONVERGING TO WITHIN THE
! SPECIFIED TOLERANCE IN THE MAXIMUM NUMBER OF ITERATIONS.
            write(run_log_unit,911)zone_idx,(extrap_err(species_idx),species_idx=1,3)
  911       format(1x,'***LIBURN CONVERGENCE FAILURE IN SHELL ',i4, &
            'ERRORS '/1p3e10.3)
         end if
! WRITE NEW ABUNDANCES AND EXIT.
         composition(i_li6,zone_idx)=extrap_result(1)
         if(composition(i_li6,zone_idx).lt.1.0d-24)composition(i_li6,zone_idx)=0.0d0
         composition(i_li7,zone_idx)=extrap_result(2)
         if(composition(i_li7,zone_idx).lt.1.0d-24)composition(i_li7,zone_idx)=0.0d0
         composition(i_be9,zone_idx)=extrap_result(3)
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
           li6_cz_start, li7_cz_start, be9_cz_start, &
           log_rate_li6_cz_start, log_rate_li7_cz_start, log_rate_be9_cz_start, &
           li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end, &
           log_rate_li6_cz_end, log_rate_li7_cz_end, log_rate_be9_cz_end)
! FOR THE STARTING ABUNDANCE, USE THE MIXED ABUNDANCE FOR THE
! DEEPER CZ.
      if(cz_base_zone.lt.cz_base_zone_old)then
         li6_cz_start = li6_cz_end
         li7_cz_start = li7_cz_end
         be9_cz_start = be9_cz_end
      endif
      do refine_idx = 1,11
! PERFORM BURNING IN substep_counts SMALLER TIMESTEPS.
!  DEFINE FRAC - THE FRACTION OF THE TOTAL TIMESTEP IN EACH SMALL ONE -
!  AND DELT - THE FRACTIONAL TIMESTEP IN YEARS.
         substep_frac = 1.0d0/dfloat(substep_counts(refine_idx))
         substep_dt = substep_frac*timestep
         converged = .false.
! INITIALIZE ABUNDANCES.
         li6_cz_end = li6_cz_start
         li7_cz_end = li7_cz_start
         be9_cz_end = be9_cz_start
! INITIALIZE RATES.
         cz_log_rate_li6 = log_rate_li6_cz_start-0.5d0*substep_frac* &
              (log_rate_li6_cz_end-log_rate_li6_cz_start)
         cz_log_rate_li7 = log_rate_li7_cz_start-0.5d0*substep_frac* &
              (log_rate_li7_cz_end-log_rate_li7_cz_start)
         cz_log_rate_be9 = log_rate_be9_cz_start-0.5d0*substep_frac* &
              (log_rate_be9_cz_end-log_rate_be9_cz_start)
! BURN IN substep_counts SMALL TIMESTEPS.
! MHP 05/02 ADDED THE IMPACT OF MASS ACCRETION ON LIGHT
! ELEMENT ABUNDANCES.  INITIALIZE ARRAYS; SPECIES ARE
! ASSUMED TO EXPERIENCE BURNING FOR ON AVERAGE 1/2 OF
! THE STEP WHERE THEY ARE INITIALLY ACCRETED AND ARE
! THEN FULLY BURNED.
         if(star%job%use_mass_accretion.and.star%ctrl%mass_accretion_rate.gt.0.0d0)then
            accretion_active = .true.
            li6_accreted = 0.0d0
            li7_accreted = 0.0d0
            be9_accreted = 0.0d0
         else
            accretion_active = .false.
         endif
         do substep_idx = 1,substep_counts(refine_idx)
! INCREMENT RATES BY LINEAR INTERPOLATION IN THE LOG.
            cz_log_rate_li6 = cz_log_rate_li6 + substep_frac* &
                 (log_rate_li6_cz_end-log_rate_li6_cz_start)
            cz_log_rate_li7 = cz_log_rate_li7 + substep_frac* &
                 (log_rate_li7_cz_end-log_rate_li7_cz_start)
            cz_log_rate_be9 = cz_log_rate_be9 + substep_frac* &
                 (log_rate_be9_cz_end-log_rate_be9_cz_start)
            li6_depletion = substep_dt*exp(cz_log_rate_li6)
            li7_depletion = substep_dt*exp(cz_log_rate_li7)
            be9_depletion = substep_dt*exp(cz_log_rate_be9)
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
! MHP 05/02 MASS ACCRETION.  THE NET LITHIUM ADDED IN
! A GIVEN STEP IS XXXADD WHILE THE LITHIUM ADDED IN PRIOR
! STEPS IS XXXA; THE TWO ARE ADDED AFTER BURNING.
! N.B. SINCE THIS IS A HALF-LIFE PROBLEM IT IS OK TO BURN
! THE COMPONENTS SEPARATELY.
            if(accretion_active)then
               li6_added = star%ctrl%accreted_composition(i_li6)*substep_frac/exp(0.5d0*li6_depletion)
               li6_accreted = li6_accreted/exp(li6_depletion) + li6_added
               li7_added = star%ctrl%accreted_composition(i_li7)*substep_frac/exp(0.5d0*li7_depletion)
               li7_accreted = li7_accreted/exp(li7_depletion) + li7_added
               be9_added = star%ctrl%accreted_composition(i_be9)*substep_frac/exp(0.5d0*be9_depletion)
               be9_accreted = be9_accreted/exp(be9_depletion) + be9_added
            endif
         end do
! DO A MASS-WEIGHTED AVERAGE OF THE ORIGINAL LIGHT
! ELEMENT CONTENT AND THE NET AMOUNT ADDED.
! FMASSACC = DMDT*DT/ORIGINAL CZ MASS
! NOTE: THIS FORMULATION ASSUMES THAT ALL ACCRETED MATTER
! LANDS IN A CZ WITH THE END-OF-TIMESTEP DEPTH.
! NOT STRICTLY TRUE, BUT NOT A BAD APPROXIMATION EITHER.
         if(accretion_active)then
            write(*,913)li6_cz_start,li6_cz_end,li6_accreted,li7_cz_start, &
                 li7_cz_end,li7_accreted,star%accreted_mass
 913        format(1p7e12.3)
            li6_cz_end = (li6_cz_end*cz_mass_end+li6_accreted*star%accreted_mass)/ &
                 (cz_mass_end+star%accreted_mass)
            li7_cz_end = (li7_cz_end*cz_mass_end+li7_accreted*star%accreted_mass)/ &
                 (cz_mass_end+star%accreted_mass)
            be9_cz_end = (be9_cz_end*cz_mass_end+be9_accreted*star%accreted_mass)/ &
                 (cz_mass_end+star%accreted_mass)
         endif
! STORE THE ABUNDANCE AS A FUNCTION OF THE TIME STEP.
! STORE ABUNDANCES AS A FUNCTION OF TIMESTEP IN VECTOR YEST.
         extrap_y(1) = li6_cz_end
         extrap_y(2) = li7_cz_end
         extrap_y(3) = be9_cz_end
         extrap_x = substep_frac**2
! USE A RATIONAL POLYNOMIAL EXTRAPOLATOR TO EXTRAPOLATE TO THE SOLUTION FOR
! ZERO TIMESTEP.
         call ratext(refine_idx,extrap_x,extrap_y,extrap_result,extrap_err,3,7)
! CHECK IF THE SOLUTION HAS CONVERGED.
         if(refine_idx.gt.1)then
            converged = .true.
            do species_idx=1,3
               if(extrap_y(species_idx).lt.1.0d-24)cycle
               extrap_err(species_idx) = abs(extrap_err(species_idx)/extrap_y(species_idx))
               if(extrap_err(species_idx).gt.extrap_tol(species_idx))converged=.false.
            end do
            if(converged)then
!              WRITE(ISHORT,912)JENV,J,(YEXT(K2),YERR(K2),K2=1,3)
               exit
            endif
         endif
      end do
      if (.not. converged) then
! THE REFINEMENT LOOP RAN OUT WITHOUT CONVERGING TO WITHIN THE
! SPECIFIED TOLERANCE IN THE MAXIMUM NUMBER OF ITERATIONS.
         write(run_log_unit,911)cz_base_zone,(extrap_err(species_idx),species_idx=1,3)
      end if
! WRITE NEW ABUNDANCES AND EXIT.
      li6_cz_end = extrap_result(1)
      if(li6_cz_end.lt.1.0d-24)li6_cz_end=0.0d0
      li7_cz_end = extrap_result(2)
      if(li7_cz_end.lt.1.0d-24)li7_cz_end=0.0d0
      be9_cz_end = extrap_result(3)
      if(be9_cz_end.lt.1.0d-24)be9_cz_end=0.0d0
      do zone_idx = max_zone,num_zones
         composition(i_li6,zone_idx) = li6_cz_end
         composition(i_li7,zone_idx) = li7_cz_end
         composition(i_be9,zone_idx) = be9_cz_end
      end do
! STORE ENDING RATE FOR USE AT THE BEGINNING OF THE NEXT STEP.
      star%log_rate_li6_prev = log_rate_li6_cz_end
      star%log_rate_li7_prev = log_rate_li7_cz_end
      star%log_rate_be9_prev = log_rate_be9_cz_end
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
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION: Li6 and Li7 go
! through safe_divide_by_exp (overflow-safe divide by exp); Be9 here, and all
! three species in burn_li_be_midpoint's twin of this loop, still use the plain
! division. Kept as is (changing it could change numbers).
         call safe_divide_by_exp(composition(i_li6,zone_idx),li6_depletion)
         if(composition(i_li6,zone_idx).lt.1.0d-24)composition(i_li6,zone_idx)=0.0d0
         call safe_divide_by_exp(composition(i_li7,zone_idx),li7_depletion)
         if(composition(i_li7,zone_idx).lt.1.0d-24)composition(i_li7,zone_idx)=0.0d0
         composition(i_be9,zone_idx) = composition(i_be9,zone_idx)/exp(be9_depletion)
         if(composition(i_be9,zone_idx).lt.1.0d-24)composition(i_be9,zone_idx)=0.0d0
      end do
      return
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
! shell_mass genuinely needs intent(inout) -- same reasoning as
! burn_li_be's identical fix above: the write (a save-then-restore pair
! around the CZ-base mass adjustment) is read by the following
! rate/abundance sums before being restored, so it's real, not dead
! code, even though the net effect on the caller's array is zero.
! burn_li_be_midpoint (formerly liburn2)
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
      call locate_cz_base(radius, env_cz_zone, env_cz_zone_old, num_zones, &
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
           li6_cz_start, li7_cz_start, be9_cz_start, &
           log_rate_li6_cz_start, log_rate_li7_cz_start, log_rate_be9_cz_start, &
           li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end, &
           log_rate_li6_cz_end, log_rate_li7_cz_end, log_rate_be9_cz_end)
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
! STORE ENDING RATE FOR USE AT THE BEGINNING OF THE NEXT STEP.
      star%log_rate_li6_prev = log_rate_li6_cz_end
      star%log_rate_li7_prev = log_rate_li7_cz_end
      star%log_rate_be9_prev = log_rate_be9_cz_end
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
subroutine store_start_li_be_rates()
      use star_info_lib, only: star
      call store_li_be_rates(star%xa, star%logRho, star%logT, star%nz, 2)
end subroutine store_start_li_be_rates

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
!  rate_by_zone(rr_*,k) net_lib's rate vectors, per Gyr per amu
!  shell_mass          the mass weights for a convection zone
!  composition         out: the new abundances (star%xa layout); the
!                      solved shell's previous neighbour is the first
!                      guess, so shells are solved in increasing order
!  timestep_years      the step
! The start-of-step abundances are star%xa_start. Z changes only once
! hydrogen is exhausted (X < 5e-7 at the top of the region), as
! 1 - X - Y - He3.
subroutine burn_network_region(log_temperature, zone_begin, zone_end, &
     rate_by_zone, shell_mass, composition, timestep_years, ierr)

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
      double precision, intent(in) :: timestep_years
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

      call implicit_burn_step(x_start, rate_vec, timestep_years*1.0d-9, &
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
! The burning steps of mix (2026-10: moved here from mixing/mix.f90,
! which calls them in its original order; bodies unchanged).
! =================================================================

!----------------------------------------------------------------------
! compute_burn_rates
!----------------------------------------------------------------------
! The composition rate vector of every shell above the nuclear cutoff
! (net_lib's network_rate_vector; zero above it) and its deuterium burning rate
! (light_element_lib's deuterium_rate, stored in
! star%deuterium_burning_rate and, at iteration level 1, _start). The
! rates are also kept in rot_scr%reaction_rate_by_zone for rotmix and
! rezone.
subroutine compute_burn_rates(iteration_level, rate_by_zone)
      use star_info_lib, only: star, json, i_h1, i_h2
      use rotation_scratch_lib
      use light_element_lib, only: deuterium_rate
      integer, intent(in) :: iteration_level
      double precision, intent(out) :: rate_by_zone(num_rate_rows, json)
      double precision :: deuterium_test
      integer :: zone_idx, clear_idx

! FIND BURNING RATES (HR1- HR13,HF1,HF2).
      if (star%job%use_mass_accretion .and. star%ctrl%mass_accretion_rate.gt.0.0d0) then
         deuterium_test = max(star%xa(i_h2,star%nz), &
              star%ctrl%accreted_composition(12))
      else
         deuterium_test = star%xa(i_h2,star%nz)
      end if
      do zone_idx = 1, star%nz
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
         if (star%logT(zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
! SETUP NUCLEAR ENERGY TERMS
         call network_rate_vector(star%logRho(zone_idx), star%logT(zone_idx), &
              star%xa(:,zone_idx), rate_by_zone(:, zone_idx))
! MHP 5/02 COMPUTE RATE FOR DEUTERIUM BURNING
         if (deuterium_test.gt.1.0d-11) then
            star%deuterium_burning_rate(zone_idx) = deuterium_rate( &
                 star%logRho(zone_idx), star%logT(zone_idx), star%xa(i_h1,zone_idx), &
                 zone_idx.ge.star%jcz, star%convective_turnover_timescale)
            if (iteration_level.eq.1) star%deuterium_burning_rate_start(zone_idx) = &
                 star%deuterium_burning_rate(zone_idx)
         else
            star%deuterium_burning_rate(zone_idx) = 0.0d0
         end if
      end do
! ZERO THE RATES BELOW THE CUTOFF (zone_idx = nz+1 IF NO SHELL WAS COLD).
      do clear_idx = zone_idx, star%nz
         rate_by_zone(:, clear_idx) = 0.0d0
! MHP 5/02 ZERO OUT DEUTERIUM BURNING RATE
         star%deuterium_burning_rate(clear_idx) = 0.0d0
      end do
! keep the rates for rotmix/rezone (columns 1..nz; the rows are the
! same rr_* layout, so each column is copied whole).
      do zone_idx = 1, star%nz
         rot_scr%reaction_rate_by_zone(:,zone_idx) = rate_by_zone(:,zone_idx)
      end do
end subroutine compute_burn_rates

!----------------------------------------------------------------------
! burn_network
!----------------------------------------------------------------------
! The implicit network burn over the step (burn_network_region, i.e.
! network_lib's implicit_burn_step): shell by shell in the radiative
! regions above the nuclear cutoff, and once per mixed region, which it
! also homogenizes.
subroutine burn_network(radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, rate_by_zone, timestep_years, ierr)
      use star_info_lib, only: star, json, max_radiative_zones
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
      double precision, intent(in) :: timestep_years
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
                 rate_by_zone, star%dm, star%xa, timestep_years, ierr)
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
              rate_by_zone, star%dm, star%xa, timestep_years, ierr)
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
subroutine burn_equilibrium(radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, rate_by_zone, timestep_years)
      use star_info_lib, only: star, json, max_radiative_zones, i_h1, &
           i_metals, i_c12, i_o16
      use network_lib, only: equilibrium_burn_rates
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
      double precision, intent(in) :: timestep_years
      integer :: radiative_region_idx, inner_zone_idx, mixed_zone_idx, &
           zone_begin, zone_end, zone_idx
      double precision :: dt_gyr, dc_dt, do_dt, dx_dt, dy_dt
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
         dt_gyr = timestep_years*1.0d-9
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
                       dx_dt*dt_gyr
               end if
!  USE THE HELIUM BURNING RATE FROM EQBURN AND THE CARBON,ALPHA
!  BURNING RATE
               if (dy_dt.ne.0.0d0) then
                  star%xa(i_metals,inner_zone_idx) = metal_fraction - &
                       dy_dt*dt_gyr
                  star%xa(i_c12,inner_zone_idx) = c12_fraction + &
                       dc_dt*dt_gyr
                  star%xa(i_o16,inner_zone_idx) = o16_fraction + &
                       do_dt*dt_gyr
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
               hydrogen_fraction = hydrogen_fraction + dx_dt*dt_gyr
               do zone_idx = zone_begin, zone_end
                  star%xa(i_h1,zone_idx) = hydrogen_fraction
               end do
            end if
!  USE THE HELIUM BURNING RATE FROM EQBURN AND THE CARBON,ALPHA
!  BURNING RATE
            if (dy_dt.ne.0.0d0) then
               metal_fraction = metal_fraction - dy_dt*dt_gyr
               c12_fraction = c12_fraction + dc_dt*dt_gyr
               o16_fraction = o16_fraction + do_dt*dt_gyr
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
subroutine burn_deuterium(radiative_zone_bounds, num_radiative_zones, &
     num_mixed_zones, timestep_years)
      use star_info_lib, only: star, max_radiative_zones
      use light_element_lib, only: burn_deuterium_region, min_deuterium_step
      integer, intent(in) :: radiative_zone_bounds(max_radiative_zones,2), &
           num_radiative_zones, num_mixed_zones
      double precision, intent(in) :: timestep_years
      integer :: radiative_region_idx, inner_zone_idx, mixed_zone_idx, &
           zone_begin, zone_end
      double precision :: dt_gyr

! MHP 5/02 ADDED DEUTERIUM BURNING
! RADIATIVE ZONES.
!
         dt_gyr = timestep_years*1.0d-9
         do radiative_region_idx = 1, num_radiative_zones
            do inner_zone_idx = &
                 radiative_zone_bounds(radiative_region_idx,1), &
                 radiative_zone_bounds(radiative_region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (star%logT(inner_zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
               zone_begin = inner_zone_idx
               zone_end = inner_zone_idx
               call burn_deuterium_region(zone_begin, zone_end, star%nz, star%dm, &
                    star%xa, dt_gyr, star%xa_start, star%deuterium_burning_rate_start, &
                    star%deuterium_burning_rate, star%accreted_mass, min_deuterium_step)
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
                 star%xa, dt_gyr, star%xa_start, star%deuterium_burning_rate_start, &
                 star%deuterium_burning_rate, star%accreted_mass, min_deuterium_step)
         end do
end subroutine burn_deuterium

end module burn_step
