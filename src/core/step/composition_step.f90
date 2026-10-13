!----------------------------------------------------------------------
! composition_step (formerly mix, mixing/mix.f90)
!----------------------------------------------------------------------
! One operator-split composition step over the timestep, at fixed
! structure. It changes star%xa (and the mixed-region bounds); the
! structure is solved afterwards with the new composition.
!
! evolve_step runs it in two kinds of pass (the two named options
! below):
!   - predictor_pass, once per timestep from advance_composition_and_age:
!     on the start-of-step structure, before the structure solve (the
!     composition that iteration levels 1-2 solve with);
!   - corrector_pass, at iteration levels >= 3 from solve_level: star%xa
!     is reset to star%xa_start and the step is redone with the
!     updated structure.
!
! The order of operations (with the options that switch each one):
!   1. the mixed regions: convection zones with overshoot
!      (find_convection_zones), with the ctrl%dpenv deep-mixing option
!      (deep_mixing)
!   2. the burning rates of every shell (compute_burn_rates; with
!      store_start_rates, the deuterium rates are also kept as the
!      start-of-step rates)
!   3. the implicit network burn (burn_network)
!   4. semiconvection (semiconvection, when job%use_semiconvection is set)
!   5. homogenize the mixed regions (homogenize, or always with
!      semiconvection on)
!   6. gravitational settling (when enabled), sub-stepped
!   7. clip the abundances to [0, 1] and close the mass fractions
!   8. after deep mixing, recompute the zone bounds with the ordinary
!      flags (deep_mixing)
!   9. deuterium burning (extended composition)
! (A deviation from the legacy code, 2026-10: the legacy predictor pass
! also ran an explicit burn with He3 and the CN cycle in equilibrium
! after the network, overwriting X and, in helium burning, Z, C12 and
! O16 with a cruder estimate of the same step. Both passes now burn the
! same way.)
! Li/Be are burned later, after the structure converges (burn_li_be).
module composition_step
      use star_info_lib, only: star_info
      implicit none
      private
      public :: advance_composition

! What a pass does beyond the operations every pass runs.
      type, public :: composition_step_options
! Keep this pass's deuterium burning rates as the start-of-step rates.
         logical :: store_start_rates = .false.
! Homogenize the mixed regions (always done when semiconvection is on).
         logical :: homogenize = .false.
! Run semiconvection, when job%use_semiconvection is set.
         logical :: semiconvection = .false.
! Apply the deep-mixing option, when ctrl%dpenv < 1.
         logical :: deep_mixing = .false.
      end type composition_step_options

      type(composition_step_options), parameter, public :: &
           predictor_pass = composition_step_options(store_start_rates=.true., &
                homogenize=.true., semiconvection=.false., deep_mixing=.false.), &
           corrector_pass = composition_step_options(store_start_rates=.false., &
                homogenize=.false., semiconvection=.true., deep_mixing=.true.)

contains

!----------------------------------------------------------------------
! advance_composition
!----------------------------------------------------------------------
! The composition step over timestep (s), as the pass opts selects (the
! operations in the module header). It changes star%xa and the
! mixed-region bounds, and returns the edges of the convective core and
! of the surface convection zone (core_cz_edge, envelope_cz_edge) and the
! mixed-region bounds without overshoot. ierr /= 0 on a failure in the
! zone search, the network burn, semiconvection or settling.
!   star reads: ctrl%dpenv, ctrl%helium_diffusion_min,
!     ctrl%hydrogen_diffusion_floor, ctrl%settling_timestep_fraction,
!     dm, gradT, job%use_extended_composition, job%use_semiconvection,
!     log_mass, log_Teff, log_total_mass, logP, logR, logT,
!     luminosity_lsun, nz
!   star writes: convective_flag, iov1, iov2, job%diffuse_helium_active,
!     logRho, m, mixed_zone_bounds, xa
subroutine advance_composition(s, timestep, opts, core_cz_edge, &
     envelope_cz_edge, mixed_zone_bounds_no_overshoot, ierr)
      use star_info_lib, only: i_h1, i_he3, i_he4, i_metals, json, &
           max_convective_zones, max_radiative_zones
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use net_lib, only: num_rate_rows
      use burn_step, only: compute_burn_rates, burn_network, &
           burn_deuterium
      use convection_lib, only: find_convection_zones, extend_by_semiconvection
      use element_diffusion_lib, only: diffuse_elements
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

! timestep: s
      double precision, intent(in) :: timestep
      type(composition_step_options), intent(in) :: opts
      integer, intent(out) :: core_cz_edge, envelope_cz_edge
      integer, intent(inout) :: mixed_zone_bounds_no_overshoot(max_convective_zones,2)
      integer, intent(out) :: ierr

! rate_by_zone(r_*/rv_*, zone): the composition rate vector of every
! shell, from compute_burn_rates.
      double precision :: rate_by_zone(num_rate_rows, json)
      integer :: radiative_zone_bounds(max_radiative_zones,2)
      logical :: deep_mix_flag(json)
      integer :: num_species
      integer :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      integer :: zone_idx, inner_zone_idx

      ierr = 0

      if (s%job%use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if

! 1. The mixed regions: the convection zones with and without overshoot,
! from the mixing flags of set_deep_mix_flags.
      call set_deep_mix_flags
      call find_convection_zones(s%xa, s%logRho, s%logP, s%logR, &
           s%log_mass, s%logT, deep_mix_flag, s%nz, &
           radiative_zone_bounds, s%mixed_zone_bounds, &
           mixed_zone_bounds_no_overshoot, core_cz_edge, envelope_cz_edge, &
           num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot, s%iov1, s%iov2, ierr)
      if (ierr /= 0) return

! 2-3. Burning: the rates, then the implicit network.
      call compute_burn_rates(s, opts%store_start_rates, rate_by_zone)
      call burn_network(s, radiative_zone_bounds, num_radiative_zones, &
           num_mixed_zones, rate_by_zone, timestep, ierr)
      if (ierr /= 0) return

! 4. Semiconvection, on the passes that run it.
      if (s%job%use_semiconvection) then
         if (opts%semiconvection) &
              call extend_by_semiconvection(timestep, s%xa, s%logRho, &
              s%luminosity_lsun, s%logP, s%logR, s%log_mass, &
              s%logT, s%nz, s%mixed_zone_bounds, &
              num_mixed_zones, s%log_Teff, ierr)
              if (ierr /= 0) return
      end if

! 5. Homogenize the mixed regions: always with semiconvection on, and on
! the passes that ask for it (opts%homogenize). The bounds go to the run
! log with the solver diagnostics.
      if (s%job%use_semiconvection .or. opts%homogenize) &
           call homogenize_mixed_regions
      if (num_mixed_zones.ge.1 .and. solver_diagnostics()) then
         write(run_log_unit,110) ((s%mixed_zone_bounds(zone_idx,inner_zone_idx), &
              inner_zone_idx=1,2),zone_idx=1,num_mixed_zones)
  110    format(' ZONES MIXED IN ORDER--',12('(',i5,',',i5,') ') )
      end if

! 6. Gravitational settling.
      if (s%job%diffuse_helium_active) then
         call settle
         if (ierr /= 0) return
      end if

! 7. Clip and renormalize.
      call clip_and_renormalize

! 8. After deep mixing (MHP 1/95), recompute the zone bounds and the
! edges (core_cz_edge, envelope_cz_edge; legacy JCORE, JENV) from the
! ordinary convective flags.
      if (s%ctrl%dpenv.lt.1.0d0 .and. opts%deep_mixing) then
         call find_convection_zones(s%xa, s%logRho, s%logP, s%logR, &
              s%log_mass, s%logT, s%convective_flag, s%nz, &
              radiative_zone_bounds, s%mixed_zone_bounds, &
              mixed_zone_bounds_no_overshoot, core_cz_edge, &
              envelope_cz_edge, num_radiative_zones, num_mixed_zones, &
              num_mixed_zones_no_overshoot, s%iov1, s%iov2, ierr)
         if (ierr /= 0) return
      end if

! 9. Deuterium burning over the step (extended composition only).
      if (s%job%use_extended_composition) call burn_deuterium(s,  &
           radiative_zone_bounds, num_radiative_zones, num_mixed_zones, &
           timestep)

contains

! deep_mix_flag: the mixing flag of each shell. ctrl%dpenv < 1 (with
! deep_mixing) forces mixing from the centre out to the mass fraction
! dpenv, m <= dpenv M; outside that, and otherwise, each shell keeps its
! ordinary convective flag.
subroutine set_deep_mix_flags
      double precision :: max_mixed_mass
      integer :: k, j
      if (s%ctrl%dpenv.lt.1.0d0 .and. opts%deep_mixing) then
         max_mixed_mass = s%ctrl%dpenv*exp(ln10*s%log_total_mass)
         do k = 1, s%nz
            if (s%m(k).gt.max_mixed_mass) then
               do j = k, s%nz
                  deep_mix_flag(j) = s%convective_flag(j)
               end do
               exit
            else
               deep_mix_flag(k) = .true.
            end if
         end do
      else
         do j = 1, s%nz
            deep_mix_flag(j) = s%convective_flag(j)
         end do
      end if
end subroutine set_deep_mix_flags

! The mass-weighted mean of every species, X_i = sum(X_i dm)/sum(dm),
! applied throughout each mixed region (star%mixed_zone_bounds: first
! and last shell). A one-shell region is skipped unless it starts at the
! centre.
subroutine homogenize_mixed_regions
      double precision :: species_sum(15), weight_sum
      integer :: mixed_zone_idx, mix_start, mix_end, shell_idx, species_idx
      do mixed_zone_idx = 1, num_mixed_zones
         mix_start = s%mixed_zone_bounds(mixed_zone_idx,1)
         mix_end = s%mixed_zone_bounds(mixed_zone_idx,2)
         if (mix_start.ne.1 .and. mix_start.ge.mix_end) cycle
         weight_sum = 0.0d0
         do species_idx = 1, num_species
            species_sum(species_idx) = 0.0d0
         end do
         do shell_idx = mix_start, mix_end
            weight_sum = weight_sum + s%dm(shell_idx)
            do species_idx = 1, num_species
               species_sum(species_idx) = species_sum(species_idx) + &
                    s%xa(species_idx,shell_idx)*s%dm(shell_idx)
            end do
         end do
         do species_idx = 1, num_species
            species_sum(species_idx) = species_sum(species_idx)/weight_sum
         end do
         do shell_idx = mix_start, mix_end
            do species_idx = 1, num_species
               s%xa(species_idx,shell_idx) = species_sum(species_idx)
            end do
         end do
      end do
end subroutine homogenize_mixed_regions

! Gravitational settling over the step (diffuse_elements; MHP 12/89,
! after Bahcall & Loeb 1990), sub-stepped so that no sub-step exceeds the
! fraction ctrl%settling_timestep_fraction (legacy DT_GS) of the settling
! timescale at the outer boundary (MHP 6/90). The outer boundary is the
! base of the surface convection zone or, where Y there is below
! ctrl%helium_diffusion_min (legacy YMIN), the first zone inward where Y
! exceeds it. Settling switches itself off for good
! (star%job%diffuse_helium_active) once the central X drops below
! ctrl%hydrogen_diffusion_floor.
subroutine settle
      double precision :: total_mass_unlogged, dlnp_dr(json)
      double precision :: mass_fraction_above, settling_timescale, &
           max_settling_dt, settling_dt
      integer :: k, search_idx, num_settling_steps, substep_idx
      if (s%xa(i_h1,1).lt.s%ctrl%hydrogen_diffusion_floor) then
         s%job%diffuse_helium_active = .false.
         return
      end if
      total_mass_unlogged = exp(ln10*s%log_total_mass)
! d ln P/dr = -G m rho/(r^2 P) per shell (hydrostatic equilibrium).
      do k = 1, s%nz
         dlnp_dr(k) = -exp(ln10*(s%logRho(k)+cgl+ &
              s%log_mass(k)-2.0d0*s%logR(k)- &
              s%logP(k)))
      end do
      if (.not.s%convective_flag(s%nz)) then
! Without a surface convection zone, skip settling with a warning in the
! run log (DBG 2/92: a warning, not a stop).
         write(run_log_unit,911)
  911    format(1x,'NO SURFACE CZ - DIFFUSION NOT MEANINGFUL')
         return
      end if
      do search_idx = envelope_cz_edge, 1, -1
         if (s%xa(i_he4,search_idx).gt.s%ctrl%helium_diffusion_min) exit
      end do
! Y < helium_diffusion_min everywhere from the convection zone's base
! inward: nothing to settle.
      if (search_idx < (1)) return
! The mass fraction q above the outer point, and the settling timescale
! of helium there (Michaud et al. 1984, ApJ 282, 206), stored in s:
!   t_settle = 4.348e21 yr q / T^1.5.
! The sub-step count is int(timestep/max_settling_dt), plus one when
! mod(max_settling_dt, timestep) /= 0 (the arguments reversed, as in
! the F77 mix.f), and at least one: when max_settling_dt is an exact
! multiple of timestep it would be zero, dividing by zero and skipping
! the settling (the guard rotation's burn_and_settle, after rotmix.f,
! already has).
      mass_fraction_above = (total_mass_unlogged- &
           s%m(search_idx))/total_mass_unlogged
      settling_timescale = 4.348d21*seconds_per_year*mass_fraction_above/ &
           exp(ln10*1.5d0*s%logT(search_idx))
      max_settling_dt = s%ctrl%settling_timestep_fraction*settling_timescale
      num_settling_steps = int(timestep/max_settling_dt)
      if (mod(max_settling_dt,timestep).ne.0.0d0 .or. num_settling_steps.eq.0) &
           num_settling_steps = num_settling_steps + 1
      settling_dt = timestep/dfloat(num_settling_steps)
      do substep_idx = 1, num_settling_steps
         call diffuse_elements(settling_dt, s%xa, dlnp_dr, s%logR, &
              s%logRho, s%m, s%logT, s%gradT, deep_mix_flag, &
              s%nz, total_mass_unlogged, ierr)
         if (ierr /= 0) return
      end do
end subroutine settle

! Clip every abundance to [0, 1], then close the mass fractions:
! Z <= 1 - X - He3, and Y = 1 - X - Z - He3.
subroutine clip_and_renormalize
      integer :: k, species_idx
      do k = 1, s%nz
         do species_idx = 1, num_species
            s%xa(species_idx,k) = max(s%xa(species_idx,k),0.0d0)
            s%xa(species_idx,k) = min(s%xa(species_idx,k),1.0d0)
         end do
         s%xa(i_metals,k) = min(s%xa(i_metals,k), &
              1.0d0-s%xa(i_h1,k)-s%xa(i_he3,k))
         s%xa(i_he4,k) = 1.0d0 - s%xa(i_h1,k) - &
              s%xa(i_metals,k) - s%xa(i_he3,k)
      end do
end subroutine clip_and_renormalize

end subroutine advance_composition

end module composition_step
