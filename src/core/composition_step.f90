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
! The order of operations (the options that switch each one):
!   1. mixed regions: convection zones with overshoot (find_convection_zones),
!      with the dpenv deep-mixing option (deep_mixing)
!   2. the burning rates of every shell (compute_burn_rates; with
!      store_start_rates, the deuterium rates are also kept as the
!      start-of-step rates)
!   3. the implicit network burn (burn_network)
!   4. semiconvection (semiconvection, and the control)
!   5. homogenize the mixed regions (homogenize, or with semiconvection on)
!   6. gravitational settling (when enabled), sub-stepped
!   7. clip abundances to [0, 1] and close the mass fractions
!   8. dpenv: recompute the zone bounds with the ordinary flags (deep_mixing)
!   9. deuterium burning (extended composition)
! (2026-10: the predictor pass also ran an explicit burn with He3 and the
! CN cycle in equilibrium after the network, overwriting X and, in
! helium burning, Z, C12 and O16 with a cruder estimate of the same
! step; it is gone, and both passes burn the same way.)
! Li/Be are burned later, after the structure converges (burn_li_be).
module composition_step
      implicit none
      private
      public :: advance_composition

! What a pass does beyond the operations every pass runs.
      type, public :: composition_step_options
! keep this pass's deuterium burning rates as the start-of-step rates
         logical :: store_start_rates = .false.
! homogenize the mixed regions (always done when semiconvection is on)
         logical :: homogenize = .false.
! run semiconvection, when use_semiconvection is set
         logical :: semiconvection = .false.
! apply the dpenv deep-mixing option, when dpenv < 1
         logical :: deep_mixing = .false.
      end type composition_step_options

      type(composition_step_options), parameter, public :: &
           predictor_pass = composition_step_options(store_start_rates=.true., &
                homogenize=.true., semiconvection=.false., deep_mixing=.false.), &
           corrector_pass = composition_step_options(store_start_rates=.false., &
                homogenize=.false., semiconvection=.true., deep_mixing=.true.)

contains

!   star reads: ctrl%dpenv, ctrl%helium_diffusion_min,
!     ctrl%hydrogen_diffusion_floor, ctrl%settling_timestep_fraction,
!     dm, gradT, job%use_extended_composition, job%use_semiconvection,
!     log_mass, log_Teff, log_total_mass, logP, logR, logT,
!     luminosity_lsun, nz
!   star writes: convective_flag, iov1, iov2, job%diffuse_helium_active,
!     logRho, m, mixed_zone_bounds, xa
subroutine advance_composition(timestep, opts, core_cz_edge, &
     envelope_cz_edge, mixed_zone_bounds_no_overshoot, ierr)
      use star_info_lib, only: star, i_h1, i_he3, i_he4, i_metals, json, &
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

! timestep: s
      double precision, intent(in) :: timestep
      type(composition_step_options), intent(in) :: opts
      integer, intent(out) :: core_cz_edge, envelope_cz_edge
      integer, intent(inout) :: mixed_zone_bounds_no_overshoot(max_convective_zones,2)
      integer, intent(out) :: ierr

! rate_by_zone(r_*/rv_*, zone): the composition rate vector of every shell,
! from compute_burn_rates
      double precision :: rate_by_zone(num_rate_rows, json)
      integer :: radiative_zone_bounds(max_radiative_zones,2)
      logical :: deep_mix_flag(json)
      integer :: num_species
      integer :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      integer :: zone_idx, inner_zone_idx

      ierr = 0

      if (star%job%use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if

! 1. MIXED REGIONS: CONVECTION ZONES WITH AND WITHOUT OVERSHOOT
      call set_deep_mix_flags
      call find_convection_zones(star%xa, star%logRho, star%logP, star%logR, &
           star%log_mass, star%logT, deep_mix_flag, star%nz, &
           radiative_zone_bounds, star%mixed_zone_bounds, &
           mixed_zone_bounds_no_overshoot, core_cz_edge, envelope_cz_edge, &
           num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot, star%iov1, star%iov2, ierr)
      if (ierr /= 0) return

! 2-3. BURNING: THE RATES, THEN THE IMPLICIT NETWORK
      call compute_burn_rates(opts%store_start_rates, rate_by_zone)
      call burn_network(radiative_zone_bounds, num_radiative_zones, &
           num_mixed_zones, rate_by_zone, timestep, ierr)
      if (ierr /= 0) return

! 4. SEMICONVECTION
      if (star%job%use_semiconvection) then
         if (opts%semiconvection) &
              call extend_by_semiconvection(timestep, star%xa, star%logRho, &
              star%luminosity_lsun, star%logP, star%logR, star%log_mass, &
              star%logT, star%nz, star%mixed_zone_bounds, &
              num_mixed_zones, star%log_Teff, ierr)
              if (ierr /= 0) return
      end if

! 5. HOMOGENIZE THE MIXED REGIONS. THIS NEEDS TO BE DONE IF
! SEMI-CONVECTION IS BEING CHECKED, OR IF EXPLICIT H AND HE BURNING IS
! BEING DONE.
      if (star%job%use_semiconvection .or. opts%homogenize) &
           call homogenize_mixed_regions
      if (num_mixed_zones.ge.1 .and. solver_diagnostics()) then
         write(run_log_unit,110) ((star%mixed_zone_bounds(zone_idx,inner_zone_idx), &
              inner_zone_idx=1,2),zone_idx=1,num_mixed_zones)
  110    format(' ZONES MIXED IN ORDER--',12('(',i5,',',i5,') ') )
      end if

! 6. GRAVITATIONAL SETTLING
      if (star%job%diffuse_helium_active) then
         call settle
         if (ierr /= 0) return
      end if

! 7. CLIP AND RENORMALIZE
      call clip_and_renormalize

! 8. MHP 1/95: RESET JENV, JCORE FOR DEEP MIXING
      if (star%ctrl%dpenv.lt.1.0d0 .and. opts%deep_mixing) then
         call find_convection_zones(star%xa, star%logRho, star%logP, star%logR, &
              star%log_mass, star%logT, star%convective_flag, star%nz, &
              radiative_zone_bounds, star%mixed_zone_bounds, &
              mixed_zone_bounds_no_overshoot, core_cz_edge, &
              envelope_cz_edge, num_radiative_zones, num_mixed_zones, &
              num_mixed_zones_no_overshoot, star%iov1, star%iov2, ierr)
         if (ierr /= 0) return
      end if

! 9. DEUTERIUM BURNING OVER THE STEP
      if (star%job%use_extended_composition) call burn_deuterium( &
           radiative_zone_bounds, num_radiative_zones, num_mixed_zones, &
           timestep)

contains

! deep_mix_flag: the mixing flag per shell. DPENV < 1 (with deep_mixing)
! forces mixing from the centre out to mass fraction DPENV; outside
! that, and otherwise, the ordinary convective flag.
subroutine set_deep_mix_flags
      double precision :: max_mixed_mass
      integer :: k, j
      if (star%ctrl%dpenv.lt.1.0d0 .and. opts%deep_mixing) then
         max_mixed_mass = star%ctrl%dpenv*exp(ln10*star%log_total_mass)
         do k = 1, star%nz
            if (star%m(k).gt.max_mixed_mass) then
               do j = k, star%nz
                  deep_mix_flag(j) = star%convective_flag(j)
               end do
               exit
            else
               deep_mix_flag(k) = .true.
            end if
         end do
      else
         do j = 1, star%nz
            deep_mix_flag(j) = star%convective_flag(j)
         end do
      end if
end subroutine set_deep_mix_flags

! The mass-weighted mean of every species applied throughout each
! mixed region (star%mixed_zone_bounds: first and last shell).
subroutine homogenize_mixed_regions
      double precision :: species_sum(15), weight_sum
      integer :: mixed_zone_idx, mix_start, mix_end, shell_idx, species_idx
      do mixed_zone_idx = 1, num_mixed_zones
         mix_start = star%mixed_zone_bounds(mixed_zone_idx,1)
         mix_end = star%mixed_zone_bounds(mixed_zone_idx,2)
         if (mix_start.ne.1 .and. mix_start.ge.mix_end) cycle
         weight_sum = 0.0d0
         do species_idx = 1, num_species
            species_sum(species_idx) = 0.0d0
         end do
         do shell_idx = mix_start, mix_end
            weight_sum = weight_sum + star%dm(shell_idx)
            do species_idx = 1, num_species
               species_sum(species_idx) = species_sum(species_idx) + &
                    star%xa(species_idx,shell_idx)*star%dm(shell_idx)
            end do
         end do
         do species_idx = 1, num_species
            species_sum(species_idx) = species_sum(species_idx)/weight_sum
         end do
         do shell_idx = mix_start, mix_end
            do species_idx = 1, num_species
               star%xa(species_idx,shell_idx) = species_sum(species_idx)
            end do
         end do
      end do
end subroutine homogenize_mixed_regions

! Gravitational settling over the step (MHP 12/89, Bahcall & Loeb),
! sub-stepped so that no sub-step exceeds a fraction (DT_GS) of the
! settling timescale at the outer boundary (MHP 6/90). The outer
! boundary is the base of the surface convection zone, or the first
! point inward where Y exceeds a minimum (YMIN). Settling switches
! itself off for good (star%job%diffuse_helium_active) once the
! central X drops below hydrogen_diffusion_floor.
subroutine settle
      double precision :: total_mass_unlogged, dlnp_dr(json)
      double precision :: mass_fraction_above, settling_timescale, &
           max_settling_dt, settling_dt
      integer :: k, search_idx, num_settling_steps, substep_idx
      if (star%xa(i_h1,1).lt.star%ctrl%hydrogen_diffusion_floor) then
         star%job%diffuse_helium_active = .false.
         return
      end if
      total_mass_unlogged = exp(ln10*star%log_total_mass)
! d ln P / dr per shell
      do k = 1, star%nz
         dlnp_dr(k) = -exp(ln10*(star%logRho(k)+cgl+ &
              star%log_mass(k)-2.0d0*star%logR(k)- &
              star%logP(k)))
      end do
      if (.not.star%convective_flag(star%nz)) then
! DBG 2/92 CHANGED STOP TO JUST A WARNING MESSAGE, EXECUTION CONTINUES
         write(run_log_unit,911)
  911    format(1x,'NO SURFACE CZ - DIFFUSION NOT MEANINGFUL')
         return
      end if
      do search_idx = envelope_cz_edge, 1, -1
         if (star%xa(i_he4,search_idx).gt.star%ctrl%helium_diffusion_min) exit
      end do
! Y < YMIN FOR THE WHOLE STAR
      if (search_idx < (1)) return
! the mass fraction above the outer point, and the settling timescale
! of helium there (Michaud et al. 1984, ApJ 282, 206)
      mass_fraction_above = (total_mass_unlogged- &
           star%m(search_idx))/total_mass_unlogged
      settling_timescale = 4.348d21*seconds_per_year*mass_fraction_above/ &
           exp(ln10*1.5d0*star%logT(search_idx))
      max_settling_dt = star%ctrl%settling_timestep_fraction*settling_timescale
      num_settling_steps = int(timestep/max_settling_dt)
      if (mod(max_settling_dt,timestep).ne.0.0d0) num_settling_steps = &
           num_settling_steps + 1
      settling_dt = timestep/dfloat(num_settling_steps)
      do substep_idx = 1, num_settling_steps
         call diffuse_elements(settling_dt, star%xa, dlnp_dr, star%logR, &
              star%logRho, star%m, star%logT, star%gradT, deep_mix_flag, &
              star%nz, total_mass_unlogged, ierr)
         if (ierr /= 0) return
      end do
end subroutine settle

! Guard against small negative or >1 abundances, then close the mass
! fractions: Z <= 1 - X - He3, and Y4 = 1 - X - Z - He3.
subroutine clip_and_renormalize
      integer :: k, species_idx
      do k = 1, star%nz
         do species_idx = 1, num_species
            star%xa(species_idx,k) = max(star%xa(species_idx,k),0.0d0)
            star%xa(species_idx,k) = min(star%xa(species_idx,k),1.0d0)
         end do
         star%xa(i_metals,k) = min(star%xa(i_metals,k), &
              1.0d0-star%xa(i_h1,k)-star%xa(i_he3,k))
         star%xa(i_he4,k) = 1.0d0 - star%xa(i_h1,k) - &
              star%xa(i_metals,k) - star%xa(i_he3,k)
      end do
end subroutine clip_and_renormalize

end subroutine advance_composition

end module composition_step
