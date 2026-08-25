!----------------------------------------------------------------------
! bursmix
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original bursmix.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This subroutine is designed to alternate between nuclear burning
! plus gravitational settling and mixing in a series of progressively
! smaller timesteps, extrapolating to the answer for zero timestep.
! It assumes that the diffusion coefficients, structure, and nuclear
! reaction rates are held constant during the time step. The procedure
! is as follows:
! Store the results for the whole timestep in entry one.
! Then call the nuclear burning routines to burn and the rotational
! mixing routine to mix for a series of smaller time steps, storing the
! answer as a function of time step.
subroutine bursmix(diffusion_coeff, timestep, composition, log_density, &
     log_luminosity, log_pressure, log_radius, log_mass, enclosed_mass, &
     shell_mass, log_total_mass, log_temperature, velocity, zone_max, &
     zone_min, env_cz_zone_old, env_cz_zone, final_iteration_flag, &
     convective_flag, num_zones, radiative_zone_bounds, mixed_zone_bounds, &
     num_radiative_zones, num_zones_mixed, ierr)
      use net_lib
      use star_info_lib, only: star, json
      use burn_lib
      implicit none






! INPUT VARIABLES
      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: log_density(json)
      double precision, intent(inout) :: log_luminosity(json)
      double precision, intent(inout) :: log_pressure(json)
      double precision, intent(inout) :: log_radius(json)
      double precision, intent(inout) :: log_mass(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(inout) :: log_temperature(json)
      double precision, intent(inout) :: velocity(json)
      integer, intent(inout) :: zone_max, zone_min
      integer, intent(inout) :: env_cz_zone_old, env_cz_zone
      logical, intent(in) :: final_iteration_flag
      logical, intent(inout) :: convective_flag(json)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: radiative_zone_bounds(13,2)
      integer, intent(inout) :: mixed_zone_bounds(12,2)
      integer, intent(inout) :: num_radiative_zones, num_zones_mixed
! LOCAL VARIABLES
      double precision :: composition_saved(15,json)
      integer :: burn_rate_sequence(11)
      double precision :: composition_kept(3,json)
      data burn_rate_sequence/1,2,4,6,8,10,12,14,16,18,20/
! equally_spaced_diffusion_coeff/equally_spaced_mass (originally ECOD2/
! EM) are scratch arrays for ndifcom's internal equally-spaced-grid
! solve; zone_begin/zone_end (originally IBEG/IEND) are scratch
! outputs from ndifcom that are not used again here.
      double precision :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      integer :: zone_begin, zone_end
      integer :: num_species, species_begin, species_end, num_substeps
      integer :: zone_idx, species_idx, extrap_idx, substep_idx, num_steps
      double precision :: substep_dt
      logical :: converged

      integer, intent(out) :: ierr

      ierr = 0

      if (star%job%use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if
      species_begin = 4
      species_end = num_species
      num_substeps = 11
! STORE ALL ELEMENTS HERE - HCOMPP IS USED AS THE START OF
! TIMESTEP VARIABLE IN KEMCOM, SO IT NEEDS TO BE TEMPORARILY
! OVERWRITTEN IN THE SOLUTION PROCESS
      do zone_idx = 1, num_zones
         do species_idx = 1, num_species
            composition_saved(species_idx,zone_idx) = &
                 star%xa_start(species_idx,zone_idx)
         end do
         do species_idx = 1, 3
            composition_kept(species_idx,zone_idx) = &
                 composition(species_idx,zone_idx)
         end do
      end do
! STORE RESULTS OF INITIAL CALCULATION AS FIRST ENTRY IN BS
! ROUTINE
      extrap_idx = 1
      call bsrotmix(timestep, composition, extrap_idx, num_zones, &
           species_begin, num_species, burn_rate_sequence, converged)
      do extrap_idx = 2, num_substeps
! RESET COMPOSITION TO INITIAL VALUE.
         do zone_idx = 1, num_zones
            do species_idx = 1, num_species
               composition(species_idx,zone_idx) = &
                    composition_saved(species_idx,zone_idx)
               star%xa_start(species_idx,zone_idx) = &
                    composition_saved(species_idx,zone_idx)
            end do
         end do
         num_steps = burn_rate_sequence(extrap_idx)
         substep_dt = timestep/dfloat(num_steps)
         do substep_idx = 1, num_steps
! PERFORM NUCLEAR BURNING
            call rotmix(substep_dt, composition, shell_mass, &
                 log_temperature, num_zones, radiative_zone_bounds, &
                 mixed_zone_bounds, num_radiative_zones, num_zones_mixed, &
                 log_total_mass, log_density, log_mass, log_radius, &
                 log_pressure, convective_flag, enclosed_mass, ierr)
            if (ierr /= 0) return
! PERFORM MIXING
            call ndifcom(substep_dt, diffusion_coeff, &
                 equally_spaced_diffusion_coeff, equally_spaced_mass, &
                 log_density, log_luminosity, log_pressure, log_radius, &
                 log_mass, enclosed_mass, shell_mass, log_total_mass, &
                 velocity, zone_begin, zone_end, zone_max, zone_min, &
                 convective_flag, final_iteration_flag, num_zones, &
                 composition, species_begin, species_end)
            if (star%job%use_extended_composition) then
               call liburn2(substep_dt, composition, log_radius, &
                    enclosed_mass, shell_mass, log_temperature, &
                    env_cz_zone, env_cz_zone_old, num_zones)
            end if
            do zone_idx = 1, num_zones
               do species_idx = 1, num_species
                  star%xa_start(species_idx,zone_idx) = &
                       composition(species_idx,zone_idx)
               end do
            end do
         end do
         call bsrotmix(timestep, composition, extrap_idx, num_zones, &
              species_begin, num_species, burn_rate_sequence, converged)
         if (converged) then
            do zone_idx = 1, num_zones
               do species_idx = 1, num_species
                  star%xa_start(species_idx,zone_idx) = &
                       composition(species_idx,zone_idx)
               end do
            end do
            exit
         end if
      end do
      if (extrap_idx > num_substeps) then
! FAILED TO CONVERGE; PRINT WARNING
! IN THIS CASE THE UNEXTRAPOLATED FINAL COMPOSITION IS USED.
      write(*,5)
    5 format(' WARNING - EXTRAPOLATION IN BSBURN DID NOT CONVERGE')
      end if
      do zone_idx = 1, num_zones
         composition(1,zone_idx) = composition_kept(1,zone_idx)
         composition(3,zone_idx) = composition_kept(3,zone_idx)
         composition(2,zone_idx) = 1.0d0 - composition(1,zone_idx) - &
              composition(3,zone_idx) - composition(4,zone_idx)
      end do
      if (star%job%use_extended_composition) then
         do zone_idx = 1, num_zones
            star%rate_li6_start(zone_idx) = star%rate_li6(zone_idx)
            star%rate_li7_start(zone_idx) = star%rate_li7(zone_idx)
            star%rate_be9_start(zone_idx) = star%rate_be9(zone_idx)
         end do
         env_cz_zone_old = env_cz_zone
         star%pressure_scale_height_start = star%pressure_scale_height_end
      end if
      return
end subroutine bursmix
