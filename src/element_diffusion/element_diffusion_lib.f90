!----------------------------------------------------------------------
! element_diffusion_lib
!----------------------------------------------------------------------
! The public interface of element diffusion: gravitational settling and
! the microscopic diffusion of helium, the metals and the light
! elements over one timestep (README.md).
!   diffuse_elements        one settling step: the settling region
!                           (settling_region), then the newer routines
!                           (microdiff) when use_new_diffusion_routines
!                           is set, else the original ones
!                           (gravitational_settling);
!   plan_settling_substeps  the sub-steps of a step's settling
!                           (settling_substeps), which its callers take;
!   set_diffusion_settings  installs the settings record
!                           (element_diffusion_configuration), as the
!                           star layer does at the start of every run
!                           card.
module element_diffusion_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use original_settling, only: gravitational_settling
      use microdiff_driver, only: microdiff
      use element_diffusion_configuration, only: diffusion_settings_t, &
           diffusion_settings, set_diffusion_settings
      use element_diffusion_scratch, only: diffusion_scratch
      use settling_region, only: find_settling_region, region_found
      use settling_substeps, only: settling_substep_plan, plan_settling_substeps
      use element_diffusion_messages, only: note_settling_region
      implicit none
      private
      public :: diffuse_elements, diffusion_settings_t, set_diffusion_settings, &
           settling_substep_plan, plan_settling_substeps

contains

! ---------------------------------------------------------------
! One settling step of timestep seconds over a model of num_zones zones
! (centre first): the mass fractions composition(species, zone) change
! in place. The structure: log10 radius (cm), density (g/cm^3) and
! temperature (K), the enclosed mass (g), d ln P/dr (1/cm), the
! temperature gradient del_grad = d ln T/d ln P, the convective flags and
! the total mass (g). microdiff leaves them as they are; the original
! routines convert timestep, dlnp_dr, enclosed_mass and total_mass to
! Bahcall & Loeb units in place and back, which can change their last
! bit. Hydrogen and helium settle always (the callers call only while
! the star's settling_active holds), the metals with use_diffusion_z and
! the light elements with diffuse_lithium (microdiff only).
! suspension_reported: the run log's latch (the star's), set when a
! model that cannot settle is noted, cleared when settling resumes
! (element_diffusion_messages' note_settling_region). ierr: a numerics
! routine's status (a singular matrix), after its message; nonzero stops
! the run.
subroutine diffuse_elements(timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, &
     convective_flag, num_zones, total_mass, &
     suspension_reported, ierr)
      real(dp), intent(inout) :: timestep
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(inout) :: dlnp_dr(:)
      real(dp), intent(in) :: log_radius(:), log_density(:)
      real(dp), intent(inout) :: enclosed_mass(:)
      real(dp), intent(in) :: log_temperature(:), del_grad(:)
      logical, intent(in) :: convective_flag(:)
      integer, intent(in) :: num_zones
      real(dp), intent(inout) :: total_mass
      logical, intent(inout) :: suspension_reported
      integer, intent(out) :: ierr
! the original routines' working arrays for this step
! (element_diffusion_scratch), passed down as an argument
      type(diffusion_scratch), allocatable :: diff_scr
! the zones that settle
      integer :: zone_begin, zone_end, region

      ierr = 0
      call find_settling_region(convective_flag, composition, num_zones, &
           diffusion_settings%hydrogen_diffusion_floor, &
           diffusion_settings%helium_diffusion_min, zone_begin, zone_end, region)
      call note_settling_region(region, diffusion_settings%hydrogen_diffusion_floor, &
           diffusion_settings%helium_diffusion_min, suspension_reported)
      if (region /= region_found) return
      if (diffusion_settings%use_new_diffusion_routines) then
         call microdiff(timestep, composition, dlnp_dr, log_radius, &
              log_density, enclosed_mass, log_temperature, del_grad, &
              num_zones, total_mass, zone_begin, zone_end, ierr)
      else
         allocate(diff_scr)
         call gravitational_settling(diff_scr, timestep, composition, dlnp_dr, log_radius, &
              log_density, enclosed_mass, log_temperature, del_grad, &
              num_zones, total_mass, zone_begin, zone_end, ierr)
      end if
end subroutine diffuse_elements

end module element_diffusion_lib
