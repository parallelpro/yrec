!----------------------------------------------------------------------
! element_diffusion_lib
!----------------------------------------------------------------------
! The public interface of element diffusion: gravitational settling and
! the microscopic diffusion of helium, the metals and the light
! elements over one timestep. diffuse_elements runs the newer routines
! (microdiff) when use_new_diffusion_routines is set, else the original
! ones (gravitational_settling). Its settings come from
! set_diffusion_settings (element_diffusion_configuration), which the
! star layer calls at the start of every run card.
module element_diffusion_lib
      use original_settling, only: gravitational_settling
      use microdiff_driver, only: microdiff
      use element_diffusion_configuration, only: diffusion_settings_t, &
           diffusion_settings, set_diffusion_settings
      use element_diffusion_scratch, only: diffusion_scratch
      implicit none
      private
      public :: diffuse_elements, diffusion_settings_t, set_diffusion_settings

contains

! ---------------------------------------------------------------
! One settling step of length timestep (s) over num_zones zones (center
! first): the mass fractions in composition change in place. The other
! arrays are the structure (log radius, density and temperature, the
! enclosed mass in grams, d ln P/dr, the temperature gradient del_grad
! and the convective flags); the routines convert some of them to their
! own units in place and back. diffuse_helium: settle hydrogen and
! helium (the newer routines; the metals and light elements settle
! regardless). suspension_reported: the run log's latch (the star's), set
! when a model too convective or exhausted to settle is reported, cleared
! when settling resumes.
subroutine diffuse_elements(timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, &
     convective_flag, num_zones, total_mass, diffuse_helium, &
     suspension_reported, ierr)
      use composition_layout, only: json
      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json), del_grad(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      logical, intent(in) :: diffuse_helium
      logical, intent(inout) :: suspension_reported
      integer, intent(out) :: ierr
! the settling routines' working arrays for this step
! (element_diffusion_scratch), passed down as an argument
      type(diffusion_scratch), allocatable :: diff_scr

      allocate(diff_scr)
      if (diffusion_settings%use_new_diffusion_routines) then
         call microdiff(diff_scr, timestep, composition, dlnp_dr, log_radius, &
              log_density, enclosed_mass, log_temperature, del_grad, &
              convective_flag, num_zones, total_mass, diffuse_helium, &
              suspension_reported, ierr)
      else
         call gravitational_settling(diff_scr, timestep, composition, dlnp_dr, log_radius, &
              log_density, enclosed_mass, log_temperature, del_grad, &
              convective_flag, num_zones, total_mass, suspension_reported, ierr)
      end if
end subroutine diffuse_elements

end module element_diffusion_lib
