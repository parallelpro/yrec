!----------------------------------------------------------------------
! element_diffusion_lib
!----------------------------------------------------------------------
! The public interface of element diffusion: gravitational settling and
! the microscopic diffusion of helium, the metals and the light
! elements over one timestep. diffuse_elements runs the newer routines
! (microdiff) when use_new_diffusion_routines is set, else the original
! ones (gravitational_settling).
module element_diffusion_lib
      use original_settling, only: gravitational_settling
      use microdiff_driver, only: microdiff
      implicit none
      private
      public :: diffuse_elements

contains

! ---------------------------------------------------------------
! One settling step of length timestep (s) over num_zones zones (center
! first): the mass fractions in composition change in place. The other
! arrays are the structure (log radius, density and temperature, the
! enclosed mass in grams, d ln P/dr, the temperature gradient del_grad
! and the convective flags); the routines convert some of them to their
! own units in place and back.
subroutine diffuse_elements(timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, &
     convective_flag, num_zones, total_mass, ierr)
      use star_info_lib, only: star, json
      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json), del_grad(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      integer, intent(out) :: ierr

      if (star%ctrl%use_new_diffusion_routines) then
         call microdiff(timestep, composition, dlnp_dr, log_radius, &
              log_density, enclosed_mass, log_temperature, del_grad, &
              convective_flag, num_zones, total_mass, ierr)
      else
         call gravitational_settling(timestep, composition, dlnp_dr, log_radius, &
              log_density, enclosed_mass, log_temperature, del_grad, &
              convective_flag, num_zones, total_mass, ierr)
      end if
end subroutine diffuse_elements

end module element_diffusion_lib
