!----------------------------------------------------------------------
! element_diffusion_lib
!----------------------------------------------------------------------
! The public interface of element diffusion: gravitational settling and
! the microscopic diffusion of helium, the metals and the light
! elements over one timestep. diffuse_elements runs the newer routines
! (microdiff) when use_new_diffusion_routines is set, else the original
! ones (gravitational_settling). reset_element_diffusion restores the
! domain's working state at the start of each run.
module element_diffusion_lib
      use original_settling, only: gravitational_settling
      use microdiff_driver, only: microdiff
      use element_diffusion_scratch, only: diffusion_scratch, diff_scr
      implicit none
      private
      public :: diffuse_elements, reset_element_diffusion

      type(diffusion_scratch), save :: pristine_scratch
      logical, save :: have_pristine_scratch = .false.

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

! ---------------------------------------------------------------
! At the start of every run (yrec_reset): the first time, record the
! pristine working state; afterwards, restore it, so an in-process
! rerun settles exactly as a fresh run does.
subroutine reset_element_diffusion()
      if (.not. have_pristine_scratch) then
         pristine_scratch = diff_scr
         have_pristine_scratch = .true.
      else
         diff_scr = pristine_scratch
      end if
end subroutine reset_element_diffusion

end module element_diffusion_lib
