!----------------------------------------------------------------------
! envelope_refit
!----------------------------------------------------------------------
! Appending envelope points to the model as shells, the step the two
! envelope refits share. Before the call, atm_get has integrated a
! fresh envelope from the surface down to the current fitting point,
! and the caller has shifted env_points so that its point 1 is the
! first point outside the last interior shell.
!
! Callers: core/envelope/envelope_rebuild.f90 (the fitting point moved
! outward during the run) and core/model/starting_model.f90 (a
! requested envelope thinner than the starting model's). What they do
! before the call (the surface L, Teff and g, the env_points shift)
! and after it (shell masses, rotation, messages) stays with them. They
! differ in one test, carried by accept_point_at_fit_mass: whether an
! envelope point exactly at the fitting mass is taken as a shell
! (read_starting_model) or interpolated (rebuild_envelope). Each
! reports its own message when ierr comes back nonzero.
!
! The arithmetic is envelope_refit_kernel's append_envelope_shells,
! star-free; this driver hands it the star's arrays.
module envelope_refit
      use star_info_lib, only: star_info
      implicit none
      private
      public :: append_envelope_points

contains

! ---------------------------------------------------------------
! Appends the star's env_points inside the new fitting mass
! star%log_fit_point_mass_fraction to the model as zones nz+1, nz+2,
! ..., the first point outside it interpolated onto the fitting mass,
! and advances star%nz (append_envelope_shells; species_end_index is its
! last_copied_species). ierr = 1 when the mass interval of the
! interpolation is degenerate.
!   star reads: env_points, log_fit_point_mass_fraction, log_total_mass
!   star writes: convective_flag, log_mass, logP, logR, logRho, logT,
!     luminosity_lsun, nz, xa
subroutine append_envelope_points(s, accept_point_at_fit_mass, species_end_index, &
     ierr)
      use envelope_refit_kernel, only: append_envelope_shells
      type(star_info), intent(inout) :: s
      logical, intent(in) :: accept_point_at_fit_mass
      integer, intent(in) :: species_end_index
      integer, intent(out) :: ierr

      call append_envelope_shells(s%env_points, s%log_fit_point_mass_fraction, &
           s%log_total_mass, accept_point_at_fit_mass, species_end_index, s%nz, &
           s%log_mass, s%logP, s%logT, s%logR, s%logRho, s%luminosity_lsun, &
           s%convective_flag, s%xa, ierr)
end subroutine append_envelope_points

end module envelope_refit
