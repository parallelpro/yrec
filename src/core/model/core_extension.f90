!----------------------------------------------------------------------
! core_extension
!----------------------------------------------------------------------
! The points a model gains inward of its innermost one when the deck
! asks for a deeper central fitting point: a star-blind kernel
! (docs/GUIDELINES.md, star rule 6). extend_core_toward_center
! (starting_model.f90) reads the model from the star, calls these and
! stores the points; nothing here reads the star, prints or keeps state.
!
! The new points follow what the model assumes inside its central
! fitting point: constant density, constant energy generation rate and
! constant d log T/d log P (MHP 9/98). test_model_kernels checks them.
module core_extension
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      use physical_constants, only: c4pi3l, c4pil, cc13, ln10
      implicit none
      private
      public :: core_extension_count, extend_core_points

contains

! ---------------------------------------------------------------
! The number of points to add so that the innermost mass drops by
! mass_reduction_factor, in whole steps of log_mass_spacing in log m, so
! that the mesh needs no reshuffling: one more than the whole steps in
! log10 of the factor.
pure integer function core_extension_count(mass_reduction_factor, log_mass_spacing)
      real(dp), intent(in) :: mass_reduction_factor, log_mass_spacing

      core_extension_count = int(log10(mass_reduction_factor)/log_mass_spacing) + 1
end function core_extension_count

! ---------------------------------------------------------------
! The structure of the points added inward of a model's innermost point:
! log_mass(1:N), log_radius, luminosity, log_pressure, log_temperature
! and log_density, point N the one next to the old innermost point,
! whose log m, log P, log T and L are the inner_* arguments. Marching
! inward from it in steps of log_mass_spacing:
!   log m_i = log m_(i+1) - spacing;
!   constant density rho_c (central_log_density), m = (4 pi/3) rho_c r^3:
!     log r = (1/3)(log m - log(4 pi/3) - log rho_c);
!   constant epsilon: L = L_inner m/m_inner;
!   hydrostatic equilibrium, d log P/d log m = -G m^2/(4 pi r^4 P), taken
!   in two stages (the slope at log P half a step in from point i+1, then
!   the full step with it):
!     log P_i = log P_(i+1) + (G m_i^2/(4 pi r_i^4 P_half)) spacing;
!   constant del (temperature_gradient):
!     log T_i = log T_(i+1) + del (log P_i - log P_(i+1));
!   the density estimate log rho = log P - log T - density_offset, at
!   fixed mean molecular weight, in place of an EOS call (MHP 4/12).
! log10_newton_g is log10 G; the luminosity is in the unit of
! inner_luminosity.
pure subroutine extend_core_points(inner_log_mass, inner_log_pressure, &
     inner_log_temperature, inner_luminosity, central_log_density, temperature_gradient, &
     density_offset, log_mass_spacing, log10_newton_g, log_mass, log_radius, luminosity, &
     log_pressure, log_temperature, log_density)
      real(dp), intent(in) :: inner_log_mass, inner_log_pressure, inner_log_temperature, &
           inner_luminosity, central_log_density, temperature_gradient, density_offset, &
           log_mass_spacing, log10_newton_g
      real(dp), intent(out) :: log_mass(:), log_radius(:), luminosity(:), log_pressure(:), &
           log_temperature(:), log_density(:)
      real(dp) :: outer_log_mass, outer_log_pressure, outer_log_temperature, &
           trial_log_pressure, dlogp_dlogm
      integer :: i

      outer_log_mass = inner_log_mass
      outer_log_pressure = inner_log_pressure
      outer_log_temperature = inner_log_temperature
      do i = size(log_mass), 1, -1
         log_mass(i) = outer_log_mass-log_mass_spacing
         log_radius(i) = cc13*(log_mass(i)-c4pi3l-central_log_density)
         luminosity(i) = exp(ln10*(log_mass(i)-inner_log_mass))*inner_luminosity
         trial_log_pressure = outer_log_pressure
         dlogp_dlogm = exp(ln10*(log10_newton_g+2.0d0*log_mass(i)-c4pil- &
              trial_log_pressure-4.0d0*log_radius(i)))
         trial_log_pressure = trial_log_pressure+0.5d0*dlogp_dlogm* &
              log_mass_spacing
         dlogp_dlogm = exp(ln10*(log10_newton_g+2.0d0*log_mass(i)-c4pil- &
              trial_log_pressure-4.0d0*log_radius(i)))
         log_pressure(i) = outer_log_pressure+dlogp_dlogm* &
              log_mass_spacing
         log_temperature(i) = outer_log_temperature+dlogp_dlogm* &
              temperature_gradient*log_mass_spacing
         log_density(i) = log_pressure(i) - log_temperature(i) - &
              density_offset
         outer_log_mass = log_mass(i)
         outer_log_pressure = log_pressure(i)
         outer_log_temperature = log_temperature(i)
      end do
end subroutine extend_core_points

end module core_extension
