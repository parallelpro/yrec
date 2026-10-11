!----------------------------------------------------------------------
! convection_lib
!----------------------------------------------------------------------
! The convection domain's one public module: the temperature gradient
! (mixing-length theory), the convection zones with overshoot and
! semiconvection, and their homogenization. It reads nothing from the
! star; its settings arrive through set_convection_settings (every run
! card) and set_convection_envelope_metals (the starting model).
module convection_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use convection_configuration, only: convection_settings_t, convection_settings
      use temperature_gradients, only: compute_temperature_gradients
      use convection_zones, only: find_convection_zones, homogenize_convection_zones
      use semiconvection, only: extend_by_semiconvection
      implicit none
      private
      public :: convection_settings_t, set_convection_settings, set_convection_envelope_metals
      public :: compute_temperature_gradients, find_convection_zones, &
           homogenize_convection_zones, extend_by_semiconvection

contains

! Install the settings (the star layer's convection_settings_from_star).
subroutine set_convection_settings(settings)
      type(convection_settings_t), intent(in) :: settings
      convection_settings = settings
end subroutine set_convection_settings

! The envelope's metal fraction outside C, N and O, once the starting
! model's mixture is known.
subroutine set_convection_envelope_metals(envelope_metals_no_cno)
      real(dp), intent(in) :: envelope_metals_no_cno
      convection_settings%envelope_metals_no_cno = envelope_metals_no_cno
end subroutine set_convection_envelope_metals

end module convection_lib
