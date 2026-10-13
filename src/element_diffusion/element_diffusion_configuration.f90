!----------------------------------------------------------------------
! element_diffusion_configuration
!----------------------------------------------------------------------
! Element diffusion's settings: every control it reads, installed by the
! star layer through set_diffusion_settings at the start of every run
! card (setup/microphysics_settings builds the record from star%ctrl and
! star%job). Element diffusion reads no control directly and never the
! star (tools/check_boundaries.py, STAR_FREE_DIRS).
module element_diffusion_configuration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: diffusion_settings_t, diffusion_settings, set_diffusion_settings

      type :: diffusion_settings_t
! the newer routines (microdiff), else the original ones
! (gravitational_settling)
         logical :: use_new_diffusion_routines = .false.
! settle the metals; diffuse the light elements (newer routines)
         logical :: use_diffusion_z = .false., diffuse_lithium = .false.
! the scale factors of the helium and the metal diffusion coefficients
! (fgrz: a Monte-Carlo run rescales it)
         real(dp) :: fgry = 1.0d0, fgrz = 1.0d0
! no settling once X falls below hydrogen_diffusion_floor, or Y below
! helium_diffusion_min, in the whole model
         real(dp) :: hydrogen_diffusion_floor = 0.0d0, helium_diffusion_min = 0.0d0
! the iterations of the implicit solve and their tolerance
         integer :: settling_num_iterations = 0
         real(dp) :: settling_tolerance = 0.0d0
! the original routines: the Thoul et al. (1994) coefficients (or their
! fit) with this Coulomb logarithm, or a constant mixing coefficient and
! its settling reduction
         integer :: coulomb_log_choice = 0
         logical :: use_thoul_fit = .false., use_thoul_diffusion = .false.
         real(dp) :: constant_mixing_coeff = 1.0d0, constant_settling_reduction = 1.0d0
      end type diffusion_settings_t

      type(diffusion_settings_t), save :: diffusion_settings

contains

! Installs the run card's diffusion settings.
subroutine set_diffusion_settings(settings)
      type(diffusion_settings_t), intent(in) :: settings
      diffusion_settings = settings
end subroutine set_diffusion_settings

end module element_diffusion_configuration
