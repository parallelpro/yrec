!----------------------------------------------------------------------
! element_diffusion_configuration
!----------------------------------------------------------------------
! Element diffusion's settings: every control it reads, installed by the
! star layer through set_diffusion_settings at the start of every run
! card (setup/diffusion_settings_builder builds the record from the star's
! controls and settling_z_factor_in_use). Element diffusion reads no
! control directly and never the star (tools/check_boundaries.py,
! STAR_FREE_DIRS). README.md section 3 says which pipeline each setting
! acts in.
module element_diffusion_configuration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: diffusion_settings_t, diffusion_settings, set_diffusion_settings

      type :: diffusion_settings_t
! the newer routines (microdiff), else the original ones
! (gravitational_settling)
         logical :: use_new_diffusion_routines = .false.
! settle the metals; settle the light elements (microdiff only)
         logical :: use_diffusion_z = .false., diffuse_lithium = .false.
! the factors on the coefficients of hydrogen and helium (fgry) and of
! the metals (fgrz; in microdiff the light elements' too); a Monte Carlo
! run varies fgrz
         real(dp) :: fgry = 1.0d0, fgrz = 1.0d0
! no settling once X falls below hydrogen_diffusion_floor, or Y below
! helium_diffusion_min, in the whole model
         real(dp) :: hydrogen_diffusion_floor = 0.0d0, helium_diffusion_min = 0.0d0
! the passes of the implicit solve and their tolerance: the original
! routines iterate on the coefficients; microdiff solves once, and these
! decide only its messages (docs/bug_sweep.md M253)
         integer :: settling_num_iterations = 0
         real(dp) :: settling_tolerance = 0.0d0
! the original routines only: the coefficients (Thoul et al. 1994 with
! this Coulomb logarithm, or their fit, or the older analytic ones), and
! with Thoul's a factor on every coefficient (constant_mixing_coeff) and
! one more on the metals' drift (constant_settling_reduction)
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
