!----------------------------------------------------------------------
! element_diffusion_scratch
!----------------------------------------------------------------------
! The original routines' working state (diff_scr): the metals'
! coefficients and abundance, which gravitational_settling_setup,
! model_to_equal and equal_to_model pass between them. diffuse_elements
! owns the step's instance and passes it down as an argument; this
! module holds only the type. Nothing in it outlives a step, but a
! field the step does not write keeps whatever the memory held: with
! use_diffusion_z on and use_thoul_diffusion off the metals'
! coefficients are never set (docs/bug_sweep.md M33). (The run log's
! suspension latch is the star's, settling_suspended_reported.)
module element_diffusion_scratch
      use iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      implicit none
      private

      type, public :: diffusion_scratch
! the metals' drift and diffusion coefficients (D1, D2) and their
! derivatives in Z, on the model's zones
           real(dp) :: zone_metal_coeff1(json), zone_metal_coeff2(json), &
                zone_metal_coeff1_dz(json), zone_metal_coeff2_dz(json)
! on the equal grids: D1 at the centres; D1, D2 and their derivatives in
! Z at the midpoints
           real(dp) :: metal_coeff1(json), metal_coeff1_mid(json), &
                metal_coeff2_mid(json), metal_coeff1_dz_mid(json), &
                metal_coeff2_dz_mid(json)
! Z on the equal grids: at the centres, then its change over the step;
! at the midpoints, then the first Lax-Wendroff step's change
           real(dp) :: eq_metal_fraction(json), eq_metal_fraction_mid(json)
      end type diffusion_scratch

end module element_diffusion_scratch
