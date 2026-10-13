!----------------------------------------------------------------------
! element_diffusion_scratch
!----------------------------------------------------------------------
! The settling routines' own working state, one module-level instance
! (diff_scr): the metal diffusion coefficients on the model and the
! equally spaced grids, the change of the metal abundance, and the
! Bahcall & Loeb (1990) unit scales. Nothing in it outlives a step: each
! settling call writes what it reads (diff_scr poisoned with NaN at
! every step start, 2026-10-13: every pinned deck byte-identical), so an
! in-process rerun needs no restore. (The run-log latch that reports a
! suspension of settling once is the star's,
! star%settling_suspended_reported.)
module element_diffusion_scratch
      use star_info_lib, only: json
      implicit none
      private

      type, public :: diffusion_scratch
! the metal diffusion coefficients: on the model grid, at the grid
! mid-points, and on the equally spaced grid's mid-points
           double precision :: metal_diffusion_coeff1(json), &
                metal_diffusion_coeff1_mid(json), &
                metal_diffusion_coeff2_mid(json), &
                eq_metal_diffusion_coeff1_mid(json), &
                eq_metal_diffusion_coeff2_mid(json), &
                metal_abundance_change_mid(json)
! the original routines' coefficients on the source grid, and their
! derivatives with respect to Z
           double precision :: src_grid_metal_diffusion_coeff1(json), &
                src_grid_metal_diffusion_coeff2(json), &
                src_grid_metal_diffusion_coeff1_dz(json), &
                src_grid_metal_diffusion_coeff2_dz(json)
           double precision :: metal_abundance_change(json)
! the Bahcall & Loeb units (set_bahcall_loeb_scales)
           double precision :: bl_radius_scale, bl_mass_scale, &
                bl_temp_scale, bl_time_scale
      end type diffusion_scratch

      type(diffusion_scratch), public, save :: diff_scr

end module element_diffusion_scratch
