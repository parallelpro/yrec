!----------------------------------------------------------------------
! control_slots
!----------------------------------------------------------------------
! The slot names of the controls that are vectors: atime(itime_*), the
! timestep-control factors, and chi_grid_scale(ichi_*), the zoning
! tolerances (star%ctrl%atime, star%ctrl%chi_grid_scale). Constants
! only, so that any module may use them: the controls buffer
! (io/controls_buffer.f90) re-exports them for the controls read, and
! the star layer and setup/'s builders use this module. (They were the
! buffer's own until 2026-10-29.)
module control_slots
      implicit none
      private

! Slot names for atime (the per-phase timestep-control factors). The
! meaning of each slot is fixed by setup/user_inputs.f90, which
! fills it from the like-named time_* namelist inputs when lnewvars
! is set; readers subscript with these names.
      integer, parameter, public :: &
           itime_core_min = 1, &      ! time_core_min: min central abundance for the core-burning criterion
           itime_dx_core_tot = 2, &   ! time_dx_core_tot: max absolute delta-X burned in core per step
           itime_dx_core_frac = 3, &  ! time_dx_core_frac: max fraction of core X burned per step
           itime_dy_core_tot = 4, &   ! time_dy_core_tot: max absolute delta-Y burned in core per step
           itime_dy_core_frac = 5, &  ! time_dy_core_frac: max fraction of core Y burned per step
           itime_dx_total = 6, &      ! time_dx_total: max Msun of X burned per step
           itime_dx_shell = 7, &      ! time_dx_shell: max fraction of X burned in the shell per step
           itime_dt = 8, &            ! time_dt: max change in T at any shell per step
           itime_dp = 9, &            ! time_dp: max change in P at any shell per step
           itime_dr = 10, &           ! time_dr: max change in R at any shell per step
           itime_dl = 11, &           ! time_dl: max change in L at any shell per step
           itime_dy_shell = 12, &     ! time_dy_shell: max fraction of Y burned in the shell per step
           itime_max_dt_frac = 13, &  ! time_max_dt_frac: max timestep change factor model to model
           itime_dy_total = 14        ! time_dy_total: max Msun of Y burned in the shell per step

! Slot names for chi_grid_scale, the zoning tolerances: the deck's
! hpttol (zoning_tolerances), or, when lnewvars is set, the like-named
! tol_*/flag_* controls (setup/user_inputs.f90); readers subscript with
! these names.
      integer, parameter, public :: &
           ichi_dm_min = 1, &            ! tol_dm_min: min log-m spacing
           ichi_dm_max = 2, &            ! tol_dm_max: max log-m spacing
           ichi_flag_dx = 3, &           ! flag_dx: X jump that flags a shell
           ichi_flag_dz = 4, &           ! flag_dz: Z jump that flags a shell
           ichi_dx_max = 5, &            ! tol_dx_max: max X spacing
           ichi_dz_max = 6, &            ! tol_dz_max: max Z spacing
           ichi_czbase_fine_width = 7, & ! tol_czbase_fine_width: log-P width of fine zoning at CZ base
           ichi_dp_env_max = 8, &        ! tol_dp_env_max: max log-P spacing, convective envelope
           ichi_dl_max = 9, &            ! tol_dl_max: max L/Lsurf spacing
           ichi_dp_czbase_max = 10, &    ! tol_dp_czbase_max: max log-P spacing near CZ base
           ichi_dp_core_max = 11, &      ! tol_dp_core_max: max log-P spacing, radiative interior
           ichi_flag_dw = 12             ! flag_dw: max log-omega spacing (rotating models)

end module control_slots
