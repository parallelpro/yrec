!----------------------------------------------------------------------
! control_slots
!----------------------------------------------------------------------
! The slot names of the controls that are vectors, and the run cards'
! capacity:
!  - timestep_factors(itime_*) (legacy atime), the timestep-control
!    factors, and chi_grid_scale(ichi_*), the zoning tolerances
!    (star%physics%timestep_factors, star%physics%chi_grid_scale);
!  - the rows of a rescale card's targets, i_target_*
!    (star%control%rescale_params and its run copy, the run plan's
!    card_rescale_targets);
!  - max_runs, the length of the run-card arrays.
! Constants only, so that any module may use them; each user imports
! what it takes from here (the controls buffer, io/controls_buffer.f90,
! imports only max_runs, and nothing re-exports them). (They were the
! buffer's own until 2026-10-29, max_runs until 2026-10-30.)
module control_slots
      implicit none
      private

! Maximum number of runs (kind cards) in one job: the dimension of
! every per-run array of the controls and of the star's run plan (the
! calibration protocols in core/run/calibration.f90 fill them up to
! this limit). The one source of the value: the controls buffer, the
! star's run plan and the calibrations use this parameter; the
! generated declarations of the controls (tools/gen_controls_state.py,
! tools/gen_inlist.py) spell it by name, and tools/controls_registry.py
! reads its value from here (MAX_RUNS, and the registry's per-run type
! cells, real(max_runs)). tools/check_boundaries.py (check_max_runs)
! fails on a registry cell or a generated per-run size that is not it.
      integer, parameter, public :: max_runs = 50

! Slot names for atime (the per-phase timestep-control factors). The
! meaning of each slot is fixed by io/user_inputs.f90, which
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
! tol_*/flag_* controls (io/user_inputs.f90); readers subscript with
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

! Slot names for the rows of a rescale card's targets,
! control%rescale_params(:,card) and its run copy, the run plan's
! card_rescale_targets(:,card), which io/run_cards.f90 fills from the
! like-named rescale_* controls; a target of 0 keeps the model's value
! (core/model/model_rescaling.f90). Readers subscript with these names.
      integer, parameter, public :: &
           i_target_mass = 1, &       ! rescale_mass: the total mass [Msun]
           i_target_x = 2, &          ! rescale_x: the surface X
           i_target_z = 3, &          ! rescale_z: the surface Z
           i_target_core_mass = 4, &  ! rescale_core_mass: the core mass [Msun]
           n_rescale_targets = 4      ! the number of rows

end module control_slots
