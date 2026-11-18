!----------------------------------------------------------------------
! rezoning
!----------------------------------------------------------------------
! The rezone of each step (rezone): a new mesh in log m, its points
! spaced evenly in a coordinate chi built from the changes of log m, L
! and log P, with points added at X and Z jumps; the model interpolated
! onto it and, for rotating models, the rotation re-spun on it. The
! arithmetic is two star-blind kernels: the new mesh (chi, the points
! along it, the points at jumps, the crowded points) in mesh_spacing.f90,
! and the regridding (the shell masses of a mesh, the interpolation of
! the profiles, the angular momentum totals) in regridding.f90; rezone
! gathers their inputs from the star, reports their failures and stores
! their results.
module rezoning
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use shell_physics, only: compute_shell_physics
      use star_info_lib, only: star_info
      use am_regions, only: am_convective_regions
      use omega_profile, only: omega_from_j
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use star_info_lib, only: i_eps_grav, i_eps_neu, i_h1, i_metals, json, n_species_basic, &
           n_species_extended, max_convective_zones, max_radiative_zones
      use control_slots, only: ichi_czbase_fine_width, ichi_dl_max, ichi_dm_max, &
           ichi_dm_min, ichi_dp_core_max, ichi_dp_czbase_max, ichi_dp_env_max, ichi_dx_max, &
           ichi_dz_max, ichi_flag_dw, ichi_flag_dx, ichi_flag_dz
      use settings_install, only: install_surface_x_settings
      use run_messages, only: solver_diagnostics
      use core_messages, only: core_log, core_report, message_length, unwritten_line
      use physical_constants, only: cc23, ln10
      use numerics_lib, only: numerics_status_text
      use regridding, only: structure_profiles, rotation_profiles, &
           compute_shell_masses, regrid_profile, regrid_composition, &
           regrid_structure, regrid_rotation, compute_angular_momentum_totals
      use mesh_spacing, only: chi_scales, accumulate_chi, sort_kept_points, &
           space_points_along_chi, spline_onto_points, insert_points_at_jumps, &
           remove_crowded_points
      implicit none
      private
      public :: rezone

! The most points flag_fixed_points keeps (append_flag_point guards
! every append; one slot stays for the star%nz terminator).
      integer, parameter :: max_flag_points = 100

! The point added outside a surface hotter than the envelope window:
! its log T rises by added_point_dlogt_dlogp of its log P's rise [-], and
! its log rho is log P - log T - added_point_log_density_offset, a rough
! ideal-gas estimate [dex]. The surface opacity tables are rebuilt when
! the surface X has moved by more than surface_x_refresh_tol since they
! were built [mass fraction].
      real(dp), parameter :: added_point_dlogt_dlogp = 0.250D0, &
           added_point_log_density_offset = 8.0D0, surface_x_refresh_tol = 1.0D-8

! The rezone's state, which its steps pass on to each other:
!   flag_point(1:flag_count)  the points the rezone keeps and does not
!                             smooth across (flag_fixed_points), ending
!                             with star%nz;
!   flag_overflow_warned      append_flag_point has reported a dropped
!                             point in this rezone;
!   max_luminosity_step, max_x_step, max_z_step
!                             the largest steps in L, X and Z between
!                             new points (seed_new_point_run);
!   chi                       the running sum of the normalized changes of
!                             log m, L and log P from shell to shell (MHP
!                             7/02), which sets the spacing;
!   new_log_mass(1:new_num_zones)  the new mesh in log m (legacy HSS),
!                             adopted as star%log_mass by
!                             transfer_new_points;
!   ft_splined, fp_splined    the shape factors splined onto the new mesh,
!                             which transfer_new_points adopts;
!   old_point_count, new_point_count  the old and the new number of points;
!   num_species_tracked       the species carried: 11, or 15 with the
!                             extended composition.
      type :: rezone_state
         integer :: flag_point(max_flag_points) = 0
         integer :: flag_count = 0
         logical :: flag_overflow_warned = .false.
         real(dp) :: max_luminosity_step = 0.0d0
         real(dp) :: max_x_step = 0.0d0
         real(dp) :: max_z_step = 0.0d0
         real(dp) :: chi(json) = 0.0d0
         real(dp) :: new_log_mass(json) = 0.0d0
         integer :: new_num_zones = 0
         real(dp) :: ft_splined(json) = 0.0d0
         real(dp) :: fp_splined(json) = 0.0d0
         integer :: old_point_count = 0
         integer :: new_point_count = 0
         integer :: num_species_tracked = 0
      end type rezone_state

contains

!----------------------------------------------------------------------
! rezone (legacy hpoint.f)
!----------------------------------------------------------------------
! Rezones the model so that no step between successive points exceeds
! the user's limits on log m, L, log P, X and Z (chi_grid_scale, legacy
! HPTTOL), without smoothing across real discontinuities:
!   1. check_envelope_temperature_range: the outermost point must lie in
!      the envelope temperature window (tenv0, tenv1); outer points are
!      deleted, or one is added, otherwise.
!   2. flag_fixed_points: the points the rezone keeps and does not
!      smooth across (convection-zone edges, the H-shell base, X, Z and
!      omega jumps, the overshoot and fine-zoning bases).
!   3. assign_new_points: the new mesh, evenly spaced in chi between the
!      flagged points, with points added at X and Z jumps and crowded
!      points dropped.
!   4. locate_new_cz_edges: the convection-zone edges on the new mesh.
!   5. interpolate_onto_new_grid: every variable splined onto the new
!      mesh and its shell masses set (the regridding kernels); for
!      rotating models omega re-derived from j and the total angular
!      momentum and rotational kinetic energy re-totalled.
! From the star: stored_vertex_index, the envelope-triangle vertex whose
! integration stored its state where T reached
! envelope_store_log_temperature (0 for none), used when a point is
! added at the surface; h_shell_bot_index and has_h_shell, the base of
! the H-burning shell and whether there is one. To the star:
! reset_triangle, set when outer points were deleted or added;
! total_angular_momentum and total_rotational_ke, re-totalled on the new
! mesh (rotating models). ierr /= 0 on failure.
!   star reads: physics%use_extended_composition
!   star writes: reset_triangle
subroutine rezone(s, ierr)
      implicit none
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
! the rezone's state, from the flagged points to the new mesh
      type(rezone_state), allocatable :: rz

      ierr = 0
      allocate(rz)

      s%reset_triangle = .false.
! The number of species carried, with or without the extended composition.
      rz%num_species_tracked = n_species_basic
      if (s%physics%use_extended_composition) rz%num_species_tracked = n_species_extended
      call check_envelope_temperature_range(s, rz, ierr)
      if (ierr /= 0) return
      call flag_fixed_points(s, rz, ierr)
      if (ierr /= 0) return
      call assign_new_points(s, rz, ierr)
      if (ierr /= 0) return
      call locate_new_cz_edges(s, rz)
      if (ierr /= 0) return
      call interpolate_onto_new_grid(s, rz, ierr)
      if (ierr /= 0) return

      return

end subroutine rezone



! ---------------------------------------------------------------
! Keep the outermost Henyey point inside the envelope temperature window
! tenv0 < log T < tenv1. If it is cooler than tenv0, the outer points
! down to the first one hotter than tenv0 are deleted; if it is hotter
! than tenv1 and a stored envelope state exists, one point is added
! outside it. Either sets star%reset_triangle (and star%log_fit_point_mass_fraction); ierr if the
! whole model is cooler than tenv0.
!   star reads: fit_point_logP, fit_point_logR, log_total_mass,
!     physics%rezone_env_logt0, physics%rezone_env_logt1,
!     stored_envelope_state, stored_vertex_index
!   star writes: log_fit_point_mass_fraction, log_mass, logP, logR,
!     logRho, logT, luminosity_lsun, nz, reset_triangle, xa
subroutine check_envelope_temperature_range(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      integer :: i
      integer :: j
      character(len=message_length) :: lines(4)
! The outer point is cooler than tenv0: delete points down to the first
! one hotter than tenv0.
      if (s%logT(s%nz) < s%physics%rezone_env_logt0) then
         do i = s%nz-1,1,-1
            if (s%logT(i) > s%physics%rezone_env_logt0) then
               lines = unwritten_line
               write(lines,10) s%nz,i
               call core_log(lines)
   10          format(' OUTER POINTS DELETED OLD M =',I5,'  NEW M =',I5)
               s%nz = i
               s%log_fit_point_mass_fraction = s%log_mass(s%nz) - s%log_total_mass
               s%reset_triangle = .true.
               exit
            end if
         end do
         if (i < (1)) then
! The whole model is cooler than tenv0 (unlikely): stop with ierr, which
! run_yrec returns and the command-line wrapper stops on.
            lines = unwritten_line
            lines(1) = ' ERROR: rezone: the hottest shell is cooler than the minimum envelope' &
                 //' temperature tenv0; run stopped'
            call core_report(lines)
            ierr = 1
            return
         end if
! The outer point is hotter than tenv1: add a point outside it at the
! stored envelope state (log m/M = stored_envelope_state(4)), shifting
! log P and log R by the stored state's offsets from the fitting point,
! with Delta log T = Delta log P/4 (added_point_dlogt_dlogp), L and the
! composition copied, and log rho = log P - log T - 8
! (added_point_log_density_offset).
      else if (s%logT(s%nz) > s%physics%rezone_env_logt1.and. &
           s%stored_vertex_index /= 0) then
         s%nz = s%nz + 1
         s%log_fit_point_mass_fraction = s%stored_envelope_state(4)
         s%log_mass(s%nz) = s%log_total_mass + s%log_fit_point_mass_fraction
         s%logP(s%nz) = s%logP(s%nz-1) + &
              (s%stored_envelope_state(1) - s%fit_point_logP(s%stored_vertex_index))
         s%logT(s%nz) = s%logT(s%nz-1) + &
              (s%stored_envelope_state(1) - &
              s%fit_point_logP(s%stored_vertex_index))*added_point_dlogt_dlogp
         s%logR(s%nz) = s%logR(s%nz-1) + &
              (s%stored_envelope_state(3) - s%fit_point_logR(s%stored_vertex_index))
         s%luminosity_lsun(s%nz) = s%luminosity_lsun(s%nz-1)
         do i = 1,rz%num_species_tracked
            s%xa(i,s%nz) = s%xa(i,s%nz-1)
         end do
         s%logRho(s%nz) = s%logP(s%nz) - &
              s%logT(s%nz) - added_point_log_density_offset
         j = s%nz - 1
         lines = unwritten_line
         write(lines,60) j,s%log_mass(j),s%logP(j), &
              s%logT(j),s%logR(j),s%nz,s%log_mass(s%nz), &
              s%logP(s%nz),s%logT(s%nz), &
              s%logR(s%nz)
         call core_log(lines)
   60    format(' OUTER POINT ADDED',I5,F15.10,'  PTR',3F10.6)
         s%reset_triangle = .true.
      end if
end subroutine check_envelope_temperature_range


! ---------------------------------------------------------------
! Flag the points the rezone keeps and does not smooth across, one
! criterion per step, each adding its points with append_flag_point:
!   flag_convection_zone_edges  the edges of the convective core and of
!                               the surface convection zone;
!   flag_h_shell_base           the point below the H-burning shell;
!   flag_jumps                  shells where X or Z jumps by more than
!                               flag_dx or flag_dz, or (rotating models)
!                               log omega by more than flag_dw;
!   flag_envelope_bases         the overshoot base and the base of the
!                               finely zoned region below the surface
!                               convection zone, when flag_jumps left
!                               room in the list.
! sort_flag_points then ends the list with star%nz, sorted and free of
! duplicates. A new criterion is a step of the same form, called before
! sort_flag_points (core/README.md, "Add a mesh criterion"). The list
! holds max_flag_points; append_flag_point drops what does not fit.
!   star reads: none
!   star writes: none
subroutine flag_fixed_points(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      logical :: list_full

      rz%flag_count = 1
      rz%flag_overflow_warned = .false.
      call flag_convection_zone_edges(s, rz)
      call flag_h_shell_base(s, rz)
      call flag_jumps(s, rz, list_full)
      if (.not. list_full) then
         call flag_envelope_bases(s, rz, ierr)
         if (ierr /= 0) return
      end if
      call sort_flag_points(s, rz)
end subroutine flag_fixed_points

! ---------------------------------------------------------------
! The edges of the central and surface convection zones.
!   star reads: core_cz_top_index, envelope_cz_bottom_index, nz
!   star writes: none
subroutine flag_convection_zone_edges(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      if (s%core_cz_top_index > 1) then
         call append_flag_point(rz, s%core_cz_top_index)
      end if
      if (s%envelope_cz_bottom_index < s%nz .and. &
           s%envelope_cz_bottom_index > 1) then
         call append_flag_point(rz, s%envelope_cz_bottom_index)
      end if
end subroutine flag_convection_zone_edges

! ---------------------------------------------------------------
! The point just below the H-burning shell.
!   star reads: h_shell_bot_index, has_h_shell
!   star writes: none
subroutine flag_h_shell_base(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      if (s%has_h_shell.and.s%h_shell_bot_index > 1) then
         call append_flag_point(rz, s%h_shell_bot_index - 1)
      end if
end subroutine flag_h_shell_base

! ---------------------------------------------------------------
! A jump in X between neighbouring shells, else a jump in Z, else
! (rotating models) a jump in log omega. list_full: the scan stopped
! because the list filled up.
!   star reads: nz, omega, physics%chi_grid_scale,
!     physics%rotation_active, xa
!   star writes: none
subroutine flag_jumps(s, rz, list_full)
      use math_lib
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      logical, intent(out) :: list_full
      integer :: i
      real(dp) :: log_omega_top, log_omega_bot
      do i = 2,s%nz
! A jump in X between neighbouring shells,
         if (abs(s%xa(i_h1,i)-s%xa(i_h1,i-1)) > s%physics%chi_grid_scale(ichi_flag_dx)) then
            call append_flag_point(rz, i)
! else a jump in Z,
         else if (abs(s%xa(i_metals,i)-s%xa(i_metals,i-1)) > &
              s%physics%chi_grid_scale(ichi_flag_dz)) then
            call append_flag_point(rz, i)
! else (rotating models) a jump in log omega.
         else if (s%physics%rotation_active) then
            log_omega_top = log10(s%omega(i))
            log_omega_bot = log10(s%omega(i-1))
            if (abs(log_omega_top-log_omega_bot) > s%physics%chi_grid_scale(ichi_flag_dw)) then
               call append_flag_point(rz, i)
            end if
         end if
! Stop scanning once the list is full (append_flag_point warns on the
! first rejected append; one slot stays reserved for the star%nz
! terminator).
         if (rz%flag_count >= max_flag_points) exit
      end do
      list_full = i <= s%nz
end subroutine flag_jumps

! ---------------------------------------------------------------
! The overshoot base and the base of the finely zoned region below the
! surface convection zone. Without a surface convection zone, or in a
! fully convective star, there is nothing to flag. (The legacy
! per-region log P limits PMAX1 to PMAX5 computed here were never read;
! the live ones are in accumulate_chi_spacing.) ierr from
! envelope_overshoot_depth.
!   star reads: convective_flag, envelope_cz_bottom_index, log_mass,
!     logP, logR, logRho, logT, nz, physics%chi_grid_scale,
!     physics%envelope_overshoot_active, xa
!   star writes: none
subroutine flag_envelope_bases(s, rz, ierr)
      use math_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use convection_lib, only: envelope_overshoot_depth
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      integer :: overshoot_base_zone, fine_zone_base
      real(dp) :: overshoot_depth, overshoot_base_radius
      if (.not.s%convective_flag(s%nz)) then
         overshoot_base_zone = s%nz
         fine_zone_base = s%nz
      else if (s%envelope_cz_bottom_index == 1) then
         overshoot_base_zone = 1
         fine_zone_base = 1
      else
! The overshoot base: the convection-zone base without envelope
! overshoot; with it, the deepest point less than the envelope overshoot
! depth alpha_env H_P (convection_lib's envelope_overshoot_depth) below
! the base, as the overshoot moves the zone's lower edge.
         if (.not.s%physics%envelope_overshoot_active) then
            overshoot_base_zone = s%envelope_cz_bottom_index
         else
            call envelope_overshoot_depth(s%xa, s%logRho, s%logP, s%logR, s%log_mass, &
                 s%logT, s%envelope_cz_bottom_index, overshoot_depth, ierr)
            if (ierr /= 0) return
            overshoot_base_radius = exp(ln10*s%logR(s%envelope_cz_bottom_index)) - overshoot_depth
            do overshoot_base_zone = s%envelope_cz_bottom_index-1,1,-1
               if (exp(ln10*s%logR(overshoot_base_zone)) < overshoot_base_radius) exit
            end do
            overshoot_base_zone = overshoot_base_zone + 1
            if (overshoot_base_zone > 1) then
               call append_flag_point(rz, overshoot_base_zone)
            end if
         end if
         if (s%physics%chi_grid_scale(ichi_czbase_fine_width) == 0.0D0) then
            fine_zone_base = overshoot_base_zone
         else
! The base of the finely zoned region: the deepest point with
! log P - log P(czb) <= tol_czbase_fine_width, flagged when it differs
! from the overshoot base.
            if (overshoot_base_zone == 1) then
               fine_zone_base = 1
            else
               do fine_zone_base = overshoot_base_zone-1,1,-1
                  if (s%logP(fine_zone_base) - &
                       s%logP(s%envelope_cz_bottom_index) > s%physics%chi_grid_scale(ichi_czbase_fine_width)) &
                       exit
               end do
               fine_zone_base = fine_zone_base + 1
               if (fine_zone_base /= overshoot_base_zone) then
                  if (fine_zone_base > 1) then
                     call append_flag_point(rz, fine_zone_base)
                  end if
               end if
            end if
         end if
      end if
end subroutine flag_envelope_bases

! ---------------------------------------------------------------
! End the flagged points with star%nz, sort them into ascending order
! and remove the duplicates (sort_kept_points); with the solver
! diagnostics on, list them in the run log.
!   star reads: nz
!   star writes: none
subroutine sort_flag_points(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer :: j
      character(len=message_length), allocatable :: flag_lines(:)
! append_flag_point keeps flag_count <= max_flag_points, so the
! terminator is always in bounds.
      rz%flag_point(rz%flag_count) = s%nz
      call sort_kept_points(rz%flag_point, rz%flag_count)
      if (solver_diagnostics()) then
         allocate(flag_lines(rz%flag_count/20 + 1))
         flag_lines = unwritten_line
         write(flag_lines,185) (rz%flag_point(j),j=1,rz%flag_count)
         call core_log(flag_lines)
      end if
  185 format(1X,'FLAG-POINTS',20I4)
end subroutine sort_flag_points


! ---------------------------------------------------------------
! Append zone to the flagged points, unless that would leave no room for
! the star%nz terminator; the first dropped append of a rezone is
! reported in the run log.
subroutine append_flag_point(rz, zone)
      type(rezone_state), intent(inout) :: rz
      character(len=message_length) :: lines(4)
      integer, intent(in) :: zone
      if (rz%flag_count > max_flag_points - 1) then
         if (.not. rz%flag_overflow_warned) then
            lines = unwritten_line
            write(lines,110) max_flag_points
            call core_log(lines)
  110       format(1X,'MORE THAN ',i4,' FLAG POINTS - REST DROPPED')
            rz%flag_overflow_warned = .true.
         end if
         return
      end if
      rz%flag_point(rz%flag_count) = zone
      rz%flag_count = rz%flag_count + 1
end subroutine append_flag_point


! ---------------------------------------------------------------
! Build the new mesh new_log_mass(1:new_num_zones): points evenly spaced
! in chi between the flagged points, more points where X rises or Z
! falls outward too steeply, and points closer than tol_dm_min in log m
! dropped. ierr past json points or on a spline failure.
!   star reads: none
!   star writes: none
subroutine assign_new_points(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      call seed_new_point_run(s, rz)
      call accumulate_chi_spacing(s, rz)
      call place_points_along_chi(s, rz, ierr)
      if (ierr /= 0) return
! X rising outward by more than tol_dx_max, then Z falling outward by
! more than tol_dz_max.
      call add_points_for_jumps(s, rz, i_h1, 'X', .false., rz%max_x_step, ierr)
      if (ierr /= 0) return
      call add_points_for_jumps(s, rz, i_metals, 'Z', .true., rz%max_z_step, ierr)
      if (ierr /= 0) return
      call drop_crowded_points(s, rz)
end subroutine assign_new_points


! ---------------------------------------------------------------
! Save the start-of-step centre values, Teff and surface L (for the
! HR-motion timestep limit), and set the largest steps between new points
! in L, X and Z (legacy HPMAX), the L one scaled by the largest L in the
! model.
!   star reads: log_Teff, logP, luminosity_lsun, nz,
!     physics%chi_grid_scale
!   star writes: log_Teff_start, logP_start, luminosity_lsun_start,
!     surface_luminosity_lsun_start
subroutine seed_new_point_run(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer :: i
      real(dp) :: luminosity_max
      s%logP_start(1) = s%logP(1)
      s%luminosity_lsun_start(1) = s%luminosity_lsun(1)
      luminosity_max = s%luminosity_lsun(s%nz)
! The start-of-step log Teff (JVS 04/14) and surface L, saved before the
! rezone moves the points.
      s%log_Teff_start = s%log_Teff
      s%surface_luminosity_lsun_start = s%luminosity_lsun(s%nz)
! luminosity_max: the largest L in the model.
      do i = s%nz-1,1,-1
         if (s%luminosity_lsun(i) > luminosity_max) then
            luminosity_max = s%luminosity_lsun(i)
         end if
      end do
      rz%max_luminosity_step = s%physics%chi_grid_scale(ichi_dl_max)*luminosity_max
      rz%max_x_step = s%physics%chi_grid_scale(ichi_dx_max)
      rz%max_z_step = s%physics%chi_grid_scale(ichi_dz_max)
end subroutine seed_new_point_run


! ---------------------------------------------------------------
! Accumulate chi, the running sum of the normalized changes from shell
! to shell (accumulate_chi):
!   chi(1) = 1,
!   chi(j) = chi(j-1) + Delta log m/tol_dm_max
!            + Delta L/(tol_dl_max L_max) - Delta log P/dp_scale,
! the L term only where L increases outward. dp_scale depends on the
! region: tol_dp_czbase_max within tol_czbase_fine_width in log P of the
! surface convection zone's base (czb), tol_dp_core_max deeper, and
! tol_dp_env_max above (in the convection zone). A unit step in chi is
! the largest step allowed between points.
!   star reads: envelope_cz_bottom_index, log_mass, logP,
!     luminosity_lsun, nz, physics%chi_grid_scale
!   star writes: none
subroutine accumulate_chi_spacing(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      type(chi_scales) :: scales

      scales%log_mass = s%physics%chi_grid_scale(ichi_dm_max)
      scales%luminosity = rz%max_luminosity_step
      scales%log_pressure_core = s%physics%chi_grid_scale(ichi_dp_core_max)
      scales%log_pressure_base = s%physics%chi_grid_scale(ichi_dp_czbase_max)
      scales%log_pressure_envelope = s%physics%chi_grid_scale(ichi_dp_env_max)
      scales%base_fine_width = s%physics%chi_grid_scale(ichi_czbase_fine_width)
      call accumulate_chi(s%log_mass(1:s%nz), s%luminosity_lsun(1:s%nz), s%logP(1:s%nz), &
           s%logP(s%envelope_cz_bottom_index), scales, rz%chi(1:s%nz))
end subroutine accumulate_chi_spacing


! ---------------------------------------------------------------
! Place the new points evenly in chi between successive flagged points,
! int(Delta chi) + 1 intervals per segment so that every step is below
! one, reading log m off a natural cubic spline of log m against chi
! (space_points_along_chi). The flagged points themselves are kept (the
! spline passes through its knots). ierr past json points or on a spline
! failure.
!   star reads: log_mass, nz
!   star writes: none
subroutine place_points_along_chi(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      integer :: needed_points, spline_ierr
      character(len=message_length) :: lines(4)

      call space_points_along_chi(rz%chi(1:s%nz), s%log_mass(1:s%nz), &
           rz%flag_point(1:rz%flag_count), rz%new_log_mass, rz%new_num_zones, &
           needed_points, spline_ierr)
      if (needed_points > 0) then
         call report_mesh_overflow(needed_points)
         ierr = 1
      else if (spline_ierr /= 0) then
         lines = unwritten_line
         write(lines,'(a)') ' ERROR: rezone: the spline of log m along the mesh' &
              //' criterion failed: '//trim(numerics_status_text(spline_ierr))
         call core_report(lines)
         ierr = spline_ierr
      end if
end subroutine place_points_along_chi


! ---------------------------------------------------------------
! The error of a mesh that needs more points than the arrays hold; the
! caller returns ierr = 1, and the run stops.
subroutine report_mesh_overflow(needed_points)
      integer, intent(in) :: needed_points
      character(len=message_length) :: lines(4)

      lines = unwritten_line
      write(lines,101) needed_points, json
      call core_report(lines)
  101 format(' ERROR: rezone: the new mesh needs ',I0,' points, more than json (', &
           I0,'); run stopped')
end subroutine report_mesh_overflow


! ---------------------------------------------------------------
! Interpolate one species (slot, called profile_name in the messages: X
! or Z) onto the new points and, where it steps outward by more than
! max_step between two new points (its rise, or with falling its fall),
! insert int(step/max_step) points evenly in log m between them
! (spline_onto_points, insert_points_at_jumps). Nothing is inserted when
! max_step is zero. ierr past json points or on a spline failure.
!   star reads: log_mass, nz, xa
!   star writes: none
subroutine add_points_for_jumps(s, rz, slot, profile_name, falling, max_step, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(in) :: slot
      character(len=*), intent(in) :: profile_name
      logical, intent(in) :: falling
      real(dp), intent(in) :: max_step
      integer, intent(inout) :: ierr
      real(dp) :: values(rz%new_num_zones)
      integer :: needed_points, spline_ierr
      character(len=message_length) :: lines(4)

      call spline_onto_points(s%log_mass(1:s%nz), s%xa(slot,1:s%nz), &
           rz%new_log_mass(1:rz%new_num_zones), values, spline_ierr)
      if (spline_ierr /= 0) then
         lines = unwritten_line
         write(lines,'(a)') ' ERROR: rezone: the spline of '//profile_name// &
              ' on the new mesh failed: '//trim(numerics_status_text(spline_ierr))
         call core_report(lines)
         ierr = spline_ierr
         return
      end if
      call insert_points_at_jumps(values, falling, max_step, rz%new_log_mass, &
           rz%new_num_zones, needed_points)
      if (needed_points > 0) then
         call report_mesh_overflow(needed_points)
         ierr = 1
      end if
end subroutine add_points_for_jumps


! ---------------------------------------------------------------
! Drop the new points that lie within tol_dm_min in log m of the last
! point kept; the first and last points are always kept
! (remove_crowded_points).
!   star reads: physics%chi_grid_scale
!   star writes: none
subroutine drop_crowded_points(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz

      call remove_crowded_points(rz%new_log_mass, rz%new_num_zones, &
           s%physics%chi_grid_scale(ichi_dm_min))
end subroutine drop_crowded_points


! ---------------------------------------------------------------
! Locate the edges of the convection zones on the new mesh: the
! convective core's outer edge becomes the last new point at or inside
! the old edge, and the surface convection zone's base the first new
! point at or outside the old base. Then set convective_flag on the new
! points: convective in the core and the envelope zone, radiative
! between.
!   star reads: log_mass, nz
!   star writes: convective_flag, core_cz_top_index,
!     envelope_cz_bottom_index
subroutine locate_new_cz_edges(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer :: i
      integer :: j
      integer :: radiative_zone_begin
      integer :: radiative_zone_end
      if (s%core_cz_top_index > 1) then
         do j = 2,rz%new_num_zones
            if (rz%new_log_mass(j) > s%log_mass(s%core_cz_top_index)) exit
         end do
         s%core_cz_top_index = j - 1
      else
         s%core_cz_top_index = 1
      end if
      if (s%envelope_cz_bottom_index < s%nz) then
         do j = rz%new_num_zones-1,1,-1
            if (rz%new_log_mass(j) < s%log_mass(s%envelope_cz_bottom_index)) &
                 exit
         end do
         s%envelope_cz_bottom_index = j + 1
      else
         s%envelope_cz_bottom_index = rz%new_num_zones
      end if
      if (s%core_cz_top_index > 1) then
         do i = 1,s%core_cz_top_index
            s%convective_flag(i) = .true.
         end do
         radiative_zone_begin = s%core_cz_top_index + 1
      else
         radiative_zone_begin = 1
      end if
      if (s%envelope_cz_bottom_index < rz%new_num_zones) then
         do i = s%envelope_cz_bottom_index,rz%new_num_zones
            s%convective_flag(i) = .true.
         end do
         radiative_zone_end = s%envelope_cz_bottom_index - 1
      else
         radiative_zone_end = rz%new_num_zones
      end if
      if (radiative_zone_end >= 1 .and. radiative_zone_begin <= rz%new_num_zones) then
         do j = radiative_zone_begin,radiative_zone_end
            s%convective_flag(j) = .false.
         end do
      end if
end subroutine locate_new_cz_edges


! ---------------------------------------------------------------
! Move the model onto the new mesh: spline the composition, structure
! and rotation variables and the start-of-step energy generation onto it
! with the shape-preserving osculatory quadratic spline (osplin, through
! the regridding kernels), adopt the mesh and set its shell masses; for
! rotating models re-spin and re-total the angular momentum; re-run the
! shell physics where the transport needs it; and rebuild the surface
! opacity tables if the surface X moved.
!   star reads: physics%rotation_active
!   star writes: none
subroutine interpolate_onto_new_grid(s, rz, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      character(len=message_length) :: lines(4)
      call regrid_composition_terms(s, rz, ierr)
      if (ierr == 0) call regrid_structure_variables(s, rz, ierr)
      if (ierr == 0) call regrid_rotation_variables(s, rz, ierr)
      if (ierr == 0) call regrid_start_of_step_eps(s, rz, ierr)
      if (ierr /= 0) then
! osplin could not regrid a profile (the regridding kernels' status):
! the run stops with ierr, as for the other rezoning failures.
         lines = unwritten_line
         write(lines,'(a)') ' rezone: osplin: '//numerics_status_text(ierr)
         call core_report(lines)
         return
      end if
      call transfer_new_points(s, rz)
      call set_new_shell_masses(s)
      if (s%physics%rotation_active) then
         call respin_on_new_grid(s, ierr)
         if (ierr /= 0) return
         call retally_angular_momentum(s)
      end if
      call rerun_physics_on_new_grid(s, ierr)
      if (ierr /= 0) return
      call refresh_surface_opacity(s, ierr)
end subroutine interpolate_onto_new_grid


! ---------------------------------------------------------------
! Set the old and new point counts, and spline the entropy terms, every
! species (now and at the start of the step) and, with the extended
! composition, the start-of-step deuterium burning rate onto the new
! points.
!   star reads: log_mass, nz, physics%use_extended_composition
!   star writes: deuterium_burning_rate_start, dln_luminosity_dt,
!     dlog_pressure_dt, dlog_radius_dt, dlog_temperature_dt, xa,
!     xa_start
subroutine regrid_composition_terms(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
! regridded_xa: the species tracked at the new points, before they
! replace star%xa and star%xa_start.
      real(dp), allocatable :: regridded_xa(:,:)
! Each field is splined from the old points (star%log_mass,
! old_point_count of them) to the new ones (new_log_mass,
! new_point_count); regrid_in_place overwrites the field it is given.
      rz%old_point_count = s%nz
      rz%new_point_count = rz%new_num_zones
! The entropy terms (7/91).
      call regrid_in_place(s, rz, ierr, s%dlog_temperature_dt)
      call regrid_in_place(s, rz, ierr, s%dlog_pressure_dt)
      call regrid_in_place(s, rz, ierr, s%dln_luminosity_dt)
      call regrid_in_place(s, rz, ierr, s%dlog_radius_dt)
      if (ierr /= 0) return
      allocate(regridded_xa(rz%num_species_tracked,rz%new_point_count))
      call regrid_composition(s%log_mass(1:rz%old_point_count), &
           s%xa(1:rz%num_species_tracked,1:rz%old_point_count), &
           rz%new_log_mass(1:rz%new_point_count), regridded_xa, ierr)
      if (ierr /= 0) return
      s%xa(1:rz%num_species_tracked,1:rz%new_point_count) = regridded_xa
! xa_start (legacy HCOMPP), the composition at the start of the step, is
! needed for composition diffusion in rotating models.
      call regrid_composition(s%log_mass(1:rz%old_point_count), &
           s%xa_start(1:rz%num_species_tracked,1:rz%old_point_count), &
           rz%new_log_mass(1:rz%new_point_count), regridded_xa, ierr)
      if (ierr /= 0) return
      s%xa_start(1:rz%num_species_tracked,1:rz%new_point_count) = regridded_xa
! The start-of-step deuterium burning rate, as every per-zone array that
! outlives the rezone (the rotation sub-steps read it zone by zone on the
! new mesh).
      if (s%physics%use_extended_composition) then
         call regrid_in_place(s, rz, ierr, s%deuterium_burning_rate_start)
      end if
end subroutine regrid_composition_terms


! ---------------------------------------------------------------
! Spline log P, log R, L, log T and log rho onto the new points, into the
! start-of-step arrays star%*_start. Once transfer_new_points has copied
! them into star%logP etc., the *_start arrays hold the start-of-step
! model on the new mesh (what rezone_or_snapshot stores explicitly when
! it does not rezone), which timestep_limit_structure,
! mid_timestep_model, burn_step and convection_geometry read.
!   star reads: log_mass, logP, logR, logRho, logT, luminosity_lsun
!   star writes: logP_start, logR_start, logRho_start, logT_start,
!     luminosity_lsun_start
subroutine regrid_structure_variables(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      type(structure_profiles) :: regridded
      call regrid_structure(s%log_mass(1:rz%old_point_count), &
           structure_profiles(log_pressure=s%logP(1:rz%old_point_count), &
           log_radius=s%logR(1:rz%old_point_count), &
           luminosity_lsun=s%luminosity_lsun(1:rz%old_point_count), &
           log_temperature=s%logT(1:rz%old_point_count), &
           log_density=s%logRho(1:rz%old_point_count)), &
           rz%new_log_mass(1:rz%new_point_count), regridded, ierr)
      if (ierr /= 0) return
      s%logP_start(1:rz%new_point_count) = regridded%log_pressure
      s%logR_start(1:rz%new_point_count) = regridded%log_radius
      s%luminosity_lsun_start(1:rz%new_point_count) = regridded%luminosity_lsun
      s%logT_start(1:rz%new_point_count) = regridded%log_temperature
      s%logRho_start(1:rz%new_point_count) = regridded%log_density
end subroutine regrid_structure_variables


! ---------------------------------------------------------------
! For rotating models, spline omega, the specific angular momentum j,
! the shape factors fp and ft, r0 and eta2 onto the new points, into the
! *_start arrays (and fp_splined, ft_splined) that transfer_new_points adopts.
!   star reads: eta_squared, fp_rot, ft_rot, j_rot, log_mass,
!     mean_radius, omega, physics%rotation_active
!   star writes: eta_squared_start, j_rot_start, mean_radius_start,
!     omega_start
subroutine regrid_rotation_variables(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      type(rotation_profiles) :: regridded
      if (s%physics%rotation_active) then
         call regrid_rotation(s%log_mass(1:rz%old_point_count), &
              rotation_profiles(omega=s%omega(1:rz%old_point_count), &
              specific_angular_momentum=s%j_rot(1:rz%old_point_count), &
              fp=s%fp_rot(1:rz%old_point_count), &
              ft=s%ft_rot(1:rz%old_point_count), &
              mean_radius=s%mean_radius(1:rz%old_point_count), &
              eta_squared=s%eta_squared(1:rz%old_point_count)), &
              rz%new_log_mass(1:rz%new_point_count), regridded, ierr)
         if (ierr /= 0) return
         s%omega_start(1:rz%new_point_count) = regridded%omega
         s%j_rot_start(1:rz%new_point_count) = &
              regridded%specific_angular_momentum
         rz%fp_splined(1:rz%new_point_count) = regridded%fp
         rz%ft_splined(1:rz%new_point_count) = regridded%ft
         s%mean_radius_start(1:rz%new_point_count) = regridded%mean_radius
         s%eta_squared_start(1:rz%new_point_count) = regridded%eta_squared
      end if

end subroutine regrid_rotation_variables


! ---------------------------------------------------------------
! For rotating models (or the extended composition with envelope
! overshoot), spline the start-of-step energy generation onto the new
! points (MHP 6/00): eps_total into eps_total_start and
! eps_total + eps_neu + eps_grav into eps_plus_neu_grav_start, for the
! mid-step Eddington-Sweet velocities. This runs before
! transfer_new_points, while star%log_mass is still the old mesh. (2026
! fix, bugsweep sec-11: hpoint.f did this after the transfer, against
! the new mesh, so the values were misplaced on every rotating rezone.)
!   star reads: eps_channels, eps_total, log_mass, nz,
!     physics%envelope_overshoot_active, physics%rotation_active,
!     physics%use_extended_composition
!   star writes: eps_plus_neu_grav_start, eps_total_start
subroutine regrid_start_of_step_eps(s, rz, ierr)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      real(dp) :: spline_y(json)
      integer :: zone_index
      if (s%physics%rotation_active .or. (s%physics%use_extended_composition .and. &
           s%physics%envelope_overshoot_active)) then
         call regrid_profile(s%log_mass(1:rz%old_point_count), &
              s%eps_total(1:rz%old_point_count), rz%new_log_mass(1:rz%new_point_count), &
              s%eps_total_start(1:rz%new_point_count), ierr)
         if (ierr /= 0) return
         do zone_index = 1,s%nz
            spline_y(zone_index) = s%eps_total(zone_index)+s%eps_channels(i_eps_neu,zone_index)+ &
                 s%eps_channels(i_eps_grav,zone_index)
         end do
         call regrid_profile(s%log_mass(1:rz%old_point_count), &
              spline_y(1:rz%old_point_count), rz%new_log_mass(1:rz%new_point_count), &
              s%eps_plus_neu_grav_start(1:rz%new_point_count), ierr)
      end if
end subroutine regrid_start_of_step_eps


! ---------------------------------------------------------------
! Adopt the new points: the new mesh and the splined structure and
! rotation variables become the model, with the new point count.
!   star reads: eta_squared_start, j_rot_start, logP_start, logR_start,
!     logRho_start, logT_start, luminosity_lsun_start,
!     mean_radius_start, omega_start, physics%rotation_active
!   star writes: eta_squared, fp_rot, ft_rot, j_rot, log_mass, logP,
!     logR, logRho, logT, luminosity_lsun, mean_radius, nz, omega
subroutine transfer_new_points(s, rz)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer :: j
      character(len=message_length) :: lines(4)
      do j = 1,rz%new_num_zones
         s%log_mass(j) = rz%new_log_mass(j)
         s%logP(j) = s%logP_start(j)
         s%logT(j) = s%logT_start(j)
         s%logR(j) = s%logR_start(j)
         s%luminosity_lsun(j) = s%luminosity_lsun_start(j)
         s%logRho(j) = s%logRho_start(j)
      end do
      if (s%physics%rotation_active) then
         do j = 1, rz%new_num_zones
            s%j_rot(j) = s%j_rot_start(j)
            s%omega(j) = s%omega_start(j)
            s%fp_rot(j) = rz%fp_splined(j)
            s%ft_rot(j) = rz%ft_splined(j)
            s%eta_squared(j) = s%eta_squared_start(j)
            s%mean_radius(j) = s%mean_radius_start(j)
         end do
      end if
      if (solver_diagnostics()) then
         lines = unwritten_line
         write(lines,1020) s%nz,rz%new_num_zones
         call core_log(lines)
      end if
 1020 format(' POINTS  OLD',I5,'   NEW',I5)
      s%nz = rz%new_num_zones
end subroutine transfer_new_points


! ---------------------------------------------------------------
! Set the point masses m_i = 10^(log m_i) and the mass weights
!   dm_i  = (m_(i+1) - m_(i-1))/2,  with m_0 = -m_1,
!   dm_nz = M - (m_(nz-1) + m_nz)/2,
! of the new mesh (compute_shell_masses), so dm_i is the mass between
! the midpoints on either side of point i (from the centre for i = 1),
! and dm_nz includes the envelope above the outer point
! (M = 10^log_total_mass).
!   star reads: log_mass, log_total_mass, nz
!   star writes: dm, m
subroutine set_new_shell_masses(s)
      type(star_info), intent(inout) :: s
      call compute_shell_masses(s%log_mass(1:s%nz), s%log_total_mass, &
           s%m(1:s%nz), s%dm(1:s%nz))
end subroutine set_new_shell_masses


! ---------------------------------------------------------------
! Rotating models: on the new mesh, rebuild the moments of inertia, the
! convection regions for angular momentum transport, omega from the
! splined specific angular momentum j, and the shape factors; ierr from
! the regions or the shape factors.
!   star reads: disk_gate_age_gyr, dm, log_mass, logP, logR, logRho,
!     logT, nz, xa
!   star writes: am_transport_convective_flag, convective_flag,
!     envelope_overshoot_bottom_index, envelope_overshoot_top_index,
!     eta_squared, fp_rot, ft_rot, i_rot, j_rot, mean_gravity,
!     mean_radius, omega
subroutine respin_on_new_grid(s, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      integer, intent(inout) :: ierr
      integer :: radiative_zone_bounds(max_radiative_zones,2)
      integer :: convective_zone_bounds(max_convective_zones,2)
      integer :: i
      integer :: num_radiative_zones
      integer :: num_convective_zones
! di_domega of omega_from_j, which this caller does not use
      real(dp) :: unused_di_domega(json)
! First guess at the moments of inertia, I_i = (2/3) dm_i r_i^2
! (spherical shells).
      do i=1,s%nz
         s%i_rot(i) = cc23*s%dm(i)* &
              exp(ln10*2.0D0*s%logR(i))
      end do
! The convection zones, extended by overshoot, for angular momentum
! transport (am_transport_convective_flag).
      call am_convective_regions(s%xa,s%logRho,s%logP,s%logR, &
           s%log_mass,s%logT,s%convective_flag,s%nz, &
           s%am_transport_convective_flag,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones, &
           s%envelope_overshoot_bottom_index, s%envelope_overshoot_top_index, ierr)
      if (ierr /= 0) return
! The moments of inertia of the distorted shells (with eta2 and r0)
! before omega is re-derived, as the 2015/04/06 version did (JNT
! 2025/09/03). It may be redundant, since the convection-zone rotation
! law (walpcz) runs in this version, but it is harmless.
      call zone_moments_of_inertia(s%eta_squared,s%logR,s%log_mass,s%dm,1,s%nz, &
           s%omega,s%mean_radius,s%i_rot,unused_di_domega)

! Omega from the specific angular momentum j (legacy HJM).
      call omega_from_j(s%logRho,s%j_rot,s%logR, &
           s%log_mass,s%dm,s%am_transport_convective_flag,s%nz, &
           s%eta_squared,s%i_rot,s%omega,unused_di_domega,s%mean_radius,s%disk_gate_age_gyr)
! The shape factors fp and ft, r0 and eta2 for this omega.
      call rotation_shape_factors(s%logRho,s%logR,s%log_mass,s%nz,s%omega, &
           s%eta_squared,s%fp_rot,s%ft_rot,s%mean_gravity,s%mean_radius,ierr)
      if (ierr /= 0) return
! The moments of inertia for the new shape.
      call zone_moments_of_inertia(s%eta_squared,s%logR,s%log_mass,s%dm,1,s%nz, &
           s%omega,s%mean_radius,s%i_rot,unused_di_domega)
end subroutine respin_on_new_grid


! ---------------------------------------------------------------
! Rotating models: re-total the angular momentum and the rotational
! kinetic energy on the new mesh (compute_angular_momentum_totals),
!   J = sum_i j_i dm_i,   E_rot = sum_i (1/2) omega_i j_i dm_i,
! and store the rotation state as the start of the step.
!   star reads: convective_flag, dm, eta_squared, j_rot, mean_gravity,
!     mean_radius, nz, omega
!   star writes: convective_flag_start, eta_squared_start, j_rot_start,
!     kinetic_energy_rot, mean_gravity_start, mean_radius_start,
!     omega_start, total_angular_momentum, total_rotational_ke
subroutine retally_angular_momentum(s)
      type(star_info), intent(inout) :: s
      integer :: i
      real(dp) :: sum_angular_momentum
      real(dp) :: sum_rotational_ke
      character(len=message_length) :: lines(4)
      call compute_angular_momentum_totals(s%j_rot(1:s%nz), &
           s%omega(1:s%nz), s%dm(1:s%nz), &
           s%kinetic_energy_rot(1:s%nz), sum_angular_momentum, &
           sum_rotational_ke)
! The totals before and after the rezone, a solver diagnostic.
      if (solver_diagnostics()) then
         lines = unwritten_line
         write(lines,1120)s%total_angular_momentum, &
              sum_angular_momentum,s%total_rotational_ke,sum_rotational_ke
         call core_log(lines)
      end if
 1120 format(1X,'TOTAL J OF STAR - PREVIOUS ',1PE21.13,' NEW ', &
           1PE21.13/' TOTAL ROTATIONAL K.E. OF STAR-PREVIOUS ',1PE21.13, &
           ' NEW ',1PE21.13)
      s%total_angular_momentum = sum_angular_momentum
      s%total_rotational_ke = sum_rotational_ke
! The start-of-step rotation state, for the angular momentum transport.
      do i = 1,s%nz
         s%omega_start(i) = s%omega(i)
         s%mean_gravity_start(i) = s%mean_gravity(i)
         s%eta_squared_start(i) = s%eta_squared(i)
         s%mean_radius_start(i) = s%mean_radius(i)
         s%convective_flag_start(i) = s%convective_flag(i)
! j also goes to j_rot_start (legacy HJX), for I/O
! as well (MHP 10/91).
         s%j_rot_start(i) = s%j_rot(i)
      end do
end subroutine retally_angular_momentum


! ---------------------------------------------------------------
! For rotating models, or with the extended composition and envelope
! overshoot, re-run the shell physics on the new points and store the
! start-of-step quantities the transport needs; ierr if the physics
! fails. The temperature gradients are needed with overshoot so that the
! base of the overshoot region of pre-main-sequence models can be located
! accurately (MHP 9/91); even without diffusion, the new physics locates
! the convection-zone edges accurately.
!   star reads: conv_vel, cp, dlnrho_dlnt, grada, gradr, gradT,
!     kinematic_viscosity, mu, nz, opacity_zone,
!     physics%envelope_overshoot_active, physics%rotation_active,
!     physics%use_extended_composition, thermal_diffusivity
!   star writes: conv_vel_start, cp_start, dlnrho_dlnt_start,
!     grada_start, gradr_start, gradT_start, kinematic_viscosity_start,
!     mu_start, opacity_zone_start, thermal_diffusivity_start
subroutine rerun_physics_on_new_grid(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(inout) :: ierr
      integer :: zone_index
      integer :: jerr
      if (s%physics%rotation_active .or. (s%physics%use_extended_composition .and. &
           s%physics%envelope_overshoot_active)) then
         call compute_shell_physics(s, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
! The quantities for the velocities at the start of the step on the new
! mesh, so that several short diffusion sub-steps can be taken within
! one evolutionary step.
         do zone_index = 1,s%nz
            s%gradr_start(zone_index) = s%gradr(zone_index)
            s%gradT_start(zone_index) = s%gradT(zone_index)
            s%grada_start(zone_index) = s%grada(zone_index)
            s%mu_start(zone_index) = s%mu(zone_index)
            s%opacity_zone_start(zone_index) = s%opacity_zone(zone_index)
            s%cp_start(zone_index) = s%cp(zone_index)
            s%dlnrho_dlnt_start(zone_index) = s%dlnrho_dlnt(zone_index)
            s%conv_vel_start(zone_index) = s%conv_vel(zone_index)
            s%kinematic_viscosity_start(zone_index) = s%kinematic_viscosity(zone_index)
            s%thermal_diffusivity_start(zone_index) = s%thermal_diffusivity(zone_index)
! The start-of-step eps sum (legacy ESUMO) is splined in
! regrid_start_of_step_eps (MHP 6/00). (Until 2026-10 this loop also
! stored grad_ad - grad (MHP 06/02) in the rotation scratch's
! del_grad_diff_interface,
! which rotation_stability_setup overwrites before anything reads it.)
         end do
! (The advection-diffusion transport's theta term on the new mesh is in
! model/retired/rezoning_theta_term.f90.txt.)
      end if
end subroutine rerun_physics_on_new_grid


! ---------------------------------------------------------------
! Rebuild the surface opacity tables (for the new surface X and Z) when
! the surface X has moved by more than surface_x_refresh_tol since they
! were built (settings_install's install_surface_x_settings).
!   star reads: nz, xa
!   star writes: surface_hydrogen_fraction, surface_metal_fraction
subroutine refresh_surface_opacity(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(inout) :: ierr
      if (abs(s%surface_hydrogen_fraction-s%xa(i_h1,s%nz)) > surface_x_refresh_tol) then
         s%surface_hydrogen_fraction = s%xa(i_h1,s%nz)
         s%surface_metal_fraction = s%xa(i_metals,s%nz)
         call install_surface_x_settings(s%surface_hydrogen_fraction, ierr)
         if (ierr /= 0) return

      end if
end subroutine refresh_surface_opacity


!----------------------------------------------------------------------
! regrid_in_place
!----------------------------------------------------------------------
! Spline one field from the old mesh (star%log_mass, old_point_count
! points) onto the new one (new_log_mass, new_point_count points) with
! regrid_profile, overwriting field(1:new_point_count). regrid_out holds
! the result until the old values have been read (the legacy HTO). After
! a failure (rezone's ierr set) it does nothing, so that a run of calls
! stops at the first.
!   star reads: log_mass
!   star writes: none
subroutine regrid_in_place(s, rz, ierr, field)
      type(star_info), intent(inout) :: s
      type(rezone_state), intent(inout) :: rz
      integer, intent(inout) :: ierr
      real(dp) :: regrid_out(json)
      real(dp), intent(inout) :: field(:)

      if (ierr /= 0) return
      call regrid_profile(s%log_mass(1:rz%old_point_count), field(1:rz%old_point_count), &
           rz%new_log_mass(1:rz%new_point_count), regrid_out(1:rz%new_point_count), ierr)
      if (ierr == 0) field(1:rz%new_point_count) = regrid_out(1:rz%new_point_count)
end subroutine regrid_in_place


end module rezoning
