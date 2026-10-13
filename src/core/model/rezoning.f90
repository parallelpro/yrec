!----------------------------------------------------------------------
! rezoning
!----------------------------------------------------------------------
! The rezone of each step (rezone): a new mesh in log m, its points
! spaced evenly in a coordinate chi built from the changes of log m, L
! and log P, with points added at X and Z jumps; the model interpolated
! onto it and, for rotating models, the rotation re-spun on it. The
! regridding arithmetic (the shell masses of a mesh, the interpolation
! of the profiles, the angular momentum totals) is the star-blind kernel
! in regridding.f90; rezone gathers its inputs from the star and stores
! its results.
module rezoning
      use shell_physics, only: compute_shell_physics
      use star_info_lib, only: star_info
      implicit none
      private
      public :: rezone

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
! envelope_store_index: the envelope-triangle vertex whose integration
! stored its state where T reached tenv (0 for none), used when a point
! is added at the surface. point_reset_flag: set when outer points were
! deleted or added. h_shell_zone_begin, h_shell_active: the base of the
! H-burning shell and whether there is one. total_angular_momentum,
! total_rotational_ke: re-totalled on the new mesh (rotating models).
! ierr /= 0 on failure.
!   star reads: conv_vel, cp, ctrl%chi_grid_scale,
!     ctrl%overshoot_alpha_envelope, ctrl%tenv0, ctrl%tenv1,
!     disk_gate_age_gyr, eps_channels, eps_total, fit_point_pressure,
!     fit_point_radius, grada, gradr, gradT,
!     job%envelope_overshoot_active, job%rotation_active,
!     job%use_extended_composition, log_Teff, log_total_mass, mu,
!     opacity_zone, qdt, stored_envelope_state, thdif, visc
!   star writes: am_transport_convective_flag, convective_flag,
!     convective_flag_start, core_cz_top_index,
!     deuterium_burning_rate_start, dln_luminosity_dt, dlog_pressure_dt,
!     dlog_radius_dt, dlog_temperature_dt, dm, envelope_cz_bottom_index,
!     eta_squared, fp_rot, ft_rot, i_rot, iov1, iov2, j_rot,
!     kinetic_energy_rot, log_mass, log_Teff_start, logP, logP_start,
!     logR, logR_start, logRho, logRho_start, logT, logT_start,
!     luminosity_lsun, luminosity_lsun_start, m, mean_gravity,
!     mean_radius, nz, old_conv_vel, old_cp, old_eps_plus_neu_grav,
!     old_eps_total, old_eta_squared, old_grada, old_gradr, old_gradT,
!     old_hg, old_mean_radius, old_mu, old_omega, old_opacity, old_qdt,
!     old_specific_angular_momentum, old_thdif, old_visc, omega, qiw,
!     senv, surface_luminosity_lsun_start, xa, xa_start, xnew, znew
subroutine rezone(s, envelope_store_index, point_reset_flag, &
     h_shell_zone_begin, h_shell_active, total_angular_momentum, &
     total_rotational_ke, ierr)
      use am_regions, only: am_convective_regions
      use omega_profile, only: omega_from_j
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use star_info_lib, only: i_eps_grav, i_eps_neu, i_h1, i_h2, i_metals, json
      use controls_buffer, only: ichi_czbase_fine_width, ichi_dl_max, ichi_dm_max, &
           ichi_dm_min, ichi_dp_core_max, ichi_dp_czbase_max, ichi_dp_env_max, ichi_dx_max, &
           ichi_dz_max, ichi_flag_dw, ichi_flag_dx, ichi_flag_dz
      use kap_lib
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use numerics_lib, only: natural_spline, splint
      use regridding, only: structure_profiles, rotation_profiles, &
           compute_shell_masses, regrid_profile, regrid_composition, &
           regrid_structure, regrid_rotation, compute_angular_momentum_totals

      implicit none
      type(star_info), intent(inout) :: s
      integer :: jerr_gate

      integer, intent(in) :: envelope_store_index
      logical, intent(out) :: point_reset_flag
      integer, intent(in) :: h_shell_zone_begin
      logical, intent(in) :: h_shell_active
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke

      double precision :: point_spacing_max(4)
! flag_point(1:flag_count): the flagged points, at most max_flag_points.
! append_flag_point guards every append (2026 audit fix: the legacy code
! guarded only the gradient scan, so the appends outside it could write
! past the array).
      integer, parameter :: max_flag_points = 100
      integer :: flag_point(max_flag_points)
      logical :: flag_overflow_warned
! ft_old, fp_old: the shape factors splined onto the new mesh, which
! transfer_new_points adopts.
      double precision :: ft_old(json), fp_old(json)
      integer :: radiative_zone_bounds(13,2), convective_zone_bounds(12,2)
      double precision :: z_new(json), x_new(json)
! spline_x, spline_y, spline_second_deriv: a spline table and its second
! derivatives. chi (MHP 7/02): the running sum of the normalized changes
! of log m, L and log P from shell to shell, which sets the spacing.
      double precision :: spline_x(json), spline_y(json), chi(json), &
           spline_second_deriv(json)
! new_log_mass: the new mesh in log m (legacy HSS), built in
! assign_new_points and adopted as star%log_mass in transfer_new_points.
! compact_log_mass: the buffer of drop_crowded_points. regrid_out: the
! regridded field of regrid_in_place. (The legacy code kept these, and
! the copy of the field it splined, in a star member since removed,
! logRho_start, logP_start and logT_start; every element is written here
! before it is read.)
      double precision :: new_log_mass(json), compact_log_mass(json), &
           regrid_out(json)
      integer :: gradient_flag_index(json)
      integer :: num_species_tracked, i, j, k
      integer :: flag_count
      double precision :: pressure_test, dp_scale
      integer :: new_num_zones, working_num_zones, num_new_points
      double precision :: delta_x_over_max, delta_z_over_max, &
           point_insert_spacing
      double precision :: mass_scale, luminosity_scale
      double precision :: dchi
      double precision :: chi_prev
      integer :: point_prev_index
      double precision :: delta_chi
      integer :: segment_point_count
      double precision :: spline_eval_x, spline_eval_y
      integer :: gradient_flag_count
      double precision :: luminosity_max
      double precision :: log_omega_top, log_omega_bot
      integer :: overshoot_base_zone, fine_zone_base
      logical :: sort_done
      integer :: temp_swap
      integer :: old_point_count, new_point_count
      integer :: zone_index
      double precision :: sum_angular_momentum, sum_rotational_ke
! The theta term's locals, commented out with the advection-diffusion
! transport (see rerun_physics_on_new_grid).
!      double precision :: delta_radius, delta_omega, omega_mid, log_factor
      integer :: radiative_zone_begin, radiative_zone_end
      integer :: num_radiative_zones, num_convective_zones

      ! jerr: the error code of compute_shell_physics, passed on in ierr.
      integer :: jerr

      integer, intent(out) :: ierr

      ierr = 0

      point_reset_flag = .false.
! The number of species carried: 11, or 15 with the extended composition.
      num_species_tracked = 11
      if (s%job%use_extended_composition) num_species_tracked = 15
      call check_envelope_temperature_range
      if (ierr /= 0) return
      call flag_fixed_points
      if (ierr /= 0) return
      call assign_new_points
      if (ierr /= 0) return
      call locate_new_cz_edges
      if (ierr /= 0) return
      call interpolate_onto_new_grid
      if (ierr /= 0) return

      return

contains

! ---------------------------------------------------------------
! Keep the outermost Henyey point inside the envelope temperature window
! tenv0 < log T < tenv1. If it is cooler than tenv0, the outer points
! down to the first one hotter than tenv0 are deleted; if it is hotter
! than tenv1 and a stored envelope state exists, one point is added
! outside it. Either sets point_reset_flag (and star%senv); ierr if the
! whole model is cooler than tenv0.
subroutine check_envelope_temperature_range
! The outer point is cooler than tenv0: delete points down to the first
! one hotter than tenv0.
      if (s%logT(s%nz).lt.s%ctrl%tenv0) then
       do i = s%nz-1,1,-1
          if (s%logT(i).gt.s%ctrl%tenv0) then
             write(run_log_unit,10) s%nz,i
   10     format(' OUTER POINTS DELETED OLD M =',I5,'  NEW M =',I5)
             s%nz = i
             s%senv = s%log_mass(s%nz) - s%log_total_mass
             point_reset_flag = .true.
             exit
          endif
       end do
       if (i < (1)) then
! The whole model is cooler than tenv0 (unlikely): stop with ierr, which
! run_yrec returns and the command-line wrapper stops on.
       write(run_log_unit,30)
       write(terminal_unit,30)
   30    format(1X,39('>'),40('<')/1X,'ERROR IN HPOINT'/ &
     1X,'MAX. STAR T LESS THAN MINIMUM ENVELOPE T.RUN STOPPED')
       ierr = 1
       return
       end if
! The outer point is hotter than tenv1: add a point outside it at the
! stored envelope state (log m/M = stored_envelope_state(4)), shifting
! log P and log R by the stored state's offsets from the fitting point,
! with Delta log T = Delta log P/4, L and the composition copied, and
! log rho = log P - log T - 8 (a rough ideal-gas estimate).
      else if (s%logT(s%nz).gt.s%ctrl%tenv1.and. &
           envelope_store_index.ne.0) then
       s%nz = s%nz + 1
       s%senv = s%stored_envelope_state(4)
       s%log_mass(s%nz) = s%log_total_mass + s%senv
       s%logP(s%nz) = s%logP(s%nz-1) + &
            (s%stored_envelope_state(1) - s%fit_point_pressure(envelope_store_index))
       s%logT(s%nz) = s%logT(s%nz-1) + &
            (s%stored_envelope_state(1) - &
            s%fit_point_pressure(envelope_store_index))*0.250D0
       s%logR(s%nz) = s%logR(s%nz-1) + &
            (s%stored_envelope_state(3) - s%fit_point_radius(envelope_store_index))
       s%luminosity_lsun(s%nz) = s%luminosity_lsun(s%nz-1)
       do i = 1,num_species_tracked
          s%xa(i,s%nz) = s%xa(i,s%nz-1)
       end do
       s%logRho(s%nz) = s%logP(s%nz) - &
            s%logT(s%nz) - 8.0D0
       j = s%nz - 1
       write(run_log_unit,60) j,s%log_mass(j),s%logP(j), &
            s%logT(j),s%logR(j),s%nz,s%log_mass(s%nz), &
            s%logP(s%nz),s%logT(s%nz), &
            s%logR(s%nz)
   60    format(' OUTER POINT ADDED',I5,F15.10,'  PTR',3F10.6)
       point_reset_flag = .true.
      endif
end subroutine check_envelope_temperature_range

! ---------------------------------------------------------------
! Flag the points the rezone keeps and does not smooth across: the edges
! of the convective core and of the surface convection zone, the point
! below the H-burning shell, shells where X or Z jumps by more than
! flag_dx or flag_dz (or, rotating models, log omega by more than
! flag_dw), the overshoot base and the base of the finely zoned region
! below the surface convection zone. The list ends with star%nz and is
! sorted and free of duplicates.
subroutine flag_fixed_points
      use math_lib
      flag_count = 1
      flag_overflow_warned = .false.
! The edges of the central and surface convection zones.
      if (s%core_cz_top_index.gt.1) then
       call append_flag_point(s%core_cz_top_index)
      endif
      if (s%envelope_cz_bottom_index.lt.s%nz .and. &
           s%envelope_cz_bottom_index.gt.1) then
       call append_flag_point(s%envelope_cz_bottom_index)
      endif
! The point just below the H-burning shell.
      if (h_shell_active.and.h_shell_zone_begin.gt.1) then
       call append_flag_point(h_shell_zone_begin - 1)
      endif
      do i = 2,s%nz
! A jump in X between neighbouring shells,
       if (dabs(s%xa(i_h1,i)-s%xa(i_h1,i-1)).gt.s%ctrl%chi_grid_scale(ichi_flag_dx)) then
          call append_flag_point(i)
! else a jump in Z,
       else if (dabs(s%xa(i_metals,i)-s%xa(i_metals,i-1)).gt. &
            s%ctrl%chi_grid_scale(ichi_flag_dz)) then
          call append_flag_point(i)
! else (rotating models) a jump in log omega.
       else if (s%job%rotation_active) then
          log_omega_top = log10(s%omega(i))
          log_omega_bot = log10(s%omega(i-1))
          if (dabs(log_omega_top-log_omega_bot).gt.s%ctrl%chi_grid_scale(ichi_flag_dw)) then
             call append_flag_point(i)
          endif
       endif
! Stop scanning once the list is full (append_flag_point warns on the
! first rejected append; one slot stays reserved for the star%nz
! terminator).
       if (flag_count.ge.max_flag_points) exit
      end do
! The overshoot-base and fine-zone-base flags are added only when the
! scan above ran to completion (i = nz+1), i.e. the list is not full.
! Without a surface convection zone, or in a fully convective star,
! there is nothing to flag. (The legacy per-region log P limits PMAX1 to
! PMAX5 computed here were never read; the live ones are in
! accumulate_chi_spacing.)
      if (i > (s%nz)) then
      if (.not.s%convective_flag(s%nz)) then
       overshoot_base_zone = s%nz
       fine_zone_base = s%nz
      else if (s%envelope_cz_bottom_index.eq.1) then
       overshoot_base_zone = 1
       fine_zone_base = 1
      else
! The overshoot base: the convection-zone base without envelope
! overshoot; with it, the deepest point with
! log P - log P(czb) <= overshoot_alpha_envelope (alpha used here as a
! depth in log P).
       if (.not.s%job%envelope_overshoot_active) then
          overshoot_base_zone = s%envelope_cz_bottom_index
       else
          do overshoot_base_zone = s%envelope_cz_bottom_index-1,1,-1
             if (s%logP(overshoot_base_zone)- &
                  s%logP(s%envelope_cz_bottom_index).gt.s%ctrl%overshoot_alpha_envelope) exit
          end do
            overshoot_base_zone = overshoot_base_zone + 1
            if (overshoot_base_zone.gt.1) then
               call append_flag_point(overshoot_base_zone)
            endif
       endif
       if (s%ctrl%chi_grid_scale(ichi_czbase_fine_width).eq.0.0D0) then
          fine_zone_base = overshoot_base_zone
       else
! The base of the finely zoned region: the deepest point with
! log P - log P(czb) <= tol_czbase_fine_width, flagged when it differs
! from the overshoot base.
          if (overshoot_base_zone.eq.1) then
             fine_zone_base = 1
          else
          do fine_zone_base = overshoot_base_zone-1,1,-1
             if (s%logP(fine_zone_base) - &
                  s%logP(s%envelope_cz_bottom_index).gt.s%ctrl%chi_grid_scale(ichi_czbase_fine_width)) &
                  exit
          end do
            fine_zone_base = fine_zone_base + 1
          if (fine_zone_base.ne.overshoot_base_zone) then
            if (fine_zone_base.gt.1) then
               call append_flag_point(fine_zone_base)
            endif
          end if
          end if
       endif
      endif
      end if
! append_flag_point keeps flag_count <= max_flag_points, so the
! terminator is always in bounds.
      flag_point(flag_count) = s%nz
! Sort the flagged points into ascending order (bubble sort).
      if (flag_count.ne.1) then
      do
      sort_done = .true.
      do i = 1,flag_count-1
       if (flag_point(i+1).lt.flag_point(i)) then
          temp_swap = flag_point(i)
          flag_point(i) = flag_point(i+1)
          flag_point(i+1) = temp_swap
          sort_done = .false.
       endif
      end do
      if (sort_done) exit
      end do
! Remove the duplicates.
      i = 2
      do
      if (flag_point(i).eq.flag_point(i-1)) then
       if (i.lt.flag_count) then
          do j = i,flag_count-1
             flag_point(j) = flag_point(j+1)
          end do
       endif
       flag_count = flag_count - 1
      endif
      i = i + 1
      if (i.gt.flag_count) exit
      end do
      end if
      if (solver_diagnostics()) then
         write(run_log_unit,185) (flag_point(j),j=1,flag_count)
      end if
  185 format(1X,'FLAG-POINTS',20I4)
end subroutine flag_fixed_points

! ---------------------------------------------------------------
! Append zone to the flagged points, unless that would leave no room for
! the star%nz terminator; the first dropped append of a rezone is
! reported in the run log.
subroutine append_flag_point(zone)
      integer, intent(in) :: zone
      if (flag_count .gt. max_flag_points - 1) then
         if (.not. flag_overflow_warned) then
            write(run_log_unit,110) max_flag_points
  110       format(1X,'MORE THAN ',i4,' FLAG POINTS - REST DROPPED')
            flag_overflow_warned = .true.
         end if
         return
      end if
      flag_point(flag_count) = zone
      flag_count = flag_count + 1
end subroutine append_flag_point

! ---------------------------------------------------------------
! Build the new mesh new_log_mass(1:new_num_zones): points evenly spaced
! in chi between the flagged points, more points where X rises or Z
! falls outward too steeply, and points closer than tol_dm_min in log m
! dropped. ierr past json points or on a spline failure.
subroutine assign_new_points
      use math_lib
      call seed_new_point_run
      call accumulate_chi_spacing
      call place_points_along_chi
      if (ierr /= 0) return
      call add_points_for_x_jumps
      if (ierr /= 0) return
      call add_points_for_z_jumps
      if (ierr /= 0) return
      call drop_crowded_points
end subroutine assign_new_points

! ---------------------------------------------------------------
! Start the new mesh at the centre, save the start-of-step Teff and
! surface L (for the HR-motion timestep limit), and set the maximum spacings
! point_spacing_max in log P, L, X and Z, the L one scaled by the largest
! L in the model.
subroutine seed_new_point_run
! The centre point stays where it is.
      new_log_mass(1) = s%log_mass(1)
      s%logP_start(1) = s%logP(1)
      s%luminosity_lsun_start(1) = s%luminosity_lsun(1)
      x_new(1) = s%xa(i_h1,1)
      z_new(1) = s%xa(i_metals,1)
      luminosity_max = s%luminosity_lsun(s%nz)
! The start-of-step log Teff (JVS 04/14) and surface L, saved before the
! rezone moves the points.
      s%log_Teff_start = s%log_Teff
      s%surface_luminosity_lsun_start = s%luminosity_lsun(s%nz)
! luminosity_max: the largest L in the model.
      do i = s%nz-1,1,-1
         if (s%luminosity_lsun(i).gt.luminosity_max) then
            luminosity_max = s%luminosity_lsun(i)
         endif
      end do
      point_spacing_max(1) = s%ctrl%chi_grid_scale(ichi_dp_env_max)
      point_spacing_max(2) = s%ctrl%chi_grid_scale(ichi_dl_max)*luminosity_max
      point_spacing_max(3) = s%ctrl%chi_grid_scale(ichi_dx_max)
      point_spacing_max(4) = s%ctrl%chi_grid_scale(ichi_dz_max)
end subroutine seed_new_point_run

! ---------------------------------------------------------------
! Accumulate chi, the running sum of the normalized changes from shell
! to shell:
!   chi(1) = 1,
!   chi(j) = chi(j-1) + Delta log m/tol_dm_max
!            + Delta L/(tol_dl_max L_max) - Delta log P/dp_scale,
! the L term only where L increases outward. dp_scale depends on the
! region: tol_dp_czbase_max within tol_czbase_fine_width in log P of the
! surface convection zone's base (czb), tol_dp_core_max deeper, and
! tol_dp_env_max above (in the convection zone). A unit step in chi is
! the largest step allowed between points.
subroutine accumulate_chi_spacing
      mass_scale = s%ctrl%chi_grid_scale(ichi_dm_max)
      luminosity_scale = point_spacing_max(2)
      chi(1) = 1.0D0
      do j = 2, s%nz
         pressure_test = s%logP(j) - s%logP(s%envelope_cz_bottom_index)
         if (abs(pressure_test).lt.s%ctrl%chi_grid_scale(ichi_czbase_fine_width)) then
! The finely zoned region around the convection-zone base.
            dp_scale = s%ctrl%chi_grid_scale(ichi_dp_czbase_max)
         else if (pressure_test.gt.s%ctrl%chi_grid_scale(ichi_czbase_fine_width)) then
! Below it, in the radiative interior.
            dp_scale = s%ctrl%chi_grid_scale(ichi_dp_core_max)
         else
! Above it, in the surface convection zone.
            dp_scale = s%ctrl%chi_grid_scale(ichi_dp_env_max)
         endif
         if (s%luminosity_lsun(j).gt.s%luminosity_lsun(j-1)) then
            dchi = (s%log_mass(j)-s%log_mass(j-1))/mass_scale + &
                 (s%luminosity_lsun(j)-s%luminosity_lsun(j-1))/luminosity_scale - &
                 (s%logP(j)-s%logP(j-1))/dp_scale
         else
            dchi = (s%log_mass(j)-s%log_mass(j-1))/mass_scale - &
                 (s%logP(j)-s%logP(j-1))/dp_scale
         endif
         chi(j) = chi(j-1)+dchi
      end do
end subroutine accumulate_chi_spacing

! ---------------------------------------------------------------
! Place the new points evenly in chi between successive flagged points,
! int(Delta chi) + 1 intervals per segment so that every step is below
! one, reading log m off a natural cubic spline of log m against chi.
! The flagged points themselves are kept (the spline passes through its
! knots). ierr past json points or on a spline failure.
subroutine place_points_along_chi
      do j = 1,s%nz
         spline_x(j) = chi(j)
         spline_y(j) = s%log_mass(j)
      end do
! The spline of log m against chi.
      call natural_spline(spline_x,spline_y,s%nz,spline_second_deriv)
      chi_prev = chi(1)
      point_prev_index = 1
      new_num_zones = 1
      do i = 1,flag_count
         delta_chi = chi(flag_point(i))-chi_prev
         segment_point_count = int(delta_chi)+1
         if (new_num_zones+segment_point_count.gt.json) then
! More points than the arrays hold: stop with ierr.
            write(*,101) new_num_zones+segment_point_count,json
 101        format(' DESIRED NUMBER OF POINTS ',I6,'EXCEEEDS JSON ', &
     I6/' RUN STOPPED')
            ierr = 1
            return
         endif
         dchi = delta_chi/dfloat(segment_point_count)
! The new points of this segment.
         do j = point_prev_index+1,point_prev_index+segment_point_count
            spline_eval_x = chi_prev + dchi
            call splint(spline_x, spline_y, s%nz, &
                 spline_second_deriv, spline_eval_x, spline_eval_y, jerr_gate)
            ! A spline failure returns through ierr (splint prints the
            ! diagnostic).
            if (jerr_gate /= 0) then
               ierr = jerr_gate
               return
            end if
            new_log_mass(j) = spline_eval_y
            chi_prev = spline_eval_x
         end do
         new_num_zones = new_num_zones + segment_point_count
         point_prev_index = new_num_zones
      end do
end subroutine place_points_along_chi

! ---------------------------------------------------------------
! Interpolate X onto the new points (x_new) and, where X rises outward
! by more than tol_dx_max between two new points, insert
! int(Delta X/tol_dx_max) points evenly in log m between them. Skipped
! when tol_dx_max is zero; ierr on a spline failure.
subroutine add_points_for_x_jumps
      do j = 1,s%nz
         spline_x(j) = s%log_mass(j)
         spline_y(j) = s%xa(i_h1,j)
      end do
! The spline of X against log m, evaluated at the new points.
      call natural_spline(spline_x,spline_y,s%nz,spline_second_deriv)
      do i = 2,new_num_zones
         spline_eval_x = new_log_mass(i)
         call splint(spline_x, spline_y, s%nz, spline_second_deriv, &
              spline_eval_x, spline_eval_y, jerr_gate)
         if (jerr_gate /= 0) then
            ierr = jerr_gate
            return
         end if
         x_new(i) = spline_eval_y
      end do
! Skip when tol_dx_max (legacy HPMAX(3)) is zero.
      if (point_spacing_max(3).gt.1.0D-15) then
! The intervals where X rises outward by more than tol_dx_max, collected
! outermost first, so that each insertion below shifts only points
! outside the intervals still to be filled.
      gradient_flag_count = 0
      do j = new_num_zones,2,-1
         if (x_new(j)-x_new(j-1).gt.point_spacing_max(3)) then
            gradient_flag_count = gradient_flag_count + 1
            gradient_flag_index(gradient_flag_count) = j
         endif
      end do
      if (gradient_flag_count.ne.0) then
      working_num_zones = new_num_zones
      do i = 1,gradient_flag_count
         j = gradient_flag_index(i)
         delta_x_over_max = (x_new(j) - x_new(j-1))/point_spacing_max(3)
! The number of points to insert. Shift the points from j outward up by
! num_new_points, then fill the interval below j evenly in log m.
         num_new_points = int(delta_x_over_max)
         point_insert_spacing = (new_log_mass(j)-new_log_mass(j-1))/ &
              dfloat(num_new_points+1)
! More points than the arrays hold: stop with ierr, as
! place_points_along_chi does (no check here before 2026-10).
         if (working_num_zones+num_new_points.gt.json) then
            write(*,'(a,i6,a,i6)') ' DESIRED NUMBER OF POINTS ', &
                 working_num_zones+num_new_points, ' EXCEEDS JSON ', json
            ierr = 1
            return
         end if
         do k = working_num_zones+num_new_points,j+num_new_points,-1
            new_log_mass(k) = new_log_mass(k-num_new_points)
         end do
         do k = j + num_new_points -1, j -1, -1
            new_log_mass(k) = new_log_mass(k+1) - point_insert_spacing
         end do
         working_num_zones = working_num_zones + num_new_points
      end do
      new_num_zones = working_num_zones
      end if
      end if
end subroutine add_points_for_x_jumps

! ---------------------------------------------------------------
! Interpolate Z onto the new points (z_new) and, where Z falls outward
! by more than tol_dz_max between two new points, insert
! int(Delta Z/tol_dz_max) points evenly in log m between them, as
! add_points_for_x_jumps does for X. Skipped when tol_dz_max is zero;
! ierr on a spline failure.
subroutine add_points_for_z_jumps
      do j = 1,s%nz
         spline_x(j) = s%log_mass(j)
         spline_y(j) = s%xa(i_metals,j)
      end do
! The spline of Z against log m, evaluated at the new points.
      call natural_spline(spline_x,spline_y,s%nz,spline_second_deriv)
      do i = 2,new_num_zones
         spline_eval_x = new_log_mass(i)
         call splint(spline_x, spline_y, s%nz, spline_second_deriv, &
              spline_eval_x, spline_eval_y, jerr_gate)
         if (jerr_gate /= 0) then
            ierr = jerr_gate
            return
         end if
         z_new(i) = spline_eval_y
      end do
! Skip when tol_dz_max (legacy HPMAX(4)) is zero; otherwise collect the
! intervals where Z falls outward by more than tol_dz_max, outermost
! first.
      if (point_spacing_max(4).gt.1.0D-15) then
      gradient_flag_count = 0
      do j = new_num_zones,2,-1
         if (z_new(j-1)-z_new(j).gt.point_spacing_max(4)) then
            gradient_flag_count = gradient_flag_count + 1
            gradient_flag_index(gradient_flag_count) = j
         endif
      end do
      if (gradient_flag_count.ne.0) then
      working_num_zones = new_num_zones
      do i = 1,gradient_flag_count
         j = gradient_flag_index(i)
         delta_z_over_max = (z_new(j-1) - z_new(j))/point_spacing_max(4)
! The number of points to insert; shift and fill as for X.
         num_new_points = int(delta_z_over_max)
         point_insert_spacing = (new_log_mass(j)-new_log_mass(j-1))/ &
              dfloat(num_new_points+1)
! More points than the arrays hold: stop with ierr, as
! place_points_along_chi does (no check here before 2026-10).
         if (working_num_zones+num_new_points.gt.json) then
            write(*,'(a,i6,a,i6)') ' DESIRED NUMBER OF POINTS ', &
                 working_num_zones+num_new_points, ' EXCEEDS JSON ', json
            ierr = 1
            return
         end if
         do k = working_num_zones+num_new_points,j+num_new_points,-1
            new_log_mass(k) = new_log_mass(k-num_new_points)
         end do
         do k = j + num_new_points -1, j -1, -1
            new_log_mass(k) = new_log_mass(k+1) - point_insert_spacing
         end do
         working_num_zones = working_num_zones + num_new_points
      end do
      new_num_zones = working_num_zones
      end if
      end if
end subroutine add_points_for_z_jumps

! ---------------------------------------------------------------
! Drop the new points that lie within tol_dm_min in log m of the last
! point kept; the first and last points are always kept. The mesh is
! compacted through compact_log_mass, then copied back.
subroutine drop_crowded_points
      j = 1
      compact_log_mass(j) = new_log_mass(j)
      do k = 2,new_num_zones-1
       if (new_log_mass(k) - compact_log_mass(j).gt.s%ctrl%chi_grid_scale(ichi_dm_min)) then
          j = j + 1
          compact_log_mass(j) = new_log_mass(k)
       endif
      end do
      j = j + 1
      compact_log_mass(j) = new_log_mass(new_num_zones)
      new_num_zones = j
      do j = 2,new_num_zones
       new_log_mass(j) = compact_log_mass(j)
      end do
end subroutine drop_crowded_points

! ---------------------------------------------------------------
! Locate the edges of the convection zones on the new mesh: the
! convective core's outer edge becomes the last new point at or inside
! the old edge, and the surface convection zone's base the first new
! point at or outside the old base. Then set convective_flag on the new
! points: convective in the core and the envelope zone, radiative
! between.
subroutine locate_new_cz_edges
      if (s%core_cz_top_index.gt.1) then
       do j = 2,new_num_zones
          if (new_log_mass(j).gt.s%log_mass(s%core_cz_top_index)) exit
       end do
         s%core_cz_top_index = j - 1
      else
       s%core_cz_top_index = 1
      endif
      if (s%envelope_cz_bottom_index.lt.s%nz) then
       do j = new_num_zones-1,1,-1
          if (new_log_mass(j).lt.s%log_mass(s%envelope_cz_bottom_index)) &
               exit
       end do
         s%envelope_cz_bottom_index = j + 1
      else
       s%envelope_cz_bottom_index = new_num_zones
      endif
      if (s%core_cz_top_index.gt.1) then
       do i = 1,s%core_cz_top_index
          s%convective_flag(i) = .true.
       end do
       radiative_zone_begin = s%core_cz_top_index + 1
      else
       radiative_zone_begin = 1
      endif
      if (s%envelope_cz_bottom_index.lt.new_num_zones) then
       do i = s%envelope_cz_bottom_index,new_num_zones
          s%convective_flag(i) = .true.
       end do
       radiative_zone_end = s%envelope_cz_bottom_index - 1
      else
       radiative_zone_end = new_num_zones
      endif
      if (radiative_zone_end.ge.1 .and. radiative_zone_begin.le.new_num_zones) then
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
subroutine interpolate_onto_new_grid
      use math_lib
      call regrid_composition_terms
      call regrid_structure_variables
      call regrid_rotation_variables
      call regrid_start_of_step_eps
      call transfer_new_points
      call set_new_shell_masses
      if (s%job%rotation_active) then
       call respin_on_new_grid
       if (ierr /= 0) return
       call retally_angular_momentum
      endif
      call rerun_physics_on_new_grid
      if (ierr /= 0) return
      call refresh_surface_opacity
end subroutine interpolate_onto_new_grid

! ---------------------------------------------------------------
! Set the old and new point counts, and spline the entropy terms, every
! species (now and at the start of the step) and, when the surface
! deuterium is above 1e-14, the start-of-step deuterium burning rate
! onto the new points.
subroutine regrid_composition_terms
! regridded_xa: the species tracked at the new points, before they
! replace star%xa and star%xa_start.
      double precision, allocatable :: regridded_xa(:,:)
! Each field is splined from the old points (star%log_mass,
! old_point_count of them) to the new ones (new_log_mass,
! new_point_count); regrid_in_place overwrites the field it is given.
      old_point_count = s%nz
      new_point_count = new_num_zones
! The entropy terms (7/91).
      call regrid_in_place(s%dlog_temperature_dt)
      call regrid_in_place(s%dlog_pressure_dt)
      call regrid_in_place(s%dln_luminosity_dt)
      call regrid_in_place(s%dlog_radius_dt)
      allocate(regridded_xa(num_species_tracked,new_point_count))
      call regrid_composition(s%log_mass(1:old_point_count), &
           s%xa(1:num_species_tracked,1:old_point_count), &
           new_log_mass(1:new_point_count), regridded_xa)
      s%xa(1:num_species_tracked,1:new_point_count) = regridded_xa
! xa_start (legacy HCOMPP), the composition at the start of the step, is
! needed for composition diffusion in rotating models.
      call regrid_composition(s%log_mass(1:old_point_count), &
           s%xa_start(1:num_species_tracked,1:old_point_count), &
           new_log_mass(1:new_point_count), regridded_xa)
      s%xa_start(1:num_species_tracked,1:new_point_count) = regridded_xa
! The start-of-step deuterium burning rate, when the surface deuterium
! is above 1e-14 (MHP 05/02). (star%nz is still the old point count
! here: when the mesh grew, this is the regridded deuterium at that
! index, below the new surface.)
      if (s%job%use_extended_composition .and. &
           s%xa(i_h2,s%nz).ge.1.0D-14) then
         call regrid_in_place(s%deuterium_burning_rate_start)
      endif
end subroutine regrid_composition_terms

! ---------------------------------------------------------------
! Spline log P, log R, L, log T and log rho onto the new points, into the
! start-of-step arrays star%*_start. Once transfer_new_points has copied
! them into star%logP etc., the *_start arrays hold the start-of-step
! model on the new mesh (what rezone_or_snapshot stores explicitly when
! it does not rezone), which timestep_limit_structure,
! mid_timestep_model, burn_step and convection_geometry read.
subroutine regrid_structure_variables
      type(structure_profiles) :: regridded
      call regrid_structure(s%log_mass(1:old_point_count), &
           structure_profiles(log_pressure=s%logP(1:old_point_count), &
                log_radius=s%logR(1:old_point_count), &
                luminosity_lsun=s%luminosity_lsun(1:old_point_count), &
                log_temperature=s%logT(1:old_point_count), &
                log_density=s%logRho(1:old_point_count)), &
           new_log_mass(1:new_point_count), regridded)
      s%logP_start(1:new_point_count) = regridded%log_pressure
      s%logR_start(1:new_point_count) = regridded%log_radius
      s%luminosity_lsun_start(1:new_point_count) = regridded%luminosity_lsun
      s%logT_start(1:new_point_count) = regridded%log_temperature
      s%logRho_start(1:new_point_count) = regridded%log_density
end subroutine regrid_structure_variables

! ---------------------------------------------------------------
! For rotating models, spline omega, the specific angular momentum j,
! the shape factors fp and ft, r0 and eta2 onto the new points, into the
! old_* arrays (and fp_old, ft_old) that transfer_new_points adopts.
subroutine regrid_rotation_variables
      type(rotation_profiles) :: regridded
      if (s%job%rotation_active) then
         call regrid_rotation(s%log_mass(1:old_point_count), &
              rotation_profiles(omega=s%omega(1:old_point_count), &
                   specific_angular_momentum=s%j_rot(1:old_point_count), &
                   fp=s%fp_rot(1:old_point_count), &
                   ft=s%ft_rot(1:old_point_count), &
                   mean_radius=s%mean_radius(1:old_point_count), &
                   eta_squared=s%eta_squared(1:old_point_count)), &
              new_log_mass(1:new_point_count), regridded)
         s%old_omega(1:new_point_count) = regridded%omega
         s%old_specific_angular_momentum(1:new_point_count) = &
              regridded%specific_angular_momentum
         fp_old(1:new_point_count) = regridded%fp
         ft_old(1:new_point_count) = regridded%ft
         s%old_mean_radius(1:new_point_count) = regridded%mean_radius
         s%old_eta_squared(1:new_point_count) = regridded%eta_squared
      endif

end subroutine regrid_rotation_variables

! ---------------------------------------------------------------
! For rotating models (or the extended composition with envelope
! overshoot), spline the start-of-step energy generation onto the new
! points (MHP 6/00): eps_total into old_eps_total and
! eps_total + eps_neu + eps_grav into old_eps_plus_neu_grav, for the
! mid-step Eddington-Sweet velocities. This runs before
! transfer_new_points, while star%log_mass is still the old mesh. (2026
! fix, bugsweep sec-11: hpoint.f did this after the transfer, against
! the new mesh, so the values were misplaced on every rotating rezone.)
subroutine regrid_start_of_step_eps
      if (s%job%rotation_active .or. (s%job%use_extended_composition .and. &
           s%job%envelope_overshoot_active)) then
         call regrid_profile(s%log_mass(1:old_point_count), &
              s%eps_total(1:old_point_count), new_log_mass(1:new_point_count), &
              s%old_eps_total(1:new_point_count))
         do zone_index = 1,s%nz
            spline_y(zone_index) = s%eps_total(zone_index)+s%eps_channels(i_eps_neu,zone_index)+ &
                 s%eps_channels(i_eps_grav,zone_index)
         end do
         call regrid_profile(s%log_mass(1:old_point_count), &
              spline_y(1:old_point_count), new_log_mass(1:new_point_count), &
              s%old_eps_plus_neu_grav(1:new_point_count))
      endif
end subroutine regrid_start_of_step_eps

! ---------------------------------------------------------------
! Adopt the new points: the new mesh and the splined structure and
! rotation variables become the model, with the new point count.
subroutine transfer_new_points
      do j = 1,new_num_zones
       s%log_mass(j) = new_log_mass(j)
       s%logP(j) = s%logP_start(j)
       s%logT(j) = s%logT_start(j)
       s%logR(j) = s%logR_start(j)
       s%luminosity_lsun(j) = s%luminosity_lsun_start(j)
       s%logRho(j) = s%logRho_start(j)
      end do
      if (s%job%rotation_active) then
       do j = 1, new_num_zones
          s%j_rot(j) = s%old_specific_angular_momentum(j)
          s%omega(j) = s%old_omega(j)
          s%fp_rot(j) = fp_old(j)
          s%ft_rot(j) = ft_old(j)
          s%eta_squared(j) = s%old_eta_squared(j)
          s%mean_radius(j) = s%old_mean_radius(j)
       end do
      endif
      if (solver_diagnostics()) then
         write(run_log_unit,1020) s%nz,new_num_zones
      end if
 1020 format(' POINTS  OLD',I5,'   NEW',I5)
      s%nz = new_num_zones
end subroutine transfer_new_points

! ---------------------------------------------------------------
! Set the point masses m_i = 10^(log m_i) and the mass weights
!   dm_i  = (m_(i+1) - m_(i-1))/2,  with m_0 = -m_1,
!   dm_nz = M - (m_(nz-1) + m_nz)/2,
! of the new mesh (compute_shell_masses), so dm_i is the mass between
! the midpoints on either side of point i (from the centre for i = 1),
! and dm_nz includes the envelope above the outer point
! (M = 10^log_total_mass).
subroutine set_new_shell_masses
      call compute_shell_masses(s%log_mass(1:s%nz), s%log_total_mass, &
           s%m(1:s%nz), s%dm(1:s%nz))
end subroutine set_new_shell_masses

! ---------------------------------------------------------------
! Rotating models: on the new mesh, rebuild the moments of inertia, the
! convection regions for angular momentum transport, omega from the
! splined specific angular momentum j, and the shape factors; ierr from
! the regions or the shape factors.
subroutine respin_on_new_grid
      use math_lib
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
            convective_zone_bounds,num_radiative_zones,num_convective_zones, s%iov1, s%iov2, ierr)
       if (ierr /= 0) return
! The moments of inertia of the distorted shells (with eta2 and r0)
! before omega is re-derived, as the 2015/04/06 version did (JNT
! 2025/09/03). It may be redundant, since the convection-zone rotation
! law (walpcz) runs in this version, but it is harmless.
       call zone_moments_of_inertia(s%eta_squared,s%logR,s%log_mass,s%dm,1,s%nz, &
            s%omega,s%mean_radius,s%i_rot,s%qiw)

! Omega from the specific angular momentum j (legacy HJM).
       call omega_from_j(s%logRho,s%j_rot,s%logR, &
            s%log_mass,s%dm,s%am_transport_convective_flag,s%nz, &
            s%eta_squared,s%i_rot,s%omega,s%qiw,s%mean_radius,s%disk_gate_age_gyr)
! The shape factors fp and ft, r0 and eta2 for this omega.
       call rotation_shape_factors(s%logRho,s%logR,s%log_mass,s%nz,s%omega, &
            s%eta_squared,s%fp_rot,s%ft_rot,s%mean_gravity,s%mean_radius,ierr)
       if (ierr /= 0) return
! The moments of inertia for the new shape.
       call zone_moments_of_inertia(s%eta_squared,s%logR,s%log_mass,s%dm,1,s%nz, &
            s%omega,s%mean_radius,s%i_rot,s%qiw)
end subroutine respin_on_new_grid

! ---------------------------------------------------------------
! Rotating models: re-total the angular momentum and the rotational
! kinetic energy on the new mesh (compute_angular_momentum_totals),
!   J = sum_i j_i dm_i,   E_rot = sum_i (1/2) omega_i j_i dm_i,
! and store the rotation state as the start of the step.
subroutine retally_angular_momentum
       call compute_angular_momentum_totals(s%j_rot(1:s%nz), &
            s%omega(1:s%nz), s%dm(1:s%nz), &
            s%kinetic_energy_rot(1:s%nz), sum_angular_momentum, &
            sum_rotational_ke)
! The totals before and after the rezone, a solver diagnostic.
       if (solver_diagnostics()) &
            write(run_log_unit,1120)total_angular_momentum, &
            sum_angular_momentum,total_rotational_ke,sum_rotational_ke
 1120    format(1X,'TOTAL J OF STAR - PREVIOUS ',1PE21.13,' NEW ', &
     1PE21.13/' TOTAL ROTATIONAL K.E. OF STAR-PREVIOUS ',1PE21.13, &
     ' NEW ',1PE21.13)
       total_angular_momentum = sum_angular_momentum
       total_rotational_ke = sum_rotational_ke
! The start-of-step rotation state, for the angular momentum transport.
       do i = 1,s%nz
          s%old_omega(i) = s%omega(i)
          s%old_hg(i) = s%mean_gravity(i)
          s%old_eta_squared(i) = s%eta_squared(i)
          s%old_mean_radius(i) = s%mean_radius(i)
          s%convective_flag_start(i) = s%convective_flag(i)
! j also goes to old_specific_angular_momentum (legacy HJX), for I/O
! as well (MHP 10/91).
            s%old_specific_angular_momentum(i) = s%j_rot(i)
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
subroutine rerun_physics_on_new_grid
      if (s%job%rotation_active .or. (s%job%use_extended_composition .and. &
           s%job%envelope_overshoot_active)) then
       call compute_shell_physics(s, s%fp_rot,s%ft_rot,s%xa,s%logRho,s%luminosity_lsun, &
            s%logP,s%logR,s%log_mass,s%logT, &
            s%convective_flag,s%nz,s%log_Teff, jerr)
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
! The quantities for the velocities at the start of the step on the new
! mesh, so that several short diffusion sub-steps can be taken within
! one evolutionary step.
         do zone_index = 1,s%nz
            s%old_gradr(zone_index) = s%gradr(zone_index)
            s%old_gradT(zone_index) = s%gradT(zone_index)
            s%old_grada(zone_index) = s%grada(zone_index)
            s%old_mu(zone_index) = s%mu(zone_index)
            s%old_opacity(zone_index) = s%opacity_zone(zone_index)
            s%old_cp(zone_index) = s%cp(zone_index)
            s%old_qdt(zone_index) = s%qdt(zone_index)
            s%old_conv_vel(zone_index) = s%conv_vel(zone_index)
            s%old_visc(zone_index) = s%visc(zone_index)
            s%old_thdif(zone_index) = s%thdif(zone_index)
! The start-of-step eps sum (legacy ESUMO) is splined in
! regrid_start_of_step_eps (MHP 6/00). (Until 2026-10 this loop also
! stored grad_ad - grad (MHP 06/02) in rot_scr%del_grad_diff_interface,
! which rotation_stability_setup overwrites before anything reads it.)
         end do
! The theta term of the advection-diffusion transport on the new mesh,
!   theta_i = r^4/(G m) omega d omega/dr   (r, m, omega: means of i-1, i)
! (MHP 6/02), commented out with that transport (2026-10-12).
!         do i = 2,star%nz
!            delta_radius = exp(ln10*star%logR(i))-exp(ln10*star%logR(i-1))
!            delta_omega = star%omega(i) - star%omega(i-1)
!            omega_mid = 0.5D0*(star%omega(i)+star%omega(i-1))
!            log_factor = 2.0D0*(star%logR(i)+star%logR(i-1))-0.5D0* &
!     (star%log_mass(i)+star%log_mass(i-1))-cgl
!            rot_scr%theta_rezoned(i) = exp(ln10*log_factor)*omega_mid*delta_omega/delta_radius
!         end do
      endif
end subroutine rerun_physics_on_new_grid

! ---------------------------------------------------------------
! Rebuild the surface opacity tables (for the new surface X and Z) when
! the surface X has moved by more than 1e-8 since they were built.
subroutine refresh_surface_opacity
      if (dabs(s%xnew-s%xa(i_h1,s%nz)).gt.1.0D-8) then
               s%xnew = s%xa(i_h1,s%nz)
               s%znew = s%xa(i_metals,s%nz)
               call kap_update_surface_tables(s%xnew, ierr=ierr)
               if (ierr /= 0) return

      end if
end subroutine refresh_surface_opacity

!----------------------------------------------------------------------
! regrid_in_place
!----------------------------------------------------------------------
! Spline one field from the old mesh (star%log_mass, old_point_count
! points) onto the new one (new_log_mass, new_point_count points) with
! regrid_profile, overwriting field(1:new_point_count). regrid_out holds
! the result until the old values have been read (the legacy HTO).
subroutine regrid_in_place(field)
      double precision, intent(inout) :: field(:)

      call regrid_profile(s%log_mass(1:old_point_count), field(1:old_point_count), &
           new_log_mass(1:new_point_count), regrid_out(1:new_point_count))
      field(1:new_point_count) = regrid_out(1:new_point_count)
end subroutine regrid_in_place

end subroutine rezone

end module rezoning
