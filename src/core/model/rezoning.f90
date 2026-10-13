!----------------------------------------------------------------------
! rezoning
!----------------------------------------------------------------------
! The rezone of each step (rezone): a new mesh in log m, its points
! spaced evenly in a coordinate chi built from the changes of log m, L
! and log P, with points added at X and Z jumps; the model interpolated
! onto it and, for rotating models, the rotation re-spun on it.
module rezoning
      use shell_physics, only: compute_shell_physics
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
!      mesh; for rotating models omega re-derived from j and the total
!      angular momentum and rotational kinetic energy re-totalled.
! envelope_store_index: the envelope-triangle vertex whose integration
! stored its state where T reached tenv (0 for none), used when a point
! is added at the surface. point_reset_flag: set when outer points were
! deleted or added. h_shell_zone_begin, h_shell_active: the base of the
! H-burning shell and whether there is one. total_angular_momentum,
! total_rotational_ke: re-totalled on the new mesh (rotating models).
! ierr /= 0 on failure.
subroutine rezone(envelope_store_index, point_reset_flag, &
     h_shell_zone_begin, h_shell_active, total_angular_momentum, &
     total_rotational_ke, ierr)
      use am_regions, only: am_convective_regions
      use omega_profile, only: omega_from_j
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, i_h1, i_h2, i_metals, json
      use controls_buffer, only: ichi_czbase_fine_width, ichi_dl_max, ichi_dm_max, &
           ichi_dm_min, ichi_dp_core_max, ichi_dp_czbase_max, ichi_dp_env_max, ichi_dx_max, &
           ichi_dz_max, ichi_flag_dw, ichi_flag_dx, ichi_flag_dz
      use kap_lib
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use numerics_lib

      implicit none
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
! compact_log_mass: the buffer of drop_crowded_points. regrid_in,
! regrid_out: the table copy and the osplin output of regrid_in_place.
! (The legacy code kept these four in star%old_shell_mass, logRho_start,
! logP_start and logT_start; every element is written here before it is
! read.)
      double precision :: new_log_mass(json), compact_log_mass(json), &
           regrid_in(json), regrid_out(json)
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
      double precision :: angular_momentum_shell
! The theta term's locals, commented out with the advection-diffusion
! transport (see rerun_physics_on_new_grid).
!      double precision :: delta_radius, delta_omega, omega_mid, log_factor
      double precision :: mass_two_back, mass_prev, mass_curr
      integer :: radiative_zone_begin, radiative_zone_end
      integer :: num_radiative_zones, num_convective_zones

      ! jerr: the error code of compute_shell_physics, passed on in ierr.
      integer :: jerr

      integer, intent(out) :: ierr

      ierr = 0

      point_reset_flag = .false.
! The number of species carried: 11, or 15 with the extended composition.
      num_species_tracked = 11
      if (star%job%use_extended_composition) num_species_tracked = 15
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
      if (star%logT(star%nz).lt.star%ctrl%tenv0) then
       do i = star%nz-1,1,-1
          if (star%logT(i).gt.star%ctrl%tenv0) then
             write(run_log_unit,10) star%nz,i
   10     format(' OUTER POINTS DELETED OLD M =',I5,'  NEW M =',I5)
             star%nz = i
             star%senv = star%log_mass(star%nz) - star%log_total_mass
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
      else if (star%logT(star%nz).gt.star%ctrl%tenv1.and. &
           envelope_store_index.ne.0) then
       star%nz = star%nz + 1
       star%senv = star%stored_envelope_state(4)
       star%log_mass(star%nz) = star%log_total_mass + star%senv
       star%logP(star%nz) = star%logP(star%nz-1) + &
            (star%stored_envelope_state(1) - star%fit_point_pressure(envelope_store_index))
       star%logT(star%nz) = star%logT(star%nz-1) + &
            (star%stored_envelope_state(1) - &
            star%fit_point_pressure(envelope_store_index))*0.250D0
       star%logR(star%nz) = star%logR(star%nz-1) + &
            (star%stored_envelope_state(3) - star%fit_point_radius(envelope_store_index))
       star%luminosity_lsun(star%nz) = star%luminosity_lsun(star%nz-1)
       do i = 1,num_species_tracked
          star%xa(i,star%nz) = star%xa(i,star%nz-1)
       end do
       star%logRho(star%nz) = star%logP(star%nz) - &
            star%logT(star%nz) - 8.0D0
       j = star%nz - 1
       write(run_log_unit,60) j,star%log_mass(j),star%logP(j), &
            star%logT(j),star%logR(j),star%nz,star%log_mass(star%nz), &
            star%logP(star%nz),star%logT(star%nz), &
            star%logR(star%nz)
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
      if (star%core_cz_top_index.gt.1) then
       call append_flag_point(star%core_cz_top_index)
      endif
      if (star%envelope_cz_bottom_index.lt.star%nz .and. &
           star%envelope_cz_bottom_index.gt.1) then
       call append_flag_point(star%envelope_cz_bottom_index)
      endif
! The point just below the H-burning shell.
      if (h_shell_active.and.h_shell_zone_begin.gt.1) then
       call append_flag_point(h_shell_zone_begin - 1)
      endif
      do i = 2,star%nz
! A jump in X between neighbouring shells,
       if (dabs(star%xa(i_h1,i)-star%xa(i_h1,i-1)).gt.star%ctrl%chi_grid_scale(ichi_flag_dx)) then
          call append_flag_point(i)
! else a jump in Z,
       else if (dabs(star%xa(i_metals,i)-star%xa(i_metals,i-1)).gt. &
            star%ctrl%chi_grid_scale(ichi_flag_dz)) then
          call append_flag_point(i)
! else (rotating models) a jump in log omega.
       else if (star%job%rotation_active) then
          log_omega_top = log10(star%omega(i))
          log_omega_bot = log10(star%omega(i-1))
          if (dabs(log_omega_top-log_omega_bot).gt.star%ctrl%chi_grid_scale(ichi_flag_dw)) then
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
      if (i > (star%nz)) then
      if (.not.star%convective_flag(star%nz)) then
       overshoot_base_zone = star%nz
       fine_zone_base = star%nz
      else if (star%envelope_cz_bottom_index.eq.1) then
       overshoot_base_zone = 1
       fine_zone_base = 1
      else
! The overshoot base: the convection-zone base without envelope
! overshoot; with it, the deepest point with
! log P - log P(czb) <= overshoot_alpha_envelope (alpha used here as a
! depth in log P).
       if (.not.star%job%envelope_overshoot_active) then
          overshoot_base_zone = star%envelope_cz_bottom_index
       else
          do overshoot_base_zone = star%envelope_cz_bottom_index-1,1,-1
             if (star%logP(overshoot_base_zone)- &
                  star%logP(star%envelope_cz_bottom_index).gt.star%ctrl%overshoot_alpha_envelope) exit
          end do
            overshoot_base_zone = overshoot_base_zone + 1
            if (overshoot_base_zone.gt.1) then
               call append_flag_point(overshoot_base_zone)
            endif
       endif
       if (star%ctrl%chi_grid_scale(ichi_czbase_fine_width).eq.0.0D0) then
          fine_zone_base = overshoot_base_zone
       else
! The base of the finely zoned region: the deepest point with
! log P - log P(czb) <= tol_czbase_fine_width, flagged when it differs
! from the overshoot base.
          if (overshoot_base_zone.eq.1) then
             fine_zone_base = 1
          else
          do fine_zone_base = overshoot_base_zone-1,1,-1
             if (star%logP(fine_zone_base) - &
                  star%logP(star%envelope_cz_bottom_index).gt.star%ctrl%chi_grid_scale(ichi_czbase_fine_width)) &
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
      flag_point(flag_count) = star%nz
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
      new_log_mass(1) = star%log_mass(1)
      star%logP_start(1) = star%logP(1)
      star%luminosity_lsun_start(1) = star%luminosity_lsun(1)
      x_new(1) = star%xa(i_h1,1)
      z_new(1) = star%xa(i_metals,1)
      luminosity_max = star%luminosity_lsun(star%nz)
! The start-of-step log Teff (JVS 04/14) and surface L, saved before the
! rezone moves the points.
      star%log_Teff_start = star%log_Teff
      star%surface_luminosity_lsun_start = star%luminosity_lsun(star%nz)
! luminosity_max: the largest L in the model.
      do i = star%nz-1,1,-1
         if (star%luminosity_lsun(i).gt.luminosity_max) then
            luminosity_max = star%luminosity_lsun(i)
         endif
      end do
      point_spacing_max(1) = star%ctrl%chi_grid_scale(ichi_dp_env_max)
      point_spacing_max(2) = star%ctrl%chi_grid_scale(ichi_dl_max)*luminosity_max
      point_spacing_max(3) = star%ctrl%chi_grid_scale(ichi_dx_max)
      point_spacing_max(4) = star%ctrl%chi_grid_scale(ichi_dz_max)
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
      mass_scale = star%ctrl%chi_grid_scale(ichi_dm_max)
      luminosity_scale = point_spacing_max(2)
      chi(1) = 1.0D0
      do j = 2, star%nz
         pressure_test = star%logP(j) - star%logP(star%envelope_cz_bottom_index)
         if (abs(pressure_test).lt.star%ctrl%chi_grid_scale(ichi_czbase_fine_width)) then
! The finely zoned region around the convection-zone base.
            dp_scale = star%ctrl%chi_grid_scale(ichi_dp_czbase_max)
         else if (pressure_test.gt.star%ctrl%chi_grid_scale(ichi_czbase_fine_width)) then
! Below it, in the radiative interior.
            dp_scale = star%ctrl%chi_grid_scale(ichi_dp_core_max)
         else
! Above it, in the surface convection zone.
            dp_scale = star%ctrl%chi_grid_scale(ichi_dp_env_max)
         endif
         if (star%luminosity_lsun(j).gt.star%luminosity_lsun(j-1)) then
            dchi = (star%log_mass(j)-star%log_mass(j-1))/mass_scale + &
                 (star%luminosity_lsun(j)-star%luminosity_lsun(j-1))/luminosity_scale - &
                 (star%logP(j)-star%logP(j-1))/dp_scale
         else
            dchi = (star%log_mass(j)-star%log_mass(j-1))/mass_scale - &
                 (star%logP(j)-star%logP(j-1))/dp_scale
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
      do j = 1,star%nz
         spline_x(j) = chi(j)
         spline_y(j) = star%log_mass(j)
      end do
! The spline of log m against chi.
      call natural_spline(spline_x,spline_y,star%nz,spline_second_deriv)
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
            call splint(spline_x, spline_y, star%nz, &
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
      do j = 1,star%nz
         spline_x(j) = star%log_mass(j)
         spline_y(j) = star%xa(i_h1,j)
      end do
! The spline of X against log m, evaluated at the new points.
      call natural_spline(spline_x,spline_y,star%nz,spline_second_deriv)
      do i = 2,new_num_zones
         spline_eval_x = new_log_mass(i)
         call splint(spline_x, spline_y, star%nz, spline_second_deriv, &
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
      do j = 1,star%nz
         spline_x(j) = star%log_mass(j)
         spline_y(j) = star%xa(i_metals,j)
      end do
! The spline of Z against log m, evaluated at the new points.
      call natural_spline(spline_x,spline_y,star%nz,spline_second_deriv)
      do i = 2,new_num_zones
         spline_eval_x = new_log_mass(i)
         call splint(spline_x, spline_y, star%nz, spline_second_deriv, &
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
       if (new_log_mass(k) - compact_log_mass(j).gt.star%ctrl%chi_grid_scale(ichi_dm_min)) then
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
      if (star%core_cz_top_index.gt.1) then
       do j = 2,new_num_zones
          if (new_log_mass(j).gt.star%log_mass(star%core_cz_top_index)) exit
       end do
         star%core_cz_top_index = j - 1
      else
       star%core_cz_top_index = 1
      endif
      if (star%envelope_cz_bottom_index.lt.star%nz) then
       do j = new_num_zones-1,1,-1
          if (new_log_mass(j).lt.star%log_mass(star%envelope_cz_bottom_index)) &
               exit
       end do
         star%envelope_cz_bottom_index = j + 1
      else
       star%envelope_cz_bottom_index = new_num_zones
      endif
      if (star%core_cz_top_index.gt.1) then
       do i = 1,star%core_cz_top_index
          star%convective_flag(i) = .true.
       end do
       radiative_zone_begin = star%core_cz_top_index + 1
      else
       radiative_zone_begin = 1
      endif
      if (star%envelope_cz_bottom_index.lt.new_num_zones) then
       do i = star%envelope_cz_bottom_index,new_num_zones
          star%convective_flag(i) = .true.
       end do
       radiative_zone_end = star%envelope_cz_bottom_index - 1
      else
       radiative_zone_end = new_num_zones
      endif
      if (radiative_zone_end.ge.1 .and. radiative_zone_begin.le.new_num_zones) then
      do j = radiative_zone_begin,radiative_zone_end
       star%convective_flag(j) = .false.
      end do
      end if
end subroutine locate_new_cz_edges

! ---------------------------------------------------------------
! Move the model onto the new mesh: spline the composition, structure
! and rotation variables and the start-of-step energy generation onto it
! with the shape-preserving osculatory quadratic spline (osplin), adopt
! the mesh and set its shell masses; for rotating models re-spin and
! re-total the angular momentum; re-run the shell physics where the
! transport needs it; and rebuild the surface opacity tables if the
! surface X moved.
subroutine interpolate_onto_new_grid
      use math_lib
      use rotation_scratch_lib
      call regrid_composition_terms
      call regrid_structure_variables
      call regrid_rotation_variables
      call regrid_start_of_step_eps
      call transfer_new_points
      call set_new_shell_masses
      if (star%job%rotation_active) then
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
! osplin(xval, yval, xtab, ytab, ntab, ntot) evaluates at the ntot new
! points xval (here new_log_mass) the spline through the ntab old points
! (xtab, ytab) (here star%log_mass and a variable), into yval. Each field
! below goes through regrid_in_place, which overwrites it.
      old_point_count = star%nz
      new_point_count = new_num_zones
! The entropy terms (7/91).
      call regrid_in_place(star%temperature_entropy_term)
      call regrid_in_place(star%pressure_entropy_term)
      call regrid_in_place(star%luminosity_entropy_term)
      call regrid_in_place(star%radius_entropy_term)
      do i = 1,num_species_tracked
         call regrid_in_place(star%xa(i,:))
! xa_start (legacy HCOMPP), the composition at the start of the step, is
! needed for composition diffusion in rotating models.
         call regrid_in_place(star%xa_start(i,:))
      end do
! The start-of-step deuterium burning rate, when the surface deuterium
! is above 1e-14 (MHP 05/02).
      if (star%job%use_extended_composition .and. &
           star%xa(i_h2,star%nz).ge.1.0D-14) then
         call regrid_in_place(star%deuterium_burning_rate_start)
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

      call osplin(new_log_mass,star%logP_start,star%log_mass,star%logP, &
           old_point_count,new_point_count)
      call osplin(new_log_mass,star%logR_start,star%log_mass,star%logR, &
           old_point_count,new_point_count)
      call osplin(new_log_mass,star%luminosity_lsun_start,star%log_mass,star%luminosity_lsun, &
           old_point_count,new_point_count)
      call osplin(new_log_mass,star%logT_start,star%log_mass,star%logT, &
           old_point_count,new_point_count)
      call osplin(new_log_mass,star%logRho_start,star%log_mass,star%logRho, &
           old_point_count,new_point_count)
end subroutine regrid_structure_variables

! ---------------------------------------------------------------
! For rotating models, spline omega, the specific angular momentum j,
! the shape factors fp and ft, r0 and eta2 onto the new points, into the
! old_* arrays (and fp_old, ft_old) that transfer_new_points adopts.
subroutine regrid_rotation_variables
      if (star%job%rotation_active) then
         call osplin(new_log_mass,star%old_omega,star%log_mass,star%omega, &
              old_point_count,new_point_count)
         call osplin(new_log_mass,star%old_specific_angular_momentum,star%log_mass, &
              star%j_rot,old_point_count,new_point_count)
         call osplin(new_log_mass,fp_old,star%log_mass,star%fp_rot, &
              old_point_count,new_point_count)
         call osplin(new_log_mass,ft_old,star%log_mass,star%ft_rot, &
              old_point_count,new_point_count)
         call osplin(new_log_mass,star%old_mean_radius,star%log_mass,star%mean_radius, &
              old_point_count,new_point_count)
         call osplin(new_log_mass,star%old_eta_squared,star%log_mass,star%eta_squared, &
              old_point_count,new_point_count)
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
      use rotation_scratch_lib
      if (star%job%rotation_active .or. (star%job%use_extended_composition .and. &
           star%job%envelope_overshoot_active)) then
         call osplin(new_log_mass,rot_scr%old_eps_total,star%log_mass,star%eps_total, &
              old_point_count,new_point_count)
         do zone_index = 1,star%nz
            spline_y(zone_index) = star%eps_total(zone_index)+star%eps_channels(i_eps_neu,zone_index)+ &
                 star%eps_channels(i_eps_grav,zone_index)
         end do
         call osplin(new_log_mass,rot_scr%old_eps_plus_neu_grav,star%log_mass,spline_y, &
              old_point_count,new_point_count)
      endif
end subroutine regrid_start_of_step_eps

! ---------------------------------------------------------------
! Adopt the new points: the new mesh and the splined structure and
! rotation variables become the model, with the new point count.
subroutine transfer_new_points
      do j = 1,new_num_zones
       star%log_mass(j) = new_log_mass(j)
       star%logP(j) = star%logP_start(j)
       star%logT(j) = star%logT_start(j)
       star%logR(j) = star%logR_start(j)
       star%luminosity_lsun(j) = star%luminosity_lsun_start(j)
       star%logRho(j) = star%logRho_start(j)
      end do
      if (star%job%rotation_active) then
       do j = 1, new_num_zones
          star%j_rot(j) = star%old_specific_angular_momentum(j)
          star%omega(j) = star%old_omega(j)
          star%fp_rot(j) = fp_old(j)
          star%ft_rot(j) = ft_old(j)
          star%eta_squared(j) = star%old_eta_squared(j)
          star%mean_radius(j) = star%old_mean_radius(j)
       end do
      endif
      if (solver_diagnostics()) then
         write(run_log_unit,1020) star%nz,new_num_zones
      end if
 1020 format(' POINTS  OLD',I5,'   NEW',I5)
      star%nz = new_num_zones
end subroutine transfer_new_points

! ---------------------------------------------------------------
! Set the point masses m_i = 10^(log m_i) and the mass weights
!   dm_i  = (m_(i+1) - m_(i-1))/2,  with m_0 = -m_1,
!   dm_nz = M - (m_(nz-1) + m_nz)/2,
! so dm_i is the mass between the midpoints on either side of point i
! (from the centre for i = 1), and dm_nz includes the envelope above the
! outer point (M = 10^log_total_mass).
subroutine set_new_shell_masses
      use math_lib
      mass_curr = exp(clndp*star%log_mass(1))
      mass_prev = - mass_curr
      do i = 2,star%nz
       mass_two_back = mass_prev
       mass_prev = mass_curr
       mass_curr = exp(clndp*star%log_mass(i))
       star%m(i-1) = mass_prev
       star%dm(i-1) = 0.5D0*(mass_curr-mass_two_back)
      end do
      star%m(star%nz) = mass_curr
      star%dm(star%nz) = exp(ln10*star%log_total_mass) - &
           0.5D0*(mass_prev+mass_curr)
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
       do i=1,star%nz
          star%i_rot(i) = cc23*star%dm(i)* &
               exp(ln10*2.0D0*star%logR(i))
       end do
! The convection zones, extended by overshoot, for angular momentum
! transport (am_transport_convective_flag).
       call am_convective_regions(star%xa,star%logRho,star%logP,star%logR, &
            star%log_mass,star%logT,star%convective_flag,star%nz, &
            star%am_transport_convective_flag,radiative_zone_bounds, &
            convective_zone_bounds,num_radiative_zones,num_convective_zones, star%iov1, star%iov2, ierr)
       if (ierr /= 0) return
! The moments of inertia of the distorted shells (with eta2 and r0)
! before omega is re-derived, as the 2015/04/06 version did (JNT
! 2025/09/03). It may be redundant, since the convection-zone rotation
! law (walpcz) runs in this version, but it is harmless.
       call zone_moments_of_inertia(star%eta_squared,star%logR,star%log_mass,star%dm,1,star%nz, &
            star%omega,star%mean_radius,star%i_rot,star%qiw)

! Omega from the specific angular momentum j (legacy HJM).
       call omega_from_j(star%logRho,star%j_rot,star%logR, &
            star%log_mass,star%dm,star%am_transport_convective_flag,star%nz, &
            star%eta_squared,star%i_rot,star%omega,star%qiw,star%mean_radius)
! The shape factors fp and ft, r0 and eta2 for this omega.
       call rotation_shape_factors(star%logRho,star%logR,star%log_mass,star%nz,star%omega, &
            star%eta_squared,star%fp_rot,star%ft_rot,star%mean_gravity,star%mean_radius,ierr)
       if (ierr /= 0) return
! The moments of inertia for the new shape.
       call zone_moments_of_inertia(star%eta_squared,star%logR,star%log_mass,star%dm,1,star%nz, &
            star%omega,star%mean_radius,star%i_rot,star%qiw)
end subroutine respin_on_new_grid

! ---------------------------------------------------------------
! Rotating models: re-total the angular momentum and the rotational
! kinetic energy on the new mesh,
!   J = sum_i j_i dm_i,   E_rot = sum_i (1/2) omega_i j_i dm_i,
! and store the rotation state as the start of the step.
subroutine retally_angular_momentum
       sum_angular_momentum = 0.0D0
       sum_rotational_ke = 0.0D0
       do i = 1,star%nz
          angular_momentum_shell = star%j_rot(i)*star%dm(i)
          star%kinetic_energy_rot(i) = 0.5D0*star%omega(i)*angular_momentum_shell
          sum_angular_momentum = sum_angular_momentum + angular_momentum_shell
          sum_rotational_ke = sum_rotational_ke + star%kinetic_energy_rot(i)
       end do
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
       do i = 1,star%nz
          star%old_omega(i) = star%omega(i)
          star%old_hg(i) = star%mean_gravity(i)
          star%old_moment_of_inertia(i) = star%i_rot(i)
          star%old_eta_squared(i) = star%eta_squared(i)
          star%old_mean_radius(i) = star%mean_radius(i)
          star%convective_flag_start(i) = star%convective_flag(i)
! j also goes to old_specific_angular_momentum (legacy HJX), for I/O
! as well (MHP 10/91).
            star%old_specific_angular_momentum(i) = star%j_rot(i)
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
      use rotation_scratch_lib
      if (star%job%rotation_active .or. (star%job%use_extended_composition .and. &
           star%job%envelope_overshoot_active)) then
       call compute_shell_physics(star%fp_rot,star%ft_rot,star%xa,star%logRho,star%luminosity_lsun, &
            star%logP,star%logR,star%log_mass,star%logT, &
            star%convective_flag,star%nz,star%log_Teff, jerr)
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
! The quantities for the velocities at the start of the step on the new
! mesh, so that several short diffusion sub-steps can be taken within
! one evolutionary step.
         do zone_index = 1,star%nz
            rot_scr%old_gradr(zone_index) = star%gradr(zone_index)
            rot_scr%old_gradT(zone_index) = star%gradT(zone_index)
            rot_scr%old_grada(zone_index) = star%grada(zone_index)
            rot_scr%old_mu(zone_index) = star%mu(zone_index)
            rot_scr%old_opacity(zone_index) = star%opacity_zone(zone_index)
            rot_scr%old_cp(zone_index) = star%cp(zone_index)
            rot_scr%old_qdt(zone_index) = star%qdt(zone_index)
            rot_scr%old_conv_vel(zone_index) = star%conv_vel(zone_index)
            rot_scr%old_visc(zone_index) = star%visc(zone_index)
            rot_scr%old_thdif(zone_index) = star%thdif(zone_index)
! grad_ad - grad (MHP 06/02). The start-of-step eps sum (legacy ESUMO)
! is splined in regrid_start_of_step_eps (MHP 6/00).
            rot_scr%del_grad_diff_interface(zone_index) = &
                 rot_scr%old_grada(zone_index) - rot_scr%old_gradT(zone_index)
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
      if (dabs(star%xnew-star%xa(i_h1,star%nz)).gt.1.0D-8) then
               star%xnew = star%xa(i_h1,star%nz)
               star%znew = star%xa(i_metals,star%nz)
               call kap_update_surface_tables(star%xnew, ierr=ierr)
               if (ierr /= 0) return

      end if
end subroutine refresh_surface_opacity

!----------------------------------------------------------------------
! regrid_in_place
!----------------------------------------------------------------------
! Spline one field from the old mesh (star%log_mass, old_point_count
! points) onto the new one (new_log_mass, new_point_count points) with
! osplin, overwriting field(1:new_num_zones). regrid_in holds the copy
! of the table and regrid_out the result (the legacy HPO and HTO).
subroutine regrid_in_place(field)
      double precision, intent(inout) :: field(:)
      integer :: j

      do j = 1,star%nz
         regrid_in(j) = field(j)
      end do
      call osplin(new_log_mass,regrid_out,star%log_mass,regrid_in, &
           old_point_count,new_point_count)
      do j = 1,new_num_zones
         field(j) = regrid_out(j)
      end do
end subroutine regrid_in_place

end subroutine rezone

end module rezoning
