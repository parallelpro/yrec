!----------------------------------------------------------------------
! turnover_timescale
!----------------------------------------------------------------------
! The convective turnover timescale near the base of the surface
! convection zone, from the stitched model: compute_turnover_timescale
! gathers the points, and turnover_from_interior_new finds the zone and
! times it.
module turnover_timescale
      use star_info_lib, only: star_info
      implicit none
      private
      public :: compute_turnover_timescale

contains

!----------------------------------------------------------------------
! compute_turnover_timescale (formerly gettau; G. Somers 3/17)
!----------------------------------------------------------------------
! The convective turnover timescale star%convective_turnover_timescale
! (s), one pressure scale height above the base of the surface
! convection zone, and radius_at_bcz, the radius of that base (Rsun; 0
! without a surface convection zone). It reads the stitched model
! (core/diagnostics/stitched_model.f90), the converged interior followed
! by the fixed-step envelope integration, built once per step in
! evolve_step before compute_observables, and walks it as one structure
! with turnover_from_interior_new. It integrates nothing itself, and
! star%pphot belongs to the stitched build (its atm_get). The mid-step
! consumers (the wind saturation, the deuterium limiter) read the values
! stored at the end of the previous step (the _old and pphot0 copies).
! Called by read_starting_model and compute_observables.
subroutine compute_turnover_timescale(s, radius_at_bcz, ierr)
      use stitched_model, only: &
           ip_mass, ip_logR, ip_logRho, ip_logP, ip_conv, ip_gradr, &
           ip_grada, ip_conv_vel
      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      double precision, intent(out) :: radius_at_bcz
      integer, intent(out) :: ierr

      double precision :: combined_radius(json), combined_pressure(json), &
           combined_density(json), combined_mass(json), &
           combined_gravity(json), combined_velocity(json)
      double precision :: combined_grad1(json), combined_grad2(json)
      logical :: combined_convective_flag(json)
      integer :: j, num_points

      ierr = 0
! Start from zero; pphot is left alone (the stitched build owns it).
      s%convective_turnover_timescale = 0.0

! The span. In the interior model, zones j+1..nz-1 are the convective
! ones below the surface (j is the outermost radiative zone below nz).
! When the base of that zone lies less than one dex in pressure below
! the interior's surface, log P(j+1) - log P(nz) < 1 (or the outer
! zones are radiative), the walk takes the stitched interior and
! envelope (regions 1-2 of the stitched grid, never the atmosphere),
! clamped to the json-sized arrays the walker takes. A deeper zone, or a
! fully convective star, is walked on the nz interior points alone.
      do j = s%nz-1, 1, -1
         if (.not. s%convective_flag(j)) exit
      end do
      if (s%logP(j+1) - s%logP(s%nz) .lt. 1.0d0) then
         num_points = min(s%stitched%n_ie, json)
      else
         num_points = s%nz
      end if
! The walker's inputs: m (g), log r, log P, log rho, g = G m/r^2, the
! convective velocity, gradr, grada and the convective flag.
      do j = 1, num_points
         combined_mass(j) = s%stitched%prof(ip_mass,j)*s%solar_mass_cgs
         combined_radius(j) = s%stitched%prof(ip_logR,j)
         combined_pressure(j) = s%stitched%prof(ip_logP,j)
         combined_density(j) = s%stitched%prof(ip_logRho,j)
         combined_gravity(j) = combined_mass(j)* &
              exp(ln10*(cgl - 2.0d0*combined_radius(j)))
         combined_velocity(j) = s%stitched%prof(ip_conv_vel,j)
         combined_grad1(j) = s%stitched%prof(ip_gradr,j)
         combined_grad2(j) = s%stitched%prof(ip_grada,j)
         combined_convective_flag(j) = s%stitched%prof(ip_conv,j) .ne. 0.0d0
      end do

      call turnover_from_interior_new(s, combined_mass, combined_convective_flag, &
           combined_radius, combined_pressure, combined_density, &
           combined_gravity, num_points, s%nz, combined_velocity, &
           combined_grad1, combined_grad2, radius_at_bcz, ierr)
      if (ierr /= 0) return

      return
end subroutine compute_turnover_timescale

!----------------------------------------------------------------------
! turnover_from_interior_new (legacy tauintnew.f; G. Somers 3/17)
!----------------------------------------------------------------------
! The base of the surface convection zone and its turnover timescale,
! on the num_points points passed in (the first num_interior_points of
! them the interior). The method: search from the center outward for the
! base of the surface convection zone, then from the surface inward. If
! both find the same base, the zone is unique and the timescale is the
! top-down one, H_p/v_c one pressure scale height above the base. If
! they differ, several convection cells lie in the envelope, and the
! timescale is the crossing time of the lowest one, its width over v_c
! at its middle (also the fallback when one H_p above the base lies
! outside the star).
! Outputs: star%convective_turnover_timescale (s), radius_at_bcz (Rsun;
! 0 without a surface convection zone) and ierr.
! Legacy names: cz_base_index (IMAX), the base found from below;
! upper_cz_base_index (JMAX), the base of the uppermost zone;
! radius_at_bcz (RBCZ); the timescale (TAUCZ).
subroutine turnover_from_interior_new(s, shell_mass, convective_flag, log10_radius, &
     log10_pressure, log10_density, local_gravity, num_points, &
     num_interior_points, convective_velocity, radiative_gradient, &
     adiabatic_gradient, radius_at_bcz, ierr)
      use math_lib
      use numerics_lib
      use convection_geometry, only: cz_base_fraction
      use physical_constants, only: ln10
      use star_info_lib, only: json
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: shell_mass(json)
      logical, intent(in) :: convective_flag(json)
      double precision, intent(inout) :: log10_radius(json), &
           log10_pressure(json), log10_density(json), local_gravity(json)
      integer, intent(in) :: num_points, num_interior_points
      double precision, intent(in) :: convective_velocity(json)
      double precision, intent(in) :: radiative_gradient(json), &
           adiabatic_gradient(json)
      double precision, intent(out) :: radius_at_bcz
      integer, intent(inout) :: ierr
! --- locals ---
! The nodes and values of the four-point cubic splines.
      double precision :: spline_x_delta(4), spline_x_radius(4), &
           spline_y_velocity(4), spline_y_pscale(4), &
           spline_deriv(4)
      logical :: fully_convective_flag, surface_cz_deep_enough
      integer :: i, j, k, kk, search_start_index
      logical :: in_radiative_zone
      integer :: core_cz_top_index, cz_base_index, upper_cz_base_index, &
           index_gap, cz_top_index
      double precision :: mass_at_cz_top, mass_at_cz_bottom
      double precision :: dd2, dd1, interp_fraction
      double precision :: pressure_scale_height2, radius_test2, val_test
      double precision :: pressure_scale_height1, radius_test1
      double precision :: pressure_scale_height, convective_velocity_bcz
      double precision :: radius_top_cz, cz_width, radius_test
      logical :: spline_taucz_done

! The extent of the surface convection zone. The search runs outward
! from the center, past a convective core, to the base of the next
! convection zone; a rejected candidate restarts it above that base.
! (The search loop, the if/else chain and spline_taucz_done replace the
! legacy jumps to labels 50, 70/90/100 and 140.)
      fully_convective_flag = .false.
      search_start_index = 1
! No convective core found yet (the scan sets it only when the centre is
! convective; the guard below reads it).
      core_cz_top_index = 0
      search: do
      call scan_for_cz_base_from_below()
! A scan that runs off the top finds either a fully convective star or
! a radiative envelope (perhaps over a convective core).
      if (i .gt. num_points .and. in_radiative_zone) then ! radiative envelope
         fully_convective_flag = .false.
         surface_cz_deep_enough = .false.
         cz_base_index = num_points
         upper_cz_base_index = num_points
         core_cz_top_index = num_points
      else
      call find_cz_bases()
! Trap the pathologies of locating the right convection zone.
!
! The zone must lie well above a convective core: on a stitched walk, a
! few convective shells just above the core (less than one dex in
! pressure from the core's top) are skipped. The guard on positive
! indices avoids a runtime error (KC 2025-05-30).
      if (cz_base_index .gt. 0 .and. core_cz_top_index .gt. 0) then
         if (num_points.ne.num_interior_points .and. &
              log10_pressure(core_cz_top_index)-log10_pressure(cz_base_index).lt.1) then
            search_start_index = cz_base_index+1
            cycle search
         endif
      endif
! The zone must be at least three shells thick; a thinner one is likely
! a single convective shell deep in the interior, and the search
! restarts above it.
      if (.not.(convective_flag(cz_base_index).and.convective_flag(cz_base_index+1) &
           .and.convective_flag(cz_base_index+2))) then
         if (.not.convective_flag(cz_base_index+1)) then
            search_start_index = cz_base_index+1
         else
            search_start_index = cz_base_index+2
         endif
         cycle search
      endif
      call assess_surface_cz_depth()
      end if

! The timescale, top down, for a fully or partly convective star; zero
! for a radiative envelope.
      if (surface_cz_deep_enough) then
         call locate_cz_base_radius()
!        A unique zone (index_gap = 0) is timed one H_p above its base.
!        Several cells in a mostly radiative envelope, or no spline
!        result, take the bottom cell's crossing time.
         spline_taucz_done = .false.
         if (index_gap.eq.0) then
            call time_scale_height_above_base()
            if (ierr /= 0) return
         end if
         if (.not. spline_taucz_done) then
            call time_bottom_cell_crossing()
            if (ierr /= 0) return
         end if
!        The base radius in Rsun.
         radius_at_bcz = radius_at_bcz/s%solar_radius_cgs
      else
         s%convective_turnover_timescale = 0.0d0
         radius_at_bcz = 0.0
      endif
!     A timescale above 1e20 s means the zone found lies deep in the
!     interior: search again above its base.
      if (s%convective_turnover_timescale.gt.1.0d20) then
         search_start_index = cz_base_index + 1
         cycle search
      endif
      exit search
      end do search
      return

contains

! ---------------------------------------------------------------
! Scans outward from search_start_index, past a convective core (its top
! is core_cz_top_index), to the base i of the next convection zone;
! i > num_points when there is none. core_cz_top_index is set only when
! the scan starts in a convective region.
subroutine scan_for_cz_base_from_below()
      in_radiative_zone = .true.
      if (convective_flag(search_start_index)) in_radiative_zone = .false.
      do i = search_start_index,num_points,1
         if (.not.in_radiative_zone.and..not.convective_flag(i)) then
            core_cz_top_index = i-1
            in_radiative_zone = .true.
         endif
         if (in_radiative_zone.and.convective_flag(i)) exit
      end do
end subroutine scan_for_cz_base_from_below

! ---------------------------------------------------------------
! The base of the zone found from below (cz_base_index = i) and that of
! the uppermost zone (upper_cz_base_index = j+1, searched from the
! surface down), and the index gap between them. A fully convective star
! has both at the center.
subroutine find_cz_bases()
      if (i .gt. num_points) then ! fully convective
         fully_convective_flag = .true.
         surface_cz_deep_enough = .true.
         i = 1
         j = 0
         core_cz_top_index = 1
      else
! The base of the uppermost convection zone: the outermost j that is
! radiative with j+1 convective. Requiring j+1 convective skips outer
! points that are already radiative (JvS 02/2026).
      do j = num_points-1,1,-1
         if (convective_flag(j+1) .and. .not.convective_flag(j)) exit
      end do
      endif
      cz_base_index = i
      upper_cz_base_index = j + 1
      index_gap = upper_cz_base_index-cz_base_index
end subroutine find_cz_bases

! ---------------------------------------------------------------
! An interior-only walk takes the uppermost zone's base; then whether the
! surface convection zone holds any mass (surface_cz_deep_enough).
subroutine assess_surface_cz_depth()
! On an interior-only walk, take the uppermost zone's base and treat the
! zone as unique.
      if (index_gap.ne.0.and.num_points.eq.num_interior_points) then
         cz_base_index = upper_cz_base_index
         index_gap = 0
      endif
! The masses at the top of the uppermost convection zone (the outermost
! point) and at its base (midway between the base point and the one
! below).
      mass_at_cz_top = shell_mass(num_points)
      if (upper_cz_base_index.gt.1) then
         mass_at_cz_bottom = 0.5d0*(shell_mass(upper_cz_base_index)+shell_mass(upper_cz_base_index-1))
      else
         mass_at_cz_bottom = 0.0d0
      endif
! surface_cz_deep_enough (legacy LCZSUR): a surface convection zone deep
! enough for angular-momentum loss exists; the threshold is any mass.
      if ((mass_at_cz_top-mass_at_cz_bottom)/s%solar_mass_cgs.gt.0.0d0) then
         surface_cz_deep_enough = .true.
      else
         surface_cz_deep_enough= .false.
      endif
end subroutine assess_surface_cz_depth

! ---------------------------------------------------------------
! The radius at the base of the surface convection zone, radius_at_bcz
! in cm (10 cm when fully convective).
subroutine locate_cz_base_radius()
! radius_at_bcz is found as log r, then turned into cm.
      if (fully_convective_flag) then
!           The center: log r = 1, i.e. 10 cm, negligible next to the
!           radii of stars.
         radius_at_bcz = 1.
      else
!           Linear interpolation in log r to the zero of gradr - grada
!           between the last radiative point and the base,
!           f = d_below/(d_below - d) (cz_base_fraction).
         dd2 = radiative_gradient(cz_base_index-1)-adiabatic_gradient(cz_base_index-1)
         dd1 = radiative_gradient(cz_base_index)-adiabatic_gradient(cz_base_index)
         interp_fraction = cz_base_fraction(dd2, dd1)
         radius_at_bcz = log10_radius(cz_base_index-1)+interp_fraction* &
              (log10_radius(cz_base_index)-log10_radius(cz_base_index-1))
      endif
      radius_at_bcz = exp(ln10*radius_at_bcz)
end subroutine locate_cz_base_radius

! ---------------------------------------------------------------
! The turnover time one pressure scale height above the base,
! tau = H_p/v_c, with H_p = P/(rho g) and v_c splined at the point where
! r - r_bcz = H_p.
subroutine time_scale_height_above_base()
! H_p at the base point, and its distance above the base minus H_p.
      pressure_scale_height2 = exp(ln10*(log10_pressure(cz_base_index)- &
           log10_density(cz_base_index)))/local_gravity(cz_base_index)
      radius_test2 = exp(ln10*log10_radius(cz_base_index))-radius_at_bcz
      val_test = radius_test2 - pressure_scale_height2
      if (val_test.gt.0.0) then
! The innermost convective point already lies more than H_p above the
! base: take H_p/v_c there. This leaves spline_taucz_done false, so
! time_bottom_cell_crossing runs next and overwrites the result.
         k = cz_base_index
         convective_velocity_bcz = convective_velocity(cz_base_index)
         pressure_scale_height = pressure_scale_height2
         s%convective_turnover_timescale = pressure_scale_height/convective_velocity_bcz
      else
         do k = cz_base_index+1,num_points-1,1
            pressure_scale_height1 = pressure_scale_height2
            radius_test1 = radius_test2
            pressure_scale_height2 = exp(ln10*(log10_pressure(k)-log10_density(k)))/local_gravity(k)
            radius_test2 = exp(ln10*log10_radius(k))-radius_at_bcz
            val_test = radius_test2 - pressure_scale_height2
!              The first point where r - r_bcz exceeds the local H_p.
            if (val_test.gt.0.0) then
!                 Shift the stencil outward when k is the second
!                 convective point, so that it stays above the base.
               kk = k
               if (k.eq.cz_base_index+1) kk = k+1
!                 Keep the four points kk-2..kk+1 inside the model
!                 (a zone in the top shells would reach past it).
               kk = max(3, min(kk, num_points-1))
!                 v_c and H_p by natural cubic splines through four
!                 points in x = r - r_bcz - H_p, evaluated at x = 0.
               spline_x_delta(1) = exp(ln10*log10_radius(kk-2))-radius_at_bcz- &
                         exp(ln10*(log10_pressure(kk-2)-log10_density(kk-2)))/local_gravity(kk-2)
               spline_x_delta(2) = exp(ln10*log10_radius(kk-1))-radius_at_bcz- &
                         exp(ln10*(log10_pressure(kk-1)-log10_density(kk-1)))/local_gravity(kk-1)
               spline_x_delta(3) = exp(ln10*log10_radius(kk))-radius_at_bcz- &
                         exp(ln10*(log10_pressure(kk)-log10_density(kk)))/local_gravity(kk)
               spline_x_delta(4) = exp(ln10*log10_radius(kk+1))-radius_at_bcz- &
                         exp(ln10*(log10_pressure(kk+1)-log10_density(kk+1)))/local_gravity(kk+1)
               spline_y_velocity(1) = convective_velocity(kk-2)
               spline_y_velocity(2) = convective_velocity(kk-1)
               spline_y_velocity(3) = convective_velocity(kk)
               spline_y_velocity(4) = convective_velocity(kk+1)
               call natural_spline(spline_x_delta, spline_y_velocity, 4, spline_deriv)
               call splint(spline_x_delta, spline_y_velocity, 4, spline_deriv,0.0d0,convective_velocity_bcz, ierr)
               if (ierr /= 0) return
               spline_y_pscale(1) = exp(ln10*(log10_pressure(kk-2)-log10_density(kk-2)))/local_gravity(kk-2)
               spline_y_pscale(2) = exp(ln10*(log10_pressure(kk-1)-log10_density(kk-1)))/local_gravity(kk-1)
               spline_y_pscale(3) = exp(ln10*(log10_pressure(kk)-log10_density(kk)))/local_gravity(kk)
               spline_y_pscale(4) = exp(ln10*(log10_pressure(kk+1)-log10_density(kk+1)))/local_gravity(kk+1)
               call natural_spline(spline_x_delta, spline_y_pscale, 4, spline_deriv)
               call splint(spline_x_delta, spline_y_pscale, 4, spline_deriv,0.0d0,pressure_scale_height, ierr)
               if (ierr /= 0) return
! tau = H_p/v_c.
               s%convective_turnover_timescale = pressure_scale_height/convective_velocity_bcz
               spline_taucz_done = .true.
               exit
            endif
         end do
      endif
end subroutine time_scale_height_above_base

! ---------------------------------------------------------------
! Without that spline: the crossing time of the bottom convection cell,
! tau = width/v_c, with v_c splined at the cell's middle.
subroutine time_bottom_cell_crossing()
! Reached when one H_p above the base lies outside the star, when
! several convection cells sit in the radiative envelope, or when the
! innermost point already lay above one H_p. First the cell's upper
! edge: the last convective point before a radiative one.
      do k = cz_base_index,num_points-1,1
         if (.not.convective_flag(k))exit
      enddo
      cz_top_index = k-1
! The top of the cell by linear interpolation in log r to the zero of
! gradr - grada.
      dd2 = radiative_gradient(cz_top_index)-adiabatic_gradient(cz_top_index)
      dd1 = radiative_gradient(cz_top_index+1)-adiabatic_gradient(cz_top_index+1)
      interp_fraction = cz_base_fraction(dd2, dd1)
      radius_top_cz = log10_radius(cz_top_index)+interp_fraction*(log10_radius(cz_top_index+1)-log10_radius(cz_top_index))
      radius_top_cz = exp(ln10*radius_top_cz)
! The width of the cell.
      cz_width = radius_top_cz - radius_at_bcz
! The first point above the cell's middle, r_bcz + width/2.
      radius_test = log10(radius_at_bcz + 0.5*cz_width)
      do k = cz_base_index,cz_top_index+1,1
         if (log10_radius(k).gt.radius_test) exit
      enddo
      radius_test = exp(ln10*radius_test)
!           Shift the stencil outward when k is the second convective
!           point, so that it stays above the base.
      kk = k
      if (k.eq.cz_base_index+1) kk = k+1
! Keep the four points kk-2..kk+1 inside the model.
      kk = max(3, min(kk, num_points-1))
! v_c at the middle by a natural cubic spline in r through four points.
      spline_x_radius(1) = exp(ln10*log10_radius(kk-2))
      spline_x_radius(2) = exp(ln10*log10_radius(kk-1))
      spline_x_radius(3) = exp(ln10*log10_radius(kk))
      spline_x_radius(4) = exp(ln10*log10_radius(kk+1))
      spline_y_velocity(1) = convective_velocity(kk-2)
      spline_y_velocity(2) = convective_velocity(kk-1)
      spline_y_velocity(3) = convective_velocity(kk)
      spline_y_velocity(4) = convective_velocity(kk+1)
      call natural_spline(spline_x_radius, spline_y_velocity, 4, spline_deriv)
      call splint(spline_x_radius, spline_y_velocity, 4, spline_deriv,radius_test,convective_velocity_bcz, ierr)
      if (ierr /= 0) return
! tau = width/v_c.
      s%convective_turnover_timescale = cz_width/convective_velocity_bcz
end subroutine time_bottom_cell_crossing

end subroutine turnover_from_interior_new

end module turnover_timescale
