!----------------------------------------------------------------------
! am_grid
!----------------------------------------------------------------------
! The equally spaced chi grid on which the angular momentum transport is
! solved: the model's unstable region mapped onto it
! (am_transport_grid), and the solution's change of angular momentum
! mapped back to the model (equal_grid_to_model).
module am_grid
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: am_transport_grid, equal_grid_to_model

contains

! ---------------------------------------------------------------
! Maps the unstable region zone_begin..zone_end onto the equally spaced
! chi grid of the angular momentum solve (am_diffusion_coeffs, or
! am_advection_diffusion_coeffs under use_diffusion_advection_transport);
! composition_grid is its companion for the composition.
! equal_spaced_grid lays down the grid: rot_scr%chi at the region's
! model zones, rot_scr%eq_chi at the rot_scr%num_eq_points grid points
! (here as many as the region has zones) and rot_scr%eq_enclosed_mass
! there. Then:
!   interpolate_to_grid_centers: log r, I/(m r^2), j and omega at the
!     grid points;
!   sum_zone_contents: each grid zone's mass, moment of inertia I and
!     angular momentum J, the end zones taking in any convection zone
!     beyond them;
!   edge_coefficients: the diffusion coefficients at the grid zone
!     edges (and the advection-diffusion ones);
!   apply_jacobians: the factors that make them coefficients in chi.
!
! Inputs at the model zones: log_density, log_radius, radius (cm),
! log_mass, enclosed_mass and shell_mass (g), luminosity_lsun,
! log_pressure, moment_of_inertia, specific_angular_momentum and omega;
! am_diffusion_coeff and mixing_diffusion_coeff at the interfaces
! (element i between zones i-1 and i; 4 pi times a diffusivity, from
! diffusion_velocity_scales); am_transport_convective_flag;
! log_total_mass; num_zones.
! Outputs on the grid: eq_mass (g), eq_moment_of_inertia and
! eq_angular_momentum (the zones' I and J), eq_omega; at the grid edges
! (element i between points i-1 and i) eq_am_diffusion_coeff and
! eq_mixing_diffusion_coeff; grid_spacing, the chi spacing. A region with
! a single unstable interface (zone_end - zone_begin <= 1) is not
! mapped: single_interface_flag is set, the other outputs are not.
subroutine am_transport_grid(am_diffusion_coeff, mixing_diffusion_coeff, log_density, &
     moment_of_inertia, specific_angular_momentum, luminosity_lsun, &
     log_pressure, log_radius, radius, log_mass, enclosed_mass, shell_mass, &
     log_total_mass, zone_begin, zone_end, am_transport_convective_flag, &
     num_zones, omega, grid_spacing, eq_am_diffusion_coeff, &
     eq_mixing_diffusion_coeff, eq_moment_of_inertia, eq_angular_momentum, &
     eq_mass, eq_omega, single_interface_flag)
      use rotation_scratch_lib
      use equal_grid
      use star_info_lib, only: json
      use physical_constants
      use numerics_lib
      use math_lib
      implicit none

! Inputs
      double precision, intent(in) :: am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json), log_density(json), &
           moment_of_inertia(json), specific_angular_momentum(json), &
           luminosity_lsun(json), log_pressure(json), log_radius(json), &
           radius(json), log_mass(json), enclosed_mass(json), &
           shell_mass(json), log_total_mass
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: omega(json)
! Outputs
      double precision, intent(out) :: grid_spacing, &
           eq_am_diffusion_coeff(json), eq_mixing_diffusion_coeff(json), &
           eq_moment_of_inertia(json), eq_angular_momentum(json), &
           eq_mass(json), eq_omega(json)
      logical, intent(out) :: single_interface_flag

      double precision :: eq_reduced_moment_of_inertia(json)
      integer :: ntab, i, ii, ntabb
      double precision :: emtop, embot, scale_factor
! The spline tables: table_x/table_y at the model zones (or their
! interfaces), eval_x/eval_y on the grid (points or edges).
      double precision :: table_x(json), table_y(json), eval_x(json), eval_y(json)

! A single unstable interface: flag it and return (the caller skips
! the region).
      if (zone_end-zone_begin.le.1) then
         single_interface_flag = .true.
         return
      else
         single_interface_flag = .false.
      end if
! The grid equally spaced in chi.
      call equal_spaced_grid(luminosity_lsun,log_pressure,log_mass,zone_begin, &
           zone_end,num_zones)
      call interpolate_to_grid_centers()
      call sum_zone_contents()
      call edge_coefficients()
      call apply_jacobians()
! The solve's grid spacing is the spacing in chi.
      grid_spacing = rot_scr%eq_chi_spacing

contains

! ---------------------------------------------------------------
! The model quantities at the grid points, each splined in chi from the
! region's ntab model zones: log r into eval_y (kept for
! sum_zone_contents), I/(m r^2) into eq_reduced_moment_of_inertia, j
! into eq_angular_momentum and omega into eq_omega. sum_zone_contents
! turns I/(m r^2) and j into each grid zone's I and J.
subroutine interpolate_to_grid_centers()
! log r.
      ntab = zone_end - zone_begin + 1
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_x(i) = rot_scr%chi(i)
         table_y(i) = log_radius(ii)
      end do
      do i = 1,rot_scr%num_eq_points
         eval_x(i) = rot_scr%eq_chi(i)
      end do
      call osplin(eval_x,eval_y,table_x,table_y,ntab,rot_scr%num_eq_points)
! I/(m r^2) (2/3 for a thin spherical shell).
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_y(i) = moment_of_inertia(ii)/(shell_mass(ii)*radius(ii)**2)
      end do
      call osplin(eval_x,eq_reduced_moment_of_inertia,table_x,table_y,ntab,rot_scr%num_eq_points)
! j.
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_y(i) = specific_angular_momentum(ii)
      end do
      call osplin(eval_x,eq_angular_momentum,table_x,table_y,ntab,rot_scr%num_eq_points)
! omega.
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_y(i) = omega(ii)
      end do
      call osplin(eval_x,eq_omega,table_x,table_y,ntab,rot_scr%num_eq_points)
end subroutine interpolate_to_grid_centers

! ---------------------------------------------------------------
! Each grid zone's mass, moment of inertia and angular momentum. An
! interior point i spans the mass between the midpoints to its
! neighbours; with M = rot_scr%eq_enclosed_mass and r = 10**eval_y,
!   eq_mass(i) = (M(i+1) - M(i-1))/2,
!   eq_moment_of_inertia(i) = (I/m r^2)(i) eq_mass(i) r(i)^2,
!   eq_angular_momentum(i) = j(i) eq_mass(i).
! The end points are treated below.
subroutine sum_zone_contents()
      do i = 2, rot_scr%num_eq_points-1
         eq_mass(i) = 0.5d0*(rot_scr%eq_enclosed_mass(i+1)-rot_scr%eq_enclosed_mass(i-1))
         eq_moment_of_inertia(i) = eq_reduced_moment_of_inertia(i)* &
              eq_mass(i)*exp(ln10*2.0d0*eval_y(i))
         eq_angular_momentum(i) = eq_angular_momentum(i)*eq_mass(i)
      end do
! The end points, which can border convection zones. Inward, each spans
! the mass to the midpoint to its grid neighbour (not its model
! neighbour); outward, to the midpoint between the region's end zone and
! the next model zone beyond (or to m = 0 at the centre, the total mass
! at the surface). Each then takes in whole every zone of a convection
! zone directly beyond it: its mass, I and j dm. The end points sit on
! the end zones in chi, so their r is the model radius of those zones.
!
! Centre end.
      emtop = 0.5d0*(rot_scr%eq_enclosed_mass(2)+rot_scr%eq_enclosed_mass(1))
      if (zone_begin.gt.1) then
         embot = 0.5d0*(enclosed_mass(zone_begin)+enclosed_mass(zone_begin-1))
      else
         embot = 0.0d0
      end if
      eq_mass(1) = emtop - embot
      eq_moment_of_inertia(1) = eq_reduced_moment_of_inertia(1)*eq_mass(1)* &
           radius(zone_begin)**2
      eq_angular_momentum(1) = eq_angular_momentum(1)*eq_mass(1)
! The convection zone directly below zone_begin, if any (whether or not
! zone_begin itself is convective).
      if (zone_begin.gt.1) then
         do ii = zone_begin-1,1,-1
            if (.not.am_transport_convective_flag(ii)) exit
            eq_mass(1) = eq_mass(1)+shell_mass(ii)
            eq_moment_of_inertia(1) = eq_moment_of_inertia(1)+ &
                 moment_of_inertia(ii)
            eq_angular_momentum(1) = eq_angular_momentum(1)+ &
                 specific_angular_momentum(ii)*shell_mass(ii)
         end do
      end if
! Surface end.
      embot = 0.5d0*(rot_scr%eq_enclosed_mass(rot_scr%num_eq_points)+rot_scr%eq_enclosed_mass(rot_scr%num_eq_points-1))
      if (zone_end.lt.num_zones) then
         emtop = 0.5d0*(enclosed_mass(zone_end)+enclosed_mass(zone_end+1))
      else
         emtop = exp(ln10*log_total_mass)
      end if
      eq_mass(rot_scr%num_eq_points) = emtop - embot
      eq_moment_of_inertia(rot_scr%num_eq_points) = eq_reduced_moment_of_inertia(rot_scr%num_eq_points)* &
           eq_mass(rot_scr%num_eq_points)*radius(zone_end)**2
      eq_angular_momentum(rot_scr%num_eq_points) = eq_angular_momentum(rot_scr%num_eq_points)*eq_mass(rot_scr%num_eq_points)
! The convection zone directly above zone_end, if any.
      if (zone_end.lt.num_zones) then
         do ii = zone_end+1,num_zones
            if (.not.am_transport_convective_flag(ii)) exit
            eq_mass(rot_scr%num_eq_points) = eq_mass(rot_scr%num_eq_points)+shell_mass(ii)
            eq_moment_of_inertia(rot_scr%num_eq_points) = eq_moment_of_inertia(rot_scr%num_eq_points)+ &
                 moment_of_inertia(ii)
            eq_angular_momentum(rot_scr%num_eq_points) = eq_angular_momentum(rot_scr%num_eq_points)+ &
                 specific_angular_momentum(ii)*shell_mass(ii)
         end do
      end if
end subroutine sum_zone_contents

! ---------------------------------------------------------------
! The transport coefficients at the grid zone edges (element i between
! grid points i-1 and i), splined in chi from the model's interfaces.
! edge_grid_abscissae sets the table abscissae (table_x: the region's
! end zones and the interfaces, midway in chi between zones; ntabb
! entries) and the grid edges (eval_x). Below the first and above the
! last unstable interface a coefficient is held at its value there.
subroutine edge_coefficients()
      call edge_grid_abscissae(rot_scr%chi, rot_scr%eq_chi, rot_scr%eq_chi_spacing, ntab, &
           rot_scr%num_eq_points, table_x, eval_x, ntabb)
! The angular momentum diffusion coefficient.
      call interp_edge_coeff(am_diffusion_coeff, zone_begin, zone_end, ntab, ntabb, &
           rot_scr%num_eq_points, eval_x, table_x, table_y, eq_am_diffusion_coeff)
! The composition (mixing) diffusion coefficient.
      call interp_edge_coeff(mixing_diffusion_coeff, zone_begin, zone_end, ntab, ntabb, &
           rot_scr%num_eq_points, eval_x, table_x, table_y, eq_mixing_diffusion_coeff)
! Under use_diffusion_advection_transport, the advective and diffusive
! coefficients rot_scr%am_advective_coeff and am_diffusive_coeff: the
! Eddington-Sweet velocities rot_scr%es_advective_velocity and
! es_diffusive_velocity times 0.2 * 4 pi * difad_velocity_scale,
! splined the same way.
      if (rotation_settings%use_diffusion_advection_transport) then
         scale_factor = 0.2d0*c4pi*rotation_settings%difad_velocity_scale
         table_y(1) = scale_factor*rot_scr%es_advective_velocity(zone_begin + 1)
         do i = 2,ntab
            ii = zone_begin + i - 1
            table_y(i) = scale_factor*rot_scr%es_advective_velocity(ii)
         end do
         table_y(ntabb) = scale_factor*rot_scr%es_advective_velocity(zone_end)
         call osplin(eval_x,rot_scr%am_advective_coeff,table_x,table_y,ntabb,rot_scr%num_eq_points)
         table_y(1) = scale_factor*rot_scr%es_diffusive_velocity(zone_begin + 1)
         do i = 2,ntab
            ii = zone_begin + i - 1
            table_y(i) = scale_factor*rot_scr%es_diffusive_velocity(ii)
         end do
         table_y(ntabb) = scale_factor*rot_scr%es_diffusive_velocity(zone_end)
         call osplin(eval_x,rot_scr%am_diffusive_coeff,table_x,table_y,ntabb,rot_scr%num_eq_points)
! The composition diffusion coefficient is added to the diffusive part,
! cf. Zahn 1992, A&A 265, 115: the Eddington-Sweet term
! (es_diffusive_velocity) is the diffusion from horizontal transport,
! the mixing coefficient that from vertical transport (MHP 5/02).
         do i = 1, rot_scr%num_eq_points
            rot_scr%am_diffusive_coeff(i) = rot_scr%am_diffusive_coeff(i)+ &
                 eq_mixing_diffusion_coeff(i)
         end do
      end if
end subroutine edge_coefficients

! ---------------------------------------------------------------
! Turns the edge coefficients into coefficients in chi. Unlike a grid
! equally spaced in r, a grid in chi needs the Jacobian dchi/dr; with
! the geometry of each flux (the 4 pi is already in the coefficients,
! from diffusion_velocity_scales),
!   eq_mixing_diffusion_coeff = 4 pi D rho r^2 dchi/dr,
!   eq_am_diffusion_coeff     = 4 pi nu rho r^4 (I/m r^2) dchi/dr,
! and under use_diffusion_advection_transport the advective and
! diffusive coefficients take the angular momentum factor too.
! I/(m r^2) of the model zone stands in for the 2/3 of a thin spherical
! shell in the flux of angular momentum (MHP 5/02).
subroutine apply_jacobians()
! log10(rho r^2 dchi/dr): at the model zones in table_y (table_x reset
! to chi), splined to the grid edges in eval_y.
      call dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
           enclosed_mass, mid_scr%eps_plus_neu_grav, luminosity_lsun(num_zones), zone_begin, &
           ntab, rot_scr%num_eq_points, rot_scr%chi, eval_x, table_x, &
           table_y, eval_y)
! The composition coefficient times rho r^2 dchi/dr.
      call multiply_by_exp10(eq_mixing_diffusion_coeff, eval_y, rot_scr%num_eq_points)
! log10(rho r^4 (I/m r^2) dchi/dr) = log10(rho r^2 dchi/dr) + log10(I/m)
! at the model zones, splined to the grid edges in eval_y.
      do i = 1, ntab
         ii = zone_begin + i - 1
         table_y(i) = table_y(i)  + log10(moment_of_inertia(ii)/shell_mass(ii))
      end do
      call osplin(eval_x,eval_y,table_x,table_y,ntab,rot_scr%num_eq_points)
! The angular momentum coefficient times rho r^4 (I/m r^2) dchi/dr.
      call multiply_by_exp10(eq_am_diffusion_coeff, eval_y, rot_scr%num_eq_points)
! The advection-diffusion coefficients by the same factor.
      if (rotation_settings%use_diffusion_advection_transport) then
         call multiply_by_exp10(rot_scr%am_advective_coeff, eval_y, rot_scr%num_eq_points)
         call multiply_by_exp10(rot_scr%am_diffusive_coeff, eval_y, rot_scr%num_eq_points)
      end if
end subroutine apply_jacobians

end subroutine am_transport_grid

! ---------------------------------------------------------------
! Maps the angular momentum solve's result on the equally spaced grid
! back to the model (the legacy EQ2MOD). The fractional change of each
! grid zone's angular momentum, delta_angular_momentum/angular_momentum,
! is splined in chi to the region's model zones and applied there,
! j -> j (1 + dJ/J). If zone_begin lies in a convection zone, that
! convection zone's zones below it take zone_begin's change; likewise
! above zone_end. The changes must add up to
! total_delta_angular_momentum; if they miss it by more than 1e-5 of it,
! j over all the zones changed is scaled by one factor to restore it.
!
! Inputs on the grid (rot_scr%num_eq_points points):
!   delta_angular_momentum : the change of each grid zone's angular
!      momentum J (g cm^2/s) over the sub-step.
!   angular_momentum : each grid zone's J before it
!      (am_transport_grid's eq_angular_momentum).
!   total_delta_angular_momentum : the change of the region's total J
!      on the grid (the solve's sum of delta_angular_momentum): the
!      loss through the surface, zero without one.
! Inputs on the model:
!   shell_mass : the zone masses (g).
!   zone_begin, zone_end : the first and last zones of the region.
!   convective_flag : true where a zone is convective for angular
!      momentum transport.
!   num_points : the number of zones.
! rot_scr%chi and eq_chi (from equal_spaced_grid) place the model zones
! and the grid points in chi.
! In/out:
!   specific_angular_momentum : j at the model zones (cm^2/s).
subroutine equal_grid_to_model(delta_angular_momentum, angular_momentum, shell_mass, &
     zone_begin, zone_end, convective_flag, num_points, &
     total_delta_angular_momentum, specific_angular_momentum)
      use rotation_scratch_lib
      use star_info_lib, only: json
      use numerics_lib
      implicit none

      double precision, intent(in) :: delta_angular_momentum(json), &
           angular_momentum(json), shell_mass(json)
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_points
      double precision, intent(in) :: total_delta_angular_momentum
      double precision, intent(inout) :: specific_angular_momentum(json)
      integer :: i0, i1, i, nmod, ii
      double precision :: sumjmod, sumdjmod, test, val, ratio
! The spline tables: table_x/table_y at the grid points, eval_x/eval_y
! at the model zones.
      double precision :: table_x(json), table_y(json), eval_x(json), eval_y(json)

! i0..i1: the zones changed, the region widened through the convection
! zones that zone_begin and zone_end lie in, if they do.
      if (.not.convective_flag(zone_begin) .or. zone_begin.eq.1) then
         i0 = zone_begin
      else
         do i = zone_begin - 1, 1, -1
            if (.not.convective_flag(i)) then
               i0 = i + 1
               exit
            end if
         end do
         if (i < (1)) then
         i0 = 1
         end if
      end if
      if (.not.convective_flag(zone_end) .or. zone_end.eq.num_points) then
         i1 = zone_end
      else
         do i = zone_end+1, num_points
            if (.not.convective_flag(i)) then
               i1 = i - 1
               exit
            end if
         end do
         if (i > num_points) then
         i1 = num_points
         end if
      end if
! dJ/J, splined in chi from the grid points to the region's model zones
! (eval_y(1:nmod)).
      nmod = zone_end - zone_begin + 1
      do i = 1, nmod
         eval_x(i) = rot_scr%chi(i)
      end do
      do i = 1, rot_scr%num_eq_points
         table_x(i) = rot_scr%eq_chi(i)
         table_y(i) = delta_angular_momentum(i)/angular_momentum(i)
      end do
      call osplin(eval_x, eval_y, table_x, table_y, rot_scr%num_eq_points, nmod)
! j -> j (1 + dJ/J) over i0..i1, summing the zones' J (sumjmod) and its
! change (sumdjmod) on the way; the zones beyond the region take the
! change at its end zone.
      sumjmod = 0.0d0
      sumdjmod = 0.0d0
      if (i0.lt.zone_begin) then
         do i = i0, zone_begin-1
            sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
            sumdjmod = sumdjmod + eval_y(1)*specific_angular_momentum(i)* &
                 shell_mass(i)
            specific_angular_momentum(i) = specific_angular_momentum(i)* &
                 (1.0d0+eval_y(1))
         end do
      end if
      do i = zone_begin, zone_end
         ii = i - zone_begin + 1
         sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
         sumdjmod = sumdjmod + eval_y(ii)*specific_angular_momentum(i)* &
              shell_mass(i)
         specific_angular_momentum(i) = specific_angular_momentum(i)* &
              (1.0d0+eval_y(ii))
      end do
      if (i1.gt.zone_end) then
         do i = zone_end+1, i1
            sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
            sumdjmod = sumdjmod + eval_y(nmod)*specific_angular_momentum(i)* &
                 shell_mass(i)
            specific_angular_momentum(i) = specific_angular_momentum(i)* &
                 (1.0d0+eval_y(nmod))
         end do
      end if
! Conservation: the model's change of J, sumdjmod, should equal the
! grid's, total_delta_angular_momentum (the surface loss). The
! tolerance is 1e-5 of that total, so without a loss it is zero and the
! scaling below always runs, removing the interpolation's small error.
      test = abs(sumdjmod - total_delta_angular_momentum)
      val = 1.0d-5*abs(total_delta_angular_momentum)
! Scale j over i0..i1 by one factor so that the zones' total J becomes
! sumjmod + total_delta_angular_momentum.
      if (test.gt.val) then
         ratio = (sumjmod + total_delta_angular_momentum)/ &
              (sumjmod + sumdjmod)
         do i = i0, i1
            specific_angular_momentum(i) = specific_angular_momentum(i)*ratio
         end do
      end if
      return
end subroutine equal_grid_to_model

end module am_grid
