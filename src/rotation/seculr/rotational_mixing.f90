!----------------------------------------------------------------------
! rotational_mixing
!----------------------------------------------------------------------
! The rotational mixing of the composition. Each transport region (a
! run of model interfaces with positive circulation velocity) is mapped
! onto an equally spaced grid in chi (composition_grid, with the shared
! pieces in equal_grid); each species that varies over the region is
! diffused there by a fully implicit difference scheme
! (composition_diffusion_coeffs, solved by solve_tridiagonal), and the
! change is mapped back to the model (diffuse_composition).
! diffuse_composition_driver finds the regions and runs these steps.
module rotational_mixing
      implicit none
      private
      public :: diffuse_composition_driver, composition_diffusion_coeffs

contains

! ---------------------------------------------------------------
! The rotational mixing of species species_begin..species_end over a
! sub-step of length timestep (s); composition is updated in place.
! Two modes, chosen by final_iteration_flag (the caller's
! diffusion_solve_ok):
!  - false, during the angular momentum iterations: mix the region
!    zone_begin..zone_end on the grid that am_transport_grid has just
!    built for it (rot_scr's chi grid, equally_spaced_diffusion_coeff
!    and equally_spaced_mass). The caller passes species 1..4 here (X,
!    Z, He3, with He4 set from them), the species that set the
!    gradient of the mean molecular weight.
!  - true, once the iterations are done: scan the interfaces from
!    zone_min outward for every transport region, build its grid with
!    composition_grid, and mix it. The caller passes the remaining
!    species (5 on) here.
!
! diffusion_coeff: the mixing coefficient 4 pi D_mix at the model
!    interfaces (from diffusion_velocity_scales; legacy COD2).
! equally_spaced_diffusion_coeff, equally_spaced_mass: the coefficient
!    at the grid edges (with the Jacobian factor; legacy ECOD2) and the
!    grid's zone masses (g; legacy EM). Inputs in the first mode,
!    rebuilt for each region in the second. At a convective boundary
!    the end point's mass includes the adjacent convection zone.
! log_density, luminosity_lsun, log_pressure, log_radius, log_mass,
!    enclosed_mass (g), shell_mass (g), log_total_mass: the model, for
!    composition_grid.
! velocity: the total circulation velocity (cm/s) at the interfaces,
!    element k between zones k-1 and k; transport regions are where it
!    is positive.
! zone_begin, zone_end: the zone below a region's first unstable
!    interface and the zone above its last (set here in the second
!    mode).
! zone_min, zone_max: the first and last interfaces scanned (the
!    innermost radiative zone outside zone 1 and the outermost
!    radiative zone).
! convective_flag: true where a zone is convective for angular
!    momentum transport (overshoot regions included).
! num_zones: the number of model zones.
! composition(15,:): the mass fractions (legacy HCOMP).
subroutine diffuse_composition_driver(rot_scr, mid_scr, &
     timestep, diffusion_coeff, equally_spaced_diffusion_coeff, &
     equally_spaced_mass, log_density, luminosity_lsun, log_pressure, &
     log_radius, log_mass, enclosed_mass, shell_mass, log_total_mass, &
     velocity, zone_begin, zone_end, zone_max, zone_min, convective_flag, &
     final_iteration_flag, num_zones, composition, species_begin, species_end, ierr)
      use rotation_scratch_lib, only: rotation_diffusion_state, midstep_physics_state
      use star_info_lib, only: json
      implicit none
      type(rotation_diffusion_state), intent(inout) :: rot_scr
      type(midstep_physics_state), intent(inout) :: mid_scr

      double precision, intent(in) :: timestep
      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(inout) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      double precision, intent(in) :: log_density(json), luminosity_lsun(json), &
           log_pressure(json), log_radius(json), log_mass(json), &
           enclosed_mass(json), shell_mass(json), log_total_mass
      double precision, intent(in) :: velocity(json)
      integer, intent(inout) :: zone_begin, zone_end
      integer, intent(in) :: zone_max, zone_min
      logical, intent(in) :: convective_flag(json), final_iteration_flag
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: species_begin, species_end
      integer, intent(out) :: ierr
      logical :: unstable_region_active, unstable_zone_found, &
           two_zone_region
      integer :: search_start, zone_idx

      ierr = 0
      if (.not.final_iteration_flag) then
         call diffuse_composition(rot_scr, timestep, equally_spaced_diffusion_coeff, &
              equally_spaced_mass, zone_begin, zone_end, &
              convective_flag, final_iteration_flag, num_zones, composition, &
              species_begin, species_end, ierr)
         if (ierr /= 0) return
      else
! Scan the interfaces from zone_min outward and mix each transport
! region in turn. unstable_zone_found is set when the scan meets a
! positive velocity; zone_begin is the zone below the region's first
! positive velocity, zone_end the zone above its last.
         unstable_region_active = .false.
         search_start = zone_min
         region_loop: do
         unstable_zone_found = .false.
         do zone_idx = search_start,zone_max
            if (velocity(zone_idx).gt.0.0d0) then
               unstable_zone_found = .true.
               if (.not.unstable_region_active) then
! The start of a region.
                  unstable_region_active = .true.
                  zone_begin = zone_idx - 1
               end if
            else if (unstable_region_active) then
! The end of a region.
               zone_end = zone_idx - 1
               unstable_region_active = .false.
               search_start = zone_idx + 1
               exit
            end if
         end do
         if (zone_idx > zone_max) then
! The scan reached zone_max: a region still open ends there.
         if (unstable_region_active) zone_end = zone_max
         search_start = zone_max + 1
         end if
! Stop when the scan finds no further positive velocity.
         if (unstable_zone_found) then
! The equally spaced chi grid of the region.
         call composition_grid(rot_scr, mid_scr, diffusion_coeff, log_density, luminosity_lsun, &
              log_pressure, log_radius, log_mass, enclosed_mass, shell_mass, &
              log_total_mass, zone_begin, zone_end, convective_flag, &
              num_zones, equally_spaced_diffusion_coeff, equally_spaced_mass, &
              two_zone_region)
! A region of two zones (a single unstable interface) is skipped.
         if (two_zone_region) then
            if (search_start.le.zone_max) then
               cycle region_loop
            else
               exit region_loop
            end if
         end if
! Mix species species_begin..species_end over the region.
         call diffuse_composition(rot_scr, timestep, equally_spaced_diffusion_coeff, &
              equally_spaced_mass, zone_begin, zone_end, &
              convective_flag, final_iteration_flag, num_zones, composition, &
              species_begin, species_end, ierr)
         if (ierr /= 0) return
! On to the next region, unless the scan has reached zone_max.
         if (search_start.gt.zone_max) exit region_loop
         else
         exit region_loop
         end if
         end do region_loop
      end if
      return
end subroutine diffuse_composition_driver

! ---------------------------------------------------------------
! Mixes species species_begin..species_end over one transport region
! (zone_begin..zone_end) for a time timestep (s), on the equally spaced
! chi grid already in rot_scr (built by am_transport_grid or
! composition_grid). For each species that varies over the region:
!  1. interpolate X from the model points to the grid points;
!  2. solve the implicit diffusion equation on the grid
!     (composition_diffusion_coeffs, solve_tridiagonal);
!  3. interpolate the change X_new - X_old back to the model points and
!     add it to composition, then give the convection zones adjacent to
!     the region the new composition of the region's end zones.
! Before the last iteration (final_iteration_flag false) the caller
! passes species 1..4; He4 is not solved for but reset to
! 1 - X - Z - He3.
!
! equally_spaced_diffusion_coeff: the coefficient at the grid edges,
!    element i at the edge between grid points i-1 and i (g/s, the
!    Jacobian factor included).
! equally_spaced_mass: the grid's zone masses (g). At a convective
!    boundary the end point's mass includes the adjacent convection
!    zone.
! zone_begin, zone_end: the zone below the region's first unstable
!    interface and the zone above its last.
! convective_flag: true where a zone is convective for angular
!    momentum transport (overshoot regions included).
! num_zones: the number of model zones.
! composition(15,:): the mass fractions, updated in place.
subroutine diffuse_composition(rot_scr, timestep, equally_spaced_diffusion_coeff, &
     equally_spaced_mass, zone_begin, zone_end, &
     convective_flag, final_iteration_flag, num_zones, composition, &
     species_begin, species_end, ierr)
      use rotation_scratch_lib, only: rotation_diffusion_state
      use star_info_lib, only: json, i_h1, i_he4, i_metals, i_he3
      use numerics_lib
      implicit none
      type(rotation_diffusion_state), intent(inout) :: rot_scr

      double precision, intent(in) :: timestep
      double precision, intent(in) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json), final_iteration_flag
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: species_begin, species_end
      integer, intent(out) :: ierr

! The tridiagonal system (filled by composition_diffusion_coeffs) and
! its solution (by solve_tridiagonal, the Thomas algorithm).
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json)

      double precision :: equally_spaced_composition(json)
! The splines' tables and evaluation points (model to grid and back).
      double precision :: table_x(json), table_y(json), eval_x(json), eval_y(json)
      integer :: varying_species_id(15)
      integer :: num_varying_species, j_idx, zone_idx, ntab, species_num, &
           orig_zone_idx, i0, i1
      double precision :: test_value

      ierr = 0

! Find the species whose mass fraction varies over the region by more
! than 1e-14.
      num_varying_species = 0
      do j_idx = species_begin, species_end
         test_value = composition(j_idx,zone_begin)
! He4 (slot 2) is 1 - X - Z - He3, so it is not solved for.
         if (j_idx.ne.2) then
            do zone_idx = zone_begin+1, zone_end
               if (dabs(composition(j_idx,zone_idx)-test_value).gt.1.0d-14) &
                    then
                  num_varying_species = num_varying_species + 1
                  varying_species_id(num_varying_species) = j_idx
                  exit
               end if
            end do
         end if
      end do
      if (num_varying_species.ne.0) then
! Diffuse each varying species, all with the same coefficients.
      ntab = zone_end - zone_begin + 1
! X at the grid points at the start of the sub-step, interpolated from
! the model points with the osculatory spline osplin.
      do species_num = 1, num_varying_species
         do zone_idx = 1, ntab
            orig_zone_idx = zone_idx + zone_begin - 1
            table_x(zone_idx) = rot_scr%chi(zone_idx)
            table_y(zone_idx) = composition(varying_species_id(species_num), &
                 orig_zone_idx)
         end do
         call osplin(rot_scr%eq_chi, equally_spaced_composition, table_x, table_y, ntab, &
              rot_scr%num_eq_points)
! The implicit difference equations for X at the end of the sub-step.
         call composition_diffusion_coeffs(equally_spaced_diffusion_coeff, rot_scr%eq_chi_spacing, timestep, &
              equally_spaced_composition, equally_spaced_mass, rot_scr%num_eq_points, &
              sub_diag, diag, super_diag, rhs)
! Solve for X at the end of the sub-step on the grid.
         call solve_tridiagonal(sub_diag, diag, super_diag, rhs, rot_scr%num_eq_points, solution, ierr)
         if (ierr /= 0) return
! Map back to the model points. Extrapolating the new X past the end
! grid points can give large errors at the region's boundary points, so
! the change X_new - X_old at the grid points is interpolated and added
! to composition, rather than composition being replaced by X_new.
         do zone_idx = 1, rot_scr%num_eq_points
            table_x(zone_idx) = rot_scr%eq_chi(zone_idx)
            table_y(zone_idx) = solution(zone_idx) - &
                 equally_spaced_composition(zone_idx)
         end do
         do zone_idx = 1, ntab
            eval_x(zone_idx) = rot_scr%chi(zone_idx)
         end do
! eval_y: the change in X at the model points of the region.
         call osplin(eval_x, eval_y, table_x, table_y, rot_scr%num_eq_points, ntab)
! i0: the bottom of the convection zone that contains zone_begin
! (zone_begin itself when it is radiative or zone 1).
         if (convective_flag(zone_begin).and.zone_begin.gt.1) then
            do zone_idx = zone_begin-1, 1, -1
               if (.not.convective_flag(zone_idx)) then
                  i0 = zone_idx + 1
                  exit
               end if
            end do
            if (zone_idx .lt. 1) i0 = 1
         else
            i0 = zone_begin
         end if
! i1: the top of the convection zone that contains zone_end (zone_end
! itself when it is radiative or the surface zone).
         if (convective_flag(zone_end).and.zone_end.lt.num_zones) then
            do zone_idx = zone_end+1, num_zones
               if (.not.convective_flag(zone_idx)) then
                  i1 = zone_idx - 1
                  exit
               end if
            end do
            if (zone_idx .gt. num_zones) i1 = num_zones
         else
! (The F77 mixcom.f left i1 unset on this branch, so the He4 reset
! below covered no zones there.)
            i1 = zone_end
         end if
! Add the change to the model points of the region.
         do zone_idx = 1, ntab
            j_idx = zone_begin + zone_idx - 1
            composition(varying_species_id(species_num),j_idx) = &
                 composition(varying_species_id(species_num),j_idx) + &
                 eval_y(zone_idx)
         end do
! Give the convection zone below the region the new X of zone_begin.
         if (i0.lt.zone_begin) then
            do zone_idx = zone_begin-1, i0, -1
               composition(varying_species_id(species_num),zone_idx) = &
                    composition(varying_species_id(species_num),zone_begin)
            end do
         end if
! Give the convection zone above the region the new X of zone_end.
         if (zone_end.lt.i1) then
            do zone_idx = zone_end+1, i1
               composition(varying_species_id(species_num),zone_idx) = &
                    composition(varying_species_id(species_num),zone_end)
            end do
         end if
      end do
! He4 = 1 - X - Z - He3 wherever X, Z and He3 were mixed.
      if (.not.final_iteration_flag) then
         do zone_idx = i0, i1
            composition(i_he4,zone_idx) = 1.0d0 - composition(i_h1,zone_idx) - &
                 composition(i_metals,zone_idx) - composition(i_he3,zone_idx)
         end do
      end if
      end if
      return
end subroutine diffuse_composition

! ---------------------------------------------------------------
! The implicit difference equations for the rotational mixing of one
! species on the equally spaced chi grid. The diffusion equation
!   dX/dt = (1/(4 pi r^2 rho)) d/dr (4 pi r^2 rho D dX/dr),
! with X the mass fraction and D the mixing coefficient, becomes, with
! the grid zone masses m_i = 4 pi r^2 rho dr and the edge coefficient
! C = 4 pi r^2 rho D d chi/dr (diffusion_coeff, g/s),
!   X_i(n+1) - X_i(n) = dt/(dchi m_i) [ C_{i+1/2} (X_{i+1} - X_i)
!                                     - C_{i-1/2} (X_i - X_{i-1}) ](n+1),
! where n and n+1 are the start and end of the sub-step, dchi is the
! grid spacing, C_{i-1/2} = diffusion_coeff(i) and C_{i+1/2} =
! diffusion_coeff(i+1). No material crosses the bottom of the first
! point or the top of the last. The equations form the tridiagonal
! system
!   sub_diag(i) X_{i-1} + diag(i) X_i + super_diag(i) X_{i+1} = rhs(i),
! rhs = X(n), which diffuse_composition solves with solve_tridiagonal.
! am_diffusion_coeffs sets up the angular momentum analogue.
!
! diffusion_coeff: C at the grid edges (g/s).
! grid_spacing: the chi spacing of the grid, dchi (dimensionless).
! timestep: the sub-step (s).
! eq_composition: X at the grid points at the start of the sub-step.
! eq_mass: the grid's zone masses m_i (g).
! num_eq_points: the number of grid points.
subroutine composition_diffusion_coeffs(diffusion_coeff, grid_spacing, timestep, eq_composition, &
     eq_mass, num_eq_points, sub_diag, diag, super_diag, rhs)
      use star_info_lib, only: json
      implicit none

      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: grid_spacing, timestep
      double precision, intent(in) :: eq_composition(json), eq_mass(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      integer :: i
      double precision :: fact0, fact

! The right-hand side: X at the start of the sub-step.
      do i = 1,num_eq_points
         rhs(i) = eq_composition(i)
      end do
! First point: no flux through its lower boundary.
      fact0 = timestep/grid_spacing
      fact = fact0/eq_mass(1)
      sub_diag(1) = 0.0d0
      diag(1) = 1.0d0 + fact*diffusion_coeff(2)
      super_diag(1) = -fact*diffusion_coeff(2)
      do i = 2,num_eq_points-1
         fact = fact0/eq_mass(i)
         sub_diag(i) = -fact*diffusion_coeff(i)
         diag(i) = 1.0d0 + fact*(diffusion_coeff(i)+diffusion_coeff(i+1))
         super_diag(i) = -fact*diffusion_coeff(i+1)
      end do
! Last point: no flux through its upper boundary.
      fact = fact0/eq_mass(num_eq_points)
      sub_diag(num_eq_points) = -fact*diffusion_coeff(num_eq_points)
      diag(num_eq_points) = 1.0d0 + fact*diffusion_coeff(num_eq_points)
      super_diag(num_eq_points) = 0.0d0

      return
end subroutine composition_diffusion_coeffs

! ---------------------------------------------------------------
! The equally spaced chi grid of one transport region for
! diffuse_composition. Calls equal_spaced_grid for the grid itself (in
! rot_scr), then forms the grid's zone masses (equally_spaced_mass, g)
! and the mixing coefficient at the grid edges,
!   equally_spaced_diffusion_coeff = diffusion_coeff rho r^2 d chi/dr,
! with diffusion_coeff = 4 pi D at the model interfaces (from
! diffusion_velocity_scales), so that the edge coefficient is
! 4 pi r^2 rho D d chi/dr (g/s). When the region has only two zones (a
! single unstable interface) it sets single_interface_flag and does
! nothing else. am_transport_grid is the angular momentum counterpart.
subroutine composition_grid(rot_scr, mid_scr, diffusion_coeff, log_density, luminosity_lsun, &
     log_pressure, log_radius, log_mass, enclosed_mass, shell_mass, &
     log_total_mass, zone_begin, zone_end, convective_flag, num_zones, &
     equally_spaced_diffusion_coeff, equally_spaced_mass, &
     single_interface_flag)
      use rotation_scratch_lib, only: rotation_diffusion_state, midstep_physics_state
      use equal_grid
      use star_info_lib, only: json
      use physical_constants, only: ln10
      use math_lib
      implicit none
      type(rotation_diffusion_state), intent(inout) :: rot_scr
      type(midstep_physics_state), intent(inout) :: mid_scr

      double precision, intent(in) :: diffusion_coeff(json), &
           log_density(json), luminosity_lsun(json), log_pressure(json), &
           log_radius(json), log_mass(json), enclosed_mass(json), &
           shell_mass(json)
      double precision, intent(in) :: log_total_mass
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      logical, intent(out) :: single_interface_flag
      integer :: idx, search_idx, ntab, ntabb
      double precision :: em_top, em_bot
! the spline's table (the model points) and evaluation points (the grid)
      double precision :: table_x(json), table_y(json), eval_x(json), eval_y(json)

! A region with a single unstable interface is flagged and skipped.
      if (zone_end - zone_begin.le.1) then
         single_interface_flag = .true.
         return
      else
         single_interface_flag = .false.
      end if
! The equally spaced grid: chi, eq_chi and eq_enclosed_mass.
      call equal_spaced_grid(rot_scr, luminosity_lsun, log_pressure, log_mass, zone_begin, &
           zone_end, num_zones)
! The zone mass of each grid point. An interior point's zone runs from
! halfway to the point below to halfway to the point above.
      do idx = 2, rot_scr%num_eq_points-1
         equally_spaced_mass(idx) = 0.5d0*(rot_scr%eq_enclosed_mass(idx+1) - rot_scr%eq_enclosed_mass(idx-1))
      end do
! An end point's zone runs from halfway to its neighbouring grid point
! out to the model zone boundary beyond the region (halfway in enclosed
! mass between the end zone and the next model zone, or the centre or
! the surface), and takes in any convection zone directly beyond.
!
! The bottom point.
      em_top = 0.5d0*(rot_scr%eq_enclosed_mass(2) + rot_scr%eq_enclosed_mass(1))
      if (zone_begin.gt.1) then
         em_bot = 0.5d0*(enclosed_mass(zone_begin) + &
              enclosed_mass(zone_begin-1))
      else
         em_bot = 0.0d0
      end if
      equally_spaced_mass(1) = em_top - em_bot
! Add the convective shells directly below zone_begin, down to the first
! radiative one.
      if (zone_begin.gt.1) then
         do search_idx = zone_begin-1, 1, -1
            if (.not.convective_flag(search_idx)) exit
            equally_spaced_mass(1) = equally_spaced_mass(1) + &
                 shell_mass(search_idx)
         end do
      end if
! The top point.
      em_bot = 0.5d0*(rot_scr%eq_enclosed_mass(rot_scr%num_eq_points) + rot_scr%eq_enclosed_mass(rot_scr%num_eq_points-1))
      if (zone_end.lt.num_zones) then
         em_top = 0.5d0*(enclosed_mass(zone_end) + &
              enclosed_mass(zone_end+1))
      else
         em_top = exp(ln10*log_total_mass)
      end if
      equally_spaced_mass(rot_scr%num_eq_points) = em_top - em_bot
! Add the convective shells directly above zone_end, up to the first
! radiative one.
      if (zone_end.lt.num_zones) then
         do search_idx = zone_end+1, num_zones
            if (.not.convective_flag(search_idx)) exit
            equally_spaced_mass(rot_scr%num_eq_points) = equally_spaced_mass(rot_scr%num_eq_points) + &
                 shell_mass(search_idx)
         end do
      end if
! The mixing coefficient at the grid edges. The grid is equally spaced
! in chi, not in radius, so the coefficient carries the Jacobian factor
! of that change of variable.
      ntab = zone_end - zone_begin + 1
      call edge_grid_abscissae(rot_scr%chi, rot_scr%eq_chi, rot_scr%eq_chi_spacing, ntab, &
           rot_scr%num_eq_points, table_x, eval_x, ntabb)
! The model coefficient at the edges, held constant below the bottom
! unstable interface and above the top one.
      call interp_edge_coeff(diffusion_coeff, zone_begin, zone_end, ntab, ntabb, &
           rot_scr%num_eq_points, eval_x, table_x, table_y, &
           equally_spaced_diffusion_coeff)
! log10(rho r^2 d chi/dr) at the edges, in eval_y.
      call dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
           enclosed_mass, mid_scr%eps_plus_neu_grav, luminosity_lsun(num_zones), zone_begin, &
           ntab, rot_scr%num_eq_points, rot_scr%chi, eval_x, table_x, &
           table_y, eval_y)
! Apply the factor (diffusion_coeff already carries the 4 pi).
      call multiply_by_exp10(equally_spaced_diffusion_coeff, eval_y, rot_scr%num_eq_points)
      return
end subroutine composition_grid

end module rotational_mixing
