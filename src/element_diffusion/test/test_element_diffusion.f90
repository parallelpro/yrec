!----------------------------------------------------------------------
! test_element_diffusion
!----------------------------------------------------------------------
! Standalone test of element diffusion, printed for byte-comparison
! against expected_test_element_diffusion.out. Both settling routines run one
! step on a synthetic 200-zone solar-like model (radiative interior,
! convective envelope above 0.71 R):
!   gravitational_settling  the original routines (use_new_diffusion_
!                           routines false): helium only, helium and
!                           metals, each Coulomb-logarithm choice, the
!                           constant mixing and settling factors, the
!                           helium and metal scale factors, the
!                           Thoul et al. fit and the older analytic
!                           coefficients instead of the Thoul et al.
!                           ones, and a solve stopped short of
!                           convergence
!   microdiff               the newer routines (Bahcall & Loeb 1990
!                           units): helium only, helium and metals, the
!                           light elements (they ignore use_thoul_fit and
!                           use_thoul_diffusion)
! and, for both, the abundance floors: the outer edge of settling where
! Y first exceeds helium_diffusion_min, and no settling when Y or X is
! below its floor everywhere.
! Each case prints the change of H, He4, Z, Li7 and Be9 at nine zones
! and the mass-weighted change of each species over the star, all at
! full precision. No deck pins the original routines.
program test_element_diffusion
      use star_info_lib, only: star_info, json, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, i_li6, &
           i_li7, i_be9
      use run_messages, only: run_log_unit, terminal_unit
      use element_diffusion_lib, only: diffuse_elements, set_diffusion_settings
      use diffusion_settings_builder, only: diffusion_settings_from_star
      implicit none
! the test's star (star_info_lib declares no instance); static, as it
! is too large for the stack
      type(star_info), save :: star
      integer, parameter :: nz = 200
      double precision, parameter :: msun = 1.989d33, rsun = 6.9598d10, &
           year = 3.1558d7, cgs_g = 6.6726d-8, k_over_mu_mh = 1.3806d-16/(0.62d0*1.6605d-24)
      double precision :: xa0(15,json), xa(15,json), dlnp_dr(json), log_r(json), &
           log_rho(json), m(json), log_t(json), del(json), dm(json)
      logical :: convective(json)
      double precision :: dt, total_mass
      integer :: i
! the checks that printed FAIL
      integer :: failures = 0
! each zone's mass for the conservation check: the cell centred on the
! zone, from (m(k-1) + m(k))/2 to (m(k) + m(k+1))/2, the first from the
! centre and the last to the surface
      double precision :: dm_centred(json)

      run_log_unit = 6
      terminal_unit = 6
      call build_model()

      call reset_controls()
      star%physics%diffuse_helium_active = .true.
      call settle_case('old: helium', .false.)
      star%physics%use_diffusion_z = .true.
      call settle_case('old: helium and metals', .false.)
      star%physics%coulomb_log_choice = 2
      call settle_case('old: Coulomb logarithm choice 2', .false.)
      star%physics%coulomb_log_choice = 3
      call settle_case('old: Coulomb logarithm choice 3', .false.)
      star%physics%coulomb_log_choice = 4
      call settle_case('old: Coulomb logarithm choice 4', .false.)
      star%physics%coulomb_log_choice = 1
      star%physics%constant_mixing_coeff = 0.5d0
      star%physics%constant_settling_reduction = 0.8d0
      call settle_case('old: constant mixing 0.5, settling reduction 0.8', .false.)
      star%physics%constant_mixing_coeff = 1.0d0
      star%physics%constant_settling_reduction = 1.0d0
      star%physics%settling_y_factor = 0.7d0
      star%settling_z_factor_in_use = 0.9d0
      call settle_case('old: helium factor 0.7, metal factor 0.9', .false.)
      star%physics%settling_y_factor = 1.0d0
      star%settling_z_factor_in_use = 1.0d0
      star%physics%use_thoul_fit = .true.
      call settle_case('old: the Thoul et al. fit for metals', .false.)
      star%physics%use_thoul_fit = .false.
      star%physics%use_thoul_diffusion = .false.
      call settle_case('old: without the Thoul et al. coefficients', .false.)
      star%physics%use_thoul_diffusion = .true.
      star%physics%settling_num_iterations = 1
      call settle_case('old: one iteration, short of convergence', .false.)
      star%physics%settling_num_iterations = 10
      call floors(.false.)

      call reset_controls()
      star%physics%diffuse_helium_active = .true.
      call settle_case('new: helium', .true.)
      star%physics%use_diffusion_z = .true.
      call settle_case('new: helium and metals', .true.)
      star%physics%diffuse_lithium = .true.
      star%physics%coulomb_log_choice = 4
      call settle_case('new: helium, metals and light elements', .true.)
      call floors(.true.)

      write(*,'(/,a)') '# checks'
      call pipeline_checks()
      call region_and_substep_checks()
      call grid_checks()
      call scheme_checks()
      call coulomb_checks()
      call thoul_checks()
      if (failures > 0) then
         write(*,'(i0,a)') failures, ' check(s) failed'
         stop 1
      end if

contains

! A solar-like model: center (zone 1) to surface (zone nz), the mass
! and pressure gradient consistent with an ideal gas in hydrostatic
! equilibrium, the composition with a burned core and the light
! elements present only outside it.
subroutine build_model()
      double precision :: f, r(json), rho(json), p
      xa0 = 0.0d0
      do i = 1, nz
         f = dble(i - 1)/dble(nz - 1)
         r(i) = rsun*(0.002d0 + 0.998d0*f)
         log_r(i) = log10(r(i))
         log_t(i) = 7.196d0 - 0.85d0*f**1.3d0 - 2.6d0*f**25
         log_rho(i) = 2.18d0 - 3.8d0*f**1.1d0 - 5.0d0*f**30
         rho(i) = 10.0d0**log_rho(i)
         convective(i) = f > 0.71d0
         del(i) = merge(0.40d0, 0.25d0 + 0.1d0*f, convective(i))
         xa0(i_h1,i) = 0.71d0 - 0.35d0*exp(-(f/0.12d0)**2)
         xa0(i_metals,i) = 0.018d0
         xa0(i_he3,i) = 1.0d-5*(1.0d0 + 4.0d0*f)
         xa0(i_c12,i) = 3.0d-3
         xa0(i_c13,i) = 4.0d-5
         xa0(i_n14,i) = 1.0d-3
         xa0(i_n15,i) = 4.0d-6
         xa0(i_o16,i) = 9.0d-3
         xa0(i_o17,i) = 4.0d-6
         xa0(i_o18,i) = 2.0d-5
         if (f > 0.6d0) then
            xa0(i_h2,i) = 2.8d-5
            xa0(i_li6,i) = 7.0d-10
            xa0(i_li7,i) = 1.0d-8
            xa0(i_be9,i) = 1.6d-10
         end if
         xa0(i_he4,i) = 1.0d0 - xa0(i_h1,i) - xa0(i_metals,i) - xa0(i_he3,i)
      end do
! the mass from the density (trapezoid rule), then the pressure
! gradient of an ideal gas in hydrostatic equilibrium
      m(1) = 4.0d0/3.0d0*acos(-1.0d0)*r(1)**3*rho(1)
      do i = 2, nz
         m(i) = m(i-1) + 2.0d0*acos(-1.0d0)*(r(i)**2*rho(i) + r(i-1)**2*rho(i-1))*(r(i) - r(i-1))
      end do
      total_mass = m(nz)
      do i = 1, nz
         p = rho(i)*k_over_mu_mh*10.0d0**log_t(i)
         dlnp_dr(i) = -cgs_g*m(i)*rho(i)/(p*r(i)*r(i))
      end do
      dm(1) = m(1)
      do i = 2, nz
         dm(i) = m(i) - m(i-1)
      end do
      dm_centred(1) = 0.5d0*(m(1) + m(2))
      do i = 2, nz-1
         dm_centred(i) = 0.5d0*(m(i+1) - m(i-1))
      end do
      dm_centred(nz) = total_mass - 0.5d0*(m(nz-1) + m(nz))
      star%gradT(1:nz) = del(1:nz)
      dt = 1.0d8*year
end subroutine build_model

! The settling controls at their defaults (the tolerance 1e-7: the original
! solve stalls near 4e-8 at this model's centre), the switches off.
subroutine reset_controls()
      star%physics%diffuse_helium_active = .false.
      star%physics%use_diffusion_z = .false.
      star%physics%diffuse_lithium = .false.
      star%physics%use_thoul_fit = .false.
      star%physics%use_thoul_diffusion = .true.
      star%physics%coulomb_log_choice = 1
      star%physics%settling_num_iterations = 10
      star%physics%settling_tolerance = 1.0d-7
      star%physics%hydrogen_diffusion_floor = 1.0d-3
      star%physics%helium_diffusion_min = 1.0d-3
      star%physics%constant_mixing_coeff = 1.0d0
      star%physics%constant_settling_reduction = 1.0d0
      star%physics%settling_y_factor = 1.0d0
      star%settling_z_factor_in_use = 1.0d0
      star%settling_suspended_reported = .false.
end subroutine reset_controls

! The abundance floors: settling only below the outermost zone with
! Y > helium_diffusion_min, and none at all (a run-log note) when Y or X
! is below its floor everywhere.
subroutine floors(new_routines)
      logical, intent(in) :: new_routines
      character(len=4) :: tag
      tag = merge('new:', 'old:', new_routines)
      star%physics%helium_diffusion_min = 0.5d0
      call settle_case(tag//' outer edge where Y > 0.5', new_routines)
      star%physics%helium_diffusion_min = 0.99d0
      star%settling_suspended_reported = .false.
      call settle_case(tag//' Y below its floor everywhere', new_routines)
      star%physics%helium_diffusion_min = 1.0d-3
      star%physics%hydrogen_diffusion_floor = 0.9d0
      star%settling_suspended_reported = .false.
      call settle_case(tag//' X below its floor everywhere', new_routines)
      star%physics%hydrogen_diffusion_floor = 1.0d-3
end subroutine floors

! One settling step from the pristine model; prints the changes.
subroutine settle_case(label, new_routines)
      character(len=*), intent(in) :: label
      logical, intent(in) :: new_routines
      double precision :: step, mass, local_dlnp_dr(json), local_m(json)
      integer :: ierr, k, s
      integer, parameter :: shown(9) = [1, 40, 80, 120, 140, 143, 160, 190, 200]
      integer, parameter :: species(5) = [i_h1, i_he4, i_metals, i_li7, i_be9]

      xa = xa0
      step = dt
      mass = total_mass
      local_dlnp_dr = dlnp_dr
      local_m = m
      ierr = 0
      star%physics%use_new_diffusion_routines = new_routines
! the controls set above, as a run card installs them
      call set_diffusion_settings(diffusion_settings_from_star(star))
      call diffuse_elements(step, xa, local_dlnp_dr, log_r, log_rho, local_m, &
           log_t, del, convective, nz, mass, &
           star%settling_suspended_reported, ierr)
      write(*,'(/,a)') '# '//label
      write(*,'(a,i0)') 'ierr = ', ierr
      write(*,'(a)') ' zone   dX  dY  dZ  dLi7  dBe9'
      do k = 1, size(shown)
         write(*,'(i5,5es25.16e3)') shown(k), (xa(species(s),shown(k)) - &
              xa0(species(s),shown(k)), s = 1, size(species))
      end do
      write(*,'(a)') ' species  mass-weighted change'
      do s = 1, 15
         write(*,'(i5,es25.16e3)') s, sum((xa(s,1:nz) - xa0(s,1:nz))*dm(1:nz))/total_mass
      end do
end subroutine settle_case

! ======================================================================
! The checks: each prints ok or FAIL (README.md section 4 lists what
! each checks); the program stops with status 1 after a FAIL.
! ======================================================================

subroutine report_check(label, holds)
      character(len=*), intent(in) :: label
      logical, intent(in) :: holds
      if (holds) then
         write(*,'(a,t76,a)') label, 'ok'
      else
         write(*,'(a,t76,a)') label, 'FAIL'
         failures = failures + 1
      end if
end subroutine report_check

! One microdiff step of step_dt seconds from the pristine model with the
! convective flags given (the controls as last installed); no printout.
subroutine settle_quietly(step_dt, flags, xa_out, ierr)
      double precision, intent(in) :: step_dt
      logical, intent(in) :: flags(:)
      double precision, intent(out) :: xa_out(15,json)
      integer, intent(out) :: ierr
      double precision :: step, mass, local_dlnp_dr(json), local_m(json)
      xa_out = xa0
      step = step_dt
      mass = total_mass
      local_dlnp_dr = dlnp_dr
      local_m = m
      star%settling_suspended_reported = .false.
      call diffuse_elements(step, xa_out, local_dlnp_dr, log_r, log_rho, local_m, &
           log_t, del, flags, nz, mass, star%settling_suspended_reported, ierr)
end subroutine settle_quietly

! The pipeline decks run (microdiff with helium, the metals and the
! light elements) on the synthetic model.
subroutine pipeline_checks()
      use settling_region, only: find_settling_region
      integer, parameter :: settled(6) = [i_h1, i_he4, i_metals, i_li6, i_li7, i_be9]
      double precision :: xa1(15,json), xa2(15,json), net, moved, ratio, d1, d2
      logical :: holds, core_flags(json)
      integer :: ierr, k, s, zone_begin, zone_end, region

      call reset_controls()
      star%physics%diffuse_helium_active = .true.
      star%physics%use_diffusion_z = .true.
      star%physics%diffuse_lithium = .true.
      star%physics%coulomb_log_choice = 4
      star%physics%use_new_diffusion_routines = .true.
      call set_diffusion_settings(diffusion_settings_from_star(star))

      call settle_quietly(dt, convective, xa1, ierr)
      call report_check('microdiff: X + Y + Z + He3 = 1 in every zone after a step', &
           ierr == 0 .and. maxval(abs(xa1(i_h1,1:nz) + xa1(i_he4,1:nz) + &
           xa1(i_metals,1:nz) + xa1(i_he3,1:nz) - 1.0d0)) <= 1.0d-15)
! Each species' net change of mass against the mass it moves (the
! scheme conserves exactly on its grid; the map back by interpolation
! keeps 1e-4 of it here).
      holds = ierr == 0
      do k = 1, size(settled)
         s = settled(k)
         net = sum((xa1(s,1:nz) - xa0(s,1:nz))*dm_centred(1:nz))
         moved = sum(abs(xa1(s,1:nz) - xa0(s,1:nz))*dm_centred(1:nz))
         holds = holds .and. moved > 0.0d0 .and. abs(net) <= 3.0d-4*moved
      end do
      call report_check('microdiff: each species'' mass conserved (3e-4 of the mass moved)', holds)
      call report_check('microdiff: helium, the metals and Li6, Li7, Be9 sink; hydrogen rises', &
           ierr == 0 .and. xa1(i_h1,nz) > xa0(i_h1,nz) .and. xa1(i_he4,nz) < xa0(i_he4,nz) &
           .and. xa1(i_metals,nz) < xa0(i_metals,nz) .and. xa1(i_li6,nz) < xa0(i_li6,nz) &
           .and. xa1(i_li7,nz) < xa0(i_li7,nz) .and. xa1(i_be9,nz) < xa0(i_be9,nz) &
           .and. xa1(i_he4,1) > xa0(i_he4,1))

      call settle_quietly(0.0d0, convective, xa1, ierr)
      call report_check('microdiff: a zero timestep changes nothing', &
           ierr == 0 .and. all(xa1(:,1:nz) == xa0(:,1:nz)))

! At 1e4 yr the surface changes of a step twice as long are twice as
! large to 1.5e-5 (the departure falls with the step).
      call settle_quietly(1.0d4*year, convective, xa1, ierr)
      call settle_quietly(2.0d4*year, convective, xa2, ierr)
      holds = ierr == 0
      do k = 1, size(settled)
         s = settled(k)
         d1 = xa1(s,nz) - xa0(s,nz)
         d2 = xa2(s,nz) - xa0(s,nz)
         ratio = 0.0d0
         if (d1 /= 0.0d0) ratio = d2/d1
         holds = holds .and. abs(ratio - 2.0d0) <= 1.0d-4
      end do
      call report_check('microdiff: the change is linear in a small timestep', holds)

! A convective core (zones 1-30): the region starts at its top, and the
! core's zones all take the region's first change.
      core_flags = convective
      core_flags(1:30) = .true.
      call find_settling_region(core_flags, xa0, nz, 1.0d-3, 1.0d-3, zone_begin, &
           zone_end, region)
      call settle_quietly(dt, core_flags, xa1, ierr)
      holds = ierr == 0 .and. zone_begin == 30 .and. xa1(i_he4,1) > xa0(i_he4,1)
      do k = 2, 30
         holds = holds .and. abs((xa1(i_h1,k) - xa0(i_h1,k)) - (xa1(i_h1,1) - xa0(i_h1,1))) <= 1.0d-15 &
              .and. abs((xa1(i_metals,k) - xa0(i_metals,k)) - (xa1(i_metals,1) - xa0(i_metals,1))) <= 1.0d-15
      end do
      holds = holds .and. maxval(abs(xa1(i_h1,1:nz) + xa1(i_he4,1:nz) + &
           xa1(i_metals,1:nz) + xa1(i_he3,1:nz) - 1.0d0)) <= 1.0d-15
      call report_check('microdiff: a convective core''s zones change together', holds)
end subroutine pipeline_checks

! The settling region and the sub-step plan.
subroutine region_and_substep_checks()
      use settling_region, only: find_settling_region, region_found, &
           region_fully_convective, region_no_hydrogen, region_no_helium
      use settling_substeps, only: settling_substep_plan, plan_settling_substeps
      use physical_constants, only: seconds_per_year
      integer, parameter :: n = 100
      double precision :: comp(15,n), mass_star, boundary_mass, log_tb, longest
      logical :: flags(n), holds
      integer :: zone_begin, zone_end, region, r1, r2, r3
      type(settling_substep_plan) :: plan, plan_q, plan_t

      comp = 0.0d0
      comp(i_h1,:) = 0.7d0
      comp(i_he4,:) = 0.28d0
      flags = .false.
      flags(1:10) = .true.
      flags(81:n) = .true.
      call find_settling_region(flags, comp, n, 1.0d-3, 1.0d-3, zone_begin, zone_end, region)
      holds = region == region_found .and. zone_begin == 10 .and. zone_end == 81
      comp(i_he4,71:n) = 0.0d0
      call find_settling_region(flags, comp, n, 1.0d-3, 1.0d-3, zone_begin, zone_end, region)
      holds = holds .and. region == region_found .and. zone_end == 70
      call report_check('find_settling_region: from a core''s top to the envelope''s base', holds)
      flags = .true.
      call find_settling_region(flags, comp, n, 1.0d-3, 1.0d-3, zone_begin, zone_end, r1)
      flags = .false.
      call find_settling_region(flags, comp, n, 0.8d0, 1.0d-3, zone_begin, zone_end, r2)
      call find_settling_region(flags, comp, n, 1.0d-3, 0.5d0, zone_begin, zone_end, r3)
      call report_check('find_settling_region: fully convective, X or Y below its floor: no region', &
           r1 == region_fully_convective .and. r2 == region_no_hydrogen .and. &
           r3 == region_no_helium)

! Michaud et al.'s timescale, 4.348e21 yr q/T^1.5, and its scalings.
      mass_star = 2.0d33
      boundary_mass = 0.98d0*mass_star
      log_tb = log10(2.0d6)
      plan = plan_settling_substeps(1.0d8*year, 0.1d0, mass_star, boundary_mass, log_tb)
      plan_q = plan_settling_substeps(1.0d8*year, 0.1d0, mass_star, 0.96d0*mass_star, log_tb)
      plan_t = plan_settling_substeps(1.0d8*year, 0.1d0, mass_star, boundary_mass, &
           log_tb + log10(4.0d0))
      call report_check('plan_settling_substeps: tau = 4.348e21 yr q/T^1.5, linear in q, T^-1.5', &
           abs(plan%timescale/(4.348d21*seconds_per_year*0.02d0/2.0d6**1.5d0) - 1.0d0) <= 1.0d-12 &
           .and. abs(plan_q%timescale/plan%timescale - 2.0d0) <= 1.0d-12 &
           .and. abs(plan_t%timescale/plan%timescale - 0.125d0) <= 1.0d-12)
! The count: an exact multiple of the longest sub-step takes that many
! (docs/bug_sweep.md O16), anything more one more, a short step one.
      longest = plan%longest_substep
      plan = plan_settling_substeps(4.0d0*longest, 0.1d0, mass_star, boundary_mass, log_tb)
      holds = plan%num_substeps == 4 .and. plan%substep == longest
      plan = plan_settling_substeps(4.5d0*longest, 0.1d0, mass_star, boundary_mass, log_tb)
      holds = holds .and. plan%num_substeps == 5
      plan = plan_settling_substeps(0.3d0*longest, 0.1d0, mass_star, boundary_mass, log_tb)
      holds = holds .and. plan%num_substeps == 1
      call report_check('plan_settling_substeps: the count is max(1, ceiling(dt/longest))', holds)
end subroutine region_and_substep_checks

! The equal grids and their map (diffusion_grid), on a model whose zones
! widen outwards.
subroutine grid_checks()
      use diffusion_grid, only: equal_grid_map, size_equal_grid, place_equal_grids, &
           build_equal_grid_map, map_field, locate_stencil
      integer, parameter :: n = 60, zb = 5, ze = 50
      double precision :: r(n), linear(n), cubic(n), centre(json), mid(json), &
           centre_value(json), mid_value(json), spacing, narrowest, target
      type(equal_grid_map) :: map
      integer :: k, num_points, k0, search_from
      logical :: holds, found

      do k = 1, n
         r(k) = 0.01d0 + (dble(k)/dble(n))**1.5d0
         linear(k) = 3.0d0 - 2.0d0*r(k)
         cubic(k) = 1.0d0 + r(k) - 2.0d0*r(k)**2 + 0.5d0*r(k)**3
      end do
      call size_equal_grid(r, zb, ze, zb+1, num_points, spacing)
      narrowest = minval(r(zb+1:ze) - r(zb:ze-1))
      call report_check('size_equal_grid: ceiling(span/narrowest zone) points, evenly spaced', &
           num_points == ceiling((r(ze) - r(zb))/narrowest) .and. &
           abs(spacing*dble(num_points-1) - (r(ze) - r(zb))) <= 1.0d-14)
      call place_equal_grids(r(zb), spacing, num_points, centre, mid)
      call build_equal_grid_map(r, zb, ze, n, num_points, centre, mid, map)
      call map_field(map, linear, centre_value, mid_value)
      holds = maxval(abs(centre_value(1:num_points) - (3.0d0 - 2.0d0*centre(1:num_points)))) <= 1.0d-13 &
           .and. maxval(abs(mid_value(1:num_points-1) - (3.0d0 - 2.0d0*mid(1:num_points-1)))) <= 1.0d-13
      call report_check('map_field: exact on a linear field, every point', holds)
      call map_field(map, cubic, centre_value, mid_value)
      holds = .true.
      do k = 2, num_points-1
         holds = holds .and. abs(centre_value(k) - (1.0d0 + centre(k) - 2.0d0*centre(k)**2 + &
              0.5d0*centre(k)**3)) <= 1.0d-13 .and. abs(mid_value(k) - (1.0d0 + mid(k) - &
              2.0d0*mid(k)**2 + 0.5d0*mid(k)**3)) <= 1.0d-13
      end do
      call report_check('map_field: exact on a cubic at the interior points (4-point Lagrange)', holds)
! The stencil: two table points below the target and two above, within
! the table.
      holds = .true.
      search_from = 1
      do k = 1, 40
         target = r(3) + (r(n-3) - r(3))*dble(k)/41.0d0
         call locate_stencil(r, n, target, search_from, n, k0, found)
         holds = holds .and. found .and. r(k0+1) < target .and. target <= r(k0+2)
      end do
      search_from = 1
      call locate_stencil(r, n, 0.0d0, search_from, n, k0, found)
      holds = holds .and. found .and. k0 == 1
      call locate_stencil(r, n, 2.0d0*r(n), search_from, n, k0, found)
      holds = holds .and. .not. found .and. k0 == n-3
      call report_check('locate_stencil: two points each side, kept within the table', holds)
end subroutine grid_checks

! The scheme's steps conserve the species' mass on the grid (flux form).
subroutine scheme_checks()
      use diffusion_scheme, only: lax_wendroff_step1, lax_wendroff_step2, &
           implicit_diffusion_coeffs
      use numerics_lib, only: solve_tridiagonal
      integer, parameter :: n = 40
      double precision :: centre_mass(n), mid_mass(n), coeff(n), coeff_mid(n), x(n), &
           x_drift(n), change_mid(n), cell(n), mid_cell(n), alpha(n), zero(n), &
           sub_diag(n), diag(n), super_diag(n), x_new(n), net1, net2, net3, scale
      double precision, parameter :: total = 1.0d0, step = 1.0d-3
      integer :: i, ierr

      do i = 1, n
         centre_mass(i) = 0.05d0 + 0.9d0*(dble(i-1)/dble(n-1))**2
         coeff(i) = 0.3d0 + 0.2d0*sin(0.37d0*dble(i))
         x(i) = 0.7d0 - 0.2d0*cos(0.21d0*dble(i))
      end do
      do i = 1, n-1
         mid_mass(i) = 0.5d0*(centre_mass(i) + centre_mass(i+1))
         coeff_mid(i) = 0.25d0 + 0.1d0*cos(0.53d0*dble(i))
      end do
! step 1's cells: the midpoints' (m_2, m_(i+1) - m_i, M - m_(n-1))
      call lax_wendroff_step1(step, coeff, centre_mass, n, total, change_mid)
      mid_cell(1) = centre_mass(2)
      do i = 2, n-2
         mid_cell(i) = centre_mass(i+1) - centre_mass(i)
      end do
      mid_cell(n-1) = total - centre_mass(n-1)
      net1 = sum(change_mid(1:n-1)*mid_cell(1:n-1))
      scale = sum(abs(change_mid(1:n-1))*mid_cell(1:n-1))
! step 2's and the implicit step's cells: the centres'
      cell(1) = mid_mass(1)
      do i = 2, n-1
         cell(i) = mid_mass(i) - mid_mass(i-1)
      end do
      cell(n) = total - mid_mass(n-1)
      x_drift = x
      call lax_wendroff_step2(step, coeff_mid, mid_mass, x_drift, n, total)
      net2 = sum((x_drift - x)*cell)
      do i = 1, n
         alpha(i) = 50.0d0*step/cell(i)
      end do
      zero = 0.0d0
      change_mid = 0.0d0
      call implicit_diffusion_coeffs(alpha, coeff_mid, change_mid, zero, sub_diag, diag, &
           super_diag, n)
      call solve_tridiagonal(sub_diag, diag, super_diag, x_drift, n, x_new, ierr)
      net3 = sum((x_new - x_drift)*cell)
      call report_check('diffusion_scheme: both Lax-Wendroff steps and the implicit step conserve', &
           ierr == 0 .and. abs(net1) <= 1.0d-12*scale .and. &
           abs(net2) <= 1.0d-10*sum(abs(x_drift - x)*cell) .and. &
           abs(net3) <= 1.0d-10*sum(abs(x_new - x_drift)*cell) .and. &
           maxval(abs(x_new - x_drift)) > 1.0d-6)
end subroutine scheme_checks

! The Coulomb logarithms: symmetric; the Debye length from k and e; the
! fit's two limits (x >> 1: 0.97494 ln x - 1.3592, so the H-H and He-He
! logarithms differ by 0.97494 ln 4; x << 1: 0.15249 x^1.2, so their
! ratio is 4^1.2).
subroutine coulomb_checks()
      use coulomb_logarithm, only: compute_concentrations, compute_coulomb_logarithms
      use species_table, only: thoul_h1, thoul_he4, thoul_fe, thoul_electron
      use physical_constants, only: m_proton_cgs
      double precision, parameter :: k_boltzmann = 1.380649d-16, e_charge = 4.80320471d-10
      double precision :: a(4), z(4), x(4), c(4), cl(4,4), electron_density, debye, xhh, expected
      double precision, parameter :: pi = acos(-1.0d0)

      a = [thoul_h1%weight, thoul_he4%weight, thoul_fe%weight, thoul_electron%weight]
      z = [thoul_h1%charge, thoul_he4%charge, thoul_fe%charge, thoul_electron%charge]
      x = [0.7d0, 0.28d0, 0.02d0, 0.0d0]
      call compute_concentrations(a, z, x, c)
! a hot, thin plasma: Debye screening, x of order 1e5
      call compute_coulomb_logarithms(a, z, c, 1.0d-6, 1.0d7, cl)
      electron_density = 1.0d-6/(m_proton_cgs*sum(a*c))
      debye = sqrt(k_boltzmann*1.0d7/(4.0d0*pi*e_charge**2*electron_density*sum(c*z**2)))
      xhh = 4.0d0*k_boltzmann*1.0d7*debye/e_charge**2
      expected = 0.81245d0*log(1.0d0 + 0.18769d0*xhh**1.2d0)
      call report_check('compute_coulomb_logarithms: symmetric, and the H-H one from k, e (1e-4)', &
           all(cl == transpose(cl)) .and. abs(cl(1,1)/expected - 1.0d0) <= 1.0d-4)
      call report_check('compute_coulomb_logarithms: weak coupling, ln L(H) - ln L(He) = 0.975 ln 4', &
           abs((cl(1,1) - cl(2,2))/(1.2d0*0.81245d0*log(4.0d0)) - 1.0d0) <= 1.0d-5)
! a cool, dense plasma: the interionic distance screens, x of order 1e-3
      call compute_coulomb_logarithms(a, z, c, 1.0d6, 1.0d4, cl)
      call report_check('compute_coulomb_logarithms: strong coupling, ln L(H)/ln L(He) = 4^1.2', &
           abs(cl(1,1)/cl(2,2)/4.0d0**1.2d0 - 1.0d0) <= 1.0d-3)
end subroutine coulomb_checks

! The Burgers solution (Thoul et al. 1994) for H, He, Fe and electrons:
! the velocities carry no current and no net mass for any gradient; the
! coefficients scale as 1/ln Lambda; hydrogen rises, helium and iron
! sink (A_P < 0 for H, > 0 for He and Fe, d ln P/dr < 0).
subroutine thoul_checks()
      use thoul_solver, only: thoul_diffusion
      use coulomb_logarithm, only: compute_concentrations
      use species_table, only: thoul_h1, thoul_he4, thoul_fe, thoul_electron
      double precision :: a(4), z(4), x(4), c(4), cl(4,4), ap(4), at(4), ac(4,4), &
           ap1(4), at1(4), ac1(4,4)
      logical :: holds
      integer :: ierr, ierr1, j

      a = [thoul_h1%weight, thoul_he4%weight, thoul_fe%weight, thoul_electron%weight]
      z = [thoul_h1%charge, thoul_he4%charge, thoul_fe%charge, thoul_electron%charge]
      x = [0.7d0, 0.28d0, 0.02d0, 0.0d0]
      call compute_concentrations(a, z, x, c)
      cl = 2.2d0
      call thoul_diffusion(a, z, x, cl, ap, at, ac, ierr)
      holds = ierr == 0 .and. abs(sum(a*c*ap)) <= 1.0d-12*sum(abs(a*c*ap)) &
           .and. abs(sum(z*c*ap)) <= 1.0d-12*sum(abs(z*c*ap)) &
           .and. abs(sum(a*c*at)) <= 1.0d-12*sum(abs(a*c*at)) &
           .and. abs(sum(z*c*at)) <= 1.0d-12*sum(abs(z*c*at))
      do j = 1, 4
         if (j == 2 .or. j == 4) cycle
         holds = holds .and. abs(sum(a*c*ac(:,j))) <= 1.0d-12*sum(abs(a*c*ac(:,j))) &
              .and. abs(sum(z*c*ac(:,j))) <= 1.0d-12*sum(abs(z*c*ac(:,j)))
      end do
      call report_check('thoul_diffusion: no current and no net mass flow', holds)
      cl = 1.0d0
      call thoul_diffusion(a, z, x, cl, ap1, at1, ac1, ierr1)
      call report_check('thoul_diffusion: the coefficients scale as 1/ln Lambda', ierr1 == 0 &
           .and. maxval(abs(2.2d0*ap - ap1)) <= 1.0d-12*maxval(abs(ap1)) &
           .and. maxval(abs(2.2d0*at - at1)) <= 1.0d-12*maxval(abs(at1)) &
           .and. maxval(abs(2.2d0*ac - ac1)) <= 1.0d-12*maxval(abs(ac1)))
      call report_check('thoul_diffusion: hydrogen rises, helium and iron sink', &
           ap(1) < 0.0d0 .and. ap(2) > 0.0d0 .and. ap(3) > 0.0d0)
end subroutine thoul_checks

end program test_element_diffusion
