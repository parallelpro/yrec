!----------------------------------------------------------------------
! equal_grid
!----------------------------------------------------------------------
! The pieces of the equally spaced chi grid that am_transport_grid
! (am_grid) and composition_grid (rotational_mixing) share: the grid
! itself (equal_spaced_grid), the abscissae for interpolating onto the
! grid's zone edges (edge_grid_abscissae), the edge interpolation of a
! diffusion coefficient (interp_edge_coeff), and the Jacobian factor
! rho r^2 d chi/dr that a diffusion coefficient needs on the chi grid
! (dchi_dr_jacobian, multiply_by_exp10). The spline tables (table_x,
! table_y) and evaluation points (eval_x, eval_y) are the callers' own
! arrays, passed in as arguments. All interpolation uses osplin, the
! shape-preserving (osculatory) quadratic spline of numerics_lib.
!
! The blocks that differ between the two callers stay with them: the
! zone masses (am_transport_grid builds the moment of inertia and
! angular momentum in the same loops) and, under
! use_diffusion_advection_transport, the interpolation of
! es_advective_velocity and es_diffusive_velocity (a scale factor on
! each table entry).
module equal_grid
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: edge_grid_abscissae, interp_edge_coeff, dchi_dr_jacobian, &
           multiply_by_exp10, equal_spaced_grid
contains

! The abscissae for interpolating quantities at the model interfaces
! onto the zone edges of the equal grid:
!   table_x(1) = chi(1), table_x(ntabb) = chi(ntab): the region's end
!     points;
!   table_x(i) = (chi(i-1) + chi(i))/2, i = 2..ntab: the interfaces
!     between consecutive model points;
!   eval_x(1) = chi(1), and eval_x(i) = eq_chi(i) - eq_chi_spacing/2,
!     i = 2..num_eq_points: the edge between grid points i-1 and i.
! Returns ntabb = ntab + 1, the number of table entries.
subroutine edge_grid_abscissae(chi, eq_chi, eq_chi_spacing, ntab, num_eq_points, table_x, eval_x, ntabb)
      use star_info_lib, only: json
      implicit none
      double precision, intent(in) :: chi(json), eq_chi(json), eq_chi_spacing
      integer, intent(in) :: ntab, num_eq_points
      double precision, intent(inout) :: table_x(json), eval_x(json)
      integer, intent(out) :: ntabb
      integer :: i

      table_x(1) = chi(1)
      do i = 2,ntab
         table_x(i) = 0.5d0*(chi(i)+chi(i-1))
      end do
      ntabb = ntab + 1
      table_x(ntabb) = chi(ntab)
      eval_x(1) = chi(1)
      do i = 2, num_eq_points
         eval_x(i) = eq_chi(i)-0.5d0*eq_chi_spacing
      end do
end subroutine edge_grid_abscissae

! Interpolates a coefficient given at the model interfaces (coeff(k) at
! the interface between zones k-1 and k) onto the zone edges of the
! equal grid, with the abscissae from edge_grid_abscissae. The
! coefficient is held constant below the region's first unstable
! interface and above its last:
!   table_y(1) = coeff(zone_begin+1),
!   table_y(i) = coeff(zone_begin+i-1), i = 2..ntab,
!   table_y(ntabb) = coeff(zone_end).
! eq_coeff(i) is the value at the edge between grid points i-1 and i.
subroutine interp_edge_coeff(coeff, zone_begin, zone_end, ntab, ntabb, num_eq_points, &
     eval_x, table_x, table_y, eq_coeff)
      use star_info_lib, only: json
      use numerics_lib, only: osplin
      implicit none
      double precision, intent(in) :: coeff(json), eval_x(json), table_x(json)
      integer, intent(in) :: zone_begin, zone_end, ntab, ntabb, num_eq_points
      double precision, intent(inout) :: table_y(json)
      double precision, intent(out) :: eq_coeff(json)
      integer :: i, ii

      table_y(1) = coeff(zone_begin+1)
      do i = 2,ntab
         ii = zone_begin + i - 1
         table_y(i) = coeff(ii)
      end do
      table_y(ntabb) = coeff(zone_end)
      call osplin(eval_x,eq_coeff,table_x,table_y,ntabb,num_eq_points)
end subroutine interp_edge_coeff

! The Jacobian factor of a diffusion coefficient on the chi grid,
!   log10(rho r^2 d chi/dr),
! at the model points zone_begin..zone_begin+ntab-1, left in
! table_y(1:ntab) with table_x(1:ntab) = chi(1:ntab) (am_transport_grid
! reuses that table), and interpolated onto eval_x into eval_y.
! epsm is the energy generation rate per unit mass (erg/g/s), and
! surface_luminosity_lsun the caller's luminosity_lsun(num_zones), the
! L_tot of chi.
subroutine dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
     enclosed_mass, epsm, surface_luminosity_lsun, zone_begin, ntab, num_eq_points, &
     chi, eval_x, table_x, table_y, eval_y)
      use star_info_lib, only: json
      use controls_buffer, only: ichi_dl_max, ichi_dm_max, ichi_dp_core_max
      use physical_constants
      use numerics_lib, only: osplin
      use math_lib
      implicit none
      double precision, intent(in) :: log_density(json), log_radius(json), &
           log_mass(json), log_pressure(json), enclosed_mass(json), &
           epsm(json), surface_luminosity_lsun, chi(json), eval_x(json)
      integer, intent(in) :: zone_begin, ntab, num_eq_points
      double precision, intent(inout) :: table_x(json), table_y(json)
      double precision, intent(out) :: eval_y(json)
      integer :: i, ii
      double precision :: mass_scale, luminosity_scale, pressure_scale, &
           four_pi_rho_r2, dchi_dr

! The scales of chi (see equal_spaced_grid), with the luminosity scale
! dl L_tot in erg/s.
      mass_scale = rotation_settings%chi_grid_scale(ichi_dm_max)
      luminosity_scale = rotation_settings%chi_grid_scale(ichi_dl_max)*surface_luminosity_lsun* &
           rotation_settings%solar_luminosity_cgs
      pressure_scale = rotation_settings%chi_grid_scale(ichi_dp_core_max)
      do i = 1, ntab
         ii = zone_begin + i - 1
         table_x(i) = chi(i)
! With chi = log M/dm + L/(L_tot dl) - log P/dp and the structure
! equations dM/dr = 4 pi r^2 rho, dL/dr = 4 pi r^2 rho eps and
! dP/dr = -G M rho/r^2,
!   d chi/dr = 4 pi rho r^2/(ln10 dm M) + 4 pi rho r^2 eps/(dl L_tot)
!              + G M rho/(ln10 dp P r^2),
! and table_y = log10(rho r^2 d chi/dr).
         four_pi_rho_r2 = c4pi*exp(ln10*(log_density(ii)+2.0d0*log_radius(ii)))
         dchi_dr = four_pi_rho_r2/(ln10*mass_scale*enclosed_mass(ii))+ &
              four_pi_rho_r2*epsm(ii)/luminosity_scale+ &
              exp(ln10*(cgl+log_density(ii)+log_mass(ii)-log_pressure(ii)- &
              2.0d0*log_radius(ii)))/(ln10*pressure_scale)
         table_y(i) = log_density(ii) + log10(dchi_dr) + 2.0d0*log_radius(ii)
      end do
      call osplin(eval_x,eval_y,table_x,table_y,ntab,num_eq_points)
end subroutine dchi_dr_jacobian

! coeff(i) = coeff(i) 10**log_factor(i) at the num_eq_points grid points:
! applies a Jacobian factor held as a log10 (the eval_y of
! dchi_dr_jacobian) to a diffusion coefficient.
subroutine multiply_by_exp10(coeff, log_factor, num_eq_points)
      use star_info_lib, only: json
      use physical_constants, only: ln10
      use math_lib
      implicit none
      double precision, intent(inout) :: coeff(json)
      double precision, intent(in) :: log_factor(json)
      integer, intent(in) :: num_eq_points
      integer :: i

      do i = 1, num_eq_points
         coeff(i) = coeff(i)*exp(ln10*log_factor(i))
      end do
end subroutine multiply_by_exp10

! ---------------------------------------------------------------
! The equally spaced grid of the transport region zone_begin..zone_end.
! Its coordinate is
!   chi = log M/dm + L/(L_tot dl) - log P/dp,
! with L_tot the surface luminosity and dm, dl, dp the largest point
! spacings in log M, L/L_tot and log P that the rezoner allows
! (chi_grid_scale slots ichi_dm_max, ichi_dl_max, ichi_dp_core_max),
! so that the transport grid follows the structural grid. Sets in
! rot_scr: chi at the region's model points; num_eq_points values
! eq_chi, equally spaced (eq_chi_spacing) from chi(1) to the last chi;
! and eq_enclosed_mass, the enclosed mass (g) at the eq_chi points.
subroutine equal_spaced_grid(luminosity_lsun, log_pressure, log_mass, zone_begin, &
     zone_end, num_zones)
      use rotation_scratch_lib
      use star_info_lib, only: json
      use controls_buffer, only: ichi_dl_max, ichi_dm_max, ichi_dp_core_max
      use physical_constants
      use numerics_lib
      use math_lib
      implicit none

! luminosity_lsun (legacy HL) is the luminosity itself in solar units,
! not its log: both callers pass rot_scr%luminosity_lsun_mid,
! star%luminosity_lsun interpolated to the sub-step. The luminosity
! term of chi is therefore linear, the mass and pressure terms
! logarithmic.
      double precision, intent(in) :: luminosity_lsun(json), log_pressure(json), &
           log_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones

      double precision :: log_mass_table(json)
      integer :: num_points_in_range, zone_index, range_index
      double precision :: mass_scale, luminosity_scale, pressure_scale

! The scales of chi: dm, dl L_tot (in solar units) and dp.
      mass_scale = rotation_settings%chi_grid_scale(ichi_dm_max)
      luminosity_scale = rotation_settings%chi_grid_scale(ichi_dl_max)*luminosity_lsun(num_zones)
      pressure_scale = rotation_settings%chi_grid_scale(ichi_dp_core_max)
      num_points_in_range = zone_end - zone_begin + 1
      do range_index = 1, num_points_in_range
         zone_index = zone_begin + range_index - 1
         rot_scr%chi(range_index) = log_mass(zone_index)/mass_scale + &
              luminosity_lsun(zone_index)/luminosity_scale - &
              log_pressure(zone_index)/pressure_scale
      end do
! As many grid points as model points in the region. (The alternative,
! int(chi(last) - chi(1)) + 1 points, one per unit of chi, is not used.)
      rot_scr%num_eq_points = num_points_in_range
! At least three grid points.
      rot_scr%num_eq_points = max(rot_scr%num_eq_points,3)
! The grid spacing in chi.
      rot_scr%eq_chi_spacing = (rot_scr%chi(num_points_in_range)-rot_scr%chi(1))/dfloat(rot_scr%num_eq_points-1)
! The equally spaced chi values, from chi(1).
      rot_scr%eq_chi(1) = rot_scr%chi(1)
      do zone_index = 2, rot_scr%num_eq_points
         rot_scr%eq_chi(zone_index) = rot_scr%eq_chi(zone_index-1)+rot_scr%eq_chi_spacing
      end do
! The enclosed mass at the grid points, interpolated in log M, as chi is
! built from log M and log P.
      do zone_index = 1, num_points_in_range
         log_mass_table(zone_index) = log_mass(zone_begin+zone_index-1)
      end do
      call osplin(rot_scr%eq_chi,rot_scr%eq_enclosed_mass,rot_scr%chi,log_mass_table,num_points_in_range,rot_scr%num_eq_points)
! From log10 M to M (g).
      do zone_index = 1, rot_scr%num_eq_points
         rot_scr%eq_enclosed_mass(zone_index) = exp(ln10*rot_scr%eq_enclosed_mass(zone_index))
      end do
      return
end subroutine equal_spaced_grid

end module equal_grid
