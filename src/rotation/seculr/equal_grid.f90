!----------------------------------------------------------------------
! equal_grid
!----------------------------------------------------------------------
! The equally spaced grids' shared pieces: the chi grid itself
! (equal_spaced_grid) and the blocks am_transport_grid and
! composition_grid both need, which take the spline tables
! (rot_scr%xtab/table_y/eval_x/eval_y) as arguments.
!
! The blocks that differ between the two callers stay inline there:
! the zone-mass construction (am_transport_grid interleaves the moment
! of inertia and angular momentum with the mass in the same loops) and
! the two es_advective/es_diffusive interpolations under
! use_diffusion_advection_transport (a scale factor on each table entry).
module equal_grid
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: edge_grid_abscissae, interp_edge_coeff, dchi_dr_jacobian, &
           multiply_by_exp10, equal_spaced_grid
contains

! Abscissae for the zone-edge interpolations: table_x(1) and table_x(ntabb)
! are the region's end points, table_x(2..ntab) the midpoints between
! consecutive chi values; eval_x(1) is the first point and eval_x(2..num_eq_points)
! the equal-grid points shifted back by half a spacing.  Sets
! ntabb = ntab + 1, the number of table entries.
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

! Interpolates a model-point coefficient onto the zone edges of the
! equal grid through the tables built by edge_grid_abscissae, assuming
! it constant below the bottom interface and above the top interface:
! table_y(1) = coeff(zone_begin+1), table_y(ntabb) = coeff(zone_end).
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

! log10 of the product rho r^2 eq_chi_spacing/dr at the model points
! zone_begin..zone_begin+ntab-1 (left in table_y(1:ntab), with
! table_x(1:ntab) = chi(1:ntab)), interpolated onto eval_x into eval_y.
! surface_luminosity_lsun is luminosity_lsun(num_zones) of the caller.
subroutine dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
     enclosed_mass, epsm, surface_luminosity_lsun, zone_begin, ntab, num_eq_points, &
     chi, eval_x, table_x, table_y, eval_y)
      use star_info_lib, only: star, json
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

! PRODUCT OF RHO R^2 BY D CHI/DR
      mass_scale = rotation_settings%chi_grid_scale(ichi_dm_max)
      luminosity_scale = rotation_settings%chi_grid_scale(ichi_dl_max)*surface_luminosity_lsun* &
           star%solar_luminosity_cgs
      pressure_scale = rotation_settings%chi_grid_scale(ichi_dp_core_max)
      do i = 1, ntab
         ii = zone_begin + i - 1
         table_x(i) = chi(i)
! D CHI/DR = 1/DM*( D LOG M/DR) + 1/DL*(DL/DR) - 1/DP*(D LOG P/DR)
! OR, USING FAC = 4*PI*RHO*R**2
! D CHI/DR = FAC/(LN 10 * DM * M) + FAC*EPSILON/DL + RHO*GM/(LN10*DP*R**2)
! STORED IN YVAL
         four_pi_rho_r2 = c4pi*exp(ln10*(log_density(ii)+2.0d0*log_radius(ii)))
         dchi_dr = four_pi_rho_r2/(ln10*mass_scale*enclosed_mass(ii))+ &
              four_pi_rho_r2*epsm(ii)/luminosity_scale+ &
              exp(ln10*(cgl+log_density(ii)+log_mass(ii)-log_pressure(ii)- &
              2.0d0*log_radius(ii)))/(ln10*pressure_scale)
         table_y(i) = log_density(ii) + log10(dchi_dr) + 2.0d0*log_radius(ii)
      end do
      call osplin(eval_x,eval_y,table_x,table_y,ntab,num_eq_points)
end subroutine dchi_dr_jacobian

! coeff(i) = coeff(i) * 10**log_factor(i) over the num_eq_points equal-grid
! points (the multiplicative Jacobian factor on a diffusion coefficient).
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
! This routine forms the sum of rot_scr%chi, defined as
! rot_scr%chi = log(M)/dm + L/(Ltot*dl) - log(P)/dp, where dm, dp, and dl are
! the point spacings for log mass, luminosity, and pressure
! respectively. This routine transforms these variables to an equally
! spaced grid from the first point to the last point. It returns a
! set of equally spaced rot_scr%chi values and their location in mass.
subroutine equal_spaced_grid(luminosity_lsun, log_pressure, log_mass, zone_begin, &
     zone_end, num_zones)
      use rotation_scratch_lib
      use star_info_lib, only: json
      use controls_buffer, only: ichi_dl_max, ichi_dm_max, ichi_dp_core_max
      use physical_constants
      use numerics_lib
      use math_lib
      implicit none

! luminosity_lsun (originally HL) is the LINEAR luminosity in solar
! units, not a log: both callers (composition_grid.f90,
! am_transport_grid.f90) pass rot_scr%luminosity_lsun_mid, the
! mid-sub-step interpolate of star%luminosity_lsun, and the chi term
! below is L/(LTOT*DL) (unlike the LOG(M) and LOG(P) terms).
      double precision, intent(in) :: luminosity_lsun(json), log_pressure(json), &
           log_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones

      double precision :: log_mass_table(json)
      integer :: num_points_in_range, zone_index, range_index
      double precision :: mass_scale, luminosity_scale, pressure_scale

! USE THE MODEL CRITERIA FOR ASSIGNING THE SPACING BETWEEN GRID POINTS.
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
! NUMBER OF EQUALLY SPACED POINTS = NUMBER OF MODEL POINTS IN RANGE
! (THE ORIGINAL INT(CHI(NTAB)-CHI(1))+1 RULE IS NOT USED).
      rot_scr%num_eq_points = num_points_in_range
! FOR ROTATION PURPOSES, DEFINE THE MINIMUM NUMBER OF
! EQUALLY SPACED SHELLS AS 3.
      rot_scr%num_eq_points = max(rot_scr%num_eq_points,3)
! EQUALLY SPACED INCREMENT IN CHI
      rot_scr%eq_chi_spacing = (rot_scr%chi(num_points_in_range)-rot_scr%chi(1))/dfloat(rot_scr%num_eq_points-1)
! ASSIGN VECTOR OF EQUALLY SPACED CHI
      rot_scr%eq_chi(1) = rot_scr%chi(1)
      do zone_index = 2, rot_scr%num_eq_points
         rot_scr%eq_chi(zone_index) = rot_scr%eq_chi(zone_index-1)+rot_scr%eq_chi_spacing
      end do
! NOW ASSIGN MASSES TO THE NEW EQUALLY SPACED GRID POINTS.
! PERFORM INTERPOLATION IN LOG MASS - CONSISTENT WITH CHI DEFINED
! IN LOG M, LOG P, ETC.
      do zone_index = 1, num_points_in_range
         log_mass_table(zone_index) = log_mass(zone_begin+zone_index-1)
      end do
      call osplin(rot_scr%eq_chi,rot_scr%eq_enclosed_mass,rot_scr%chi,log_mass_table,num_points_in_range,rot_scr%num_eq_points)
! TRANSFORM TO PHYSICAL MASS (GM)
      do zone_index = 1, rot_scr%num_eq_points
         rot_scr%eq_enclosed_mass(zone_index) = exp(ln10*rot_scr%eq_enclosed_mass(zone_index))
      end do
      return
end subroutine equal_spaced_grid

end module equal_grid
