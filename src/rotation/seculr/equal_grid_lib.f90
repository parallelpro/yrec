!----------------------------------------------------------------------
! equal_grid_lib
!----------------------------------------------------------------------
! 2026 readability sweep (W3, rotation): the blocks that am_transport_grid
! (rotgrid.f) and composition_grid (mixgrid.f) carried as verbatim copies,
! each taking the osplin scratch tables (rot_scr%xtab/ytab/xval/yval) as
! explicit arguments instead of reading rot_scr%.  Every statement is
! token-for-token the former inline code; the callers invoke each helper
! at the position the block occupied.
!
! Rotation-internal (seculr/); not a public rotation entry.
!
! The blocks that differ between the two callers stay inline there:
! the zone-mass construction (am_transport_grid interleaves the moment
! of inertia and angular momentum with the mass in the same loops) and
! the two es_advective/es_diffusive interpolations under
! use_diffusion_advection_transport (a scale factor on each table entry).
module equal_grid_lib
      implicit none
      private
      public :: edge_grid_abscissae, interp_edge_coeff, dchi_dr_jacobian, &
           multiply_by_exp10
contains

!----------------------------------------------------------------------
! edge_grid_abscissae
!----------------------------------------------------------------------
! Abscissae for the zone-edge interpolations: xtab(1) and xtab(ntabb)
! are the region's end points, xtab(2..ntab) the midpoints between
! consecutive chi values; xval(1) is the first point and xval(2..ntot)
! the equal-grid points shifted back by half a spacing.  Sets
! ntabb = ntab + 1, the number of table entries.
subroutine edge_grid_abscissae(chi, echi, dchi, ntab, ntot, xtab, xval, ntabb)
      use star_info_lib, only: json
      implicit none
      double precision, intent(in) :: chi(json), echi(json), dchi
      integer, intent(in) :: ntab, ntot
      double precision, intent(inout) :: xtab(json), xval(json)
      integer, intent(out) :: ntabb
      integer :: i

      xtab(1) = chi(1)
      do i = 2,ntab
         xtab(i) = 0.5d0*(chi(i)+chi(i-1))
      end do
      ntabb = ntab + 1
      xtab(ntabb) = chi(ntab)
      xval(1) = chi(1)
      do i = 2, ntot
         xval(i) = echi(i)-0.5d0*dchi
      end do
end subroutine edge_grid_abscissae

!----------------------------------------------------------------------
! interp_edge_coeff
!----------------------------------------------------------------------
! Interpolates a model-point coefficient onto the zone edges of the
! equal grid through the tables built by edge_grid_abscissae, assuming
! it constant below the bottom interface and above the top interface:
! ytab(1) = coeff(zone_begin+1), ytab(ntabb) = coeff(zone_end).
subroutine interp_edge_coeff(coeff, zone_begin, zone_end, ntab, ntabb, ntot, &
     xval, xtab, ytab, eq_coeff)
      use star_info_lib, only: json
      use numerics_lib, only: osplin
      implicit none
      double precision, intent(in) :: coeff(json), xval(json), xtab(json)
      integer, intent(in) :: zone_begin, zone_end, ntab, ntabb, ntot
      double precision, intent(inout) :: ytab(json)
      double precision, intent(out) :: eq_coeff(json)
      integer :: i, ii

      ytab(1) = coeff(zone_begin+1)
      do i = 2,ntab
         ii = zone_begin + i - 1
         ytab(i) = coeff(ii)
      end do
      ytab(ntabb) = coeff(zone_end)
      call osplin(xval,eq_coeff,xtab,ytab,ntabb,ntot)
end subroutine interp_edge_coeff

!----------------------------------------------------------------------
! dchi_dr_jacobian
!----------------------------------------------------------------------
! log10 of the product rho r^2 dchi/dr at the model points
! zone_begin..zone_begin+ntab-1 (left in ytab(1:ntab), with
! xtab(1:ntab) = chi(1:ntab)), interpolated onto xval into yval.
! surface_luminosity_lsun is luminosity_lsun(num_zones) of the caller.
subroutine dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
     enclosed_mass, epsm, surface_luminosity_lsun, zone_begin, ntab, ntot, &
     chi, xval, xtab, ytab, yval)
      use star_info_lib, only: star, json
      use controls_lib, only: ichi_dl_max, ichi_dm_max, ichi_dp_core_max
      use phys_const_lib
      use numerics_lib, only: osplin
      use math_lib
      implicit none
      double precision, intent(in) :: log_density(json), log_radius(json), &
           log_mass(json), log_pressure(json), enclosed_mass(json), &
           epsm(json), surface_luminosity_lsun, chi(json), xval(json)
      integer, intent(in) :: zone_begin, ntab, ntot
      double precision, intent(inout) :: xtab(json), ytab(json)
      double precision, intent(out) :: yval(json)
      integer :: i, ii
      double precision :: mass_scale, luminosity_scale, pressure_scale, &
           four_pi_rho_r2, dchi_dr

! PRODUCT OF RHO R^2 BY D CHI/DR
      mass_scale = star%ctrl%chi_grid_scale(ichi_dm_max)
      luminosity_scale = star%ctrl%chi_grid_scale(ichi_dl_max)*surface_luminosity_lsun* &
           star%solar_luminosity_cgs
      pressure_scale = star%ctrl%chi_grid_scale(ichi_dp_core_max)
      do i = 1, ntab
         ii = zone_begin + i - 1
         xtab(i) = chi(i)
! D CHI/DR = 1/DM*( D LOG M/DR) + 1/DL*(DL/DR) - 1/DP*(D LOG P/DR)
! OR, USING FAC = 4*PI*RHO*R**2
! D CHI/DR = FAC/(LN 10 * DM * M) + FAC*EPSILON/DL + RHO*GM/(LN10*DP*R**2)
! STORED IN YVAL
         four_pi_rho_r2 = c4pi*exp(ln10*(log_density(ii)+2.0d0*log_radius(ii)))
         dchi_dr = four_pi_rho_r2/(ln10*mass_scale*enclosed_mass(ii))+ &
              four_pi_rho_r2*epsm(ii)/luminosity_scale+ &
              exp(ln10*(cgl+log_density(ii)+log_mass(ii)-log_pressure(ii)- &
              2.0d0*log_radius(ii)))/(ln10*pressure_scale)
         ytab(i) = log_density(ii) + log10(dchi_dr) + 2.0d0*log_radius(ii)
      end do
      call osplin(xval,yval,xtab,ytab,ntab,ntot)
end subroutine dchi_dr_jacobian

!----------------------------------------------------------------------
! multiply_by_exp10
!----------------------------------------------------------------------
! coeff(i) = coeff(i) * 10**log_factor(i) over the ntot equal-grid
! points (the multiplicative Jacobian factor on a diffusion coefficient).
subroutine multiply_by_exp10(coeff, log_factor, ntot)
      use star_info_lib, only: json
      use phys_const_lib, only: ln10
      use math_lib
      implicit none
      double precision, intent(inout) :: coeff(json)
      double precision, intent(in) :: log_factor(json)
      integer, intent(in) :: ntot
      integer :: i

      do i = 1, ntot
         coeff(i) = coeff(i)*exp(ln10*log_factor(i))
      end do
end subroutine multiply_by_exp10

end module equal_grid_lib
