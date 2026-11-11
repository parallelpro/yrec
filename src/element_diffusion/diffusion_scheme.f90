!----------------------------------------------------------------------
! diffusion_scheme
!----------------------------------------------------------------------
! The numerical scheme both settling pipelines share, for one species
! on an equal grid of n centres (enclosed mass m_i) and n - 1 midpoints
! (m_(i+1/2)), in Bahcall & Loeb units, for
!     dX/dt = (4 pi/dm) d/dr [ D1 + D2 dX/dr ]:
!   lax_wendroff_step1, lax_wendroff_step2  the drift term D1 by the
!       two-step Lax-Wendroff scheme (Press et al. 1986, Numerical
!       Recipes, p. 633): the change over half the step at the
!       midpoints, then over the step at the centres from the midpoints'
!       coefficients;
!   implicit_diffusion_coeffs  the tridiagonal system of the diffusion
!       term D2, solved implicitly by the caller.
! Each step is in flux form: the centre changes, weighted by the cells'
! masses (m_(3/2), m_(i+1/2) - m_(i-1/2), M - m_(n-1/2)), sum to zero,
! so the scheme conserves the species' mass on the grid.
module diffusion_scheme
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: c4pi
      implicit none
      private
      public :: lax_wendroff_step1, lax_wendroff_step2, implicit_diffusion_coeffs

contains

! ---------------------------------------------------------------
! The first Lax-Wendroff step: the drift's change of the abundance over
! half the step, timestep/2, at the n - 1 midpoints (n = num_eq_points),
!     dX_(i+1/2) = 4 pi (timestep/2) (D1_(i+1) - D1_i)/(m_(i+1) - m_i),
! from the coefficients diffusion_coeff1 at the centres (enclosed mass
! eq_mass), and at the ends
!     dX_(3/2) = 4 pi (timestep/2) D1_2/m_2,
!     dX_(n-1/2) = -4 pi (timestep/2) D1_(n-1)/(M - m_(n-1))
! (total_mass M).
pure subroutine lax_wendroff_step1(timestep, diffusion_coeff1, eq_mass, num_eq_points, &
     total_mass, diffused_abundance_mid)
      real(dp), intent(in) :: timestep
      real(dp), intent(in) :: diffusion_coeff1(:), eq_mass(:)
      integer, intent(in) :: num_eq_points
      real(dp), intent(in) :: total_mass
      real(dp), intent(inout) :: diffused_abundance_mid(:)
      real(dp) :: dt_half, zone_mass
      integer :: i

      dt_half = 0.5d0*timestep*c4pi
      zone_mass = eq_mass(2)
      diffused_abundance_mid(1) = dt_half*diffusion_coeff1(2)/zone_mass
      do i = 2, num_eq_points-2
         zone_mass = eq_mass(i+1) - eq_mass(i)
         diffused_abundance_mid(i) = dt_half*(diffusion_coeff1(i+1) - diffusion_coeff1(i))/zone_mass
      end do
      zone_mass = total_mass - eq_mass(num_eq_points-1)
      diffused_abundance_mid(num_eq_points-1) = &
           -dt_half*diffusion_coeff1(num_eq_points-1)/zone_mass
end subroutine lax_wendroff_step1

! ---------------------------------------------------------------
! The second Lax-Wendroff step: the drift's change over the step at the
! n centres from the coefficients diffusion_coeff1_mid at the midpoints
! (enclosed mass eq_mass_mid),
!     dX_i = 4 pi timestep (D1_(i+1/2) - D1_(i-1/2))/(m_(i+1/2) - m_(i-1/2)),
! added to diffused_abundance; zero flux below the first midpoint (over
! m_(3/2)) and above the last (over M - m_(n-1/2)).
pure subroutine lax_wendroff_step2(timestep, diffusion_coeff1_mid, eq_mass_mid, &
     diffused_abundance, num_eq_points, total_mass)
      real(dp), intent(in) :: timestep
      real(dp), intent(in) :: diffusion_coeff1_mid(:), eq_mass_mid(:)
      real(dp), intent(inout) :: diffused_abundance(:)
      integer, intent(in) :: num_eq_points
      real(dp), intent(in) :: total_mass
      real(dp) :: dt_full, zone_mass, delta_abundance
      integer :: i

      dt_full = timestep*c4pi
      zone_mass = eq_mass_mid(1)
      delta_abundance = dt_full*diffusion_coeff1_mid(1)/zone_mass
      diffused_abundance(1) = diffused_abundance(1) + delta_abundance
      do i = 2, num_eq_points-1
         zone_mass = eq_mass_mid(i) - eq_mass_mid(i-1)
         delta_abundance = dt_full*(diffusion_coeff1_mid(i) - diffusion_coeff1_mid(i-1))/zone_mass
         diffused_abundance(i) = diffused_abundance(i) + delta_abundance
      end do
      zone_mass = total_mass - eq_mass_mid(num_eq_points-1)
      delta_abundance = -dt_full*diffusion_coeff1_mid(num_eq_points-1)/zone_mass
      diffused_abundance(num_eq_points) = diffused_abundance(num_eq_points) + delta_abundance
end subroutine lax_wendroff_step2

! ---------------------------------------------------------------
! The tridiagonal system of the implicit diffusion step over npt
! centres,
!     X_i - alpha_i [D_(i+1/2) (X_(i+1) - X_i) - D_(i-1/2) (X_i - X_(i-1))] = X'_i,
! no flux below the first centre or above the last: its diagonals
! sub_diag, diag and super_diag from alpha (4 pi timestep/(spacing dm_i),
! the caller's) and the midpoints' coefficients diffusion_coeff_mid,
! which are first corrected to first order for the abundance's change
! at the midpoints, D += (dD/dX) dX (diffusion_coeff_deriv_mid times
! delta_abundance_mid; microdiff's derivative is zero).
pure subroutine implicit_diffusion_coeffs(alpha, diffusion_coeff_mid, delta_abundance_mid, &
     diffusion_coeff_deriv_mid, sub_diag, diag, super_diag, npt)
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(inout) :: diffusion_coeff_mid(:)
      real(dp), intent(in) :: delta_abundance_mid(:)
      real(dp), intent(in) :: diffusion_coeff_deriv_mid(:)
      real(dp), intent(inout) :: sub_diag(:), diag(:), super_diag(:)
      integer, intent(in) :: npt
      integer :: i

      do i = 1, npt-1
         diffusion_coeff_mid(i) = diffusion_coeff_mid(i) + &
              diffusion_coeff_deriv_mid(i)*delta_abundance_mid(i)
      end do
      sub_diag(1) = 0.0d0
      diag(1) = 1.0d0 + alpha(1)*diffusion_coeff_mid(1)
      super_diag(1) = -alpha(1)*diffusion_coeff_mid(1)
      do i = 2, npt-1
         sub_diag(i) = -alpha(i)*diffusion_coeff_mid(i-1)
         diag(i) = 1.0d0 + alpha(i)*(diffusion_coeff_mid(i-1) + diffusion_coeff_mid(i))
         super_diag(i) = -alpha(i)*diffusion_coeff_mid(i)
      end do
      sub_diag(npt) = -alpha(npt)*diffusion_coeff_mid(npt-1)
      diag(npt) = 1.0d0 + alpha(npt)*diffusion_coeff_mid(npt-1)
      super_diag(npt) = 0.0d0
end subroutine implicit_diffusion_coeffs

end module diffusion_scheme

! Sources
! - The scheme: MHP (12/89, after Bahcall & Loeb 1990; 3/94, the metals).
! - Press, W. H. et al. 1986, Numerical Recipes, Cambridge University
!   Press, p. 633 (the two-step Lax-Wendroff scheme).
