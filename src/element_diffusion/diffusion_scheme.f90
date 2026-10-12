!----------------------------------------------------------------------
! diffusion_scheme
!----------------------------------------------------------------------
! The numerical scheme both settling pipelines share: the two-step
! Lax-Wendroff advection of the first term and the coefficients of the
! implicit solve of the second.
module diffusion_scheme
      use element_diffusion_scratch, only: diff_scr
      implicit none
      private
      public :: lax_wendroff_step1, lax_wendroff_step2, implicit_diffusion_coeffs

contains

! ---------------------------------------------------------------
! First step of the two-step Lax-Wendroff diffusion method. Computes
! new abundances at zone midpoints using the Lax scheme:
! X(n+1/2,j+1/2) = 1/2(X(n,j+1/2) - 1/2(dt/dr)(cod1(n,j+1)-cod1(n,j))
! where n is the time variable, j is the spatial one, and cod1 is the
! diffusion coefficient.
! Variables: diffusion_coeff1 is the run of diffusion coefficients,
! eq_mass is the run of masses, and diffused_abundance_mid is the run
! of changes in the diffused species' abundance at the zone midpoints.
! MHP 3/94 added metal diffusion.
subroutine lax_wendroff_step1(timestep, diffusion_coeff1, eq_mass, num_eq_points, &
     total_mass, diffused_abundance_mid, use_generic_diffusion_vectors)

      use star_info_lib, only: star, json
      use physical_constants
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(in) :: diffusion_coeff1(json), eq_mass(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: total_mass
      double precision, intent(out) :: diffused_abundance_mid(json)
      logical, intent(in) :: use_generic_diffusion_vectors
      double precision :: dt_half, zone_mass, delta_abundance_mid, &
           delta_metal_abundance_mid
      integer :: i

      dt_half = 0.5d0*timestep*c4pi
! central boundary condition
      zone_mass = eq_mass(2)
      delta_abundance_mid = dt_half*diffusion_coeff1(2)/zone_mass
      diffused_abundance_mid(1) = delta_abundance_mid
! general case
      do i = 2,num_eq_points-2
         zone_mass = eq_mass(i+1)-eq_mass(i)
         delta_abundance_mid = dt_half*(diffusion_coeff1(i+1)-diffusion_coeff1(i))/zone_mass
         diffused_abundance_mid(i) = delta_abundance_mid
      end do
! surface boundary condition.
      zone_mass = total_mass-eq_mass(num_eq_points-1)
      delta_abundance_mid = -dt_half*diffusion_coeff1(num_eq_points-1)/zone_mass
      diffused_abundance_mid(num_eq_points-1) = delta_abundance_mid
! metal diffusion
      if(star%job%use_diffusion_z.and..not.use_generic_diffusion_vectors)then
! central boundary condition
         zone_mass = eq_mass(2)
         delta_metal_abundance_mid = dt_half*diff_scr%metal_diffusion_coeff1(2)/zone_mass
         diff_scr%metal_abundance_change_mid(1) = delta_metal_abundance_mid
! general case
         do i = 2,num_eq_points-2
            zone_mass = eq_mass(i+1)-eq_mass(i)
            delta_metal_abundance_mid = dt_half*(diff_scr%metal_diffusion_coeff1(i+1)-diff_scr%metal_diffusion_coeff1(i))/zone_mass
            diff_scr%metal_abundance_change_mid(i) = delta_metal_abundance_mid
         end do
! surface boundary condition.
         zone_mass = total_mass-eq_mass(num_eq_points-1)
         delta_metal_abundance_mid = -dt_half*diff_scr%metal_diffusion_coeff1(num_eq_points-1)/zone_mass
         diff_scr%metal_abundance_change_mid(num_eq_points-1) = delta_metal_abundance_mid
      endif
      return
end subroutine lax_wendroff_step1

! ---------------------------------------------------------------
! Second (corrector) step of the two-step Lax-Wendroff diffusion
! method. Using the zone-midpoint diffusion coefficients from
! lax_wendroff_step1's provisional solution, updates diffused_abundance in
! place with the full-step change in the diffused species' abundance.
! MHP 3/94 added metal diffusion.
subroutine lax_wendroff_step2(timestep, diffusion_coeff1_mid, eq_mass_mid, &
     diffused_abundance, num_eq_points, total_mass, &
     use_generic_diffusion_vectors)

      use star_info_lib, only: star, json
      use physical_constants
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(in) :: diffusion_coeff1_mid(json), &
           eq_mass_mid(json)
      double precision, intent(inout) :: diffused_abundance(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: total_mass
      logical, intent(in) :: use_generic_diffusion_vectors
      double precision :: dt_full, zone_mass, delta_abundance, &
           delta_metal_abundance
      integer :: i

      dt_full = timestep*c4pi
! central boundary condition
      zone_mass = eq_mass_mid(1)
      delta_abundance = dt_full*diffusion_coeff1_mid(1)/zone_mass
      diffused_abundance(1) = diffused_abundance(1)+delta_abundance
! general case
      do i = 2,num_eq_points-1
         zone_mass = eq_mass_mid(i)-eq_mass_mid(i-1)
         delta_abundance = dt_full*(diffusion_coeff1_mid(i)-diffusion_coeff1_mid(i-1))/zone_mass
         diffused_abundance(i) = diffused_abundance(i)+delta_abundance
      end do
! surface boundary condition.
      zone_mass = total_mass-eq_mass_mid(num_eq_points-1)
      delta_abundance = -dt_full*diffusion_coeff1_mid(num_eq_points-1)/zone_mass
      diffused_abundance(num_eq_points) = diffused_abundance(num_eq_points)+delta_abundance
! mhp 3/94 added metal diffusion.
      if(star%job%use_diffusion_z.and..not.use_generic_diffusion_vectors)then
         zone_mass = eq_mass_mid(1)
         delta_metal_abundance = dt_full*diff_scr%metal_diffusion_coeff1_mid(1)/zone_mass
         diff_scr%metal_abundance_change(1) = diff_scr%metal_abundance_change(1)+delta_metal_abundance
! general case
         do i = 2,num_eq_points-1
            zone_mass = eq_mass_mid(i)-eq_mass_mid(i-1)
            delta_metal_abundance = dt_full*(diff_scr%metal_diffusion_coeff1_mid(i)-diff_scr%metal_diffusion_coeff1_mid(i-1))/zone_mass
            diff_scr%metal_abundance_change(i) = diff_scr%metal_abundance_change(i)+delta_metal_abundance
         end do
! surface boundary condition.
         zone_mass = total_mass-eq_mass_mid(num_eq_points-1)
         delta_metal_abundance = -dt_full*diff_scr%metal_diffusion_coeff1_mid(num_eq_points-1)/zone_mass
         diff_scr%metal_abundance_change(num_eq_points) = diff_scr%metal_abundance_change(num_eq_points)+delta_metal_abundance
      endif
      return
end subroutine lax_wendroff_step2

! ---------------------------------------------------------------
! Rebuilds the tridiagonal-matrix coefficients for the implicit
! (second-derivative) term of the microscopic-diffusion equation
! (microdiff_run and gravitational_settling), correcting the midpoint
! diffusion coefficient for the change in abundance found on the
! previous iteration.
subroutine implicit_diffusion_coeffs(alpha, diffusion_coeff_mid, delta_abundance_mid, &
     diffusion_coeff_deriv_mid, sub_diag, diag, super_diag, npt)
      use star_info_lib, only: json

      implicit none

      double precision, intent(in) :: alpha(json)
      double precision, intent(inout) :: diffusion_coeff_mid(json)
      double precision, intent(in) :: delta_abundance_mid(json)
      double precision, intent(in) :: diffusion_coeff_deriv_mid(json)
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json)
      integer, intent(in) :: npt
      integer :: i

! CORRECT DIFFUSION COEFFICEINTS FOR CHANGE IN X IN THE PREVIOUS ITERATION.
! (microdiff_run PASSES diffusion_coeff_deriv_mid = 0, SO THE CORRECTION
!  IS A NO-OP THERE; gravitational_settling PASSES REAL DERIVATIVES.)
      do i=1,npt-1
         diffusion_coeff_mid(i)=diffusion_coeff_mid(i)+ &
              diffusion_coeff_deriv_mid(i)*delta_abundance_mid(i)
      end do
! NOW RECOMPUTE ELEMENTS OF THE TRIDIAGONAL MATRIX SYSTEM.
      sub_diag(1) = 0.0d0
      diag(1) = 1.0d0 + alpha(1)*diffusion_coeff_mid(1)
      super_diag(1) = -alpha(1)*diffusion_coeff_mid(1)
      do i=2,npt-1
         sub_diag(i) = -alpha(i)*diffusion_coeff_mid(i-1)
         diag(i) = 1.0d0 + alpha(i)*(diffusion_coeff_mid(i-1)+ &
              diffusion_coeff_mid(i))
         super_diag(i) = -alpha(i)*diffusion_coeff_mid(i)
      end do
      sub_diag(npt) = -alpha(npt)*diffusion_coeff_mid(npt-1)
      diag(npt) = 1.0d0 + alpha(npt)*diffusion_coeff_mid(npt-1)
      super_diag(npt) = 0.0d0
      return
end subroutine implicit_diffusion_coeffs

end module diffusion_scheme
