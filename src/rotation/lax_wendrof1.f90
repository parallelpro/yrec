!----------------------------------------------------------------------
! lax_wendrof1
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lax_wendrof1.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0
! regression suite (examples/run_standard_solar_model).
!
! First step of the two-step Lax-Wendroff diffusion method. Computes
! new abundances at zone midpoints using the Lax scheme:
! X(n+1/2,j+1/2) = 1/2(X(n,j+1/2) - 1/2(dt/dr)(cod1(n,j+1)-cod1(n,j))
! where n is the time variable, j is the spatial one, and cod1 is the
! diffusion coefficient.
! Variables: diffusion_coeff1 is the run of diffusion coefficients,
! eq_mass is the run of masses, and diffused_abundance_mid is the run
! of changes in the diffused species' abundance at the zone midpoints.
! MHP 3/94 added metal diffusion.
! Dummy-argument names match the actual arguments used at the call
! site in microdiff_run.f90.
subroutine lax_wendrof1(timestep, diffusion_coeff1, eq_mass, num_eq_points, &
     total_mass, diffused_abundance_mid, use_generic_diffusion_vectors)

      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: timestep
      double precision, intent(in) :: diffusion_coeff1(json), eq_mass(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: total_mass
      double precision, intent(out) :: diffused_abundance_mid(json)
      logical, intent(in) :: use_generic_diffusion_vectors

! common/gravez/: metal (Z) diffusion work arrays, parallel to the
! diffusion_coeff1/diffused_abundance_mid dummy arguments above but
! for the metal species specifically. Only metal_diffusion_coeff1 (in)
! and metal_abundance_change_mid (out) are used here. The "_mid"
! members mirror the "_mid" convention used for equally-spaced-grid
! quantities in microdiff_run.f90; the "eq_" prefixed pair's exact
! role is not exercised by this file and is preserved only for
! common-block layout (naming is a conservative guess, not confirmed
! usage -- possibly the equally-spaced-grid analogue of the metal
! diffusion coefficients, by analogy with eq_mass/eq_radius/etc. in
! microdiff_run.f90).
      double precision :: metal_diffusion_coeff1(json), &
           metal_diffusion_coeff1_mid(json), metal_diffusion_coeff2_mid(json), &
           eq_metal_diffusion_coeff1_mid(json), &
           eq_metal_diffusion_coeff2_mid(json), metal_abundance_change(json), &
           metal_abundance_change_mid(json)
      common/gravez/ metal_diffusion_coeff1, metal_diffusion_coeff1_mid, &
           metal_diffusion_coeff2_mid, eq_metal_diffusion_coeff1_mid, &
           eq_metal_diffusion_coeff2_mid, metal_abundance_change, &
           metal_abundance_change_mid
      save

      double precision :: dt_half, zone_mass, delta_abundance_mid, &
           delta_metal_abundance_mid
      integer :: i

      dt_half = 0.5d0*timestep*c4pi
! central boundary condition
      zone_mass = eq_mass(2)
      delta_abundance_mid = dt_half*diffusion_coeff1(2)/zone_mass
      diffused_abundance_mid(1) = delta_abundance_mid
! general case
      do 10 i = 2,num_eq_points-2
         zone_mass = eq_mass(i+1)-eq_mass(i)
         delta_abundance_mid = dt_half*(diffusion_coeff1(i+1)-diffusion_coeff1(i))/zone_mass
         diffused_abundance_mid(i) = delta_abundance_mid
   10 continue
! surface boundary condition.
      zone_mass = total_mass-eq_mass(num_eq_points-1)
      delta_abundance_mid = -dt_half*diffusion_coeff1(num_eq_points-1)/zone_mass
      diffused_abundance_mid(num_eq_points-1) = delta_abundance_mid
! metal diffusion
      if(use_diffusion_z.and..not.use_generic_diffusion_vectors)then
! central boundary condition
         zone_mass = eq_mass(2)
         delta_metal_abundance_mid = dt_half*metal_diffusion_coeff1(2)/zone_mass
         metal_abundance_change_mid(1) = delta_metal_abundance_mid
! general case
         do i = 2,num_eq_points-2
            zone_mass = eq_mass(i+1)-eq_mass(i)
            delta_metal_abundance_mid = dt_half*(metal_diffusion_coeff1(i+1)-metal_diffusion_coeff1(i))/zone_mass
            metal_abundance_change_mid(i) = delta_metal_abundance_mid
         end do
! surface boundary condition.
         zone_mass = total_mass-eq_mass(num_eq_points-1)
         delta_metal_abundance_mid = -dt_half*metal_diffusion_coeff1(num_eq_points-1)/zone_mass
         metal_abundance_change_mid(num_eq_points-1) = delta_metal_abundance_mid
      endif
      return
end subroutine lax_wendrof1
