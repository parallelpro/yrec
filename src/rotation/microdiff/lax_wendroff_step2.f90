!----------------------------------------------------------------------
! lax_wendrof2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lax_wendrof2.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0
! regression suite (examples/run_standard_solar_model).
!
! Second (corrector) step of the two-step Lax-Wendroff diffusion
! method. Using the zone-midpoint diffusion coefficients from
! lax_wendrof1's provisional solution, updates diffused_abundance in
! place with the full-step change in the diffused species' abundance.
! MHP 3/94 added metal diffusion.
! Dummy-argument names match the actual arguments used at the call
! site in microdiff_run.f90.
subroutine lax_wendroff_step2(timestep, diffusion_coeff1_mid, eq_mass_mid, &
     diffused_abundance, num_eq_points, total_mass, &
     use_generic_diffusion_vectors)

      use star_info_lib, only: star, json
      use phys_const_lib
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
!    99 format(5x,1p2e15.7)
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
         delta_metal_abundance = dt_full*star%rot%metal_diffusion_coeff1_mid(1)/zone_mass
         star%rot%metal_abundance_change(1) = star%rot%metal_abundance_change(1)+delta_metal_abundance
! general case
         do i = 2,num_eq_points-1
            zone_mass = eq_mass_mid(i)-eq_mass_mid(i-1)
            delta_metal_abundance = dt_full*(star%rot%metal_diffusion_coeff1_mid(i)-star%rot%metal_diffusion_coeff1_mid(i-1))/zone_mass
            star%rot%metal_abundance_change(i) = star%rot%metal_abundance_change(i)+delta_metal_abundance
         end do
! surface boundary condition.
         zone_mass = total_mass-eq_mass_mid(num_eq_points-1)
         delta_metal_abundance = -dt_full*star%rot%metal_diffusion_coeff1_mid(num_eq_points-1)/zone_mass
         star%rot%metal_abundance_change(num_eq_points) = star%rot%metal_abundance_change(num_eq_points)+delta_metal_abundance
      endif
      return
end subroutine lax_wendroff_step2
