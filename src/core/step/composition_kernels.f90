!----------------------------------------------------------------------
! composition_kernels
!----------------------------------------------------------------------
! The arithmetic of the composition step (composition_step.f90) that
! does not burn or settle, as star-blind kernels: which shells the
! deep-mixing option mixes, the mass-weighted mean of the mixed regions,
! and the clip and closure of the mass fractions. Each takes the arrays
! it works on as arguments and changes only its intent(inout) and
! intent(out) ones; nothing is read from the star or kept between calls.
! advance_composition gathers the arrays from the star and calls them;
! core/test/test_composition_kernels checks them.
!
! composition is the star's xa(species, shell); shell_mass is dm (g) and
! enclosed_mass m (g), shells 1..num_zones from the centre out.
module composition_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: i_h1, i_he3, i_he4, i_metals
      implicit none
      private
      public :: deep_mixing_flags, homogenize_regions, &
           clip_and_close_abundances

contains

!----------------------------------------------------------------------
! deep_mixing_flags
!----------------------------------------------------------------------
! The mixing flag of each shell under the deep-mixing option (MHP 1/95):
! every shell from the centre out is mixed while its enclosed mass is at
! most max_mixed_mass (g); from the first shell beyond it, each shell
! keeps its ordinary convective flag.
pure subroutine deep_mixing_flags(convective_flag, enclosed_mass, &
     max_mixed_mass, num_zones, mix_flag)
      logical, intent(in) :: convective_flag(:)
      real(dp), intent(in) :: enclosed_mass(:), max_mixed_mass
      integer, intent(in) :: num_zones
      logical, intent(inout) :: mix_flag(:)
      integer :: k

      do k = 1, num_zones
         if (enclosed_mass(k) > max_mixed_mass) then
            mix_flag(k:num_zones) = convective_flag(k:num_zones)
            exit
         else
            mix_flag(k) = .true.
         end if
      end do
end subroutine deep_mixing_flags

!----------------------------------------------------------------------
! homogenize_regions
!----------------------------------------------------------------------
! The mass-weighted mean of each of the first num_species species,
!   X_i = sum(X_i dm)/sum(dm),
! applied throughout each of the num_regions mixed regions
! (region_bounds(r,1) and (r,2): the region's first and last shell), so
! that each species' mass in a region is kept. A one-shell region is
! skipped unless it starts at the centre. (The same mean as convection's
! homogenize_region, which convection_lib does not export.)
pure subroutine homogenize_regions(composition, shell_mass, region_bounds, &
     num_regions, num_species)
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(in) :: shell_mass(:)
      integer, intent(in) :: region_bounds(:,:), num_regions, num_species
      real(dp) :: species_sum(num_species), weight_sum
      integer :: region, first, last, k, i

      do region = 1, num_regions
         first = region_bounds(region,1)
         last = region_bounds(region,2)
         if (first /= 1 .and. first >= last) cycle
         weight_sum = 0.0d0
         species_sum = 0.0d0
         do k = first, last
            weight_sum = weight_sum + shell_mass(k)
            do i = 1, num_species
               species_sum(i) = species_sum(i) + composition(i,k)*shell_mass(k)
            end do
         end do
         do i = 1, num_species
            species_sum(i) = species_sum(i)/weight_sum
         end do
         do k = first, last
            composition(1:num_species,k) = species_sum
         end do
      end do
end subroutine homogenize_regions

!----------------------------------------------------------------------
! clip_and_close_abundances
!----------------------------------------------------------------------
! In shells 1..num_zones: every one of the first num_species abundances
! clipped to [0, 1], then the mass fractions closed,
!   Z <= 1 - X - He3,   Y = 1 - X - Z - He3
! (which leaves Z and Y in [0, 1] as long as the clipped X + He3 is at
! most 1).
pure subroutine clip_and_close_abundances(composition, num_zones, num_species)
      real(dp), intent(inout) :: composition(:,:)
      integer, intent(in) :: num_zones, num_species
      integer :: k, i

      do k = 1, num_zones
         do i = 1, num_species
            composition(i,k) = max(composition(i,k), 0.0d0)
            composition(i,k) = min(composition(i,k), 1.0d0)
         end do
         composition(i_metals,k) = min(composition(i_metals,k), &
              1.0d0 - composition(i_h1,k) - composition(i_he3,k))
         composition(i_he4,k) = 1.0d0 - composition(i_h1,k) - &
              composition(i_metals,k) - composition(i_he3,k)
      end do
end subroutine clip_and_close_abundances

end module composition_kernels
