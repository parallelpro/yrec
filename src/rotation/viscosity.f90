!----------------------------------------------------------------------
! viscosity
!----------------------------------------------------------------------
! The microscopic kinematic viscosity (radiative plus ionic) and the
! radiative thermal diffusivity of each zone (viscos).
module viscosity
      implicit none
      private
      public :: viscos

contains

! ---------------------------------------------------------------
! viscos: for every zone, convective ones included, the kinematic
! viscosity visc = nu_rad + nu_ion and the radiative thermal
! diffusivity thdif = K (both cm^2/s), for the rotational instability
! and mixing coefficients, from the composition, log rho, log T, the
! opacity (cm^2/g) and c_P (erg/g/K) of each zone.
subroutine viscos(composition, log_density, log_temperature, opacity, cp, &
     num_zones, visc, thdif)
      use star_info_lib, only: json, i_metals
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_temperature(json), opacity(json), cp(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: visc(json), thdif(json)

      double precision :: amu
! Atomic weights and charges for composition rows 1..11 (i_h1..i_o18 of
! star_info_lib). Row i_metals is the lumped Z (weight 1, charge 0): it
! is skipped in the ion-viscosity sums below, but it does enter the
! mean charge, as a neutral particle. These values differ from the
! Thoul-solver descriptors in species_table and stay local.
      double precision :: weight(11), z(11)
      data amu/amu_cgs_legacy/
      data weight/1.007825d0,4.0026d0,1.0d0,3.01603d0,12.0d0, &
           13.00335d0,14.00307d0,15.0011d0,15.99491d0,16.99491d0, &
           17.99491d0/
      data z/1.0d0,2.0d0,0.0d0,2.0d0,6.0d0,6.0d0,7.0d0,7.0d0,8.0d0, &
           8.0d0,8.0d0/
! --- locals ---
      integer :: shell_idx, species_idx, species_idx2
      double precision :: opacity_local, mean_charge, number_density_sum
      double precision :: temperature_cgs, temperature_sq, density_cgs, &
           electron_number_density
      double precision :: viscosity_radiative
      double precision :: number_density(11)
      double precision :: viscosity_molecular_species(11)
      double precision :: coulomb_log_factor
      double precision :: mfp_temperature_factor, molecular_coeff, &
           viscosity_molecular
      double precision :: species_coeff, species_sum

      do shell_idx = 1,num_zones
! Abundances by number, n_i = X_i/A_i (per unit mass, in units of
! 1/m_u), the mean ion charge Z_bar = sum n_i Z_i / sum n_i, and
! electron_number_density = Z_bar rho/m_u. (The electron density is
! Z_bar (rho/m_u) sum n_i; the two agree only when sum n_i = 1.)
         opacity_local = opacity(shell_idx)
         mean_charge = 0.0d0
         number_density_sum = 0.0d0
         do species_idx = 1,11
            number_density(species_idx) = composition(species_idx,shell_idx)/ &
                 weight(species_idx)
            number_density_sum = number_density_sum+number_density(species_idx)
            mean_charge = mean_charge+number_density(species_idx)*z(species_idx)
         end do
         mean_charge = mean_charge/number_density_sum
         temperature_cgs = exp(ln10*log_temperature(shell_idx))
         temperature_sq = temperature_cgs**2
         density_cgs = exp(ln10*log_density(shell_idx))
         electron_number_density = mean_charge*density_cgs/amu
! The radiative viscosity, after Ledoux (1958, Handbuch der Physik,
! vol. 51, p. 445):
!   nu_rad = 4 a T^4 / (15 c kappa rho^2)   (kinematic),
! with the zone's opacity kappa (star%opacity_zone). The alternative
! forms of Endal & Sofia and of Thomas (1930) are not computed.
         viscosity_radiative = rad_viscosity_coeff_cgs*temperature_sq*temperature_sq/ &
              (opacity_local*density_cgs**2)
! The ion (molecular) viscosity, after Spitzer (1962, Physics of Fully
! Ionized Gases). The dynamic viscosity of species i is
!   eta_i = 3.125e-15 T^(5/2) n_i sqrt(A_i)
!           / (Z_i^2 lnL_i sum_j n_j Z_j^2 sqrt((A_i + A_j)/A_j))
! (for one species, 2.21e-15 T^(5/2) sqrt(A) / (Z^4 lnL) g/(cm s)),
! with the Coulomb logarithm
!   lnL_i = 9.4245 + 0.5 ln(T^3 min(1, 4.2e5/T) / (n_e Z_bar^2)) - ln Z_i,
! where min(1, 4.2e5/T) is the quantum correction above 4.2e5 K.
         mfp_temperature_factor = dmin1(1.0d0,4.2d5/temperature_cgs)
         coulomb_log_factor = 9.424536845d0+0.5d0*log(temperature_sq* &
              temperature_cgs*mfp_temperature_factor/ &
              (electron_number_density*mean_charge**2))
         molecular_coeff = 3.125d-15*dsqrt(temperature_cgs)*temperature_sq
         viscosity_molecular = 0.0d0
! viscosity_molecular_species(i) is eta_i. A species without collision
! partners, or with eta_i <= 0 (lnL_i <= 0), adds nothing.
         do species_idx = 1,11
            if(species_idx.eq.i_metals) cycle
            species_coeff = molecular_coeff*number_density(species_idx)* &
                 dsqrt(weight(species_idx))/ &
                 ((coulomb_log_factor-log(z(species_idx)))*z(species_idx)**2)
            species_sum = 0.0d0
            do species_idx2 = 1,11
               if(species_idx2.eq.i_metals) cycle
               species_sum = species_sum+number_density(species_idx2)* &
                    z(species_idx2)**2* &
                    dsqrt((weight(species_idx)+weight(species_idx2))/ &
                    weight(species_idx2))
            end do
            if(dabs(species_sum).lt.1.0d-38) cycle
            viscosity_molecular_species(species_idx) = species_coeff/species_sum
            if(viscosity_molecular_species(species_idx).gt.0.0d0) &
                 viscosity_molecular = viscosity_molecular+ &
                 viscosity_molecular_species(species_idx)
         end do
! nu_ion = sum_i eta_i / rho.
         viscosity_molecular = viscosity_molecular/density_cgs
         visc(shell_idx) = viscosity_radiative+viscosity_molecular
! The radiative thermal diffusivity (conduction neglected):
!   K = 16 sigma T^3 / (3 kappa rho^2 c_P).
         thdif(shell_idx) = 1.6d1*cc13*sigma_sb_cgs_legacy*temperature_cgs* &
              temperature_sq/(opacity_local*density_cgs**2*cp(shell_idx))
      end do

      return
end subroutine viscos

end module viscosity
