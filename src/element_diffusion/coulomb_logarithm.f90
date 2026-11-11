!----------------------------------------------------------------------
! coulomb_logarithm
!----------------------------------------------------------------------
! The species' concentrations and the pairwise Coulomb logarithms that
! the Thoul et al. (1994) coefficients take, for both pipelines
! (microdiff_coefficients, and original_settling with
! coulomb_log_choice 4). CGS: density in g/cm^3, temperature in K.
!
! The species are the four Thoul columns (species_table): hydrogen,
! helium, the diffused heavy species and the electrons, the electrons
! last. Their concentrations are per electron,
!     c_i = X_i / (A_i sum_j Z_j X_j/A_j)  (ions),  c_e = 1.
! The screening length is the larger of the Debye length and the mean
! interionic distance,
!     lambda_D = (k T / (4 pi e^2 n_e sum_i c_i Z_i^2))^(1/2),
!     a_0      = (3 / (4 pi n_i))^(1/3),
! with n_e = rho/(m_p sum_i A_i c_i) and n_i = n_e sum_ions c_i, and
! the logarithm of each pair is the fit
!     ln Lambda_ij = 0.81245 ln(1 + 0.18769 x_ij^1.2),
!     x_ij = 4 k T lambda / (|Z_i Z_j| e^2),
! which goes as 0.152 x^1.2 for weak screening lengths (x << 1) and as
! 0.97494 ln x - 1.3592 for long ones (x >> 1). The fit is the F77's
! (setup_grsett.f, microdiff_cod.f), which names no source; the
! evaluations of 2026-11 call it Paquette's (Paquette et al. 1986, ApJS
! 61, 177, tabulate the screened collision integrals).
module coulomb_logarithm
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: m_proton_cgs, cc13
      use math_lib
      implicit none
      private
      public :: compute_concentrations, compute_coulomb_logarithms

! 3/(4 pi), to the F77's five figures: a_0 = (3/(4 pi n_i))^(1/3)
      real(dp), parameter :: three_over_four_pi = 0.23873d0
! (k/(4 pi e^2))^(1/2), cm K^-1/2: lambda_D = 6.9010 (T/(n_e sum c Z^2))^(1/2)
      real(dp), parameter :: debye_coefficient = 6.9010d0
! 4 k/e^2, 1/(cm K): x_ij = 2393.9 T lambda/|Z_i Z_j|
      real(dp), parameter :: coulomb_x_coefficient = 2.3939d3
! the fit's amplitude, scale and power
      real(dp), parameter :: coulomb_fit_amplitude = 0.81245d0, &
           coulomb_fit_scale = 0.18769d0, coulomb_fit_power = 1.2d0

contains

! The concentrations per electron of the species (the last one the
! electrons, c = 1) from their mass fractions (the electrons' entry is
! not read), atomic weights (amu) and charges (e).
pure subroutine compute_concentrations(atomic_weight, atomic_charge, mass_fraction, &
     concentration)
      real(dp), intent(in) :: atomic_weight(:), atomic_charge(:), mass_fraction(:)
      real(dp), intent(out) :: concentration(:)
      real(dp) :: charge_per_mass
      integer :: i, num_species

      num_species = size(atomic_weight)
      charge_per_mass = 0.0d0
      do i = 1, num_species-1
         charge_per_mass = charge_per_mass + atomic_charge(i)*mass_fraction(i)/atomic_weight(i)
      end do
      do i = 1, num_species-1
         concentration(i) = mass_fraction(i)/(atomic_weight(i)*charge_per_mass)
      end do
      concentration(num_species) = 1.0d0
end subroutine compute_concentrations

! The Coulomb logarithm of every pair of species at density rho (g/cm^3)
! and temperature t (K), from their concentrations (compute_concentrations),
! atomic weights and charges.
pure subroutine compute_coulomb_logarithms(atomic_weight, atomic_charge, concentration, &
     rho, t, coulomb_log)
      real(dp), intent(in) :: atomic_weight(:), atomic_charge(:), concentration(:)
      real(dp), intent(in) :: rho, t
      real(dp), intent(out) :: coulomb_log(:,:)
      real(dp) :: weight_per_electron, electron_density, ion_density, &
           interion_distance, charge_sq_per_electron, debye_length, &
           screening_length, x
      integer :: i, j, num_species

      num_species = size(atomic_weight)
! the electron density, from the mean weight per electron
      weight_per_electron = 0.0d0
      do i = 1, num_species
         weight_per_electron = weight_per_electron + atomic_weight(i)*concentration(i)
      end do
      electron_density = rho/(m_proton_cgs*weight_per_electron)
! the mean interionic distance
      ion_density = 0.0d0
      do i = 1, num_species-1
         ion_density = ion_density + concentration(i)*electron_density
      end do
      interion_distance = pow((three_over_four_pi/ion_density), cc13)
! the Debye length
      charge_sq_per_electron = 0.0d0
      do i = 1, num_species
         charge_sq_per_electron = charge_sq_per_electron + concentration(i)*atomic_charge(i)**2
      end do
      debye_length = debye_coefficient*sqrt(t/(electron_density*charge_sq_per_electron))
      screening_length = max(debye_length, interion_distance)
      do i = 1, num_species
         do j = 1, num_species
            x = coulomb_x_coefficient*t*screening_length/abs(atomic_charge(i)*atomic_charge(j))
            coulomb_log(i,j) = coulomb_fit_amplitude &
                 *log(1.0d0 + coulomb_fit_scale*pow(x, coulomb_fit_power))
         end do
      end do
end subroutine compute_coulomb_logarithms

end module coulomb_logarithm
