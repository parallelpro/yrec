!----------------------------------------------------------------------
! screening
!----------------------------------------------------------------------
! The electron screening of the reaction rates: the composition the
! screening sees (the isotopes iso_*, their bare nuclear masses and
! charges, screening_mass_fractions), the plasma quantities of a zone
! (compute_composition_moments, compute_screening_plasma: Salpeter's
! degeneracy correction, Graboske et al.'s Lambda_0, DeWitt et al.'s
! z_tilde, the mean ionic charge), and the screening exponent of one
! reaction with its derivatives (compute_reaction_screening: weak,
! intermediate or strong; weak_intermediate_screening for Be7 + p and
! hep), from its constants (screening_t). reaction_rates applies it.
! Sources: the end of reaction_rates.f90 (Screening).
module screening
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use net_configuration, only: net_settings
      implicit none
      private

! The isotopes of the screening composition and of the mass and charge
! tables.
      integer, parameter, public :: &
           iso_n = 1, iso_h1 = 2, iso_h2 = 3, iso_h3 = 4, iso_he3 = 5, &
           iso_he4 = 6, iso_c12 = 7, iso_c13 = 8, iso_n14 = 9, &
           iso_o16 = 10, iso_o18 = 11, iso_ne20 = 12, iso_mg24 = 13, &
           num_isotopes = 13

! bare nuclear masses (amu) and charges of the iso_* isotopes
      real(dp), parameter, public :: &
           isotope_mass(num_isotopes) = [1.008665d0, 1.007276d0, 2.013553d0, 3.015501d0, 3.014933d0, 4.001506d0, &
                11.996709d0, 13.000064d0, 13.999233d0, 15.990526d0, 17.994772d0, 19.986954d0, 23.978458d0], &
           isotope_charge(num_isotopes) = [0.d0, 1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 6.d0, 6.d0, 7.d0, 8.d0, 8.d0, 10.d0, 12.d0]

! A reaction's screening constants: the charge product Z1 Z2 and the
! Graboske et al. (1973, table 4) factors of the intermediate (z86) and
! strong (z53, z43, z23) forms.
      type, public :: screening_t
         real(dp) :: zprod = 0.0d0, z86 = 0.0d0, z53 = 0.0d0, z43 = 0.0d0, &
              z23 = 0.0d0
      end type screening_t

! The screening plasma's constants (compute_screening_plasma):
!   (p_F/m_e c)^2 per (rho/mu_e)^(2/3), rho/mu_e in g/cm^3:
!     (3 pi^2 N_A)^(2/3) (hbar/m_e c)^2;
!   m_e c^2/k in 1e9 K, for E_F/kT;
!   Graboske et al.'s (1973) Lambda_0 = 1.88e8 (rho/mu_i)^(1/2) T^(-3/2)
!     with T in 1e9 K (1.88e8/(1e9)^1.5);
!   E_F/kT at or below which the electrons count as nondegenerate
!     (Salpeter's f'/f = 1), and log10(E_F/kT) from which they are fully
!     degenerate (f'/f = 0).
      real(dp), parameter :: fermi_momentum_sq_coeff = 1.017677d-4, &
           electron_rest_energy_t9 = 5.92986d0, lambda0_coeff = 5.9426d-6, &
           nondegenerate_ef_over_kt = 1.0d-2, degenerate_log_ef_over_kt = 1.5d0
! The weak-screening exponent Lambda_0 z_tilde Z1 Z2 up to which the
! intermediate form applies (above the weak threshold of the settings),
! and from which the strong form always does (compute_reaction_screening).
      real(dp), parameter :: max_intermediate_weak_exponent = 2.0d0, &
           min_strong_weak_exponent = 5.0d0

! The plasma quantities the screening formulas share: the composition
! moments (1/mu_i, 1/mu_e, sum Z^1.58 X/A, sum Z^2 X/A), Graboske et
! al.'s Lambda_0 and its powers, DeWitt et al.'s z_tilde, the mean ionic
! charge and its powers, and the weak-screening exponent per unit charge
! product.
      type :: screening_plasma
         real(dp) :: inv_mu_i = 0.0d0, inv_mu_e = 0.0d0, &
              charge_moment_158 = 0.0d0, charge_moment_2 = 0.0d0
         real(dp) :: lambda0 = 0.0d0, lambda0_23 = 0.0d0, &
              lambda0_086 = 0.0d0, z_tilde = 0.0d0, mean_charge = 0.0d0, &
              z_tilde_058 = 0.0d0, mean_charge_028 = 0.0d0, &
              mean_charge_13 = 0.0d0, weak_unit = 0.0d0
      end type screening_plasma

      public :: screening_plasma, screening_mass_fractions, &
           compute_composition_moments, compute_screening_plasma, &
           compute_reaction_screening, weak_intermediate_screening

contains

!----------------------------------------------------------------------
! The composition moments the screening needs: 1/mu_i = sum X/A,
! 1/mu_e = sum Z X/A, sum Z^1.58 X/A (intermediate screening, Graboske et
! al.: the mean of Z^(3b-1)) and sum Z^2 X/A (Salpeter's zeta).
pure subroutine compute_composition_moments(mass_frac, plasma)
      use math_lib
      real(dp), intent(in) :: mass_frac(num_isotopes)
      type(screening_plasma), intent(out) :: plasma
      real(dp) :: per_amu
      integer :: i

      plasma%inv_mu_i = 0.d0
      plasma%inv_mu_e = 0.d0
      plasma%charge_moment_158 = 0.d0
      plasma%charge_moment_2 = 0.d0
      do i = 1,num_isotopes
         per_amu = mass_frac(i)/isotope_mass(i)
         plasma%inv_mu_i = plasma%inv_mu_i+per_amu
         plasma%inv_mu_e = plasma%inv_mu_e+ &
              per_amu*isotope_charge(i)
         plasma%charge_moment_158 = plasma%charge_moment_158+per_amu*pow(isotope_charge(i), 1.58d0)
         plasma%charge_moment_2 = plasma%charge_moment_2+per_amu*isotope_charge(i)**2
      end do
end subroutine compute_composition_moments

!----------------------------------------------------------------------
! The plasma quantities of a zone's screening: Salpeter's degeneracy
! correction f'/f, Graboske et al.'s Lambda_0, DeWitt et al.'s z_tilde,
! the mean ionic charge, and the powers of them the screening forms use
! (compute_reaction_screening, weak_intermediate_screening).
pure subroutine compute_screening_plasma(density, electron_density_na, &
     t9_m1, t9_m32, plasma)
      use physical_constants, only: cc13
      use math_lib
! t9_m1, t9_m32: T9^-1 and T9^-3/2
      real(dp), intent(in) :: density, electron_density_na, t9_m1, t9_m32
      type(screening_plasma), intent(inout) :: plasma
! (p_F/m_e c)^2, E_F/kT, Salpeter's f'/f and log10(E_F/kT)
      real(dp) :: fermi_momentum_sq, fermi_energy_over_kt, &
           degeneracy_factor, log_degeneracy

      fermi_momentum_sq=fermi_momentum_sq_coeff*pow(electron_density_na, (2.0d0/3.0d0))
      fermi_energy_over_kt=electron_rest_energy_t9*t9_m1*(sqrt(1.d0+fermi_momentum_sq)-1.d0)
      if (fermi_energy_over_kt <= nondegenerate_ef_over_kt) then
         degeneracy_factor=1.0d0
      else
! Salpeter's figure 1, as a function of the degeneracy
         log_degeneracy=log10(fermi_energy_over_kt)
         if (log_degeneracy >= degenerate_log_ef_over_kt) then
            degeneracy_factor=0.0d0
         else
            degeneracy_factor=0.75793d0-0.54621d0*log_degeneracy-0.30964d0*log_degeneracy**2+ &
            0.12535d0*log_degeneracy**3+0.1203d0*log_degeneracy**4- &
            0.012857d0*log_degeneracy**5-0.014768d0*log_degeneracy**6
         end if
      end if
! Lambda_0 of Graboske et al. and its 2/3 and 0.86 powers; z_tilde of
! DeWitt et al. (Salpeter's zeta without the 1/mu_i), the mean ionic
! charge, and the powers of them the intermediate and strong forms use
      plasma%lambda0=lambda0_coeff*t9_m32*sqrt(density*plasma%inv_mu_i)
      plasma%lambda0_23=pow(plasma%lambda0, (2.0d0/3.0d0))
      plasma%lambda0_086=pow(plasma%lambda0, 0.86d0)
      plasma%z_tilde=sqrt((plasma%charge_moment_2+degeneracy_factor*plasma%inv_mu_e)/ &
           plasma%inv_mu_i)
      plasma%mean_charge=plasma%inv_mu_e/plasma%inv_mu_i
      plasma%z_tilde_058=pow(plasma%z_tilde, 0.58d0)
      plasma%mean_charge_028=pow(plasma%mean_charge, 0.28d0)
      plasma%mean_charge_13=pow(plasma%mean_charge, cc13)
! weak screening per unit charge product: Z1 Z2 Lambda_0 z_tilde
      plasma%weak_unit=plasma%lambda0*plasma%z_tilde
end subroutine compute_screening_plasma

!----------------------------------------------------------------------
! The screening exponent of one reaction and its derivatives: weak
! (Salpeter 1954, with his degeneracy correction) below
! net_settings%weak_screening_threshold (0.03 reproduces Graboske et al.
! and Salpeter; a value above 1 forces weak screening in the Sun),
! intermediate (Graboske et al. 1973, table 4) up to a weak-screening
! value of max_intermediate_weak_exponent (2), and beyond that strong or
! intermediate, whichever is smaller (strong always from
! min_strong_weak_exponent, 5). (Strong screening was never checked:
! it does not occur in the Sun.)
pure subroutine compute_reaction_screening(plasma, s, screen, &
     dscreen_dlnrho, dscreen_dlnt)
      type(screening_plasma), intent(in) :: plasma
      type(screening_t), intent(in) :: s
      real(dp), intent(out) :: screen, dscreen_dlnrho, dscreen_dlnt
      real(dp) :: weak, intermediate, strong

      weak=plasma%weak_unit*s%zprod
      if (weak <= net_settings%weak_screening_threshold) then
         screen=weak
         dscreen_dlnrho=0.5d0*weak
         dscreen_dlnt=-1.5d0*weak
      else
         intermediate=0.38d0*plasma%lambda0_086*plasma%charge_moment_158*s%z86/ &
              (plasma%inv_mu_i*plasma%z_tilde_058*plasma%mean_charge_028)
         if (weak <= max_intermediate_weak_exponent) then
            screen=intermediate
            dscreen_dlnrho=0.43d0*intermediate
            dscreen_dlnt=-1.29d0*intermediate
         else
            strong=0.624d0*plasma%mean_charge_13*plasma%lambda0_23*(s%z53+ &
                 0.316d0*plasma%mean_charge_13*s%z43+0.737d0* &
                 s%z23/(plasma%mean_charge*plasma%lambda0_23))
            if (strong < intermediate.or.weak >= min_strong_weak_exponent) then
               screen=strong
               dscreen_dlnrho=0.208d0*plasma%mean_charge_13*(s%z53+ &
                    0.316d0*plasma%mean_charge_13*s%z43)*plasma%lambda0_23
               dscreen_dlnt=-3.d0*dscreen_dlnrho
            else
               screen=intermediate
               dscreen_dlnrho=0.43d0*intermediate
               dscreen_dlnt=-1.29d0*intermediate
            end if
         end if
      end if
end subroutine compute_reaction_screening

!----------------------------------------------------------------------
! Weak screening, or intermediate above the weak threshold, for one
! reaction of charge product zprod and Graboske z86 (Be7 + p, hep).
pure real(dp) function weak_intermediate_screening(plasma, zprod, z86)
      type(screening_plasma), intent(in) :: plasma
      real(dp), intent(in) :: zprod, z86
      real(dp) :: weak

      weak = plasma%weak_unit*zprod
      if (weak <= net_settings%weak_screening_threshold) then
         weak_intermediate_screening = weak
      else
         weak_intermediate_screening = 0.38d0*plasma%lambda0_086*plasma%charge_moment_158*z86/ &
              (plasma%inv_mu_i*plasma%z_tilde_058*plasma%mean_charge_028)
      end if
end function weak_intermediate_screening

!----------------------------------------------------------------------
! screening_mass_fractions: the isotope mass fractions (iso_* slots) that
! compute_reaction_rates screens with, from a zone's composition in the
! star%xa layout. The isotopes YREC does not track (n, H3, Ne20, Mg24)
! are zero; deuterium is included (2026-10: the composition path used to
! leave it out).
pure function screening_mass_fractions(xa) result(mass_frac)
      use composition_layout, only: i_h1, i_h2, i_he4, i_he3, i_c12, i_c13, &
           i_n14, i_o16, i_o18
      real(dp), intent(in) :: xa(:)
      real(dp) :: mass_frac(num_isotopes)

      mass_frac = 0.0d0
      mass_frac(iso_h1) = xa(i_h1)
      mass_frac(iso_h2) = xa(i_h2)
      mass_frac(iso_he3) = xa(i_he3)
      mass_frac(iso_he4) = xa(i_he4)
      mass_frac(iso_c12) = xa(i_c12)
      mass_frac(iso_c13) = xa(i_c13)
      mass_frac(iso_n14) = xa(i_n14)
      mass_frac(iso_o16) = xa(i_o16)
      mass_frac(iso_o18) = xa(i_o18)
end function screening_mass_fractions

end module screening
