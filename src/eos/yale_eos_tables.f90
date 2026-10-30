!----------------------------------------------------------------------
! yale_eos_tables
!----------------------------------------------------------------------
! The state of the Yale EOS: the Fermi-Dirac table of the degenerate
! electron gas (yale: 43 abscissas, 43 degeneracies, the integrals, and
! the binned interval index fermi_table_x_lookup), read by
! read_fermi_table at eos_init; and the Debye-Hueckel correction's
! constants and its composition record (no state).
module yale_eos_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: cpi, boltzmann_constant_cgs, electron_charge_esu, amu_cgs
      implicit none

      private
      public :: debye_huckel_coefficient, debye_huckel_nu, yale
! the slots of debye_huckel_composition%metals: the metals the
! correction weights
      integer, parameter, public :: idh_c = 1, idh_n = 2, idh_o = 3, idh_ne = 4, &
           idh_na = 5, idh_mg = 6, idh_al = 7, idh_si = 8, idh_p = 9, idh_s = 10, &
           idh_cl = 11, idh_ar = 12, idh_ca = 13, idh_ti = 14, idh_cr = 15, &
           idh_mn = 16, idh_fe = 17, idh_ni = 18, n_debye_huckel_metals = 18
      type :: yale_eos_state
           real(dp) :: fermi_table_x_grid(43), fermi_table_eta(43), &
                fermi_table_data(5,43,20)
           integer :: fermi_table_x_lookup(261)
      end type yale_eos_state

      type(yale_eos_state) :: yale


! The composition the Debye-Hueckel correction weights at a zone: X, Y,
! Z and the mass fractions of its 18 metals (slots idh_*; eos_get builds
! it for each call, eos_lib's debye_huckel_composition_of).
      type, public :: debye_huckel_composition
         real(dp) :: hydrogen_fraction = 0.0d0, helium_fraction = 0.0d0, &
              metal_fraction = 0.0d0
         real(dp) :: metals(n_debye_huckel_metals) = 0.0d0
      end type debye_huckel_composition
! The Debye-Hueckel coefficient -sqrt(pi/(k m_u**3)) e**3/3 (cgs), and
! the charges of the 18 metals it weights (C, N, O, Ne, Na, Mg, Al, Si,
! P, S, Cl, Ar, Ca, Ti, Cr, Mn, Fe, Ni). The mass unit is the atomic
! mass unit, 1/N_A g (the F77 code wrote 1/6.0222137e23, a mistyped
! 6.0221367e23).
      real(dp), parameter :: debye_huckel_mass_unit = amu_cgs
      real(dp), parameter :: debye_huckel_coefficient = -sqrt(cpi/(boltzmann_constant_cgs* &
           debye_huckel_mass_unit*debye_huckel_mass_unit*debye_huckel_mass_unit))* &
           electron_charge_esu*electron_charge_esu*electron_charge_esu/3.0d0
      real(dp), parameter :: debye_huckel_nu(n_debye_huckel_metals) = [6.0d0, 7.0d0, 8.0d0, 10.0d0, &
           11.0d0, 12.0d0, 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0, &
           20.0d0, 22.0d0, 24.0d0, 25.0d0, 26.0d0, 28.0d0]

!======================================================================
! Sources: the Fermi-Dirac table and the Debye-Hueckel correction of
! the Yale code (Prather 1976; DBG 1992). This state was the COMMON
! blocks ccr (the table) and debhu (the correction: cdh, dhnue, zdh,
! xxdh, yydh, zzdh).

end module yale_eos_tables
