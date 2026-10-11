!----------------------------------------------------------------------
! yale_eos_tables
!----------------------------------------------------------------------
! The state of the Yale EOS: the Fermi-Dirac table of the degenerate
! electron gas (yale: 43 abscissas, 43 degeneracies, the integrals, and
! the binned interval index fermi_table_x_lookup), read by
! read_fermi_table at eos_init; and the Debye-Hueckel correction's
! constants and composition.
module yale_eos_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: cpi
      implicit none

      private
      public :: debye_huckel_coefficient, debye_huckel_nu, debye_huckel_x, &
           debye_huckel_y, debye_huckel_z, debye_huckel_z_total, yale
      type :: yale_eos_state
           real(dp) :: fermi_table_x_grid(43), fermi_table_eta(43), &
                fermi_table_data(5,43,20)
           integer :: fermi_table_x_lookup(261)
      end type yale_eos_state

      type(yale_eos_state) :: yale


! The Debye-Hueckel composition: X, Y, Z and the 18 metal fractions,
! set per zone by eos_set_debye_huckel_composition (elements 1-3) and
! once by eos_set_debye_huckel_z.
      real(dp) :: debye_huckel_z(18)
      real(dp) :: debye_huckel_x, debye_huckel_y, &
           debye_huckel_z_total
! The Debye-Hueckel coefficient -sqrt(pi/(k m_H**3)) e**3/3 (cgs), and
! the charges of the 18 metals it weights (C, N, O, Ne, Na, Mg, Al, Si,
! P, S, Cl, Ar, Ca, Ti, Cr, Mn, Fe, Ni).
      real(dp), parameter :: boltzmann_constant = 1.380658d-16, &
           electron_charge_esu = 4.802d-10, hydrogen_atom_mass = 1.0d0/6.0222137d23
      real(dp), parameter :: debye_huckel_coefficient = -sqrt(cpi/(boltzmann_constant* &
           hydrogen_atom_mass*hydrogen_atom_mass*hydrogen_atom_mass))*electron_charge_esu* &
           electron_charge_esu*electron_charge_esu/3.0d0
      real(dp), parameter :: debye_huckel_nu(18) = [6.0d0, 7.0d0, 8.0d0, 10.0d0, &
           11.0d0, 12.0d0, 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0, &
           20.0d0, 22.0d0, 24.0d0, 25.0d0, 26.0d0, 28.0d0]

!======================================================================
! Sources: the Fermi-Dirac table and the Debye-Hueckel correction of
! the Yale code (Prather 1976; DBG 1992). This state was the COMMON
! blocks ccr (the table) and debhu (the correction: cdh, dhnue, zdh,
! xxdh, yydh, zzdh).

end module yale_eos_tables
