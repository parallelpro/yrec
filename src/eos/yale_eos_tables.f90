!----------------------------------------------------------------------
! yale_eos_tables
!----------------------------------------------------------------------
! New (2026, phase three -- docs/ROADMAP.md stage 1). Holds the Yale/
! Prather EOS's tabulated state: the degenerate-electron (Fermi-Dirac)
! table, former common/ccr/, read once at startup (now by eos_lib's
! eos_init; historically by an inline block in setup/setups.f90) and
! consumed only by compute_fully_ionized_eos's fully-ionized-gas solve.
!
! These members lived in atm_tables through phases one and two
! purely as an accident of the original COMMON grouping -- eos-domain
! state in atm's state module, the same misplacement pattern as the
! OPAL-2006 common/eeos06/ members relocated to opal_eos_tables earlier
! in stage 1. Member names are unchanged from their atm_tables
! spelling (fermi_table_*); only the instance prefix changed
! (atm_table% -> yale%).
!
! fermi_table_x_grid/fermi_table_eta: the 43-point table abscissas.
! fermi_table_data(5,43,20): the tabulated Fermi-Dirac integral
! solutions. fermi_table_x_lookup(261): the binned index built at
! load time mapping a scaled abscissa to its table interval (see
! eos_init's bin-construction loop).
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


! 2026 (phase six, step 3 -- docs/ROADMAP.md): evicted here from
! const_lib, where this table/working data had landed during the
! phase-one COMMON conversion; it belongs with this domain's state.
! former common/debhu/: Debye-Huckel EOS correction data. No readable
! rename had ever been established for this block (every file used
! the original cryptic COMMON member spelling, unlike most other
! blocks), so these canonical names are new (2026), not picked up from
! an existing convention. use_debye_huckel_correction/
! debye_huckel_eta_min/debye_huckel_eta_max (originally ldh/etadh0/
! etadh1) are NAMELIST /physics/ values, kept local in
! core/read_input.f90 and copy-assigned. debye_huckel_coefficient/
! debye_huckel_nu (originally cdh/dhnue) are not namelist values --
! setup/setups.f90 computes them once at startup -- so they have no
! declaration-time default. debye_huckel_x/debye_huckel_y/
! debye_huckel_z_total/debye_huckel_z (originally xxdh(or xxdy)/yydh/
! zzdh/zdh) are likewise not namelist values -- mixing/compute_scale_height.f90
! recomputes them per shell from the local composition -- so they too
! have no declaration-time default.
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

end module yale_eos_tables
