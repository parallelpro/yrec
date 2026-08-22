!----------------------------------------------------------------------
! yale_eos_lib
!----------------------------------------------------------------------
! New (2026, phase three -- ROADMAP.md stage 1). Holds the Yale/
! Prather EOS's tabulated state: the degenerate-electron (Fermi-Dirac)
! table, former common/ccr/, read once at startup (now by eos_lib's
! eos_init; historically by an inline block in setup/setups.f90) and
! consumed only by eos/yale/eqrelv.f90's fully-ionized-gas solve.
!
! These members lived in atm_table_lib through phases one and two
! purely as an accident of the original COMMON grouping -- eos-domain
! state in atm's state module, the same misplacement pattern as the
! OPAL-2006 common/eeos06/ members relocated to opal_eos_lib earlier
! in stage 1. Member names are unchanged from their atm_table_lib
! spelling (fermi_table_*); only the instance prefix changed
! (atm_table% -> yale_eos%).
!
! fermi_table_x_grid/fermi_table_eta: the 43-point table abscissas.
! fermi_table_data(5,43,20): the tabulated Fermi-Dirac integral
! solutions. fermi_table_x_lookup(261): the binned index built at
! load time mapping a scaled abscissa to its table interval (see
! eos_init's bin-construction loop).
module yale_eos_lib
      implicit none

      type, public :: yale_eos_state
           double precision :: fermi_table_x_grid(43), fermi_table_eta(43), &
                fermi_table_data(5,43,20)
           integer :: fermi_table_x_lookup(261)
      end type yale_eos_state

      type(yale_eos_state), public :: yale_eos

end module yale_eos_lib
