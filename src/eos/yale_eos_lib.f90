!----------------------------------------------------------------------
! yale_eos_lib
!----------------------------------------------------------------------
! New (2026, phase three -- ROADMAP.md stage 1). Holds the Yale/
! Prather EOS's tabulated state: the degenerate-electron (Fermi-Dirac)
! table, former common/ccr/, read once at startup (now by eos_lib's
! eos_init; historically by an inline block in setup/setups.f90) and
! consumed only by eos/yale/fully_ionized_eos.f90's fully-ionized-gas solve.
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


! 2026 (phase six, step 3 -- ROADMAP.md): evicted here from
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
      double precision :: debye_huckel_coefficient
      double precision :: debye_huckel_eta_min, debye_huckel_eta_max
      double precision :: debye_huckel_z(18)
      double precision :: debye_huckel_x, debye_huckel_y, &
           debye_huckel_z_total
      double precision :: debye_huckel_nu(18)
      logical :: use_debye_huckel_correction


end module yale_eos_lib
