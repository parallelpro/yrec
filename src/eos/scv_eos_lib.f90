!----------------------------------------------------------------------
! scv_eos_lib
!----------------------------------------------------------------------
! Added 2026 (phase six, step 3 -- ROADMAP.md). The SCV
! (Saumon-Chabrier-Van Horn) EOS table state, evicted from const_lib
! where it had landed as former common/scveos/ during the phase-one
! COMMON conversion. The pure-species tables are loaded by eos_lib's
! eos_init; the envelope-mixture table tablenv is built from them by
! eos/scv/scv_envelope_table.f90 (called from core/read_starting_model);
! read by eos/scv/eqscve and eqscvg.
module scv_eos_lib
      implicit none

! former common/scveos/: the SCV EOS tables. tlogx/tablex/tabley/smix/
! tablez/tablenv/nptsx/idtt/idp are spelled identically to their
! canonical names in the majority of files that declare this block
! (eos/scv/eqscve.f90 and eos/scv/eqscvg.f90 used a distinct,
! more-descriptive naming scheme -- table_log10_temperature/
! hydrogen_table/helium_table/entropy_of_mixing_table/metal_table/
! envelope_table/num_pressure_points/scv_temp_index/scv_pressure_index
! -- renamed to the majority spelling here). idtt is itself a rename
! (from idt) to avoid colliding with the unrelated const_lib idt
!
! Table dimensions: scv_nt temperature rows, up to scv_np pressure
! points per row (nptsx holds the actual count), scv_nvar variables per
! point in the H, He and envelope tables and scv_nvar_z in the Z table.
      integer, parameter :: scv_nt = 63
      integer, parameter :: scv_np = 76
      integer, parameter :: scv_nvar = 12
      integer, parameter :: scv_nvar_z = 13
      double precision :: tlogx(scv_nt), tablex(scv_nt,scv_np,scv_nvar), &
           tabley(scv_nt,scv_np,scv_nvar), smix(scv_nt,scv_np), &
           tablez(scv_nt,scv_np,scv_nvar_z), tablenv(scv_nt,scv_np,scv_nvar)
      integer :: nptsx(scv_nt), idtt, idp
      logical :: use_scv_eos

! Readability W3 (2026): named table columns. Only the columns some
! reader or writer in eos/ addresses are named; their meanings follow
! from those uses (eos/scv/eqscvg.f90, eqscve.f90,
! scv_envelope_table.f90). eos_init reads columns 1..11 of the H and
! He tables and 1..13 of the Z table from the files verbatim.
!
! H and He tables (tablex, tabley), same layout:
      integer, parameter :: iscv_log10_p = 1        ! log10 gas pressure (row key)
      integer, parameter :: iscv_frac_neutral = 2   ! particle fraction: H2 (H table), neutral He (He table)
      integer, parameter :: iscv_frac_atom_ion1 = 3 ! particle fraction: neutral H (H table), He+ (He table)
      integer, parameter :: iscv_log10_rho = 4      ! log10 density
      integer, parameter :: iscv_log10_s = 5        ! log10 entropy
      integer, parameter :: iscv_dlnrho_dlnt = 7    ! dln(rho)/dln(T) at constant P
      integer, parameter :: iscv_dlnrho_dlnp = 8    ! dln(rho)/dln(P) at constant T
      integer, parameter :: iscv_dlns_dlnt = 9      ! dln(S)/dln(T) at constant P
      integer, parameter :: iscv_du_dt = 12         ! du/dT, built by scv_envelope_table (not read)
! Z (metal) table (tablez):
      integer, parameter :: iscvz_log10_rho = 4     ! log10 density
      integer, parameter :: iscvz_log10_du_dt = 7   ! log10 du/dT
      integer, parameter :: iscvz_dlnrho_dlnt = 10  ! dln(rho)/dln(T) at constant P
      integer, parameter :: iscvz_dlnrho_dlnp = 13  ! dln(rho)/dln(P) at constant T
! Envelope-mixture table (tablenv), built by scv_envelope_table and
! read by eqscve (columns 2..6 as its five interpolated quantities, in
! this order):
      integer, parameter :: iscvenv_log10_p = 1       ! log10 gas pressure (copied from tablex)
      integer, parameter :: iscvenv_log10_rho = 2     ! log10 mixture density
      integer, parameter :: iscvenv_dlnrho_dlnt = 3   ! dln(rho)/dln(T), gas
      integer, parameter :: iscvenv_dlnrho_dlnp = 4   ! dln(rho)/dln(P), gas
      integer, parameter :: iscvenv_cp = 5            ! cp, gas
      integer, parameter :: iscvenv_du_dt = 6         ! du/dT
      integer, parameter :: iscvenv_dqdt_dlnp = 7     ! d(dlnrho/dlnT)/dln(P)
      integer, parameter :: iscvenv_dqdt_dlnt = 8     ! d(dlnrho/dlnT)/dln(T)
      integer, parameter :: iscvenv_dlnp_dlnrho = 9   ! d(lnP)/d(ln rho) along the row
      integer, parameter :: iscvenv_dlncp_dlnt = 10   ! dln(cp)/dln(T)
      integer, parameter :: iscvenv_dqut_dlnp = 11    ! d(du/dT)/dln(P)
      integer, parameter :: iscvenv_dqut_dlnt = 12    ! d(du/dT)/dln(T)

contains

! Four-point Lagrange sum used by eqscve/eqscvg to interpolate a table
! column along one axis: w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4),
! evaluated left to right exactly as the former inline expressions.
double precision function scv_weighted_sum4(w, y)
      double precision, intent(in) :: w(4), y(4)
      scv_weighted_sum4 = w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4)
end function scv_weighted_sum4

end module scv_eos_lib
