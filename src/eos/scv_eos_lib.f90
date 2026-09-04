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

contains

! Four-point Lagrange sum used by eqscve/eqscvg to interpolate a table
! column along one axis: w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4),
! evaluated left to right exactly as the former inline expressions.
double precision function scv_weighted_sum4(w, y)
      double precision, intent(in) :: w(4), y(4)
      scv_weighted_sum4 = w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4)
end function scv_weighted_sum4

end module scv_eos_lib
