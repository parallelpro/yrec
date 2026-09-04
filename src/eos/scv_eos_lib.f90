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
      double precision :: tlogx(63), tablex(63,76,12), tabley(63,76,12), &
           smix(63,76), tablez(63,76,13), tablenv(63,76,12)
      integer :: nptsx(63), idtt, idp
      logical :: use_scv_eos

end module scv_eos_lib
