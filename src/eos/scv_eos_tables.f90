!----------------------------------------------------------------------
! scv_eos_tables
!----------------------------------------------------------------------
! The SCV EOS tables: pure hydrogen (tablex), pure helium (tabley) and
! pure metals (tablez), read by read_scv_tables at eos_init; the
! entropy of mixing (smix) and the envelope-mixture table (tablenv),
! built from them by build_scv_envelope_table; and the named columns
! of each.
module scv_eos_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none

      private
      public :: idp, idtt, iscv_dlnrho_dlnp, iscv_dlnrho_dlnt, &
           iscv_dlns_dlnt, iscv_frac_atom_ion1, iscv_frac_neutral, &
           iscv_log10_p, iscv_log10_rho, iscv_log10_s, iscvenv_cp, &
           iscvenv_dlnrho_dlnp, iscvenv_dlnrho_dlnt, iscvenv_du_dt, &
           iscvenv_log10_p, iscvenv_log10_rho, iscvz_dlnrho_dlnp, &
           iscvz_dlnrho_dlnt, iscvz_log10_du_dt, iscvz_log10_rho, nptsx, &
           scv_nt, smix, tablenv, tablex, tabley, tablez, &
           tlogx, tablenv_hydrogen_fraction, tablenv_metal_fraction
! Table dimensions: scv_nt temperature rows, up to scv_np pressure
! points per row (nptsx holds the actual count), scv_nvar variables per
! point in the H and He tables, scv_nvar_z in the Z table and
! scv_nvar_env in the envelope table.
      integer, parameter :: scv_nt = 63
      integer, parameter :: scv_np = 76
      integer, parameter :: scv_nvar = 11
      integer, parameter :: scv_nvar_z = 13
      integer, parameter :: scv_nvar_env = 6
      real(dp) :: tlogx(scv_nt), tablex(scv_nt,scv_np,scv_nvar), &
           tabley(scv_nt,scv_np,scv_nvar), smix(scv_nt,scv_np), &
           tablez(scv_nt,scv_np,scv_nvar_z), tablenv(scv_nt,scv_np,scv_nvar_env)
! the cell (first row, first column) the last SCV lookup used, where the
! next one starts its search
      integer :: nptsx(scv_nt), idtt, idp
! the envelope X and Z tablenv was built for (-1 until it is built)
      real(dp) :: tablenv_hydrogen_fraction = -1.0d0, tablenv_metal_fraction = -1.0d0
! Named table columns: only the columns some routine addresses are
! named. read_scv_tables reads columns 1-11 of the H and He tables and
! 1-13 of the Z table.
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
! Z (metal) table (tablez):
      integer, parameter :: iscvz_log10_rho = 4     ! log10 density
      integer, parameter :: iscvz_log10_du_dt = 7   ! log10 du/dT
      integer, parameter :: iscvz_dlnrho_dlnt = 10  ! dln(rho)/dln(T) at constant P
      integer, parameter :: iscvz_dlnrho_dlnp = 13  ! dln(rho)/dln(P) at constant T
! Envelope-mixture table (tablenv), built by build_scv_envelope_table and
! read by compute_scv_envelope_eos (columns 2..6 as its five interpolated quantities, in
! this order):
      integer, parameter :: iscvenv_log10_p = 1       ! log10 gas pressure (copied from tablex)
      integer, parameter :: iscvenv_log10_rho = 2     ! log10 mixture density
      integer, parameter :: iscvenv_dlnrho_dlnt = 3   ! dln(rho)/dln(T), gas
      integer, parameter :: iscvenv_dlnrho_dlnp = 4   ! dln(rho)/dln(P), gas
      integer, parameter :: iscvenv_cp = 5            ! cp, gas
      integer, parameter :: iscvenv_du_dt = 6         ! du/dT


!======================================================================
! Sources: Saumon, Chabrier & Van Horn (1995); the metal table was
! computed with the Prather EOS of the old Yale code. This state was
! the COMMON block scveos.

end module scv_eos_tables
