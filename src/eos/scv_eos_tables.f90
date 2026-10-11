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
           iscv_dlns_dlnt, iscv_du_dt, iscv_frac_atom_ion1, iscv_frac_neutral, &
           iscv_log10_p, iscv_log10_rho, iscv_log10_s, iscvenv_cp, &
           iscvenv_dlncp_dlnt, iscvenv_dlnp_dlnrho, iscvenv_dlnrho_dlnp, &
           iscvenv_dlnrho_dlnt, iscvenv_dqdt_dlnp, iscvenv_dqdt_dlnt, &
           iscvenv_dqut_dlnp, iscvenv_dqut_dlnt, iscvenv_du_dt, &
           iscvenv_log10_p, iscvenv_log10_rho, iscvz_dlnrho_dlnp, &
           iscvz_dlnrho_dlnt, iscvz_log10_du_dt, iscvz_log10_rho, nptsx, &
           scv_nt, scv_weighted_sum4, smix, tablenv, tablex, tabley, tablez, &
           tlogx
! Table dimensions: scv_nt temperature rows, up to scv_np pressure
! points per row (nptsx holds the actual count), scv_nvar variables per
! point in the H, He and envelope tables and scv_nvar_z in the Z table.
      integer, parameter :: scv_nt = 63
      integer, parameter :: scv_np = 76
      integer, parameter :: scv_nvar = 12
      integer, parameter :: scv_nvar_z = 13
      real(dp) :: tlogx(scv_nt), tablex(scv_nt,scv_np,scv_nvar), &
           tabley(scv_nt,scv_np,scv_nvar), smix(scv_nt,scv_np), &
           tablez(scv_nt,scv_np,scv_nvar_z), tablenv(scv_nt,scv_np,scv_nvar)
      integer :: nptsx(scv_nt), idtt, idp
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
      integer, parameter :: iscv_du_dt = 12         ! du/dT, built by build_scv_envelope_table (not read)
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
      integer, parameter :: iscvenv_dqdt_dlnp = 7     ! d(dlnrho/dlnT)/dln(P)
      integer, parameter :: iscvenv_dqdt_dlnt = 8     ! d(dlnrho/dlnT)/dln(T)
      integer, parameter :: iscvenv_dlnp_dlnrho = 9   ! d(lnP)/d(ln rho) along the row
      integer, parameter :: iscvenv_dlncp_dlnt = 10   ! dln(cp)/dln(T)
      integer, parameter :: iscvenv_dqut_dlnp = 11    ! d(du/dT)/dln(P)
      integer, parameter :: iscvenv_dqut_dlnt = 12    ! d(du/dT)/dln(T)

contains

! Four-point Lagrange sum used by compute_scv_envelope_eos/compute_scv_eos to interpolate a table
! column along one axis: w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4),
! evaluated left to right exactly as the former inline expressions.
pure real(dp) function scv_weighted_sum4(w, y)
      real(dp), intent(in) :: w(:), y(:)
      scv_weighted_sum4 = w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4)
end function scv_weighted_sum4

!======================================================================
! Sources: Saumon, Chabrier & Van Horn (1995); the metal table was
! computed with the Prather EOS of the old Yale code. This state was
! the COMMON block scveos.

end module scv_eos_tables
