!----------------------------------------------------------------------
! mhd_eos_tables
!----------------------------------------------------------------------
! The MHD EOS tables (mhd): three ZAMS-type tables for compositions A,
! B and C, each in a lower and an upper temperature part, and up to five
! centre-type tables, each indexed by composition in its last
! dimension (imhd_zams_a..c, imhd_centre_1..5); their temperature
! limits; and mhd_output, the result of the last interpolation, whose
! slots are i_mhd_*. read_mhd_tables (mhd_eos) fills them at eos_init.
module mhd_eos_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: i_mhd_chi_rho, i_mhd_chi_t, i_mhd_dgrad_ad_dlog10_rho, &
           i_mhd_dgrad_ad_dlog10_t, i_mhd_dgrad_ad_dx, i_mhd_dlog10cp_dx, &
           i_mhd_dlog10p_dx, i_mhd_dlog10u_dx, i_mhd_dlogcp_dlogrho, &
           i_mhd_dlogcp_dlogt, i_mhd_eta, i_mhd_grad_ad, i_mhd_ion_frac_1, &
           i_mhd_log10_cp, i_mhd_log10_p, i_mhd_log10_pgas, i_mhd_log10_rho, &
           i_mhd_log10_u, i_mhd_placeholder, imhd_centre_1, imhd_centre_2, &
           imhd_centre_3, imhd_centre_4, imhd_centre_5, imhd_zams_a, &
           imhd_zams_b, imhd_zams_c, mhd, mhd_ivarc, mhd_ivarx, mhd_n_centre, &
           mhd_n_zams, mhd_nchem0, mhd_nr1m, mhd_nr2m, mhd_nrxm, mhd_nt1m, &
           mhd_nt2m, mhd_ntxm
      integer, parameter :: mhd_ivarc = 20
      integer, parameter :: mhd_ivarx = 25
      integer, parameter :: mhd_nchem0 = 6
      integer, parameter :: mhd_nt1m = 16
      integer, parameter :: mhd_nt2m = 79
      integer, parameter :: mhd_ntxm = 10
      integer, parameter :: mhd_nr1m = 87
      integer, parameter :: mhd_nr2m = 21
      integer, parameter :: mhd_nrxm = 21
! Composition index of the ZAMS-type tables (three compositions A, B,
! C read from unit_zams_a/b/c) and of the centre-type tables (up to
! five, read from unit_centre1..5). mhdpx2's table_selector is
! -ic (lower ZAMS), +ic (upper ZAMS) or 3+ic (centre) for these.
      integer, parameter :: mhd_n_zams = 3
      integer, parameter :: mhd_n_centre = 5
      integer, parameter :: imhd_zams_a = 1, imhd_zams_b = 2, imhd_zams_c = 3
      integer, parameter :: imhd_centre_1 = 1, imhd_centre_2 = 2, &
           imhd_centre_3 = 3, imhd_centre_4 = 4, imhd_centre_5 = 5
! Slots of mhd_output (the interpolated MHD table row) that compute_mhd_eos
! reads. The MHD table documentation is not part of this repository;
! these meanings follow from how compute_mhd_eos consumes each slot. Slots 21-24
! are the numerical X-derivatives of slots 2, 3, 8 and 9 built by
! mhdst1 (named i_mhd_*_dx below, from that writer); slot 25 is a
! placeholder. Slot 3 is named only through mhdst1's derivative
! comment (log10 U); nothing reads it directly.
      integer, parameter :: i_mhd_log10_rho = 1     ! log10 density
      integer, parameter :: i_mhd_log10_p = 2       ! log10 total pressure
      integer, parameter :: i_mhd_log10_u = 3       ! log10 internal energy (mhdst1 derivative source)
      integer, parameter :: i_mhd_chi_rho = 4       ! dlnP/dlnRho at constant T
      integer, parameter :: i_mhd_chi_t = 5         ! dlnP/dlnT at constant Rho
      integer, parameter :: i_mhd_grad_ad = 8       ! adiabatic gradient
      integer, parameter :: i_mhd_log10_cp = 9      ! log10 specific heat cp
      integer, parameter :: i_mhd_dgrad_ad_dlog10_rho = 10 ! d(grad_ad)/dlog10(Rho)
      integer, parameter :: i_mhd_dgrad_ad_dlog10_t = 11   ! d(grad_ad)/dlog10(T)
      integer, parameter :: i_mhd_dlogcp_dlogrho = 12      ! dlog(cp)/dlog(Rho)
      integer, parameter :: i_mhd_dlogcp_dlogt = 13        ! dlog(cp)/dlog(T)
      integer, parameter :: i_mhd_ion_frac_1 = 14   ! first of three ionization fractions (14-16)
      integer, parameter :: i_mhd_eta = 18          ! electron degeneracy parameter
      integer, parameter :: i_mhd_log10_pgas = 20   ! log10 gas pressure
! Centre-table extension slots (mhd_ivarc+1..mhd_ivarx), built by
! mhdst1 from the X+dX and X-dX companion tables:
      integer, parameter :: i_mhd_dlog10p_dx = 21   ! d(log10 P)/dX
      integer, parameter :: i_mhd_dlog10u_dx = 22   ! d(log10 U)/dX
      integer, parameter :: i_mhd_dgrad_ad_dx = 23  ! d(grad_ad)/dX
      integer, parameter :: i_mhd_dlog10cp_dx = 24  ! d(log10 cp)/dX
      integer, parameter :: i_mhd_placeholder = 25  ! space-holder (8888844444.d0)

      type :: mhd_eos_state
! the result of the last interpolation (slots i_mhd_*)
           real(dp) :: mhd_output(mhd_ivarx)
! the log10 T boundaries between the table families, and their range
           real(dp) :: zams_lower_upper_boundary_log10t, &
                zams_centre_boundary_log10t, table_log10t_min, table_log10t_max
! the lower-T ZAMS tables, one per composition, on one log10 T grid
           real(dp) :: zams_lower_table(mhd_nt1m,mhd_nr1m,mhd_ivarc,mhd_n_zams), &
                zams_lower_log10t(mhd_nt1m)
           integer :: zams_lower_num_t, zams_lower_num_r
           real(dp) :: zams_lower_drho
! the upper-T ZAMS tables, likewise
           real(dp) :: zams_upper_table(mhd_nt2m,mhd_nr2m,mhd_ivarc,mhd_n_zams), &
                zams_upper_log10t(mhd_nt2m)
           integer :: zams_upper_num_t, zams_upper_num_r
           real(dp) :: zams_upper_drho
! the composition of each ZAMS table
           real(dp) :: zams_atomic_weight(mhd_nchem0,mhd_n_zams), &
                zams_number_abundance(mhd_nchem0,mhd_n_zams), &
                zams_mass_fraction(mhd_nchem0,mhd_n_zams), &
                zams_mean_molecular_weight(mhd_n_zams)
! the centre-type tables, one per composition, on one log10 T grid
           real(dp) :: centre_table(mhd_ntxm,mhd_nrxm,mhd_ivarx,mhd_n_centre), &
                centre_log10t(mhd_ntxm)
           integer :: centre_num_t, centre_num_r
           real(dp) :: centre_drho
! the composition of each centre table
           real(dp) :: centre_atomic_weight(mhd_nchem0,mhd_n_centre), &
                centre_number_abundance(mhd_nchem0,mhd_n_centre), &
                centre_mass_fraction(mhd_nchem0,mhd_n_centre), &
                centre_mean_molecular_weight(mhd_n_centre)
      end type mhd_eos_state

      type(mhd_eos_state) :: mhd

!======================================================================
! Sources: the MHD EOS of Hummer & Mihalas (1988), Mihalas, Daeppen &
! Hummer (1988) and Daeppen et al. (1988); the tables and their
! interpolation package are Daeppen's. This state was the COMMON blocks
! mhdout, tttt, tab1a-c, tab2a-c, chea-c, tabx1-5 and che1-5.

end module mhd_eos_tables
