!----------------------------------------------------------------------
! mhd_eos_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles the MHD equation-of-state ZAMS-type
! (A/B/C composition) and centre-type (1-5) lookup tables
! -- former common/mhdout/, tttt/, tab1a/, tab2a/, tab1b/, tab2b/,
! tab1c/, tab2c/, chea/, cheb/, chec/, tabx1/, tabx2/, tabx3/, tabx4/,
! tabx5/, che1/, che2/, che3/, che4/, che5/ -- into one derived type.
! These are Category 1 (loaded once at startup, read broadly
! thereafter), matching const_lib's treatment, but are kept in their
! own module rather than folded into const_lib to keep the (large)
! MHD table arrays out of the globally-`use`d const_lib. Only two files
! address the tables (eos/mhd/mhdpx2.f90, eos/mhd/mhdst.f90); no DATA
! statements target any member.
!
! Readability W3 (2026): the per-composition members (former
! zams_lower_a_table/_b_/_c_, zams_a_atomic_weight/_b_/_c_, ...,
! centre1_table..centre5_table, centre1_atomic_weight..centre5_...)
! became one array each with a trailing composition index, so mhdst's
! eight reads and mhdpx2's table dispatch are loops/index arithmetic
! over imhd_zams_a..c / imhd_centre_1..5 instead of copies. Each
! composition's slice (:,:,:,ic) has the same shape and element order
! as the member it replaces.
module mhd_eos_lib
      implicit none
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
! Slots of mhd_output (the interpolated MHD table row) that meqos.f90
! reads. The MHD table documentation is not part of this repository;
! these meanings follow from how meqos consumes each slot. Slots 21-24
! are the numerical X-derivatives of slots 2, 3, 8 and 9 built by
! mhdst1.f90; slot 25 is a placeholder.
      integer, parameter :: i_mhd_log10_rho = 1     ! log10 density
      integer, parameter :: i_mhd_log10_p = 2       ! log10 total pressure
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

      type, public :: mhd_eos_state
! former common/mhdout/
           double precision :: mhd_output(mhd_ivarx)
! former common/tttt/
           double precision :: zams_lower_upper_boundary_log10t, &
                zams_centre_boundary_log10t, table_log10t_min, table_log10t_max
! former common/tab1a/ (composition A's table plus the lower-T grid
! shared by all three ZAMS compositions)
           double precision :: zams_lower_table(mhd_nt1m,mhd_nr1m,mhd_ivarc,mhd_n_zams), &
                zams_lower_log10t(mhd_nt1m)
           integer :: zams_lower_num_t, zams_lower_num_r
           double precision :: zams_lower_drho
! former common/tab2a/ (likewise for the upper-T ZAMS tables)
           double precision :: zams_upper_table(mhd_nt2m,mhd_nr2m,mhd_ivarc,mhd_n_zams), &
                zams_upper_log10t(mhd_nt2m)
           integer :: zams_upper_num_t, zams_upper_num_r
           double precision :: zams_upper_drho
! former common/chea/, cheb/, chec/ (per ZAMS composition)
           double precision :: zams_atomic_weight(mhd_nchem0,mhd_n_zams), &
                zams_number_abundance(mhd_nchem0,mhd_n_zams), &
                zams_mass_fraction(mhd_nchem0,mhd_n_zams), &
                zams_mean_molecular_weight(mhd_n_zams)
! former common/tabx1/ .. tabx5/ (centre-type tables, one per
! composition, sharing one log10 T grid)
           double precision :: centre_table(mhd_ntxm,mhd_nrxm,mhd_ivarx,mhd_n_centre), &
                centre_log10t(mhd_ntxm)
           integer :: centre_num_t, centre_num_r
           double precision :: centre_drho
! former common/che1/ .. che5/ (per centre composition)
           double precision :: centre_atomic_weight(mhd_nchem0,mhd_n_centre), &
                centre_number_abundance(mhd_nchem0,mhd_n_centre), &
                centre_mass_fraction(mhd_nchem0,mhd_n_centre), &
                centre_mean_molecular_weight(mhd_n_centre)
      end type mhd_eos_state

      type(mhd_eos_state), public :: mhd_eos

end module mhd_eos_lib
