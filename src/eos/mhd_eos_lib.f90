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
! MHD table arrays out of the globally-`use`d const_lib. All members
! use byte-identical names/order across their only two declaring files
! (eos/mhd/mhdpx2.f90, eos/mhd/mhdst.f90); no DATA statements target any
! member.
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

      type, public :: mhd_eos_state
! former common/mhdout/
           double precision :: mhd_output(mhd_ivarx)
! former common/tttt/
           double precision :: zams_lower_upper_boundary_log10t, &
                zams_centre_boundary_log10t, table_log10t_min, table_log10t_max
! former common/tab1a/
           double precision :: zams_lower_a_table(mhd_nt1m,mhd_nr1m,mhd_ivarc), &
                zams_lower_log10t(mhd_nt1m)
           integer :: zams_lower_num_t, zams_lower_num_r
           double precision :: zams_lower_drho
! former common/tab2a/
           double precision :: zams_upper_a_table(mhd_nt2m,mhd_nr2m,mhd_ivarc), &
                zams_upper_log10t(mhd_nt2m)
           integer :: zams_upper_num_t, zams_upper_num_r
           double precision :: zams_upper_drho
! former common/tab1b/, tab2b/, tab1c/, tab2c/
           double precision :: zams_lower_b_table(mhd_nt1m,mhd_nr1m,mhd_ivarc)
           double precision :: zams_upper_b_table(mhd_nt2m,mhd_nr2m,mhd_ivarc)
           double precision :: zams_lower_c_table(mhd_nt1m,mhd_nr1m,mhd_ivarc)
           double precision :: zams_upper_c_table(mhd_nt2m,mhd_nr2m,mhd_ivarc)
! former common/chea/, cheb/, chec/
           double precision :: zams_a_atomic_weight(mhd_nchem0), &
                zams_a_number_abundance(mhd_nchem0), &
                zams_a_mass_fraction(mhd_nchem0), &
                zams_a_mean_molecular_weight
           double precision :: zams_b_atomic_weight(mhd_nchem0), &
                zams_b_number_abundance(mhd_nchem0), &
                zams_b_mass_fraction(mhd_nchem0), &
                zams_b_mean_molecular_weight
           double precision :: zams_c_atomic_weight(mhd_nchem0), &
                zams_c_number_abundance(mhd_nchem0), &
                zams_c_mass_fraction(mhd_nchem0), &
                zams_c_mean_molecular_weight
! former common/tabx1/
           double precision :: centre1_table(mhd_ntxm,mhd_nrxm,mhd_ivarx), &
                centre_log10t(mhd_ntxm)
           integer :: centre_num_t, centre_num_r
           double precision :: centre_drho
! former common/tabx2/, tabx3/, tabx4/, tabx5/
           double precision :: centre2_table(mhd_ntxm,mhd_nrxm,mhd_ivarx)
           double precision :: centre3_table(mhd_ntxm,mhd_nrxm,mhd_ivarx)
           double precision :: centre4_table(mhd_ntxm,mhd_nrxm,mhd_ivarx)
           double precision :: centre5_table(mhd_ntxm,mhd_nrxm,mhd_ivarx)
! former common/che1/, che2/, che3/, che4/, che5/
           double precision :: centre1_atomic_weight(mhd_nchem0), &
                centre1_number_abundance(mhd_nchem0), &
                centre1_mass_fraction(mhd_nchem0), &
                centre1_mean_molecular_weight
           double precision :: centre2_atomic_weight(mhd_nchem0), &
                centre2_number_abundance(mhd_nchem0), &
                centre2_mass_fraction(mhd_nchem0), &
                centre2_mean_molecular_weight
           double precision :: centre3_atomic_weight(mhd_nchem0), &
                centre3_number_abundance(mhd_nchem0), &
                centre3_mass_fraction(mhd_nchem0), &
                centre3_mean_molecular_weight
           double precision :: centre4_atomic_weight(mhd_nchem0), &
                centre4_number_abundance(mhd_nchem0), &
                centre4_mass_fraction(mhd_nchem0), &
                centre4_mean_molecular_weight
           double precision :: centre5_atomic_weight(mhd_nchem0), &
                centre5_number_abundance(mhd_nchem0), &
                centre5_mass_fraction(mhd_nchem0), &
                centre5_mean_molecular_weight
      end type mhd_eos_state

      type(mhd_eos_state), public :: mhd_eos

end module mhd_eos_lib
