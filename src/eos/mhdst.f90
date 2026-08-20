!---------------------------  GROUP: SR_X  -------------------------------
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MHDST
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!----------------------------------------------------------------------
! mhdst
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdst.f; only variable names, source form, and comment
! style were updated. The ZAMS-table and centre-table common blocks
! (TAB1A/TAB2A/.../CHE1.../TABX1...) match the names established in
! mhdpx2.f90.
!
! Reads the 8 MHD equation-of-state table files (3 ZAMS-type unit
! numbers covering compositions A, B, C, plus up to 5 optional
! centre-type unit numbers) and sets the shared temperature-range
! commons used by mhdpx1.
subroutine mhdst(unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
     unit_centre2, unit_centre3, unit_centre4, unit_centre5)
      implicit none

      integer, parameter :: ivarc = 20
      integer, parameter :: ivarx = 25
      integer, parameter :: nchem0 = 6
      integer, parameter :: nt1m = 16
      integer, parameter :: nt2m = 79
      integer, parameter :: ntxm = 10
      integer, parameter :: nr1m = 87
      integer, parameter :: nr2m = 21
      integer, parameter :: nrxm = 21

      integer, intent(in) :: unit_zams_a, unit_zams_b, unit_zams_c, &
           unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

!     ZAMS TABLES (LABELLED BY A,B,C)
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch
      double precision :: zams_lower_a_table(nt1m,nr1m,ivarc), &
           zams_lower_log10t(nt1m)
      integer :: zams_lower_num_t, zams_lower_num_r
      double precision :: zams_lower_drho
      common/tab1a/zams_lower_a_table, zams_lower_log10t, &
           zams_lower_num_t, zams_lower_num_r, zams_lower_drho
      double precision :: zams_upper_a_table(nt2m,nr2m,ivarc), &
           zams_upper_log10t(nt2m)
      integer :: zams_upper_num_t, zams_upper_num_r
      double precision :: zams_upper_drho
      common/tab2a/zams_upper_a_table, zams_upper_log10t, &
           zams_upper_num_t, zams_upper_num_r, zams_upper_drho
      double precision :: zams_lower_b_table(nt1m,nr1m,ivarc)
      common/tab1b/zams_lower_b_table
      double precision :: zams_upper_b_table(nt2m,nr2m,ivarc)
      common/tab2b/zams_upper_b_table
      double precision :: zams_lower_c_table(nt1m,nr1m,ivarc)
      common/tab1c/zams_lower_c_table
      double precision :: zams_upper_c_table(nt2m,nr2m,ivarc)
      common/tab2c/zams_upper_c_table
      double precision :: zams_a_atomic_weight(nchem0), &
           zams_a_number_abundance(nchem0), zams_a_mass_fraction(nchem0), &
           zams_a_mean_molecular_weight
      common/chea/zams_a_atomic_weight, zams_a_number_abundance, &
           zams_a_mass_fraction, zams_a_mean_molecular_weight
      double precision :: zams_b_atomic_weight(nchem0), &
           zams_b_number_abundance(nchem0), zams_b_mass_fraction(nchem0), &
           zams_b_mean_molecular_weight
      common/cheb/zams_b_atomic_weight, zams_b_number_abundance, &
           zams_b_mass_fraction, zams_b_mean_molecular_weight
      double precision :: zams_c_atomic_weight(nchem0), &
           zams_c_number_abundance(nchem0), zams_c_mass_fraction(nchem0), &
           zams_c_mean_molecular_weight
      common/chec/zams_c_atomic_weight, zams_c_number_abundance, &
           zams_c_mass_fraction, zams_c_mean_molecular_weight
!     CENTRE TABLES (LABELLED BY 1,2,3,4,5)
      double precision :: centre1_table(ntxm,nrxm,ivarx), centre_log10t(ntxm)
      integer :: centre_num_t, centre_num_r
      double precision :: centre_drho
      common/tabx1/centre1_table, centre_log10t, centre_num_t, &
           centre_num_r, centre_drho
      double precision :: centre2_table(ntxm,nrxm,ivarx)
      common/tabx2/centre2_table
      double precision :: centre3_table(ntxm,nrxm,ivarx)
      common/tabx3/centre3_table
      double precision :: centre4_table(ntxm,nrxm,ivarx)
      common/tabx4/centre4_table
      double precision :: centre5_table(ntxm,nrxm,ivarx)
      common/tabx5/centre5_table
      double precision :: centre1_atomic_weight(nchem0), &
           centre1_number_abundance(nchem0), centre1_mass_fraction(nchem0), &
           centre1_mean_molecular_weight
      common/che1/centre1_atomic_weight, centre1_number_abundance, &
           centre1_mass_fraction, centre1_mean_molecular_weight
      double precision :: centre2_atomic_weight(nchem0), &
           centre2_number_abundance(nchem0), centre2_mass_fraction(nchem0), &
           centre2_mean_molecular_weight
      common/che2/centre2_atomic_weight, centre2_number_abundance, &
           centre2_mass_fraction, centre2_mean_molecular_weight
      double precision :: centre3_atomic_weight(nchem0), &
           centre3_number_abundance(nchem0), centre3_mass_fraction(nchem0), &
           centre3_mean_molecular_weight
      common/che3/centre3_atomic_weight, centre3_number_abundance, &
           centre3_mass_fraction, centre3_mean_molecular_weight
      double precision :: centre4_atomic_weight(nchem0), &
           centre4_number_abundance(nchem0), centre4_mass_fraction(nchem0), &
           centre4_mean_molecular_weight
      common/che4/centre4_atomic_weight, centre4_number_abundance, &
           centre4_mass_fraction, centre4_mean_molecular_weight
      double precision :: centre5_atomic_weight(nchem0), &
           centre5_number_abundance(nchem0), centre5_mass_fraction(nchem0), &
           centre5_mean_molecular_weight
      common/che5/centre5_atomic_weight, centre5_number_abundance, &
           centre5_mass_fraction, centre5_mean_molecular_weight
!     TL<TLIM1:       LOWER PART OF ZAMS TABLES
!     TLIM1<TL<TLIM2: UPPER PART OF ZAMS TABLES
!     TL>TLIM2:       VARIABLE X TABLES
      double precision :: zams_lower_upper_boundary_log10t, &
           zams_centre_boundary_log10t, table_log10t_min, table_log10t_max
      common/tttt/zams_lower_upper_boundary_log10t, &
           zams_centre_boundary_log10t, table_log10t_min, table_log10t_max
      double precision :: unused_zams_table(nt1m,nr1m,ivarc), unused_zams_log10t(nt1m)
!     WORKING STORAGE FOR NUMERICAL X-DERIVATIVES. COMMON TO ALL
!     X-TABLES
      double precision :: log10t_down(nt2m), log10t_up(nt2m)
      double precision :: table_centre_vars(nt2m,nr2m,ivarc)
      double precision :: table_down_vars(nt2m,nr2m,ivarc)
      double precision :: table_up_vars(nt2m,nr2m,ivarc)
      double precision :: atomic_weight_down(nchem0), number_abundance_down(nchem0), &
           mass_fraction_down(nchem0)
      double precision :: atomic_weight_up(nchem0), number_abundance_up(nchem0), &
           mass_fraction_up(nchem0)
      save
!
!     DEFINE, WITH UNUSED STATEMENTS, STORAGE FOR VARIABLES
!     THAT WOULD OTHERWISE ONLY APPEAR AS FORMAL PARAMETERS
!     IN THIS SUBROUTINE. VICIOUS BUGS CAN BE THE RESULT IF
!     THIS STORAGE WERE NOT PROVIDED (REMEMBER: FORTRAN IS NOT
!     A RECURSIVE LANGUAGE.)
      integer :: num_chem_species, unused_num_t, unused_num_r, table_index
      double precision :: unused_drho

      num_chem_species = 0
      unused_drho = 0.d0
      unused_num_t = 0
      unused_num_r = 0
!     READ ZAMS TABLES
      if (unit_zams_a.gt.0) then
         table_index = 0
         call mhdst1(unit_zams_a,table_index,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0, &
                     zams_lower_num_t,zams_lower_num_r,zams_upper_num_t,zams_upper_num_r, &
                     zams_lower_log10t,zams_upper_log10t,zams_lower_a_table,zams_upper_a_table, &
                     zams_lower_drho,zams_upper_drho,num_chem_species,zams_a_atomic_weight, &
                     zams_a_number_abundance,zams_a_mass_fraction,zams_a_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
      if (unit_zams_b.gt.0) then
         table_index = 0
         call mhdst1(unit_zams_b,table_index,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0, &
                     zams_lower_num_t,zams_lower_num_r,zams_upper_num_t,zams_upper_num_r, &
                     zams_lower_log10t,zams_upper_log10t,zams_lower_b_table,zams_upper_b_table, &
                     zams_lower_drho,zams_upper_drho,num_chem_species,zams_b_atomic_weight, &
                     zams_b_number_abundance,zams_b_mass_fraction,zams_b_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
      if (unit_zams_c.gt.0) then
         table_index = 0
         call mhdst1(unit_zams_c,table_index,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0, &
                     zams_lower_num_t,zams_lower_num_r,zams_upper_num_t,zams_upper_num_r, &
                     zams_lower_log10t,zams_upper_log10t,zams_lower_c_table,zams_upper_c_table, &
                     zams_lower_drho,zams_upper_drho,num_chem_species,zams_c_atomic_weight, &
                     zams_c_number_abundance,zams_c_mass_fraction,zams_c_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
!     READ CENTRE TABLES
      if (unit_centre1.gt.0) then
         table_index = 1
         call mhdst1(unit_centre1,table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                     unused_num_t,unused_num_r,centre_num_t,centre_num_r,unused_zams_log10t, &
                     centre_log10t,unused_zams_table,centre1_table, &
                     unused_drho,centre_drho,num_chem_species,centre1_atomic_weight, &
                     centre1_number_abundance,centre1_mass_fraction,centre1_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
      if (unit_centre2.gt.0) then
         table_index = 1
         call mhdst1(unit_centre2,table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                     unused_num_t,unused_num_r,centre_num_t,centre_num_r,unused_zams_log10t, &
                     centre_log10t,unused_zams_table,centre2_table, &
                     unused_drho,centre_drho,num_chem_species,centre2_atomic_weight, &
                     centre2_number_abundance,centre2_mass_fraction,centre2_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
      if (unit_centre3.gt.0) then
         table_index = 1
         call mhdst1(unit_centre3,table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                     unused_num_t,unused_num_r,centre_num_t,centre_num_r,unused_zams_log10t, &
                     centre_log10t,unused_zams_table,centre3_table, &
                     unused_drho,centre_drho,num_chem_species,centre3_atomic_weight, &
                     centre3_number_abundance,centre3_mass_fraction,centre3_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
      if (unit_centre4.gt.0) then
         table_index = 1
         call mhdst1(unit_centre4,table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                     unused_num_t,unused_num_r,centre_num_t,centre_num_r,unused_zams_log10t, &
                     centre_log10t,unused_zams_table,centre4_table, &
                     unused_drho,centre_drho,num_chem_species,centre4_atomic_weight, &
                     centre4_number_abundance,centre4_mass_fraction,centre4_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
      if (unit_centre5.gt.0) then
         table_index = 1
         call mhdst1(unit_centre5,table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                     unused_num_t,unused_num_r,centre_num_t,centre_num_r,unused_zams_log10t, &
                     centre_log10t,unused_zams_table,centre5_table, &
                     unused_drho,centre_drho,num_chem_species,centre5_atomic_weight, &
                     centre5_number_abundance,centre5_mass_fraction,centre5_mean_molecular_weight, &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up)
      end if
!     TEMPERATURE LIMITS
      table_log10t_min = zams_lower_log10t(  1)
      zams_lower_upper_boundary_log10t = zams_lower_log10t(zams_lower_num_t)
      if (unit_centre1.le.0) then
         zams_centre_boundary_log10t = zams_upper_log10t(zams_upper_num_t)
         table_log10t_max = zams_centre_boundary_log10t
      else
         zams_centre_boundary_log10t = centre_log10t(  1)
         table_log10t_max = centre_log10t(centre_num_t)
      end if
      if (unit_centre1.le.0) goto 500
  500 continue
! 8002  FORMAT('      AT. WEIGHT     NUMBER ',
!      & 'ABUNDANCE  MASS FRACTION',(/1X,1P3G16.7))
! 8003  FORMAT(/' MEAN MOLECULAR WEIGHT = ',F12.7//)
      return
end subroutine mhdst
