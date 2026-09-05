!---------------------------  GROUP: SR_X  -------------------------------
!
!----------------------------------------------------------------------
! mhdst
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdst.f; only variable names, source form, and comment
! style were updated.
!
! Reads the 8 MHD equation-of-state table files (3 ZAMS-type unit
! numbers covering compositions A, B, C, plus up to 5 optional
! centre-type unit numbers) into mhd_eos_lib's mhd_eos state and sets
! the table temperature limits used by mhdpx1.
subroutine mhdst(unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
     unit_centre2, unit_centre3, unit_centre4, unit_centre5, ierr)
! `use const_lib` removed (2026): unused (nothing in this file's body
! reaches a const_lib symbol), and its dummy arguments above now share
! names with the const_lib lunum members, which would otherwise be an
! ambiguous reference.
      use mhd_eos_lib
      implicit none

      integer, parameter :: ivarc = mhd_ivarc
      integer, parameter :: ivarx = mhd_ivarx
      integer, parameter :: nchem0 = mhd_nchem0
      integer, parameter :: nt1m = mhd_nt1m
      integer, parameter :: nt2m = mhd_nt2m
      integer, parameter :: ntxm = mhd_ntxm
      integer, parameter :: nr1m = mhd_nr1m
      integer, parameter :: nr2m = mhd_nr2m
      integer, parameter :: nrxm = mhd_nrxm

      integer, intent(in) :: unit_zams_a, unit_zams_b, unit_zams_c, &
           unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

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
!     Placeholder actuals for the mhdst1 arguments that the centre-table
!     calls do not need (the ZAMS lower-table slots).
      integer :: num_chem_species, unused_num_t, unused_num_r, table_index
      double precision :: unused_drho
!     Readability W3 (2026): the eight unit numbers gathered per table
!     family so the eight mhdst1 calls (formerly written out one per
!     composition, in this same order) are two loops over the
!     composition index of mhd_eos's arrays.
      integer :: unit_zams(mhd_n_zams), unit_centre(mhd_n_centre), ic

      integer, intent(out) :: ierr

      ierr = 0

      num_chem_species = 0
      unused_drho = 0.d0
      unused_num_t = 0
      unused_num_r = 0
      unit_zams(imhd_zams_a) = unit_zams_a
      unit_zams(imhd_zams_b) = unit_zams_b
      unit_zams(imhd_zams_c) = unit_zams_c
      unit_centre(imhd_centre_1) = unit_centre1
      unit_centre(imhd_centre_2) = unit_centre2
      unit_centre(imhd_centre_3) = unit_centre3
      unit_centre(imhd_centre_4) = unit_centre4
      unit_centre(imhd_centre_5) = unit_centre5
!     READ ZAMS TABLES
      do ic = 1, mhd_n_zams
         if (unit_zams(ic).gt.0) then
            table_index = 0
            call mhdst1(unit_zams(ic),table_index,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0, &
                        mhd_eos%zams_lower_num_t,mhd_eos%zams_lower_num_r,mhd_eos%zams_upper_num_t,mhd_eos%zams_upper_num_r, &
                        mhd_eos%zams_lower_log10t,mhd_eos%zams_upper_log10t, &
                        mhd_eos%zams_lower_table(:,:,:,ic),mhd_eos%zams_upper_table(:,:,:,ic), &
                        mhd_eos%zams_lower_drho,mhd_eos%zams_upper_drho,num_chem_species, &
                        mhd_eos%zams_atomic_weight(:,ic),mhd_eos%zams_number_abundance(:,ic), &
                        mhd_eos%zams_mass_fraction(:,ic),mhd_eos%zams_mean_molecular_weight(ic), &
                        log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                        atomic_weight_down,atomic_weight_up, &
                        number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
            if (ierr /= 0) return
         end if
      end do
!     READ CENTRE TABLES
      do ic = 1, mhd_n_centre
         if (unit_centre(ic).gt.0) then
            table_index = 1
            call mhdst1(unit_centre(ic),table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                        unused_num_t,unused_num_r,mhd_eos%centre_num_t,mhd_eos%centre_num_r,unused_zams_log10t, &
                        mhd_eos%centre_log10t,unused_zams_table,mhd_eos%centre_table(:,:,:,ic), &
                        unused_drho,mhd_eos%centre_drho,num_chem_species, &
                        mhd_eos%centre_atomic_weight(:,ic),mhd_eos%centre_number_abundance(:,ic), &
                        mhd_eos%centre_mass_fraction(:,ic),mhd_eos%centre_mean_molecular_weight(ic), &
                        log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                        atomic_weight_down,atomic_weight_up, &
                        number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
            if (ierr /= 0) return
         end if
      end do
!     TEMPERATURE LIMITS
      mhd_eos%table_log10t_min = mhd_eos%zams_lower_log10t(  1)
      mhd_eos%zams_lower_upper_boundary_log10t = mhd_eos%zams_lower_log10t(mhd_eos%zams_lower_num_t)
      if (unit_centre1.le.0) then
         mhd_eos%zams_centre_boundary_log10t = mhd_eos%zams_upper_log10t(mhd_eos%zams_upper_num_t)
         mhd_eos%table_log10t_max = mhd_eos%zams_centre_boundary_log10t
      else
         mhd_eos%zams_centre_boundary_log10t = mhd_eos%centre_log10t(  1)
         mhd_eos%table_log10t_max = mhd_eos%centre_log10t(mhd_eos%centre_num_t)
      end if
      return
end subroutine mhdst
