!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MHDPX2
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!----------------------------------------------------------------------
! mhdpx2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdpx2.f; only variable names, source form, and comment
! style were updated. The ZAMS-type and centre-type tables live in
! mhd_eos_lib's mhd_eos (loaded by mhdst.f90).
!
! Interpolates one table (selected by table_selector: -1/-2/-3 lower
! ZAMS A/B/C, 1/2/3 upper ZAMS A/B/C, 4..8 centre 1..5) at
! (log10 P, log10 T) and stores the result as row |table_selector|
! (ZAMS) or table_selector (centre) of table_vars, with that table's
! hydrogen mass fraction in table_hfrac. The tables are the
! composition-indexed arrays of mhd_eos (imhd_zams_a..c,
! imhd_centre_1..5 in mhd_eos_lib).
subroutine mhdpx2(log10_pressure, log10_temperature, table_selector, &
     table_vars, table_hfrac, ndimt)
      use mhd_eos_lib
      use numerics_lib
      implicit none
      integer, parameter :: ivarc = mhd_ivarc
      integer, parameter :: ivarx = mhd_ivarx
      integer, parameter :: nt1m = mhd_nt1m
      integer, parameter :: nt2m = mhd_nt2m
      integer, parameter :: ntxm = mhd_ntxm
      integer, parameter :: nr1m = mhd_nr1m
      integer, parameter :: nr2m = mhd_nr2m
      integer, parameter :: nrxm = mhd_nrxm
!     ZAMS TABLES (LABELLED BY A,B,C)
      integer, intent(in) :: table_selector, ndimt
      double precision, intent(in) :: log10_pressure, log10_temperature
      double precision, intent(inout) :: table_vars(ndimt,ivarx), &
           table_hfrac(ndimt)
!     NOMENCLATURE FOR ACCESSING THE TABLES
!     ZAMS TABLES
!     CONTRAL (VARIABLE X) TABLES
!     NOTE: THE SAME OUTPUT VARIABLES VAROUT(ITBL,.) IS USED
!           FOR ITBL=1 AND -1, 2 AND -2, 3 AND -3, RESPECTIVELY.
      double precision :: intpt_work1(ivarx,4), intpt_work2(ivarx,4), &
           intpt_y(ivarx), interpolated_vars(ivarx)

      integer :: i, ic

      call zero(interpolated_vars, ivarx)
!     MAIN SELECTION OF TABLES
!     Readability W3 (2026): the former 11-way ladder over the named
!     per-composition members is now index arithmetic into mhd_eos's
!     composition-indexed arrays: selector -ic is lower ZAMS composition
!     ic (row ic), +ic is upper ZAMS composition ic (row ic), and
!     mhd_n_zams+ic is centre composition ic (row mhd_n_zams+ic). Any
!     other selector still does nothing, as before.
!     ZAMS TABLES
      if (table_selector.le.-1 .and. table_selector.ge.-mhd_n_zams) then
         ic = -table_selector
         call interpolate_and_store(mhd_eos%zams_lower_table(:,:,:,ic), &
              nt1m, nr1m, ivarc, mhd_eos%zams_lower_log10t, mhd_eos%zams_lower_num_t, &
              mhd_eos%zams_lower_num_r, ic, mhd_eos%zams_mass_fraction(1,ic))
      else if (table_selector.ge.1 .and. table_selector.le.mhd_n_zams) then
         ic = table_selector
         call interpolate_and_store(mhd_eos%zams_upper_table(:,:,:,ic), &
              nt2m, nr2m, ivarc, mhd_eos%zams_upper_log10t, mhd_eos%zams_upper_num_t, &
              mhd_eos%zams_upper_num_r, ic, mhd_eos%zams_mass_fraction(1,ic))
!     CENTER TABLES
      else if (table_selector.ge.mhd_n_zams+1 .and. &
           table_selector.le.mhd_n_zams+mhd_n_centre) then
         ic = table_selector - mhd_n_zams
         call interpolate_and_store(mhd_eos%centre_table(:,:,:,ic), &
              ntxm, nrxm, ivarx, mhd_eos%centre_log10t, mhd_eos%centre_num_t, &
              mhd_eos%centre_num_r, table_selector, mhd_eos%centre_mass_fraction(1,ic))
      end if
!     END SELECTION OF TABLES
      return

contains

! Interpolate one table at (log10 P, log10 T) and store its num_vars
! results as row `row` of table_vars, with the table's hydrogen mass
! fraction in table_hfrac(row). The work arrays and the (P,T) point are
! host-associated from mhdpx2.
subroutine interpolate_and_store(table_data, table_dim_t, table_dim_r, &
     num_vars, table_log10t, num_t, num_r, row, hydrogen_mass_fraction)
      integer, intent(in) :: table_dim_t, table_dim_r, num_vars, num_t, num_r, row
      double precision, intent(in) :: table_data(table_dim_t,table_dim_r,num_vars)
      double precision, intent(in) :: table_log10t(table_dim_t)
      double precision, intent(in) :: hydrogen_mass_fraction

      call intpt (log10_pressure, log10_temperature, table_data, &
           table_dim_t, table_dim_r, num_vars, table_log10t, num_t, &
           num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
      do i=1,num_vars
         table_vars(row,i) = interpolated_vars(i)
      end do
      table_hfrac(row)    = hydrogen_mass_fraction
end subroutine interpolate_and_store

end subroutine mhdpx2
