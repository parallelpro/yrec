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
! style were updated. common/ccout2/ member names match those
! established in meqos.f90. The ZAMS-table and centre-table common
! blocks (TAB1A/TAB2A/.../CHE1.../TABX1...) are private to this
! mhdpx/mhdst family and are named consistently across mhdpx2.f90 and
! mhdst.f90.
subroutine mhdpx2(log10_pressure, log10_temperature, table_selector, &
     table_vars, table_hfrac, ndimt)
      use mhd_eos_lib
      use const_lib
      use numerics_lib
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
!     ZAMS TABLES (LABELLED BY A,B,C)
      integer, intent(in) :: table_selector, ndimt
      double precision, intent(in) :: log10_pressure, log10_temperature
      double precision, intent(inout) :: table_vars(ndimt,ivarx), &
           table_hfrac(ndimt)

      save
!     NOMENCLATURE FOR ACCESSING THE TABLES
!     ZAMS TABLES
!     ITBL = -1   : TDVR1A
!     ITBL =  1   : TDVR2A
!     ITBL = -2   : TDVR1B
!     ITBL =  2   : TDVR2B
!     ITBL = -3   : TDVR1C
!     ITBL =  3   : TDVR2C
!     CONTRAL (VARIABLE X) TABLES
!     ITBL =  4   : TDVRX1
!     ITBL =  5   : TDVRX2
!     ITBL =  6   : TDVRX3
!     ITBL =  7   : TDVRX4
!     ITBL =  8   : TDVRX5
!     NOTE: THE SAME OUTPUT VARIABLES VAROUT(ITBL,.) IS USED
!           FOR ITBL=1 AND -1, 2 AND -2, 3 AND -3, RESPECTIVELY.
      double precision :: intpt_work1(ivarx,4), intpt_work2(ivarx,4), &
           intpt_y(ivarx), interpolated_vars(ivarx)

      integer :: i

!     IRANGE = 1
      call zero(interpolated_vars, ivarx)
!     MAIN SELECTION OF TABLES
!     ZAMS TABLES
      if  (table_selector.eq.-1) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%zams_lower_a_table, &
              nt1m, nr1m, ivarc, mhd_eos%zams_lower_log10t, mhd_eos%zams_lower_num_t, &
              mhd_eos%zams_lower_num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarc
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 20       VAROUT(1,I) = VARO(I)
            table_vars(1,i) = interpolated_vars(i)
20       continue
         end do
         table_hfrac(1)    = mhd_eos%zams_a_mass_fraction(1)
      else if (table_selector.eq. 1) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%zams_upper_a_table, &
              nt2m, nr2m, ivarc, mhd_eos%zams_upper_log10t, mhd_eos%zams_upper_num_t, &
              mhd_eos%zams_upper_num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarc
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 30       VAROUT(1,I) = VARO(I)
            table_vars(1,i) = interpolated_vars(i)
30       continue
         end do
         table_hfrac(1)    = mhd_eos%zams_a_mass_fraction(1)
      else if (table_selector.eq.-2) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%zams_lower_b_table, &
              nt1m, nr1m, ivarc, mhd_eos%zams_lower_log10t, mhd_eos%zams_lower_num_t, &
              mhd_eos%zams_lower_num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarc
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 40       VAROUT(2,I) = VARO(I)
            table_vars(2,i) = interpolated_vars(i)
40       continue
         end do
         table_hfrac(2)    = mhd_eos%zams_b_mass_fraction(1)
      else if (table_selector.eq. 2) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%zams_upper_b_table, &
              nt2m, nr2m, ivarc, mhd_eos%zams_upper_log10t, mhd_eos%zams_upper_num_t, &
              mhd_eos%zams_upper_num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarc
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 50       VAROUT(2,I) = VARO(I)
            table_vars(2,i) = interpolated_vars(i)
50       continue
         end do
         table_hfrac(2)    = mhd_eos%zams_b_mass_fraction(1)
      else if (table_selector.eq.-3) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%zams_lower_c_table, &
              nt1m, nr1m, ivarc, mhd_eos%zams_lower_log10t, mhd_eos%zams_lower_num_t, &
              mhd_eos%zams_lower_num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarc
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 60       VAROUT(3,I) = VARO(I)
            table_vars(3,i) = interpolated_vars(i)
60       continue
         end do
         table_hfrac(3)    = mhd_eos%zams_c_mass_fraction(1)
      else if (table_selector.eq. 3) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%zams_upper_c_table, &
              nt2m, nr2m, ivarc, mhd_eos%zams_upper_log10t, mhd_eos%zams_upper_num_t, &
              mhd_eos%zams_upper_num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarc
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 70       VAROUT(3,I) = VARO(I)
            table_vars(3,i) = interpolated_vars(i)
70       continue
         end do
         table_hfrac(3)    = mhd_eos%zams_c_mass_fraction(1)
!     CENTER TABLES
      else if (table_selector.eq. 4) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%centre1_table, &
              ntxm, nrxm, ivarx, mhd_eos%centre_log10t, mhd_eos%centre_num_t, mhd_eos%centre_num_r, &
              intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarx
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 100      VAROUT(4,I) = VARO(I)
            table_vars(4,i) = interpolated_vars(i)
100      continue
         end do
         table_hfrac(4)    = mhd_eos%centre1_mass_fraction(1)
      else if (table_selector.eq. 5) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%centre2_table, &
              ntxm, nrxm, ivarx, mhd_eos%centre_log10t, mhd_eos%centre_num_t, mhd_eos%centre_num_r, &
              intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarx
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 110      VAROUT(5,I) = VARO(I)
            table_vars(5,i) = interpolated_vars(i)
110      continue
         end do
         table_hfrac(5)    = mhd_eos%centre2_mass_fraction(1)
      else if (table_selector.eq. 6) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%centre3_table, &
              ntxm, nrxm, ivarx, mhd_eos%centre_log10t, mhd_eos%centre_num_t, mhd_eos%centre_num_r, &
              intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarx
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 120      VAROUT(6,I) = VARO(I)
            table_vars(6,i) = interpolated_vars(i)
120      continue
         end do
         table_hfrac(6)    = mhd_eos%centre3_mass_fraction(1)
      else if (table_selector.eq. 7) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%centre4_table, &
              ntxm, nrxm, ivarx, mhd_eos%centre_log10t, mhd_eos%centre_num_t, mhd_eos%centre_num_r, &
              intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarx
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 130      VAROUT(7,I) = VARO(I)
            table_vars(7,i) = interpolated_vars(i)
130      continue
         end do
         table_hfrac(7)    = mhd_eos%centre4_mass_fraction(1)
      else if (table_selector.eq. 8) then
         call intpt (log10_pressure, log10_temperature, mhd_eos%centre5_table, &
              ntxm, nrxm, ivarx, mhd_eos%centre_log10t, mhd_eos%centre_num_t, mhd_eos%centre_num_r, &
              intpt_work1, intpt_work2, intpt_y, interpolated_vars)
         do i=1,ivarx
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 140      VAROUT(8,I) = VARO(I)
            table_vars(8,i) = interpolated_vars(i)
140      continue
         end do
         table_hfrac(8)    = mhd_eos%centre5_mass_fraction(1)
      end if
!     END SELECTION OF TABLES
      return
!   999 RETURN
end subroutine mhdpx2
