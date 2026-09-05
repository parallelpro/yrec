!----------------------------------------------------------------------
! rdlaol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rdlaol.f; only variable names, source form, and comment
! style were updated. FORMAT strings and the READ/OPEN statements
! that must byte-for-byte match the on-disk table file are preserved
! verbatim.
!
! DBG 4/94 Modified to read in second table for ZRAMP core stuff.
! Reads the LAOL89 opacity table(s) (and, if a second Z table is
! requested, a second LAOL89 table) into opacity_table%laol(1)/(2).
!
! 2026 wave 3 (R5): the two verbatim read sequences became the
! internal read_laol_table below, called once per table; the second
! table's header mixture was read into work_array2(6),(9),(8),(11),
! (1),(3),(2),(5),(10),(4), which is exactly (ix_c, ix_n, ix_o, ix_ne,
! ix_na, ix_mg, ix_al, ix_si, ix_ar, ix_fe), so one read list serves
! both. The run-log message "SECOND OPACITY ARRAY TOO LARGE." is
! kept verbatim for table 2.
subroutine rdlaol(laol_work_array, laol_debye_huckel_z, laol_table_path, &
     laol_table2_path, ierr)
      use star_info_lib, only: star, ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, &
           ix_o, ix_n, ix_ar, ix_ne

      use opacity_table_lib
      use luout_lib
      implicit none
      integer, intent(out) :: ierr
! the 18-element metal mixture from the table header (see the read
! below); the caller hands it to the eos domain (eos_set_debye_huckel_z)
      double precision, intent(out) :: laol_debye_huckel_z(18)
      double precision, intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: laol_table_path, laol_table2_path

!     THE work_array2/zdh2 TWINS OF laol_work_array/laol_debye_huckel_z
!     ARE READ-LIST TARGETS FOR THE SECOND TABLE THAT ARE NEVER USED
!     AFTERWARDS.
      double precision :: work_array2(12), zdh2(18)
      integer :: ii, ix, ir, it

      ierr = 0
      opacity_table%laol(2)%table_number = 2
      call read_laol_table(laol_table_path, opacity_table%laol(1), &
           laol_work_array, laol_debye_huckel_z, ierr)
      if (ierr /= 0) return

! DBG 4/94 New stuff follows
      if (star%use_two_z_tables) then
         call read_laol_table(laol_table2_path, opacity_table%laol(2), &
              work_array2, zdh2, ierr)
         if (ierr /= 0) return
      end if
      return

contains

! read_laol_table: read one LAOL89 table file at `path` into tbl. The
! header's Cox & Stewart metal mixture goes to work_array (the
! caller's laol_work_array for table 1, a discarded local for table
! 2) and its 18-element Debye-Huckel mixture to debye_huckel_z.
subroutine read_laol_table(path, tbl, work_array, debye_huckel_z, ierr)
      character(len=256), intent(in) :: path
      type(laol_table_set), intent(inout) :: tbl
      double precision, intent(inout) :: work_array(12)
      double precision, intent(out) :: debye_huckel_z(18)
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
! runtime-allocated unit for the table (formerly the fixed iolaol =
! 61 for table 1, luout_lib's fixed iolaol2 = 63 for table 2)
      integer :: laol_unit
!     zlot/zhit ARE READ-LIST TARGETS THAT ARE NEVER USED AFTERWARDS.
      double precision :: zlot, zhit

      open(newunit=laol_unit,file=path, form='FORMATTED', &
           status='OLD')
!     READ IN ARRAY SIZES
      read(laol_unit,100) tbl%num_x,tbl%num_rho,tbl%num_t
  100 format(/,18x,i2,9x,i3,14x,i3)
      if (tbl%num_x.gt.n_laol_x-1.or.tbl%num_rho.gt.n_laol_rho.or.tbl%num_t.gt.n_laol_t) then
         if (tbl%table_number == 2) then
         write(run_log_unit,*)' SECOND OPACITY ARRAY TOO LARGE.'
         else
         write(run_log_unit,*)' OPACITY ARRAY TOO LARGE.'
         end if
         ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
         ierr = 1
         return
      end if
!     READ IN RELATIVE ABUNDANCES BY WEIGTH OF THE METALS C,N,O,NE,NA,
!     MG,SI,AR,FE (THE COX&STEWART MIX).  MIXTURE IS FOR THE ZLOT PART
!     OF THE OPACITY TABLE WITH LAOL RELATIVE ABUNDANCES SCALED
!     TO SUM TO ZLOT.
      read(laol_unit,120) zlot, work_array(ix_c), work_array(ix_n), &
           work_array(ix_o), work_array(ix_ne), work_array(ix_na), &
           work_array(ix_mg), work_array(ix_al), work_array(ix_si), &
           work_array(ix_ar), work_array(ix_fe)
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
!     READ ZHIT THE METAL ABUNDANCE FOR THE HIT PART OF TABLE.
! DBG 7/92 NEED RELATIVE ABUNDANCES OF METALS FOR DEBYE-HUCKEL
!     CORRECTION. 18 ELEMENTS, C,N,O,Ne,Na,Mg,Al,Si,P,
!     S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni scaled to sum to ZHIT
      read(laol_unit,130) zhit, debye_huckel_z
! 130 FORMAT(54X,F8.5,///////)
  130 format(54x,f8.5,/////,1p6e12.5,/,1p6e12.5,/1p6e12.5)
!     READ IN H MASS FRACTIONS OF TABLE
      read(laol_unit,140) (tbl%grid_x(ii),ii=1,tbl%num_x)
  140 format(/,(1p6e12.5))
!     READ IN DENSITY GRID OF TABLE
      read(laol_unit,150) (tbl%grid_rho(ii),ii=1,tbl%num_rho)
  150 format(/,(1p6e12.5))
!     READ IN TEMPERATURE GRID OF TABLE
      read(laol_unit,160) (tbl%grid_t(ii),ii=1,tbl%num_t)
  160 format(/,(1p6e12.5))
!     READ IN OPACITIES
      read(laol_unit,170)
  170 format(1x)
      do ix=1,tbl%num_x
         do ir=1,tbl%num_rho
            read(laol_unit,200)
  200       format(1x)
            read(laol_unit,210) (tbl%opacity(ix,ir,it),it=1,tbl%num_t)
  210       format(1p6e12.5)
         end do
      end do
      close(laol_unit)
end subroutine read_laol_table

end subroutine rdlaol
