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
! requested, a second LAOL89 table).
subroutine rdlaol(laol_work_array, laol_table_path, laol_table2_path)

      implicit none
      double precision, intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: laol_table_path, laol_table2_path

! MHP 8/25 Removed character file names from common block
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez

! DBG 4/94 New common block for second opacity table
      double precision :: olaol2(12,104,52), oxa2(12), ot2(52), orho2(104)
      integer :: nxyz2, nrho2, nt2
      common/nwlaol2/ olaol2, oxa2, ot2, orho2, nxyz2, nrho2, nt2

! MHP 8/25 Removed character file names from common block
! common/zramp/: only iolaol2 is used here.
      double precision :: rsclzc(50), rsclzm1(50), rsclzm2(50)
      integer :: iolaol2, ioopal2, nk
      logical :: use_z_ramp
      common/zramp/ rsclzc, rsclzm1, rsclzm2, iolaol2, ioopal2, nk, &
           use_z_ramp

      integer :: ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr

! DBG 7/92 COMMON BLOCK ADDED TO COMPUTE DEBYE-HUCKEL CORRECTION.
! common/debhu/: only zdh is used here. Naming matches eqstat2.f90.
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/ cdh, etadh0, etadh1, zdh, xxdy, yydh, zzdh, dhnue, ldh

! common/newopac/: only use_two_z_tables is used here.
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max
      logical :: use_alex06_tables, use_laol89_tables, use_opal92_tables, &
           use_opal95_tables, use_kurucz90_tables, use_alex95_tables, &
           use_two_z_tables
      common /newopac/ laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max, use_alex06_tables, &
           use_laol89_tables, use_opal92_tables, use_opal95_tables, &
           use_kurucz90_tables, use_alex95_tables, use_two_z_tables

      double precision :: work_array2(12), zdh2(18)
      save

      integer :: ii, ix, ir, it
      double precision :: zlot, zhit, zlot2, zhit2

      open(unit=iolaol,file=laol_table_path, form='FORMATTED', &
           status='OLD')
!     READ IN ARRAY SIZES
      read(iolaol,100) numofxyz,numrho,numt
  100 format(/,18x,i2,9x,i3,14x,i3)
      if (numofxyz.gt.11.or.numrho.gt.104.or.numt.gt.52) then
         write(short_file_unit,*)' OPACITY ARRAY TOO LARGE.'
         stop
      end if
!     READ IN RELATIVE ABUNDANCES BY WEIGTH OF THE METALS C,N,O,NE,NA,
!     MG,SI,AR,FE (THE COX&STEWART MIX).  MIXTURE IS FOR THE ZLOT PART
!     OF THE OPACITY TABLE WITH LAOL RELATIVE ABUNDANCES SCALED
!     TO SUM TO ZLOT.
      read(iolaol,120) zlot, laol_work_array(6), laol_work_array(9), &
           laol_work_array(8), laol_work_array(11), laol_work_array(1), &
           laol_work_array(3), laol_work_array(2), laol_work_array(5), &
           laol_work_array(10), laol_work_array(4)
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
!     READ ZHIT THE METAL ABUNDANCE FOR THE HIT PART OF TABLE.
! DBG 7/92 NEED RELATIVE ABUNDANCES OF METALS FOR DEBYE-HUCKEL
!     CORRECTION. 18 ELEMENTS, C,N,O,Ne,Na,Mg,Al,Si,P,
!     S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni scaled to sum to ZHIT
      read(iolaol,130) zhit, zdh
! 130 FORMAT(54X,F8.5,///////)
  130 format(54x,f8.5,/////,1p6e12.5,/,1p6e12.5,/1p6e12.5)
!     READ IN H MASS FRACTIONS OF TABLE
      read(iolaol,140) (oxa(ii),ii=1,numofxyz)
  140 format(/,(1p6e12.5))
!     READ IN DENSITY GRID OF TABLE
      read(iolaol,150) (orho(ii),ii=1,numrho)
  150 format(/,(1p6e12.5))
!     READ IN TEMPERATURE GRID OF TABLE
      read(iolaol,160) (ot(ii),ii=1,numt)
  160 format(/,(1p6e12.5))
!     READ IN OPACITIES
      read(iolaol,170)
  170 format(1x)
      do ix=1,numofxyz
         do ir=1,numrho
            read(iolaol,200)
  200       format(1x)
            read(iolaol,210) (olaol(ix,ir,it),it=1,numt)
  210       format(1p6e12.5)
         end do
      end do
      close(iolaol)


! DBG 4/94 New stuff follows
      if (use_two_z_tables) then
         open(unit=iolaol2,file=laol_table2_path, form='FORMATTED', &
              status='OLD')
!        READ IN ARRAY SIZES
         read(iolaol2,100) nxyz2,nrho2,nt2
         if (nxyz2.gt.11.or.nrho2.gt.104.or.nt2.gt.52) then
            write(short_file_unit,*)' SECOND OPACITY ARRAY TOO LARGE.'
            stop
         end if
!        READ IN RELATIVE ABUNDANCES BY WEIGTH OF THE METALS C,N,O,NE,NA,
!        MG,SI,AR,FE (THE COX&STEWART MIX).  MIXTURE IS FOR THE ZLOT PART
!        OF THE OPACITY TABLE WITH LAOL RELATIVE ABUNDANCES SCALED
!        TO SUM TO ZLOT2.
         read(iolaol2,120) zlot2, work_array2(6), work_array2(9), &
              work_array2(8), work_array2(11), &
              work_array2(1), work_array2(3), work_array2(2), work_array2(5), &
              work_array2(10), work_array2(4)
!        READ ZHIT2 THE METAL ABUNDANCE FOR THE HIT PART OF TABLE.
! DBG 7/92 NEED RELATIVE ABUNDANCES OF METALS FOR DEBYE-HUCKEL
!        CORRECTION. 18 ELEMENTS, C,N,O,Ne,Na,Mg,Al,Si,P,
!        S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni scaled to sum to ZHIT
         read(iolaol2,130) zhit2, zdh2
!        READ IN H MASS FRACTIONS OF TABLE
         read(iolaol2,140) (oxa2(ii),ii=1,nxyz2)
!        READ IN DENSITY GRID OF TABLE
         read(iolaol2,150) (orho2(ii),ii=1,nrho2)
!        READ IN TEMPERATURE GRID OF TABLE
         read(iolaol2,160) (ot2(ii),ii=1,nt2)
!        READ IN OPACITIES
         read(iolaol2,170)
         do ix=1,nxyz2
            do ir=1,nrho2
               read(iolaol2,200)
               read(iolaol2,210) (olaol2(ix,ir,it),it=1,nt2)
            end do
         end do
         close(iolaol2)
      end if
      return
end subroutine rdlaol
