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

      use opacity_table_lib
      use const_lib
      use luout_lib
      implicit none
      double precision, intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: laol_table_path, laol_table2_path







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
      read(iolaol,130) zhit, debye_huckel_z
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
         read(iolaol2,100) opacity_table%nxyz2,opacity_table%nrho2,opacity_table%nt2
         if (opacity_table%nxyz2.gt.11.or.opacity_table%nrho2.gt.104.or.opacity_table%nt2.gt.52) then
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
         read(iolaol2,140) (opacity_table%oxa2(ii),ii=1,opacity_table%nxyz2)
!        READ IN DENSITY GRID OF TABLE
         read(iolaol2,150) (opacity_table%orho2(ii),ii=1,opacity_table%nrho2)
!        READ IN TEMPERATURE GRID OF TABLE
         read(iolaol2,160) (opacity_table%ot2(ii),ii=1,opacity_table%nt2)
!        READ IN OPACITIES
         read(iolaol2,170)
         do ix=1,opacity_table%nxyz2
            do ir=1,opacity_table%nrho2
               read(iolaol2,200)
               read(iolaol2,210) (opacity_table%olaol2(ix,ir,it),it=1,opacity_table%nt2)
            end do
         end do
         close(iolaol2)
      end if
      return
end subroutine rdlaol
