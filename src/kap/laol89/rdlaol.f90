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
! runtime-allocated units for the LAOL table (formerly the fixed
! iolaol = 61) and the second (Z-ramp) LAOL table (formerly
! luout_lib's fixed iolaol2 = 63)
      integer :: laol_unit, laol2_unit
      double precision, intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: laol_table_path, laol_table2_path

!     zlot/zhit AND THE work_array2/zdh2/zlot2/zhit2 TWINS ARE
!     READ-LIST TARGETS THAT ARE NEVER USED AFTERWARDS.
      double precision :: work_array2(12), zdh2(18)
      integer :: ii, ix, ir, it
      double precision :: zlot, zhit, zlot2, zhit2

      ierr = 0
      open(newunit=laol_unit,file=laol_table_path, form='FORMATTED', &
           status='OLD')
!     READ IN ARRAY SIZES
      read(laol_unit,100) opacity_table%laol_num_x,opacity_table%laol_num_rho,opacity_table%laol_num_t
  100 format(/,18x,i2,9x,i3,14x,i3)
      if (opacity_table%laol_num_x.gt.n_laol_x-1.or.opacity_table%laol_num_rho.gt.n_laol_rho.or.opacity_table%laol_num_t.gt.n_laol_t) then
         write(run_log_unit,*)' OPACITY ARRAY TOO LARGE.'
         ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
         ierr = 1
         return
      end if
!     READ IN RELATIVE ABUNDANCES BY WEIGTH OF THE METALS C,N,O,NE,NA,
!     MG,SI,AR,FE (THE COX&STEWART MIX).  MIXTURE IS FOR THE ZLOT PART
!     OF THE OPACITY TABLE WITH LAOL RELATIVE ABUNDANCES SCALED
!     TO SUM TO ZLOT.
      read(laol_unit,120) zlot, laol_work_array(ix_c), laol_work_array(ix_n), &
           laol_work_array(ix_o), laol_work_array(ix_ne), laol_work_array(ix_na), &
           laol_work_array(ix_mg), laol_work_array(ix_al), laol_work_array(ix_si), &
           laol_work_array(ix_ar), laol_work_array(ix_fe)
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
!     READ ZHIT THE METAL ABUNDANCE FOR THE HIT PART OF TABLE.
! DBG 7/92 NEED RELATIVE ABUNDANCES OF METALS FOR DEBYE-HUCKEL
!     CORRECTION. 18 ELEMENTS, C,N,O,Ne,Na,Mg,Al,Si,P,
!     S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni scaled to sum to ZHIT
      read(laol_unit,130) zhit, laol_debye_huckel_z
! 130 FORMAT(54X,F8.5,///////)
  130 format(54x,f8.5,/////,1p6e12.5,/,1p6e12.5,/1p6e12.5)
!     READ IN H MASS FRACTIONS OF TABLE
      read(laol_unit,140) (opacity_table%laol_grid_x(ii),ii=1,opacity_table%laol_num_x)
  140 format(/,(1p6e12.5))
!     READ IN DENSITY GRID OF TABLE
      read(laol_unit,150) (opacity_table%laol_grid_rho(ii),ii=1,opacity_table%laol_num_rho)
  150 format(/,(1p6e12.5))
!     READ IN TEMPERATURE GRID OF TABLE
      read(laol_unit,160) (opacity_table%laol_grid_t(ii),ii=1,opacity_table%laol_num_t)
  160 format(/,(1p6e12.5))
!     READ IN OPACITIES
      read(laol_unit,170)
  170 format(1x)
      do ix=1,opacity_table%laol_num_x
         do ir=1,opacity_table%laol_num_rho
            read(laol_unit,200)
  200       format(1x)
            read(laol_unit,210) (opacity_table%laol_opacity(ix,ir,it),it=1,opacity_table%laol_num_t)
  210       format(1p6e12.5)
         end do
      end do
      close(laol_unit)


! DBG 4/94 New stuff follows
      if (star%use_two_z_tables) then
         open(newunit=laol2_unit,file=laol_table2_path, form='FORMATTED', &
              status='OLD')
!        READ IN ARRAY SIZES
         read(laol2_unit,100) opacity_table%laol2_num_x,opacity_table%laol2_num_rho,opacity_table%laol2_num_t
         if (opacity_table%laol2_num_x.gt.n_laol_x-1.or.opacity_table%laol2_num_rho.gt.n_laol_rho.or.opacity_table%laol2_num_t.gt.n_laol_t) then
            write(run_log_unit,*)' SECOND OPACITY ARRAY TOO LARGE.'
            ! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
            ierr = 1
            return
         end if
!        READ IN RELATIVE ABUNDANCES BY WEIGTH OF THE METALS C,N,O,NE,NA,
!        MG,SI,AR,FE (THE COX&STEWART MIX).  MIXTURE IS FOR THE ZLOT PART
!        OF THE OPACITY TABLE WITH LAOL RELATIVE ABUNDANCES SCALED
!        TO SUM TO ZLOT2.
         read(laol2_unit,120) zlot2, work_array2(6), work_array2(9), &
              work_array2(8), work_array2(11), &
              work_array2(1), work_array2(3), work_array2(2), work_array2(5), &
              work_array2(10), work_array2(4)
!        READ ZHIT2 THE METAL ABUNDANCE FOR THE HIT PART OF TABLE.
! DBG 7/92 NEED RELATIVE ABUNDANCES OF METALS FOR DEBYE-HUCKEL
!        CORRECTION. 18 ELEMENTS, C,N,O,Ne,Na,Mg,Al,Si,P,
!        S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni scaled to sum to ZHIT
         read(laol2_unit,130) zhit2, zdh2
!        READ IN H MASS FRACTIONS OF TABLE
         read(laol2_unit,140) (opacity_table%laol2_grid_x(ii),ii=1,opacity_table%laol2_num_x)
!        READ IN DENSITY GRID OF TABLE
         read(laol2_unit,150) (opacity_table%laol2_grid_rho(ii),ii=1,opacity_table%laol2_num_rho)
!        READ IN TEMPERATURE GRID OF TABLE
         read(laol2_unit,160) (opacity_table%laol2_grid_t(ii),ii=1,opacity_table%laol2_num_t)
!        READ IN OPACITIES
         read(laol2_unit,170)
         do ix=1,opacity_table%laol2_num_x
            do ir=1,opacity_table%laol2_num_rho
               read(laol2_unit,200)
               read(laol2_unit,210) (opacity_table%laol2_opacity(ix,ir,it),it=1,opacity_table%laol2_num_t)
            end do
         end do
         close(laol2_unit)
      end if
      return
end subroutine rdlaol
