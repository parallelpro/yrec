!----------------------------------------------------------------------
! rdzlaol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rdzlaol.f; only variable names, source form, and comment
! style were updated. FORMAT strings and the READ/OPEN statements
! that must byte-for-byte match the on-disk table file are preserved
! verbatim.
!
! Reads the pure-Z (Z=1) LAOL89 opacity table.
! MHP 10/02 vector v not used
subroutine rdzlaol(pure_z_table_path)

      use const_lib
      use luout_lib
      implicit none
      character(len=256), intent(in) :: pure_z_table_path


! DBG 12/95 ARRAYS FOR PURE Z TABLE
      double precision :: zlaol_opacity(104,52), zlaol_logt_grid(52), &
           zlaol_logrho_grid(104)
      integer :: zlaol_num_rho, zlaol_num_t
      common/zlaol/ zlaol_opacity, zlaol_logt_grid, zlaol_logrho_grid, &
           zlaol_num_rho, zlaol_num_t


! common/newopac/: not used here; declared only to preserve the shared
! storage layout (see getopac.f90 for these names).
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

      save

      integer :: n, i, ii, ir, it
      double precision :: dummy(104)

      open(unit=iopurez, file=pure_z_table_path,form='FORMATTED', &
              status='OLD')
!     READ IN ARRAY SIZES
      read(iopurez,100) n,zlaol_num_rho,zlaol_num_t
  100 format(/,18x,i2,9x,i3,14x,i3)
      if (n.ne.1.or.zlaol_num_rho.gt.104.or.zlaol_num_t.gt.52) then
         write(short_file_unit,*)' Z OPACITY INPUT ERROR.'
         stop
      end if
      read(iopurez,120) (dummy(i),i=1,11)
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
      read(iopurez,131) dummy(1)
  131 format(54x,f8.5,///////)
!     READ IN H MASS FRACTIONS OF TABLE
      read(iopurez,140) dummy(1)
  140 format(/,(1p6e12.5))
!     READ IN DENSITY GRID OF TABLE
      read(iopurez,150) (zlaol_logrho_grid(ii),ii=1,zlaol_num_rho)
  150 format(/,(1p6e12.5))
!     READ IN TEMPERATURE GRID OF TABLE
      read(iopurez,160) (zlaol_logt_grid(ii),ii=1,zlaol_num_t)
  160 format(/,(1p6e12.5))
!     READ IN PURE Z OPACITIES
      read(iopurez,170)
  170 format(1x)
      do ir=1,zlaol_num_rho
         read(iopurez,200)
  200    format(1x)
         read(iopurez,210) (zlaol_opacity(ir,it),it=1,zlaol_num_t)
  210    format(1p6e12.5)
         end do
      close(iopurez)

      return
end subroutine rdzlaol
