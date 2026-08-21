!----------------------------------------------------------------------
! sulaol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original sulaol.f; only variable names, source form, and comment
! style were updated.
!
! DBG 4/94 Modified to do ZRAMP stuff.
! Builds the natural cubic spline coefficients (in log rho, at each
! tabulated X and T) for the LAOL89 opacity table(s) read by
! rdlaol.f90, for use by gtlaol.f90/gtlaol2.f90.
subroutine sulaol

      use const_lib
      use numerics_lib
      implicit none

! MHP 8/25 Removed unused variables
      double precision :: row_log10_opacity(104), row_log_rho(104), &
           row_d2opacity(104)


      double precision :: olaol2(12,104,52), oxa2(12), ot2(52), orho2(104)
      integer :: nxyz2, nrho2, nt2
      common/nwlaol2/ olaol2, oxa2, ot2, orho2, nxyz2, nrho2, nt2

      double precision :: slaol_opacity(12,104,52), slaol_log_rho(12,104,52), &
           slaol_d2opacity(12,104,52)
      integer :: slaol_num_points(12,52)
      common/slaol/ slaol_opacity, slaol_log_rho, slaol_d2opacity, &
           slaol_num_points

      double precision :: slaol2_opacity(12,104,52), slaol2_log_rho(12,104,52), &
           slaol2_d2opacity(12,104,52)
      integer :: slaol2_num_points(12,52)
      common/slaol2/ slaol2_opacity, slaol2_log_rho, slaol2_d2opacity, &
           slaol2_num_points


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

      save

      integer :: it, ix, ir, num_valid_rho

      do it=1, numt
         ot(it) = log10(ot(it))
      end do
      do ix=1, numofxyz
         do it=1, numt
            num_valid_rho=0
            do ir=1, numrho
                slaol_opacity(ix,ir,it) = 0.0d0
                slaol_log_rho(ix,ir,it) = 0.0d0
                slaol_d2opacity(ix,ir,it) = 0.0d0
                if (olaol(ix,ir,it) .ne. 0.0d0) then
                   num_valid_rho = num_valid_rho+1
                   row_log10_opacity(num_valid_rho) = log10(olaol(ix,ir,it))
                   row_log_rho(num_valid_rho) = log10(orho(ir))
                end if
            end do
            if (num_valid_rho .ge. 4) then
               slaol_num_points(ix,it)=num_valid_rho
               call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    1.0d30, 1.0d30, row_d2opacity)
               do ir=1,num_valid_rho
                   slaol_opacity(ix,ir,it) = row_log10_opacity(ir)
                   slaol_log_rho(ix,ir,it) = row_log_rho(ir)
                   slaol_d2opacity(ix,ir,it) = row_d2opacity(ir)
               end do
            else
               slaol_num_points(ix,it) = 0
            end if
         end do
      end do
! DBG 4/94 Do SPLINE on second opacity table if ZRAMP
      if (use_two_z_tables) then
       do it=1, numt
          ot2(it) = log10(ot2(it))
       end do
         do ix=1, nxyz2
            do it=1, nt2
               num_valid_rho=0
               do ir=1, nrho2
                  slaol2_opacity(ix,ir,it) = 0.0d0
                  slaol2_log_rho(ix,ir,it) = 0.0d0
                  slaol2_d2opacity(ix,ir,it) = 0.0d0
                  if (olaol2(ix,ir,it) .ne. 0.0d0) then
                     num_valid_rho = num_valid_rho+1
                     row_log10_opacity(num_valid_rho) = log10(olaol2(ix,ir,it))
                     row_log_rho(num_valid_rho) = log10(orho2(ir))
                  end if
               end do
               if (num_valid_rho .ge. 4) then
                  slaol2_num_points(ix,it)=num_valid_rho
                  call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                       1.0d30, 1.0d30, row_d2opacity)
                  do ir=1,num_valid_rho
                     slaol2_opacity(ix,ir,it) = row_log10_opacity(ir)
                     slaol2_log_rho(ix,ir,it) = row_log_rho(ir)
                     slaol2_d2opacity(ix,ir,it) = row_d2opacity(ir)
                  end do
               else
                  slaol2_num_points(ix,it) = 0
               end if
            end do
         end do
      end if
      return
end subroutine sulaol
