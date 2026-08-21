!----------------------------------------------------------------------
! readcoeos06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original readcoeos06.f; only variable names, source form, and
! comment style were updated.
!
! Reads the OPAL 2006 EOS data tables into common/aeos06/ and the
! other 2006-EOS common blocks (block names preserved verbatim -- see
! readco.f90's header for why). Called once (guarded by
! common/lreadco/) from esac06.f90.
!
! NOTE: preserved verbatim from the original -- var_idx (the "j" loop
! variable used to zero opal_eos%eos_table_06's variable dimension and to reset
! atm_table%eos_output to 1.0 during the one-time init block below) is reused,
! stale, as the column index into log10_ne_grid in the per-row READ
! statement further down. Since var_idx is SAVE'd and the init block
! only runs once, var_idx is left sitting at mv+1 (its value just
! after the "do var_idx=1,mv" loops complete normally) for every
! subsequent call to this routine, so every row's read of
! log10_ne_grid(density_row, var_idx) always overwrites the SAME
! column (var_idx == mv+1) instead of tracking t6_row as the analogous
! opal_eos%t6_list_06/amu_grid reads do. This looks like an original bug (using
! "j" instead of "i" in that READ's variable list) but is
! transliterated exactly rather than "fixed".
subroutine readcoeos06

      use opal_eos_lib
      use atm_table_lib
      use const_lib
      use luout_lib
      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197

      character(len=1) :: blank_line

! former common/eeeeos06/: this batch's own block, only used by
! readcoeos06.f90, now use-associated from opal_eos_lib.
! amu_grid_06/log10_ne_grid_06 are read per (density_row, t6_row) but
! not used elsewhere in this batch (log10_ne_grid_06 is also subject
! to the stale-index bug noted above).



! --- locals ---
      integer :: x_loop_index_06
      integer :: table_init_flag
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      integer :: species_read_idx, table_var_idx
      double precision :: unused_field

      save

      blank_line = ' '

      if (table_init_flag.ne.12345678) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     opal_eos%eos_table_06(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         do var_idx = 1, mv
            atm_table%eos_output(var_idx) = 1.0d0
         end do
         table_init_flag = 12345678
      end if

      close (2)
! .....read  tables
! MHP 8/25 Moved opening of file to parmin
      do 3 x_loop_index_06 = 1, mx
         read (iopale,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              opal_eos%hydrogen_fraction_header_06(x_loop_index_06), &
              opal_eos%z_table_06(x_loop_index_06), &
              opal_eos%moles_per_gram_table_06(x_loop_index_06), &
              opal_eos%mean_molecular_weight_header_06(x_loop_index_06)
         read (iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (opal_eos%species_fraction_header_06(x_loop_index_06,species_read_idx), &
              species_read_idx=1,6)
         read (iopale,'(A)') blank_line
         do 2 density_row = 1, nr
            read (iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 opal_eos%temperature_count_used_06(x_loop_index_06,density_row), &
                 unused_field, unused_field, &
                 opal_eos%density_grid_table_06(x_loop_index_06,density_row)
            if (record_number.ne.density_row) then
               write (short_file_unit,'(" OEOS06 Data file incorrect: numtot,jcs= " &
                    &,2I5)') record_number, density_row
               stop
            end if
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
            if (opal_eos%temperature_count_used_06(x_loop_index_06,density_row).lt. &
                 opal_eos%t6_index_lo_06(density_row)) then
               write (short_file_unit,'("Problem with OEOS96 data files: X=",F6.4, &
                    &" density=",E14.4)') opal_eos%hydrogen_fraction_header_06(x_loop_index_06), &
                    opal_eos%density_grid_table_06(x_loop_index_06,density_row)
               stop
            end if
            do t6_row = 1, opal_eos%temperature_count_used_06(x_loop_index_06,density_row)
               if (t6_row.gt.opal_eos%t6_index_lo_06(density_row)) then
                  read (iopale,'(A)') blank_line
                  go to 4
               end if
               read (iopale,'(F11.6,1X,F6.4,E11.4,2E13.6,2E11.3,5F10.6)') &
                    opal_eos%t6_list_06(density_row,t6_row), opal_eos%amu_grid_06(density_row,t6_row), &
                    opal_eos%log10_ne_grid_06(density_row,var_idx), &
                    (opal_eos%eos_table_06(x_loop_index_06,opal_eos%eos_var_order_06(table_var_idx), &
                    t6_row,density_row), table_var_idx=1,9)
    4          continue
            end do
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
    2    continue
         read(iopale,'(A)') blank_line
    3 continue

      do t6_scan_idx = 1, nt
         if (opal_eos%t6_list_06(1,t6_scan_idx).eq.0.0d0) then
            write(short_file_unit,'("READCOEOS06: Error:",I4, &
                 &"-th T6 value is zero")') t6_scan_idx
            stop
         end if
         opal_eos%t6_grid_06(t6_scan_idx) = opal_eos%t6_list_06(1,t6_scan_idx)
      end do
      do 12 t6_idx = 2, nt
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!    12 dfs(i)=1D0/(t6a(i)-t6a(i-1))
         opal_eos%t6_grid_spacing_inv_06(t6_idx) = 1.0d0/(opal_eos%t6_grid_06(t6_idx) - &
              opal_eos%t6_grid_06(t6_idx-1))
   12 continue
      opal_eos%density_grid_06(1) = opal_eos%density_grid_table_06(1,1)
      do 13 r_idx = 2, nr
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       rho(i)=rhogr(1,i)
!    13 dfsr(i)=1D0/(rho(i)-rho(i-1))
         opal_eos%density_grid_06(r_idx) = opal_eos%density_grid_table_06(1,r_idx)
         opal_eos%density_grid_spacing_inv_06(r_idx) = 1.0d0/(opal_eos%density_grid_06(r_idx) - &
              opal_eos%density_grid_06(r_idx-1))
   13 continue
      do x_idx = 2, mx
         opal_eos%x_grid_spacing_inv_06(x_idx) = 1.0d0/(opal_eos%x_grid_copy_06(x_idx) - opal_eos%x_grid_copy_06(x_idx-1))
      end do
      close (iopale)

      return
end subroutine readcoeos06
