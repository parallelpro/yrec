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
! variable used to zero eos_table's variable dimension and to reset
! eos_output to 1.0 during the one-time init block below) is reused,
! stale, as the column index into log10_ne_grid in the per-row READ
! statement further down. Since var_idx is SAVE'd and the init block
! only runs once, var_idx is left sitting at mv+1 (its value just
! after the "do var_idx=1,mv" loops complete normally) for every
! subsequent call to this routine, so every row's read of
! log10_ne_grid(density_row, var_idx) always overwrites the SAME
! column (var_idx == mv+1) instead of tracking t6_row as the analogous
! t6_list/amu_grid reads do. This looks like an original bug (using
! "j" instead of "i" in that READ's variable list) but is
! transliterated exactly rather than "fixed".
subroutine readcoeos06

      use const_lib
      use luout_lib
      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197

      double precision :: moles_per_gram_table(mx)
      character(len=1) :: blank_line

! common/aaeos06/: not used in this file; placeholders (see esac06.f90).
      double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
      common/aaeos06/ rho_interp_hi, rho_interp_lo, xxh

! common/aeos06/: see esac06.f90 for the full description.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr)
      integer :: x_loop_index, x_index_lo
      double precision :: x_grid(mx)
      common/aeos06/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_loop_index, &
           x_index_lo, x_grid

! common/beos06/: z_table and t6_index_lo used here; eos_index_inverse
! and eos_var_order used here too (eos_var_order to order the table
! columns being read); density_index_edge is a placeholder here.
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), &
           t6_index_lo(nr), density_index_edge(nt)
      common/beos06/ z_table, eos_index_inverse, eos_var_order, &
           t6_index_lo, density_index_edge

! common/eeos06/: eos_output is zeroed (to 1.0) during the one-time
! init block below.
      double precision :: esact, eos_output(mv)
      common/eeos06/ esact, eos_output

! common/eeeos06/: x_grid_copy used below (see readco.f90's note on
! why the dfsx-equivalent computation reads the copy, not x_grid
! directly). x_interp_workspace not used.
      double precision :: x_interp_workspace(mx,nt,nr), x_grid_copy(mx)
      common/eeeos06/ x_interp_workspace, x_grid_copy

! common/eeeeos06/: this batch's own block, only used by
! readcoeos06.f90. amu_grid/log10_ne_grid are read per (density_row,
! t6_row) but not used elsewhere in this batch (log10_ne_grid is also
! subject to the stale-index bug noted above).
      double precision :: hydrogen_fraction_header(mx), &
           mean_molecular_weight_header(mx), amu_grid(nr,nt), &
           log10_ne_grid(nr,nt), density_grid_table(mx,nr), &
           species_fraction_header(mx,6), log10_r_value(nr,nt)
      integer :: temperature_count_used(mx,nr)
      common/eeeeos06/ moles_per_gram_table, hydrogen_fraction_header, &
           mean_molecular_weight_header, amu_grid, log10_ne_grid, &
           density_grid_table, species_fraction_header, log10_r_value, &
           temperature_count_used



! --- locals ---
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
                     eos_table(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         do var_idx = 1, mv
            eos_output(var_idx) = 1.0d0
         end do
         table_init_flag = 12345678
      end if

      close (2)
! .....read  tables
! MHP 8/25 Moved opening of file to parmin
      do 3 x_loop_index = 1, mx
         read (iopale,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              hydrogen_fraction_header(x_loop_index), &
              z_table(x_loop_index), &
              moles_per_gram_table(x_loop_index), &
              mean_molecular_weight_header(x_loop_index)
         read (iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (species_fraction_header(x_loop_index,species_read_idx), &
              species_read_idx=1,6)
         read (iopale,'(A)') blank_line
         do 2 density_row = 1, nr
            read (iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 density_grid_table(x_loop_index,density_row)
            if (record_number.ne.density_row) then
               write (short_file_unit,'(" OEOS06 Data file incorrect: numtot,jcs= " &
                    &,2I5)') record_number, density_row
               stop
            end if
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
            if (temperature_count_used(x_loop_index,density_row).lt. &
                 t6_index_lo(density_row)) then
               write (short_file_unit,'("Problem with OEOS96 data files: X=",F6.4, &
                    &" density=",E14.4)') hydrogen_fraction_header(x_loop_index), &
                    density_grid_table(x_loop_index,density_row)
               stop
            end if
            do t6_row = 1, temperature_count_used(x_loop_index,density_row)
               if (t6_row.gt.t6_index_lo(density_row)) then
                  read (iopale,'(A)') blank_line
                  go to 4
               end if
               read (iopale,'(F11.6,1X,F6.4,E11.4,2E13.6,2E11.3,5F10.6)') &
                    t6_list(density_row,t6_row), amu_grid(density_row,t6_row), &
                    log10_ne_grid(density_row,var_idx), &
                    (eos_table(x_loop_index,eos_var_order(table_var_idx), &
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
         if (t6_list(1,t6_scan_idx).eq.0.0d0) then
            write(short_file_unit,'("READCOEOS06: Error:",I4, &
                 &"-th T6 value is zero")') t6_scan_idx
            stop
         end if
         t6_grid(t6_scan_idx) = t6_list(1,t6_scan_idx)
      end do
      do 12 t6_idx = 2, nt
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!    12 dfs(i)=1D0/(t6a(i)-t6a(i-1))
         t6_grid_spacing_inv(t6_idx) = 1.0d0/(t6_grid(t6_idx) - &
              t6_grid(t6_idx-1))
   12 continue
      density_grid(1) = density_grid_table(1,1)
      do 13 r_idx = 2, nr
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       rho(i)=rhogr(1,i)
!    13 dfsr(i)=1D0/(rho(i)-rho(i-1))
         density_grid(r_idx) = density_grid_table(1,r_idx)
         density_grid_spacing_inv(r_idx) = 1.0d0/(density_grid(r_idx) - &
              density_grid(r_idx-1))
   13 continue
      do x_idx = 2, mx
         x_grid_spacing_inv(x_idx) = 1.0d0/(x_grid_copy(x_idx) - x_grid_copy(x_idx-1))
      end do
      close (iopale)

      return
end subroutine readcoeos06
