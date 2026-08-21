!----------------------------------------------------------------------
! readcoeos01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original readcoeos01.f; only variable names, source form, and
! comment style were updated.
!
! Reads the OPAL 2001 EOS data tables into common/aeos/ and the other
! 2001-EOS common blocks (block names preserved verbatim -- see
! readco.f90's header for why). Called once (guarded by
! common/lreadco/) from esac01.f90.
subroutine readcoeos01

      use luout_lib
      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 191

      double precision :: moles_per_gram_table(mx)
      character(len=1) :: blank_line

! common/aaeos/: not used in this file; placeholders (see esac01.f90).
      double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
      common/aaeos/ rho_interp_hi, rho_interp_lo, xxh

! common/aeos/: see esac01.f90 for the full description.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr)
      integer :: x_loop_index, x_index_lo
      double precision :: x_grid(mx)
      common/aeos/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_loop_index, &
           x_index_lo, x_grid

! common/beos/: z_table and t6_index_lo used here; eos_index_inverse
! and eos_var_order are used here too (eos_var_order to order the
! table columns being read).
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), t6_index_lo(nr)
      common/beos/ z_table, eos_index_inverse, eos_var_order, t6_index_lo

! common/eeos/: not used in this file.
      double precision :: esact, eos_output(mv)
      common/eeos/ esact, eos_output

! common/eeeos/: x_grid_copy used below (see readco.f90 for why
! readco.f90's dfsx computation reads the copy, not x_grid directly --
! the analogous logic applies here). x_interp_workspace not used.
      double precision :: x_interp_workspace(mx,nt,nr), x_grid_copy(mx)
      common/eeeos/ x_interp_workspace, x_grid_copy

! common/eeeeos/: this batch's own block, only used by readcoeos01.f90.
      double precision :: hydrogen_fraction_header(mx), &
           mean_molecular_weight_header(mx), density_grid_table(mx,nr), &
           species_fraction_header(mx,6), log10_r_value(nr,nt)
      integer :: temperature_count_used(mx,nr)
      common/eeeeos/ moles_per_gram_table, hydrogen_fraction_header, &
           mean_molecular_weight_header, density_grid_table, &
           species_fraction_header, log10_r_value, temperature_count_used


! common/opaleos/: reused member names from eqstat2.f90 (order/type
! match the original LOPALE,IOPALE,lopale01,lopale06,lNumDeriv).
! iopale doubles as the Fortran unit number the tables are read from.
      logical :: use_opal95_eos, use_opal2001_eos, use_opal2006_eos, &
           lnumderiv
      integer :: iopale
      common/opaleos/ use_opal95_eos, iopale, use_opal2001_eos, &
           use_opal2006_eos, lnumderiv

! common/rmpopeos01/: edge-of-table ramp data (consumed elsewhere,
! e.g. eqbound01.f, out of scope for this batch).
      double precision :: density_edge_at_t(nt)
      integer :: density_index_edge_at_t(nt), t_row_index
      common/rmpopeos01/ density_edge_at_t, density_index_edge_at_t, &
           t_row_index

! --- locals ---
      integer :: table_init_flag
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx, fill_idx
      integer :: density_row, t6_row, record_number
      double precision :: unused_field

      data (density_index_edge_at_t(fill_idx), fill_idx=1,nt) &
           /16*169, 168, 167, 166, 165, 2*164, 163, 2*162, 161, 160, 2*159, &
           4*143, 5*137, 6*134, 2*125, 5*123, 2*122, 6*121, 4*119, 8*116, &
           9*115, 5*113, 7*111, 6*110, 34*109, 107, 104, 40*100, 10*99, &
           98, 97, 96, 95, 94, 93, 92/

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
         table_init_flag = 12345678
      end if

      close (2)
! ..... read  tables
! MHP 8/25 Moved opening of file to parmin
      do 3 x_loop_index = 1, mx

         read (iopale,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              hydrogen_fraction_header(x_loop_index), &
              z_table(x_loop_index), &
              moles_per_gram_table(x_loop_index), &
              mean_molecular_weight_header(x_loop_index)
         read (iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (species_fraction_header(x_loop_index,var_idx), &
              var_idx=1,6)
         read (iopale,'(A)') blank_line
         do 2 density_row = 1, nr
            read (iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 density_grid_table(x_loop_index,density_row)
            if (record_number.ne.density_row) then
               write (short_file_unit,'(" Data file incorrect: numtot,jcs= ",2I5)') &
                    record_number, density_row
               stop
            end if
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
            if (temperature_count_used(x_loop_index,density_row).lt. &
                 t6_index_lo(density_row)) then
               write (short_file_unit,'("problem with data files: X=",F6.4, &
                    &" density=",E14.4)') hydrogen_fraction_header(x_loop_index), &
                    density_grid_table(x_loop_index,density_row)
               stop
            end if
            do t6_row = 1, temperature_count_used(x_loop_index,density_row)
               if (t6_row.gt.t6_index_lo(density_row)) then
                  read (iopale,'(A)') blank_line
                  go to 4
               end if
               read (iopale,'(F9.5,1X,F6.2,3E13.5,6F8.4)') &
                    t6_list(density_row,t6_row), &
                    log10_r_value(density_row,t6_row), &
                    (eos_table(x_loop_index,eos_var_order(var_idx), &
                    t6_row,density_row), var_idx=1,9)
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
            write(short_file_unit,'("READCOEOS01: Error:",I4, &
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
! MHP 7/2003 ADDED RAMP BETWEEN OPAL AND OTHER EOS
! NEED EDGE OF TABLE AT HIGH RHO, FIXED T.
      t_row_index = 1
      do t6_idx = 1, nt
         density_edge_at_t(t6_idx) = density_grid(density_index_edge_at_t(t6_idx))
      end do

      return
end subroutine readcoeos01
