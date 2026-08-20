!----------------------------------------------------------------------
! readco
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original readco.f; only variable names, source form, and comment
! style were updated.
!
! Reads the OPAL 1995 EOS data tables (one block per tabulated X
! value) into common/a/ and the other 1995-EOS common blocks (block
! names are preserved verbatim from the original -- Fortran COMMON is
! linked by name, and other translation units outside this batch,
! e.g. eqbound.f, still reference these blocks under their original
! names), and derives the auxiliary grid-spacing arrays used by
! esac.f90/t6rinterp.f90 for interpolation. Called once (guarded by
! common/lreadco/) from esac.f90.
subroutine readco

      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 77, nt = 56

! common/opaleos/: reused member names from eqstat2.f90 (order/type
! match the original LOPALE,IOPALE,lopale01,lopale06,lNumDeriv).
! iopale doubles as the Fortran unit number the tables are read from.
      logical :: use_opal95_eos, use_opal2001_eos, use_opal2006_eos, &
           lnumderiv
      integer :: iopale
      common/opaleos/ use_opal95_eos, iopale, use_opal2001_eos, &
           use_opal2006_eos, lnumderiv

! common/luout/: only short_file_unit (the .short log unit) is used
! here; the rest are placeholders. Naming matches getopac.f90.
      integer :: ilast, idebug, itrack, short_file_unit, imilne, imodpt, &
           istor, main_output_unit
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, main_output_unit

      double precision :: moles_per_gram_table(mx)
      character(len=1) :: blank_line

! common/aa/: rho_interp_hi, rho_interp_lo, xxh are not used in this
! file; declared only to preserve the storage layout shared with
! esac.f90/t6rinterp.f90 (this batch's own common block, not one of
! the reference files' -- names chosen here and reused consistently
! across every file in this batch that shares it).
      double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
      common/aa/ rho_interp_hi, rho_interp_lo, xxh

! common/ee/: x_grid_copy is used below (esac.f90 copies x_grid into
! it before calling readco, so by the time we read it here it equals
! x_grid -- preserved verbatim rather than reading x_grid directly).
! x_interp_workspace is not used in this file.
      double precision :: x_interp_workspace(mx,nt,nr), x_grid_copy(mx)
      common/ee/ x_interp_workspace, x_grid_copy

! common/a/: the main OPAL 1995 EOS table and its interpolation
! grids/scratch arrays. Names documented in esac.f90.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr), &
           x_grid(mx)
      integer :: x_loop_index, x_index_lo
      common/a/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_grid, &
           x_loop_index, x_index_lo

! common/b/: z_table and t6_index_lo are used here; eos_index_inverse
! and eos_var_order are used here too (eos_var_order to order the
! table columns being read; eos_index_inverse is a placeholder here).
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), &
           t6_index_lo(nr)
      common/b/ z_table, eos_index_inverse, eos_var_order, t6_index_lo

! common/e/: esact, eos_output not used in this file.
      double precision :: esact, eos_output(mv)
      common/e/ esact, eos_output

! common/eee/: this batch's own block, only used by readco.f90.
      double precision :: hydrogen_fraction_header(mx), &
           mean_molecular_weight_header(mx), density_grid_table(mx,nr), &
           species_fraction_header(mx,6), log10_r_value(nr,nt)
      integer :: temperature_count_used(mx,nr)
      common/eee/ moles_per_gram_table, hydrogen_fraction_header, &
           mean_molecular_weight_header, density_grid_table, &
           species_fraction_header, log10_r_value, temperature_count_used

! common/rmpopeos/: edge-of-table ramp data (consumed elsewhere,
! e.g. eqbound.f, out of scope for this batch).
      double precision :: density_edge_at_t(nt)
      integer :: density_index_edge_at_t(nt), t_row_index
      common/rmpopeos/ density_edge_at_t, density_index_edge_at_t, &
           t_row_index

! --- locals ---
      integer :: table_init_flag
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number, t6_count_used
      double precision :: unused_field

      data (density_index_edge_at_t(t6_scan_idx), t6_scan_idx=1,nt) &
           /7*77, 2*76, 2*74, 2*72, 2*70, 68, 67, 66, 65, 64, 63, 61, &
           60, 59, 58, 57, 55, 54, 53, 52, 51, 2*49, 48, 2*47, 46, &
           2*45, 15*44, 2*37/

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

! ..... read  tables
! MHP 8/25 Moved opening of file to parmin
      do 3 x_loop_index = 1, mx

         read (iopale,'(3X,F6.4,3X,F6.4,11X,F10.7,17X,F10.7)') &
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
               write(short_file_unit,'(" DATA FILE INCORRECT: NUMTOT,JCS= ",2I5)') &
                    record_number, density_row
               stop
            end if
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
            if (temperature_count_used(x_loop_index,density_row).lt. &
                 t6_index_lo(density_row)) then
               write(short_file_unit,'("PROBLEM WITH DATA FILES: X=",F6.4," DENSITY=", &
                    &E14.4)') hydrogen_fraction_header(x_loop_index), &
                    density_grid_table(x_loop_index,density_row)
               stop
            end if
            do t6_row = 1, temperature_count_used(x_loop_index,density_row)
               if (t6_row.gt.t6_index_lo(density_row)) then
                  read (iopale,'(A)') blank_line
                  go to 4
               end if
               read (iopale,'(F8.4,1X,F6.2,3E13.5,E11.3,6F8.4)') &
                    t6_list(density_row,t6_row), &
                    log10_r_value(density_row,t6_row), &
                    (eos_table(x_loop_index,eos_var_order(var_idx), &
                    t6_row,density_row), var_idx=1,10)
    4          continue
            end do
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
            read(iopale,'(A)') blank_line
    2    continue
         read(iopale,'(A)') blank_line
    3 continue

      do 11 t6_scan_idx = 1, nt
         if (t6_list(1,t6_scan_idx).eq.0.0d0) then
            t6_count_used = t6_scan_idx
            go to 14
         end if
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!    11 T6A(I)=T6LIST(1,I)
         t6_grid(t6_scan_idx) = t6_list(1,t6_scan_idx)
   11 continue
   14 do 12 t6_idx = 2, nt
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!    12 DFS(I)=1.D0/(T6A(I)-T6A(I-1))
         t6_grid_spacing_inv(t6_idx) = 1.0d0/(t6_grid(t6_idx) - &
              t6_grid(t6_idx-1))
   12 continue
      density_grid(1) = density_grid_table(1,1)
      do 13 r_idx = 2, nr
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       RHO(I)=RHOGR(1,I)
!    13 DFSR(I)=1.D0/(RHO(I)-RHO(I-1))
         density_grid(r_idx) = density_grid_table(1,r_idx)
         density_grid_spacing_inv(r_idx) = 1.0d0/(density_grid(r_idx) - &
              density_grid(r_idx-1))
   13 continue
      do x_idx = 2, mx
         x_grid_spacing_inv(x_idx) = 1.0d0/(x_grid_copy(x_idx) - x_grid_copy(x_idx-1))
      end do

      close (iopale)
!  MHP 8/98 ADDED RAMP BETWEEN OPAL AND OTHER EOS
!  NEED EDGE OF TABLE AT HIGH RHO, FIXED T.
      t_row_index = 1
      do t6_idx = 1, nt
         density_edge_at_t(t6_idx) = density_grid(density_index_edge_at_t(t6_idx))
      end do
      return
end subroutine readco
