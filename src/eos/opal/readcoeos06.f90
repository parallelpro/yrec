!----------------------------------------------------------------------
! readcoeos06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original readcoeos06.f; only variable names, source form, and
! comment style were updated.
!
! Reads the OPAL 2006 EOS data tables into opal_eos_lib's opal_eos state
! and derives the auxiliary grid-spacing arrays used by esac06.f90/
! t6rinteos06.f90 for interpolation. Called once from esac06.f90
! (guarded by opal_eos%table_loaded_flag).
!
! NOTE: preserved verbatim from the original -- var_idx, the "j" loop
! variable of the one-time init block below (which zeroes
! v%eos_table and resets v%eos_output to 1.0), is
! reused, stale, as the column index into log10_ne_grid_06 in the
! per-row READ statement further down. Because the init loops leave
! var_idx at mv+1, every row's read of
! log10_ne_grid_06(density_row, var_idx) overwrites the SAME column
! (var_idx == mv+1) instead of tracking t6_row as the analogous
! v%t6_list/amu_grid_06 reads do. This looks like an original
! bug (using "j" instead of "i" in that READ's variable list) but is
! transliterated exactly rather than "fixed"; log10_ne_grid_06 is never
! read anywhere.
subroutine readcoeos06(v, ierr)
      use star_info_lib, only: star

      use opal_eos_lib
      use luout_lib
      implicit none

      type(opal_eos_vintage), intent(inout) :: v

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv, nr = n_eos06_nr, nt = n_eos06_nt

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index_06
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      integer :: species_read_idx, table_var_idx
      double precision :: unused_field
      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '

      if (opal_eos%readco_init_flag(iv_opal06).ne.opal_flag_set) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     v%eos_table(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         do var_idx = 1, mv
            v%eos_output(var_idx) = 1.0d0
         end do
         opal_eos%readco_init_flag(iv_opal06) = opal_flag_set
      end if

      close (2)
! .....read  tables (the file is opened on star%ctrl%iopale by read_controls)
      do x_loop_index_06 = 1, mx
         read (star%ctrl%iopale,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              v%hydrogen_fraction_header(x_loop_index_06), &
              v%z_table(x_loop_index_06), &
              v%moles_per_gram_table(x_loop_index_06), &
              v%mean_molecular_weight_header(x_loop_index_06)
         read (star%ctrl%iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (v%species_fraction_header(x_loop_index_06,species_read_idx), &
              species_read_idx=1,6)
         read (star%ctrl%iopale,'(A)') blank_line
         do density_row = 1, nr
            read (star%ctrl%iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 v%temperature_count_used(x_loop_index_06,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index_06,density_row)
            if (record_number.ne.density_row) then
               write (run_log_unit,'(" OEOS06 Data file incorrect: numtot,jcs= " &
                    &,2I5)') record_number, density_row
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            if (v%temperature_count_used(x_loop_index_06,density_row).lt. &
                 opal06_t6_index_lo(density_row)) then
               write (run_log_unit,'("Problem with OEOS96 data files: X=",F6.4, &
                    &" density=",E14.4)') v%hydrogen_fraction_header(x_loop_index_06), &
                    v%density_grid_table(x_loop_index_06,density_row)
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index_06,density_row)
               if (t6_row.gt.opal06_t6_index_lo(density_row)) then
                  read (star%ctrl%iopale,'(A)') blank_line
                  cycle
               end if
               read (star%ctrl%iopale,'(F11.6,1X,F6.4,E11.4,2E13.6,2E11.3,5F10.6)') &
                    v%t6_list(density_row,t6_row), v%amu_grid(density_row,t6_row), &
                    v%log10_ne_grid(density_row,var_idx), &
                    (v%eos_table(x_loop_index_06,opal_eos_var_order(table_var_idx), &
                    t6_row,density_row), table_var_idx=1,9)
            end do
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
         end do
         read(star%ctrl%iopale,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx).eq.0.0d0) then
            write(run_log_unit,'("READCOEOS06: Error:",I4, &
                 &"-th T6 value is zero")') t6_scan_idx
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
            ! facades stop when their caller passes no ierr.
            ierr = 1
            return
         end if
         v%t6_grid(t6_scan_idx) = v%t6_list(1,t6_scan_idx)
      end do
      do t6_idx = 2, nt
         v%t6_grid_spacing_inv(t6_idx) = 1.0d0/(v%t6_grid(t6_idx) - &
              v%t6_grid(t6_idx-1))
      end do
      v%density_grid(1) = v%density_grid_table(1,1)
      do r_idx = 2, nr
         v%density_grid(r_idx) = v%density_grid_table(1,r_idx)
         v%density_grid_spacing_inv(r_idx) = 1.0d0/(v%density_grid(r_idx) - &
              v%density_grid(r_idx-1))
      end do
      do x_idx = 2, mx
         v%x_grid_spacing_inv(x_idx) = 1.0d0/(v%x_grid_copy(x_idx) - v%x_grid_copy(x_idx-1))
      end do
      close (star%ctrl%iopale)

      return
end subroutine readcoeos06
