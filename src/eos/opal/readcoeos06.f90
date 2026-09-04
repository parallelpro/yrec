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
! opal_eos%eos_table_06 and resets opal_eos%eos_output_06 to 1.0), is
! reused, stale, as the column index into log10_ne_grid_06 in the
! per-row READ statement further down. Because the init loops leave
! var_idx at mv+1, every row's read of
! log10_ne_grid_06(density_row, var_idx) overwrites the SAME column
! (var_idx == mv+1) instead of tracking t6_row as the analogous
! opal_eos%t6_list_06/amu_grid_06 reads do. This looks like an original
! bug (using "j" instead of "i" in that READ's variable list) but is
! transliterated exactly rather than "fixed"; log10_ne_grid_06 is never
! read anywhere.
subroutine readcoeos06(ierr)
      use star_info_lib, only: star

      use opal_eos_lib
      use luout_lib
      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197

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

      if (opal_eos%readcoeos06_init_flag.ne.12345678) then
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
            opal_eos%eos_output_06(var_idx) = 1.0d0
         end do
         opal_eos%readcoeos06_init_flag = 12345678
      end if

      close (2)
! .....read  tables (the file is opened on star%ctrl%iopale by read_controls)
      do x_loop_index_06 = 1, mx
         read (star%ctrl%iopale,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              opal_eos%hydrogen_fraction_header_06(x_loop_index_06), &
              opal_eos%z_table_06(x_loop_index_06), &
              opal_eos%moles_per_gram_table_06(x_loop_index_06), &
              opal_eos%mean_molecular_weight_header_06(x_loop_index_06)
         read (star%ctrl%iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (opal_eos%species_fraction_header_06(x_loop_index_06,species_read_idx), &
              species_read_idx=1,6)
         read (star%ctrl%iopale,'(A)') blank_line
         do density_row = 1, nr
            read (star%ctrl%iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 opal_eos%temperature_count_used_06(x_loop_index_06,density_row), &
                 unused_field, unused_field, &
                 opal_eos%density_grid_table_06(x_loop_index_06,density_row)
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
            if (opal_eos%temperature_count_used_06(x_loop_index_06,density_row).lt. &
                 opal_eos%t6_index_lo_06(density_row)) then
               write (run_log_unit,'("Problem with OEOS96 data files: X=",F6.4, &
                    &" density=",E14.4)') opal_eos%hydrogen_fraction_header_06(x_loop_index_06), &
                    opal_eos%density_grid_table_06(x_loop_index_06,density_row)
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            do t6_row = 1, opal_eos%temperature_count_used_06(x_loop_index_06,density_row)
               if (t6_row.gt.opal_eos%t6_index_lo_06(density_row)) then
                  read (star%ctrl%iopale,'(A)') blank_line
                  cycle
               end if
               read (star%ctrl%iopale,'(F11.6,1X,F6.4,E11.4,2E13.6,2E11.3,5F10.6)') &
                    opal_eos%t6_list_06(density_row,t6_row), opal_eos%amu_grid_06(density_row,t6_row), &
                    opal_eos%log10_ne_grid_06(density_row,var_idx), &
                    (opal_eos%eos_table_06(x_loop_index_06,opal_eos%eos_var_order_06(table_var_idx), &
                    t6_row,density_row), table_var_idx=1,9)
            end do
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
         end do
         read(star%ctrl%iopale,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (opal_eos%t6_list_06(1,t6_scan_idx).eq.0.0d0) then
            write(run_log_unit,'("READCOEOS06: Error:",I4, &
                 &"-th T6 value is zero")') t6_scan_idx
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
            ! facades stop when their caller passes no ierr.
            ierr = 1
            return
         end if
         opal_eos%t6_grid_06(t6_scan_idx) = opal_eos%t6_list_06(1,t6_scan_idx)
      end do
      do t6_idx = 2, nt
         opal_eos%t6_grid_spacing_inv_06(t6_idx) = 1.0d0/(opal_eos%t6_grid_06(t6_idx) - &
              opal_eos%t6_grid_06(t6_idx-1))
      end do
      opal_eos%density_grid_06(1) = opal_eos%density_grid_table_06(1,1)
      do r_idx = 2, nr
         opal_eos%density_grid_06(r_idx) = opal_eos%density_grid_table_06(1,r_idx)
         opal_eos%density_grid_spacing_inv_06(r_idx) = 1.0d0/(opal_eos%density_grid_06(r_idx) - &
              opal_eos%density_grid_06(r_idx-1))
      end do
      do x_idx = 2, mx
         opal_eos%x_grid_spacing_inv_06(x_idx) = 1.0d0/(opal_eos%x_grid_copy_06(x_idx) - opal_eos%x_grid_copy_06(x_idx-1))
      end do
      close (star%ctrl%iopale)

      return
end subroutine readcoeos06
