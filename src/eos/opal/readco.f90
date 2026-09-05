!----------------------------------------------------------------------
! readco
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original readco.f; only variable names, source form, and comment
! style were updated.
!
! Reads the OPAL 1995 EOS data tables (one block per tabulated X
! value) into opal_eos_lib's opal_eos state and derives the auxiliary
! grid-spacing arrays used by esac.f90/t6rinterp.f90 for
! interpolation. Called once from esac.f90 (guarded by
! opal_eos%table_loaded_flag).
subroutine readco(v, ierr)
      use star_info_lib, only: star

      use opal_eos_lib
      use luout_lib
      implicit none

      type(opal_eos_vintage95), intent(inout) :: v

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv, nr = n_eos95_nr, nt = n_eos95_nt

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      double precision :: unused_field

      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '

      if (opal_eos%readco_init_flag(iv_opal95).ne.opal_flag_set) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     v%eos_table(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         opal_eos%readco_init_flag(iv_opal95) = opal_flag_set
      end if

! ..... read  tables (the file is opened on star%ctrl%iopale by read_controls)
      do x_loop_index = 1, mx

         read (star%ctrl%iopale,'(3X,F6.4,3X,F6.4,11X,F10.7,17X,F10.7)') &
              v%hydrogen_fraction_header(x_loop_index), &
              v%z_table(x_loop_index), &
              v%moles_per_gram_table(x_loop_index), &
              v%mean_molecular_weight_header(x_loop_index)
         read (star%ctrl%iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (v%species_fraction_header(x_loop_index,var_idx), &
              var_idx=1,6)
         read (star%ctrl%iopale,'(A)') blank_line
         do density_row = 1, nr
            read (star%ctrl%iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 v%temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index,density_row)
            if (record_number.ne.density_row) then
               write(run_log_unit,'(" DATA FILE INCORRECT: NUMTOT,JCS= ",2I5)') &
                    record_number, density_row
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            if (v%temperature_count_used(x_loop_index,density_row).lt. &
                 opal95_t6_index_lo(density_row)) then
               write(run_log_unit,'("PROBLEM WITH DATA FILES: X=",F6.4," DENSITY=", &
                    &E14.4)') v%hydrogen_fraction_header(x_loop_index), &
                    v%density_grid_table(x_loop_index,density_row)
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index,density_row)
               if (t6_row.gt.opal95_t6_index_lo(density_row)) then
                  read (star%ctrl%iopale,'(A)') blank_line
                  cycle
               end if
               read (star%ctrl%iopale,'(F8.4,1X,F6.2,3E13.5,E11.3,6F8.4)') &
                    v%t6_list(density_row,t6_row), &
                    v%log10_r_value(density_row,t6_row), &
                    (v%eos_table(x_loop_index,opal_eos_var_order(var_idx), &
                    t6_row,density_row), var_idx=1,10)
            end do
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
         end do
         read(star%ctrl%iopale,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx).eq.0.0d0) exit
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
!  MHP 8/98 ADDED RAMP BETWEEN OPAL AND OTHER EOS
!  NEED EDGE OF TABLE AT HIGH RHO, FIXED T.
      v%t_row_index = 1
      do t6_idx = 1, nt
         v%density_edge_at_t(t6_idx) = v%density_grid(opal95_density_index_edge(t6_idx))
      end do
      return
end subroutine readco
