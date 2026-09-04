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
subroutine readco(ierr)
      use star_info_lib, only: star

      use opal_eos_lib
      use luout_lib
      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 77, nt = 56

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      double precision :: unused_field

      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '

      if (opal_eos%readco_init_flag.ne.12345678) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     opal_eos%eos_table(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         opal_eos%readco_init_flag = 12345678
      end if

! ..... read  tables (the file is opened on star%ctrl%iopale by read_controls)
      do x_loop_index = 1, mx

         read (star%ctrl%iopale,'(3X,F6.4,3X,F6.4,11X,F10.7,17X,F10.7)') &
              opal_eos%hydrogen_fraction_header(x_loop_index), &
              opal_eos%z_table(x_loop_index), &
              opal_eos%moles_per_gram_table(x_loop_index), &
              opal_eos%mean_molecular_weight_header(x_loop_index)
         read (star%ctrl%iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (opal_eos%species_fraction_header(x_loop_index,var_idx), &
              var_idx=1,6)
         read (star%ctrl%iopale,'(A)') blank_line
         do density_row = 1, nr
            read (star%ctrl%iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 opal_eos%temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 opal_eos%density_grid_table(x_loop_index,density_row)
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
            if (opal_eos%temperature_count_used(x_loop_index,density_row).lt. &
                 opal_eos%t6_index_lo(density_row)) then
               write(run_log_unit,'("PROBLEM WITH DATA FILES: X=",F6.4," DENSITY=", &
                    &E14.4)') opal_eos%hydrogen_fraction_header(x_loop_index), &
                    opal_eos%density_grid_table(x_loop_index,density_row)
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            do t6_row = 1, opal_eos%temperature_count_used(x_loop_index,density_row)
               if (t6_row.gt.opal_eos%t6_index_lo(density_row)) then
                  read (star%ctrl%iopale,'(A)') blank_line
                  cycle
               end if
               read (star%ctrl%iopale,'(F8.4,1X,F6.2,3E13.5,E11.3,6F8.4)') &
                    opal_eos%t6_list(density_row,t6_row), &
                    opal_eos%log10_r_value(density_row,t6_row), &
                    (opal_eos%eos_table(x_loop_index,opal_eos%eos_var_order(var_idx), &
                    t6_row,density_row), var_idx=1,10)
            end do
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
         end do
         read(star%ctrl%iopale,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (opal_eos%t6_list(1,t6_scan_idx).eq.0.0d0) exit
         opal_eos%t6_grid(t6_scan_idx) = opal_eos%t6_list(1,t6_scan_idx)
      end do
      do t6_idx = 2, nt
         opal_eos%t6_grid_spacing_inv(t6_idx) = 1.0d0/(opal_eos%t6_grid(t6_idx) - &
              opal_eos%t6_grid(t6_idx-1))
      end do
      opal_eos%density_grid(1) = opal_eos%density_grid_table(1,1)
      do r_idx = 2, nr
         opal_eos%density_grid(r_idx) = opal_eos%density_grid_table(1,r_idx)
         opal_eos%density_grid_spacing_inv(r_idx) = 1.0d0/(opal_eos%density_grid(r_idx) - &
              opal_eos%density_grid(r_idx-1))
      end do
      do x_idx = 2, mx
         opal_eos%x_grid_spacing_inv(x_idx) = 1.0d0/(opal_eos%x_grid_copy(x_idx) - opal_eos%x_grid_copy(x_idx-1))
      end do

      close (star%ctrl%iopale)
!  MHP 8/98 ADDED RAMP BETWEEN OPAL AND OTHER EOS
!  NEED EDGE OF TABLE AT HIGH RHO, FIXED T.
      opal_eos%t_row_index = 1
      do t6_idx = 1, nt
         opal_eos%density_edge_at_t(t6_idx) = opal_eos%density_grid(opal_eos%density_index_edge_at_t(t6_idx))
      end do
      return
end subroutine readco
