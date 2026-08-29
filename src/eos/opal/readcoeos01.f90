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
subroutine readcoeos01(ierr)
      use star_info_lib, only: star

      use opal_eos_lib
      use luout_lib
      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 191

      character(len=1) :: blank_line

! former common/eeeeos/: this batch's own block, only used by
! readcoeos01.f90, now use-associated from opal_eos_lib.




! --- locals ---
      integer :: x_loop_index_01
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx, fill_idx
      integer :: density_row, t6_row, record_number
      double precision :: unused_field

! density_index_edge_at_t_01's default moved to opal_eos_lib.f90: DATA
! can no longer target it here now that it's use-associated.
      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '

      if (opal_eos%readcoeos01_init_flag.ne.12345678) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     opal_eos%eos_table_01(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         opal_eos%readcoeos01_init_flag = 12345678
      end if

      close (2)
! ..... read  tables
! MHP 8/25 Moved opening of file to parmin
      do x_loop_index_01 = 1, mx

         read (star%ctrl%iopale,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              opal_eos%hydrogen_fraction_header_01(x_loop_index_01), &
              opal_eos%z_table_01(x_loop_index_01), &
              opal_eos%moles_per_gram_table_01(x_loop_index_01), &
              opal_eos%mean_molecular_weight_header_01(x_loop_index_01)
         read (star%ctrl%iopale,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (opal_eos%species_fraction_header_01(x_loop_index_01,var_idx), &
              var_idx=1,6)
         read (star%ctrl%iopale,'(A)') blank_line
         do density_row = 1, nr
            read (star%ctrl%iopale,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 opal_eos%temperature_count_used_01(x_loop_index_01,density_row), &
                 unused_field, unused_field, &
                 opal_eos%density_grid_table_01(x_loop_index_01,density_row)
            if (record_number.ne.density_row) then
               write (run_log_unit,'(" Data file incorrect: numtot,jcs= ",2I5)') &
                    record_number, density_row
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            if (opal_eos%temperature_count_used_01(x_loop_index_01,density_row).lt. &
                 opal_eos%t6_index_lo_01(density_row)) then
               write (run_log_unit,'("problem with data files: X=",F6.4, &
                    &" density=",E14.4)') opal_eos%hydrogen_fraction_header_01(x_loop_index_01), &
                    opal_eos%density_grid_table_01(x_loop_index_01,density_row)
               ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            do t6_row = 1, opal_eos%temperature_count_used_01(x_loop_index_01,density_row)
               if (t6_row.gt.opal_eos%t6_index_lo_01(density_row)) then
                  read (star%ctrl%iopale,'(A)') blank_line
                  cycle
               end if
               read (star%ctrl%iopale,'(F9.5,1X,F6.2,3E13.5,6F8.4)') &
                    opal_eos%t6_list_01(density_row,t6_row), &
                    opal_eos%log10_r_value_01(density_row,t6_row), &
                    (opal_eos%eos_table_01(x_loop_index_01,opal_eos%eos_var_order_01(var_idx), &
                    t6_row,density_row), var_idx=1,9)
            end do
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
            read(star%ctrl%iopale,'(A)') blank_line
         end do
         read(star%ctrl%iopale,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (opal_eos%t6_list_01(1,t6_scan_idx).eq.0.0d0) then
            write(run_log_unit,'("READCOEOS01: Error:",I4, &
                 &"-th T6 value is zero")') t6_scan_idx
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
            ! facades stop when their caller passes no ierr.
            ierr = 1
            return
         end if
         opal_eos%t6_grid_01(t6_scan_idx) = opal_eos%t6_list_01(1,t6_scan_idx)
      end do
      do t6_idx = 2, nt
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!    12 dfs(i)=1D0/(t6a(i)-t6a(i-1))
         opal_eos%t6_grid_spacing_inv_01(t6_idx) = 1.0d0/(opal_eos%t6_grid_01(t6_idx) - &
              opal_eos%t6_grid_01(t6_idx-1))
      end do
      opal_eos%density_grid_01(1) = opal_eos%density_grid_table_01(1,1)
      do r_idx = 2, nr
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       rho(i)=rhogr(1,i)
!    13 dfsr(i)=1D0/(rho(i)-rho(i-1))
         opal_eos%density_grid_01(r_idx) = opal_eos%density_grid_table_01(1,r_idx)
         opal_eos%density_grid_spacing_inv_01(r_idx) = 1.0d0/(opal_eos%density_grid_01(r_idx) - &
              opal_eos%density_grid_01(r_idx-1))
      end do
      do x_idx = 2, mx
         opal_eos%x_grid_spacing_inv_01(x_idx) = 1.0d0/(opal_eos%x_grid_copy_01(x_idx) - opal_eos%x_grid_copy_01(x_idx-1))
      end do

      close (star%ctrl%iopale)
! MHP 7/2003 ADDED RAMP BETWEEN OPAL AND OTHER EOS
! NEED EDGE OF TABLE AT HIGH RHO, FIXED T.
      opal_eos%t_row_index_01 = 1
      do t6_idx = 1, nt
         opal_eos%density_edge_at_t_01(t6_idx) = opal_eos%density_grid_01(opal_eos%density_index_edge_at_t_01(t6_idx))
      end do

      return
end subroutine readcoeos01
