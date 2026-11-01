!----------------------------------------------------------------------
! run_files
!----------------------------------------------------------------------
! The run's files: their fixed unit numbers, the {YREC_...}
! placeholders in the path controls expanded (expand_value), the output
! directory created, the run log opened first, a verbatim copy of the
! inlists (inlist_used), and the other files opened.
module run_files
      use io_units, only: terminal_unit, last_model_unit, run_log_unit, check_io, &
           starting_model_unit, dynamics_unit, mc_summary_unit, mc_model_unit, &
           isochrone_unit, flush_pending_log
      use legacy_controls, only: alex06_table_path, alex95_table_paths, &
           allard_table_path, atm_table_path, centre1_table_path, &
           centre2_table_path, centre3_table_path, centre4_table_path, &
           centre5_table_path, conductive_table_path, control_nml_file, fdyn, &
           fermi_table_path, ffirst, fiso, flast, fshort, history_columns_file, &
           inlist_used_file, isochrone_output_active, kurucz_table2_path, &
           kurucz_table_path, laol_table2_path, laol_table_path, &
           monte_carlo_active, monte_carlo_file1_path, monte_carlo_file2_path, &
           opal2001_eos_table_path, opal2006_eos_table_path, opal92_table2_path, &
           opal92_table_path, opal95_eos_table_path, opal95_table_path, &
           physics_nml_file, profile_columns_file, pure_z_table_path, &
           scv_h_table_path, scv_he_table_path, scv_z_table_path, &
           zams_a_table_path, zams_b_table_path, zams_c_table_path
      implicit none
      private
      public :: assign_unit_numbers, expand_paths, open_run_files, expand_value, &
           output_directory, directory_of, copy_inlists_used, open_run_log_after_error

! The run's output directory, the run log's (open_run_files sets it at
! every read; output_directory returns it).
      character(len=256) :: output_dir = ' '
! Whether this read has expanded the deck's paths (expand_paths), so
! that the run log's path is the deck's (open_run_log_after_error);
! each read starts without (assign_unit_numbers).
      logical :: paths_expanded = .false.

contains

! ---------------------------------------------------------------
! The fixed units written from everywhere: the terminal, the final
! model and the run log (io_units); the read calls this first, before
! its paths are expanded.
subroutine assign_unit_numbers
      terminal_unit = 6
      last_model_unit = 11
      run_log_unit = 20
      paths_expanded = .false.
end subroutine assign_unit_numbers

! ---------------------------------------------------------------
! The {YREC_INPUT}, {YREC_START} and {YREC_OUTPUT} placeholders in
! every path control, expanded (expand_value): every string control
! that names a file, by its canonical name where it has one (after
! adopt_canonical_names), else by its legacy name (test_generated_files.py
! checks the list against the registry).
subroutine expand_paths
      integer :: i

      call expand_value(alex06_table_path)
      call expand_value(allard_table_path)
      call expand_value(atm_table_path)
      call expand_value(conductive_table_path)
      call expand_value(fdyn)
      call expand_value(fermi_table_path)
      call expand_value(ffirst)
      call expand_value(fiso)
      call expand_value(kurucz_table_path)
      call expand_value(kurucz_table2_path)
      call expand_value(laol_table_path)
      call expand_value(laol_table2_path)
      call expand_value(flast)
      call expand_value(opal95_table_path)
      call expand_value(opal92_table_path)
      call expand_value(zams_a_table_path)
      call expand_value(zams_b_table_path)
      call expand_value(zams_c_table_path)
      call expand_value(centre1_table_path)
      call expand_value(centre2_table_path)
      call expand_value(centre3_table_path)
      call expand_value(centre4_table_path)
      call expand_value(centre5_table_path)
      call expand_value(monte_carlo_file1_path)
      call expand_value(monte_carlo_file2_path)
      call expand_value(opal92_table2_path)
      call expand_value(opal95_eos_table_path)
      call expand_value(opal2001_eos_table_path)
      call expand_value(opal2006_eos_table_path)
      call expand_value(pure_z_table_path)
      call expand_value(scv_h_table_path)
      call expand_value(scv_he_table_path)
      call expand_value(scv_z_table_path)
      call expand_value(fshort)
      call expand_value(history_columns_file)
      call expand_value(inlist_used_file)
      call expand_value(profile_columns_file)
      do i = 1, size(alex95_table_paths)
         call expand_value(alex95_table_paths(i))
      end do
      paths_expanded = .true.
end subroutine expand_paths

! ---------------------------------------------------------------
! The files the controls read named, which run_yrec opens after the
! read, in this order: the output directory (the run log's) created and
! the run log opened, so that it holds the read's lines (its header,
! the run plan, its warnings) and any error that follows; the copy of
! the inlists (inlist_used, or inlist_used_file when the deck names
! one); the starting model; the final model (the restart file, and the
! solver's recovery source when a step diverges); and the Monte Carlo
! and isochrone files when those are on. (The history and profile files
! are the output's: yrec_output's output_open.)
subroutine open_run_files(ierr)
      integer, intent(inout) :: ierr
      character(len=256) :: msg, inlist_copy
      integer :: ios

      call start_run_log(ierr)
      if (ierr /= 0) return
      write(terminal_unit,*) "OUTPUT placed in :  ",trim(output_dir)
      write(terminal_unit,*) ''
      inlist_copy = inlist_used_file
      if (len_trim(inlist_copy) == 0) inlist_copy = trim(output_dir) // 'inlist_used'
      call copy_inlists_used([control_nml_file, physics_nml_file], trim(inlist_copy), ierr)
      if (ierr /= 0) return

      call close_run_files
      open(newunit=starting_model_unit, file=ffirst, form='FORMATTED', status='OLD', &
           iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the starting model', ffirst, ierr)
      if (ierr /= 0) return
      open(unit=last_model_unit,file=flast,form='FORMATTED',status='UNKNOWN', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the final model file', flast, ierr)
      if (ierr /= 0) return
      if (monte_carlo_active) then
         open(newunit=dynamics_unit, file=fdyn, form='FORMATTED', status='OLD', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo dynamics file', fdyn, ierr)
         if (ierr /= 0) return
         open(newunit=mc_summary_unit, file=monte_carlo_file1_path, status='UNKNOWN', form='FORMATTED', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo summary file', monte_carlo_file1_path, ierr)
         if (ierr /= 0) return
         open(newunit=mc_model_unit, file=monte_carlo_file2_path, status='UNKNOWN', form='FORMATTED', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo model file', monte_carlo_file2_path, ierr)
         if (ierr /= 0) return
      endif
      if (isochrone_output_active) then
         open(newunit=isochrone_unit, file=fiso, status='UNKNOWN', form='FORMATTED', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the isochrone file', fiso, ierr)
         if (ierr /= 0) return
      endif
end subroutine open_run_files

! ---------------------------------------------------------------
! After a controls error (read_controls' ierr): the run log alone, so
! that it holds the read's lines and the error, once the read has
! expanded the deck's paths. A read that stopped before (the deck
! missing or unreadable) has no run log of the deck's: its messages
! stay on the terminal.
subroutine open_run_log_after_error
      integer :: ierr

      if (.not. paths_expanded) return
      ierr = 0
      call start_run_log(ierr)
end subroutine open_run_log_after_error

! ---------------------------------------------------------------
! The output directory, the run log's, created (with its parents), and
! the run log opened.
subroutine start_run_log(ierr)
      integer, intent(inout) :: ierr

      output_dir = directory_of(fshort)
      call make_directory(output_dir)
      call open_run_log(fshort, ierr)
end subroutine start_run_log

! ---------------------------------------------------------------
! The files an earlier read opened (an in-process rerun), closed before
! they are opened again.
subroutine close_run_files
      call close_unit(starting_model_unit)
      call close_unit(dynamics_unit)
      call close_unit(mc_summary_unit)
      call close_unit(mc_model_unit)
      call close_unit(isochrone_unit)
contains
subroutine close_unit(unit)
      integer, intent(inout) :: unit
      logical :: is_open
      integer :: ios
      if (unit == -1) return
      inquire(unit=unit, opened=is_open, iostat=ios)
      if (ios == 0 .and. is_open) close(unit)
      unit = -1
end subroutine close_unit
end subroutine close_run_files

! ---------------------------------------------------------------
! The run's output directory: the run log's (blank for the working
! directory).
function output_directory() result(dir)
      character(len=256) :: dir
      dir = output_dir
end function output_directory

! ---------------------------------------------------------------
! The run log on run_log_unit, from log_output_file: an explicit .log
! name is kept, a legacy .short name becomes .log, and .log is appended
! to any other (so the default run.log stays run.log).
subroutine open_run_log(log_output_file, ierr)
      character(len=*), intent(in) :: log_output_file
      integer, intent(inout) :: ierr
      character(len=256) :: log_path, msg
      integer :: n, ios

      n = len_trim(log_output_file)
      if (n >= 4 .and. log_output_file(max(1,n-3):n) == '.log') then
         log_path = log_output_file
      else
         n = index(log_output_file, '.short')
         if (n > 0) then
            log_path = log_output_file(1:n) // 'log'
         else
            log_path = trim(log_output_file) // '.log'
         end if
      end if
      open(run_log_unit, file=log_path, form='FORMATTED', &
           status='REPLACE', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the run log', log_path, ierr)
      if (ierr == 0) call flush_pending_log
end subroutine open_run_log

! ---------------------------------------------------------------
! The directory part of a path, up to and including its last '/';
! blank for a bare file name.
pure function directory_of(path) result(dir)
      character(len=*), intent(in) :: path
      character(len=len(path)) :: dir
      integer :: islash

      islash = index(path, '/', back=.true.)
      dir = ' '
      if (islash > 0) dir = path(1:islash)
end function directory_of

! ---------------------------------------------------------------
! A directory and its parents created, as mkdir -p does, through the C
! library's mkdir (no shell). A directory that exists already is left
! alone, and a failure is left to the first open in it to report.
subroutine make_directory(dir)
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
      character(len=*), intent(in) :: dir
      interface
         integer(c_int) function c_mkdir(path, mode) bind(c, name='mkdir')
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
         end function c_mkdir
      end interface
! rwx for all, less the umask (what mkdir -p gives)
      integer(c_int), parameter :: mode_rwx_all = int(o'777', c_int)
      integer :: i, n
      integer(c_int) :: rc

      n = len_trim(dir)
      do i = 2, n
         if (dir(i:i) == '/') rc = c_mkdir(dir(1:i-1) // c_null_char, mode_rwx_all)
      end do
      if (n > 0) then
         if (dir(n:n) /= '/') rc = c_mkdir(dir(1:n) // c_null_char, mode_rwx_all)
      end if
end subroutine make_directory

! ---------------------------------------------------------------
! A verbatim copy of the deck's files (sources: the CONTROL and PHYSICS
! files, the same file twice for an inlist, copied once) at dest, in the
! output directory: the run's record of its settings. A source that
! cannot be opened is skipped; a failed write is an error
! (output_files).
subroutine copy_inlists_used(sources, dest, ierr)
      use output_files, only: open_output_file, close_output_file
      character(len=*), intent(in) :: sources(2), dest
      integer, intent(inout) :: ierr
      integer :: src_unit, dst_unit
      character(len=4096) :: line
      character(len=256) :: msg
      integer :: ios, read_ios, nfile

      call open_output_file(dest, dst_unit, ierr)
      if (ierr /= 0) return
      ios = 0
      do nfile = 1, 2
         if (nfile == 2 .and. trim(sources(2)) == trim(sources(1))) exit
         open(newunit=src_unit, file=sources(nfile), form='FORMATTED', &
              status='OLD', iostat=read_ios)
         if (read_ios /= 0) cycle
         write(dst_unit,'(2a)',iostat=ios,iomsg=msg) '! ==== copied verbatim from: ', &
              trim(sources(nfile))
         do while (ios == 0)
            read(src_unit,'(a)',iostat=read_ios) line
            if (read_ios /= 0) exit
            write(dst_unit,'(a)',iostat=ios,iomsg=msg) trim(line)
         end do
         close(src_unit)
         if (ios /= 0) exit
      end do
      call close_output_file(dst_unit, dest, ios, msg, ierr)
end subroutine copy_inlists_used

! ---------------------------------------------------------------
! A leading {YREC_INPUT}, {YREC_START} or {YREC_OUTPUT} in path_value,
! replaced by that environment variable's value, or by its default
! (relative to the deck's directory) when the variable is unset or
! empty. A value that itself starts with a later placeholder is
! expanded again.
subroutine expand_value(path_value)
      character(len=*), intent(inout) :: path_value
      character(len=*), parameter :: placeholders(3) = [character(len=13) :: &
           '{YREC_INPUT}', '{YREC_START}', '{YREC_OUTPUT}']
      character(len=*), parameter :: defaults(3) = [character(len=17) :: &
           '../../input', '../../startmodels', 'output']
      character(len=1024) :: value
      integer :: i, n, value_len, status

      do i = 1, size(placeholders)
         n = len_trim(placeholders(i))
         if (len(path_value) < n) cycle
         if (path_value(1:n) /= placeholders(i)(1:n)) cycle
         call get_environment_variable(placeholders(i)(2:n-1), value, value_len, status)
         if (status == 1 .or. len_trim(value) == 0) value = defaults(i)
         path_value = trim(value) // path_value(n+1:)
      end do
end subroutine expand_value

end module run_files
