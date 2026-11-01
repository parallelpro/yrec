!----------------------------------------------------------------------
! run_files
!----------------------------------------------------------------------
! The run's files: their fixed unit numbers, the {YREC_...}
! placeholders a path control may start with (expand_value, which the
! controls read applies to every path control), the output directory
! created, the run log opened first, a verbatim copy of the deck
! (inlist_used), and the other files opened. The files' paths are
! controls the star holds (star%control's first_model_file, ...,
! star%physics%log_output_file), so nothing here keeps state from the
! read.
module run_files
      use run_messages, only: terminal_unit, run_log_unit, check_io, flush_pending_log
      use star_info_lib, only: star_info
      implicit none
      private
      public :: assign_unit_numbers, open_run_files, expand_value, &
           output_directory, directory_of, copy_inlists_used, open_run_log_after_error
      public :: last_model_unit, starting_model_unit, dynamics_unit, mc_summary_unit, &
           mc_model_unit, isochrone_unit

! The units of the run's files besides the run log and the terminal
! (run_messages'): the final model, written from everywhere (assigned
! with them by assign_unit_numbers), and the starting model, the Monte
! Carlo samples, summary and model files and the isochrone file, which
! open_run_files opens with newunit (-1 until then):
!   last_model_unit     -- the final .mod model (last_model_file)
!   starting_model_unit -- the starting model (first_model_file)
!   dynamics_unit       -- the Monte Carlo samples (dynamics_file)
!   mc_summary_unit, mc_model_unit -- the Monte Carlo summary and model
!                          files (monte_carlo_file1_path, _file2_path)
!   isochrone_unit      -- the isochrone file (isochrone_file)
! The table readers open their tables themselves, with newunit.
      integer :: last_model_unit
      integer :: starting_model_unit = -1, dynamics_unit = -1, &
           mc_summary_unit = -1, mc_model_unit = -1, isochrone_unit = -1

contains

! ---------------------------------------------------------------
! The fixed units written from everywhere: the terminal and the run log
! (run_messages'), and the final model; the read calls this first,
! before its paths are expanded.
subroutine assign_unit_numbers
      terminal_unit = 6
      last_model_unit = 11
      run_log_unit = 20
end subroutine assign_unit_numbers

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
!   star reads: control%control_nml_file, control%dynamics_file,
!     control%first_model_file, control%inlist_used_file,
!     control%isochrone_file, control%isochrone_output_active,
!     control%last_model_file, control%monte_carlo_file1_path,
!     control%monte_carlo_file2_path, control%physics_nml_file,
!     physics%log_output_file, physics%monte_carlo_active
!   star writes: none
subroutine open_run_files(s, ierr)
      type(star_info), intent(in) :: s
      integer, intent(inout) :: ierr
      character(len=256) :: msg, inlist_copy, output_dir
      integer :: ios

      call start_run_log(s%physics%log_output_file, ierr)
      if (ierr /= 0) return
      output_dir = output_directory(s)
      write(terminal_unit,*) "OUTPUT placed in :  ",trim(output_dir)
      write(terminal_unit,*) ''
      inlist_copy = s%control%inlist_used_file
      if (len_trim(inlist_copy) == 0) inlist_copy = trim(output_dir) // 'inlist_used'
      call copy_inlists_used([s%control%control_nml_file, s%control%physics_nml_file], &
           trim(inlist_copy), ierr)
      if (ierr /= 0) return

      call close_run_files
      open(newunit=starting_model_unit, file=s%control%first_model_file, form='FORMATTED', status='OLD', &
           iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the starting model', s%control%first_model_file, ierr)
      if (ierr /= 0) return
      open(unit=last_model_unit,file=s%control%last_model_file,form='FORMATTED',status='UNKNOWN', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the final model file', s%control%last_model_file, ierr)
      if (ierr /= 0) return
      if (s%physics%monte_carlo_active) then
         open(newunit=dynamics_unit, file=s%control%dynamics_file, form='FORMATTED', status='OLD', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo dynamics file', s%control%dynamics_file, ierr)
         if (ierr /= 0) return
         open(newunit=mc_summary_unit, file=s%control%monte_carlo_file1_path, status='UNKNOWN', form='FORMATTED', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo summary file', s%control%monte_carlo_file1_path, ierr)
         if (ierr /= 0) return
         open(newunit=mc_model_unit, file=s%control%monte_carlo_file2_path, status='UNKNOWN', form='FORMATTED', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo model file', s%control%monte_carlo_file2_path, ierr)
         if (ierr /= 0) return
      endif
      if (s%control%isochrone_output_active) then
         open(newunit=isochrone_unit, file=s%control%isochrone_file, status='UNKNOWN', form='FORMATTED', &
              iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the isochrone file', s%control%isochrone_file, ierr)
         if (ierr /= 0) return
      endif
end subroutine open_run_files

! ---------------------------------------------------------------
! After a controls error: the run log alone, at the deck's
! log_output_file (expanded), so that it holds the read's lines and the
! error. The controls read calls it when check_controls stops it, the
! one error after the deck's paths are expanded (the records are not
! stored after an error, so the read passes its own value). A read that
! stopped before (the deck missing or unreadable) opens no run log: its
! messages stay on the terminal.
subroutine open_run_log_after_error(log_output_file)
      character(len=*), intent(in) :: log_output_file
      integer :: ierr

      ierr = 0
      call start_run_log(log_output_file, ierr)
end subroutine open_run_log_after_error

! ---------------------------------------------------------------
! The output directory, the run log's, created (with its parents), and
! the run log opened.
subroutine start_run_log(log_output_file, ierr)
      character(len=*), intent(in) :: log_output_file
      integer, intent(inout) :: ierr

      call make_directory(directory_of(log_output_file))
      call open_run_log(log_output_file, ierr)
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
!   star reads: physics%log_output_file
!   star writes: none
function output_directory(s) result(dir)
      type(star_info), intent(in) :: s
      character(len=256) :: dir
      dir = directory_of(s%physics%log_output_file)
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
