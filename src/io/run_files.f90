!----------------------------------------------------------------------
! run_files
!----------------------------------------------------------------------
! The run's files: their fixed unit numbers, the {YREC_...}
! placeholders in the path controls expanded (expand_value), the output
! directory created with a verbatim copy of the inlists (inlist_used),
! and the files opened.
module run_files
      use io_units, only: terminal_unit, last_model_unit, run_log_unit, check_io
      use legacy_controls
      implicit none
      private
      public :: assign_unit_numbers, expand_paths, open_run_files, expand_value, &
           output_directory

! The run's output directory, the run log's (open_run_files sets it at
! every read; output_directory returns it).
      character(len=256) :: output_dir = ' '

contains

! ---------------------------------------------------------------
! The fixed unit numbers of the run's files and tables (the newer
! readers and writers use newunit instead).
subroutine assign_unit_numbers
      terminal_unit = 6
      last_model_unit = 11
      first_unit = 12
! the physics and control namelists
      run_unit = 13
      standard_unit = 14
      fermi_unit = 15
      run_log_unit = 20
      dynamics_unit = 30
      opal92_table_unit = 32
      opal95_table_unit = 48
      composition_unit = 34
      kurucz_table_unit = 36
      isochrone_file_unit = 37
      alex94_table_unit = 39
! the MHD EOS tables (mhdst treats a unit <= 0 as a missing table)
      unit_zams_a = 40
      unit_zams_b = 41
      unit_zams_c = 42
      unit_centre1 = 43
      unit_centre2 = 44
      unit_centre3 = 45
      unit_centre4 = 46
      unit_centre5 = 47
! the OPAL EOS
      iopale = 49
      kurucz_table2_unit = 65
      allard_table_unit = 66
      monte_carlo_unit1 = 70
      monte_carlo_unit2 = 71
      scv_h_unit = 72
      scv_he_unit = 73
      scv_z_unit = 74
! the Potekhin conductive opacities
      icondopacp = 75
      alex06_table_unit = 90
end subroutine assign_unit_numbers

! ---------------------------------------------------------------
! The {YREC_INPUT}, {YREC_START} and {YREC_OUTPUT} placeholders in
! every path control, expanded (expand_value): every string control
! that names a file, by its legacy name (test_generated_files.py checks
! the list against the registry).
subroutine expand_paths
      integer :: i

      call expand_value(falex06)
      call expand_value(fallard)
      call expand_value(fatm)
      call expand_value(fcondopacp)
      call expand_value(fdyn)
      call expand_value(ffermi)
      call expand_value(ffirst)
      call expand_value(fiso)
      call expand_value(fkur)
      call expand_value(fkur2)
      call expand_value(flaol)
      call expand_value(flaol2)
      call expand_value(flast)
      call expand_value(fliv95)
      call expand_value(flldat)
      call expand_value(fmhd1)
      call expand_value(fmhd2)
      call expand_value(fmhd3)
      call expand_value(fmhd4)
      call expand_value(fmhd5)
      call expand_value(fmhd6)
      call expand_value(fmhd7)
      call expand_value(fmhd8)
      call expand_value(fmonte1)
      call expand_value(fmonte2)
      call expand_value(fopal2)
      call expand_value(fopale)
      call expand_value(fopale01)
      call expand_value(fopale06)
      call expand_value(fpurez)
      call expand_value(fscomp)
      call expand_value(fscvh)
      call expand_value(fscvhe)
      call expand_value(fscvz)
      call expand_value(fshort)
      call expand_value(fstch)
      call expand_value(history_columns_file)
      call expand_value(inlist_used_file)
      call expand_value(profile_columns_file)
      do i = 1, size(opecalex)
         call expand_value(opecalex(i))
      end do
end subroutine expand_paths

! ---------------------------------------------------------------
! The output directory (the run log's directory) created and the
! inlists copied into it, then the run's files opened: the starting
! model, the final model (the restart file, and the solver's recovery
! source when a step diverges), the run log, and the Monte Carlo and
! isochrone files when those are on. (The history and profile files
! are the output's: yrec_output's output_open, after the read. s is
! unused.)
subroutine open_run_files(s, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
      integer, intent(inout) :: ierr
      character(len=256) :: msg
      integer :: ios

      output_dir = directory_of(fshort)
      print *,"OUTPUT placed in :  ",trim(output_dir)
      print *, ''
      call make_directory(output_dir)
      call copy_inlists_used(trim(inlist_used_file), ierr)
      if (ierr /= 0) return

      open(unit=first_unit,file=ffirst,form='FORMATTED',status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the starting model', ffirst, ierr)
      if (ierr /= 0) return
      open(unit=last_model_unit,file=flast,form='FORMATTED',status='UNKNOWN', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the final model file', flast, ierr)
      if (ierr /= 0) return
      call open_run_log(fshort, ierr)
      if (ierr /= 0) return
      if (monte_carlo_active) then
         open(unit=dynamics_unit,file=fdyn,form='FORMATTED',status='OLD', iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo dynamics file', fdyn, ierr)
         if (ierr /= 0) return
         open(monte_carlo_unit1, file=fmonte1,status='UNKNOWN',form='FORMATTED', iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo summary file', fmonte1, ierr)
         if (ierr /= 0) return
         open(monte_carlo_unit2, file=fmonte2,status='UNKNOWN',form='FORMATTED', iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the Monte Carlo model file', fmonte2, ierr)
         if (ierr /= 0) return
      endif
      if (isochrone_output_active) then
         open(isochrone_file_unit, file=fiso,status='UNKNOWN', form='FORMATTED', iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'open the isochrone file', fiso, ierr)
         if (ierr /= 0) return
      endif
end subroutine open_run_files

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
! A verbatim copy of the deck's file(s) in the output directory: the
! run's record of its settings.
subroutine copy_inlists_used(dest, ierr)
      character(len=*), intent(in) :: dest
      integer, intent(inout) :: ierr
      integer :: src_unit, dst_unit
      character(len=4096) :: line
      character(len=256) :: msg
      integer :: ios, nfile
      character(len=256) :: sources(2)
      sources(1) = control_nml_file
      sources(2) = physics_nml_file
      open(newunit=dst_unit, file=dest, &
           form='FORMATTED', status='UNKNOWN', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'write the copy of the inlist', dest, ierr)
      if (ierr /= 0) return
      do nfile = 1, 2
         if (nfile == 2 .and. trim(sources(2)) == trim(sources(1))) exit
         open(newunit=src_unit, file=sources(nfile), form='FORMATTED', &
              status='OLD', iostat=ios)
         if (ios /= 0) cycle
         write(dst_unit,'(2a)') '! ==== copied verbatim from: ', &
              trim(sources(nfile))
         do
            read(src_unit,'(a)',iostat=ios) line
            if (ios /= 0) exit
            write(dst_unit,'(a)') trim(line)
         end do
         close(src_unit)
      end do
      close(dst_unit)
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
