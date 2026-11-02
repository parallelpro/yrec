!----------------------------------------------------------------------
! test_io
!----------------------------------------------------------------------
! io's own procedures against their contracts, with no deck and no
! tables: every file a check reads is written here, in the working
! directory. Checks:
!   paths        expand_value: a leading {YREC_INPUT}, {YREC_START} or
!                {YREC_OUTPUT} becomes the environment variable's value
!                (set here with the C library's setenv), or its default
!                (../../input, ../../startmodels, output) when the
!                variable is unset or empty; a placeholder that is not
!                leading, a value without one and a value shorter than
!                any placeholder stay as they are; a value that expands
!                to a later placeholder is expanded again;
!   directories  directory_of: the path up to and including its last
!                '/', blank for a bare file name;
!   pulse format pulse_format_code: GYRE, FGONG and GSM in any case and
!                with blanks around them; anything else is unknown;
!   fields       integer_field: a number that fits is right-justified in
!                its width, as an I<width> edit writes it; a wider one
!                (a million models in a field of 6) is written in full
!                after a blank; interior_gyre_path: the model number
!                zero-padded to five digits, in full from 100000 on;
!   columns      parse_columns on a four-column table: a blank file name
!                selects the default names, in their order; a columns
!                file selects its names in file order, past '!' comments
!                (whole-line and trailing) and blank lines; an unknown
!                name, a file that selects nothing, a missing file and
!                more than max_cols names are errors (ierr = 1) with
!                their messages, the unknown name's with the valid names
!                in the run log;
!   groups       nml_file_has_group: a group is found in any case, with
!                a blank, a tab or '/' after its name, and not as the
!                prefix of a longer name or inside a comment; a missing
!                file has none;
!   .mod         write_mod_model then read_mod_model gives back every
!                value it holds bit for bit (the globals, the flags and
!                every shell's columns, at values from 1e-310 to 1e120),
!                and the file has the YMOD 1 columns line with a row in
!                its order;
!                an unknown global is noted in the run log and skipped;
!                replace keeps one model in the file, append adds one;
!                another columns line, a file cut inside its shells and
!                a zone count above json are errors with their messages;
!   format tag   read_model_format_tag reads the first four characters;
!                an empty file is an error with its message;
!   YREC7        read_yrec7 on a six-shell file written here in the
!                layout model_files' header gives: the header values
!                (MODEL2's at their defaults); X
!                and Z in units of 1e-6; He3-O18 per radiative shell and
!                once per convection zone, copied through the zone; He4
!                = 1 - X - Z - He3; H2-Be9 two shells per line, the odd
!                one on its own line; omega from -log10(omega), zero at
!                the 58.9 sentinel; the triangle reflected by a negative
!                record number, and not without one; luminosities in
!                erg/s converted to Lsun, in Lsun kept;
!                a rotating file without its omega records reads with
!                omega zero; a file cut inside its shells, a zone count
!                of 0, a shell mass above the star's and an envelope
!                convection zone outside the model are errors;
!   MODEL2       read_model2 on a four-shell file: the header values,
!                the codes and flags, the luminosities, the envelope
!                records, every shell's columns; a file cut inside its
!                shells and a zone count above json are errors;
!   tables       history_columns and profile_columns: every name non-
!                blank, shorter than column_name_len and unique, every
!                kind real_column or integer_column; each default
!                selection (blank columns file) resolves, in the .list
!                order;
!   writes       a failed write is an error (ierr = 1) with its message
!                ("cannot write ..."): a write to a unit open for
!                reading only (which gfortran reports) for the files
!                kept open, the .mod, a history row and the isochrone
!                file's header; a file whose directory cannot be made;
!                and a full disk, as a file-size limit (RLIMIT_FSIZE)
!                stands in for one, for the files io closes, the GYRE,
!                FGONG and profile files and the copy of the deck (at
!                their close), and the history file (at its row), which
!                end up shorter on disk than what was written (gfortran's
!                runtime reports no error on the writes themselves,
!                io/output_files.f90; a runtime that does reports it at
!                the write, with the same message head), and for the
!                files open for the run, the .mod model and the isochrone
!                file at their write and the run log at check_run_log (a
!                run log rewound and written again, shorter, is not cut
!                short);
!                two runs that write one path (M174) are not a full
!                disk: a history file that another run truncates in
!                place after this run's last row (the C library's fopen
!                "w", the open gfortran's status='replace' makes) or
!                replaces by a new file, and a file io closes whose path
!                another file has taken, close without error (the file
!                measured is the unit's own); and the size check keeps
!                the unit's place (a line written after another writer
!                lengthened the file follows the unit's last line);
!   FGONG        write_fgong_pulse's globals 4-6 are the first card's
!                starting Z and X and the mixing length (MESA's initial
!                Z, X and alpha), whatever the current card's are.
! The expected values are formed here: through the text of the layouts'
! formats (expected = what the text written here reads back as, by an
! internal read), never by the reader under test. Each check prints ok
! or FAIL; the program stops with status 1 on any FAIL. Each check was
! confirmed to fail under a mutation of the code it covers: the first
! ones on 2026-10-29, the writes' and FGONG's when they came (2026-10-31),
! M174's on 2026-11-01, and the fields' and the full-disk checks of the
! .mod, the isochrone file and the run log on 2026-11-01.
! The messages io reports go to the files terminal.txt and run_log.txt
! (run_messages' terminal_unit and run_log_unit), where the error checks
! look for them.
program test_io
      use, intrinsic :: iso_fortran_env, only: dp => real64, int64
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char, c_int64_t, &
           c_funptr, c_intptr_t, c_null_funptr, c_ptr, c_associated, c_new_line
      use star_info_lib, only: star_info, json, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, i_li6, i_li7, &
           i_be9
      use run_messages, only: run_log_unit, terminal_unit, log_line
      use run_files, only: expand_value, directory_of, copy_inlists_used, isochrone_unit, &
           check_run_log
      use yrec_output, only: pulse_format_code, unknown_pulse_format, &
           pulse_gyre, pulse_fgong, pulse_gsm, output_run_header, interior_gyre_path
      use run_cards, only: card_evolve
      use output_columns, only: output_column, real_column, integer_column, &
           parse_columns, max_cols, column_name_len, hist_default_names, &
           prof_default_names, integer_field
      use controls_reader, only: nml_file_has_group
      use model_files, only: model_file_header, write_mod_model, read_mod_model, &
           read_yrec7, read_model2, read_model_format_tag
      use history_output, only: history_columns, history_file, history_output_init, &
           write_history_row, close_history
      use profile_output, only: profile_selection, profile_output_init, write_profile
      use pulsation_output, only: write_gyre_pulse, write_fgong_pulse
      use output_files, only: open_output_file, flush_output_file, close_output_file
      use stitched_model, only: profile_columns, n_prof_cols, n_pulse_cols
      implicit none

! A process's file-size limits (struct rlimit; rlim_t is 64 bits on
! 64-bit Linux and macOS), RLIMIT_FSIZE and SIGXFSZ (the same numbers
! on both), and SIG_IGN, the handler that ignores a signal.
      type, bind(c) :: file_size_limits
         integer(c_int64_t) :: soft, hard
      end type file_size_limits
      integer(c_int), parameter :: rlimit_fsize = 1, sigxfsz = 25
      integer(c_intptr_t), parameter :: sig_ign = 1

      interface
         integer(c_int) function c_setenv(name, value, overwrite) bind(c, name='setenv')
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: name(*), value(*)
            integer(c_int), value :: overwrite
         end function c_setenv
         integer(c_int) function c_unsetenv(name) bind(c, name='unsetenv')
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: name(*)
         end function c_unsetenv
         integer(c_int) function c_getrlimit(resource, limits) bind(c, name='getrlimit')
            import :: c_int, file_size_limits
            integer(c_int), value :: resource
            type(file_size_limits), intent(out) :: limits
         end function c_getrlimit
         integer(c_int) function c_setrlimit(resource, limits) bind(c, name='setrlimit')
            import :: c_int, file_size_limits
            integer(c_int), value :: resource
            type(file_size_limits), intent(in) :: limits
         end function c_setrlimit
         type(c_funptr) function c_signal(signal, handler) bind(c, name='signal')
            import :: c_int, c_funptr
            integer(c_int), value :: signal
            type(c_funptr), value :: handler
         end function c_signal
         type(c_ptr) function c_fopen(path, mode) bind(c, name='fopen')
            import :: c_ptr, c_char
            character(kind=c_char), intent(in) :: path(*), mode(*)
         end function c_fopen
         integer(c_int) function c_fputs(text, stream) bind(c, name='fputs')
            import :: c_int, c_char, c_ptr
            character(kind=c_char), intent(in) :: text(*)
            type(c_ptr), value :: stream
         end function c_fputs
         integer(c_int) function c_fclose(stream) bind(c, name='fclose')
            import :: c_int, c_ptr
            type(c_ptr), value :: stream
         end function c_fclose
         integer(c_int) function c_rename(old_path, new_path) bind(c, name='rename')
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: old_path(*), new_path(*)
         end function c_rename
      end interface

! the files the checks write and read, in the working directory
      character(len=*), parameter :: terminal_file = 'terminal.txt', &
           run_log_file = 'run_log.txt'

! The synthetic starting models: a YREC7 file of six shells (core
! convection zone top y7_jcore, envelope convection zone bottom y7_jenv)
! and a MODEL2 file of four.
      integer, parameter :: y7_nz = 6, y7_jcore = 2, y7_jenv = 4, m2_nz = 4
      type :: yrec7_fixture
         real(dp) :: log_mass(y7_nz), log_r(y7_nz), lum(y7_nz), log_p(y7_nz), &
              log_t(y7_nz), log_rho(y7_nz), cno(4:11,y7_nz), light(12:15,y7_nz), &
              omega_log10(y7_nz), lums(7), erg_lums(7), env(8,3)
         integer :: ix(y7_nz), iz(y7_nz)
         logical :: conv(y7_nz)
      end type yrec7_fixture
      type :: model2_fixture
         real(dp) :: col(6,m2_nz), omega(m2_nz), xa(15,m2_nz), lums(7), env(8,3)
         logical :: conv(m2_nz), flags(12)
      end type model2_fixture

      type(star_info), save :: s, s2
      integer :: failures
! the file-size limit the program started with (limit_file_size)
      type(file_size_limits) :: starting_limits

      failures = 0
      call start_capture()

      call paths_check()
      call directories_check()
      call pulse_format_check()
      call fields_check()
      call columns_check()
      call groups_check()
      call mod_round_trip_check()
      call mod_files_check()
      call format_tag_check()
      call yrec7_check()
      call model2_check()
      call tables_check()
      call writes_check()
      call fgong_check()

      call end_capture()
      if (failures > 0) then
         write(*,'(i0,a)') failures, " check(s) failed"
         stop 1
      end if

contains

! ---------------------------------------------------------------
! one check: its label, ok or FAIL
subroutine report(label, passed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: passed
      if (passed) then
         write(*,'(a,t72,a)') label, "ok"
      else
         write(*,'(a,t72,a)') label, "FAIL"
         failures = failures + 1
      end if
end subroutine report

! ---------------------------------------------------------------
! io's messages into fresh files terminal.txt and run_log.txt
subroutine start_capture()
      open(newunit=terminal_unit, file=terminal_file, status='replace', action='write')
      open(newunit=run_log_unit, file=run_log_file, status='replace', action='write')
end subroutine start_capture

subroutine end_capture()
      close(terminal_unit)
      close(run_log_unit)
end subroutine end_capture

! Does the file hold a line containing text? (The capture is closed
! first and restarted after, so each check sees only its own messages.)
logical function captured(file, text)
      character(len=*), intent(in) :: file, text
      character(len=1024) :: line
      integer :: u, ios
      call end_capture()
      captured = .false.
      open(newunit=u, file=file, status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) exit
            if (index(line, text) > 0) captured = .true.
         end do
         close(u)
      end if
      call start_capture()
end function captured

! ---------------------------------------------------------------
! a text file of the given lines, replacing any file of that name
subroutine write_text_file(file, lines)
      character(len=*), intent(in) :: file
      character(len=*), intent(in) :: lines(:)
      integer :: u, i
      open(newunit=u, file=file, status='replace', action='write')
      do i = 1, size(lines)
         write(u,'(a)') trim(lines(i))
      end do
      close(u)
end subroutine write_text_file

! the first n lines of a file, copied to another (a truncated file)
subroutine copy_head(source, dest, n)
      character(len=*), intent(in) :: source, dest
      integer, intent(in) :: n
      character(len=1024) :: line
      integer :: u, v, i, ios
      open(newunit=u, file=source, status='old', action='read')
      open(newunit=v, file=dest, status='replace', action='write')
      do i = 1, n
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         write(v,'(a)') trim(line)
      end do
      close(u)
      close(v)
end subroutine copy_head

! the number of lines of a file that start with text
integer function count_lines_starting(file, text)
      character(len=*), intent(in) :: file, text
      character(len=1024) :: line
      integer :: u, ios
      count_lines_starting = 0
      open(newunit=u, file=file, status='old', action='read')
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (index(line, text) == 1) count_lines_starting = count_lines_starting + 1
      end do
      close(u)
end function count_lines_starting

! v as the text of format fmt reads it back: the value a reader of
! that field gets (an internal write and read, not the reader under test)
real(dp) function through(fmt, v)
      character(len=*), intent(in) :: fmt
      real(dp), intent(in) :: v
      character(len=64) :: text
      write(text, fmt) v
      read(text, fmt) through
end function through

! the same bits
elemental logical function same_bits(a, b)
      real(dp), intent(in) :: a, b
      same_bits = transfer(a, 0_int64) == transfer(b, 0_int64)
end function same_bits

! ---------------------------------------------------------------
subroutine set_env(name, value)
      character(len=*), intent(in) :: name, value
      integer(c_int) :: rc
      rc = c_setenv(name // c_null_char, value // c_null_char, 1_c_int)
end subroutine set_env

subroutine unset_env(name)
      character(len=*), intent(in) :: name
      integer(c_int) :: rc
      rc = c_unsetenv(name // c_null_char)
end subroutine unset_env

! expand_value of value, compared with the expected result
logical function expands_to(value, expected)
      character(len=*), intent(in) :: value, expected
      character(len=256) :: v
      v = value
      call expand_value(v)
      expands_to = v == expected
end function expands_to

subroutine paths_check()
      character(len=2) :: short_value
      logical :: ok

      call set_env('YREC_INPUT', '/data/in')
      call unset_env('YREC_START')
      call set_env('YREC_OUTPUT', '')
      call report("paths: {YREC_INPUT} set becomes its value", &
           expands_to('{YREC_INPUT}/eos/x.dat', '/data/in/eos/x.dat') .and. &
           expands_to('{YREC_INPUT}', '/data/in'))
      call report("paths: {YREC_START} unset becomes ../../startmodels", &
           expands_to('{YREC_START}/zahb/m.last', '../../startmodels/zahb/m.last'))
      call report("paths: {YREC_OUTPUT} empty becomes output", &
           expands_to('{YREC_OUTPUT}/run.log', 'output/run.log'))
      call unset_env('YREC_INPUT')
      call report("paths: {YREC_INPUT} unset becomes ../../input", &
           expands_to('{YREC_INPUT}/kap', '../../input/kap'))
      short_value = 'ab'
      call expand_value(short_value)
      ok = short_value == 'ab'
      call report("paths: no leading placeholder, the value unchanged", ok .and. &
           expands_to('plain/{YREC_INPUT}/x', 'plain/{YREC_INPUT}/x') .and. &
           expands_to('rel/path.dat', 'rel/path.dat') .and. &
           expands_to('{YREC_INPU}/x', '{YREC_INPU}/x'))
      call set_env('YREC_INPUT', '{YREC_OUTPUT}/in')
      call unset_env('YREC_OUTPUT')
      call report("paths: an expansion to a later placeholder is expanded again", &
           expands_to('{YREC_INPUT}/x', 'output/in/x'))
      call unset_env('YREC_INPUT')
end subroutine paths_check

! ---------------------------------------------------------------
subroutine directories_check()
      call report("directories: the path up to its last '/'", &
           directory_of('output/samplecases/run.log') == 'output/samplecases/' .and. &
           directory_of('/abs/path/x.log') == '/abs/path/' .and. &
           directory_of('dir/') == 'dir/')
      call report("directories: blank for a bare file name", &
           directory_of('run.log') == ' ' .and. directory_of('') == '')
end subroutine directories_check

! ---------------------------------------------------------------
subroutine fields_check()
      character(len=6) :: i6
      character(len=3) :: i3

      write(i6,'(i6)') 400
      write(i3,'(i3)') -12
      call report("fields: a number that fits is an I<width> field", &
           integer_field(400, 6) == i6 .and. integer_field(-12, 3) == i3 .and. &
           integer_field(999999, 6) == '999999')
      call report("fields: a wider one in full after a blank, no asterisks", &
           integer_field(1000000, 6) == ' 1000000' .and. integer_field(-123, 3) == ' -123')
      call report("fields: the interior GYRE file's model number, 5 digits or all", &
           interior_gyre_path('out/', 'profile', 42) == 'out/profile00042.data.GYRE' .and. &
           interior_gyre_path('', 'p', 99999) == 'p99999.data.GYRE' .and. &
           interior_gyre_path('out/', 'profile', 123456) == 'out/profile123456.data.GYRE' .and. &
           interior_gyre_path('d/', 'x', -3) == 'd/x-3.data.GYRE')
end subroutine fields_check

! ---------------------------------------------------------------
subroutine pulse_format_check()
      call report("pulse format: GYRE, FGONG, GSM in any case", &
           pulse_format_code('GYRE') == pulse_gyre .and. &
           pulse_format_code('gyre') == pulse_gyre .and. &
           pulse_format_code('Gyre') == pulse_gyre .and. &
           pulse_format_code('FGONG') == pulse_fgong .and. &
           pulse_format_code('Fgong') == pulse_fgong .and. &
           pulse_format_code('GSM') == pulse_gsm .and. &
           pulse_format_code('gsm') == pulse_gsm)
      call report("pulse format: leading and trailing blanks ignored", &
           pulse_format_code('  fgong  ') == pulse_fgong)
      call report("pulse format: any other name unknown", &
           pulse_format_code('GYR') == unknown_pulse_format .and. &
           pulse_format_code('GYRE2') == unknown_pulse_format .and. &
           pulse_format_code('') == unknown_pulse_format .and. &
           pulse_format_code('MESA') == unknown_pulse_format)
end subroutine pulse_format_check

! ---------------------------------------------------------------
subroutine columns_check()
      type(output_column) :: table(4)
      character(len=column_name_len) :: defaults(2)
      character(len=16) :: many(max_cols + 1)
      integer :: sel(max_cols), nsel, ierr, i

      table(1) = output_column('alpha', real_column)
      table(2) = output_column('beta', integer_column)
      table(3) = output_column('gamma', real_column)
      table(4) = output_column('delta_x', real_column)
      defaults = [character(len=column_name_len) :: 'gamma', 'alpha']

      call parse_columns(' ', table, defaults, sel, nsel, 'test', ierr)
      call report("columns: a blank file name selects the defaults, in order", &
           ierr == 0 .and. nsel == 2 .and. all(sel(1:2) == [3, 1]))

      call write_text_file('columns_good.list', [character(len=40) :: &
           '! the test columns', '', '  beta    ! the id column', 'delta_x', &
           '   ', 'alpha'])
      call parse_columns('columns_good.list', table, defaults, sel, nsel, 'test', ierr)
      call report("columns: a file selects its names in order, past comments", &
           ierr == 0 .and. nsel == 3 .and. all(sel(1:3) == [2, 4, 1]))

      call write_text_file('columns_unknown.list', [character(len=40) :: 'beta', 'epsilon'])
      call parse_columns('columns_unknown.list', table, defaults, sel, nsel, 'test', ierr)
      call report("columns: an unknown name is an error", ierr == 1 .and. &
           captured(terminal_file, 'unknown test column: epsilon'))
      call parse_columns('columns_unknown.list', table, defaults, sel, nsel, 'test', ierr)
      call report("columns: the run log lists the valid names", ierr == 1 .and. &
           captured(run_log_file, 'delta_x'))

      call write_text_file('columns_empty.list', [character(len=40) :: '! nothing', ''])
      call parse_columns('columns_empty.list', table, defaults, sel, nsel, 'test', ierr)
      call report("columns: a file that selects nothing is an error", ierr == 1 .and. &
           captured(terminal_file, 'selects no columns'))

      call parse_columns('columns_missing.list', table, defaults, sel, nsel, 'test', ierr)
      call report("columns: a missing file is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot open the test_columns_file columns_missing.list'))

      do i = 1, size(many)
         many(i) = 'alpha'
      end do
      call write_text_file('columns_many.list', many)
      call parse_columns('columns_many.list', table, defaults, sel, nsel, 'test', ierr)
      call report("columns: more than max_cols names is an error", ierr == 1 .and. &
           captured(terminal_file, 'more than 128 test columns'))
end subroutine columns_check

! ---------------------------------------------------------------
subroutine groups_check()
      call write_text_file('groups.nml', [character(len=40) :: &
           '! &physics in a comment', '&control  ', '  x = 1', '/', &
           '   &PHYSICS', '/', '&controls_extra', '/', '&empty/', &
           '&tabbed' // achar(9) // 'y = 2 /'])
      call report("groups: found in any case, after blanks", &
           nml_file_has_group('groups.nml', 'control') .and. &
           nml_file_has_group('groups.nml', 'physics') .and. &
           nml_file_has_group('groups.nml', 'Control'))
      call report("groups: a blank, a tab or '/' may end the name", &
           nml_file_has_group('groups.nml', 'empty') .and. &
           nml_file_has_group('groups.nml', 'tabbed') .and. &
           nml_file_has_group('groups.nml', 'controls_extra'))
      call report("groups: not a prefix of a longer name, nor a comment", &
           .not. nml_file_has_group('groups.nml', 'controls') .and. &
           .not. nml_file_has_group('groups.nml', 'contro') .and. &
           .not. nml_file_has_group('groups.nml', 'star_job'))
      call write_text_file('groups_comment.nml', [character(len=40) :: &
           '! &physics in a comment', '&control', '/'])
      call report("groups: none in a comment, none in a missing file", &
           .not. nml_file_has_group('groups_comment.nml', 'physics') .and. &
           .not. nml_file_has_group('groups_missing.nml', 'control'))
end subroutine groups_check

! ---------------------------------------------------------------
! A synthetic seven-shell model on star s, its values spread over the
! range a double holds: the .mod round trip's test model.
subroutine fill_mod_model(star)
      type(star_info), intent(inout) :: star
      integer :: k, i
      real(dp) :: f

      star%nz = 7
      star%model_number = 4217
      star%star_mass = 1.0d0/3.0d0
      star%log_total_mass = 33.29863242798d0
      star%age_gyr = 4.5677777777777779d0
      star%timestep_yr = 1.234567890123456789d7
      star%log_L_lsun = -1.0d-310
      star%log_Teff = 3.7617d0
      star%mixing_length_alpha = 1.8518851851851851d0
      star%physics%use_extended_composition = .true.
      star%physics%rotation_active = .false.
      do k = 1, star%nz
         f = sin(dble(k))
         star%log_mass(k) = 30.0d0 + f
         star%logR(k) = 9.0d0 + f/7.0d0
         star%luminosity_lsun(k) = f*10.0d0**(40*k - 160)
         star%logP(k) = 17.0d0 - f*acos(-1.0d0)
         star%logT(k) = 7.0d0 - exp(-dble(k))
         star%logRho(k) = 2.0d0 - dble(k)/3.0d0
         star%convective_flag(k) = mod(k, 3) == 0
         star%omega(k) = 2.5d-6*f**3
         do i = 1, size(star%xa, 1)
            star%xa(i,k) = 0.1d0*abs(cos(dble(i*k)))*10.0d0**(-i)
         end do
      end do
      star%omega(1) = 0.0d0
      star%xa(i_be9,2) = 1.0d-300
end subroutine fill_mod_model

! the .mod's slots, by slot (model_files' mod_slots)
subroutine mod_round_trip_check()
      integer, parameter :: slots(15) = [i_h1, i_he4, i_metals, i_he3, i_c12, &
           i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, i_li6, i_li7, i_be9]
      type(model_file_header) :: h
      integer :: u, ierr, n, izone
      real(dp) :: row(22)
      logical :: ok, conv
      character(len=1024) :: line

      call fill_mod_model(s)
      open(newunit=u, file='round_trip.mod', status='replace')
      ierr = 0
      call write_mod_model(s, u, replace=.true., ierr=ierr)
      s2%nz = 0
      s2%xa = -1.0d0
      call read_mod_model(s2, u, h, ierr)
      close(u)
      n = s%nz
      call report(".mod: read back without error", ierr == 0)
      ok = s2%nz == n .and. s2%model_number == s%model_number .and. &
           same_bits(s2%star_mass, s%star_mass) .and. &
           same_bits(s2%log_total_mass, s%log_total_mass) .and. &
           same_bits(s2%age_gyr, s%age_gyr) .and. &
           same_bits(h%timestep_yr, s%timestep_yr) .and. &
           same_bits(s2%log_L_lsun, s%log_L_lsun) .and. &
           same_bits(s2%log_Teff, s%log_Teff) .and. &
           same_bits(h%mixing_length, s%mixing_length_alpha) .and. &
           h%use_extended_composition .and. .not. h%rotation_active
      call report(".mod: every global value bit for bit", ok)
      ok = all(same_bits(s2%log_mass(1:n), s%log_mass(1:n))) .and. &
           all(same_bits(s2%logR(1:n), s%logR(1:n))) .and. &
           all(same_bits(s2%luminosity_lsun(1:n), s%luminosity_lsun(1:n))) .and. &
           all(same_bits(s2%logP(1:n), s%logP(1:n))) .and. &
           all(same_bits(s2%logT(1:n), s%logT(1:n))) .and. &
           all(same_bits(s2%logRho(1:n), s%logRho(1:n))) .and. &
           all(s2%convective_flag(1:n) .eqv. s%convective_flag(1:n)) .and. &
           all(same_bits(s2%omega(1:n), s%omega(1:n)))
      call report(".mod: every shell's structure bit for bit", ok)
      ok = all(same_bits(s2%xa(slots,1:n), s%xa(slots,1:n)))
      call report(".mod: every shell's composition bit for bit, in its slot", ok)

! the file itself: its columns line, and the first row in that order
      open(newunit=u, file='round_trip.mod', status='old', action='read')
      do
         read(u,'(a)') line
         if (index(line, 'columns:') == 1) exit
      end do
      read(u,*) izone, row(1:6), conv, row(7:22)
      close(u)
      ok = trim(line) == 'columns: zone log_mass logR luminosity logP logT logRho conv' &
           // ' omega h1 he4 z he3 c12 c13 n14 n15 o16 o17 o18 h2 li6 li7 be9' .and. &
           izone == 1 .and. all(same_bits(row(1:6), [s%log_mass(1), s%logR(1), &
           s%luminosity_lsun(1), s%logP(1), s%logT(1), s%logRho(1)])) .and. &
           (conv .eqv. s%convective_flag(1)) .and. same_bits(row(7), s%omega(1)) .and. &
           all(same_bits(row(8:22), s%xa(slots,1)))
      call report(".mod: the columns line, and a row in its order", ok)
end subroutine mod_round_trip_check

! the .mod's files: unknown globals, replace and append, errors
subroutine mod_files_check()
      character(len=1024) :: lines(40)
      type(model_file_header) :: h
      integer :: u, ierr, n, i

      call fill_mod_model(s)
! one model replaced, then two appended
      ierr = 0
      open(newunit=u, file='replace.mod', status='replace')
      call write_mod_model(s, u, replace=.true., ierr=ierr)
      call write_mod_model(s, u, replace=.true., ierr=ierr)
      close(u)
      open(newunit=u, file='append.mod', status='replace')
      call write_mod_model(s, u, replace=.false., ierr=ierr)
      call write_mod_model(s, u, replace=.false., ierr=ierr)
      close(u)
      call report(".mod: replace keeps one model, append adds one", ierr == 0 .and. &
           count_lines_starting('replace.mod', 'YMOD 1') == 1 .and. &
           count_lines_starting('append.mod', 'YMOD 1') == 2)

! an unknown global, inserted before the columns line
      n = line_count('round_trip.mod')
      lines(1:n) = read_back('round_trip.mod', n)
      call write_text_file('unknown_global.mod', [lines(1:2), &
           [character(len=1024) :: 'future_global        1.5'], lines(3:n)])
      open(newunit=u, file='unknown_global.mod', status='old')
      call read_mod_model(s2, u, h, ierr)
      close(u)
      call report(".mod: an unknown global is noted and skipped", ierr == 0 .and. &
           s2%nz == s%nz .and. same_bits(s2%logT(s%nz), s%logT(s%nz)) .and. &
           captured(run_log_file, 'ignoring the unknown global future_global'))

! errors: another columns line, a cut, too many zones
      do i = 1, n
         if (index(lines(i), 'columns:') == 1) lines(i) = 'columns: zone log_mass logR'
      end do
      call write_text_file('other_columns.mod', lines(1:n))
      call read_mod_error('other_columns.mod', 'has the columns', &
           ".mod: another columns line is an error")
      call copy_head('round_trip.mod', 'cut.mod', 17)
      call read_mod_error('cut.mod', 'cannot read the shells of the starting model', &
           ".mod: a file cut inside its shells is an error")
      lines(1:n) = read_back('round_trip.mod', n)
      write(lines(3), '(a,1x,i12)') 'num_zones           ', json + 1
      call write_text_file('many_zones.mod', lines(1:n))
      call read_mod_error('many_zones.mod', 'zones, not 1 to', &
           ".mod: a zone count above json is an error")
end subroutine mod_files_check

! the number of lines of a file
integer function line_count(file)
      character(len=*), intent(in) :: file
      integer :: u, ios
      line_count = 0
      open(newunit=u, file=file, status='old', action='read')
      do
         read(u,*,iostat=ios)
         if (ios /= 0) exit
         line_count = line_count + 1
      end do
      close(u)
end function line_count

! the first n lines of a file
function read_back(file, n) result(lines)
      character(len=*), intent(in) :: file
      integer, intent(in) :: n
      character(len=1024) :: lines(n)
      integer :: u, i
      open(newunit=u, file=file, status='old', action='read')
      do i = 1, n
         read(u,'(a)') lines(i)
      end do
      close(u)
end function read_back

! read_mod_model on a damaged file: ierr = 1 and the message
subroutine read_mod_error(file, message, label)
      character(len=*), intent(in) :: file, message, label
      type(model_file_header) :: h
      integer :: u, ierr
      open(newunit=u, file=file, status='old')
      call read_mod_model(s2, u, h, ierr)
      close(u)
      call report(label, ierr == 1 .and. captured(terminal_file, message))
end subroutine read_mod_error

! ---------------------------------------------------------------
subroutine format_tag_check()
      character(len=4) :: tag
      integer :: u, ierr

      open(newunit=u, file='round_trip.mod', status='old')
      read(u,*)
      call read_model_format_tag(u, tag, ierr)
      close(u)
      call report("format tag: the first four characters, from the start", &
           ierr == 0 .and. tag == 'YMOD')
      open(newunit=u, file='empty.model', status='replace')
      call read_model_format_tag(u, tag, ierr)
      close(u)
      call report("format tag: an empty file is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot read the format tag of the starting model'))
end subroutine format_tag_check

! ---------------------------------------------------------------
! A YREC7 file of six shells: a convective core (shells 1-2), the
! radiative shells between the zones' edges (jcore = 2, jenv = 4), a
! convective envelope (4-6), rotating, with the extended composition:
! three shells from jcore to jenv, so the light elements take a pair
! line and a single one. Its values are the fixture's (fill_yrec7).
subroutine fill_yrec7(f)
      type(yrec7_fixture), intent(out) :: f
      integer :: k, i

      do k = 1, y7_nz
         f%log_mass(k) = through('(F13.10)', 32.0d0 + 0.15d0*k + 1.0d-11/3.0d0)
         f%log_r(k) = through('(F10.7)', 8.0d0 + 0.4d0*k)
         f%lum(k) = through('(1PE14.7)', 0.1d0*k*k/3.0d0)
         f%log_p(k) = through('(F10.7)', 17.0d0 - 1.7d0*k)
         f%log_t(k) = through('(F10.7)', 7.2d0 - 0.4d0*k)
         f%log_rho(k) = through('(F10.7)', 2.0d0 - 1.1d0*k)
         f%conv(k) = k <= y7_jcore .or. k >= y7_jenv
         f%ix(k) = 700000 - 11111*k
         f%iz(k) = 18804 + k
         do i = 4, 11
            f%cno(i,k) = through('(1PE9.3)', 1.0d-3*i/(k + 0.7d0))
         end do
         do i = 12, 15
            f%light(i,k) = through('(1PE9.3)', 1.0d-9*i*k/7.0d0)
         end do
         f%omega_log10(k) = through('(F10.7)', 5.3d0 + 0.1d0*k/3.0d0)
      end do
      f%omega_log10(1) = 59.0d0
      do i = 1, 7
         f%lums(i) = through('(E10.3)', 0.1d0*i/3.0d0)
         f%erg_lums(i) = through('(E10.3)', f%lums(i)*1.0d33)
      end do
      do k = 1, 3
         f%env(1,k) = through('(F7.5)', 3.76d0 + 0.01d0*k/3.0d0)
         do i = 2, 5
            f%env(i,k) = through('(F8.5)', 0.1d0*i + 0.001d0*k/3.0d0)
         end do
         do i = 6, 8
            f%env(i,k) = through('(E12.5)', 1.0d0*i*k/7.0d0)
         end do
      end do
end subroutine fill_yrec7

! The fixture as a YREC7 file, in the layout of model_files' header:
! nz_header is the header's zone count, envelope_edge the flags line's
! envelope convection zone bottom; erg writes the luminosities in erg/s;
! without with_omega the file ends before the rotation rates; reflected
! gives the triangle's second record a negative record number.
subroutine write_yrec7(file, f, nz_header, envelope_edge, with_omega, erg, reflected)
      character(len=*), intent(in) :: file
      type(yrec7_fixture), intent(in) :: f
      integer, intent(in) :: nz_header, envelope_edge
      logical, intent(in) :: with_omega, erg, reflected
      character(len=*), parameter :: f_header = '(A4,1X,2I5,3F8.5,1X,F14.10,F13.9,F13.0)', &
           f_flags = '(2(L1,1X),5X,2I4,F7.3)', f_lum = '(10X,7E10.3)', &
           f_env = '(3X,I2,F7.5,4F8.5,3E12.5)', &
           f_shell = '(0PF13.10,F10.7,1PE14.7,0PF10.7,2F10.7,L1,2I6)', &
           f_cno = '(8(1PE9.3,1X))', f_light = '(4(1PE9.3,1X))', f_omega = '(0P8F10.7)'
      integer :: v, j, k

      open(newunit=v, file=file, status='replace', action='write')
      write(v, f_header) 'NMOD', 12, nz_header, 0.95d0, 3.7617d0, 0.0123d0, 33.2d0, &
           4.567d0, 1.0d6
      write(v, f_flags) .true., .true., y7_jcore, envelope_edge, 1.875d0
      if (erg) then
         write(v, f_lum) f%erg_lums
      else
         write(v, f_lum) f%lums
      end if
      do k = 1, 3
         j = k
         if (k == 2 .and. reflected) j = -k
         write(v, f_env) j, f%env(:,k)
      end do
      do k = 1, y7_nz
         write(v, f_shell) f%log_mass(k), f%log_r(k), f%lum(k), f%log_p(k), f%log_t(k), &
              f%log_rho(k), f%conv(k), f%ix(k), f%iz(k)
      end do
      do k = y7_jcore, y7_jenv
         write(v, f_cno) f%cno(:,k)
      end do
! H2-Be9, jcore to jenv: two shells per line, the odd one on its own
      do k = y7_jcore, y7_jenv - 1, 2
         write(v, f_cno) f%light(:,k), f%light(:,k+1)
      end do
      if (mod(y7_jenv - 1 - y7_jcore, 2) /= 0) write(v, f_light) f%light(:,y7_jenv)
      if (with_omega) write(v, f_omega) f%omega_log10
      close(v)
end subroutine write_yrec7

! read_yrec7 on a file; ierr out
subroutine read_yrec7_file(file, header, ierr)
      character(len=*), intent(in) :: file
      type(model_file_header), intent(out) :: header
      integer, intent(out) :: ierr
      integer :: u
      open(newunit=u, file=file, status='old')
      call read_yrec7(s, u, header, ierr)
      close(u)
end subroutine read_yrec7_file

subroutine yrec7_check()
      type(yrec7_fixture) :: f
      type(model_file_header) :: h
      integer :: k, i, ierr
      logical :: ok
      real(dp) :: saved_log_mass

      call fill_yrec7(f)
      s%solar_luminosity_cgs = 3.828d33
      call write_yrec7('good.yrec7', f, y7_nz, y7_jenv, &
           with_omega=.true., erg=.false., reflected=.true.)
      call read_yrec7_file('good.yrec7', h, ierr)
      call report("YREC7: read without error", ierr == 0)
      ok = s%model_number == 12 .and. s%nz == y7_nz .and. &
           s%star_mass == through('(F8.5)', 0.95d0) .and. &
           s%log_Teff == through('(F8.5)', 3.7617d0) .and. &
           s%log_L_lsun == through('(F8.5)', 0.0123d0) .and. &
           s%log_total_mass == through('(F14.10)', 33.2d0) .and. &
           s%age_gyr == through('(F13.9)', 4.567d0) .and. h%timestep_yr == 1.0d6 .and. &
           h%rotation_active .and. h%use_extended_composition .and. h%core_cz_top_index == y7_jcore .and. h%envelope_cz_bottom_index == y7_jenv &
           .and. h%mixing_length == through('(F7.3)', 1.875d0)
      call report("YREC7: the header and the flags", ok)
      ok = h%eos_code == '    ? ' .and. h%atm_code == '  ? ' .and. h%alok_code == '  ? ' &
           .and. h%hik_code == '  ? ' .and. .not. (h%use_diffusion_y .or. h%use_diffusion_z &
           .or. h%disk_locking_active .or. h%instability_transport_active .or. &
           h%use_wind_torque .or. h%core_overshoot_active .or. h%envelope_overshoot_active &
           .or. h%intermediate_overshoot_active .or. h%use_pure_z_table .or. &
           h%use_semiconvection) .and. &
           h%disk_omega_rad_s == 7.2722d-6 .and. h%disk_locking_age_gyr == 0.0d0 .and. &
           h%wind_saturation_omega == 3.0d-4
      call report("YREC7: MODEL2's header values at their defaults", ok)
      ok = all(s%log_mass(1:y7_nz) == f%log_mass) .and. all(s%logR(1:y7_nz) == f%log_r) &
           .and. all(s%luminosity_lsun(1:y7_nz) == f%lum) .and. &
           all(s%logP(1:y7_nz) == f%log_p) .and. all(s%logT(1:y7_nz) == f%log_t) .and. &
           all(s%logRho(1:y7_nz) == f%log_rho) .and. &
           all(s%convective_flag(1:y7_nz) .eqv. f%conv)
      call report("YREC7: the shells' structure", ok)
      ok = all(s%xa(i_h1,1:y7_nz) == 1.0d-6*dble(f%ix)) .and. &
           all(s%xa(i_metals,1:y7_nz) == 1.0d-6*dble(f%iz))
      call report("YREC7: X and Z in units of 1e-6", ok)
      ok = .true.
      do k = 1, y7_nz
         do i = 4, 11
            ok = ok .and. s%xa(i,k) == f%cno(i, min(max(k, y7_jcore), y7_jenv))
         end do
      end do
      call report("YREC7: He3-O18 per shell, copied through the convection zones", ok)
      ok = all(s%xa(i_he4,1:y7_nz) == 1.0d0 - s%xa(i_h1,1:y7_nz) - s%xa(i_metals,1:y7_nz) &
           - s%xa(i_he3,1:y7_nz))
      call report("YREC7: He4 = 1 - X - Z - He3", ok)
      ok = .true.
      do k = 1, y7_nz
         do i = 12, 15
            ok = ok .and. s%xa(i,k) == f%light(i, min(max(k, y7_jcore), y7_jenv))
         end do
      end do
      call report("YREC7: H2-Be9 two shells per line, the odd one alone", ok)
      ok = s%omega(1) == 0.0d0 .and. &
           all(abs(s%omega(2:y7_nz)/10.0d0**(-f%omega_log10(2:y7_nz)) - 1.0d0) < 1.0d-14)
      call report("YREC7: omega from -log10(omega), zero at the sentinel", ok)
      ok = all(s%triangle_log_teff == f%env(1,:)) .and. &
           all(s%triangle_log_l_lsun == f%env(2,:)) .and. &
           all(s%fit_point_logP == f%env(3,:)) .and. &
           all(s%fit_point_logT == f%env(4,:)) .and. &
           all(s%fit_point_logR == f%env(5,:)) .and. &
           all(s%envelope_fit_coeffs == reshape(f%env(6:8,:), [9])) .and. &
           h%triangle_orientation == -1.0d0
      call report("YREC7: the envelope records, the triangle reflected", ok)
      call report("YREC7: luminosities in Lsun kept", &
           all(s%luminosity_breakdown_lsun(1:7) == f%lums))

      call write_yrec7('erg.yrec7', f, y7_nz, y7_jenv, with_omega=.true., erg=.true., &
           reflected=.false.)
      call read_yrec7_file('erg.yrec7', h, ierr)
      call report("YREC7: luminosities in erg/s converted to Lsun", ierr == 0 .and. &
           all(s%luminosity_breakdown_lsun(1:7) == f%erg_lums/s%solar_luminosity_cgs))
      call report("YREC7: no negative record number, the triangle as it is", &
           ierr == 0 .and. h%triangle_orientation == 1.0d0)

      call write_yrec7('no_omega.yrec7', f, y7_nz, y7_jenv, &
           with_omega=.false., erg=.false., reflected=.true.)
      call read_yrec7_file('no_omega.yrec7', h, ierr)
      call report("YREC7: a rotating file without omega reads, omega zero", ierr == 0 .and. &
           all(s%omega(1:y7_nz) == 0.0d0) .and. captured(run_log_file, 'has no omega'))

      call copy_head('good.yrec7', 'cut.yrec7', 8)
      call yrec7_error('cut.yrec7', 'cannot read the shells of the starting model', &
           "YREC7: a file cut inside its shells is an error")
      call write_yrec7('zero_zones.yrec7', f, 0, y7_jenv, &
           with_omega=.true., erg=.false., reflected=.true.)
      call yrec7_error('zero_zones.yrec7', 'has 0 zones', &
           "YREC7: a zone count of 0 is an error")
      saved_log_mass = f%log_mass(3)
      f%log_mass(3) = 33.5d0
      call write_yrec7('heavy_shell.yrec7', f, y7_nz, y7_jenv, &
           with_omega=.true., erg=.false., reflected=.true.)
      call yrec7_error('heavy_shell.yrec7', 'has a mass below zero or above', &
           "YREC7: a shell mass above the star's is an error")
      f%log_mass(3) = saved_log_mass
      call write_yrec7('bad_edge.yrec7', f, y7_nz, y7_nz + 3, &
           with_omega=.true., erg=.false., reflected=.true.)
      call yrec7_error('bad_edge.yrec7', 'envelope convection zone bottom at zone 9', &
           "YREC7: an envelope convection zone outside the model is an error")
end subroutine yrec7_check

! read_yrec7 on a damaged file: ierr = 1 and the message
subroutine yrec7_error(file, message, label)
      character(len=*), intent(in) :: file, message, label
      type(model_file_header) :: h
      integer :: ierr
      call read_yrec7_file(file, h, ierr)
      call report(label, ierr == 1 .and. captured(terminal_file, message))
end subroutine yrec7_error

! ---------------------------------------------------------------
! A MODEL2 file of four shells, every column distinct.
subroutine fill_model2(f)
      type(model2_fixture), intent(out) :: f
      integer :: k, i

      do k = 1, m2_nz
         f%col(1,k) = through('(F18.14)', 32.0d0 + k/3.0d0)
         f%col(2,k) = through('(F18.14)', 8.0d0 + k/7.0d0)
         f%col(3,k) = through('(1PE24.16)', k/3.0d0)
         do i = 4, 6
            f%col(i,k) = through('(F18.14)', i + k/9.0d0)
         end do
         f%omega(k) = through('(1PE24.16)', 1.0d-6*k/3.0d0)
         f%conv(k) = mod(k, 2) == 1
         do i = 1, 3
            f%xa(i,k) = through('(F12.9)', 0.1d0*i + k/97.0d0)
         end do
         do i = 4, 15
            f%xa(i,k) = through('(E16.8)', 1.0d-4*i/(k + 0.3d0))
         end do
      end do
      do i = 1, 7
         f%lums(i) = through('(1PE17.9)', 0.1d0*i/7.0d0)
      end do
      do k = 1, 3
         do i = 1, 5
            f%env(i,k) = through('(F16.12)', 0.5d0*i + k/11.0d0)
         end do
         do i = 6, 8
            f%env(i,k) = through('(1PE20.12)', 1.0d0*i*k/13.0d0)
         end do
      end do
! the flags line's logicals, in file order: the pure-Z table, then the
! extended composition, the diffusion of Y and of Z, semiconvection,
! core and envelope overshoot, intermediate overshoot, rotation, instability transport,
! the wind torque, disk locking
      f%flags = [.false., .true., .false., .true., .false., .false., .true., .false., &
           .true., .false., .true., .true.]
end subroutine fill_model2

! The fixture as a MODEL2 file (model_files' header gives the layout);
! nz_header is the header's zone count.
subroutine write_model2(file, f, nz_header)
      character(len=*), intent(in) :: file
      type(model2_fixture), intent(in) :: f
      integer, intent(in) :: nz_header
      character(len=*), parameter :: &
           f_header = '(A5,2I8,5F16.12,1PE18.10,0P2F16.12)', &
           f_flags = '(2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,11(L1,1X),3(1PE18.10))', &
           f_lum = '(10X,1P7E17.9)', f_env = '(4X,I2,5F16.12,1P3E20.12)', &
           f_shell = '(I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1x,L1,3(0PF12.9),12(0PE16.8))'
      integer :: v, k, j

      open(newunit=v, file=file, status='replace', action='write')
      write(v, f_header) 'MOD2 ', 77, nz_header, 1.0d0/3.0d0, 3.7d0, 0.1d0, 33.3d0, &
           4.5d0, 1.0d5/3.0d0, f%col(1,1), f%col(1,m2_nz)
      write(v, f_flags) 1, 3, 1.9d0, 'OPAL06', 'GRAY', 'ALEX', 'OP17', f%flags(1), 'GS98', &
           f%flags(2:12), 1.0d0/3.0d0, 7.2722d-6, 3.0d-4
      write(v, f_lum) f%lums
! the triangle's third record reflected
      do k = 1, 3
         j = k
         if (k == 3) j = -k
         write(v, f_env) j, f%env(:,k)
      end do
      write(v,'(a)') 'heading 1', 'heading 2', 'heading 3'
      do k = 1, m2_nz
         write(v, f_shell) k, f%col(1:6,k), f%omega(k), f%conv(k), f%xa(:,k)
      end do
      close(v)
end subroutine write_model2

! read_model2 on a file; ierr out
subroutine read_model2_file(file, header, ierr)
      character(len=*), intent(in) :: file
      type(model_file_header), intent(out) :: header
      integer, intent(out) :: ierr
      integer :: u
      open(newunit=u, file=file, status='old')
      call read_model2(s, u, header, ierr)
      close(u)
end subroutine read_model2_file

subroutine model2_check()
      type(model2_fixture) :: f
      type(model_file_header) :: h
      integer :: ierr
      logical :: ok

      call fill_model2(f)
      s%solar_luminosity_cgs = 3.828d33
      call write_model2('good.model2', f, m2_nz)
      call read_model2_file('good.model2', h, ierr)
      call report("MODEL2: read without error", ierr == 0)
      ok = s%model_number == 77 .and. s%nz == m2_nz .and. &
           s%star_mass == through('(F16.12)', 1.0d0/3.0d0) .and. &
           s%log_Teff == through('(F16.12)', 3.7d0) .and. &
           s%log_L_lsun == through('(F16.12)', 0.1d0) .and. &
           s%log_total_mass == through('(F16.12)', 33.3d0) .and. &
           s%age_gyr == through('(F16.12)', 4.5d0) .and. &
           h%timestep_yr == through('(1PE18.10)', 1.0d5/3.0d0) .and. &
           h%core_cz_top_index == 1 .and. h%envelope_cz_bottom_index == 3 .and. &
           h%mixing_length == through('(F16.10)', 1.9d0) .and. h%triangle_orientation == -1.0d0
      call report("MODEL2: the header", ok)
      ok = h%eos_code == 'OPAL06' .and. h%atm_code == 'GRAY' .and. &
           h%alok_code == 'ALEX' .and. h%hik_code == 'OP17' .and. &
           (h%use_pure_z_table .eqv. f%flags(1)) .and. &
           (h%use_extended_composition .eqv. f%flags(2)) .and. &
           (h%use_diffusion_y .eqv. f%flags(3)) .and. (h%use_diffusion_z .eqv. f%flags(4)) .and. &
           (h%use_semiconvection .eqv. f%flags(5)) .and. &
           (h%core_overshoot_active .eqv. f%flags(6)) .and. &
           (h%envelope_overshoot_active .eqv. f%flags(7)) .and. &
           (h%intermediate_overshoot_active .eqv. f%flags(8)) .and. &
           (h%rotation_active .eqv. f%flags(9)) .and. &
           (h%instability_transport_active .eqv. f%flags(10)) .and. &
           (h%use_wind_torque .eqv. f%flags(11)) .and. (h%disk_locking_active .eqv. f%flags(12)) &
           .and. h%disk_locking_age_gyr == through('(1PE18.10)', 1.0d0/3.0d0) .and. &
           h%disk_omega_rad_s == through('(1PE18.10)', 7.2722d-6) .and. &
           h%wind_saturation_omega == through('(1PE18.10)', 3.0d-4)
      call report("MODEL2: the codes, the flags and the disk constants", ok)
      ok = all(s%luminosity_breakdown_lsun(1:7) == f%lums) .and. &
           all(s%triangle_log_teff == f%env(1,:)) .and. &
           all(s%triangle_log_l_lsun == f%env(2,:)) .and. &
           all(s%fit_point_logP == f%env(3,:)) .and. &
           all(s%fit_point_logT == f%env(4,:)) .and. &
           all(s%fit_point_logR == f%env(5,:)) .and. &
           all(s%envelope_fit_coeffs == reshape(f%env(6:8,:), [9]))
      call report("MODEL2: the luminosities and the envelope records", ok)
      ok = all(s%log_mass(1:m2_nz) == f%col(1,:)) .and. all(s%logR(1:m2_nz) == f%col(2,:)) &
           .and. all(s%luminosity_lsun(1:m2_nz) == f%col(3,:)) .and. &
           all(s%logP(1:m2_nz) == f%col(4,:)) .and. all(s%logT(1:m2_nz) == f%col(5,:)) &
           .and. all(s%logRho(1:m2_nz) == f%col(6,:)) .and. &
           all(s%omega(1:m2_nz) == f%omega) .and. &
           all(s%convective_flag(1:m2_nz) .eqv. f%conv) .and. all(s%xa(1:15,1:m2_nz) == f%xa)
      call report("MODEL2: every shell's columns", ok)

      call copy_head('good.model2', 'cut.model2', 10)
      call model2_error('cut.model2', 'cannot read the shells of the starting model', &
           "MODEL2: a file cut inside its shells is an error")
      call write_model2('many_zones.model2', f, json + 1)
      call model2_error('many_zones.model2', 'zones, not 1 to', &
           "MODEL2: a zone count above json is an error")
end subroutine model2_check

! read_model2 on a damaged file: ierr = 1 and the message
subroutine model2_error(file, message, label)
      character(len=*), intent(in) :: file, message, label
      type(model_file_header) :: h
      integer :: ierr
      call read_model2_file(file, h, ierr)
      call report(label, ierr == 1 .and. captured(terminal_file, message))
end subroutine model2_error

! ---------------------------------------------------------------
subroutine tables_check()
      call table_check("tables: history", history_columns(), hist_default_names)
      call table_check("tables: profile", profile_columns(), prof_default_names)
end subroutine tables_check

! one column table: its names and kinds, and its default selection
subroutine table_check(label, table, defaults)
      character(len=*), intent(in) :: label
      type(output_column), intent(in) :: table(:)
      character(len=column_name_len), intent(in) :: defaults(:)
      integer :: sel(max_cols), nsel, ierr, i, j
      logical :: ok

      ok = .true.
      do i = 1, size(table)
         ok = ok .and. len_trim(table(i)%name) > 0 .and. &
              len_trim(table(i)%name) < column_name_len
         do j = 1, i - 1
            ok = ok .and. table(i)%name /= table(j)%name
         end do
      end do
      call report(label // " names nonblank, short enough and unique", ok)
      call report(label // " kinds real_column or integer_column", &
           all(table%kind == real_column .or. table%kind == integer_column))
      call parse_columns(' ', table, defaults, sel, nsel, 'test', ierr)
      ok = ierr == 0 .and. nsel == size(defaults)
      if (ok) ok = all(table(sel(1:nsel))%name == defaults)
      call report(label // " default selection resolves, in order", ok)
end subroutine table_check

! ---------------------------------------------------------------
! A full disk, as a writer sees it: files no longer than max_bytes
! (RLIMIT_FSIZE, with the signal a write past it raises ignored, so the
! write fails with EFBIG, as one to a full disk fails with ENOSPC);
! max_bytes < 0 restores the limit the program started with.
subroutine limit_file_size(max_bytes)
      integer, intent(in) :: max_bytes
      type(file_size_limits) :: limits
      type(c_funptr) :: previous_handler
      integer(c_int) :: rc

      if (max_bytes >= 0) then
         previous_handler = c_signal(sigxfsz, transfer(sig_ign, c_null_funptr))
         rc = c_getrlimit(rlimit_fsize, starting_limits)
         limits = starting_limits
         limits%soft = max_bytes
         rc = c_setrlimit(rlimit_fsize, limits)
      else
         rc = c_setrlimit(rlimit_fsize, starting_limits)
      end if
end subroutine limit_file_size

! ---------------------------------------------------------------
! Another run that writes the same path, through the C library (the
! file is open on a unit here, which gfortran would not open twice):
! mode 'w' opens it as gfortran's status='replace' does, the same file
! truncated in place; 'a' appends to it. The line text is written and
! the file closed.
subroutine other_run_writes(path, mode, text)
      character(len=*), intent(in) :: path, mode, text
      type(c_ptr) :: stream
      integer(c_int) :: rc
      stream = c_fopen(path // c_null_char, mode // c_null_char)
      if (.not. c_associated(stream)) error stop 'other_run_writes: fopen failed'
      rc = c_fputs(text // c_new_line // c_null_char, stream)
      rc = c_fclose(stream)
end subroutine other_run_writes

! Another run's file put at path in place of the one there: a new file
! of one line renamed over it.
subroutine other_run_replaces(path, text)
      character(len=*), intent(in) :: path, text
      integer(c_int) :: rc
      call write_text_file('other_run.tmp', [text])
      rc = c_rename('other_run.tmp' // c_null_char, path // c_null_char)
      if (rc /= 0) error stop 'other_run_replaces: rename failed'
end subroutine other_run_replaces

! ---------------------------------------------------------------
subroutine writes_check()
      integer, parameter :: n_points = 20, limit_bytes = 2000
      type(history_file) :: history
      type(profile_selection) :: profiles
      real(dp) :: pts(n_pulse_cols, n_points)
      character(len=200) :: deck_lines(100)
      character(len=256) :: msg
      character(len=8) :: lines(3)
      integer :: u, ierr, i, j, ios, row_ierr, test_log_unit

      call fill_mod_model(s)
      do j = 1, n_points
         do i = 1, n_pulse_cols
            pts(i,j) = dble(i) + 1.0d-3*dble(j)
         end do
      end do
      call write_text_file('readonly.txt', [character(len=8) :: 'locked'])
      call write_text_file('not_a_dir', [character(len=8) :: 'a file'])

! the files kept open for the run, on a unit open for reading only
      open(newunit=u, file='readonly.txt', status='old', action='read')
      ierr = 0
      call write_mod_model(s, u, replace=.true., ierr=ierr)
      call report("writes: a failed write of the .mod model is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write the model file'))
      history%path = 'readonly.txt'
      history%unit = u
      history%is_open = .true.
      history%nsel = 2
      history%sel(1:2) = [1, 4]
      ierr = 0
      call write_history_row(s, history, 0, ierr)
      call report("writes: a failed history row is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write readonly.txt'))
      isochrone_unit = u
      s%control%isochrone_output_active = .true.
      s%current_card = 1
      s%card_kind(1) = card_evolve
      ierr = 0
      call output_run_header(s, 1.0d0, ierr)
      call report("writes: a failed write of the isochrone file is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write the isochrone file'))
      close(u)
      isochrone_unit = -1
      s%control%isochrone_output_active = .false.

! a file whose directory cannot be made (a file of that name is in the
! way)
      ierr = 0
      call write_gyre_pulse(n_points, pts, 1.0d33, 7.0d10, 3.8d33, 'not_a_dir/x.GYRE', ierr)
! (the open's own reason, the operating system's)
      call report("writes: a file whose directory cannot be made is an error", ierr == 1 .and. &
           captured(terminal_file, "not_a_dir/x.GYRE': Not a directory"))

! a full disk, for the files io closes
      ierr = 0
      call limit_file_size(limit_bytes)
      call write_gyre_pulse(n_points, pts, 1.0d33, 7.0d10, 3.8d33, 'full.GYRE', ierr)
      call limit_file_size(-1)
      call report("writes: a GYRE file cut short by a full disk is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write full.GYRE'))
      ierr = 0
      call limit_file_size(limit_bytes)
      call write_fgong_pulse(s, n_points, pts, 1.0d33, 7.0d10, 3.8d33, 'full.FGONG', ierr)
      call limit_file_size(-1)
      call report("writes: an FGONG file cut short by a full disk is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write full.FGONG'))
      allocate(s%stitched%prof(n_prof_cols, n_points))
      s%stitched%prof = 1.0d0
      s%stitched%n_ext = n_points
      s%physics%profile_data_prefix = 'full_profile'
      call profile_output_init(s, profiles, ierr)
      call limit_file_size(limit_bytes)
      call write_profile(s, profiles, '', 1, ierr)
      call limit_file_size(-1)
      call report("writes: a profile cut short by a full disk is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write full_profile1.data'))
      deallocate(s%stitched%prof)
      s%stitched%n_ext = 0
      history = history_file()
      call history_output_init(s, history, 'full_', ierr)
      call limit_file_size(limit_bytes)
      call write_history_row(s, history, 0, ierr)
      row_ierr = ierr
      call close_history(history, ierr)
      call limit_file_size(-1)
      call report("writes: a history row cut short by a full disk is an error", row_ierr == 1 .and. &
           ierr == 1 .and. captured(terminal_file, 'cannot write full_history.data'))

! a full disk, for the files kept open for the run: the .mod model and
! the isochrone file at each write, the run log at each model
! (check_run_log)
      open(newunit=u, file='full.mod', status='replace', action='write')
      ierr = 0
      call limit_file_size(limit_bytes)
      call write_mod_model(s, u, replace=.true., ierr=ierr)
      call limit_file_size(-1)
      close(u)
      call report("writes: a .mod model cut short by a full disk is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write full.mod'))
! (card kind 0: output_run_header writes the isochrone header and no
! terminal line)
      open(newunit=isochrone_unit, file='full.iso', status='replace', action='write')
      s%control%isochrone_output_active = .true.
      s%current_card = 1
      s%card_kind(1) = 0
      ierr = 0
      call limit_file_size(limit_bytes)
      do i = 1, 30
         if (ierr == 0) call output_run_header(s, 1.0d0, ierr)
      end do
      call limit_file_size(-1)
      close(isochrone_unit)
      isochrone_unit = -1
      s%control%isochrone_output_active = .false.
      call report("writes: an isochrone file cut short by a full disk is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write full.iso'))
      test_log_unit = run_log_unit
      open(newunit=run_log_unit, file='full_run.log', status='replace', action='write')
      ierr = 0
      call limit_file_size(limit_bytes)
      do i = 1, size(deck_lines)
         call log_line(trim(deck_lines(i)))
      end do
      call check_run_log(ierr)
      call limit_file_size(-1)
      close(run_log_unit)
      run_log_unit = test_log_unit
      call report("writes: a run log cut short by a full disk is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write full_run.log'))
! a run log rewound and written again, shorter (rewind_short_file), is
! not cut short
      test_log_unit = run_log_unit
      open(newunit=run_log_unit, file='rewound_run.log', status='replace', action='write')
      do i = 1, size(deck_lines)
         call log_line(trim(deck_lines(i)))
      end do
      ierr = 0
      call check_run_log(ierr)
      rewind(run_log_unit)
      call log_line('one line')
      call check_run_log(ierr)
      close(run_log_unit)
      run_log_unit = test_log_unit
      call report("writes: a run log rewound and written shorter is not cut short", ierr == 0 .and. &
           .not. captured(terminal_file, 'rewound_run.log'))

! the copy of the deck (an inlist: one file, copied once)
      do i = 1, size(deck_lines)
         write(deck_lines(i), '(a,i0,a)') ' line_', i, ' = 1  ! a line of a deck of some kilobytes'
      end do
      call write_text_file('big_deck.nml', deck_lines)
      ierr = 0
      call limit_file_size(limit_bytes)
      call copy_inlists_used([character(len=12) :: 'big_deck.nml', 'big_deck.nml'], &
           'full_inlist_used', ierr)
      call limit_file_size(-1)
      call report("writes: a copy of the deck cut short by a full disk is an error", ierr == 1 .and. &
           captured(terminal_file, 'cannot write full_inlist_used'))

! Two runs that write one path (M174: several decks in one directory
! all write output/history.data). The other run's open, after this
! run's last row, truncates the file in place (gfortran's
! status='replace' keeps the file, O_TRUNC) or puts a new file at the
! path; this run's close is not a full disk.
      history = history_file()
      call history_output_init(s, history, 'shared_', ierr)
      call write_history_row(s, history, 0, ierr)
      call other_run_writes('shared_history.data', 'w', 'another run''s header')
      call close_history(history, ierr)
      call report("writes: a history file another run truncated is not a full disk", ierr == 0 .and. &
           .not. captured(terminal_file, 'shared_history.data'))
      history = history_file()
      call history_output_init(s, history, 'replaced_', ierr)
      call write_history_row(s, history, 0, ierr)
      call other_run_replaces('replaced_history.data', 'another run''s header')
      call close_history(history, ierr)
      call report("writes: a history file another run replaced is not a full disk", ierr == 0 .and. &
           .not. captured(terminal_file, 'replaced_history.data'))
! a file io closes, its path taken by another file before the close:
! the file measured is the one this unit wrote
      ierr = 0
      call open_output_file('replaced.txt', u, ierr)
      ios = 0
      do i = 1, size(deck_lines)
         if (ios == 0) write(u, '(a)', iostat=ios, iomsg=msg) trim(deck_lines(i))
      end do
      call other_run_replaces('replaced.txt', 'another run''s line')
      call close_output_file(u, 'replaced.txt', ios, msg, ierr)
      call report("writes: a closed file is measured, not the one now at its path", ierr == 0 .and. &
           .not. captured(terminal_file, 'replaced.txt'))
! the size check keeps the unit's place: another writer has made the
! file longer than this unit's writes when they are checked, and the
! next line still follows this unit's last one
      ierr = 0
      call open_output_file('place.txt', u, ierr)
      write(u, '(a)') 'line 1'
      call flush_output_file(u, 'place.txt', 0, ' ', ierr)
      call other_run_writes('place.txt', 'a', 'a longer line that another run appended')
      write(u, '(a)') 'line 2'
      call flush_output_file(u, 'place.txt', 0, ' ', ierr)
      write(u, '(a)') 'line 3'
      call close_output_file(u, 'place.txt', 0, ' ', ierr)
      open(newunit=u, file='place.txt', status='old', action='read')
      read(u, '(a)') lines
      close(u)
      call report("writes: the size check keeps the unit's place in its file", ierr == 0 .and. &
           all(lines == [character(len=8) :: 'line 1', 'line 2', 'line 3']))
end subroutine writes_check

! ---------------------------------------------------------------
subroutine fgong_check()
      integer, parameter :: n_points = 4
      character(len=*), parameter :: fgong_value = '(1P,X,E26.18E3)'
      real(dp) :: pts(n_pulse_cols, n_points), glob(15)
      integer :: u, ierr, i

      call fill_mod_model(s)
      pts = 0.5d0
      s%card_initial_x(1) = 0.7123456789012345d0
      s%card_initial_z(1) = 0.0187654321098765d0
      s%current_card_envelope_x = 0.0d0
      s%current_card_envelope_z = 0.0d0
      s%mixing_length_alpha = 1.8518851851851851d0
      ierr = 0
      call write_fgong_pulse(s, n_points, pts, 1.0d33, 7.0d10, 3.8d33, 'globals.FGONG', ierr)
      open(newunit=u, file='globals.FGONG', status='old', action='read')
      do i = 1, 5
         read(u,*)
      end do
      read(u,*) glob
      close(u)
      call report("FGONG: globals 4-6 the run's initial Z, X and the mixing length", &
           ierr == 0 .and. glob(4) == through(fgong_value, s%card_initial_z(1)) .and. &
           glob(5) == through(fgong_value, s%card_initial_x(1)) .and. &
           glob(6) == through(fgong_value, s%mixing_length_alpha))
end subroutine fgong_check

end program test_io
