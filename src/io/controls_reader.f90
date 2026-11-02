!----------------------------------------------------------------------
! controls_reader
!----------------------------------------------------------------------
! Reads a deck into star%control and star%physics, the records of the
! inlist's two namelist groups. Each read starts from the pristine
! defaults (control_group_controls() and physics_group_controls(),
! generated from controls_buffer's declarations), seeded into the buffer
! (seed_controls_buffer); read_input reads the deck over the buffer,
! straight into the controls' code names (an inlist through read_inlist,
! a legacy pair through legacy_controls' read_legacy_deck), the steps
! below turn it into the controls the run uses, store_controls_to_star
! copies the result to the star, and arm_run_copies sets the run's
! working copies of the controls it changes (the run plan, the latches;
! the registry's run_copy column) from them. The read never reads the
! star's run state. The read opens no file of the run (run_yrec's
! open_run_files does, after it), but for the run log after a controls
! error (open_run_log_after_error, so that the log holds the error); its
! messages reach the run log when that opens (run_messages). The read's
! other steps live in user_inputs (map_user_inputs), controls_check (the
! errors and warnings), derived_controls and run_cards. Errors return via
! ierr.
module controls_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use run_messages, only: check_io, report_warning, discard_pending_log, &
           set_solver_diagnostics
      use user_inputs, only: map_user_inputs
      use star_info_lib, only: star_info, control_nml_override, physics_nml_override
      use controls_sync, only: seed_controls_buffer, store_controls_to_star, &
           arm_run_copies
      use yrec_version, only: yrec_version_string, git_hash_string
      use controls_buffer, only: control_nml_file, physics_nml_file, report_solver_diagnostics, &
           log_output_file
      use legacy_controls, only: read_legacy_deck
      use run_files, only: assign_unit_numbers, expand_value, open_run_log_after_error
      use derived_controls, only: derive_controls
      use controls_check, only: check_controls, warn_inert_controls
      use run_cards, only: interpret_kind_cards
      implicit none
      private
      public :: read_controls, nml_file_has_group

! A value other than the one given, of its type (read_inlist's first read
! of the aliased names, inlist_read.inc): for the types of today's aliases
! (gen_inlist asks for another when an alias of a new type appears).
      interface other_than
         module procedure other_than_real, other_than_text
      end interface other_than

contains

!   star reads: none
!   star writes: none
subroutine read_controls(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

      ierr = 0
      call discard_pending_log
! the declared defaults, whatever an earlier read left on s
      call seed_controls_buffer(s)
      call set_solver_diagnostics(report_solver_diagnostics)
! the deck, under the code names, and the readable inputs mapped
      call read_input(ierr)
      if (ierr /= 0) return
      call map_user_inputs
! (a path that does not fit its control stops the read; the run log's
! path may be the one, so no run log is opened for it)
      call expand_paths(ierr)
      if (ierr /= 0) return
! the errors, the inert legacy controls, the derived controls, the run
! plan
      call check_controls(ierr)
      if (ierr /= 0) then
! a controls error: the run log still gets the read's lines and the
! error, at the deck's path (expanded above)
         call open_run_log_after_error(log_output_file)
         return
      end if
      call warn_inert_controls
      call derive_controls
      call interpret_kind_cards
      call store_controls_to_star(s)
      call arm_run_copies(s)
      call set_solver_diagnostics(report_solver_diagnostics)
end subroutine read_controls

!----------------------------------------------------------------------
! read_input
!----------------------------------------------------------------------
! Reads the deck over the defaults: a single-file inlist (&control and
! &physics) or a legacy CONTROL/PHYSICS namelist pair (yrec8.nml1 and
! yrec8.nml2, or the two files named on the command line).
subroutine read_input(ierr)
      integer, intent(out) :: ierr

      ierr = 0
      call assign_unit_numbers

      print *,''
      print *,'Yale Rotating Evolution Code - YREC, v',yrec_version_string(1:len_trim(yrec_version_string)),' (',git_hash_string(1:len_trim(git_hash_string)),')'

! The deck: the first command-line argument (default yrec8.nml1), or
! the override an embedding host (pyyrec) sets, since argv is the host's.
      if (control_nml_override /= ' ') then
         control_nml_file = control_nml_override
      else
         call getarg(1, control_nml_file)
         if (control_nml_file(1:2) == ' ') control_nml_file = 'yrec8.nml1'
      end if
      print *, ' '
      call read_namelist_files(ierr)
end subroutine read_input


! ---------------------------------------------------------------
! Read the deck. A file that holds both a &control and a &physics group
! is an inlist, read by read_inlist; anything else is a legacy pair, the
! CONTROL namelist in the first file and the PHYSICS namelist in the
! second, read by read_legacy_deck. An inlist from before 2026-10, with
! &star_job and &controls groups, is refused with the command that
! renames its groups.
subroutine read_namelist_files(ierr)
      integer, intent(inout) :: ierr
      logical :: single_inlist

      if (nml_file_has_group(control_nml_file, 'star_job')) then
         write(*,'(a)') ' '//trim(control_nml_file)//' names its groups &star_job and &controls;'
         write(*,'(a)') ' they are &control and &physics now. Rename them with'
         write(*,'(a)') '   python3 src/tools/upgrade_inlist.py '//trim(control_nml_file)
         ierr = 1
         return
      end if
      single_inlist = nml_file_has_group(control_nml_file, 'control') .and. &
           nml_file_has_group(control_nml_file, 'physics')

      if (physics_nml_override /= ' ') then
         physics_nml_file = physics_nml_override
      else
         call getarg(2, physics_nml_file)
         if (physics_nml_file(1:2) == ' ') physics_nml_file = 'yrec8.nml2'
      end if
! a single-file inlist is both files (copy_inlists_used copies it once)
      if (single_inlist) physics_nml_file = control_nml_file

      if (single_inlist) then
         write(*,*) 'inlist           :  ',control_nml_file(1:len_trim(control_nml_file))
         call read_inlist(control_nml_file, ierr)
      else
         write(*,*) 'CONTROL namelist :  ',control_nml_file(1:len_trim(control_nml_file))
         write(*,*) 'PHYSICS namelist :  ',physics_nml_file(1:len_trim(physics_nml_file))
         call read_legacy_deck(control_nml_file, physics_nml_file, ierr)
      end if
end subroutine read_namelist_files

! ---------------------------------------------------------------
! Read a single-file inlist, inlist_file: its &control and &physics
! groups, under the inlist names (defaults/control.defaults,
! defaults/physics.defaults), into the controls' code names. Both
! includes are generated by tools/gen_inlist.py: inlist_decl.inc binds
! the groups to the code names (and to a local, copied into its code
! name after the read, for an inlist name that is not its code name),
! and inlist_read.inc reads them (twice, for the aliases).
subroutine read_inlist(inlist_file, ierr)
      include 'inlist_decl.inc'
      character(len=*), intent(in) :: inlist_file
      integer, intent(inout) :: ierr
      integer :: ios, namelist_unit
      character(len=256) :: msg

      include 'inlist_read.inc'
end subroutine read_inlist

! ---------------------------------------------------------------
! The {YREC_INPUT}, {YREC_START} and {YREC_OUTPUT} placeholders in
! every path control, expanded: the generated expand_paths.inc calls
! run_files' expand_value on each control whose registry type is path,
! naming it by its inlist name; one that does not fit returns ierr = 1.
subroutine expand_paths(ierr)
      include 'expand_paths.inc'
end subroutine expand_paths

! ---------------------------------------------------------------
! A value other than x, of its type, element by element: read_inlist
! seeds the aliased names with it, so that an element the deck names
! differs from its seed in one of the two reads.
elemental function other_than_real(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = merge(1.0_dp, 0.0_dp, x == 0.0_dp)
end function other_than_real

elemental function other_than_text(x) result(y)
      character(len=*), intent(in) :: x
      character(len=len(x)) :: y
      y = merge('b', 'a', x(1:1) == 'a')
end function other_than_text

! Does this namelist file contain the group &<gname> (in any case, and
! not merely a group whose name begins with gname)?
logical function nml_file_has_group(fname, gname)
      character(len=*), intent(in) :: fname, gname
      character(len=512) :: probe_line
      character(len=1) :: next_char
      integer :: probe_unit, ios, n

      nml_file_has_group = .false.
      n = len_trim(gname) + 1
      open(newunit=probe_unit, file=fname, status='OLD', &
           action='READ', iostat=ios)
      if (ios /= 0) return
      do
         read(probe_unit, '(a)', iostat=ios) probe_line
         if (ios /= 0) exit
         probe_line = adjustl(probe_line)
         next_char = probe_line(n+1:n+1)
         if (lowercase(probe_line(1:n)) == '&' // lowercase(trim(gname)) &
              .and. (next_char == ' ' .or. next_char == '/' .or. &
              next_char == achar(9))) then
            nml_file_has_group = .true.
            exit
         end if
      end do
      close(probe_unit)
end function nml_file_has_group

! text with its ASCII capitals lowered
pure function lowercase(text) result(lowered)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lowered
      integer :: i, code

      lowered = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) &
              lowered(i:i) = achar(code + 32)
      end do
end function lowercase

end module controls_reader
