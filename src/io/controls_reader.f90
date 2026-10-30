!----------------------------------------------------------------------
! controls_reader
!----------------------------------------------------------------------
! Reads a deck into star%ctrl and star%job. Each read starts from the
! pristine defaults: star%ctrl = controls_state() and star%job =
! star_job() (generated from controls_buffer's declarations), seeded
! into the buffer (seed_controls_buffer); read_input reads the deck over
! the buffer, the steps below turn it into the controls the run uses,
! and store_controls_to_star copies the result back. The read opens no
! file of the run (run_yrec's open_run_files does, after it); its
! messages reach the run log when that opens (io_units). The read's steps
! live in legacy_controls (the namelist variables and their defaults),
! run_files, controls_check (the errors and warnings), derived_controls
! and run_cards; the copies between the legacy and the canonical names
! are generated from the registry (adopt_canonical.inc,
! legacy_mirror.inc). Errors return via ierr.
module controls_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use run_log, only: set_solver_diagnostics
      use user_inputs, only: map_user_inputs
      use net_lib, only: cross_section_scales_from_sstandard
      use star_info_lib, only: star_info, controls_state, star_job, &
           control_nml_override, physics_nml_override
      use controls_sync, only: seed_controls_buffer, store_controls_to_star
      use yrec_version, only: yrec_version_string, git_hash_string
      use io_units, only: check_io, discard_pending_log
! the generated includes name every namelist variable
      use legacy_controls
      use run_files, only: assign_unit_numbers, expand_paths
      use derived_controls, only: derive_controls
      use controls_check, only: check_controls, warn_inert_controls
      use run_cards, only: interpret_kind_cards
      implicit none
      private
      public :: read_controls

contains

subroutine read_controls(s, ierr)
      type(star_info), intent(inout) :: s

      integer, intent(out) :: ierr

      ierr = 0
      call discard_pending_log
! the declared defaults, whatever an earlier read or run left on s
      s%ctrl = controls_state()
      s%job = star_job()
      call seed_controls_buffer(s)
      call set_solver_diagnostics(report_solver_diagnostics)
! the deck, under the canonical names, and the readable inputs mapped
      call read_input(ierr)
      if (ierr /= 0) return
      call adopt_canonical_names(s)
      call map_user_inputs(s)
      call expand_paths
! the errors, the inert legacy controls, the derived controls, the run
! plan
      call check_controls(ierr)
      if (ierr /= 0) return
      call warn_inert_controls
      call derive_controls(s)
      call interpret_kind_cards(s)
      call store_controls_to_star(s)
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
      call set_legacy_defaults
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
! second. An inlist from before 2026-10, with &star_job and &controls
! groups, is refused with the command that renames its groups.
subroutine read_namelist_files(ierr)
      integer, intent(inout) :: ierr
      logical :: single_inlist
      integer :: ios, namelist_unit
      character(len=256) :: msg

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
      else
      write(*,*) 'CONTROL namelist :  ',control_nml_file(1:len_trim(control_nml_file))
      write(*,*) 'PHYSICS namelist :  ',physics_nml_file(1:len_trim(physics_nml_file))
      end if

      if (single_inlist) then
      call read_inlist(ierr)
      else
      open(newunit=namelist_unit, file=control_nml_file, status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the CONTROL namelist file', control_nml_file, ierr)
      if (ierr /= 0) return
      read(unit=namelist_unit, nml=control, iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'read the CONTROL namelist from', control_nml_file, ierr, namelist_unit)
      if (ierr /= 0) return
      close(namelist_unit)
      open(newunit=namelist_unit, file=physics_nml_file, status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the PHYSICS namelist file', physics_nml_file, ierr)
      if (ierr /= 0) return
      read(unit=namelist_unit, nml=physics, iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'read the PHYSICS namelist from', physics_nml_file, ierr, namelist_unit)
      if (ierr /= 0) return
      close(namelist_unit)
      end if
end subroutine read_namelist_files

! ---------------------------------------------------------------
! Read a single-file inlist: its &control and &physics groups, under the
! new names (defaults/control.defaults, defaults/physics.defaults).
! Both files are generated by tools/gen_inlist.py: inlist_new_decl.inc
! declares the groups -- named like legacy_controls' CONTROL and PHYSICS
! groups, which they hide here -- and the new-name variables that have
! no legacy twin; inlist_new_read.inc seeds the new names from the
! legacy defaults, reads the groups and copies the values back onto the
! legacy names.
subroutine read_inlist(ierr)
      integer, intent(inout) :: ierr
      integer :: ios, namelist_unit
      character(len=256) :: msg
      include 'inlist_new_decl.inc'

      include 'inlist_new_read.inc'
end subroutine read_inlist

! ---------------------------------------------------------------
! Copy each control read under a legacy name into its canonical home
! (controls_buffer or star). map_user_inputs,
! which runs next, reads use_new_nuclear_rates and the s0_* values.
subroutine adopt_canonical_names(s)
      type(star_info), intent(inout) :: s
      include 'adopt_canonical.inc'
! What the copies do not give.
! lnewnuc = .false. takes the cross sections from sstandard (ratios to
! Bahcall & Ulrich 1988, mapped onto the reaction and extra scales by
! net); map_user_inputs sets them from the s0_* inputs otherwise. An
! inlist gets the default ratios.
      if (.not. lnewnuc) call cross_section_scales_from_sstandard(sstandard, &
           s%cross_section_scale, s%extra_cross_section_scale)
! the atmosphere choice as read; the run may change star%job%atm_choice
      s%atm_choice_initial = kttau
      s%use_ttau_relation = .false.
end subroutine adopt_canonical_names

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
