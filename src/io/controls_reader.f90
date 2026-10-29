!----------------------------------------------------------------------
! controls_reader
!----------------------------------------------------------------------
! Reads a deck into star%ctrl and star%job. Each read starts from the
! pristine defaults: star%ctrl = controls_state() (generated from
! controls_buffer's declarations), seeded into the buffer
! (seed_controls_buffer); read_input reads the deck over the buffer,
! and store_controls_to_star copies the result back. The read's steps
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
      use star_info_lib, only: star_info, controls_state, control_nml_override, &
           physics_nml_override
      use controls_sync, only: seed_controls_buffer, store_controls_to_star
      use yrec_version, only: yrec_version_string, git_hash_string
      use io_units, only: check_io
! the generated includes name every namelist variable
      use legacy_controls
      use run_files, only: assign_unit_numbers, expand_paths, open_run_files
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
      s%ctrl = controls_state()
      call set_solver_diagnostics(s%ctrl%report_solver_diagnostics)
      call seed_controls_buffer(s)
      call read_input(s, ierr)
      if (ierr /= 0) return
      call store_controls_to_star(s)
      call set_solver_diagnostics(s%ctrl%report_solver_diagnostics)

      return
end subroutine read_controls

!----------------------------------------------------------------------
! read_input
!----------------------------------------------------------------------
! Reads the deck: a single-file inlist (&control and &physics) or a
! legacy CONTROL/PHYSICS namelist pair (yrec8.nml1 and yrec8.nml2, or
! the two files named on the command line). Then copies the values to
! their canonical names, opens the run's files, derives the remaining
! controls and turns the kind cards into the run plan.
subroutine read_input(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

      ierr = 0
      call set_legacy_defaults
      call assign_unit_numbers
      call load_table_paths(s)

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
      if (ierr /= 0) return
      call adopt_canonical_names(s)
      call map_user_inputs(s)
! check_controls, below, reads the controls from star%ctrl and star%job
      call store_controls_to_star(s)
      call set_solver_diagnostics(s%ctrl%report_solver_diagnostics)
      call expand_paths
      call open_run_files(s, ierr)
      call store_table_paths(s)
      if (ierr /= 0) return
      call check_controls(s, ierr)
      if (ierr /= 0) return
      call warn_inert_controls
      call derive_controls(s)
      call interpret_kind_cards(s)

end subroutine read_input


! ---------------------------------------------------------------
! The table paths between star%job and their namelist variables. The
! load comes first, so a deck that names no table keeps the previous
! read's path; the store follows the read and the path expansion.
subroutine load_table_paths(s)
      type(star_info), intent(inout) :: s
      falex06 = s%job%alex06_table_path
      fallard = s%job%allard_table_path
      fatm = s%job%atm_table_path
      ffermi = s%job%fermi_table_path
      fkur = s%job%kurucz_table_path
      fkur2 = s%job%kurucz_table2_path
      flaol = s%job%laol_table_path
      flaol2 = s%job%laol_table2_path
      fliv95 = s%job%opal95_table_path
      flldat = s%job%opal92_table_path
      fmhd1 = s%job%zams_a_table_path
      fmhd2 = s%job%zams_b_table_path
      fmhd3 = s%job%zams_c_table_path
      fmhd4 = s%job%centre1_table_path
      fmhd5 = s%job%centre2_table_path
      fmhd6 = s%job%centre3_table_path
      fmhd7 = s%job%centre4_table_path
      fmhd8 = s%job%centre5_table_path
      fopal2 = s%job%opal92_table2_path
      fpurez = s%job%pure_z_table_path
      fscvh = s%job%scv_h_table_path
      fscvhe = s%job%scv_he_table_path
      fscvz = s%job%scv_z_table_path
      opecalex = s%job%alex95_table_paths
end subroutine load_table_paths

subroutine store_table_paths(s)
      type(star_info), intent(inout) :: s
      s%job%alex06_table_path = falex06
      s%job%allard_table_path = fallard
      s%job%atm_table_path = fatm
      s%job%fermi_table_path = ffermi
      s%job%kurucz_table_path = fkur
      s%job%kurucz_table2_path = fkur2
      s%job%laol_table_path = flaol
      s%job%laol_table2_path = flaol2
      s%job%opal95_table_path = fliv95
      s%job%opal92_table_path = flldat
      s%job%zams_a_table_path = fmhd1
      s%job%zams_b_table_path = fmhd2
      s%job%zams_c_table_path = fmhd3
      s%job%centre1_table_path = fmhd4
      s%job%centre2_table_path = fmhd5
      s%job%centre3_table_path = fmhd6
      s%job%centre4_table_path = fmhd7
      s%job%centre5_table_path = fmhd8
      s%job%opal92_table2_path = fopal2
      s%job%pure_z_table_path = fpurez
      s%job%scv_h_table_path = fscvh
      s%job%scv_he_table_path = fscvhe
      s%job%scv_z_table_path = fscvz
      s%job%alex95_table_paths = opecalex
      s%job%conductive_table_path = fcondopacp
end subroutine store_table_paths

! ---------------------------------------------------------------
! Read the deck. A file that holds both a &control and a &physics group
! is an inlist, read by read_inlist; anything else is a legacy pair, the
! CONTROL namelist in the first file and the PHYSICS namelist in the
! second. An inlist from before 2026-10, with &star_job and &controls
! groups, is refused with the command that renames its groups.
subroutine read_namelist_files(ierr)
      integer, intent(inout) :: ierr
      logical :: single_inlist
      integer :: ios
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
      open(unit=standard_unit, file=control_nml_file, status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the CONTROL namelist file', control_nml_file, ierr)
      if (ierr /= 0) return
      read(unit=standard_unit, nml=control, iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'read the CONTROL namelist from', control_nml_file, ierr, standard_unit)
      if (ierr /= 0) return
      close(standard_unit)
      open(unit=run_unit, file=physics_nml_file, status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the PHYSICS namelist file', physics_nml_file, ierr)
      if (ierr /= 0) return
      read(unit=run_unit, nml=physics, iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'read the PHYSICS namelist from', physics_nml_file, ierr, run_unit)
      if (ierr /= 0) return
      close(run_unit)
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
      integer :: ios
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
! The legacy names whose home is in the star, and what the copies
! do not give.
      s%solar_luminosity_cgs = clsun
      s%solar_radius_cgs = crsun
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
