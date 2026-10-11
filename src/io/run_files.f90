!----------------------------------------------------------------------
! run_files
!----------------------------------------------------------------------
! The run's files: the {YREC_...} placeholders in path controls
! expanded (expand_value), the output directory created and the output
! files opened, and the inlists copied verbatim into the output
! directory as inlist_used (formerly part of parmin).
module run_files
      use star_info_lib, only: star, control_nml_override, physics_nml_override
      use io_units
      use yrec_output, only: output_init_mesa
      use legacy_controls
      implicit none
      private
      public :: resolve_output_mode_and_paths, derive_options_and_open_files, expand_value

contains

! ---------------------------------------------------------------
! Expand environment variables in every configured file path.
subroutine resolve_output_mode_and_paths
! 2026 log redesign: the full namelist echo into the run log is
! replaced by the verbatim inlist copy written to the output
! directory (inlist_used -- see copy_inlists_used below).

! Post-process all CONTROL namelist vars that hold path values.
! Expand any placeholders found in the string with the value taken from a
! corresponding environment variable, if one is defined.
      call expand_value(falex06)
      call expand_value(fallard)
      call expand_value(fatm)
      call expand_value(fcondopacp)
      call expand_value(ffermi)
      call expand_value(ffirst)
      call expand_value(flast)
      call expand_value(fliv95)
      call expand_value(fopale06)
      call expand_value(fpurez)
      call expand_value(fscomp)
      call expand_value(fscvh)
      call expand_value(fscvhe)
      call expand_value(fscvz)
      call expand_value(fshort)
      call expand_value(inlist_used_file)
end subroutine resolve_output_mode_and_paths

! ---------------------------------------------------------------
! Output directory creation, derived options (track-name parsing,
! isotope switches, EOS table selection) and the output-unit opens
! (one open path for every run since the 2026 use_legacy_output
! retirement). Sets ierr nonzero (config error) if semiconvection
! and overshoot are both enabled; the caller returns immediately in
! that case.
subroutine derive_options_and_open_files(ierr)
      integer, intent(inout) :: ierr
! Create the output directory if it doesn't already exist. It is
! taken from the directory part of the .short log path (FSHORT) --
! every deck keeps its outputs together. (Historically this came
! from FTRACK, retired with the .track file.)
      shell_cmd = 'mkdir -p '
! find index of last '/' char. Use that to snip out the directory name.
      do i = len_trim(fshort), 1, -1
          if (fshort(i:i) .eq. '/') then
              last_slash_idx = i
              exit
          endif
      end do
      shell_cmd(len_trim(shell_cmd)+2:len_trim(fshort(1:last_slash_idx))+len_trim(shell_cmd)) = fshort(1:last_slash_idx)
      print *,"OUTPUT placed in :  ",fshort(1:last_slash_idx)
      print *, ''
      call system(shell_cmd)
      call copy_inlists_used(trim(inlist_used_file))

! G Somers 6/14, DEFINE SCALING COEFFICIENT FOR LI/BE CROSS SECTIONS
!          DEFAULT LI6 = 5.5 MeV b FROM FOWLER ET AL. 1967
!          DEFAULT LI7 = 52 keV b FROM ROLFS & KAVANAGH 1986
!           DEFAULT BE9(P,G)B10 (BE91) = 1.1 keV b INFERRED FROM FOWLER 1988
!           DEFAULT BE9(P,D)2HE4 (BE92) = 15,000 keV b FROM FOWLER ET AL. 1967
!           DEFAULT BE9(P,A)LI6 (BE93) = 15,000 keV b FROM FOWLER ET AL. 1967
      if (lxli6) then
          sli6 = xsli6/5.5d3
      endif
      if (lxli7) then
          sli7 = xsli7/5.2d1
      endif
      if (lxbe91) then
          sbe91 = xsbe91/1.1d0
      endif
      if (lxbe92) then
          sbe92 = xsbe92/1.5d4
      endif
      if (lxbe93) then
          sbe93 = xsbe93/1.5d4
      endif
! sli6/sli7/sbe91/sbe92/sbe93 must keep their NAMELIST spelling (see
! declaration above), so copy into const_lib's canonical names here.
      li6_rate_scale = sli6
      li7_rate_scale = sli7
      be9_pg_rate_scale = sbe91
      be9_pd_rate_scale = sbe92
      be9_palpha_rate_scale = sbe93
! G Somers END
! MHP 8/25 open relevant table as well as ensuring that only one is selected
!  Disable Older OPAL EOS's if a newer one is specified
!  (2026-10: eos_init opens the selected table, star%job%opal_eos_table_path)
      star%job%opal_eos_table_path = ""
      if (lopale06) then
          lopale01 = .false.
          lopale = .false.
         star%job%opal_eos_table_path = fopale06
      endif
      if (lopale01) then
         star%job%opal_eos_table_path = fopale01
         lopale = .false.
      else if(lopale) then
         star%job%opal_eos_table_path = fopale
      endif
! 3/09 Disable older Alexander opacities if a newer one is specified
      if(lalex06)then
         lalex95 = .false.
         lkur90 = .false.
      endif
      if(lalex95)then
         lkur90 = .false.
      endif
! lopale/lopale01/lalex95/lkur90 can all be overridden above, after
! they were already copy-assigned into their const_lib canonical names
! (use_opal95_eos/use_opal2001_eos/use_alex95_tables/
! use_kurucz90_tables) earlier in this subroutine -- re-sync now so the
! const_lib values reflect any override from this block, not just the
! raw namelist read.
      use_opal95_eos = lopale
      use_opal2001_eos = lopale01
      use_alex95_tables = lalex95
      use_kurucz90_tables = lkur90

! 2026 use_legacy_output retirement: one open path for every run.
! (the LDEBUG/.debug stream open that sat here is retired)
      open(unit=first_unit,file=ffirst,form='FORMATTED',status='OLD')
! last_model_unit (the .mod model file) is written every model in BOTH output
! modes: it is the restart file AND the solver's divergence/timestep-
! cutting recovery source, so it must always be connected (2026: the
! legacy-only guard here left MESA-mode runs with an unconnected
! unit -- no restart file and broken divergence recovery).
      open(unit=last_model_unit,file=flast,form='FORMATTED',status='UNKNOWN')
! open the run log (fshort's dir; .short renamed .log) and the
! history stream, parse the column selections
      call output_init_mesa(fshort, ierr)
      if (ierr /= 0) return
! MHP 6/98
! MHP 8/25 Moved call from main to here for opening dynamics_unit
      if(lmonte)then
         open(unit=dynamics_unit,file=fdyn,form='FORMATTED',status='OLD')
         open(monte_carlo_unit1, file=fmonte1,status='UNKNOWN',form='FORMATTED')
         open(monte_carlo_unit2, file=fmonte2,status='UNKNOWN',form='FORMATTED')
      endif
!     MHP 8/25 Moved opening of the EoS tables here, to avoid complicated passages of declared variables.
!     (2026-10: the conductive opacity table is read by kap_init from
!     star%job%conductive_table_path.)
      if(liso) then
         open(isochrone_file_unit, file=fiso,status='UNKNOWN', form='FORMATTED')
      endif
      if(lsemic)then
         if(lovstc.or.lovste.or.lovstm)then
            write(run_log_unit,2)lsemic,lovste,lovstc,lovstm
      2       format(1x,'ERROR IN SUBROUTINE PARMIN'/'SEMI-CONVECTION', &
           &  ' AND OVERSHOOT FLAGS BOTH TURNED ON'/'FLAGS LSEMIC',l2, &
           &  ' OVERSHOOT - CORE,ENVELOPE,INTERMEDIATE-',3l2/'RUN STOPPED')
            ! 2026 (phase five, step B): stop converted to ierr; run_yrec
            ! returns the error and the CLI wrapper (main) stops.
            ierr = 1
            return
         endif
      endif
end subroutine derive_options_and_open_files

! ---------------------------------------------------------------
! Write a verbatim copy of the namelist input file(s) into the
! output directory as "inlist_used" -- run provenance (2026 log
! redesign: replaces the STANDARD/CURRENT settings tables and the
! full namelist echoes the run log used to carry).
subroutine copy_inlists_used(dest)
      character(len=*), intent(in) :: dest
      integer :: src_unit, dst_unit
      character(len=4096) :: line
      integer :: ios, nfile
      character(len=256) :: sources(2)
      sources(1) = control_nml_file
      sources(2) = physics_nml_file
      open(newunit=dst_unit, file=dest, &
           form='FORMATTED', status='UNKNOWN')
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

!----------------------------------------------------------------------
! expand_value
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original parmin.f (this subroutine lived in the same source file);
! only variable names, source form, and comment style were updated.
!
! Replaces any defined "{YREC_XXX}" placeholder string in the passed
! namelist path variable with the value of the corresponding
! environment variable (falling back to a built-in default path if
! that environment variable is not set).
subroutine expand_value(path_value)

      implicit none

      integer, parameter :: n_env_vars = 3  ! Number of possible env. vars.
      character(len=256) :: path_value
      character(len=256) :: temp_value
      character(len=256) :: placeholder_names(n_env_vars)
      character(len=256) :: default_values(n_env_vars)
      character(len=256) :: placeholder
      integer :: env_val_len, orig_val_len, placeholder_len, default_len
      integer :: i
! Each placeholder (env var name enclosed in curly braces) to be
! supported for expansion in namelists, along with default value in the
! case where the var is not defined in the execution environment.
! The number of assignment pairs here must match the value of the
! parameter n_env_vars, above.
      placeholder_names(1) = "{YREC_INPUT}"
      default_values(1) = "../../input"

      placeholder_names(2) = "{YREC_START}"
      default_values(2) = "../../startmodels"

      placeholder_names(3) = "{YREC_OUTPUT}"
      default_values(3) = "output"

      do i=1, n_env_vars
        placeholder_len = len_trim(placeholder_names(i))
        default_len = len_trim(default_values(i))
        placeholder = placeholder_names(i)(1:placeholder_len)

        if (path_value(1:placeholder_len) .eq. placeholder) then    ! If placeholder
            call getenv(placeholder(2:placeholder_len-1), temp_value)
            env_val_len = len_trim(temp_value)
            orig_val_len = len_trim(path_value)
            if (env_val_len .ne. 0) then       ! If env var set
                temp_value(env_val_len+1:env_val_len+orig_val_len-placeholder_len) = path_value(placeholder_len+1:orig_val_len)
                  !print *,"Override: ",TEMP(1:LEN_TRIM(TEMP))
            else                          ! No env var. Use default path.
                temp_value(1:default_len) = default_values(i)(1:len_trim(default_values(i)))
                temp_value(default_len+1:default_len+orig_val_len) = path_value(placeholder_len+1:orig_val_len)
                  !print *,"default: ",TEMP(1:LEN_TRIM(TEMP))
            end if
            path_value = temp_value
        end if
      end do

end subroutine expand_value

end module run_files
