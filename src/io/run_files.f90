!----------------------------------------------------------------------
! run_files
!----------------------------------------------------------------------
! The run's files: their fixed unit numbers, the {YREC_...}
! placeholders in the path controls expanded (expand_value), the output
! directory created with a verbatim copy of the inlists (inlist_used),
! and the files opened.
module run_files
      use io_units
      use yrec_output, only: output_init_mesa
      use legacy_controls
      implicit none
      private
      public :: assign_unit_numbers, expand_paths, open_run_files, expand_value

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
! The {YREC_INPUT}, {YREC_START} and {YREC_OUTPUT} placeholders in the
! path controls, expanded (expand_value).
subroutine expand_paths
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
end subroutine expand_paths

! ---------------------------------------------------------------
! The output directory (the run log's directory) created and the
! inlists copied into it, then the run's files opened: the starting
! model, the final model (the restart file, and the solver's recovery
! source when a step diverges), the run log and history (output_init_mesa),
! and the Monte Carlo and isochrone files when those are on.
subroutine open_run_files(ierr)
      integer, intent(inout) :: ierr
      character(len=256) :: shell_cmd, msg
      integer :: i, last_slash_idx, ios

      last_slash_idx = 0
      do i = len_trim(fshort), 1, -1
          if (fshort(i:i) .eq. '/') then
              last_slash_idx = i
              exit
          endif
      end do
      shell_cmd = 'mkdir -p '
      shell_cmd(len_trim(shell_cmd)+2:len_trim(fshort(1:last_slash_idx))+len_trim(shell_cmd)) = fshort(1:last_slash_idx)
      print *,"OUTPUT placed in :  ",fshort(1:last_slash_idx)
      print *, ''
      call system(shell_cmd)
      call copy_inlists_used(trim(inlist_used_file), ierr)
      if (ierr /= 0) return

      open(unit=first_unit,file=ffirst,form='FORMATTED',status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the starting model', ffirst, ierr)
      if (ierr /= 0) return
      open(unit=last_model_unit,file=flast,form='FORMATTED',status='UNKNOWN', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the final model file', flast, ierr)
      if (ierr /= 0) return
      call output_init_mesa(fshort, ierr)
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
! replaced by that environment variable's value, or by the default
! below when it is unset.
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
