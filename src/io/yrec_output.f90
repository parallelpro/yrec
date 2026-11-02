!----------------------------------------------------------------------
! yrec_output
!----------------------------------------------------------------------
! The run's output, per job, in the output directory (the run log's):
!   <star_history_name>        one row per converged model, MESA's
!                              history.data layout (history_output); its
!                              profile_number column maps models to
!                              profiles.
!   <profile_data_prefix><N>.data  every profile_interval models, MESA's
!                              profile layout, zone 1 the surface
!                              (profile_output), with the pulsation file
!                              in pulse_format beside it
!                              (pulsation_output).
!   the .mod model             every model (model_files).
! The run log (run_files opens it) gets the progress line every
! terminal_interval models (run_log).
!
! The job's output state is one record, output_state, which run_yrec
! owns: output_open sets it up after the controls are read, evolve_step
! passes it to output_write_model every model, and output_close closes
! the history file at the end. The writers only read the star:
! compute_observables and the stitched model (core) compute every
! quantity first. Every write is checked (output_files): a failed one
! returns ierr = 1, which stops the run.
module yrec_output
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: ln10, solar_bolometric_magnitude
      use history_output, only: history_file, history_output_init, &
           write_history_row, close_history
      use profile_output, only: profile_selection, profile_output_init, &
           write_profile
      use run_log, only: progress_lines, log_model_line, log_final_model_line
      use pulsation_output, only: write_fgong_pulse, write_gsm_pulse, &
           write_gyre_pulse, gsm_supported
      use model_files, only: write_mod_model
      use output_columns, only: integer_field
      use math_lib
      implicit none
      private
      public :: output_state, output_open, output_close, output_run_header, &
           output_write_model, output_final_model_line
      public :: pulse_format_code, unknown_pulse_format, pulse_gyre, pulse_fgong, &
           pulse_gsm, interior_gyre_path

! The pulsation file formats pulse_format names (case-insensitive);
! unknown_pulse_format for any other name.
      integer, parameter :: unknown_pulse_format = 0, pulse_gyre = 1, &
           pulse_fgong = 2, pulse_gsm = 3

! One job's output state: the output directory, the profile counter
! (the N of profile<N>.data), the pulse format, the history file, the
! profile columns and the progress lines.
      type :: output_state
         character(len=256) :: out_dir = ' '
         integer :: profile_counter = 0
         integer :: pulse_format = pulse_gyre
         type(history_file) :: history
         type(profile_selection) :: profiles
         type(progress_lines) :: progress
      end type output_state

contains

! ---------------------------------------------------------------
! The job's output, after the controls are read: a history file a
! previous job left open is closed, the column selections are read,
! and the pulse format is checked, so that a bad one stops the run
! before its first model (an unknown name, or GSM on a build without
! HDF5).
!   star reads: physics%pulse_format, physics%write_pulse_flag
!   star writes: none
subroutine output_open(s, output, ierr)
      use star_info_lib, only: star_info
      use run_files, only: output_directory
      use run_messages, only: report_error
      type(star_info), intent(in) :: s
      type(output_state), intent(inout) :: output
      integer, intent(out) :: ierr

      ierr = 0
      output%out_dir = output_directory(s)
      output%profile_counter = 0
      output%progress = progress_lines()
      call history_output_init(s, output%history, output%out_dir, ierr)
      if (ierr /= 0) return
      call profile_output_init(s, output%profiles, ierr)
      if (ierr /= 0) return
      output%pulse_format = pulse_format_code(s%physics%pulse_format)
      if (s%physics%write_pulse_flag) then
         select case (output%pulse_format)
         case (unknown_pulse_format)
            call report_error("pulse_format = '" // trim(s%physics%pulse_format) // &
                 "' is not GYRE, FGONG or GSM")
            ierr = 1
            return
         case (pulse_gsm)
            if (.not. gsm_supported()) then
               call report_error('pulse_format = GSM requires an HDF5-enabled ' // &
                    'build: make clean && make USE_HDF5=1')
               ierr = 1
               return
            end if
         end select
      end if
end subroutine output_open

! ---------------------------------------------------------------
! The end of the job's output: the history file closed, and the run
! log's last lines checked on disk.
subroutine output_close(output, ierr)
      use run_files, only: check_run_log
      type(output_state), intent(inout) :: output
      integer, intent(inout) :: ierr
      call close_history(output%history, ierr)
      if (ierr == 0) call check_run_log(ierr)
end subroutine output_close

! ---------------------------------------------------------------
! The format a pulse_format value names: pulse_gyre, pulse_fgong or
! pulse_gsm, whatever its case; unknown_pulse_format otherwise.
pure integer function pulse_format_code(name)
      character(len=*), intent(in) :: name
      character(len=len(name)) :: upper
      integer :: i, c

      upper = name
      do i = 1, len(upper)
         c = iachar(upper(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) upper(i:i) = achar(c - 32)
      end do
      select case (trim(adjustl(upper)))
      case ('GYRE');  pulse_format_code = pulse_gyre
      case ('FGONG'); pulse_format_code = pulse_fgong
      case ('GSM');   pulse_format_code = pulse_gsm
      case default;   pulse_format_code = unknown_pulse_format
      end select
end function pulse_format_code

! ---------------------------------------------------------------
! One converged model's output: the progress line, the isochrone
! record, the profile and pulse files when due, the history row, the
! .mod model, and the interval GYRE file; then the run log checked on
! disk (run_files' check_run_log).
!   star reads: age_gyr, control%isochrone_output_active,
!     control%pulse_gyre_interval, h_shell_bot_index, has_h_shell,
!     log_g_surface, log_L_lsun, log_R_surface_rsun, log_Teff, m,
!     model_number, physics%write_profile_flag,
!     physics%write_pulse_flag, solar_luminosity_cgs, solar_radius_cgs,
!     xa
!   star writes: none
subroutine output_write_model(s, output, ierr)
      use star_info_lib, only: star_info, i_he4
      use run_files, only: last_model_unit, isochrone_unit, check_run_log
      type(star_info), intent(in) :: s
      type(output_state), intent(inout) :: output
      integer, intent(inout) :: ierr
      integer :: iprof, ios
      character(len=256) :: msg
      real(dp) :: age_yr, luminosity_erg_s, radius_cm, teff_k, &
           gravity_cgs, ycenter_local, he_core_mass_grams

! the progress line (every terminal_interval models)
      call log_model_line(s, output%progress)

! isochrone record: model no. (a field of 5, wider from model 100000 on),
! age (yr), L (erg/s), R (cm), Teff (K), g (cm/s**2), Ycenter, He-core
! mass (g)
      if (s%control%isochrone_output_active) then
         age_yr = s%age_gyr*1.0D9
         luminosity_erg_s = exp10(s%log_L_lsun)*s%solar_luminosity_cgs
         radius_cm = exp10(s%log_R_surface_rsun)*s%solar_radius_cgs
         teff_k = exp10(s%log_Teff)
         gravity_cgs = exp10(s%log_g_surface)
         ycenter_local = s%xa(i_he4,1)
         if (s%has_h_shell) then
            he_core_mass_grams = s%m(s%h_shell_bot_index-1)
         else
            he_core_mass_grams = 0.0D0
         end if
         write(isochrone_unit,'(1X, A, 1P7E17.8)',iostat=ios,iomsg=msg) integer_field(s%model_number, 5), &
              age_yr,luminosity_erg_s,radius_cm,teff_k,gravity_cgs, &
              ycenter_local,he_core_mass_grams
         call check_isochrone_write(ios, msg, ierr)
         if (ierr /= 0) return
      end if

! Profiles and pulse files share the trigger and the number (MESA's
! write_pulse_data_with_profile): every profile_interval models the
! counter advances and each enabled file is written, profile<N>.data
! and profile<N>.data.<format>. The history's profile_number column
! records N (0 for a model without one).
      iprof = 0
      if (profile_write_due(s)) then
         output%profile_counter = output%profile_counter + 1
         iprof = output%profile_counter
         if (s%physics%write_profile_flag) call write_profile(s, output%profiles, &
              output%out_dir, iprof, ierr)
         if (ierr /= 0) return
         if (s%physics%write_pulse_flag) call write_pulse(s, output, iprof, ierr)
         if (ierr /= 0) return
      end if
      call write_history_row(s, output%history, iprof, ierr)
      if (ierr /= 0) return
! the .mod model: the restart file, and the recovery source of a step
! that fails
      call write_mod_model(s, last_model_unit, replace=.true., ierr=ierr)
      if (ierr /= 0) return

! the interval GYRE file (independent of write_pulse_flag)
      if (s%control%pulse_gyre_interval > 0 .and. &
          mod(s%model_number,s%control%pulse_gyre_interval) == 0) &
           call write_interior_gyre(s, output, ierr)
      if (ierr /= 0) return
      call check_run_log(ierr)
end subroutine output_write_model

! ---------------------------------------------------------------
! The interior-only GYRE file, interior_gyre_path's
! profile<model>.data.GYRE: the stitched model's interior points, which
! are the first nz of its point set (centre to the fitting point), with
! the outermost interior shell's mass, radius and luminosity as M_star,
! R_star and L_star.
!   star reads: model_number, nz, physics%profile_data_prefix, stitched
!   star writes: none
subroutine write_interior_gyre(s, output, ierr)
      use star_info_lib, only: star_info
      use stitched_model, only: ipul_m, ipul_r, ipul_L
      type(star_info), intent(in) :: s
      type(output_state), intent(in) :: output
      integer, intent(inout) :: ierr
      character(len=:), allocatable :: path

      path = interior_gyre_path(output%out_dir, s%physics%profile_data_prefix, s%model_number)
      call write_gyre_pulse(s%nz, s%stitched%pulse(:,1:s%nz), s%stitched%pulse(ipul_m,s%nz), &
           s%stitched%pulse(ipul_r,s%nz), s%stitched%pulse(ipul_L,s%nz), path, ierr)
end subroutine write_interior_gyre

! ---------------------------------------------------------------
! The interior GYRE file's path: the output directory, the profiles'
! prefix, the model number zero-padded to five digits (as it is written
! in full from model 100000 on, and for a negative one), '.data.GYRE'.
pure function interior_gyre_path(out_dir, prefix, model_number) result(path)
      character(len=*), intent(in) :: out_dir, prefix
      integer, intent(in) :: model_number
      character(len=:), allocatable :: path
      character(len=24) :: model_text

      write(model_text,'(i0)') model_number
      if (model_number >= 0 .and. len_trim(model_text) < 5) write(model_text,'(i5.5)') model_number
      path = trim(out_dir) // trim(prefix) // trim(model_text) // '.data.GYRE'
end function interior_gyre_path

! ---------------------------------------------------------------
! The end of a card: the progress line of its last converged model,
! whatever the terminal_interval phase.
!   star reads: none
!   star writes: none
subroutine output_final_model_line(s, output)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(output_state), intent(inout) :: output
      call log_final_model_line(s, output%progress)
end subroutine output_final_model_line

! ---------------------------------------------------------------
! Is a profile/pulse write due this model? A positive cadence, at
! least one of the two enabled, and the model number on the cadence.
!   star reads: model_number, physics%profile_interval,
!     physics%write_profile_flag, physics%write_pulse_flag
!   star writes: none
logical function profile_write_due(s)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      profile_write_due = .false.
      if (s%physics%profile_interval <= 0) return
      if (.not. (s%physics%write_profile_flag .or. &
           s%physics%write_pulse_flag)) return
      profile_write_due = mod(s%model_number, s%physics%profile_interval) == 0
end function profile_write_due

! ---------------------------------------------------------------
! The pulsation file of profile iprof, profile<N>.data.GYRE / .FGONG /
! .GSM, over the whole extended model (s%stitched's pulse points).
! M_star, R_star and L_star refer to the photosphere (atmosphere points
! sit above R_star, as in MESA's add_atmosphere).
!   star reads: log10_solar_radius, log_L_lsun, log_R_surface_rsun,
!     log_total_mass, physics%profile_data_prefix, solar_luminosity_cgs,
!     stitched
!   star writes: none
subroutine write_pulse(s, output, iprof, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(output_state), intent(in) :: output
      integer, intent(in) :: iprof
      integer, intent(inout) :: ierr
      character(len=:), allocatable :: path
      character(len=12) :: numstr
      real(dp) :: mstar_g, rstar_cm, lstar_cgs
      integer :: n

      mstar_g = exp(ln10*s%log_total_mass)
      rstar_cm = exp(ln10*(s%log_R_surface_rsun + s%log10_solar_radius))
      lstar_cgs = exp(ln10*s%log_L_lsun)*s%solar_luminosity_cgs

      write(numstr, '(i0)') iprof
      path = trim(output%out_dir) // trim(s%physics%profile_data_prefix) // trim(numstr) // '.data'
      n = s%stitched%n_ext
! (output_open has rejected any other format)
      select case (output%pulse_format)
      case (pulse_fgong)
         call write_fgong_pulse(s, n, s%stitched%pulse(:,1:n), mstar_g, rstar_cm, &
              lstar_cgs, path // '.FGONG', ierr)
      case (pulse_gsm)
         call write_gsm_pulse(n, s%stitched%pulse(:,1:n), mstar_g, rstar_cm, &
              lstar_cgs, path // '.GSM', ierr)
      case default
         call write_gyre_pulse(n, s%stitched%pulse(:,1:n), mstar_g, rstar_cm, &
              lstar_cgs, path // '.GYRE', ierr)
      end select
end subroutine write_pulse

! ---------------------------------------------------------------
! At the start of each run card: the card's line on the terminal, and
! the isochrone file's header.
!   star reads: card_kind, card_num_models,
!     control%isochrone_output_active, current_card,
!     current_card_envelope_x, current_card_envelope_z,
!     mixing_length_alpha, solar_mass_cgs
!   star writes: none
subroutine output_run_header(s, total_mass_msun, ierr)
      use star_info_lib, only: star_info
      use run_messages, only: terminal_unit
      use run_files, only: isochrone_unit
      use run_cards, only: card_evolve, card_rescale, card_rescale_and_evolve
      type(star_info), intent(in) :: s
      real(dp), intent(in) :: total_mass_msun
      integer, intent(inout) :: ierr
      real(dp) :: total_mass_grams
      integer :: ios
      character(len=256) :: msg

! the card's line: its kind, the rescale targets, the mixing length and
! the number of models
! (the number of models a field of 6, wider when it does not fit)
      character(len=*), parameter :: evolve_format = &
           "(/,1x,'card',i3,' (evolve):  mixing length =',f8.6,'  models =',a)", &
           rescale_format = "(/,1x,'card',i3,' (rescale):         X =',f8.6," // &
           "'  Z =',f8.6,'  mixing length =',f8.6,'  models =',a)", &
           rescale_and_evolve_format = "(/,1x,'card',i3,' (rescale+evolve):  X =',f8.6," // &
           "'  Z =',f8.6,'  mixing length =',f8.6,'  models =',a)"

      select case (s%card_kind(s%current_card))
      case (card_evolve)
         write(terminal_unit, evolve_format) s%current_card, &
              s%mixing_length_alpha, integer_field(s%card_num_models(s%current_card), 6)
      case (card_rescale)
         write(terminal_unit, rescale_format) s%current_card, s%current_card_envelope_x, &
              s%current_card_envelope_z, s%mixing_length_alpha, &
              integer_field(s%card_num_models(s%current_card), 6)
      case (card_rescale_and_evolve)
         write(terminal_unit, rescale_and_evolve_format) s%current_card, s%current_card_envelope_x, &
              s%current_card_envelope_z, s%mixing_length_alpha, &
              integer_field(s%card_num_models(s%current_card), 6)
      end select

! the isochrone file's card header: the mass (g), X, Z, the mixing
! length and the solar bolometric magnitude
      if (s%control%isochrone_output_active) then
         total_mass_grams = total_mass_msun*s%solar_mass_cgs
         write(isochrone_unit, '(7X, 1P5E16.8)', iostat=ios, iomsg=msg) total_mass_grams, &
              s%current_card_envelope_x, s%current_card_envelope_z, s%mixing_length_alpha, &
              solar_bolometric_magnitude
         call check_isochrone_write(ios, msg, ierr)
      end if
end subroutine output_run_header

! ---------------------------------------------------------------
! After a write to the isochrone file (open for the run) with
! iostat=ios and iomsg=msg: a failure reported, ierr = 1; otherwise the
! record handed to the system and the file's size on disk checked
! (output_files' flush_output_file: a full disk).
subroutine check_isochrone_write(ios, msg, ierr)
      use run_messages, only: check_io
      use run_files, only: isochrone_unit
      use output_files, only: flush_output_file
      integer, intent(in) :: ios
      character(len=*), intent(in) :: msg
      integer, intent(inout) :: ierr
      character(len=256) :: path

      inquire(unit=isochrone_unit, name=path)
      if (ios /= 0) then
         call check_io(ios, msg, 'write the isochrone file', path, ierr)
      else
         call flush_output_file(isochrone_unit, path, 0, ' ', ierr)
      end if
end subroutine check_isochrone_write

end module yrec_output
