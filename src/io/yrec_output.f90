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
! quantity first.
module yrec_output
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: ln10
      use history_output, only: history_file, history_output_init, &
           write_history_row, close_history
      use profile_output, only: profile_selection, profile_output_init, &
           write_profile
      use run_log, only: progress_lines, log_model_line, log_final_model_line
      use pulsation_output, only: write_fgong_pulse, write_gsm_pulse, &
           write_gyre_pulse, write_gyre_ext, gsm_supported
      use model_files, only: write_mod_model
      use math_lib
      implicit none
      private
      public :: output_state, output_open, output_close, output_run_header, &
           output_write_model, output_final_model_line

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
subroutine output_open(s, output, ierr)
      use star_info_lib, only: star_info
      use run_files, only: output_directory
      use io_units, only: report_error
      type(star_info), intent(in) :: s
      type(output_state), intent(inout) :: output
      integer, intent(out) :: ierr

      ierr = 0
      call close_history(output%history)
      output%out_dir = output_directory()
      output%profile_counter = 0
      output%progress = progress_lines()
      call history_output_init(s, output%history, output%out_dir, ierr)
      if (ierr /= 0) return
      call profile_output_init(s, output%profiles, ierr)
      if (ierr /= 0) return
      output%pulse_format = pulse_format_code(s%ctrl%pulse_format)
      if (s%ctrl%write_pulse_flag) then
         select case (output%pulse_format)
         case (unknown_pulse_format)
            call report_error("pulse_format = '" // trim(s%ctrl%pulse_format) // &
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
! The end of the job's output: the history file closed.
subroutine output_close(output)
      type(output_state), intent(inout) :: output
      call close_history(output%history)
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
! .mod model, and the interval GYRE file.
subroutine output_write_model(s, output, ierr)
      use star_info_lib, only: star_info, i_he4
      use io_units, only: last_model_unit
      type(star_info), intent(in) :: s
      type(output_state), intent(inout) :: output
      integer, intent(inout) :: ierr
      integer :: iprof
      real(dp) :: age_yr, luminosity_erg_s, radius_cm, teff_k, &
           gravity_cgs, ycenter_local, he_core_mass_grams
      character(len=5) :: gyre_suffix
      character(len=320) :: gyre_path

! the progress line (every terminal_interval models)
      call log_model_line(s, output%progress)

! isochrone record: model no., age (yr), L (erg/s), R (cm), Teff (K),
! g (cm/s**2), Ycenter, He-core mass (g)
      if (s%ctrl%isochrone_output_active) then
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
        write(s%ctrl%isochrone_file_unit,1005) s%model_number, &
             age_yr,luminosity_erg_s,radius_cm,teff_k,gravity_cgs, &
             ycenter_local,he_core_mass_grams
 1005   format(1X, I5, 1P7E17.8)
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
         if (s%ctrl%write_profile_flag) call write_profile(s, output%profiles, &
              output%out_dir, iprof, s%stitched, ierr)
         if (ierr /= 0) return
         if (s%ctrl%write_pulse_flag) call write_pulse(s, output, iprof, s%stitched, ierr)
         if (ierr /= 0) return
      end if
      call write_history_row(s, output%history, iprof, ierr)
      if (ierr /= 0) return
! the .mod model: the restart file, and the recovery source of a step
! that fails
      call write_mod_model(s, last_model_unit, replace=.true.)

! the interval GYRE file (independent of write_pulse_flag): numbered by
! the model (zero-padded), with the profiles' prefix and directory
      if (s%ctrl%pulse_gyre_interval > 0 .and. &
          mod(s%model_number,s%ctrl%pulse_gyre_interval) == 0) then
         write(gyre_suffix,'(I5.5)') s%model_number
         gyre_path = trim(output%out_dir)//trim(s%ctrl%profile_data_prefix)//gyre_suffix//'.data.GYRE'
         call write_gyre_pulse(s, s%nz,s%m, &
              s%logRho,s%luminosity_lsun,s%logP,s%logR, &
              s%logT,s%omega, gyre_path, ierr)
      endif
end subroutine output_write_model

! ---------------------------------------------------------------
! The end of a card: the progress line of its last converged model,
! whatever the terminal_interval phase.
subroutine output_final_model_line(s, output)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(output_state), intent(inout) :: output
      call log_final_model_line(s, output%progress)
end subroutine output_final_model_line

! ---------------------------------------------------------------
! Is a profile/pulse write due this model? A positive cadence, at
! least one of the two enabled, and the model number on the cadence.
logical function profile_write_due(s)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      profile_write_due = .false.
      if (s%ctrl%profile_interval <= 0) return
      if (.not. (s%ctrl%write_profile_flag .or. &
           s%ctrl%write_pulse_flag)) return
      profile_write_due = mod(s%model_number, s%ctrl%profile_interval) == 0
end function profile_write_due

! ---------------------------------------------------------------
! The pulsation file of profile iprof, profile<N>.data.GYRE / .FGONG /
! .GSM, over the whole extended model. M_star, R_star and L_star refer
! to the photosphere (atmosphere points sit above R_star, as in MESA's
! add_atmosphere).
subroutine write_pulse(s, output, iprof, stitched, ierr)
      use star_info_lib, only: star_info
      use stitched_model_state, only: stitched_model_record
      type(star_info), intent(in) :: s
      type(output_state), intent(in) :: output
      integer, intent(in) :: iprof
! stitched: the stitched model whose pulse columns are written.
      type(stitched_model_record), intent(in) :: stitched
      integer, intent(inout) :: ierr
      character(len=256) :: path
      character(len=12) :: numstr
      real(dp) :: mstar_g, rstar_cm, lstar_cgs

      mstar_g = exp(ln10*s%log_total_mass)
      rstar_cm = exp(ln10*(s%log_R_surface_rsun + s%log10_solar_radius))
      lstar_cgs = exp(ln10*s%log_L_lsun)*s%solar_luminosity_cgs

      write(numstr, '(i0)') iprof
      path = trim(output%out_dir) // trim(s%ctrl%profile_data_prefix) // trim(numstr) // '.data'
! (output_open has rejected any other format)
      select case (output%pulse_format)
      case (pulse_fgong)
         call write_fgong_pulse(s, stitched%n_ext, stitched%pulse(:,1:stitched%n_ext), mstar_g, rstar_cm, &
              lstar_cgs, trim(path) // '.FGONG', ierr)
      case (pulse_gsm)
         call write_gsm_pulse(stitched%n_ext, stitched%pulse(:,1:stitched%n_ext), mstar_g, rstar_cm, &
              lstar_cgs, trim(path) // '.GSM', ierr)
      case default
         call write_gyre_ext(stitched%n_ext, stitched%pulse(:,1:stitched%n_ext), mstar_g, rstar_cm, &
              lstar_cgs, trim(path) // '.GYRE', ierr)
      end select
end subroutine write_pulse

! ---------------------------------------------------------------
! At the start of each run card: the card's line on the terminal, and
! the isochrone file's header.
subroutine output_run_header(s, total_mass_msun)
      use star_info_lib, only: star_info
      use io_units, only: terminal_unit
      use run_cards, only: card_evolve, card_rescale, card_rescale_and_evolve
      type(star_info), intent(in) :: s
      real(dp), intent(in) :: total_mass_msun
      real(dp) :: total_mass_grams

      select case (s%job%rescale_kind(s%job%nk))
      case (card_evolve)
         write(terminal_unit, 47) s%job%nk, &
              s%mixing_length_alpha, s%job%num_models(s%job%nk)
      case (card_rescale)
         write(terminal_unit, 48) s%job%nk, s%job%initial_envelope_x, s%job%initial_envelope_z, &
              s%mixing_length_alpha, s%job%num_models(s%job%nk)
      case (card_rescale_and_evolve)
         write(terminal_unit, 49) s%job%nk, s%job%initial_envelope_x, s%job%initial_envelope_z, &
              s%mixing_length_alpha, s%job%num_models(s%job%nk)
      end select
  47  format(/,1x,'card',i3,' (evolve):  mixing length =',f8.6, &
             '  models =',i6)
  48  format(/,1x,'card',i3,' (rescale):         X =',f8.6, &
             '  Z =',f8.6,'  mixing length =',f8.6,'  models =',i6)
  49  format(/,1x,'card',i3,' (rescale+evolve):  X =',f8.6, &
             '  Z =',f8.6,'  mixing length =',f8.6,'  models =',i6)

      if (s%ctrl%isochrone_output_active) then
         total_mass_grams = total_mass_msun*s%solar_mass_cgs
         write(s%ctrl%isochrone_file_unit, 1495) total_mass_grams, &
              s%job%initial_envelope_x,s%job%initial_envelope_z,s%mixing_length_alpha,s%solar_bolometric_magnitude
 1495    format(7X, 1P5E16.8)
      end if
end subroutine output_run_header

end module yrec_output
