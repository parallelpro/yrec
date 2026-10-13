!----------------------------------------------------------------------
! yrec_output
!----------------------------------------------------------------------
! The run's output, per job:
!   <outdir>/<star_history_name>  one row per converged model, MESA's
!                                 history.data layout (history_output);
!                                 the profile_number column maps models
!                                 to profiles.
!   <outdir>/profile<N>.data      every profile_interval models, MESA's
!                                 profile layout, zone 1 the surface
!                                 (profile_output), with the pulsation
!                                 file in pulse_format beside it
!                                 (pulsation_output).
!   the run log                   everything written to run_log_unit.
!   the .mod model                every model (model_files).
! Columns are selectable through history_columns_file and
! profile_columns_file (output_columns). This module opens the job's
! files, writes the run card's header, and decides each model what is
! due. The writers only read: compute_observables (core) and the
! stitched model (core/stitched_model) compute every quantity first.
! The output directory and the profile counter are per job, reset by
! output_init_mesa together with the writers' own state.
module yrec_output
      use physical_constants
      use history_output, only: history_output_init, write_history_row
      use profile_output, only: profile_output_init, write_profile
      use pulsation_output, only: write_fgong_pulse, write_gsm_pulse, write_gyre_pulse, &
           write_gyre_ext, gsm_supported
      use model_files, only: write_mod_model
      use math_lib
      implicit none
      private
      public :: output_init_mesa, output_run_header, output_write_model

      character(len=256) :: out_dir = ' '
      integer :: profile_counter = 0

contains

! ---------------------------------------------------------------
subroutine output_init_mesa(s, log_output_file, ierr)
      use star_info_lib, only: star_info
      use io_units
      use run_log, only: log_reset
      type(star_info), intent(inout) :: s
      character(len=*), intent(in) :: log_output_file
      integer, intent(out) :: ierr
      character(len=256) :: log_path, msg
      integer :: n, islash, ios

      ierr = 0
! run-log name: keep an explicit .log name as-is, swap a legacy
! .short name to .log, append .log otherwise (so the default
! run.log stays run.log, not run.log.log).
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
      if (ierr /= 0) return
      call log_reset()

      islash = index(log_output_file, '/', back=.true.)
      if (islash > 0) then
         out_dir = log_output_file(1:islash)
      else
         out_dir = ' '
      end if
      profile_counter = 0

      call history_output_init(s, out_dir, ierr)
      if (ierr /= 0) return
      call profile_output_init(s, out_dir, ierr)
      if (ierr /= 0) return
! 2026 io-writer stops -> ierr: fail GSM-without-HDF5 at config time
! (the stub's stop at first write remains only as a last resort).
      if (s%ctrl%write_pulse_flag .and. &
          trim(s%job%pulse_format) == 'GSM') then
         if (.not. gsm_supported()) then
            write(*,*) 'pulse_format = GSM requires an HDF5-enabled build:'
            write(*,*) '  make clean && make USE_HDF5=1'
            write(run_log_unit,*) &
                 'pulse_format = GSM requires USE_HDF5=1 build'
            ierr = 1
            return
         end if
      end if
end subroutine output_init_mesa

! ---------------------------------------------------------------
subroutine output_write_model(s, ierr)
      use star_info_lib, only: star_info
! 2026 use_legacy_output retirement: ONE output path for every run.
! Per converged model: the run-log progress line, the isochrone
! record (if enabled), profiles/pulse files when due, the history
! row, the .mod restart model, and the interval GYRE pulse. The
! remnants of the deleted write_legacy_output (wrtout) live here.
      use star_info_lib, only: i_he4
      use io_units
      use run_log, only: log_model_line
      type(star_info), intent(inout) :: s
      integer, intent(inout) :: ierr
      integer :: iprof
      double precision :: age_yr, luminosity_erg_s, radius_cm, teff_k, &
           gravity_cgs, ycenter_local, he_core_mass_grams
      character(len=5) :: gyre_suffix
      character(len=320) :: gyre_path

! the compact progress line (throttled by terminal_interval)
      call log_model_line(s)

! isochrone record: model no., age (yr), L (erg/s), R (cm), Teff (K),
! g (cm/s**2), Ycenter, He-core mass (g)
      if (s%ctrl%isochrone_output_active) then
        age_yr = s%age_gyr*1.0D9
        luminosity_erg_s = exp10(s%log_L)*s%solar_luminosity_cgs
        radius_cm = exp10(s%log_R_surface)*s%solar_radius_cgs
        teff_k = exp10(s%log_Teff)
        gravity_cgs = exp10(s%log_g_surface)
        ycenter_local = s%xa(i_he4,1)
        if (s%has_h_shell) then
           he_core_mass_grams = s%m(s%h_shell_zone_begin-1)
        else
           he_core_mass_grams = 0.0D0
        end if
        write(s%ctrl%isochrone_file_unit,1005) s%model_number, &
             age_yr,luminosity_erg_s,radius_cm,teff_k,gravity_cgs, &
             ycenter_local,he_core_mass_grams
 1005   format(1X, I5, 1P7E17.8)
      end if

! Profiles and pulse files share the trigger and the number (MESA's
! write_pulse_data_with_profile coupling): every profile_interval
! models, the counter advances and each enabled product is written --
! profile<N>.data, and/or profile<N>.data.GYRE / .data.FGONG per
! pulse_format. The history profile_number column records N. The
! stitched arrays are built every step by evolve_step (ahead of
! compute_observables); this block only decides whether THIS model
! gets profile/pulse files, then numbers and writes.
      iprof = 0
      if (profile_write_due(s)) then
         profile_counter = profile_counter + 1
         iprof = profile_counter
         if (s%ctrl%write_profile_flag) call write_profile(s, iprof, s%stitched, ierr)
         if (ierr /= 0) return
         if (s%ctrl%write_pulse_flag) call write_pulse(s, iprof, s%stitched, ierr)
         if (ierr /= 0) return
      end if
      call write_history_row(s, iprof, ierr)
      if (ierr /= 0) return
! The .mod model file: restart + in-run divergence recovery.
      call write_mod_model(s, last_model_unit)

! interval GYRE pulse output (independent of write_pulse_flag);
! model-numbered (zero-padded), same prefix and directory as the
! counter-numbered profile/pulse stream
      if (s%ctrl%pulse_gyre_interval.gt.0 .and. &
          mod(s%model_number,s%ctrl%pulse_gyre_interval).eq.0) then
         write(gyre_suffix,'(I5.5)') s%model_number
         gyre_path = trim(out_dir)//trim(s%ctrl%profile_data_prefix)//gyre_suffix//'.data.GYRE'
         call write_gyre_pulse(s, s%nz,s%m, &
              s%logRho,s%luminosity_lsun,s%logP,s%logR, &
              s%logT,s%omega, gyre_path, ierr)
      endif
end subroutine output_write_model

! ---------------------------------------------------------------
! Is a profile/pulse write due this model? The historical trigger:
! MESA-style output, a positive cadence, at least one product
! enabled, and the model number on the cadence.
logical function profile_write_due(s)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
      profile_write_due = .false.
      if (s%ctrl%profile_interval <= 0) return
      if (.not. (s%ctrl%write_profile_flag .or. &
           s%ctrl%write_pulse_flag)) return
      profile_write_due = mod(s%model_number, s%ctrl%profile_interval) == 0
end function profile_write_due

! Pulse file alongside profile <iprof>: profile<N>.data.GYRE /
! .data.FGONG / .data.GSM in the output directory, covering the FULL
! extended model. Global M_star/R_star/L_star refer to the
! photosphere (atmosphere points extend above R_star, as in MESA's
! add_atmosphere).
subroutine write_pulse(s, iprof, stitched, ierr)
      use star_info_lib, only: star_info
      use math_lib
      use stitched_model_state, only: stitched_model_record
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iprof
! stitched: the stitched model whose pulse columns are written.
      type(stitched_model_record), intent(in) :: stitched
      integer, intent(inout) :: ierr
      character(len=256) :: path
      character(len=12) :: numstr
      double precision :: mstar_g, rstar_cm, lstar_cgs

      mstar_g = exp(ln10*s%log_total_mass)
      rstar_cm = exp(ln10*(s%log_R_surface + s%log10_solar_radius))
      lstar_cgs = exp(ln10*s%log_L)*s%solar_luminosity_cgs

      write(numstr, '(i0)') iprof
      if (s%job%pulse_format(1:5) == 'FGONG' .or. &
          s%job%pulse_format(1:5) == 'fgong') then
         path = trim(out_dir) // trim(s%ctrl%profile_data_prefix) // trim(numstr) // '.data.FGONG'
         call write_fgong_pulse(s, stitched%n_ext, stitched%pulse(:,1:stitched%n_ext), mstar_g, rstar_cm, &
              lstar_cgs, path, ierr)
      else if (s%job%pulse_format(1:3) == 'GSM' .or. &
               s%job%pulse_format(1:3) == 'gsm') then
         path = trim(out_dir) // trim(s%ctrl%profile_data_prefix) // trim(numstr) // '.data.GSM'
         call write_gsm_pulse(stitched%n_ext, stitched%pulse(:,1:stitched%n_ext), mstar_g, rstar_cm, &
              lstar_cgs, path, ierr)
      else
         path = trim(out_dir) // trim(s%ctrl%profile_data_prefix) // trim(numstr) // '.data.GYRE'
         call write_gyre_ext(stitched%n_ext, stitched%pulse(:,1:stitched%n_ext), mstar_g, rstar_cm, &
              lstar_cgs, path, ierr)
      end if
end subroutine write_pulse

! ---------------------------------------------------------------
! At the start of each run card: the card's line on the terminal, and
! the isochrone file's header.
subroutine output_run_header(s, total_mass_msun)
      use star_info_lib, only: star_info
      use io_units
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: total_mass_msun
      double precision :: total_mass_grams

      if (s%job%rescale_kind(s%job%nk) .eq. 1) then
         write(terminal_unit, 47) s%job%nk, &
              s%mixing_length_alpha, s%job%num_models(s%job%nk)
      else if (s%job%rescale_kind(s%job%nk) .eq. 2) then
         write(terminal_unit, 48) s%job%nk, s%job%initial_envelope_x, s%job%initial_envelope_z, &
              s%mixing_length_alpha, s%job%num_models(s%job%nk)
      else if (s%job%rescale_kind(s%job%nk) .eq. 3) then
         write(terminal_unit, 49) s%job%nk, s%job%initial_envelope_x, s%job%initial_envelope_z, &
              s%mixing_length_alpha, s%job%num_models(s%job%nk)
      end if
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

      return
end subroutine output_run_header

end module yrec_output
