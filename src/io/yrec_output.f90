!----------------------------------------------------------------------
! yrec_output
!----------------------------------------------------------------------
! Added 2026 (MESA-style output centralization); use_legacy_output
! retired 2026-08-28 -- there is ONE output path for every run.
! Split 2026-08-29: this module is the coordinator (job init, the
! per-model dispatch, the pulse-format dispatch); the column
! machinery and the per-format tables/writers live in
!   io/output_columns.f90   parse_columns + the compiled-in defaults
!   io/history_output.f90   history column tables + row writer
!   io/profile_output.f90   profile column table + file writer
!
! Output per job:
!   <outdir>/<star_history_name>   MESA history.data layout, one row
!                                  per converged model; columns
!                                  selectable via history_columns_file
!                                  (one name per line, ! comments);
!                                  the profile_number column maps
!                                  models to profiles (YREC's
!                                  replacement for profiles.index).
!   <outdir>/profile<N>.data       every profile_interval models;
!                                  MESA profile layout (zone 1 = the
!                                  surface); columns selectable via
!                                  profile_columns_file.
!   CASE.log                       everything the code writes to
!                                  run_log_unit (diagnostics).
!   GYRE/FGONG pulse files         every pulse_gyre_interval models,
!                                  format per pulse_format.
!
! Writers are PURE READERS: every quantity is computed by
! compute_observables (core/observables_lib.f90, after the step's
! physics) into star_info.
! Module state below is per-job configuration (output directory, the
! profile counter), reset by output_init_mesa on every job --
! re-entrant runs re-enter through it, and it resets the writer
! modules' per-job state too.
module yrec_output
      use phys_const_lib
! The stitched full-star model (interior + envelope + atmosphere)
! is assembled and materialized by core/stitched_model.f90; the
! writers below are pure readers of its stx_* arrays.
      use stitched_model_lib, only: n_ext, stx_pulse, n_pulse_cols, ipul_gyre_last
      use history_output, only: history_output_init, write_history_row
      use profile_output, only: profile_output_init, write_profile
      use math_lib
      implicit none
      private
      public :: output_init_mesa, output_run_header, output_write_model

      character(len=256) :: out_dir = ' '
      integer :: profile_counter = 0

contains

! ---------------------------------------------------------------
subroutine output_init_mesa(log_output_file, ierr)
      use star_info_lib, only: star
      use luout_lib
      use run_log_lib, only: log_reset
      character(len=*), intent(in) :: log_output_file
      integer, intent(out) :: ierr
      logical :: gsm_supported
      external gsm_supported
      character(len=256) :: log_path
      integer :: n, islash

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
           status='REPLACE')
      call log_reset()

      islash = index(log_output_file, '/', back=.true.)
      if (islash > 0) then
         out_dir = log_output_file(1:islash)
      else
         out_dir = ' '
      end if
      profile_counter = 0

      call history_output_init(out_dir, ierr)
      if (ierr /= 0) return
      call profile_output_init(out_dir, ierr)
      if (ierr /= 0) return
! 2026 io-writer stops -> ierr: fail GSM-without-HDF5 at config time
! (the stub's stop at first write remains only as a last resort).
      if (star%ctrl%write_pulse_flag .and. &
          trim(star%job%pulse_format) == 'GSM') then
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
subroutine output_run_header(star_mass_msun)
      use star_info_lib, only: star
      double precision, intent(in) :: star_mass_msun

      call write_output_headers(star_mass_msun)
end subroutine output_run_header

! ---------------------------------------------------------------
subroutine output_write_model()
! 2026 use_legacy_output retirement: ONE output path for every run.
! Per converged model: the run-log progress line, the isochrone
! record (if enabled), profiles/pulse files when due, the history
! row, the .mod restart model, and the interval GYRE pulse. The
! remnants of the deleted write_legacy_output (wrtout) live here.
      use star_info_lib, only: star, i_he4
      use luout_lib
      use run_log_lib, only: log_model_line
      integer :: iprof
      double precision :: age_yr, luminosity_erg_s, radius_cm, teff_k, &
           gravity_cgs, ycenter_local, he_core_mass_grams
      character(len=5) :: gyre_suffix
      character(len=320) :: gyre_path

! the compact progress line (throttled by terminal_interval)
      call log_model_line()

! isochrone record: model no., age (yr), L (erg/s), R (cm), Teff (K),
! g (cm/s**2), Ycenter, He-core mass (g)
      if (star%ctrl%isochrone_output_active) then
        age_yr = star%dage*1.0D9
        luminosity_erg_s = exp10(star%log_L)*star%solar_luminosity_cgs
        radius_cm = exp10(star%log_R_surface)*star%solar_radius_cgs
        teff_k = exp10(star%log_Teff)
        gravity_cgs = exp10(star%log_g_surface)
        ycenter_local = star%xa(i_he4,1)
        if (star%has_h_shell) then
           he_core_mass_grams = star%m(star%h_shell_zone_begin-1)
        else
           he_core_mass_grams = 0.0D0
        end if
        write(star%ctrl%isochrone_file_unit,1005) star%model_number, &
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
      if (profile_write_due()) then
         profile_counter = profile_counter + 1
         iprof = profile_counter
         if (star%ctrl%write_profile_flag) call write_profile(iprof)
         if (star%ctrl%write_pulse_flag) call write_pulse(iprof)
      end if
      call write_history_row(iprof)
! The .mod model file: restart + in-run divergence recovery.
      call write_mod_model(last_model_unit)

! interval GYRE pulse output (independent of write_pulse_flag);
! model-numbered (zero-padded), same prefix and directory as the
! counter-numbered profile/pulse stream
      if (star%ctrl%pulse_gyre_interval.gt.0 .and. &
          mod(star%model_number,star%ctrl%pulse_gyre_interval).eq.0) then
         write(gyre_suffix,'(I5.5)') star%model_number
         gyre_path = trim(out_dir)//trim(star%ctrl%profile_data_prefix)//gyre_suffix//'.data.GYRE'
         call write_gyre_pulse(star%nz,star%model_number,star%m, &
              star%logRho,star%luminosity_lsun,star%logP,star%logR, &
              star%logT,star%omega, gyre_path)
      endif
end subroutine output_write_model

! ---------------------------------------------------------------
! Is a profile/pulse write due this model? The historical trigger:
! MESA-style output, a positive cadence, at least one product
! enabled, and the model number on the cadence.
logical function profile_write_due()
      use star_info_lib, only: star
      profile_write_due = .false.
      if (star%ctrl%profile_interval <= 0) return
      if (.not. (star%ctrl%write_profile_flag .or. &
           star%ctrl%write_pulse_flag)) return
      profile_write_due = mod(star%model_number, star%ctrl%profile_interval) == 0
end function profile_write_due

! Pulse file alongside profile <iprof>: profile<N>.data.GYRE /
! .data.FGONG / .data.GSM in the output directory, covering the FULL
! extended model. Global M_star/R_star/L_star refer to the
! photosphere (atmosphere points extend above R_star, as in MESA's
! add_atmosphere).
subroutine write_pulse(iprof)
      use math_lib
      use star_info_lib, only: star
      integer, intent(in) :: iprof
      character(len=256) :: path
      character(len=12) :: numstr
      double precision :: mstar_g, rstar_cm, lstar_cgs

      mstar_g = exp(ln10*star%log_total_mass)
      rstar_cm = exp(ln10*(star%log_R_surface + star%log10_solar_radius))
      lstar_cgs = exp(ln10*star%log_L)*star%solar_luminosity_cgs

      write(numstr, '(i0)') iprof
      if (star%job%pulse_format(1:5) == 'FGONG' .or. &
          star%job%pulse_format(1:5) == 'fgong') then
         path = trim(out_dir) // trim(star%ctrl%profile_data_prefix) // trim(numstr) // '.data.FGONG'
         call write_fgong_pulse(n_ext, stx_pulse(:,1:n_ext), mstar_g, rstar_cm, &
              lstar_cgs, path)
      else if (star%job%pulse_format(1:3) == 'GSM' .or. &
               star%job%pulse_format(1:3) == 'gsm') then
         path = trim(out_dir) // trim(star%ctrl%profile_data_prefix) // trim(numstr) // '.data.GSM'
         call write_gsm_pulse(n_ext, stx_pulse(:,1:n_ext), mstar_g, rstar_cm, &
              lstar_cgs, path)
      else
         path = trim(out_dir) // trim(star%ctrl%profile_data_prefix) // trim(numstr) // '.data.GYRE'
         call write_gyre_ext(n_ext, stx_pulse(:,1:n_ext), mstar_g, rstar_cm, &
              lstar_cgs, path)
      end if
end subroutine write_pulse

! GYRE text (MESA schema 101) over the extended point set. The
! legacy-path writer io/write_gyre_pulse.f90 keeps its historical
! interior-only behavior (byte-pinned); this is the MESA-mode writer.
subroutine write_gyre_ext(n, pts, mstar_g, rstar_cm, lstar_cgs, path)
      use star_info_lib, only: star
      integer, intent(in) :: n
      double precision, intent(in) :: pts(n_pulse_cols, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: path
      integer, parameter :: gyre_schema = 101
      integer :: u, j, i

      open(newunit=u, file=path, status='REPLACE', form='FORMATTED')
      write(u,'(I6,3(1X,1PE26.16),1X,I6)') n, mstar_g, rstar_cm, &
           lstar_cgs, gyre_schema
      do j = 1, n
         write(u,'(I6,99(1X,1PE26.16))') j, (pts(i,j), i = 1, ipul_gyre_last)
      end do
      close(u)
end subroutine write_gyre_ext

end module yrec_output
