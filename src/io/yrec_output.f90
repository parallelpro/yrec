!----------------------------------------------------------------------
! yrec_output
!----------------------------------------------------------------------
! Added 2026 (MESA-style output centralization); use_legacy_output
! retired 2026-08-28 -- there is ONE output path for every run.
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
! Module state below is per-job configuration (paths, parsed column
! selections, the profile counter), reset by output_init_mesa on
! every job -- re-entrant runs re-enter through it.
module yrec_output
      use phys_const_lib
! The stitched full-star model (interior + envelope + atmosphere)
! is assembled and materialized by core/stitched_model.f90; the
! writers below are pure readers of its stx_* arrays.
      use stitched_model_lib, only: n_ext, &
           stx_prof, stx_pulse, n_prof_cols, n_pulse_cols
      implicit none
      private
      public :: output_init_mesa, output_run_header, output_write_model

      integer, parameter :: max_cols = 128
      integer, parameter :: n_hist_cols = 86
! the default column selections (blank history_columns_file /
! profile_columns_file) are compiled in from
! defaults/{history,profile}_columns.list: uncommented names, file
! order. Opt-in columns (2026: the seismic pair) ship commented out
! there; uncomment and rebuild to write them by default.
      include 'default_columns.inc'

      character(len=256) :: out_dir = ' '
      character(len=256) :: hist_path = ' '
      integer :: profile_counter = 0
      integer :: hist_nsel = 0, prof_nsel = 0
      integer :: hist_sel(max_cols), prof_sel(max_cols)


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
      character(len=24) :: hist_names(n_hist_cols), prof_names(n_prof_cols)
      integer :: n, islash, prior_unit
      logical :: prior_open

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
      hist_path = trim(out_dir) // trim(star%ctrl%star_history_name)
! in-process re-entry: a history unit still connected from a
! previous job would silently accumulate; close it so the lazy
! open in write_history_row starts the file over.
      inquire(file=hist_path, opened=prior_open, number=prior_unit)
      if (prior_open) close(prior_unit)
      profile_counter = 0

      ierr = 0
      call history_column_names(hist_names)
      call parse_columns(star%ctrl%history_columns_file, hist_names, n_hist_cols, &
           hist_default_names, n_hist_default, hist_sel, hist_nsel, 'history', ierr)
      if (ierr /= 0) return
! model_number / profile_number / num_zones lead the file whenever
! they are selected, in that fixed order (they are columns 1-3 of the
! built-in table, so a stable hoist of indices 1..3 does it).
      call hoist_id_columns(hist_sel, hist_nsel)
      call profile_column_names(prof_names)
      call parse_columns(star%ctrl%profile_columns_file, prof_names, n_prof_cols, &
           prof_default_names, n_prof_default, prof_sel, prof_nsel, 'profile', ierr)
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
        luminosity_erg_s = 10.0D0**star%log_L*star%solar_luminosity_cgs
        radius_cm = 10.0D0**star%log_R_surface*star%solar_radius_cgs
        teff_k = 10.0D0**star%log_Teff
        gravity_cgs = 10.0D0**star%log_g_surface
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
! Column-selection files: one column name per line, '!' comments,
! blank lines ignored. Blank/absent control -> all columns in the
! built-in order. Unknown names are fatal (config error), with the
! valid names listed in the log.
subroutine parse_columns(fname, names, ncol, default_names, ndefault, &
           sel, nsel, label, ierr)
      use luout_lib
      character(len=*), intent(in) :: fname, label
      integer, intent(out) :: ierr
      integer, intent(in) :: ncol, ndefault
      character(len=24), intent(in) :: names(ncol), default_names(ndefault)
      integer, intent(out) :: sel(max_cols), nsel
      character(len=256) :: line
      integer :: u, ios, i

      ierr = 0
      nsel = 0
      if (len_trim(fname) == 0) then
! blank columns file: the compiled-in default selection (generated
! from defaults/<label>_columns.list -- see default_columns.inc)
         do i = 1, ndefault
            call append_column(default_names(i), names, ncol, sel, nsel, &
                 label, ierr)
            if (ierr /= 0) return
         end do
         return
      end if
      open(newunit=u, file=fname, status='OLD', action='READ', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'cannot open ', trim(label), '_columns_file: ', &
              trim(fname)
         write(run_log_unit,*) 'cannot open ', trim(label), &
              '_columns_file: ', trim(fname)
! 2026 io-writer stops -> ierr (stage-3 pattern): config error
! returned to output_init_mesa -> read_input -> read_controls.
         ierr = 1
         return
      end if
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         i = index(line, '!')
         if (i > 0) line = line(1:i-1)
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         call append_column(line, names, ncol, sel, nsel, label, ierr)
         if (ierr /= 0) return
      end do
      close(u)
      if (nsel == 0) then
         write(*,*) trim(label), '_columns_file selected no columns: ', &
              trim(fname)
         ierr = 1
         return
      end if
end subroutine parse_columns

! ---------------------------------------------------------------
! Look a column name up in the writer's table and append its index
! to the selection; unknown names are fatal, with the valid list
! written to the run log.
subroutine append_column(name, names, ncol, sel, nsel, label, ierr)
      use luout_lib
      character(len=*), intent(in) :: name, label
      integer, intent(in) :: ncol
      character(len=24), intent(in) :: names(ncol)
      integer, intent(inout) :: sel(max_cols), nsel
      integer, intent(out) :: ierr
      integer :: j

      ierr = 0
      do j = 1, ncol
         if (trim(name) == trim(names(j))) then
            nsel = nsel + 1
            sel(nsel) = j
            return
         end if
      end do
      write(*,*) 'unknown ', trim(label), ' column: ', trim(name)
      write(run_log_unit,*) 'unknown ', trim(label), &
           ' column: ', trim(name)
      write(run_log_unit,*) 'valid ', trim(label), ' columns:'
      do j = 1, ncol
         write(run_log_unit,'(2x,a)') trim(names(j))
      end do
      ierr = 1
end subroutine append_column

! ---------------------------------------------------------------
subroutine history_column_names(names)
      character(len=24), intent(out) :: names(n_hist_cols)

      names(1) = 'model_number'
      names(2) = 'profile_number'
      names(3) = 'num_zones'
      names(4) = 'log_L'
      names(5) = 'log_R'
      names(6) = 'log_g'
      names(7) = 'log_Teff'
      names(8) = 'mass_conv_core'
      names(9) = 'cz_bot_mass'
      names(10) = 'cz_bot_radius'
      names(11) = 'cz_bot_temperature'
      names(12) = 'cz_bot_density'
      names(13) = 'cz_bot_pressure'
      names(14) = 'cz_bot_opacity'
      names(15) = 'log_center_T'
      names(16) = 'log_center_Rho'
      names(17) = 'log_center_P'
      names(18) = 'center_beta'
      names(19) = 'center_degeneracy'
      names(20) = 'center_h1'
      names(21) = 'center_he4'
      names(22) = 'center_z'
      names(23) = 'lum_ppI'
      names(24) = 'lum_ppII'
      names(25) = 'lum_ppIII'
      names(26) = 'lum_cno'
      names(27) = 'lum_3alpha'
      names(28) = 'lum_heC'
      names(29) = 'lum_grav'
      names(30) = 'lum_neu'
      names(31) = 'snu_cl37'
      names(32) = 'snu_ga71'
      names(33) = 'neut_pp'
      names(34) = 'neut_pep'
      names(35) = 'neut_hep'
      names(36) = 'neut_be7'
      names(37) = 'neut_b8'
      names(38) = 'neut_n13'
      names(39) = 'neut_o15'
      names(40) = 'neut_f17'
      names(41) = 'neut_diag_1'
      names(42) = 'neut_diag_2'
      names(43) = 'center_he3'
      names(44) = 'center_c12'
      names(45) = 'center_c13'
      names(46) = 'center_n14'
      names(47) = 'center_n15'
      names(48) = 'center_o16'
      names(49) = 'center_o17'
      names(50) = 'center_o18'
      names(51) = 'surface_he3'
      names(52) = 'surface_c12'
      names(53) = 'surface_c13'
      names(54) = 'surface_n14'
      names(55) = 'surface_n15'
      names(56) = 'surface_o16'
      names(57) = 'surface_o17'
      names(58) = 'surface_o18'
      names(59) = 'surface_h2'
      names(60) = 'surface_li6'
      names(61) = 'surface_li7'
      names(62) = 'surface_be9'
      names(63) = 'surface_h1'
      names(64) = 'surface_he4'
      names(65) = 'surface_z'
      names(66) = 'surface_z_div_x'
      names(67) = 'total_angular_momentum'
      names(68) = 'total_rotational_ke'
      names(69) = 'total_moment_of_inertia'
      names(70) = 'cz_moment_of_inertia'
      names(71) = 'surf_avg_omega'
      names(72) = 'center_omega'
      names(73) = 'rotation_period_days'
      names(74) = 'surf_avg_v_rot'
      names(75) = 'tau_conv_cz'
      names(76) = 'h_shell_bot_mass'
      names(77) = 'h_shell_mid_mass'
      names(78) = 'h_shell_top_mass'
      names(79) = 'h_shell_bot_radius'
      names(80) = 'h_shell_mid_radius'
      names(81) = 'h_shell_top_radius'
      names(82) = 'log_P_photosphere'
      names(83) = 'star_mass'
      names(84) = 'star_age'
! opt-in columns (commented out in defaults/history_columns.list)
      names(85) = 'nu_max'
      names(86) = 'delta_nu'
end subroutine history_column_names

subroutine history_values(vals, iprof)
      use star_info_lib, only: star, i_h1, i_he4, i_lum_grav, i_lum_he_c, i_lum_neu, i_metals
      double precision, intent(out) :: vals(n_hist_cols)
      integer, intent(in) :: iprof
      integer :: i

      vals(1) = dble(star%model_number)
      vals(2) = dble(iprof)
      vals(3) = dble(star%nz)
      vals(4) = star%log_L
      vals(5) = star%log_R_surface
      vals(6) = star%log_g_surface
      vals(7) = star%log_Teff
      vals(8) = star%core_cz_mass
      vals(9) = star%envelope_mass
      vals(10) = star%envelope_radius
      vals(11) = star%envelope_cz_temperature
      vals(12) = star%envelope_cz_density
      vals(13) = star%envelope_cz_pressure
      vals(14) = star%envelope_cz_opacity
      vals(15) = star%central_log10_temperature
      vals(16) = star%central_log10_density
      vals(17) = star%central_log10_pressure
      vals(18) = star%central_beta
      vals(19) = star%central_degeneracy_eta
      vals(20) = star%xa(i_h1,1)
      vals(21) = star%xa(i_he4,1)
      vals(22) = star%xa(i_metals,1)
      do i = 1, 5
         vals(22+i) = star%luminosity_breakdown(i)
      end do
      vals(28) = star%luminosity_breakdown(i_lum_he_c)
      vals(29) = star%luminosity_breakdown(i_lum_grav)
      vals(30) = star%luminosity_breakdown(i_lum_neu)
      vals(31) = star%cl37_snu_rate
      vals(32) = star%ga71_snu_rate
      do i = 1, 10
         vals(32+i) = star%neutrino_flux_total(i)
      end do
      do i = 4, 11
         vals(39+i) = star%xa(i,1)
      end do
      do i = 4, 15
         vals(47+i) = star%xa(i,star%nz)
      end do
      vals(63) = star%xa(i_h1,star%nz)
      vals(64) = star%xa(i_he4,star%nz)
      vals(65) = star%xa(i_metals,star%nz)
      vals(66) = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
      vals(67) = star%total_angular_momentum
      vals(68) = star%total_rotational_ke
      vals(69) = star%total_moment_of_inertia
      vals(70) = star%cz_moment_of_inertia
      vals(71) = star%omega(star%nz)
      vals(72) = star%omega(1)
      vals(73) = star%rotation_period_days
      vals(74) = star%surf_velocity_kms
      vals(75) = star%convective_turnover_timescale
      vals(76) = star%h_shell_bot_mass
      vals(77) = star%h_shell_mid_mass
      vals(78) = star%h_shell_top_mass
      vals(79) = star%h_shell_bot_radius
      vals(80) = star%h_shell_mid_radius
      vals(81) = star%h_shell_top_radius
      vals(82) = star%pphot
      vals(83) = star%star_mass
      vals(84) = star%dage*1.0d9
      vals(85) = star%nu_max
      vals(86) = star%delta_nu
end subroutine history_values

subroutine write_history_row(iprof)
      use star_info_lib, only: star, i_metals
      integer, intent(in) :: iprof
      character(len=24) :: names(n_hist_cols)
      double precision :: vals(n_hist_cols)
      integer :: hist_unit, i, k
      logical :: hist_open

      call history_values(vals, iprof)
      inquire(file=hist_path, opened=hist_open, number=hist_unit)
      if (.not. hist_open) then
         open(newunit=hist_unit, file=hist_path, status='REPLACE', &
              action='WRITE')
         call history_column_names(names)
! Global block: run metadata, matching what the legacy .track header
! carried (initial composition + mixing length). The initial_* values
! are the kind card's starting envelope composition (initial_y by
! closure), not the current surface values.
         write(hist_unit, '(6(1x,i40))') 1, 2, 3, 4, 5, 6
         write(hist_unit, '(6(1x,a40))') adjustr('version_number'), &
              adjustr('initial_mass'), adjustr('initial_x'), &
              adjustr('initial_y'), adjustr('initial_z'), &
              adjustr('mixing_length_alpha')
         write(hist_unit, '(1x,a40,5(1x,es40.16e3))') &
              adjustr('"' // trim(yrec_version_string) // '"'), &
              star%star_mass, star%job%initial_envelope_x, &
              1.0d0 - star%job%initial_envelope_x - star%job%initial_envelope_z, &
              star%job%initial_envelope_z, star%mixing_length_alpha
         write(hist_unit, '(a)') ''
         write(hist_unit, '(999(1x,i40))') (k, k = 1, hist_nsel)
         write(hist_unit, '(999(1x,a40))') &
              (adjustr(trim(names(hist_sel(k)))), k = 1, hist_nsel)
      end if
      do i = 1, hist_nsel
         if (is_int_hist_col(hist_sel(i))) then
            write(hist_unit, '(1x,i40)', advance='no') &
                 nint(vals(hist_sel(i)))
         else
            write(hist_unit, '(1x,es40.16e3)', advance='no') &
                 vals(hist_sel(i))
         end if
      end do
      write(hist_unit, '(a)') ''
! Flush so the file is live during the run (tail -f / editor view),
! matching MESA's behavior.
      flush(hist_unit)
end subroutine write_history_row

! ---------------------------------------------------------------
subroutine profile_column_names(names)
      character(len=24), intent(out) :: names(n_prof_cols)

      names(1) = 'zone'
      names(2) = 'mass'
      names(3) = 'logR'
      names(4) = 'logT'
      names(5) = 'logRho'
      names(6) = 'logP'
      names(7) = 'luminosity'
      names(8) = 'dm'
      names(9) = 'mixing_type'
      names(10) = 'gamma1'
      names(11) = 'opacity'
      names(12) = 'gradr'
      names(13) = 'gradT'
      names(14) = 'grada'
      names(15) = 'conv_vel'
      names(16) = 'beta'
      names(17) = 'eta'
      names(18) = 'grav'
      names(19) = 'eps_ppI'
      names(20) = 'eps_ppII'
      names(21) = 'eps_ppIII'
      names(22) = 'eps_cno'
      names(23) = 'eps_he3'
      names(24) = 'eps_nuc'
      names(25) = 'eps_neu'
      names(26) = 'eps_grav'
      names(27) = 'h1'
      names(28) = 'h2'
      names(29) = 'he3'
      names(30) = 'he4'
      names(31) = 'li6'
      names(32) = 'li7'
      names(33) = 'be9'
      names(34) = 'c12'
      names(35) = 'c13'
      names(36) = 'n14'
      names(37) = 'n15'
      names(38) = 'o16'
      names(39) = 'o17'
      names(40) = 'o18'
      names(41) = 'z'
      names(42) = 'omega'
      names(43) = 'j_rot'
      names(44) = 'i_rot'
      names(45) = 'fp_rot'
      names(46) = 'ft_rot'
      names(47) = 'v_es'
      names(48) = 'v_gsf'
      names(49) = 'v_ss'
      names(50) = 'cp'
      names(51) = 'delta'
      names(52) = 'mu'
      names(53) = 'mu_e_inv'
      names(54) = 'brunt_N2'
      names(55) = 'lamb_S2'
      names(56) = 'gradL'
      names(57) = 'gradr_div_grada'
      names(58) = 'csound'
      names(59) = 'D_omega'
      names(60) = 'D_mix_rot'
end subroutine profile_column_names

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

! ---------------------------------------------------------------
subroutine write_profile(iprof)
      use star_info_lib, only: star
      integer, intent(in) :: iprof
      character(len=24) :: names(n_prof_cols)
      character(len=256) :: path
      character(len=12) :: numstr
      integer :: u, i, k, kz
      double precision :: v

      call profile_column_names(names)
      write(numstr, '(i0)') iprof
      path = trim(out_dir) // trim(star%ctrl%profile_data_prefix) // trim(numstr) // '.data'
      open(newunit=u, file=path, status='REPLACE', action='WRITE')
! global block (MESA profile shape: numbers / names / values)
      write(u, '(8(1x,i40))') 1, 2, 3, 4, 5, 6, 7, 8
      write(u, '(8(1x,a40))') adjustr('model_number'), &
           adjustr('num_zones'), adjustr('initial_mass'), &
           adjustr('initial_z'), adjustr('star_age'), &
           adjustr('time_step'), adjustr('log_Teff'), &
           adjustr('star_mass')
! initial_mass/initial_z are the kind card's starting values
! (pulsation_mass_msun is stamped by begin_kind_card); star_age and
! time_step are in years, matching MESA's profile header.
      write(u, '(2(1x,i40),6(1x,es40.16e3))') star%model_number, &
           n_ext, star%pulsation_mass_msun, &
           star%job%initial_envelope_z, star%dage*1.0d9, &
           star%timestep_yr, star%log_Teff, star%star_mass
      write(u, '(a)') ''
      write(u, '(999(1x,i40))') (k, k = 1, prof_nsel)
      write(u, '(999(1x,a40))') &
           (adjustr(trim(names(prof_sel(k)))), k = 1, prof_nsel)
! rows: zone 1 = the outermost point (MESA convention). The extended
! model (interior + envelope + atmosphere, from stitched_model_lib)
! is stored inward-to-outward, so walk it in reverse.
      do k = 1, n_ext
         kz = n_ext - k + 1
         do i = 1, prof_nsel
            if (prof_sel(i) == 1) then
               v = dble(k)
            else
               v = stx_prof(prof_sel(i), kz)
            end if
            if (prof_sel(i) == 1 .or. prof_sel(i) == 9) then
               write(u, '(1x,i40)', advance='no') nint(v)
            else
               write(u, '(1x,es40.16e3)', advance='no') v
            end if
         end do
         write(u, '(a)') ''
      end do
      close(u)
end subroutine write_profile

! Pulse file alongside profile <iprof>: profile<N>.data.GYRE /
! .data.FGONG / .data.GSM in the output directory, covering the FULL
! extended model. Global M_star/R_star/L_star refer to the
! photosphere (atmosphere points extend above R_star, as in MESA's
! add_atmosphere).
subroutine write_pulse(iprof)
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
      double precision, intent(in) :: pts(35, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: path
      integer, parameter :: gyre_schema = 101
      integer :: u, j, i

      open(newunit=u, file=path, status='REPLACE', form='FORMATTED')
      write(u,'(I6,3(1X,1PE26.16),1X,I6)') n, mstar_g, rstar_cm, &
           lstar_cgs, gyre_schema
      do j = 1, n
         write(u,'(I6,99(1X,1PE26.16))') j, (pts(i,j), i = 1, 18)
      end do
      close(u)
end subroutine write_gyre_ext

! History columns written as true integers.
logical function is_int_hist_col(icol)
      integer, intent(in) :: icol
      is_int_hist_col = (icol >= 1 .and. icol <= 3)
end function is_int_hist_col

! Stable hoist of the id columns (built-in indices 1..3, i.e.
! model_number/profile_number/num_zones) to the front of a selection.
subroutine hoist_id_columns(sel, nsel)
      integer, intent(inout) :: sel(max_cols), nsel
      integer :: tmp(max_cols), n, id, i

      n = 0
      do id = 1, 3
         do i = 1, nsel
            if (sel(i) == id) then
               n = n + 1
               tmp(n) = id
            end if
         end do
      end do
      do i = 1, nsel
         if (sel(i) > 3) then
            n = n + 1
            tmp(n) = sel(i)
         end if
      end do
      sel(1:n) = tmp(1:n)
      nsel = n
end subroutine hoist_id_columns

end module yrec_output
