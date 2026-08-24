!----------------------------------------------------------------------
! yrec_output
!----------------------------------------------------------------------
! Added 2026 (MESA-style output centralization). The output-mode
! branch (controls_lib's use_legacy_output) has exactly three
! decision points: parmin's legacy open-block (wrapped there; its
! MESA else-branch calls output_init_mesa), output_run_header, and
! output_write_model. run_yrec's final-model store and calibration
! chain are legacy-only blocks wrapped in place.
!
! MESA-mode output per job:
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
!                                  short_file_unit (diagnostics).
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
      use const_lib
      implicit none
      private
      public :: output_init_mesa, output_run_header, output_write_model

      integer, parameter :: max_cols = 128
      integer, parameter :: n_hist_cols = 84
      integer, parameter :: n_prof_cols = 53

      character(len=256) :: out_dir = ' '
      character(len=256) :: hist_path = ' '
      integer :: profile_counter = 0
      integer :: hist_nsel = 0, prof_nsel = 0
      integer :: hist_sel(max_cols), prof_sel(max_cols)

! The extended model: interior (center -> fitting point) + envelope
! (fitting point -> photosphere) + atmosphere (photosphere -> tau~0),
! assembled inward-to-outward, exactly the regions io/stitch.f90
! splices for the legacy .store format. Profiles and pulse files both
! cover the full star: truncating at the fitting point would drop the
! superadiabatic layer and photosphere, which dominate p-mode
! frequencies.
      integer, parameter :: max_ext = 3*5000
      integer :: n_ext = 0
      integer :: ext_region(max_ext)   ! 1 interior, 2 envelope, 3 atmosphere
      integer :: ext_index(max_ext)    ! index within that region

contains

! ---------------------------------------------------------------
subroutine output_init_mesa(fshort)
      use luout_lib
      character(len=*), intent(in) :: fshort
      character(len=256) :: log_path
      character(len=24) :: hist_names(n_hist_cols), prof_names(n_prof_cols)
      integer :: n, islash

      n = index(fshort, '.short')
      if (n > 0) then
         log_path = fshort(1:n) // 'log'
      else
         log_path = trim(fshort) // '.log'
      end if
      open(short_file_unit, file=log_path, form='FORMATTED', &
           status='REPLACE')

      islash = index(fshort, '/', back=.true.)
      if (islash > 0) then
         out_dir = fshort(1:islash)
      else
         out_dir = ' '
      end if
      hist_path = trim(out_dir) // trim(star_history_name)
      profile_counter = 0

      call history_column_names(hist_names)
      call parse_columns(history_columns_file, hist_names, n_hist_cols, &
           hist_sel, hist_nsel, 'history')
! model_number / profile_number / num_zones lead the file whenever
! they are selected, in that fixed order (they are columns 1-3 of the
! built-in table, so a stable hoist of indices 1..3 does it).
      call hoist_id_columns(hist_sel, hist_nsel)
      call profile_column_names(prof_names)
      call parse_columns(profile_columns_file, prof_names, n_prof_cols, &
           prof_sel, prof_nsel, 'profile')
end subroutine output_init_mesa

! ---------------------------------------------------------------
subroutine output_run_header(star_mass_msun)
      double precision, intent(in) :: star_mass_msun

      if (use_legacy_output) then
         call wrthead(star_mass_msun)
      end if
end subroutine output_run_header

! ---------------------------------------------------------------
subroutine output_write_model(timestep_yr, log_gravity, has_h_shell, &
     h_shell_begin_index, h_shell_end_index, h_shell_mid_index, &
     trial_sign_flag, punch_pending_flag, total_angular_momentum, &
     total_rotational_kinetic_energy)
      use star_info_lib, only: star
      use luout_lib
      double precision, intent(in) :: timestep_yr
      double precision, intent(out) :: log_gravity
      logical, intent(in) :: has_h_shell
      integer, intent(in) :: h_shell_begin_index, h_shell_end_index, &
           h_shell_mid_index
      double precision, intent(in) :: trial_sign_flag
      logical, intent(inout) :: punch_pending_flag
      double precision, intent(in) :: total_angular_momentum, &
           total_rotational_kinetic_energy
      integer :: iprof

      if (use_legacy_output) then
         call wrtout(timestep_yr, log_gravity, has_h_shell, &
              h_shell_begin_index, h_shell_end_index, h_shell_mid_index, &
              trial_sign_flag, punch_pending_flag, total_angular_momentum, &
              total_rotational_kinetic_energy)
      else
! Every quantity below was computed by compute_observables and
! stored in star_info; the writers are pure readers. wrtout computes
! log_gravity as an output on the legacy path; hand back the stored
! value here.
         log_gravity = star%run%log_g_surface
! Profiles and pulse files share the trigger and the number (MESA's
! write_pulse_data_with_profile coupling): every profile_interval
! models, the counter advances and each enabled product is written --
! profile<N>.data, and/or profile<N>.data.GYRE / .data.FGONG per
! pulse_format. The history profile_number column records N.
         iprof = 0
         if (profile_interval > 0 .and. &
             (write_profile_flag .or. write_pulse_flag)) then
            if (mod(star%model_number, profile_interval) == 0) then
               profile_counter = profile_counter + 1
               iprof = profile_counter
               call build_extended
               if (write_profile_flag) call write_profile(iprof)
               if (write_pulse_flag) call write_pulse(iprof)
            end if
         end if
         call write_history_row(iprof)
! Keep the log live during the run, like the history file.
         flush(short_file_unit)
      end if
end subroutine output_write_model

! ---------------------------------------------------------------
! Column-selection files: one column name per line, '!' comments,
! blank lines ignored. Blank/absent control -> all columns in the
! built-in order. Unknown names are fatal (config error), with the
! valid names listed in the log.
subroutine parse_columns(fname, names, ncol, sel, nsel, label)
      use luout_lib
      character(len=*), intent(in) :: fname, label
      integer, intent(in) :: ncol
      character(len=24), intent(in) :: names(ncol)
      integer, intent(out) :: sel(max_cols), nsel
      character(len=256) :: line
      integer :: u, ios, i, j
      logical :: found

      if (len_trim(fname) == 0) then
         nsel = ncol
         do i = 1, ncol
            sel(i) = i
         end do
         return
      end if
      open(newunit=u, file=fname, status='OLD', action='READ', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'cannot open ', trim(label), '_columns_file: ', &
              trim(fname)
         write(short_file_unit,*) 'cannot open ', trim(label), &
              '_columns_file: ', trim(fname)
         stop 1
      end if
      nsel = 0
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         i = index(line, '!')
         if (i > 0) line = line(1:i-1)
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         found = .false.
         do j = 1, ncol
            if (trim(line) == trim(names(j))) then
               nsel = nsel + 1
               sel(nsel) = j
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            write(*,*) 'unknown ', trim(label), ' column: ', trim(line)
            write(short_file_unit,*) 'unknown ', trim(label), &
                 ' column: ', trim(line)
            write(short_file_unit,*) 'valid ', trim(label), ' columns:'
            do j = 1, ncol
               write(short_file_unit,'(2x,a)') trim(names(j))
            end do
            stop 1
         end if
      end do
      close(u)
      if (nsel == 0) then
         write(*,*) trim(label), '_columns_file selected no columns: ', &
              trim(fname)
         stop 1
      end if
end subroutine parse_columns

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
      vals(5) = star%run%log_R_surface
      vals(6) = star%run%log_g_surface
      vals(7) = star%log_Teff
      vals(8) = star%run%core_cz_mass
      vals(9) = star%run%envelope_mass
      vals(10) = star%run%envelope_radius
      vals(11) = star%run%envelope_cz_temperature
      vals(12) = star%run%envelope_cz_density
      vals(13) = star%run%envelope_cz_pressure
      vals(14) = star%run%envelope_cz_o16
      vals(15) = star%run%central_log10_temperature
      vals(16) = star%run%central_log10_density
      vals(17) = star%run%central_log10_pressure
      vals(18) = star%run%central_beta
      vals(19) = star%run%central_degeneracy_eta
      vals(20) = star%xa(i_h1,1)
      vals(21) = star%xa(i_he4,1)
      vals(22) = star%xa(i_metals,1)
      do i = 1, 5
         vals(22+i) = star%luminosity_breakdown(i)
      end do
      vals(28) = star%luminosity_breakdown(i_lum_he_c)
      vals(29) = star%luminosity_breakdown(i_lum_grav)
      vals(30) = star%luminosity_breakdown(i_lum_neu)
      vals(31) = star%flux%cl37_snu_rate
      vals(32) = star%flux%ga71_snu_rate
      do i = 1, 10
         vals(32+i) = star%flux%neutrino_flux_total(i)
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
      vals(67) = star%evo%total_angular_momentum
      vals(68) = star%evo%total_rotational_ke
      vals(69) = star%run%total_moment_of_inertia
      vals(70) = star%run%cz_moment_of_inertia
      vals(71) = star%omega(star%nz)
      vals(72) = star%omega(1)
      vals(73) = star%run%rotation_period_days
      vals(74) = star%run%surf_velocity_kms
      vals(75) = star%turnover%convective_turnover_timescale
      vals(76) = star%run%h_shell_bot_mass
      vals(77) = star%run%h_shell_mid_mass
      vals(78) = star%run%h_shell_top_mass
      vals(79) = star%run%h_shell_bot_radius
      vals(80) = star%run%h_shell_mid_radius
      vals(81) = star%run%h_shell_top_radius
      vals(82) = star%turnover%pphot
      vals(83) = star%star_mass
      vals(84) = star%run%dage*1.0d9
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
         write(hist_unit, '(3(1x,i40))') 1, 2, 3
         write(hist_unit, '(3(1x,a40))') adjustr('version_number'), &
              adjustr('initial_mass'), adjustr('initial_z')
         write(hist_unit, '(1x,a40,2(1x,es40.16e3))') &
              adjustr('"' // trim(yrec_version_string) // '"'), &
              star%star_mass, star%xa(i_metals,star%nz)
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
end subroutine profile_column_names

! Per-zone value for profile column icol at YREC zone index k
! (1 = center .. nz = surface). Sources match putstore's per-shell
! block and the pulse arrays coefft fills every model.
double precision function profile_value(icol, k)
      use star_info_lib, only: star, i_be9, i_c12, i_c13, i_eps_cno, i_eps_grav, i_eps_he3, i_eps_neu, i_eps_pp1, i_eps_pp2, i_eps_pp3, i_grad_actual, i_grad_ad, i_grad_rad, i_h1, i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18
      integer, intent(in) :: icol, k

      select case (icol)
      case (1);  profile_value = 0.0d0   ! zone number set by the writer
      case (2);  profile_value = star%m(k)/solar_mass_cgs
      case (3);  profile_value = star%logR(k)
      case (4);  profile_value = star%logT(k)
      case (5);  profile_value = star%logRho(k)
      case (6);  profile_value = star%logP(k)
      case (7);  profile_value = star%luminosity_lsun(k)
      case (8);  profile_value = star%dm(k)/solar_mass_cgs
      case (9)
         if (star%convective_flag(k)) then
            profile_value = 1.0d0
         else
            profile_value = 0.0d0
         end if
      case (10); profile_value = star%run%adiabatic_index_gamma1(k)
      case (11); profile_value = star%diag%so(k)
      case (12); profile_value = star%diag%del_grad(i_grad_rad,k)
      case (13); profile_value = star%diag%del_grad(i_grad_actual,k)
      case (14); profile_value = star%diag%del_grad(i_grad_ad,k)
      case (15); profile_value = star%diag%svel(k)
      case (16); profile_value = star%diag%sbeta(k)
      case (17); profile_value = star%diag%seta(k)
      case (18)
         profile_value = exp(ln10*(cgl - 2.0d0*star%logR(k)))*star%m(k)
      case (19); profile_value = star%diag%seg(i_eps_pp1,k)
      case (20); profile_value = star%diag%seg(i_eps_pp2,k)
      case (21); profile_value = star%diag%seg(i_eps_pp3,k)
      case (22); profile_value = star%diag%seg(i_eps_cno,k)
      case (23); profile_value = star%diag%seg(i_eps_he3,k)
      case (24); profile_value = star%diag%sesum(k)
      case (25); profile_value = star%diag%seg(i_eps_neu,k)
      case (26); profile_value = star%diag%seg(i_eps_grav,k)
      case (27); profile_value = star%xa(i_h1,k)
      case (28); profile_value = star%xa(i_h2,k)
      case (29); profile_value = star%xa(i_he3,k)
      case (30); profile_value = star%xa(i_he4,k)
      case (31); profile_value = star%xa(i_li6,k)
      case (32); profile_value = star%xa(i_li7,k)
      case (33); profile_value = star%xa(i_be9,k)
      case (34); profile_value = star%xa(i_c12,k)
      case (35); profile_value = star%xa(i_c13,k)
      case (36); profile_value = star%xa(i_n14,k)
      case (37); profile_value = star%xa(i_n15,k)
      case (38); profile_value = star%xa(i_o16,k)
      case (39); profile_value = star%xa(i_o17,k)
      case (40); profile_value = star%xa(i_o18,k)
      case (41); profile_value = star%xa(i_metals,k)
      case (42); profile_value = star%omega(k)
      case (43); profile_value = star%j_rot(k)
      case (44); profile_value = star%i_rot(k)
      case (45); profile_value = star%fp_rot(k)
      case (46); profile_value = star%ft_rot(k)
      case (47); profile_value = star%circ%es_circulation_velocity(k)
      case (48); profile_value = star%circ%gsf_circulation_velocity(k)
      case (49); profile_value = star%circ%secular_shear_velocity(k)
      case (50); profile_value = star%pulse%pulse_specific_heat(k)
      case (51); profile_value = -star%pulse%pulse_dlnrho_dlnt(k)
      case (52); profile_value = star%pulse%pulse_mean_molecular_weight(k)
      case (53)
         if (star%pulse%pulse_electron_mean_molecular_weight(k) > 0.0d0) then
            profile_value = &
                 1.0d0/star%pulse%pulse_electron_mean_molecular_weight(k)
         else
            profile_value = 0.0d0
         end if
      case default
         profile_value = 0.0d0
      end select
end function profile_value

! Regenerate the envelope/atmosphere structures for the CONVERGED
! model and build the inward-to-outward index map. The envelope the
! solver last integrated belongs to some trial (Teff, L) -- or was
! never integrated at all, when the envelope triangle interpolated --
! so a fresh atm_get at the converged values is required (this is
! what io/stitch.f90 does for the legacy .store profiles, with the
! same fixed output step sizes).
subroutine build_extended
      use star_info_lib, only: star, i_h1, i_metals
      use envstruct_lib
      use atmstruct_lib
      use envint_lib, only: atm_get
      integer :: j, i, jerr
      double precision :: atm_beg0, atm_min0, atm_max0
      double precision :: env_beg0, env_min0, env_max0
      double precision :: b, gl, rl, ateffl, plim, dum1(4), dum2(3), &
           dum3(3), dum4(3)
      integer :: ixx, idum, katm, kenv, ksaha
      logical :: lprt, lsbc0, lpulpt

! interior always present
      n_ext = 0
      do j = 1, star%nz
         n_ext = n_ext + 1
         ext_region(n_ext) = 1
         ext_index(n_ext) = j
      end do
      if (.not. calc_envelope_flag) return

! ---- re-integrate at the converged model (stitch's recipe) ----
      atm_beg0 = atm_step_begin
      atm_min0 = atm_step_min
      atm_max0 = atm_step_max
      env_beg0 = env_step_begin
      env_min0 = env_step_min
      env_max0 = env_step_max
      atm_step_begin = atm_step_size
      atm_step_min = atm_step_size
      atm_step_max = atm_step_size
      env_step_begin = envelope_step_size
      env_step_min = envelope_step_size
      env_step_max = envelope_step_size

      idum = 0
      ixx = 0
      katm = 0
      kenv = 0
      ksaha = 0
      lprt = .false.
      lsbc0 = .false.
      lpulpt = .false.
      b = exp(ln10*star%log_L)
      rl = 0.5d0*(star%log_L + log10_solar_luminosity - 4.0d0*star%log_Teff &
           - c4pil - csigl)
      gl = cgl + star%log_total_mass - rl - rl
      plim = star%logP(star%nz)
      if (star%convective_flag(star%nz) .and. spot_filling_factor /= 0.0d0 &
          .and. spot_temp_contrast /= 1.0d0) then
         ateffl = star%log_Teff - 0.25d0*log10(spot_filling_factor* &
              spot_temp_contrast**4.0d0 + 1.0d0 - spot_filling_factor)
      else
         ateffl = star%log_Teff
      end if
      jerr = 0
      call atm_get(b, star%fp_rot(star%nz), star%ft_rot(star%nz), gl, &
           star%log_total_mass, ixx, lprt, lsbc0, plim, rl, ateffl, &
           star%xa(i_h1,star%nz), star%xa(i_metals,star%nz), dum1, idum, katm, &
           kenv, ksaha, dum2, dum3, dum4, lpulpt, jerr)

      atm_step_begin = atm_beg0
      atm_step_min = atm_min0
      atm_step_max = atm_max0
      env_step_begin = env_beg0
      env_step_min = env_min0
      env_step_max = env_max0
      if (jerr /= 0) return

! envelope: env_struct runs fitting point -> photosphere (envint
! inverts it), so it appends directly. Its innermost point repeats
! the fitting point (to re-integration roundoff, so its radius can
! even land marginally BELOW the last interior point's) -- skip any
! envelope point not strictly outside the interior, keeping the
! extended radius strictly monotonic for GYRE.
      do i = 1, env_struct%num_env_points
         if (n_ext >= max_ext) exit
         if (env_struct%env_log10_radius(i) <= star%logR(star%nz)) cycle
         n_ext = n_ext + 1
         ext_region(n_ext) = 2
         ext_index(n_ext) = i
      end do
! atmosphere: atmo_struct runs outward-in, so walk it in reverse.
      do i = atmo_struct%num_atm_points, 1, -1
         if (n_ext >= max_ext) exit
         n_ext = n_ext + 1
         ext_region(n_ext) = 3
         ext_index(n_ext) = i
      end do
end subroutine build_extended

! Profile column value at extended point j. Envelope/atmosphere points
! carry what those integrations store; quantities they do not track
! (per-species abundances beyond X/Z, burning terms, rotation
! internals) are zero, as io/stitch.f90 also writes them.
double precision function ext_profile_value(icol, j)
      use star_info_lib, only: star, i_h1, i_metals
      use envstruct_lib
      use atmstruct_lib
      integer, intent(in) :: icol, j
      integer :: i

      i = ext_index(j)
      select case (ext_region(j))
      case (1)
         ext_profile_value = profile_value(icol, i)
      case (2)
         select case (icol)
         case (2);  ext_profile_value = exp(ln10*(env_struct%env_log10_mass(i) &
                       + star%log_total_mass))/solar_mass_cgs
         case (3);  ext_profile_value = env_struct%env_log10_radius(i)
         case (4);  ext_profile_value = env_struct%env_log10_temperature(i)
         case (5);  ext_profile_value = env_struct%env_log10_density(i)
         case (6);  ext_profile_value = env_struct%env_log10_pressure(i)
         case (7);  ext_profile_value = env_struct%env_luminosity(i)
         case (9)
            if (env_struct%env_convective_flag(i)) then
               ext_profile_value = 1.0d0
            else
               ext_profile_value = 0.0d0
            end if
         case (10); ext_profile_value = env_struct%env_gamma1(i)
         case (11); ext_profile_value = env_struct%env_opacity(i)
         case (12); ext_profile_value = env_struct%env_gradients(1,i)
         case (13); ext_profile_value = env_struct%env_gradients(2,i)
         case (14); ext_profile_value = env_struct%env_gradients(3,i)
         case (27); ext_profile_value = env_struct%env_hydrogen_fraction(i)
         case (41); ext_profile_value = env_struct%env_metal_fraction(i)
         case (42); ext_profile_value = star%omega(star%nz)
         case (50); ext_profile_value = env_struct%env_specific_heat_cp(i)
         case (51); ext_profile_value = -env_struct%env_dlnrho_dlnt(i)
         case default; ext_profile_value = 0.0d0
         end select
      case (3)
         select case (icol)
         case (2);  ext_profile_value = star%star_mass
         case (3);  ext_profile_value = log10(exp(ln10* &
                       env_struct%env_log10_radius(env_struct%num_env_points)) &
                       + atmo_struct%atmo_delta_depth(i))
         case (4);  ext_profile_value = atmo_struct%atmo_log10_temperature(i)
         case (5);  ext_profile_value = atmo_struct%atmo_log10_density(i)
         case (6);  ext_profile_value = atmo_struct%atmo_log10_pressure(i)
         case (7);  ext_profile_value = exp(ln10*star%log_L)
         case (10); ext_profile_value = atmo_struct%atmo_gamma1(i)
         case (11); ext_profile_value = atmo_struct%atmo_opacity(i)
         case (12); ext_profile_value = atmo_struct%atmo_gradients(1,i)
         case (13); ext_profile_value = atmo_struct%atmo_gradients(2,i)
         case (14); ext_profile_value = atmo_struct%atmo_gradients(3,i)
         case (16); ext_profile_value = atmo_struct%atmo_beta(i)
         case (27); ext_profile_value = star%xa(i_h1,star%nz)
         case (41); ext_profile_value = star%xa(i_metals,star%nz)
         case (42); ext_profile_value = star%omega(star%nz)
         case (50); ext_profile_value = atmo_struct%atmo_specific_heat_cp(i)
         case (51); ext_profile_value = -atmo_struct%atmo_dlnrho_dlnt(i)
         case default; ext_profile_value = 0.0d0
         end select
      case default
         ext_profile_value = 0.0d0
      end select
end function ext_profile_value

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
      path = trim(out_dir) // 'profile' // trim(numstr) // '.data'
      open(newunit=u, file=path, status='REPLACE', action='WRITE')
! global block (MESA profile shape: numbers / names / values)
      write(u, '(5(1x,i40))') 1, 2, 3, 4, 5
      write(u, '(5(1x,a40))') adjustr('model_number'), &
           adjustr('num_zones'), adjustr('star_age'), &
           adjustr('log_Teff'), adjustr('star_mass')
      write(u, '(2(1x,i40),3(1x,es40.16e3))') star%model_number, &
           n_ext, star%run%dage*1.0d9, star%log_Teff, star%star_mass
      write(u, '(a)') ''
      write(u, '(999(1x,i40))') (k, k = 1, prof_nsel)
      write(u, '(999(1x,a40))') &
           (adjustr(trim(names(prof_sel(k)))), k = 1, prof_nsel)
! rows: zone 1 = the outermost point (MESA convention). The extended
! model (interior + envelope + atmosphere, built by build_extended)
! is stored inward-to-outward, so walk it in reverse.
      do k = 1, n_ext
         kz = n_ext - k + 1
         do i = 1, prof_nsel
            if (prof_sel(i) == 1) then
               v = dble(k)
            else
               v = ext_profile_value(prof_sel(i), kz)
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

! Assemble the per-point pulse data for the FULL extended model
! (interior + envelope + atmosphere; build_extended must have run).
! Columns 1-18 are the GYRE schema-101 set; 19-34 the extras FGONG
! needs. Envelope/atmosphere points carry what those integrations
! store; opacity/energy derivative columns they do not track are
! zero there (they matter only for nonadiabatic work), and the
! composition above the fitting point is the surface composition, as
! io/stitch.f90 also writes.
subroutine build_pulse_points(pts)
      use star_info_lib, only: star, i_eps_grav, i_grad_actual, i_grad_ad, i_h1, i_metals
      use envstruct_lib
      use atmstruct_lib
      double precision, intent(out) :: pts(35, n_ext)
      integer :: j, i, k
      double precision :: r, m, P, T, rho, delta, nab, nab_ad, grav

      do j = 1, n_ext
         i = ext_index(j)
         pts(:,j) = 0.0d0
         select case (ext_region(j))
         case (1)
            r = exp(ln10*star%logR(i))
            m = star%m(i)
            P = exp(ln10*star%logP(i))
            T = exp(ln10*star%logT(i))
            rho = exp(ln10*star%logRho(i))
            delta = -star%pulse%pulse_dlnrho_dlnt(i)
            nab = star%diag%del_grad(i_grad_actual,i)
            nab_ad = star%diag%del_grad(i_grad_ad,i)
            pts(3,j) = star%luminosity_lsun(i)*solar_luminosity_cgs
            pts(9,j) = star%run%adiabatic_index_gamma1(i)
            pts(12,j) = star%diag%so(i)
            pts(13,j) = star%pulse%pulse_dlnkap_dlnt(i)
            pts(14,j) = star%pulse%pulse_dlnkap_dlnrho(i)
            pts(15,j) = star%diag%sesum(i)
            pts(16,j) = star%pulse%pulse_dlneps_dlnt(i)
            pts(17,j) = star%pulse%pulse_dlneps_dlnrho(i)
            pts(18,j) = star%omega(i)
            pts(19,j) = star%pulse%pulse_specific_heat(i)
            if (star%pulse%pulse_electron_mean_molecular_weight(i) &
                 > 0.0d0) pts(20,j) = &
                 1.0d0/star%pulse%pulse_electron_mean_molecular_weight(i)
            pts(21,j) = star%xa(i_h1,i)
            pts(22,j) = star%xa(i_metals,i)
            do k = 1, 11
               pts(22+k,j) = star%xa(species_slot(k),i)
            end do
            pts(34,j) = star%diag%seg(i_eps_grav,i)
         case (2)
            r = exp(ln10*env_struct%env_log10_radius(i))
            m = exp(ln10*(env_struct%env_log10_mass(i) + &
                 star%log_total_mass))
            P = exp(ln10*env_struct%env_log10_pressure(i))
            T = exp(ln10*env_struct%env_log10_temperature(i))
            rho = exp(ln10*env_struct%env_log10_density(i))
            delta = -env_struct%env_dlnrho_dlnt(i)
            nab = env_struct%env_gradients(2,i)
            nab_ad = env_struct%env_gradients(3,i)
            pts(3,j) = env_struct%env_luminosity(i)*solar_luminosity_cgs
            pts(9,j) = env_struct%env_gamma1(i)
            pts(12,j) = env_struct%env_opacity(i)
            pts(18,j) = star%omega(star%nz)
            pts(19,j) = env_struct%env_specific_heat_cp(i)
            pts(21,j) = env_struct%env_hydrogen_fraction(i)
            pts(22,j) = env_struct%env_metal_fraction(i)
            do k = 1, 11
               pts(22+k,j) = star%xa(species_slot(k),star%nz)
            end do
         case default   ! atmosphere
            r = exp(ln10* &
                 env_struct%env_log10_radius(env_struct%num_env_points)) &
                 + atmo_struct%atmo_delta_depth(i)
            m = exp(ln10*star%log_total_mass)
            P = exp(ln10*atmo_struct%atmo_log10_pressure(i))
            T = exp(ln10*atmo_struct%atmo_log10_temperature(i))
            rho = exp(ln10*atmo_struct%atmo_log10_density(i))
            delta = -atmo_struct%atmo_dlnrho_dlnt(i)
            nab = atmo_struct%atmo_gradients(2,i)
            nab_ad = atmo_struct%atmo_gradients(3,i)
            pts(3,j) = exp(ln10*star%log_L)*solar_luminosity_cgs
            pts(9,j) = atmo_struct%atmo_gamma1(i)
            pts(12,j) = atmo_struct%atmo_opacity(i)
            pts(18,j) = star%omega(star%nz)
            pts(19,j) = atmo_struct%atmo_specific_heat_cp(i)
            pts(21,j) = star%xa(i_h1,star%nz)
            pts(22,j) = star%xa(i_metals,star%nz)
            do k = 1, 11
               pts(22+k,j) = star%xa(species_slot(k),star%nz)
            end do
         end select
         pts(1,j) = r
         pts(2,j) = m
         pts(4,j) = P
         pts(5,j) = T
         pts(6,j) = rho
         pts(7,j) = nab
         pts(10,j) = nab_ad
         pts(11,j) = delta
         if (r > 0.0d0) then
            grav = exp(ln10*cgl)*m/(r*r)
            pts(8,j) = grav*grav*(rho/P)*delta*(nab_ad - nab)
         end if
      end do
end subroutine build_pulse_points

! star%xa slot for pulse column 22+k (k = 1..11), in FGONG species
! order (column 34 is eps_grav, not a species; be9 has no FGONG
! column).
integer function species_slot(k)
      use star_info_lib, only: i_he3, i_c12, i_c13, i_n14, i_o16, &
           i_h2, i_he4, i_li7, i_n15, i_o17, i_o18
      integer, intent(in) :: k
      integer, parameter :: slots(11) = [i_he3, i_c12, i_c13, i_n14, &
           i_o16, i_h2, i_he4, i_li7, i_n15, i_o17, i_o18]
      species_slot = slots(k)
end function species_slot

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
      double precision, allocatable :: pts(:,:)
      double precision :: mstar_g, rstar_cm, lstar_cgs

      allocate(pts(35, n_ext))
      call build_pulse_points(pts)
      mstar_g = exp(ln10*star%log_total_mass)
      rstar_cm = exp(ln10*(star%run%log_R_surface + log10_solar_radius))
      lstar_cgs = exp(ln10*star%log_L)*solar_luminosity_cgs

      write(numstr, '(i0)') iprof
      if (pulse_format(1:5) == 'FGONG' .or. &
          pulse_format(1:5) == 'fgong') then
         path = trim(out_dir) // 'profile' // trim(numstr) // '.data.FGONG'
         call write_fgong_pulse(n_ext, pts, mstar_g, rstar_cm, &
              lstar_cgs, path)
      else if (pulse_format(1:3) == 'GSM' .or. &
               pulse_format(1:3) == 'gsm') then
         path = trim(out_dir) // 'profile' // trim(numstr) // '.data.GSM'
         call write_gsm_pulse(n_ext, pts, mstar_g, rstar_cm, &
              lstar_cgs, path)
      else
         path = trim(out_dir) // 'profile' // trim(numstr) // '.data.GYRE'
         call write_gyre_ext(n_ext, pts, mstar_g, rstar_cm, &
              lstar_cgs, path)
      end if
      deallocate(pts)
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
