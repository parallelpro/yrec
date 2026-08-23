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
! update_output_diagnostics (after the step's physics) into star_info.
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
! Every quantity below was computed by update_output_diagnostics and
! stored in star_info; the writers are pure readers. wrtout computes
! log_gravity as an output on the legacy path; hand back the stored
! value here.
         log_gravity = star%run%log_g_surface
         iprof = 0
         if (profile_interval > 0) then
            if (mod(star%model_number, profile_interval) == 0) then
               profile_counter = profile_counter + 1
               iprof = profile_counter
               call write_profile(iprof)
            end if
         end if
         call write_history_row(iprof)
         if (pulse_gyre_interval > 0) then
            if (mod(star%model_number, pulse_gyre_interval) == 0) then
               if (pulse_format(1:5) == 'FGONG' .or. &
                   pulse_format(1:5) == 'fgong') then
                  call write_fgong_pulse(star%nz, star%model_number)
               else
                  call write_gyre_pulse(star%nz, star%model_number, &
                       star%m, star%logRho, star%luminosity_lsun, &
                       star%logP, star%logR, star%logT, star%omega)
               end if
            end if
         end if
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
      names(2) = 'num_zones'
      names(3) = 'star_age'
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
      names(84) = 'profile_number'
end subroutine history_column_names

subroutine history_values(vals, iprof)
      use star_info_lib, only: star
      double precision, intent(out) :: vals(n_hist_cols)
      integer, intent(in) :: iprof
      integer :: i

      vals(1) = dble(star%model_number)
      vals(2) = dble(star%nz)
      vals(3) = star%run%dage*1.0d9
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
      vals(20) = star%xa(1,1)
      vals(21) = star%xa(2,1)
      vals(22) = star%xa(3,1)
      do i = 1, 5
         vals(22+i) = star%luminosity_breakdown(i)
      end do
      vals(28) = star%luminosity_breakdown(8)
      vals(29) = star%luminosity_breakdown(7)
      vals(30) = star%luminosity_breakdown(6)
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
      vals(63) = star%xa(1,star%nz)
      vals(64) = star%xa(2,star%nz)
      vals(65) = star%xa(3,star%nz)
      vals(66) = star%xa(3,star%nz)/star%xa(1,star%nz)
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
      vals(84) = dble(iprof)
end subroutine history_values

subroutine write_history_row(iprof)
      use star_info_lib, only: star
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
         write(hist_unit, '(1x,3i28)') 1, 2, 3
         write(hist_unit, '(1x,3a28)') adjustr('version_number'), &
              adjustr('initial_mass'), adjustr('initial_z')
         write(hist_unit, '(1x,a28,2es28.16e3)') &
              adjustr(yrec_version_string(1:min(28, &
              len_trim(yrec_version_string)))), &
              star%star_mass, star%xa(3,star%nz)
         write(hist_unit, '(a)') ''
         write(hist_unit, '(1x,999i28)') (k, k = 1, hist_nsel)
         write(hist_unit, '(1x,999a28)') &
              (adjustr(trim(names(hist_sel(k)))), k = 1, hist_nsel)
      end if
      write(hist_unit, '(1x,999es28.16e3)') &
           (vals(hist_sel(i)), i = 1, hist_nsel)
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
      use star_info_lib, only: star
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
      case (12); profile_value = star%diag%del_grad(1,k)
      case (13); profile_value = star%diag%del_grad(2,k)
      case (14); profile_value = star%diag%del_grad(3,k)
      case (15); profile_value = star%diag%svel(k)
      case (16); profile_value = star%diag%sbeta(k)
      case (17); profile_value = star%diag%seta(k)
      case (18)
         profile_value = exp(ln10*(cgl - 2.0d0*star%logR(k)))*star%m(k)
      case (19); profile_value = star%diag%seg(1,k)
      case (20); profile_value = star%diag%seg(2,k)
      case (21); profile_value = star%diag%seg(3,k)
      case (22); profile_value = star%diag%seg(4,k)
      case (23); profile_value = star%diag%seg(5,k)
      case (24); profile_value = star%diag%sesum(k)
      case (25); profile_value = star%diag%seg(6,k)
      case (26); profile_value = star%diag%seg(7,k)
      case (27); profile_value = star%xa(1,k)
      case (28); profile_value = star%xa(12,k)
      case (29); profile_value = star%xa(4,k)
      case (30); profile_value = star%xa(2,k)
      case (31); profile_value = star%xa(13,k)
      case (32); profile_value = star%xa(14,k)
      case (33); profile_value = star%xa(15,k)
      case (34); profile_value = star%xa(5,k)
      case (35); profile_value = star%xa(6,k)
      case (36); profile_value = star%xa(7,k)
      case (37); profile_value = star%xa(8,k)
      case (38); profile_value = star%xa(9,k)
      case (39); profile_value = star%xa(10,k)
      case (40); profile_value = star%xa(11,k)
      case (41); profile_value = star%xa(3,k)
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
      write(u, '(1x,5i28)') 1, 2, 3, 4, 5
      write(u, '(1x,5a28)') adjustr('model_number'), &
           adjustr('num_zones'), adjustr('star_age'), &
           adjustr('log_Teff'), adjustr('star_mass')
      write(u, '(1x,2i28,3es28.16e3)') star%model_number, star%nz, &
           star%run%dage*1.0d9, star%log_Teff, star%star_mass
      write(u, '(a)') ''
      write(u, '(1x,999i28)') (k, k = 1, prof_nsel)
      write(u, '(1x,999a28)') &
           (adjustr(trim(names(prof_sel(k)))), k = 1, prof_nsel)
! rows: zone 1 = the surface (MESA convention); YREC stores
! center (1) .. surface (nz), so walk the arrays in reverse.
      do k = 1, star%nz
         kz = star%nz - k + 1
         do i = 1, prof_nsel
            if (prof_sel(i) == 1) then
               v = dble(k)
            else
               v = profile_value(prof_sel(i), kz)
            end if
            write(u, '(es28.16e3)', advance='no') v
         end do
         write(u, '(a)') ''
      end do
      close(u)
end subroutine write_profile

end module yrec_output
