!----------------------------------------------------------------------
! write_history
!----------------------------------------------------------------------
! New (2026, MESA-style output): writes one history row per converged
! model in MESA's history.data layout -- 3 header rows of run-global
! data, a blank row, a column-number row, a column-name row, then one
! fixed-width row per model. pandas.read_fwf(path, skiprows=5) and
! mesa_reader both parse it directly.
!
! Called from io/wrtout.f90 at exactly the point the legacy .track row
! would be written (track_file_version 0 branch), when
! use_legacy_output is .false.; the caller passes the same locals that
! feed the legacy write, so the two streams stay numerically
! interchangeable column for column. The column set is the .track v0
! set under MESA-vocabulary names.
!
! The file is <track path with .track replaced by .history>. It is
! opened (REPLACE) on first use in a process/run and found again by
! INQUIRE thereafter -- no SAVEd state, so re-entrant runs (which
! close all units in the prologue) start a fresh file exactly like a
! fresh process.
subroutine write_history(itrack, radius_log_surface, log_gravity, &
     total_angular_momentum, total_rotational_kinetic_energy, &
     total_moment_of_inertia, cz_moment_of_inertia, &
     rotation_period_days, equatorial_velocity_kms, &
     h_shell_begin_mass, h_shell_mid_mass2, h_shell_end_mass, &
     h_shell_begin_radius, h_shell_mid_radius, h_shell_end_radius)
      use star_info_lib, only: star
      use const_lib
      implicit none

      integer, intent(in) :: itrack
      double precision, intent(in) :: radius_log_surface, log_gravity, &
           total_angular_momentum, total_rotational_kinetic_energy, &
           total_moment_of_inertia, cz_moment_of_inertia, &
           rotation_period_days, equatorial_velocity_kms, &
           h_shell_begin_mass, h_shell_mid_mass2, h_shell_end_mass, &
           h_shell_begin_radius, h_shell_mid_radius, h_shell_end_radius

      integer, parameter :: ncol = 83
      character(len=24) :: names(ncol)
      double precision :: vals(ncol)
      character(len=512) :: track_path, history_path
      integer :: hist_unit, i, n
      logical :: track_named, hist_open

! ---- column names (MESA vocabulary where an equivalent exists) ----
      names(1) = 'model_number'
      names(2) = 'num_zones'
      names(3) = 'star_age'          ! years (MESA convention; .track is Gyr)
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

! ---- values: same sources as the legacy .track v0 row ----
      vals(1) = dble(star%model_number)
      vals(2) = dble(star%nz)
      vals(3) = star%run%dage*1.0d9
      vals(4) = star%log_L
      vals(5) = radius_log_surface
      vals(6) = log_gravity
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
      vals(67) = total_angular_momentum
      vals(68) = total_rotational_kinetic_energy
      vals(69) = total_moment_of_inertia
      vals(70) = cz_moment_of_inertia
      vals(71) = star%omega(star%nz)
      vals(72) = star%omega(1)
      vals(73) = rotation_period_days
      vals(74) = equatorial_velocity_kms
      vals(75) = star%turnover%convective_turnover_timescale
      vals(76) = h_shell_begin_mass
      vals(77) = h_shell_mid_mass2
      vals(78) = h_shell_end_mass
      vals(79) = h_shell_begin_radius
      vals(80) = h_shell_mid_radius
      vals(81) = h_shell_end_radius
      vals(82) = star%turnover%pphot
      vals(83) = star%star_mass

! ---- locate/open the history file next to the .track ----
      inquire(unit=itrack, named=track_named, name=track_path)
      if (.not. track_named) track_path = 'output/yrec'
      n = len_trim(track_path)
      if (n > 6 .and. track_path(n-5:n) == '.track') then
         history_path = track_path(1:n-6) // '.history'
      else
         history_path = track_path(1:n) // '.history'
      end if
      inquire(file=history_path, opened=hist_open, number=hist_unit)
      if (.not. hist_open) then
         open(newunit=hist_unit, file=history_path, status='REPLACE', &
              action='WRITE')
         call write_header(hist_unit)
      end if

      write(hist_unit, '(1x,i28,i28,81es28.16e3)') &
           star%model_number, star%nz, (vals(i), i = 3, ncol)
      return

contains

subroutine write_header(u)
      integer, intent(in) :: u
      integer :: k

! run-global block (3 columns), MESA history.data shape:
! numbers / names / values, then a blank line, then the per-model
! column-number and column-name rows. Data rows start at line 7, so
! pandas.read_fwf(..., skiprows=5) sees the name row as the header.
      write(u, '(1x,3i28)') 1, 2, 3
      write(u, '(1x,3a28)') adjustr('version_number'), &
           adjustr('initial_mass'), adjustr('initial_z')
      write(u, '(1x,a28,2es28.16e3)') &
           adjustr(yrec_version_string(1:min(28,len_trim(yrec_version_string)))), &
           star%star_mass, star%xa(3,star%nz)
      write(u, '(a)') ''
      write(u, '(1x,83i28)') (k, k = 1, ncol)
      write(u, '(1x,83a28)') (adjustr(trim(names(k))), k = 1, ncol)
end subroutine write_header

end subroutine write_history
