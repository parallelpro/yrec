!----------------------------------------------------------------------
! history_output
!----------------------------------------------------------------------
! The MESA-layout history file:
! the column name/value tables -- the part that grows with every new
! column -- and the row writer. One row per converged model; columns
! selectable via the history_columns_file control (default selection
! from defaults/history_columns.list -- see output_columns).
!
! Writers are PURE READERS: every quantity is computed by
! compute_observables (core/diagnostics/observables.f90, after the step's
! physics) into star_info. Module state below is per-job (path,
! parsed selection), reset by history_output_init on every job.
module history_output
      use physical_constants
      use output_columns, only: max_cols, parse_columns, &
           hist_default_names, n_hist_default
      implicit none
      private
      public :: history_output_init, write_history_row

      integer, parameter :: n_hist_cols = 89

      character(len=256) :: hist_path = ' '
      integer :: hist_nsel = 0
      integer :: hist_sel(max_cols)

contains

! ---------------------------------------------------------------
subroutine history_output_init(s, out_dir, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
      character(len=*), intent(in) :: out_dir
      integer, intent(out) :: ierr
      character(len=24) :: hist_names(n_hist_cols)
      integer :: prior_unit
      logical :: prior_open

      hist_path = trim(out_dir) // trim(s%ctrl%star_history_name)
! in-process re-entry: a history unit still connected from a
! previous job would silently accumulate; close it so the lazy
! open in write_history_row starts the file over.
      inquire(file=hist_path, opened=prior_open, number=prior_unit)
      if (prior_open) close(prior_unit)

      ierr = 0
      call history_column_names(hist_names)
      call parse_columns(s%ctrl%history_columns_file, hist_names, &
           n_hist_cols, hist_default_names, n_hist_default, &
           hist_sel, hist_nsel, 'history', ierr)
      if (ierr /= 0) return
! model_number / profile_number / num_zones lead the file whenever
! they are selected, in that fixed order (they are columns 1-3 of the
! built-in table, so a stable hoist of indices 1..3 does it).
      call hoist_id_columns(hist_sel, hist_nsel)
end subroutine history_output_init

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
      names(86) = 'delta_nu_rho'
      names(87) = 'delta_nu'
      names(88) = 'delta_Pg'
      names(89) = 'lum_extra'
end subroutine history_column_names

subroutine history_values(s, vals, iprof)
      use star_info_lib, only: star_info, i_h1, i_he4, i_lum_grav, &
           i_lum_he_c, i_lum_neu, i_lum_extra, i_metals
      type(star_info), intent(inout) :: s
      double precision, intent(out) :: vals(n_hist_cols)
      integer, intent(in) :: iprof
      integer :: i

      vals(1) = dble(s%model_number)
      vals(2) = dble(iprof)
      vals(3) = dble(s%nz)
      vals(4) = s%log_L
      vals(5) = s%log_R_surface
      vals(6) = s%log_g_surface
      vals(7) = s%log_Teff
      vals(8) = s%core_cz_mass
      vals(9) = s%envelope_mass
      vals(10) = s%envelope_radius
      vals(11) = s%envelope_cz_temperature
      vals(12) = s%envelope_cz_density
      vals(13) = s%envelope_cz_pressure
      vals(14) = s%envelope_cz_opacity
      vals(15) = s%central_log10_temperature
      vals(16) = s%central_log10_density
      vals(17) = s%central_log10_pressure
      vals(18) = s%central_beta
      vals(19) = s%central_degeneracy_eta
      vals(20) = s%xa(i_h1,1)
      vals(21) = s%xa(i_he4,1)
      vals(22) = s%xa(i_metals,1)
      do i = 1, 5
         vals(22+i) = s%luminosity_breakdown(i)
      end do
      vals(28) = s%luminosity_breakdown(i_lum_he_c)
      vals(29) = s%luminosity_breakdown(i_lum_grav)
      vals(30) = s%luminosity_breakdown(i_lum_neu)
      vals(31) = s%cl37_snu_rate
      vals(32) = s%ga71_snu_rate
      do i = 1, 10
         vals(32+i) = s%neutrino_flux_total(i)
      end do
      do i = 4, 11
         vals(39+i) = s%xa(i,1)
      end do
      do i = 4, 15
         vals(47+i) = s%xa(i,s%nz)
      end do
      vals(63) = s%xa(i_h1,s%nz)
      vals(64) = s%xa(i_he4,s%nz)
      vals(65) = s%xa(i_metals,s%nz)
      vals(66) = s%xa(i_metals,s%nz)/s%xa(i_h1,s%nz)
      vals(67) = s%total_angular_momentum
      vals(68) = s%total_rotational_ke
      vals(69) = s%total_moment_of_inertia
      vals(70) = s%cz_moment_of_inertia
      vals(71) = s%omega(s%nz)
      vals(72) = s%omega(1)
      vals(73) = s%rotation_period_days
      vals(74) = s%surf_velocity_kms
      vals(75) = s%convective_turnover_timescale
      vals(76) = s%h_shell_bot_mass
      vals(77) = s%h_shell_mid_mass
      vals(78) = s%h_shell_top_mass
      vals(79) = s%h_shell_bot_radius
      vals(80) = s%h_shell_mid_radius
      vals(81) = s%h_shell_top_radius
      vals(82) = s%pphot
      vals(83) = s%star_mass
      vals(84) = s%age_gyr*1.0d9
      vals(85) = s%nu_max
      vals(86) = s%delta_nu_rho
      vals(87) = s%delta_nu
      vals(88) = s%delta_Pg
      vals(89) = s%luminosity_breakdown(i_lum_extra)
end subroutine history_values

subroutine write_history_row(s, iprof, ierr)
      use star_info_lib, only: star_info, i_metals
      use io_units, only: check_io
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iprof
      integer, intent(inout) :: ierr
      character(len=24) :: names(n_hist_cols)
      double precision :: vals(n_hist_cols)
      integer :: hist_unit, i, k, ios
      character(len=256) :: msg
      logical :: hist_open

      call history_values(s, vals, iprof)
      inquire(file=hist_path, opened=hist_open, number=hist_unit)
      if (.not. hist_open) then
         open(newunit=hist_unit, file=hist_path, status='REPLACE', &
              action='WRITE', iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'write', hist_path, ierr)
         if (ierr /= 0) return
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
              s%star_mass, s%job%initial_envelope_x, &
              1.0d0 - s%job%initial_envelope_x - s%job%initial_envelope_z, &
              s%job%initial_envelope_z, s%mixing_length_alpha
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

end module history_output
