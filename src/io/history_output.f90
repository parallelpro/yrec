!----------------------------------------------------------------------
! history_output
!----------------------------------------------------------------------
! The history file, in MESA's history.data layout: its column table and
! its row writer. One row per converged model; the columns are selected
! by the history_columns_file control (the default selection is
! defaults/history_columns.list, see output_columns).
!
! Each column has a named index ih_<name>, a row in history_columns
! keyed by it (its name and kind), and a value in history_values. The
! writer only reads the star: compute_observables
! (core/diagnostics/observables.f90) computes every quantity first.
! The file's state (its path, unit and selection) is a history_file
! record, part of yrec_output's output_state, which run_yrec owns.
module history_output
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use yrec_version, only: yrec_version_string
      use output_columns, only: max_cols, parse_columns, output_column, &
           real_column, integer_column, hist_default_names
      implicit none
      private
      public :: history_file, history_output_init, write_history_row, &
           close_history, history_columns

! The history columns' named indices, in the order of the column table.
      integer, parameter :: &
           ih_model_number = 1, ih_profile_number = 2, ih_num_zones = 3, &
           ih_log_L = 4, ih_log_R = 5, ih_log_g = 6, &
           ih_log_Teff = 7, ih_mass_conv_core = 8, ih_cz_bot_mass = 9, &
           ih_cz_bot_radius = 10, ih_cz_bot_temperature = 11, ih_cz_bot_density = 12, &
           ih_cz_bot_pressure = 13, ih_cz_bot_opacity = 14, ih_log_center_T = 15, &
           ih_log_center_Rho = 16, ih_log_center_P = 17, ih_center_beta = 18, &
           ih_center_degeneracy = 19, ih_center_h1 = 20, ih_center_he4 = 21, &
           ih_center_z = 22, ih_lum_ppI = 23, ih_lum_ppII = 24, &
           ih_lum_ppIII = 25, ih_lum_cno = 26, ih_lum_3alpha = 27, &
           ih_lum_heC = 28, ih_lum_grav = 29, ih_lum_neu = 30, &
           ih_snu_cl37 = 31, ih_snu_ga71 = 32, ih_neut_pp = 33, &
           ih_neut_pep = 34, ih_neut_hep = 35, ih_neut_be7 = 36, &
           ih_neut_b8 = 37, ih_neut_n13 = 38, ih_neut_o15 = 39, &
           ih_neut_f17 = 40, ih_neut_diag_1 = 41, ih_neut_diag_2 = 42, &
           ih_center_he3 = 43, ih_center_c12 = 44, ih_center_c13 = 45, &
           ih_center_n14 = 46, ih_center_n15 = 47, ih_center_o16 = 48, &
           ih_center_o17 = 49, ih_center_o18 = 50, ih_surface_he3 = 51, &
           ih_surface_c12 = 52, ih_surface_c13 = 53, ih_surface_n14 = 54, &
           ih_surface_n15 = 55, ih_surface_o16 = 56, ih_surface_o17 = 57, &
           ih_surface_o18 = 58, ih_surface_h2 = 59, ih_surface_li6 = 60, &
           ih_surface_li7 = 61, ih_surface_be9 = 62, ih_surface_h1 = 63, &
           ih_surface_he4 = 64, ih_surface_z = 65, ih_surface_z_div_x = 66, &
           ih_total_angular_momentum = 67, ih_total_rotational_ke = 68, ih_total_moment_of_inertia = 69, &
           ih_cz_moment_of_inertia = 70, ih_surf_avg_omega = 71, ih_center_omega = 72, &
           ih_rotation_period_days = 73, ih_surf_avg_v_rot = 74, ih_tau_conv_cz = 75, &
           ih_h_shell_bot_mass = 76, ih_h_shell_mid_mass = 77, ih_h_shell_top_mass = 78, &
           ih_h_shell_bot_radius = 79, ih_h_shell_mid_radius = 80, ih_h_shell_top_radius = 81, &
           ih_log_P_photosphere = 82, ih_star_mass = 83, ih_star_age = 84, &
           ih_nu_max = 85, ih_delta_nu_rho = 86, ih_delta_nu = 87, &
           ih_delta_Pg = 88, ih_lum_extra = 89
      integer, parameter :: n_hist_cols = ih_lum_extra
! The id columns, which lead the file in this order whenever selected.
      integer, parameter :: id_columns(3) = &
           [ih_model_number, ih_profile_number, ih_num_zones]

! One job's history file: its path, its unit once the first row has
! opened it, and the selected columns (indices into history_columns).
      type :: history_file
         character(len=256) :: path = ' '
         integer :: unit = 0
         logical :: is_open = .false.
         integer :: nsel = 0
         integer :: sel(max_cols) = 0
      end type history_file

contains

! ---------------------------------------------------------------
! The job's history file: a file a previous job left open is closed (an
! in-process rerun), and the path and the column selection are set; the
! first row opens the file.
subroutine history_output_init(s, history, out_dir, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(history_file), intent(inout) :: history
      character(len=*), intent(in) :: out_dir
      integer, intent(out) :: ierr

      call close_history(history)
      history%path = trim(out_dir) // trim(s%physics%star_history_name)
      call parse_columns(s%control%history_columns_file, history_columns(), &
           hist_default_names, history%sel, history%nsel, 'history', ierr)
      if (ierr /= 0) return
      call hoist_id_columns(history%sel, history%nsel)
end subroutine history_output_init

! The history file closed, if a row has opened it.
subroutine close_history(history)
      type(history_file), intent(inout) :: history
      if (history%is_open) close(history%unit)
      history%is_open = .false.
end subroutine close_history

! ---------------------------------------------------------------
! The history column table: each column's name and kind, at its named
! index.
function history_columns() result(cols)
      type(output_column) :: cols(n_hist_cols)

      cols(ih_model_number)            = output_column('model_number', integer_column)
      cols(ih_profile_number)          = output_column('profile_number', integer_column)
      cols(ih_num_zones)               = output_column('num_zones', integer_column)
      cols(ih_log_L)                   = output_column('log_L', real_column)
      cols(ih_log_R)                   = output_column('log_R', real_column)
      cols(ih_log_g)                   = output_column('log_g', real_column)
      cols(ih_log_Teff)                = output_column('log_Teff', real_column)
      cols(ih_mass_conv_core)          = output_column('mass_conv_core', real_column)
      cols(ih_cz_bot_mass)             = output_column('cz_bot_mass', real_column)
      cols(ih_cz_bot_radius)           = output_column('cz_bot_radius', real_column)
      cols(ih_cz_bot_temperature)      = output_column('cz_bot_temperature', real_column)
      cols(ih_cz_bot_density)          = output_column('cz_bot_density', real_column)
      cols(ih_cz_bot_pressure)         = output_column('cz_bot_pressure', real_column)
      cols(ih_cz_bot_opacity)          = output_column('cz_bot_opacity', real_column)
      cols(ih_log_center_T)            = output_column('log_center_T', real_column)
      cols(ih_log_center_Rho)          = output_column('log_center_Rho', real_column)
      cols(ih_log_center_P)            = output_column('log_center_P', real_column)
      cols(ih_center_beta)             = output_column('center_beta', real_column)
      cols(ih_center_degeneracy)       = output_column('center_degeneracy', real_column)
      cols(ih_center_h1)               = output_column('center_h1', real_column)
      cols(ih_center_he4)              = output_column('center_he4', real_column)
      cols(ih_center_z)                = output_column('center_z', real_column)
      cols(ih_lum_ppI)                 = output_column('lum_ppI', real_column)
      cols(ih_lum_ppII)                = output_column('lum_ppII', real_column)
      cols(ih_lum_ppIII)               = output_column('lum_ppIII', real_column)
      cols(ih_lum_cno)                 = output_column('lum_cno', real_column)
      cols(ih_lum_3alpha)              = output_column('lum_3alpha', real_column)
      cols(ih_lum_heC)                 = output_column('lum_heC', real_column)
      cols(ih_lum_grav)                = output_column('lum_grav', real_column)
      cols(ih_lum_neu)                 = output_column('lum_neu', real_column)
      cols(ih_snu_cl37)                = output_column('snu_cl37', real_column)
      cols(ih_snu_ga71)                = output_column('snu_ga71', real_column)
      cols(ih_neut_pp)                 = output_column('neut_pp', real_column)
      cols(ih_neut_pep)                = output_column('neut_pep', real_column)
      cols(ih_neut_hep)                = output_column('neut_hep', real_column)
      cols(ih_neut_be7)                = output_column('neut_be7', real_column)
      cols(ih_neut_b8)                 = output_column('neut_b8', real_column)
      cols(ih_neut_n13)                = output_column('neut_n13', real_column)
      cols(ih_neut_o15)                = output_column('neut_o15', real_column)
      cols(ih_neut_f17)                = output_column('neut_f17', real_column)
      cols(ih_neut_diag_1)             = output_column('neut_diag_1', real_column)
      cols(ih_neut_diag_2)             = output_column('neut_diag_2', real_column)
      cols(ih_center_he3)              = output_column('center_he3', real_column)
      cols(ih_center_c12)              = output_column('center_c12', real_column)
      cols(ih_center_c13)              = output_column('center_c13', real_column)
      cols(ih_center_n14)              = output_column('center_n14', real_column)
      cols(ih_center_n15)              = output_column('center_n15', real_column)
      cols(ih_center_o16)              = output_column('center_o16', real_column)
      cols(ih_center_o17)              = output_column('center_o17', real_column)
      cols(ih_center_o18)              = output_column('center_o18', real_column)
      cols(ih_surface_he3)             = output_column('surface_he3', real_column)
      cols(ih_surface_c12)             = output_column('surface_c12', real_column)
      cols(ih_surface_c13)             = output_column('surface_c13', real_column)
      cols(ih_surface_n14)             = output_column('surface_n14', real_column)
      cols(ih_surface_n15)             = output_column('surface_n15', real_column)
      cols(ih_surface_o16)             = output_column('surface_o16', real_column)
      cols(ih_surface_o17)             = output_column('surface_o17', real_column)
      cols(ih_surface_o18)             = output_column('surface_o18', real_column)
      cols(ih_surface_h2)              = output_column('surface_h2', real_column)
      cols(ih_surface_li6)             = output_column('surface_li6', real_column)
      cols(ih_surface_li7)             = output_column('surface_li7', real_column)
      cols(ih_surface_be9)             = output_column('surface_be9', real_column)
      cols(ih_surface_h1)              = output_column('surface_h1', real_column)
      cols(ih_surface_he4)             = output_column('surface_he4', real_column)
      cols(ih_surface_z)               = output_column('surface_z', real_column)
      cols(ih_surface_z_div_x)         = output_column('surface_z_div_x', real_column)
      cols(ih_total_angular_momentum)  = output_column('total_angular_momentum', real_column)
      cols(ih_total_rotational_ke)     = output_column('total_rotational_ke', real_column)
      cols(ih_total_moment_of_inertia) = output_column('total_moment_of_inertia', real_column)
      cols(ih_cz_moment_of_inertia)    = output_column('cz_moment_of_inertia', real_column)
      cols(ih_surf_avg_omega)          = output_column('surf_avg_omega', real_column)
      cols(ih_center_omega)            = output_column('center_omega', real_column)
      cols(ih_rotation_period_days)    = output_column('rotation_period_days', real_column)
      cols(ih_surf_avg_v_rot)          = output_column('surf_avg_v_rot', real_column)
      cols(ih_tau_conv_cz)             = output_column('tau_conv_cz', real_column)
      cols(ih_h_shell_bot_mass)        = output_column('h_shell_bot_mass', real_column)
      cols(ih_h_shell_mid_mass)        = output_column('h_shell_mid_mass', real_column)
      cols(ih_h_shell_top_mass)        = output_column('h_shell_top_mass', real_column)
      cols(ih_h_shell_bot_radius)      = output_column('h_shell_bot_radius', real_column)
      cols(ih_h_shell_mid_radius)      = output_column('h_shell_mid_radius', real_column)
      cols(ih_h_shell_top_radius)      = output_column('h_shell_top_radius', real_column)
      cols(ih_log_P_photosphere)       = output_column('log_P_photosphere', real_column)
      cols(ih_star_mass)               = output_column('star_mass', real_column)
      cols(ih_star_age)                = output_column('star_age', real_column)
! opt-in columns (commented out in defaults/history_columns.list)
      cols(ih_nu_max)                  = output_column('nu_max', real_column)
      cols(ih_delta_nu_rho)            = output_column('delta_nu_rho', real_column)
      cols(ih_delta_nu)                = output_column('delta_nu', real_column)
      cols(ih_delta_Pg)                = output_column('delta_Pg', real_column)
      cols(ih_lum_extra)               = output_column('lum_extra', real_column)
end function history_columns

! ---------------------------------------------------------------
! Every column's value for this model, at its named index (iprof: the
! profile written with this model, 0 for none).
subroutine history_values(s, vals, iprof)
      use star_info_lib, only: star_info, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, i_li6, &
           i_li7, i_be9, i_lum_pp1, i_lum_pp2, i_lum_pp3, i_lum_cno, &
           i_lum_3alpha, i_lum_he_c, i_lum_grav, i_lum_neu, i_lum_extra
      use neutrino_flux_layout, only: i_nu_pp, i_nu_pep, i_nu_hep, &
           i_nu_be7, i_nu_b8, i_nu_n13, i_nu_o15, i_nu_f17, i_nu_he3he3, &
           i_nu_he3he4
      type(star_info), intent(in) :: s
      real(dp), intent(out) :: vals(n_hist_cols)
      integer, intent(in) :: iprof

      vals(ih_model_number)            = dble(s%model_number)
      vals(ih_profile_number)          = dble(iprof)
      vals(ih_num_zones)               = dble(s%nz)
      vals(ih_log_L)                   = s%log_L_lsun
      vals(ih_log_R)                   = s%log_R_surface_rsun
      vals(ih_log_g)                   = s%log_g_surface
      vals(ih_log_Teff)                = s%log_Teff
      vals(ih_mass_conv_core)          = s%core_cz_mass_msun
      vals(ih_cz_bot_mass)             = s%envelope_mass_msun
      vals(ih_cz_bot_radius)           = s%envelope_radius_rsun
      vals(ih_cz_bot_temperature)      = s%envelope_cz_temperature
      vals(ih_cz_bot_density)          = s%envelope_cz_density
      vals(ih_cz_bot_pressure)         = s%envelope_cz_pressure
      vals(ih_cz_bot_opacity)          = s%envelope_cz_opacity
      vals(ih_log_center_T)            = s%central_log10_temperature
      vals(ih_log_center_Rho)          = s%central_log10_density
      vals(ih_log_center_P)            = s%central_log10_pressure
      vals(ih_center_beta)             = s%central_beta
      vals(ih_center_degeneracy)       = s%central_degeneracy_eta
      vals(ih_center_h1)               = s%xa(i_h1,1)
      vals(ih_center_he4)              = s%xa(i_he4,1)
      vals(ih_center_z)                = s%xa(i_metals,1)
      vals(ih_lum_ppI)                 = s%luminosity_breakdown_lsun(i_lum_pp1)
      vals(ih_lum_ppII)                = s%luminosity_breakdown_lsun(i_lum_pp2)
      vals(ih_lum_ppIII)               = s%luminosity_breakdown_lsun(i_lum_pp3)
      vals(ih_lum_cno)                 = s%luminosity_breakdown_lsun(i_lum_cno)
      vals(ih_lum_3alpha)              = s%luminosity_breakdown_lsun(i_lum_3alpha)
      vals(ih_lum_heC)                 = s%luminosity_breakdown_lsun(i_lum_he_c)
      vals(ih_lum_grav)                = s%luminosity_breakdown_lsun(i_lum_grav)
      vals(ih_lum_neu)                 = s%luminosity_breakdown_lsun(i_lum_neu)
      vals(ih_snu_cl37)                = s%cl37_snu_rate
      vals(ih_snu_ga71)                = s%ga71_snu_rate
      vals(ih_neut_pp)                 = s%neutrino_flux_total_e10_per_cm2_s(i_nu_pp)
      vals(ih_neut_pep)                = s%neutrino_flux_total_e10_per_cm2_s(i_nu_pep)
      vals(ih_neut_hep)                = s%neutrino_flux_total_e10_per_cm2_s(i_nu_hep)
      vals(ih_neut_be7)                = s%neutrino_flux_total_e10_per_cm2_s(i_nu_be7)
      vals(ih_neut_b8)                 = s%neutrino_flux_total_e10_per_cm2_s(i_nu_b8)
      vals(ih_neut_n13)                = s%neutrino_flux_total_e10_per_cm2_s(i_nu_n13)
      vals(ih_neut_o15)                = s%neutrino_flux_total_e10_per_cm2_s(i_nu_o15)
      vals(ih_neut_f17)                = s%neutrino_flux_total_e10_per_cm2_s(i_nu_f17)
      vals(ih_neut_diag_1)             = s%neutrino_flux_total_e10_per_cm2_s(i_nu_he3he3)
      vals(ih_neut_diag_2)             = s%neutrino_flux_total_e10_per_cm2_s(i_nu_he3he4)
      vals(ih_center_he3)              = s%xa(i_he3,1)
      vals(ih_center_c12)              = s%xa(i_c12,1)
      vals(ih_center_c13)              = s%xa(i_c13,1)
      vals(ih_center_n14)              = s%xa(i_n14,1)
      vals(ih_center_n15)              = s%xa(i_n15,1)
      vals(ih_center_o16)              = s%xa(i_o16,1)
      vals(ih_center_o17)              = s%xa(i_o17,1)
      vals(ih_center_o18)              = s%xa(i_o18,1)
      vals(ih_surface_he3)             = s%xa(i_he3,s%nz)
      vals(ih_surface_c12)             = s%xa(i_c12,s%nz)
      vals(ih_surface_c13)             = s%xa(i_c13,s%nz)
      vals(ih_surface_n14)             = s%xa(i_n14,s%nz)
      vals(ih_surface_n15)             = s%xa(i_n15,s%nz)
      vals(ih_surface_o16)             = s%xa(i_o16,s%nz)
      vals(ih_surface_o17)             = s%xa(i_o17,s%nz)
      vals(ih_surface_o18)             = s%xa(i_o18,s%nz)
      vals(ih_surface_h2)              = s%xa(i_h2,s%nz)
      vals(ih_surface_li6)             = s%xa(i_li6,s%nz)
      vals(ih_surface_li7)             = s%xa(i_li7,s%nz)
      vals(ih_surface_be9)             = s%xa(i_be9,s%nz)
      vals(ih_surface_h1)              = s%xa(i_h1,s%nz)
      vals(ih_surface_he4)             = s%xa(i_he4,s%nz)
      vals(ih_surface_z)               = s%xa(i_metals,s%nz)
      vals(ih_surface_z_div_x)         = s%xa(i_metals,s%nz)/s%xa(i_h1,s%nz)
      vals(ih_total_angular_momentum)  = s%total_angular_momentum
      vals(ih_total_rotational_ke)     = s%total_rotational_ke
      vals(ih_total_moment_of_inertia) = s%total_moment_of_inertia
      vals(ih_cz_moment_of_inertia)    = s%cz_moment_of_inertia
      vals(ih_surf_avg_omega)          = s%omega(s%nz)
      vals(ih_center_omega)            = s%omega(1)
      vals(ih_rotation_period_days)    = s%rotation_period_days
      vals(ih_surf_avg_v_rot)          = s%surf_velocity_kms
      vals(ih_tau_conv_cz)             = s%convective_turnover_timescale
      vals(ih_h_shell_bot_mass)        = s%h_shell_bot_mass_msun
      vals(ih_h_shell_mid_mass)        = s%h_shell_mid_mass_msun
      vals(ih_h_shell_top_mass)        = s%h_shell_top_mass_msun
      vals(ih_h_shell_bot_radius)      = s%h_shell_bot_radius_fraction
      vals(ih_h_shell_mid_radius)      = s%h_shell_mid_radius_fraction
      vals(ih_h_shell_top_radius)      = s%h_shell_top_radius_fraction
      vals(ih_log_P_photosphere)       = s%photosphere_logP
      vals(ih_star_mass)               = s%star_mass
      vals(ih_star_age)                = s%age_gyr*1.0d9
      vals(ih_nu_max)                  = s%nu_max_uhz
      vals(ih_delta_nu_rho)            = s%delta_nu_rho_uhz
      vals(ih_delta_nu)                = s%delta_nu_uhz
      vals(ih_delta_Pg)                = s%delta_Pg
      vals(ih_lum_extra)               = s%luminosity_breakdown_lsun(i_lum_extra)
end subroutine history_values

! This model's row; the first row opens the file and writes its header
! (MESA's layout: the global block, then the column numbers and names).
subroutine write_history_row(s, history, iprof, ierr)
      use star_info_lib, only: star_info
      use io_units, only: check_io
      type(star_info), intent(in) :: s
      type(history_file), intent(inout) :: history
      integer, intent(in) :: iprof
      integer, intent(inout) :: ierr
      type(output_column) :: cols(n_hist_cols)
      real(dp) :: vals(n_hist_cols)
      integer :: hist_unit, i, k, ios
      character(len=256) :: msg

      call history_values(s, vals, iprof)
      cols = history_columns()
      if (.not. history%is_open) then
         open(newunit=history%unit, file=history%path, status='REPLACE', &
              action='WRITE', iostat=ios, iomsg=msg)
         call check_io(ios, msg, 'write', history%path, ierr)
         if (ierr /= 0) return
         history%is_open = .true.
         hist_unit = history%unit
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
              s%star_mass, s%current_card_envelope_x, &
              1.0d0 - s%current_card_envelope_x - s%current_card_envelope_z, &
              s%current_card_envelope_z, s%mixing_length_alpha
         write(hist_unit, '(a)') ''
         write(hist_unit, '(999(1x,i40))') (k, k = 1, history%nsel)
         write(hist_unit, '(999(1x,a40))') &
              (adjustr(trim(cols(history%sel(k))%name)), k = 1, history%nsel)
      end if
      hist_unit = history%unit
      do i = 1, history%nsel
         if (cols(history%sel(i))%kind == integer_column) then
            write(hist_unit, '(1x,i40)', advance='no') &
                 nint(vals(history%sel(i)))
         else
            write(hist_unit, '(1x,es40.16e3)', advance='no') &
                 vals(history%sel(i))
         end if
      end do
      write(hist_unit, '(a)') ''
! Flush so the file is live during the run (tail -f / editor view),
! matching MESA's behavior.
      flush(hist_unit)
end subroutine write_history_row

! ---------------------------------------------------------------
! The id columns moved to the front of a selection, in id_columns'
! order; the other columns keep theirs.
subroutine hoist_id_columns(sel, nsel)
      integer, intent(inout) :: sel(max_cols), nsel
      integer :: tmp(max_cols), n, id, i

      n = 0
      do id = 1, size(id_columns)
         do i = 1, nsel
            if (sel(i) == id_columns(id)) then
               n = n + 1
               tmp(n) = sel(i)
            end if
         end do
      end do
      do i = 1, nsel
         if (all(sel(i) /= id_columns)) then
            n = n + 1
            tmp(n) = sel(i)
         end if
      end do
      sel(1:n) = tmp(1:n)
      nsel = n
end subroutine hoist_id_columns

end module history_output
