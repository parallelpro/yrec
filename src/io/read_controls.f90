!----------------------------------------------------------------------
! read_controls
!----------------------------------------------------------------------
! Reads a deck into star%ctrl and star%job. Each read starts from the
! pristine defaults: star%ctrl = controls_state() (generated from
! controls_buffer's declarations), seeded into the buffer
! (seed_controls_buffer); read_input reads the deck over the buffer,
! and store_controls_to_star copies the result back. The read's steps
! live in legacy_controls (the namelist variables and their defaults),
! run_files, derived_controls and run_cards. Errors return via ierr.
module controls_reader
      use user_inputs, only: map_user_inputs
      use star_info_lib, only: star, control_nml_override, physics_nml_override
      use controls_sync, only: store_controls_to_star
      use physical_constants
      use io_units
      use legacy_controls
      use run_files, only: assign_unit_numbers, expand_paths, open_run_files
      use derived_controls, only: derive_controls
      use run_cards, only: interpret_kind_cards
      implicit none
      private
      public :: read_controls

contains

subroutine read_controls(ierr)

      use star_info_lib, only: star, controls_state
      use controls_sync, only: seed_controls_buffer, &
           store_controls_to_star
      implicit none

      integer, intent(out) :: ierr

      ierr = 0
      star%ctrl = controls_state()
      call seed_controls_buffer
      call read_input(ierr)
      if (ierr /= 0) return
      call store_controls_to_star

      return
end subroutine read_controls

!----------------------------------------------------------------------
! read_input
!----------------------------------------------------------------------
! Reads the deck: a single-file inlist (&control and &physics) or a
! legacy CONTROL/PHYSICS namelist pair (yrec8.nml1 and yrec8.nml2, or
! the two files named on the command line). Then copies the values to
! their canonical names, opens the run's files, derives the remaining
! controls and turns the kind cards into the run plan.
subroutine read_input(ierr)
      integer, intent(out) :: ierr

      ierr = 0
      call set_legacy_defaults
      call assign_unit_numbers
      call load_table_paths

      print *,''
      print *,'Yale Rotating Evolution Code - YREC, v',yrec_version_string(1:len_trim(yrec_version_string)),' (',git_hash_string(1:len_trim(git_hash_string)),')'

! the run log banner's format, fitted to the version strings (a zero
! width is illegal)
      write(version_fmt, 315) max(1, len_trim(yrec_version_string)), &
           max(1, len_trim(git_hash_string))
      315 format('(''# YREC v'', A', i2.2, ', '' ('', A', i2.2, &
           &        ', '')'')')

! The deck: the first command-line argument (default yrec8.nml1), or
! the override an embedding host (pyyrec) sets, since argv is the host's.
      if (control_nml_override .ne. ' ') then
         control_nml_file = control_nml_override
      else
      call getarg(1, control_nml_file)
      if (control_nml_file(1:2) .eq. ' ') control_nml_file = 'yrec8.nml1'
      end if
      print *, ' '
      call read_namelist_files(ierr)
      if (ierr /= 0) return
      call adopt_canonical_names
      call map_user_inputs(star)
! open_run_files' output_init_mesa reads its controls from star%ctrl
      call store_controls_to_star
      call expand_paths
      call open_run_files(ierr)
      call store_table_paths
      if (ierr /= 0) return
      call derive_controls(ierr)
      if (ierr /= 0) return
      call interpret_kind_cards

end subroutine read_input


! ---------------------------------------------------------------
! The table paths between star%job and their namelist variables. The
! load comes first, so a deck that names no table keeps the previous
! read's path; the store follows the read and the path expansion.
subroutine load_table_paths
      falex06 = star%job%alex06_table_path
      fallard = star%job%allard_table_path
      fatm = star%job%atm_table_path
      ffermi = star%job%fermi_table_path
      fkur = star%job%kurucz_table_path
      fkur2 = star%job%kurucz_table2_path
      flaol = star%job%laol_table_path
      flaol2 = star%job%laol_table2_path
      fliv95 = star%job%opal95_table_path
      flldat = star%job%opal92_table_path
      fmhd1 = star%job%zams_a_table_path
      fmhd2 = star%job%zams_b_table_path
      fmhd3 = star%job%zams_c_table_path
      fmhd4 = star%job%centre1_table_path
      fmhd5 = star%job%centre2_table_path
      fmhd6 = star%job%centre3_table_path
      fmhd7 = star%job%centre4_table_path
      fmhd8 = star%job%centre5_table_path
      fopal2 = star%job%opal92_table2_path
      fpurez = star%job%pure_z_table_path
      fscvh = star%job%scv_h_table_path
      fscvhe = star%job%scv_he_table_path
      fscvz = star%job%scv_z_table_path
      opecalex = star%job%alex95_table_paths
end subroutine load_table_paths

subroutine store_table_paths
      star%job%alex06_table_path = falex06
      star%job%allard_table_path = fallard
      star%job%atm_table_path = fatm
      star%job%fermi_table_path = ffermi
      star%job%kurucz_table_path = fkur
      star%job%kurucz_table2_path = fkur2
      star%job%laol_table_path = flaol
      star%job%laol_table2_path = flaol2
      star%job%opal95_table_path = fliv95
      star%job%opal92_table_path = flldat
      star%job%zams_a_table_path = fmhd1
      star%job%zams_b_table_path = fmhd2
      star%job%zams_c_table_path = fmhd3
      star%job%centre1_table_path = fmhd4
      star%job%centre2_table_path = fmhd5
      star%job%centre3_table_path = fmhd6
      star%job%centre4_table_path = fmhd7
      star%job%centre5_table_path = fmhd8
      star%job%opal92_table2_path = fopal2
      star%job%pure_z_table_path = fpurez
      star%job%scv_h_table_path = fscvh
      star%job%scv_he_table_path = fscvhe
      star%job%scv_z_table_path = fscvz
      star%job%alex95_table_paths = opecalex
      star%job%conductive_table_path = fcondopacp
end subroutine store_table_paths

! ---------------------------------------------------------------
! Read the deck. A file that holds both a &control and a &physics group
! is an inlist, read by read_inlist; anything else is a legacy pair, the
! CONTROL namelist in the first file and the PHYSICS namelist in the
! second. An inlist from before 2026-10, with &star_job and &controls
! groups, is refused with the command that renames its groups.
subroutine read_namelist_files(ierr)
      integer, intent(inout) :: ierr
      logical :: single_inlist
      integer :: ios
      character(len=256) :: msg

      if (nml_file_has_group(control_nml_file, 'star_job')) then
         write(*,'(a)') ' '//trim(control_nml_file)//' names its groups &star_job and &controls;'
         write(*,'(a)') ' they are &control and &physics now. Rename them with'
         write(*,'(a)') '   python3 src/tools/upgrade_inlist.py '//trim(control_nml_file)
         ierr = 1
         return
      end if
      single_inlist = nml_file_has_group(control_nml_file, 'control') .and. &
           nml_file_has_group(control_nml_file, 'physics')

      if (physics_nml_override .ne. ' ') then
         physics_nml_file = physics_nml_override
      else
      call getarg(2, physics_nml_file)
      if (physics_nml_file(1:2) .eq. ' ') physics_nml_file = 'yrec8.nml2'
      end if
! a single-file inlist is both files (copy_inlists_used copies it once)
      if (single_inlist) physics_nml_file = control_nml_file

      if (single_inlist) then
         write(*,*) 'inlist           :  ',control_nml_file(1:len_trim(control_nml_file))
      else
      write(*,*) 'CONTROL namelist :  ',control_nml_file(1:len_trim(control_nml_file))
      write(*,*) 'PHYSICS namelist :  ',physics_nml_file(1:len_trim(physics_nml_file))
      end if

      if (single_inlist) then
      call read_inlist(ierr)
      else
      open(unit=standard_unit, file=control_nml_file, status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the CONTROL namelist file', control_nml_file, ierr)
      if (ierr /= 0) return
      read(unit=standard_unit, nml=control, iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'read the CONTROL namelist from', control_nml_file, ierr, standard_unit)
      if (ierr /= 0) return
      close(standard_unit)
      open(unit=run_unit, file=physics_nml_file, status='OLD', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'open the PHYSICS namelist file', physics_nml_file, ierr)
      if (ierr /= 0) return
      read(unit=run_unit, nml=physics, iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'read the PHYSICS namelist from', physics_nml_file, ierr, run_unit)
      if (ierr /= 0) return
      close(run_unit)
      end if
end subroutine read_namelist_files

! ---------------------------------------------------------------
! Read a single-file inlist: its &control and &physics groups, under the
! new names (defaults/control.defaults, defaults/physics.defaults).
! Both files are generated by tools/gen_inlist.py: inlist_new_decl.inc
! declares the groups -- named like legacy_controls' CONTROL and PHYSICS
! groups, which they hide here -- and the new-name variables that have
! no legacy twin; inlist_new_read.inc seeds the new names from the
! legacy defaults, reads the groups and copies the values back onto the
! legacy names.
subroutine read_inlist(ierr)
      integer, intent(inout) :: ierr
      integer :: ios
      character(len=256) :: msg
      include 'inlist_new_decl.inc'

      include 'inlist_new_read.inc'
end subroutine read_inlist

! ---------------------------------------------------------------
! Copy each control read under a legacy name into its canonical home
! (controls_buffer, integration_parameters or star). map_user_inputs,
! which runs next, reads use_new_nuclear_rates and the s0_* values.
subroutine adopt_canonical_names
      tolerance_fraction = stolr0
      max_stage_index = imax
      extrap_order = nuse
      saha_log10t_cutoff = tscut
      use_new_nuclear_rates = lnewnuc
      weak_screening_threshold = weakscreening
      overshoot_alpha_core = alphac
      overshoot_alpha_envelope = alphae
      diffusion_timestep_factor = dtdif
      num_rotation_structure_iters = itdif1
      max_diffusion_iters = itdif2
      diffuse_lithium = ldifli
      calc_neutrinos = lsnu
      core_overshoot_active = lovstc
      set_initial_omega = lwnew
      use_wind_torque = ljdot0
      use_thoul_diffusion = lthoul
      use_semiconvection = lsemic
      improved_first_guess_flag = lnews
      max_iter_level1 = niter1
      max_iter_level2 = niter2
      max_iter_level3 = niter3
      max_iter_level4 = niter4
      nuclear_logT_cutoffs = tcut
      mixture_weights_seed = vnew
      initial_omega = wnew
      monte_carlo_active = lmonte
      mc_run_start = imbeg
      mc_run_end = imend
      reimers_scaling_factor = creim
      use_reimers_wind = lreimer
      constant_mixing_coeff = cstmixing
      constant_settling_reduction = cstdiffmix
      envelope_overshoot_active = lovste
      rotation_active = lrot
      instability_transport_active = linstb
      mass_accretion_rate = dmdt0
      accreted_composition = compacc
      use_mass_accretion = lmdot
      use_itoh_neutrino_loss = lnulos1
      calc_envelope_flag = lcalcenv
      star%solar_luminosity_cgs = clsun
      star%solar_radius_cgs = crsun
      use_extended_composition = lexcom
      spot_filling_factor = spotf
      spot_temp_contrast = spotx
      spot_depth_varies = lsdepth
      disk_locking_age_gyr = tdisk
      disk_omega_rad_s = pdisk
      disk_locking_active = ldisk
      target_end_age = endage
      timestep_override = setdt
      central_deuterium_stop = end_dcen
      central_hydrogen_stop = end_xcen
      central_helium_stop = end_ycen
      use_mhd_eos = lmhd
      metal_fraction_match_tolerance = optol
      convergence_tolerance = djok
      settling_tolerance = grtol
      coulomb_log_choice = ilambda
      settling_num_iterations = niter_gs
      diffuse_helium_active = ldify
      settling_timestep_fraction = dt_gs
      hydrogen_diffusion_floor = xmin
      helium_diffusion_min = ymin
      use_thoul_fit = lthoulfit
      use_diffusion_z = ldifz
      use_new_diffusion_routines = lnewdif
      atm_error_tol = atmerr
      atm_step_initial = atmd0
      atm_step_begin = atmbeg
      atm_step_min = atmmin
      atm_step_max = atmmax
      env_error_tol = enverr
      env_step_begin = envbeg
      env_step_min = envmin
      env_step_max = envmax
      atm_choice = kttau
      debye_huckel_eta_min = etadh0
      debye_huckel_eta_max = etadh1
      use_debye_huckel_correction = ldh
      monte_carlo_file1_path = fmonte1
      monte_carlo_file2_path = fmonte2
      num_models = nmodls
      first_call_flag = lfirst
      num_runs = numrun
      max_domega_global = dtwind
      use_structure_dt_limits = lptime
      tri_delta_teffl = tridt
      tri_delta_logl = tridl
      atm_step_size = atmstp
      envelope_step_size = envstp
      envelope_generation_flag = lenvg
      helium_flash_active = lkuthe
      new_species_value = xnewcp
      rescale_species_active = lnewcp
      mixture_change_mode = isetmix
      isotope_change_mode = isetiso
      target_carbon_cno_fraction = frac_c
      target_nitrogen_cno_fraction = frac_n
      target_oxygen_cno_fraction = frac_o
      c12_to_c13_ratio = r12_13
      n14_to_n15_ratio = r14_15
      o16_to_o17_ratio = r16_17
      o16_to_o18_ratio = r16_18
      target_metal_fraction = zxmix
      initial_h2_fraction = xh2_ini
      initial_he3_fraction = xhe3_ini
      initial_li6_fraction = xli6_ini
      initial_li7_fraction = xli7_ini
      initial_be9_fraction = xbe9_ini
      initial_b10_fraction = xb10_ini
      initial_b11_fraction = xb11_ini
      difad_velocity_scale = fw
      mixing_velocity_scale = fc
      es_velocity_scale = fes
      gsf_velocity_scale = fgsf
      mu_gradient_scale = fmu
      secular_shear_velocity_scale = fss
      critical_reynolds = rcrit
      es_mixing_scale = fesc
      secular_shear_mixing_scale = fssc
      gsf_mixing_scale = fgsfc
      gsf_inhibition_mode = igsf
      min_abundance = cmin
      absolute_tolerance = abstol
      relative_tolerance = reltol
      max_burn_iterations = kemmax
      extend_core_inward = lcore
      num_core_shells_added = mcore
      core_mass_reduction_factor = fcore
      use_pure_z_table = lpurez
      rewind_short_file = lrwsh
      isochrone_output_active = liso
      initial_x_array = xenv0a
      initial_z_array = zenv0a
      mixing_length_array = cmixla
      has_senv0_array = lsenv0a
      senv0_array = senv0a
      luminosity_tolerance = toll
      radius_tolerance = tolr
      zx_tolerance = tolz
      calibrate_solar_model = lcals
      calibrate_solar_zx = lcalsolzx
      target_solar_zx = calsolzx
      target_solar_age = calsolage
      use_z_ramp = lzramp
      target_luminosity_lsun = xls
      target_star_luminosity_tolerance = xlstol
      target_teff = steff
      target_radius_rsun = sr
      specify_teff_flag = lteff
      calibrate_star_flag = lcalst
      use_opal95_eos = lopale
      use_opal2001_eos = lopale01
      use_opal2006_eos = lopale06
      use_numerical_derivatives = lnumderiv
      laol_table_z1 = zlaol1
      laol_table_z2 = zlaol2
      opal_table_z1 = zopal1
      opal_table_z2 = zopal2
      opal95_single_table_z = zopal951
      alex_table_z1 = zalex1
      kurucz_table_z1 = zkur1
      kurucz_table_z2 = zkur2
      molecular_opacity_logt_min = tmolmin
      molecular_opacity_logt_max = tmolmax
      use_alex06_tables = lalex06
      use_laol89_tables = llaol89
      use_opal92_tables = lopal92
      use_opal95_tables = lopal95
      use_kurucz90_tables = lkur90
      use_alex95_tables = lalex95
      use_conductive_opacity = lcondopacp
      use_diffusion_advection_transport = ldifad
      no_am_transport_in_core = lnoj
      allard_target_feh = alatm_feh
      allard_target_alpha = alatm_alpha
      allard_use_tau100 = laltptau100
      force_solid_body_rotation = lsolid
      solid_body_mode_flag = impjmod
      use_envelope_triangle_dt = ltrist
      pmm_exponent_a = pmma
      pmm_exponent_b = pmmb
      pmm_exponent_c = pmmc
      pmm_exponent_d = pmmd
      pmm_exponent_m = pmmm
      pmm_norm_jdot = pmmjd
      pmm_norm_mdot = pmmmd
      pmm_solar_pressure = pmmsolp
      pmm_solar_omega = pmmsolw
      pmm_solar_turnover_timescale = pmmsoltau
      use_pmm_wind_law = lmwind
      scale_by_rossby_number = lrossby
      scale_by_b_field = lbscale
      wind_law_name = awind
      wind_saturation_omega = wmax
      constant_background_diffusion_coeff = codm
      use_constant_background_diffusion = lcodm
      s0_pp = s0_1_1
      s0_he3he3 = s0_3_3
      s0_he3he4 = s0_3_4
      s0_p_c12 = s0_1_12
      s0_p_c13 = s0_1_13
      s0_p_n14 = s0_1_14
      s0_p_o16 = s0_1_16
      s0_be7_electron = s0_1_be7e
      s0_be7_p = s0_1_be7p
      s0_n15_p_c12_branch = s0_1_15_c12alp
      s0_n15_p_o16_branch = s0_1_15_o16
      s0p_pp = s0p_1_1
      s0p_he3he3 = s0p_3_3
      s0p_he3he4 = s0p_3_4
      s0p_p_c12 = s0p_1_12
      s0p_p_c13 = s0p_1_13
      s0p_p_n14 = s0p_1_14
      s0p_p_o16 = s0p_1_16
      s0pp_p_c12 = s0pp_1_12
      s0pp_p_c13 = s0pp_1_13
      s0pp_p_o16 = s0pp_1_16
      s0p_be7_p = s0p_1_be7p
      s0pp_be7_p = s0pp_1_be7p
      use_scv_eos = lscv
! lnewnuc = .false. takes the cross sections from sstandard (ratios to
! Bahcall & Ulrich 1988); map_user_inputs sets them from the s0_* inputs
! otherwise. An inlist gets the default ratios.
      if (.not. lnewnuc) star%cross_section_scale = sstandard
! the atmosphere choice as read; the run may change star%job%atm_choice
      star%atm_choice_initial = kttau
      star%use_ttau_relation = .false.
end subroutine adopt_canonical_names

! Does this namelist file contain the group &<gname> (in any case, and
! not merely a group whose name begins with gname)?
logical function nml_file_has_group(fname, gname)
      character(len=*), intent(in) :: fname, gname
      character(len=512) :: probe_line
      character(len=1) :: next_char
      integer :: probe_unit, ios, n

      nml_file_has_group = .false.
      n = len_trim(gname) + 1
      open(newunit=probe_unit, file=fname, status='OLD', &
           action='READ', iostat=ios)
      if (ios /= 0) return
      do
         read(probe_unit, '(a)', iostat=ios) probe_line
         if (ios /= 0) exit
         probe_line = adjustl(probe_line)
         next_char = probe_line(n+1:n+1)
         if (lowercase(probe_line(1:n)) == '&' // lowercase(trim(gname)) &
              .and. (next_char == ' ' .or. next_char == '/' .or. &
              next_char == achar(9))) then
            nml_file_has_group = .true.
            exit
         end if
      end do
      close(probe_unit)
end function nml_file_has_group

! text with its ASCII capitals lowered
pure function lowercase(text) result(lowered)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lowered
      integer :: i, code

      lowered = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) &
              lowered(i:i) = achar(code + 32)
      end do
end function lowercase

end module controls_reader
