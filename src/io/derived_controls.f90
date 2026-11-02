!----------------------------------------------------------------------
! derived_controls
!----------------------------------------------------------------------
! The controls computed from the ones read: the light-element rate
! scales, one OPAL EOS and one low-temperature opacity source, a fresh
! envelope triangle for rotating runs, the starting model's mixture
! changes, the switches that diffusion options imply and the envelope's
! storage temperature; and it writes the run log's header. Its messages go to the terminal and the run log
! through run_messages (report_warning, report_note); the errors are
! controls_check's. Everything here works on the controls' code names,
! in the read's buffer (controls_buffer).
module derived_controls
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use yrec_version, only: yrec_version_string, git_hash_string
      use run_messages, only: report_warning, report_note, log_line
      use controls_buffer, only: atm_choice, envelope_store_log_temperature, &
           be9_palpha_rate_scale, be9_palpha_reference_xs, be9_pd_rate_scale, &
           be9_pd_reference_xs, be9_pg_rate_scale, be9_pg_reference_xs, &
           be9_rate_value1, be9_rate_value2, be9_rate_value3, &
           c12_to_c13_ratio, change_cno_mixture_active, &
           change_isotope_ratios_active, change_species_name, cno_mixture_name, &
           coulomb_log_choice, diffuse_helium_active, diffuse_lithium, &
           initial_b10_fraction, initial_b11_fraction, &
           initial_be9_fraction, initial_h2_fraction, initial_he3_fraction, &
           initial_li6_fraction, initial_li7_fraction, isotope_change_mode, &
           isotope_mixture_name, li6_rate_scale, li6_rate_value, li6_reference_xs, &
           li7_rate_scale, li7_rate_value, li7_reference_xs, mixture_change_mode, &
           n14_to_n15_ratio, new_species_index, new_species_value_mode, &
           new_triangle_init_flag, o16_to_o17_ratio, o16_to_o18_ratio, &
           opal2001_eos_table_path, opal2006_eos_table_path, opal95_eos_table_path, &
           opal_eos_table_path, opal_eos_vintage, pulse_gyre_interval, &
           rescale_species_active, rezone_env_logt0, rezone_env_logt1, &
           rotation_active, run_description, scale_be9_rate1, scale_be9_rate2, &
           scale_be9_rate3, scale_li6_rate, scale_li7_rate, &
           target_carbon_cno_fraction, target_metal_fraction, &
           target_nitrogen_cno_fraction, target_oxygen_cno_fraction, &
           use_alex06_tables, use_alex95_tables, use_diffusion_z, &
           use_kurucz90_tables, use_new_diffusion_routines, use_opal2001_eos, &
           use_opal2006_eos, use_opal95_eos, use_pure_z_table, use_thoul_diffusion, &
           use_thoul_fit, value_relative_to_h
      use eos_lib, only: opal_vintage_none, opal_vintage_1995, opal_vintage_2001, &
           opal_vintage_2006
      implicit none
      private
      public :: derive_controls

      integer, parameter :: message_length = 512

! The messages' formats: a species that is not one, the CNO mixture's
! and the isotope mixture's (rejected, unknown, applied), and the
! format of the run log's version banner.
      character(len=*), parameter :: &
           unknown_species_format = "('change_species_name (ANEWCP) ',a3," // &
           "' is not a species of the starting model; no species rescaled')", &
           negative_cno_format = "('negative CNO fraction (',3es12.4,'); mixture not modified')", &
           cno_sum_format = "('the CNO fractions (',3es12.4,') sum to 1 or more; mixture not modified')", &
           unknown_cno_format = "('CNO mixture ',a8,' not recognized; mixture not altered')", &
           cno_applied_format = "('CNO mixture ',a8,': C =',es12.4,'  N =',es12.4," // &
           "'  O =',es12.4,' (fractions of Z); applied to the starting model')", &
           negative_isotope_format = "('negative isotope ratio or light-element'," // &
           "' mass fraction (',11es12.4,'); mixture not modified')", &
           light_element_sum_format = "('light-element mass fractions sum above 1'," // &
           "' (',7es12.4,'); mixture not modified')", &
           unknown_isotope_format = "('isotope mixture ',a8,' not recognized; mixture not altered')", &
           isotope_applied_format = "('isotope/light-element mixture ',a8,': C12/C13 ='," // &
           "es12.4,'  O16/O18 =',es12.4,'  H2 =',es12.4,'  He3 =',es12.4," // &
           "'  Li6 =',es12.4,'  Li7 =',es12.4,'  Be9 =',es12.4,'; applied to the starting model')", &
           banner_format_format = "('(''# YREC v'', A', i2.2, ', '' ('', A', i2.2, ', '')'')')"

! The species change_species_name may name, from He3 on in
! composition-slot order (He3 is slot i_he3).
      character(len=3), parameter :: element_id(12) = (/ 'HE3', 'C12', 'C13', 'N14', &
           'N15', 'O16', 'O17', 'O18', 'H2 ', 'LI6', 'LI7', 'BE9' /)
! The CNO mixtures cno_mixture_name may name: Z/X and the C, N and O
! fractions of Z of GS98 (Grevesse & Sauval 1998), AAG21 (Asplund et al.
! 2021), M22P and M22M (Magg et al. 2022, photospheric and meteoritic).
      character(len=8), parameter :: mixture_id_table(4) = &
           (/ character(len=8) :: 'GS98', 'AAG21', 'M22P', 'M22M' /)
      real(dp), parameter :: zx_mix_table(4) = (/ 0.02292d0, 0.0187d0, 0.0225d0, 0.0226d0 /)
      real(dp), parameter :: frac_c_table(4) = (/ 0.172148d0, 0.184156d0, 0.19239d0, 0.191425d0 /)
      real(dp), parameter :: frac_n_table(4) = (/ 0.050426d0, 0.050344d0, 0.059012d0, 0.058716d0 /)
      real(dp), parameter :: frac_o_table(4) = (/ 0.468195d0, 0.416592d0, 0.415631d0, 0.413545d0 /)

contains

! ---------------------------------------------------------------
! The derivations, in order, and the run log's header. (controls_check's
! check_controls has rejected the controls the run cannot use.)
subroutine derive_controls
      call derive_rate_scales
      call select_eos_and_opacity_tables
      envelope_store_log_temperature = 0.5d0*(rezone_env_logt0 + rezone_env_logt1)
      if (rotation_active) call derive_rotation_controls
      if (rescale_species_active) call find_rescaled_species
      call select_cno_mixture
      call select_isotope_mixture
      call derive_diffusion_switches
      if (pulse_gyre_interval < 0) pulse_gyre_interval = 0
      call write_run_log_header
end subroutine derive_controls

! ---------------------------------------------------------------
! A light-element rate scale given as a cross section (keV b,
! li6_rate_value etc., with scale_li6_rate etc. on): the cross section
! over its reference value (controls_buffer's *_reference_xs, which are
! also the cross sections' defaults).
subroutine derive_rate_scales
      if (scale_li6_rate) li6_rate_scale = li6_rate_value/li6_reference_xs
      if (scale_li7_rate) li7_rate_scale = li7_rate_value/li7_reference_xs
      if (scale_be9_rate1) be9_pg_rate_scale = be9_rate_value1/be9_pg_reference_xs
      if (scale_be9_rate2) be9_pd_rate_scale = be9_rate_value2/be9_pd_reference_xs
      if (scale_be9_rate3) be9_palpha_rate_scale = be9_rate_value3/be9_palpha_reference_xs
end subroutine derive_rate_scales

! ---------------------------------------------------------------
! One OPAL EOS vintage and one family of molecular opacity tables. The
! newest OPAL EOS requested wins: opal_eos_vintage and its table path,
! opal_eos_table_path, which the eos settings builder copies. Of the
! molecular tables, Alexander 2006 wins over Alexander 1994/95, which wins
! over Kurucz 1990: the losers are switched off (kap reads one family, and
! rebuilds the Alexander 1994/95 surface table whenever that switch is on).
! Selecting more than one of either is a warning.
subroutine select_eos_and_opacity_tables
      character(len=message_length) :: line

      if (count([use_opal95_eos, use_opal2001_eos, use_opal2006_eos]) > 1) then
         write(line,'(a,3l2,a)') 'use_opal95_eos, use_opal2001_eos, use_opal2006_eos =', &
              use_opal95_eos, use_opal2001_eos, use_opal2006_eos, &
              ': more than one OPAL EOS; the newest is used'
         call report_warning(line)
      end if
      opal_eos_vintage = opal_vintage_none
      opal_eos_table_path = ""
      if (use_opal2006_eos) then
         opal_eos_vintage = opal_vintage_2006
         opal_eos_table_path = opal2006_eos_table_path
      else if (use_opal2001_eos) then
         opal_eos_vintage = opal_vintage_2001
         opal_eos_table_path = opal2001_eos_table_path
      else if (use_opal95_eos) then
         opal_eos_vintage = opal_vintage_1995
         opal_eos_table_path = opal95_eos_table_path
      end if
      if (count([use_alex06_tables, use_alex95_tables, use_kurucz90_tables]) > 1) then
         write(line,'(a,3l2,a)') 'use_alex06_tables, use_alex95_tables, use_kurucz90_tables =', &
              use_alex06_tables, use_alex95_tables, use_kurucz90_tables, &
              ': more than one molecular opacity family; the first is used'
         call report_warning(line)
      end if
      if (use_alex06_tables) then
         use_alex95_tables = .false.
         use_kurucz90_tables = .false.
      end if
      if (use_alex95_tables) use_kurucz90_tables = .false.
end subroutine select_eos_and_opacity_tables

! ---------------------------------------------------------------
! Rotating runs start with a fresh envelope-fit triangle. (The wind
! law's coefficients, the cz_rotation_law_exponent clamp and the solid-body flag are the
! rotation settings builder's, setup/rotation_settings_builder.f90.)
subroutine derive_rotation_controls
      new_triangle_init_flag = .true.
end subroutine derive_rotation_controls

! ---------------------------------------------------------------
! The species change_species_name names (an index into the composition
! array), and whether its new value is relative to hydrogen
! (new_species_value_mode /= 'ABS').
subroutine find_rescaled_species
      use composition_layout, only: i_he3
      integer :: i
      character(len=message_length) :: line

      value_relative_to_h = .true.
      if (new_species_value_mode == 'ABS') value_relative_to_h = .false.
      do i = 1, size(element_id)
         if (change_species_name == element_id(i)) then
            new_species_index = i_he3 - 1 + i
            exit
         end if
      end do
      if (i > size(element_id)) then
         rescale_species_active = .false.
         write(line,unknown_species_format) change_species_name
         call report_warning(line)
      end if
end subroutine find_rescaled_species

! ---------------------------------------------------------------
! mixture_change_mode = 1: the C, N and O fractions of Z and Z/X from
! the cno_mixture_name table row, or as given for 'CUS' unless negative
! or summing to 1 or more.
subroutine select_cno_mixture
      real(dp) :: sum_frac
      integer :: i
      character(len=message_length) :: line

      change_cno_mixture_active = .false.
      if (mixture_change_mode == 1) then
         if (cno_mixture_name == 'CUS') then
            if (target_carbon_cno_fraction < 0.0d0 .or. target_nitrogen_cno_fraction < 0.0d0 &
                 .or. target_oxygen_cno_fraction < 0.0d0) then
               write(line,negative_cno_format) target_carbon_cno_fraction, &
                    target_nitrogen_cno_fraction, target_oxygen_cno_fraction
               call report_warning(line)
            else
               sum_frac = target_carbon_cno_fraction + target_nitrogen_cno_fraction &
                    + target_oxygen_cno_fraction
               if (sum_frac >= 1.0d0) then
                  write(line,cno_sum_format) target_carbon_cno_fraction, &
                       target_nitrogen_cno_fraction, target_oxygen_cno_fraction
                  call report_warning(line)
               else
                  change_cno_mixture_active = .true.
               end if
            end if
         else
            do i = 1, 4
               if (cno_mixture_name == mixture_id_table(i)) then
                  target_metal_fraction = zx_mix_table(i)
                  target_carbon_cno_fraction = frac_c_table(i)
                  target_nitrogen_cno_fraction = frac_n_table(i)
                  target_oxygen_cno_fraction = frac_o_table(i)
                  change_cno_mixture_active = .true.
                  exit
               end if
            end do
            if (.not. change_cno_mixture_active) then
               write(line,unknown_cno_format) cno_mixture_name
               call report_warning(line)
            end if
         end if
      end if
      if (change_cno_mixture_active) then
         write(line,cno_applied_format) cno_mixture_name, target_carbon_cno_fraction, &
              target_nitrogen_cno_fraction, target_oxygen_cno_fraction
         call report_note(line)
      end if
end subroutine select_cno_mixture

! ---------------------------------------------------------------
! isotope_change_mode = 1: the CNO isotope ratios and light-element
! abundances, the defaults for isotope_mixture_name 'L21' or as given for
! 'CUS' unless negative or summing to 1 or more.
subroutine select_isotope_mixture
      real(dp) :: sum_frac
      character(len=message_length) :: line

      change_isotope_ratios_active = .false.
      if (isotope_change_mode == 1) then
         if (isotope_mixture_name == 'CUS') then
            if (c12_to_c13_ratio < 0.0d0 .or. n14_to_n15_ratio < 0.0d0 &
                 .or. o16_to_o17_ratio < 0.0d0 .or. o16_to_o18_ratio < 0.0d0 &
                 .or. initial_h2_fraction < 0.0d0 .or. initial_he3_fraction < 0.0d0 &
                 .or. initial_li6_fraction < 0.0d0 .or. initial_li7_fraction < 0.0d0 &
                 .or. initial_be9_fraction < 0.0d0 .or. initial_b10_fraction < 0.0d0 &
                 .or. initial_b11_fraction < 0.0d0) then
               write(line,negative_isotope_format) c12_to_c13_ratio, n14_to_n15_ratio, o16_to_o17_ratio, &
                    o16_to_o18_ratio, initial_h2_fraction, initial_he3_fraction, &
                    initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                    initial_b10_fraction, initial_b11_fraction
               call report_warning(line)
            else
               sum_frac = initial_h2_fraction + initial_he3_fraction + initial_li6_fraction &
                    + initial_li7_fraction + initial_be9_fraction + initial_b10_fraction &
                    + initial_b11_fraction
               if (sum_frac >= 1.0d0) then
                  write(line,light_element_sum_format) initial_h2_fraction, initial_he3_fraction, &
                       initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                       initial_b10_fraction, initial_b11_fraction
                  call report_warning(line)
               else
                  change_isotope_ratios_active = .true.
               end if
            end if
         else if (isotope_mixture_name == 'L21') then
            change_isotope_ratios_active = .true.
         else
            write(line,unknown_isotope_format) isotope_mixture_name
            call report_warning(line)
         end if
      end if
      if (change_isotope_ratios_active) then
         write(line,isotope_applied_format) isotope_mixture_name, c12_to_c13_ratio, &
              o16_to_o18_ratio, initial_h2_fraction, initial_he3_fraction, &
              initial_li6_fraction, initial_li7_fraction, initial_be9_fraction
         call report_note(line)
      end if
end subroutine select_isotope_mixture

! ---------------------------------------------------------------
! Metal diffusion needs helium diffusion and the Thoul et al. routines;
! lithium diffusion needs all of the new diffusion routines.
subroutine derive_diffusion_switches
      if (use_diffusion_z) then
         diffuse_helium_active = .true.
         use_thoul_diffusion = .true.
      end if
      if (diffuse_lithium) then
         use_new_diffusion_routines = .true.
         diffuse_helium_active = .true.
         use_diffusion_z = .true.
         use_thoul_diffusion = .true.
         use_thoul_fit = .false.
         coulomb_log_choice = 4
      end if
end subroutine derive_diffusion_switches

! ---------------------------------------------------------------
! The surface boundary condition, the pure C and N opacity tables, and
! the version and description banner.
subroutine write_run_log_header
! the banner's format, fitted to the version strings (a zero width is
! illegal)
      character(len=64) :: version_fmt
      character(len=message_length) :: line

      write(version_fmt, banner_format_format) max(1, len_trim(yrec_version_string)), &
           max(1, len_trim(git_hash_string))
      select case (atm_choice)
      case (0)
         call log_line(' surface boundary: Eddington gray T(tau) relation')
      case (1)
         call log_line(' surface boundary: Krishna-Swamy T(tau) relation')
      case (2)
         call log_line(' surface boundary: Harvard-Smithsonian reference atmosphere')
      case (3)
         call log_line(' surface boundary: Kurucz atmosphere tables')
      case (4)
         call log_line(' surface boundary: Allard atmosphere tables')
      case (5)
         call log_line(' surface boundary: Kurucz/Castelli atmosphere tables')
      end select
      if (use_pure_z_table) call log_line(' opacity: pure C and N tables enabled')
      call log_line('')
      call log_line(' ' // repeat('=', 100))
      write(line, version_fmt) yrec_version_string, git_hash_string
      call log_line(trim(line))
      call log_line(' description: ' // trim(run_description(1)))
      if (len_trim(run_description(2)) > 0) call log_line(repeat(' ', 14) // trim(run_description(2)))
      call log_line(' ' // repeat('=', 100))
end subroutine write_run_log_header

end module derived_controls
