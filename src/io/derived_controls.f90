!----------------------------------------------------------------------
! derived_controls
!----------------------------------------------------------------------
! The controls computed from the ones read: the light-element rate
! scales, one OPAL EOS and one low-temperature opacity source, a fresh
! envelope triangle for rotating runs, the starting model's mixture
! changes, and the switches that diffusion options imply; and it writes
! the run log's header. Its messages go to the terminal and the run log
! through io_units (report_warning, report_note); the errors are
! controls_check's. Everything here works on the canonical names, after
! adopt_canonical_names; the legacy names it reads (amix, aiso,
! anewcp, atmp, descrip, the Li/Be inputs) have no canonical home.
module derived_controls
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use yrec_version, only: yrec_version_string, git_hash_string
      use io_units, only: run_log_unit, report_warning, report_note
      use legacy_controls
      use eos_lib, only: opal_vintage_none, opal_vintage_1995, opal_vintage_2001, &
           opal_vintage_2006
      use star_info_lib, only: star_info
      implicit none
      private
      public :: derive_controls

      integer, parameter :: message_length = 512

contains

! ---------------------------------------------------------------
! The derivations, in order, and the run log's header. (controls_check's
! check_controls has rejected the controls the run cannot use.)
subroutine derive_controls(s)
      type(star_info), intent(inout) :: s

      call derive_rate_scales
      call select_eos_and_opacity_tables
      s%envelope_store_log_temperature = 0.5d0*(tenv0 + tenv1)
      if (rotation_active) call derive_rotation_controls
      if (rescale_species_active) call find_rescaled_species
      call select_cno_mixture
      call select_isotope_mixture
      call derive_diffusion_switches
      if (npenv <= 0) npenv = 9999
      if (pulse_gyre_interval < 0) pulse_gyre_interval = 0
      call write_run_log_header
end subroutine derive_controls

! ---------------------------------------------------------------
! The light-element cross sections (keV b, the xs* controls) as scales
! of their reference values (legacy_controls' *_reference_xs), which are
! also their defaults.
subroutine derive_rate_scales
      if (lxli6) sli6 = xsli6/li6_reference_xs
      if (lxli7) sli7 = xsli7/li7_reference_xs
      if (lxbe91) sbe91 = xsbe91/be9_pg_reference_xs
      if (lxbe92) sbe92 = xsbe92/be9_pd_reference_xs
      if (lxbe93) sbe93 = xsbe93/be9_palpha_reference_xs
      li6_rate_scale = sli6
      li7_rate_scale = sli7
      be9_pg_rate_scale = sbe91
      be9_pd_rate_scale = sbe92
      be9_palpha_rate_scale = sbe93
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
! law's coefficients, the walpcz clamp and the solid-body flag are the
! rotation settings builder's, setup/rotation_settings_setup.f90.)
subroutine derive_rotation_controls
      lnew0 = .true.
end subroutine derive_rotation_controls

! ---------------------------------------------------------------
! The species anewcp names (an index into the composition array), and
! whether its new value is relative to hydrogen (atmp /= 'ABS').
! element_id lists the species from He3 on, in composition-slot order
! (He3 is slot i_he3).
subroutine find_rescaled_species
      use composition_layout, only: i_he3
      integer :: i
      character(len=message_length) :: line

      value_relative_to_h = .true.
      if (atmp == 'ABS') value_relative_to_h = .false.
      do i = 1, size(element_id)
         if (anewcp == element_id(i)) then
            new_species_index = i_he3 - 1 + i
            exit
         end if
      end do
      if (i > size(element_id)) then
         rescale_species_active = .false.
         write(line,20) anewcp
      20    format('change_species_name (ANEWCP) ',a3, &
           &    ' is not a species of the starting model; no species rescaled')
         call report_warning(line)
      end if
end subroutine find_rescaled_species

! ---------------------------------------------------------------
! mixture_change_mode = 1: the C, N and O fractions of Z and Z/X from
! the amix table row, or as given for 'CUS' unless negative or summing
! to 1 or more.
subroutine select_cno_mixture
      real(dp) :: sum_frac
      integer :: i
      character(len=message_length) :: line

      change_cno_mixture_active = .false.
      if (mixture_change_mode == 1) then
         if (amix == 'CUS') then
            if (target_carbon_cno_fraction < 0.0d0 .or. target_nitrogen_cno_fraction < 0.0d0 &
                 .or. target_oxygen_cno_fraction < 0.0d0) then
               write(line,591) target_carbon_cno_fraction, target_nitrogen_cno_fraction, &
                    target_oxygen_cno_fraction
      591          format('negative CNO fraction (',3es12.4, &
           &              '); mixture not modified')
               call report_warning(line)
            else
               sum_frac = target_carbon_cno_fraction + target_nitrogen_cno_fraction &
                    + target_oxygen_cno_fraction
               if (sum_frac >= 1.0d0) then
                  write(line,598) target_carbon_cno_fraction, target_nitrogen_cno_fraction, &
                       target_oxygen_cno_fraction
      598             format('the CNO fractions (',3es12.4, &
           &                 ') sum to 1 or more; mixture not modified')
                  call report_warning(line)
               else
                  change_cno_mixture_active = .true.
               end if
            end if
         else
            do i = 1, 4
               if (amix == mixture_id_table(i)) then
                  target_metal_fraction = zx_mix_table(i)
                  target_carbon_cno_fraction = frac_c_table(i)
                  target_nitrogen_cno_fraction = frac_n_table(i)
                  target_oxygen_cno_fraction = frac_o_table(i)
                  change_cno_mixture_active = .true.
                  exit
               end if
            end do
            if (.not. change_cno_mixture_active) then
               write(line,589) amix
      589       format('CNO mixture ',a8,' not recognized; mixture not altered')
               call report_warning(line)
            end if
         end if
      end if
      if (change_cno_mixture_active) then
         write(line,604) amix, target_carbon_cno_fraction, target_nitrogen_cno_fraction, &
              target_oxygen_cno_fraction
      604    format('CNO mixture ',a8,': C =',es12.4,'  N =',es12.4, &
           &         '  O =',es12.4,' (fractions of Z); applied to the starting model')
         call report_note(line)
      end if
end subroutine select_cno_mixture

! ---------------------------------------------------------------
! isotope_change_mode = 1: the CNO isotope ratios and light-element
! abundances, the defaults for aiso 'L21' or as given for 'CUS' unless
! negative or summing to 1 or more.
subroutine select_isotope_mixture
      real(dp) :: sum_frac
      character(len=message_length) :: line

      change_isotope_ratios_active = .false.
      if (isotope_change_mode == 1) then
         if (aiso == 'CUS') then
            if (c12_to_c13_ratio < 0.0d0 .or. n14_to_n15_ratio < 0.0d0 &
                 .or. o16_to_o17_ratio < 0.0d0 .or. o16_to_o18_ratio < 0.0d0 &
                 .or. initial_h2_fraction < 0.0d0 .or. initial_he3_fraction < 0.0d0 &
                 .or. initial_li6_fraction < 0.0d0 .or. initial_li7_fraction < 0.0d0 &
                 .or. initial_be9_fraction < 0.0d0 .or. initial_b10_fraction < 0.0d0 &
                 .or. initial_b11_fraction < 0.0d0) then
               write(line,596) c12_to_c13_ratio, n14_to_n15_ratio, o16_to_o17_ratio, &
                    o16_to_o18_ratio, initial_h2_fraction, initial_he3_fraction, &
                    initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                    initial_b10_fraction, initial_b11_fraction
      596          format('negative isotope ratio or light-element', &
           &   ' mass fraction (',11es12.4,'); mixture not modified')
               call report_warning(line)
            else
               sum_frac = initial_h2_fraction + initial_he3_fraction + initial_li6_fraction &
                    + initial_li7_fraction + initial_be9_fraction + initial_b10_fraction &
                    + initial_b11_fraction
               if (sum_frac >= 1.0d0) then
                  write(line,595) initial_h2_fraction, initial_he3_fraction, &
                       initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                       initial_b10_fraction, initial_b11_fraction
      595             format('light-element mass fractions sum above 1', &
           &     ' (',7es12.4,'); mixture not modified')
                  call report_warning(line)
               else
                  change_isotope_ratios_active = .true.
               end if
            end if
         else if (aiso == 'L21') then
            change_isotope_ratios_active = .true.
         else
            write(line,597) aiso
      597    format('isotope mixture ',a8,' not recognized; mixture not altered')
            call report_warning(line)
         end if
      end if
      if (change_isotope_ratios_active) then
         write(line,605) aiso, c12_to_c13_ratio, o16_to_o18_ratio, initial_h2_fraction, &
              initial_he3_fraction, initial_li6_fraction, initial_li7_fraction, &
              initial_be9_fraction
      605    format('isotope/light-element mixture ',a8,': C12/C13 =', &
           &    es12.4,'  O16/O18 =',es12.4,'  H2 =',es12.4,'  He3 =',es12.4, &
           &    '  Li6 =',es12.4,'  Li7 =',es12.4,'  Be9 =',es12.4, &
           &    '; applied to the starting model')
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

      write(version_fmt, 315) max(1, len_trim(yrec_version_string)), &
           max(1, len_trim(git_hash_string))
      315 format('(''# YREC v'', A', i2.2, ', '' ('', A', i2.2, &
           &        ', '')'')')
      select case (atm_choice)
      case (0)
         write(run_log_unit,'(1x,a)') 'surface boundary: Eddington gray T(tau) relation'
      case (1)
         write(run_log_unit,'(1x,a)') 'surface boundary: Krishna-Swamy T(tau) relation'
      case (2)
         write(run_log_unit,'(1x,a)') 'surface boundary: Harvard-Smithsonian reference atmosphere'
      case (3)
         write(run_log_unit,'(1x,a)') 'surface boundary: Kurucz atmosphere tables'
      case (4)
         write(run_log_unit,'(1x,a)') 'surface boundary: Allard atmosphere tables'
      case (5)
         write(run_log_unit,'(1x,a)') 'surface boundary: Kurucz/Castelli atmosphere tables'
      end select
      if (use_pure_z_table) &
           write(run_log_unit,'(1x,a)') 'opacity: pure C and N tables enabled'
      write(run_log_unit,'(/,1x,100(''=''))')
      write(run_log_unit,version_fmt) yrec_version_string, git_hash_string
      write(run_log_unit,'(1x,2a)') 'description: ', trim(descrip(1))
      if (len_trim(descrip(2)) > 0) &
           write(run_log_unit,'(14x,a)') trim(descrip(2))
      write(run_log_unit,'(1x,100(''=''))')
end subroutine write_run_log_header

end module derived_controls
