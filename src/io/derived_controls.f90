!----------------------------------------------------------------------
! derived_controls
!----------------------------------------------------------------------
! The controls computed from the ones read: the light-element rate
! scales, one OPAL EOS and one low-temperature opacity source, the wind
! law's coefficients and exponents, the starting model's mixture
! changes, and the switches that diffusion options imply. It also
! rejects semiconvection with overshoot and writes the run log's
! header. Everything here works on the canonical names, after
! adopt_canonical_names; the legacy names it reads (alfa, fk, wmax_sun,
! amix, aiso, anewcp, atmp, descrip, the Li/Be inputs) have no
! canonical home.
module derived_controls
      use star_info_lib, only: star
      use physical_constants
      use io_units
      use legacy_controls
      implicit none
      private
      public :: derive_controls

contains

! ---------------------------------------------------------------
! ierr = 1: semiconvection and overshoot are both on.
subroutine derive_controls(ierr)
      integer, intent(inout) :: ierr

      call derive_rate_scales
      call select_eos_and_opacity_tables
      if (use_semiconvection) then
         if (core_overshoot_active .or. envelope_overshoot_active .or. lovstm) then
            write(run_log_unit,2) use_semiconvection, envelope_overshoot_active, &
                 core_overshoot_active, lovstm
      2       format(1x,'ERROR IN SUBROUTINE PARMIN'/'SEMI-CONVECTION', &
           &  ' AND OVERSHOOT FLAGS BOTH TURNED ON'/'FLAGS LSEMIC',l2, &
           &  ' OVERSHOOT - CORE,ENVELOPE,INTERMEDIATE-',3l2/'RUN STOPPED')
            ierr = 1
            return
         end if
      end if
      star%tenv = 0.5d0*(tenv0 + tenv1)
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
! The light-element cross sections (keV b) as scales of their defaults
! (Li6 5500, Li7 52, Be9(p,g) 1.1, Be9(p,d) and Be9(p,a) 15000; see
! set_legacy_defaults).
subroutine derive_rate_scales
      if (lxli6) sli6 = xsli6/5.5d3
      if (lxli7) sli7 = xsli7/5.2d1
      if (lxbe91) sbe91 = xsbe91/1.1d0
      if (lxbe92) sbe92 = xsbe92/1.5d4
      if (lxbe93) sbe93 = xsbe93/1.5d4
      li6_rate_scale = sli6
      li7_rate_scale = sli7
      be9_pg_rate_scale = sbe91
      be9_pd_rate_scale = sbe92
      be9_palpha_rate_scale = sbe93
end subroutine derive_rate_scales

! ---------------------------------------------------------------
! The newest OPAL EOS and Alexander opacities requested win; eos_init
! opens star%job%opal_eos_table_path.
subroutine select_eos_and_opacity_tables
      star%job%opal_eos_table_path = ""
      if (use_opal2006_eos) then
         use_opal2001_eos = .false.
         use_opal95_eos = .false.
         star%job%opal_eos_table_path = fopale06
      end if
      if (use_opal2001_eos) then
         star%job%opal_eos_table_path = fopale01
         use_opal95_eos = .false.
      else if (use_opal95_eos) then
         star%job%opal_eos_table_path = fopale
      end if
      if (use_alex06_tables) then
         use_alex95_tables = .false.
         use_kurucz90_tables = .false.
      end if
      if (use_alex95_tables) use_kurucz90_tables = .false.
end subroutine select_eos_and_opacity_tables

! ---------------------------------------------------------------
! Rotating runs: a fresh envelope-fit triangle, walpcz clamped to
! [-2, 0], solid-body rotation, and the wind law's coefficients.
subroutine derive_rotation_controls
      use math_lib
      double precision :: one_third, two_thirds

      lnew0 = .true.
      if (walpcz < -2.0d0) walpcz = -2.0d0
      if (walpcz > 0.0d0) walpcz = 0.0d0
      if (force_solid_body_rotation) solid_body_mode_flag = 1
      if (.not. use_pmm_wind_law) then
! The Kawaler (1988) law. For an index alfa,
!   Jdot = fk 2.036e33 1.452e9**alfa (Mdot/1e-14)**(1-2alfa/3)
!          Omega**(1+4alfa/3) (R/Rsun)**(2-alfa) (M/Msun)**(-alfa/3);
! exmd, wind_law_omega_exponent, exr and exm are the exponents of Mdot,
! Omega, R and M. fk = 1 reproduces the observed solar angular momentum
! loss.
         one_third = 1.0d0/3.0d0
         two_thirds = 2.0d0/3.0d0
         constfactor = fk*2.036d33*pow(1.452d9, alfa)
         exmd = 1.0d0 - two_thirds*alfa
         wind_law_omega_exponent = 1.0d0 + 2.0d0*two_thirds*alfa
         exr = 2.0d0 - alfa
         exm = -one_third*alfa
      else
! The PMM law: wind_law_name 'K97' is Krishnamurthi et al. (1997), 'V13'
! van Saders & Pinsonneault (2013, after Matt et al. 2008, 2012), and
! 'CUS' keeps the pmm_exponent_* given. With every term scaled to solar,
!   Jdot = w Bmag**4m Mdot**(1-2m) R**(5m+2) M**(-m) Fcen**md
!   Bmag = Pphot**0.5 (w R**c)**b Tcz**d
!   Mdot = R M**-1 Lbol w**a  (from GM/R Mdot = Lx, Lx/Lbol ~ w**a)
! so Jdot ~ w**(1+a-2ma+4mb) Tcz**(4md) Lbol**(1-2m) Pphot**(2m) Fcen**md,
! and, ignoring the GM/R term in Mdot, R**(2+5m+4mbc) M**(-m). excen is
! the exponent of the centrifugal term {K2**2/(K2**2 + 0.5 w**2 R**3/GM)};
! d switches the Tcz and centrifugal terms off. The Tcz term is
! computed in the braking itself, so extau is 0.
         if (wind_law_name == 'K97') then
            pmm_exponent_m = alfa/3.0
            pmm_exponent_a = 2.0
            pmm_exponent_b = 1.0
            pmm_exponent_c = 2.0
            pmm_exponent_d = 0.0
            scale_by_b_field = .false.
         else if (wind_law_name == 'V13') then
            pmm_exponent_m = 0.22
            pmm_exponent_a = 2.0
            pmm_exponent_b = 1.0
            pmm_exponent_c = 0.0
            pmm_exponent_d = 1.0
            scale_by_b_field = .true.
         end if
         associate (a => pmm_exponent_a, b => pmm_exponent_b, c => pmm_exponent_c, &
                    d => pmm_exponent_d, m => pmm_exponent_m)
         wind_law_omega_exponent = 1.0d0 + a - 2.0d0*a*m + 4.0d0*m*b
         extau = 0.0
         exr = 2.0d0 + 5.0d0*m - 4.0d0*m*c
         exm = -m
         exl = 1.0d0 - 2.0d0*m
         expr = 2.0d0*m*d
         excen = m*d
         end associate
         c_2 = 0.0506
         constfactor = fk*pmm_norm_jdot/pow(pmm_solar_omega, wind_law_omega_exponent)
! wmax_sun < 1000 gives the saturation rate in units of the solar rate
         if (wmax_sun < 1.0e3) wind_saturation_omega = wmax_sun*pmm_solar_omega
      end if
end subroutine derive_rotation_controls

! ---------------------------------------------------------------
! The species anewcp names (an index into the composition array), and
! whether its new value is relative to hydrogen (atmp /= 'ABS').
subroutine find_rescaled_species
      integer :: i

      value_relative_to_h = .true.
      if (atmp == 'ABS') value_relative_to_h = .false.
      do i = 1, 12
         if (anewcp == element_id(i)) then
            new_species_index = i + 3
            exit
         end if
      end do
      if (i > 12) then
         rescale_species_active = .false.
         write(run_log_unit,20) anewcp
      20    format(1x,'VARIABLE',a4,1x,'NOT A RECOGNIZED ELEMENT'/1x, &
           &    'RESCALING NOT PERFORMED')
      end if
end subroutine find_rescaled_species

! ---------------------------------------------------------------
! mixture_change_mode = 1: the C, N and O fractions of Z and Z/X from
! the amix table row, or as given for 'CUS' unless negative or summing
! to 1 or more.
subroutine select_cno_mixture
      double precision :: sum_frac
      integer :: i

      change_cno_mixture_active = .false.
      if (mixture_change_mode == 1) then
         if (amix == 'CUS') then
            if (target_carbon_cno_fraction < 0.0d0 .or. target_nitrogen_cno_fraction < 0.0d0 &
                 .or. target_oxygen_cno_fraction < 0.0d0) then
               write(*,591) target_carbon_cno_fraction, target_nitrogen_cno_fraction, &
                    target_oxygen_cno_fraction
               write(run_log_unit,591) target_carbon_cno_fraction, target_nitrogen_cno_fraction, &
                    target_oxygen_cno_fraction
      591          format('NEGATIVE INPUT CNO FRACTION ',3e12.4, &
           &              ' MIX NOT MODIFIED')
            else
               sum_frac = target_carbon_cno_fraction + target_nitrogen_cno_fraction &
                    + target_oxygen_cno_fraction
               if (sum_frac >= 1.0d0) then
                  write(*,598) target_carbon_cno_fraction, target_nitrogen_cno_fraction, &
                       target_oxygen_cno_fraction
                  write(run_log_unit,598) target_carbon_cno_fraction, &
                       target_nitrogen_cno_fraction, target_oxygen_cno_fraction
      598             format('INPUT CNO FRACTION ',3e12.4, &
           &                 ' EXCEEDS 1. MIX NOT MODIFIED')
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
               write(*,589) amix
               write(run_log_unit,589) amix
      589       format(1x,'warning: CNO mixture ',a8,' not recognized; mixture not altered')
            end if
         end if
      end if
      if (change_cno_mixture_active) then
         write(*,604) amix, target_carbon_cno_fraction, target_nitrogen_cno_fraction, &
              target_oxygen_cno_fraction
         write(run_log_unit,604) amix, target_carbon_cno_fraction, &
              target_nitrogen_cno_fraction, target_oxygen_cno_fraction
      604    format(1x,'CNO mixture ',a8,': C =',es12.4,'  N =',es12.4, &
           &         '  O =',es12.4,' (fractions of Z); applied to the starting model')
      end if
end subroutine select_cno_mixture

! ---------------------------------------------------------------
! isotope_change_mode = 1: the CNO isotope ratios and light-element
! abundances, the defaults for aiso 'L21' or as given for 'CUS' unless
! negative or summing to 1 or more.
subroutine select_isotope_mixture
      double precision :: sum_frac

      change_isotope_ratios_active = .false.
      if (isotope_change_mode == 1) then
         if (aiso == 'CUS') then
            if (c12_to_c13_ratio < 0.0d0 .or. n14_to_n15_ratio < 0.0d0 &
                 .or. o16_to_o17_ratio < 0.0d0 .or. o16_to_o18_ratio < 0.0d0 &
                 .or. initial_h2_fraction < 0.0d0 .or. initial_he3_fraction < 0.0d0 &
                 .or. initial_li6_fraction < 0.0d0 .or. initial_li7_fraction < 0.0d0 &
                 .or. initial_be9_fraction < 0.0d0 .or. initial_b10_fraction < 0.0d0 &
                 .or. initial_b11_fraction < 0.0d0) then
               write(*,596) c12_to_c13_ratio, n14_to_n15_ratio, o16_to_o17_ratio, &
                    o16_to_o18_ratio, initial_h2_fraction, initial_he3_fraction, &
                    initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                    initial_b10_fraction, initial_b11_fraction
               write(run_log_unit,596) c12_to_c13_ratio, n14_to_n15_ratio, o16_to_o17_ratio, &
                    o16_to_o18_ratio, initial_h2_fraction, initial_he3_fraction, &
                    initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                    initial_b10_fraction, initial_b11_fraction
      596          format(1x,'warning: negative isotope ratio or light-element', &
           &   ' mass fraction (',11es12.4,'); mixture not modified')
            else
               sum_frac = initial_h2_fraction + initial_he3_fraction + initial_li6_fraction &
                    + initial_li7_fraction + initial_be9_fraction + initial_b10_fraction &
                    + initial_b11_fraction
               if (sum_frac >= 1.0d0) then
                  write(*,595) initial_h2_fraction, initial_he3_fraction, &
                       initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                       initial_b10_fraction, initial_b11_fraction
                  write(run_log_unit,595) initial_h2_fraction, initial_he3_fraction, &
                       initial_li6_fraction, initial_li7_fraction, initial_be9_fraction, &
                       initial_b10_fraction, initial_b11_fraction
      595             format(1x,'warning: light-element mass fractions sum above 1', &
           &     ' (',11es12.4,'); mixture not modified')
               else
                  change_isotope_ratios_active = .true.
               end if
            end if
         else if (aiso == 'L21') then
            change_isotope_ratios_active = .true.
         end if
      end if
      if (change_isotope_ratios_active) then
         write(*,605) aiso, c12_to_c13_ratio, o16_to_o18_ratio, initial_h2_fraction, &
              initial_he3_fraction, initial_li6_fraction, initial_li7_fraction, &
              initial_be9_fraction
         write(run_log_unit,605) aiso, c12_to_c13_ratio, o16_to_o18_ratio, &
              initial_h2_fraction, initial_he3_fraction, initial_li6_fraction, &
              initial_li7_fraction, initial_be9_fraction
      605    format(1x,'isotope/light-element mixture ',a8,': C12/C13 =', &
           &    es12.4,'  O16/O18 =',es12.4,'  H2 =',es12.4,'  He3 =',es12.4, &
           &    '  Li6 =',es12.4,'  Li7 =',es12.4,'  Be9 =',es12.4, &
           &    '; applied to the starting model')
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
