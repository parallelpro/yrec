!----------------------------------------------------------------------
! neutrino_flux_table
!----------------------------------------------------------------------
! Resurrected 2026 as an opt-in diagnostic (controls_lib's
! compute_neutrino_fluxes, default .false.). This is the former
! LNUTAB neutrino-table block of program main (MHP 2/04), which had
! been hardcoded off since 2004 and had bit-rotted while dead: the
! engeb call passed uninitialized locals in the deuterium-fraction
! and shell-index slots, and the block overwrote star%m/star%dm with
! its own shell-center weights. Fixed on resurrection:
!   - engeb receives star%xa(i_h2,k) and the loop index;
!   - the integration weights are LOCAL arrays (shell_center_mass,
!     shell_mass), star%m/star%dm are not touched;
!   - output goes to short_file_unit (the .short / CASE.log stream)
!     instead of the historical bare unit 76 (a stray fort.76 file)
!     and stdout.
! Called from run_yrec once per kind card, right after starin, so the
! table describes the starting model of each run.
!
! Content: per-shell neutrino production from engeb (erg/g/s scaled to
! erg/s by the shell mass), the totals for the ten flux slots (units
! of 1e10; slots 1-8 = pp, pep, hep, Be7, B8, N13, O15, F17), then the
! normalized production profiles in the Bahcall & Pinsonneault layout
! (r/Rsun, T6, log n_e, shell mass fraction, X_Be7, then pp, B8, N13,
! O15, F17, Be7, pep, hep).
subroutine neutrino_flux_table
      use star_info_lib, only: star, json, i_h1, i_h2, i_he3, i_he4, &
           i_c12, i_c13, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18, &
           i_nu_pp, i_nu_pep, i_nu_hep, i_nu_be7, i_nu_b8, i_nu_n13, &
           i_nu_o15, i_nu_f17
      use luout_lib
      use const_lib
      use burn_lib
      implicit none

      double precision :: shell_center_mass(json), shell_mass(json)
      double precision :: prev_mass_bound, curr_mass_bound, next_mass_bound
      double precision :: shell_log_density, shell_log_temperature, &
           hydrogen_fraction, helium_fraction, he3_fraction, &
           c12_fraction, c13_fraction, n14_fraction, &
           o16_fraction, o18_fraction
      double precision :: pp_chain_energy_gen, he3he4_be7_electron_energy_gen, &
           he3he4_be7_proton_energy_gen, cno_cycle_energy_gen, &
           triple_alpha_energy_gen, dlnepsilon_dlnrho, dlnepsilon_dlnt, &
           total_energy_gen_rate
      double precision :: t6_million_k, log_electron_density, &
           zone_mass_fraction, zone_radius_fraction
      integer :: i, j

! SET UP WEIGHTS AND MASSES (local; the star%m/star%dm model arrays
! are left untouched).
! shell_center_mass = LOCATION IN GM (UNLOGGED) OF SHELL CENTERS.
! shell_mass = MASS IN GM OF EACH SHELL.
      curr_mass_bound = exp(ln10*star%log_mass(1))
      prev_mass_bound = -curr_mass_bound
      do i = 2,star%nz
         next_mass_bound = prev_mass_bound
         prev_mass_bound = curr_mass_bound
         curr_mass_bound = exp(ln10*star%log_mass(i))
         shell_center_mass(i-1) = prev_mass_bound
         shell_mass(i-1) = 0.5d0*(curr_mass_bound-next_mass_bound)
      end do
      shell_center_mass(star%nz) = curr_mass_bound
      shell_mass(star%nz) = exp(ln10*star%log_total_mass) &
           - 0.5d0*(prev_mass_bound+curr_mass_bound)

      do j = 1,10
         star%flux%neutrino_flux_total(j) = 0.0d0
         do i = 1,star%nz
            star%neutrino_flux_zone(j,i) = 0.0d0
         end do
      end do

! PER-SHELL PRODUCTION FROM ENGEB.
      do i = 1,star%nz
         shell_log_density = star%logRho(i)
         shell_log_temperature = star%logT(i)
! SKIP CALCULATIONS FOR LOW TEMPERATURES.
         if (shell_log_temperature.lt.6.0d0) exit
         hydrogen_fraction = star%xa(i_h1,i)
         helium_fraction = star%xa(i_he4,i)
         he3_fraction = star%xa(i_he3,i)
         c12_fraction = star%xa(i_c12,i)
         c13_fraction = star%xa(i_c13,i)
         n14_fraction = star%xa(i_n14,i)
         o16_fraction = star%xa(i_o16,i)
         o18_fraction = star%xa(i_o18,i)
         call engeb(pp_chain_energy_gen,he3he4_be7_electron_energy_gen, &
              he3he4_be7_proton_energy_gen,cno_cycle_energy_gen, &
              triple_alpha_energy_gen,dlnepsilon_dlnrho,dlnepsilon_dlnt, &
              total_energy_gen_rate,shell_log_density, &
              shell_log_temperature,hydrogen_fraction,helium_fraction, &
              he3_fraction,c12_fraction,c13_fraction,n14_fraction, &
              o16_fraction,o18_fraction,star%xa(i_h2,i),i, &
              star%reaction_rate_1,star%reaction_rate_2, &
              star%reaction_rate_3,star%reaction_rate_4, &
              star%reaction_rate_5,star%reaction_rate_6, &
              star%reaction_rate_7,star%reaction_rate_8, &
              star%reaction_rate_9,star%reaction_rate_10, &
              star%reaction_rate_11,star%reaction_rate_12, &
              star%reaction_rate_13,star%n15_alpha_branch_fraction, &
              star%be7_electron_capture_fraction)
! BE7 MASS FRACTION.
         star%be7_mass_fraction_zone(i) = star%engeb%be7_mass_fraction
! CONVERT FROM ERG/GM/S TO ERG/S FOR EACH SHELL BY MULTIPLYING
! BY THE MASS OF EACH SHELL IN GM.
         do j = 1,10
            star%neutrino_flux_zone(j,i) = star%flux%neutrino_flux(j)*shell_mass(i)
            star%flux%neutrino_flux_total(j) = star%flux%neutrino_flux_total(j) &
                 + star%neutrino_flux_zone(j,i)
         end do
         write(short_file_unit,911)i,shell_mass(i), &
              (star%neutrino_flux_zone(j,i),j=1,10)
 911     format(I5,1P11E10.3)
      end do

! WRITE OUT TOTAL NEUTRINO FLUXES.
! ***NOTE THAT THESE ARE IN UNITS OF 10**10. ***
      write(short_file_unit,222)(star%flux%neutrino_flux_total(i),i=1,10)
 222  format(1P10E10.3)

! NORMALIZE FLUXES.
      do j = 1,10
         if (star%flux%neutrino_flux_total(j) .ne. 0.0d0) then
            do i = 1,star%nz
               star%neutrino_flux_zone(j,i) = star%neutrino_flux_zone(j,i) &
                    /star%flux%neutrino_flux_total(j)
            end do
         end if
      end do
      do i = 1,star%nz
! TEMPERATURE IN UNITS OF 10**6 K.
         t6_million_k = exp(ln10*(star%logT(i)-6.0d0))
         if (t6_million_k.lt.5.0d0) exit
! ELECTRON DENSITY.
         log_electron_density = star%logRho(i)+log10((1.0d0+star%xa(i_h1,i))/2.0d0)
! MASS FRACTION.
         zone_mass_fraction = shell_mass(i)/1.9891d33
! RADIUS FRACTION.
         zone_radius_fraction = exp(ln10*star%logR(i))/star%solar_radius_cgs
! FLUXES ARE PRINTED IN THE SAME ORDER AS BAHCALL AND PINSONNEAULT.
         write(short_file_unit,145)zone_radius_fraction,t6_million_k, &
              log_electron_density,zone_mass_fraction, &
              star%be7_mass_fraction_zone(i), &
              star%neutrino_flux_zone(i_nu_pp,i), &
              star%neutrino_flux_zone(i_nu_b8,i), &
              star%neutrino_flux_zone(i_nu_n13,i), &
              star%neutrino_flux_zone(i_nu_o15,i), &
              star%neutrino_flux_zone(i_nu_f17,i), &
              star%neutrino_flux_zone(i_nu_be7,i), &
              star%neutrino_flux_zone(i_nu_pep,i), &
              star%neutrino_flux_zone(i_nu_hep,i)
 145     format(F9.5,F7.3,F6.3,1P10E10.3)
      end do
      return
end subroutine neutrino_flux_table
