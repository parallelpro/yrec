!----------------------------------------------------------------------
! neutrino_flux_table
!----------------------------------------------------------------------
! Resurrected 2026 as an opt-in diagnostic (controls_buffer's
! compute_neutrino_fluxes, default .false.). This is the former
! LNUTAB neutrino-table block of program main (MHP 2/04), which had
! been hardcoded off since 2004 and had bit-rotted while dead: the
! compute_zone_energy call passed uninitialized locals in the deuterium-fraction
! and shell-index slots, and the block overwrote star%m/star%dm with
! its own shell-center weights. Fixed on resurrection:
!   - compute_zone_energy receives star%xa(i_h2,k) and the loop index;
!   - the integration weights are LOCAL arrays (shell_center_mass,
!     shell_mass), star%m/star%dm are not touched;
!   - output goes to run_log_unit (the .short / CASE.log stream)
!     instead of the historical bare unit 76 (a stray fort.76 file)
!     and stdout.
! Called from run_yrec once per kind card, right after
! read_starting_model, so the
! table describes the starting model of each run.
!
! Content: per-shell neutrino production from compute_zone_energy (erg/g/s scaled to
! erg/s by the shell mass), the totals for the ten flux slots (units
! of 1e10; slots 1-8 = pp, pep, hep, Be7, B8, N13, O15, F17), then the
! normalized production profiles in the Bahcall & Pinsonneault layout
! (r/Rsun, T6, log n_e, shell mass fraction, X_Be7, then pp, B8, N13,
! O15, F17, Be7, pep, hep).
subroutine neutrino_flux_table
      use star_info_lib, only: star, json, i_h1, i_h2, i_he3, i_he4, &
           i_c12, i_c13, i_n14, i_o16, i_o18, &
           i_nu_pp, i_nu_pep, i_nu_hep, i_nu_be7, i_nu_b8, i_nu_n13, &
           i_nu_o15, i_nu_f17, n_nu_fluxes
      use io_units
      use physical_constants
      use net_lib, only: nuclear_energy_result
      use burn_step, only: compute_zone_energy
      use math_lib
      use observables, only: shell_masses_from_log_mass
      implicit none

      double precision :: shell_center_mass(json), shell_mass(json)
      double precision :: shell_log_density, shell_log_temperature
! burn: the shell's energy record; only its solar neutrino fluxes and
! Be7 abundance are used here.
      type(nuclear_energy_result) :: burn
      double precision :: t6_million_k, log_electron_density, &
           zone_mass_fraction, zone_radius_fraction
      integer :: i, j

! SET UP WEIGHTS AND MASSES (local; the star%m/star%dm model arrays
! are left untouched).
! shell_center_mass = LOCATION IN GM (UNLOGGED) OF SHELL CENTERS.
! shell_mass = MASS IN GM OF EACH SHELL.
      call shell_masses_from_log_mass(star%log_mass, star%log_total_mass, &
           star%nz, shell_center_mass, shell_mass)

      do j = 1,n_nu_fluxes
         star%neutrino_flux_total(j) = 0.0d0
         do i = 1,star%nz
            star%neutrino_flux_zone(j,i) = 0.0d0
         end do
      end do

! PER-SHELL PRODUCTION FROM ENGEB.
      do i = 1,star%nz
         shell_log_density = star%logRho(i)
         shell_log_temperature = star%logT(i)
! SKIP CALCULATIONS FOR LOW TEMPERATURES -- below log T = 6, and at or
! below the nuclear cutoff, where nothing burns.
         if (shell_log_temperature.lt.6.0d0 .or. &
              shell_log_temperature.le.star%ctrl%nuclear_logT_cutoffs(1)) exit
         call compute_zone_energy(shell_log_density, &
              shell_log_temperature, star%xa(:,i), i, burn)
! BE7 MASS FRACTION (2026-10: zero where X <= 1e-6, where none is
! computed; this used to repeat the last shell's value that had one)
         star%be7_mass_fraction_zone(i) = burn%be7_mass_fraction
! CONVERT FROM ERG/GM/S TO ERG/S FOR EACH SHELL BY MULTIPLYING
! BY THE MASS OF EACH SHELL IN GM.
         do j = 1,n_nu_fluxes
            star%neutrino_flux_zone(j,i) = burn%neutrino_flux(j)*shell_mass(i)
            star%neutrino_flux_total(j) = star%neutrino_flux_total(j) &
                 + star%neutrino_flux_zone(j,i)
         end do
         write(run_log_unit,911)i,shell_mass(i), &
              (star%neutrino_flux_zone(j,i),j=1,n_nu_fluxes)
 911     format(I5,1P11E10.3)
      end do

! WRITE OUT TOTAL NEUTRINO FLUXES.
! ***NOTE THAT THESE ARE IN UNITS OF 10**10. ***
      write(run_log_unit,222)(star%neutrino_flux_total(i),i=1,n_nu_fluxes)
 222  format(1P10E10.3)

! NORMALIZE FLUXES.
      do j = 1,n_nu_fluxes
         if (star%neutrino_flux_total(j) .ne. 0.0d0) then
            do i = 1,star%nz
               star%neutrino_flux_zone(j,i) = star%neutrino_flux_zone(j,i) &
                    /star%neutrino_flux_total(j)
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
         zone_mass_fraction = shell_mass(i)/star%solar_mass_cgs
! RADIUS FRACTION.
         zone_radius_fraction = exp(ln10*star%logR(i))/star%solar_radius_cgs
! FLUXES ARE PRINTED IN THE SAME ORDER AS BAHCALL AND PINSONNEAULT.
         write(run_log_unit,145)zone_radius_fraction,t6_million_k, &
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
