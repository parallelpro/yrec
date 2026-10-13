!----------------------------------------------------------------------
! neutrino_flux_output
!----------------------------------------------------------------------
! The optional per-shell neutrino-flux table of each kind card's
! starting model (write_neutrino_flux_table, switched on by the
! compute_neutrino_fluxes control).
module neutrino_flux_output
      use star_info_lib, only: star_info
      implicit none
      private
      public :: write_neutrino_flux_table

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! write_neutrino_flux_table (formerly the LNUTAB block of program
! main, MHP 2/04)
!----------------------------------------------------------------------
! An optional diagnostic (compute_neutrino_fluxes, default .false.),
! called by run_yrec once per kind card right after
! read_starting_model, so the table describes each card's starting
! model. It writes to run_log_unit:
!  1. for each shell from the center out, as long as log T >= 6 and
!     log T is above the nuclear cutoff: the shell index, the shell
!     mass (g) and the shell's contributions to the ten flux slots (the
!     per-gram fluxes from compute_zone_energy times the shell mass);
!  2. the totals of the ten slots over those shells;
!  3. for each shell from the center out while T >= 5e6 K, the
!     normalized production profiles in the order Bahcall &
!     Pinsonneault use: r/Rsun, T6 = T/1e6 K, log(n_e/N_A), the shell
!     mass (Msun), X(Be7), then the pp, B8, N13, O15, F17, Be7, pep and
!     hep fractions.
! The fluxes are those at the Earth, in units of 10^10 cm^-2 s^-1 (see
! compute_solar_neutrino_fluxes). Slots 1-8 are pp, pep, hep, Be7, B8,
! N13, O15 and F17; slots 9-10 are the He3+He3 and He3+He4 reaction
! rates in the same units. The shell masses are local arrays, so the
! model's star%m and star%dm are not touched.
subroutine write_neutrino_flux_table(s)
      use star_info_lib, only: json, i_h1, i_h2, i_he3, i_he4, i_c12, &
           i_c13, i_n14, i_o16, i_o18, i_nu_pp, i_nu_pep, i_nu_hep, &
           i_nu_be7, i_nu_b8, i_nu_n13, i_nu_o15, i_nu_f17, n_nu_fluxes
      use io_units
      use physical_constants
      use net_lib, only: nuclear_energy_result
      use burn_step, only: compute_zone_energy
      use math_lib
      use observables, only: shell_masses_from_log_mass
      implicit none
      type(star_info), intent(inout) :: s

      double precision :: shell_center_mass(json), shell_mass(json)
      double precision :: shell_log_density, shell_log_temperature
! burn: the shell's energy record; only its solar-neutrino fluxes and
! Be7 abundance are used here.
      type(nuclear_energy_result) :: burn
      double precision :: t6_million_k, log_electron_density, &
           zone_mass_fraction, zone_radius_fraction
      integer :: i, j

! The shell masses from the log m grid: shell_center_mass is the mass
! (g) inside each shell center, and shell_mass the mass (g) of each
! shell.
      call shell_masses_from_log_mass(s%log_mass, s%log_total_mass, &
           s%nz, shell_center_mass, shell_mass)

      do j = 1,n_nu_fluxes
         s%neutrino_flux_total(j) = 0.0d0
         do i = 1,s%nz
            s%neutrino_flux_zone(j,i) = 0.0d0
         end do
      end do

! The production in each shell, from compute_zone_energy.
      do i = 1,s%nz
         shell_log_density = s%logRho(i)
         shell_log_temperature = s%logT(i)
! Stop at the first shell cooler than log T = 6, or at or below the
! nuclear cutoff, where nothing burns.
         if (shell_log_temperature.lt.6.0d0 .or. &
              shell_log_temperature.le.s%ctrl%nuclear_logT_cutoffs(1)) exit
         call compute_zone_energy(s, shell_log_density, &
              shell_log_temperature, s%xa(:,i), i, burn)
! The Be7 mass fraction (zero where X <= 1e-6, where none is computed).
         s%be7_mass_fraction_zone(i) = burn%be7_mass_fraction
! The shell's fluxes: the per-gram fluxes times the shell mass (g).
         do j = 1,n_nu_fluxes
            s%neutrino_flux_zone(j,i) = burn%neutrino_flux(j)*shell_mass(i)
            s%neutrino_flux_total(j) = s%neutrino_flux_total(j) &
                 + s%neutrino_flux_zone(j,i)
         end do
         write(run_log_unit,911)i,shell_mass(i), &
              (s%neutrino_flux_zone(j,i),j=1,n_nu_fluxes)
 911     format(I5,1P11E10.3)
      end do

! The total fluxes, in units of 10^10 cm^-2 s^-1.
      write(run_log_unit,222)(s%neutrino_flux_total(i),i=1,n_nu_fluxes)
 222  format(1P10E10.3)

! Normalize each shell's flux to the total of its slot.
      do j = 1,n_nu_fluxes
         if (s%neutrino_flux_total(j) .ne. 0.0d0) then
            do i = 1,s%nz
               s%neutrino_flux_zone(j,i) = s%neutrino_flux_zone(j,i) &
                    /s%neutrino_flux_total(j)
            end do
         end if
      end do
      do i = 1,s%nz
! T6 = T/1e6 K; stop below 5e6 K.
         t6_million_k = exp(ln10*(s%logT(i)-6.0d0))
         if (t6_million_k.lt.5.0d0) exit
! The electron density for full ionization,
! log(n_e/N_A) = log[rho (1+X)/2].
         log_electron_density = s%logRho(i)+log10((1.0d0+s%xa(i_h1,i))/2.0d0)
! The shell mass in solar masses.
         zone_mass_fraction = shell_mass(i)/s%solar_mass_cgs
! The radius in solar radii.
         zone_radius_fraction = exp(ln10*s%logR(i))/s%solar_radius_cgs
! The fluxes in the order Bahcall & Pinsonneault use.
         write(run_log_unit,145)zone_radius_fraction,t6_million_k, &
              log_electron_density,zone_mass_fraction, &
              s%be7_mass_fraction_zone(i), &
              s%neutrino_flux_zone(i_nu_pp,i), &
              s%neutrino_flux_zone(i_nu_b8,i), &
              s%neutrino_flux_zone(i_nu_n13,i), &
              s%neutrino_flux_zone(i_nu_o15,i), &
              s%neutrino_flux_zone(i_nu_f17,i), &
              s%neutrino_flux_zone(i_nu_be7,i), &
              s%neutrino_flux_zone(i_nu_pep,i), &
              s%neutrino_flux_zone(i_nu_hep,i)
 145     format(F9.5,F7.3,F6.3,1P10E10.3)
      end do
      return
end subroutine write_neutrino_flux_table

end module neutrino_flux_output
