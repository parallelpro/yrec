!----------------------------------------------------------------------
! nuclear_energy
!----------------------------------------------------------------------
! The nuclear energy generation of one zone. compute_nuclear_energy
! returns, in a nuclear_energy_result:
!   - eps (erg/g/s) and its absolute derivatives d eps/d ln rho and
!     d eps/d ln T, for the Henyey Jacobian;
!   - the energy of each chain (eps_by_chain, network's chain_* slots:
!     pp I with d(p,g), the two Be7 branches, CNO, the alpha captures,
!     triple alpha);
!   - the thermal neutrino loss (included in eps) and the He3 energy
!     diagnostics;
!   - when net_settings%calc_neutrinos, the solar neutrino fluxes and the
!     equilibrium Be7 abundance.
! It reads no model state; its configuration is net_settings. The energy
! of each network reaction comes from the reaction table (network's
! compute_reaction_energies); deuterium burning, outside the table, is
! added here. Sources and history: the end of this file.
module nuclear_energy
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use neutrino_flux_layout, only: n_nu_fluxes
      use net_configuration, only: net_settings, num_reactions, r_pp, &
           r_he3he3, r_he3he4, r_pc12, r_pn14, r_po16, xs_pep
      use reaction_rates, only: nuclear_rates, screening_mass_fractions, &
           compute_deuterium_proton_capture, cap_deuterium_rate
      use thermal_neutrinos, only: compute_thermal_neutrino_losses
      use network, only: compute_reaction_rates, reaction_counts, &
           compute_reaction_energies, n_chains, n_network, &
           network_species_slot, chain_pp1
      implicit none
      private

      type, public :: nuclear_energy_result
! .false. at or below nuclear_logT_cutoffs(1): nothing else is set
         logical :: burning = .false.
! erg/g/s; derivatives absolute (d eps/dln rho, d eps/dln T)
         real(dp) :: total_energy_gen_rate = 0.0d0
         real(dp) :: deps_dlnrho = 0.0d0, deps_dlnt = 0.0d0
! each chain's energy generation, erg/g/s (network's chain_* slots)
         real(dp) :: eps_by_chain(n_chains) = 0.0d0
         real(dp) :: eps_neutrino_loss = 0.0d0
! set only when net_settings%calc_neutrinos (fluxes_set) and X > 1e-6 (be7_set)
         logical :: fluxes_set = .false., be7_set = .false.
         real(dp) :: neutrino_flux(n_nu_fluxes) = 0.0d0
         real(dp) :: be7_mass_fraction = 0.0d0
      end type nuclear_energy_result

      public :: compute_nuclear_energy

! MeV to erg
      real(dp), parameter :: erg_per_mev = 1.602177d-6
! deuterium is burned (and its energy counted) above this mass fraction
      real(dp), parameter :: min_deuterium = 1.0d-11
! energy released by d(p,g)3He, MeV
      real(dp), parameter :: q_deuterium_mev = 5.494d0
! the energy path drops a rate of the S-factor reactions (r_pp to r_po16)
! below this (per gram per s, abundance factors omitted), with its
! derivatives (the composition path drops rates at or below
! network's min_composition_rate instead)
      real(dp), parameter :: min_energy_path_rate = 1.0d-30
! eps (erg/g/s, before the thermal neutrino loss) at or below which its
! derivatives are zero, and the solar neutrino fluxes' reaction counts
      real(dp), parameter :: min_eps_with_derivatives = 1.0d-12
! X(H) at or below which the solar neutrino fluxes are zero (MHP 9/91)
      real(dp), parameter :: min_flux_hydrogen = 1.0d-6

contains

!----------------------------------------------------------------------
! compute_nuclear_energy: the energy generation of one zone.
!   xa                  the zone's composition, star%xa layout (i_* slots)
!   in_surface_cz       the zone lies in the surface convection zone: there
!   turnover_timescale  deuterium burns no faster than the convective
!                       overturn (when turnover_timescale > 1 s)
pure subroutine compute_nuclear_energy(log_density, log_temperature, xa, &
     in_surface_cz, turnover_timescale, e)
      real(dp), intent(in) :: log_density, log_temperature, xa(:)
      logical, intent(in) :: in_surface_cz
      real(dp), intent(in) :: turnover_timescale
      type(nuclear_energy_result), intent(out) :: e

      type(nuclear_rates) :: nr
! the energy path's rates (energy_path_rates)
      real(dp) :: rate(num_reactions), dlnrate_dlnrho(num_reactions), &
           dlnrate_dlnt(num_reactions)
! each slot's rate times its reactants' mass fractions, per gram per s
      real(dp) :: reactions_per_gram(num_reactions)
      real(dp) :: energy_by_slot(num_reactions), energy_by_chain(n_chains)
      real(dp) :: deuterium_energy, deuterium_dlneps_dlnrho, deuterium_dlneps_dlnt
! the network species' mass fractions (network's s_* order)
      real(dp) :: x(n_network)

      if (log_temperature <= net_settings%logT_cutoffs(1)) return
      e%burning = .true.

      call compute_reaction_rates(log_density, log_temperature, &
           screening_mass_fractions(xa), nr)
      call energy_path_rates(nr, rate, dlnrate_dlnrho, dlnrate_dlnt)
      call compute_deuterium_burning(nr, xa, in_surface_cz, turnover_timescale, &
           deuterium_energy, deuterium_dlneps_dlnrho, deuterium_dlneps_dlnt)
      x = network_mass_fractions(xa)
      reactions_per_gram = reaction_counts(rate, x)
      call compute_reaction_energies(rate, x, nr%frac_n15_pa, nr%be7%frac_ecap, &
           energy_by_slot, energy_by_chain)
      call sum_energy_generation(energy_by_slot, dlnrate_dlnrho, dlnrate_dlnt, &
           deuterium_energy, deuterium_dlneps_dlnrho, deuterium_dlneps_dlnt, e)
      call store_chain_energies(energy_by_chain, deuterium_energy, e)
      if (net_settings%calc_neutrinos) &
           call compute_solar_neutrino_fluxes(nr, xa, reactions_per_gram, e)
      if (log_temperature > net_settings%logT_cutoffs(5)) &
           call add_thermal_neutrino_loss(log_density, log_temperature, nr%zone%t%t9, xa, e)
end subroutine compute_nuclear_energy

!----------------------------------------------------------------------
! The energy path's copy of the rates: those of the S-factor reactions
! below min_energy_path_rate are set to zero with their derivatives.
pure subroutine energy_path_rates(nr, rate, dlnrate_dlnrho, dlnrate_dlnt)
      type(nuclear_rates), intent(in) :: nr
      real(dp), intent(out) :: rate(num_reactions), &
           dlnrate_dlnrho(num_reactions), dlnrate_dlnt(num_reactions)
      integer :: i

      rate = nr%rate
      dlnrate_dlnrho = nr%dlnrate_dlnrho
      dlnrate_dlnt = nr%dlnrate_dlnt
      do i = r_pp, r_po16
         if (rate(i) < min_energy_path_rate) then
            rate(i) = 0.d0
            dlnrate_dlnrho(i) = 0.d0
            dlnrate_dlnt(i) = 0.d0
         end if
      end do
end subroutine energy_path_rates

!----------------------------------------------------------------------
! Deuterium burning, d(p,g)3He (MHP 5/02): its energy generation (erg/g/s)
! and d ln eps / d ln rho, d ln eps / d ln T. Unscreened -- reasonable at
! the low temperatures involved -- and, in the surface convection zone,
! no faster than the convective overturn: there the capped rate,
! N_A/(m_d tau) per unit X, depends on neither rho nor T, and both
! derivatives are zero. Zero at or below min_deuterium.
pure subroutine compute_deuterium_burning(nr, xa, in_surface_cz, turnover_timescale, &
     energy, dlneps_dlnrho, dlneps_dlnt)
      use composition_layout, only: i_h1, i_h2
      type(nuclear_rates), intent(in) :: nr
      real(dp), intent(in) :: xa(:), turnover_timescale
      logical, intent(in) :: in_surface_cz
      real(dp), intent(out) :: energy, dlneps_dlnrho, dlneps_dlnt
      real(dp) :: rate
      logical :: capped

      energy = 0.0d0
      dlneps_dlnrho = 0.0d0
      dlneps_dlnt = 0.0d0
      if (xa(i_h2) <= min_deuterium) return
      call compute_deuterium_proton_capture(nr%zone%density, nr%zone%t, rate, &
           dlneps_dlnt)
      dlneps_dlnrho = 1.0d0
      call cap_deuterium_rate(xa(i_h1), in_surface_cz, turnover_timescale, &
           rate, capped)
      if (capped) then
         dlneps_dlnrho = 0.0d0
         dlneps_dlnt = 0.0d0
      end if
      energy = rate*xa(i_h1)*xa(i_h2)*q_deuterium_mev*erg_per_mev
end subroutine compute_deuterium_burning


!----------------------------------------------------------------------
! The network species' mass fractions (network's s_* order) from a
! zone's composition in the star%xa layout.
pure function network_mass_fractions(xa) result(x)
      real(dp), intent(in) :: xa(:)
      real(dp) :: x(n_network)

      x = xa(network_species_slot)
end function network_mass_fractions

!----------------------------------------------------------------------
! eps and its absolute derivatives: the sum over reactions of eps_j and
! eps_j * d ln rate_j / d ln (rho, T), plus deuterium burning's. At or
! below eps = min_eps_with_derivatives the derivatives are zero (and so
! are the solar neutrino fluxes' reaction counts:
! compute_solar_neutrino_fluxes). These are absolute derivatives, not
! d ln eps: the Henyey Jacobian and the pulsation output (GYRE's
! eps*eps_T) both want this form.
pure subroutine sum_energy_generation(energy_by_slot, dlnrate_dlnrho, dlnrate_dlnt, &
     deuterium_energy, deuterium_dlneps_dlnrho, deuterium_dlneps_dlnt, e)
      real(dp), intent(in) :: energy_by_slot(num_reactions), &
           deuterium_energy, deuterium_dlneps_dlnrho, deuterium_dlneps_dlnt
      real(dp), intent(in) :: dlnrate_dlnrho(num_reactions), &
           dlnrate_dlnt(num_reactions)
      type(nuclear_energy_result), intent(inout) :: e
      real(dp) :: sum_rho, sum_t
      integer :: i

      e%total_energy_gen_rate = 0.0d0
      sum_rho = 0.0d0
      sum_t = 0.0d0
      do i = 1, num_reactions
         e%total_energy_gen_rate = e%total_energy_gen_rate+energy_by_slot(i)
         sum_rho = sum_rho+energy_by_slot(i)*dlnrate_dlnrho(i)
         sum_t = sum_t+energy_by_slot(i)*dlnrate_dlnt(i)
      end do
      e%total_energy_gen_rate = e%total_energy_gen_rate + deuterium_energy
      sum_rho = sum_rho + deuterium_energy*deuterium_dlneps_dlnrho
      sum_t = sum_t + deuterium_energy*deuterium_dlneps_dlnt
      if (e%total_energy_gen_rate <= min_eps_with_derivatives) then
         e%deps_dlnrho = 0.d0
         e%deps_dlnt = 0.d0
      else
         e%deps_dlnrho = sum_rho
         e%deps_dlnt = sum_t
      end if
end subroutine sum_energy_generation

!----------------------------------------------------------------------
! The energy of each chain (Neutrino Astrophysics, table 21): each table
! row counts under its chain; pp I also gets the d(p,g) energy (deuterium
! is burned outside the table). pp II and pp III are the Be7 electron-
! and proton-capture branches, each with its share of the He3 + He4
! release. (The original Yale routine carried two branching factors on
! those two instead of one.) Triple alpha is negligible in the Sun.
pure subroutine store_chain_energies(energy_by_chain, deuterium_energy, e)
      real(dp), intent(in) :: energy_by_chain(n_chains), deuterium_energy
      type(nuclear_energy_result), intent(inout) :: e

      e%eps_by_chain = energy_by_chain
      e%eps_by_chain(chain_pp1) = e%eps_by_chain(chain_pp1) + deuterium_energy
end subroutine store_chain_energies

!----------------------------------------------------------------------
! The solar neutrino fluxes at the Earth from one gram of this zone, in
! 10^10 per cm^2 per s per g (Neutrino Astrophysics, tables 3.1 and 3.3,
! eqs. 6.1-6.8), assuming nothing happens to the neutrinos after they are
! made; slots i_nu_he3he3/he3he4 hold the He3+He3 and He3+He4 reaction
! counts in the same units (diagnostic "fictional" fluxes). Also the
! equilibrium Be7 mass fraction (Be7 is always in equilibrium with the
! slower He3 + He4 production; computed because some applications want
! it). Off in hydrogen-exhausted zones (X <= min_flux_hydrogen, MHP
! 9/91). reactions_per_gram is each slot's reactions per gram per s; it
! counts as zero where eps is at or below min_eps_with_derivatives (e,
! from sum_energy_generation, has the zone's eps), so that only the hep
! flux, from its own rate factor, is left there.
pure subroutine compute_solar_neutrino_fluxes(nr, xa, reactions_per_gram, e)
      use composition_layout, only: i_h1, i_he3
      use neutrino_flux_layout, only: i_nu_b8, i_nu_be7, i_nu_f17, i_nu_hep, &
           i_nu_n13, i_nu_o15, i_nu_pep, i_nu_pp, i_nu_he3he3, i_nu_he3he4
      type(nuclear_rates), intent(in) :: nr
      real(dp), intent(in) :: xa(:), reactions_per_gram(num_reactions)
      type(nuclear_energy_result), intent(inout) :: e
! 4 pi (1 AU)^2, cm^2
      real(dp), parameter :: four_pi_au2 = 2.812295d+27
! the fluxes' unit, 1e10 per cm^2 per s per g, inverted
      real(dp), parameter :: flux_unit_inverse = 1.0d-10
! the reactions per gram per s that make the neutrinos
      real(dp) :: counts(num_reactions)
      integer :: k

      e%fluxes_set = .true.
      if (xa(i_h1) <= min_flux_hydrogen) then
         e%neutrino_flux = 0.0d0
         return
      end if
      if (e%total_energy_gen_rate <= min_eps_with_derivatives) then
         counts = 0.0d0
      else
         counts = reactions_per_gram
      end if
      e%neutrino_flux(i_nu_pp) = counts(r_pp)/four_pi_au2
! pep from pp (Neutrino Astrophysics eq. 3.17); xs_pep scales the pep/pp
! ratio, so leave it at 1 unless that ratio is meant to change
      e%neutrino_flux(i_nu_pep) = (3.4848d-6)*nr%zone%electron_density_na*nr%zone%t%m12* &
           (1.0d0 + 20.d0*nr%zone%t%t9)*counts(r_pp)
      e%neutrino_flux(i_nu_pep) = e%neutrino_flux(i_nu_pep)*net_settings%extra_cross_section_scale(xs_pep)/four_pi_au2
! hep (eq. 3.12), from the screened rate factor of compute_reaction_rates
      e%neutrino_flux(i_nu_hep) = nr%hep
      e%neutrino_flux(i_nu_hep) = e%neutrino_flux(i_nu_hep)*xa(i_h1)*xa(i_he3)/ &
           four_pi_au2
      e%be7_set = .true.
      e%be7_mass_fraction = counts(r_he3he4)/(nr%be7%proton + nr%be7%electron)
      e%neutrino_flux(i_nu_be7) = counts(r_he3he4)*nr%be7%frac_ecap/four_pi_au2
      e%neutrino_flux(i_nu_b8) = counts(r_he3he4)*nr%be7%frac_pcap/four_pi_au2
      e%neutrino_flux(i_nu_n13) = counts(r_pc12)/four_pi_au2
      e%neutrino_flux(i_nu_o15) = counts(r_pn14)/four_pi_au2
      e%neutrino_flux(i_nu_f17) = counts(r_po16)/four_pi_au2
      e%neutrino_flux(i_nu_he3he3) = counts(r_he3he3)/four_pi_au2
      e%neutrino_flux(i_nu_he3he4) = counts(r_he3he4)/four_pi_au2
! to 10^10 per cm^2 per s per g; non-positive values to zero (KC
! 2025-05-30: was a 1e-50 floor, which underflowed)
      do k = 1, n_nu_fluxes
         e%neutrino_flux(k) = flux_unit_inverse*e%neutrino_flux(k)
         if (e%neutrino_flux(k) <= 0.0d0) e%neutrino_flux(k) = 0.0d0
      end do
end subroutine compute_solar_neutrino_fluxes

!----------------------------------------------------------------------
! The thermal neutrino losses (MHP 5/91; for neutrino-cooled cores of
! evolved stars), subtracted from eps with their absolute derivatives:
! Itoh et al. (1996) when net_settings%use_itoh_neutrino_loss (Grant
! Newsham 9/06), otherwise the Beaudet, Petrosian & Salpeter (1967) fits
! for the pair, photo and plasma losses.
pure subroutine add_thermal_neutrino_loss(log_density, log_temperature, t9, xa, e)
      use composition_layout, only: i_h1, i_he4, i_c12, i_c13, i_o16, i_o18
      use physical_constants, only: cc13, ln10
      use math_lib
      real(dp), intent(in) :: log_density, log_temperature, t9, xa(:)
      type(nuclear_energy_result), intent(inout) :: e
! lambda = kT/m_e c^2 = T9/bps_electron_rest_energy_t9 (m_e c^2/k in
! 1e9 K, the fits' value; the screening's is 5.92986); the pair losses
! only from T9 = min_pair_t9
      real(dp), parameter :: bps_electron_rest_energy_t9 = 5.9302d0, min_pair_t9 = 0.2d0
! Beaudet, Petrosian & Salpeter fit coefficients: pair, photo, plasma
      real(dp), parameter :: pair_fit(7) = [6.002d19, 2.084d20, &
           1.872d21, 9.383d-1, -4.141d-1, 5.829d-2, 5.5924d0]
      real(dp), parameter :: photo_fit(7) = [4.886d10, 7.580d10, &
           6.023d10, 6.290d-3, 7.483d-3, 3.061d-4, 1.5654d0]
      real(dp), parameter :: plasma_fit(7) = [2.320d-7, 8.449d-8, &
           1.787d-8, 2.581d-2, 1.734d-2, 6.990d-4, 0.56457d0]
      real(dp) :: temperature, density, snu, dsnu_dt, dsnu_drho
      real(dp) :: el, eli, ez, ez3, emue, ex1, ex2, ex3, el2
      real(dp) :: polx10, polx11, polx12, polx21, polx22, polx31, &
           polx32, qedn, qetn, qetnx, qednx

      temperature = exp10(log_temperature)
      density = exp10(log_density)
      if (net_settings%use_itoh_neutrino_loss) then
! snu >= 0 and its ABSOLUTE derivatives d snu/dT, d snu/d rho, so the
! loss contributes d(-snu)/d ln T = -T d snu/dT (2026 bugsweep Batch 2:
! the original converted to log-derivatives, the wrong units here, and
! built the T line from the updated rho line)
         call compute_thermal_neutrino_losses(temperature, density, xa(i_h1), &
              xa(i_he4), xa(i_c12)+xa(i_c13), xa(i_o16)+xa(i_o18), &
              snu, dsnu_dt, dsnu_drho)
         e%eps_neutrino_loss = -snu
         e%total_energy_gen_rate = e%total_energy_gen_rate + e%eps_neutrino_loss
         e%deps_dlnrho = e%deps_dlnrho - density*dsnu_drho
         e%deps_dlnt = e%deps_dlnt - temperature*dsnu_dt
         return
      end if
      el = t9/bps_electron_rest_energy_t9
      eli = 1.0d0/el
      ez = exp(cc13*ln10*(log_density-9.0d0))*eli*(0.7937d0+0.2063d0*xa(i_h1))
      ez3 = ez**3
      emue = 0.5d0*(1.0d0+xa(i_h1))
      ex1 = 0.0d0
      if (t9 >= min_pair_t9) then
! pair neutrinos
         el2 = el*el
         polx10=(1.d0+el2*(-13.04d0+el2*(133.5d0+el2*(1534.d0+el2*918.6d0))))
         polx11 = pair_fit(1) + ez*(pair_fit(2) + ez*pair_fit(3))
         polx12 = ez3 + eli*(pair_fit(4) + eli*(pair_fit(5) + eli*pair_fit(6)))
         ex1 = exp(-ez*pair_fit(7)-eli-eli-ln10*log_density)*polx10*polx11/polx12
      end if
! photo neutrinos
      polx21 = photo_fit(1) + ez*(photo_fit(2) + ez*photo_fit(3))
      polx22 = ez3 + eli*(photo_fit(4) + eli*(photo_fit(5) + eli*photo_fit(6)))
      ex2 = emue*el**5*exp(-ez*photo_fit(7))*polx21/polx22
! plasma neutrinos
      polx31 = plasma_fit(1) + ez*(plasma_fit(2) + ez*plasma_fit(3))
      polx32 = ez3 + eli*(plasma_fit(4) + eli*(plasma_fit(5) + eli*plasma_fit(6)))
      ex3 = emue**3*exp(-ez*plasma_fit(7)+ln10*(log_density+log_density))*polx31/polx32
      e%eps_neutrino_loss = -(ex1 + ex2 + ex3)
      e%total_energy_gen_rate = e%total_energy_gen_rate + e%eps_neutrino_loss
! their d/d ln T (qetnx) and d/d ln rho (qednx)
      qetnx = 0.0d0
      qednx = 0.0d0
      if (t9 >= min_pair_t9) then
         qedn= ez*(pair_fit(2)+2.d0*ez*pair_fit(3))/polx11 - ez*pair_fit(7) &
               - 3.0d0*ez3/polx12
         qetn= -qedn+ eli*(pair_fit(4)+eli*(2.d0*pair_fit(5)+3.d0*eli*pair_fit(6)))/polx12
         qetnx= el2*(-26.08d0+el2*(534.d0+el2*(9204.d0+el2*7348.8d0)))/polx10
         qetnx = (qetnx + qetn + eli+eli)*ex1
         qednx = (-1.0d0 +cc13*qedn)*ex1
      end if
      qedn= ez*(photo_fit(2)+2.d0*ez*photo_fit(3))/polx21 - ez*photo_fit(7) - 3.d0*ez3/polx22
      qetn= -qedn+ eli*(photo_fit(4)+eli*(2.d0*photo_fit(5)+3.d0*eli*photo_fit(6)))/polx22
      qetnx = qetnx + (5.0d0 + qetn)*ex2
      qednx = qednx +cc13*qedn*ex2
      qedn= ez*(plasma_fit(2)+2.d0*ez*plasma_fit(3))/polx31 - ez*plasma_fit(7) - 3.d0*ez3/polx32
      qetn=-qedn + eli*(plasma_fit(4)+eli*(2.d0*plasma_fit(5)+3.d0*eli*plasma_fit(6)))/polx32
      qetnx = qetnx + qetn*ex3
      qednx = qednx + (2.0d0 +cc13*qedn)*ex3
      e%deps_dlnt = e%deps_dlnt - qetnx
      e%deps_dlnrho = e%deps_dlnrho - qednx
end subroutine add_thermal_neutrino_loss

end module nuclear_energy

!======================================================================
! Sources and history (engeb.f: J. N. Bahcall 11/6/1990, M. H.
! Pinsonneault, D. B. Guenther). The F77 comments are kept verbatim in
! net/retired/f77_rate_comments.txt.
!
! Legacy names: EPP1, EPP2, EPP3, ECN, E3AL and EALPCA -> eps_by_chain
! (chain_pp1, chain_pp2, chain_pp3, chain_cno, chain_triple_alpha,
! chain_alpha_capture); PEP and PET -> deps_dlnrho and deps_dlnt; SUM1 ->
! total_energy_gen_rate; EG(I) -> reactions_per_gram (each reaction's
! count per gram per s); FLUX -> neutrino_flux; X, Y, XHE3, ... -> the xa
! slots.
!
! Energy releases. From Bahcall & Ulrich (1988, RMP 60, 297), allowing
! for the neutrino losses. JNB recomputed all of them on 9/23-25/97 for
! the Solar Fusion workshop paper, from the neutral atomic mass
! differences of the Table of Isotopes (8th ed., 1996) and the neutrino
! losses of Bahcall (1997, Phys. Rev. C, the gallium experiments); his
! notes, vol. 19, pp. 132-141. They are the q_mev of network's table
! rows; the d(p,g) release is q_deuterium_mev here.
!
! Solar neutrino fluxes. The fluxes at the Earth per gram of the zone
! (Bahcall's section, 1990; NA tables 3.1 and 3.3, eqs. 6.1-6.8), assuming
! nothing happens to the neutrinos on the way; times a shell's mass they
! are its flux. The He3 + He3 and He3 + He4 entries are fictional fluxes,
! the reaction counts in the same units, useful diagnostics of a solar
! model. Zero where X <= 1e-6 (MHP 9/91). pep from pp: NA eq. 3.17.
!
! Thermal neutrino losses. Beaudet, Petrosian & Salpeter (1967) for the
! pair, photo and plasma losses (MHP 5/91; DBG 3/92 added new routines);
! Itoh et al. (1996) through thermal_neutrinos (G. Newsham 9/06).
!
! Deuterium burning. d(p,g)3He (MHP 5/02), from CF88's fit (the
! T9^(-2/3) form; 2026 bug sweep: it was T9^(+2/3)), capped in the surface
! convection zone by the convective overturn.
!
! Earlier Yale code took the rates mostly from Fowler, Caughlan &
! Zimmerman (1975, ARA&A), the Be7 and B8 rates from Bahcall & Sears
! (1972, ARA&A) and some from Harris et al. (1983); Bahcall corrected its
! numerical errors and physically wrong statements and made the cross
! sections explicit factors (SSTANDARD = 1 for NA's value).
!======================================================================
