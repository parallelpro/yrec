!----------------------------------------------------------------------
! nuclear_energy
!----------------------------------------------------------------------
! The nuclear energy generation of one zone. compute_nuclear_energy
! returns, in a nuclear_energy_result:
!   - eps (erg/g/s) and its absolute derivatives d eps/d ln rho and
!     d eps/d ln T, for the Henyey Jacobian;
!   - the energy of each chain (pp I with d(p,g), the two Be7 branches,
!     CNO, triple alpha) and of the alpha captures;
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
      use reaction_rates, only: nuclear_rates, &
           net_settings, screening_mass_fractions, &
           deuterium_proton_capture_rate, avogadro_dpg, deuteron_mass_amu, &
           num_reactions, r_pp, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, &
           r_po16, r_ac13, r_ao16_unused, r_ac12, r_an14, r_3alpha, &
           r_c12c12_unused, xs_pep
      use thermal_neutrinos, only: compute_thermal_neutrino_losses
      use network, only: compute_reaction_rates, reaction_counts, &
           compute_reaction_energies, n_chains, n_network, &
           network_species_slot, chain_pp1, chain_pp2, chain_pp3, chain_cno, chain_alpha_capture, &
           chain_triple_alpha
      implicit none
      private

      type, public :: nuclear_energy_result
! .false. at or below nuclear_logT_cutoffs(1): nothing else is set
         logical :: burning = .false.
! erg/g/s; derivatives absolute (d eps/dln rho, d eps/dln T)
         real(dp) :: total_energy_gen_rate = 0.0d0
         real(dp) :: deps_dlnrho = 0.0d0, deps_dlnt = 0.0d0
         real(dp) :: pp_chain_energy_gen = 0.0d0
         real(dp) :: he3he4_be7_electron_energy_gen = 0.0d0
         real(dp) :: he3he4_be7_proton_energy_gen = 0.0d0
         real(dp) :: cno_cycle_energy_gen = 0.0d0
         real(dp) :: triple_alpha_energy_gen = 0.0d0
         real(dp) :: eps_alpha_capture = 0.0d0
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

contains

!----------------------------------------------------------------------
! compute_nuclear_energy: the energy generation of one zone.
!   xa                  the zone's composition, star%xa layout (i_* slots)
!   in_surface_cz       the zone lies in the surface convection zone: there
!   turnover_timescale  deuterium burns no faster than the convective
!                       overturn (when turnover_timescale > 1 s)
subroutine compute_nuclear_energy(log_density, log_temperature, xa, &
     in_surface_cz, turnover_timescale, e)
      real(dp), intent(in) :: log_density, log_temperature, xa(:)
      logical, intent(in) :: in_surface_cz
      real(dp), intent(in) :: turnover_timescale
      type(nuclear_energy_result), intent(out) :: e

      type(nuclear_rates) :: nr
! the energy path's rates: reactions 1-7 below 1e-30 are dropped (with
! their derivatives) before eps is summed (the composition path drops
! rates at or below 1e-5 instead)
      real(dp) :: rate(num_reactions), dlnrate_dlnrho(num_reactions), &
           dlnrate_dlnt(num_reactions)
! each slot's rate times its reactants' mass fractions, per gram per s
      real(dp) :: reactions_per_gram(num_reactions)
      real(dp) :: energy_by_slot(num_reactions), energy_by_chain(n_chains)
      real(dp) :: deuterium_energy, deuterium_dlneps_dlnt

      if (log_temperature <= net_settings%logT_cutoffs(1)) return
      e%burning = .true.

      call compute_reaction_rates(log_density, log_temperature, &
           screening_mass_fractions(xa), nr)
      call energy_path_rates(nr, rate, dlnrate_dlnrho, dlnrate_dlnt)
      call compute_deuterium_burning(nr, xa, in_surface_cz, turnover_timescale, &
           deuterium_energy, deuterium_dlneps_dlnt)
      reactions_per_gram = reaction_counts(rate, network_mass_fractions(xa))
      call compute_reaction_energies(rate, network_mass_fractions(xa), &
           nr%frac_n15_pa, nr%frac_be7_ecap, energy_by_slot, energy_by_chain)
      call sum_energy_generation(energy_by_slot, dlnrate_dlnrho, dlnrate_dlnt, &
           deuterium_energy, deuterium_dlneps_dlnt, reactions_per_gram, e)
      call store_chain_energies(energy_by_chain, deuterium_energy, e)
      if (net_settings%calc_neutrinos) &
           call compute_solar_neutrino_fluxes(nr, xa, reactions_per_gram, e)
      if (log_temperature > net_settings%logT_cutoffs(5)) &
           call add_thermal_neutrino_loss(log_density, log_temperature, nr%t9, xa, e)
end subroutine compute_nuclear_energy

!----------------------------------------------------------------------
! The energy path's copy of the rates: reactions 1-7 below 1e-30 are
! set to zero together with their derivatives.
pure subroutine energy_path_rates(nr, rate, dlnrate_dlnrho, dlnrate_dlnt)
      type(nuclear_rates), intent(in) :: nr
      real(dp), intent(out) :: rate(num_reactions), &
           dlnrate_dlnrho(num_reactions), dlnrate_dlnt(num_reactions)
      integer :: i

      rate = nr%rate
      dlnrate_dlnrho = nr%dlnrate_dlnrho
      dlnrate_dlnt = nr%dlnrate_dlnt
      do i = r_pp, r_po16
         if (rate(i) < 1.d-30) then
            rate(i) = 0.d0
            dlnrate_dlnrho(i) = 0.d0
            dlnrate_dlnt(i) = 0.d0
         end if
      end do
end subroutine energy_path_rates

!----------------------------------------------------------------------
! Deuterium burning, d(p,g)3He (MHP 5/02): its energy generation (erg/g/s)
! and d ln eps / d ln T. Unscreened -- reasonable at the low temperatures
! involved -- and, in the surface convection zone, no faster than the
! convective overturn. Zero at or below min_deuterium.
pure subroutine compute_deuterium_burning(nr, xa, in_surface_cz, turnover_timescale, &
     energy, dlneps_dlnt)
      use composition_layout, only: i_h1, i_h2
      use physical_constants, only: cc13
      type(nuclear_rates), intent(in) :: nr
      real(dp), intent(in) :: xa(:), turnover_timescale
      logical, intent(in) :: in_surface_cz
      real(dp), intent(out) :: energy, dlneps_dlnt
      real(dp) :: rate, max_rate, zz, tfac, tfac_deriv

      energy = 0.0d0
      dlneps_dlnt = 0.0d0
      if (xa(i_h2) <= min_deuterium) return
! CF88's fit; zz and tfac are its exponent and polynomial, for the
! T-derivative (2026 bugsweep Batch 2: T9**(-2/3), as in CF88)
      zz = -3.72d0*nr%t9_m13
      tfac = 1.0d0+0.112d0*nr%t9_p13+3.38d0*nr%t9_p23+2.65d0*nr%t9
      rate = deuterium_proton_capture_rate(nr%density, nr%t9, nr%t9_p13, &
           nr%t9_p23, nr%t9_m13, nr%t9_m23)
      tfac_deriv = 0.112d0*nr%t9_p13+6.76d0*nr%t9_p23+7.95d0*nr%t9
      dlneps_dlnt = cc13*((tfac_deriv/tfac) -2.0d0 - zz)
      if (in_surface_cz .and. turnover_timescale > 1.0d0) then
         max_rate = avogadro_dpg/deuteron_mass_amu/turnover_timescale
         if (rate*xa(i_h1) > max_rate) then
            if (xa(i_h1) > 1.0d-6) rate = max_rate/xa(i_h1)
         end if
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
! eps_j * d ln rate_j / d ln (rho, T), plus deuterium burning (whose rate
! depends on rho linearly). Below eps = 1e-12 the derivatives are zero,
! and so are the reaction counts the neutrino fluxes use. These are
! absolute derivatives, not d ln eps: the Henyey Jacobian and the
! pulsation output (GYRE's eps*eps_T) both want this form.
pure subroutine sum_energy_generation(energy_by_slot, dlnrate_dlnrho, dlnrate_dlnt, &
     deuterium_energy, deuterium_dlneps_dlnt, reactions_per_gram, e)
      real(dp), intent(in) :: energy_by_slot(num_reactions), &
           deuterium_energy, deuterium_dlneps_dlnt
      real(dp), intent(inout) :: dlnrate_dlnrho(num_reactions), &
           dlnrate_dlnt(num_reactions), reactions_per_gram(num_reactions)
      type(nuclear_energy_result), intent(inout) :: e
      real(dp) :: sum_rho, sum_t
      integer :: i

! reactions 9 (O16 + alpha) and 13 (C12 + C12) are never computed
      dlnrate_dlnrho(r_ao16_unused) = 0.0d0
      dlnrate_dlnt(r_ao16_unused) = 0.0d0
      dlnrate_dlnrho(r_c12c12_unused) = 0.0d0
      dlnrate_dlnt(r_c12c12_unused) = 0.0d0
      e%total_energy_gen_rate = 0.0d0
      sum_rho = 0.0d0
      sum_t = 0.0d0
      do i = 1, num_reactions
         e%total_energy_gen_rate = e%total_energy_gen_rate+energy_by_slot(i)
         sum_rho = sum_rho+energy_by_slot(i)*dlnrate_dlnrho(i)
         sum_t = sum_t+energy_by_slot(i)*dlnrate_dlnt(i)
      end do
      e%total_energy_gen_rate = e%total_energy_gen_rate + deuterium_energy
      sum_rho = sum_rho + deuterium_energy
      sum_t = sum_t + deuterium_energy*deuterium_dlneps_dlnt
      if (e%total_energy_gen_rate <= 1.d-12) then
         e%deps_dlnrho = 0.d0
         e%deps_dlnt = 0.d0
         reactions_per_gram = 0.d0
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

      e%pp_chain_energy_gen = energy_by_chain(chain_pp1) + deuterium_energy
      e%he3he4_be7_electron_energy_gen = energy_by_chain(chain_pp2)
      e%he3he4_be7_proton_energy_gen = energy_by_chain(chain_pp3)
      e%cno_cycle_energy_gen = energy_by_chain(chain_cno)
      e%triple_alpha_energy_gen = energy_by_chain(chain_triple_alpha)
      e%eps_alpha_capture = energy_by_chain(chain_alpha_capture)
end subroutine store_chain_energies

!----------------------------------------------------------------------
! The solar neutrino fluxes at the Earth from one gram of this zone, in
! 10^10 per cm^2 per s per g (Neutrino Astrophysics, tables 3.1 and 3.3,
! eqs. 6.1-6.8), assuming nothing happens to the neutrinos after they are
! made; slots i_nu_he3he3/he3he4 hold the He3+He3 and He3+He4 reaction
! counts in the same units (diagnostic "fictional" fluxes). Also the
! equilibrium Be7 mass fraction (Be7 is always in equilibrium with the
! slower He3 + He4 production; computed because some applications want
! it). Off in hydrogen-exhausted zones (X <= 1e-6, MHP 9/91).
pure subroutine compute_solar_neutrino_fluxes(nr, xa, reactions_per_gram, e)
      use composition_layout, only: i_h1, i_he3
      use neutrino_flux_layout, only: i_nu_b8, i_nu_be7, i_nu_f17, i_nu_hep, &
           i_nu_n13, i_nu_o15, i_nu_pep, i_nu_pp, i_nu_he3he3, i_nu_he3he4
      type(nuclear_rates), intent(in) :: nr
      real(dp), intent(in) :: xa(:), reactions_per_gram(num_reactions)
      type(nuclear_energy_result), intent(inout) :: e
! 4 pi (1 AU)^2, cm^2
      real(dp), parameter :: four_pi_au2 = 2.812295d+27
      real(dp) :: eg(num_reactions)
      integer :: k

      e%fluxes_set = .true.
      if (xa(i_h1) <= 1.0d-6) then
         e%neutrino_flux = 0.0d0
         return
      end if
      eg = reactions_per_gram
      e%neutrino_flux(i_nu_pp) = eg(r_pp)/four_pi_au2
! pep from pp (Neutrino Astrophysics eq. 3.17); xs_pep scales the pep/pp
! ratio, so leave it at 1 unless that ratio is meant to change
      e%neutrino_flux(i_nu_pep) = (3.4848d-6)*nr%electron_density_na*nr%t9_m12* &
           (1.0d0 + 20.d0*nr%t9)*eg(r_pp)
      e%neutrino_flux(i_nu_pep) = e%neutrino_flux(i_nu_pep)*net_settings%cross_section_scale(xs_pep)/four_pi_au2
! hep (eq. 3.12), from the screened rate factor of compute_reaction_rates
      e%neutrino_flux(i_nu_hep) = nr%hep
      e%neutrino_flux(i_nu_hep) = e%neutrino_flux(i_nu_hep)*xa(i_h1)*xa(i_he3)/ &
           four_pi_au2
      e%be7_set = .true.
      e%be7_mass_fraction = eg(r_he3he4)/(nr%be7_proton + nr%be7_electron)
      e%neutrino_flux(i_nu_be7) = eg(r_he3he4)*nr%frac_be7_ecap/four_pi_au2
      e%neutrino_flux(i_nu_b8) = eg(r_he3he4)*nr%frac_be7_pcap/four_pi_au2
      e%neutrino_flux(i_nu_n13) = eg(r_pc12)/four_pi_au2
      e%neutrino_flux(i_nu_o15) = eg(r_pn14)/four_pi_au2
      e%neutrino_flux(i_nu_f17) = eg(r_po16)/four_pi_au2
      e%neutrino_flux(i_nu_he3he3) = eg(r_he3he3)/four_pi_au2
      e%neutrino_flux(i_nu_he3he4) = eg(r_he3he4)/four_pi_au2
! to 10^10 per cm^2 per s per g; non-positive values to zero (KC
! 2025-05-30: was a 1e-50 floor, which underflowed)
      do k = 1,10
         e%neutrino_flux(k) = (1.0d-10)*e%neutrino_flux(k)
         if (e%neutrino_flux(k) <= 0.0d0) e%neutrino_flux(k) = 0.0d0
      end do
end subroutine compute_solar_neutrino_fluxes

!----------------------------------------------------------------------
! The thermal neutrino losses (MHP 5/91; for neutrino-cooled cores of
! evolved stars), subtracted from eps with their absolute derivatives:
! Itoh et al. (1996) when net_settings%use_itoh_neutrino_loss (Grant
! Newsham 9/06), otherwise the Beaudet, Petrosian & Salpeter (1967) fits
! for the pair, photo and plasma losses.
subroutine add_thermal_neutrino_loss(log_density, log_temperature, t9, xa, e)
      use composition_layout, only: i_h1, i_he4, i_c12, i_c13, i_o16, i_o18
      use physical_constants, only: cc13, ln10
      use math_lib
      real(dp), intent(in) :: log_density, log_temperature, t9, xa(:)
      type(nuclear_energy_result), intent(inout) :: e
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
      el = t9/5.9302d0
      eli = 1.0d0/el
      ez = exp(cc13*ln10*(log_density-9.0d0))*eli*(0.7937d0+0.2063d0*xa(i_h1))
      ez3 = ez**3
      emue = 0.5d0*(1.0d0+xa(i_h1))
      ex1 = 0.0d0
      if (t9 >= 0.2d0) then
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
      if (t9 >= 0.2d0) then
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
! Sources and history (the header of engeb.f, kept with the physics it
! describes). Legacy names: EPP1/EPP2/EPP3/ECN/E3AL -> the chain energies
! (pp_chain_energy_gen, he3he4_be7_electron_energy_gen,
! he3he4_be7_proton_energy_gen, cno_cycle_energy_gen,
! triple_alpha_energy_gen); PEP/PET -> deps_dlnrho/deps_dlnt; SUM1 ->
! total_energy_gen_rate; EG(I) -> rates_times_reactants; X, Y, XHE3, ...
! -> the xa slots.
!======================================================================
! 3/92 DBG Added new neutrino loss calculation routines
! NOVEMBER 6, 1990 (JNB)
! THE FOLLOWING SUBROUTINE CALCULATES THE NUCLEAR REACTION RATES WITH
!   SPECIAL ATTENTION TO DETAIL REQUIRED FOR CALCULATING SOLAR NEUTRINO
!   FLUXES. THE NEUTRINO FLUXES ARE EVALUATED ALSO IN THIS SUBROUTINE.
! ALL PHYSICAL QUANTITIES ARE CGS. THE REACTION RATES ARE IN
!  GM^{-1}S^{-1}.
!  MHP 10/97
!
!   On 10/13/97, JNB converted the nuclear masses from neutral nuclear
!   masses to bare nuclear masses by subtracting Z(I)*(m_e)*c^2 from
!   the neutral nuclear masses. This caused changes in a number of
!   places: in ANUC(I), in Q1(I)-Q7(I), in the calculation of the Be7electron
!   Be7proton rates, and in the calculation of the N15p branching ratio.
!
!   JNB made some purely cosmetic changes on 1/20/96. Revised all of the
!   input Q's on 9/23-25/97 to agree with submitted version of Solar
!   Fusion Workshop paper. The SStandard are fixed to agree with
!   the Workshop paper. JNB recalculated all of the EG(I) to determine
!   the best values for the energy generation for all the reactions,
!   taking account of my improved calculations of neutrino energy loss.
!   The calculations are documented in Vol. 19, 132-141, 1997 of my notes.
!   The neutral atom mass differences are taken from Table of Isotopes,
!   8th Ed, 1996 and the neutrino energy losses from Bahcall, Gallium
!   solar neutrino experiments, Phys Rev C, in press, 1997.
!
!  All numbers in this subroutine have been calculated by John Bahcall
!  so that they agree with the modern numbers in Neutrino
!  Astrophysics (1989) or much more recent results, as indicated in
!  comments lines.
!
!  PREVIOUS VERSIONS OF THE YALE CODE
!  USED THE NUCLEAR REACTION RATES PRIMARILY FROM FOWLER, CAUGHLAN, AND
!  ZIMMERMAN IN ANN. REV. ASTR. AND ASTROPHYS. (1975) OR, FOR THE BE7 AND
!  B8 RATES, THE VALUES FROM BAHCALL AND SEARS (1972), ALSO. ANN. REV. ASTR.
!  AND ASTROPHYS. SOME VALUES WERE FROM HARRIS ET AL (1983, ANN. REV. ASTR.
!  IN A FEW CASES, NUMERICAL ERRORS (NOT DUE TO REVISIONS IN THE NUCLEAR DATA)
!  HAVE BEEN CORRECTED.  IN A SMALL NUMBER OF OTHER PLACES, THE PREVIOUS
!  PROGRAM CONTAINED PHYSICALLY INCORRECT STATEMENTS, WHICH HAVE BEEN CORRECTED
!  IN THE PRESENT VERSION. IN THE COURSE OF THE CHECKING AND REVISIONS, I
!  HAVE ADDED MANY COMMENT STATEMENTS IN ORDER TO MAKE THE PROGRAM MORE
!  TRANSPARENT AND SIMPLER TO REVISE.  SOME STATEMENTS HAVE BEEN ADDED
!  TO SIMPLIFY THE PROCESS OF REVISING THE SUBROUTINE AS IMPROVED DATA
!  BECOME AVAILABLE.  FOR EXAMPLE, I HAVE INSERTED EXPLICITLY
!  THE STANDARD CROSS SECTION FACTORS IN A WAY THAT IS EASY TO REVISE;
!  SSTANDARD(I) = 1.0 WHEN THE ITH CROSS SECTION HAS THE STANDARD VALUE
!  MEASUREMENTS USED IN NEUTRINO ASTROPHYSICS.
!  MHP 10/97
!
!  The PREVIOUS standard values of SStandard(I) used were those that
!   were given in Table 1 of Bahcall and Pinsonneault (1992).  This work
!   was  published in the Reviews of Modern Physics, 64,885, 1992. The
!   current version of SStandard(I) refers to Bahcall and
!   Pinsonneault(1995) , Rev. Mod. Phys. 67, 781 (1995), Table 1, if
!   changes have occured.  Otherwise, they are the same as in the
!   Bahcall-Pinsonnealt (1992) Rev. Mod. Phys. article.
!
!  Weakscreening is a parameter set in this subroutine. To obtain
!  the Graboske et al. and Salpeter standard results, use:
!  weakscreening = 0.03.  To investigate the effect of always using
!  weak screening, use a value for weakscreening greater than unity,
!  e. g., 30. The transition region between weak and strong screening
!  is defined by Graboske et al. as $\Lambda_{12} = 0.1$, see page
!  ApJ 181, 465 (1973) [5/15/97].
!
!  The value of SStandard(17) for hep is taken from Carlson et al (1991)
!  Phys. Rev. C 44, 619. It corresponds to an S sub 0 =
!  1.3 E-20 keV-barns, a factor of 0.1625  smaller than indicated by the
!  older measurements and calculations used in Neutrino Astrophysics.
!
! THE NUCLEAR ENERGY RELEASE TO THE STAR FROM EACH REACTION
!  IS TAKEN FROM BAHCALL AND ULRICH, RMP 60, 297 (1988) AND TAKES ACCOUNT
!  ACCURATELY OF NEUTRINO ENERGY LOSS.
!
! I HAVE ADDED A NEW SECTION AT THE END OF THE SUBROUTINE THAT
!  CALCULATES THE SOLAR NEUTRINO FLUXES AT THE EARTH.  THESE FLUXES
!  ARE IN THE UNITS OF CM^-2 SEC^-2 PER GM. TO GET THE FLUX
!  FROM A SHELL, MULTIPLY BY THE MASS OF THE SHELL IN UNITS OF
!  GRAMS.  THE FLUXES ARE STORED IN star%neutrino_flux.  I ALSO CALCULATE
!  THE FICTIONAL NEUTRINO FLUXES ASSOCIATED WITH THE HE3 + HE3 AND WITH
!  THE HE3 + HE4 RECTIONS; THESE FICTIONAL FLUXES ARE USEFUL DIAGNOSTICS
!  OF THE SOLAR MODEL.
!
! *****************************************
! CHANGING NUCLEAR REACTION CROSS SECTIONS.
! *****************************************
! NUCLEAR CROSS SECTIONS CAN BE CHANGED SIMPLY BY INSERTING NEW NUMBERS
!  FOR THE DATA VALUES OF SSTANDARD(I), WHICH ARE THE RATIOS OF THE
!  DESIRED CROSS SECTION FACTORS TO THE VALUES GIVEN IN TABLES 3.2 AND
!  3.4 OF NEUTRINO ASTROPHYSICS (1989). IF THE VALUE GIVEN IN NEUTRINO
!  ASTROPHYSICS IS USED, THEN SSTANDARD(I) = 1.0 .  TO INCREASE THE CROSS
!  SECTION FOR REACTION K BY A FACTOR OF TWO COMPARED TO THE STANDARD
!  VALUE, SET SSTANDARD(K) = 2.0 .  TO DETERMINE WHICH VALUE OF I GOES
!  WITH WHICH REACTION, SEE THE SECTION JUST BELOW.
!  THE ENERGY DERIVATIVES ENTER IN A FORM IN WHICH THEY ARE DIVIDED BY
!  THE ABSOLUTE VALUES OF THE CROSS SECTIONS AT ZERO ENERGY.  THUS IF
!  THE SHAPE OF THE CROSS SECTION EXTRAPOLATION IS UNCHANGED AND ONLY
!  THE INTERCEPT OF S(E) AT ZERO ENERGY IS CHANGED, THEN NO CORRECTION
!  NEED BE MADE FOR THE DERIVATIVES.  THEY ARE AUTOMATICALLY SCALED
!  CORRECTLY.  THE EXACT WAY THAT THE DERIVATIVES ENTER THE
!  CALCULATIONS IS DESCRIBED IN THE SECTION LABELED ``DEFINING THE
!  Q(I)'' THAT IS PRESENTED BELOW.
! ************************************
! IDENTIFYING THE REACTIONS.
! ************************************
!  THE VALUE OF J DENOTES WHICH OF THE REACTIONS THE COEFFICIENTS
!  REFER TO:
!  J = 1, PP; J = 2, HE3+HE3; J = 3, HE3+ HE4; J =4, P + C12;  J = 5, P+C13;
!  J = 6. P + N14; J = 7, P + O16.
!  REACTIONS J = 8, 13 ARE NOT RELEVANT FOR THE SOLAR INTERIOR; THEY ARE
!   HOLDOVERS FROM THE EARLIER YALE CODE.
!  REACTION 14 IS PEP; REACTION 15 IS BE7 ELECTRON CAPTURE; REACTION 16 IS
!   BE7 PROTON CAPTURE; REACTION 17 IS THE HEP REACTION.
!  do not change sstandard(14) unless you want to change the ratio of
!  pep to pp.
!   REACTIONS 14-17 WERE NOT EXPLICITLY INCLUDED IN THE YALE
!    PREVIOUS VERSION OF THE CODE, BUT THEY ARE IMPORTANT FOR
!    SOLAR NEUTRINO CALCULATIONS.
!   THE BRANCHING OF THE N15 + P REACTIONS IS TREATED IN A SERIES OF
!    SEPARATE STATEMENTS FOLLOWING THE CALCULATION OF THE BE7 + P
!    REACTION. SEE THE DEFINITIONS OF F3 AND F4.  IF THE CROSS-SECTION
!    FACTORS OF THE N15 + P REACTIONS ARE REVISED, THEN THE NUMERICAL
!    COEFFICIENTS MUST BE CHANGED IN THE DEFINITION OF C12ALPHA AND
!    O16GAMMA.
! FOR Q1(I), ...,Q(5(I), I = 8 CORRESPONDS TO THE BE7 +  P REACTION.
!  THIS ASSIGNMENT FOR I = 8 IS ONLY VALID FOR THE LISTED Q'S AND NOT
!  FOR OTHER ARRAYS IN THE PROGRAM.
!  IU IS THE SHELL NUMBER.
