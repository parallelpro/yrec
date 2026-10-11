!----------------------------------------------------------------------
! reaction_rates (formerly the body of net_lib)
!----------------------------------------------------------------------
! The rate physics of the network: pure functions of (log rho, log T,
! composition) and the run's settings.
!     compute_reaction_rates   the one rate evaluation of a zone: 13
!                              reaction rates, their T/rho derivatives,
!                              the branchings
!     compute_network_rates    those rates repacked as the composition
!                              rate vector (per s per amu) for the
!                              network burn
!     deuterium_proton_capture_rate   the shared d(p,g) fit
!     set_net_settings         installs the run's settings (net_settings)
! It also holds the index families (r_* reaction slots, rv_* fraction
! slots, xs_* cross-section slots, iso_* isotopes), the fit and
! screening data, and the settings record. The rest of net/ builds on
! it; code outside net/ reaches it through net_lib.
!
! (History: these routines were once separate files, nuclear/<name>.f90,
! concatenated into net_lib with their per-file conversion notes kept
! below. lir and ratext moved out to numerics_lib in 2026 -- generic
! numerical utilities filed here only because the Li/Be burn called
! them.)
! ---------------------------------------------------------------------
module reaction_rates
      implicit none
      private

! The reaction slots: the index of every per-reaction array (rates,
! derivatives, screening, energies, cross_section_scale(1:13)). Slots 9
! (O16 + alpha) and 13 (C12 + C12) are carried but never computed.
      integer, parameter, public :: &
           r_pp = 1, r_he3he3 = 2, r_he3he4 = 3, r_pc12 = 4, r_pc13 = 5, &
           r_pn14 = 6, r_po16 = 7, r_ac13 = 8, r_ao16_unused = 9, &
           r_ac12 = 10, r_an14 = 11, r_3alpha = 12, r_c12c12_unused = 13, &
           num_reactions = 13
! The composition rate vector (compute_network_rates): the reactions at
! their r_* slots, per s per amu, then the N15(p,alpha) and Be7
! electron-capture fractions.
      integer, parameter, public :: num_rate_rows = 15
      integer, parameter, public :: rv_frac_n15_pa = 14, rv_frac_be7_ecap = 15
! cross_section_scale(14:17): the four cross-section factors after the
! reactions -- pep, the Be7 electron and proton captures, hep.
      integer, parameter, public :: &
           xs_pep = 14, xs_be7e = 15, xs_be7p = 16, xs_hep = 17
! In fit_q1..fit_q5, qs0e_scale and qqs0ee_scale, slot 8 holds Be7 + p
! (in these arrays only; reaction slot 8 is He4 + C13).
      integer, parameter, public :: iq_be7p = 8
! The isotopes of the screening composition and of the mass and charge
! tables.
      integer, parameter, public :: &
           iso_n = 1, iso_h1 = 2, iso_h2 = 3, iso_h3 = 4, iso_he3 = 5, &
           iso_he4 = 6, iso_c12 = 7, iso_c13 = 8, iso_n14 = 9, &
           iso_o16 = 10, iso_o18 = 11, iso_ne20 = 12, iso_mg24 = 13, &
           num_isotopes = 13
! S0 (keV-b) of the two N15 + p branches (compute_n15_branching); integers,
! as the original literals were (both exact when promoted).
      integer, parameter, public :: s0_n15pg_kevb = 64, s0_n15pa_kevb = 67500

! d(p,g)3He: Avogadro's number over the deuteron mass (amu) turns CF88's
! N_A<sigma v> into the per-gram rate of both the deuterium depletion and
! its energy release, and sets the convective-overturn cap. (The two
! former copies differed by 0.68%: one divided by A = 2.)
      double precision, parameter, public :: avogadro_dpg = 6.023d23
      double precision, parameter, public :: deuteron_mass_amu = 2.013553d0
! S-factor fits of reactions 1-7 and of Be7(p,g) (slot iq_be7p of q1..q5):
! FCZ67 eq. 53 / Neutrino Astrophysics eq. 3.14 times T9^(-2/3); q6 the
! tau coefficient, q7 the log prefactor, q8 the Harris et al. T9^2 term
! (history of the values: the end of this file).
      double precision, parameter, private :: &
           fit_q1(iq_be7p) = [0.12317d0, .03392d0, .0325d0, .0304d0, .03035d0, .0273d0, .02494d0, .040572d0], &
           fit_q2(iq_be7p) = [1.08749d0, -.273d0, -.2085d0, .7630d0, -0.4044d0, -1.60d0, -1.224d0, -0.2095d0], &
           fit_q3(iq_be7p) = [.93833d0, -.0648d0, -.0474d0, .1626d0, -.08598d0, -.3064d0, -.2139d0, -0.0595d0], &
           fit_q4(iq_be7p) = [0.d0, 0.d0, 0.d0, 4.79d0, 7.456d0, 0.0d0, .69703d0, 0.16762d0], &
           fit_q5(iq_be7p) = [0.d0, 0.d0, 0.d0, 2.595d0, 4.032d0, 0.0d0, 0.3097d0, 0.12114d0], &
           fit_q6(r_po16) = [-3.3804d0, -12.2757d0, -12.826d0, -13.6899d0, -13.7173d0, -15.2281d0, -16.6925d0], &
           fit_q7(r_po16) = [20.8964d0, 76.6003d0, 67.8036d0, 69.130d0, 70.3809d0, 69.8517d0, 70.8012d0], &
           fit_q8(r_po16) = [0.d0, 0.d0, 0.d0, 0.0d0, 0.0d0, 0.0d0, 0.d0]
! screening: charge products and the Graboske et al. (1973) table-4
! factors of each reaction slot
      double precision, parameter, private :: &
           screen_zprod(num_reactions) = [1.d0, 4.d0, 4.d0, 6.d0, 6.d0, 7.d0, 8.d0, 12.d0, 16.d0, 12.d0, 14.d0, 12.d0, 36.d0], &
           screen_z53(num_reactions) = [1.175d0, 3.73d0, 3.73d0, 4.804d0, 4.804d0, 5.385d0, 5.941d0, 9.014d0, 11.24d0, &
                9.014d0, 10.15d0, 9.104d0, 23.28d0], &
           screen_z43(num_reactions) = [0.52d0, 1.31d0, 1.31d0, 1.488d0, 1.488d0, 1.61d0, 1.721d0, 2.577d0, 3.025d0, &
                2.577d0, 2.81d0, 2.577d0, 5.668d0], &
           screen_z23(num_reactions) = [-0.413d0, -0.655d0, -0.655d0, -0.643d0, -0.643d0, -0.659d0, -0.673d0, -0.889d0, &
                -0.946d0, -0.889d0, -0.92d0, -0.889d0, -1.36d0], &
           screen_z86(num_reactions) = [1.630d0, 5.917d0, 5.917d0, 8.302d0, 8.302d0, 9.520d0, 10.716d0, 16.192d0, &
                20.978d0, 16.192d0, 18.606d0, 16.192d0, 45.6635d0]
! bare nuclear masses (amu) and charges of the iso_* isotopes
      double precision, parameter, public :: &
           isotope_mass(num_isotopes) = [1.008665d0, 1.007276d0, 2.013553d0, 3.015501d0, 3.014933d0, 4.001506d0, &
                11.996709d0, 13.000064d0, 13.999233d0, 15.990526d0, 17.994772d0, 19.986954d0, 23.978458d0], &
           isotope_charge(num_isotopes) = [0.d0, 1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 6.d0, 6.d0, 7.d0, 8.d0, 8.d0, 10.d0, 12.d0]

! The rates of one zone, from compute_reaction_rates: rate(r_*) per gram
! per second with the abundance factors omitted (screened, scaled), and
! its logarithmic derivatives; the Be7 and N15 branchings; the two Be7
! channel rates and the screened hep factor (flux diagnostics); and the
! zone values the energy and flux consumers reuse. All zero below
! net_settings%logT_cutoffs(1).
      type, public :: nuclear_rates
         double precision :: rate(num_reactions) = 0.0d0
         double precision :: dlnrate_dlnt(num_reactions) = 0.0d0
         double precision :: dlnrate_dlnrho(num_reactions) = 0.0d0
         double precision :: frac_be7_ecap = 0.0d0, frac_be7_pcap = 0.0d0
         double precision :: frac_n15_pa = 0.0d0, frac_n15_pg = 0.0d0
         double precision :: be7_electron = 0.0d0, be7_proton = 0.0d0
         double precision :: hep = 0.0d0
         double precision :: density = 0.0d0, t9 = 0.0d0, t9_p13 = 0.0d0, &
              t9_p23 = 0.0d0, t9_m13 = 0.0d0, t9_m23 = 0.0d0, t9_m12 = 0.0d0
         double precision :: electron_density_na = 0.0d0
      end type nuclear_rates

! T9 = T/10^9 K and the powers of it the rate formulas use (p13 = T9^1/3,
! m23 = T9^-2/3, m12 = T9^-1/2, ...)
      type :: t9_powers
         double precision :: t9 = 0.0d0, p13 = 0.0d0, p23 = 0.0d0, &
              m13 = 0.0d0, m23 = 0.0d0, m1 = 0.0d0, m2 = 0.0d0, &
              m12 = 0.0d0, m32 = 0.0d0
      end type t9_powers

! The plasma quantities the screening formulas share: the composition
! moments (1/mu_i, 1/mu_e, sum Z^1.58 X/A, sum Z^2 X/A), Graboske et
! al.'s Lambda_0 and its powers, DeWitt et al.'s z_tilde, the mean ionic
! charge and its powers, and the weak-screening exponent per unit charge
! product.
      type :: screening_plasma
         double precision :: inv_mu_i = 0.0d0, inv_mu_e = 0.0d0, &
              charge_moment_158 = 0.0d0, charge_moment_2 = 0.0d0
         double precision :: lambda0 = 0.0d0, lambda0_23 = 0.0d0, &
              lambda0_086 = 0.0d0, z_tilde = 0.0d0, mean_charge = 0.0d0, &
              z_tilde_058 = 0.0d0, mean_charge_028 = 0.0d0, &
              mean_charge_13 = 0.0d0, weak_unit = 0.0d0
      end type screening_plasma

! ---- the network's settings (2026-10) ----
! Everything net/ needs from the run's configuration, in one record. The
! star layer fills it at the start of every run, after the controls are
! read and the Monte Carlo rate factors applied (configure_network,
! core/burn_step.f90), by calling set_net_settings; nothing else writes it.
      integer, parameter, public :: num_cross_section_scales = xs_hep
      type, public :: net_settings_t
! log10 T at or below which: (1) nothing burns; (2), (3) parts of the
! equilibrium burn are off; (4) helium burning is off; (5) the solar
! neutrino fluxes and thermal neutrino losses are skipped
         double precision :: logT_cutoffs(5) = 0.0d0
! the screening parameter below which screening is weak
         double precision :: weak_screening_threshold = 0.0d0
! multipliers on the S-factors (r_* slots 1-13, then xs_pep..xs_hep), on
! their first and second energy derivatives (slots 1-7 and iq_be7p),
! and on the two N15 + p branches
         double precision :: cross_section_scale(num_cross_section_scales) = 1.0d0
         double precision :: qs0e_scale(iq_be7p) = 1.0d0, qqs0ee_scale(iq_be7p) = 1.0d0
         double precision :: o16_gamma_scale = 1.0d0, c12_alpha_scale = 1.0d0
! Li6(p,a), Li7(p,a), Be9(p,g), Be9(p,d), Be9(p,a) rate multipliers
         double precision :: li6_rate_scale = 1.0d0, li7_rate_scale = 1.0d0, &
              be9_pg_rate_scale = 1.0d0, be9_pd_rate_scale = 1.0d0, &
              be9_palpha_rate_scale = 1.0d0
! compute the solar neutrino fluxes; use the Itoh et al. (1996) thermal
! neutrino losses
         logical :: calc_neutrinos = .false., use_itoh_neutrino_loss = .false.
      end type net_settings_t
      type(net_settings_t), public, protected :: net_settings
      public :: set_net_settings, screening_mass_fractions, &
           deuterium_proton_capture_rate, compute_reaction_rates, &
           compute_network_rates

contains

! set_net_settings: install the run's settings (net_settings_t above).
subroutine set_net_settings(settings)
      type(net_settings_t), intent(in) :: settings
      net_settings = settings
end subroutine set_net_settings


!----------------------------------------------------------------------
! deuterium_proton_capture_rate
!----------------------------------------------------------------------
! The d(p,g)3He rate excluding the abundance factors, rho N_A<sigma v> /
! m_d (Caughlan & Fowler 1988, the T9**(-2/3) form), for the deuterium
! depletion (light_elements) and its energy (nuclear_energy) alike.
! Unscreened: the burning temperatures are low.
double precision function deuterium_proton_capture_rate(density, t9, t9_p13, t9_p23, t9_m13, t9_m23)
      use math_lib
      implicit none
      double precision, intent(in) :: density, t9, t9_p13, t9_p23, &
           t9_m13, t9_m23
      double precision :: zz, tfacdeut

      zz = -3.72d0*t9_m13
      tfacdeut = 1.0d0+0.112d0*t9_p13+3.38d0*t9_p23+2.65d0*t9
      deuterium_proton_capture_rate = density*2.240d3*t9_m23*exp(zz)*tfacdeut*avogadro_dpg/ &
           deuteron_mass_amu
end function deuterium_proton_capture_rate

!----------------------------------------------------------------------
! compute_reaction_rates: the H/He reaction rates of one zone and their
! logarithmic derivatives with respect to T and rho -- the one rate
! evaluation behind both the composition rate vector
! (compute_network_rates) and the energy generation (nuclear_energy).
! mass_frac(iso_*) is the composition the screening sees
! (screening_mass_fractions). nr%rate(r_*) is per gram per second with the
! abundance factors omitted, screened and scaled by net_settings. Each
! consumer applies its own floor afterwards. Below logT_cutoffs(1)
! everything is zero; without hydrogen, reactions 1-7 and the branchings
! are; the alpha captures need logT_cutoffs(4). Sources: the end of this
! file.
subroutine compute_reaction_rates(log_density, log_temperature, mass_frac, nr)
      use physical_constants, only: ln10
      use math_lib
      double precision, intent(in) :: log_density, log_temperature, &
           mass_frac(num_isotopes)
      type(nuclear_rates), intent(out) :: nr
      type(t9_powers) :: t
      type(screening_plasma) :: plasma
! the screening exponent of each reaction and its d/d ln rho, d/d ln T
      double precision :: screen(num_reactions), dscreen_dlnrho(num_reactions), &
           dscreen_dlnt(num_reactions)
      integer :: i, first_zeroed_rxn

      call compute_composition_moments(mass_frac, plasma)
      nr%electron_density_na = ( exp(ln10*log_density) )*plasma%inv_mu_e
      if (log_temperature.le.net_settings%logT_cutoffs(1)) return
      nr%density = exp(ln10*log_density)
      t = powers_of_t9(log_temperature)
      nr%t9 = t%t9
      nr%t9_p13 = t%p13
      nr%t9_p23 = t%p23
      nr%t9_m13 = t%m13
      nr%t9_m23 = t%m23
      nr%t9_m12 = t%m12
      call compute_screening(nr%density, nr%electron_density_na, t, plasma, &
           screen, dscreen_dlnrho, dscreen_dlnt)
! reactions 1-7 and the two branchings need hydrogen; the slots from
! first_zeroed_rxn on are zeroed (only the alpha captures when hydrogen
! burns, all 13 when it is gone)
      first_zeroed_rxn = r_pp
      if (mass_frac(iso_h1).ne.0.0d0) then
         first_zeroed_rxn = r_ac13
         call compute_hydrogen_burning_rates(nr%density, t, screen, &
              dscreen_dlnrho, dscreen_dlnt, nr)
         call compute_be7_branching(nr%density, t, mass_frac(iso_h1), plasma, nr)
         call compute_n15_branching(t, nr)
      end if
      do i=first_zeroed_rxn,num_reactions
         nr%rate(i)=0.d0
         nr%dlnrate_dlnrho(i)=0.d0
         nr%dlnrate_dlnt(i)=0.d0
      end do
      if (log_temperature.ge.net_settings%logT_cutoffs(4)) &
           call compute_helium_burning_rates(nr%density, t, screen, &
           dscreen_dlnrho, dscreen_dlnt, nr)
      call compute_hep_rate(nr%density, t, plasma, nr)
end subroutine compute_reaction_rates

!----------------------------------------------------------------------
! The composition moments the screening needs: 1/mu_i = sum X/A,
! 1/mu_e = sum Z X/A, sum Z^1.58 X/A (intermediate screening, Graboske et
! al.: the mean of Z^(3b-1)) and sum Z^2 X/A (Salpeter's zeta).
subroutine compute_composition_moments(mass_frac, plasma)
      use math_lib
      double precision, intent(in) :: mass_frac(num_isotopes)
      type(screening_plasma), intent(out) :: plasma
      double precision :: per_amu
      integer :: i

      plasma%inv_mu_i = 0.d0
      plasma%inv_mu_e = 0.d0
      plasma%charge_moment_158 = 0.d0
      plasma%charge_moment_2 = 0.d0
      do i = 1,num_isotopes
         per_amu = mass_frac(i)/isotope_mass(i)
         plasma%inv_mu_i = plasma%inv_mu_i+per_amu
         plasma%inv_mu_e = plasma%inv_mu_e+ &
              per_amu*isotope_charge(i)
         plasma%charge_moment_158 = plasma%charge_moment_158+per_amu*pow(isotope_charge(i), 1.58d0)
         plasma%charge_moment_2 = plasma%charge_moment_2+per_amu*isotope_charge(i)**2
      end do
end subroutine compute_composition_moments

!----------------------------------------------------------------------
! T9 = T/10^9 K and the powers of it the rate formulas use.
function powers_of_t9(log_temperature) result(t)
      use physical_constants, only: ln10, cc13
      use math_lib
      double precision, intent(in) :: log_temperature
      type(t9_powers) :: t

      t%t9 = exp(ln10*(log_temperature - 9.0d0))
      t%p13 = pow(t%t9, cc13)
      t%p23 = t%p13**2
      t%m13=1.d0/t%p13
      t%m23=t%m13**2
      t%m1=1.d0/t%t9
      t%m2=t%m1**2
      t%m12=1.d0/sqrt(t%t9)
      t%m32=t%m1*t%m12
end function powers_of_t9

!----------------------------------------------------------------------
! The screening exponent of each reaction slot and its derivatives:
! weak (Salpeter 1954, with his degeneracy correction f'/f) below
! net_settings%weak_screening_threshold (0.03 reproduces Graboske et al.
! and Salpeter; a value above 1 forces weak screening in the Sun),
! intermediate (Graboske et al. 1973, table 4) up to a weak-screening
! value of 2, and beyond that strong or intermediate, whichever is
! smaller (strong always from 5). Also fills the plasma quantities the
! single-reaction screening of Be7 + p and hep reuses. (Strong screening
! was never checked: it does not occur in the Sun.)
subroutine compute_screening(density, electron_density_na, t, plasma, &
     screen, dscreen_dlnrho, dscreen_dlnt)
      use physical_constants, only: cc13
      use math_lib
      double precision, intent(in) :: density, electron_density_na
      type(t9_powers), intent(in) :: t
      type(screening_plasma), intent(inout) :: plasma
      double precision, intent(out) :: screen(num_reactions), &
           dscreen_dlnrho(num_reactions), dscreen_dlnt(num_reactions)
! (p_F/m_e c)^2, E_F/kT, Salpeter's f'/f and log10(E_F/kT)
      double precision :: fermi_momentum_sq, fermi_energy_over_kt, &
           degeneracy_factor, log_degeneracy
      double precision :: weak, intermediate, strong
      integer :: i

      fermi_momentum_sq=1.017677d-4*pow(electron_density_na, (2.0d0/3.0d0))
      fermi_energy_over_kt=5.92986d0*t%m1*(sqrt(1.d0+fermi_momentum_sq)-1.d0)
      if (fermi_energy_over_kt.le.1.d-2) then
         degeneracy_factor=1.0d0
      else
! Salpeter's figure 1, as a function of the degeneracy
         log_degeneracy=log10(fermi_energy_over_kt)
         if (log_degeneracy.ge.1.5d0) then
            degeneracy_factor=0.0d0
         else
            degeneracy_factor=0.75793d0-0.54621d0*log_degeneracy-0.30964d0*log_degeneracy**2+ &
            0.12535d0*log_degeneracy**3+0.1203d0*log_degeneracy**4- &
            0.012857d0*log_degeneracy**5-0.014768d0*log_degeneracy**6
         end if
      end if
! Lambda_0 of Graboske et al. and its 2/3 and 0.86 powers; z_tilde of
! DeWitt et al. (Salpeter's zeta without the 1/mu_i), the mean ionic
! charge, and the powers of them the intermediate and strong forms use
      plasma%lambda0=5.9426d-6*t%m32*sqrt(density*plasma%inv_mu_i)
      plasma%lambda0_23=pow(plasma%lambda0, (2.0d0/3.0d0))
      plasma%lambda0_086=pow(plasma%lambda0, 0.86d0)
      plasma%z_tilde=sqrt((plasma%charge_moment_2+degeneracy_factor*plasma%inv_mu_e)/ &
           plasma%inv_mu_i)
      plasma%mean_charge=plasma%inv_mu_e/plasma%inv_mu_i
      plasma%z_tilde_058=pow(plasma%z_tilde, 0.58d0)
      plasma%mean_charge_028=pow(plasma%mean_charge, 0.28d0)
      plasma%mean_charge_13=pow(plasma%mean_charge, cc13)
! weak screening per unit charge product: Z1 Z2 Lambda_0 z_tilde
      plasma%weak_unit=plasma%lambda0*plasma%z_tilde
      do i=1,num_reactions
         weak=plasma%weak_unit*screen_zprod(i)
         if (weak.le.net_settings%weak_screening_threshold) then
            screen(i)=weak
            dscreen_dlnrho(i)=0.5d0*weak
            dscreen_dlnt(i)=-1.5d0*weak
         else
            intermediate=0.38d0*plasma%lambda0_086*plasma%charge_moment_158*screen_z86(i)/ &
                 (plasma%inv_mu_i*plasma%z_tilde_058*plasma%mean_charge_028)
            if (weak.le.2.d0) then
               screen(i)=intermediate
               dscreen_dlnrho(i)=0.43d0*intermediate
               dscreen_dlnt(i)=-1.29d0*intermediate
            else
               strong=0.624d0*plasma%mean_charge_13*plasma%lambda0_23*(screen_z53(i)+ &
                    0.316d0*plasma%mean_charge_13*screen_z43(i)+0.737d0* &
                    screen_z23(i)/(plasma%mean_charge*plasma%lambda0_23))
               if (strong.lt.intermediate.or.weak.ge.5.d0) then
                  screen(i)=strong
                  dscreen_dlnrho(i)=0.208d0*plasma%mean_charge_13*(screen_z53(i)+ &
                       0.316d0*plasma%mean_charge_13*screen_z43(i))*plasma%lambda0_23
                  dscreen_dlnt(i)=-3.d0*dscreen_dlnrho(i)
               else
                  screen(i)=intermediate
                  dscreen_dlnrho(i)=0.43d0*intermediate
                  dscreen_dlnt(i)=-1.29d0*intermediate
               end if
            end if
         end if
      end do
end subroutine compute_screening

!----------------------------------------------------------------------
! Weak screening, or intermediate above the weak threshold, for one
! reaction of charge product zprod and Graboske z86 (Be7 + p, hep).
double precision function weak_intermediate_screening(plasma, zprod, z86)
      type(screening_plasma), intent(in) :: plasma
      double precision, intent(in) :: zprod, z86
      double precision :: weak

      weak = plasma%weak_unit*zprod
      if (weak.le.net_settings%weak_screening_threshold) then
         weak_intermediate_screening = weak
      else
         weak_intermediate_screening = 0.38d0*plasma%lambda0_086*plasma%charge_moment_158*z86/ &
              (plasma%inv_mu_i*plasma%z_tilde_058*plasma%mean_charge_028)
      end if
end function weak_intermediate_screening

!----------------------------------------------------------------------
! Reactions 1-7 (pp, He3+He3, He3+He4, p+C12, p+C13, p+N14, p+O16):
! rate = rho * [T9^(-2/3) S_eff/S0] * exp(tau + q7 + (q8 T9)^2 + screening)
! (Neutrino Astrophysics eq. 3.14, Fowler, Caughlan & Zimmerman 1967 eq.
! 52), scaled by cross_section_scale, with qs0e_scale and qqs0ee_scale on
! the first and second energy derivatives of S, and the exact
! T-derivative of that form (2026-10: the q3 term lacked qs0e_scale and
! the q4, q5 terms qqs0ee_scale).
subroutine compute_hydrogen_burning_rates(density, t, screen, dscreen_dlnrho, &
     dscreen_dlnt, nr)
      use math_lib
      double precision, intent(in) :: density, screen(num_reactions), &
           dscreen_dlnrho(num_reactions), dscreen_dlnt(num_reactions)
      type(t9_powers), intent(in) :: t
      type(nuclear_rates), intent(inout) :: nr
! T9^(-2/3) S_eff/S0
      double precision :: s_ratio
      integer :: i

      do i=r_pp,r_po16
         s_ratio=t%m23+fit_q1(i)*t%m13+net_settings%qs0e_scale(i)*(fit_q2(i)+fit_q3(i)*t%p13)+ &
              net_settings%qqs0ee_scale(i)*(fit_q4(i)*t%p23+fit_q5(i)*t%t9)
         nr%rate(i)=density*s_ratio*exp(fit_q6(i)*t%m13+fit_q7(i)+(fit_q8(i)*t%t9)**2+ &
              screen(i))
         nr%rate(i) = nr%rate(i)*net_settings%cross_section_scale(i)
         nr%dlnrate_dlnrho(i)=1.d0+dscreen_dlnrho(i)
         nr%dlnrate_dlnt(i)=dscreen_dlnt(i)-(fit_q6(i)*t%m13+(2.d0*t%m23+ &
         fit_q1(i)*t%m13-net_settings%qs0e_scale(i)*fit_q3(i)*t%p13- &
         net_settings%qqs0ee_scale(i)*(2.d0*fit_q4(i)*t%p23+3.d0*fit_q5(i)*t%t9))/s_ratio)/3.d0+ &
         2.d0*(fit_q8(i)*t%t9)**2
      end do
end subroutine compute_hydrogen_burning_rates

!----------------------------------------------------------------------
! The fate of Be7: electron capture (Neutrino Astrophysics eq. 3.18) or
! proton capture (table 3.2, eq. 3.12, with the S_eff/S0 corrections of
! slot iq_be7p and weak/intermediate screening), both per gram (rho over
! the bare Be7 mass, JNB 10/14/97), and the fraction burned each way.
subroutine compute_be7_branching(density, t, hydrogen_fraction, plasma, nr)
      use math_lib
      double precision, intent(in) :: density, hydrogen_fraction
      type(t9_powers), intent(in) :: t
      type(screening_plasma), intent(in) :: plasma
      type(nuclear_rates), intent(inout) :: nr
      double precision :: tau_be7p, s_ratio_be7p, per_gram_be7

      nr%be7_electron = (1.752d-10)*t%m12*(1.0d0 + 0.004d0*(1000.d0*t%t9 - 16.d0))
      nr%be7_electron = nr%be7_electron*plasma%inv_mu_e*net_settings%cross_section_scale(xs_be7e)
      tau_be7p = (-10.2625d0*t%m13)
      nr%be7_proton = (3.128813d+5)*hydrogen_fraction*net_settings%cross_section_scale(xs_be7p)*exp(tau_be7p)
      s_ratio_be7p = t%m23 + fit_q1(iq_be7p)*t%m13 + net_settings%qs0e_scale(iq_be7p)*(fit_q2(iq_be7p)+ fit_q3(iq_be7p)*t%p13) &
           + net_settings%qqs0ee_scale(iq_be7p)*(fit_q4(iq_be7p)*t%p23 + fit_q5(iq_be7p)*t%t9)
      nr%be7_proton = nr%be7_proton*s_ratio_be7p
      nr%be7_proton = nr%be7_proton*exp(weak_intermediate_screening(plasma, 4.0d0, 5.7790d0))
      per_gram_be7 = density*8.584981d+22
      nr%be7_proton = per_gram_be7*nr%be7_proton
      nr%be7_electron = per_gram_be7*nr%be7_electron
      nr%frac_be7_ecap = nr%be7_electron/(nr%be7_electron + nr%be7_proton)
      nr%frac_be7_pcap = nr%be7_proton/(nr%be7_electron + nr%be7_proton)
end subroutine compute_be7_branching

!----------------------------------------------------------------------
! The fate of N15 + p: (p,alpha)C12 or (p,gamma)O16, in the ratio of
! their effective S-factors (S0 in keV-b times the S_eff/S0 polynomial;
! Neutrino Astrophysics table 3.4, S0(p,alpha) = 67500 keV-b since the
! Solar Fusion workshop, JNB 9/25/97), each with its scale (MHP 9/14).
subroutine compute_n15_branching(t, nr)
      type(t9_powers), intent(in) :: t
      type(nuclear_rates), intent(inout) :: nr
      double precision :: pg_strength, pa_strength

      pg_strength = t%m23 + 0.0273016d0*t%m13 + 0.14374d0 + 0.027490d0*t%p13 &
                + 6.14685d0*t%p23 + 2.98940d0*t%t9
      pg_strength = pg_strength*s0_n15pg_kevb*net_settings%o16_gamma_scale
      pa_strength = t%m23 + 0.0273016d0*t%m13 + 2.01186d0 + 0.384763d0*t%p13 &
                + 17.0579d0*t%p23 + 8.29580d0*t%t9
      pa_strength = pa_strength*s0_n15pa_kevb*net_settings%c12_alpha_scale
      nr%frac_n15_pa = pa_strength/(pa_strength + pg_strength)
      nr%frac_n15_pg = 1.0d0 - nr%frac_n15_pa
end subroutine compute_n15_branching

!----------------------------------------------------------------------
! The alpha captures C13(a,n)O16, C12(a,g)O16, N14(a,g)F18(e+ nu)O18 and
! triple alpha (Caughlan & Fowler 1988, MHP 3/91; the prefactors from
! VandenBerg's notes), with their exact derivatives. In each, a1..a5 are
! the terms of the CF88 fit and r1, r2 its polynomial factors.
subroutine compute_helium_burning_rates(density, t, screen, dscreen_dlnrho, &
     dscreen_dlnt, nr)
      use physical_constants, only: cc13, cc23
      use math_lib
      double precision, intent(in) :: density, screen(num_reactions), &
           dscreen_dlnrho(num_reactions), dscreen_dlnt(num_reactions)
      type(t9_powers), intent(in) :: t
      type(nuclear_rates), intent(inout) :: nr
      double precision :: r1, r2, a1, a2, a3, a4, a5, dr1, da1, s_sum

! C13(alpha,n)O16
      r1=t%m23+0.0129d0*t%m13+2.04d0+0.184d0*t%p13
      a1 = 6.77d15*exp(-32.329d0*t%m13-(t%t9/1.284d0)**2)
      a2 = 3.82d5*exp(-9.373d0*t%m1)
      a3 = 1.41d6*exp(-11.873d0*t%m1)
      a4 = 2.00d9*exp(-20.409d0*t%m1)
      a5 = 2.92d9*exp(-29.283d0*t%m1)
      s_sum = a1*r1+t%m32*(a2+a3+a4+a5)
      nr%rate(r_ac13) = 1.157126d22*density*exp(screen(r_ac13))*s_sum
      nr%dlnrate_dlnrho(r_ac13)=1.0d0+dscreen_dlnrho(r_ac13)
      dr1 = cc13*(-2.0d0*t%m23-0.0129d0*t%m13+0.184d0*t%p13)
      da1 = a1*(cc13*32.329d0*t%m13 - 2.0d0*(t%t9/1.284d0)**2)
      nr%dlnrate_dlnt(r_ac13) = dscreen_dlnt(r_ac13)+(dr1*a1 + r1*da1 + &
           t%m32*(a2*(9.373d0*t%m1-1.5d0)+a3*(11.873d0*t%m1-1.5d0)+ &
           a4*(20.409d0*t%m1-1.5d0)+a5*(29.283d0*t%m1-1.5d0)))/s_sum
! C12(alpha,gamma)O16
      r1 = 1.0d0/(1.0d0+0.0489d0*t%m23)
      r2 = 1.0d0/(1.0d0+0.2654d0*t%m23)
      a1 = t%m2*exp(-32.120d0*t%m13)
      a2 = 1.04d8*r1**2*exp(-(t%t9/3.496d0)**2)
      a3 = 1.76d8*r2**2
      a4 = 1.25d3*t%m32*exp(-27.499d0*t%m1)
      a5 = 1.43d-2*t%t9**5*exp(-15.541d0*t%m1)
      s_sum = a1*(a2+a3)+a4+a5
      nr%rate(r_ac12) = 1.25388d22*density*exp(screen(r_ac12))*s_sum
      nr%dlnrate_dlnrho(r_ac12) = 1.0d0+dscreen_dlnrho(r_ac12)
      nr%dlnrate_dlnt(r_ac12) = dscreen_dlnt(r_ac12)+ &
           (a1*((cc13*32.120d0*t%m13-2.0d0)*(a2+a3)+ &
           a2*(cc13*0.1956d0*t%m23*r1-2.0d0*(t%t9/3.496d0)**2)+ &
           a3*(cc13*1.0616d0*t%m23*r2))+ &
           a4*(27.499d0*t%m1-1.5d0)+a5*(5.0d0+15.541d0*t%m1))/s_sum
! N14(alpha,gamma)F18(e+ nu)O18; linear coefficient 0.117 (2026 bugsweep
! Batch 3: both former copies carried 0.177)
      r1 = t%m23+0.012d0*t%m13+1.45d0+0.117d0*t%p13+1.97d0*t%p23 &
           +0.406d0*t%t9
      a1 = 7.78d9*exp(-36.031d0*t%m13-(t%t9/0.881d0)**2)
      a2 = t%m32*2.36d-10*exp(-2.798d0*t%m1)
      a3 = t%m32*2.03d0*exp(-5.054d0*t%m1)
      a4 = t%m23*1.15d4*exp(-12.310d0*t%m1)
      s_sum = a1*r1+a2+a3+a4
      nr%rate(r_an14)= 1.07452d22*density*exp(screen(r_an14))*s_sum
      nr%dlnrate_dlnrho(r_an14)=1.d0+dscreen_dlnrho(r_an14)
      dr1 = cc13*(-2.0d0*t%m23-0.012d0*t%m13+0.117d0*t%p13+ &
            3.94d0*t%p23)+0.406d0*t%t9
      da1 = a1*(cc13*36.031d0*t%m13-2.0d0*(t%t9/0.881d0)**2)
      nr%dlnrate_dlnt(r_an14) = dscreen_dlnt(r_an14)+(dr1*a1+r1*da1+ &
              a2*(2.798d0*t%m1-1.5d0)+a3*(5.054d0*t%m1-1.5d0)+ &
              a4*(12.310d0*t%m1-cc23))/s_sum
! triple alpha
      nr%rate(r_3alpha) = 1.565315d21*density**2*t%m1*t%m2*2.79d-8* &
                 exp(-4.4027d0*t%m1+screen(r_3alpha))
      nr%dlnrate_dlnrho(r_3alpha) = 2.0d0+dscreen_dlnrho(r_3alpha)
      nr%dlnrate_dlnt(r_3alpha) = -3.0d0+dscreen_dlnt(r_3alpha)+4.4027d0*t%m1
end subroutine compute_helium_burning_rates

!----------------------------------------------------------------------
! hep, He3(p,e+ nu)He4, for the neutrino flux only: the screened, scaled
! rate factor, which the consumer multiplies by X(H) X(He3).
subroutine compute_hep_rate(density, t, plasma, nr)
      use math_lib
      double precision, intent(in) :: density
      type(t9_powers), intent(in) :: t
      type(screening_plasma), intent(in) :: plasma
      type(nuclear_rates), intent(inout) :: nr

      nr%hep = (1.71724d+11)*density*t%m23*exp(-6.1399d0*t%m13)
      nr%hep = (1.0d0 + 0.067862d0*t%p13)*net_settings%cross_section_scale(xs_hep)*nr%hep
      nr%hep = nr%hep*exp(weak_intermediate_screening(plasma, 2.0d0, 3.08687d0))
end subroutine compute_hep_rate

!----------------------------------------------------------------------
! screening_mass_fractions: the isotope mass fractions (iso_* slots) that
! compute_reaction_rates screens with, from a zone's composition in the
! star%xa layout. The isotopes YREC does not track (n, H3, Ne20, Mg24)
! are zero; deuterium is included (2026-10: the composition path used to
! leave it out).
function screening_mass_fractions(xa) result(mass_frac)
      use star_info_lib, only: i_h1, i_h2, i_he4, i_he3, i_c12, i_c13, &
           i_n14, i_o16, i_o18
      double precision, intent(in) :: xa(:)
      double precision :: mass_frac(num_isotopes)

      mass_frac = 0.0d0
      mass_frac(iso_h1) = xa(i_h1)
      mass_frac(iso_h2) = xa(i_h2)
      mass_frac(iso_he3) = xa(i_he3)
      mass_frac(iso_he4) = xa(i_he4)
      mass_frac(iso_c12) = xa(i_c12)
      mass_frac(iso_c13) = xa(i_c13)
      mass_frac(iso_n14) = xa(i_n14)
      mass_frac(iso_o16) = xa(i_o16)
      mass_frac(iso_o18) = xa(i_o18)
end function screening_mass_fractions


!----------------------------------------------------------------------
! compute_network_rates: the composition rate vector of one zone
! (rate_vec, num_rate_rows): the reaction rates at their r_* slots,
! zeroed at or below 1e-5 and turned into per second per amu (times
! amu_cgs), then the N15(p,alpha) and Be7 electron-capture fractions at
! rv_frac_n15_pa and rv_frac_be7_ecap. Slots 9 and 13 are zero. xa is
! the zone's composition in the star%xa layout.
subroutine compute_network_rates(log_density, log_temperature, xa, rate_vec)
      use physical_constants, only: amu_cgs
      double precision, intent(in) :: log_density, log_temperature, xa(:)
      double precision, intent(out) :: rate_vec(num_rate_rows)

      double precision :: mass_frac(num_isotopes)
      type(nuclear_rates) :: nr
      integer :: i

      mass_frac = screening_mass_fractions(xa)
      call compute_reaction_rates(log_density, log_temperature, mass_frac, nr)
      do i=1,num_reactions
         if(nr%rate(i).le.1.d-5) nr%rate(i) = 0.0d0
      end do
! per s per amu (amu_cgs)
      rate_vec(r_pp)=nr%rate(r_pp)*amu_cgs
      rate_vec(r_he3he3)=nr%rate(r_he3he3)*amu_cgs
      rate_vec(r_he3he4)=nr%rate(r_he3he4)*amu_cgs
      rate_vec(r_pc12)=nr%rate(r_pc12)*amu_cgs
      rate_vec(r_pc13)=nr%rate(r_pc13)*amu_cgs
      rate_vec(r_pn14)=nr%rate(r_pn14)*amu_cgs
      rate_vec(r_po16)=nr%rate(r_po16)*amu_cgs
      rate_vec(r_ac13)=nr%rate(r_ac13)*amu_cgs
      rate_vec(r_ao16_unused)=nr%rate(r_ao16_unused)*amu_cgs
      rate_vec(r_ac12)=nr%rate(r_ac12)*amu_cgs
      rate_vec(r_an14)=nr%rate(r_an14)*amu_cgs
      rate_vec(r_3alpha)=nr%rate(r_3alpha)*amu_cgs
      rate_vec(r_c12c12_unused)=nr%rate(r_c12c12_unused)*amu_cgs
      rate_vec(rv_frac_n15_pa)=nr%frac_n15_pa
      rate_vec(rv_frac_be7_ecap)=nr%frac_be7_ecap
end subroutine compute_network_rates



end module reaction_rates

!======================================================================
! Sources and history of the rate formulas and coefficients, kept
! verbatim from the bodies of engeb.f and rates.f (JNB, MHP, DBG). The
! Q1..Q8 they describe are fit_q1..fit_q8; ZPRD, Z86, Z53, Z43, Z23 are
! screen_zprod, screen_z86, ...; SSTANDARD is cross_section_scale; XXL,
! XXL6, XXL8 are lambda0, lambda0_23, lambda0_086; ZCURL is z_tilde, ZBAR
! mean_charge, AMU/EMU inv_mu_i/inv_mu_e, XTR/ZET the charge moments;
! PFMC2, EFMKT, FPRF, DEGD are fermi_momentum_sq, fermi_energy_over_kt,
! degeneracy_factor, log_degeneracy; UWK/UINT/USTR the weak, intermediate
! and strong exponents; R1 s_ratio; TEMP3 tau_be7p, QRBE7 s_ratio_be7p,
! CAMUBE7 per_gram_be7; O16GAMMA/C12ALPHA pg_strength/pa_strength.
!======================================================================
! Provenance of the formulas and coefficients (moved verbatim from the
! body of rates, 2026-10, when this routine became the one copy):
! ANUC ARE ATOMIC MASS UNITS.
! On 10/13/97, JNB converted ANUC(I) from neutral atomic masses to
! bare nuclear masses.
! The scale is the mass of C12 divided by 12 or 931.49432 MeV,
! which is 1.6605402 times 10^{-24} gm. The values are obtained
! by dividing the mass excess (expressed in MeV) by 931.49432 MeV and
! adding to this the atomic mass number, A.  The value for Be7, which
! is used implicitly in this subroutine was, 7.016930, until 13/13/97.
! We now use the bare nuclear mass of 7.014735 .
! THE ISOTOPES ARE NEUTRON, H1, D, H3, HE3, HE4, C12, C13, N14, O16, O18,
! NE20, MG24, RESPECTIVELY. ALL OF THESE NUMBERS WERE CHECKED.
! THE QUANTITIES Q1(J), Q2(J), ...,Q5(J) ARE THE TERMS IN EQUATION 3.14 OF
! NEUTRINO ASTROPHYSICS AND IN EQUATION 53 OF FOWLER, CAUGHLAN, AND
! ZIMMERMAN (1967) EQ. 53, IN BOTH CASES MULTIPLIED BY T SUB 9 ^(-2/3).
! THE REACTIONS CORRESPONDING TO EACH J ARE LISTED ABOVE, UNDER:
! IDENTIFYING THE RECTIONS.
! FOR THIS SET OF PARAMETERS, AND ONLY FOR THIS SET OF PARAMETERS,
! J = 8 CORRESPONDS TO THE BE7 + P REACTION.
! GENERAL EXPRESSION FOR Q'S:
! T_9^(-2/3)[S_EFF/S(0)] =
! [T_9^(-2/3) + Q1(I)T_9^(-1/3) + Q2(I) + Q3(I)T_9^(1/3) +
! Q4(I)T_9^(2/3) + Q(5)T_9 ]
! BY COMPARISON WITH EQUATION 3.14 OF NEUTRINO ASTROPHYSICS, WE SEE THAT
! EACH OF THE Q'S IS INDEPENDENT OF TEMPERATURE (T), AS CAN BE SEEN FROM
! EQUATIONS 3.10 AND 3.11 .
! ALL OF THE VALUES OF THE Q1, ...., Q5 HAVE BEEN RECALCULATED, USING
! WHERE NEEDED NUCLEAR CROSS SECTIONS GIVEN IN TABLES 3.2 AND 3.4 OF NEUTRINO
! ASTROPHYSICS.
! Q6 IS THE COEFFICIENT OF THE TEMPERATURE TERM IN THE DEFINITION OF
! TAU, EQUATION 3.10 OF NEUTRINO ASTROPHYSICS.
! TAU = Q6*(T SUB 9 TO THE (-1/3) POWER ).
! SLIGHT CHANGES HAVE BEEN MADE IN THE PREVIOUS VALUES OF Q6 TO MAKE
! THE DATA MORE ACCURATE.
! NOTE Q6 IS NEGATIVE.
! Q7 IS THE CONSTANT IN FRONT OF THE REACTION RATE. THE LOGARITHMS ARE
! RELATED TO S FACTORS BY A LOGARITHM AND VARIOUS NUMERICAL FACTORS.
! THE GENERAL RELATION IS:
! Q7 = 70.62860 + ((LN(Z0*Z1/A))/3) -LN(A0*A1) - LN(S SUB 0)
! -LN(1 + DELTA_01)
! HERE S SUB 0 IS THE CROSS SECTION FACTOR IN UNITS OF KEV-BARNS.
! THE NUMERICAL VALUES USED HERE ARE TAKEN FROM TABLES 3.4 AND 3.2
! OF NUCLEAR ASTROPHYSICS.
! THE QUANTITY DELTA_01 IS NON-ZERO (EQUAL TO UNITY) ONLY WHEN THE
! TWO REACTING NUCLEI, 0 AND 1, ARE IDENTICAL.
! Q8 REFLECTS A TERM IN THE EXPONENTIAL THAT OCCURS IN THE RATES, THE
! TERM BEING PROPORTIONAL TO E^( CONSTANT*T_9^2). SEE HARRIS, ET AL.
! 1983, ANNUAL REV. ASTRON. ASTROPHYS.21, 165 (1983) FOR THE MEANING
! OF THIS OBSCURE TERM.
! THE VALUES OF Q1(I),...Q7(I) GIVEN BELOW WERE OBTAINED USING THE DATA
! IN NEUTRINO ASTROPHYSICS WITH THE AID OF AUXILARY COMPUTER CODES
! THAT GENERATED ACCURATE EVALUATIONS. THE NUMBERS GIVEN HERE ARE
! IN MANY CASES COMPLETELY DIFFERENT FROM THE VALUES IN THE ORIGINAL
! YALE CODE.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! 10/13/97. Changed Q(I) so that are now calculated for the bare nuclear
! masses. The calculates were made cues.f . Last date on the calculations
! in my notes is 10/10/97.  Note Q6 is -tau(T_9 = 1); in cues.f, I calculate
! tau(T_6 = 1). The connection is:  Q6 = 0.1*tau(cues.f).
! BP00 values
! DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
! 1 Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, -0.2095/,
! 2 Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/,
! 3 Q4/0.,0.,0.,4.79,7.456,0.0,.69703,.16762/,
! 4 Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,.12114/,
! 5 Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,
! 5    -16.6925/,
! 6 Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
! 7 Q8/0.,0.,0.,0.0,0.0,0.0,0.0/,
! 8NRXNS/13/
! Revised values 7/21/03
! DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
! 1Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, 0.0/,
! 2Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,0.0/,
! 3Q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.0/,
! 4Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.0/,
! 5Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/,
! 6Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
! 7Q8/0.,0.,0.,0.0,0.0,0.0,0./,
! 8NRXNS/13/
! MHP 9/14 RESTORED BE7+P DERIVATIVES THAT WERE ZEROED OUT
! For different values of SSTANDARD, check the comments in ENGEB
! THE VALUES OF SSTANDARD(I) ARE TO BE CHANGED FROM UNITY IF THE CROSS
! SECTION FACTORS ARE NOT THE ONES GIVEN IN NEUTRINO ASTROPHYSICS,
! TABLE 3.2 AND TABLE 3.4 .
! THE CURRENT VALUES ARE CHANGED TO BE THE PREFERRED VALUES LISTED IN
! THE LAST COLUMN OF TABLE 1 OF BAHCALL AND PINSONNEAULT (1991),
! 6/14/91.
! THE VALUE OF SSTANDARD(17) FOR HEP IS TAKEN FROM WOLFS ET AL. , PHYS.
! REV. LETTERS, 63, 2721 (1989).  IT CORRESPONDS TO AN S SUB 0 =
! 15.3E-20 KEV-BARNS, A FACTOR OF 1.913 LARGER THAN FOR THE OLDER
! MEASUREMENTS USED IN NEUTRINO ASTROPHYSICS.
! ZPRD IS USED IN SCREENING CALCULATIONS. IT IS THE PRODUCT OF THE
! CHARGES OF THE INTERACTING IONS. ZPRD WAS CHECKED. Z86 IS USED
! IN CALCULATING INTERMEDIATE SCREENING AND IS DEFINED BY GRABOSKE ET
! AL, AP. J. 181, PAGE 465 (1973), IN TABLE 4.  Z86 WAS CHECKED AND
! SOME NUMERICAL VALUES WERE MADE SLIGHTLY MORE ACCURATE.  Z53, Z43,
! AND Z23 ARE ALSO DEFINED IN TABLE 4 (SEE ABOVE).  SINCE THEY ARE
! ONLY USED IN STRONG SCREENING, THE VALUES OF Z53, Z43, AND Z23
! WERE NOT CHECKED.
! DEFINE NEXT THE FRACTIONAL ABUNDANCES BY MASS OF THE IMPORTANT
! ISOTOPES.
! X, Y, Z, XHE3,..., XBE9 ARE THE MASS FRACTIONS OF THE ISOTOPES.
! THE ABUNDANCES OF NEUTRONS, H2, H3, NE20,AND MG24, WHICH ARE,
! RESPECTIVELY, XFRAC(I) FOR I = 1,3,4,12,13, ARE NO LONGER USED.
! BEGIN CALCULATION OF SCREENING CORRECTION.
! THE BASIC REFERENCES ARE SALPETER, AUSTRALIAN JOURNAL OF PHYSICS,
! VOL. 7, 373 (1954). THE FORMULA FOR WEAK SCREENING THAT IS BEING
! PROGRAMMED IS EQUATION (25) OF THIS PAPER. THE OTHER IMPORTANT
! REFERENCES ARE: DEWITT, GRABOSKE, AND COOPER, AP. J. 181, 439 (1973)
! AND GRABOSKE ET AL., AP. J. 181, 457 (1973). THE VALUES OF EMU AND
! ZET ARE ESSENTIAL FOR COMPUTING WEAK SCREENING; THE VALUE OF AMU IS
! USED IN AN NON-ESSENTIAL WAY IN THIS COMPUTATION. XTR IS USED IN
! COMPUTING INTERMEDIATE SCREENING.
! AMU IS ONE OVER THE MEAN MOLECULAR WEIGHT OF THE IONS, MU SUB I .
! EMU IS ONE OVER THE ELECTRON MEAN MOLECULAR WEIGHT, MU SUB E.
! EMU IS USED HERE AS THE NAME FOR THE SECOND PART OF THE ZETA FUNCTION
! IN THE SALPETER EXPRESSION FOR WEAK SCREENING.
! XTR IS USED LATER IN THE INTERMEDIATE SCREENING CALCULATION.  THE AVERAGE
! OF THE QUANTITY Z**(3B -1) IS EQUAL TO XTR/AMU.
! ZET IS THE FIRST PART OF THE SALPETER SCREENING ZETA VARIABLE.
! MU = SUM OVER I OF [X(I)/A(I)].
! MU SUB E = SUM OVER I [ Z(I)*X(I)/A(I)].
! DL AND DT ARE THE THE LOG10 OF THE DENSITY AND TEMPERATURE.
! THE UNIT OF TEMPERATURE IS 10^9 K AND THE UNIT OF DENSITY IS
! GM PER CM^3 .
! DD = LOG RHO TO THE BASE 10.
! CLN = LN10.  CLN IS CONVERSION BETWEEN LOG10 AND LN.
! CONVERT DENSITY TO UNLOGGED FORM.
! RWE = RHO/(MU SUB E), I. E., THE NUMBER OF ELECTRONS DIVIDED BY
! AVOGADRO'S NUMBER.
! THE EXPRESSION FOR RWE WAS INCORRECT IN THE ORIGINAL YALE SUBROUTINE.
! THE ORIGINAL VERSION HAD ( EXP(CLN*DL) ) DIVIDED BY EMU INSTEAD OF
! MULTIPLIED BY EMU.  RWE IS USED LATER IN COMPUTING THE SCREENING
! CORRECTION.
! SET RATES EQUAL TO ZERO FOR THE LOG_10(T) < TCUT(1)
! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
! POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
! OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
! CGS UNITS.
! F PRIME/F
! THE NEXT PIECE OF CODE COMPUTES FIRST THE FERMI ENERGY DIVIDED BY KT, WHERE
! PFMC2 IS THE FERMI MOMENTUM DIVIDED BY MASS OF ELECTRON TIMES C, ALL
! SQUARED AND EFMKT IS THE FERMI ENERGY DIVIDED BY KT.  THE QUANTITY
! FPRF IS THE RATIO OF F PRIME TO F IN SALPETER'S SCREENING CORRECTION.
! THE VALUE OF FPRF IS DETERMINED IN THE INTERMEDIATE CASE BY AN
! INTERPOLATION FORMULA DEPENDING UPON THE DEGREE OF DEGENERACY, AS
! MEASURED BY DEGD = LOG10(E_F/KT). THE NUMERICAL VALUES FOR THE FIT
! WERE TAKEN FROM SALPETER'S ORIGINAL PAPER, FIGURE 1. THE ONLY CHANGES
! IN THIS PART OF THE SUBROUTINE WERE THE CORRECTION OF THE ERROR IN
! THE DEFINITION OF RWE (SEE ABOVE) AND REFINEMENTS OF THE COEFFICIENTS
! IN THE EXPRESSIONS FOR PFMC2 AND EFMKT.
! END OF CALCULATION OF FERMI ENERGY DIVIDED BY KT AND THE INTERPOLATION
! FORMULA FOR (F PRIME/F), WHICH APPEARS IN SALPETER'S SCREENING
! FORMULA.
! NOW WE GET TO THE COMPUTATION OF WEAK SCREENING (SEE ALSO PAGE 61 OF
! NEUTRINO ASTROPHYSICS, WHICH GIVES ONLY A SIMPLIFIED FORMULA) AND
! ALSO THE MORE COMPLICATED INTERMEDIATE AND STRONG SCREENING CASES
! (SEE REFERENCES TO AP. J. 181 ABOVE FOR INTERMEDIATE
! AND STRONG SCREENING), ESPECIALLY PAGE 465.
! XXL IS USED IN ALL THE SCREENING FORMULAE.  THIS QUANTITY IS THE
! IS THE FUNCTION CALLED LAMBDA SUB ZERO BY GRABOSKE ET AL. THE
! QUANTITY XXL**0.86 = XXL8 IS USED IN CALCULATING INTERMEDIATE
! SCREENING.
! ZCURL IS THE QUANTITY DEFINED BY EQUATION (4) OF DEWITT ET AL; IT
! IS THEIR Z WITH A CURLY SYMBOL ON ITS TOP. ZCURL IS USED IN WEAK
! AND IN INTERMEDIATE SCREENING AND WAS FIRST DEFINED BY SALPETER.
! ZCURL IS THE SAME AS SALPETER'S ZETA EXCEPT FOR THE FACTOR OF 1/AMU.
! (Z SUB 1 TIMES Z SUB 2)*XXL*ZCURL GIVES THE WEAK SCREENING FACTOR,
! THE SAME AS SALPETER OR AS EQUATION (19) OF GRABOSKE ET AL.
! Z BAR IS THE SAME AS Z BAR OF DE WITT ET AL AND IS THE AVERAGE CHARGE
! OF THE IONS.  IT IS EQUAL TO EMU/AMU.  THE QUANTITY (Z BAR)**0.28 =
! Z28 OCCURS IN THE COMPUTATION OF INTERMEDIATE SCREENING.
! XXL6 IS USED FOR COMPUTING STRONG STRONG SCREENING.
! THE NOTATION USED HERE IS EXPLAINED IN LARGE PART BY
! EQUATION (4) OF DEWITT ET AL., AP. J. 181, PAGE 439.
! THE FINAL EXPRESSION FOR WEAK SCREENING IS EXACTLY EQUAL TO SALPETER'S
! FORMULA, WHICH INCLUDES A DEGENERACY CORRECTION. THE MORE GENERAL
! EXPRESSIONS ARE GIVEN IN TABLE 4 AND EQUATION (19) OF GRABOSKE ET AL.
! COMPUTE SCREENING FOR EACH OF THE REACTIONS.
! net_settings%weak_screening_threshold (WEAKSCREENING) IS A NUMERICAL
! PARAMETER READ FROM THE INPUT. TO OBTAIN THE GRABOSKE ET AL. AND SALPETER STANDARD RESULTS,
! USE: WEAKSCREENING = 0.03.  FOR THE STANDARD SOLAR MODEL, THIS IS THE
! VALUE THAT SHOULD BE ADOPTED. TO INVESTIGATE THE EFFECT OF ALWAYS USING
! WEAK SCREENING, USE A LARGE VALUE FOR WEAKSCREENING, E. G., 30.  AS
! LONG AS WEAKSCREENING IS ASSUMED TO BE BIGGER THAN ONE, THE PROGRAM
! WILL ALWAYS CALCULATE FOR THE SUN WITH THE WEAK SCREENING
! APPROXIMATION.
! UTOT IS THE FINAL SCREENING CORRECTION WHICH APPEARS IN THE
! RATE EXPRESSION AS: EXP(UTOT) .
! END OF SCREENING CALCULATION. WEAK AND INTERMEDIATE SCREENING FORMS
! ARE GIVEN CORRECTLY.  STRONG SCREENING WAS NOT CHECKED BECAUSE IT IS
! NOT RELEVANT FOR THE SUN.
! Reactions first_zeroed_rxn..num_reactions are zeroed after this block:
! only the alpha captures when hydrogen burns, all 13 when it is gone.
! CALCULATE REACTION RATES FOR THE THREE PRINCIPAL REACTIONS OF
! THE PP CHAIN: PP, HE3+HE3, HE3 +HE4, AND THE FOUR PROTON
! BURNING REACTIONS ON C12: C13, N14, AND O16.
! R1 IS (T SUB 9)^(-3/2) TIMES (S SUB EFF)/(S SUB 0). THE CORRECT
! EXPRESSION FOR (S SUB EFF)/(S SUB 0) IS GIVEN IN EQUATION 3.14 OF
! NEUTRINO ASTROPHYSICS.  THE NUMERICAL FORM THAT IS USED IS EQUATION
! 52 OF FOWLER, CAUHLAN, AND ZIMMERMAN, VOL. 5, 1967.
! RATE(I) IS THE RATE OF THE DIFFERENT REACTIONS PER SECOND PER GRAM,
! EXCEPT THAT THE MASS FRACTIONS ARE OMITTED.
! THE PREVIOUS YALE VERSION HAD THE RATE FOR THE O16 + P REACTION
! MULTIPLIED BY T9**(-1/7). THIS FACTOR IS INCORRECT AND HAS BEEN
! REMOVED; IT APPEARED BEFORE AS AN IF STATEMENT REFERRING ONLY TO
! RATE(7).
! R1=T9M23+Q1(I)*T9M13+Q2(I)+Q3(I)*T9P13+Q4(I)*T9P23+Q5(I)*T9
! MHP 8/14 RATES CORRECTED TO PERMIT USER MODIFICATION OF REACTION
! RATE DERIVATIVES
! END OF CALCULATION OF REACTION RATES FOR FIRST 7 REACTIONS.
! BE7 BURNING
! CALCULATE THE BURNING OF BE7 BY PROTONS (WHICH PRODUCES THE MOST
! EXPERIMENTALLY ACCESSIBLE SOLAR NEUTRINOS) AND THE BURNING OF BE7
! BY ELECTRON CAPTURE (THE DOMINANT PROCESS). THE PREVIOUSLY USED
! EXPRESSION CONTAINED A TERM THAT WAS PHYSICALLY INCORRECT
! AND ALSO USED CROSS SECTION FACTORS THAT WERE OUT OF DATE.
! THE ELECTRON CAPTURE RATE IN SEC^(-1) IS GIVEN BY EQUATION (3.18) OF
! NEUTRINO ASTROPHYSICS.  CAN OMIT THE FACTOR OF RHO (IN CGS UNITS)
! WHICH ALSO APPEARS IN THE BE7PROTON RATE. THE BE7 PROTON CAPTURE
! RATE IS GIVEN BY TABLE 3.2 AND EQUATION (3.12).  NOTE THAT 1 OVER MU
! SUB E IS EQUAL TO EMU.
! HERE WE USE THE NOTATION BE7 + PROTON = BE7PROTON AND
! BE7 + E = BE7ELECTRON.
! F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY ELECTRON CAPTURE.
! 1-F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY PROTON CAPTURE.
! F3 IS THE FRACTION OF THE N14 THAT IS BURNED BY P, ALPHA REACTION.
! 1-F3 IS THE FRACTION OF THE N15 THAT IS BURNED BY P, GAMMA REACTION.
! SEE TABLE 21 OF BAHCALL AND ULRICH (1988), REV. MOD. PHYS. 60.
! 10/13/97. I changed Temp3 (i.e., tau) by 5/10^5 as a result of
! using bare nuclear masses. Previously coefficient was -10.26202.
! I did not change Be7electron for the slightly different S0 for
! electron capture, since that is done in SStandard. Previous coefficient
! in Be7electron expression was (3.126571E+5). 10/14/97.
! INCLUDE FOR BE7PROTON THE T9M23 FACTOR AND ALL CORRECTIONS PROPORTIONAL TO
! Q1,...,Q5 FROM EQUATION 3.14 OF NEUTRINO ASTROPHYSICS. THESE
! CORRECTIONS ARE DEFINED EARLIER IN THIS SUBROUTINE.
! QRBE7 = T9M23 + Q1(8)*T9M13 + Q2(8)+ Q3(8)*T9P13
! $          + Q4(8)*T9P23 + Q5(8)*T9
! MHP 9/14 ADDED THE ABILITY TO ALTER DERIVATIVES INDEPENDENTLY
! CALCULATE THE SCREENING CORRECTION FOR BE7 + P REACTION.  USE WEAK AND
! INTERMEDIATE SCREENING FORMULAE.
! END OF CALCULATION OF SCREENING CORRECTION FOR BE7 + P REACTION.
! MULTIPLY RATES BY FACTOR OF RHO/[(ATOMIC MASS UNIT)*A(BE7)] TO GET
! IN UNITS OF GM^{-1}.  CALL FACTOR CAMUBE7. Corrected for Be7 bare
! mass number. Previously used 8.582295E+22 in the expression below.
! These corrections were made on 10/14/97.
! END OF MULTIPLICATION INSERTED NOVEMBER 6, 1990.
! END OF CALCULATION OF CRUCIAL BE7 ELECTRON CAPTURE AND PROTON CAPTURE
! RATES AND THEIR RATIO.
! N15 + P BRANCHING.
! THE FOLLOWING STATEMENTS COMPUTE THE BRANCHING OF N15(P,ALPHA)C12
! AND N15(P,GAMMA)O16. THESE STATEMENTS REPLACE OUTDATED STATEMENTS
! IN THE YALE CODE. THE CNO CROSS-SECTION FACTORS ARE FROM TABLE 3.4
! OF NEUTRINO ASTROPHYSICS.  THE RATIO OF THE REACTIONS DEPENDS ONLY
! UPON THE EFFECTIVE ZERO ENERGY S-FACTOR, WHICH IS S0(ZERO ENERGY)
! TIMES THE COMBINATION OF TEMPERATURE AND S-FACTOR DERIVATIVES/S0
! THAT WAS USED PREVIOUSLY AS R1 IN THE RATE CALCULATIONS. THE
! NUMERICAL COEFFICIENTS THAT APPEAR IN THE RATE WERE REPRESENTED
! BY THE Q1(J),...Q5(J) FOR THE OTHER REACTIONS.
! THE QVALUES FOR THE N15 REACTIONS HAVE BEEN COMPUTED SEPARATELY.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! To agree with Solar Fusion Workshop paper, the value of S0 in keV-b
! has been changed from 78000. to 675000. on 9/25/97. On 10/14/97,
! JNB changed cues.f so as to compute the Q-coefficients for the N15 + p
! reactions.  Also, checked that the coefficients are the same as the
! ones given earlier in energy.f when we use the older CNO data.
! MULTIPLY BY THE VALUE OF S0 IN KEV-B.
! MHP 9/14 ADDED THE OPTION TO MODIFY THE RELATIVE CROSS SECTIONS
! FOR N15+P -> C12+ALPHA AND O16+GAMMA
! END OF NEW ROUTINE FOR THE BRANCHING OF N15 + P .
! ***MHP 3/91 ALPHA CAPTURE REACTIONS UPDATED TO CAUGHLAN AND FOWLER(1988)
! RATES.  THE RATES ARE EXPRESSED IN THE SAME TERMS USED BY CZ, WITH
! THE CONVERSION FACTOR IN THE FRONT OBTAINED FROM VANDENBERG'S
! NOTES ON THE REACTION RATES.
! RATE(8)  HE4+C13
! RATE(10) HE4+C12=>O16
! RATE(11) HE4+N14=>O18
! RATE(12) TRIPLE ALPHA
! C13(ALPHA,N) O16
! C12(ALPHA,GAMMA)O16
! N14(ALPHA,GAMMA)F18 + F18=>O18+EPLUS+NU
! 2026 (bugsweep Batch 3): the linear coefficient of the CF88
! N14(a,g)F18 polynomial is 0.117 (1 + 0.012 T9^1/3 + 1.45 T9^2/3
! + 0.117 T9 + 1.97 T9^4/3 + 0.406 T9^5/3); both copies of the rate
! carried 0.177 -- a digit transposition inherited from the F77
! source (~0.5% of the rate at T9 = 0.15).
! TRIPLE ALPHA
! END OF XEROING OUT OF REACTIONS 9 AND 13.
! RATES PER SECOND PER ATOMIC MASS UNIT (2026-10; per 10^9 sidereal
! years per amu, the C21 factor, before): HRK(IU)
! HR1, ..., HR13 ARE THE RATES OF THE INDIVIDUAL REACTIONS.
! THE INTERPRETATION OF WHICH REACTION GOES WITH WHICH SYMBOL CAN BE
! MADE EASILY BY LOOKING AT THE DEFINITIONS OF THE EG(I)'S.
! THE ABUNDANCES ARE UPDATED BY THE NETWORK BURN USING THESE.
! rows r_ao16_unused/r_c12c12_unused: zero-filled placeholders (rate(9), rate(13)
! are always 0 -- see the XEROING block above).
! END OF COMPUTATION OF HRK(IU).
!
! JULY 3, 1991 (MHP)
! THIS SUBROUTINE COMPUTES THE NUCLEAR BURNING RATES FOR USE IN KEMCOM.
! IT IS A STRIPPED-DOWN VERSION OF ENGEB; SEE ENGEB FOR DETAILED NOTES
! ON THE REACTIONS.
! ALL PHYSICAL QUANTITIES ARE CGS. THE REACTION RATES ARE IN
!  GM^{-1}S^{-1}.
!  MHP 10/97
!
! On 10/13/97, JNB converted the nuclear masses from neutral nuclear
! masses to bare nuclear masses by subtracting Z(I)*(m_e)*c^2 from
! the neutral nuclear masses. This caused changes in a number of
! places: in ANUC(I), in Q1(I)-Q7(I), in the calculation of the Be7electron
! Be7proton rates, and in the calculation of the N15p branching ratio.
!
! JNB made some purely cosmetic changes on 1/20/96. Revised all of the
! input Q's on 9/23-25/97 to agree with submitted version of Solar
! Fusion Workshop paper. The SStandard are fixed to agree with
! the Workshop paper. JNB recalculated all of the EG(I) to determine
! the best values for the energy generation for all the reactions,
! taking account of my improved calculations of neutrino energy loss.
! The calculations are documented in Vol. 19, 132-141, 1997 of my notes.
! The neutral atom mass differences are taken from Table of Isotopes,
! 8th Ed, 1996 and the neutrino energy losses from Bahcall, Gallium
! solar neutrino experiments, Phys Rev C, in press, 1997.
!
! ALL NUMBERS IN THIS SUBROUTINE HAVE BEEN CHECKED AND REVISED, WHERE
!  NECESSARY, BY JOHN BAHCALL SO THAT THEY AGREE WITH THE MODERN NUMBERS
!  IN NEUTRINO ASTROPHYSICS (1989).
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
! THE ENERGY DERIVATIVES ENTER IN A FORM IN WHICH THEY ARE DIVIDED BY
!  THE ABSOLUTE VALUES OF THE CROSS SECTIONS AT ZERO ENERGY.  THUS IF
!  THE SHAPE OF THE CROSS SECTION EXTRAPOLATION IS UNCHANGED AND ONLY
!  THE INTERCEPT OF S(E) AT ZERO ENERGY IS CHANGED, THEN NO CORRECTION
!  NEED BE MADE FOR THE DERIVATIVES.  THEY ARE AUTOMATICALLY SCALED
!  CORRECTLY.  THE EXACT WAY THAT THE DERIVATIVES ENETER THE
!  CALCUALTIONS IS DESCRIBED IN THE SECTION LABELED ``DEFINING THE
!  Q(I)' THAT IS PRESENTED BELOW.
! ************************************
! IDENTIFYING THE REACTIONS.
! ************************************
!  THE VALUE OF J DENOTES WHICH OF THE REACTIONS THE COEFFICIENTS
!  REFER TO:
!  J = 1, PP; J = 2, HE3+HE3; J = 3, HE3+ HE4; J =4, P + C12;  J = 5, P+C13;
!  J = 6. P + N14; J = 7, P + O16 ; J = 8, HE4+C13; J = 10, HE4+C12;
!  J = 11, HE4 + N14; J = 12, TRIPLE ALPHA.
!  REACTION 14 IS PEP; REACTION 15 IS BE7 ELECTRON CAPTURE; REACTION 16 IS
!   BE7 PROTON CAPTURE; REACTION 17 IS THE HEP REACTION.
!  do not change sstandard(14) unless you want to change the ratio of
!  pep to pp.
!   REACTIONS 14-17 WERE NOT EXPLICITLY INCLUDED IN THE YALE
!    PREVIOUS VERSION OF THE CODE, BUT THEY ARE MOST OF THE STORY FOR
!    THE SOLAR NEUTRINO PROBLEM.
!   THE BRANCHING OF THE N15 + P REACTIONS IS TREATED IN A SERIES OF
!    SEPARATE STATEMENTS FOLLOWING THE CALCULATION OF THE BE7 + P
!    REACTION. SEE THE DEFINITIONS OF F3 AND F4.  IF THE CROSS-SECTION
!    FACTORS OF THE N15 + P REACTIONS ARE REVISED, THEN THE NUMERICAL
!    COEFFICIENTS MUST BE CHANGED IN THE DEFINITION OF C12ALPHA AND
!    O16GAMMA.
! FOR Q1(I), ...,Q(5(I), I = 8 CORRESPONDS TO THE BE7 +  P REACTION.
!  THIS ASSIGNMENT FOR I = 8 IS ONLY VALID FOR THE LISTED Q'S AND NOT
!  FOR OTHER ARRAYS IN THE PROGRAM.
