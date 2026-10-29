!----------------------------------------------------------------------
! reaction_rates (formerly the body of net_lib)
!----------------------------------------------------------------------
! The rate physics of the network: pure functions of (log rho, log T,
! composition) and the run's settings. network's compute_reaction_rates
! evaluates each reaction of its table with these:
!     compute_zone_conditions  the density, T9 powers and screening
!                              plasma of a zone
!     compute_screened_rate    one reaction's screened rate and its T/rho
!                              derivatives, from the rate law, fit,
!                              slot and screening constants of its
!                              table row
!     the rate laws            the S-factor law of the seven proton and
!                              He3 captures and the Caughlan & Fowler
!                              alpha captures and triple alpha
!     compute_be7_branching, n15_alpha_fraction, hep_rate_factor
!                              the branchings and the hep factor
!     compute_deuterium_proton_capture, cap_deuterium_rate
!                              the d(p,g) fit and its convective-overturn
!                              cap, for the deuterium depletion
!                              (light_elements) and its energy
!                              (nuclear_energy) alike
! It also holds the screening isotopes (iso_*) and their masses and
! charges, the rate-law names a table row can give, and the screening
! constants' type. The reaction slots and the settings it reads are
! net_configuration's.
! Sources and history: the end of this file.
! ---------------------------------------------------------------------
module reaction_rates
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use net_configuration, only: num_reactions, xs_be7e, xs_be7p, xs_hep, &
           net_settings
      implicit none
      private

! The isotopes of the screening composition and of the mass and charge
! tables.
      integer, parameter, public :: &
           iso_n = 1, iso_h1 = 2, iso_h2 = 3, iso_h3 = 4, iso_he3 = 5, &
           iso_he4 = 6, iso_c12 = 7, iso_c13 = 8, iso_n14 = 9, &
           iso_o16 = 10, iso_o18 = 11, iso_ne20 = 12, iso_mg24 = 13, &
           num_isotopes = 13
! S0 (keV-b) of the two N15 + p branches (n15_alpha_fraction); integers,
! as the original literals were (both exact when promoted).
      integer, parameter :: s0_n15pg_kevb = 64, s0_n15pa_kevb = 67500

! d(p,g)3He: Avogadro's number over the deuteron mass (amu) turns CF88's
! N_A<sigma v> into the per-gram rate of both the deuterium depletion and
! its energy release, and sets the convective-overturn cap. (The two
! former copies differed by 0.68%: one divided by A = 2.)
      real(dp), parameter :: avogadro_dpg = 6.023d23
      real(dp), parameter :: deuteron_mass_amu = 2.013553d0
! The overturn cap applies in the surface convection zone when the
! turnover time exceeds min_cap_turnover_time (s) and X(H) exceeds
! min_cap_hydrogen (cap_deuterium_rate).
      real(dp), parameter :: min_cap_turnover_time = 1.0d0, &
           min_cap_hydrogen = 1.0d-6
! The Be7(p,g) S_eff/S0 fit, q1..q5 (FCZ67 eq. 53 / Neutrino Astrophysics
! eq. 3.14; history: the end of this file). The fits of the network's
! S-factor reactions are in network's reaction table, one row each.
      real(dp), parameter, private :: &
           be7_pcapture_fit(5) = [.040572d0, -0.2095d0, -0.0595d0, 0.16762d0, 0.12114d0]
! The rate laws a reaction-table row can name (compute_screened_rate),
! 1 to num_rate_laws. shares_rate marks a second outcome of the reaction
! in the row above.
      integer, parameter, public :: shares_rate = 0, s_factor_fit = 1, &
           cf88_c13_alpha_n = 2, cf88_c12_alpha_gamma = 3, &
           cf88_n14_alpha_gamma = 4, cf88_triple_alpha = 5, num_rate_laws = 5
! A reaction's screening constants: the charge product Z1 Z2 and the
! Graboske et al. (1973, table 4) factors of the intermediate (z86) and
! strong (z53, z43, z23) forms.
      type, public :: screening_t
         real(dp) :: zprod = 0.0d0, z86 = 0.0d0, z53 = 0.0d0, z43 = 0.0d0, &
              z23 = 0.0d0
      end type screening_t
! The screening plasma's constants (compute_screening_plasma):
!   (p_F/m_e c)^2 per (rho/mu_e)^(2/3), rho/mu_e in g/cm^3:
!     (3 pi^2 N_A)^(2/3) (hbar/m_e c)^2;
!   m_e c^2/k in 1e9 K, for E_F/kT;
!   Graboske et al.'s (1973) Lambda_0 = 1.88e8 (rho/mu_i)^(1/2) T^(-3/2)
!     with T in 1e9 K (1.88e8/(1e9)^1.5);
!   E_F/kT at or below which the electrons count as nondegenerate
!     (Salpeter's f'/f = 1), and log10(E_F/kT) from which they are fully
!     degenerate (f'/f = 0).
      real(dp), parameter :: fermi_momentum_sq_coeff = 1.017677d-4, &
           electron_rest_energy_t9 = 5.92986d0, lambda0_coeff = 5.9426d-6, &
           nondegenerate_ef_over_kt = 1.0d-2, degenerate_log_ef_over_kt = 1.5d0
! The weak-screening exponent Lambda_0 z_tilde Z1 Z2 up to which the
! intermediate form applies (above the weak threshold of the settings),
! and from which the strong form always does (compute_reaction_screening).
      real(dp), parameter :: max_intermediate_weak_exponent = 2.0d0, &
           min_strong_weak_exponent = 5.0d0
! The per-gram prefactors of the Caughlan & Fowler (1988) rates: their
! N_A<sigma v> (cm^3/mol/s; N_A^2<sigma v> for triple alpha) becomes the
! rate per gram per second with the abundance factors omitted through
! rho N_A/(A1 A2), and rho^2 N_A/(3! A_He4^3) for triple alpha, with the
! neutral atomic masses (4.002603, 12, 13.003355, 14.003074; VandenBerg's
! notes) -- unlike the bare nuclear masses of the screening and Be7.
      real(dp), parameter :: per_gram_c13_alpha = 1.157126d22, &
           per_gram_c12_alpha = 1.25388d22, per_gram_n14_alpha = 1.07452d22, &
           per_gram_triple_alpha = 1.565315d21
! N_A/m(Be7) per gram, the bare nuclear mass 7.014735 amu (JNB 10/14/97):
! the Be7 capture rates per unit Be7 abundance, per gram.
      real(dp), parameter :: per_gram_be7 = 8.584981d+22
! bare nuclear masses (amu) and charges of the iso_* isotopes
      real(dp), parameter, public :: &
           isotope_mass(num_isotopes) = [1.008665d0, 1.007276d0, 2.013553d0, 3.015501d0, 3.014933d0, 4.001506d0, &
                11.996709d0, 13.000064d0, 13.999233d0, 15.990526d0, 17.994772d0, 19.986954d0, 23.978458d0], &
           isotope_charge(num_isotopes) = [0.d0, 1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 6.d0, 6.d0, 7.d0, 8.d0, 8.d0, 10.d0, 12.d0]

! T9 = T/10^9 K and the powers of it the rate formulas use (p13 = T9^1/3,
! m23 = T9^-2/3, m12 = T9^-1/2, ...)
      type, public :: t9_powers
         real(dp) :: t9 = 0.0d0, p13 = 0.0d0, p23 = 0.0d0, &
              m13 = 0.0d0, m23 = 0.0d0, m1 = 0.0d0, m2 = 0.0d0, &
              m12 = 0.0d0, m32 = 0.0d0
      end type t9_powers

! The plasma quantities the screening formulas share: the composition
! moments (1/mu_i, 1/mu_e, sum Z^1.58 X/A, sum Z^2 X/A), Graboske et
! al.'s Lambda_0 and its powers, DeWitt et al.'s z_tilde, the mean ionic
! charge and its powers, and the weak-screening exponent per unit charge
! product.
      type :: screening_plasma
         real(dp) :: inv_mu_i = 0.0d0, inv_mu_e = 0.0d0, &
              charge_moment_158 = 0.0d0, charge_moment_2 = 0.0d0
         real(dp) :: lambda0 = 0.0d0, lambda0_23 = 0.0d0, &
              lambda0_086 = 0.0d0, z_tilde = 0.0d0, mean_charge = 0.0d0, &
              z_tilde_058 = 0.0d0, mean_charge_028 = 0.0d0, &
              mean_charge_13 = 0.0d0, weak_unit = 0.0d0
      end type screening_plasma

! What every rate of one zone needs (compute_zone_conditions): the
! density (g/cm^3), the electron density over N_A (rho/mu_e), T9 and its
! powers, and the screening plasma.
      type, public :: zone_conditions
         real(dp) :: density = 0.0d0, electron_density_na = 0.0d0
         type(t9_powers) :: t
         type(screening_plasma) :: plasma
      end type zone_conditions

! The fate of Be7 (compute_be7_branching): the electron- and
! proton-capture rates per gram per second, per unit Be7 abundance (the
! flux diagnostics), and the fraction burned each way.
      type, public :: be7_branching
         real(dp) :: electron = 0.0d0, proton = 0.0d0
         real(dp) :: frac_ecap = 0.0d0, frac_pcap = 0.0d0
      end type be7_branching

! The rates of one zone, from network's compute_reaction_rates: rate(r_*)
! per gram per second with the abundance factors omitted (screened,
! scaled), and its logarithmic derivatives; the Be7 branching and the
! N15(p,alpha) fraction; the screened hep factor (the flux diagnostics);
! and the zone conditions they came from, which the energy and flux
! consumers reuse. All zero below net_settings%logT_cutoffs(1).
      type, public :: nuclear_rates
         real(dp) :: rate(num_reactions) = 0.0d0
         real(dp) :: dlnrate_dlnt(num_reactions) = 0.0d0
         real(dp) :: dlnrate_dlnrho(num_reactions) = 0.0d0
         type(be7_branching) :: be7
         real(dp) :: frac_n15_pa = 0.0d0
         real(dp) :: hep = 0.0d0
         type(zone_conditions) :: zone
      end type nuclear_rates

      public :: screening_mass_fractions, compute_deuterium_proton_capture, &
           cap_deuterium_rate, powers_of_t9, compute_zone_conditions, &
           compute_screened_rate, compute_be7_branching, &
           n15_alpha_fraction, hep_rate_factor

contains

!----------------------------------------------------------------------
! compute_deuterium_proton_capture: the d(p,g)3He rate per gram per
! second excluding the abundance factors, rho N_A<sigma v> N_A/m_d
! (Caughlan & Fowler 1988, the T9^(-2/3) form: N_A<sigma v> = 2240
! T9^(-2/3) exp(-3.72 T9^(-1/3)) (1 + 0.112 T9^(1/3) + 3.38 T9^(2/3) +
! 2.65 T9)), and its exact d ln rate/d ln T. Unscreened: the burning
! temperatures are low.
pure subroutine compute_deuterium_proton_capture(density, t, rate, dlnrate_dlnt)
      use physical_constants, only: cc13
      use math_lib
      real(dp), intent(in) :: density
      type(t9_powers), intent(in) :: t
      real(dp), intent(out) :: rate, dlnrate_dlnt
! the fit's exponent and polynomial, and T9 d(polynomial)/dT9 times 3
      real(dp) :: zz, tfac, tfac_deriv

      zz = -3.72d0*t%m13
      tfac = 1.0d0+0.112d0*t%p13+3.38d0*t%p23+2.65d0*t%t9
      rate = density*2.240d3*t%m23*exp(zz)*tfac*avogadro_dpg/ &
           deuteron_mass_amu
      tfac_deriv = 0.112d0*t%p13+6.76d0*t%p23+7.95d0*t%t9
      dlnrate_dlnt = cc13*((tfac_deriv/tfac) -2.0d0 - zz)
end subroutine compute_deuterium_proton_capture

!----------------------------------------------------------------------
! cap_deuterium_rate: in the surface convection zone (in_surface_cz),
! deuterium burns no faster than the convection overturns: when the
! turnover time exceeds min_cap_turnover_time and X(H) exceeds
! min_cap_hydrogen, rate (per gram per second, abundance factors
! omitted) is capped so that rate X(H) <= N_A/(m_d tau). capped says
! whether the cap bound; a capped rate depends on neither rho nor T.
pure subroutine cap_deuterium_rate(hydrogen_fraction, in_surface_cz, &
     turnover_timescale, rate, capped)
      real(dp), intent(in) :: hydrogen_fraction, turnover_timescale
      logical, intent(in) :: in_surface_cz
      real(dp), intent(inout) :: rate
      logical, intent(out) :: capped
      real(dp) :: max_rate

      capped = .false.
      if (.not. in_surface_cz .or. turnover_timescale <= min_cap_turnover_time) return
      max_rate = avogadro_dpg/deuteron_mass_amu/turnover_timescale
      if (rate*hydrogen_fraction > max_rate .and. hydrogen_fraction > min_cap_hydrogen) then
         rate = max_rate/hydrogen_fraction
         capped = .true.
      end if
end subroutine cap_deuterium_rate

!----------------------------------------------------------------------
! compute_zone_conditions: what every rate of one zone needs -- the
! density and the electron density, T9 and its powers, and the plasma
! quantities of the screening -- for network's compute_reaction_rates.
! mass_frac(iso_*) is the composition the screening sees
! (screening_mass_fractions).
pure subroutine compute_zone_conditions(log_density, log_temperature, &
     mass_frac, zone)
      use physical_constants, only: ln10
      use math_lib
      real(dp), intent(in) :: log_density, log_temperature, &
           mass_frac(num_isotopes)
      type(zone_conditions), intent(out) :: zone

      call compute_composition_moments(mass_frac, zone%plasma)
      zone%density = exp(ln10*log_density)
      zone%electron_density_na = zone%density*zone%plasma%inv_mu_e
      zone%t = powers_of_t9(log_temperature)
      call compute_screening_plasma(zone%density, zone%electron_density_na, &
           zone%t, zone%plasma)
end subroutine compute_zone_conditions

!----------------------------------------------------------------------
! The composition moments the screening needs: 1/mu_i = sum X/A,
! 1/mu_e = sum Z X/A, sum Z^1.58 X/A (intermediate screening, Graboske et
! al.: the mean of Z^(3b-1)) and sum Z^2 X/A (Salpeter's zeta).
pure subroutine compute_composition_moments(mass_frac, plasma)
      use math_lib
      real(dp), intent(in) :: mass_frac(num_isotopes)
      type(screening_plasma), intent(out) :: plasma
      real(dp) :: per_amu
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
pure function powers_of_t9(log_temperature) result(t)
      use physical_constants, only: ln10, cc13
      use math_lib
      real(dp), intent(in) :: log_temperature
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
! The plasma quantities of a zone's screening: Salpeter's degeneracy
! correction f'/f, Graboske et al.'s Lambda_0, DeWitt et al.'s z_tilde,
! the mean ionic charge, and the powers of them the screening forms use
! (compute_reaction_screening, weak_intermediate_screening).
pure subroutine compute_screening_plasma(density, electron_density_na, t, &
     plasma)
      use physical_constants, only: cc13
      use math_lib
      real(dp), intent(in) :: density, electron_density_na
      type(t9_powers), intent(in) :: t
      type(screening_plasma), intent(inout) :: plasma
! (p_F/m_e c)^2, E_F/kT, Salpeter's f'/f and log10(E_F/kT)
      real(dp) :: fermi_momentum_sq, fermi_energy_over_kt, &
           degeneracy_factor, log_degeneracy

      fermi_momentum_sq=fermi_momentum_sq_coeff*pow(electron_density_na, (2.0d0/3.0d0))
      fermi_energy_over_kt=electron_rest_energy_t9*t%m1*(sqrt(1.d0+fermi_momentum_sq)-1.d0)
      if (fermi_energy_over_kt <= nondegenerate_ef_over_kt) then
         degeneracy_factor=1.0d0
      else
! Salpeter's figure 1, as a function of the degeneracy
         log_degeneracy=log10(fermi_energy_over_kt)
         if (log_degeneracy >= degenerate_log_ef_over_kt) then
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
      plasma%lambda0=lambda0_coeff*t%m32*sqrt(density*plasma%inv_mu_i)
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
end subroutine compute_screening_plasma

!----------------------------------------------------------------------
! The screening exponent of one reaction and its derivatives: weak
! (Salpeter 1954, with his degeneracy correction) below
! net_settings%weak_screening_threshold (0.03 reproduces Graboske et al.
! and Salpeter; a value above 1 forces weak screening in the Sun),
! intermediate (Graboske et al. 1973, table 4) up to a weak-screening
! value of max_intermediate_weak_exponent (2), and beyond that strong or
! intermediate, whichever is smaller (strong always from
! min_strong_weak_exponent, 5). (Strong screening was never checked:
! it does not occur in the Sun.)
pure subroutine compute_reaction_screening(plasma, s, screen, &
     dscreen_dlnrho, dscreen_dlnt)
      type(screening_plasma), intent(in) :: plasma
      type(screening_t), intent(in) :: s
      real(dp), intent(out) :: screen, dscreen_dlnrho, dscreen_dlnt
      real(dp) :: weak, intermediate, strong

      weak=plasma%weak_unit*s%zprod
      if (weak <= net_settings%weak_screening_threshold) then
         screen=weak
         dscreen_dlnrho=0.5d0*weak
         dscreen_dlnt=-1.5d0*weak
      else
         intermediate=0.38d0*plasma%lambda0_086*plasma%charge_moment_158*s%z86/ &
              (plasma%inv_mu_i*plasma%z_tilde_058*plasma%mean_charge_028)
         if (weak <= max_intermediate_weak_exponent) then
            screen=intermediate
            dscreen_dlnrho=0.43d0*intermediate
            dscreen_dlnt=-1.29d0*intermediate
         else
            strong=0.624d0*plasma%mean_charge_13*plasma%lambda0_23*(s%z53+ &
                 0.316d0*plasma%mean_charge_13*s%z43+0.737d0* &
                 s%z23/(plasma%mean_charge*plasma%lambda0_23))
            if (strong < intermediate.or.weak >= min_strong_weak_exponent) then
               screen=strong
               dscreen_dlnrho=0.208d0*plasma%mean_charge_13*(s%z53+ &
                    0.316d0*plasma%mean_charge_13*s%z43)*plasma%lambda0_23
               dscreen_dlnt=-3.d0*dscreen_dlnrho
            else
               screen=intermediate
               dscreen_dlnrho=0.43d0*intermediate
               dscreen_dlnt=-1.29d0*intermediate
            end if
         end if
      end if
end subroutine compute_reaction_screening

!----------------------------------------------------------------------
! Weak screening, or intermediate above the weak threshold, for one
! reaction of charge product zprod and Graboske z86 (Be7 + p, hep).
pure real(dp) function weak_intermediate_screening(plasma, zprod, z86)
      type(screening_plasma), intent(in) :: plasma
      real(dp), intent(in) :: zprod, z86
      real(dp) :: weak

      weak = plasma%weak_unit*zprod
      if (weak <= net_settings%weak_screening_threshold) then
         weak_intermediate_screening = weak
      else
         weak_intermediate_screening = 0.38d0*plasma%lambda0_086*plasma%charge_moment_158*z86/ &
              (plasma%inv_mu_i*plasma%z_tilde_058*plasma%mean_charge_028)
      end if
end function weak_intermediate_screening

!----------------------------------------------------------------------
! compute_screened_rate: the screened rate of one reaction of network's
! table, per gram per second with the abundance factors omitted, and
! its logarithmic derivatives, from its row's rate law, fit
! coefficients, reaction slot (the S-factor law's scales are the
! settings' at that slot) and screening constants. Each law applies the
! screening exponent where its formula has always put it (inside the
! fit's exponential for the S-factor and triple-alpha laws).
pure subroutine compute_screened_rate(rate_law, fit, slot, s, zone, &
     rate, dlnrate_dlnrho, dlnrate_dlnt)
      integer, intent(in) :: rate_law, slot
      real(dp), intent(in) :: fit(:)
      type(screening_t), intent(in) :: s
      type(zone_conditions), intent(in) :: zone
      real(dp), intent(out) :: rate, dlnrate_dlnrho, dlnrate_dlnt
      real(dp) :: screen, dscreen_dlnrho, dscreen_dlnt

      call compute_reaction_screening(zone%plasma, s, screen, dscreen_dlnrho, &
           dscreen_dlnt)
      select case (rate_law)
      case (s_factor_fit)
         call compute_s_factor_rate(fit, slot, zone%density, zone%t, &
              screen, dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, &
              dlnrate_dlnt)
      case (cf88_c13_alpha_n)
         call compute_c13_alpha_n_rate(zone%density, zone%t, screen, &
              dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      case (cf88_c12_alpha_gamma)
         call compute_c12_alpha_gamma_rate(zone%density, zone%t, screen, &
              dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      case (cf88_n14_alpha_gamma)
         call compute_n14_alpha_gamma_rate(zone%density, zone%t, screen, &
              dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      case (cf88_triple_alpha)
         call compute_triple_alpha_rate(zone%density, zone%t, screen, &
              dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      case default
         rate = 0.0d0
         dlnrate_dlnrho = 0.0d0
         dlnrate_dlnt = 0.0d0
      end select
end subroutine compute_screened_rate

!----------------------------------------------------------------------
! The S-factor law (reactions 1-7: pp, He3+He3, He3+He4, p+C12, p+C13,
! p+N14, p+O16):
! rate = rho * [T9^(-2/3) S_eff/S0] * exp(tau + q7 + (q8 T9)^2 + screening)
! (Neutrino Astrophysics eq. 3.14, Fowler, Caughlan & Zimmerman 1967 eq.
! 52), with S_eff/S0 = 1 + q1 T9^(1/3) + q2.. q5 terms and tau = q6
! T9^(-1/3). The reaction's multipliers in net_settings, at its slot:
! cross_section_scale on the rate, s_prime_scale and
! s_double_prime_scale on the first and second energy derivatives of S
! (the q2, q3 and the q4, q5 terms). The T-derivative is exact for this
! form. (2026-10: the q3 term lacked the first-derivative scale and the
! q4, q5 terms the second.)
pure subroutine compute_s_factor_rate(q, slot, density, t, screen, &
     dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      use math_lib
      real(dp), intent(in) :: q(:), density, screen, dscreen_dlnrho, &
           dscreen_dlnt
      integer, intent(in) :: slot
      type(t9_powers), intent(in) :: t
      real(dp), intent(out) :: rate, dlnrate_dlnrho, dlnrate_dlnt
! T9^(-2/3) S_eff/S0, and the three multipliers
      real(dp) :: s_ratio, qs0e, qqs0ee, xs

      qs0e = net_settings%s_prime_scale(slot)
      qqs0ee = net_settings%s_double_prime_scale(slot)
      xs = net_settings%cross_section_scale(slot)
      s_ratio=t%m23+q(1)*t%m13+qs0e*(q(2)+q(3)*t%p13)+ &
           qqs0ee*(q(4)*t%p23+q(5)*t%t9)
      rate=density*s_ratio*exp(q(6)*t%m13+q(7)+(q(8)*t%t9)**2+screen)
      rate = rate*xs
      dlnrate_dlnrho=1.d0+dscreen_dlnrho
      dlnrate_dlnt=dscreen_dlnt-(q(6)*t%m13+(2.d0*t%m23+ &
           q(1)*t%m13-qs0e*q(3)*t%p13- &
           qqs0ee*(2.d0*q(4)*t%p23+3.d0*q(5)*t%t9))/s_ratio)/3.d0+ &
           2.d0*(q(8)*t%t9)**2
end subroutine compute_s_factor_rate

!----------------------------------------------------------------------
! The fate of Be7: electron capture (Neutrino Astrophysics eq. 3.18) or
! proton capture (table 3.2, eq. 3.12, with the S_eff/S0 corrections,
! their scales be7p_s_prime_scale and be7p_s_double_prime_scale, and
! weak/intermediate screening), both per gram (rho over
! the bare Be7 mass, JNB 10/14/97), and the fraction burned each way.
pure subroutine compute_be7_branching(zone, hydrogen_fraction, be7)
      use math_lib
      type(zone_conditions), intent(in) :: zone
      real(dp), intent(in) :: hydrogen_fraction
      type(be7_branching), intent(out) :: be7
      real(dp) :: tau_be7p, s_ratio_be7p, per_gram

      associate (density => zone%density, t => zone%t, plasma => zone%plasma)

      be7%electron = (1.752d-10)*t%m12*(1.0d0 + 0.004d0*(1000.d0*t%t9 - 16.d0))
      be7%electron = be7%electron*plasma%inv_mu_e*net_settings%extra_cross_section_scale(xs_be7e)
      tau_be7p = (-10.2625d0*t%m13)
      be7%proton = (3.128813d+5)*hydrogen_fraction*net_settings%extra_cross_section_scale(xs_be7p)*exp(tau_be7p)
      s_ratio_be7p = t%m23 + be7_pcapture_fit(1)*t%m13 + &
           net_settings%be7p_s_prime_scale*(be7_pcapture_fit(2)+ be7_pcapture_fit(3)*t%p13) &
           + net_settings%be7p_s_double_prime_scale*(be7_pcapture_fit(4)*t%p23 + be7_pcapture_fit(5)*t%t9)
      be7%proton = be7%proton*s_ratio_be7p
      be7%proton = be7%proton*exp(weak_intermediate_screening(plasma, 4.0d0, 5.7790d0))
      per_gram = density*per_gram_be7
      be7%proton = per_gram*be7%proton
      be7%electron = per_gram*be7%electron
      be7%frac_ecap = be7%electron/(be7%electron + be7%proton)
      be7%frac_pcap = be7%proton/(be7%electron + be7%proton)
      end associate
end subroutine compute_be7_branching

!----------------------------------------------------------------------
! The fate of N15 + p: the fraction that goes (p,alpha)C12 rather than
! (p,gamma)O16, in the ratio of their effective S-factors (S0 in keV-b
! times the S_eff/S0 polynomial; Neutrino Astrophysics table 3.4,
! S0(p,alpha) = 67500 keV-b since the Solar Fusion workshop, JNB
! 9/25/97), each with its scale (MHP 9/14).
pure real(dp) function n15_alpha_fraction(zone)
      type(zone_conditions), intent(in) :: zone
      real(dp) :: pg_strength, pa_strength

      associate (t => zone%t)

      pg_strength = t%m23 + 0.0273016d0*t%m13 + 0.14374d0 + 0.027490d0*t%p13 &
                + 6.14685d0*t%p23 + 2.98940d0*t%t9
      pg_strength = pg_strength*s0_n15pg_kevb*net_settings%n15_pg_scale
      pa_strength = t%m23 + 0.0273016d0*t%m13 + 2.01186d0 + 0.384763d0*t%p13 &
                + 17.0579d0*t%p23 + 8.29580d0*t%t9
      pa_strength = pa_strength*s0_n15pa_kevb*net_settings%n15_pa_scale
      n15_alpha_fraction = pa_strength/(pa_strength + pg_strength)
      end associate
end function n15_alpha_fraction

!----------------------------------------------------------------------
! C13(alpha,n)O16 (Caughlan & Fowler 1988, MHP 3/91; the prefactor from
! VandenBerg's notes), with its exact derivatives: a1..a5 are the terms
! of the CF88 fit and r1 its polynomial factor.
pure subroutine compute_c13_alpha_n_rate(density, t, screen, &
     dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      use physical_constants, only: cc13
      use math_lib
      real(dp), intent(in) :: density, screen, dscreen_dlnrho, dscreen_dlnt
      type(t9_powers), intent(in) :: t
      real(dp), intent(out) :: rate, dlnrate_dlnrho, dlnrate_dlnt
      real(dp) :: r1, a1, a2, a3, a4, a5, dr1, da1, s_sum

      r1=t%m23+0.0129d0*t%m13+2.04d0+0.184d0*t%p13
      a1 = 6.77d15*exp(-32.329d0*t%m13-(t%t9/1.284d0)**2)
      a2 = 3.82d5*exp(-9.373d0*t%m1)
      a3 = 1.41d6*exp(-11.873d0*t%m1)
      a4 = 2.00d9*exp(-20.409d0*t%m1)
      a5 = 2.92d9*exp(-29.283d0*t%m1)
      s_sum = a1*r1+t%m32*(a2+a3+a4+a5)
      rate = per_gram_c13_alpha*density*exp(screen)*s_sum
      dlnrate_dlnrho=1.0d0+dscreen_dlnrho
      dr1 = cc13*(-2.0d0*t%m23-0.0129d0*t%m13+0.184d0*t%p13)
      da1 = a1*(cc13*32.329d0*t%m13 - 2.0d0*(t%t9/1.284d0)**2)
      dlnrate_dlnt = dscreen_dlnt+(dr1*a1 + r1*da1 + &
           t%m32*(a2*(9.373d0*t%m1-1.5d0)+a3*(11.873d0*t%m1-1.5d0)+ &
           a4*(20.409d0*t%m1-1.5d0)+a5*(29.283d0*t%m1-1.5d0)))/s_sum
end subroutine compute_c13_alpha_n_rate

!----------------------------------------------------------------------
! C12(alpha,gamma)O16 (Caughlan & Fowler 1988), with its exact
! derivatives: a1..a5 are the terms of the CF88 fit and r1, r2 its
! polynomial factors.
pure subroutine compute_c12_alpha_gamma_rate(density, t, screen, &
     dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      use physical_constants, only: cc13
      use math_lib
      real(dp), intent(in) :: density, screen, dscreen_dlnrho, dscreen_dlnt
      type(t9_powers), intent(in) :: t
      real(dp), intent(out) :: rate, dlnrate_dlnrho, dlnrate_dlnt
      real(dp) :: r1, r2, a1, a2, a3, a4, a5, s_sum

      r1 = 1.0d0/(1.0d0+0.0489d0*t%m23)
      r2 = 1.0d0/(1.0d0+0.2654d0*t%m23)
      a1 = t%m2*exp(-32.120d0*t%m13)
      a2 = 1.04d8*r1**2*exp(-(t%t9/3.496d0)**2)
      a3 = 1.76d8*r2**2
      a4 = 1.25d3*t%m32*exp(-27.499d0*t%m1)
      a5 = 1.43d-2*t%t9**5*exp(-15.541d0*t%m1)
      s_sum = a1*(a2+a3)+a4+a5
      rate = per_gram_c12_alpha*density*exp(screen)*s_sum
      dlnrate_dlnrho = 1.0d0+dscreen_dlnrho
      dlnrate_dlnt = dscreen_dlnt+ &
           (a1*((cc13*32.120d0*t%m13-2.0d0)*(a2+a3)+ &
           a2*(cc13*0.1956d0*t%m23*r1-2.0d0*(t%t9/3.496d0)**2)+ &
           a3*(cc13*1.0616d0*t%m23*r2))+ &
           a4*(27.499d0*t%m1-1.5d0)+a5*(5.0d0+15.541d0*t%m1))/s_sum
end subroutine compute_c12_alpha_gamma_rate

!----------------------------------------------------------------------
! N14(alpha,gamma)F18(e+ nu)O18 (Caughlan & Fowler 1988), with its exact
! derivatives: a1..a4 are the terms of the CF88 fit and r1 its
! polynomial factor. Linear coefficient 0.117 (2026 bugsweep Batch 3:
! both former copies carried 0.177).
pure subroutine compute_n14_alpha_gamma_rate(density, t, screen, &
     dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      use physical_constants, only: cc13, cc23
      use math_lib
      real(dp), intent(in) :: density, screen, dscreen_dlnrho, dscreen_dlnt
      type(t9_powers), intent(in) :: t
      real(dp), intent(out) :: rate, dlnrate_dlnrho, dlnrate_dlnt
      real(dp) :: r1, a1, a2, a3, a4, dr1, da1, s_sum

      r1 = t%m23+0.012d0*t%m13+1.45d0+0.117d0*t%p13+1.97d0*t%p23 &
           +0.406d0*t%t9
      a1 = 7.78d9*exp(-36.031d0*t%m13-(t%t9/0.881d0)**2)
      a2 = t%m32*2.36d-10*exp(-2.798d0*t%m1)
      a3 = t%m32*2.03d0*exp(-5.054d0*t%m1)
      a4 = t%m23*1.15d4*exp(-12.310d0*t%m1)
      s_sum = a1*r1+a2+a3+a4
      rate = per_gram_n14_alpha*density*exp(screen)*s_sum
      dlnrate_dlnrho=1.d0+dscreen_dlnrho
      dr1 = cc13*(-2.0d0*t%m23-0.012d0*t%m13+0.117d0*t%p13+ &
            3.94d0*t%p23)+0.406d0*t%t9
      da1 = a1*(cc13*36.031d0*t%m13-2.0d0*(t%t9/0.881d0)**2)
      dlnrate_dlnt = dscreen_dlnt+(dr1*a1+r1*da1+ &
              a2*(2.798d0*t%m1-1.5d0)+a3*(5.054d0*t%m1-1.5d0)+ &
              a4*(12.310d0*t%m1-cc23))/s_sum
end subroutine compute_n14_alpha_gamma_rate

!----------------------------------------------------------------------
! Triple alpha (Caughlan & Fowler 1988), with its exact derivatives.
pure subroutine compute_triple_alpha_rate(density, t, screen, &
     dscreen_dlnrho, dscreen_dlnt, rate, dlnrate_dlnrho, dlnrate_dlnt)
      use math_lib
      real(dp), intent(in) :: density, screen, dscreen_dlnrho, dscreen_dlnt
      type(t9_powers), intent(in) :: t
      real(dp), intent(out) :: rate, dlnrate_dlnrho, dlnrate_dlnt

      rate = per_gram_triple_alpha*density**2*t%m1*t%m2*2.79d-8* &
                 exp(-4.4027d0*t%m1+screen)
      dlnrate_dlnrho = 2.0d0+dscreen_dlnrho
      dlnrate_dlnt = -3.0d0+dscreen_dlnt+4.4027d0*t%m1
end subroutine compute_triple_alpha_rate

!----------------------------------------------------------------------
! hep, He3(p,e+ nu)He4, for the neutrino flux only: the screened, scaled
! rate factor, which the consumer multiplies by X(H) X(He3).
pure real(dp) function hep_rate_factor(zone)
      use math_lib
      type(zone_conditions), intent(in) :: zone

      associate (density => zone%density, t => zone%t, plasma => zone%plasma)

      hep_rate_factor = (1.71724d+11)*density*t%m23*exp(-6.1399d0*t%m13)
      hep_rate_factor = (1.0d0 + 0.067862d0*t%p13)*net_settings%extra_cross_section_scale(xs_hep)*hep_rate_factor
      hep_rate_factor = hep_rate_factor*exp(weak_intermediate_screening(plasma, 2.0d0, 3.08687d0))
      end associate
end function hep_rate_factor

!----------------------------------------------------------------------
! screening_mass_fractions: the isotope mass fractions (iso_* slots) that
! compute_reaction_rates screens with, from a zone's composition in the
! star%xa layout. The isotopes YREC does not track (n, H3, Ne20, Mg24)
! are zero; deuterium is included (2026-10: the composition path used to
! leave it out).
pure function screening_mass_fractions(xa) result(mass_frac)
      use composition_layout, only: i_h1, i_h2, i_he4, i_he3, i_c12, i_c13, &
           i_n14, i_o16, i_o18
      real(dp), intent(in) :: xa(:)
      real(dp) :: mass_frac(num_isotopes)

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

end module reaction_rates

!======================================================================
! Sources and history of the rate formulas (rates.f and engeb.f, J. N.
! Bahcall, M. H. Pinsonneault, D. B. Guenther, 1990-2014). The comments
! of the F77 code are kept verbatim in net/retired/f77_rate_comments.txt.
!
! Legacy names. Q1..Q8 are the fit coefficients of network's table rows
! (fit) and be7_pcapture_fit; ZPRD, Z86, Z53, Z43, Z23 the rows'
! screening constants (screening_t); SSTANDARD the cross-section scales
! (net_configuration, which maps its 17 entries); QS0E, QQS0EE the
! derivative scales; XXL, XXL6, XXL8 lambda0, lambda0_23, lambda0_086;
! ZCURL z_tilde, ZBAR mean_charge, AMU and EMU inv_mu_i and inv_mu_e, XTR
! and ZET the charge moments; PFMC2, EFMKT, FPRF, DEGD fermi_momentum_sq,
! fermi_energy_over_kt, degeneracy_factor, log_degeneracy; UWK, UINT,
! USTR the weak, intermediate and strong exponents; R1 s_ratio; TEMP3
! tau_be7p, QRBE7 s_ratio_be7p, CAMUBE7 per_gram_be7; O16GAMMA and
! C12ALPHA pg_strength and pa_strength; F1 frac_be7_ecap, F3 frac_n15_pa.
!
! Every number was checked and revised by Bahcall to agree with Neutrino
! Astrophysics (Bahcall 1989; NA below) or later results; many differ
! from the original Yale code's.
!
! The S-factor law. T9^(-2/3) S_eff/S0 = T9^(-2/3) + q1 T9^(-1/3) + q2 +
! q3 T9^(1/3) + q4 T9^(2/3) + q5 T9 (NA eq. 3.14; Fowler, Caughlan &
! Zimmerman 1967, eqs. 52-53), the q's independent of T (NA eqs. 3.10,
! 3.11) and recomputed from NA tables 3.2 and 3.4. q6 = -tau(T9 = 1), tau
! = q6 T9^(-1/3) (NA eq. 3.10; negative). q7 = 70.62860 + ln(Z0 Z1/A)/3 -
! ln(A0 A1) - ln S0 - ln(1 + delta_01), S0 in keV b, delta_01 = 1 for
! identical nuclei. q8 multiplies T9^2 in the exponent (Harris et al.
! 1983, ARA&A 21, 165). The values are BP00's: DBG applied MHP's update
! (8/94); JNB recomputed them for bare nuclear masses with cues.f
! (10/13/97; there q6 = 0.1 times cues.f's tau(T6 = 1)). The S'/S0 and S''/S0 terms are
! scaled apart from S0 (MHP 8/14 and 9/14, which also restored the
! Be7 + p derivative terms that the 7/21/03 revision had zeroed).
!
! Masses. Bare nuclear masses since 10/13/97 (JNB: the neutral atomic
! mass less Z m_e c^2, from the mass excess over 931.49432 MeV plus A),
! which changed the q's, the Be7 rates and the N15 + p branching. The
! Caughlan & Fowler prefactors keep the neutral masses.
!
! The cross-section scales. SSTANDARD(i), the ratio of the wanted S0 to
! NA's (tables 3.2 and 3.4), multiplies the rate; the energy derivatives
! enter divided by S0, so changing only S(0) needs nothing else. Its
! preferred values were those of Bahcall & Pinsonneault (1992, RMP 64,
! 885, table 1), then (1995, RMP 67, 781). Entry 14, pep, scales the
! ratio of pep to pp: leave it at 1 unless that ratio is meant to change.
! Entry 17, hep: the code's fit has NA's S0 = 8.0e-20 keV b (s0_hep_ref
! in setup/user_inputs.f90); earlier choices of SSTANDARD(17) were 1.913
! (Wolfs et al. 1989, PRL 63, 2721: 15.3e-20 keV b) and 0.1625 (Carlson et
! al. 1991, PRC 44, 619: 1.3e-20 keV b). Entries 8 and 10-12 (the
! helium-burning rates) and 9 and 13 (O16 + alpha and C12 + C12, which
! were never computed) have no effect.
!
! Screening. Weak screening is Salpeter's (1954, Aust. J. Phys. 7, 373,
! eq. 25) with his degeneracy correction f'/f, fitted to his figure 1 in
! log10(E_F/kT); Bahcall refined the coefficients of fermi_momentum_sq
! and fermi_energy_over_kt and corrected the Yale code's electron density
! (it divided by 1/mu_e instead of multiplying). Intermediate and strong
! screening are those of Graboske, DeWitt, Grossman & Cooper (1973, ApJ
! 181, 457, table 4 and eq. 19), with z_tilde from DeWitt, Graboske &
! Cooper (1973, ApJ 181, 439, eq. 4). weak_screening_threshold = 0.03
! gives the standard results of Graboske et al. and Salpeter (the
! standard solar model's value); above 1 (30, say) the Sun is screened
! weakly throughout. Graboske et al. put the transition at Lambda_12 =
! 0.1 (p. 465; 5/15/97). The weak and intermediate forms were checked;
! the strong form and its constants Z53, Z43 and Z23 were not, as the Sun
! has no strong screening (test_net now checks the strong form's
! derivatives).
!
! The O16 + p rate of the Yale code carried a factor T9^(-1/7), which was
! wrong and was removed.
!
! Be7. Electron capture per NA eq. 3.18 (1/mu_e is EMU), proton capture
! per NA table 3.2 and eq. 3.12, with the S_eff/S0 corrections of
! be7_pcapture_fit and weak or intermediate screening; both per gram
! through rho/(m_u A_Be7) (inserted 11/6/1990), with the bare mass
! 7.014735 since 10/14/97 (7.016930 and the factor 8.582295e22 before;
! tau_be7p's coefficient was -10.26202). The S0 change of the electron
! capture is carried by its scale, not its coefficient (3.126571e5
! before). The Yale code's Be7 expression had a physically wrong term and
! old cross sections.
!
! N15 + p. The branching between (p,alpha)C12 and (p,gamma)O16 is the
! ratio of their effective S-factors, S0 times the S_eff/S0 polynomial
! (NA table 3.4; Bahcall & Ulrich 1988, RMP 60, table 21), whose
! coefficients JNB computed with cues.f (10/14/97, checked against
! energy.f with the older CNO data); new cross sections need new
! coefficients there. S0(p,alpha) = 67500 keV b since 9/25/97 (the Solar
! Fusion workshop; 78000 before; the F77 comment's "675000" is a typo),
! S0(p,gamma) = 64 keV b; their scales are MHP's (9/14).
!
! Alpha captures. The Caughlan & Fowler (1988) rates (MHP 3/91), in the
! form Caughlan & Zimmerman used, with the prefactors from VandenBerg's
! notes. The linear coefficient of the N14(a,g)F18 polynomial is 0.117
! (1 + 0.012 T9^(1/3) + 1.45 T9^(2/3) + 0.117 T9 + 1.97 T9^(4/3) + 0.406
! T9^(5/3)); both former copies had 0.177, a digit transposition from the
! F77 source, about 0.5% of the rate at T9 = 0.15 (2026 bug sweep).
!
! Units. The composition rate vector is per s per amu (2026-10; per
! 1e9 sidereal years per amu, the C21 factor, before). rates.f (MHP
! 7/91) was a stripped-down copy of engeb's rates for kemcom; YREC now
! has one rate evaluation.
!======================================================================
