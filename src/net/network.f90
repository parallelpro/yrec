!----------------------------------------------------------------------
! network
!----------------------------------------------------------------------
! The hydrogen- and helium-burning network as data, and the one implicit
! integrator that reads it (2026-10; replaces the hand-written
! 7-species equations and Jacobian of kemcom, mixing/burn_network_region).
!
! Species. The implicit solve carries H, He3, He4, C12, C13, N14 and
! O16 (in that order -- the order of the old system, which fixes the
! pivoting). O18 is a product only (of N14(alpha,g)) and is updated
! explicitly from the converged abundances. N15, O17, Be7 and the
! intermediate nuclei of each chain are taken to be in equilibrium and
! do not appear; deuterium and the light elements are burned elsewhere.
!
! Rates. Each row names how its rate is computed: a rate law of
! reaction_rates (the S-factor fit, with the row's coefficients, or one
! of the Caughlan & Fowler alpha-capture forms), its cross-section scale
! slot, its screening constants and whether it burns with hydrogen or in
! helium burning. compute_reaction_rates evaluates every row that owns a
! rate; a second outcome of the same reaction (shares_rate) uses the
! rate of the row above. The rates land at the rows' slots r_*.
!
! Reactions. Reaction j proceeds at
!     R_j = b_j * g_j * prod_m X(reactants_j(m))    [per s per amu]
! where g_j is row rate_row of the composition rate vector
! (compute_network_rates: N_A<sigma v>, density, screening and the 1/A
! factors are already in it), X are mass fractions, and b_j is 1, or
! the N15(p,alpha) branching fraction f (or 1-f) carried in row
! rv_frac_n15_pa, for the two outcomes of the N14(p,g) chain. Each
! reaction changes the per-amu abundance of species i by nu_ij, so
!     dX_i/dt = A_i * sum_j nu_ij R_j
! with A_i the integer mass number, the convention the rate rows use.
!
! Energy. Each row also carries the energy it deposits per reaction,
! q_mev, net of the neutrinos it emits (Bahcall & Ulrich 1988, table 21,
! as revised by Bahcall 1997), and the chain it belongs to. The nuclear
! energy generation sums R_j q_j over the table (compute_reaction_energies), so
! a reaction's composition change and its energy come from one row.
! (JNB 9/25/97 changed the pp energy release by 0.002 MeV, from a better
! estimate of the neutrino energy loss; his notes, vol. 19, pp. 139-141,
! document other keV-sized changes, all but the rare 8B. The C12 + p
! sequence and C13 + p are broken down separately.)
! He3(He4,g)Be7 appears twice, once per Be7 branch (electron capture:
! pp II; proton capture: pp III): the same change in abundances, with
! the two branch fractions and their own energies.
!
! Implicit step. Over dt the solve finds X with
!     F_i(X) = X_i - X0_i - dt A_i sum_j nu_ij R_j(X) = 0
! by Newton iteration on J_ik = dF_i/dX_k.
!
! Adding a reaction means adding a row to the table, and naming its rate
! slot r_* in reaction_rates (raising num_reactions); a rate of a new
! functional form also needs its rate law there. Adding a species means
! adding it to the species list and to the rows that touch it.
module network
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: i_h1, i_he3, i_he4, i_c12, i_c13, i_n14, &
           i_o16, i_o18
      use reaction_rates, only: r_pp, r_he3he3, r_he3he4, &
           r_pc12, r_pc13, r_pn14, r_po16, r_ac13, &
           r_ac12, r_an14, r_3alpha, rv_frac_n15_pa, rv_frac_be7_ecap, &
           num_reactions, num_rate_rows, iso_h1, num_isotopes, &
           nuclear_rates, zone_conditions, screening_t, shares_rate, &
           s_factor_fit, cf88_c13_alpha_n, cf88_c12_alpha_gamma, &
           cf88_n14_alpha_gamma, cf88_triple_alpha, net_settings, &
           compute_zone_conditions, compute_screened_rate, &
           compute_be7_branching, compute_n15_branching, compute_hep_rate, &
           screening_mass_fractions
      implicit none
      private

! network species: positions in the abundance vector of the solve
      integer, parameter, public :: s_h1 = 1, s_he3 = 2, s_he4 = 3, &
           s_c12 = 4, s_c13 = 5, s_n14 = 6, s_o16 = 7, s_o18 = 8
! the first n_solved species are solved implicitly; the rest are
! products updated explicitly after convergence
      integer, parameter, public :: n_solved = 7, n_network = 8
! the composition-vector (xa) slot of each network species, s_* order
      integer, parameter, public :: network_species_slot(n_network) = &
           [i_h1, i_he3, i_he4, i_c12, i_c13, i_n14, i_o16, i_o18]
      real(dp), parameter, public :: mass_number(n_network) = &
           [1.0d0, 3.0d0, 4.0d0, 12.0d0, 13.0d0, 14.0d0, 16.0d0, 18.0d0]
! atomic mass excesses, MeV (AME2020: Wang et al. 2021), of the network
! species and of the neutron, the one untracked particle a row emits;
! they check the rows' energies (q_from_masses)
      real(dp), parameter :: mass_excess_mev(n_network) = &
           [7.288971d0, 14.931219d0, 2.424916d0, 0.0d0, 3.125009d0, &
           2.863417d0, -4.737001d0, -0.782816d0]
      real(dp), parameter :: neutron_mass_excess_mev = 8.071318d0

! branch factor b_j: 1, the N15(p,alpha) fraction f or 1 - f, or the Be7
! electron-capture fraction f_e or 1 - f_e
      integer, parameter :: unbranched = 0, n15_to_c12 = 1, n15_to_o16 = 2, &
           be7_electron_capture = 3, be7_proton_capture = 4

! when a reaction burns: with hydrogen present, or at helium-burning
! temperatures (log T >= logT_cutoffs(4))
      integer, parameter :: hydrogen_burning = 1, helium_burning = 2

! the chain a reaction's energy is counted under
      integer, parameter, public :: chain_pp1 = 1, chain_pp2 = 2, &
           chain_pp3 = 3, chain_cno = 4, chain_alpha_capture = 5, &
           chain_triple_alpha = 6, n_chains = 6

! One reaction. reactants lists the species whose mass fractions multiply
! the rate (repeat a species for a square or cube; 0 = unused). change lists
! (species, nu) pairs: the change in per-amu abundance per reaction
! (0 = unused). emitted_mass is the mass number carried off by particles
! the network does not track (the neutron of C13(a,n)O16), so that every
! row conserves baryon number exactly: sum_i A_i nu_i + emitted_mass = 0.
! q_mev is the energy released (from the atomic masses, so positron
! annihilation included) less the mean energy of the row's neutrinos:
! never more than the mass-energy released. network_table_errors checks
! both; net/test/test_net pins the table. The rate: rate_law (a rate law
! of reaction_rates, or shares_rate), burns_in (hydrogen_burning or
! helium_burning), fit (the S-factor law's q1..q8), scale_slot (its
! multipliers in net_settings' cross_section_scale, qs0e_scale and
! qqs0ee_scale; 0: none) and screening (Z1 Z2 and the Graboske factors).
      integer, parameter :: max_reactants = 3, max_change = 5
      type :: reaction_t
         character(len=24) :: name
         integer :: rate_row
         integer :: branch
         integer :: reactants(max_reactants)
         integer :: change(2*max_change)
         integer :: emitted_mass = 0
! energy deposited per reaction, MeV (net of neutrino losses), and the
! chain it is counted under
         real(dp) :: q_mev = 0.0d0
         integer :: chain = 0
! how the rate is computed
         integer :: rate_law = shares_rate
         integer :: burns_in = hydrogen_burning
         real(dp) :: fit(8) = 0.0d0
         integer :: scale_slot = 0
         type(screening_t) :: screening = screening_t()
      end type reaction_t

! rows of the reaction table
      integer, parameter :: n_rows = 13
      type(reaction_t), parameter :: reactions(n_rows) = [ &
!                     name                  rate slot  branch                reactants                 changes (species, nu)
           reaction_t('pp chain to He3',     r_pp,      unbranched,           [s_h1, s_h1, 0],          &
                      [s_h1, -3, s_he3, 1, 0, 0, 0, 0, 0, 0], &
                      q_mev=6.664d0, chain=chain_pp1, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, scale_slot=r_pp, &
                      fit=[0.12317d0, 1.08749d0, .93833d0, 0.d0, 0.d0, -3.3804d0, 20.8964d0, 0.d0], &
                      screening=screening_t(1.d0, 1.630d0, 1.175d0, 0.52d0, -0.413d0)), &
           reaction_t('He3(He3,2p)He4',      r_he3he3,  unbranched,           [s_he3, s_he3, 0],        &
                      [s_h1, 2, s_he3, -2, s_he4, 1, 0, 0, 0, 0], &
                      q_mev=12.860d0, chain=chain_pp1, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, scale_slot=r_he3he3, &
                      fit=[.03392d0, -.273d0, -.0648d0, 0.d0, 0.d0, -12.2757d0, 76.6003d0, 0.d0], &
                      screening=screening_t(4.d0, 5.917d0, 3.73d0, 1.31d0, -0.655d0)), &
           reaction_t('He3(He4,g) Be7(e) ppII', r_he3he4, be7_electron_capture, [s_he3, s_he4, 0],      &
                      [s_h1, -1, s_he3, -1, s_he4, 1, 0, 0, 0, 0], &
                      q_mev=1.586d0+17.394d0, chain=chain_pp2, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, scale_slot=r_he3he4, &
                      fit=[.0325d0, -.2085d0, -.0474d0, 0.d0, 0.d0, -12.826d0, 67.8036d0, 0.d0], &
                      screening=screening_t(4.d0, 5.917d0, 3.73d0, 1.31d0, -0.655d0)), &
           reaction_t('He3(He4,g) Be7(p) ppIII', r_he3he4, be7_proton_capture, [s_he3, s_he4, 0],      &
                      [s_h1, -1, s_he3, -1, s_he4, 1, 0, 0, 0, 0], &
                      q_mev=1.586d0+11.499d0, chain=chain_pp3), &
           reaction_t('C12(p,g) to C13',     r_pc12,    unbranched,           [s_h1, s_c12, 0],         &
                      [s_h1, -1, s_c12, -1, s_c13, 1, 0, 0, 0, 0], &
                      q_mev=3.457372d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, scale_slot=r_pc12, &
                      fit=[.0304d0, .7630d0, .1626d0, 4.79d0, 2.595d0, -13.6899d0, 69.130d0, 0.0d0], &
                      screening=screening_t(6.d0, 8.302d0, 4.804d0, 1.488d0, -0.643d0)), &
           reaction_t('C13(p,g)N14',         r_pc13,    unbranched,           [s_h1, s_c13, 0],         &
                      [s_h1, -1, s_c13, -1, s_n14, 1, 0, 0, 0, 0], &
                      q_mev=7.550628d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, scale_slot=r_pc13, &
                      fit=[.03035d0, -0.4044d0, -.08598d0, 7.456d0, 4.032d0, -13.7173d0, 70.3809d0, 0.0d0], &
                      screening=screening_t(6.d0, 8.302d0, 4.804d0, 1.488d0, -0.643d0)), &
           reaction_t('N14(p,g) N15(p,a)C12', r_pn14,   n15_to_c12,           [s_h1, s_n14, 0],         &
                      [s_h1, -2, s_n14, -1, s_c12, 1, s_he4, 1, 0, 0], &
                      q_mev=9.054d0+4.966d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, scale_slot=r_pn14, &
                      fit=[.0273d0, -1.60d0, -.3064d0, 0.0d0, 0.0d0, -15.2281d0, 69.8517d0, 0.0d0], &
                      screening=screening_t(7.d0, 9.520d0, 5.385d0, 1.61d0, -0.659d0)), &
           reaction_t('N14(p,g) N15(p,g)O16', r_pn14,   n15_to_o16,           [s_h1, s_n14, 0],         &
                      [s_h1, -2, s_n14, -1, s_o16, 1, 0, 0, 0, 0], &
                      q_mev=9.054d0+12.128d0, chain=chain_cno), &
           reaction_t('O16(p,g) O17(p,a)N14', r_po16,   unbranched,           [s_h1, s_o16, 0],         &
                      [s_h1, -2, s_o16, -1, s_n14, 1, s_he4, 1, 0, 0], &
                      q_mev=3.553d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, scale_slot=r_po16, &
                      fit=[.02494d0, -1.224d0, -.2139d0, .69703d0, 0.3097d0, -16.6925d0, 70.8012d0, 0.d0], &
                      screening=screening_t(8.d0, 10.716d0, 5.941d0, 1.721d0, -0.673d0)), &
           reaction_t('C12(a,g)O16',         r_ac12,    unbranched,           [s_he4, s_c12, 0],        &
                      [s_he4, -1, s_c12, -1, s_o16, 1, 0, 0, 0, 0], &
                      q_mev=7.162d0, chain=chain_alpha_capture, &
                      rate_law=cf88_c12_alpha_gamma, burns_in=helium_burning, &
                      screening=screening_t(12.d0, 16.192d0, 9.014d0, 2.577d0, -0.889d0)), &
           reaction_t('C13(a,n)O16',         r_ac13,    unbranched,           [s_he4, s_c13, 0],        &
                      [s_he4, -1, s_c13, -1, s_o16, 1, 0, 0, 0, 0], emitted_mass=1, &
                      q_mev=2.216d0, chain=chain_alpha_capture, &
                      rate_law=cf88_c13_alpha_n, burns_in=helium_burning, &
                      screening=screening_t(12.d0, 16.192d0, 9.014d0, 2.577d0, -0.889d0)), &
           reaction_t('N14(a,g) to O18',     r_an14,    unbranched,           [s_he4, s_n14, 0],        &
                      [s_he4, -1, s_n14, -1, s_o18, 1, 0, 0, 0, 0], &
                      q_mev=5.815d0, chain=chain_alpha_capture, &
                      rate_law=cf88_n14_alpha_gamma, burns_in=helium_burning, &
                      screening=screening_t(14.d0, 18.606d0, 10.15d0, 2.81d0, -0.92d0)), &
           reaction_t('triple alpha',        r_3alpha,  unbranched,           [s_he4, s_he4, s_he4],    &
                      [s_he4, -3, s_c12, 1, 0, 0, 0, 0, 0, 0], &
                      q_mev=7.275d0, chain=chain_triple_alpha, &
                      rate_law=cf88_triple_alpha, burns_in=helium_burning, &
                      screening=screening_t(12.d0, 16.192d0, 9.104d0, 2.577d0, -0.889d0)) ]

      public :: burn_implicitly, network_table_errors, write_network_table, &
           compute_equilibrium_burn_rates, compute_reaction_energies, &
           compute_reaction_rates, compute_network_rates, reaction_counts

contains

!----------------------------------------------------------------------
! compute_reaction_rates: the H/He reaction rates of one zone and their
! logarithmic derivatives with respect to T and rho -- the one rate
! evaluation behind both the composition rate vector
! (compute_network_rates) and the energy generation (nuclear_energy).
! Every row that owns a rate gets it from reaction_rates'
! compute_screened_rate when the row burns: hydrogen-burning rows need
! hydrogen, helium-burning rows log T >= logT_cutoffs(4). nr%rate(r_*)
! is per gram per second with the abundance factors omitted, screened
! and scaled by net_settings; each consumer applies its own floor. With
! hydrogen, also the Be7 and N15 branchings; always the hep factor.
! Below logT_cutoffs(1) everything is zero. mass_frac(iso_*) is the
! composition the screening sees (screening_mass_fractions).
pure subroutine compute_reaction_rates(log_density, log_temperature, &
     mass_frac, nr)
      real(dp), intent(in) :: log_density, log_temperature, &
           mass_frac(num_isotopes)
      type(nuclear_rates), intent(out) :: nr
      type(zone_conditions) :: zone
      logical :: hydrogen, helium_burning_temperature
      integer :: j, k

      call compute_zone_conditions(log_density, log_temperature, mass_frac, &
           zone, nr)
      if (log_temperature <= net_settings%logT_cutoffs(1)) return
      hydrogen = mass_frac(iso_h1) /= 0.0d0
      helium_burning_temperature = &
           log_temperature >= net_settings%logT_cutoffs(4)
      do j = 1, n_rows
         if (reactions(j)%rate_law == shares_rate) cycle
         if (reactions(j)%burns_in == hydrogen_burning .and. .not. hydrogen) cycle
         if (reactions(j)%burns_in == helium_burning .and. &
              .not. helium_burning_temperature) cycle
         k = reactions(j)%rate_row
         call compute_screened_rate(reactions(j)%rate_law, reactions(j)%fit, &
              reactions(j)%scale_slot, reactions(j)%screening, zone, &
              nr%rate(k), nr%dlnrate_dlnrho(k), nr%dlnrate_dlnt(k))
      end do
      if (hydrogen) then
         call compute_be7_branching(zone, mass_frac(iso_h1), nr)
         call compute_n15_branching(zone, nr)
      end if
      call compute_hep_rate(zone, nr)
end subroutine compute_reaction_rates

!----------------------------------------------------------------------
! compute_network_rates: the composition rate vector of one zone
! (rate_vec, num_rate_rows): the reaction rates at their r_* slots,
! zeroed at or below 1e-5 and turned into per second per amu (times
! amu_cgs), then the N15(p,alpha) and Be7 electron-capture fractions at
! rv_frac_n15_pa and rv_frac_be7_ecap. Slots without a table row are
! zero. xa is the zone's composition in the star%xa layout.
pure subroutine compute_network_rates(log_density, log_temperature, xa, &
     rate_vec)
      use physical_constants, only: amu_cgs
      real(dp), intent(in) :: log_density, log_temperature, xa(:)
      real(dp), intent(out) :: rate_vec(num_rate_rows)
      type(nuclear_rates) :: nr
      integer :: i

      call compute_reaction_rates(log_density, log_temperature, &
           screening_mass_fractions(xa), nr)
      do i = 1, num_reactions
         if (nr%rate(i) <= 1.d-5) nr%rate(i) = 0.0d0
      end do
! per s per amu (amu_cgs)
      rate_vec(1:num_reactions) = nr%rate*amu_cgs
      rate_vec(rv_frac_n15_pa) = nr%frac_n15_pa
      rate_vec(rv_frac_be7_ecap) = nr%frac_be7_ecap
end subroutine compute_network_rates

!----------------------------------------------------------------------
! reaction_counts: each rate slot's reactions per gram per second, its
! rate (abundance factors omitted) times its reactants' mass fractions,
! from the row that owns the rate; slots without a row are zero. x holds
! the network species (s_* order). For the solar neutrino fluxes.
pure function reaction_counts(rate, x) result(counts)
      real(dp), intent(in) :: rate(num_reactions), x(n_network)
      real(dp) :: counts(num_reactions)
      integer :: j, k, m

      counts = 0.0d0
      do j = 1, n_rows
         if (reactions(j)%rate_law == shares_rate) cycle
         k = reactions(j)%rate_row
         counts(k) = rate(k)
         do m = 1, max_reactants
            if (reactions(j)%reactants(m) /= 0) &
                 counts(k) = counts(k)*x(reactions(j)%reactants(m))
         end do
      end do
end function reaction_counts

! ---------------------------------------------------------------
! Self-check of the table: the number of rows that do not conserve
! baryon number, or whose reactants include a species outside the
! implicit solve (an explicit product cannot feed back).
pure integer function network_table_errors()
      integer :: j

      network_table_errors = 0
      do j = 1, n_rows
         if (baryon_imbalance(j) /= 0 .or. &
              any(reactions(j)%reactants > n_solved) .or. &
              neutrino_energy(j) < -1.0d-3) &
              network_table_errors = network_table_errors + 1
      end do
end function network_table_errors

! The energy row j releases, MeV, from the atomic mass excesses of what
! it consumes and makes (neutrinos included; positron annihilation
! included, through the atomic masses).
pure real(dp) function q_from_masses(j)
      integer, intent(in) :: j
      integer :: c

      q_from_masses = -reactions(j)%emitted_mass*neutron_mass_excess_mev
      do c = 1, 2*max_change - 1, 2
         if (reactions(j)%change(c) /= 0) q_from_masses = q_from_masses - &
              reactions(j)%change(c+1)*mass_excess_mev(reactions(j)%change(c))
      end do
end function q_from_masses

! The mean neutrino energy per reaction the table implies, MeV: zero,
! to the 1 keV the table's Q values are rounded to, for a row with no
! beta+ decay or electron capture.
pure real(dp) function neutrino_energy(j)
      integer, intent(in) :: j
      neutrino_energy = q_from_masses(j) - reactions(j)%q_mev
end function neutrino_energy

pure integer function baryon_imbalance(j)
      integer, intent(in) :: j
      integer :: c

      baryon_imbalance = reactions(j)%emitted_mass
      do c = 1, 2*max_change - 1, 2
         if (reactions(j)%change(c) /= 0) baryon_imbalance = baryon_imbalance + &
              nint(mass_number(reactions(j)%change(c)))*reactions(j)%change(c+1)
      end do
end function baryon_imbalance

! ---------------------------------------------------------------
! The table, one row per reaction, for tests and documentation.
subroutine write_network_table(unit)
      integer, intent(in) :: unit
      integer :: j

      write(unit, '(a)') '# reaction                 rate row  branch  baryon imbalance     Q (MeV)  chain' // &
           '  Q masses (MeV)  neutrinos (MeV)'
      do j = 1, n_rows
         write(unit, '(a24, 2x, i8, 2x, i6, 2x, i16, 2x, f10.6, 2x, i5, 2x, f14.6, 2x, f15.6)') &
              reactions(j)%name, reactions(j)%rate_row, reactions(j)%branch, &
              baryon_imbalance(j), reactions(j)%q_mev, reactions(j)%chain, &
              q_from_masses(j), neutrino_energy(j)
      end do
      write(unit, '(a, i0)') '# table errors: ', network_table_errors()
end subroutine write_network_table

! ---------------------------------------------------------------
! Rate of reaction j at abundances x: R_j = b_j g_j prod X(law).
pure real(dp) function reaction_rate(j, rate_vec, x)
      integer, intent(in) :: j
      real(dp), intent(in) :: rate_vec(:), x(n_network)
      integer :: m

      reaction_rate = branch_factor(j, rate_vec)*rate_vec(reactions(j)%rate_row)
      do m = 1, max_reactants
         if (reactions(j)%reactants(m) /= 0) &
              reaction_rate = reaction_rate*x(reactions(j)%reactants(m))
      end do
end function reaction_rate

! ---------------------------------------------------------------
! dR_j/dX_k: the product rule over the factors of the rate law.
pure real(dp) function reaction_rate_derivative(j, k, rate_vec, x)
      integer, intent(in) :: j, k
      real(dp), intent(in) :: rate_vec(:), x(n_network)
      integer :: m, m2
      real(dp) :: term

      reaction_rate_derivative = 0.0d0
      do m = 1, max_reactants
         if (reactions(j)%reactants(m) /= k) cycle
         term = branch_factor(j, rate_vec)*rate_vec(reactions(j)%rate_row)
         do m2 = 1, max_reactants
            if (m2 /= m .and. reactions(j)%reactants(m2) /= 0) &
                 term = term*x(reactions(j)%reactants(m2))
         end do
         reaction_rate_derivative = reaction_rate_derivative + term
      end do
end function reaction_rate_derivative

pure real(dp) function branch_factor(j, rate_vec)
      integer, intent(in) :: j
      real(dp), intent(in) :: rate_vec(:)

      branch_factor = branch_fraction(j, rate_vec(rv_frac_n15_pa), &
           rate_vec(rv_frac_be7_ecap))
end function branch_factor

! b_j from the two branching fractions
pure real(dp) function branch_fraction(j, frac_n15_pa, frac_be7_ecap)
      integer, intent(in) :: j
      real(dp), intent(in) :: frac_n15_pa, frac_be7_ecap

      select case (reactions(j)%branch)
      case (n15_to_c12)
         branch_fraction = frac_n15_pa
      case (n15_to_o16)
         branch_fraction = 1.0d0 - frac_n15_pa
      case (be7_electron_capture)
         branch_fraction = frac_be7_ecap
      case (be7_proton_capture)
         branch_fraction = 1.0d0 - frac_be7_ecap
      case default
         branch_fraction = 1.0d0
      end select
end function branch_fraction

! ---------------------------------------------------------------
! The energy generation of the table's reactions (erg/g/s), by rate
! slot (energy_by_slot, r_*) and by chain (energy_by_chain, chain_*):
! the sum over rows of b_j * rate(slot) * prod X(law_j) * q_j, with
! rate the per-gram rates of compute_reaction_rates (abundance factors
! omitted) and x the network species' mass fractions.
pure subroutine compute_reaction_energies(rate, x, frac_n15_pa, frac_be7_ecap, &
     energy_by_slot, energy_by_chain)
      use reaction_rates, only: num_reactions
      real(dp), intent(in) :: rate(num_reactions), x(n_network), &
           frac_n15_pa, frac_be7_ecap
      real(dp), intent(out) :: energy_by_slot(num_reactions), &
           energy_by_chain(n_chains)
! MeV to erg
      real(dp), parameter :: erg_per_mev = 1.602177d-6
      real(dp) :: e
      integer :: j, m

      energy_by_slot = 0.0d0
      energy_by_chain = 0.0d0
      do j = 1, n_rows
         e = rate(reactions(j)%rate_row)
         do m = 1, max_reactants
            if (reactions(j)%reactants(m) /= 0) e = e*x(reactions(j)%reactants(m))
         end do
         e = e*branch_fraction(j, frac_n15_pa, frac_be7_ecap)* &
              reactions(j)%q_mev*erg_per_mev
         energy_by_slot(reactions(j)%rate_row) = &
              energy_by_slot(reactions(j)%rate_row) + e
         energy_by_chain(reactions(j)%chain) = &
              energy_by_chain(reactions(j)%chain) + e
      end do
end subroutine compute_reaction_energies

! ---------------------------------------------------------------
! dX_i/dt / A_i = sum_j nu_ij R_j for every network species.
pure subroutine compute_abundance_rates(rate_vec, x, dxdt_per_a)
      real(dp), intent(in) :: rate_vec(:), x(n_network)
      real(dp), intent(out) :: dxdt_per_a(n_network)
      integer :: j, c
      real(dp) :: r

      dxdt_per_a = 0.0d0
      do j = 1, n_rows
         r = reaction_rate(j, rate_vec, x)
         do c = 1, 2*max_change - 1, 2
            if (reactions(j)%change(c) == 0) cycle
            dxdt_per_a(reactions(j)%change(c)) = dxdt_per_a(reactions(j)%change(c)) + &
                 reactions(j)%change(c+1)*r
         end do
      end do
end subroutine compute_abundance_rates

! ---------------------------------------------------------------
! One implicit (backward-Euler) burning step of length dt (s).
!   x_start   abundances at the start of the step (all network species)
!   rate_vec  the rate vector (slots r_*, rv_*) for this zone or zone average
!   x         in: the first guess for the solved species; out: the
!             solution, with the explicit products (O18) updated
! Newton converges when every |correction| < absolute_tolerance and every
! |correction/X| < relative_tolerance (species below min_abundance are
! set to zero and skip the relative test); ierr = 1 after
! max_iterations attempts. As in kemcom, Jacobian entries below 1e-8 in
! magnitude (off the diagonal) are dropped, and if none remain the
! residual itself is taken as the correction.
pure subroutine burn_implicitly(x_start, rate_vec, dt, min_abundance, &
     absolute_tolerance, relative_tolerance, max_iterations, x, ierr)
      use numerics_lib, only: simeqc
      real(dp), intent(in) :: x_start(n_network), rate_vec(:), dt
      real(dp), intent(in) :: min_abundance, absolute_tolerance, &
           relative_tolerance
      integer, intent(in) :: max_iterations
      real(dp), intent(inout) :: x(n_network)
      integer, intent(out) :: ierr

! system(:, 1:n) is the Jacobian, system(:, n+1) the residual; simeqc
! leaves the Newton correction in system(:, n+1)
      real(dp) :: system(n_solved, n_solved+1), dt_a(n_network)
      real(dp) :: dxdt_per_a(n_network), correction(n_solved)
      real(dp) :: max_abs_change, max_rel_change
      integer :: i, k, j, c, iterations
      logical :: off_diagonal

      ierr = 0
      dt_a = mass_number*dt
      iterations = 0
      do
         call compute_abundance_rates(rate_vec, x, dxdt_per_a)
         do i = 1, n_solved
            system(i, n_solved+1) = x(i) - x_start(i) - dt_a(i)*dxdt_per_a(i)
         end do
! J_ik = delta_ik - dt A_i sum_j nu_ij dR_j/dX_k
         system(:, 1:n_solved) = 0.0d0
         do j = 1, n_rows
            do c = 1, 2*max_change - 1, 2
               i = reactions(j)%change(c)
               if (i == 0 .or. i > n_solved) cycle
               do k = 1, n_solved
                  system(i, k) = system(i, k) - dt_a(i)*reactions(j)%change(c+1)* &
                       reaction_rate_derivative(j, k, rate_vec, x)
               end do
            end do
         end do
         off_diagonal = .false.
         do k = 1, n_solved
            system(k, k) = 1.0d0 + system(k, k)
            do i = 1, n_solved
               if (i == k) cycle
               if (abs(system(i, k)) < 1.0d-8) then
                  system(i, k) = 0.0d0
               else
                  off_diagonal = .true.
               end if
            end do
         end do
         if (off_diagonal) then
            call simeqc(system, n_solved+1, n_solved, ierr)
            if (ierr /= 0) return
         end if
         correction = system(:, n_solved+1)

         max_abs_change = 0.0d0
         max_rel_change = 0.0d0
         do i = 1, n_solved
            x(i) = x(i) - correction(i)
            max_abs_change = max(max_abs_change, abs(correction(i)))
            if (x(i) < min_abundance) then
               x(i) = 0.0d0
            else
               max_rel_change = max(max_rel_change, abs(correction(i)/x(i)))
            end if
         end do
         if (max_abs_change < absolute_tolerance .and. &
              max_rel_change < relative_tolerance) exit
         iterations = iterations + 1
         if (iterations >= max_iterations) then
            ierr = 1
            return
         end if
      end do

! explicit products, from the converged abundances
      call compute_abundance_rates(rate_vec, x, dxdt_per_a)
      do i = n_solved+1, n_network
         x(i) = x_start(i) + dt_a(i)*dxdt_per_a(i)
      end do
end subroutine burn_implicitly

!----------------------------------------------------------------------
! compute_equilibrium_burn_rates (formerly eqburn)
!----------------------------------------------------------------------
! The explicit hydrogen- and helium-burning rates of one radiative shell
! (zone_begin = zone_end) or one mixed region, with He3 and the CN cycle
! in equilibrium: the classic approximation that follows only X, Y, C12
! and O16. burn_equilibrium (core/burn_step) advances X and, under
! helium burning, Z, C12 and O16 with them at Newton level 1; the
! hydrogen-shell timestep limiter reads dX/dt.
!
! A mixed region takes the mass-weighted averages of the start-of-step
! abundances and of the rate vectors. With R the rate-vector rows
! (per s per amu, the 1/A factors included) and X, Y the abundances:
!   He3: the root of 2 R33 X3^2 + R34 Y X3 = R11 X^2, or the current X3
!     if smaller (in outer layers He3 takes longer than a step to reach
!     equilibrium);
!   CN: N14 from C12/12 + C13/13 + N14/14 conserved, with
!     C12 = N14 R14/R12 and C13 = N14 R14/R13;
!   dX/dt = -3 R11 X^2 + 2 R33 X3^2 - R34 Y X3
!           - R12 X C12 - R13 X C13 - 2 R14 X N14 - 2 R16 X O16;
!   dY/dt = -4 (Rag Y C12 + 3 R3a Y^3), dC/dt = 12 (R3a Y^3 - Rag Y C12),
!   dO/dt = 16 Rag Y C12.
! The pp chains burn when X > 1e-10 and log T > logT_cutoffs(2), the CN
! cycle above logT_cutoffs(3), helium when Y > 1e-10 and log T >
! logT_cutoffs(4); every rate is zero at or below logT_cutoffs(1), and
! the region's temperature is that of its first zone. Also returned:
! the region's start-of-step X and Z, its C12 (the CN-equilibrium value
! when the cycle runs) and its O16.
pure subroutine compute_equilibrium_burn_rates(rate_by_zone, xa_start, &
     shell_mass, shell_temperature, zone_begin, zone_end, dc_dt, do_dt, &
     dx_dt, dy_dt, equilibrium_xc12, equilibrium_xo16, hydrogen_fraction, &
     metal_fraction)
      use composition_layout, only: n_species_basic, i_h1, i_he4, i_metals, &
           i_he3, i_c12, i_c13, i_n14, i_o16
      use reaction_rates, only: net_settings
! rate_by_zone(r_*/rv_*, zone) from compute_network_rates; xa_start the
! start-of-step composition (star%xa layout); shell_temperature log10 T
      real(dp), intent(in) :: rate_by_zone(:,:), xa_start(:,:), &
           shell_mass(:), shell_temperature(:)
      integer, intent(in) :: zone_begin, zone_end
      real(dp), intent(out) :: dc_dt, do_dt, dx_dt, dy_dt, &
           equilibrium_xc12, equilibrium_xo16, hydrogen_fraction, &
           metal_fraction

      real(dp) :: x(n_species_basic), rate(size(rate_by_zone,1)), region_mass
      real(dp) :: helium_fraction, he3, c13, n14, log_temperature
      real(dp) :: pp_term, he3he4_term, he3_equilibrium
      real(dp) :: cno_sum, c12_per_n14, c13_per_n14
      real(dp) :: triple_alpha_term, c12_alpha_term
      integer :: k

! the region's start-of-step abundances and rates
      if (zone_begin == zone_end) then
         x = xa_start(1:n_species_basic, zone_begin)
         rate = rate_by_zone(:, zone_begin)
      else
         region_mass = 0.0d0
         x = 0.0d0
         rate = 0.0d0
         do k = zone_begin, zone_end
            region_mass = region_mass + shell_mass(k)
            x = x + xa_start(1:n_species_basic, k)*shell_mass(k)
            rate = rate + shell_mass(k)*rate_by_zone(:, k)
         end do
         x = x/region_mass
         rate = rate/region_mass
      end if
      hydrogen_fraction = x(i_h1)
      helium_fraction = x(i_he4)
      metal_fraction = x(i_metals)
      he3 = x(i_he3)
      equilibrium_xc12 = x(i_c12)
      c13 = x(i_c13)
      n14 = x(i_n14)
      equilibrium_xo16 = x(i_o16)

      dx_dt = 0.0d0
      dy_dt = 0.0d0
      dc_dt = 0.0d0
      do_dt = 0.0d0
      log_temperature = shell_temperature(zone_begin)
      if (log_temperature < net_settings%logT_cutoffs(1)) return

! hydrogen burning: the pp chains, with He3 in equilibrium
      if (hydrogen_fraction > 1.0d-10 .and. &
           log_temperature > net_settings%logT_cutoffs(2)) then
         pp_term = rate(r_pp)*hydrogen_fraction**2
         he3he4_term = rate(r_he3he4)*helium_fraction
         he3_equilibrium = 0.25d0*(sqrt(he3he4_term**2 + &
              8.0d0*rate(r_he3he3)*pp_term) - he3he4_term)/rate(r_he3he3)
         he3 = min(he3, he3_equilibrium)
         dx_dt = -3.0d0*pp_term + 2.0d0*rate(r_he3he3)*he3**2 &
              - he3he4_term*he3
! the CN cycle, with C12, C13 and N14 in equilibrium
         if (log_temperature > net_settings%logT_cutoffs(3)) then
            cno_sum = equilibrium_xc12/1.2d1 + c13/1.3d1 + n14/1.4d1
            c12_per_n14 = rate(r_pn14)/rate(r_pc12)
            c13_per_n14 = rate(r_pn14)/rate(r_pc13)
            n14 = cno_sum/(c12_per_n14/1.2d1 + c13_per_n14/1.3d1 + 1.0d0/1.4d1)
            equilibrium_xc12 = n14*c12_per_n14
            c13 = n14*c13_per_n14
            dx_dt = dx_dt - rate(r_pc12)*hydrogen_fraction*equilibrium_xc12 - &
                 rate(r_pc13)*hydrogen_fraction*c13 - &
                 2.0d0*rate(r_pn14)*hydrogen_fraction*n14 - &
                 2.0d0*rate(r_po16)*hydrogen_fraction*equilibrium_xo16
         end if
      end if

! helium burning: triple alpha and C12(a,g)
      if (helium_fraction > 1.0d-10 .and. &
           log_temperature > net_settings%logT_cutoffs(4)) then
         triple_alpha_term = rate(r_3alpha)*helium_fraction**3
         c12_alpha_term = rate(r_ac12)*helium_fraction*equilibrium_xc12
         dy_dt = -4.0d0*(c12_alpha_term + 3.0d0*triple_alpha_term)
         dc_dt = 1.2d1*(triple_alpha_term - c12_alpha_term)
         do_dt = 1.6d1*c12_alpha_term
      end if
end subroutine compute_equilibrium_burn_rates

end module network
