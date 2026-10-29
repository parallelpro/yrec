!----------------------------------------------------------------------
! network
!----------------------------------------------------------------------
! The hydrogen- and helium-burning network as data, and the one implicit
! integrator that reads it (2026-10; it replaced the hand-written
! 7-species equations and Jacobian of kemcom).
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
! of the Caughlan & Fowler alpha-capture forms), its screening constants
! and whether it burns with hydrogen or in helium burning.
! compute_reaction_rates evaluates every row that owns a rate; a second
! outcome of the same reaction (shares_rate) uses the rate of the row
! above. The rates land at the rows' slots r_* (net_configuration),
! where the settings' scales of each reaction are too.
!
! Reactions. Reaction j proceeds at
!     R_j = b_j * g_j * prod_m X(reactants_j(m))    [per s per amu]
! where g_j is the composition rate vector at the row's slot
! (compute_network_rates: N_A<sigma v>, density, screening and the 1/A
! factors are already in it), X are mass fractions, and b_j is 1; or
! the N15(p,alpha) branching fraction f (or 1-f), carried at
! rv_frac_n15_pa, for the two outcomes of N14(p,g); or the Be7
! electron-capture fraction f_e (or 1-f_e), at rv_frac_be7_ecap, for the
! two of He3(He4,g). Each
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
! N14(a,g)F18(e+ nu)O18 deposits 5.646 MeV: Q = 6.071 MeV from the masses
! less the mean neutrino energy of the F18 decay, 0.425 MeV (96.73% beta+
! with <E_nu> = 633.5 - 249.8 keV, 3.27% electron capture at 1655 keV;
! ENSDF, Tilley et al. 1995). (2026-10: it was 5.815, implying 0.256 MeV,
! about the mean positron energy.)
! He3(He4,g)Be7 appears twice, once per Be7 branch (electron capture:
! pp II; proton capture: pp III): the same change in abundances, with
! the two branch fractions and their own energies.
!
! Implicit step. Over dt the solve finds X with
!     F_i(X) = X_i - X0_i - dt A_i sum_j nu_ij R_j(X) = 0
! by Newton iteration on J_ik = dF_i/dX_k.
!
! The module is the two tables and the code that reads them; the
! tables' self-check and set_net_settings, which runs it, are
! network_checks'. The table's columns are public for that module (code
! outside net/ sees only what net_lib re-exports).
!
! Adding a reaction means adding a row to the table, and naming its rate
! slot r_* in net_configuration (raising num_reactions); a rate of a new
! functional form also needs its rate law in reaction_rates. Adding a
! species means adding it to the species table, before s_o18 if it is
! solved, and to the rows that make or burn it; net/README.md's recipe
! lists the files outside net/ that follow.
module network
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: i_h1, i_h2, i_he3, i_he4, i_c12, i_c13, &
           i_n14, i_o16, i_o18
      use net_configuration, only: r_pp, r_he3he3, r_he3he4, &
           r_pc12, r_pc13, r_pn14, r_po16, r_ac13, &
           r_ac12, r_an14, r_3alpha, rv_frac_n15_pa, rv_frac_be7_ecap, &
           num_reactions, num_rate_rows, net_settings
      use screening, only: num_isotopes, screening_t, iso_h1, iso_h2, &
           iso_he3, iso_he4, iso_c12, iso_c13, iso_n14, iso_o16, iso_o18
      use reaction_rates, only: nuclear_rates, shares_rate, &
           s_factor_fit, cf88_c13_alpha_n, cf88_c12_alpha_gamma, &
           cf88_n14_alpha_gamma, cf88_triple_alpha, &
           compute_zone_conditions, compute_screened_rate, &
           compute_be7_branching, n15_alpha_fraction, erg_per_mev
      implicit none
      private

! network species: positions in the abundance vector of the solve
      integer, parameter, public :: s_h1 = 1, s_he3 = 2, s_he4 = 3, &
           s_c12 = 4, s_c13 = 5, s_n14 = 6, s_o16 = 7, s_o18 = 8
! the first n_solved species are solved implicitly; the rest are
! products updated explicitly after convergence
      integer, parameter, public :: n_solved = 7, n_network = 8
! The species table, in s_* order: each species' name, its slot in the
! composition vector xa, its isotope in the screening composition
! (screening's iso_*), its mass number (the integer A the rate
! convention uses) and its atomic mass excess (MeV, AME2020: Wang et al.
! 2021), which checks the rows' energies (network_checks).
      type :: species_t
         character(len=4) :: name
         integer :: xa_slot, screening_isotope
         real(dp) :: mass_number, mass_excess_mev
      end type species_t
      type(species_t), parameter :: species(n_network) = [ &
           species_t('H1  ', i_h1,  iso_h1,  1.0d0,  7.288971d0), &
           species_t('He3 ', i_he3, iso_he3, 3.0d0, 14.931219d0), &
           species_t('He4 ', i_he4, iso_he4, 4.0d0,  2.424916d0), &
           species_t('C12 ', i_c12, iso_c12, 12.0d0, 0.0d0), &
           species_t('C13 ', i_c13, iso_c13, 13.0d0, 3.125009d0), &
           species_t('N14 ', i_n14, iso_n14, 14.0d0, 2.863417d0), &
           species_t('O16 ', i_o16, iso_o16, 16.0d0, -4.737001d0), &
           species_t('O18 ', i_o18, iso_o18, 18.0d0, -0.782816d0)]
! the columns of the table the rest of the code reads
      integer, parameter, public :: network_species_slot(n_network) = species%xa_slot
      integer, parameter :: screening_isotope(n_network) = species%screening_isotope
      real(dp), parameter, public :: mass_number(n_network) = species%mass_number
      real(dp), parameter, public :: mass_excess_mev(n_network) = species%mass_excess_mev
! the atomic mass excess of the neutron, the one untracked particle a
! row emits, MeV
      real(dp), parameter, public :: neutron_mass_excess_mev = 8.071318d0
! the rounding of the rows' Q values, MeV: a row whose Q exceeds the
! energy its masses release by more is wrong (network_checks)
      real(dp), parameter, public :: q_rounding_mev = 1.0d-3

! branch factor b_j: 1, the N15(p,alpha) fraction f or 1 - f, or the Be7
! electron-capture fraction f_e or 1 - f_e
      integer, parameter, public :: unbranched = 0, n15_to_c12 = 1, n15_to_o16 = 2, &
           be7_electron_capture = 3, be7_proton_capture = 4

! burn_implicitly's ierr: 0, or the Newton iteration ran out of
! iterations, or its linear system was singular (simeqc)
      integer, parameter, public :: burn_not_converged = 1, &
           burn_singular_jacobian = 2

! when a reaction burns: with hydrogen present, or at helium-burning
! temperatures (log T >= logT_cutoffs(4))
      integer, parameter, public :: hydrogen_burning = 1, helium_burning = 2

! the chain a reaction's energy is counted under, and its name
! (nuclear_energy_result%eps_by_chain is indexed by these)
      integer, parameter, public :: chain_pp1 = 1, chain_pp2 = 2, &
           chain_pp3 = 3, chain_cno = 4, chain_alpha_capture = 5, &
           chain_triple_alpha = 6, n_chains = 6
      character(len=*), parameter, public :: chain_name(n_chains) = [ &
           'pp I          ', 'pp II         ', 'pp III        ', &
           'CNO           ', 'alpha capture ', 'triple alpha  ']

! One reaction. reactants lists the species whose mass fractions multiply
! the rate (repeat a species for a square or cube; 0 = unused). change lists
! (species, nu) pairs: the change in per-amu abundance per reaction
! (0 = unused). emitted_mass is the mass number carried off by particles
! the network does not track (the neutron of C13(a,n)O16), so that every
! row conserves baryon number exactly: sum_i A_i nu_i + emitted_mass = 0.
! q_mev is the energy released (from the atomic masses, so positron
! annihilation included) less the mean energy of the row's neutrinos:
! never more than the mass-energy released. network_checks checks both;
! net/test/test_net pins the table. The rate: rate_law (a rate law
! of reaction_rates, or shares_rate), burns_in (hydrogen_burning or
! helium_burning), fit (the S-factor law's q1..q8; its multipliers are
! the settings' scales at the row's slot), s0_kevb (the S0 in keV b the
! fit is normalized to, Bahcall & Ulrich 1988: the reference the input
! S0 is scaled against, reference_s0) and screening (the charges of its
! two nuclei; triple alpha takes C12 + alpha's, as engeb did).
      integer, parameter, public :: max_reactants = 3, max_change = 5
! the composition path drops a rate (per gram per s, abundance factors
! omitted) at or below min_composition_rate (compute_network_rates); the
! implicit burn drops off-diagonal Jacobian entries below
! min_jacobian_entry in magnitude, as kemcom did (burn_implicitly)
      real(dp), parameter :: min_composition_rate = 1.0d-5, &
           min_jacobian_entry = 1.0d-8
      type, public :: reaction_t
         character(len=24) :: name
         integer :: slot
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
         real(dp) :: s0_kevb = 0.0d0
         type(screening_t) :: screening = screening_t()
      end type reaction_t

! rows of the reaction table
      integer, parameter, public :: n_rows = 13
      type(reaction_t), parameter, public :: reactions(n_rows) = [ &
!                     name                  rate slot  branch                reactants                 changes (species, nu)
           reaction_t('pp chain to He3',     r_pp,      unbranched,           [s_h1, s_h1, 0],          &
                      [s_h1, -3, s_he3, 1, 0, 0, 0, 0, 0, 0], &
                      q_mev=6.664d0, chain=chain_pp1, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, &
                      fit=[0.12317d0, 1.08749d0, .93833d0, 0.d0, 0.d0, -3.3804d0, 20.8964d0, 0.d0], s0_kevb=4.07d-22, &
                      screening=screening_t(1, 1)), &
           reaction_t('He3(He3,2p)He4',      r_he3he3,  unbranched,           [s_he3, s_he3, 0],        &
                      [s_h1, 2, s_he3, -2, s_he4, 1, 0, 0, 0, 0], &
                      q_mev=12.860d0, chain=chain_pp1, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, &
                      fit=[.03392d0, -.273d0, -.0648d0, 0.d0, 0.d0, -12.2757d0, 76.6003d0, 0.d0], s0_kevb=5.15d3, &
                      screening=screening_t(2, 2)), &
           reaction_t('He3(He4,g) Be7(e) ppII', r_he3he4, be7_electron_capture, [s_he3, s_he4, 0],      &
                      [s_h1, -1, s_he3, -1, s_he4, 1, 0, 0, 0, 0], &
                      q_mev=1.586d0+17.394d0, chain=chain_pp2, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, &
                      fit=[.0325d0, -.2085d0, -.0474d0, 0.d0, 0.d0, -12.826d0, 67.8036d0, 0.d0], s0_kevb=5.4d-1, &
                      screening=screening_t(2, 2)), &
           reaction_t('He3(He4,g) Be7(p) ppIII', r_he3he4, be7_proton_capture, [s_he3, s_he4, 0],      &
                      [s_h1, -1, s_he3, -1, s_he4, 1, 0, 0, 0, 0], &
                      q_mev=1.586d0+11.499d0, chain=chain_pp3), &
           reaction_t('C12(p,g) to C13',     r_pc12,    unbranched,           [s_h1, s_c12, 0],         &
                      [s_h1, -1, s_c12, -1, s_c13, 1, 0, 0, 0, 0], &
                      q_mev=3.457372d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, &
                      fit=[.0304d0, .7630d0, .1626d0, 4.79d0, 2.595d0, -13.6899d0, 69.130d0, 0.0d0], s0_kevb=1.45d0, &
                      screening=screening_t(1, 6)), &
           reaction_t('C13(p,g)N14',         r_pc13,    unbranched,           [s_h1, s_c13, 0],         &
                      [s_h1, -1, s_c13, -1, s_n14, 1, 0, 0, 0, 0], &
                      q_mev=7.550628d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, &
                      fit=[.03035d0, -0.4044d0, -.08598d0, 7.456d0, 4.032d0, -13.7173d0, 70.3809d0, 0.0d0], s0_kevb=5.5d0, &
                      screening=screening_t(1, 6)), &
           reaction_t('N14(p,g) N15(p,a)C12', r_pn14,   n15_to_c12,           [s_h1, s_n14, 0],         &
                      [s_h1, -2, s_n14, -1, s_c12, 1, s_he4, 1, 0, 0], &
                      q_mev=9.054d0+4.966d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, &
                      fit=[.0273d0, -1.60d0, -.3064d0, 0.0d0, 0.0d0, -15.2281d0, 69.8517d0, 0.0d0], s0_kevb=3.32d0, &
                      screening=screening_t(1, 7)), &
           reaction_t('N14(p,g) N15(p,g)O16', r_pn14,   n15_to_o16,           [s_h1, s_n14, 0],         &
                      [s_h1, -2, s_n14, -1, s_o16, 1, 0, 0, 0, 0], &
                      q_mev=9.054d0+12.128d0, chain=chain_cno), &
           reaction_t('O16(p,g) O17(p,a)N14', r_po16,   unbranched,           [s_h1, s_o16, 0],         &
                      [s_h1, -2, s_o16, -1, s_n14, 1, s_he4, 1, 0, 0], &
                      q_mev=3.553d0, chain=chain_cno, &
                      rate_law=s_factor_fit, burns_in=hydrogen_burning, &
                      fit=[.02494d0, -1.224d0, -.2139d0, .69703d0, 0.3097d0, -16.6925d0, 70.8012d0, 0.d0], s0_kevb=9.4d0, &
                      screening=screening_t(1, 8)), &
           reaction_t('C12(a,g)O16',         r_ac12,    unbranched,           [s_he4, s_c12, 0],        &
                      [s_he4, -1, s_c12, -1, s_o16, 1, 0, 0, 0, 0], &
                      q_mev=7.162d0, chain=chain_alpha_capture, &
                      rate_law=cf88_c12_alpha_gamma, burns_in=helium_burning, &
                      screening=screening_t(2, 6)), &
           reaction_t('C13(a,n)O16',         r_ac13,    unbranched,           [s_he4, s_c13, 0],        &
                      [s_he4, -1, s_c13, -1, s_o16, 1, 0, 0, 0, 0], emitted_mass=1, &
                      q_mev=2.216d0, chain=chain_alpha_capture, &
                      rate_law=cf88_c13_alpha_n, burns_in=helium_burning, &
                      screening=screening_t(2, 6)), &
           reaction_t('N14(a,g) to O18',     r_an14,    unbranched,           [s_he4, s_n14, 0],        &
                      [s_he4, -1, s_n14, -1, s_o18, 1, 0, 0, 0, 0], &
                      q_mev=5.646d0, chain=chain_alpha_capture, &
                      rate_law=cf88_n14_alpha_gamma, burns_in=helium_burning, &
                      screening=screening_t(2, 7)), &
           reaction_t('triple alpha',        r_3alpha,  unbranched,           [s_he4, s_he4, s_he4],    &
                      [s_he4, -3, s_c12, 1, 0, 0, 0, 0, 0, 0], &
                      q_mev=7.275d0, chain=chain_triple_alpha, &
                      rate_law=cf88_triple_alpha, burns_in=helium_burning, &
                      screening=screening_t(2, 6)) ]

      public :: burn_implicitly, compute_abundance_rates, &
           compute_reaction_energies, compute_reaction_rates, &
           compute_network_rates, reaction_counts, s_factor_slots, &
           reference_s0, branch_fraction, screening_mass_fractions

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
! hydrogen, also the Be7 and N15 branchings.
! Below logT_cutoffs(1) everything is zero. mass_frac(iso_*) is the
! composition the screening sees (screening_mass_fractions).
pure subroutine compute_reaction_rates(log_density, log_temperature, &
     mass_frac, nr)
      real(dp), intent(in) :: log_density, log_temperature, &
           mass_frac(num_isotopes)
      type(nuclear_rates), intent(out) :: nr
      logical :: hydrogen, helium_burning_temperature
      integer :: j, k

      if (log_temperature <= net_settings%logT_cutoffs(1)) return
      call compute_zone_conditions(log_density, log_temperature, mass_frac, &
           nr%zone)
      hydrogen = mass_frac(iso_h1) /= 0.0d0
      helium_burning_temperature = &
           log_temperature >= net_settings%logT_cutoffs(4)
      do j = 1, n_rows
         if (reactions(j)%rate_law == shares_rate) cycle
         if (reactions(j)%burns_in == hydrogen_burning .and. .not. hydrogen) cycle
         if (reactions(j)%burns_in == helium_burning .and. &
              .not. helium_burning_temperature) cycle
         k = reactions(j)%slot
         call compute_screened_rate(reactions(j)%rate_law, reactions(j)%fit, &
              k, reactions(j)%screening, nr%zone, &
              nr%rate(k), nr%dlnrate_dlnrho(k), nr%dlnrate_dlnt(k))
      end do
      if (hydrogen) then
         call compute_be7_branching(nr%zone, mass_frac(iso_h1), nr%be7)
         nr%frac_n15_pa = n15_alpha_fraction(nr%zone)
      end if
end subroutine compute_reaction_rates

!----------------------------------------------------------------------
! compute_network_rates: the composition rate vector of one zone
! (rate_vec, num_rate_rows): the reaction rates at their r_* slots,
! zeroed at or below min_composition_rate and turned into per second per amu (times
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
         if (nr%rate(i) <= min_composition_rate) nr%rate(i) = 0.0d0
      end do
! per s per amu (amu_cgs)
      rate_vec(1:num_reactions) = nr%rate*amu_cgs
      rate_vec(rv_frac_n15_pa) = nr%frac_n15_pa
      rate_vec(rv_frac_be7_ecap) = nr%be7%frac_ecap
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
         k = reactions(j)%slot
         counts(k) = rate(k)
         do m = 1, max_reactants
            if (reactions(j)%reactants(m) /= 0) &
                 counts(k) = counts(k)*x(reactions(j)%reactants(m))
         end do
      end do
end function reaction_counts

!----------------------------------------------------------------------
! screening_mass_fractions: the isotope mass fractions (screening's
! iso_* slots) that compute_reaction_rates screens with, from a zone's
! composition in the star%xa layout: each network species at its
! screening isotope (the species table), and deuterium, which burns
! outside the network but screens too (2026-10: the composition path used
! to leave it out). The isotopes nothing here names (n, H3, Ne20, Mg24)
! are zero.
pure function screening_mass_fractions(xa) result(mass_frac)
      real(dp), intent(in) :: xa(:)
      real(dp) :: mass_frac(num_isotopes)

      mass_frac = 0.0d0
      mass_frac(screening_isotope) = xa(network_species_slot)
      mass_frac(iso_h2) = xa(i_h2)
end function screening_mass_fractions

!----------------------------------------------------------------------
! reference_s0: the S0 (keV b) the fit of the row owning slot k is
! normalized to (its s0_kevb; 0 for a slot whose rate law takes none):
! map_user_inputs scales the input S0 of the reaction against it.
pure real(dp) function reference_s0(k)
      integer, intent(in) :: k
      integer :: j

      reference_s0 = 0.0d0
      do j = 1, n_rows
         if (reactions(j)%slot == k .and. reactions(j)%rate_law /= shares_rate) &
              reference_s0 = reactions(j)%s0_kevb
      end do
end function reference_s0

!----------------------------------------------------------------------
! s_factor_slots: which rate slots the S-factor law computes (the rows
! with rate_law s_factor_fit), the rates the energy path floors
! (nuclear_energy's energy_path_rates).
pure function s_factor_slots() result(mask)
      logical :: mask(num_reactions)
      integer :: j

      mask = .false.
      do j = 1, n_rows
         if (reactions(j)%rate_law == s_factor_fit) mask(reactions(j)%slot) = .true.
      end do
end function s_factor_slots

! ---------------------------------------------------------------
! Rate of reaction j at abundances x: R_j = b_j g_j prod X(law).
pure real(dp) function reaction_rate(j, rate_vec, x)
      integer, intent(in) :: j
      real(dp), intent(in) :: rate_vec(:), x(n_network)
      integer :: m

      reaction_rate = branch_factor(j, rate_vec)*rate_vec(reactions(j)%slot)
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
         term = branch_factor(j, rate_vec)*rate_vec(reactions(j)%slot)
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
      real(dp), intent(in) :: rate(num_reactions), x(n_network), &
           frac_n15_pa, frac_be7_ecap
      real(dp), intent(out) :: energy_by_slot(num_reactions), &
           energy_by_chain(n_chains)
      real(dp) :: e
      integer :: j, m

      energy_by_slot = 0.0d0
      energy_by_chain = 0.0d0
      do j = 1, n_rows
         e = rate(reactions(j)%slot)
         do m = 1, max_reactants
            if (reactions(j)%reactants(m) /= 0) e = e*x(reactions(j)%reactants(m))
         end do
         e = e*branch_fraction(j, frac_n15_pa, frac_be7_ecap)* &
              reactions(j)%q_mev*erg_per_mev
         energy_by_slot(reactions(j)%slot) = &
              energy_by_slot(reactions(j)%slot) + e
         energy_by_chain(reactions(j)%chain) = &
              energy_by_chain(reactions(j)%chain) + e
      end do
end subroutine compute_reaction_energies

! ---------------------------------------------------------------
! dX_i/dt / A_i = sum_j nu_ij R_j for every network species, at the
! abundances x (s_* order) and the composition rate vector rate_vec.
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
! set to zero and skip the relative test). ierr = burn_not_converged
! after max_iterations attempts, burn_singular_jacobian when the Newton
! system is singular. As in kemcom, Jacobian entries below
! min_jacobian_entry in magnitude (off the diagonal) are dropped, and
! if none remain the residual itself is taken as the correction.
pure subroutine burn_implicitly(x_start, rate_vec, dt, min_abundance, &
     absolute_tolerance, relative_tolerance, max_iterations, x, ierr)
      use numerics_lib, only: simeqc, numerics_ok
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
      integer :: i, k, j, c, iterations, solve_status
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
               if (abs(system(i, k)) < min_jacobian_entry) then
                  system(i, k) = 0.0d0
               else
                  off_diagonal = .true.
               end if
            end do
         end do
         if (off_diagonal) then
            call simeqc(system, n_solved+1, n_solved, solve_status)
            if (solve_status /= numerics_ok) then
               ierr = burn_singular_jacobian
               return
            end if
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
            ierr = burn_not_converged
            return
         end if
      end do

! explicit products, from the converged abundances
      call compute_abundance_rates(rate_vec, x, dxdt_per_a)
      do i = n_solved+1, n_network
         x(i) = x_start(i) + dt_a(i)*dxdt_per_a(i)
      end do
end subroutine burn_implicitly


end module network
