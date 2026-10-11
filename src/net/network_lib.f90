!----------------------------------------------------------------------
! network_lib
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
! Reactions. Reaction j proceeds at
!     R_j = b_j * g_j * prod_m X(law_j(m))          [per Gyr per amu]
! where g_j is row rate_row of the rate vector filled by net_lib's
! rates (rows rr_*: N_A<sigma v>, density, screening and the 1/A
! factors are already in it), X are mass fractions, and b_j is 1, or
! the N15(p,alpha) branching fraction f (or 1-f) carried in row
! rr_frac_c12_alpha, for the two outcomes of the N14(p,g) chain. Each
! reaction changes the per-amu abundance of species i by nu_ij, so
!     dX_i/dt = A_i * sum_j nu_ij R_j
! with A_i the integer mass number, the convention the rate rows use.
!
! Implicit step. Over dt the solve finds X with
!     F_i(X) = X_i - X0_i - dt A_i sum_j nu_ij R_j(X) = 0
! by Newton iteration on J_ik = dF_i/dX_k. Adding a reaction means
! adding a row to the table; adding a species means adding it to the
! species list and to the rows that touch it.
module network_lib
      use net_lib, only: rr_pp, rr_he3_he3, rr_he3_he4, &
           rr_c12_p, rr_c13_p, rr_n14_p, rr_o16_p, rr_c13_alpha, &
           rr_c12_alpha, rr_n14_alpha, rr_triple_alpha, rr_frac_c12_alpha
      implicit none
      private

! network species: positions in the abundance vector of the solve
      integer, parameter, public :: s_h1 = 1, s_he3 = 2, s_he4 = 3, &
           s_c12 = 4, s_c13 = 5, s_n14 = 6, s_o16 = 7, s_o18 = 8
! the first n_solved species are solved implicitly; the rest are
! products updated explicitly after convergence
      integer, parameter, public :: n_solved = 7, n_network = 8
      double precision, parameter, public :: mass_number(n_network) = &
           [1.0d0, 3.0d0, 4.0d0, 12.0d0, 13.0d0, 14.0d0, 16.0d0, 18.0d0]

! branch factor b_j
      integer, parameter :: unbranched = 0, n15_to_c12 = 1, n15_to_o16 = 2

! One reaction. law lists the species whose mass fractions multiply the
! rate (repeat a species for a square or cube; 0 = unused). change lists
! (species, nu) pairs: the change in per-amu abundance per reaction
! (0 = unused). emitted_mass is the mass number carried off by particles
! the network does not track (the neutron of C13(a,n)O16), so that every
! row conserves baryon number exactly: sum_i A_i nu_i + emitted_mass = 0
! (checked by network_table_errors, pinned by net/test/test_net).
      integer, parameter :: max_law = 3, max_change = 5
      type :: reaction_t
         character(len=24) :: name
         integer :: rate_row
         integer :: branch
         integer :: law(max_law)
         integer :: change(2*max_change)
         integer :: emitted_mass = 0
      end type reaction_t

      integer, parameter :: n_reactions = 12
      type(reaction_t), parameter :: reactions(n_reactions) = [ &
!                     name                  rate row         branch      rate law                  changes (species, nu)
           reaction_t('pp chain to He3',     rr_pp,           unbranched, [s_h1, s_h1, 0],          &
                      [s_h1, -3, s_he3, 1, 0, 0, 0, 0, 0, 0]), &
           reaction_t('He3(He3,2p)He4',      rr_he3_he3,      unbranched, [s_he3, s_he3, 0],        &
                      [s_h1, 2, s_he3, -2, s_he4, 1, 0, 0, 0, 0]), &
           reaction_t('He3(He4,g) pp II/III', rr_he3_he4,     unbranched, [s_he3, s_he4, 0],        &
                      [s_h1, -1, s_he3, -1, s_he4, 1, 0, 0, 0, 0]), &
           reaction_t('C12(p,g) to C13',     rr_c12_p,        unbranched, [s_h1, s_c12, 0],         &
                      [s_h1, -1, s_c12, -1, s_c13, 1, 0, 0, 0, 0]), &
           reaction_t('C13(p,g)N14',         rr_c13_p,        unbranched, [s_h1, s_c13, 0],         &
                      [s_h1, -1, s_c13, -1, s_n14, 1, 0, 0, 0, 0]), &
           reaction_t('N14(p,g) N15(p,a)C12', rr_n14_p,       n15_to_c12, [s_h1, s_n14, 0],         &
                      [s_h1, -2, s_n14, -1, s_c12, 1, s_he4, 1, 0, 0]), &
           reaction_t('N14(p,g) N15(p,g)O16', rr_n14_p,       n15_to_o16, [s_h1, s_n14, 0],         &
                      [s_h1, -2, s_n14, -1, s_o16, 1, 0, 0, 0, 0]), &
           reaction_t('O16(p,g) O17(p,a)N14', rr_o16_p,       unbranched, [s_h1, s_o16, 0],         &
                      [s_h1, -2, s_o16, -1, s_n14, 1, s_he4, 1, 0, 0]), &
           reaction_t('C12(a,g)O16',         rr_c12_alpha,    unbranched, [s_he4, s_c12, 0],        &
                      [s_he4, -1, s_c12, -1, s_o16, 1, 0, 0, 0, 0]), &
           reaction_t('C13(a,n)O16',         rr_c13_alpha,    unbranched, [s_he4, s_c13, 0],        &
                      [s_he4, -1, s_c13, -1, s_o16, 1, 0, 0, 0, 0], emitted_mass=1), &
           reaction_t('N14(a,g) to O18',     rr_n14_alpha,    unbranched, [s_he4, s_n14, 0],        &
                      [s_he4, -1, s_n14, -1, s_o18, 1, 0, 0, 0, 0]), &
           reaction_t('triple alpha',        rr_triple_alpha, unbranched, [s_he4, s_he4, s_he4],    &
                      [s_he4, -3, s_c12, 1, 0, 0, 0, 0, 0, 0]) ]

      public :: implicit_burn_step, network_table_errors, write_network_table, &
           equilibrium_burn_rates

contains

! ---------------------------------------------------------------
! Self-check of the table: the number of rows that do not conserve
! baryon number, or whose rate law names a species outside the
! implicit solve (an explicit product cannot feed back).
integer function network_table_errors()
      integer :: j

      network_table_errors = 0
      do j = 1, n_reactions
         if (baryon_imbalance(j) /= 0 .or. &
              any(reactions(j)%law > n_solved)) &
              network_table_errors = network_table_errors + 1
      end do
end function network_table_errors

integer function baryon_imbalance(j)
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

      write(unit, '(a)') '# reaction                 rate row  branch  baryon imbalance'
      do j = 1, n_reactions
         write(unit, '(a24, 2x, i8, 2x, i6, 2x, i16)') reactions(j)%name, &
              reactions(j)%rate_row, reactions(j)%branch, baryon_imbalance(j)
      end do
      write(unit, '(a, i0)') '# table errors: ', network_table_errors()
end subroutine write_network_table

! ---------------------------------------------------------------
! Rate of reaction j at abundances x: R_j = b_j g_j prod X(law).
double precision function reaction_rate(j, rate_vec, x)
      integer, intent(in) :: j
      double precision, intent(in) :: rate_vec(:), x(n_network)
      integer :: m

      reaction_rate = branch_factor(j, rate_vec)*rate_vec(reactions(j)%rate_row)
      do m = 1, max_law
         if (reactions(j)%law(m) /= 0) &
              reaction_rate = reaction_rate*x(reactions(j)%law(m))
      end do
end function reaction_rate

! ---------------------------------------------------------------
! dR_j/dX_k: the product rule over the factors of the rate law.
double precision function reaction_rate_derivative(j, k, rate_vec, x)
      integer, intent(in) :: j, k
      double precision, intent(in) :: rate_vec(:), x(n_network)
      integer :: m, m2
      double precision :: term

      reaction_rate_derivative = 0.0d0
      do m = 1, max_law
         if (reactions(j)%law(m) /= k) cycle
         term = branch_factor(j, rate_vec)*rate_vec(reactions(j)%rate_row)
         do m2 = 1, max_law
            if (m2 /= m .and. reactions(j)%law(m2) /= 0) &
                 term = term*x(reactions(j)%law(m2))
         end do
         reaction_rate_derivative = reaction_rate_derivative + term
      end do
end function reaction_rate_derivative

double precision function branch_factor(j, rate_vec)
      integer, intent(in) :: j
      double precision, intent(in) :: rate_vec(:)

      select case (reactions(j)%branch)
      case (n15_to_c12)
         branch_factor = rate_vec(rr_frac_c12_alpha)
      case (n15_to_o16)
         branch_factor = 1.0d0 - rate_vec(rr_frac_c12_alpha)
      case default
         branch_factor = 1.0d0
      end select
end function branch_factor

! ---------------------------------------------------------------
! dX_i/dt / A_i = sum_j nu_ij R_j for every network species.
subroutine abundance_rates(rate_vec, x, dxdt_per_a)
      double precision, intent(in) :: rate_vec(:), x(n_network)
      double precision, intent(out) :: dxdt_per_a(n_network)
      integer :: j, c
      double precision :: r

      dxdt_per_a = 0.0d0
      do j = 1, n_reactions
         r = reaction_rate(j, rate_vec, x)
         do c = 1, 2*max_change - 1, 2
            if (reactions(j)%change(c) == 0) cycle
            dxdt_per_a(reactions(j)%change(c)) = dxdt_per_a(reactions(j)%change(c)) + &
                 reactions(j)%change(c+1)*r
         end do
      end do
end subroutine abundance_rates

! ---------------------------------------------------------------
! One implicit (backward-Euler) burning step of length dt_gyr.
!   x_start   abundances at the start of the step (all network species)
!   rate_vec  the rate vector (rows rr_*) for this zone or zone average
!   x         in: the first guess for the solved species; out: the
!             solution, with the explicit products (O18) updated
! Newton converges when every |correction| < absolute_tolerance and every
! |correction/X| < relative_tolerance (species below min_abundance are
! set to zero and skip the relative test); ierr = 1 after
! max_iterations attempts. As in kemcom, Jacobian entries below 1e-8 in
! magnitude (off the diagonal) are dropped, and if none remain the
! residual itself is taken as the correction.
subroutine implicit_burn_step(x_start, rate_vec, dt_gyr, min_abundance, &
     absolute_tolerance, relative_tolerance, max_iterations, x, ierr)
      use numerics_lib, only: simeqc
      double precision, intent(in) :: x_start(n_network), rate_vec(:), dt_gyr
      double precision, intent(in) :: min_abundance, absolute_tolerance, &
           relative_tolerance
      integer, intent(in) :: max_iterations
      double precision, intent(inout) :: x(n_network)
      integer, intent(out) :: ierr

! system(:, 1:n) is the Jacobian, system(:, n+1) the residual; simeqc
! leaves the Newton correction in system(:, n+1)
      double precision :: system(n_solved, n_solved+1), dt_a(n_network)
      double precision :: dxdt_per_a(n_network), correction(n_solved)
      double precision :: max_abs_change, max_rel_change
      integer :: i, k, j, c, iterations
      logical :: off_diagonal

      ierr = 0
      dt_a = mass_number*dt_gyr
      iterations = 0
      do
         call abundance_rates(rate_vec, x, dxdt_per_a)
         do i = 1, n_solved
            system(i, n_solved+1) = x(i) - x_start(i) - dt_a(i)*dxdt_per_a(i)
         end do
! J_ik = delta_ik - dt A_i sum_j nu_ij dR_j/dX_k
         system(:, 1:n_solved) = 0.0d0
         do j = 1, n_reactions
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
      call abundance_rates(rate_vec, x, dxdt_per_a)
      do i = n_solved+1, n_network
         x(i) = x_start(i) + dt_a(i)*dxdt_per_a(i)
      end do
end subroutine implicit_burn_step

!----------------------------------------------------------------------
! equilibrium_burn_rates
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqburn.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). 2026-10: moved here from
! core/burn_lib, with the start-of-step composition as an argument.
!
! Computes explicit hydrogen- and helium-burning rates (dX/dt, dY/dt,
! dC/dt, dO/dt) assuming equilibrium He3 and CN-cycle abundances --
! the classic approximation that lets YREC track only a handful of
! "slow" abundances instead of every CNO isotope individually (see
! dburn below for the analogous deuterium-burning treatment). Used as
! the initial guess for the burning rates at the start of a timestep;
! may be supplemented elsewhere by a fully implicit non-equilibrium
! calculation.
! equilibrium_burn_rates (formerly eqburn)
subroutine equilibrium_burn_rates(rate_by_zone, xa_start, shell_mass, shell_temperature, zone_begin, &
     zone_end, dc_dt, do_dt, dx_dt, dy_dt, equilibrium_xc12, &
     equilibrium_xo16, hydrogen_fraction, metal_fraction)

      use star_info_lib, only: json, n_species_extended, i_h1, i_he4, &
           i_metals, i_he3, i_c12, i_c13, i_n14, i_o16
      use net_lib, only: num_rate_rows, net_settings
      implicit none

! rate_by_zone(rr_*, zone): the per-zone rate vectors written by
! net_lib's network_rate_vector (2026 W3 -- formerly 9 separate rate_*(json)
! arguments); rows are the rr_* indices of net_lib. Only
! the nine rows named in the use statement are read here.
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
! xa_start: the start-of-step composition (star%xa_start at every caller;
! an argument since the move to net/, 2026-10)
      double precision, intent(in) :: xa_start(n_species_extended, json)
      double precision, intent(in) :: shell_mass(json), &
           shell_temperature(json)
      integer, intent(in) :: zone_begin, zone_end
      double precision, intent(out) :: dc_dt, do_dt, dx_dt, dy_dt, &
           equilibrium_xc12, equilibrium_xo16, hydrogen_fraction, &
           metal_fraction



      double precision :: zone_avg_abundance(11)
      double precision :: total_shell_mass
      integer :: species_idx, zone_idx
      double precision :: pp_rate, he3_he3_rate, he3_he4_rate, c12_p_rate, &
           c13_p_rate, n14_p_rate, o16_p_rate, c12_alpha_rate, &
           triple_alpha_rate
      double precision :: local_xhe3, helium_fraction, local_xc13, &
           local_xn14
      double precision :: pp_reaction_term, he3_he4_reaction_term, &
           equilibrium_xhe3
      double precision :: cno_sum, cn_ratio_c12, cn_ratio_c13
      double precision :: triple_alpha_term, c12_alpha_term

!     SKIP BURNING CALCULATIONS IF STARTING SHELL BELOW T CUTOFF FOR REACTIONS.
      dx_dt = 0.0d0
      dy_dt = 0.0d0
      dc_dt = 0.0d0
      do_dt = 0.0d0
      if (shell_temperature(zone_begin).lt.net_settings%logT_cutoffs(1)) return
!     COMPUTE EXPLICIT HYDROGEN AND HELIUM BURNING RATES ASSUMING EQUILIBRIUM
!     HE3 AND CN CYCLE ABUNDANCES.
!     THIS IS USED FOR AN INITIAL GUESS AT THE BURNING RATES AT THE START
!     OF THE TIME STEP, AND CAN BE SUPPLEMENTED BY A FULLY IMPLICIT
!     NON-EQUILIBRIUM CALCULATION IN THE 3RD AND 4TH LEVELS OF ITERATION.
      if (zone_begin.ne.zone_end) then
!        HOMOGENIZE CONVECTION ZONES.
!        zone_avg_abundance IS THE MASS-WEIGHTED AVERAGE ABUNDANCE FOR
!        THE CZ. total_shell_mass IS THE TOTAL MASS OF THE CZ.
!        INITIALIZE SUMS.
         total_shell_mass = 0.0d0
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = 0.0d0
         end do
         do zone_idx = zone_begin, zone_end
            total_shell_mass = total_shell_mass + shell_mass(zone_idx)
            do species_idx = 1, 11
               zone_avg_abundance(species_idx) = &
                    zone_avg_abundance(species_idx) + &
                    xa_start(species_idx,zone_idx)*shell_mass(zone_idx)
            end do
         end do
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = &
                 zone_avg_abundance(species_idx)/total_shell_mass
         end do
      else
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = &
                 xa_start(species_idx,zone_begin)
         end do
      end if
      if (zone_begin.eq.zone_end) then
!        PP
         pp_rate = rate_by_zone(rr_pp, zone_begin)
!        HE3,HE3
         he3_he3_rate = rate_by_zone(rr_he3_he3, zone_begin)
!        HE3,HE4
         he3_he4_rate = rate_by_zone(rr_he3_he4, zone_begin)
!        C12,P
         c12_p_rate = rate_by_zone(rr_c12_p, zone_begin)
!        C13,P
         c13_p_rate = rate_by_zone(rr_c13_p, zone_begin)
!        N14,P + N15,P
         n14_p_rate = rate_by_zone(rr_n14_p, zone_begin)
!        O16,P + O17,P.
         o16_p_rate = rate_by_zone(rr_o16_p, zone_begin)
!        C13,ALPHA -- not used (row rr_c13_alpha not read here)
!        O16,ALPHA (NOT USED)
!        C12,ALPHA
         c12_alpha_rate = rate_by_zone(rr_c12_alpha, zone_begin)
!        N14,ALPHA -- not used (row rr_n14_alpha not read here)
!        TRIPLE ALPHA
         triple_alpha_rate = rate_by_zone(rr_triple_alpha, zone_begin)
!        C12,C12 (NOT USED)
!        BRANCHING RATIO FOR N15,P :
!        f = FRACTION GOING TO C12+ALPHA, 1-f = FRACTION GOING TO O16
!        -- not used (row rr_frac_c12_alpha not read here)
      else
!        USE THE MASS-WEIGHTED AVERAGE RATES FOR THE CZ
         pp_rate = 0.0d0
         he3_he3_rate = 0.0d0
         he3_he4_rate = 0.0d0
         c12_p_rate = 0.0d0
         c13_p_rate = 0.0d0
         n14_p_rate = 0.0d0
         o16_p_rate = 0.0d0
         c12_alpha_rate = 0.0d0
         triple_alpha_rate = 0.0d0
         do zone_idx = zone_begin, zone_end
            pp_rate = pp_rate + shell_mass(zone_idx)*rate_by_zone(rr_pp, zone_idx)
            he3_he3_rate = he3_he3_rate + &
                 shell_mass(zone_idx)*rate_by_zone(rr_he3_he3, zone_idx)
            he3_he4_rate = he3_he4_rate + &
                 shell_mass(zone_idx)*rate_by_zone(rr_he3_he4, zone_idx)
            c12_p_rate = c12_p_rate + shell_mass(zone_idx)*rate_by_zone(rr_c12_p, zone_idx)
            c13_p_rate = c13_p_rate + shell_mass(zone_idx)*rate_by_zone(rr_c13_p, zone_idx)
            n14_p_rate = n14_p_rate + shell_mass(zone_idx)*rate_by_zone(rr_n14_p, zone_idx)
            o16_p_rate = o16_p_rate + shell_mass(zone_idx)*rate_by_zone(rr_o16_p, zone_idx)
            c12_alpha_rate = c12_alpha_rate + &
                 shell_mass(zone_idx)*rate_by_zone(rr_c12_alpha, zone_idx)
            triple_alpha_rate = triple_alpha_rate + &
                 shell_mass(zone_idx)*rate_by_zone(rr_triple_alpha, zone_idx)
         end do
         pp_rate = pp_rate/total_shell_mass
         he3_he3_rate = he3_he3_rate/total_shell_mass
         he3_he4_rate = he3_he4_rate/total_shell_mass
         c12_p_rate = c12_p_rate/total_shell_mass
         c13_p_rate = c13_p_rate/total_shell_mass
         n14_p_rate = n14_p_rate/total_shell_mass
         o16_p_rate = o16_p_rate/total_shell_mass
         c12_alpha_rate = c12_alpha_rate/total_shell_mass
         triple_alpha_rate = triple_alpha_rate/total_shell_mass
      end if

!     INITIAL ABUNDANCES OF SPECIES.

      hydrogen_fraction = zone_avg_abundance(i_h1)
      local_xhe3 = zone_avg_abundance(i_he3)
      helium_fraction = zone_avg_abundance(i_he4)
      metal_fraction = zone_avg_abundance(i_metals)
      equilibrium_xc12 = zone_avg_abundance(i_c12)
      local_xc13 = zone_avg_abundance(i_c13)
      local_xn14 = zone_avg_abundance(i_n14)
      equilibrium_xo16 = zone_avg_abundance(i_o16)
      if (hydrogen_fraction.gt.1.0d-10 .and. &
           shell_temperature(zone_begin).gt.net_settings%logT_cutoffs(2)) then
!        FIND EQUILIBRIUM HELIUM-3 ABUNDANCE, USING THE QUADRATIC FORMULA.
!        THE EQUATION IS -2*R(3,3)*XHE3**2 - R(3,4)*XHE3*XHE4 + R(1,1)X**2 = 0.
         pp_reaction_term = pp_rate*hydrogen_fraction**2
         he3_he4_reaction_term = he3_he4_rate*helium_fraction
         equilibrium_xhe3 = 0.25d0*(sqrt(he3_he4_reaction_term**2 + &
              8.0d0*he3_he3_rate*pp_reaction_term) - he3_he4_reaction_term)/ &
              he3_he3_rate
!        USE THE MINIMUM OF THE CURRENT AND EQUILIBRIUM HE3 ABUNDANCE -
!        NEEDED TO AVOID SPURIOUSLY HIGH BURNING RATES IN OUTER LAYERS,
!        WHERE THE TIMESCALE FOR HELIUM-3 TO REACH EQUILIBRIUM IS LONG.
         local_xhe3 = min(local_xhe3, equilibrium_xhe3)
!        PP BURNING.
         dx_dt = -3.0d0*pp_reaction_term + 2.0d0*he3_he3_rate*local_xhe3**2 &
              - he3_he4_reaction_term*local_xhe3
!        CN CYCLE (only on the hydrogen-burning branch).
         if (shell_temperature(zone_begin).gt.net_settings%logT_cutoffs(3)) then
!           FIND EQUILIBRIUM C12,C13,N14 ABUNDANCES TREATING CN PROCESSING AS
!           A CLOSED LOOP.
            cno_sum = equilibrium_xc12/1.2d1 + local_xc13/1.3d1 + &
                 local_xn14/1.4d1
!           FIND EXPLICIT RATES FOR P+C12,P+C13,P+N14(==R121,R131,R141)
            cn_ratio_c12 = n14_p_rate/c12_p_rate
            cn_ratio_c13 = n14_p_rate/c13_p_rate
!           IN EQUILIBRIUM WE HAVE
!           XC12 = XN14{R(1,14)/R(1,12)}
!           XC13 = XN14{R(1,14)/R(1,13)}
!           AND FROM CONSERVATION OF SPECIES XC12/12 +XC13/13 + XN14/14 IS
!           UNCHANGED.
!           EXPRESSING C12 AND C13 AS FUNCTIONS OF N14, SOLVE FOR THE NEW N14
!           EQUILIBRIUM ABUNDANCE.
            local_xn14 = cno_sum/(cn_ratio_c12/1.2d1 + cn_ratio_c13/1.3d1 + &
                 1.0d0/1.4d1)
            equilibrium_xc12 = local_xn14*cn_ratio_c12
            local_xc13 = local_xn14*cn_ratio_c13
            dx_dt = dx_dt - c12_p_rate*hydrogen_fraction*equilibrium_xc12 - &
                 c13_p_rate*hydrogen_fraction*local_xc13 - &
                 2.0d0*n14_p_rate*hydrogen_fraction*local_xn14 - &
                 2.0d0*o16_p_rate*hydrogen_fraction*equilibrium_xo16
         end if
      end if
!     HELIUM BURNING REACTIONS (dy_dt, dc_dt, do_dt stay at the zero set
!     above when helium does not burn).
      if (helium_fraction.gt.1.0d-10 .and. &
           shell_temperature(zone_begin).gt.net_settings%logT_cutoffs(4)) then
         triple_alpha_term = triple_alpha_rate*helium_fraction**3
         c12_alpha_term = c12_alpha_rate*helium_fraction*equilibrium_xc12
         dy_dt = -4.0d0*(c12_alpha_term + 3.0d0*triple_alpha_term)
         dc_dt = 1.2d1*(triple_alpha_term - c12_alpha_term)
         do_dt = 1.6d1*c12_alpha_term
      end if

end subroutine equilibrium_burn_rates

end module network_lib
