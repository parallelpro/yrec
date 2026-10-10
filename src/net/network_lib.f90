!----------------------------------------------------------------------
! network_lib
!----------------------------------------------------------------------
! The hydrogen- and helium-burning network as data, and the one implicit
! integrator that reads it (2026-10; replaces the hand-written
! 7-species equations and Jacobian of kemcom, mixing/solve_composition).
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
      use rotation_scratch_lib, only: rr_pp, rr_he3_he3, rr_he3_he4, &
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
! (0 = unused).
      integer, parameter :: max_law = 3, max_change = 5
      type :: reaction_t
         character(len=24) :: name
         integer :: rate_row
         integer :: branch
         integer :: law(max_law)
         integer :: change(2*max_change)
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
                      [s_he4, -1, s_c13, -1, s_o16, 1, 0, 0, 0, 0]), &
           reaction_t('N14(a,g) to O18',     rr_n14_alpha,    unbranched, [s_he4, s_n14, 0],        &
                      [s_he4, -1, s_n14, -1, s_o18, 1, 0, 0, 0, 0]), &
           reaction_t('triple alpha',        rr_triple_alpha, unbranched, [s_he4, s_he4, s_he4],    &
                      [s_he4, -3, s_c12, 1, 0, 0, 0, 0, 0, 0]) ]

      public :: implicit_burn_step

contains

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

end module network_lib
