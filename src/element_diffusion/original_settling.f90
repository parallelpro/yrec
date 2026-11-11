!----------------------------------------------------------------------
! original_settling
!----------------------------------------------------------------------
! The original settling routines (use_new_diffusion_routines false; no
! shipped deck selects them): gravitational_settling settles hydrogen
! and helium, and the metals (use_diffusion_z) as fully ionized iron,
! over one step, in the formalism of Bahcall & Loeb (1990). Over the
! settling region diffuse_elements found:
!   1. gravitational_settling_setup: the timestep, d ln P/dr, the
!      enclosed mass and the total mass in Bahcall & Loeb units, in
!      place; at every zone the coefficients D1 and D2 of
!      dX/dt = (1/(rho r^2)) d/dr [D1 + D2 dX/dr] and their derivatives
!      in X (in Z for the metals): the older analytic ones, or Thoul et
!      al.'s (1994) from Burgers' equations with coulomb_log_choice's
!      logarithm, or their fit (use_thoul_fit);
!   2. model_to_equal (settling_grid): the equal grids and the
!      coefficients and abundances on them;
!   3. the drift by the two Lax-Wendroff steps (diffusion_scheme), for X
!      and for Z, the midpoints' D1 corrected by its derivative for the
!      half step's change;
!   4. iterate_hydrogen_diffusion: hydrogen's diffusion term, the implicit
!      step iterated with D2 corrected by its derivative, until the
!      largest change is below settling_tolerance (settling_num_iterations
!      passes at most); the metals have no diffusion term here;
!   5. equal_to_model (settling_grid): the changes back on the model, and
!      the units back.
! Four defects found in 2026 are recorded, not fixed, until the authors
! decide to fix and test these routines or retire them
! (docs/bug_sweep.md M243-M246: the Noerdlinger term's sign, ln Lambda
! twice with choices 2 and 3, the iteration's midpoint change one slot
! off, the convergence on the largest signed change), and M33 (the
! metals' coefficients unset without use_thoul_diffusion).
module original_settling
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use element_diffusion_configuration, only: diffusion_settings
      use element_diffusion_scratch, only: diffusion_scratch
      use settling_grid, only: model_to_equal, equal_to_model
      use diffusion_scheme, only: lax_wendroff_step1, lax_wendroff_step2, &
           implicit_diffusion_coeffs
      use thoul_solver, only: thoul_diffusion
      use coulomb_logarithm, only: compute_concentrations, compute_coulomb_logarithms
      use element_diffusion_messages, only: note_solve_pass, report_unconverged, &
           report_numerics_failure
      implicit none
      private
      public :: gravitational_settling

contains

! ---------------------------------------------------------------
! One step of timestep seconds over zone_begin..zone_end of a model of
! num_zones zones: the mass fractions composition(species, zone) change
! in place. The structure as diffuse_elements describes it; timestep,
! dlnp_dr, enclosed_mass and total_mass are converted to Bahcall & Loeb
! units and back in place. ierr: a numerics routine's status, after its
! message.
subroutine gravitational_settling(diff_scr, timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, num_zones, total_mass, &
     zone_begin, zone_end, ierr)
      type(diffusion_scratch), intent(inout) :: diff_scr
      real(dp), intent(inout) :: timestep, total_mass
      real(dp), intent(inout) :: composition(:,:), dlnp_dr(:), enclosed_mass(:)
      real(dp), intent(in) :: log_radius(:), log_density(:), log_temperature(:), del_grad(:)
      integer, intent(in) :: num_zones, zone_begin, zone_end
      integer, intent(out) :: ierr
! on the model's zones: the radius and the temperature in Bahcall & Loeb
! units, hydrogen's D1 and D2 and their derivatives in X
      real(dp) :: radius_bl(json), temperature_bl(json), coeff1(json), coeff2(json), &
           coeff1_dx(json), coeff2_dx(json)
! on the equal grids (_mid: at the midpoints): the centres' radius, the
! enclosed mass, hydrogen's coefficients and derivatives, X (then its
! change), and X and Z at the step's start
      real(dp) :: eq_radius(json), eq_mass(json), eq_mass_mid(json), eq_coeff1(json), &
           eq_coeff1_mid(json), eq_coeff2_mid(json), eq_coeff1_dx_mid(json), &
           eq_coeff2_dx_mid(json), eq_hydrogen(json), eq_hydrogen_mid(json), &
           hydrogen_start(json), metal_start(json)
      real(dp) :: grid_spacing
      integer :: n, i

      ierr = 0
      call gravitational_settling_setup(diff_scr, timestep, dlnp_dr, log_radius, &
           log_density, enclosed_mass, log_temperature, del_grad, num_zones, total_mass, &
           coeff1, coeff2, composition, radius_bl, temperature_bl, coeff1_dx, coeff2_dx, &
           ierr)
      if (ierr /= 0) return
      call model_to_equal(diff_scr, coeff1, coeff2, composition, radius_bl, enclosed_mass, &
           zone_begin, zone_end, num_zones, coeff1_dx, coeff2_dx, grid_spacing, &
           eq_coeff1, eq_coeff1_mid, eq_coeff2_mid, eq_mass, eq_mass_mid, &
           eq_coeff1_dx_mid, eq_coeff2_dx_mid, eq_radius, eq_hydrogen, eq_hydrogen_mid, n)
! The abundances at the step's start: the changes, not the new runs, go
! back to the model.
      hydrogen_start(1:n) = eq_hydrogen(1:n)
      if (diffusion_settings%use_diffusion_z) metal_start(1:n) = diff_scr%eq_metal_fraction(1:n)
! Lax-Wendroff step 1, X and Z: the half step's change at the midpoints,
! and the midpoints' D1 corrected for it, D1 += (dD1/dX) dX.
      call lax_wendroff_step1(timestep, eq_coeff1, eq_mass, n, total_mass, eq_hydrogen_mid)
      if (diffusion_settings%use_diffusion_z) then
         call lax_wendroff_step1(timestep, diff_scr%metal_coeff1, eq_mass, n, total_mass, &
              diff_scr%eq_metal_fraction_mid)
      end if
      do i = 1, n-1
         eq_coeff1_mid(i) = eq_coeff1_mid(i) + eq_hydrogen_mid(i)*eq_coeff1_dx_mid(i)
      end do
      if (diffusion_settings%use_diffusion_z) then
         do i = 1, n-1
            diff_scr%metal_coeff1_mid(i) = diff_scr%metal_coeff1_mid(i) + &
                 diff_scr%eq_metal_fraction_mid(i)*diff_scr%metal_coeff1_dz_mid(i)
         end do
      end if
! Lax-Wendroff step 2, X and Z.
      call lax_wendroff_step2(timestep, eq_coeff1_mid, eq_mass_mid, eq_hydrogen, n, total_mass)
      if (diffusion_settings%use_diffusion_z) then
         call lax_wendroff_step2(timestep, diff_scr%metal_coeff1_mid, eq_mass_mid, &
              diff_scr%eq_metal_fraction, n, total_mass)
      end if
! Hydrogen's diffusion term.
      call iterate_hydrogen_diffusion(timestep, total_mass, grid_spacing, n, eq_mass_mid, &
           eq_coeff2_mid, eq_coeff2_dx_mid, hydrogen_start, eq_hydrogen, eq_hydrogen_mid, ierr)
      if (ierr /= 0) return
! The changes, back on the model.
      eq_hydrogen(1:n) = eq_hydrogen(1:n) - hydrogen_start(1:n)
      if (diffusion_settings%use_diffusion_z) then
         diff_scr%eq_metal_fraction(1:n) = diff_scr%eq_metal_fraction(1:n) - metal_start(1:n)
      end if
      call equal_to_model(diff_scr, timestep, eq_radius, eq_hydrogen, zone_begin, zone_end, &
           n, composition, dlnp_dr, radius_bl, enclosed_mass, num_zones, total_mass)
end subroutine gravitational_settling

! ---------------------------------------------------------------
! Hydrogen's diffusion term over timestep on the n centres (cells as in
! diffusion_scheme; alpha_i = 4 pi timestep/(grid_spacing dm_i)): the
! implicit step, X - alpha [D (X_(i+1) - X_i) - ...] = X' (X' the
! drift's result, eq_hydrogen on entry), passed after pass with each
! midpoint's D2 (coeff2_mid) corrected by its derivative (coeff2_dx_mid)
! times the change of X there since the step's start (hydrogen_start),
! until the largest signed change of a pass is below settling_tolerance,
! at most settling_num_iterations passes (each pass under solver
! diagnostics, and a solve that misses the tolerance, reported).
! (The midpoint change of pass k is stored one slot off from the
! midpoint it belongs to, and the test takes the signed change:
! docs/bug_sweep.md M245, M246.) mid_change: work space.
subroutine iterate_hydrogen_diffusion(timestep, total_mass, grid_spacing, n, mass_mid, &
     coeff2_mid, coeff2_dx_mid, hydrogen_start, eq_hydrogen, mid_change, ierr)
      use physical_constants, only: c4pi
      use numerics_lib, only: solve_tridiagonal
      real(dp), intent(in) :: timestep, total_mass, grid_spacing
      integer, intent(in) :: n
      real(dp), intent(in) :: mass_mid(:), coeff2_dx_mid(:), hydrogen_start(:)
      real(dp), intent(inout) :: coeff2_mid(:), eq_hydrogen(:), mid_change(:)
      integer, intent(out) :: ierr
! alpha_i, the diagonals, the drift's result (the right-hand side) and
! the previous pass's X
      real(dp) :: alpha(json), sub_diag(json), diag(json), super_diag(json), &
           hydrogen_drift(json), hydrogen_prev(json)
      real(dp) :: alpha_scale, change, dx
      integer :: i, pass, change_shell

      ierr = 0
      hydrogen_drift(1:n) = eq_hydrogen(1:n)
      alpha_scale = c4pi*timestep/grid_spacing
      alpha(1) = alpha_scale/mass_mid(1)
      do i = 2, n-1
         alpha(i) = alpha_scale/(mass_mid(i) - mass_mid(i-1))
      end do
      alpha(n) = alpha_scale/(total_mass - mass_mid(n-1))
      hydrogen_prev(1:n) = hydrogen_start(1:n)
      change = 0.0d0
      change_shell = 0
      do pass = 1, diffusion_settings%settling_num_iterations
         do i = 2, n
            mid_change(i) = 0.5d0*(eq_hydrogen(i) + eq_hydrogen(i-1) - hydrogen_prev(i) - &
                 hydrogen_prev(i-1))
         end do
         hydrogen_prev(1:n) = eq_hydrogen(1:n)
         call implicit_diffusion_coeffs(alpha, coeff2_mid, mid_change, coeff2_dx_mid, &
              sub_diag, diag, super_diag, n)
         call solve_tridiagonal(sub_diag, diag, super_diag, hydrogen_drift, n, eq_hydrogen, ierr)
         if (ierr /= 0) then
            call report_numerics_failure('solve_tridiagonal', ierr)
            return
         end if
         change = 0.0d0
         change_shell = 0
         do i = 1, n
            dx = eq_hydrogen(i) - hydrogen_prev(i)
            if (dx > change) then
               change = dx
               change_shell = i
            end if
         end do
         call note_solve_pass(pass, change, change_shell)
         if (change < diffusion_settings%settling_tolerance) exit
      end do
      if (pass > diffusion_settings%settling_num_iterations) then
         call report_unconverged('GRSETT', 'X', diffusion_settings%settling_tolerance, &
              diffusion_settings%settling_num_iterations, change, change_shell)
      end if
end subroutine iterate_hydrogen_diffusion

! ---------------------------------------------------------------
! The setup of gravitational_settling: timestep_seconds, dlnp_dr,
! enclosed_mass and total_mass converted to Bahcall & Loeb units in
! place, radius_bl and temperature_bl the unlogged radius and
! temperature in those units, and at every zone hydrogen's coefficients
! (coeff1 = D1, coeff2 = D2) and their derivatives in X (coeff1_dx,
! coeff2_dx), with the metals' in diff_scr (use_diffusion_z). With
! f = fgry r^2 T^(5/2)/ln Lambda (fgrz for the metals), X, Z + He3 and
! del = d ln T/d ln P:
!   the older analytic coefficients (use_thoul_diffusion off; Bahcall &
!   Loeb, BC 5/90),
!     D1 = f d ln P/dr X (1 - X - Z - He3) [5/4 + 6 del (X + 0.32)/
!          (5.4 + 6.3 X - 4.5 X^2)],
!     D2 = f (X + 3)/(5 X^2 + 8 X + 3);
!   Thoul et al.'s (use_thoul_diffusion on), from Burgers' equations
!   (thoul_diffusion) or their fit in X (use_thoul_fit), D1 = c f
!   d ln P/dr X (A_P + A_T del), D2 = c f A_H, c = constant_mixing_coeff,
!   and for iron D1 = r c f d ln P/dr Z (A_P + A_T del) and
!   D2 = c |f A_H|, r = constant_settling_reduction.
! ln Lambda is coulomb_log_choice's (coulomb_log_of_choice); with choice 4
! and Burgers' equations every pair has its own (coulomb_logarithm).
! ierr: ludcmp's status, after its message.
subroutine gravitational_settling_setup(diff_scr, timestep_seconds, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, num_zones, total_mass, &
     coeff1, coeff2, composition, radius_bl, temperature_bl, coeff1_dx, coeff2_dx, ierr)
      use composition_layout, only: i_h1, i_he4, i_metals, i_he3
      use species_table, only: thoul_h1, thoul_he4, thoul_fe, thoul_electron, &
           thoul_col_h, thoul_col_he, thoul_col_metal
      use bahcall_loeb_units, only: bl_radius_scale, bl_mass_scale, bl_temp_scale, &
           bl_time_scale
      use physical_constants, only: ln10
      use math_lib
      type(diffusion_scratch), intent(inout) :: diff_scr
      real(dp), intent(inout) :: timestep_seconds, total_mass
      real(dp), intent(inout) :: dlnp_dr(:), enclosed_mass(:)
      real(dp), intent(in) :: log_radius(:), log_density(:), log_temperature(:), del_grad(:)
      integer, intent(in) :: num_zones
      real(dp), intent(out) :: coeff1(:), coeff2(:), coeff1_dx(:), coeff2_dx(:)
      real(dp), intent(in) :: composition(:,:)
      real(dp), intent(out) :: radius_bl(:), temperature_bl(:)
      integer, intent(out) :: ierr
! the four Thoul columns: hydrogen, helium, iron, the electrons (MHP 8/94)
      integer, parameter :: num_species = 4
      real(dp), parameter :: atomic_weight(num_species) = [thoul_h1%weight, &
           thoul_he4%weight, thoul_fe%weight, thoul_electron%weight]
      real(dp), parameter :: atomic_charge(num_species) = [thoul_h1%charge, &
           thoul_he4%charge, thoul_fe%charge, thoul_electron%charge]
      real(dp) :: settling_ap(num_species), settling_at(num_species), &
           settling_ac(num_species,num_species), coulomb_log(num_species,num_species), &
           species_mass_fraction(num_species), concentration(num_species)
      real(dp) :: ln_lambda, x, x_sq, z_plus_he3, metal_fraction, x_times_z, denominator, &
           prefactor, coeff_p, coeff_t, helium_ah, iron_ah, dap_dx, dat_dx, dac_dx, &
           rho_cgs, t_cgs
      integer :: k

      ierr = 0
      do k = 1, num_zones
         radius_bl(k) = exp(ln10*log_radius(k))*bl_radius_scale
         temperature_bl(k) = exp(ln10*log_temperature(k))*bl_temp_scale
         enclosed_mass(k) = enclosed_mass(k)*bl_mass_scale
         dlnp_dr(k) = dlnp_dr(k)/bl_radius_scale
      end do
      timestep_seconds = timestep_seconds/bl_time_scale
      total_mass = total_mass*bl_mass_scale

      do k = 1, num_zones
         x = composition(i_h1,k)
         ln_lambda = coulomb_log_of_choice(diffusion_settings%coulomb_log_choice, x, &
              log_density(k), log_temperature(k))
         prefactor = diffusion_settings%fgry*radius_bl(k)**2*pow(temperature_bl(k), 2.5d0)/ln_lambda
         z_plus_he3 = composition(i_metals,k) + composition(i_he3,k)
         metal_fraction = composition(i_metals,k)
         x_sq = x*x
         if (.not. diffusion_settings%use_thoul_diffusion) then
! the older analytic coefficients
            x_times_z = x*z_plus_he3
            denominator = 5.4d0 + 6.3d0*x - 4.5d0*x_sq
            coeff1(k) = prefactor*dlnp_dr(k)* &
                 (x - x_sq - x_times_z)*(1.25d0 + &
                 del_grad(k)*6.0d0*(x + 0.32d0)/denominator)
            coeff2(k) = prefactor*(x + 3.0d0)/ &
                 (5.0d0*x_sq + 8.0d0*x + 3.0d0)
            coeff1_dx(k) = prefactor*dlnp_dr(k)* &
                 ( (1.0d0 - 2.0d0*x - z_plus_he3)*(1.25d0 + &
                 (6.0d0*del_grad(k)*(x + 0.32d0))/denominator) + &
                 (x - x_sq - x_times_z)*6.0d0* &
                 del_grad(k)*(3.384d0 + 2.88d0*x + 4.5d0*x_sq)/ &
                 denominator**2 )
            coeff2_dx(k) = -prefactor*(5.0d0*x_sq + &
                 3.0d1*x + 2.1d1)/ &
                 (5.0d0*x_sq + 8.0d0*x + 3.0d0)**2
         else
! Thoul et al.'s coefficients: Burgers' equations, or their fit
            species_mass_fraction(thoul_col_h) = composition(i_h1,k)
            species_mass_fraction(thoul_col_he) = composition(i_he4,k)
            species_mass_fraction(thoul_col_metal) = composition(i_metals,k)
            if (.not. diffusion_settings%use_thoul_fit) then
               if (diffusion_settings%coulomb_log_choice == 4) then
                  rho_cgs = exp(ln10*log_density(k))
                  t_cgs = exp(ln10*log_temperature(k))
                  call compute_concentrations(atomic_weight, atomic_charge, &
                       species_mass_fraction, concentration)
                  call compute_coulomb_logarithms(atomic_weight, atomic_charge, &
                       concentration, rho_cgs, t_cgs, coulomb_log)
               else
                  coulomb_log = ln_lambda
               end if
               call thoul_diffusion(atomic_weight, atomic_charge, species_mass_fraction, &
                    coulomb_log, settling_ap, settling_at, settling_ac, ierr)
               if (ierr /= 0) then
                  call report_numerics_failure('ludcmp', ierr)
                  return
               end if
               coeff_p = -settling_ap(thoul_col_h)
               coeff_t = -del_grad(k)*settling_at(thoul_col_h)
            else
               coeff_p = 1.58d0 - 2.42d0*x + 0.844d0*x_sq
               coeff_t = del_grad(k)*(1.90d0 - 2.69d0*x + 0.805d0*x_sq)
            end if
            helium_ah = 1.15d0 - 1.42d0*x + 0.647d0*x_sq
            dap_dx = -2.42d0 + 1.688d0*x
            dat_dx = del_grad(k)*(-2.69d0 + 1.61d0*x)
            dac_dx = -1.42d0 + 1.294d0*x
            coeff1(k) = diffusion_settings%constant_mixing_coeff*prefactor*dlnp_dr(k)* &
                 x*(coeff_p + coeff_t)
            coeff2(k) = diffusion_settings%constant_mixing_coeff*prefactor*helium_ah
            coeff1_dx(k) = prefactor*dlnp_dr(k)* &
                 (coeff_p + coeff_t + x*(dap_dx + dat_dx))
            coeff2_dx(k) = prefactor*dac_dx
         end if
! The metals, as fully ionized iron (with Thoul et al.'s coefficients).
         if (diffusion_settings%use_diffusion_z) then
            prefactor = diffusion_settings%fgrz*radius_bl(k)**2*pow(temperature_bl(k), 2.5d0)/ln_lambda
            if (diffusion_settings%use_thoul_diffusion) then
               if (diffusion_settings%use_thoul_fit) then
                  coeff_p = -0.157d0 - 0.511d0*x + 0.389d0*x_sq
                  coeff_t = del_grad(k)*(-1.36d0 - 1.42d0*x + &
                       0.549d0*x_sq)
               else
                  coeff_p = -settling_ap(thoul_col_metal)
                  coeff_t = -del_grad(k)*settling_at(thoul_col_metal)
               end if
               iron_ah = -0.0375d0 - 0.193d0*x + 0.107d0*x_sq
               diff_scr%zone_metal_coeff1(k) = diffusion_settings%constant_settling_reduction* &
                    diffusion_settings%constant_mixing_coeff* &
                    prefactor*dlnp_dr(k)*metal_fraction*(coeff_p + coeff_t)
! (D2 must be positive)
               diff_scr%zone_metal_coeff2(k) = diffusion_settings%constant_mixing_coeff* &
                    abs(prefactor*iron_ah)
               diff_scr%zone_metal_coeff1_dz(k) = diffusion_settings%constant_mixing_coeff* &
                    prefactor*dlnp_dr(k)*(coeff_p + coeff_t)
               diff_scr%zone_metal_coeff2_dz(k) = 0.0d0
            end if
         end if
! Noerdlinger's logarithm depends on X: its term in dD1/dX (with the
! wrong sign, docs/bug_sweep.md M243).
         if (diffusion_settings%coulomb_log_choice == 2) then
            coeff1_dx(k) = coeff1_dx(k) + coeff1(k)*1.5d0/ &
                 (ln_lambda*(3.0d0*x + 1.0d0))
         end if
      end do
end subroutine gravitational_settling_setup

! ---------------------------------------------------------------
! ln Lambda of coulomb_log_choice at hydrogen fraction x, log10 density
! and log10 temperature: 2 Noerdlinger's (1977, A&A 57, 407),
! -19.7 - ln[4/(3X + 1)]/2 - ln(rho)/2 + 1.5 ln T; 3 Loeb's (1989, Phys.
! Rev. D 39, 1009), -19.105747 - ln(rho)/2 + 1.5 ln T; otherwise 2.2.
pure function coulomb_log_of_choice(choice, x, log_density, log_temperature) &
     result(ln_lambda)
      use physical_constants, only: ln10
      use math_lib
      integer, intent(in) :: choice
      real(dp), intent(in) :: x, log_density, log_temperature
      real(dp) :: ln_lambda

      if (choice == 2) then
         ln_lambda = -1.97d1 - 0.5d0*log(1.0d0/(0.75d0*x + 0.25d0)) &
              - 0.5d0*log_density*ln10 + 1.5d0*log_temperature*ln10
      else if (choice == 3) then
         ln_lambda = -1.9105747d1 - 0.50d0*log_density*ln10 + 1.5d0*log_temperature*ln10
      else
         ln_lambda = 2.2d0
      end if
end function coulomb_log_of_choice

end module original_settling

! Sources
! - grsett and setup_grsett: MHP (12/89, after Bahcall & Loeb 1990; 6/90
!   the hydrogen-exhausted core; 3/94 the metals; 8/94 Thoul et al.'s
!   coefficients; 10/02 X set before shell 1), BC (5/90, the analytic
!   coefficients valid for all X with a variable ln Lambda), CFD (10/09,
!   constant_mixing_coeff and constant_settling_reduction).
! - Bahcall, J. N. & Loeb, A. 1990, ApJ 360, 267; Thoul, A. A., Bahcall,
!   J. N. & Loeb, A. 1994, ApJ 421, 828.
