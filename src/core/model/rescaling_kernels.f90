!----------------------------------------------------------------------
! rescaling_kernels
!----------------------------------------------------------------------
! The arithmetic of rescaling a model, on plain arrays: a star-blind
! kernel (docs/GUIDELINES.md, star rule 6). rescale_model
! (model_rescaling.f90) decides what a card rescales, gathers the model
! from the star, calls these, reports and stores; nothing here reads the
! star, prints or keeps state between calls.
!   scale_hydrogen             every shell's X for a new surface X;
!   shift_metals               every shell's Z, and its CNO isotopes, for
!                              a new surface Z;
!   ramp_metals                Z ramped from a core value to the surface's;
!   set_helium_by_difference   Y = 1 - X - Z - X(He3);
!   hb_region_masses           a horizontal-branch model's core, shell and
!                              envelope masses;
!   mass_stretch_factor, stretch_masses_above
!                              a new total mass by stretching the mass
!                              above a base point;
!   move_core_mass             a new core mass at fixed shell and total
!                              mass;
!   compute_envelope_mixture   the envelope mixture from X, Z and He3.
! A composition array xa is (species, shell), in composition_layout's
! slots; a mesh log_mass is log10 of the mass inside each point [g],
! increasing outward. test_model_kernels checks each.
module rescaling_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      use physical_constants, only: ln10
      use composition_layout, only: i_h1, i_he4, i_metals, i_he3, i_c12, i_o18, &
           n_mix_species, ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, ix_h, ix_o, ix_n, &
           ix_ar, ix_ne, ix_he
      implicit none
      private
      public :: hb_masses
      public :: scale_hydrogen, shift_metals, ramp_metals, set_helium_by_difference, &
           hb_region_masses, mass_stretch_factor, stretch_masses_above, move_core_mass, &
           compute_envelope_mixture

! The floors that keep a ratio of mass fractions finite: of the old
! surface X, and of a shell's old Z [mass fraction].
      real(dp), parameter :: x_floor = 1.0d-20, z_floor = 1.d-30

! Atomic weights of the surface-mixture elements, in ix_* slot order
! (Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He) [amu] (Sources).
      real(dp), parameter :: mixture_atomic_weight(n_mix_species) = [23.0d0, 26.99d0, &
           24.32d0, 55.86d0, 28.1d0, 12.015d0, 1.008d0, 16.0d0, 14.01d0, 39.96d0, &
           20.19d0, 4.004d0]

! The masses of a horizontal-branch model's regions [Msun]: the helium
! core, inside the hydrogen shell; the shell; the envelope, from the
! shell's end to the total mass; and the part of it inside the outermost
! point (the fitting point).
      type :: hb_masses
         real(dp) :: core = 0.0d0, shell = 0.0d0, envelope = 0.0d0, &
              inner_envelope = 0.0d0
      end type hb_masses

contains

! ---------------------------------------------------------------
! Every shell's X for a new surface X: scaled by X_new/X_surf, with
! X_surf the old surface X, and capped at X_new,
!   X_i -> min(X_new, X_i X_new/max(X_surf, x_floor)).
pure subroutine scale_hydrogen(hydrogen, new_surface_x, old_surface_x)
      real(dp), intent(inout) :: hydrogen(:)
      real(dp), intent(in) :: new_surface_x, old_surface_x
      real(dp) :: x_rescale_factor
      integer :: zone_idx

      x_rescale_factor = new_surface_x/max(old_surface_x, x_floor)
      do zone_idx = 1, size(hydrogen)
         hydrogen(zone_idx) = min(new_surface_x, hydrogen(zone_idx)*x_rescale_factor)
      end do
end subroutine scale_hydrogen

! ---------------------------------------------------------------
! Every shell's Z shifted by delta_z (the new surface Z less the old) and
! its CNO isotopes (slots i_c12 to i_o18) scaled by the shell's new over
! old Z:
!   Z_i -> max(0, Z_i + dZ),
!   X_k,i -> X_k,i max(0, Z_i + dZ)/(Z_i + z_floor).
! He3 and the light elements are left alone.
pure subroutine shift_metals(xa, delta_z)
      real(dp), intent(inout) :: xa(:,:)
      real(dp), intent(in) :: delta_z
      real(dp) :: z_rescale_factor
      integer :: zone_idx, species_idx

      do zone_idx = 1, size(xa, 2)
         z_rescale_factor = max(0.d0, xa(i_metals,zone_idx)+delta_z)/ &
              (xa(i_metals,zone_idx)+z_floor)
         xa(i_metals,zone_idx) = max(0.0d0, xa(i_metals,zone_idx)+delta_z)
         do species_idx = i_c12, i_o18
            xa(species_idx,zone_idx) = xa(species_idx,zone_idx)*z_rescale_factor
         end do
      end do
end subroutine shift_metals

! ---------------------------------------------------------------
! Z as a function of the mass fraction q = m/M (total_mass, g): core_z
! for q < inner_q, then linear up to surface_z at outer_q,
!   Z(q) = core_z + (q - inner_q)(surface_z - core_z)/(outer_q - inner_q);
! Z beyond outer_q is left alone. X takes up the change, Y being kept:
! X = 1 - Z - X(He3) - Y. The CNO isotopes (i_c12 to i_o18) are
! multiplied by Z_new/(Z_old + z_floor), as in shift_metals.
pure subroutine ramp_metals(xa, log_mass, total_mass, core_z, inner_q, outer_q, &
     surface_z)
      real(dp), intent(inout) :: xa(:,:)
      real(dp), intent(in) :: log_mass(:), total_mass, core_z, inner_q, outer_q, &
           surface_z
      real(dp) :: z_ramp_slope, mass_fraction_local, z_ramp_value, z_rescale_factor
      integer :: zone_idx, species_idx

      z_ramp_slope = (surface_z - core_z)/(outer_q - inner_q)
      do zone_idx = 1, size(log_mass)
         mass_fraction_local = exp10(log_mass(zone_idx))/total_mass
         if (mass_fraction_local < inner_q) then
            z_ramp_value = core_z
         else if (mass_fraction_local < outer_q) then
            z_ramp_value = (mass_fraction_local-inner_q)*z_ramp_slope+core_z
         else
            cycle
         end if
         z_rescale_factor = z_ramp_value/(xa(i_metals,zone_idx)+z_floor)
         xa(i_metals,zone_idx) = z_ramp_value
         xa(i_h1,zone_idx) = 1.0d0-xa(i_metals,zone_idx)-xa(i_he3,zone_idx)- &
              xa(i_he4,zone_idx)
         do species_idx = i_c12, i_o18
            xa(species_idx,zone_idx) = xa(species_idx,zone_idx)*z_rescale_factor
         end do
      end do
end subroutine ramp_metals

! ---------------------------------------------------------------
! Every shell's Y from its X, Z and He3: Y = 1 - X - Z - X(He3).
pure subroutine set_helium_by_difference(xa)
      real(dp), intent(inout) :: xa(:,:)
      integer :: zone_idx

      do zone_idx = 1, size(xa, 2)
         xa(i_he4,zone_idx) = 1.0d0-xa(i_h1,zone_idx)-xa(i_metals,zone_idx)- &
              xa(i_he3,zone_idx)
      end do
end subroutine set_helium_by_difference

! ---------------------------------------------------------------
! The masses of a horizontal-branch model whose hydrogen shell is the
! points shell_begin to shell_end, in units of solar_mass [g]:
!   core            m(shell_begin-1),
!   shell           m(shell_end) - m(shell_begin-1),
!   envelope        M - m(shell_end),
!   inner_envelope  m(nz) - m(shell_end), the envelope inside the
!                   outermost point.
pure function hb_region_masses(log_mass, log_total_mass, shell_begin, shell_end, &
     solar_mass) result(masses)
      real(dp), intent(in) :: log_mass(:), log_total_mass, solar_mass
      integer, intent(in) :: shell_begin, shell_end
      type(hb_masses) :: masses

      masses%core = exp(ln10*log_mass(shell_begin-1))/solar_mass
      masses%envelope = (exp(ln10*log_total_mass)-exp(ln10*log_mass(shell_end)))/ &
           solar_mass
      masses%inner_envelope = (exp(ln10*log_mass(size(log_mass)))- &
           exp(ln10*log_mass(shell_end)))/solar_mass
      masses%shell = (exp(ln10*log_mass(shell_end))-exp(ln10*log_mass(shell_begin-1)))/ &
           solar_mass
end function hb_region_masses

! ---------------------------------------------------------------
! The factor s by which the mass above a base point is stretched to take
! the total mass from old_total_msun to new_total_msun,
!   s = (M_new - m_b)/(M_old - m_b),
! with m_b = 10^log_base_mass [g] and the masses in units of solar_mass.
pure real(dp) function mass_stretch_factor(log_base_mass, new_total_msun, &
     old_total_msun, solar_mass)
      real(dp), intent(in) :: log_base_mass, new_total_msun, old_total_msun, solar_mass

      mass_stretch_factor = (new_total_msun-exp(ln10*log_base_mass)/solar_mass)/ &
           (old_total_msun-exp(ln10*log_base_mass)/solar_mass)
end function mass_stretch_factor

! ---------------------------------------------------------------
! The mass above the base point stretched by the factor s:
!   m_i -> m_b + s (m_i - m_b)  for every point i above base,
! and likewise the total mass M; the base point and those below it stay.
pure subroutine stretch_masses_above(log_mass, log_total_mass, base, stretch_factor)
      real(dp), intent(inout) :: log_mass(:), log_total_mass
      integer, intent(in) :: base
      real(dp), intent(in) :: stretch_factor
      real(dp) :: base_mass
      integer :: zone_idx

      base_mass = exp10(log_mass(base))
      log_total_mass = log10(base_mass+stretch_factor*(exp10(log_total_mass)-base_mass))
      do zone_idx = base+1, size(log_mass)
         log_mass(zone_idx) = log10(base_mass+ &
              stretch_factor*(exp10(log_mass(zone_idx))-base_mass))
      end do
end subroutine stretch_masses_above

! ---------------------------------------------------------------
! A new core mass at fixed hydrogen-shell mass and fixed m(nz) (and so
! fixed total mass), moving mass between the core and the envelope:
! - the core's points (below shell_begin) shift in log m by
!   log10(M_core,new/M_core,old);
! - the shell's points (shell_begin to shell_end) shift in m by
!   dM_core = M_core,new - M_core,old, which keeps the shell's mass;
! - the points above the shell move toward m(nz):
!     m_i -> m_nz - s (m_nz - m_i),
!     s = (m_nz - m_new(shell_end))/(m_nz - m_old(shell_end)).
! The core masses are in units of solar_mass [g].
pure subroutine move_core_mass(log_mass, shell_begin, shell_end, new_core_msun, &
     old_core_msun, solar_mass)
      real(dp), intent(inout) :: log_mass(:)
      integer, intent(in) :: shell_begin, shell_end
      real(dp), intent(in) :: new_core_msun, old_core_msun, solar_mass
      real(dp) :: old_shell_end_mass, log_mass_shift, mass_shift_grams, &
           mass_scale_factor, outer_mass
      integer :: zone_idx, nz

      nz = size(log_mass)
      old_shell_end_mass = exp(ln10*log_mass(shell_end))
      log_mass_shift = log10(new_core_msun/old_core_msun)
      do zone_idx = 1, shell_begin-1
         log_mass(zone_idx) = log_mass(zone_idx) + log_mass_shift
      end do
      mass_shift_grams = (new_core_msun - old_core_msun)*solar_mass
      do zone_idx = shell_begin, shell_end
         log_mass(zone_idx) = log10(exp(ln10*log_mass(zone_idx)) + mass_shift_grams)
      end do
      mass_scale_factor = (exp(ln10*log_mass(nz))-exp(ln10*log_mass(shell_end)))/ &
           (exp(ln10*log_mass(nz))-old_shell_end_mass)
      outer_mass = exp10(log_mass(nz))
      do zone_idx = shell_end+1, nz
         log_mass(zone_idx) = log10(outer_mass- &
              mass_scale_factor*(outer_mass-exp10(log_mass(zone_idx))))
      end do
end subroutine move_core_mass

! ---------------------------------------------------------------
! The envelope mixture from the envelope's X, Z and He3 mass fractions
! (x_env, z_env, he3_env) and the controls' relative mass fractions w_i
! of the ten metals (seed_weights, mixture_weights_seed, n_mix_species
! of them in ix_* slot order; H's and He's slots are not read), into
! number_fractions (n_mix_species):
!   y_env = 1 - x_env - z_env - he3_env;
!   with S the sum of w_i over the ten metals, Z outside C, N and O,
!     metals_no_cno = z_env (S - w_C - w_O - w_N)/S;
!   with f = z_env/S, the metals' mass fractions X_i = f w_i, H's x_env
!   and He's 1 - x_env - z_env (He3 with He4); inverse_atomic_weight is
!   sum_i X_i/A_i, 1/mu_ion of the mixture, and number fractions are
!   the (X_i/A_i)/inverse_atomic_weight (the EOS's, compute_blended_eos
!   and compute_saha_eos; DBG 1/96).
pure subroutine compute_envelope_mixture(x_env, z_env, he3_env, seed_weights, &
     y_env, metals_no_cno, inverse_atomic_weight, number_fractions)
      real(dp), intent(in) :: x_env, z_env, he3_env, seed_weights(:)
      real(dp), intent(out) :: y_env, metals_no_cno, inverse_atomic_weight, &
           number_fractions(:)
      real(dp) :: weights(n_mix_species), weight_sum, scale_factor
      integer :: i

      do i = 1, n_mix_species
         weights(i) = seed_weights(i)
      end do
      y_env = 1.0d0 - x_env - z_env - he3_env
      weight_sum = weights(ix_na)+weights(ix_al)+ &
           weights(ix_mg)+weights(ix_fe)+ &
           weights(ix_si)+weights(ix_c)+ &
           weights(ix_o)+weights(ix_n)+ &
           weights(ix_ar)+weights(ix_ne)
      metals_no_cno = z_env*(weight_sum - &
           weights(ix_c)-weights(ix_o)- &
           weights(ix_n))/weight_sum
      scale_factor = z_env/weight_sum
      weights(ix_h) = x_env/scale_factor
      weights(ix_he) = (1.0d0-x_env-z_env)/scale_factor
      weight_sum = 0.0d0
      do i = 1, n_mix_species
         weights(i) = scale_factor*weights(i)/mixture_atomic_weight(i)
         weight_sum = weight_sum + weights(i)
      end do
      inverse_atomic_weight = weight_sum
      scale_factor = 1.0d0/inverse_atomic_weight
      do i = 1, n_mix_species
         number_fractions(i) = weights(i)*scale_factor
      end do
end subroutine compute_envelope_mixture

! ---------------------------------------------------------------
! Sources
! 1. The rescalings are the F77 RSCALE's (git show 6cd56739:src/rscale.f);
!    the Z ramp is DBG 5/94. Its CNO factor was (Z_old - Z_new)/Z_old,
!    the fraction of Z removed: right only for a halved Z, zero for an
!    unchanged one, negative for a raised one; it is Z_new/Z_old since
!    2026-10.
! 2. compute_envelope_mixture is STARIN's AMUENV arithmetic, in its
!    order. The atomic weights are STARIN's ATOMWT
!    (git show 6cd56739:src/starin.f, line 118), the values of the Saha
!    EOS's ATOMW2 (eqstat.f). They are not the weights
!    step/shell_physics uses for the mean molecular weight (1.007825 and
!    4.002603 for H and He, against 1.008 and 4.004 here): merging the
!    two tables moves results.
end module rescaling_kernels
