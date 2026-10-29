!----------------------------------------------------------------------
! blended_eos
!----------------------------------------------------------------------
! The equation of state of every run that does not use MHD: the
! Yale EOS, with the SCV table at low temperature and an OPAL table
! blended in where it covers the point.
!     blend_eos            one evaluation at (log T, log P, X, Z):
!                          Saha (or SCV) below the ionization cutoff,
!                          fully ionized above it, interpolated across
!                          it, then the OPAL result adopted inside the
!                          table or blended in near its edge
!     compute_blended_eos  blend_eos, with the T and P derivatives of
!                          the gradients taken numerically (four more
!                          evaluations) when use_numerical_derivatives
!                          is set
module blended_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use yale_eos, only: compute_saha_eos, compute_fully_ionized_eos
      use scv_eos, only: compute_scv_envelope_eos
      use opal_eos, only: compute_opal_eos, find_opal_edge_ramp, opal_layout
      use physical_constants, only: gas_constant, ln10, &
           radiation_constant_over_3
      use io_units, only: run_log_unit
      use eos_mixture, only: eos_mix, ion_mean_weight_excess, ix_c, ix_h, &
           ix_he, ix_na, ix_ne, ix_o, n_mix_species
      use eos_configuration, only: eos_settings, opal_vintage_none
      use math_lib

      implicit none
      private
      public :: compute_blended_eos, blend_eos

! 1/A of H and He (the H and He excesses of 1/mu_ion, and 1/mu_e's H
! term), 1/A of C, and the electrons per nucleon of He and the metals,
! Z/A = 1/2 (0.4995 = 2/4.004)
      real(dp), parameter :: inverse_atomic_weights(4) = [0.9921d0, 0.24975d0, &
           0.08322d0, 0.4995d0]

contains

!----------------------------------------------------------------------
! compute_blended_eos (formerly eqstat)
!----------------------------------------------------------------------
! The blended EOS at (log T, log P, X, Z): blend_eos once. With
! want_derivatives and use_numerical_derivatives, the six second
! derivatives are instead taken by centred differences of blend_eos's
! first derivatives, four extra evaluations at log T +- dtl and
! log P +- dpl (LLP 10-22-06, 8-5-07):
!   dlnrho_dlnt_dt = d(dln rho/dln T)/dln T, dlnrho_dlnt_dp the same
!   in P; adiabatic_gradient_dt/_dp and specific_heat_cp_dt/_dp the log
!   derivatives of grad_ad and cp.
! The shifted evaluations start from a copy of r, carry their own
! density guess from one to the next (the fully ionized iteration's
! starting point), count their Saha iterations in a copy of
! saha_iterations, and leave r as it was.
subroutine compute_blended_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
     metal_fraction, r, want_derivatives, in_atmosphere, saha_iterations, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_iterations
      integer, intent(out) :: ierr
! the steps in log T and log P: a hundredth of the tables' approximate
! intervals (0.03, 0.15), which worked best empirically
      real(dp), parameter :: dtl = 0.030d0*0.01d0, dpl = 0.150d0*0.01d0
      real(dp), parameter :: dtl2 = 2.0d0*dtl, dpl2 = 2.0d0*dpl
! the shifted evaluations: their record, their density guess, and the
! first derivatives (dln rho/dln T, dln rho/dln P, cp, grad_ad) above
! and below
      type(eos_result) :: shifted
      real(dp) :: shifted_log10_density, shifted_density
      real(dp) :: above(4), below(4)
      real(dp) :: kept(6)
      integer :: saha_iterations_local

      ierr = 0
      if (want_derivatives .and. eos_settings%use_numerical_derivatives) then
         saha_iterations_local = saha_iterations
         shifted_log10_density = r%log10_density
         shifted_density = r%density
         call evaluate_shifted(log10_temperature + dtl, log10_pressure, above)
         if (ierr /= 0) return
         call evaluate_shifted(log10_temperature - dtl, log10_pressure, below)
         if (ierr /= 0) return
         r%dlnrho_dlnt_dt = (above(1) - below(1))/dtl2/ln10
         r%specific_heat_cp_dt = (log10(above(3)) - log10(below(3)))/dtl2
         r%adiabatic_gradient_dt = (log10(above(4)) - log10(below(4)))/dtl2
         call evaluate_shifted(log10_temperature, log10_pressure + dpl, above)
         if (ierr /= 0) return
         call evaluate_shifted(log10_temperature, log10_pressure - dpl, below)
         if (ierr /= 0) return
         r%dlnrho_dlnt_dp = (above(1) - below(1))/dpl2/ln10
         r%specific_heat_cp_dp = (log10(above(3)) - log10(below(3)))/dpl2
         r%adiabatic_gradient_dp = (log10(above(4)) - log10(below(4)))/dpl2
      end if

      if (want_derivatives .and. .not. eos_settings%use_numerical_derivatives) then
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, r, .true., in_atmosphere, saha_iterations, ierr)
      else
! no analytic second derivatives: r keeps the ones it has (the
! numerical ones, or the caller's)
         kept = [r%dlnrho_dlnt_dt, r%dlnrho_dlnt_dp, r%adiabatic_gradient_dt, &
              r%adiabatic_gradient_dp, r%specific_heat_cp_dt, r%specific_heat_cp_dp]
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, r, .false., in_atmosphere, saha_iterations, ierr)
         r%dlnrho_dlnt_dt = kept(1)
         r%dlnrho_dlnt_dp = kept(2)
         r%adiabatic_gradient_dt = kept(3)
         r%adiabatic_gradient_dp = kept(4)
         r%specific_heat_cp_dt = kept(5)
         r%specific_heat_cp_dp = kept(6)
      end if

contains

! One shifted evaluation (no second derivatives): this side's first
! derivatives, in first.
subroutine evaluate_shifted(lt, lp, first)
      real(dp), intent(in) :: lt, lp
      real(dp), intent(out) :: first(4)
      shifted = r
      shifted%log10_density = shifted_log10_density
      shifted%density = shifted_density
      call blend_eos(lt, lp, hydrogen_fraction, metal_fraction, shifted, &
           .false., in_atmosphere, saha_iterations_local, ierr)
      shifted_log10_density = shifted%log10_density
      shifted_density = shifted%density
      first = [shifted%dlnrho_dlnt, shifted%dlnrho_dlnp, shifted%specific_heat_cp, &
           shifted%adiabatic_gradient]
end subroutine evaluate_shifted

end subroutine compute_blended_eos

!----------------------------------------------------------------------
! blend_eos (formerly eqstat2)
!----------------------------------------------------------------------
! One evaluation of the EOS at (log T, log P, X, Z), in named steps:
!   zone_mixture           the zone's 1/mu_ion, and the Saha equations'
!                          abundances when they are solved
!   radiation              beta = P_gas/P (beta <= 0, radiation pressure
!                          at or above P, is an error)
!   the gas                above the ionization cutoff, a fully ionized
!                          gas (fully_ionized_gas); below it the Saha
!                          equations (compute_saha_eos), with SCV's
!                          density and first derivatives where SCV is on
!                          and covers the point, and within
!                          saha_ramp_width of the cutoff blended into the
!                          fully ionized solution (blend_saha_into)
!   OPAL                   inside the selected vintage's (T, rho) window,
!                          its result where its table covers the point
!                          (adopt_table_result), blended in near the
!                          table's edge (blend_table_result)
! metal_fraction is intent(in): no callee writes it.
subroutine blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
     metal_fraction, r, want_derivatives, in_atmosphere, saha_iterations, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_iterations
      integer, intent(out) :: ierr
! the Saha ramp: its width below the cutoff in log T, and the scale that
! maps it onto a weight from 1 (its foot) to 0 (the cutoff)
      real(dp), parameter :: saha_ramp_width = 0.500d0, saha_ramp_scale = 1.0d0/saha_ramp_width
      logical :: need_saha_solution, skip_fully_ionized, valid_table_point
      real(dp) :: saha_mass_fractions(n_mix_species), radiation_fraction
! the Saha (or SCV) result, blended into the fully ionized one
      type(eos_result) :: saha
! the OPAL result, adopted or blended in
      type(eos_result) :: table_r
      real(dp) :: ramp_factor
      logical :: in_opal_table, needs_ramp, opal_covers_point
      integer :: vintage

      ierr = 0
! the Saha equations are solved below the cutoff, alone below its ramp
      need_saha_solution = log10_temperature < eos_settings%saha_log10t_cutoff
      skip_fully_ionized = log10_temperature <= &
           (eos_settings%saha_log10t_cutoff - saha_ramp_width)
      call zone_mixture(hydrogen_fraction, metal_fraction, need_saha_solution, &
           r%ion_mean_weight_inverse, saha_mass_fractions)

! radiation: its share of P, beta = P_gas/P, and beta14 = 4(1 - beta)/beta
      r%temperature = exp(ln10*log10_temperature)
      r%pressure = exp(ln10*log10_pressure)
      radiation_fraction = radiation_constant_over_3*(r%temperature**2)**2/r%pressure
      r%beta = 1 - radiation_fraction
      if (r%beta <= 0d0) then
         write(run_log_unit,*) &
              'blend_eos: beta = P_gas/P is not positive: log T, log P, beta', &
              log10_temperature, log10_pressure, r%beta
         ierr = 1
         return
      end if
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 4.0d0*r%beta_inverse*min(1.0d0, radiation_fraction)
      r%electron_degeneracy_parameter = 0.0d0

! the gas
      if (.not. need_saha_solution) then
         call fully_ionized_gas(log10_temperature, log10_pressure, hydrogen_fraction, r)
      else
         call compute_saha_eos(saha_mass_fractions, log10_temperature, &
              log10_pressure, r, want_derivatives, in_atmosphere, saha_iterations)
! SCV replaces Saha's density and first derivatives (Saha's second
! derivatives stay) where its table covers the point; elsewhere Saha's
! result stands
         if (eos_settings%use_scv_eos) then
            call compute_scv_envelope_eos(log10_temperature, hydrogen_fraction, &
                 metal_fraction, r, valid_table_point, ierr)
            if (ierr /= 0) return
         end if
         if (.not. skip_fully_ionized) then
            saha = r
            call fully_ionized_gas(log10_temperature, log10_pressure, hydrogen_fraction, r)
            call blend_saha_into(saha, &
                 saha_ramp_scale*(eos_settings%saha_log10t_cutoff - log10_temperature), &
                 want_derivatives, r)
         end if
      end if

! OPAL: inside the selected vintage's (T, rho) window, its result
! replaces the Yale/SCV one, or is blended in with weight ramp_factor
! near the table's edge; a point the table does not cover keeps the
! result so far
      vintage = eos_settings%opal_vintage
      if (vintage == opal_vintage_none) return
      if (.not. inside_opal_window(vintage, r%temperature, log10_temperature, &
           r%log10_density)) return
      call compute_opal_eos(vintage, log10_temperature, log10_pressure, &
           hydrogen_fraction, metal_fraction, table_r, ierr, opal_covers_point)
      r%temperature = table_r%temperature
      r%pressure = table_r%pressure
      if (ierr /= 0 .or. .not. opal_covers_point) return
      call find_opal_edge_ramp(vintage, r%temperature, table_r%log10_density, &
           ramp_factor, in_opal_table, needs_ramp, ierr)
      if (ierr /= 0) return
      if (.not. in_opal_table) return
      if (.not. needs_ramp) then
         call adopt_table_result(table_r, r)
      else
         call blend_table_result(table_r, ramp_factor, r)
      end if
end subroutine blend_eos

!----------------------------------------------------------------------
! zone_mixture
!----------------------------------------------------------------------
! The zone's 1/mu_ion and, when the Saha equations are solved, the
! abundances they take (X_i/A_i over 1/mu_ion, the 12 species in the
! ix_* order), from the envelope mixture eos_mix (fxenv = X_i/A_i/amuenv,
! set by eos_set_mixture). With use_diffusion_z all the metals scale
! with Z/Z_env (MHP 3/94), and H and He add their excess over the
! envelope; otherwise the zone differs from the envelope by its excess
! H, He and metals, the metals counted as carbon (ion_mean_weight_excess,
! which compute_mean_weights shares).
subroutine zone_mixture(hydrogen_fraction, metal_fraction, need_saha_solution, &
     ion_mean_weight_inverse, saha_mass_fractions)
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      logical, intent(in) :: need_saha_solution
      real(dp), intent(out) :: ion_mean_weight_inverse
      real(dp), intent(inout) :: saha_mass_fractions(n_mix_species)
      integer :: species_idx
      real(dp) :: metal_ratio, amu_correction, h_excess, y_excess
      real(dp) :: envelope_amu_over_amu
      real(dp) :: dfx1, dfx12, dfx4, amu_inverse, envelope_amu_frac
      logical :: use_envelope_mixture

      if (eos_settings%use_diffusion_z) then
         metal_ratio = metal_fraction/eos_mix%envelope_metal_fraction
         amu_correction = (metal_ratio - 1.0d0)*eos_mix%amuenv
         ion_mean_weight_inverse = eos_mix%amuenv
         do species_idx = ix_na, ix_c
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 amu_correction*eos_mix%fxenv(species_idx)
         end do
         h_excess = (hydrogen_fraction - eos_mix%envelope_hydrogen_fraction)* &
              inverse_atomic_weights(1)
         ion_mean_weight_inverse = ion_mean_weight_inverse + h_excess
         do species_idx = ix_o, ix_ne
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 amu_correction*eos_mix%fxenv(species_idx)
         end do
         y_excess = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              hydrogen_fraction - metal_fraction)*inverse_atomic_weights(2)
         ion_mean_weight_inverse = ion_mean_weight_inverse + y_excess
         if (.not. need_saha_solution) return
         envelope_amu_over_amu = metal_ratio*eos_mix%amuenv/ion_mean_weight_inverse
         do species_idx = ix_na, ix_c
            saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                 eos_mix%fxenv(species_idx)
         end do
         saha_mass_fractions(ix_h) = (eos_mix%fxenv(ix_h)* &
              eos_mix%amuenv + h_excess)/ion_mean_weight_inverse
         do species_idx = ix_o, ix_ne
            saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                 eos_mix%fxenv(species_idx)
         end do
         saha_mass_fractions(ix_he) = (eos_mix%fxenv(ix_he)* &
              eos_mix%amuenv + y_excess)/ion_mean_weight_inverse
      else
         call ion_mean_weight_excess(hydrogen_fraction, metal_fraction, &
              ion_mean_weight_inverse, dfx1, dfx12, dfx4, use_envelope_mixture)
         if (.not. need_saha_solution) return
         if (use_envelope_mixture) then
            do species_idx = 1, n_mix_species
               saha_mass_fractions(species_idx) = eos_mix%fxenv(species_idx)
            end do
         else
            amu_inverse = 1.0d0/ion_mean_weight_inverse
            envelope_amu_frac = eos_mix%amuenv*amu_inverse
            do species_idx = 1, n_mix_species
               saha_mass_fractions(species_idx) = envelope_amu_frac* &
                    eos_mix%fxenv(species_idx)
            end do
            saha_mass_fractions(ix_c) = saha_mass_fractions(ix_c) + &
                 dfx12*amu_inverse
            saha_mass_fractions(ix_h) = saha_mass_fractions(ix_h) + &
                 dfx1*amu_inverse
            saha_mass_fractions(ix_he) = saha_mass_fractions(ix_he) + &
                 dfx4*amu_inverse
         end if
      end if
end subroutine zone_mixture

!----------------------------------------------------------------------
! fully_ionized_gas
!----------------------------------------------------------------------
! A fully ionized gas: 1/mu_e = X/A_H + (1 - X)(Z/A of He and the
! metals, 1/2), R/mu, the ionization fractions (H+, He+, He++ = 1, 0, 1),
! then compute_fully_ionized_eos's density and derivatives.
subroutine fully_ionized_gas(log10_temperature, log10_pressure, hydrogen_fraction, r)
      real(dp), intent(in) :: log10_temperature, log10_pressure, hydrogen_fraction
      type(eos_result), intent(inout) :: r
      r%electron_mean_weight_inverse = hydrogen_fraction*inverse_atomic_weights(1) &
           + (1.0d0 - hydrogen_fraction)*inverse_atomic_weights(4)
      r%specific_gas_constant = gas_constant*(r%ion_mean_weight_inverse + &
           r%electron_mean_weight_inverse)
      r%ion_fraction(1) = 1.0d0
      r%ion_fraction(2) = 0.0d0
      r%ion_fraction(3) = 1.0d0
      call compute_fully_ionized_eos(log10_temperature, log10_pressure, r)
end subroutine fully_ionized_gas

!----------------------------------------------------------------------
! blend_saha_into
!----------------------------------------------------------------------
! The ramp from the Saha (or SCV) result to the fully ionized one r: with
! ramp_weight w from 1 at the ramp's foot to 0 at the cutoff, the
! smoothstep f = w^3 (6w^2 - 15w + 10) (its first and second derivatives
! vanish at both ends) takes r = r + f (saha - r) for log rho (then rho,
! and R/mu and 1/mu_e from them), the ionization fractions, dln rho/dln T,
! dln rho/dln P, cp, grad_ad and, with want_derivatives, the six second
! derivatives.
subroutine blend_saha_into(saha, ramp_weight, want_derivatives, r)
      type(eos_result), intent(in) :: saha
      real(dp), intent(in) :: ramp_weight
      logical, intent(in) :: want_derivatives
      type(eos_result), intent(inout) :: r
      real(dp) :: ramp_weight_sq, f
      ramp_weight_sq = ramp_weight*ramp_weight
      f = ramp_weight_sq*ramp_weight*(6.0d0*ramp_weight_sq - 15.0d0*ramp_weight + 10.0d0)
      r%log10_density = r%log10_density + f*(saha%log10_density - r%log10_density)
      r%density = exp(ln10*r%log10_density)
      r%specific_gas_constant = r%beta*r%pressure/(r%density*r%temperature)
      r%electron_mean_weight_inverse = r%specific_gas_constant/gas_constant - &
           r%ion_mean_weight_inverse
      r%ion_fraction = r%ion_fraction + f*(saha%ion_fraction - r%ion_fraction)
      r%dlnrho_dlnt = r%dlnrho_dlnt + f*(saha%dlnrho_dlnt - r%dlnrho_dlnt)
      r%dlnrho_dlnp = r%dlnrho_dlnp + f*(saha%dlnrho_dlnp - r%dlnrho_dlnp)
      r%specific_heat_cp = r%specific_heat_cp + f*(saha%specific_heat_cp - r%specific_heat_cp)
      r%adiabatic_gradient = r%adiabatic_gradient + &
           f*(saha%adiabatic_gradient - r%adiabatic_gradient)
      if (.not. want_derivatives) return
      r%dlnrho_dlnt_dt = r%dlnrho_dlnt_dt + f*(saha%dlnrho_dlnt_dt - r%dlnrho_dlnt_dt)
      r%dlnrho_dlnt_dp = r%dlnrho_dlnt_dp + f*(saha%dlnrho_dlnt_dp - r%dlnrho_dlnt_dp)
      r%specific_heat_cp_dt = r%specific_heat_cp_dt + &
           f*(saha%specific_heat_cp_dt - r%specific_heat_cp_dt)
      r%specific_heat_cp_dp = r%specific_heat_cp_dp + &
           f*(saha%specific_heat_cp_dp - r%specific_heat_cp_dp)
      r%adiabatic_gradient_dt = r%adiabatic_gradient_dt + &
           f*(saha%adiabatic_gradient_dt - r%adiabatic_gradient_dt)
      r%adiabatic_gradient_dp = r%adiabatic_gradient_dp + &
           f*(saha%adiabatic_gradient_dp - r%adiabatic_gradient_dp)
end subroutine blend_saha_into

! the (T, rho) window inside which a vintage is consulted at all (the
! table's edge is then find_opal_edge_ramp's): its row of opal_layout,
! where 1995 bounds log10 T and the others T
logical function inside_opal_window(vintage, temperature, log10_temperature, log10_density)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: temperature, log10_temperature, log10_density
      inside_opal_window = temperature >= opal_layout(vintage)%window_t_min .and. &
           temperature <= opal_layout(vintage)%window_t_max .and. &
           log10_temperature <= opal_layout(vintage)%window_log10_t_max .and. &
           log10_density <= opal_layout(vintage)%window_log10_rho_max
end function inside_opal_window

! Inside a table source's coverage (OPAL's): its density, beta, mean
! weights and first derivatives replace r's (the second derivatives stay
! the Yale/SCV ones).
subroutine adopt_table_result(table_r, r)
      type(eos_result), intent(in) :: table_r
      type(eos_result), intent(inout) :: r
      r%log10_density = table_r%log10_density
      r%density = exp10(r%log10_density)
      r%beta = table_r%beta
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 1.0d0 - r%beta
      r%specific_gas_constant = table_r%specific_gas_constant
      r%ion_mean_weight_inverse = table_r%ion_mean_weight_inverse
      r%electron_mean_weight_inverse = table_r%electron_mean_weight_inverse
      r%dlnrho_dlnt = table_r%dlnrho_dlnt
      r%dlnrho_dlnp = table_r%dlnrho_dlnp
      r%specific_heat_cp = table_r%specific_heat_cp
      r%adiabatic_gradient = table_r%adiabatic_gradient
end subroutine adopt_table_result

! Near a table source's edge: the same quantities blended with its edge
! weight, r + w (table - r).
subroutine blend_table_result(table_r, ramp_factor, r)
      type(eos_result), intent(in) :: table_r
      real(dp), intent(in) :: ramp_factor
      type(eos_result), intent(inout) :: r
      r%log10_density = r%log10_density + &
           ramp_factor*(table_r%log10_density - r%log10_density)
      r%density = exp10(r%log10_density)
      r%beta = r%beta + ramp_factor*(table_r%beta - r%beta)
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 1.0d0 - r%beta
      r%specific_gas_constant = r%specific_gas_constant + &
           ramp_factor*(table_r%specific_gas_constant - r%specific_gas_constant)
      r%ion_mean_weight_inverse = r%ion_mean_weight_inverse + &
           ramp_factor*(table_r%ion_mean_weight_inverse - r%ion_mean_weight_inverse)
      r%electron_mean_weight_inverse = r%electron_mean_weight_inverse + &
           ramp_factor*(table_r%electron_mean_weight_inverse - r%electron_mean_weight_inverse)
      r%dlnrho_dlnt = r%dlnrho_dlnt + &
           ramp_factor*(table_r%dlnrho_dlnt - r%dlnrho_dlnt)
      r%dlnrho_dlnp = r%dlnrho_dlnp + &
           ramp_factor*(table_r%dlnrho_dlnp - r%dlnrho_dlnp)
      r%specific_heat_cp = r%specific_heat_cp + &
           ramp_factor*(table_r%specific_heat_cp - r%specific_heat_cp)
      r%adiabatic_gradient = r%adiabatic_gradient + &
           ramp_factor*(table_r%adiabatic_gradient - r%adiabatic_gradient)
end subroutine blend_table_result

end module blended_eos
