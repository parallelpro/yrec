!----------------------------------------------------------------------
! blended_eos
!----------------------------------------------------------------------
! The equation of state of every run that does not use a whole-EOS
! source (MHD): the Yale EOS, with the table sources of eos_source_list
! (SCV at low temperature, an OPAL table) applied where they cover the
! point.
!     blend_eos            one evaluation at (log T, log P, X, Z):
!                          Saha below the ionization cutoff (SCV in its
!                          place where it covers the point), fully
!                          ionized above it, interpolated across it,
!                          then the OPAL result adopted inside the table
!                          or blended in near its edge
!     compute_blended_eos  blend_eos, with the T and P derivatives of
!                          the gradients taken numerically (four more
!                          evaluations) when use_numerical_derivatives
!                          is set
module blended_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use yale_eos, only: compute_saha_eos, compute_fully_ionized_eos
      use yale_eos_tables, only: debye_huckel_composition
      use eos_table_source, only: eos_table_source_t, stage_in_place_of_saha, stage_over_the_gas
      use eos_source_list, only: table_sources, n_table_sources
      use opal_eos, only: esac06
      use opal_eos_tables, only: i_opal_gamma1, i_opal_gamma2_ratio, opal_tables, iv_opal06
      use physical_constants, only: gas_constant, ln10, &
           radiation_constant_over_3
      use eos_messages, only: eos_log, eos_report, eos_log_once, once_fermi_extrapolation, &
           message_length, unwritten_line
      use eos_mixture, only: eos_mix, ion_mean_weight_excess, mean_weights_from_density, &
           debye_huckel_composition_of, &
           inverse_a_hydrogen, inverse_a_helium, electrons_per_nucleon_metals, &
           ix_c, ix_h, ix_he, ix_na, ix_ne, ix_o, n_mix_species
      use eos_configuration, only: eos_settings, opal_vintage_2006, saha_ramp_width
      use math_lib

      implicit none
      private
      public :: compute_blended_eos, blend_eos, gamma1_at_density


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
     metal_fraction, dh, r, want_derivatives, in_atmosphere, saha_iterations, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(debye_huckel_composition), intent(in) :: dh
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
              metal_fraction, dh, r, .true., in_atmosphere, saha_iterations, ierr)
      else
! no analytic second derivatives: r keeps the ones it has (the
! numerical ones, or the caller's)
         kept = [r%dlnrho_dlnt_dt, r%dlnrho_dlnt_dp, r%adiabatic_gradient_dt, &
              r%adiabatic_gradient_dp, r%specific_heat_cp_dt, r%specific_heat_cp_dp]
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, dh, r, .false., in_atmosphere, saha_iterations, ierr)
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
      call blend_eos(lt, lp, hydrogen_fraction, metal_fraction, dh, shifted, &
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
! One evaluation of the EOS at (log T, log P, X, Z), in named steps,
! each taking the result so far:
!   zone_mixture           the zone's 1/mu_ion, and the Saha equations'
!                          abundances when they are solved
!   radiation_share        beta = P_gas/P (beta <= 0, radiation pressure
!                          at or above P, is an error)
!   gas_eos                above the ionization cutoff, a fully ionized
!                          gas; below it the Saha equations, with the
!                          selected sources of stage_in_place_of_saha
!                          (SCV) applied to them, and within
!                          saha_ramp_width of the cutoff that blended
!                          into the fully ionized one
!   the sources of stage_over_the_gas (OPAL), each applied in the list's
!                          order (apply_table_source): its result where
!                          its table covers the point, blended in near
!                          the table's edge
! A table source added later is one more entry of eos_source_list.
subroutine blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
     metal_fraction, dh, r, want_derivatives, in_atmosphere, saha_iterations, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
! the zone's composition for the Debye-Hueckel correction (eos_get's)
      type(debye_huckel_composition), intent(in) :: dh
      type(eos_result), intent(inout) :: r
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_iterations
      integer, intent(out) :: ierr
      real(dp) :: saha_number_fractions(n_mix_species)
      type(eos_table_source_t) :: sources(n_table_sources)
      integer :: k

      ierr = 0
      call zone_mixture(hydrogen_fraction, metal_fraction, &
           log10_temperature < eos_settings%saha_log10t_cutoff, &
           r%ion_mean_weight_inverse, saha_number_fractions)
      call radiation_share(log10_temperature, log10_pressure, r, ierr)
      if (ierr /= 0) return
      call gas_eos(log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction, &
           saha_number_fractions, dh, r, want_derivatives, in_atmosphere, saha_iterations, ierr)
      if (ierr /= 0) return
      sources = table_sources()
      do k = 1, n_table_sources
         if (sources(k)%stage /= stage_over_the_gas) cycle
         if (.not. sources(k)%selected()) cycle
         call apply_table_source(sources(k), log10_temperature, log10_pressure, &
              hydrogen_fraction, metal_fraction, r, ierr)
         if (ierr /= 0) return
      end do
end subroutine blend_eos

!----------------------------------------------------------------------
! radiation_share
!----------------------------------------------------------------------
! T and P from their logs, the radiation pressure's share of P, beta =
! P_gas/P and beta14 = 4(1 - beta)/beta; beta <= 0 (radiation pressure
! at or above P: no gas pressure to solve for) is ierr = 1.
subroutine radiation_share(log10_temperature, log10_pressure, r, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(eos_result), intent(inout) :: r
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)
      real(dp) :: radiation_fraction

      ierr = 0
      r%temperature = exp(ln10*log10_temperature)
      r%pressure = exp(ln10*log10_pressure)
      radiation_fraction = radiation_constant_over_3*(r%temperature**2)**2/r%pressure
      r%beta = 1 - radiation_fraction
      if (r%beta <= 0d0) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'blend_eos: beta = P_gas/P is not positive: log T, log P, beta', &
              log10_temperature, log10_pressure, r%beta
         call eos_report(eos_lines)
         ierr = 1
         return
      end if
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 4.0d0*r%beta_inverse*radiation_fraction
      r%electron_degeneracy_parameter = 0.0d0
end subroutine radiation_share

!----------------------------------------------------------------------
! gas_eos
!----------------------------------------------------------------------
! The gas, by log T against the Saha cutoff T_c: at and above it a
! fully ionized gas (fully_ionized_gas); below it the Saha equations
! (compute_saha_eos), and the selected sources of stage_in_place_of_saha
! applied to that result (apply_table_source; SCV's lookup refines a copy
! of it: its density, the mean weights that density implies, the
! ionization fractions and the first derivatives, the Saha second
! derivatives staying; outside the table the Saha result stands);
! within saha_ramp_width of T_c that blended into the fully ionized gas
! (blend_saha_into).
subroutine gas_eos(log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction, &
     saha_number_fractions, dh, r, want_derivatives, in_atmosphere, saha_iterations, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(in) :: saha_number_fractions(:)
      type(debye_huckel_composition), intent(in) :: dh
      type(eos_result), intent(inout) :: r
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_iterations
      integer, intent(out) :: ierr
! the scale that maps the Saha ramp (saha_ramp_width below the cutoff)
! onto a weight from 1 (its foot) to 0 (the cutoff)
      real(dp), parameter :: saha_ramp_scale = 1.0d0/saha_ramp_width
! the Saha (or SCV) result, blended into the fully ionized one
      type(eos_result) :: saha
      type(eos_table_source_t) :: sources(n_table_sources)
      integer :: k

      ierr = 0
      if (log10_temperature >= eos_settings%saha_log10t_cutoff) then
         call fully_ionized_gas(log10_temperature, log10_pressure, hydrogen_fraction, dh, r, ierr)
         return
      end if
      call compute_saha_eos(saha_number_fractions, log10_temperature, &
           log10_pressure, r, want_derivatives, in_atmosphere, saha_iterations)
      sources = table_sources()
      do k = 1, n_table_sources
         if (sources(k)%stage /= stage_in_place_of_saha) cycle
         if (.not. sources(k)%selected()) cycle
         call apply_table_source(sources(k), log10_temperature, log10_pressure, &
              hydrogen_fraction, metal_fraction, r, ierr)
         if (ierr /= 0) return
      end do
      if (log10_temperature <= eos_settings%saha_log10t_cutoff - saha_ramp_width) return
      saha = r
      call fully_ionized_gas(log10_temperature, log10_pressure, hydrogen_fraction, dh, r, ierr)
      if (ierr /= 0) return
      call blend_saha_into(saha, &
           saha_ramp_scale*(eos_settings%saha_log10t_cutoff - log10_temperature), &
           want_derivatives, r)
end subroutine gas_eos

!----------------------------------------------------------------------
! apply_table_source
!----------------------------------------------------------------------
! One table source at the point (eos_table_source): inside its window
! (in_window), its lookup from a copy of r; T and P become the lookup's
! (OPAL's come from its own exp10; SCV's are r's); where its table covers
! the point its result replaces r's (adopt) or, near the table's edge,
! is blended in with its edge weight (edge_weight, blend); elsewhere r
! stands. ierr is the lookup's or the edge weight's.
subroutine apply_table_source(src, log10_temperature, log10_pressure, hydrogen_fraction, &
     metal_fraction, r, ierr)
      type(eos_table_source_t), intent(in) :: src
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      integer, intent(out) :: ierr
! the source's result, adopted or blended in, and its edge weight
      type(eos_result) :: table_r
      real(dp) :: weight
      logical :: covers, in_table, near_edge

      ierr = 0
      if (.not. src%in_window(log10_temperature, r)) return
      table_r = r
      call src%lookup(log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction, &
           table_r, covers, ierr)
      r%temperature = table_r%temperature
      r%pressure = table_r%pressure
      if (ierr /= 0 .or. .not. covers) return
      call src%edge_weight(r%temperature, table_r%log10_density, weight, in_table, near_edge, ierr)
      if (ierr /= 0) return
      if (.not. in_table) return
      if (.not. near_edge) then
         call src%adopt(table_r, r)
      else
         call src%blend(table_r, weight, r)
      end if
end subroutine apply_table_source

!----------------------------------------------------------------------
! zone_mixture
!----------------------------------------------------------------------
! The zone's 1/mu_ion and, when the Saha equations are solved, the
! abundances they take (X_i/A_i over 1/mu_ion, the 12 species in the
! ix_* order), from the envelope mixture eos_mix (number_fractions =
! X_i/A_i over inverse_atomic_weight = sum_j X_j/A_j, set by
! eos_set_mixture). With use_diffusion_z all the metals scale
! with Z/Z_env (MHP 3/94), and H and He add their excess over the
! envelope; otherwise, and when the envelope has no metals to scale
! (Z_env = 0), the zone differs from the envelope by its excess H, He
! and metals, the metals counted as carbon (ion_mean_weight_excess,
! which compute_mean_weights shares).
subroutine zone_mixture(hydrogen_fraction, metal_fraction, need_saha_solution, &
     ion_mean_weight_inverse, saha_number_fractions)
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      logical, intent(in) :: need_saha_solution
      real(dp), intent(out) :: ion_mean_weight_inverse
      real(dp), intent(inout) :: saha_number_fractions(n_mix_species)
      integer :: species_idx
      real(dp) :: metal_ratio, metal_weight_excess, h_excess, y_excess
      real(dp) :: metal_number_scale
      real(dp) :: dfx1, dfx12, dfx4, ion_mean_weight, envelope_number_scale
      logical :: use_envelope_mixture

      if (eos_settings%use_diffusion_z .and. eos_mix%envelope_metal_fraction > 0.0d0) then
         metal_ratio = metal_fraction/eos_mix%envelope_metal_fraction
         metal_weight_excess = (metal_ratio - 1.0d0)*eos_mix%inverse_atomic_weight
         ion_mean_weight_inverse = eos_mix%inverse_atomic_weight
         do species_idx = ix_na, ix_c
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 metal_weight_excess*eos_mix%number_fractions(species_idx)
         end do
         h_excess = (hydrogen_fraction - eos_mix%envelope_hydrogen_fraction)* &
              inverse_a_hydrogen
         ion_mean_weight_inverse = ion_mean_weight_inverse + h_excess
         do species_idx = ix_o, ix_ne
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 metal_weight_excess*eos_mix%number_fractions(species_idx)
         end do
         y_excess = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              hydrogen_fraction - metal_fraction)*inverse_a_helium
         ion_mean_weight_inverse = ion_mean_weight_inverse + y_excess
         if (.not. need_saha_solution) return
         metal_number_scale = metal_ratio*eos_mix%inverse_atomic_weight/ion_mean_weight_inverse
         do species_idx = ix_na, ix_c
            saha_number_fractions(species_idx) = metal_number_scale* &
                 eos_mix%number_fractions(species_idx)
         end do
         saha_number_fractions(ix_h) = (eos_mix%number_fractions(ix_h)* &
              eos_mix%inverse_atomic_weight + h_excess)/ion_mean_weight_inverse
         do species_idx = ix_o, ix_ne
            saha_number_fractions(species_idx) = metal_number_scale* &
                 eos_mix%number_fractions(species_idx)
         end do
         saha_number_fractions(ix_he) = (eos_mix%number_fractions(ix_he)* &
              eos_mix%inverse_atomic_weight + y_excess)/ion_mean_weight_inverse
      else
         call ion_mean_weight_excess(hydrogen_fraction, metal_fraction, &
              ion_mean_weight_inverse, dfx1, dfx12, dfx4, use_envelope_mixture)
         if (.not. need_saha_solution) return
         if (use_envelope_mixture) then
            do species_idx = 1, n_mix_species
               saha_number_fractions(species_idx) = eos_mix%number_fractions(species_idx)
            end do
         else
            ion_mean_weight = 1.0d0/ion_mean_weight_inverse
            envelope_number_scale = eos_mix%inverse_atomic_weight*ion_mean_weight
            do species_idx = 1, n_mix_species
               saha_number_fractions(species_idx) = envelope_number_scale* &
                    eos_mix%number_fractions(species_idx)
            end do
            saha_number_fractions(ix_c) = saha_number_fractions(ix_c) + &
                 dfx12*ion_mean_weight
            saha_number_fractions(ix_h) = saha_number_fractions(ix_h) + &
                 dfx1*ion_mean_weight
            saha_number_fractions(ix_he) = saha_number_fractions(ix_he) + &
                 dfx4*ion_mean_weight
         end if
      end if
end subroutine zone_mixture

!----------------------------------------------------------------------
! fully_ionized_gas
!----------------------------------------------------------------------
! A fully ionized gas: 1/mu_e = X/A_H + (1 - X)(Z/A of He and the
! metals, 1/2), R/mu, the ionization fractions (H+, He+, He++ = 1, 0, 1),
! then compute_fully_ionized_eos's density and derivatives (ierr = 1 when
! its density iteration does not converge). The first point of a run
! where they extrapolate the Fermi-Dirac table gets a warning in the run
! log (eos_log_once).
subroutine fully_ionized_gas(log10_temperature, log10_pressure, hydrogen_fraction, dh, r, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure, hydrogen_fraction
      type(debye_huckel_composition), intent(in) :: dh
      type(eos_result), intent(inout) :: r
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)
      logical :: extrapolated
      r%electron_mean_weight_inverse = hydrogen_fraction*inverse_a_hydrogen &
           + (1.0d0 - hydrogen_fraction)*electrons_per_nucleon_metals
      r%specific_gas_constant = gas_constant*(r%ion_mean_weight_inverse + &
           r%electron_mean_weight_inverse)
      r%ion_fraction(1) = 1.0d0
      r%ion_fraction(2) = 0.0d0
      r%ion_fraction(3) = 1.0d0
      call compute_fully_ionized_eos(log10_temperature, log10_pressure, dh, r, extrapolated, ierr)
      if (ierr /= 0 .or. .not. extrapolated) return
      eos_lines = unwritten_line
      write(eos_lines,'(a/a,f0.4,a,f0.4,a)') 'eos: the fully ionized gas extrapolates the Fermi-Dirac ' // &
           'table (log T 5.2 to 9, log rho - log mu_e - 1.5 log T from -16 to -3)', &
           '     at log T ', log10_temperature, ', log P ', log10_pressure, &
           '; the run''s later points outside it are not reported'
      call eos_log_once(once_fermi_extrapolation, eos_lines)
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
      call mean_weights_from_density(r)
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

!----------------------------------------------------------------------
! gamma1_at_density (formerly eos_lib's gamma1_at_density)
!----------------------------------------------------------------------
! gamma1 and the adiabatic gradient at (T6, rho, P) and (X, Z): from
! the OPAL 2006 table (esac06, with the radiation correction) when
! the 2006 vintage is selected and the table covers the point, else from
! blend_eos, with gamma1 built from chi_rho, chi_T and cv (Cox & Giuli
! ch. 9). test_eos is its only caller (eos_lib does not export it).
subroutine gamma1_at_density(hydrogen_fraction, metal_fraction, &
     temperature_1e6k, density, pressure, gamma1, adiabatic_gradient, &
     saha_iterations, ierr)

      real(dp), intent(in) :: hydrogen_fraction, metal_fraction, &
           temperature_1e6k, density, pressure
      real(dp), intent(out) :: gamma1, adiabatic_gradient
      integer, intent(inout) :: saha_iterations
! ierr as for eos_get; an out-of-table esac06 lookup is not an error.
      integer, intent(out) :: ierr
! esac06's full evaluation (its nine variables) with the radiation
! correction
      integer, parameter :: opal06_variables = 9, with_radiation = 1

      integer :: jerr
      real(dp) :: log10_temperature, log10_pressure
      type(eos_result) :: g
      type(debye_huckel_composition) :: dh
      real(dp) :: chi_rho, chi_t, specific_heat_cv
      logical :: in_table

      ierr = 0
      jerr = 0
      in_table = .false.

      if (eos_settings%opal_vintage == opal_vintage_2006) then
         call esac06(opal_tables(iv_opal06), hydrogen_fraction, temperature_1e6k, density, &
              opal06_variables, with_radiation, jerr, in_table)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      end if
      if (eos_settings%opal_vintage == opal_vintage_2006 .and. in_table) then
         gamma1 = opal_tables(iv_opal06)%eos_output(i_opal_gamma1)
         adiabatic_gradient = 1.0d0/opal_tables(iv_opal06)%eos_output(i_opal_gamma2_ratio)
      else
! (rho is the density guess; P and T are set from their logs)
         g%density = density
         g%log10_density = log10(density)
         log10_temperature = log10(temperature_1e6k*1.0d6)
         log10_pressure = log10(pressure)
         if (eos_settings%use_debye_huckel_correction) &
              dh = debye_huckel_composition_of(hydrogen_fraction, metal_fraction)
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, dh, g, .false., .false., saha_iterations, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         chi_rho = 1.0d0/g%dlnrho_dlnp
         chi_t = -chi_rho*g%dlnrho_dlnt
         specific_heat_cv = g%specific_heat_cp - exp(ln10*(log10_pressure - &
              g%log10_density - log10_temperature))*chi_t**2/chi_rho
         gamma1 = chi_rho*g%specific_heat_cp/specific_heat_cv
         adiabatic_gradient = g%adiabatic_gradient
      end if

end subroutine gamma1_at_density

end module blended_eos
