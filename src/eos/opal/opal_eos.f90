!----------------------------------------------------------------------
! opal_eos
!----------------------------------------------------------------------
! The OPAL equation of state of the three vintages, 1995, 2001 and 2006,
! which share one state layout (opal_eos_vintage, their instances
! opal_tables(iv_opal95:iv_opal06)).
!     compute_opal_eos     density and the thermodynamic derivatives of
!                          one vintage at (log T, log P, X, Z)
!     find_opal_edge_ramp  whether the vintage's table holds the point,
!                          and the OPAL weight near its edge
! Under these: opal_density_search (OPAL's rhoofp, rhoofp01 and
! rhoofp06), interpolate_opal_eos (OPAL's esac, esac01 and esac06, whose
! names are its entry points in eos_test_entries), t6rinteos, the radiation corrections (add_radiation;
! 1995: add_radiation_to_gammas), readcoeos (the table reader, run once
! per vintage), and the edge ramps find_table_edge_ramp (1995, 2001) and
! find_ragged_edge_ramp (2006). This code is written once: what the
! vintages differ in (the table's size and grids, the slots, the
! radiation constants and correction, the mass per particle, the
! reader's formats, how the density search
! treats radiation, which edge the table has, the name its messages
! start with) is their row of opal_layout (opal_eos_tables); the
! report_* procedures write each message once. The 1995 code follows Rogers,
! Swenson & Iglesias (1996), the 2001 code Rogers & Nayfonov (2002); the
! 2006 tables are the 2005 update.
module opal_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use opal_table_reader, only: open_opal_eos_table, readcoeos
      use opal_eos_tables, only: find_opal_ramp, i_opal_p, &
           interpolate_opal_t6_rho, iv_opal95, iv_opal01, iv_opal06, &
           n_eos_mv, n_eos_mx, opal_tables, opal_layout, opal_shared, &
           opal_eos_vintage, opal_flag_set, opal_rho_not_found, opal_vintage_layout, quad, &
           opal_radiation_constant_over_3, opal_empty_entry, opal_empty_sum, &
           opal_usable_corner_sum, opal_beyond_stencil, opal_x_interp_overflow, opal_range_overflow
      use physical_constants, only: gas_constant
      use eos_messages, only: eos_log, eos_report, message_length, unwritten_line
      use eos_mixture, only: compute_mean_weights
      use eos_configuration, only: eos_settings, eos_table_paths_t, opal_vintage_none
      use eos_table_source, only: eos_table_source_t, stage_over_the_gas, adopt_table_result, &
           blend_table_result
      use math_lib
      implicit none
      private
      public :: compute_opal_eos, find_opal_edge_ramp, interpolate_opal_eos, opal_table_source
      public :: iv_opal95, iv_opal01, iv_opal06, opal_layout

contains

!----------------------------------------------------------------------
! opal_table_source
!----------------------------------------------------------------------
! OPAL as a table source (eos_table_source): the vintage the settings
! select (opal_vintage), its table opened at eos_init and read on first
! use, consulted over the gas's result inside the vintage's (T, rho)
! window (inside_opal_window), its density, beta, mean weights and first
! derivatives adopted where its table holds the point and blended in
! with its edge weight near the table's edge (find_opal_edge_ramp).
function opal_table_source() result(src)
      type(eos_table_source_t) :: src
      src%name = 'OPAL'
      src%stage = stage_over_the_gas
      src%selected => opal_selected
      src%read_tables => open_opal_source
      src%in_window => opal_in_window
      src%lookup => opal_lookup
      src%edge_weight => opal_edge_weight
      src%adopt => adopt_table_result
      src%blend => blend_table_result
end function opal_table_source

logical function opal_selected()
      opal_selected = eos_settings%opal_vintage /= opal_vintage_none
end function opal_selected

subroutine open_opal_source(paths, ierr)
      type(eos_table_paths_t), intent(in) :: paths
      integer, intent(out) :: ierr
      call open_opal_eos_table(eos_settings%opal_vintage, trim(paths%opal), ierr)
end subroutine open_opal_source

logical function opal_in_window(log10_temperature, r)
      real(dp), intent(in) :: log10_temperature
      type(eos_result), intent(in) :: r
      opal_in_window = inside_opal_window(eos_settings%opal_vintage, r%temperature, &
           log10_temperature, r%log10_density)
end function opal_in_window

subroutine opal_lookup(log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction, &
     table_r, covers, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: table_r
      logical, intent(out) :: covers
      integer, intent(out) :: ierr
      call compute_opal_eos(eos_settings%opal_vintage, log10_temperature, log10_pressure, &
           hydrogen_fraction, metal_fraction, table_r, ierr, covers)
end subroutine opal_lookup

! the edge weight at the lookup's T and density (table_r's)
subroutine opal_edge_weight(log10_temperature, table_r, weight, in_table, near_edge, ierr)
      real(dp), intent(in) :: log10_temperature
      type(eos_result), intent(in) :: table_r
      real(dp), intent(out) :: weight
      logical, intent(out) :: in_table, near_edge
      integer, intent(out) :: ierr
      call find_opal_edge_ramp(eos_settings%opal_vintage, table_r%temperature, &
           table_r%log10_density, weight, in_table, near_edge, ierr)
end subroutine opal_edge_weight

! the (T, rho) window inside which a vintage is consulted at all (the
! table's edge is then find_opal_edge_ramp's): its row of opal_layout,
! where 1995 bounds log10 T and the others T
pure logical function inside_opal_window(vintage, temperature, log10_temperature, log10_density)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: temperature, log10_temperature, log10_density
      inside_opal_window = temperature >= opal_layout(vintage)%window_t_min .and. &
           temperature <= opal_layout(vintage)%window_t_max .and. &
           log10_temperature <= opal_layout(vintage)%window_log10_t_max .and. &
           log10_density <= opal_layout(vintage)%window_log10_rho_max
end function inside_opal_window

!----------------------------------------------------------------------
! compute_opal_eos (formerly oeqos, oeqos01, oeqos06)
!----------------------------------------------------------------------
! The OPAL EOS of one vintage (iv_opal95, iv_opal01 or iv_opal06) at
! (log T, log P, X, Z): the density that gives P (rhoofp*), the
! thermodynamic quantities there (esac*), beta and the mean weights.
! in_table is false outside the vintage's T range, when no density
! matches P, or when esac* finds the point outside its table. The
! vintage's row of opal_layout gives its T6 range, how many variables
! esac fills (10 for 1995, 9 for the others), their slots, and whether
! to check that the interpolated pressure reproduces P (1995).
subroutine compute_opal_eos(vintage, log10_temperature, log10_pressure, &
     hydrogen_fraction, metal_fraction, r, ierr, in_table)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
      character(len=message_length) :: eos_lines(4)
      real(dp) :: output(n_eos_mv)
      real(dp) :: t_million_k, p_e12, hydrogen_fraction_work, density_cgs, radiation_fraction
      real(dp) :: specific_gas_constant_check
      real(dp) :: chi_rho, chi_t, cv, gamma1, gamma2_ratio
      integer :: rad_flag, deriv_order
! how far the 1995 table's interpolated pressure may be from P
! (relative) at the density the search found
      real(dp), parameter :: pressure_check_tolerance = 0.5d-6

      ierr = 0
      in_table = .true.
      rad_flag = 1     ! adds the radiation correction
      deriv_order = opal_layout(vintage)%n_vars

      r%temperature = exp10(log10_temperature)
      r%pressure = exp10(log10_pressure)
      t_million_k = r%temperature/1.0d6
      p_e12 = r%pressure/1.0d12
      if (t_million_k < opal_layout(vintage)%t6_min .or. &
           t_million_k > opal_layout(vintage)%t6_max) then
         in_table = .false.
         return
      end if
      hydrogen_fraction_work = hydrogen_fraction

! the density that gives P
      density_cgs = opal_density_search(opal_tables(vintage), vintage, hydrogen_fraction_work, &
           t_million_k, p_e12, rad_flag, ierr)
      if (ierr /= 0) return
      if (density_cgs == opal_rho_not_found) then
         in_table = .false.
         return
      end if
      r%density = density_cgs
      r%log10_density = log10(r%density)

! the thermodynamic quantities there, into the vintage's eos_output
      call interpolate_opal_eos(opal_tables(vintage), vintage, hydrogen_fraction_work, &
           t_million_k, density_cgs, deriv_order, rad_flag, ierr, in_table)
      if (.not. in_table) then
         eos_lines = unwritten_line
         write(eos_lines,'(a,a)') trim(opal_layout(vintage)%eos_name), ': the table cannot ' // &
              'interpolate at the density the search found; the result without OPAL stands'
         call eos_log(eos_lines)
         return
      end if
      if (ierr /= 0) return

! dln P/dln rho, dln P/dln T6, cv, gamma1 and gamma2/(gamma2-1) from the
! vintage's results, in its slots
      output = opal_tables(vintage)%eos_output
      if (opal_layout(vintage)%checks_pressure .and. &
           abs((p_e12-output(opal_layout(vintage)%i_p))/p_e12) > pressure_check_tolerance) then
         eos_lines = unwritten_line
         eos_lines(1) = 'OPAL 1995 EOS: the table''s pressure at the density found differs from P ' // &
              'by more than 5e-7 (P12, table):'
         write(eos_lines(2),*) p_e12, output(opal_layout(vintage)%i_p)
         call eos_report(eos_lines)
         ierr = 1
         return
      end if
      chi_rho = output(opal_layout(vintage)%i_chi_rho)
      chi_t = output(opal_layout(vintage)%i_chi_t)
      cv = output(opal_layout(vintage)%i_cv)
      gamma1 = output(opal_layout(vintage)%i_gamma1)
      gamma2_ratio = output(opal_layout(vintage)%i_gamma2_ratio)
      r%dlnrho_dlnp = 1.0d0/chi_rho
      r%dlnrho_dlnt = -chi_t/chi_rho
      r%specific_heat_cp = 1.0d6*cv*gamma1/chi_rho
      r%adiabatic_gradient = 1.0d0/gamma2_ratio

! the radiation pressure's share of P ((a/3) T^4 in OPAL's units), and
! beta
      radiation_fraction = (opal_radiation_constant_over_3*t_million_k*t_million_k)* &
           (t_million_k*t_million_k/p_e12)
      r%beta = 1.0d0 - radiation_fraction
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 4.0d0*r%beta_inverse*radiation_fraction
      r%specific_gas_constant = r%pressure*r%beta/(r%density*r%temperature)
      call compute_mean_weights(r%temperature, r%pressure, r%density, hydrogen_fraction, &
           metal_fraction, eos_settings%use_diffusion_z, specific_gas_constant_check, &
           r%ion_mean_weight_inverse, r%electron_mean_weight_inverse, r%beta)
      if (r%electron_mean_weight_inverse <= 0.0d0) then
         r%electron_mean_weight_inverse = 0.0d0
         r%ion_mean_weight_inverse = r%specific_gas_constant/gas_constant - &
              r%electron_mean_weight_inverse
      end if
end subroutine compute_opal_eos

!----------------------------------------------------------------------
! find_opal_edge_ramp
!----------------------------------------------------------------------
! Whether the vintage's table holds (T, log rho), and the weight of the
! OPAL result near the table's edge: the 2006 table's ragged edge
! (find_ragged_edge_ramp), or the 1995 and 2001 tables' one edge per T6
! row (find_table_edge_ramp).
subroutine find_opal_edge_ramp(vintage, temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
      integer, intent(out) :: ierr

      ierr = 0
      if (opal_layout(vintage)%ragged_edge) then
         call find_ragged_edge_ramp(opal_tables(vintage), opal_layout(vintage), temperature, &
              log10_density, ramp_factor, in_opal_table, needs_ramp)
      else
         call find_table_edge_ramp(opal_tables(vintage), opal_layout(vintage), temperature, &
              log10_density, ramp_factor, in_opal_table, needs_ramp, ierr)
      end if
end subroutine find_opal_edge_ramp

!----------------------------------------------------------------------
! find_table_edge_ramp (formerly eqbound01)
!----------------------------------------------------------------------
! MHP 8/98 and 7/2003: whether a (T, log10 rho) point falls within the
! OPAL 1995 or 2001 EOS table, and near its edge the ramp factor that
! blends the OPAL
! result with the Yale/SCV one (blend_eos's opal_step calls it, through
! find_opal_edge_ramp, right after compute_opal_eos): the vintage's
! density range (log10 rho from ramp_log10_rho_min to the window's
! top), then find_opal_ramp on the instance's grids.
subroutine find_table_edge_ramp(v, lay, temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)
      type(opal_eos_vintage), intent(inout) :: v
      type(opal_vintage_layout), intent(in) :: lay
      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
      integer, intent(out) :: ierr
      real(dp) :: t6

      ierr = 0
      t6 = temperature*1.0d-6
!     exit if outside table in rho
      if (log10_density < lay%ramp_log10_rho_min .or. log10_density > lay%window_log10_rho_max) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if
      call find_opal_ramp(lay%nt, v%t6_grid, v%t_row_index, v%density_edge_at_t, &
           lay%density_index_edge, v%density_grid, t6, log10_density, ramp_factor, &
           in_opal_table, needs_ramp, ierr)
end subroutine find_table_edge_ramp

!----------------------------------------------------------------------
! opal_density_search (formerly rhoofp, rhoofp01 and rhoofp06)
!----------------------------------------------------------------------
! The density at which the vintage's table gives the pressure P (in
! Mb) at (X, T6), by regula falsi through the interpolator's pressure
! alone (trial calls, deriv_order = 1), from a bracket: a first density
! from the table's highest one scaled by P over the pressure there, then
! the lower density shrunk by 0.2 until its pressure is below the target
! (at most 20 times, down to the table's 1e-14), or the upper one grown
! by 5 (to the table's edge); at most 30 steps to |dP/P| < 0.5e-7. The
! vintage's row of opal_layout says how radiation enters (2006 subtracts
! its pressure from the target and asks for the gas alone; 1995 and 2001
! add it to the table's bounds and ask for gas plus radiation), the highest
! density at each T6 row, a/c, and the name in messages. Returns
! opal_rho_not_found (-999) when the point is not in the table, and,
! with a run-log message, when the search does not converge.
function opal_density_search(v, vintage, hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr) result(found_density)
      type(opal_eos_vintage), intent(inout) :: v
      integer, intent(in) :: vintage
      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, pressure_e12
      integer, intent(in) :: rad_flag
      integer, intent(out) :: ierr
      real(dp) :: found_density
      character(len=message_length) :: eos_lines(4)
      integer, parameter :: mx = n_eos_mx, max_refine = 30, max_shrink = 20
! the table-priming call's point (its result is discarded)
      real(dp), parameter :: priming_x = 0.5d0, priming_t6 = 1.0d0, priming_density = 0.001d0
      real(dp) :: rat, radiation_pressure, bound_radiation, pressure_no_rad
      real(dp) :: pressure_max, pressure_min
      real(dp) :: density_trial1, density_trial2, density_trial3
      real(dp) :: pressure_trial1, pressure_trial2, pressure_trial3
      integer :: lo_idx, hi_idx, mid_idx, x_bisect_idx, t6_bisect_idx, nt, edge
      integer :: refine_count, shrink_count, trial_rad
      logical :: in_table

      ierr = 0
      nt = opal_layout(vintage)%nt
      rat = opal_layout(vintage)%search_rad_const_over_c
      radiation_pressure = 0.0d0
      bound_radiation = 0.0d0
      if (opal_layout(vintage)%search_subtracts_radiation) then
         if (rad_flag == 1) radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
         trial_rad = 0
      else
         bound_radiation = rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
         trial_rad = rad_flag
      end if
      pressure_no_rad = pressure_e12 - radiation_pressure

! load the table on first use (the priming call's output is discarded)
      if (opal_shared%table_loaded_flag(vintage) /= opal_flag_set) then
         call interpolate_opal_eos(v, vintage, priming_x, priming_t6, priming_density, 1, &
              trial_rad, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            found_density = opal_rho_not_found
            return
         end if
      end if

! the X interval and the T6 row of the point
      lo_idx = 2
      hi_idx = mx
      do while (hi_idx-lo_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction <= opal_layout(vintage)%x_grid(mid_idx)+1.0d-7) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
      end do
      x_bisect_idx = lo_idx
      lo_idx = nt
      hi_idx = 2
      do while (lo_idx-hi_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (t6_temperature == v%t6_list(1,mid_idx)) then
            lo_idx = mid_idx
            exit
         end if
         if (t6_temperature <= v%t6_list(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
      end do
      t6_bisect_idx = lo_idx
      edge = opal_layout(vintage)%search_density_index_edge(t6_bisect_idx)

! the table's pressure range at that row; outside it the point is not
! in the table
      pressure_max = v%eos_table(x_bisect_idx,1,t6_bisect_idx,edge)*t6_temperature* &
           v%density_grid(edge) + bound_radiation
      pressure_min = v%eos_table(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           v%density_grid(1) + bound_radiation
      if ((pressure_no_rad > 1.25d0*pressure_max) .or. &
           (pressure_no_rad < pressure_min)) then
         found_density = opal_rho_not_found
         return
      end if

! the bracket
      density_trial1 = v%density_grid(edge)*pressure_no_rad/pressure_max
      if (.not. trial_pressure(density_trial1, pressure_trial1)) return
      if (pressure_trial1 > pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         shrink_count = 0
         shrink: do
            shrink_count = shrink_count + 1
            density_trial1 = 0.2d0*density_trial1
            if (density_trial1 < 1.0d-14) density_trial1 = 1.0d-14
            if (.not. trial_pressure(density_trial1, pressure_trial1)) return
            if (pressure_trial1 <= pressure_no_rad) exit shrink
            if (density_trial1 <= 1.0d-14 .or. shrink_count >= max_shrink) exit shrink
         end do shrink
      else
         density_trial2 = 5.0d0*density_trial1
         if (density_trial2 > v%density_grid(edge)) density_trial2 = v%density_grid(edge)
         if (.not. trial_pressure(density_trial2, pressure_trial2)) return
      end if

! regula falsi
      refine_count = 0
      refine: do
         refine_count = refine_count + 1
         density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
              (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
         if (.not. trial_pressure(density_trial3, pressure_trial3)) return
         if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad) < 0.5d-7) then
            found_density = density_trial3
            return
         end if
         if (pressure_trial3 > pressure_no_rad) then
            density_trial2 = density_trial3
            pressure_trial2 = pressure_trial3
         else
            density_trial1 = density_trial3
            pressure_trial1 = pressure_trial3
         end if
         if (refine_count < max_refine) cycle refine
         eos_lines = unwritten_line
         write(eos_lines,'(A,": no convergence after ",I0," tries;", " X, T6, P12, |dP/P| =",4ES12.4," -- falling back")') &
              trim(opal_layout(vintage)%search_name), max_refine, &
              hydrogen_fraction, t6_temperature, pressure_e12, abs((pressure_trial3-pressure_no_rad)/pressure_no_rad)
         call eos_log(eos_lines)
         found_density = opal_rho_not_found
         return
      end do refine

contains

! The table's pressure at density (a trial call); false (with the
! search's result opal_rho_not_found) when the interpolator fails or the
! point is outside the table.
logical function trial_pressure(density, pressure)
      real(dp), intent(in) :: density
      real(dp), intent(out) :: pressure
      call interpolate_opal_eos(v, vintage, hydrogen_fraction, t6_temperature, density, 1, &
           trial_rad, ierr, in_table)
      trial_pressure = ierr == 0 .and. in_table
      if (trial_pressure) then
         pressure = v%eos_output(opal_layout(vintage)%i_p)
      else
         found_density = opal_rho_not_found
      end if
end function trial_pressure

end function opal_density_search

!----------------------------------------------------------------------
! interpolate_opal_eos (formerly esac, esac01 and esac06)
!----------------------------------------------------------------------
! The OPAL interpolation (OPAL's esac, esac01 and esac06, entered so by
! eos_test_entries): the EOS variables at (X, T6, rho) by quadratic
! interpolation in X at each node of a 3 x 3 (or 4 x 4) stencil in T6 and
! rho, blending the two quadratics that overlap where X lies inside the
! grid, then in T6 and rho (t6rinteos), blending the overlapping
! quadratics of a 4 x 4 stencil. Where the stencil would take an empty
! entry (the reader's fill, opal_empty_entry) its coverage sums say so,
! and it shrinks or shifts; a region with no usable stencil is outside
! the table. An X-interpolated value beyond 1e20 either way is an error
! (ierr = 1, named): an entry no real table has. What the vintages
! differ in (the table's size, the slots,
! the radiation correction and constants) is their row of opal_layout, and their messages are written by the report_*
! procedures below. The 1995 table carries 10 variables, the others 9.
!
! With the radiation correction (rad_flag = 1) a full call (deriv_order
! = the table's n_vars) makes every slot that of gas plus radiation
! (add_radiation, 1995: add_radiation_to_gammas); a density search's
! trial call (deriv_order = 1) only adds the radiation pressure to its
! pressure slot, since its other slots are an earlier call's. Without it
! (rad_flag = 0, the 2006 density search's trial calls) the slots are the
! gas's, the cv slot in the table's units.
!
! The slots (the 2001 table lacks slot 4, so its later ones are one
! lower; only 1995's has slot 10):
!   v%eos_output(1) pressure, in megabars (1e12 dyne/cm**2)
!   v%eos_output(2) energy, in 1e12 erg/gm (zero at T6=0)
!   v%eos_output(3) entropy, in units of energy/T6
!   v%eos_output(4) dE/dRho at constant T6
!   v%eos_output(5) specific heat, dE/dT6 at constant volume
!   v%eos_output(6) dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   v%eos_output(7) dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   v%eos_output(8) gamma1 (Cox & Guili eq 9.88)
!   v%eos_output(9) gamma2/(gamma2-1) (Cox & Guili eq 9.88)
!   v%eos_output(10) gamma3-1 (Cox & Guili eq 9.88)
subroutine interpolate_opal_eos(v, vintage, hydrogen_fraction, t6_temperature, &
     density, deriv_order, rad_flag, ierr, in_table)
      integer, intent(in) :: vintage
      type(opal_eos_vintage), intent(inout) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag
      integer, intent(out) :: ierr
! in_table is false when (T6, rho) is outside the table or in an
! empty region of it
      logical, intent(out) :: in_table

      integer, parameter :: mx = n_eos_mx
! the vintage's table: its size and its X grid
      integer :: nr, nt
      real(dp) :: x_grid(n_eos_mx)

      real(dp) :: species_number_fraction(7)

! --- locals ---
      integer :: x_loop_index
      real(dp) :: t6_value, density_value   ! working copies (slt/slr)
      integer :: species_idx
      integer :: lo_idx, hi_idx, mid_idx, result_idx
      integer :: x_index_2, x_index_3, x_index_4, x_index_hi
      integer :: t6_order_hi, density_order_hi   ! ipu/iqu
      real(dp) :: table_sum_2x2, table_sum_3x3, table_sum_3x4, &
           table_sum_4x4
      integer :: eos_var_idx
      integer :: density_scan_idx, t6_scan_idx
      integer :: recompute_flag, cache_slot
      real(dp) :: x_interp_weight
      real(dp) :: pressure_scale
      real(dp) :: total_moles, ground_state_energy, &
           metal_mole_fraction, mean_molecular_weight
      type(opal_vintage_layout) :: lay

      lay = opal_layout(vintage)
      nr = lay%nr
      nt = lay%nt
      x_grid = lay%x_grid
      ierr = 0
      in_table = .true.

      if (deriv_order > lay%n_vars) call report_order_too_high(vintage, deriv_order)
      if ((rad_flag /= 0) .and. (rad_flag /= 1)) then
         call report_bad_rad_flag(vintage, rad_flag)
         ierr = 1
         return
      end if
      t6_value = t6_temperature
      density_value = density
!
      if (opal_shared%table_loaded_flag(vintage) /= opal_flag_set) then
         do species_idx = 1, mx
            v%x_grid_copy(species_idx) = x_grid(species_idx)
         end do
!
! ..... read the data files (the table counts as loaded once it is read)
         call readcoeos(v, vintage, ierr)
         if (ierr /= 0) return
         opal_shared%table_loaded_flag(vintage) = opal_flag_set
         opal_shared%table_metal_fraction(vintage) = v%z_table(1)

         if (opal_shared%table_metal_fraction(vintage)+hydrogen_fraction-1.0d-6 > 1) then
            call report_mass_fractions(vintage, opal_shared%table_metal_fraction(vintage), &
                 hydrogen_fraction)
            ierr = 1
            return
         end if
      end if
!
!
! ..... Determine T6,rho grid points to use in the
!       interpolation.
      if ((t6_value > v%t6_grid(1)) .or. (t6_value < v%t6_grid(nt)) .or. &
           (density_value < v%density_grid(1)) .or. (density_value > v%density_grid(nr))) then
         call report_outside_table(vintage, t6_value, density_value, v%t6_grid(1), &
              v%t6_grid(nt), v%density_grid(1), v%density_grid(nr))
         in_table = .false.
         return
      end if
!
!
!
      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction <= x_grid(mid_idx)+1.0d-7) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
    end do
      result_idx = hi_idx
      v%x_index_lo = result_idx - 2
      x_index_2 = result_idx - 1
      x_index_3 = result_idx
      x_index_4 = result_idx + 1
      x_index_hi = x_index_4
      if (hydrogen_fraction < 1.0d-6) then
         v%x_index_lo = 1
         x_index_2 = 1
         x_index_3 = 1
         x_index_4 = 2
         x_index_hi = 1
      end if
      if ((hydrogen_fraction <= x_grid(2)+1.0d-7) .or. &
           (hydrogen_fraction >= x_grid(mx-2)-1.0d-7)) x_index_hi = x_index_3
!
      lo_idx = 2
      hi_idx = nr
   do while (hi_idx-lo_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (density_value == v%density_grid(mid_idx)) then
            hi_idx = mid_idx
            exit
         end if
         if (density_value <= v%density_grid(mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
   end do
      result_idx = hi_idx
      v%density_index_1 = result_idx - 2
      v%density_index_2 = result_idx - 1
      v%density_index_3 = result_idx
      v%density_index_4 = v%density_index_3 + 1
      density_order_hi = 3
      if (v%density_index_4 > nr) density_order_hi = 2
!
      lo_idx = nt
      hi_idx = 2
   do while (lo_idx-hi_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (t6_temperature == v%t6_list(1,mid_idx)) then
            lo_idx = mid_idx
            exit
         end if
         if (t6_temperature <= v%t6_list(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
   end do
      result_idx = lo_idx
      v%t6_index_1 = result_idx - 2
      v%t6_index_2 = result_idx - 1
      v%t6_index_3 = result_idx
      v%t6_index_4 = v%t6_index_3 + 1
      t6_order_hi = 3
      if (v%t6_index_4 > nt) t6_order_hi = 2
      if (v%t6_index_3 == 0) call report_t6_search(vintage, hi_idx, lo_idx, mid_idx)

!     check to determine if interpolation indices fall within
!     table boundaries.  choose largest allowed size.
      call clamped_coverage(v, x_index_hi, t6_order_hi, density_order_hi, table_sum_2x2, &
           table_sum_3x3, table_sum_3x4, table_sum_4x4)
      v%density_interp_order = 2
      v%t6_interp_order = 2
      if (table_sum_3x3 > opal_empty_sum) then
         if (table_sum_2x2 < opal_usable_corner_sum) then
            v%t6_index_1 = v%t6_index_3 - 3
            v%t6_index_2 = v%t6_index_1 + 1
            v%t6_index_3 = v%t6_index_2 + 1
            v%density_index_1 = v%density_index_3 - 3
            v%density_index_2 = v%density_index_1 + 1
            v%density_index_3 = v%density_index_2 + 1
         else
            call report_empty_region(vintage, hydrogen_fraction, t6_temperature, density)
            in_table = .false.
            return
         end if
      else
      if (table_sum_3x4 < opal_empty_sum) v%t6_interp_order = 3
      if (table_sum_4x4 < opal_empty_sum) v%density_interp_order = 3

      if (t6_temperature >= v%t6_list(1,2)+1.0d-7) v%t6_interp_order = 2
      if (density_value <= v%density_grid(2)+1.0d-15) v%density_interp_order = 2

      if ((v%density_index_3 == nr) .or. (v%t6_index_3 == nt)) then
         v%density_interp_order = 2
         v%t6_interp_order = 2
      end if

      end if
      do eos_var_idx = 1, deriv_order
! the stencil's entries of this variable, for each X table
      do x_loop_index = v%x_index_lo, x_index_hi
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            v%x_interp_workspace(x_loop_index,t6_scan_idx,density_scan_idx) = &
                 v%eos_table(x_loop_index,eos_var_idx,t6_scan_idx,density_scan_idx)
         end do
      end do
      end do
      if ((v%z_table(x_index_2) /= v%z_table(v%x_index_lo)) .or. &
           (v%z_table(x_index_3) /= v%z_table(v%x_index_lo))) then
         call report_tables_z(vintage)
         ierr = 1
         return
      end if
      if (opal_shared%table_metal_fraction(vintage) /= v%z_table(v%x_index_lo)) exit
! The X-interpolated value at each node of the stencil. A value beyond
! opal_x_interp_overflow either way came from an empty or corrupt entry
! (a huge entry with a negative Lagrange weight comes back hugely
! negative, M292): it is named and the call returns ierr = 1. The check
! is on the magnitude, not the sign: dE/dRho is negative in most entries.
      recompute_flag = 0
      cache_slot = 1
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            if (x_index_hi == 1) then
               v%x_interp_result(t6_scan_idx,density_scan_idx) = &
                    v%x_interp_workspace(v%x_index_lo,t6_scan_idx,density_scan_idx)
               if (abs(v%x_interp_result(t6_scan_idx,density_scan_idx)) > opal_x_interp_overflow) then
                  call report_x_interpolation(vintage, v, t6_scan_idx, density_scan_idx, &
                       v%x_index_lo, v%x_index_lo)
                  ierr = 1
                  return
               end if
               cycle
            end if
            v%x_interp_result(t6_scan_idx,density_scan_idx) = &
                 quad(v%quad, recompute_flag, cache_slot, hydrogen_fraction, &
                 v%x_interp_workspace(v%x_index_lo,t6_scan_idx,density_scan_idx), &
                 v%x_interp_workspace(x_index_2,t6_scan_idx,density_scan_idx), &
                 v%x_interp_workspace(x_index_3,t6_scan_idx,density_scan_idx), &
                 v%x_grid_copy(v%x_index_lo), v%x_grid_copy(x_index_2), &
                 v%x_grid_copy(x_index_3))
            if (abs(v%x_interp_result(t6_scan_idx,density_scan_idx)) > opal_x_interp_overflow) then
               call report_x_interpolation(vintage, v, t6_scan_idx, density_scan_idx, &
                    v%x_index_lo, x_index_3)
               ierr = 1
               return
            end if
            recompute_flag = 1
         end do
      end do

      if (x_index_4 == x_index_hi) then ! interpolate between quadratics
      recompute_flag = 0
      cache_slot = 1
      x_interp_weight = (v%x_grid_copy(x_index_3) - hydrogen_fraction)* &
           v%x_grid_spacing_inv(x_index_3)
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            v%x_interp_result_alt(t6_scan_idx,density_scan_idx) = &
                 quad(v%quad, recompute_flag, cache_slot, hydrogen_fraction, &
                 v%x_interp_workspace(x_index_2,t6_scan_idx,density_scan_idx), &
                 v%x_interp_workspace(x_index_3,t6_scan_idx,density_scan_idx), &
                 v%x_interp_workspace(x_index_4,t6_scan_idx,density_scan_idx), &
                 v%x_grid_copy(x_index_2), v%x_grid_copy(x_index_3), &
                 v%x_grid_copy(x_index_4))
            if (abs(v%x_interp_result_alt(t6_scan_idx,density_scan_idx)) > opal_x_interp_overflow) then
               call report_x_interpolation(vintage, v, t6_scan_idx, density_scan_idx, &
                    x_index_2, x_index_4)
               ierr = 1
               return
            end if
            v%x_interp_result(t6_scan_idx,density_scan_idx) = &
                 v%x_interp_result(t6_scan_idx,density_scan_idx)*x_interp_weight + &
                 v%x_interp_result_alt(t6_scan_idx,density_scan_idx)* &
                 (1.0d0 - x_interp_weight)
            recompute_flag = 1
         end do
      end do


      end if

      recompute_flag = 0
!
! ..... completed X interpolation. Now interpolate T6 and rho on a
!       4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure
!       mixes overlapping quadratics to obtain smoothed derivatives.
      call t6rinteos(v, vintage, density_value, t6_value, ierr)
      if (ierr /= 0) return
      v%eos_output(eos_var_idx) = v%esact
      end do
      if (eos_var_idx > deriv_order) then
      pressure_scale = t6_temperature*density
      v%eos_output(lay%i_p) = v%eos_output(lay%i_p)*pressure_scale ! interpolated in p/po
! only the slots this call interpolated (index <= deriv_order) are
! rescaled
      if (lay%i_e <= deriv_order) &
           v%eos_output(lay%i_e) = v%eos_output(lay%i_e)*t6_temperature ! interpolated in E/T6
      if (rad_flag == 1 .and. deriv_order >= lay%n_vars) then
         mean_molecular_weight = opal_mass_per_mole(lay, hydrogen_fraction, &
              opal_shared%table_metal_fraction(vintage), total_moles, ground_state_energy, &
              metal_mole_fraction, species_number_fraction)
         if (lay%radiation_gammas_direct) then
            call add_radiation_to_gammas(v, lay, t6_temperature, density, total_moles, &
                 mean_molecular_weight)
         else
            call add_radiation(v, lay, t6_temperature, density, total_moles, &
                 mean_molecular_weight)
         end if
      else if (rad_flag == 1) then
! a density search's trial call (only the pressure slot is fresh, the
! others hold an earlier call's values): the radiation pressure is added
! to the pressure alone, as the full correction would
         if (lay%radiation_gammas_direct) then
            v%eos_output(lay%i_p) = (v%eos_output(lay%i_p) + &
                 opal_radiation_pressure(lay, t6_temperature))*lay%energy_unit_ratio
         else
            v%eos_output(lay%i_p) = v%eos_output(lay%i_p) + &
                 opal_radiation_pressure(lay, t6_temperature)
         end if
      end if
      return
      end if
      call report_run_z(vintage, v, hydrogen_fraction, t6_temperature, density)
      ierr = 1

end subroutine interpolate_opal_eos

! The coverage sums of a stencil, which stay inside the table (the 1995
! ones too since 2026-11, M94): the pressure entries of the X tables x_index_lo..x_index_hi over
! the 2 x 2 and 3 x 3 cells from the stencil's corner, the 3 x 4 (T6 by
! density) one where the stencil has a fourth T6 row (else 2e30, past
! the empty-entry threshold) and the 4 x 4 one where it has a fourth
! density column (else 2e30).
pure subroutine clamped_coverage(v, x_index_hi, t6_order_hi, density_order_hi, sum_2x2, &
     sum_3x3, sum_3x4, sum_4x4)
      type(opal_eos_vintage), intent(in) :: v
      integer, intent(in) :: x_index_hi, t6_order_hi, density_order_hi
      real(dp), intent(out) :: sum_2x2, sum_3x3, sum_3x4, sum_4x4
      integer :: x_loop_index, density_scan_idx, t6_scan_idx
      sum_2x2 = 0.0d0
      sum_3x3 = 0.0d0
      sum_3x4 = 0.0d0
      sum_4x4 = 0.0d0
      do x_loop_index = v%x_index_lo, x_index_hi
         do density_scan_idx = v%density_index_1, v%density_index_1+1
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+1
               sum_2x2 = sum_2x2 + v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         do density_scan_idx = v%density_index_1, v%density_index_1+2
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+2
               sum_3x3 = sum_3x3 + v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         if (t6_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+2
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  sum_3x4 = sum_3x4 + v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            sum_3x4 = opal_beyond_stencil
         end if
         if (density_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+3
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  sum_4x4 = sum_4x4 + v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            sum_4x4 = opal_beyond_stencil
         end if
      end do
end subroutine clamped_coverage

!----------------------------------------------------------------------
! t6rinteos (formerly t6rinteos01 and t6rinteos06)
!----------------------------------------------------------------------
! interpolate_opal_t6_rho, with each vintage's out-of-range message.
subroutine t6rinteos(v, vintage, slr, slt, ierr)
      type(opal_eos_vintage), intent(inout) :: v
      integer, intent(in) :: vintage
      real(dp), intent(in) :: slr, slt
      integer, intent(out) :: ierr

      ierr = 0
      call interpolate_opal_t6_rho(v, slr, slt)
      if (v%esact > opal_range_overflow) then
         call report_interpolation_indices(vintage)
         ierr = 1
         return
      end if
end subroutine t6rinteos


!----------------------------------------------------------------------
! find_ragged_edge_ramp (formerly eqbound06)
!----------------------------------------------------------------------
! LLP 10/16/06: checks whether a (T, log10 rho) point falls within the
! OPAL 2006 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (blend_eos's opal_step calls it, through find_opal_edge_ramp,
! right after compute_opal_eos). Unlike find_table_edge_ramp (1995 and
! 2001), it tracks the table's ragged edge through the vintage's
! per-column and per-row index arrays (lay%t6_index_lo,
! lay%density_index_edge) rather than one edge value per T6 row, from
! the stencil the last interpolation left in v.
pure subroutine find_ragged_edge_ramp(v, lay, temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp)
      type(opal_eos_vintage), intent(in) :: v
      type(opal_vintage_layout), intent(in) :: lay

      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp



      integer :: t6_row, density_row
! --- locals ---
      real(dp) :: t6, density, density_ramp_factor, &
           t6_ramp_factor

      t6 = temperature*1.0d-6
      density = pow(10d0, log10_density)

!     Exit if outside table in rho
      if ((density < v%density_grid(1)) .or. (density >= v%density_grid(lay%nr))) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if

!     Exit if outside table in T6
      if ((t6 > v%t6_grid(1)) .or. (t6 <= v%t6_grid(lay%nt))) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if

!     Initialize
      density_ramp_factor = 1d0   ! Initialize density ramp factor to 1
      t6_ramp_factor = 1d0        ! Initialize temperature ramp factor to 1
      in_opal_table = .true.      ! Presume that T6,D point is in OPAL table
      needs_ramp = .false.        ! presume ramping is not needed,
                                   ! i.e. we are not in the border of the table.

!     Check for easy borders
      if (t6 > v%t6_grid(2)) then
          t6_ramp_factor = (t6-v%t6_grid(1))/(v%t6_grid(2)-v%t6_grid(1)) ! Ramp in temperature
          needs_ramp = .true.
      end if
      if (density < v%density_grid(2)) then
          density_ramp_factor = (density - v%density_grid(1))/ &
               (v%density_grid(2)-v%density_grid(1)) ! Ramp in density
          needs_ramp = .true.
      end if

!     Find a t6_row such that v%t6_grid(t6_row-1) >= T6 > v%t6_grid(t6_row)
!     and a density_row such that v%density_grid(density_row-1) < D <=
!     v%density_grid(density_row)
!     The ESAC06 variables (k=)v%t6_index_3 and (l-)v%density_index_2 are
!     close to what we need. A linear search will work fine.
      t6_row = v%t6_index_3
   do while (t6 <= v%t6_grid(t6_row))
        t6_row = t6_row+1
   end do
   do while (t6 > v%t6_grid(t6_row-1))
        t6_row = t6_row-1
   end do
!     We now have: v%t6_grid(t6_row-1) >= T6 > v%t6_grid(t6_row)

      density_row = v%density_index_2
   do while (density > v%density_grid(density_row))
        density_row = density_row+1
   end do
   do while (density <= v%density_grid(density_row-1))
        density_row = density_row-1
   end do
!     We now have: v%density_grid(density_row-1) < D <= v%density_grid(density_row)

!     For a given temperature in array v%t6_grid with index t6_row, e.g.,
!     v%t6_grid(t6_row), element lay%density_index_edge(t6_row) of array
!     lay%density_index_edge contains the index to the max allowed density
!     in array v%density_grid. So, for a given v%t6_grid(t6_row), the
!     associated D must be less than v%density_grid(lay%density_index_edge(t6_row)).
!     A valid index in array v%density_grid must be less than or equal
!     to lay%density_index_edge(t6_row)

      if (density_row > lay%density_index_edge(t6_row)) then
        in_opal_table = .false.
        needs_ramp = .true.
        ramp_factor = 0d0
        return
      end if
      if (density_row == lay%density_index_edge(t6_row)) then
        needs_ramp = .true.
        density_ramp_factor = (v%density_grid(density_row)-density)/ &
             (v%density_grid(density_row)-v%density_grid(density_row-1))
      end if

!     For a given density in array v%density_grid with index density_row,
!     i.e., v%density_grid(density_row), element lay%t6_index_lo(density_row)
!     of array lay%t6_index_lo contains the index to the min allowed
!     temperature in array v%t6_grid. So for a given
!     v%density_grid(density_row), the associated T6 must be greater
!     than v%t6_grid(lay%t6_index_lo(density_row)). (Note that v%t6_grid is a
!     decreasing array.) A valid index in v%t6_grid must be less than or
!     equal to lay%t6_index_lo(density_row).

      if (t6_row > lay%t6_index_lo(density_row)) then
        in_opal_table = .false.
        needs_ramp = .true.
        ramp_factor = 0d0
        return
      end if
      if (t6_row == lay%t6_index_lo(density_row)) then
        needs_ramp = .true.
        t6_ramp_factor = (t6-v%t6_grid(t6_row))/(v%t6_grid(t6_row-1)-v%t6_grid(t6_row))
      end if

      ramp_factor = density_ramp_factor * t6_ramp_factor
end subroutine find_ragged_edge_ramp



!----------------------------------------------------------------------
! add_radiation (formerly radsub01 and radsub06)
!----------------------------------------------------------------------
! The radiation correction of the 2001 and 2006 tables, on the values
! already in v%eos_output (pressure in Mb, energies in Mb cm^3/g, T6
! units): P_rad = (4/3)(a/c) T6^4, E_rad = 3 P_rad/rho and S_rad =
! (4/3) E_rad/T6 are added to P, E and S (and -E_rad/rho to dE/dRho,
! where the table has it), and chi_rho, chi_T and cv become those of
! gas plus radiation. The gammas are not recomputed from these (the cv
! is not accurate enough, Landau & Lifshitz eq. 16.16): each gets the
! difference between its gas-plus-radiation and its gas-only value from
! the same formulas (gamma3 - 1 only where the table has it, 2001). The
! vintage's row of opal_layout gives the slots and the constants.
pure subroutine add_radiation(v, lay, t6_temperature, density, total_moles, &
     mean_molecular_weight)
      type(opal_eos_vintage), intent(inout) :: v
      type(opal_vintage_layout), intent(in) :: lay
      real(dp), intent(in) :: t6_temperature, density
      real(dp), intent(in) :: total_moles
      real(dp), intent(in) :: mean_molecular_weight
      real(dp) :: moles_per_ev
      real(dp) :: radiation_pressure, radiation_energy, radiation_entropy
      real(dp) :: total_pressure, chi_rho, chi_t6, molar_specific_heat
      real(dp) :: gamma3_minus1_norad, gamma1_norad, &
           gamma2_over_gamma2_minus1_norad
      real(dp) :: gamma3_minus1, gamma1, gamma2_over_gamma2_minus1

      moles_per_ev = total_moles*lay%molar_gas_constant_mbcc  ! Mb-cc/unit T6
      radiation_pressure = opal_radiation_pressure(lay, t6_temperature)   ! Mb
      radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
      radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)

! the gammas of the gas alone
      total_pressure = v%eos_output(lay%i_p)
      chi_rho = v%eos_output(lay%i_chi_rho)*v%eos_output(lay%i_p)/ &
           total_pressure
      chi_t6 = (v%eos_output(lay%i_p)*v%eos_output(lay%i_chi_t))/ &
           total_pressure
      molar_specific_heat = (v%eos_output(lay%i_cv)*moles_per_ev/ &
           mean_molecular_weight)
      gamma3_minus1_norad = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1_norad = chi_rho + chi_t6*gamma3_minus1_norad
      gamma2_over_gamma2_minus1_norad = gamma1_norad/gamma3_minus1_norad

! gas plus radiation
      total_pressure = v%eos_output(lay%i_p) + radiation_pressure
      chi_rho = v%eos_output(lay%i_chi_rho)*v%eos_output(lay%i_p)/ &
           total_pressure
      chi_t6 = (v%eos_output(lay%i_p)*v%eos_output(lay%i_chi_t) &
           + 4.0d0*radiation_pressure)/total_pressure
      molar_specific_heat = (v%eos_output(lay%i_cv)*moles_per_ev/ &
           mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
      gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1 = chi_rho + chi_t6*gamma3_minus1
      gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1

      v%eos_output(lay%i_gamma1) = v%eos_output(lay%i_gamma1) + &
           gamma1 - gamma1_norad
      v%eos_output(lay%i_gamma2_ratio) = v%eos_output(lay%i_gamma2_ratio) + &
           gamma2_over_gamma2_minus1 - gamma2_over_gamma2_minus1_norad
      if (lay%i_gamma3m1 > 0) v%eos_output(lay%i_gamma3m1) = v%eos_output(lay%i_gamma3m1) + &
           gamma3_minus1 - gamma3_minus1_norad
      if (lay%i_dedrho > 0) v%eos_output(lay%i_dedrho) = v%eos_output(lay%i_dedrho) - &
           radiation_energy/density
      v%eos_output(lay%i_p) = total_pressure
      v%eos_output(lay%i_e) = v%eos_output(lay%i_e) + radiation_energy
      v%eos_output(lay%i_s) = v%eos_output(lay%i_s) + radiation_entropy
      v%eos_output(lay%i_cv) = molar_specific_heat
      v%eos_output(lay%i_chi_rho) = chi_rho
      v%eos_output(lay%i_chi_t) = chi_t6
end subroutine add_radiation

! The radiation pressure (4/3)(a/c) T6^4 in Mb, with the vintage's a/c.
pure real(dp) function opal_radiation_pressure(lay, t6_temperature)
      type(opal_vintage_layout), intent(in) :: lay
      real(dp), intent(in) :: t6_temperature
      opal_radiation_pressure = 4.0d0/3.0d0*lay%rad_const_over_c*t6_temperature**4
end function opal_radiation_pressure

!----------------------------------------------------------------------
! add_radiation_to_gammas (formerly radsub)
!----------------------------------------------------------------------
! The radiation correction of the 1995 table (lay%radiation_gammas_direct),
! on the values already in v%eos_output: P_rad = (4/3)(a/c) T6^4, E_rad =
! 3 P_rad/rho and S_rad = (4/3) E_rad/T6 are added to P, E and S, and
! -E_rad/rho to dE/dRho; chi_rho, chi_T and cv become those of gas plus
! radiation, and the three gammas are computed from them (not corrected
! by a difference, as add_radiation's are). P and S are then rescaled by
! lay%energy_unit_ratio, the ratio of the latest N_A k in erg/eV
! (0.9648575e12: 6.0221367e23 x 1.380658e-16 erg/K x 11604.5 K/eV) to the
! one the 1995 tables were computed with (0.9652e12); E and dE/dRho are
! not (as in OPAL's code).
pure subroutine add_radiation_to_gammas(v, lay, t6_temperature, density, total_moles, &
     mean_molecular_weight)
      type(opal_eos_vintage), intent(inout) :: v
      type(opal_vintage_layout), intent(in) :: lay
      real(dp), intent(in) :: t6_temperature, density
      real(dp), intent(in) :: total_moles
      real(dp), intent(in) :: mean_molecular_weight
      real(dp) :: moles_per_ev, radiation_pressure, &
           radiation_energy, radiation_entropy
      real(dp) :: total_pressure, total_energy, total_entropy
      real(dp) :: chi_rho, chi_t6
      real(dp) :: molar_specific_heat, gamma3_minus1, gamma1, &
           gamma2_over_gamma2_minus1
      real(dp) :: de_drho_at_t

      moles_per_ev = total_moles*lay%molar_gas_constant_mbcc  ! Mb-cc/unit T6
      radiation_pressure = opal_radiation_pressure(lay, t6_temperature)   ! Mb
      radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
      radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)
      total_pressure = v%eos_output(lay%i_p) + radiation_pressure
      total_energy = v%eos_output(lay%i_e) + radiation_energy
      total_entropy = v%eos_output(lay%i_s) + radiation_entropy
      chi_rho = v%eos_output(lay%i_chi_rho)*v%eos_output(lay%i_p)/ &
           total_pressure
      chi_t6 = (v%eos_output(lay%i_p)*v%eos_output(lay%i_chi_t) &
           + 4.0d0*radiation_pressure)/total_pressure
      molar_specific_heat = (v%eos_output(lay%i_cv)*moles_per_ev/ &
           mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
      gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1 = chi_rho + chi_t6*gamma3_minus1
      gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1
      de_drho_at_t = v%eos_output(lay%i_dedrho) - radiation_energy/density
      total_pressure = total_pressure*lay%energy_unit_ratio
      total_entropy = total_entropy*lay%energy_unit_ratio
      v%eos_output(lay%i_p) = total_pressure
      v%eos_output(lay%i_e) = total_energy
      v%eos_output(lay%i_s) = total_entropy
      v%eos_output(lay%i_dedrho) = de_drho_at_t
      v%eos_output(lay%i_cv) = molar_specific_heat
      v%eos_output(lay%i_chi_rho) = chi_rho
      v%eos_output(lay%i_chi_t) = chi_t6
      v%eos_output(lay%i_gamma1) = gamma1
      v%eos_output(lay%i_gamma2_ratio) = gamma2_over_gamma2_minus1
      v%eos_output(lay%i_gamma3m1) = gamma3_minus1
end subroutine add_radiation_to_gammas

!----------------------------------------------------------------------
! opal_mass_per_mole (formerly gmass and gmass01)
!----------------------------------------------------------------------
! The mean mass per particle (amu) of fully ionized H, He, C, N, O and Ne
! at (X, Z), the metals in a fixed solar-relative C:N:O:Ne mix, for the
! radiation correction's specific heat; also the number (mole) fractions
! of the species (Ne, O, N, C, He, H), their total number with their
! electrons (total_moles), the ground-state ionization energy and the
! metals' moles per unit mass. The vintage's row gives the ionization
! energies, and whether the electrons' mass counts (1995: the sum then
! starts from (total_moles - 1) m_e, the electrons per ion; 2001 and
! 2006 start from 0).
real(dp) function opal_mass_per_mole(lay, hydrogen_fraction, metal_fraction, &
     total_moles, ground_state_energy, metal_mole_fraction, species_number_fraction)
      type(opal_vintage_layout), intent(in) :: lay
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: total_moles, ground_state_energy, &
           metal_mole_fraction
! species_number_fraction(1:6): number (mole) fractions of Ne, O, N, C,
! He, H (same order as atomic_number); slot 7 is never assigned.
      real(dp), intent(out) :: species_number_fraction(:)
! atomic_number(i): atomic number (= electrons per fully-ionized atom)
! of species i, for i=1,6 = Ne, O, N, C, He, H.
      real(dp), parameter :: atomic_number(6) = [10.0d0, 8.0d0, &
           7.0d0, 6.0d0, 2.0d0, 1.0d0]
! atomic_weight(i): atomic weight (amu), i=1 electron, 2 Ne, 3 O, 4 N,
! 5 C, 6 He, 7 H.
      real(dp), parameter :: atomic_weight(7) = [0.00054858d0, 20.179d0, 15.9994d0, &
           14.0067d0, 12.011d0, 4.0026d0, 1.0079d0]
! solar-relative mole-fraction-per-unit-Z of C, N, O, Ne in the metal
! mix (fixed composition assumed for all Z).
      real(dp), parameter :: carbon_mix_fraction = 0.247137766d0, &
           nitrogen_mix_fraction = 0.0620782d0, &
           oxygen_mix_fraction = 0.52837118d0, &
           neon_mix_fraction = 0.1624188d0
      real(dp) :: carbon_moles, nitrogen_moles, oxygen_moles, &
           neon_moles
      real(dp) :: hydrogen_moles, helium_moles, total_moles_raw
      integer :: species_idx

      metal_mole_fraction = metal_fraction/(carbon_mix_fraction* &
           atomic_weight(5) + nitrogen_mix_fraction*atomic_weight(4) + &
           oxygen_mix_fraction*atomic_weight(3) + neon_mix_fraction* &
           atomic_weight(2))
      carbon_moles = metal_mole_fraction*carbon_mix_fraction
      nitrogen_moles = metal_mole_fraction*nitrogen_mix_fraction
      oxygen_moles = metal_mole_fraction*oxygen_mix_fraction
      neon_moles = metal_mole_fraction*neon_mix_fraction
      hydrogen_moles = hydrogen_fraction/atomic_weight(7)
      helium_moles = (1.0d0 - hydrogen_fraction - metal_fraction)/ &
           atomic_weight(6)
      total_moles_raw = hydrogen_moles + helium_moles + carbon_moles + &
           nitrogen_moles + oxygen_moles + neon_moles
      species_number_fraction(6) = hydrogen_moles/total_moles_raw
      species_number_fraction(5) = helium_moles/total_moles_raw
      species_number_fraction(4) = carbon_moles/total_moles_raw
      species_number_fraction(3) = nitrogen_moles/total_moles_raw
      species_number_fraction(2) = oxygen_moles/total_moles_raw
      species_number_fraction(1) = neon_moles/total_moles_raw
      ground_state_energy = 0.0d0
      total_moles = 0.0d0
      do species_idx = 1, 6
         ground_state_energy = ground_state_energy + &
              lay%ionization_energy(species_idx)*species_number_fraction(species_idx)
         total_moles = total_moles + (1.0d0 + atomic_number(species_idx))* &
              species_number_fraction(species_idx)
      end do
      opal_mass_per_mole = 0.0d0
      if (lay%counts_electron_mass) opal_mass_per_mole = (total_moles - 1.0d0)*atomic_weight(1)
      do species_idx = 2, 7
         opal_mass_per_mole = opal_mass_per_mole + atomic_weight(species_idx)* &
              species_number_fraction(species_idx - 1)
      end do
end function opal_mass_per_mole

!----------------------------------------------------------------------
! The interpolator's messages, one text for every vintage after its name
! (opal_layout's eos_name), to the run log (eos_log) or, for an error, to
! the terminal and the run log (eos_report). The numbers in parentheses
! are OPAL's own error numbers.
!----------------------------------------------------------------------
! deriv_order above the table's variable count (the call goes on)
subroutine report_order_too_high(vintage, deriv_order)
      integer, intent(in) :: vintage, deriv_order
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0,a,i0,a)') trim(opal_layout(vintage)%eos_name), ': deriv_order ', &
           deriv_order, ' is above the table''s ', opal_layout(vintage)%n_vars, ' variables'
      call eos_log(eos_lines)
end subroutine report_order_too_high

! rad_flag neither 0 nor 1
subroutine report_bad_rad_flag(vintage, rad_flag)
      integer, intent(in) :: vintage, rad_flag
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0)') trim(opal_layout(vintage)%eos_name), &
           ': rad_flag must be 0 or 1, not ', rad_flag
      call eos_report(eos_lines)
end subroutine report_bad_rad_flag

! the table's Z and the first call's X add up to more than 1
subroutine report_mass_fractions(vintage, table_z, hydrogen_fraction)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: table_z, hydrogen_fraction
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,2es14.6)') trim(opal_layout(vintage)%eos_name), &
           ': the table''s Z and X add up to more than 1 (61): Z, X =', table_z, hydrogen_fraction
      call eos_report(eos_lines)
end subroutine report_mass_fractions

! (T6, rho) outside the table
subroutine report_outside_table(vintage, t6, density, t6_first, t6_last, &
     density_first, density_last)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: t6, density, t6_first, t6_last, density_first, density_last
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a/a,es14.6,a,es14.6,a,es14.6,a/a,es14.6,a,es14.6,a,es14.6,a)') &
           trim(opal_layout(vintage)%eos_name), ': (T6, rho) outside the table (62):', &
           '     T6 ', t6, ' (the table: ', t6_first, ' to ', t6_last, ')', &
           '     rho', density, ' (the table: ', density_first, ' to ', density_last, ')'
      call eos_log(eos_lines)
end subroutine report_outside_table

! the T6 search found no row
subroutine report_t6_search(vintage, hi_idx, lo_idx, mid_idx)
      integer, intent(in) :: vintage, hi_idx, lo_idx, mid_idx
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,3i5)') trim(opal_layout(vintage)%eos_name), &
           ': the T6 search found no row: ihi, ilo, imd =', hi_idx, lo_idx, mid_idx
      call eos_log(eos_lines)
end subroutine report_t6_search

! (T6, rho) in an empty region of the table
subroutine report_empty_region(vintage, hydrogen_fraction, t6, density)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: hydrogen_fraction, t6, density
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,3es12.4)') trim(opal_layout(vintage)%eos_name), &
           ': (T6, rho) in an empty region of the table (65): X, T6, rho =', &
           hydrogen_fraction, t6, density
      call eos_log(eos_lines)
end subroutine report_empty_region

! the X tables of the stencil are not all for one Z
subroutine report_tables_z(vintage)
      integer, intent(in) :: vintage
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a)') trim(opal_layout(vintage)%eos_name), &
           ': the stencil''s X tables are not all for one Z'
      call eos_report(eos_lines)
end subroutine report_tables_z

! an X-interpolated value beyond 1e20 either way (an empty or corrupt
! entry in the stencil, an error): the stencil and the table values the
! interpolation read, X tables x_first to x_last
subroutine report_x_interpolation(vintage, v, t6_idx, density_idx, x_first, x_last)
      integer, intent(in) :: vintage, t6_idx, density_idx, x_first, x_last
      type(opal_eos_vintage), intent(in) :: v
      character(len=message_length) :: eos_lines(4)
      integer :: k
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0,a,i0,a,4i5)') trim(opal_layout(vintage)%eos_name), &
           ': an X-interpolated value beyond 1e20 (an empty or corrupt entry) at T6 row ', &
           t6_idx, ', density column ', density_idx, '; l3, k3, iq, ip =', v%density_index_3, &
           v%t6_index_3, v%density_interp_order, v%t6_interp_order
      write(eos_lines(2),'(a,3es12.4)') '     the table''s values there:', &
           (v%x_interp_workspace(k,t6_idx,density_idx), k=x_first,x_last)
      call eos_report(eos_lines)
end subroutine report_x_interpolation

! an X table's Z differs from the first X table's (the table's Z)
subroutine report_run_z(vintage, v, hydrogen_fraction, t6, density)
      integer, intent(in) :: vintage
      type(opal_eos_vintage), intent(in) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6, density
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0,a,es12.4,a,es12.4,a/a,4i5,a,3es12.4)') &
           trim(opal_layout(vintage)%eos_name), ': X table ', v%x_index_lo, '''s Z,', &
           v%z_table(v%x_index_lo), ', is not the table''s,', &
           opal_shared%table_metal_fraction(vintage), ' (66)', &
           '     iq, ip, k3, l3 =', v%t6_interp_order, v%density_interp_order, v%t6_index_3, &
           v%density_index_3, '; X, T6, rho =', hydrogen_fraction, t6, density
      call eos_report(eos_lines)
end subroutine report_run_z

! the T6-rho interpolation's value out of range (an empty entry)
subroutine report_interpolation_indices(vintage)
      integer, intent(in) :: vintage
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a)') trim(opal_layout(vintage)%eos_name), &
           ': the T6-rho interpolation''s value is out of range (an empty entry in its ' // &
           'stencil); please report the conditions'
      call eos_report(eos_lines)
end subroutine report_interpolation_indices
!======================================================================
! Sources: the OPAL EOS interpolation codes of Rogers, Swenson & Iglesias
! (1996; the 1995 tables: esac, rhoofp, radsub, gmass, readco) and Rogers
! & Nayfonov (2002; the 2001 tables: esac01, rhoofp01, radsub01;
! the 2006 tables, the 2005 update: esac06, rhoofp06, radsub06), with
! the YREC edge ramps EQBOUND (MHP 8/98, 1995), EQBOUND01 (MHP 7/2003)
! and EQBOUND06 (LLP 10/16/06). In 2026: the three vintages' code was
! written once, driven by opal_layout (the 2001/2006 codes first, the
! 1995 code with them later); the density searches' bracket and their
! tolerance, once loosened to 1e-5 to hide a cv slot rescaled from stale
! values on every trial call, were fixed (the trial calls rescale only
! the slots they interpolate, as the original did not); the 2001 search
! the same (M41).

end module opal_eos
