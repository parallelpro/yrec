!----------------------------------------------------------------------
! opal_eos
!----------------------------------------------------------------------
! The OPAL equation of state: the entry points for all three vintages,
! and the 2001 and 2006 tables, which share one state layout
! (opal_eos_vintage, their instances opal_tables(iv_opal01:iv_opal06);
! the 1995 tables, with their own layout, are in opal95_eos).
!     compute_opal_eos     density and the thermodynamic derivatives of
!                          one vintage at (log T, log P, X, Z)
!     find_opal_edge_ramp  whether the vintage's table holds the point,
!                          and the OPAL weight near its edge
! Under these, for 2001 and 2006: opal_density_search (OPAL's rhoofp01
! and rhoofp06), interpolate_opal_eos (OPAL's esac01 and esac06, which
! remain its entry points, public for blended_eos's gamma1_at_density and test_eos),
! t6rinteos, add_radiation (the radiation correction), readcoeos (the
! table reader, run once per vintage), and the edge ramps
! find_table_edge_ramp (2001) and find_ragged_edge_ramp (2006). This
! code is written once: what the vintages differ in (the table's size
! and grids, the slots, the radiation constants, the reader's row
! format, how the density search treats radiation, which edge the table
! has) is their row of opal_layout (opal_eos_tables); compute_opal_eos
! branches only between the 1995 code and the rest. The report_*
! procedures keep each vintage's messages in its own words. The 2001
! code follows Rogers & Nayfonov (2002); the 2006 tables are the 2005
! update.
module opal_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use opal95_eos, only: esac, rhoofp, find_opal95_ramp
      use opal_eos_tables, only: find_opal_ramp, gmass01, i_opal_p, open_opal_eos_table, &
           interpolate_opal_t6_rho, iv_opal95, iv_opal01, iv_opal06, &
           n_eos_mv, n_eos_mx, opal95, opal_tables, opal_layout, opal_shared, &
           opal_eos_vintage, opal_flag_set, opal_rho_not_found, opal_vintage_layout, quad, &
           opal_radiation_constant_over_3
      use physical_constants, only: gas_constant
      use eos_messages, only: eos_log, eos_report, message_length, unwritten_line
      use eos_mixture, only: compute_mean_weights
      use eos_configuration, only: eos_settings, eos_table_paths_t, opal_vintage_none
      use eos_table_source, only: eos_table_source_t, stage_over_the_gas, adopt_table_result, &
           blend_table_result
      use math_lib
      implicit none
      private
      public :: compute_opal_eos, find_opal_edge_ramp, esac06, opal_table_source
! for test_eos: the 2001 interpolator's entry point
      public :: esac01
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

subroutine opal_edge_weight(temperature, log10_density, weight, in_table, near_edge, ierr)
      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: weight
      logical, intent(out) :: in_table, near_edge
      integer, intent(out) :: ierr
      call find_opal_edge_ramp(eos_settings%opal_vintage, temperature, log10_density, weight, &
           in_table, near_edge, ierr)
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

! the density that gives P (the 1995 tables have their own state layout
! and code; 2001 and 2006 one search, driven by the vintage's row)
      if (vintage == iv_opal95) then
         density_cgs = rhoofp(opal95, hydrogen_fraction_work, t_million_k, p_e12, rad_flag, ierr)
      else
         density_cgs = opal_density_search(opal_tables(vintage), vintage, hydrogen_fraction_work, &
              t_million_k, p_e12, rad_flag, ierr)
      end if
      if (ierr /= 0) return
      if (density_cgs <= -998.0d0) then
         in_table = .false.
         return
      end if
      r%density = density_cgs
      r%log10_density = log10(r%density)

! the thermodynamic quantities there, into the vintage's eos_output
      if (vintage == iv_opal95) then
         call esac(opal95, hydrogen_fraction_work, t_million_k, density_cgs, &
              deriv_order, rad_flag, ierr, in_table)
      else
         call interpolate_opal_eos(opal_tables(vintage), vintage, hydrogen_fraction_work, &
              t_million_k, density_cgs, deriv_order, rad_flag, ierr, in_table)
      end if
      if (.not. in_table) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'WARNNING... OPAL TBL FAIL'
         call eos_log(eos_lines)
         return
      end if
      if (ierr /= 0) return

! dln P/dln rho, dln P/dln T6, cv, gamma1 and gamma2/(gamma2-1) from the
! vintage's results, in its slots
      if (vintage == iv_opal95) then
         output = opal95%eos_output
      else
         output = opal_tables(vintage)%eos_output
      end if
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
           metal_fraction, specific_gas_constant_check, &
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
! OPAL result near the table's edge: the 1995 table (find_opal95_ramp),
! the 2006 table's ragged edge (find_ragged_edge_ramp), or the 2001
! table's one edge per T6 row (find_table_edge_ramp).
subroutine find_opal_edge_ramp(vintage, temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
      integer, intent(out) :: ierr

      ierr = 0
      if (vintage == iv_opal95) then
         call find_opal95_ramp(temperature, log10_density, ramp_factor, &
              in_opal_table, needs_ramp, ierr)
      else if (opal_layout(vintage)%ragged_edge) then
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
! MHP 7/2003: whether a (T, log10 rho) point falls within the OPAL 2001
! EOS table, and near its edge the ramp factor that blends the OPAL
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
! opal_density_search (formerly rhoofp01 and rhoofp06)
!----------------------------------------------------------------------
! The density at which the 2001 or 2006 table gives the pressure P (in
! Mb) at (X, T6), by regula falsi through the interpolator's pressure
! alone (trial calls, deriv_order = 1), from a bracket: a first density
! from the table's highest one scaled by P over the pressure there, then
! the lower density shrunk by 0.2 until its pressure is below the target
! (at most 20 times, down to the table's 1e-14), or the upper one grown
! by 5 (to the table's edge); at most 30 steps to |dP/P| < 0.5e-7. The
! vintage's row of opal_layout says how radiation enters (2006 subtracts
! its pressure from the target and asks for the gas alone; 2001 adds it
! to the table's bounds and asks for gas plus radiation), the highest
! density at each T6 row, a/c, and the name in messages. Returns
! opal_rho_not_found (-999) when the point is not in the table, and,
! with a run-log message, when the search does not converge.
! (2026 bug sweep, batch 2: the bracket and the tolerance, once loosened
! to 1e-5 to hide a cv slot rescaled from stale values on every trial
! call, fixed in batch 0; 2026-10: the 2001 search the same, M41.)
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
! interpolate_opal_eos (formerly esac01 and esac06)
!----------------------------------------------------------------------
! The OPAL 2001 and 2006 interpolation (OPAL's esac01 and esac06, which
! are its entry points): the EOS variables at (X, T6, rho) by quadratic
! interpolation in rho at fixed T6, then in T6 (t6rinteos), blended
! across overlapping quadratics, then in X; see esac in opal95_eos for
! the scheme. Both tables carry 9 variables. Unlike the 1995 code, every
! stencil is clamped to the table. What the vintages differ in (the
! table's size, the slots, the radiation constants) is their row of
! opal_layout, and their messages are written by the report_*
! procedures below.
!
! With the radiation correction (rad_flag = 1) a full call (deriv_order
! = the table's n_vars) makes every slot that of gas plus radiation
! (add_radiation); a density search's trial call (deriv_order = 1) only
! adds the radiation pressure to its pressure slot, since its other
! slots are an earlier call's. Without it (rad_flag = 0, the 2006
! density search's trial calls) the slots are the gas's, the cv slot in
! the table's units.
!
! The 2006 slots (the 2001 table lacks slot 4):
!   v%eos_output(1) pressure, in megabars (1e12 dyne/cm**2)
!   v%eos_output(2) energy, in 1e12 erg/gm (zero at T6=0)
!   v%eos_output(3) entropy, in units of energy/T6
!   v%eos_output(4) dE/dRho at constant T6
!   v%eos_output(5) specific heat, dE/dT6 at constant volume
!   v%eos_output(6) dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   v%eos_output(7) dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   v%eos_output(8) gamma1 (Cox & Guili eq 9.88)
!   v%eos_output(9) gamma2/(gamma2-1) (Cox & Guili eq 9.88)
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
! the vintage's table: its size, its X grid and its name in messages
      integer :: nr, nt
      real(dp) :: x_grid(n_eos_mx)
      character(len=15) :: routine_id

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
      character(len=message_length) :: eos_lines(4)

      nr = opal_layout(vintage)%nr
      nt = opal_layout(vintage)%nt
      x_grid = opal_layout(vintage)%x_grid
      routine_id = opal_layout(vintage)%routine_id
      ierr = 0
      in_table = .true.

      if (deriv_order > 9) then
         eos_lines = unwritten_line
         write(eos_lines,'(A, " iorder cannot exceed 9")') routine_id
         call eos_log(eos_lines)
      end if
      if ((rad_flag /= 0) .and. (rad_flag /= 1)) then
         eos_lines = unwritten_line
         write(eos_lines,'(A, " Irad must be 0 or 1")') routine_id
         call eos_report(eos_lines)
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
            eos_lines = unwritten_line
            write(eos_lines,*) routine_id, "Mass fractions exceed unity (61)"
            call eos_report(eos_lines)
            eos_lines = unwritten_line
            write(eos_lines,*) 'Z, XH', opal_shared%table_metal_fraction(vintage), hydrogen_fraction
            call eos_report(eos_lines)
            
            
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
      table_sum_2x2 = 0.0d0
      table_sum_3x3 = 0.0d0
      table_sum_3x4 = 0.0d0
      table_sum_4x4 = 0.0d0
      do x_loop_index = v%x_index_lo, x_index_hi
         do density_scan_idx = v%density_index_1, v%density_index_1+1
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+1
               table_sum_2x2 = table_sum_2x2 + &
                    v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         do density_scan_idx = v%density_index_1, v%density_index_1+2
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+2
               table_sum_3x3 = table_sum_3x3 + &
                    v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         if (t6_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+2
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  table_sum_3x4 = table_sum_3x4 + &
                       v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_3x4 = 2.0d+30
         end if
         if (density_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+3
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  table_sum_4x4 = table_sum_4x4 + &
                       v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_4x4 = 2.0d+30
         end if
      end do
      v%density_interp_order = 2
      v%t6_interp_order = 2
      if (table_sum_3x3 > 1.0d+30) then
         if (table_sum_2x2 < 1.0d+25) then
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
      if (table_sum_3x4 < 1.0d+30) v%t6_interp_order = 3
      if (table_sum_4x4 < 1.0d+30) v%density_interp_order = 3

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
      recompute_flag = 0
      cache_slot = 1
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            if (x_index_hi == 1) then
               v%x_interp_result(t6_scan_idx,density_scan_idx) = &
                    v%x_interp_workspace(v%x_index_lo,t6_scan_idx,density_scan_idx)
               cycle
            end if
            v%x_interp_result(t6_scan_idx,density_scan_idx) = &
                 quad(v%quad, recompute_flag, cache_slot, hydrogen_fraction, &
                 v%x_interp_workspace(v%x_index_lo,t6_scan_idx,density_scan_idx), &
                 v%x_interp_workspace(x_index_2,t6_scan_idx,density_scan_idx), &
                 v%x_interp_workspace(x_index_3,t6_scan_idx,density_scan_idx), &
                 v%x_grid_copy(v%x_index_lo), v%x_grid_copy(x_index_2), &
                 v%x_grid_copy(x_index_3))
            if (v%x_interp_result(t6_scan_idx,density_scan_idx) > 1.0d+20) &
                 call report_x_interpolation(vintage, v, t6_scan_idx, density_scan_idx, &
                 v%x_index_lo)
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
            if (v%x_interp_result_alt(t6_scan_idx,density_scan_idx) > 1.0d+20) &
                 call report_x_interpolation(vintage, v, t6_scan_idx, density_scan_idx, &
                 x_index_2)
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
!
!
      call t6rinteos(v, vintage, density_value, t6_value, ierr)
      if (ierr /= 0) return
      v%eos_output(eos_var_idx) = v%esact
      end do
      if (eos_var_idx > deriv_order) then
      pressure_scale = t6_temperature*density
      lay = opal_layout(vintage)
      v%eos_output(lay%i_p) = v%eos_output(lay%i_p)*pressure_scale ! interpolated in p/po
! Only slots that this call actually re-interpolated (index <= deriv_order)
! are rescaled; the original scaled E and cv unconditionally, so the
! deriv_order=1 trial calls from the rho(P,T) inversion compounded the
! stale E and cv slots by T6 and moles*R/mu on every call.
      if (lay%i_e <= deriv_order) &
           v%eos_output(lay%i_e) = v%eos_output(lay%i_e)*t6_temperature ! interpolated in E/T6
      if (rad_flag == 1 .and. deriv_order >= lay%n_vars) then
         mean_molecular_weight = gmass01(hydrogen_fraction, opal_shared%table_metal_fraction(vintage), &
              total_moles, ground_state_energy, metal_mole_fraction, &
              species_number_fraction)
         call add_radiation(v, lay, t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else if (rad_flag == 1) then
! a density search's trial call (only the pressure slot is fresh, the
! others hold an earlier call's values): the radiation pressure is added
! to the pressure alone, as add_radiation would
         v%eos_output(lay%i_p) = v%eos_output(lay%i_p) + &
              opal_radiation_pressure(lay, t6_temperature)
      end if
      return
      end if
      call report_run_z(vintage, v, hydrogen_fraction, t6_temperature, density)
      ierr = 1

end subroutine interpolate_opal_eos

!----------------------------------------------------------------------
! esac01, esac06
!----------------------------------------------------------------------
! OPAL's entry points for the 2001 and 2006 tables (interpolate_opal_eos).
subroutine esac01(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)
      type(opal_eos_vintage), intent(inout) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, density
      integer, intent(in) :: deriv_order, rad_flag
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
      call interpolate_opal_eos(v, iv_opal01, hydrogen_fraction, t6_temperature, &
           density, deriv_order, rad_flag, ierr, in_table)
end subroutine esac01

subroutine esac06(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)
      type(opal_eos_vintage), intent(inout) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, density
      integer, intent(in) :: deriv_order, rad_flag
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
      call interpolate_opal_eos(v, iv_opal06, hydrogen_fraction, t6_temperature, &
           density, deriv_order, rad_flag, ierr, in_table)
end subroutine esac06

!----------------------------------------------------------------------
! readcoeos (formerly readcoeos01 and readcoeos06)
!----------------------------------------------------------------------
! Reads a 2001 or 2006 OPAL EOS table (one block per tabulated X) into
! v, then derives the grid spacings interpolate_opal_eos uses and the
! table's edge in density at each T6 row (find_opal_ramp, 2001). Called
! once per vintage from interpolate_opal_eos, on the unit eos_init's
! open_opal_eos_table opened. The vintage's row of opal_layout gives the
! table's size and its row format; the header's moles per gram, mean
! molecular weight and species fractions, and each row's log10 R (2001)
! or amu and log10 n_e (2006), are read and not kept.
subroutine readcoeos(v, vintage, ierr)
      integer, intent(in) :: vintage

      type(opal_eos_vintage), intent(inout) :: v
      integer, intent(out) :: ierr

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv
! the vintage's table size
      integer :: nr, nt

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      integer :: column_idx, table_var_idx
      real(dp) :: unused_field, unused_header(2), unused_species(6), unused_columns(2)
      integer :: ios
      logical :: layout_error

      nr = opal_layout(vintage)%nr
      nt = opal_layout(vintage)%nt
      ierr = 0
      if (opal_shared%table_unit == 0) then
         call report_no_table(vintage)
         ierr = 1
         return
      end if

      blank_line = ' '

      if (opal_shared%readco_init_flag(vintage) /= opal_flag_set) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     v%eos_table(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         do var_idx = 1, mv
            v%eos_output(var_idx) = 1.0d0
         end do
         opal_shared%readco_init_flag(vintage) = opal_flag_set
      end if

! .....read  tables (on the unit open_opal_eos_table opened); a read
! that fails or a block out of its layout ends the read with ierr = 1,
! the unit closed
      call read_table_blocks(ios, layout_error)
      if (ios /= 0 .or. layout_error) then
         if (ios /= 0) call report_bad_table(vintage, x_loop_index, density_row)
         close (opal_shared%table_unit)
         opal_shared%table_unit = 0
         ierr = 1
         return
      end if

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx) == 0.0d0) then
            call report_zero_t6(vintage, t6_scan_idx)
            close (opal_shared%table_unit)
            opal_shared%table_unit = 0
            ierr = 1
            return
         end if
         v%t6_grid(t6_scan_idx) = v%t6_list(1,t6_scan_idx)
      end do
      do t6_idx = 2, nt
         v%t6_grid_spacing_inv(t6_idx) = 1.0d0/(v%t6_grid(t6_idx) - &
              v%t6_grid(t6_idx-1))
      end do
      v%density_grid(1) = v%density_grid_table(1,1)
      do r_idx = 2, nr
         v%density_grid(r_idx) = v%density_grid_table(1,r_idx)
         v%density_grid_spacing_inv(r_idx) = 1.0d0/(v%density_grid(r_idx) - &
              v%density_grid(r_idx-1))
      end do
      do x_idx = 2, mx
         v%x_grid_spacing_inv(x_idx) = 1.0d0/(v%x_grid_copy(x_idx) - v%x_grid_copy(x_idx-1))
      end do
      close (opal_shared%table_unit)
      opal_shared%table_unit = 0
! MHP 7/2003 added ramp between OPAL and other EOS: the edge of the
! table at high rho, fixed T (find_opal_ramp reads it for 2001)
      v%t_row_index = 1
      do t6_idx = 1, nt
         v%density_edge_at_t(t6_idx) = v%density_grid(opal_layout(vintage)%density_index_edge(t6_idx))
      end do

contains

! Each X block of the table: its header, then for each density column
! its record (which must be the column's number, with at least the
! T6 rows the table needs) and its T6 rows in the vintage's row format.
subroutine read_table_blocks(ios, layout_error)
      integer, intent(out) :: ios
      logical, intent(out) :: layout_error
      integer :: u
      u = opal_shared%table_unit
      layout_error = .false.
      density_row = 0
      do x_loop_index = 1, mx
         read (u,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)',iostat=ios) &
              v%hydrogen_fraction_header(x_loop_index), &
              v%z_table(x_loop_index), unused_header(1), unused_header(2)
         if (ios /= 0) return
         read (u,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)',iostat=ios) unused_species
         if (ios /= 0) return
         read (u,'(A)',iostat=ios) blank_line
         if (ios /= 0) return
         do density_row = 1, nr
            read (u,'(2I5,2F12.7,17X,E15.7)',iostat=ios) record_number, &
                 v%temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index,density_row)
            if (ios /= 0) return
            if (record_number /= density_row) then
               call report_record_number(vintage, record_number, density_row)
               layout_error = .true.
               return
            end if
            read(u,'(A)',iostat=ios) blank_line
            if (ios == 0) read(u,'(A)',iostat=ios) blank_line
            if (ios /= 0) return
            if (v%temperature_count_used(x_loop_index,density_row) < &
                 opal_layout(vintage)%t6_index_lo(density_row)) then
               call report_short_column(vintage, v%hydrogen_fraction_header(x_loop_index), &
                    v%density_grid_table(x_loop_index,density_row))
               layout_error = .true.
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index,density_row)
               if (t6_row > opal_layout(vintage)%t6_index_lo(density_row)) then
                  read (u,'(A)',iostat=ios) blank_line
                  if (ios /= 0) return
                  cycle
               end if
               read (u, opal_layout(vintage)%row_format,iostat=ios) &
                    v%t6_list(density_row,t6_row), &
                    (unused_columns(column_idx), column_idx=1,opal_layout(vintage)%n_skipped_columns), &
                    (v%eos_table(x_loop_index,table_var_idx,t6_row,density_row), &
                    table_var_idx=1,9)
               if (ios /= 0) return
            end do
            read(u,'(A)',iostat=ios) blank_line
            if (ios == 0) read(u,'(A)',iostat=ios) blank_line
            if (ios == 0) read(u,'(A)',iostat=ios) blank_line
            if (ios /= 0) return
         end do
         read(u,'(A)',iostat=ios) blank_line
         if (ios /= 0) return
      end do
end subroutine read_table_blocks

end subroutine readcoeos

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
      if (v%esact > 1.0d+15) then
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
! The 2001 and 2006 interpolators' messages, each vintage's in its
! original words (run log).
!----------------------------------------------------------------------
! (T6, rho) outside the table
subroutine report_outside_table(vintage, t6, density, t6_first, t6_last, &
     density_first, density_last)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: t6, density, t6_first, t6_last, density_first, density_last
      character(len=message_length) :: eos_lines(4)
      character(len=15) :: id
      id = opal_layout(vintage)%routine_id
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,*) id, " T6 or rho outside of table range (62)"
         call eos_log(eos_lines)
         eos_lines = unwritten_line
         write(eos_lines,*) "t6, t6a(1),t6a(nt):", t6, t6_first, t6_last
         call eos_log(eos_lines)
         eos_lines = unwritten_line
         write(eos_lines,*) "slr,r,rho(1),rho(nr):", density, density, density_first, density_last
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,*) id, " T6/LogR outside of table range (62)"
         call eos_log(eos_lines)
         eos_lines = unwritten_line
         write(eos_lines,*) "slt, t6a(1),t6a(nt):", t6, t6_first, t6_last
         call eos_log(eos_lines)
         eos_lines = unwritten_line
         write(eos_lines,*) "slr, rho(1), rho(nr):", density, density_first, density_last
         call eos_log(eos_lines)
      end if
end subroutine report_outside_table

! the T6 search found no row (the 2001 message never printed its indices)
subroutine report_t6_search(vintage, hi_idx, lo_idx, mid_idx)
      integer, intent(in) :: vintage, hi_idx, lo_idx, mid_idx
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'(A, " ihi,ilo,imd",3I5)') opal_layout(vintage)%routine_id
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'(A, " ihi,ilo,imd",3I5)') opal_layout(vintage)%routine_id, hi_idx, lo_idx, mid_idx
         call eos_log(eos_lines)
      end if
end subroutine report_t6_search

! (T6, rho) in an empty region of the table
subroutine report_empty_region(vintage, hydrogen_fraction, t6, density)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: hydrogen_fraction, t6, density
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,*) opal_layout(vintage)%routine_id, "T6/rho in empty region of OPAL 2001 EOS", &
              " table (65)"
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,*) opal_layout(vintage)%routine_id, "T6/log rho in empty region of table (65)"
         call eos_log(eos_lines)
      end if
      eos_lines = unwritten_line
      write(eos_lines,'("xh,t6,r=", 3E12.4)') hydrogen_fraction, t6, density
      call eos_log(eos_lines)
end subroutine report_empty_region

! the X tables of the stencil are not all for one Z
subroutine report_tables_z(vintage)
      integer, intent(in) :: vintage
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'(A,"Z does not match Z in OEOS01 files you are" ," using")') opal_layout(vintage)%routine_id
         call eos_report(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'(A,"Z does not match Z in OEOS06 files you are" ," using")') opal_layout(vintage)%routine_id
         call eos_report(eos_lines)
      end if
end subroutine report_tables_z

! an X-interpolated value above 1e20 (an empty entry in the stencil):
! the stencil and the three table values from x_first
subroutine report_x_interpolation(vintage, v, t6_idx, density_idx, x_first)
      integer, intent(in) :: vintage, t6_idx, density_idx, x_first
      type(opal_eos_vintage), intent(in) :: v
      character(len=message_length) :: eos_lines(4)
      integer :: k
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'(A," problem it ir,l3,k3,iq,ip=", 6I5)') opal_layout(vintage)%routine_id, &
              t6_idx, density_idx, v%density_index_3, v%t6_index_3, v%density_interp_order, &
              v%t6_interp_order
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') opal_layout(vintage)%routine_id, &
              t6_idx, density_idx, v%density_index_3, v%t6_index_3, v%density_interp_order, &
              v%t6_interp_order
         call eos_log(eos_lines)
      end if
      eos_lines = unwritten_line
      write(eos_lines,'(3E12.4)') (v%x_interp_workspace(k,t6_idx,density_idx), k=x_first,x_first+2)
      call eos_log(eos_lines)
end subroutine report_x_interpolation

! the run's Z differs from the table's
subroutine report_run_z(vintage, v, hydrogen_fraction, t6, density)
      integer, intent(in) :: vintage
      type(opal_eos_vintage), intent(in) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6, density
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,*) opal_layout(vintage)%routine_id, " Z does not match Z in OPAL 2001 EOS files", &
              " you are using (66)"
         call eos_report(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,*) opal_layout(vintage)%routine_id, " Z does not match Z in EOSdata* files you ", &
              "are using (66)"
         call eos_report(eos_lines)
      end if
      eos_lines = unwritten_line
      write(eos_lines,'("mf,zz(mf)=",I5,E12.4)') v%x_index_lo, v%z_table(v%x_index_lo)
      call eos_report(eos_lines)
      eos_lines = unwritten_line
      write(eos_lines,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') v%t6_interp_order, &
           v%density_interp_order, v%t6_index_3, v%density_index_3, hydrogen_fraction, &
           t6, density, opal_shared%table_metal_fraction(vintage)
      call eos_report(eos_lines)
end subroutine report_run_z

! the reader: a density block out of order
subroutine report_record_number(vintage, record_number, density_row)
      integer, intent(in) :: vintage, record_number, density_row
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'(" Data file incorrect: numtot,jcs= ",2I5)') record_number, density_row
         call eos_report(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'(" OEOS06 Data file incorrect: numtot,jcs= " ,2I5)') record_number, density_row
         call eos_report(eos_lines)
      end if
end subroutine report_record_number

! the reader: a density column with fewer T6 rows than the table needs
subroutine report_short_column(vintage, hydrogen_fraction, density)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: hydrogen_fraction, density
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'("problem with data files: X=",F6.4, " density=",E14.4)') hydrogen_fraction, density
         call eos_report(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'("Problem with OEOS96 data files: X=",F6.4, " density=",E14.4)') hydrogen_fraction, &
              density
         call eos_report(eos_lines)
      end if
end subroutine report_short_column

! the reader: no table open (eos_init opens the selected vintage's; a
! vintage installed without it, or a table read before, has none)
subroutine report_no_table(vintage)
      integer, intent(in) :: vintage
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a)') trim(opal_layout(vintage)%routine_id), &
           ' no OPAL EOS table is open to read (eos_init opens it)'
      call eos_report(eos_lines)
end subroutine report_no_table

! the reader: a read that failed (the table ends early or is corrupt)
subroutine report_bad_table(vintage, x_block, density_row)
      integer, intent(in) :: vintage, x_block, density_row
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0,a,i0)') trim(opal_layout(vintage)%routine_id), &
           ' the OPAL EOS table ends early or is corrupt: X block ', x_block, ', density column ', &
           density_row
      call eos_report(eos_lines)
end subroutine report_bad_table

! the reader: a zero T6
subroutine report_zero_t6(vintage, t6_idx)
      integer, intent(in) :: vintage, t6_idx
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'("READCOEOS01: Error:",I4, "-th T6 value is zero")') t6_idx
         call eos_report(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'("READCOEOS06: Error:",I4, "-th T6 value is zero")') t6_idx
         call eos_report(eos_lines)
      end if
end subroutine report_zero_t6

! the T6-rho interpolation's indices out of range
subroutine report_interpolation_indices(vintage)
      integer, intent(in) :: vintage
      character(len=message_length) :: eos_lines(4)
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'("Interpolation indices out of range", ";please report conditions.")')
         call eos_report(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'("T6RINTEOS06: Interpolation indices out", " of range;please report conditions.")')
         call eos_report(eos_lines)
      end if
end subroutine report_interpolation_indices

end module opal_eos
