!----------------------------------------------------------------------
! opal2006_eos
!----------------------------------------------------------------------
! The OPAL 2006 equation of state (the 2005 update of Rogers &
! Nayfonov), with the same two entries as opal95_eos:
!     compute_opal2006_eos  density and the thermodynamic derivatives
!     find_opal2006_ramp    table coverage and the edge-blend weight
! Under these, OPAL's routines under their original names (rhoofp06,
! esac06, t6rinteos06, radsub06, readcoeos06); esac06 is also public
! for eos_get_gamma1 and test_eos. The tables and the interpolation
! state live in opal_eos_tables (opal06).
module opal2006_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opal_eos_tables, only: gmass01, i_opal_chi_rho, i_opal_chi_t, &
           i_opal_cv, i_opal_dedrho, i_opal_e, i_opal_gamma1, &
           i_opal_gamma2_ratio, i_opal_p, i_opal_s, interpolate_opal_t6_rho, &
           iv_opal06, n_eos06_nr, n_eos06_nt, n_eos_mv, n_eos_mx, opal06, &
           opal06_density_index_edge, opal06_t6_index_lo, opal06_x_grid, &
           opal_eos, opal_eos_var_order, opal_eos_vintage, opal_flag_set, &
           opal_rho_not_found, quad
      use physical_constants, only: gas_constant
      use io_units, only: run_log_unit
      use eos_mixture, only: compute_mean_weights
      use math_lib
      implicit none
      private
      public :: compute_opal2006_eos, find_opal2006_ramp, esac06

contains

!----------------------------------------------------------------------
! compute_opal2006_eos (formerly oeqos06)
!----------------------------------------------------------------------
! 2006 OPAL equation of state. LLP, October 17, 2006. Given log10(T),
! log10(P), X, and Z, looks up density and the thermodynamic
! derivatives needed by the rest of the EOS machinery from the OPAL
! 2006 tables (via esac06). Called from blend_eos
! (compute_blended_eos); in_table is false when the point is
! outside the OPAL tables.
subroutine compute_opal2006_eos(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, dlnrho_dlnt, &
     dlnrho_dlnp, specific_heat_cp, adiabatic_gradient, ierr, in_table)


      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(out) :: temperature, pressure, &
           log10_density, density
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: beta, beta_inverse, beta14, &
           specific_gas_constant, ion_mean_weight_inverse, &
           electron_mean_weight_inverse, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient

      real(dp) :: t_million_k, p_e12
      real(dp) :: hydrogen_fraction_work
      real(dp) :: density_cgs
      real(dp) :: specific_gas_constant_check
      integer :: rad_flag, deriv_order

      integer, intent(out) :: ierr
! in_table is false outside the table's T range, when no density
! matches P, or when esac06 finds the point outside the table
      logical, intent(out) :: in_table

      ierr = 0
      in_table = .true.

      deriv_order = 9  ! gives all 1st and 2nd order data. See instructions
!                  in esac06.
!     NOTE: rad_flag=0 does not add radiation; rad_flag=1 adds radiation
      rad_flag = 1     ! does add radiation  corrections

      temperature = exp10(log10_temperature)
      pressure = exp10(log10_pressure)
      t_million_k = temperature/1.0d6
      p_e12 = pressure/1.0d12
      if (t_million_k < 0.001870d0 .or. t_million_k > 200.0d0) then
         in_table = .false.
         return
      end if
      hydrogen_fraction_work = hydrogen_fraction

      density_cgs = rhoofp06(opal06, hydrogen_fraction_work, t_million_k, p_e12, &
           rad_flag, ierr)
      if (ierr /= 0) return
      if (density_cgs <= -998.0d0) then
         in_table = .false.
         return
      end if
      density = density_cgs
      log10_density = log10(density)

      call esac06(opal06, hydrogen_fraction_work, t_million_k, density_cgs, &
           deriv_order, rad_flag, ierr, in_table)
      if (.not. in_table) then
         write(run_log_unit, *) 'WARNNING... OPAL TBL FAIL'
         return
      end if
      if (ierr /= 0) return


      dlnrho_dlnp = 1.0d0/opal06%eos_output(i_opal_chi_rho)               ! O2006 EOS(6) is dlogP/dlogRho at const T6
      dlnrho_dlnt = -opal06%eos_output(i_opal_chi_t)/opal06%eos_output(i_opal_chi_rho)       ! O2006 EOS(7) is dlogp/dlogT6 at const Rho
      specific_heat_cp = 1.0d6*opal06%eos_output(i_opal_cv)*opal06%eos_output(i_opal_gamma1)/opal06%eos_output(i_opal_chi_rho)
                                      ! O2006 EOS(5) is the specific heat. dE/dT6
                                      !              at const Vol
                                      ! O2006 EOS(8) is gamma1
      adiabatic_gradient = 1.0d0/opal06%eos_output(i_opal_gamma2_ratio)         ! O2006 EOS(9) is gamma2/(gamma2-1)

      beta14 = (2.521971383d-3*t_million_k*t_million_k)* &
           (t_million_k*t_million_k/p_e12)
      beta = 1.0d0 - beta14
      beta_inverse = 1.0d0/beta
      specific_gas_constant = pressure*beta/(density*temperature)
      call compute_mean_weights(temperature, pressure, density, hydrogen_fraction, &
           metal_fraction, specific_gas_constant_check, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, beta)
      if (electron_mean_weight_inverse <= 0.0d0) then
      electron_mean_weight_inverse = 0.0d0
      ion_mean_weight_inverse = specific_gas_constant/gas_constant - &
           electron_mean_weight_inverse
      end if

end subroutine compute_opal2006_eos

!----------------------------------------------------------------------
! find_opal2006_ramp (formerly eqbound06)
!----------------------------------------------------------------------
! LLP 10/16/06: checks whether a (T, log10 rho) point falls within the
! OPAL 2006 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (see blend_eos.f90's use_opal2006_eos branch, which calls this
! immediately after compute_opal2006_eos). Unlike find_opal95_ramp/find_opal2001_ramp (1995/2001),
! this version tracks the table's ragged edge via per-row/per-column
! index arrays (opal06_t6_index_lo/opal06_density_index_edge) rather than a single
! linear-interpolated edge value.
pure subroutine find_opal2006_ramp(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp)


      integer, parameter :: nr = n_eos06_nr, nt = n_eos06_nt

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
      if ((density < opal06%density_grid(1)) .or. (density >= opal06%density_grid(nr))) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if

!     Exit if outside table in T6
      if ((t6 > opal06%t6_grid(1)) .or. (t6 <= opal06%t6_grid(nt))) then
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
      if (t6 > opal06%t6_grid(2)) then
          t6_ramp_factor = (t6-opal06%t6_grid(1))/(opal06%t6_grid(2)-opal06%t6_grid(1)) ! Ramp in temperature
          needs_ramp = .true.
      end if
      if (density < opal06%density_grid(2)) then
          density_ramp_factor = (density - opal06%density_grid(1))/ &
               (opal06%density_grid(2)-opal06%density_grid(1)) ! Ramp in density
          needs_ramp = .true.
      end if

!     Find a t6_row such that opal06%t6_grid(t6_row-1) >= T6 > opal06%t6_grid(t6_row)
!     and a density_row such that opal06%density_grid(density_row-1) < D <=
!     opal06%density_grid(density_row)
!     The ESAC06 variables (k=)opal06%t6_index_3 and (l-)opal06%density_index_2 are
!     close to what we need. A linear search will work fine.
      t6_row = opal06%t6_index_3
   do while (t6 <= opal06%t6_grid(t6_row))
        t6_row = t6_row+1
   end do
   do while (t6 > opal06%t6_grid(t6_row-1))
        t6_row = t6_row-1
   end do
!     We now have: opal06%t6_grid(t6_row-1) >= T6 > opal06%t6_grid(t6_row)

      density_row = opal06%density_index_2
   do while (density > opal06%density_grid(density_row))
        density_row = density_row+1
   end do
   do while (density <= opal06%density_grid(density_row-1))
        density_row = density_row-1
   end do
!     We now have: opal06%density_grid(density_row-1) < D <= opal06%density_grid(density_row)

!     For a given temperature in array opal06%t6_grid with index t6_row, e.g.,
!     opal06%t6_grid(t6_row), element opal06_density_index_edge(t6_row) of array
!     opal06_density_index_edge contains the index to the max allowed density
!     in array opal06%density_grid. So, for a given opal06%t6_grid(t6_row), the
!     associated D must be less than opal06%density_grid(opal06_density_index_edge(t6_row)).
!     A valid index in array opal06%density_grid must be less than or equal
!     to opal06_density_index_edge(t6_row)

      if (density_row > opal06_density_index_edge(t6_row)) then
        in_opal_table = .false.
        needs_ramp = .true.
        ramp_factor = 0d0
        return
      end if
      if (density_row == opal06_density_index_edge(t6_row)) then
        needs_ramp = .true.
        density_ramp_factor = (opal06%density_grid(density_row)-density)/ &
             (opal06%density_grid(density_row)-opal06%density_grid(density_row-1))
      end if

!     For a given density in array opal06%density_grid with index density_row,
!     i.e., opal06%density_grid(density_row), element opal06_t6_index_lo(density_row)
!     of array opal06_t6_index_lo contains the index to the min allowed
!     temperature in array opal06%t6_grid. So for a given
!     opal06%density_grid(density_row), the associated T6 must be greater
!     than opal06%t6_grid(opal06_t6_index_lo(density_row)). (Note that opal06%t6_grid is a
!     decreasing array.) A valid index in opal06%t6_grid must be less than or
!     equal to opal06_t6_index_lo(density_row).

      if (t6_row > opal06_t6_index_lo(density_row)) then
        in_opal_table = .false.
        needs_ramp = .true.
        ramp_factor = 0d0
        return
      end if
      if (t6_row == opal06_t6_index_lo(density_row)) then
        needs_ramp = .true.
        t6_ramp_factor = (t6-opal06%t6_grid(t6_row))/(opal06%t6_grid(t6_row-1)-opal06%t6_grid(t6_row))
      end if

      ramp_factor = density_ramp_factor * t6_ramp_factor
end subroutine find_opal2006_ramp

!----------------------------------------------------------------------
! rhoofp06
!----------------------------------------------------------------------
! OPAL 2006 EOS analogue of rhoofp01 (see there and rhoofp
! for the general description), with two preserved differences:
! (1) pressure_max/pressure_min here do NOT add a radiation term
! (rhoofp01 both add rad_flag*4/3*rat*t6**4); (2) the
! trial esac06 calls below pass a hardcoded 0 for the radiation
! flag (not rad_flag); since Batch 3 so does the table-priming call
! (see the note at that call).
!
! 2026 (bugsweep Batch 2): the convergence tolerance is back at the
! original 0.5d-7 (same as rhoofp/rhoofp01). It had been loosened to
! 1.0d-5 in 2025 to stop v%eos_output(5) from "growing
! without bound"; the real cause was esac06's derivative tail
! rescaling the cv slot from stale data on every deriv_order=1 trial
! call, fixed in bug-sweep batch 0 (the tail is now gated on
! deriv_order). The loose tolerance only hid the symptom by making
! the inversion converge in fewer trials. Also in this change: the
! lower bracket is verified (repeat the x0.2 shrink until
! P(rho1) <= P_target or the floor is reached -- plain regula falsi
! only converges from a valid bracket), the refinement cap is raised
! from 11 to 30 (regula falsi is linear once an endpoint sticks), and
! the non-convergence message is written to the run log instead of
! silently returning -999 (which compute_blended_eos turns into a Yale/SCV
! fallback with ierr = 0).
real(dp) function rhoofp06(v, hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = n_eos_mx, nt = n_eos06_nt

      real(dp) :: rad_const_over_c
! NOTE: no D-suffix in the original (rhoofp06.f) -- preserved
! verbatim: parsed as a single-precision constant and then widened to
! double precision, which is NOT bit-identical to the correctly-
! rounded double value.
      data rad_const_over_c/1.8914785e-3/
! --- locals ---
      real(dp) :: rat, radiation_pressure, pressure_no_rad
      real(dp) :: hydrogen_fraction_dbg, t6_dbg, density_dbg
      integer :: ideriv_dbg
      integer :: lo_idx, hi_idx, mid_idx
      integer :: x_bisect_idx, t6_bisect_idx
      real(dp) :: pressure_max, pressure_min
      real(dp) :: density_trial1, density_trial2, density_trial3
      real(dp) :: pressure_trial1, pressure_trial2, pressure_trial3
      integer :: refine_count, shrink_count
      integer, parameter :: max_refine = 30, max_shrink = 20

      integer, intent(out) :: ierr
      logical :: in_table

      ierr = 0

      rat = rad_const_over_c
      radiation_pressure = 0.0d0
      if (rad_flag == 1) radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_eos%table_loaded_flag(iv_opal06) /= opal_flag_set) then
         hydrogen_fraction_dbg = 0.5d0
         t6_dbg = 1.0d0
         density_dbg = 0.001d0
         ideriv_dbg = 1
! 2026 (bugsweep Batch 3): radiation flag 0 here, not rad_flag. This
! call only loads the tables and its output is discarded; with
! ideriv=1 only the pressure slot is interpolated, so radsub06 would
! divide by a zero cv (0/0 -> NaN in the gamma slots, a trap under
! -ffpe-trap=invalid). The next full esac06 call rewrites every slot.
         call esac06(v, hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              0, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp06 = opal_rho_not_found
            return
         end if
      end if

      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction <= opal06_x_grid(mid_idx)+1.0d-7) then
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

      pressure_max = v%eos_table(x_bisect_idx,1,t6_bisect_idx, &
           opal06_density_index_edge(t6_bisect_idx))*t6_temperature* &
           v%density_grid(opal06_density_index_edge(t6_bisect_idx))
      pressure_min = v%eos_table(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           v%density_grid(1)
      if ((pressure_no_rad > 1.25d0*pressure_max) .or. &
           (pressure_no_rad < pressure_min)) then
!        requested pressure-temperature not in the OPAL 2006 EOS table
         rhoofp06 = opal_rho_not_found
         return
      end if

      density_trial1 = v%density_grid(opal06_density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac06(v, hydrogen_fraction, t6_temperature, density_trial1, 1, 0, ierr, in_table)
      if (ierr /= 0 .or. .not. in_table) then
         rhoofp06 = opal_rho_not_found
         return
      end if
      pressure_trial1 = v%eos_output(i_opal_p)
      if (pressure_trial1 > pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         shrink_count = 0
         shrink: do
            shrink_count = shrink_count + 1
            density_trial1 = 0.2d0*density_trial1
            if (density_trial1 < 1.0d-14) density_trial1 = 1.0d-14
            call esac06(v, hydrogen_fraction, t6_temperature, density_trial1, 1, 0, ierr, in_table)
            if (ierr /= 0 .or. .not. in_table) then
               rhoofp06 = opal_rho_not_found
               return
            end if
            pressure_trial1 = v%eos_output(i_opal_p)
! 2026 (bugsweep Batch 2): keep shrinking until the root is bracketed
! (the original did a single x0.2 step and let regula falsi
! extrapolate from an unbracketed pair).
            if (pressure_trial1 <= pressure_no_rad) exit shrink
            if (density_trial1 <= 1.0d-14 .or. shrink_count >= max_shrink) exit shrink
         end do shrink
      else
         density_trial2 = 5.0d0*density_trial1
         if (density_trial2 > v%density_grid(opal06_density_index_edge(t6_bisect_idx))) &
              density_trial2 = v%density_grid(opal06_density_index_edge(t6_bisect_idx))
         call esac06(v, hydrogen_fraction, t6_temperature, density_trial2, 1, 0, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp06 = opal_rho_not_found
            return
         end if
         pressure_trial2 = v%eos_output(i_opal_p)
      end if

      refine_count = 0
      refine: do
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
      call esac06(v, hydrogen_fraction, t6_temperature, density_trial3, 1, 0, ierr, in_table)
      if (ierr /= 0 .or. .not. in_table) then
         rhoofp06 = opal_rho_not_found
         return
      end if
      pressure_trial3 = v%eos_output(i_opal_p)
! 2026 (bugsweep Batch 2): tolerance restored to the original 0.5d-7
! (was 1.0d-5 from 2025-10-10 to 2026-09; see the header note).
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad) < 0.5d-7) then
         rhoofp06 = density_trial3
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
      write (run_log_unit,'("RHOOFP06: no convergence after ",I0," tries;", &
           & " X, T6, P12, |dP/P| =",4ES12.4," -- falling back")') &
           max_refine, hydrogen_fraction, t6_temperature, pressure_e12, &
           abs((pressure_trial3-pressure_no_rad)/pressure_no_rad)
      rhoofp06 = opal_rho_not_found
      return
      end do refine
end function rhoofp06

!----------------------------------------------------------------------
! esac06
!----------------------------------------------------------------------
! OPAL 2006 EOS analogue of esac01 (see there and esac for
! the general description of the interpolation scheme); the
! structure is otherwise identical to esac01 (same iorder cap of
! 9, same iqu/ipu pre-checks, same l3==nr .or. k3==nt boundary reset).
!
!   v%eos_output(1) pressure, in megabars (1e12 dyne/cm**2)
!   v%eos_output(2) energy, in 1e12 erg/gm (zero at T6=0)
!   v%eos_output(3) entropy, in units of energy/T6
!   v%eos_output(4) dE/dRho at constant T6
!   v%eos_output(5) specific heat, dE/dT6 at constant volume
!   v%eos_output(6) dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   v%eos_output(7) dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   v%eos_output(8) gamma1 (Cox & Guili eq 9.88)
!   v%eos_output(9) gamma2/(gamma2-1) (Cox & Guili eq 9.88)
subroutine esac06(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = n_eos_mx, nr = n_eos06_nr, nt = n_eos06_nt

      real(dp) :: species_number_fraction(7)
      real(dp) :: molar_gas_constant_mbcc
      character(len=15) :: routine_id

! NOTE: this literal has no D0 suffix in the original (esac06.f) --
! preserved verbatim: it is parsed as a single-precision constant and
! then widened to double precision, which is NOT bit-identical to the
! correctly-rounded double value.
      data molar_gas_constant_mbcc/83.14511/
      data routine_id/"OPALEOS/ESAC06:"/
! --- locals ---
      integer :: x_loop_index_06
      real(dp) :: t6_value, density_value   ! working copies (slt/slr)
      integer :: species_idx, index_idx
      integer :: lo_idx, hi_idx, mid_idx, result_idx
      integer :: x_index_2, x_index_3, x_index_4, x_index_hi
      integer :: t6_order_hi, density_order_hi   ! ipu/iqu
      real(dp) :: table_sum_2x2, table_sum_3x3, table_sum_3x4, &
           table_sum_4x4
      integer :: eos_var_idx, x_print_idx
      integer :: density_scan_idx, t6_scan_idx
      integer :: recompute_flag, cache_slot
      real(dp) :: x_interp_weight
      real(dp) :: pressure_scale
      real(dp) :: total_moles, ground_state_energy, &
           metal_mole_fraction, mean_molecular_weight

      integer, intent(out) :: ierr
! in_table is false when (T6, rho) is outside the table or in an
! empty region of it
      logical, intent(out) :: in_table

      ierr = 0
      in_table = .true.

      if (deriv_order > 9) then
         write (run_log_unit,'(A, " iorder cannot exceed 9")') routine_id
      end if
      if ((rad_flag /= 0) .and. (rad_flag /= 1)) then
         write (run_log_unit,'(A, " Irad must be 0 or 1")') routine_id
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      t6_value = t6_temperature
      density_value = density
!
      if (opal_eos%table_loaded_flag(iv_opal06) /= opal_flag_set) then
         opal_eos%table_loaded_flag(iv_opal06) = opal_flag_set
         do index_idx = 1, 10
            do species_idx = 1, 10
               if (opal_eos_var_order(index_idx) == species_idx) &
                    v%eos_index_inverse(index_idx) = species_idx
            end do
         end do
         do species_idx = 1, mx
            v%x_grid_copy(species_idx) = opal06_x_grid(species_idx)
         end do
!
! ..... read the data files
         call readcoeos06(v, ierr)
         if (ierr /= 0) return
         opal_eos%table_metal_fraction(iv_opal06) = v%z_table(1)

         if (opal_eos%table_metal_fraction(iv_opal06)+hydrogen_fraction-1.0d-6 > 1) then
            write(run_log_unit,*) routine_id, "Mass fractions exceed unity (61)"
            write(run_log_unit,*) 'Z, XH', opal_eos%table_metal_fraction(iv_opal06), hydrogen_fraction
            
            
            ierr = 1
            return
         end if
      end if
!
!
! ..... Determine T6,rho grid points to use in the
!       interpolation.
      if ((t6_value > v%t6_grid(1)) .or. (t6_value < v%t6_grid(nt))) then
         write(run_log_unit,*) routine_id, " T6/LogR outside of table range (62)"
         write(run_log_unit,*) "slt, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
         v%t6_grid(nt)
         write(run_log_unit,*) "slr, rho(1), rho(nr):", density_value, &
         v%density_grid(1), v%density_grid(nr)
         in_table = .false.
         return
      end if
      if ((density_value < v%density_grid(1)) .or. (density_value > v%density_grid(nr))) then
         write(run_log_unit,*) routine_id, " T6/LogR outside of table range (62)"
         write(run_log_unit,*) "slt, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
         v%t6_grid(nt)
         write(run_log_unit,*) "slr, rho(1), rho(nr):", density_value, &
         v%density_grid(1), v%density_grid(nr)
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
         if (hydrogen_fraction <= opal06_x_grid(mid_idx)+1.0d-7) then
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
      if ((hydrogen_fraction <= opal06_x_grid(2)+1.0d-7) .or. &
           (hydrogen_fraction >= opal06_x_grid(mx-2)-1.0d-7)) x_index_hi = x_index_3
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
      if (v%t6_index_3 == 0) then
         write (run_log_unit,'(A, " ihi,ilo,imd",3I5)') routine_id, &
              hi_idx, lo_idx, mid_idx
      end if
!     endif

!     check to determine if interpolation indices fall within
!     table boundaries.  choose largest allowed size.
      table_sum_2x2 = 0.0d0
      table_sum_3x3 = 0.0d0
      table_sum_3x4 = 0.0d0
      table_sum_4x4 = 0.0d0
      do x_loop_index_06 = v%x_index_lo, x_index_hi
         do density_scan_idx = v%density_index_1, v%density_index_1+1
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+1
               table_sum_2x2 = table_sum_2x2 + &
                    v%eos_table(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         do density_scan_idx = v%density_index_1, v%density_index_1+2
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+2
               table_sum_3x3 = table_sum_3x3 + &
                    v%eos_table(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         if (t6_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+2
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  table_sum_3x4 = table_sum_3x4 + &
                       v%eos_table(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_3x4 = 2.0d+30
         end if
         if (density_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+3
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  table_sum_4x4 = table_sum_4x4 + &
                       v%eos_table(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
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
            write(run_log_unit,*) routine_id, "T6/log rho in empty region of table (65)"
            write(run_log_unit,'("xh,t6,r=", 3E12.4)') hydrogen_fraction, &
            t6_temperature, density
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
      do x_loop_index_06 = v%x_index_lo, x_index_hi

      recompute_flag = 0
!__________
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            v%x_interp_workspace(x_loop_index_06,t6_scan_idx,density_scan_idx) = &
                 v%eos_table(x_loop_index_06,eos_var_idx,t6_scan_idx,density_scan_idx)
            recompute_flag = 1
         end do
      end do
      end do
      if ((v%z_table(x_index_2) /= v%z_table(v%x_index_lo)) .or. &
           (v%z_table(x_index_3) /= v%z_table(v%x_index_lo))) then
         write(run_log_unit,'(A,"Z does not match Z in OEOS06 files you are" &
              &," using")') routine_id
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (opal_eos%table_metal_fraction(iv_opal06) /= v%z_table(v%x_index_lo)) exit
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
            if (v%x_interp_result(t6_scan_idx,density_scan_idx) > 1.0d+20) then
               write(run_log_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
                    routine_id, t6_scan_idx, density_scan_idx, v%density_index_3, &
                    v%t6_index_3, v%density_interp_order, v%t6_interp_order
               write(run_log_unit,'(3E12.4)') &
                    (v%x_interp_workspace(x_print_idx,t6_scan_idx,density_scan_idx), &
                    x_print_idx=v%x_index_lo,v%x_index_lo+2)
            end if
            recompute_flag = 1
         end do
      end do

      if (x_index_4 == x_index_hi) then  ! interpolate between quadratics
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
            if (v%x_interp_result(t6_scan_idx,density_scan_idx) > 1.0d+20) then
               write(run_log_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
                    routine_id, t6_scan_idx, density_scan_idx, v%density_index_3, &
                    v%t6_index_3, v%density_interp_order, v%t6_interp_order
               write(run_log_unit,'(3E12.4)') &
                    (v%x_interp_workspace(x_print_idx,t6_scan_idx,density_scan_idx), &
                    x_print_idx=x_index_2,x_index_2+2)
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
!
!
      call t6rinteos06(v, density_value, t6_value, ierr)
      if (ierr /= 0) return
      v%eos_output(eos_var_idx) = v%esact
      end do
      if (eos_var_idx > deriv_order) then
      pressure_scale = t6_temperature*density
      v%eos_output(v%eos_index_inverse(i_opal_p)) = v%eos_output(v%eos_index_inverse(i_opal_p))* &
           pressure_scale   ! interpolated in p/po
! Only slots that this call actually re-interpolated (index <= deriv_order)
! are rescaled; the original scaled E and cv unconditionally, so the
! deriv_order=1 trial calls from the rho(P,T) inversion compounded the
! stale E and cv slots by T6 and moles*R/mu on every call.
      if (v%eos_index_inverse(i_opal_e) <= deriv_order) &
           v%eos_output(v%eos_index_inverse(i_opal_e)) = v%eos_output(v%eos_index_inverse(i_opal_e))* &
           t6_temperature   ! interpolated in E/T6
      mean_molecular_weight = gmass01(hydrogen_fraction, opal_eos%table_metal_fraction(iv_opal06), &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_number_fraction)
      if (rad_flag == 1) then
         call radsub06(v, rad_flag, t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else
         if (v%eos_index_inverse(i_opal_cv) <= deriv_order) &
              v%eos_output(v%eos_index_inverse(i_opal_cv)) = v%eos_output(v%eos_index_inverse(i_opal_cv))* &
              total_moles*molar_gas_constant_mbcc/mean_molecular_weight
      end if
      return
      end if
      write(run_log_unit,*) routine_id, " Z does not match Z in EOSdata* files you ", &
           "are using (66)"
      write(run_log_unit,'("mf,zz(mf)=",I5,E12.4)') v%x_index_lo, &
           v%z_table(v%x_index_lo)
      write(run_log_unit,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') &
           v%t6_interp_order, v%density_interp_order, v%t6_index_3, v%density_index_3, &
           hydrogen_fraction, t6_temperature, density, opal_eos%table_metal_fraction(iv_opal06)
      ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1

end subroutine esac06

!----------------------------------------------------------------------
! t6rinteos06
!----------------------------------------------------------------------
! OPAL 2006 EOS analogue of t6rinteos01 (see there and
! t6rinterp for the general description). Declared RECURSIVE in
! the original; preserved verbatim even though nothing here actually
! recurses.
!
! Readability W3 (2026): the interpolation itself is interpolate_opal_t6_rho,
! shared with t6rinteos01; this wrapper keeps the 2006 out-of-range
! message.
recursive subroutine t6rinteos06(v, slr, slt, ierr)


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: slr, slt

      integer, intent(out) :: ierr

      ierr = 0

      call interpolate_opal_t6_rho(v, slr, slt)
      if (v%esact > 1.0d+15) then
         write(run_log_unit,'("T6RINTEOS06: Interpolation indices out", &
              &" of range;please report conditions.")')
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if

end subroutine t6rinteos06

!----------------------------------------------------------------------
! radsub06
!----------------------------------------------------------------------
! OPAL 2006 EOS analogue of radsub01 (see there for the general
! "norad, then diff in the gammas" strategy), but restructured so the
! rad_flag/irad test is a real IF around the radiation-correction
! block (radsub01 always compute the radiation-added
! values; the caller, esac06, currently only calls this routine
! when rad_flag=1, so rad_flag is effectively always nonzero here in
! practice, but the guard is preserved as written). total_pressure,
! total_energy, total_entropy, de_drho_at_t, chi_rho, chi_t6, and
! molar_specific_heat are computed once "without radiation" and then,
! if rad_flag/=0, OVERWRITTEN in place with the "with radiation"
! values -- unlike radsub01 there is no separate set of "_norad"
! variables for these seven; only the three gammas get separate
! "_norad" values (needed for the later difference). Also note only
! gamma1 and gamma2_over_gamma2_minus1 get the "add the difference"
! treatment here -- there is no gamma3_minus1 slot in the OPAL 2006
! table (9 variables, not 10; see esac06's header), and the
! corresponding line is commented out in the original.
subroutine radsub06(v, rad_flag, t6_temperature, density, total_moles, &
     mean_molecular_weight)


      type(opal_eos_vintage), intent(inout) :: v

      integer, intent(in) :: rad_flag
      real(dp), intent(in) :: t6_temperature, density
      real(dp), intent(in) :: total_moles
      real(dp), intent(in) :: mean_molecular_weight

      real(dp) :: rad_const_over_c, molar_gas_constant_mbcc
      data rad_const_over_c/1.8914785d-3/, molar_gas_constant_mbcc/83.14510d0/
! --- locals ---
      real(dp) :: rat, moles_per_ev
      real(dp) :: radiation_pressure, radiation_energy, radiation_entropy
      real(dp) :: total_pressure, total_energy, total_entropy
      real(dp) :: de_drho_at_t, chi_rho, chi_t6, molar_specific_heat
      real(dp) :: gamma3_minus1_norad, gamma1_norad, &
           gamma2_over_gamma2_minus1_norad
      real(dp) :: gamma3_minus1, gamma1, gamma2_over_gamma2_minus1

      rat = rad_const_over_c

      moles_per_ev = total_moles*molar_gas_constant_mbcc  ! Mb-cc/unit T6

!-----Calculate EOS without radiation correction

      total_pressure = v%eos_output(v%eos_index_inverse(i_opal_p))
      total_energy = v%eos_output(v%eos_index_inverse(i_opal_e))
      total_entropy = v%eos_output(v%eos_index_inverse(i_opal_s))
      de_drho_at_t = v%eos_output(v%eos_index_inverse(i_opal_dedrho))
      chi_rho = v%eos_output(v%eos_index_inverse(i_opal_chi_rho))*v%eos_output(v%eos_index_inverse(i_opal_p))/ &
           total_pressure
      chi_t6 = (v%eos_output(v%eos_index_inverse(i_opal_p))*v%eos_output(v%eos_index_inverse(i_opal_chi_t)))/ &
           total_pressure
      molar_specific_heat = (v%eos_output(v%eos_index_inverse(i_opal_cv))*moles_per_ev/ &
           mean_molecular_weight)
      gamma3_minus1_norad = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1_norad = chi_rho + chi_t6*gamma3_minus1_norad
      gamma2_over_gamma2_minus1_norad = gamma1_norad/gamma3_minus1_norad
!---- End  no radiation calculation

      if (rad_flag /= 0) then
!---- Calculate EOS with radiation calculation
         radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
         radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
         radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)
         total_pressure = v%eos_output(v%eos_index_inverse(i_opal_p)) + radiation_pressure
         total_energy = v%eos_output(v%eos_index_inverse(i_opal_e)) + radiation_energy
         total_entropy = v%eos_output(v%eos_index_inverse(i_opal_s)) + radiation_entropy
         de_drho_at_t = v%eos_output(v%eos_index_inverse(i_opal_dedrho)) - radiation_energy/density
         chi_rho = v%eos_output(v%eos_index_inverse(i_opal_chi_rho))*v%eos_output(v%eos_index_inverse(i_opal_p))/ &
              total_pressure
         chi_t6 = (v%eos_output(v%eos_index_inverse(i_opal_p))*v%eos_output(v%eos_index_inverse(i_opal_chi_t)) &
              + 4.0d0*radiation_pressure)/total_pressure
         molar_specific_heat = (v%eos_output(v%eos_index_inverse(i_opal_cv))*moles_per_ev/ &
              mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
         gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
              t6_temperature)                                        ! DIRECT
         gamma1 = chi_rho + chi_t6*gamma3_minus1 !eq 16.16 Landau_Lifshitz (Stat. Mech) ! DIRECT
         gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1             ! DIRECT

!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
!-----Add difference between EOS with and without radiation.  cvtt
!       calculation is not accurate enough to give accurate results using
!       eq. 16.16 Landau&Lifshitz (SEE line labeled DIRECT)
         v%eos_output(v%eos_index_inverse(i_opal_gamma1)) = v%eos_output(v%eos_index_inverse(i_opal_gamma1)) + &
              gamma1 - gamma1_norad
         v%eos_output(v%eos_index_inverse(i_opal_gamma2_ratio)) = v%eos_output(v%eos_index_inverse(i_opal_gamma2_ratio)) + &
              gamma2_over_gamma2_minus1 - gamma2_over_gamma2_minus1_norad
!     eos(iri(10))=eos(iri(10))+gam3pt-gam3pt_norad
      end if
!-----End EOS calculations with radiation
      v%eos_output(v%eos_index_inverse(i_opal_p)) = total_pressure
      v%eos_output(v%eos_index_inverse(i_opal_e)) = total_energy
      v%eos_output(v%eos_index_inverse(i_opal_s)) = total_entropy
      v%eos_output(v%eos_index_inverse(i_opal_dedrho)) = de_drho_at_t
      v%eos_output(v%eos_index_inverse(i_opal_cv)) = molar_specific_heat
      v%eos_output(v%eos_index_inverse(i_opal_chi_rho)) = chi_rho
      v%eos_output(v%eos_index_inverse(i_opal_chi_t)) = chi_t6
end subroutine radsub06

!----------------------------------------------------------------------
! readcoeos06
!----------------------------------------------------------------------
! Reads the OPAL 2006 EOS data tables into opal_eos_tables's opal_eos state
! and derives the auxiliary grid-spacing arrays used by esac06/
! t6rinteos06 for interpolation. Called once from esac06
! (guarded by opal_eos%table_loaded_flag).
!
! NOTE: preserved verbatim from the original -- var_idx, the "j" loop
! variable of the one-time init block below (which zeroes
! v%eos_table and resets v%eos_output to 1.0), is
! reused, stale, as the column index into log10_ne_grid_06 in the
! per-row READ statement further down. Because the init loops leave
! var_idx at mv+1, every row's read of
! log10_ne_grid_06(density_row, var_idx) overwrites the SAME column
! (var_idx == mv+1) instead of tracking t6_row as the analogous
! v%t6_list/amu_grid_06 reads do. This looks like an original
! bug (using "j" instead of "i" in that READ's variable list) but is
! transliterated exactly rather than "fixed"; log10_ne_grid_06 is never
! read anywhere.
subroutine readcoeos06(v, ierr)

      type(opal_eos_vintage), intent(inout) :: v

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv, nr = n_eos06_nr, nt = n_eos06_nt

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index_06
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      integer :: species_read_idx, table_var_idx
      real(dp) :: unused_field
      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '

      if (opal_eos%readco_init_flag(iv_opal06) /= opal_flag_set) then
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
         opal_eos%readco_init_flag(iv_opal06) = opal_flag_set
      end if

      close (2)
! .....read  tables (the file is opened on opal_eos%table_unit by read_controls)
      do x_loop_index_06 = 1, mx
         read (opal_eos%table_unit,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              v%hydrogen_fraction_header(x_loop_index_06), &
              v%z_table(x_loop_index_06), &
              v%moles_per_gram_table(x_loop_index_06), &
              v%mean_molecular_weight_header(x_loop_index_06)
         read (opal_eos%table_unit,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (v%species_fraction_header(x_loop_index_06,species_read_idx), &
              species_read_idx=1,6)
         read (opal_eos%table_unit,'(A)') blank_line
         do density_row = 1, nr
            read (opal_eos%table_unit,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 v%temperature_count_used(x_loop_index_06,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index_06,density_row)
            if (record_number /= density_row) then
               write (run_log_unit,'(" OEOS06 Data file incorrect: numtot,jcs= " &
                    &,2I5)') record_number, density_row
               ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
            if (v%temperature_count_used(x_loop_index_06,density_row) < &
                 opal06_t6_index_lo(density_row)) then
               write (run_log_unit,'("Problem with OEOS96 data files: X=",F6.4, &
                    &" density=",E14.4)') v%hydrogen_fraction_header(x_loop_index_06), &
                    v%density_grid_table(x_loop_index_06,density_row)
               ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
               ! facades stop when their caller passes no ierr.
               ierr = 1
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index_06,density_row)
               if (t6_row > opal06_t6_index_lo(density_row)) then
                  read (opal_eos%table_unit,'(A)') blank_line
                  cycle
               end if
               read (opal_eos%table_unit,'(F11.6,1X,F6.4,E11.4,2E13.6,2E11.3,5F10.6)') &
                    v%t6_list(density_row,t6_row), v%amu_grid(density_row,t6_row), &
                    v%log10_ne_grid(density_row,var_idx), &
                    (v%eos_table(x_loop_index_06,opal_eos_var_order(table_var_idx), &
                    t6_row,density_row), table_var_idx=1,9)
            end do
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
         end do
         read(opal_eos%table_unit,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx) == 0.0d0) then
            write(run_log_unit,'("READCOEOS06: Error:",I4, &
                 &"-th T6 value is zero")') t6_scan_idx
            ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
            ! facades stop when their caller passes no ierr.
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
      close (opal_eos%table_unit)

end subroutine readcoeos06

end module opal2006_eos
