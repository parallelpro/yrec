!----------------------------------------------------------------------
! opal_eos
!----------------------------------------------------------------------
! The OPAL equation of state: the entry points for all three vintages,
! and the 2001 and 2006 tables, which share one state layout
! (opal_eos_vintage; the 1995 tables, with their own layout, are in
! opal95_eos).
!     compute_opal_eos     density and the thermodynamic derivatives of
!                          one vintage at (log T, log P, X, Z)
!     find_opal_edge_ramp  whether the vintage's table holds the point,
!                          and the OPAL weight near its edge
! Under these: rhoofp01/06 (OPAL's density search; the two differ),
! interpolate_opal_eos (OPAL's esac01/esac06, which remain its entry
! points; esac06 is also public for eos_get_gamma1 and test_eos),
! t6rinteos, add_radiation (the radiation correction) and readcoeos
! (the table reader, run once per vintage). The 2001 and 2006 code is
! written once: what the vintages differ in (the table's size and grids,
! the slots, the radiation constants, the reader's row format) is their
! row of opal_layout (opal_eos_tables), and the report_* procedures keep
! each vintage's messages in its own words. The 2001 code follows
! Rogers & Nayfonov (2002); the 2006 tables are the 2005 update.
module opal_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use opal95_eos, only: esac, rhoofp, find_opal95_ramp
      use opal_eos_tables, only: find_opal_ramp, gmass01, i_opal01_p, i_opal_p, &
           interpolate_opal_t6_rho, iv_opal95, iv_opal01, iv_opal06, n_eos01_nt, &
           n_eos06_nr, n_eos06_nt, n_eos_mv, n_eos_mx, opal95, opal01, opal06, &
           opal01_density_index_edge, opal01_x_grid, opal06_density_index_edge, &
           opal06_t6_index_lo, opal06_x_grid, opal_layout, opal_shared, &
           opal_eos_vintage, opal_flag_set, opal_rho_not_found, opal_vintage_layout, quad
      use physical_constants, only: gas_constant
      use eos_messages, only: eos_log, message_length, unwritten_line
      use eos_mixture, only: compute_mean_weights
      use math_lib
      implicit none
      private
      public :: compute_opal_eos, find_opal_edge_ramp, esac06
      public :: iv_opal95, iv_opal01, iv_opal06, opal_layout

contains

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
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
      real(dp) :: output(n_eos_mv)
      real(dp) :: t_million_k, p_e12, hydrogen_fraction_work, density_cgs
      real(dp) :: specific_gas_constant_check
      real(dp) :: chi_rho, chi_t, cv, gamma1, gamma2_ratio
      integer :: rad_flag, deriv_order

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

      select case (vintage)
      case (iv_opal95)
         density_cgs = rhoofp(opal95, hydrogen_fraction_work, t_million_k, p_e12, rad_flag, ierr)
      case (iv_opal01)
         density_cgs = rhoofp01(opal01, hydrogen_fraction_work, t_million_k, p_e12, &
              rad_flag, ierr)
      case default
         density_cgs = rhoofp06(opal06, hydrogen_fraction_work, t_million_k, p_e12, &
              rad_flag, ierr)
      end select
      if (ierr /= 0) return
      if (density_cgs <= -998.0d0) then
         in_table = .false.
         return
      end if
      r%density = density_cgs
      r%log10_density = log10(r%density)

      select case (vintage)
      case (iv_opal95)
         call esac(opal95, hydrogen_fraction_work, t_million_k, density_cgs, &
              deriv_order, rad_flag, ierr, in_table)
      case (iv_opal01)
         call esac01(opal01, hydrogen_fraction_work, t_million_k, density_cgs, &
              deriv_order, rad_flag, ierr, in_table)
      case default
         call esac06(opal06, hydrogen_fraction_work, t_million_k, density_cgs, &
              deriv_order, rad_flag, ierr, in_table)
      end select
      if (.not. in_table) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'WARNNING... OPAL TBL FAIL'
         call eos_log(eos_lines)
         return
      end if
      if (ierr /= 0) return

! dln P/dln rho, dln P/dln T6, cv, gamma1 and gamma2/(gamma2-1) from the
! vintage's results, in its slots
      select case (vintage)
      case (iv_opal95)
         output = opal95%eos_output
      case (iv_opal01)
         output = opal01%eos_output
      case default
         output = opal06%eos_output
      end select
      if (opal_layout(vintage)%checks_pressure .and. &
           abs((p_e12-output(opal_layout(vintage)%i_p))/p_e12) > 0.5d-6) then
         eos_lines = unwritten_line
         write(eos_lines,*) p_e12, output(opal_layout(vintage)%i_p)
         call eos_log(eos_lines)
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

      r%beta14 = (2.521971383d-3*t_million_k*t_million_k)* &
           (t_million_k*t_million_k/p_e12)
      r%beta = 1.0d0 - r%beta14
      r%beta_inverse = 1.0d0/r%beta
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
! OPAL result near the table's edge (find_opal95_ramp,
! find_opal2001_ramp, find_opal2006_ramp; the last returns no ierr).
subroutine find_opal_edge_ramp(vintage, temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
      integer, intent(out) :: ierr

      select case (vintage)
      case (iv_opal95)
         call find_opal95_ramp(temperature, log10_density, ramp_factor, &
              in_opal_table, needs_ramp, ierr)
      case (iv_opal01)
         call find_opal2001_ramp(temperature, log10_density, ramp_factor, &
              in_opal_table, needs_ramp, ierr)
      case default
         ierr = 0
         call find_opal2006_ramp(temperature, log10_density, ramp_factor, &
              in_opal_table, needs_ramp)
      end select
end subroutine find_opal_edge_ramp

!----------------------------------------------------------------------
! find_opal2001_ramp (formerly eqbound01)
!----------------------------------------------------------------------
! MHP 7/2003: checks whether a (T, log10 rho) point falls within the
! OPAL 2001 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (blend_eos calls it, through find_opal_edge_ramp, right after
! compute_opal_eos).
!
! The 2001 density check (log10 rho outside -14..7 is outside the
! table), then find_opal_ramp.
subroutine find_opal2001_ramp(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)

      integer, parameter :: nt = n_eos01_nt

      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
! --- locals ---
      real(dp) :: t6

      integer, intent(out) :: ierr

      ierr = 0

      t6 = temperature*1.0d-6
!     exit if outside table in rho
      if (log10_density < -14d0 .or. log10_density > 7.0d0) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if
      call find_opal_ramp(nt, opal01%t6_grid, opal01%t_row_index, opal01%density_edge_at_t, &
           opal01_density_index_edge, opal01%density_grid, t6, log10_density, ramp_factor, &
           in_opal_table, needs_ramp, ierr)

end subroutine find_opal2001_ramp

!----------------------------------------------------------------------
! rhoofp01
!----------------------------------------------------------------------
! OPAL 2001 EOS analogue of rhoofp (see there for the general
! description: the bracket, at most 30 regula falsi steps, the run-log
! message on failure), with the same convergence tolerance (0.5d-7) as
! rhoofp and rhoofp06.
real(dp) function rhoofp01(v, hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)


      character(len=message_length) :: eos_lines(4)
      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = n_eos_mx, nt = n_eos01_nt

! NOTE: no D-suffix in the original (rhoofp01.f) -- preserved
! verbatim: parsed as a single-precision constant and then widened to
! double precision, which is NOT bit-identical to the correctly-
! rounded double value.
      real(dp), parameter :: rad_const_over_c = 1.8914785e-3

! --- locals ---
      integer :: t6_scan_idx
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

! The highest valid density column at each T6 row. These are the
! values of opal01_density_index_edge except at rows 78-80 (113 here,
! 115 there), as in OPAL's rhoofp01.
      integer, parameter :: density_index_edge(nt) = [ &
           (169,t6_scan_idx=1,16), 168, 167, 166, 165, &
           (164,t6_scan_idx=1,2), 163, (162,t6_scan_idx=1,2), 161, 160, &
           (159,t6_scan_idx=1,2), (143,t6_scan_idx=1,4), &
           (137,t6_scan_idx=1,5), (134,t6_scan_idx=1,6), &
           (125,t6_scan_idx=1,2), (123,t6_scan_idx=1,5), &
           (122,t6_scan_idx=1,2), (121,t6_scan_idx=1,6), &
           (119,t6_scan_idx=1,4), (116,t6_scan_idx=1,8), &
           (115,t6_scan_idx=1,6), (113,t6_scan_idx=1,8), &
           (111,t6_scan_idx=1,7), (110,t6_scan_idx=1,6), &
           (109,t6_scan_idx=1,34), 107, 104, (100,t6_scan_idx=1,40), &
           (99,t6_scan_idx=1,10), 98, 97, 96, 95, 94, 93, 92]
      integer, intent(out) :: ierr
      logical :: in_table

      ierr = 0

      rat = rad_const_over_c
! The radiation pressure is deliberately not subtracted from the target
! here (the original's IF(IRAD.EQ.1) PR=4/3*RAT*T6**4 is commented out);
! rad_flag instead adds it to pressure_max/pressure_min below.
      radiation_pressure = 0.0d0
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_shared%table_loaded_flag(iv_opal01) /= opal_flag_set) then
         hydrogen_fraction_dbg = 0.5d0
         t6_dbg = 1.0d0
         density_dbg = 0.001d0
         ideriv_dbg = 1
         call esac01(v, hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              rad_flag, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp01 = opal_rho_not_found
            return
         end if
      end if

      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction <= opal01_x_grid(mid_idx)+1.0d-7) then
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
           density_index_edge(t6_bisect_idx))*t6_temperature* &
           v%density_grid(density_index_edge(t6_bisect_idx)) + &
           rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      pressure_min = v%eos_table(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           v%density_grid(1) + rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      if ((pressure_no_rad > 1.25d0*pressure_max) .or. &
           (pressure_no_rad < pressure_min)) then
!        requested pressure-temperature not in table
         rhoofp01 = opal_rho_not_found
         return
      end if

      density_trial1 = v%density_grid(density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac01(v, hydrogen_fraction, t6_temperature, density_trial1, 1, &
           rad_flag, ierr, in_table)
      if (ierr /= 0 .or. .not. in_table) then
         rhoofp01 = opal_rho_not_found
         return
      end if
      pressure_trial1 = v%eos_output(i_opal01_p)
      if (pressure_trial1 > pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
! shrink the lower density until it brackets the root (regula falsi
! converges only from a bracket), down to the table's floor
         shrink_count = 0
         shrink: do
            shrink_count = shrink_count + 1
            density_trial1 = 0.2d0*density_trial1
            if (density_trial1 < 1.0d-14) density_trial1 = 1.0d-14
            call esac01(v, hydrogen_fraction, t6_temperature, density_trial1, 1, &
                 rad_flag, ierr, in_table)
            if (ierr /= 0 .or. .not. in_table) then
               rhoofp01 = opal_rho_not_found
               return
            end if
            pressure_trial1 = v%eos_output(i_opal01_p)
            if (pressure_trial1 <= pressure_no_rad) exit shrink
            if (density_trial1 <= 1.0d-14 .or. shrink_count >= max_shrink) exit shrink
         end do shrink
      else
         density_trial2 = 5.0d0*density_trial1
         if (density_trial2 > v%density_grid(density_index_edge(t6_bisect_idx))) &
              density_trial2 = v%density_grid(density_index_edge(t6_bisect_idx))
         call esac01(v, hydrogen_fraction, t6_temperature, density_trial2, 1, &
              rad_flag, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp01 = opal_rho_not_found
            return
         end if
         pressure_trial2 = v%eos_output(i_opal01_p)
      end if

      refine_count = 0
      refine: do
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
      call esac01(v, hydrogen_fraction, t6_temperature, density_trial3, 1, &
           rad_flag, ierr, in_table)
      if (ierr /= 0 .or. .not. in_table) then
         rhoofp01 = opal_rho_not_found
         return
      end if
      pressure_trial3 = v%eos_output(i_opal01_p)
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad) < 0.5d-7) then
         rhoofp01 = density_trial3
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
      write(eos_lines,'("RHOOFP01: no convergence after ",I0," tries;", " X, T6, P12, |dP/P| =",4ES12.4," -- falling back")') max_refine, &
           hydrogen_fraction, t6_temperature, pressure_e12, abs((pressure_trial3-pressure_no_rad)/pressure_no_rad)
      call eos_log(eos_lines)
      rhoofp01 = opal_rho_not_found
      return
      end do refine
end function rhoofp01

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
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = n_eos_mx
! the vintage's table: its size, its X grid and its name in messages
      integer :: nr, nt
      real(dp) :: x_grid(n_eos_mx)
      character(len=15) :: routine_id

      real(dp) :: species_number_fraction(7)

! NOTE: this literal has no D0 suffix in the original (esac06.f) --
! preserved verbatim: it is parsed as a single-precision constant and
! then widened to double precision, which is NOT bit-identical to the
! correctly-rounded double value.
      real(dp), parameter :: molar_gas_constant_mbcc = 83.14511
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

      integer, intent(out) :: ierr
! in_table is false when (T6, rho) is outside the table or in an
! empty region of it
      logical, intent(out) :: in_table

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
         call eos_log(eos_lines)
         ierr = 1
         return
      end if
      t6_value = t6_temperature
      density_value = density
!
      if (opal_shared%table_loaded_flag(vintage) /= opal_flag_set) then
         opal_shared%table_loaded_flag(vintage) = opal_flag_set
         do species_idx = 1, mx
            v%x_grid_copy(species_idx) = x_grid(species_idx)
         end do
!
! ..... read the data files
         call readcoeos(v, vintage, ierr)
         if (ierr /= 0) return
         opal_shared%table_metal_fraction(vintage) = v%z_table(1)

         if (opal_shared%table_metal_fraction(vintage)+hydrogen_fraction-1.0d-6 > 1) then
            eos_lines = unwritten_line
            write(eos_lines,*) routine_id, "Mass fractions exceed unity (61)"
            call eos_log(eos_lines)
            eos_lines = unwritten_line
            write(eos_lines,*) 'Z, XH', opal_shared%table_metal_fraction(vintage), hydrogen_fraction
            call eos_log(eos_lines)
            
            
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
      do x_loop_index = v%x_index_lo, x_index_hi

      recompute_flag = 0
!__________
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            v%x_interp_workspace(x_loop_index,t6_scan_idx,density_scan_idx) = &
                 v%eos_table(x_loop_index,eos_var_idx,t6_scan_idx,density_scan_idx)
            recompute_flag = 1
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
            if (v%x_interp_result(t6_scan_idx,density_scan_idx) > 1.0d+20) &
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
      mean_molecular_weight = gmass01(hydrogen_fraction, opal_shared%table_metal_fraction(vintage), &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_number_fraction)
      if (rad_flag == 1) then
         call add_radiation(v, lay, t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else
         if (lay%i_cv <= deriv_order) &
              v%eos_output(lay%i_cv) = v%eos_output(lay%i_cv)* &
              total_moles*molar_gas_constant_mbcc/mean_molecular_weight
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
      integer, intent(out) :: ierr

      nr = opal_layout(vintage)%nr
      nt = opal_layout(vintage)%nt
      ierr = 0

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

      close (2)
! .....read  tables (on the unit open_opal_eos_table opened)
      do x_loop_index = 1, mx
         read (opal_shared%table_unit,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              v%hydrogen_fraction_header(x_loop_index), &
              v%z_table(x_loop_index), unused_header(1), unused_header(2)
         read (opal_shared%table_unit,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') unused_species
         read (opal_shared%table_unit,'(A)') blank_line
         do density_row = 1, nr
            read (opal_shared%table_unit,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 v%temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index,density_row)
            if (record_number /= density_row) then
               call report_record_number(vintage, record_number, density_row)
               ierr = 1
               return
            end if
            read(opal_shared%table_unit,'(A)') blank_line
            read(opal_shared%table_unit,'(A)') blank_line
            if (v%temperature_count_used(x_loop_index,density_row) < &
                 opal_layout(vintage)%t6_index_lo(density_row)) then
               call report_short_column(vintage, v%hydrogen_fraction_header(x_loop_index), &
                    v%density_grid_table(x_loop_index,density_row))
               ierr = 1
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index,density_row)
               if (t6_row > opal_layout(vintage)%t6_index_lo(density_row)) then
                  read (opal_shared%table_unit,'(A)') blank_line
                  cycle
               end if
               read (opal_shared%table_unit, opal_layout(vintage)%row_format) &
                    v%t6_list(density_row,t6_row), &
                    (unused_columns(column_idx), column_idx=1,opal_layout(vintage)%n_skipped_columns), &
                    (v%eos_table(x_loop_index,table_var_idx,t6_row,density_row), &
                    table_var_idx=1,9)
            end do
            read(opal_shared%table_unit,'(A)') blank_line
            read(opal_shared%table_unit,'(A)') blank_line
            read(opal_shared%table_unit,'(A)') blank_line
         end do
         read(opal_shared%table_unit,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx) == 0.0d0) then
            call report_zero_t6(vintage, t6_scan_idx)
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
! find_opal2006_ramp (formerly eqbound06)
!----------------------------------------------------------------------
! LLP 10/16/06: checks whether a (T, log10 rho) point falls within the
! OPAL 2006 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (blend_eos calls it, through find_opal_edge_ramp, right after
! compute_opal_eos). Unlike find_opal95_ramp/find_opal2001_ramp (1995/2001),
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


      character(len=message_length) :: eos_lines(4)
      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = n_eos_mx, nt = n_eos06_nt

! NOTE: no D-suffix in the original (rhoofp06.f) -- preserved
! verbatim: parsed as a single-precision constant and then widened to
! double precision, which is NOT bit-identical to the correctly-
! rounded double value.
      real(dp), parameter :: rad_const_over_c = 1.8914785e-3
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

      if (opal_shared%table_loaded_flag(iv_opal06) /= opal_flag_set) then
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
      eos_lines = unwritten_line
      write(eos_lines,'("RHOOFP06: no convergence after ",I0," tries;", " X, T6, P12, |dP/P| =",4ES12.4," -- falling back")') max_refine, &
           hydrogen_fraction, t6_temperature, pressure_e12, abs((pressure_trial3-pressure_no_rad)/pressure_no_rad)
      call eos_log(eos_lines)
      rhoofp06 = opal_rho_not_found
      return
      end do refine
end function rhoofp06


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
      radiation_pressure = 4.0d0/3.0d0*lay%rad_const_over_c*t6_temperature**4   ! Mb
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

!----------------------------------------------------------------------
! The 2001 and 2006 interpolators' messages, each vintage's in its
! original words (run log).
!----------------------------------------------------------------------
! (T6, rho) outside the table
subroutine report_outside_table(vintage, t6, density, t6_first, t6_last, &
     density_first, density_last)
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: t6, density, t6_first, t6_last, density_first, density_last
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
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage, hi_idx, lo_idx, mid_idx
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
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: hydrogen_fraction, t6, density
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
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'(A,"Z does not match Z in OEOS01 files you are" ," using")') opal_layout(vintage)%routine_id
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'(A,"Z does not match Z in OEOS06 files you are" ," using")') opal_layout(vintage)%routine_id
         call eos_log(eos_lines)
      end if
end subroutine report_tables_z

! an X-interpolated value above 1e20 (an empty entry in the stencil):
! the stencil and the three table values from x_first
subroutine report_x_interpolation(vintage, v, t6_idx, density_idx, x_first)
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage, t6_idx, density_idx, x_first
      type(opal_eos_vintage), intent(in) :: v
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
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage
      type(opal_eos_vintage), intent(in) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6, density
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,*) opal_layout(vintage)%routine_id, " Z does not match Z in OPAL 2001 EOS files", &
              " you are using (66)"
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,*) opal_layout(vintage)%routine_id, " Z does not match Z in EOSdata* files you ", &
              "are using (66)"
         call eos_log(eos_lines)
      end if
      eos_lines = unwritten_line
      write(eos_lines,'("mf,zz(mf)=",I5,E12.4)') v%x_index_lo, v%z_table(v%x_index_lo)
      call eos_log(eos_lines)
      eos_lines = unwritten_line
      write(eos_lines,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') v%t6_interp_order, &
           v%density_interp_order, v%t6_index_3, v%density_index_3, hydrogen_fraction, &
           t6, density, opal_shared%table_metal_fraction(vintage)
      call eos_log(eos_lines)
end subroutine report_run_z

! the reader: a density block out of order
subroutine report_record_number(vintage, record_number, density_row)
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage, record_number, density_row
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'(" Data file incorrect: numtot,jcs= ",2I5)') record_number, density_row
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'(" OEOS06 Data file incorrect: numtot,jcs= " ,2I5)') record_number, density_row
         call eos_log(eos_lines)
      end if
end subroutine report_record_number

! the reader: a density column with fewer T6 rows than the table needs
subroutine report_short_column(vintage, hydrogen_fraction, density)
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: hydrogen_fraction, density
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'("problem with data files: X=",F6.4, " density=",E14.4)') hydrogen_fraction, density
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'("Problem with OEOS96 data files: X=",F6.4, " density=",E14.4)') hydrogen_fraction, &
              density
         call eos_log(eos_lines)
      end if
end subroutine report_short_column

! the reader: a zero T6
subroutine report_zero_t6(vintage, t6_idx)
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage, t6_idx
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'("READCOEOS01: Error:",I4, "-th T6 value is zero")') t6_idx
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'("READCOEOS06: Error:",I4, "-th T6 value is zero")') t6_idx
         call eos_log(eos_lines)
      end if
end subroutine report_zero_t6

! the T6-rho interpolation's indices out of range
subroutine report_interpolation_indices(vintage)
      character(len=message_length) :: eos_lines(4)
      integer, intent(in) :: vintage
      if (vintage == iv_opal01) then
         eos_lines = unwritten_line
         write(eos_lines,'("Interpolation indices out of range", ";please report conditions.")')
         call eos_log(eos_lines)
      else
         eos_lines = unwritten_line
         write(eos_lines,'("T6RINTEOS06: Interpolation indices out", " of range;please report conditions.")')
         call eos_log(eos_lines)
      end if
end subroutine report_interpolation_indices

end module opal_eos
