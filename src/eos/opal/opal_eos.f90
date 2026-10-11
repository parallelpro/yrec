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
! t6rinteos, radsub01/06 (the radiation corrections; the two differ)
! and readcoeos (the table reader, run once per vintage). Where 2001 and
! 2006 share code it is written once, branching on the vintage only
! where they differ. The 2001 code follows
! Rogers & Nayfonov (2002); the 2006 tables are the 2005 update.
module opal_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use opal95_eos, only: esac, rhoofp, find_opal95_ramp
      use opal_eos_tables, only: find_opal_ramp, gmass01, &
           i_opal01_chi_rho, i_opal01_chi_t, i_opal01_cv, i_opal01_e, i_opal01_gamma1, &
           i_opal01_gamma2_ratio, i_opal01_gamma3m1, i_opal01_p, i_opal01_s, &
           i_opal_chi_rho, i_opal_chi_t, i_opal_cv, i_opal_dedrho, i_opal_e, i_opal_gamma1, &
           i_opal_gamma2_ratio, i_opal_p, i_opal_s, interpolate_opal_t6_rho, &
           iv_opal95, iv_opal01, iv_opal06, n_opal_vintages, n_eos01_nr, n_eos01_nt, &
           n_eos06_nr, n_eos06_nt, n_eos_mv, n_eos_mx, n_eos_nr_max, opal95, opal01, opal06, &
           opal01_density_index_edge, opal01_t6_index_lo, opal01_x_grid, &
           opal06_density_index_edge, opal06_t6_index_lo, opal06_x_grid, &
           opal_shared, opal_eos_var_order, opal_eos_vintage, opal_flag_set, &
           opal_rho_not_found, quad
      use physical_constants, only: gas_constant
      use io_units, only: run_log_unit
      use eos_mixture, only: compute_mean_weights
      use math_lib
      implicit none
      private
      public :: compute_opal_eos, find_opal_edge_ramp, esac06
      public :: iv_opal95, iv_opal01, iv_opal06

contains

!----------------------------------------------------------------------
! compute_opal_eos (formerly oeqos, oeqos01, oeqos06)
!----------------------------------------------------------------------
! The OPAL EOS of one vintage (iv_opal95, iv_opal01 or iv_opal06) at
! (log T, log P, X, Z): the density that gives P (rhoofp*), the
! thermodynamic quantities there (esac*), beta and the mean weights.
! in_table is false outside the vintage's T range, when no density
! matches P, or when esac* finds the point outside its table. The
! vintages differ in their T6 range, in how many variables esac fills
! (10 for 1995, 9 for the others; the 2001 table has no dE/drho, so its
! slots are one lower), and in 1995's check that the interpolated
! pressure reproduces P.
subroutine compute_opal_eos(vintage, log10_temperature, log10_pressure, &
     hydrogen_fraction, metal_fraction, r, ierr, in_table)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
! each vintage's T6 range and the number of variables asked of esac
      real(dp), parameter :: t6_min(n_opal_vintages) = [0.0050d0, 0.0020d0, 0.001870d0]
      real(dp), parameter :: t6_max(n_opal_vintages) = [100.0d0, 100.0d0, 200.0d0]
      integer, parameter :: deriv_orders(n_opal_vintages) = [10, 9, 9]
      real(dp) :: t_million_k, p_e12, hydrogen_fraction_work, density_cgs
      real(dp) :: specific_gas_constant_check
      real(dp) :: chi_rho, chi_t, cv, gamma1, gamma2_ratio
      integer :: rad_flag, deriv_order

      ierr = 0
      in_table = .true.
      rad_flag = 1     ! adds the radiation correction
      deriv_order = deriv_orders(vintage)

      r%temperature = exp10(log10_temperature)
      r%pressure = exp10(log10_pressure)
      t_million_k = r%temperature/1.0d6
      p_e12 = r%pressure/1.0d12
      if (t_million_k < t6_min(vintage) .or. t_million_k > t6_max(vintage)) then
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
         write(run_log_unit, *) 'WARNNING... OPAL TBL FAIL'
         return
      end if
      if (ierr /= 0) return

! dln P/dln rho, dln P/dln T6, cv, gamma1 and gamma2/(gamma2-1) from the
! vintage's results
      select case (vintage)
      case (iv_opal95)
         if (abs((p_e12-opal95%eos_output(i_opal_p))/p_e12) > 0.5d-6) then
            write(run_log_unit,*) p_e12, opal95%eos_output(i_opal_p)
            ierr = 1
            return
         end if
         chi_rho = opal95%eos_output(i_opal_chi_rho)
         chi_t = opal95%eos_output(i_opal_chi_t)
         cv = opal95%eos_output(i_opal_cv)
         gamma1 = opal95%eos_output(i_opal_gamma1)
         gamma2_ratio = opal95%eos_output(i_opal_gamma2_ratio)
      case (iv_opal01)
         chi_rho = opal01%eos_output(i_opal01_chi_rho)
         chi_t = opal01%eos_output(i_opal01_chi_t)
         cv = opal01%eos_output(i_opal01_cv)
         gamma1 = opal01%eos_output(i_opal01_gamma1)
         gamma2_ratio = opal01%eos_output(i_opal01_gamma2_ratio)
      case default
         chi_rho = opal06%eos_output(i_opal_chi_rho)
         chi_t = opal06%eos_output(i_opal_chi_t)
         cv = opal06%eos_output(i_opal_cv)
         gamma1 = opal06%eos_output(i_opal_gamma1)
         gamma2_ratio = opal06%eos_output(i_opal_gamma2_ratio)
      end select
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
! result (see blend_eos.f90's use_opal2001_eos branch, which calls this
! immediately after compute_opal2001_eos).
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
! description), with the same convergence tolerance (0.5d-7) as rhoofp
! and rhoofp06.
real(dp) function rhoofp01(v, hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)


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
      integer :: refine_count

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
         density_trial1 = 0.2d0*density_trial1
         if (density_trial1 < 1.0d-14) density_trial1 = 1.0d-14
         call esac01(v, hydrogen_fraction, t6_temperature, density_trial1, 1, &
              rad_flag, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp01 = opal_rho_not_found
            return
         end if
         pressure_trial1 = v%eos_output(i_opal01_p)
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
         if (refine_count < 11) cycle refine
!        no convergence after 10 tries
         rhoofp01 = opal_rho_not_found
         return
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
         if (refine_count < 11) cycle refine
!        no convergence after 10 tries
         rhoofp01 = opal_rho_not_found
         return
      end if
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
! stencil is clamped to the table. The vintages differ in their table
! size, their messages and the 2001 table's missing dE/dRho (its slots
! after the entropy are one lower), and the code branches on vintage
! only there.
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

      if (vintage == iv_opal01) then
         nr = n_eos01_nr
         nt = n_eos01_nt
         x_grid = opal01_x_grid
         routine_id = "OPALEOS/ESAC01:"
      else
         nr = n_eos06_nr
         nt = n_eos06_nt
         x_grid = opal06_x_grid
         routine_id = "OPALEOS/ESAC06:"
      end if
      ierr = 0
      in_table = .true.

      if (deriv_order > 9) then
         write (run_log_unit,'(A, " iorder cannot exceed 9")') routine_id
      end if
      if ((rad_flag /= 0) .and. (rad_flag /= 1)) then
         write (run_log_unit,'(A, " Irad must be 0 or 1")') routine_id
         ierr = 1
         return
      end if
      t6_value = t6_temperature
      density_value = density
!
      if (opal_shared%table_loaded_flag(vintage) /= opal_flag_set) then
         opal_shared%table_loaded_flag(vintage) = opal_flag_set
         do index_idx = 1, 10
            do species_idx = 1, 10
               if (opal_eos_var_order(index_idx) == species_idx) &
                    v%eos_index_inverse(index_idx) = species_idx
            end do
         end do
         do species_idx = 1, mx
            v%x_grid_copy(species_idx) = x_grid(species_idx)
         end do
!
! ..... read the data files
         call readcoeos(v, vintage, ierr)
         if (ierr /= 0) return
         opal_shared%table_metal_fraction(vintage) = v%z_table(1)

         if (opal_shared%table_metal_fraction(vintage)+hydrogen_fraction-1.0d-6 > 1) then
            write(run_log_unit,*) routine_id, "Mass fractions exceed unity (61)"
            write(run_log_unit,*) 'Z, XH', opal_shared%table_metal_fraction(vintage), hydrogen_fraction
            
            
            ierr = 1
            return
         end if
      end if
!
!
! ..... Determine T6,rho grid points to use in the
!       interpolation.
      if ((t6_value > v%t6_grid(1)) .or. (t6_value < v%t6_grid(nt))) then
         if (vintage == iv_opal01) then
            write(run_log_unit,*) routine_id, " T6 or rho outside of table range (62)"
            write(run_log_unit,*) "t6, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
            v%t6_grid(nt)
            write(run_log_unit,*) "slr,r,rho(1),rho(nr):", density_value, &
            density, v%density_grid(1), v%density_grid(nr)
         else
            write(run_log_unit,*) routine_id, " T6/LogR outside of table range (62)"
            write(run_log_unit,*) "slt, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
            v%t6_grid(nt)
            write(run_log_unit,*) "slr, rho(1), rho(nr):", density_value, &
            v%density_grid(1), v%density_grid(nr)
         end if
         in_table = .false.
         return
      end if
      if ((density_value < v%density_grid(1)) .or. (density_value > v%density_grid(nr))) then
         if (vintage == iv_opal01) then
            write(run_log_unit,*) routine_id, " T6 or rho outside of table range (62)"
            write(run_log_unit,*) "t6, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
            v%t6_grid(nt)
            write(run_log_unit,*) "slr,r,rho(1),rho(nr):", density_value, &
            density, v%density_grid(1), v%density_grid(nr)
         else
            write(run_log_unit,*) routine_id, " T6/LogR outside of table range (62)"
            write(run_log_unit,*) "slt, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
            v%t6_grid(nt)
            write(run_log_unit,*) "slr, rho(1), rho(nr):", density_value, &
            v%density_grid(1), v%density_grid(nr)
         end if
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
      if (v%t6_index_3 == 0) then
         if (vintage == iv_opal01) then
            write (run_log_unit,'(A, " ihi,ilo,imd",3I5)') routine_id
         else
            write (run_log_unit,'(A, " ihi,ilo,imd",3I5)') routine_id, &
                 hi_idx, lo_idx, mid_idx
         end if
      end if
!     endif

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
            if (vintage == iv_opal01) then
               write(run_log_unit,*) routine_id, "T6/rho in empty region of OPAL 2001 EOS", &
               " table (65)"
            else
               write(run_log_unit,*) routine_id, "T6/log rho in empty region of table (65)"
            end if
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
         if (vintage == iv_opal01) then
            write(run_log_unit,'(A,"Z does not match Z in OEOS01 files you are" &
                 &," using")') routine_id
         else
            write(run_log_unit,'(A,"Z does not match Z in OEOS06 files you are" &
                 &," using")') routine_id
         end if
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
            if (v%x_interp_result(t6_scan_idx,density_scan_idx) > 1.0d+20) then
               if (vintage == iv_opal01) then
                  write(run_log_unit,'(A," problem it ir,l3,k3,iq,ip=", 6I5)') &
                       routine_id, t6_scan_idx, density_scan_idx, v%density_index_3, &
                       v%t6_index_3, v%density_interp_order, v%t6_interp_order
               else
                  write(run_log_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
                       routine_id, t6_scan_idx, density_scan_idx, v%density_index_3, &
                       v%t6_index_3, v%density_interp_order, v%t6_interp_order
               end if
               write(run_log_unit,'(3E12.4)') &
                    (v%x_interp_workspace(x_print_idx,t6_scan_idx,density_scan_idx), &
                    x_print_idx=v%x_index_lo,v%x_index_lo+2)
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
            if (v%x_interp_result(t6_scan_idx,density_scan_idx) > 1.0d+20) then
               if (vintage == iv_opal01) then
                  write(run_log_unit,'(A," problem it ir,l3,k3,iq,ip=", 6I5)') &
                       routine_id, t6_scan_idx, density_scan_idx, v%density_index_3, &
                       v%t6_index_3, v%density_interp_order, v%t6_interp_order
               else
                  write(run_log_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
                       routine_id, t6_scan_idx, density_scan_idx, v%density_index_3, &
                       v%t6_index_3, v%density_interp_order, v%t6_interp_order
               end if
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
      call t6rinteos(v, vintage, density_value, t6_value, ierr)
      if (ierr /= 0) return
      v%eos_output(eos_var_idx) = v%esact
      end do
      if (eos_var_idx > deriv_order) then
      pressure_scale = t6_temperature*density
      if (vintage == iv_opal01) then
         v%eos_output(v%eos_index_inverse(i_opal01_p)) = v%eos_output(v%eos_index_inverse(i_opal01_p))* &
              pressure_scale ! interpolated in p/po
         ! Only slots that this call actually re-interpolated (index <= deriv_order)
         ! are rescaled; the original scaled E and cv unconditionally, so the
         ! deriv_order=1 trial calls from the rho(P,T) inversion compounded the
         ! stale E and cv slots by T6 and moles*R/mu on every call.
         if (v%eos_index_inverse(i_opal01_e) <= deriv_order) &
              v%eos_output(v%eos_index_inverse(i_opal01_e)) = v%eos_output(v%eos_index_inverse(i_opal01_e))* &
              t6_temperature ! interpolated in E/T6
      else
         v%eos_output(v%eos_index_inverse(i_opal_p)) = v%eos_output(v%eos_index_inverse(i_opal_p))* &
              pressure_scale ! interpolated in p/po
         ! Only slots that this call actually re-interpolated (index <= deriv_order)
         ! are rescaled; the original scaled E and cv unconditionally, so the
         ! deriv_order=1 trial calls from the rho(P,T) inversion compounded the
         ! stale E and cv slots by T6 and moles*R/mu on every call.
         if (v%eos_index_inverse(i_opal_e) <= deriv_order) &
              v%eos_output(v%eos_index_inverse(i_opal_e)) = v%eos_output(v%eos_index_inverse(i_opal_e))* &
              t6_temperature ! interpolated in E/T6
      end if
      mean_molecular_weight = gmass01(hydrogen_fraction, opal_shared%table_metal_fraction(vintage), &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_number_fraction)
      if (rad_flag == 1) then
         if (vintage == iv_opal01) then
            call radsub01(v, t6_temperature, density, total_moles, &
                 mean_molecular_weight)
         else
            call radsub06(v, rad_flag, t6_temperature, density, total_moles, &
                 mean_molecular_weight)
         end if
      else
         if (vintage == iv_opal01) then
            if (v%eos_index_inverse(i_opal01_cv) <= deriv_order) &
                 v%eos_output(v%eos_index_inverse(i_opal01_cv)) = v%eos_output(v%eos_index_inverse(i_opal01_cv))* &
                 total_moles*molar_gas_constant_mbcc/mean_molecular_weight
         else
            if (v%eos_index_inverse(i_opal_cv) <= deriv_order) &
                 v%eos_output(v%eos_index_inverse(i_opal_cv)) = v%eos_output(v%eos_index_inverse(i_opal_cv))* &
                 total_moles*molar_gas_constant_mbcc/mean_molecular_weight
         end if
      end if
      return
      end if
      if (vintage == iv_opal01) then
         write(run_log_unit,*) routine_id, " Z does not match Z in OPAL 2001 EOS files", &
              " you are using (66)"
      else
         write(run_log_unit,*) routine_id, " Z does not match Z in EOSdata* files you ", &
              "are using (66)"
      end if
      write(run_log_unit,'("mf,zz(mf)=",I5,E12.4)') v%x_index_lo, &
           v%z_table(v%x_index_lo)
      write(run_log_unit,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') &
           v%t6_interp_order, v%density_interp_order, v%t6_index_3, v%density_index_3, &
           hydrogen_fraction, t6_temperature, density, opal_shared%table_metal_fraction(vintage)
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
! v, then derives the grid spacings interpolate_opal_eos uses. Called
! once per vintage from interpolate_opal_eos. The 2006 tables also
! carry amu and log10 ne columns; the 2001 reader also builds the table
! edge at each T6 row for find_opal_ramp.
!
! NOTE (2006): preserved verbatim from the original -- var_idx, the "j" loop
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
subroutine readcoeos(v, vintage, ierr)
      integer, intent(in) :: vintage

      type(opal_eos_vintage), intent(inout) :: v

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv
! the vintage's table size and the first tabulated T6 row of each density
      integer :: nr, nt
      integer :: t6_index_lo(n_eos_nr_max)

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      integer :: species_read_idx, table_var_idx
      real(dp) :: unused_field
      integer, intent(out) :: ierr

      if (vintage == iv_opal01) then
         nr = n_eos01_nr
         nt = n_eos01_nt
         t6_index_lo = opal01_t6_index_lo
      else
         nr = n_eos06_nr
         nt = n_eos06_nt
         t6_index_lo = opal06_t6_index_lo
      end if
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
         if (vintage == iv_opal06) then
            do var_idx = 1, mv
               v%eos_output(var_idx) = 1.0d0
            end do
         end if
         opal_shared%readco_init_flag(vintage) = opal_flag_set
      end if

      close (2)
! .....read  tables (the file is opened on opal_shared%table_unit by read_controls)
      do x_loop_index = 1, mx
         read (opal_shared%table_unit,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              v%hydrogen_fraction_header(x_loop_index), &
              v%z_table(x_loop_index), &
              v%moles_per_gram_table(x_loop_index), &
              v%mean_molecular_weight_header(x_loop_index)
         if (vintage == iv_opal01) then
            read (opal_shared%table_unit,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
                 &4X,E11.4)') (v%species_fraction_header(x_loop_index,var_idx), &
                 var_idx=1,6)
         else
            read (opal_shared%table_unit,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
                 &4X,E11.4)') (v%species_fraction_header(x_loop_index,species_read_idx), &
                 species_read_idx=1,6)
         end if
         read (opal_shared%table_unit,'(A)') blank_line
         do density_row = 1, nr
            read (opal_shared%table_unit,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 v%temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index,density_row)
            if (record_number /= density_row) then
               if (vintage == iv_opal01) then
                  write (run_log_unit,'(" Data file incorrect: numtot,jcs= ",2I5)') &
                       record_number, density_row
               else
                  write (run_log_unit,'(" OEOS06 Data file incorrect: numtot,jcs= " &
                       &,2I5)') record_number, density_row
               end if
               ierr = 1
               return
            end if
            read(opal_shared%table_unit,'(A)') blank_line
            read(opal_shared%table_unit,'(A)') blank_line
            if (v%temperature_count_used(x_loop_index,density_row) < &
                 t6_index_lo(density_row)) then
               if (vintage == iv_opal01) then
                  write (run_log_unit,'("problem with data files: X=",F6.4, &
                       &" density=",E14.4)') v%hydrogen_fraction_header(x_loop_index), &
                       v%density_grid_table(x_loop_index,density_row)
               else
                  write (run_log_unit,'("Problem with OEOS96 data files: X=",F6.4, &
                       &" density=",E14.4)') v%hydrogen_fraction_header(x_loop_index), &
                       v%density_grid_table(x_loop_index,density_row)
               end if
               ierr = 1
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index,density_row)
               if (t6_row > t6_index_lo(density_row)) then
                  read (opal_shared%table_unit,'(A)') blank_line
                  cycle
               end if
               if (vintage == iv_opal01) then
                  read (opal_shared%table_unit,'(F9.5,1X,F6.2,3E13.5,6F8.4)') &
                       v%t6_list(density_row,t6_row), &
                       v%log10_r_value(density_row,t6_row), &
                       (v%eos_table(x_loop_index,opal_eos_var_order(var_idx), &
                       t6_row,density_row), var_idx=1,9)
               else
                  read (opal_shared%table_unit,'(F11.6,1X,F6.4,E11.4,2E13.6,2E11.3,5F10.6)') &
                       v%t6_list(density_row,t6_row), v%amu_grid(density_row,t6_row), &
                       v%log10_ne_grid(density_row,var_idx), &
                       (v%eos_table(x_loop_index,opal_eos_var_order(table_var_idx), &
                       t6_row,density_row), table_var_idx=1,9)
               end if
            end do
            read(opal_shared%table_unit,'(A)') blank_line
            read(opal_shared%table_unit,'(A)') blank_line
            read(opal_shared%table_unit,'(A)') blank_line
         end do
         read(opal_shared%table_unit,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx) == 0.0d0) then
            if (vintage == iv_opal01) then
               write(run_log_unit,'("READCOEOS01: Error:",I4, &
                    &"-th T6 value is zero")') t6_scan_idx
            else
               write(run_log_unit,'("READCOEOS06: Error:",I4, &
                    &"-th T6 value is zero")') t6_scan_idx
            end if
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
      if (vintage == iv_opal01) then
         ! MHP 7/2003 added ramp between OPAL and other EOS
         ! need edge of table at high rho, fixed T.
         v%t_row_index = 1
         do t6_idx = 1, nt
            v%density_edge_at_t(t6_idx) = v%density_grid(opal01_density_index_edge(t6_idx))
         end do
      end if

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
         if (vintage == iv_opal01) then
            write(run_log_unit,'("Interpolation indices out of range", &
                 &";please report conditions.")')
         else
            write(run_log_unit,'("T6RINTEOS06: Interpolation indices out", &
                 &" of range;please report conditions.")')
         end if
         ierr = 1
         return
      end if
end subroutine t6rinteos

!----------------------------------------------------------------------
! radsub01
!----------------------------------------------------------------------
! OPAL 2001 EOS analogue of radsub, but with a different
! algorithm (preserved verbatim): rather than folding the radiation
! term directly into a single set of derivatives, this computes the
! thermodynamic derivatives BOTH without and with the radiation
! correction, then stores the "with radiation" pressure/energy/
! entropy/specific-heat/chi_rho/chi_t6 directly but only ADDS THE
! DIFFERENCE (with minus without) into gamma1/gamma2_over_gamma2_
! minus1/gamma3_minus1 -- because, per the original's comment, the
! molar_specific_heat calculation isn't accurate enough to give
! accurate gammas directly via the eq. 16.16 Landau & Lifshitz route
! (marked "DIRECT" below) when radiation dominates.
pure subroutine radsub01(v, t6_temperature, density, total_moles, &
     mean_molecular_weight)


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: t6_temperature, density
      real(dp), intent(in) :: total_moles
      real(dp), intent(in) :: mean_molecular_weight

! NOTE: neither literal has a D-suffix in the original (radsub01.f) --
! preserved verbatim: both are parsed as single-precision constants
! and then widened to double precision, which is NOT bit-identical to
! the correctly-rounded double values.
      real(dp), parameter :: rad_const_over_c = 1.8914785e-3
      real(dp), parameter :: molar_gas_constant_mbcc = 83.14510
! --- locals ---
      real(dp) :: rat
      real(dp) :: moles_per_ev, radiation_pressure, &
           radiation_energy, radiation_entropy
      real(dp) :: total_pressure, total_energy, total_entropy
      real(dp) :: chi_rho, chi_t6, molar_specific_heat
      real(dp) :: gamma3_minus1_norad, gamma1_norad, &
           gamma2_over_gamma2_minus1_norad
      real(dp) :: gamma3_minus1, gamma1, gamma2_over_gamma2_minus1

      rat = rad_const_over_c

      moles_per_ev = total_moles*molar_gas_constant_mbcc  ! Mb-cc/unit T6
      radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
      radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
      radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)

!-----Calculate EOS without radiation correction
      total_pressure = v%eos_output(v%eos_index_inverse(i_opal01_p))
      total_energy = v%eos_output(v%eos_index_inverse(i_opal01_e))
      total_entropy = v%eos_output(v%eos_index_inverse(i_opal01_s))
      chi_rho = v%eos_output(v%eos_index_inverse(i_opal01_chi_rho))*v%eos_output(v%eos_index_inverse(i_opal01_p))/ &
           total_pressure
      chi_t6 = (v%eos_output(v%eos_index_inverse(i_opal01_p))*v%eos_output(v%eos_index_inverse(i_opal01_chi_t)))/ &
           total_pressure
      molar_specific_heat = (v%eos_output(v%eos_index_inverse(i_opal01_cv))*moles_per_ev/ &
           mean_molecular_weight)
      gamma3_minus1_norad = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1_norad = chi_rho + chi_t6*gamma3_minus1_norad
      gamma2_over_gamma2_minus1_norad = gamma1_norad/gamma3_minus1_norad
!---- End  no radiation calculation

!---- Calculate EOS with radiation calculation
      total_pressure = v%eos_output(v%eos_index_inverse(i_opal01_p)) + radiation_pressure
      total_energy = v%eos_output(v%eos_index_inverse(i_opal01_e)) + radiation_energy
      total_entropy = v%eos_output(v%eos_index_inverse(i_opal01_s)) + radiation_entropy
      chi_rho = v%eos_output(v%eos_index_inverse(i_opal01_chi_rho))*v%eos_output(v%eos_index_inverse(i_opal01_p))/ &
           total_pressure
      chi_t6 = (v%eos_output(v%eos_index_inverse(i_opal01_p))*v%eos_output(v%eos_index_inverse(i_opal01_chi_t)) &
           + 4.0d0*radiation_pressure)/total_pressure
      molar_specific_heat = (v%eos_output(v%eos_index_inverse(i_opal01_cv))*moles_per_ev/ &
           mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
      gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)                                        ! DIRECT
      gamma1 = chi_rho + chi_t6*gamma3_minus1 !eq 16.16 Landau_Lifshitz (Stat. Mech) ! DIRECT
      gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1             ! DIRECT
!-----End Eos calculations with radiation

!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
      v%eos_output(v%eos_index_inverse(i_opal01_p)) = total_pressure
      v%eos_output(v%eos_index_inverse(i_opal01_e)) = total_energy
      v%eos_output(v%eos_index_inverse(i_opal01_s)) = total_entropy
      v%eos_output(v%eos_index_inverse(i_opal01_cv)) = molar_specific_heat
      v%eos_output(v%eos_index_inverse(i_opal01_chi_rho)) = chi_rho
      v%eos_output(v%eos_index_inverse(i_opal01_chi_t)) = chi_t6
!-----Add difference between EOS with and without radiation.  cvtt
!       calculation is not accurate enough to give accurate results using
!       eq. 16.16 Landau&Lifshitz (SEE line labeled DIRECT)
      v%eos_output(v%eos_index_inverse(i_opal01_gamma1)) = v%eos_output(v%eos_index_inverse(i_opal01_gamma1)) + &
           gamma1 - gamma1_norad
      v%eos_output(v%eos_index_inverse(i_opal01_gamma2_ratio)) = v%eos_output(v%eos_index_inverse(i_opal01_gamma2_ratio)) + &
           gamma2_over_gamma2_minus1 - gamma2_over_gamma2_minus1_norad
      v%eos_output(v%eos_index_inverse(i_opal01_gamma3m1)) = v%eos_output(v%eos_index_inverse(i_opal01_gamma3m1)) + &
           gamma3_minus1 - gamma3_minus1_norad
end subroutine radsub01

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
      write (run_log_unit,'("RHOOFP06: no convergence after ",I0," tries;", &
           & " X, T6, P12, |dP/P| =",4ES12.4," -- falling back")') &
           max_refine, hydrogen_fraction, t6_temperature, pressure_e12, &
           abs((pressure_trial3-pressure_no_rad)/pressure_no_rad)
      rhoofp06 = opal_rho_not_found
      return
      end do refine
end function rhoofp06

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
pure subroutine radsub06(v, rad_flag, t6_temperature, density, total_moles, &
     mean_molecular_weight)


      type(opal_eos_vintage), intent(inout) :: v

      integer, intent(in) :: rad_flag
      real(dp), intent(in) :: t6_temperature, density
      real(dp), intent(in) :: total_moles
      real(dp), intent(in) :: mean_molecular_weight

      real(dp), parameter :: rad_const_over_c = 1.8914785d-3
      real(dp), parameter :: molar_gas_constant_mbcc = 83.14510d0
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

end module opal_eos
