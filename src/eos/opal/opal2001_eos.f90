!----------------------------------------------------------------------
! opal2001_eos
!----------------------------------------------------------------------
! The OPAL 2001 equation of state (Rogers & Nayfonov 2002), with
! the same two entries as opal95_eos:
!     compute_opal2001_eos  density and the thermodynamic derivatives
!     find_opal2001_ramp    table coverage and the edge-blend weight
! Under these, OPAL's routines under their original names (rhoofp01,
! esac01, t6rinteos01, radsub01, readcoeos01). The tables and the
! interpolation state live in opal_eos_tables (opal01).
module opal2001_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opal_eos_tables, only: find_opal_ramp, gmass01, i_opal01_chi_rho, &
           i_opal01_chi_t, i_opal01_cv, i_opal01_e, i_opal01_gamma1, &
           i_opal01_gamma2_ratio, i_opal01_gamma3m1, i_opal01_p, i_opal01_s, &
           interpolate_opal_t6_rho, iv_opal01, n_eos01_nr, n_eos01_nt, &
           n_eos_mv, n_eos_mx, opal01, opal01_density_index_edge, &
           opal01_t6_index_lo, opal01_x_grid, opal_eos, opal_eos_var_order, &
           opal_eos_vintage, opal_flag_set, opal_rho_not_found, quad
      use physical_constants, only: gas_constant
      use io_units, only: run_log_unit
      use eos_mixture, only: compute_mean_weights
      use math_lib
      implicit none
      private
      public :: compute_opal2001_eos, find_opal2001_ramp

contains

!----------------------------------------------------------------------
! compute_opal2001_eos (formerly oeqos01)
!----------------------------------------------------------------------
! 2001 OPAL equation of state. LLP, August 10, 2003. Given log10(T),
! log10(P), X, and Z, looks up density and the thermodynamic
! derivatives needed by the rest of the EOS machinery from the OPAL
! 2001 tables (via esac01). Called from blend_eos
! (compute_blended_eos); in_table is false when the point is
! outside the OPAL tables.
subroutine compute_opal2001_eos(log10_temperature, temperature, log10_pressure, &
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
! matches P, or when esac01 finds the point outside the table
      logical, intent(out) :: in_table

      ierr = 0
      in_table = .true.

      deriv_order = 9  ! gives all 1st and 2nd order data. See instructions
!                  in esac01.
!     NOTE: rad_flag=0 does not add radiation; rad_flag=1 adds radiation
      rad_flag = 1     ! does add radiation  corrections

      temperature = exp10(log10_temperature)
      pressure = exp10(log10_pressure)
      t_million_k = temperature/1.0d6
      p_e12 = pressure/1.0d12
      if (t_million_k < 0.0020d0 .or. t_million_k > 100.0d0) then
         in_table = .false.
         return
      end if
      hydrogen_fraction_work = hydrogen_fraction

      density_cgs = rhoofp01(opal01, hydrogen_fraction_work, t_million_k, p_e12, &
           rad_flag, ierr)
      if (ierr /= 0) return
      if (density_cgs <= -998.0d0) then
         in_table = .false.
         return
      end if
      density = density_cgs
      log10_density = log10(density)

      call esac01(opal01, hydrogen_fraction_work, t_million_k, density_cgs, &
           deriv_order, rad_flag, ierr, in_table)
      if (.not. in_table) then
         write(run_log_unit, *) 'WARNNING... OPAL TBL FAIL'
         return
      end if
      if (ierr /= 0) return

! The 2001 tables carry no dE/dRho slot, so the slot numbers are one
! lower than in compute_opal2006_eos (see esac01's header).
      dlnrho_dlnp = 1.0d0/opal01%eos_output(i_opal01_chi_rho)
      dlnrho_dlnt = -opal01%eos_output(i_opal01_chi_t)/opal01%eos_output(i_opal01_chi_rho)

      specific_heat_cp = 1.0d6*opal01%eos_output(i_opal01_cv)*opal01%eos_output(i_opal01_gamma1)/opal01%eos_output(i_opal01_chi_rho)
      adiabatic_gradient = 1.0d0/opal01%eos_output(i_opal01_gamma2_ratio)

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

end subroutine compute_opal2001_eos

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

      integer, parameter :: nt = n_eos01_nt, nr = n_eos01_nr

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
      call find_opal_ramp(nt, nr, opal01%t6_grid, opal01%t_row_index, opal01%density_edge_at_t, &
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

! density_index_edge(t6_idx): highest valid density-grid index for
! temperature-grid row t6_idx (a local copy, DATA-initialized here;
! NOTE its values differ from opal_eos_tables's
! v%density_index_edge_at_t in the 115/113 runs -- preserved
! verbatim).
      integer :: density_index_edge(nt)
      real(dp) :: rad_const_over_c
! NOTE: no D-suffix in the original (rhoofp01.f) -- preserved
! verbatim: parsed as a single-precision constant and then widened to
! double precision, which is NOT bit-identical to the correctly-
! rounded double value.
      data rad_const_over_c/1.8914785e-3/

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

      data (density_index_edge(t6_scan_idx), t6_scan_idx=1,nt) &
           /16*169, 168, 167, 166, 165, 2*164, 163, 2*162, 161, 160, 2*159, &
           4*143, 5*137, 6*134, 2*125, 5*123, 2*122, 6*121, 4*119, 8*116, &
           6*115, 8*113, 7*111, 6*110, 34*109, 107, 104, 40*100, 10*99, &
           98, 97, 96, 95, 94, 93, 92/
      integer, intent(out) :: ierr
      logical :: in_table

      ierr = 0

      rat = rad_const_over_c
! The radiation pressure is deliberately not subtracted from the target
! here (the original's IF(IRAD.EQ.1) PR=4/3*RAT*T6**4 is commented out);
! rad_flag instead adds it to pressure_max/pressure_min below.
      radiation_pressure = 0.0d0
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_eos%table_loaded_flag(iv_opal01) /= opal_flag_set) then
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
! esac01
!----------------------------------------------------------------------
! OPAL 2001 EOS analogue of esac (see there for the general
! description of the interpolation scheme). Differences from the
! 1995 version are preserved verbatim, including: iorder capped at 9
! (not 10, since the 2001 tables only carry 9 thermodynamic
! variables); the extra iqu/ipu pre-checks that gate whether the
! "3x4"/"4x4" boundary sums are even accumulated; and the boundary
! reset checking both l3==nr and k3==nt (1995's esac.f only checks
! l3==nr).
!
!   v%eos_output(1) pressure, in megabars (1e12 dyne/cm**2)
!   v%eos_output(2) energy, in 1e12 erg/gm (zero at T6=0)
!   v%eos_output(3) entropy, in units of energy/T6
!   v%eos_output(4) specific heat, dE/dT6 at constant volume
!   v%eos_output(5) dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   v%eos_output(6) dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   v%eos_output(7) gamma1 (Cox & Guili eq 9.88)
!   v%eos_output(8) gamma2/(gamma2-1) (Cox & Guili eq 9.88)
!   v%eos_output(9) gamma3-1 (Cox & Guili eq 9.88)
subroutine esac01(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = n_eos_mx, nr = n_eos01_nr, nt = n_eos01_nt

      real(dp) :: species_number_fraction(7)
      real(dp) :: molar_gas_constant_mbcc
      character(len=15) :: routine_id

! NOTE: this literal has no D0 suffix in the original (esac01.f), unlike
! esac's analogous 83.1446304d0 -- preserved verbatim: it is parsed
! as a single-precision constant and then widened to double precision,
! which is NOT bit-identical to the correctly-rounded double value.
      data molar_gas_constant_mbcc/83.14511/
      data routine_id/"OPALEOS/ESAC01:"/
! --- locals ---
      integer :: x_loop_index_01
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
         ierr = 1
         return
      end if

      t6_value = t6_temperature
      density_value = density
!
      if (opal_eos%table_loaded_flag(iv_opal01) /= opal_flag_set) then
         opal_eos%table_loaded_flag(iv_opal01) = opal_flag_set
         do index_idx = 1, 10
            do species_idx = 1, 10
               if (opal_eos_var_order(index_idx) == species_idx) &
                    v%eos_index_inverse(index_idx) = species_idx
            end do
         end do
         do species_idx = 1, mx
            v%x_grid_copy(species_idx) = opal01_x_grid(species_idx)
         end do
!
! ..... read the data files
         call readcoeos01(v, ierr)
         if (ierr /= 0) return
         opal_eos%table_metal_fraction(iv_opal01) = v%z_table(1)

         if (opal_eos%table_metal_fraction(iv_opal01)+hydrogen_fraction-1.0d-6 > 1) then
            write(run_log_unit,*) routine_id, "Mass fractions exceed unity (61)"
            write(run_log_unit,*) 'Z, XH', opal_eos%table_metal_fraction(iv_opal01), hydrogen_fraction
            ierr = 1
            return
         end if
      end if
!
!
! ..... Determine T6,rho grid points to use in the
!       interpolation.
      if ((t6_value > v%t6_grid(1)) .or. (t6_value < v%t6_grid(nt))) then
         continue
         write(run_log_unit,*) routine_id, " T6 or rho outside of table range (62)"
         write(run_log_unit,*) "t6, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
         v%t6_grid(nt)
         write(run_log_unit,*) "slr,r,rho(1),rho(nr):", density_value, &
         density, v%density_grid(1), v%density_grid(nr)
         in_table = .false.
         return
      end if
      if ((density_value < v%density_grid(1)) .or. (density_value > v%density_grid(nr))) then
         continue
         write(run_log_unit,*) routine_id, " T6 or rho outside of table range (62)"
         write(run_log_unit,*) "t6, t6a(1),t6a(nt):", t6_value, v%t6_grid(1), &
         v%t6_grid(nt)
         write(run_log_unit,*) "slr,r,rho(1),rho(nr):", density_value, &
         density, v%density_grid(1), v%density_grid(nr)
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
         if (hydrogen_fraction <= opal01_x_grid(mid_idx)+1.0d-7) then
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
      if ((hydrogen_fraction <= opal01_x_grid(2)+1.0d-7) .or. &
           (hydrogen_fraction >= opal01_x_grid(mx-2)-1.0d-7)) x_index_hi = x_index_3
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
         write (run_log_unit,'(A, " ihi,ilo,imd",3I5)') routine_id
      end if

!     check to determine if interpolation indices fall within
!     table boundaries.  choose largest allowed size.
      table_sum_2x2 = 0.0d0
      table_sum_3x3 = 0.0d0
      table_sum_3x4 = 0.0d0
      table_sum_4x4 = 0.0d0
      do x_loop_index_01 = v%x_index_lo, x_index_hi
         do density_scan_idx = v%density_index_1, v%density_index_1+1
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+1
               table_sum_2x2 = table_sum_2x2 + &
                    v%eos_table(x_loop_index_01,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         do density_scan_idx = v%density_index_1, v%density_index_1+2
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+2
               table_sum_3x3 = table_sum_3x3 + &
                    v%eos_table(x_loop_index_01,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         if (t6_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+2
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  table_sum_3x4 = table_sum_3x4 + &
                       v%eos_table(x_loop_index_01,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_3x4 = 2.0d+30
         end if
         if (density_order_hi == 3) then
            do density_scan_idx = v%density_index_1, v%density_index_1+3
               do t6_scan_idx = v%t6_index_1, v%t6_index_1+t6_order_hi
                  table_sum_4x4 = table_sum_4x4 + &
                       v%eos_table(x_loop_index_01,1,t6_scan_idx,density_scan_idx)
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
            continue
            write(run_log_unit,*) routine_id, "T6/rho in empty region of OPAL 2001 EOS", &
            " table (65)"
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
      do x_loop_index_01 = v%x_index_lo, x_index_hi

      recompute_flag = 0

!__________
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            v%x_interp_workspace(x_loop_index_01,t6_scan_idx,density_scan_idx) = &
                 v%eos_table(x_loop_index_01,eos_var_idx,t6_scan_idx,density_scan_idx)
            recompute_flag = 1
         end do
      end do
      end do
      if ((v%z_table(x_index_2) /= v%z_table(v%x_index_lo)) .or. &
           (v%z_table(x_index_3) /= v%z_table(v%x_index_lo))) then
         write(run_log_unit,'(A,"Z does not match Z in OEOS01 files you are" &
              &," using")') routine_id
         ierr = 1
         return
      end if
      if (opal_eos%table_metal_fraction(iv_opal01) /= v%z_table(v%x_index_lo)) exit
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
               write(run_log_unit,'(A," problem it ir,l3,k3,iq,ip=", 6I5)') &
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
               write(run_log_unit,'(A," problem it ir,l3,k3,iq,ip=", 6I5)') &
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
      call t6rinteos01(v, density_value, t6_value, ierr)
      if (ierr /= 0) return
      v%eos_output(eos_var_idx) = v%esact
      end do
      if (eos_var_idx > deriv_order) then

      pressure_scale = t6_temperature*density
      v%eos_output(v%eos_index_inverse(i_opal01_p)) = v%eos_output(v%eos_index_inverse(i_opal01_p))* &
           pressure_scale   ! interpolated in p/po
! Only slots that this call actually re-interpolated (index <= deriv_order)
! are rescaled; the original scaled E and cv unconditionally, so the
! deriv_order=1 trial calls from the rho(P,T) inversion compounded the
! stale E and cv slots by T6 and moles*R/mu on every call.
      if (v%eos_index_inverse(i_opal01_e) <= deriv_order) &
           v%eos_output(v%eos_index_inverse(i_opal01_e)) = v%eos_output(v%eos_index_inverse(i_opal01_e))* &
           t6_temperature   ! interpolated in E/T6
      mean_molecular_weight = gmass01(hydrogen_fraction, opal_eos%table_metal_fraction(iv_opal01), &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_number_fraction)
      if (rad_flag == 1) then
         call radsub01(v, t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else
         if (v%eos_index_inverse(i_opal01_cv) <= deriv_order) &
              v%eos_output(v%eos_index_inverse(i_opal01_cv)) = v%eos_output(v%eos_index_inverse(i_opal01_cv))* &
              total_moles*molar_gas_constant_mbcc/mean_molecular_weight
      end if
      return
      end if
      write(run_log_unit,*) routine_id, " Z does not match Z in OPAL 2001 EOS files", &
           " you are using (66)"
      write(run_log_unit,'("mf,zz(mf)=",I5,E12.4)') v%x_index_lo, &
           v%z_table(v%x_index_lo)
      write(run_log_unit,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') &
           v%t6_interp_order, v%density_interp_order, v%t6_index_3, v%density_index_3, &
           hydrogen_fraction, t6_temperature, density, opal_eos%table_metal_fraction(iv_opal01)
      ierr = 1

end subroutine esac01

!----------------------------------------------------------------------
! t6rinteos01
!----------------------------------------------------------------------
! OPAL 2001 EOS analogue of t6rinterp (see there for the general
! description). Interpolates the already X-interpolated table slice
! (v%x_interp_result) in T6 and density to produce v%esact.
!
! interpolate_opal_t6_rho, with the 2001 out-of-range message.
subroutine t6rinteos01(v, slr, slt, ierr)

      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: slr, slt

      integer, intent(out) :: ierr

      ierr = 0

      call interpolate_opal_t6_rho(v, slr, slt)
      if (v%esact > 1.0d+15) then
         write(run_log_unit,'("Interpolation indices out of range", &
              &";please report conditions.")')
         ierr = 1
         return
      end if

end subroutine t6rinteos01

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
subroutine radsub01(v, t6_temperature, density, total_moles, &
     mean_molecular_weight)


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: t6_temperature, density
      real(dp), intent(in) :: total_moles
      real(dp), intent(in) :: mean_molecular_weight

      real(dp) :: rad_const_over_c, molar_gas_constant_mbcc
! NOTE: neither literal has a D-suffix in the original (radsub01.f) --
! preserved verbatim: both are parsed as single-precision constants
! and then widened to double precision, which is NOT bit-identical to
! the correctly-rounded double values.
      data rad_const_over_c/1.8914785e-3/, molar_gas_constant_mbcc/83.14510/
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
! readcoeos01
!----------------------------------------------------------------------
! Reads the OPAL 2001 EOS data tables into opal_eos_tables's opal_eos state
! and derives the auxiliary grid-spacing arrays used by esac01/
! t6rinteos01 for interpolation. Called once from esac01
! (guarded by v%table_loaded_flag).
subroutine readcoeos01(v, ierr)

      type(opal_eos_vintage), intent(inout) :: v

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv, nr = n_eos01_nr, nt = n_eos01_nt

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index_01
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      real(dp) :: unused_field

      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '

      if (opal_eos%readco_init_flag(iv_opal01) /= opal_flag_set) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     v%eos_table(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         opal_eos%readco_init_flag(iv_opal01) = opal_flag_set
      end if

      close (2)
! ..... read  tables (the file is opened on opal_eos%table_unit by read_controls)
      do x_loop_index_01 = 1, mx

         read (opal_eos%table_unit,'(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)') &
              v%hydrogen_fraction_header(x_loop_index_01), &
              v%z_table(x_loop_index_01), &
              v%moles_per_gram_table(x_loop_index_01), &
              v%mean_molecular_weight_header(x_loop_index_01)
         read (opal_eos%table_unit,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (v%species_fraction_header(x_loop_index_01,var_idx), &
              var_idx=1,6)
         read (opal_eos%table_unit,'(A)') blank_line
         do density_row = 1, nr
            read (opal_eos%table_unit,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 v%temperature_count_used(x_loop_index_01,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index_01,density_row)
            if (record_number /= density_row) then
               write (run_log_unit,'(" Data file incorrect: numtot,jcs= ",2I5)') &
                    record_number, density_row
               ierr = 1
               return
            end if
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
            if (v%temperature_count_used(x_loop_index_01,density_row) < &
                 opal01_t6_index_lo(density_row)) then
               write (run_log_unit,'("problem with data files: X=",F6.4, &
                    &" density=",E14.4)') v%hydrogen_fraction_header(x_loop_index_01), &
                    v%density_grid_table(x_loop_index_01,density_row)
               ierr = 1
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index_01,density_row)
               if (t6_row > opal01_t6_index_lo(density_row)) then
                  read (opal_eos%table_unit,'(A)') blank_line
                  cycle
               end if
               read (opal_eos%table_unit,'(F9.5,1X,F6.2,3E13.5,6F8.4)') &
                    v%t6_list(density_row,t6_row), &
                    v%log10_r_value(density_row,t6_row), &
                    (v%eos_table(x_loop_index_01,opal_eos_var_order(var_idx), &
                    t6_row,density_row), var_idx=1,9)
            end do
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
         end do
         read(opal_eos%table_unit,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx) == 0.0d0) then
            write(run_log_unit,'("READCOEOS01: Error:",I4, &
                 &"-th T6 value is zero")') t6_scan_idx
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
! MHP 7/2003 added ramp between OPAL and other EOS
! need edge of table at high rho, fixed T.
      v%t_row_index = 1
      do t6_idx = 1, nt
         v%density_edge_at_t(t6_idx) = v%density_grid(opal01_density_index_edge(t6_idx))
      end do

end subroutine readcoeos01

end module opal2001_eos
