!----------------------------------------------------------------------
! opal95_eos
!----------------------------------------------------------------------
! The OPAL 1995 equation of state (Rogers, Swenson & Iglesias 1996),
! tabulated in X, T6 and density at one Z.
!     compute_opal95_eos   density and the thermodynamic derivatives at
!                          (log T, log P, X, Z); in_table is false when
!                          the tables do not cover the point
!     find_opal95_ramp     whether a (T, log rho) point is inside the
!                          table, and its weight in the blend with the
!                          Yale/SCV result near the table's edge
! Under these: OPAL's own routines, under their original names: rhoofp
! (the density that gives P), esac (interpolation in X, T6 and rho),
! t6rinterp, radsub (the radiation correction), gmass and readco (the
! table reader, run once on first use). The tables and the
! interpolation state live in opal_eos_tables (opal95).
module opal95_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opal_eos_tables, only: find_opal_ramp, i_opal_chi_rho, i_opal_chi_t, &
           i_opal_cv, i_opal_dedrho, i_opal_e, i_opal_gamma1, &
           i_opal_gamma2_ratio, i_opal_gamma3m1, i_opal_p, i_opal_s, &
           iv_opal95, n_eos95_nr, n_eos95_nt, n_eos_mv, n_eos_mx, opal95, &
           opal95_density_index_edge, opal95_t6_index_lo, opal95_x_grid, &
           opal_eos, opal_eos_var_order, opal_eos_vintage95, opal_flag_set, &
           opal_rho_not_found, quad
      use physical_constants, only: gas_constant
      use io_units, only: run_log_unit
      use eos_mixture, only: compute_mean_weights
      use math_lib
      implicit none
      private
      public :: compute_opal95_eos, find_opal95_ramp

contains

!----------------------------------------------------------------------
! compute_opal95_eos (formerly oeqos)
!----------------------------------------------------------------------
! OPAL (1995) equation of state. Implemented by YC (Feb. 21, 1995).
! Given log10(T), log10(P), X, and Z, looks up density and the
! thermodynamic derivatives needed by the rest of the EOS machinery
! from the OPAL 1995 tables (via esac). Called from
! blend_eos (compute_blended_eos); in_table is false when the point
! is outside the OPAL tables.
subroutine compute_opal95_eos(log10_temperature, temperature, log10_pressure, &
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
! matches P, or when esac finds the point outside the table
      logical, intent(out) :: in_table

      ierr = 0
      in_table = .true.

      temperature = exp10(log10_temperature)
      pressure = exp10(log10_pressure)
      t_million_k = temperature/1.0d6
      p_e12 = pressure/1.0d12
      if (t_million_k < 0.0050d0 .or. t_million_k > 100.0d0) then
         in_table = .false.
         return
      end if
      hydrogen_fraction_work = hydrogen_fraction

      rad_flag = 1
      deriv_order = 10

      density_cgs = rhoofp(opal95, hydrogen_fraction_work, t_million_k, p_e12, rad_flag, ierr)
      if (ierr /= 0) return
      if (density_cgs <= -998.0d0) then
         in_table = .false.
         return
      end if
      density = density_cgs
      log10_density = log10(density)

      call esac(opal95, hydrogen_fraction_work, t_million_k, density_cgs, &
           deriv_order, rad_flag, ierr, in_table)
      if (.not. in_table) then
         write(run_log_unit, *) 'WARNNING... OPAL TBL FAIL'
         return
      end if
      if (ierr /= 0) return

      if (abs((p_e12-opal95%eos_output(i_opal_p))/p_e12) > 0.5d-6) then
         write(run_log_unit,*) p_e12, opal95%eos_output(i_opal_p)
         ierr = 1
         return
      end if
      dlnrho_dlnp = 1.0d0/opal95%eos_output(i_opal_chi_rho)
      dlnrho_dlnt = -opal95%eos_output(i_opal_chi_t)/opal95%eos_output(i_opal_chi_rho)

      specific_heat_cp = 1.0d6*opal95%eos_output(i_opal_cv)*opal95%eos_output(i_opal_gamma1)/opal95%eos_output(i_opal_chi_rho)
      adiabatic_gradient = 1.0d0/opal95%eos_output(i_opal_gamma2_ratio)

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

end subroutine compute_opal95_eos

!----------------------------------------------------------------------
! find_opal95_ramp (formerly eqbound)
!----------------------------------------------------------------------
! MHP 8/98: checks whether a (T, log10 rho) point falls within the
! OPAL 1995 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (see blend_eos.f90's use_opal95_eos branch, which calls this
! immediately after compute_opal95_eos).
!
! The 1995 density check (log10 rho outside -14..5 is outside the
! table), then find_opal_ramp.
subroutine find_opal95_ramp(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)

      integer, parameter :: nt = n_eos95_nt, nr = n_eos95_nr

      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
! --- locals ---
      real(dp) :: t6

      integer, intent(out) :: ierr

      ierr = 0

      t6 = temperature*1.0d-6
!     exit if outside table in rho
      if (log10_density < -14d0 .or. log10_density > 5.0d0) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if
      call find_opal_ramp(nt, nr, opal95%t6_grid, opal95%t_row_index, opal95%density_edge_at_t, &
           opal95_density_index_edge, opal95%density_grid, t6, log10_density, ramp_factor, &
           in_opal_table, needs_ramp, ierr)

end subroutine find_opal95_ramp

!----------------------------------------------------------------------
! rhoofp
!----------------------------------------------------------------------
! Finds the OPAL 1995 EOS density that reproduces a target pressure
! at a given (hydrogen_fraction, t6_temperature), by secant-method
! refinement of trial densities through esac. Returns -999 (a
! sentinel checked by compute_opal95_eos) if the point isn't bracketed in the
! table or the iteration fails to converge in 10 tries.
real(dp) function rhoofp(v, hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)


      type(opal_eos_vintage95), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = n_eos_mx, nt = n_eos95_nt

! density_index_edge(t6_idx): highest valid density-grid index for
! temperature-grid row t6_idx (a local copy, DATA-initialized here;
! same values as opal_eos_tables's v%density_index_edge_at_t but
! separate storage).
      integer :: density_index_edge(nt)
      real(dp) :: rad_const_over_c
      data rad_const_over_c/1.8914785d-3/

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
           /7*77, 2*76, 2*74, 2*72, 2*70, 68, 67, 66, 65, 64, 63, 61, &
           60, 59, 58, 57, 55, 54, 53, 52, 51, 2*49, 48, 2*47, 46, &
           2*45, 15*44, 2*37/
      integer, intent(out) :: ierr
      logical :: in_table

      ierr = 0

      rat = rad_const_over_c
! The radiation pressure is deliberately not subtracted from the target
! here (the original's IF(IRAD.EQ.1) PR=4/3*RAT*T6**4 is commented out);
! rad_flag instead adds it to pressure_max/pressure_min below.
      radiation_pressure = 0.0d0
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_eos%table_loaded_flag(iv_opal95) /= opal_flag_set) then
         hydrogen_fraction_dbg = 0.5d0
         t6_dbg = 1.0d0
         density_dbg = 0.001d0
         ideriv_dbg = 1
         call esac(v, hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              rad_flag, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp = opal_rho_not_found
            return
         end if
      end if

      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx > 1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction <= opal95_x_grid(mid_idx)+1.0d-7) then
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
         rhoofp = opal_rho_not_found
         return
      end if

      density_trial1 = v%density_grid(density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac(v, hydrogen_fraction, t6_temperature, density_trial1, 1, &
           rad_flag, ierr, in_table)
      if (ierr /= 0 .or. .not. in_table) then
         rhoofp = opal_rho_not_found
         return
      end if
      pressure_trial1 = v%eos_output(i_opal_p)
      if (pressure_trial1 > pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         density_trial1 = 0.2d0*density_trial1
         if (density_trial1 < 1.0d-14) density_trial1 = 1.0d-14
         call esac(v, hydrogen_fraction, t6_temperature, density_trial1, 1, &
              rad_flag, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp = opal_rho_not_found
            return
         end if
         pressure_trial1 = v%eos_output(i_opal_p)
      else
         density_trial2 = 5.0d0*density_trial1
         if (density_trial2 > v%density_grid(density_index_edge(t6_bisect_idx))) &
              density_trial2 = v%density_grid(density_index_edge(t6_bisect_idx))
         call esac(v, hydrogen_fraction, t6_temperature, density_trial2, 1, &
              rad_flag, ierr, in_table)
         if (ierr /= 0 .or. .not. in_table) then
            rhoofp = opal_rho_not_found
            return
         end if
         pressure_trial2 = v%eos_output(i_opal_p)
      end if

      refine_count = 0
      refine: do
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
      call esac(v, hydrogen_fraction, t6_temperature, density_trial3, 1, &
           rad_flag, ierr, in_table)
      if (ierr /= 0 .or. .not. in_table) then
         rhoofp = opal_rho_not_found
         return
      end if
      pressure_trial3 = v%eos_output(i_opal_p)
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad) < 0.5d-7) then
         rhoofp = density_trial3
         return
      end if
      if (pressure_trial3 > pressure_no_rad) then
         density_trial2 = density_trial3
         pressure_trial2 = pressure_trial3
         if (refine_count < 11) cycle refine
!        no convergence after 10 tries
         rhoofp = opal_rho_not_found
         return
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
         if (refine_count < 11) cycle refine
!        no convergence after 10 tries
         rhoofp = opal_rho_not_found
         return
      end if
      end do refine
end function rhoofp

!----------------------------------------------------------------------
! esac
!----------------------------------------------------------------------
! Interpolates the OPAL 1995 equation of state and its derivatives in
! X (hydrogen_fraction), T6 (temperature in millions of K), and
! density. Values are obtained by a quadratic interpolation at fixed
! T6 at three (or, near table edges/holes, blended pairs of
! overlapping-quadratic) values of density, followed by a quadratic
! interpolation along T6 (done in t6rinterp); results are
! smoothed by mixing overlapping quadratics.
!
!   v%eos_output(1)  pressure, in megabars (1e12 dyne/cm**2)
!   v%eos_output(2)  energy, in 1e12 erg/gm (zero at T6=0)
!   v%eos_output(3)  entropy, in units of energy/T6
!   v%eos_output(4)  dE/dRho at constant T6
!   v%eos_output(5)  specific heat, dE/dT6 at constant volume
!   v%eos_output(6)  dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   v%eos_output(7)  dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   v%eos_output(8)  gamma1 (Cox & Guili eq 9.88)
!   v%eos_output(9)  gamma2/(gamma2-1) (Cox & Guili eq 9.88)
!   v%eos_output(10) gamma3-1 (Cox & Guili eq 9.88)
!
! deriv_order sets the maximum index of v%eos_output filled in (e.g.
! deriv_order=1 gives just the pressure). rad_flag=0 -> no radiation
! correction; rad_flag=1 -> radiation correction added (via
! radsub). opal_eos_var_order(i),i=1,10 sets the order in which the
! EOS variables are stored in v%eos_output/v%eos_table (by default
! 1,2,...,10, i.e. unpermuted).
subroutine esac(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)


      type(opal_eos_vintage95), intent(inout) :: v

      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = n_eos_mx, nr = n_eos95_nr, nt = n_eos95_nt

      real(dp) :: species_number_fraction(7)
      real(dp) :: molar_gas_constant_mbcc
      data molar_gas_constant_mbcc/83.1446304d0/
! --- locals ---
      integer :: x_loop_index
      real(dp) :: t6_value, density_value   ! working copies (SLT/SLR)
      integer :: species_idx, index_idx
      integer :: lo_idx, hi_idx, mid_idx, result_idx
      integer :: x_index_2, x_index_3, x_index_4, x_index_hi
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

      if (deriv_order > 10) then
         write(run_log_unit,'(" IORDER CANNOT EXCEED 10")')
      end if
      if ((rad_flag /= 0) .and. (rad_flag /= 1)) then
         write(run_log_unit,'(" IRAD MUST BE 0 OR 1")')
         ierr = 1
         return
      end if

      t6_value = t6_temperature
      density_value = density

      if (opal_eos%table_loaded_flag(iv_opal95) /= opal_flag_set) then
         opal_eos%table_loaded_flag(iv_opal95) = opal_flag_set
         do index_idx = 1, 10
            do species_idx = 1, 10
               if (opal_eos_var_order(index_idx) == species_idx) &
                    v%eos_index_inverse(index_idx) = species_idx
            end do
         end do
         do species_idx = 1, mx
            v%x_grid_copy(species_idx) = opal95_x_grid(species_idx)
         end do
!
! ..... read the data files
         call readco(v, ierr)
         if (ierr /= 0) return
         opal_eos%table_metal_fraction(iv_opal95) = v%z_table(1)

         if (opal_eos%table_metal_fraction(iv_opal95) + hydrogen_fraction - 1.0d-6 > 1.0d0) then
            write(run_log_unit,'(" MASS FRACTIONS EXCEED UNITY (61)")')
            write(run_log_unit,*) opal_eos%table_metal_fraction(iv_opal95), hydrogen_fraction
            ierr = 1
            return
         end if
      end if
!
!
! ..... determine T6,rho grid points to use in the
!       interpolation.
      if ((t6_value > v%t6_grid(1)) .or. (t6_value < v%t6_grid(nt))) then
         write(run_log_unit,'(" T6/LOGR OUTSIDE OF TABLE RANGE (62)")')
         write(run_log_unit,*) v%t6_grid(1), t6_value, v%t6_grid(nt)
         write(run_log_unit,*) v%density_grid(1), density_value, v%density_grid(nr)
         in_table = .false.
         return
      end if
      if ((density_value < v%density_grid(1)) .or. (density_value > v%density_grid(nr))) then
         write(run_log_unit,'(" T6/LOGR OUTSIDE OF TABLE RANGE (62)")')
         write(run_log_unit,*) v%t6_grid(1), t6_value, v%t6_grid(nt)
         write(run_log_unit,*) v%density_grid(1), density_value, v%density_grid(nr)
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
         if (hydrogen_fraction <= opal95_x_grid(mid_idx)+1.0d-7) then
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
         x_index_3 = 1
         x_index_2 = 1
         x_index_4 = 2
         x_index_hi = 1
      end if
      if ((hydrogen_fraction <= opal95_x_grid(2)+1.0d-7) .or. &
           (hydrogen_fraction >= opal95_x_grid(mx-2)-1.0d-7)) x_index_hi = x_index_3
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
      if (v%t6_index_3 == 0) then
         write(run_log_unit,'(" IHI,ILO,IMD",3I5)')
      end if

!     check to determine if interpolation indices fall within
!     table boundaries.  choose largest allowed size.
      table_sum_2x2 = 0.0d0
      table_sum_3x3 = 0.0d0
      table_sum_3x4 = 0.0d0
      table_sum_4x4 = 0.0d0
      do x_loop_index = v%x_index_lo, v%x_index_lo+3
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
         do density_scan_idx = v%density_index_1, v%density_index_1+2
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+3
               table_sum_3x4 = table_sum_3x4 + &
                    v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         do density_scan_idx = v%density_index_1, v%density_index_1+3
            do t6_scan_idx = v%t6_index_1, v%t6_index_1+3
               table_sum_4x4 = table_sum_4x4 + &
                    v%eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
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
            write(run_log_unit,'("T6/LOG RHO IN EMPTY REGION OF TABLE (65)")')
            write(run_log_unit,'("XH,T6,R=", 3E12.4)') hydrogen_fraction, &
            t6_temperature, density
            in_table = .false.
            return
         end if
      else
      if (table_sum_3x4 < 1.0d+30) v%t6_interp_order = 3
      if (table_sum_4x4 < 1.0d+30) v%density_interp_order = 3

      if (t6_temperature >= v%t6_list(1,2)+1.0d-7) v%t6_interp_order = 2
      if (density_value <= v%density_grid(2)+1.0d-15) v%density_interp_order = 2

      if (v%density_index_3 == nr) v%density_interp_order = 2

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
         write(run_log_unit,'("Z DOES NOT MATCH Z IN EOSDATA FILES YOU ARE" &
              &," USING")')
         ierr = 1
         return
      end if
      if (opal_eos%table_metal_fraction(iv_opal95) /= v%z_table(v%x_index_lo)) exit
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
               write(run_log_unit,'(" PROBLEM IT IR,L3,K3,IQ,IP=", 6I5)') &
                    t6_scan_idx, density_scan_idx, v%density_index_3, &
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
               write(run_log_unit,'(" PROBLEM IT IR,L3,K3,IQ,IP=", 6I5)') &
                    t6_scan_idx, density_scan_idx, v%density_index_3, &
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
! ..... completed X interpolation. Now interpolate T6 and RHO on a
!       4x4 grid. (T6a(I),I=I1,I1+3),rho(j),j=J1,J1+3)).Procedure
!       mixes overlapping quadratics to obtain smoothed derivatives.
!
!
      call t6rinterp(v, density_value, t6_value, ierr)
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
! dE/dRho (slot 4) is deliberately not rescaled here (YCK commented out
! EOS(IRI(4))=EOS(IRI(4))/SQRT(R*T6) in the original).
      mean_molecular_weight = gmass(hydrogen_fraction, opal_eos%table_metal_fraction(iv_opal95), &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_number_fraction)
      if (rad_flag == 1) then
         call radsub(v, t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else
         if (v%eos_index_inverse(i_opal_cv) <= deriv_order) &
              v%eos_output(v%eos_index_inverse(i_opal_cv)) = v%eos_output(v%eos_index_inverse(i_opal_cv))* &
              total_moles*molar_gas_constant_mbcc/mean_molecular_weight
      end if
      return
      end if
      write(run_log_unit,'(" Z DOES NOT MATCH Z IN EOSDATA* FILES YOU ARE", &
           &" USING (66)")')
      write(run_log_unit,'("MF,ZZ(MF)=",I5,E12.4)') v%x_index_lo, &
           v%z_table(v%x_index_lo)
      write(run_log_unit,'("  IQ,IP,K3,L3,XH,T6,R,Z= ",4I5,4E12.4)') &
           v%t6_interp_order, v%density_interp_order, v%t6_index_3, v%density_index_3, &
           hydrogen_fraction, t6_temperature, density, opal_eos%table_metal_fraction(iv_opal95)
      ierr = 1

end subroutine esac

!----------------------------------------------------------------------
! t6rinterp
!----------------------------------------------------------------------
! Second stage of the OPAL 1995 EOS interpolation: given the already
! X-interpolated table slice (v%x_interp_result, computed by esac's
! call site right before this), interpolates in T6 and density on the
! surrounding 3x3 (or 4x4, smoothed by mixing overlapping quadratics)
! grid to produce v%esact, the interpolated value of the current
! thermodynamic variable at (slt, slr) = (T6, density).
subroutine t6rinterp(v, slr, slt, ierr)


      type(opal_eos_vintage95), intent(inout) :: v

      real(dp), intent(in) :: slr, slt

! --- locals ---
      integer :: hi_loop_count, recompute_flag, cache_slot, t6_grid_idx
      real(dp) :: esactq, esact2, esactq2, dix, dix2

      integer, intent(out) :: ierr

      ierr = 0

      hi_loop_count = 0
      recompute_flag = 0

      do t6_grid_idx = v%t6_index_1, v%t6_index_1 + v%t6_interp_order
         cache_slot = 1
         hi_loop_count = hi_loop_count + 1
         v%rho_interp_lo(hi_loop_count) = quad(v%quad, recompute_flag, cache_slot, &
              slr, v%x_interp_result(t6_grid_idx,v%density_index_1), &
              v%x_interp_result(t6_grid_idx,v%density_index_2), &
              v%x_interp_result(t6_grid_idx,v%density_index_3), &
              v%density_grid(v%density_index_1), v%density_grid(v%density_index_2), &
              v%density_grid(v%density_index_3))
         if (v%density_interp_order == 3) then
            cache_slot = 2
            v%rho_interp_hi(hi_loop_count) = quad(v%quad, recompute_flag, cache_slot, &
                 slr, v%x_interp_result(t6_grid_idx,v%density_index_2), &
                 v%x_interp_result(t6_grid_idx,v%density_index_3), &
                 v%x_interp_result(t6_grid_idx,v%density_index_4), &
                 v%density_grid(v%density_index_2), v%density_grid(v%density_index_3), &
                 v%density_grid(v%density_index_4))
         end if
         recompute_flag = 1
      end do

      recompute_flag = 0
      cache_slot = 1
! ..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      v%esact = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(1), &
           v%rho_interp_lo(2), v%rho_interp_lo(3), v%t6_grid(v%t6_index_1), &
           v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      if (v%density_interp_order == 3) then
! .....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
         esactq = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_hi(1), &
              v%rho_interp_hi(2), v%rho_interp_hi(3), v%t6_grid(v%t6_index_1), &
              v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      end if
      if (v%t6_interp_order == 3) then
! .....    eos(i) in lower-left 3x3.
         esact2 = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(2), &
              v%rho_interp_lo(3), v%rho_interp_lo(4), v%t6_grid(v%t6_index_2), &
              v%t6_grid(v%t6_index_3), v%t6_grid(v%t6_index_4))
! .....    eos(i) smoothed in left 3x4
         dix = (v%t6_grid(v%t6_index_3) - slt)*v%t6_grid_spacing_inv(v%t6_index_3)
         v%esact = v%esact*dix + esact2*(1.0d0 - dix)
      end if
! NOTE: preserved verbatim from the original -- this block is a
! sibling of (not nested inside) the "v%t6_interp_order.eq.3" block
! above, so if v%density_interp_order.eq.3 but v%t6_interp_order.ne.3,
! dix is not computed in this call: it is an ordinary (non-SAVE'd)
! local, so it then holds the -finit-local-zero value 0.0d0.
      if (v%density_interp_order == 3) then
! .....     eos(i) in upper-right 3x3.
         esactq2 = quad(v%quad, recompute_flag, cache_slot, slt, &
              v%rho_interp_hi(2), v%rho_interp_hi(3), v%rho_interp_hi(4), &
              v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3), &
              v%t6_grid(v%t6_index_4))
         esactq = esactq*dix + esactq2*(1.0d0 - dix)
      end if
      if (v%density_interp_order == 3) then
         dix2 = (v%density_grid(v%density_index_3) - slr)* &
              v%density_grid_spacing_inv(v%density_index_3)
         if (v%t6_interp_order == 3) then
! .....        eos(i) smoothed in both log(T6) and log(R)
            v%esact = v%esact*dix2 + esactq*(1.0d0 - dix2)
         end if
      end if
      if (v%esact > 1.0d+15) then
         write(run_log_unit,'("INTERPOLATION INDICES OUT OF RANGE", &
              &";PLEASE REPORT CONDITIONS.")')
         ierr = 1
         return
      end if

end subroutine t6rinterp

!----------------------------------------------------------------------
! radsub
!----------------------------------------------------------------------
! Adds the radiation-pressure/energy/entropy correction to the OPAL
! 1995 EOS values already in v%eos_output, and recomputes
! the derived thermodynamic derivatives (chi_rho, chi_t6, gamma1,
! gamma2/(gamma2-1), gamma3-1) including that correction. Called from
! esac when rad_flag=1.
!
! Physical constants (see comments below): na*k = 6.0221367e23 *
! 1.380658e-16 erg/degree K = 8.314511e7 erg/degree K = 8.314511e7 *
! 11604.5 erg/eV = 0.9648575e12 erg/eV. unit_factor = 0.9648575
! (latest physical constants); unit_factor_legacy = 0.9652 (older
! units, still used elsewhere in this EOS code). In these units,
! energy/density is in Mb-cc/gm and pressure is in units of 1e12
! bars = Mb. sigma (Stefan-Boltzmann) = 5.67051e-5 erg/(s cm**2 K**4),
! c = 2.99792458e10 cm/sec; rad_const_over_c = sigma/c, converted to
! units of 1e6 K (T6) via *1e24.
subroutine radsub(v, t6_temperature, density, total_moles, &
     mean_molecular_weight)


      type(opal_eos_vintage95), intent(inout) :: v

      real(dp), intent(in) :: t6_temperature, density
      real(dp), intent(in) :: total_moles
      real(dp), intent(in) :: mean_molecular_weight

      real(dp) :: rad_const_over_c, unit_factor, unit_factor_legacy, &
           molar_gas_constant_mbcc
      data unit_factor/0.9648575d0/, unit_factor_legacy/0.9652d0/, &
           rad_const_over_c/1.8914785d-3/, molar_gas_constant_mbcc/83.1446304d0/
! --- locals ---
      real(dp) :: moles_per_ev, radiation_pressure, &
           radiation_energy, radiation_entropy
      real(dp) :: total_pressure, total_energy, total_entropy
      real(dp) :: chi_rho, chi_t6
      real(dp) :: molar_specific_heat, gamma3_minus1, gamma1, &
           gamma2_over_gamma2_minus1
      real(dp) :: de_drho_at_t
      real(dp) :: unit_ratio

      moles_per_ev = total_moles*molar_gas_constant_mbcc  ! Mb-cc/unit T6

      radiation_pressure = 4.0d0/3.0d0*rad_const_over_c*t6_temperature**4   ! Mb
      radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
      radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)
      total_pressure = v%eos_output(v%eos_index_inverse(i_opal_p)) + radiation_pressure
      total_energy = v%eos_output(v%eos_index_inverse(i_opal_e)) + radiation_energy
      total_entropy = v%eos_output(v%eos_index_inverse(i_opal_s)) + radiation_entropy
      chi_rho = v%eos_output(v%eos_index_inverse(i_opal_chi_rho))*v%eos_output(v%eos_index_inverse(i_opal_p))/ &
           total_pressure
      chi_t6 = (v%eos_output(v%eos_index_inverse(i_opal_p))*v%eos_output(v%eos_index_inverse(i_opal_chi_t)) &
           + 4.0d0*radiation_pressure)/total_pressure
      molar_specific_heat = (v%eos_output(v%eos_index_inverse(i_opal_cv))*moles_per_ev/ &
           mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
      gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1 = chi_rho + chi_t6*gamma3_minus1
      gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1

!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
      de_drho_at_t = v%eos_output(v%eos_index_inverse(i_opal_dedrho)) - radiation_energy/density
      unit_ratio = unit_factor/unit_factor_legacy
      total_pressure = total_pressure*unit_ratio
! MHP 10/02 EN is never used; should this be ET=ET*REVISE???
!      EN=EN*REVISE
      total_entropy = total_entropy*unit_ratio
! DEDRHOA is never used; should this be DEDRHOAT=DEDRHOAT*REVISE????
!      Dedrhoa=dedrhoa*revise
      v%eos_output(v%eos_index_inverse(i_opal_p)) = total_pressure
      v%eos_output(v%eos_index_inverse(i_opal_e)) = total_energy
      v%eos_output(v%eos_index_inverse(i_opal_s)) = total_entropy
      v%eos_output(v%eos_index_inverse(i_opal_dedrho)) = de_drho_at_t
      v%eos_output(v%eos_index_inverse(i_opal_cv)) = molar_specific_heat
      v%eos_output(v%eos_index_inverse(i_opal_chi_rho)) = chi_rho
      v%eos_output(v%eos_index_inverse(i_opal_chi_t)) = chi_t6
      v%eos_output(v%eos_index_inverse(i_opal_gamma1)) = gamma1
      v%eos_output(v%eos_index_inverse(i_opal_gamma2_ratio)) = gamma2_over_gamma2_minus1
      v%eos_output(v%eos_index_inverse(i_opal_gamma3m1)) = gamma3_minus1
end subroutine radsub

!----------------------------------------------------------------------
! gmass
!----------------------------------------------------------------------
! Called from esac (OPAL 1995 EOS reader) to get the total mass
! per particle (in amu) for a mix of fully-ionized H, He, C, N, O,
! Ne at composition (hydrogen_fraction, metal_fraction), assuming a
! fixed solar-relative C:N:O:Ne mix within the metals. Also returns
! the ground-state ionization energy and number (mole) fractions of each
! species, used by esac for the radiation and specific
! heat corrections.
real(dp) function gmass(hydrogen_fraction, metal_fraction, &
     total_moles, ground_state_energy, metal_mole_fraction, &
     species_number_fraction)


      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: total_moles, ground_state_energy, &
           metal_mole_fraction
! species_number_fraction(1:6): number (mole) fractions of Ne, O, N, C,
! He, H (same order as atomic_number); slot 7 is never assigned.
      real(dp), intent(out) :: species_number_fraction(:)

! atomic_number(i): atomic number (= electrons per fully-ionized atom)
! of species i, for i=1,6 = Ne, O, N, C, He, H.
      real(dp) :: atomic_number(6)
! ionization_energy(i): total ground-state ionization energy (eV) of
! species i, same ordering as atomic_number.
      real(dp) :: ionization_energy(6)
! atomic_weight(i): atomic weight (amu), i=1 electron, 2 Ne, 3 O, 4 N,
! 5 C, 6 He, 7 H.
      real(dp) :: atomic_weight(7)
! solar-relative mole-fraction-per-unit-Z of C, N, O, Ne in the metal
! mix (fixed composition assumed for all Z).
      real(dp) :: carbon_mix_fraction, nitrogen_mix_fraction, &
           oxygen_mix_fraction, neon_mix_fraction
      real(dp) :: carbon_moles, nitrogen_moles, oxygen_moles, &
           neon_moles
      real(dp) :: hydrogen_moles, helium_moles, total_moles_raw
! electron_mole_excess: electrons per mole of ions. total_moles sums
! (1 + Z_i) over the species mole fractions (which sum to 1), i.e.
! ions plus their electrons, so total_moles - 1 is the electron count.
      real(dp) :: electron_mole_excess
      integer :: species_idx

      data (ionization_energy(species_idx), species_idx=1,6) &
           /-3388.637d0, -1970.918d0, -1431.717d0, -993.2303d0, &
           -76.2315d0, -15.29409d0/
      data (atomic_number(species_idx), species_idx=1,6) &
           /10.0d0, 8.0d0, 7.0d0, 6.0d0, 2.0d0, 1.0d0/
      carbon_mix_fraction = 0.247137766d0
      nitrogen_mix_fraction = 0.0620782d0
      oxygen_mix_fraction = 0.52837118d0
      neon_mix_fraction = 0.1624188d0
      atomic_weight(7) = 1.0079d0
      atomic_weight(6) = 4.0026d0
      atomic_weight(5) = 12.011d0
      atomic_weight(4) = 14.0067d0
      atomic_weight(3) = 15.9994d0
      atomic_weight(2) = 20.179d0
      atomic_weight(1) = 0.00054858d0
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
              ionization_energy(species_idx)*species_number_fraction(species_idx)
         total_moles = total_moles + (1.0d0 + atomic_number(species_idx))* &
              species_number_fraction(species_idx)
      end do
      electron_mole_excess = total_moles - 1.0d0
      gmass = electron_mole_excess*atomic_weight(1)
      do species_idx = 2, 7
         gmass = gmass + atomic_weight(species_idx)* &
              species_number_fraction(species_idx - 1)
      end do

end function gmass

!----------------------------------------------------------------------
! readco
!----------------------------------------------------------------------
! Reads the OPAL 1995 EOS data tables (one block per tabulated X
! value) into opal_eos_tables's opal_eos state and derives the auxiliary
! grid-spacing arrays used by t6rinterp for
! interpolation. Called once from esac (guarded by
! opal_eos%table_loaded_flag).
subroutine readco(v, ierr)

      type(opal_eos_vintage95), intent(inout) :: v

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv, nr = n_eos95_nr, nt = n_eos95_nt

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      real(dp) :: unused_field

      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '

      if (opal_eos%readco_init_flag(iv_opal95) /= opal_flag_set) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     v%eos_table(x_idx, var_idx, t6_idx, r_idx) = 1.0d+35
                  end do
               end do
            end do
         end do
         opal_eos%readco_init_flag(iv_opal95) = opal_flag_set
      end if

! ..... read  tables (the file is opened on opal_eos%table_unit by read_controls)
      do x_loop_index = 1, mx

         read (opal_eos%table_unit,'(3X,F6.4,3X,F6.4,11X,F10.7,17X,F10.7)') &
              v%hydrogen_fraction_header(x_loop_index), &
              v%z_table(x_loop_index), &
              v%moles_per_gram_table(x_loop_index), &
              v%mean_molecular_weight_header(x_loop_index)
         read (opal_eos%table_unit,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)') (v%species_fraction_header(x_loop_index,var_idx), &
              var_idx=1,6)
         read (opal_eos%table_unit,'(A)') blank_line
         do density_row = 1, nr
            read (opal_eos%table_unit,'(2I5,2F12.7,17X,E15.7)') record_number, &
                 v%temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index,density_row)
            if (record_number /= density_row) then
               write(run_log_unit,'(" DATA FILE INCORRECT: NUMTOT,JCS= ",2I5)') &
                    record_number, density_row
               ierr = 1
               return
            end if
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
            if (v%temperature_count_used(x_loop_index,density_row) < &
                 opal95_t6_index_lo(density_row)) then
               write(run_log_unit,'("PROBLEM WITH DATA FILES: X=",F6.4," DENSITY=", &
                    &E14.4)') v%hydrogen_fraction_header(x_loop_index), &
                    v%density_grid_table(x_loop_index,density_row)
               ierr = 1
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index,density_row)
               if (t6_row > opal95_t6_index_lo(density_row)) then
                  read (opal_eos%table_unit,'(A)') blank_line
                  cycle
               end if
               read (opal_eos%table_unit,'(F8.4,1X,F6.2,3E13.5,E11.3,6F8.4)') &
                    v%t6_list(density_row,t6_row), &
                    v%log10_r_value(density_row,t6_row), &
                    (v%eos_table(x_loop_index,opal_eos_var_order(var_idx), &
                    t6_row,density_row), var_idx=1,10)
            end do
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
            read(opal_eos%table_unit,'(A)') blank_line
         end do
         read(opal_eos%table_unit,'(A)') blank_line
      end do

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx) == 0.0d0) exit
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
!  MHP 8/98 added ramp between OPAL and other EOS
!  need edge of table at high rho, fixed T.
      v%t_row_index = 1
      do t6_idx = 1, nt
         v%density_edge_at_t(t6_idx) = v%density_grid(opal95_density_index_edge(t6_idx))
      end do
end subroutine readco

end module opal95_eos
