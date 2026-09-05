!----------------------------------------------------------------------
! esac
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original esac.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Interpolates the OPAL 1995 equation of state and its derivatives in
! X (hydrogen_fraction), T6 (temperature in millions of K), and
! density. Values are obtained by a quadratic interpolation at fixed
! T6 at three (or, near table edges/holes, blended pairs of
! overlapping-quadratic) values of density, followed by a quadratic
! interpolation along T6 (done in t6rinterp.f90); results are
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
! radsub.f90). opal_eos_var_order(i),i=1,10 sets the order in which the
! EOS variables are stored in v%eos_output/v%eos_table (by default
! 1,2,...,10, i.e. unpermuted).
subroutine esac(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, *)

      use opal_eos_lib
      use luout_lib
      implicit none

      type(opal_eos_vintage95), intent(inout) :: v

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = n_eos_mx, nr = n_eos95_nr, nt = n_eos95_nt

      double precision :: species_number_fraction(7)
      double precision :: molar_gas_constant_mbcc
      data molar_gas_constant_mbcc/83.1446304d0/
! --- locals ---
      integer :: x_loop_index
      double precision :: t6_value, density_value   ! working copies (SLT/SLR)
      integer :: species_idx, index_idx
      integer :: lo_idx, hi_idx, mid_idx, result_idx
      integer :: x_index_2, x_index_3, x_index_4, x_index_hi
      double precision :: table_sum_2x2, table_sum_3x3, table_sum_3x4, &
           table_sum_4x4
      integer :: eos_var_idx, x_print_idx
      integer :: density_scan_idx, t6_scan_idx
      integer :: recompute_flag, cache_slot
      double precision :: x_interp_weight
      double precision :: pressure_scale
      double precision :: total_moles, ground_state_energy, &
           metal_mole_fraction, mean_molecular_weight
      double precision, external :: quad, gmass

      integer, intent(out) :: ierr

      ierr = 0

      if (deriv_order.gt.10) then
         write(run_log_unit,'(" IORDER CANNOT EXCEED 10")')
      end if
      if ((rad_flag.ne.0) .and. (rad_flag.ne.1)) then
         write(run_log_unit,'(" IRAD MUST BE 0 OR 1")')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if

      t6_value = t6_temperature
      density_value = density

      if (opal_eos%table_loaded_flag.ne.opal_flag_set) then
         opal_eos%table_loaded_flag = opal_flag_set
         do index_idx = 1, 10
            do species_idx = 1, 10
               if (opal_eos_var_order(index_idx).eq.species_idx) &
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

         if (opal_eos%table_metal_fraction(iv_opal95) + hydrogen_fraction - 1.0d-6.gt.1.0d0) then
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
      if ((t6_value.gt.v%t6_grid(1)) .or. (t6_value.lt.v%t6_grid(nt))) then
         write(run_log_unit,'(" T6/LOGR OUTSIDE OF TABLE RANGE (62)")')
         write(run_log_unit,*) v%t6_grid(1), t6_value, v%t6_grid(nt)
         write(run_log_unit,*) v%density_grid(1), density_value, v%density_grid(nr)
         return 1
      end if
      if ((density_value.lt.v%density_grid(1)) .or. (density_value.gt.v%density_grid(nr))) then
         write(run_log_unit,'(" T6/LOGR OUTSIDE OF TABLE RANGE (62)")')
         write(run_log_unit,*) v%t6_grid(1), t6_value, v%t6_grid(nt)
         write(run_log_unit,*) v%density_grid(1), density_value, v%density_grid(nr)
         return 1
      end if
!
!
!
      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx.gt.1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.opal95_x_grid(mid_idx)+1.0d-7) then
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
      if (hydrogen_fraction.lt.1.0d-6) then
         x_index_3 = 1
         x_index_2 = 1
         x_index_4 = 2
         x_index_hi = 1
      end if
      if ((hydrogen_fraction.le.opal95_x_grid(2)+1.0d-7) .or. &
           (hydrogen_fraction.ge.opal95_x_grid(mx-2)-1.0d-7)) x_index_hi = x_index_3
!
      lo_idx = 2
      hi_idx = nr
   do while (hi_idx-lo_idx.gt.1)
         mid_idx = (hi_idx+lo_idx)/2
         if (density_value.eq.v%density_grid(mid_idx)) then
            hi_idx = mid_idx
            exit
         end if
         if (density_value.le.v%density_grid(mid_idx)) then
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
   do while (lo_idx-hi_idx.gt.1)
         mid_idx = (hi_idx+lo_idx)/2
         if (t6_temperature.eq.v%t6_list(1,mid_idx)) then
            lo_idx = mid_idx
            exit
         end if
         if (t6_temperature.le.v%t6_list(1,mid_idx)) then
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
      if (v%t6_index_3.eq.0) then
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
      if (table_sum_3x3.gt.1.0d+30) then
         if (table_sum_2x2.lt.1.0d+25) then
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
            return 1
         end if
      else
      if (table_sum_3x4.lt.1.0d+30) v%t6_interp_order = 3
      if (table_sum_4x4.lt.1.0d+30) v%density_interp_order = 3

      if (t6_temperature.ge.v%t6_list(1,2)+1.0d-7) v%t6_interp_order = 2
      if (density_value.le.v%density_grid(2)+1.0d-15) v%density_interp_order = 2

      if (v%density_index_3.eq.nr) v%density_interp_order = 2

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
      if ((v%z_table(x_index_2).ne.v%z_table(v%x_index_lo)) .or. &
           (v%z_table(x_index_3).ne.v%z_table(v%x_index_lo))) then
         write(run_log_unit,'("Z DOES NOT MATCH Z IN EOSDATA FILES YOU ARE" &
              &," USING")')
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (opal_eos%table_metal_fraction(iv_opal95).ne.v%z_table(v%x_index_lo)) exit
      recompute_flag = 0
      cache_slot = 1
      do density_scan_idx = v%density_index_1, v%density_index_1+v%density_interp_order
         do t6_scan_idx = v%t6_index_1, v%t6_index_1+v%t6_interp_order
            if (x_index_hi.eq.1) then
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
            if (v%x_interp_result(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
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

      if (x_index_4.eq.x_index_hi) then  ! interpolate between quadratics
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
            if (v%x_interp_result(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
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
!       4X4 grid. (T6A(I),I=I1,I1+3),RHO(J),J=J1,J1+3)).PROCEDURE
!       MIXES OVERLAPPING QUADRATICS TO OBTAIN SMOOTHED DERIVATIVES.
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
      if (rad_flag.eq.1) then
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
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return

end subroutine esac
