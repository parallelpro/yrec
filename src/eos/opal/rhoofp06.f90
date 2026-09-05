!----------------------------------------------------------------------
! rhoofp06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rhoofp06.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2006 EOS analogue of rhoofp01.f90 (see there and rhoofp.f90
! for the general description), with two preserved differences:
! (1) pressure_max/pressure_min here do NOT add a radiation term
! (rhoofp.f90/rhoofp01.f90 both add rad_flag*4/3*rat*t6**4); (2) the
! trial esac06.f90 calls below pass a hardcoded 0 for the radiation
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
! silently returning -999 (which eqstat turns into a Yale/SCV
! fallback with ierr = 0).
double precision function rhoofp06(v, hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)

      use opal_eos_lib
      use luout_lib, only: run_log_unit
      implicit none

      type(opal_eos_vintage), intent(inout) :: v

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = n_eos_mx, nt = n_eos06_nt

      double precision :: rad_const_over_c
! NOTE: no D-suffix in the original (rhoofp06.f) -- preserved
! verbatim: parsed as a single-precision constant and then widened to
! double precision, which is NOT bit-identical to the correctly-
! rounded double value.
      data rad_const_over_c/1.8914785e-3/
! --- locals ---
      double precision :: rat, radiation_pressure, pressure_no_rad
      double precision :: hydrogen_fraction_dbg, t6_dbg, density_dbg
      integer :: ideriv_dbg
      integer :: lo_idx, hi_idx, mid_idx
      integer :: x_bisect_idx, t6_bisect_idx
      double precision :: pressure_max, pressure_min
      double precision :: density_trial1, density_trial2, density_trial3
      double precision :: pressure_trial1, pressure_trial2, pressure_trial3
      integer :: refine_count, shrink_count
      integer, parameter :: max_refine = 30, max_shrink = 20

      integer, intent(out) :: ierr

      ierr = 0

      rat = rad_const_over_c
      radiation_pressure = 0.0d0
      if (rad_flag.eq.1) radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_eos%table_loaded_flag.ne.opal_flag_set) then
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
              0, ierr, *999)
         if (ierr /= 0) then
            rhoofp06 = opal_rho_not_found
            return
         end if
      end if

      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx.gt.1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.opal06_x_grid(mid_idx)+1.0d-7) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
    end do
      x_bisect_idx = lo_idx

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
      t6_bisect_idx = lo_idx

      pressure_max = v%eos_table(x_bisect_idx,1,t6_bisect_idx, &
           opal06_density_index_edge(t6_bisect_idx))*t6_temperature* &
           v%density_grid(opal06_density_index_edge(t6_bisect_idx))
      pressure_min = v%eos_table(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           v%density_grid(1)
      if ((pressure_no_rad.gt.1.25d0*pressure_max) .or. &
           (pressure_no_rad.lt.pressure_min)) then
!        requested pressure-temperature not in the OPAL 2006 EOS table
         rhoofp06 = opal_rho_not_found
         return
      end if

      density_trial1 = v%density_grid(opal06_density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac06(v, hydrogen_fraction, t6_temperature, density_trial1, 1, 0, ierr, *999)
      if (ierr /= 0) then
         rhoofp06 = opal_rho_not_found
         return
      end if
      pressure_trial1 = v%eos_output(i_opal_p)
      if (pressure_trial1.gt.pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         shrink_count = 0
         shrink: do
            shrink_count = shrink_count + 1
            density_trial1 = 0.2d0*density_trial1
            if (density_trial1.lt.1.0d-14) density_trial1 = 1.0d-14
            call esac06(v, hydrogen_fraction, t6_temperature, density_trial1, 1, 0, ierr, *999)
            if (ierr /= 0) then
               rhoofp06 = opal_rho_not_found
               return
            end if
            pressure_trial1 = v%eos_output(i_opal_p)
! 2026 (bugsweep Batch 2): keep shrinking until the root is bracketed
! (the original did a single x0.2 step and let regula falsi
! extrapolate from an unbracketed pair).
            if (pressure_trial1.le.pressure_no_rad) exit shrink
            if (density_trial1.le.1.0d-14 .or. shrink_count.ge.max_shrink) exit shrink
         end do shrink
      else
         density_trial2 = 5.0d0*density_trial1
         if (density_trial2.gt.v%density_grid(opal06_density_index_edge(t6_bisect_idx))) &
              density_trial2 = v%density_grid(opal06_density_index_edge(t6_bisect_idx))
         call esac06(v, hydrogen_fraction, t6_temperature, density_trial2, 1, 0, ierr, *999)
         if (ierr /= 0) then
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
      call esac06(v, hydrogen_fraction, t6_temperature, density_trial3, 1, 0, ierr, *999)
      if (ierr /= 0) then
         rhoofp06 = opal_rho_not_found
         return
      end if
      pressure_trial3 = v%eos_output(i_opal_p)
! 2026 (bugsweep Batch 2): tolerance restored to the original 0.5d-7
! (was 1.0d-5 from 2025-10-10 to 2026-09; see the header note).
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad).lt.0.5d-7) then
         rhoofp06 = density_trial3
         return
      end if
      if (pressure_trial3.gt.pressure_no_rad) then
         density_trial2 = density_trial3
         pressure_trial2 = pressure_trial3
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
      end if
      if (refine_count.lt.max_refine) cycle refine
      write (run_log_unit,'("RHOOFP06: no convergence after ",I0," tries;", &
           & " X, T6, P12, |dP/P| =",4ES12.4," -- falling back")') &
           max_refine, hydrogen_fraction, t6_temperature, pressure_e12, &
           abs((pressure_trial3-pressure_no_rad)/pressure_no_rad)
      rhoofp06 = opal_rho_not_found
      return
      end do refine
  999 continue
!     esac06 took its alternate return: failed to find rho
      rhoofp06 = opal_rho_not_found
      return
end function rhoofp06
