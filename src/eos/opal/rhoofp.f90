!----------------------------------------------------------------------
! rhoofp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rhoofp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Finds the OPAL 1995 EOS density that reproduces a target pressure
! at a given (hydrogen_fraction, t6_temperature), by secant-method
! refinement of trial densities through esac.f90. Returns -999 (a
! sentinel checked by oeqos.f90) if the point isn't bracketed in the
! table or the iteration fails to converge in 10 tries.
double precision function rhoofp(hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)

      use opal_eos_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = n_eos_mx, nt = n_eos95_nt

! density_index_edge(t6_idx): highest valid density-grid index for
! temperature-grid row t6_idx (a local copy, DATA-initialized here;
! same values as opal_eos_lib's opal_eos%density_index_edge_at_t but
! separate storage).
      integer :: density_index_edge(nt)
      double precision :: rad_const_over_c
      data rad_const_over_c/1.8914785d-3/

! --- locals ---
      integer :: t6_scan_idx
      double precision :: rat, radiation_pressure, pressure_no_rad
      double precision :: hydrogen_fraction_dbg, t6_dbg, density_dbg
      integer :: ideriv_dbg
      integer :: lo_idx, hi_idx, mid_idx
      integer :: x_bisect_idx, t6_bisect_idx
      double precision :: pressure_max, pressure_min
      double precision :: density_trial1, density_trial2, density_trial3
      double precision :: pressure_trial1, pressure_trial2, pressure_trial3
      integer :: refine_count

      data (density_index_edge(t6_scan_idx), t6_scan_idx=1,nt) &
           /7*77, 2*76, 2*74, 2*72, 2*70, 68, 67, 66, 65, 64, 63, 61, &
           60, 59, 58, 57, 55, 54, 53, 52, 51, 2*49, 48, 2*47, 46, &
           2*45, 15*44, 2*37/
      integer, intent(out) :: ierr

      ierr = 0

      rat = rad_const_over_c
! The radiation pressure is deliberately not subtracted from the target
! here (the original's IF(IRAD.EQ.1) PR=4/3*RAT*T6**4 is commented out);
! rad_flag instead adds it to pressure_max/pressure_min below.
      radiation_pressure = 0.0d0
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_eos%table_loaded_flag.ne.opal_flag_set) then
         hydrogen_fraction_dbg = 0.5d0
         t6_dbg = 1.0d0
         density_dbg = 0.001d0
         ideriv_dbg = 1
         call esac(hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              rad_flag, ierr, *999)
         if (ierr /= 0) then
            rhoofp = opal_rho_not_found
            return
         end if
      end if

      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx.gt.1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.opal_eos%x_grid(mid_idx)+1.0d-7) then
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
         if (t6_temperature.eq.opal_eos%t6_list(1,mid_idx)) then
            lo_idx = mid_idx
            exit
         end if
         if (t6_temperature.le.opal_eos%t6_list(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
   end do
      t6_bisect_idx = lo_idx

      pressure_max = opal_eos%eos_table(x_bisect_idx,1,t6_bisect_idx, &
           density_index_edge(t6_bisect_idx))*t6_temperature* &
           opal_eos%density_grid(density_index_edge(t6_bisect_idx)) + &
           rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      pressure_min = opal_eos%eos_table(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           opal_eos%density_grid(1) + rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      if ((pressure_no_rad.gt.1.25d0*pressure_max) .or. &
           (pressure_no_rad.lt.pressure_min)) then
!        requested pressure-temperature not in table
         rhoofp = opal_rho_not_found
         return
      end if

      density_trial1 = opal_eos%density_grid(density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac(hydrogen_fraction, t6_temperature, density_trial1, 1, &
           rad_flag, ierr, *999)
      if (ierr /= 0) then
         rhoofp = opal_rho_not_found
         return
      end if
      pressure_trial1 = opal_eos%eos_output(i_opal_p)
      if (pressure_trial1.gt.pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         density_trial1 = 0.2d0*density_trial1
         if (density_trial1.lt.1.0d-14) density_trial1 = 1.0d-14
         call esac(hydrogen_fraction, t6_temperature, density_trial1, 1, &
              rad_flag, ierr, *999)
         if (ierr /= 0) then
            rhoofp = opal_rho_not_found
            return
         end if
         pressure_trial1 = opal_eos%eos_output(i_opal_p)
      else
         density_trial2 = 5.0d0*density_trial1
         if (density_trial2.gt.opal_eos%density_grid(density_index_edge(t6_bisect_idx))) &
              density_trial2 = opal_eos%density_grid(density_index_edge(t6_bisect_idx))
         call esac(hydrogen_fraction, t6_temperature, density_trial2, 1, &
              rad_flag, ierr, *999)
         if (ierr /= 0) then
            rhoofp = opal_rho_not_found
            return
         end if
         pressure_trial2 = opal_eos%eos_output(i_opal_p)
      end if

      refine_count = 0
      refine: do
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
      call esac(hydrogen_fraction, t6_temperature, density_trial3, 1, &
           rad_flag, ierr, *999)
      if (ierr /= 0) then
         rhoofp = opal_rho_not_found
         return
      end if
      pressure_trial3 = opal_eos%eos_output(i_opal_p)
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad).lt.0.5d-7) then
         rhoofp = density_trial3
         return
      end if
      if (pressure_trial3.gt.pressure_no_rad) then
         density_trial2 = density_trial3
         pressure_trial2 = pressure_trial3
         if (refine_count.lt.11) cycle refine
!        no convergence after 10 tries
         rhoofp = opal_rho_not_found
         return
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
         if (refine_count.lt.11) cycle refine
!        no convergence after 10 tries
         rhoofp = opal_rho_not_found
         return
      end if
      end do refine
  999 continue
!     esac took its alternate return: failed to find rho
      rhoofp = opal_rho_not_found
      return
end function rhoofp
