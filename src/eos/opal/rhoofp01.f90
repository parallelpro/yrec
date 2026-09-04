!----------------------------------------------------------------------
! rhoofp01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rhoofp01.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2001 EOS analogue of rhoofp.f90 (see there for the general
! description). The convergence tolerance below (0.5d-7) matches
! rhoofp.f90's 1995 version and, since bug-sweep Batch 2 (2026),
! rhoofp06.f90 again (it was 1.0d-5 there for a while; see its
! header note).
double precision function rhoofp01(hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)

      use opal_eos_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = n_eos_mx, nt = n_eos01_nt

! density_index_edge(t6_idx): highest valid density-grid index for
! temperature-grid row t6_idx (a local copy, DATA-initialized here;
! NOTE its values differ from opal_eos_lib's
! opal_eos%density_index_edge_at_t_01 in the 115/113 runs -- preserved
! verbatim).
      integer :: density_index_edge(nt)
      double precision :: rad_const_over_c
! NOTE: no D-suffix in the original (rhoofp01.f) -- preserved
! verbatim: parsed as a single-precision constant and then widened to
! double precision, which is NOT bit-identical to the correctly-
! rounded double value.
      data rad_const_over_c/1.8914785e-3/

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
           /16*169, 168, 167, 166, 165, 2*164, 163, 2*162, 161, 160, 2*159, &
           4*143, 5*137, 6*134, 2*125, 5*123, 2*122, 6*121, 4*119, 8*116, &
           6*115, 8*113, 7*111, 6*110, 34*109, 107, 104, 40*100, 10*99, &
           98, 97, 96, 95, 94, 93, 92/
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
         call esac01(hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              rad_flag, ierr, *999)
         if (ierr /= 0) then
            rhoofp01 = opal_rho_not_found
            return
         end if
      end if

      lo_idx = 2
      hi_idx = mx
    do while (hi_idx-lo_idx.gt.1)
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.opal_eos%x_grid_01(mid_idx)+1.0d-7) then
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
         if (t6_temperature.eq.opal_eos%t6_list_01(1,mid_idx)) then
            lo_idx = mid_idx
            exit
         end if
         if (t6_temperature.le.opal_eos%t6_list_01(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
   end do
      t6_bisect_idx = lo_idx

      pressure_max = opal_eos%eos_table_01(x_bisect_idx,1,t6_bisect_idx, &
           density_index_edge(t6_bisect_idx))*t6_temperature* &
           opal_eos%density_grid_01(density_index_edge(t6_bisect_idx)) + &
           rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      pressure_min = opal_eos%eos_table_01(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           opal_eos%density_grid_01(1) + rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      if ((pressure_no_rad.gt.1.25d0*pressure_max) .or. &
           (pressure_no_rad.lt.pressure_min)) then
!        requested pressure-temperature not in table
         rhoofp01 = opal_rho_not_found
         return
      end if

      density_trial1 = opal_eos%density_grid_01(density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac01(hydrogen_fraction, t6_temperature, density_trial1, 1, &
           rad_flag, ierr, *999)
      if (ierr /= 0) then
         rhoofp01 = opal_rho_not_found
         return
      end if
      pressure_trial1 = opal_eos%eos_output_01(i_opal01_p)
      if (pressure_trial1.gt.pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         density_trial1 = 0.2d0*density_trial1
         if (density_trial1.lt.1.0d-14) density_trial1 = 1.0d-14
         call esac01(hydrogen_fraction, t6_temperature, density_trial1, 1, &
              rad_flag, ierr, *999)
         if (ierr /= 0) then
            rhoofp01 = opal_rho_not_found
            return
         end if
         pressure_trial1 = opal_eos%eos_output_01(i_opal01_p)
      else
         density_trial2 = 5.0d0*density_trial1
         if (density_trial2.gt.opal_eos%density_grid_01(density_index_edge(t6_bisect_idx))) &
              density_trial2 = opal_eos%density_grid_01(density_index_edge(t6_bisect_idx))
         call esac01(hydrogen_fraction, t6_temperature, density_trial2, 1, &
              rad_flag, ierr, *999)
         if (ierr /= 0) then
            rhoofp01 = opal_rho_not_found
            return
         end if
         pressure_trial2 = opal_eos%eos_output_01(i_opal01_p)
      end if

      refine_count = 0
      refine: do
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
      call esac01(hydrogen_fraction, t6_temperature, density_trial3, 1, &
           rad_flag, ierr, *999)
      if (ierr /= 0) then
         rhoofp01 = opal_rho_not_found
         return
      end if
      pressure_trial3 = opal_eos%eos_output_01(i_opal01_p)
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad).lt.0.5d-7) then
         rhoofp01 = density_trial3
         return
      end if
      if (pressure_trial3.gt.pressure_no_rad) then
         density_trial2 = density_trial3
         pressure_trial2 = pressure_trial3
         if (refine_count.lt.11) cycle refine
!        no convergence after 10 tries
         rhoofp01 = opal_rho_not_found
         return
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
         if (refine_count.lt.11) cycle refine
!        no convergence after 10 tries
         rhoofp01 = opal_rho_not_found
         return
      end if
      end do refine
  999 continue
!     esac01 took its alternate return: failed to find rho
      rhoofp01 = opal_rho_not_found
      return
end function rhoofp01
