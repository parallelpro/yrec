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
      use luout_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = 5, mv = 10, nr = 77, nt = 56






! density_index_edge(t6_idx): highest valid density-grid index for
! temperature-grid row t6_idx (a local copy, DATA-initialized here;
! not the same storage as readco.f90's density_edge_at_t/
! density_index_edge_at_t in common/rmpopeos/, though it encodes the
! same edge-of-table concept).
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

      save

      integer, intent(out) :: ierr

      ierr = 0

      rat = rad_const_over_c
      radiation_pressure = 0.0d0
!      IF(IRAD .EQ. 1) PR=4.D0/3.D0*RAT*T6**4   ! MB
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_eos%table_loaded_flag.ne.12345678) then
         hydrogen_fraction_dbg = 0.5d0
         t6_dbg = 1.0d0
         density_dbg = 0.001d0
         ideriv_dbg = 1
         call esac(hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              rad_flag, ierr, *999)
         if (ierr /= 0) go to 999
      end if

      lo_idx = 2
      hi_idx = mx
    8 if (hi_idx-lo_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.opal_eos%x_grid(mid_idx)+1.0d-7) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 8
      end if
      x_bisect_idx = lo_idx

      lo_idx = nt
      hi_idx = 2
   11 if (lo_idx-hi_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (t6_temperature.eq.opal_eos%t6_list(1,mid_idx)) then
            lo_idx = mid_idx
            go to 14
         end if
         if (t6_temperature.le.opal_eos%t6_list(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 11
      end if
   14 t6_bisect_idx = lo_idx

      pressure_max = opal_eos%eos_table(x_bisect_idx,1,t6_bisect_idx, &
           density_index_edge(t6_bisect_idx))*t6_temperature* &
           opal_eos%density_grid(density_index_edge(t6_bisect_idx)) + &
           rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      pressure_min = opal_eos%eos_table(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           opal_eos%density_grid(1) + rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      if ((pressure_no_rad.gt.1.25d0*pressure_max) .or. &
           (pressure_no_rad.lt.pressure_min)) then
!      WRITE(ISHORT,'(" THE REQUESTED PRESSURE-TEMPERATURE NOT IN ",
!     *       "TABLE")')
!     STOP
!      WRITE(ISHORT,'("PNR, PMAX,PMIN=",3E14.4)') PNR,PMAX,PMIN
         go to 999
      end if

      density_trial1 = opal_eos%density_grid(density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac(hydrogen_fraction, t6_temperature, density_trial1, 1, &
           rad_flag, ierr, *999)
      if (ierr /= 0) go to 999
      pressure_trial1 = opal_eos%eos_output(1)
      if (pressure_trial1.gt.pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         density_trial1 = 0.2d0*density_trial1
         if (density_trial1.lt.1.0d-14) density_trial1 = 1.0d-14
         call esac(hydrogen_fraction, t6_temperature, density_trial1, 1, &
              rad_flag, ierr, *999)
         if (ierr /= 0) go to 999
         pressure_trial1 = opal_eos%eos_output(1)
      else
         density_trial2 = 5.0d0*density_trial1
!          IF(RHOG2 .GT. RHO(KLO)) RHOG2=RHO(KLO) ! Corrected below  llp 8/19/08
         if (density_trial2.gt.opal_eos%density_grid(density_index_edge(t6_bisect_idx))) &
              density_trial2 = opal_eos%density_grid(density_index_edge(t6_bisect_idx)) ! Had wrong pointer, see RHOG1= ten lines up
         call esac(hydrogen_fraction, t6_temperature, density_trial2, 1, &
              rad_flag, ierr, *999)
         if (ierr /= 0) go to 999
         pressure_trial2 = opal_eos%eos_output(1)
      end if

      refine_count = 0
    1 continue
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
      call esac(hydrogen_fraction, t6_temperature, density_trial3, 1, &
           rad_flag, ierr, *999)
      if (ierr /= 0) go to 999
      pressure_trial3 = opal_eos%eos_output(1)
!      IF (ABS((P3-PNR)/PNR) .LT. 1.D-5) THEN
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad).lt.0.5d-7) then
         rhoofp = density_trial3
!      WRITE(ISHORT,*)RHOG3,P,PNR,P3
!      WRITE(ISHORT,*)X,ZTAB,T6,P,RHOG3,IORDER,IRAD
         return
      end if
      if (pressure_trial3.gt.pressure_no_rad) then
         density_trial2 = density_trial3
         pressure_trial2 = pressure_trial3
         if (refine_count.lt.11) go to 1
!        WRITE(ISHORT,'("NO CONVERGENCE AFTER 10 TRIES")')
         go to 999
!        STOP
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
         if (refine_count.lt.11) go to 1
!        WRITE(ISHORT,'("NO CONVERGENCE AFTER 10 TRIES")')
         go to 999
!        STOP
      end if
  999 continue
      rhoofp = -999.0d0
!      WRITE(ISHORT,'("FAIL TO FIND RHO")')
      return
end function rhoofp
