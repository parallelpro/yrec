!----------------------------------------------------------------------
! rhoofp01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rhoofp01.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2001 EOS analogue of rhoofp.f90 (see there for the general
! description). Note the convergence tolerance below (0.5d-7) matches
! rhoofp.f90's 1995 version; rhoofp06.f90 later loosened this to
! 1.0d-5 (see its own header note).
double precision function rhoofp01(hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag)

      use luout_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 191

! common/lreadco/: shared (by COMMON block name) with esac.f90/
! esac01.f90/esac06.f90/rhoofp.f90/rhoofp06.f90 -- see esac.f90.
      integer :: table_loaded_flag
      common/lreadco/ table_loaded_flag

! common/aeos/: see esac01.f90 for the full description.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr)
      integer :: x_loop_index, x_index_lo
      double precision :: x_grid(mx)
      common/aeos/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_loop_index, &
           x_index_lo, x_grid

! common/beos/: not used in this file; placeholders (see esac01.f90).
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), t6_index_lo(nr)
      common/beos/ z_table, eos_index_inverse, eos_var_order, t6_index_lo

! common/eeos/: only eos_output(1) is used here.
      double precision :: esact, eos_output(mv)
      common/eeos/ esact, eos_output


! density_index_edge(t6_idx): highest valid density-grid index for
! temperature-grid row t6_idx (local copy; see rhoofp.f90's note on
! why this isn't shared storage with readcoeos01.f90).
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

      save

      rat = rad_const_over_c
      radiation_pressure = 0.0d0
!      IF(IRAD .EQ. 1) PR=4.D0/3.D0*RAT*T6**4   ! MB
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (table_loaded_flag.ne.12345678) then
         hydrogen_fraction_dbg = 0.5d0
         t6_dbg = 1.0d0
         density_dbg = 0.001d0
         ideriv_dbg = 1
         call esac01(hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              rad_flag, *999)
      end if

      lo_idx = 2
      hi_idx = mx
    8 if (hi_idx-lo_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.x_grid(mid_idx)+1.0d-7) then
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
         if (t6_temperature.eq.t6_list(1,mid_idx)) then
            lo_idx = mid_idx
            go to 14
         end if
         if (t6_temperature.le.t6_list(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 11
      end if
   14 t6_bisect_idx = lo_idx

      pressure_max = eos_table(x_bisect_idx,1,t6_bisect_idx, &
           density_index_edge(t6_bisect_idx))*t6_temperature* &
           density_grid(density_index_edge(t6_bisect_idx)) + &
           rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      pressure_min = eos_table(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           density_grid(1) + rad_flag*4.0d0/3.0d0*rat*t6_temperature**4
      if ((pressure_no_rad.gt.1.25d0*pressure_max) .or. &
           (pressure_no_rad.lt.pressure_min)) then
!      write (ISHORT,'(" The requested pressure-temperature not in ",
!     X       "table")')
!     stop
!      write (ISHORT,'("pnr, pmax,pmin=",3e14.4)') pnr,pmax,pmin
         go to 999
      end if

      density_trial1 = density_grid(density_index_edge(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac01(hydrogen_fraction, t6_temperature, density_trial1, 1, &
           rad_flag, *999)
      pressure_trial1 = eos_output(1)
      if (pressure_trial1.gt.pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         density_trial1 = 0.2d0*density_trial1
         if (density_trial1.lt.1.0d-14) density_trial1 = 1.0d-14
         call esac01(hydrogen_fraction, t6_temperature, density_trial1, 1, &
              rad_flag, *999)
         pressure_trial1 = eos_output(1)
      else
         density_trial2 = 5.0d0*density_trial1
!          if(rhog2 .gt. rho(klo)) rhog2=rho(klo)  ! Corrected below llp 8/19/08
         if (density_trial2.gt.density_grid(density_index_edge(t6_bisect_idx))) &
              density_trial2 = density_grid(density_index_edge(t6_bisect_idx)) ! Had wrong pointer, see rhog1= ten lines up
         call esac01(hydrogen_fraction, t6_temperature, density_trial2, 1, &
              rad_flag, *999)
         pressure_trial2 = eos_output(1)
      end if

      refine_count = 0
    1 continue
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)
      call esac01(hydrogen_fraction, t6_temperature, density_trial3, 1, &
           rad_flag, *999)
      pressure_trial3 = eos_output(1)
!      if (abs((p3-pnr)/pnr) .lt. 1.D-5) then
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad).lt.0.5d-7) then
         rhoofp01 = density_trial3

         return
      end if
      if (pressure_trial3.gt.pressure_no_rad) then
         density_trial2 = density_trial3
         pressure_trial2 = pressure_trial3
         if (refine_count.lt.11) go to 1
!        write (ISHORT,'("No convergence after 10 tries")')
         go to 999
!        stop
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
         if (refine_count.lt.11) go to 1
!        write (*,'("No convergence after 10 tries")')
         go to 999
!        stop
      end if

  999 continue
      rhoofp01 = -999.0d0
!      WRITE(ISHORT,'("FAIL TO FIND RHO")')
      return

end function rhoofp01
