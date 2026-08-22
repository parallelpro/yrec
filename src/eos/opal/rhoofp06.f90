!----------------------------------------------------------------------
! rhoofp06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rhoofp06.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2006 EOS analogue of rhoofp01.f90 (see there and rhoofp.f90
! for the general description), with three preserved differences:
! (1) pressure_max/pressure_min here do NOT add a radiation term
! (rhoofp.f90/rhoofp01.f90 both add rad_flag*4/3*rat*t6**4); (2) the
! trial esac06.f90 calls below pass a hardcoded 0 for the radiation
! flag (not rad_flag) -- only the very first, table-priming call
! passes rad_flag through; (3) the convergence tolerance was loosened
! from 0.5d-7 to 1.0d-5 (see the dated comment below) to stop
! opal_eos%eos_output_06(5) from growing without bound and crashing some model
! runs.
double precision function rhoofp06(hydrogen_fraction, t6_temperature, &
     pressure_e12, rad_flag, ierr)

      use opal_eos_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           pressure_e12
      integer, intent(in) :: rad_flag

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197





      double precision :: rad_const_over_c
! NOTE: no D-suffix in the original (rhoofp06.f) -- preserved
! verbatim: parsed as a single-precision constant and then widened to
! double precision, which is NOT bit-identical to the correctly-
! rounded double value.
      data rad_const_over_c/1.8914785e-3/
!--------------------------------------------------------------------
! --- locals ---
      double precision :: rat, radiation_pressure, pressure_no_rad
      double precision :: hydrogen_fraction_dbg, t6_dbg, density_dbg
      integer :: ideriv_dbg
      integer :: lo_idx, hi_idx, mid_idx
      integer :: x_bisect_idx, t6_bisect_idx
      double precision :: pressure_max, pressure_min
      double precision :: density_trial1, density_trial2, density_trial3
      double precision :: pressure_trial1, pressure_trial2, pressure_trial3
      integer :: refine_count

      integer, intent(out) :: ierr

      ierr = 0

      rat = rad_const_over_c
      radiation_pressure = 0.0d0
      if (rad_flag.eq.1) radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
      pressure_no_rad = pressure_e12 - radiation_pressure

      if (opal_eos%table_loaded_flag.ne.12345678) then
         hydrogen_fraction_dbg = 0.5d0
         t6_dbg = 1.0d0
         density_dbg = 0.001d0
         ideriv_dbg = 1
         call esac06(hydrogen_fraction_dbg, t6_dbg, density_dbg, ideriv_dbg, &
              rad_flag, ierr, *999)
         if (ierr /= 0) go to 999
      end if

      lo_idx = 2
      hi_idx = mx
    8 if (hi_idx-lo_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.opal_eos%x_grid_06(mid_idx)+1.0d-7) then
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
         if (t6_temperature.eq.opal_eos%t6_list_06(1,mid_idx)) then
            lo_idx = mid_idx
            go to 14
         end if
         if (t6_temperature.le.opal_eos%t6_list_06(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 11
      end if
   14 t6_bisect_idx = lo_idx

      pressure_max = opal_eos%eos_table_06(x_bisect_idx,1,t6_bisect_idx, &
           opal_eos%density_index_edge_06(t6_bisect_idx))*t6_temperature* &
           opal_eos%density_grid_06(opal_eos%density_index_edge_06(t6_bisect_idx))
      pressure_min = opal_eos%eos_table_06(x_bisect_idx,1,t6_bisect_idx,1)*t6_temperature* &
           opal_eos%density_grid_06(1)
      if ((pressure_no_rad.gt.1.25d0*pressure_max) .or. &
           (pressure_no_rad.lt.pressure_min)) then
!      write (ISHORT,'(" The requested pressure-temperature not in",
!     x   " the OPAL 2006 EOS table")')
!     stop
!      write (ISHORT,'("pnr, pmax,pmin=",3e14.4)') pnr,pmax,pmin
         go to 999     !RHOOFP06 error exit
      end if

      density_trial1 = opal_eos%density_grid_06(opal_eos%density_index_edge_06(t6_bisect_idx))* &
           pressure_no_rad/pressure_max
      call esac06(hydrogen_fraction, t6_temperature, density_trial1, 1, 0, ierr, *999)
      if (ierr /= 0) go to 999
      pressure_trial1 = opal_eos%eos_output_06(1)
      if (pressure_trial1.gt.pressure_no_rad) then
         pressure_trial2 = pressure_trial1
         density_trial2 = density_trial1
         density_trial1 = 0.2d0*density_trial1
         if (density_trial1.lt.1.0d-14) density_trial1 = 1.0d-14
         call esac06(hydrogen_fraction, t6_temperature, density_trial1, 1, 0, ierr, *999)
         if (ierr /= 0) go to 999
         pressure_trial1 = opal_eos%eos_output_06(1)
      else
         density_trial2 = 5.0d0*density_trial1
!          if(rhog2 .gt. rho(klo)) rhog2=rho(klo)  ! Corrected below   llp  8/19/08
         if (density_trial2.gt.opal_eos%density_grid_06(opal_eos%density_index_edge_06(t6_bisect_idx))) &
              density_trial2 = opal_eos%density_grid_06(opal_eos%density_index_edge_06(t6_bisect_idx)) ! Had wrong pointer, see rhog1= ten lines up
         call esac06(hydrogen_fraction, t6_temperature, density_trial2, 1, 0, ierr, *999)
         if (ierr /= 0) go to 999
         pressure_trial2 = opal_eos%eos_output_06(1)
      end if

      refine_count = 0
    1 continue
      refine_count = refine_count + 1
      density_trial3 = density_trial1 + (density_trial2-density_trial1)* &
           (pressure_no_rad-pressure_trial1)/(pressure_trial2-pressure_trial1)  ! KC 2025-05-31
      call esac06(hydrogen_fraction, t6_temperature, density_trial3, 1, 0, ierr, *999)
      if (ierr /= 0) go to 999
      pressure_trial3 = opal_eos%eos_output_06(1)
! Changed the comparison below to use the commented-out value 1.D-5
! found here to prevent array value eos(5) from growing without bound
! and crashing the program during certain model runs. - MR 2025-10-10
      if (abs((pressure_trial3-pressure_no_rad)/pressure_no_rad).lt.1.0d-5) then
!      IF (DABS((P3-PNR)/PNR) .LT. 0.5D-7) THEN
         rhoofp06 = density_trial3
         return
      end if
      if (pressure_trial3.gt.pressure_no_rad) then
         density_trial2 = density_trial3
         pressure_trial2 = pressure_trial3
         if (refine_count.lt.11) go to 1
!        write (ISHORT,'("Rhoofp06: No convergence after 10 tries")')
         go to 999     !RHOOFP06 error exit
!        stop
      else
         density_trial1 = density_trial3
         pressure_trial1 = pressure_trial3
         if (refine_count.lt.11) go to 1
!        write (ISHORT,'("RHOOFP06: No convergence after 10 tries")')
         go to 999   ! RHOOFP06 error exit
!        stop
      end if

  999 continue
      rhoofp06 = -999.0d0
!      WRITE(ISHORT,'("RHOOFP06: FAILED TO FIND RHO")')
      return


end function rhoofp06
