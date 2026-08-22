!----------------------------------------------------------------------
! eqbound
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqbound.f; only variable names, source form, and comment
! style were updated.
!
! MHP 8/98: checks whether a (T, log10 rho) point falls within the
! OPAL 1995 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (see eqstat2.f90's use_opal95_eos branch, which calls this
! immediately after oeqos).
subroutine eqbound(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)

      use opal_eos_lib
      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 77, nt = 56

      double precision, intent(in) :: temperature, log10_density
      double precision, intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
! --- locals ---
      double precision :: t6, table_edge_density, ramp_start_density
      double precision :: t6_top_of_table, t_fraction
      integer :: t6_scan_idx

      integer, intent(out) :: ierr

      ierr = 0

      t6 = temperature*1.0d-6
!     exit if outside table in rho
      if (log10_density.lt.-14d0 .or. log10_density.gt.5.0d0) then
         continue
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         
         return
      end if
!     find nearest table element in t.
      if (t6.lt.opal_eos%t6_grid(opal_eos%t_row_index)) then
         do t6_scan_idx = opal_eos%t_row_index+1, nt
            if (t6.ge.opal_eos%t6_grid(t6_scan_idx)) then
               opal_eos%t_row_index = t6_scan_idx - 1
               goto 10
            end if
         end do
         t6_top_of_table = opal_eos%t6_grid(nt)
!        sr call should have been stopped outside table bounds; stop code
         write(*,5) t6, t6_top_of_table, opal_eos%t_row_index
    5    format(' ERROR IN OPAL EOS: OUTSIDE TABLE IN T6',2F10.6,I5)
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
   10    continue
      else
         do t6_scan_idx = opal_eos%t_row_index, 1, -1
            if (t6.le.opal_eos%t6_grid(t6_scan_idx)) then
               opal_eos%t_row_index = t6_scan_idx
               goto 20
            end if
         end do
         t6_top_of_table = opal_eos%t6_grid(1)
!        sr call should have been stopped outside table bounds; stop code
         write(*,5) t6, t6_top_of_table, opal_eos%t_row_index
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
   20    continue
      end if
      t_fraction = (t6 - opal_eos%t6_grid(opal_eos%t_row_index+1))/ &
           (opal_eos%t6_grid(opal_eos%t_row_index)-opal_eos%t6_grid(opal_eos%t_row_index+1))
!     define table edge in rho by linear interpolation.
      table_edge_density = opal_eos%density_edge_at_t(opal_eos%t_row_index+1)
      table_edge_density = log10(table_edge_density)
!     define beginning of ramp in the same fashion -
!     ramp is defined as one table element wide.
      ramp_start_density = opal_eos%density_grid(opal_eos%density_index_edge_at_t(opal_eos%t_row_index+1)-1)
      ramp_start_density = log10(ramp_start_density)
!     check if within table bounds in rho
      if (log10_density.gt.table_edge_density) then
         continue
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         
         return
      end if

!     If we get here, the point is in the table.
      in_opal_table = .true.

!     Now we check if ramping is needed.
!     First we check if ramping in temperature is needed.
      if (t6.le.opal_eos%t6_grid(nt-1)) then
         needs_ramp = .true.
         ramp_factor = (t6-opal_eos%t6_grid(nt))/(opal_eos%t6_grid(nt-1)-opal_eos%t6_grid(nt))
      else if (t6.ge.opal_eos%t6_grid(2)) then
         needs_ramp = .true.
         ramp_factor = (opal_eos%t6_grid(1)-t6)/(opal_eos%t6_grid(1)-opal_eos%t6_grid(2))
      else if (log10_density.ge.ramp_start_density) then
!        If we get here, ramping in density is needed.
         needs_ramp = .true.
         ramp_factor = (table_edge_density-log10_density)/ &
              (table_edge_density-ramp_start_density)
      else
!        If we get here, we are in the middle of the table and
!        no ramping is needed.
         needs_ramp = .false.
      end if

      return

!     Error Exit.
 9999 continue
      in_opal_table = .false.       ! Not in table
      needs_ramp = .true.           ! Turn on ramping
      ramp_factor = 0d0             ! Set ramping factor to zero
!                            This way, out of table results are ramped to zero.
      return
end subroutine eqbound
