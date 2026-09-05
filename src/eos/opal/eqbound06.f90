!----------------------------------------------------------------------
! eqbound06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqbound06.f; only variable names, source form, and comment
! style were updated.
!
! LLP 10/16/06: checks whether a (T, log10 rho) point falls within the
! OPAL 2006 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (see eqstat2.f90's use_opal2006_eos branch, which calls this
! immediately after oeqos06). Unlike eqbound/eqbound01 (1995/2001),
! this version tracks the table's ragged edge via per-row/per-column
! index arrays (opal06_t6_index_lo/opal06_density_index_edge) rather than a single
! linear-interpolated edge value.
subroutine eqbound06(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp)

      use opal_eos_lib
      use math_lib
      implicit none

      integer, parameter :: nr = n_eos06_nr, nt = n_eos06_nt

      double precision, intent(in) :: temperature, log10_density
      double precision, intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp

      integer :: t6_row, density_row
! --- locals ---
      double precision :: t6, density, density_ramp_factor, &
           t6_ramp_factor

      t6 = temperature*1.0d-6
      density = pow(10d0, log10_density)

!     Exit if outside table in rho
      if ((density.lt.opal06%density_grid(1)) .or. (density.ge.opal06%density_grid(nr))) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if

!     Exit if outside table in T6
      if ((t6.gt.opal06%t6_grid(1)) .or. (t6.le.opal06%t6_grid(nt))) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if

!     Initialize
      density_ramp_factor = 1d0   ! Initialize density ramp factor to 1
      t6_ramp_factor = 1d0        ! Initialize temperature ramp factor to 1
      in_opal_table = .true.      ! Presume that T6,D point is in OPAL table
      needs_ramp = .false.        ! presume ramping is not needed,
                                   ! i.e. we are not in the border of the table.

!     Check for easy borders
      if (t6.gt.opal06%t6_grid(2)) then
          t6_ramp_factor = (t6-opal06%t6_grid(1))/(opal06%t6_grid(2)-opal06%t6_grid(1)) ! Ramp in temperature
          needs_ramp = .true.
      end if
      if (density.lt.opal06%density_grid(2)) then
          density_ramp_factor = (density - opal06%density_grid(1))/ &
               (opal06%density_grid(2)-opal06%density_grid(1)) ! Ramp in density
          needs_ramp = .true.
      end if

!     Find a t6_row such that opal06%t6_grid(t6_row-1) >= T6 > opal06%t6_grid(t6_row)
!     and a density_row such that opal06%density_grid(density_row-1) < D <=
!     opal06%density_grid(density_row)
!     The ESAC06 variables (k=)opal06%t6_index_3 and (l-)opal06%density_index_2 are
!     close to what we need. A linear search will work fine.
      t6_row = opal06%t6_index_3
   do while (t6.le.opal06%t6_grid(t6_row))
        t6_row = t6_row+1
   end do
   do while (t6.gt.opal06%t6_grid(t6_row-1))
        t6_row = t6_row-1
   end do
!     We now have: opal06%t6_grid(t6_row-1) >= T6 > opal06%t6_grid(t6_row)

      density_row = opal06%density_index_2
   do while (density.gt.opal06%density_grid(density_row))
        density_row = density_row+1
   end do
   do while (density.le.opal06%density_grid(density_row-1))
        density_row = density_row-1
   end do
!     We now have: opal06%density_grid(density_row-1) < D <= opal06%density_grid(density_row)

!     For a given temperature in array opal06%t6_grid with index t6_row, e.g.,
!     opal06%t6_grid(t6_row), element opal06_density_index_edge(t6_row) of array
!     opal06_density_index_edge contains the index to the max allowed density
!     in array opal06%density_grid. So, for a given opal06%t6_grid(t6_row), the
!     associated D must be less than opal06%density_grid(opal06_density_index_edge(t6_row)).
!     A valid index in array opal06%density_grid must be less than or equal
!     to opal06_density_index_edge(t6_row)

      if (density_row.gt.opal06_density_index_edge(t6_row)) then
        in_opal_table = .false.
        needs_ramp = .true.
        ramp_factor = 0d0
        return
      end if
      if (density_row.eq.opal06_density_index_edge(t6_row)) then
        needs_ramp = .true.
        density_ramp_factor = (opal06%density_grid(density_row)-density)/ &
             (opal06%density_grid(density_row)-opal06%density_grid(density_row-1))
      end if

!     For a given density in array opal06%density_grid with index density_row,
!     i.e., opal06%density_grid(density_row), element opal06_t6_index_lo(density_row)
!     of array opal06_t6_index_lo contains the index to the min allowed
!     temperature in array opal06%t6_grid. So for a given
!     opal06%density_grid(density_row), the associated T6 must be greater
!     than opal06%t6_grid(opal06_t6_index_lo(density_row)). (Note that opal06%t6_grid is a
!     decreasing array.) A valid index in opal06%t6_grid must be less than or
!     equal to opal06_t6_index_lo(density_row).

      if (t6_row.gt.opal06_t6_index_lo(density_row)) then
        in_opal_table = .false.
        needs_ramp = .true.
        ramp_factor = 0d0
        return
      end if
      if (t6_row.eq.opal06_t6_index_lo(density_row)) then
        needs_ramp = .true.
        t6_ramp_factor = (t6-opal06%t6_grid(t6_row))/(opal06%t6_grid(t6_row-1)-opal06%t6_grid(t6_row))
      end if

      ramp_factor = density_ramp_factor * t6_ramp_factor
      return
end subroutine eqbound06
