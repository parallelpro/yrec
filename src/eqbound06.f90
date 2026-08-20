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
! index arrays (t6_index_lo/density_index_edge) rather than a single
! linear-interpolated edge value.
subroutine eqbound06(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp)

      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197

      double precision, intent(in) :: temperature, log10_density
      double precision, intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp

      integer :: t6_row, density_row

! common/aeos06/: the main OPAL 2006 EOS table and its interpolation
! grids/scratch arrays. Only density_grid and t6_grid are used here;
! the rest are placeholders. Naming matches esac06.f90.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr)
      integer :: x_loop_index, x_index_lo
      double precision :: x_grid(mx)
      common/aeos06/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_loop_index, &
           x_index_lo, x_grid

! common/beos06/: t6_index_lo(density_row) is the index in t6_grid of
! the lowest available temperature at that density row;
! density_index_edge(t6_row) is the index in density_grid of the
! highest available density at that temperature row. Both are used
! here to detect the table's ragged edge. The rest are placeholders.
! Naming matches esac06.f90.
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), &
           t6_index_lo(nr), density_index_edge(nt)
      common/beos06/ z_table, eos_index_inverse, eos_var_order, &
           t6_index_lo, density_index_edge

! common/bbeos06/: density_index_3/t6_index_3 seed the linear search
! below for the table cell containing the point (carried over from
! the ESAC06 interpolation as a good starting guess); the rest are
! unused placeholders here. Naming matches esac06.f90.
      integer :: density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order
      common/bbeos06/ density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order

      save

! --- locals ---
      double precision :: t6, density, density_ramp_factor, &
           t6_ramp_factor

      t6 = temperature*1.0d-6
      density = 10d0**log10_density

!     Exit if outside table in rho
      if ((density.lt.density_grid(1)) .or. (density.ge.density_grid(nr))) then
         goto 9999        ! Out of Table in density. Go to out of table exit.
      end if

!     Exit if outside table in T6
      if ((t6.gt.t6_grid(1)) .or. (t6.le.t6_grid(nt))) then
         goto 9999        ! Out of Table in temperature. Go to out of table exit.
      end if

!     Initialize
      density_ramp_factor = 1d0   ! Initialize density ramp factor to 1
      t6_ramp_factor = 1d0        ! Initialize temperature ramp factor to 1
      in_opal_table = .true.      ! Presume that T6,D point is in OPAL table
      needs_ramp = .false.        ! presume ramping is not needed,
                                   ! i.e. we are not in the border of the table.

!     Check for easy borders
      if (t6.gt.t6_grid(2)) then
          t6_ramp_factor = (t6-t6_grid(1))/(t6_grid(2)-t6_grid(1)) ! Ramp in temperature
          needs_ramp = .true.
      end if
      if (density.lt.density_grid(2)) then
          density_ramp_factor = (density - density_grid(1))/ &
               (density_grid(2)-density_grid(1)) ! Ramp in density
          needs_ramp = .true.
      end if

!     Find a t6_row such that t6_grid(t6_row-1) >= T6 > t6_grid(t6_row)
!     and a density_row such that density_grid(density_row-1) < D <=
!     density_grid(density_row)
!     The ESAC06 variables (k=)t6_index_3 and (l-)density_index_2 are
!     close to what we need. A linear search will work fine.
      t6_row = t6_index_3
   10 if (t6.le.t6_grid(t6_row)) then
        t6_row = t6_row+1
        goto 10
      end if
   20 if (t6.gt.t6_grid(t6_row-1)) then
        t6_row = t6_row-1
        goto 20
      end if
!     We now have: t6_grid(t6_row-1) >= T6 > t6_grid(t6_row)

      density_row = density_index_2
   30 if (density.gt.density_grid(density_row)) then
        density_row = density_row+1
        goto 30
      end if
   40 if (density.le.density_grid(density_row-1)) then
        density_row = density_row-1
        goto 40
      end if
!     We now have: density_grid(density_row-1) < D <= density_grid(density_row)

!     For a given temperature in array t6_grid with index t6_row, e.g.,
!     t6_grid(t6_row), element density_index_edge(t6_row) of array
!     density_index_edge contains the index to the max allowed density
!     in array density_grid. So, for a given t6_grid(t6_row), the
!     associated D must be less than density_grid(density_index_edge(t6_row)).
!     A valid index in array density_grid must be less than or equal
!     to density_index_edge(t6_row)

      if (density_row.gt.density_index_edge(t6_row)) then
        goto 9999       ! Out of table exit in density.
      end if
      if (density_row.eq.density_index_edge(t6_row)) then
        needs_ramp = .true.
        density_ramp_factor = (density_grid(density_row)-density)/ &
             (density_grid(density_row)-density_grid(density_row-1))
      end if

!     For a given density in array density_grid with index density_row,
!     i.e., density_grid(density_row), element t6_index_lo(density_row)
!     of array t6_index_lo contains the index to the min allowed
!     temperature in array t6_grid. So for a given
!     density_grid(density_row), the associated T6 must be greater
!     than t6_grid(t6_index_lo(density_row)). (Note that t6_grid is a
!     decreasing array.) A valid index in t6_grid must be less than or
!     equal to t6_index_lo(density_row).

      if (t6_row.gt.t6_index_lo(density_row)) then
        goto 9999      ! Out of table exit in temperature
      end if
      if (t6_row.eq.t6_index_lo(density_row)) then
        needs_ramp = .true.
        t6_ramp_factor = (t6-t6_grid(t6_row))/(t6_grid(t6_row-1)-t6_grid(t6_row))
      end if

      ramp_factor = density_ramp_factor * t6_ramp_factor
      return

!     OUT OF TABLE EXIT
 9999 continue
      in_opal_table = .false.       ! Not in table
      needs_ramp = .true.           ! Turn on ramping
      ramp_factor = 0d0             ! Set ramping factor to zero
!                            This way, out of table results are ramped to zero.
      return

end subroutine eqbound06
