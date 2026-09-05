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
!
! Readability W3 (2026): the row search and ramp calculation are
! eqbound_core.f90, shared with eqbound01; this wrapper keeps the
! 1995 density pre-check (log10 rho > 5.0d0 is outside the table).
subroutine eqbound(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)

      use opal_eos_lib
      implicit none

      integer, parameter :: nt = n_eos95_nt, nr = n_eos95_nr

      double precision, intent(in) :: temperature, log10_density
      double precision, intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
! --- locals ---
      double precision :: t6

      integer, intent(out) :: ierr

      ierr = 0

      t6 = temperature*1.0d-6
!     exit if outside table in rho
      if (log10_density.lt.-14d0 .or. log10_density.gt.5.0d0) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if
      call eqbound_core(nt, nr, opal95%t6_grid, opal95%t_row_index, opal95%density_edge_at_t, &
           opal95_density_index_edge, opal95%density_grid, t6, log10_density, ramp_factor, &
           in_opal_table, needs_ramp, ierr)

      return
end subroutine eqbound
