!----------------------------------------------------------------------
! findex
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original findex.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Interpolation package for Cox's and for Kurucz's opacities.
! YCK 3/91
!
! Finds index such that grid_x(index) <= x_eval < grid_x(index+1),
! searching outward from the incoming value of index as an initial
! guess. On entry, index is reset to 1 if it is out of [1,n_grid]. If
! x_eval falls below grid_x(1), index is returned as -1; if it falls
! at or above grid_x(n_grid), index is returned as -n_grid.
module findex_mod
contains
subroutine findex(grid_x, n_grid, x_eval, index)
      implicit none
      integer, intent(in) :: n_grid
      double precision, intent(in) :: grid_x(n_grid)
      double precision, intent(in) :: x_eval
      integer, intent(inout) :: index

      integer :: found_index, j
      save

! find the 'index'
      if(index.lt.1.or.index.gt.n_grid)index=1
      found_index=index
      if(x_eval.lt.grid_x(found_index))then
         do 211 j=found_index-1,1,-1
            if(grid_x(j).le.x_eval)then
               found_index=j
               goto 213
            endif
 211     continue
         found_index=-1
      else
         do 212 j=found_index,n_grid-1
            if(grid_x(j+1).gt.x_eval)then
               found_index=j
               goto 213
            endif
 212     continue
         found_index = -n_grid
      endif
 213  index=found_index

      return
end subroutine findex
end module findex_mod
