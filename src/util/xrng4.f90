!----------------------------------------------------------------------
! xrng4
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original xrng4.f; only variable names, source form, and comment
! style were updated.
!
! Get four points in the x array which surround x.
subroutine xrng4(grid_index, grid_size, window_start, window_end)

      implicit none
      integer, intent(in) :: grid_index, grid_size
      integer, intent(out) :: window_start, window_end

      save

      if (grid_index.le.2) then
         window_start = 1
         window_end = 4
      else if (grid_index.ge.grid_size-2) then
         window_start = grid_size-3
         window_end = grid_size
      else
         window_start = grid_index - 1
         window_end = grid_index + 2
      end if
      return
end subroutine xrng4
