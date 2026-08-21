!----------------------------------------------------------------------
! locate
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original locate.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Given an array xx of length n, and given a value x, returns a
! value j such that x is between xx(j) and xx(j+1). j=0 or n then
! out of range.
module locate_mod
contains
subroutine locate(xx, n, x, j)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: xx(n), x
      integer, intent(out) :: j

      integer :: jl, ju, jm
      save

      jl = 0
      ju = n+1
   10 if (ju-jl.gt.1) then
         jm = (ju+jl)/2
         if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
         goto 10
      end if
      j = jl
      if ((j .eq. 0) .and. (x .gt. 0.99d0*xx(1))) then
         j = 1
      end if
      if ((j .eq. n) .and. (x .lt. 1.01d0*xx(n))) then
         j = n-1
      end if
      return
end subroutine locate
end module locate_mod
