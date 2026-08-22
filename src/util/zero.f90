!----------------------------------------------------------------------
! zero
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original zero.f; only variable names, source form, and comment
! style were updated.
!
! Sets array(n) = 0, n = 1, array_size.
subroutine zero(array, array_size)

      implicit none
      integer, intent(in) :: array_size
      double precision, intent(out) :: array(array_size)

      integer :: index
      do index=1,array_size
         array(index)=0.d0
    1 continue
      end do
      return
end subroutine zero
