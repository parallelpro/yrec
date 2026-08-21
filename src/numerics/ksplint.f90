!----------------------------------------------------------------------
! ksplint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ksplint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Fixed-size (4-point) cubic-spline evaluation, companion to kspline:
! given the table xa/ya and the second derivatives y2a from kspline,
! evaluate the spline at x.
module ksplint_mod
contains
subroutine ksplint(xa, ya, y2a, x, y)
      implicit none
      integer, parameter :: nm = 4

      double precision, intent(in) :: xa(nm), ya(nm), y2a(nm), x
      double precision, intent(out) :: y

      double precision :: h, a, b
      integer :: klo, khi, k
      save

      klo = 1
      khi = nm
    1 continue
      if (khi-klo .gt. 1) then
         k = (khi+klo)/2
         if (xa(k) .gt. x) then
            khi = k
         else
            klo = k
         end if
         goto 1
      end if
!      write(*,*) khi, klo, xa(khi), xa(klo), x
      h = xa(khi) - xa(klo)
      if (h .eq. 0d0) then
            print*, 'Ksplint failure'
            stop
      end if
!      if (h .eq. 0d0) stop 911
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+ &
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
      return
end subroutine ksplint
end module ksplint_mod
