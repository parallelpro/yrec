!----------------------------------------------------------------------
! splint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! double precision version for opacities, changed from real*4 10/02
! MHP. Taken from Numerical Recipes, Press, et al, p89. Given the
! arrays xa and ya of length n, which tabulate a function (with the
! xa(i)'s in order), and given the array y2a, which is the output of
! cspline above, and given a value of x, this routine returns a
! cubic-spline interpolated value y.
subroutine splint(xa, ya, n, y2a, x, y, klo, khi)
      use luout_lib
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: xa(n), ya(n), y2a(n), x
      double precision, intent(out) :: y
      integer, intent(out) :: klo, khi


      integer :: k
      double precision :: h, a, b
      save

      klo = 1
      khi = n
    1 if (khi-klo .gt. 1) then
         k = (khi+klo)/2
         if (xa(k) .gt. x) then
            khi = k
         else
            klo = k
         end if
         goto 1
      end if
      h = xa(khi) - xa(klo)
      if (h .eq. 0d0) then
           write(short_file_unit,*) 'ERROR IN SPLINT ROUTINE.'
         stop
      end if
      a = (xa(khi)-x)/h
      b = (x - xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+ &
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
      return
end subroutine splint
