!----------------------------------------------------------------------
! splintd2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splintd2.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! double precision version. Taken from Numerical Recipes, Press, et
! al, p89. Given the arrays xa and ya of length n, which tabulate a
! function (with the xa(i)'s in order), and given the array y2a, which
! is the output of cspline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y.
!
! Note: xa/ya/y2a are dimensioned to the json=5000 module-wide
! maximum rather than to n, exactly as in the original file.
subroutine splintd2(xa, ya, n, y2a, x, y, klo, khi)
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: xa(json), ya(json), y2a(json), x
      double precision, intent(out) :: y
      integer, intent(out) :: klo, khi

! common/luout/: only short_file_unit (the .short log unit) is used
! here. Naming matches getopac.f90.
      integer :: ilast, idebug, itrack, short_file_unit, imilne, imodpt, &
           istor, iowr
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr

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
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0
      return
end subroutine splintd2
