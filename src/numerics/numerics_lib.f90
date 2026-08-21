!----------------------------------------------------------------------
! numerics_lib
!----------------------------------------------------------------------
! Aggregates YREC's generic, COMMON-free numerical utility routines
! (splines, interpolation, linear algebra, integration) into a single
! module, following MESA's convention of one small public *_lib module
! per domain. Each routine below was originally its own file
! (numerics/<name>.f90, with per-file conversion notes preserved below);
! bodies are unchanged from that conversion, only the per-routine module
! wrapper was removed and the routines concatenated under one module.
! Callers use this via `use numerics_lib`.
!
! lir and ratext (below, at the end of this module) were moved in
! (2026) from nuclear_lib.f90, where they'd originally been filed
! (along with the rest of nuclear/) purely because nuclear/liburn.f90
! is one of their callers -- both are generic numerics with no
! nuclear-physics content, and this module's own bsstep/intpt already
! called them via an errant `use nuclear_lib` (ratext is in fact the
! textbook Numerical-Recipes companion to bsstep/mmid above). See
! GUIDELINES.md's rule that folder/module placement should track
! function, not caller.
module numerics_lib
contains

!----------------------------------------------------------------------
! boole
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original boole.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Performs 5-point Newton-Cotes (Boole's rule) integration of
! tabulated data.
!      x: vector of independent variable values (x)
!      y: vector of dependent variable values (f(x))
!      n: number of elements in x
!      n_grid: number of elements desired in the interpolated grid
!              (should be 1+4*integer)
! integral: returned, integrated value of the function (kept as a
!           length-1 array, exactly as in the original file)
!
! Data need to be evenly gridded for the integration formula to
! work. Data that are not already so are interpolated onto an even
! grid using spline interpolation (Numerical Recipes SPLINE, renamed
! splinj here, and SPLINT).
subroutine boole(x, y, n, n_grid, integral)
      implicit none
      double precision, parameter :: scalex = 1e-11
      double precision, parameter :: scaley = 1e7

      double precision, intent(in) :: x(n), y(n)
      integer, intent(in) :: n, n_grid
      double precision, intent(out) :: integral(1)

      double precision :: h_step(1), y2_deriv(n)
      double precision :: x_even(n_grid), y_even(n_grid)
      double precision :: x_scaled(n), y_scaled(n)
      integer :: i, num_quads, klo, khi

! rescale radius and cs vectors to have values ~1
      do i = 1, n
            x_scaled(i) = x(i)*scalex
            y_scaled(i) = y(i)*scaley
      end do
! deal with unevenly gridded datasets:
            call splinj(x_scaled, y_scaled, y2_deriv, n) ! get derivs of interp. fn.
            do i = 1, n_grid
                  x_even(i) = x_scaled(1)+(i-1)*(x_scaled(n)-x_scaled(1))/(n_grid-1)
                  call splint(x_scaled, y_scaled, n, y2_deriv, x_even(i), y_even(i), klo, khi)
            end do

! how many sets of four points do we have?
      num_quads = (n_grid-1)/4
      h_step = (x_scaled(n)-x_scaled(1))/(n_grid-1)
! for each set of 4 points before the last, apply formula and add up
      integral = 0.0d0
      do i = 0, num_quads-1
            integral = integral+2.0*h_step*(7*y_even(1+4*i)+ 32*y_even(2+4*i) + 12*y_even(3+4*i) + &
                  32*y_even(4+4*i) + 7*y_even(5+4*i))/45.0
      end do
! rescale result back to actual units:
      integral = integral/(scalex*scaley)
!      print*, 'Last point in vector:', n_grid
!      print*, 'Last point integrated:', 5+4*(num_quads-1)

!--------------------------------------------------------------
!                  open(unit=100,file='diagnostic3.out',status='old')
!                  do i=1,n
!                              write(100,1504) x(i), y(i), x_scaled(i), y_scaled(i),
!     *                         x_even(i), y_even(i), integral(1), h_step(1)
!                  end do
!1504                  format(1x,8e16.8)
!                  close(100)
!----------------------------------------------------------------

end subroutine boole

!----------------------------------------------------------------------
! cspline
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original cspline.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! double precision version for opacities. Taken from Numerical
! Recipes, Press et al, p88. Given arrays x and y of length n
! containing a tabulated function, i.e. y(i) = f(x(i)) with
! x(1) < x(2) < ... < x(n), and given values yp1 and ypn for the first
! derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2 of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i). If yp1 and/or ypn are equal to 1.0e30 or
! larger, the routine is signalled to set the corresponding boundary
! condition for a natural spline, with zero second derivative on that
! boundary.
subroutine cspline(x, y, n, yp1, ypn, y2)
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: x(n), y(n), yp1, ypn
      double precision, intent(out) :: y2(n)

      double precision :: u(json)
      integer :: i, k
      double precision :: sig, p, qn, un
      save

      if (yp1 .gt. 0.99d30) then
         y2(1) = 0.0d0
         u(1) = 0.0d0
      else
         y2(1) = -0.5d0
         u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      end if
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.0d0
         y2(i) = (sig-1.0d0)/p
         u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn .gt. 0.99d30) then
         qn = 0.0d0
         un = 0.0d0
      else
         qn = 0.5d0
         un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      end if
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k = n-1, 1, -1
         y2(k) = y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine cspline

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

!----------------------------------------------------------------------
! inter3
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original inter3.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! inter3 is the interpolation routine for density in the Livermore
! opacity tables, and it uses a 3-point Lagrangian interpolation
! scheme. Given the 3 abscissa points x_nodes and evaluation point
! x_eval, returns the Lagrangian weights (weight) and their
! derivatives with respect to x_eval (dweight) such that a function's
! interpolated value/derivative at x_eval are sum(weight(:)*f(:)) and
! sum(dweight(:)*f(:)) for f sampled at x_nodes.
subroutine inter3(x_nodes, weight, dweight, x_eval)

      implicit none
      double precision, intent(in) :: x_nodes(3)
      double precision, intent(out) :: weight(3), dweight(3)
      double precision, intent(in) :: x_eval

      double precision :: diff32, diff31, diff21
      double precision :: denom1, denom2, denom3
      double precision :: dx1, dx2, dx3
      save

! inter3 is the interpolation routine for density in the livermore
! opacity tables, and it uses a 3-point lagrangian interpolation scheme.
      diff32 = x_nodes(3) - x_nodes(2)
      diff31 = x_nodes(3) - x_nodes(1)
      diff21 = x_nodes(2) - x_nodes(1)
      denom1 = diff21*diff31
      denom2 = -diff21*diff32
      denom3 = diff31*diff32
      dx1 = x_eval - x_nodes(1)
      dx2 = x_eval - x_nodes(2)
      dx3 = x_eval - x_nodes(3)
      weight(1) = (dx2*dx3)/denom1
      weight(2) = (dx1*dx3)/denom2
      weight(3) = (dx1*dx2)/denom3
      dweight(1) = (dx2 + dx3)/denom1
      dweight(2) = (dx1 + dx3)/denom2
      dweight(3) = (dx1 + dx2)/denom3

      return
end subroutine inter3

!----------------------------------------------------------------------
! interp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original interp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! interp is the interpolation routine for the VandenBerg opacity
! tables (CAPPA), and it uses a 4-point Lagrangian interpolation
! scheme. Given the 4 abscissa points x_nodes and evaluation point
! x_eval, returns the Lagrangian weights (weight) and their
! derivatives with respect to x_eval (dweight) such that a function's
! interpolated value/derivative at x_eval are sum(weight(:)*f(:)) and
! sum(dweight(:)*f(:)) for f sampled at x_nodes.
subroutine interp(x_nodes, weight, dweight, x_eval)

      implicit none
      double precision, intent(in) :: x_nodes(4)
      double precision, intent(out) :: weight(4), dweight(4)
      double precision, intent(in) :: x_eval

      double precision :: diff43, diff42, diff41, diff32, diff31, diff21
      double precision :: denom1, denom2, denom3, denom4
      double precision :: dx1, dx2, dx3, dx4
      save

      diff43 = x_nodes(4) - x_nodes(3)
      diff42 = x_nodes(4) - x_nodes(2)
      diff41 = x_nodes(4) - x_nodes(1)
      diff32 = x_nodes(3) - x_nodes(2)
      diff31 = x_nodes(3) - x_nodes(1)
      diff21 = x_nodes(2) - x_nodes(1)
      denom1 = -diff21*diff31*diff41
      denom2 = diff21*diff32*diff42
      denom3 = -diff31*diff32*diff43
      denom4 = diff41*diff42*diff43
      dx1 = x_eval - x_nodes(1)
      dx2 = x_eval - x_nodes(2)
      dx3 = x_eval - x_nodes(3)
      dx4 = x_eval - x_nodes(4)
      weight(1) = (dx2*dx3*dx4)/denom1
      weight(2) = (dx1*dx3*dx4)/denom2
      weight(3) = (dx1*dx2*dx4)/denom3
      weight(4) = (dx1*dx2*dx3)/denom4
      dweight(1) = (dx3*dx4 + dx2*dx4 + dx2*dx3)/denom1
      dweight(2) = (dx3*dx4 + dx1*dx4 + dx1*dx3)/denom2
      dweight(3) = (dx2*dx4 + dx1*dx4 + dx1*dx2)/denom3
      dweight(4) = (dx2*dx3 + dx1*dx3 + dx1*dx2)/denom4

      return
end subroutine interp

!----------------------------------------------------------------------
! intrp2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intrp2.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! interp without derivatives. Interp is the interpolation routine for
! the VandenBerg opacity tables (CAPPA), and it uses a 4-point
! Lagrangian interpolation scheme. Given the 4 abscissa points
! x_nodes and evaluation point x_eval, returns the Lagrangian weights
! (weight) such that a function's interpolated value at x_eval is
! sum(weight(:)*f(:)) for f sampled at x_nodes.
subroutine intrp2(x_nodes, weight, x_eval)

      implicit none
      double precision, intent(in) :: x_nodes(4)
      double precision, intent(out) :: weight(4)
      double precision, intent(in) :: x_eval

      double precision :: diff43, diff42, diff41, diff32, diff31, diff21
      double precision :: denom1, denom2, denom3, denom4
      double precision :: dx1, dx2, dx3, dx4
      save

! interp is the interpolation routine for the VandenBerg
! opacity tables(CAPPA), and it uses a 4-point Lagrangian
! interpolation scheme.
      diff43 = x_nodes(4) - x_nodes(3)
      diff42 = x_nodes(4) - x_nodes(2)
      diff41 = x_nodes(4) - x_nodes(1)
      diff32 = x_nodes(3) - x_nodes(2)
      diff31 = x_nodes(3) - x_nodes(1)
      diff21 = x_nodes(2) - x_nodes(1)
      denom1 = -diff21*diff31*diff41
      denom2 = diff21*diff32*diff42
      denom3 = -diff31*diff32*diff43
      denom4 = diff41*diff42*diff43
      dx1 = x_eval - x_nodes(1)
      dx2 = x_eval - x_nodes(2)
      dx3 = x_eval - x_nodes(3)
      dx4 = x_eval - x_nodes(4)
      weight(1) = (dx2*dx3*dx4)/denom1
      weight(2) = (dx1*dx3*dx4)/denom2
      weight(3) = (dx1*dx2*dx4)/denom3
      weight(4) = (dx1*dx2*dx3)/denom4
      return
end subroutine intrp2

!----------------------------------------------------------------------
! kspline
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original kspline.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Fixed-size (4-point) natural cubic spline coefficient generator,
! same algorithm as cspline/splinj/splinc but hardwired to n=4 points
! and always the natural-spline boundary condition.
subroutine kspline(x, y, y2)
      implicit none
      integer, parameter :: nm = 4

      double precision, intent(in) :: x(nm), y(nm)
      double precision, intent(out) :: y2(nm)

      double precision :: u(nm), sig, qn, un, p
      integer :: n, i, k
      save

      n = nm
! natural spline
      y2(1) = 0.0d0
      u(1) = 0.0d0
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.0d0
         y2(i) = (sig-1.0d0)/p
         u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      qn = 0.0d0
      un = 0.0d0
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k = n-1, 1, -1
         y2(k) = y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine kspline

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

!----------------------------------------------------------------------
! lubksb
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lubksb.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! LU back-substitution (Numerical Recipes LUBKSB companion to
! LUDCMP): solves a*x=b for x, given the LU decomposition of a (in
! place, as produced by ludcmp) and its pivot record indx. b is
! overwritten with the solution.
subroutine lubksb(a, n, np, indx, b)
      implicit none
      integer, intent(in) :: n, np
      double precision, intent(in) :: a(np,np)
      integer, intent(in) :: indx(n)
      double precision, intent(inout) :: b(n)

      integer :: i, ii, j, ll
      double precision :: sum
      save

      ii = 0
      do i = 1, n
          ll = indx(i)
          sum = b(ll)
          b(ll) = b(i)
          if (ii.ne.0) then
              do j = ii, i - 1
                  sum = sum - a(i,j)*b(j)
              end do

          else if (sum.ne.0d0) then
              ii = i
          end if

          b(i) = sum
      end do
      do i = n, 1, -1
          sum = b(i)
          if (i.lt.n) then
              do j = i + 1, n
                  sum = sum - a(i,j)*b(j)
              end do
          end if

          b(i) = sum/a(i,i)
      end do
      return

end subroutine lubksb

!----------------------------------------------------------------------
! ludcmp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ludcmp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! LU decomposition with partial pivoting (Numerical Recipes LUDCMP):
! replaces a by the LU decomposition of a row-wise permutation of
! itself; indx records the row permutation, d returns +/-1 depending
! on whether the number of row interchanges was even or odd (used by
! the caller to get the sign of the determinant).
subroutine ludcmp(a, n, np, indx, d)
      implicit none
      integer, parameter :: nmax = 100
      double precision, parameter :: tiny = 1.0d-20

      integer, intent(in) :: n, np
      double precision, intent(inout) :: a(np,np)
      integer, intent(out) :: indx(n)
      double precision, intent(out) :: d

      integer :: i, imax, j, k
      double precision :: vv(nmax)
      double precision :: aamax, dum, sum
      save

      d = 1.d0
      do i = 1, n
          aamax = 0.d0
          do j = 1, n
              if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
          end do
          if (aamax.eq.0d0) stop 'Singular matrix.'
          vv(i) = 1.d0/aamax
      end do
      do j = 1, n
          if (j.gt.1) then
              do i = 1, j - 1
                  sum = a(i,j)
                  if (i.gt.1) then
                      do k = 1, i - 1
                          sum = sum - a(i,k)*a(k,j)
                      end do
                      a(i,j) = sum
                  end if

              end do
          end if

          aamax = 0.d0
          do i = j, n
              sum = a(i,j)
              if (j.gt.1) then
                  do k = 1, j - 1
                      sum = sum - a(i,k)*a(k,j)
                  end do
                  a(i,j) = sum
              end if

              dum = vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax = i
                  aamax = dum
              end if

          end do
          if (j.ne.imax) then
              do k = 1, n
                  dum = a(imax,k)
                  a(imax,k) = a(j,k)
                  a(j,k) = dum
              end do
              d = -d
              vv(imax) = vv(j)
          end if

          indx(j) = imax
          if (j.ne.n) then
              if (a(j,j).eq.0d0) a(j,j) = tiny
              dum = 1.d0/a(j,j)
              do i = j + 1, n
                  a(i,j) = a(i,j)*dum
              end do
          end if

      end do
      if (a(n,n).eq.0d0) a(n,n) = tiny
      return

end subroutine ludcmp

!----------------------------------------------------------------------
! mmid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mmid.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Integrates the dependent variables y from x_start to
! x_start + h_total in n_step increments (the "modified midpoint"
! substep of the Bulirsch-Stoer method; called by bsstep). Input are
! the y's and dy/dx's at x_start; values of y at x_start + h_total are
! stored in y_out. Derivatives are calculated in subroutine deriv.
! This SR from Numerical Recipes, p.562.
! b/pressure_rotation_factor/.../saha_state are opaque pass-through
! arguments forwarded unchanged to deriv (the caller-supplied
! derivative routine, e.g. qatm/qenv) -- named to match the actual
! arguments used at the bsstep call sites in atm/atm_lib.f90.
subroutine mmid(y, dydx, n_var, x_start, h_total, n_step, y_out, deriv, &
     luminosity_linear, pressure_rotation_factor, temperature_rotation_factor, &
     log10_gravity, in_atmosphere, want_derivatives, conductive_opacity_flag, &
     print_flag, log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     call_count, saha_state)
      implicit none

      double precision, intent(in) :: y(3), dydx(3)
      integer, intent(in) :: n_var
      double precision, intent(in) :: x_start, h_total
      integer, intent(in) :: n_step
      double precision, intent(out) :: y_out(3)
      external deriv
      double precision, intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(inout) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag, print_flag
      double precision, intent(inout) :: log10_radius, log10_teff, &
           hydrogen_fraction, metal_fraction
      integer, intent(inout) :: call_count, saha_state

! h_sub is the size of each small step.
      double precision :: y_mid(3), y_new(3)
      double precision :: h_sub, h_sub2, x_current, y_swap
      integer :: i, step_index
      save

      h_sub = h_total/dfloat(n_step)
! first step
      do i = 1,n_var
       y_mid(i) = y(i)
       y_new(i) = y(i) + dydx(i)*h_sub
      end do
      x_current = x_start + h_sub
! y_out temporarily used for storage of derivatives.
      call deriv(x_current, y_new, y_out, luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity, in_atmosphere, want_derivatives, &
           conductive_opacity_flag, print_flag, log10_radius, log10_teff, &
           hydrogen_fraction, metal_fraction, call_count, saha_state)
      h_sub2 = 2.0d0*h_sub
! general step.
      do step_index = 2,n_step
       do i = 1,n_var
          y_swap = y_mid(i) + h_sub2*y_out(i)
          y_mid(i) = y_new(i)
          y_new(i) = y_swap
       end do
       x_current = x_current + h_sub
       call deriv(x_current, y_new, y_out, luminosity_linear, &
            pressure_rotation_factor, temperature_rotation_factor, &
            log10_gravity, in_atmosphere, want_derivatives, &
            conductive_opacity_flag, print_flag, log10_radius, log10_teff, &
            hydrogen_fraction, metal_fraction, call_count, saha_state)
      end do
! last step.
      do i = 1,n_var
       y_out(i) = 0.5d0*(y_mid(i) + y_new(i) + h_sub*y_out(i))
      end do
      return
end subroutine mmid

!----------------------------------------------------------------------
! osplin
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original osplin.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! ACM Algorithm 574: shape-preserving osculatory quadratic splines,
! by D.F. McAllister and J.A. Roulier, ACM Transactions on
! Mathematical Software, September 1981.
!
! xtab contains the abscissas of the points of interpolation.
! ytab contains the ordinates of the points of interpolation.
! n is the number of data points.
! k is the number of points at which the spline is to be evaluated
! (the points themselves are xval; the evaluated values come back in
! yval).
!
! Upon exit from subroutine 'slopes' -- first_derivs contains the
! computed first derivatives at each data point.
subroutine osplin(xval, yval, xtab, ytab, n, k)
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: xval(json)
      double precision, intent(out) :: yval(json)
      double precision, intent(in) :: xtab(json), ytab(json)
      integer, intent(in) :: n, k

      double precision :: first_derivs(json), eps
      integer :: err
      save

! calculate the slopes at each data point.
      call slopes(xtab, ytab, first_derivs, n)

! set the error tolerance eps, which is used in subroutine 'choose'.
      eps = 1.d-04
! call meval to evaluate the spline at the run of points xval.
      call meval(xval, yval, xtab, ytab, first_derivs, n, k, eps, err)

      return
end subroutine osplin

!----------------------------------------------------------------------
! polint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original polint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Neville's algorithm polynomial interpolation/extrapolation (as in
! Numerical Recipes POLINT), returning the interpolated value y at x
! and an error estimate dy.
! MHP 10/02: dimensions changed for consistency with fpft (fixed at
! 20 rather than sized to n).
subroutine polint(xa, ya, n, x, y, dy)
      implicit none

      double precision, intent(in) :: xa(20), ya(20)
      integer, intent(in) :: n
      double precision, intent(in) :: x
      double precision, intent(out) :: y, dy

      double precision :: c(20), d(20)
      integer :: ns, i, j
      double precision :: dif, dift, ho, hp, w, den
      save

      ns = 1
      dif = dabs(x-xa(1))
      do i = 1, n
       dift = dabs(x-xa(i))
       if (dift.lt.dif) then
          ns = i
          dif = dift
       end if
       c(i) = ya(i)
       d(i) = ya(i)
      end do
      y = ya(ns)
      ns = ns - 1
      do j = 1, n-1
       do i = 1, n-j
          ho = xa(i)-x
          hp = xa(i+j)-x
          w = c(i+1) - d(i)
          den = ho - hp
          if (dabs(den).lt.1.0d-20) stop
          den = w/den
          d(i) = hp*den
          c(i) = ho*den
       end do
       if (2*ns.lt.n-j) then
          dy = c(ns+1)
       else
          dy = d(ns)
          ns = ns-1
       end if
       y = y + dy
      end do
      return
end subroutine polint

!----------------------------------------------------------------------
! quint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original quint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Quadratic interpolation for equidistant points.
! y0=y(x0), y1=y(x1), y2=y(x2); h=x1-x0=x2-x1; computes y=y(x).
subroutine quint(x, x0, h, y0, y1, y2, y)
      implicit none
      double precision, intent(in) :: x, x0, h, y0, y1, y2
      double precision, intent(out) :: y

      double precision :: d1, d2, t
      save

      d1 = y1 - y0
      d2 = y2 - 2.d0*y1 + y0
      t = (x - x0)/h
      y = y0 + t*d1 + 0.5d0*t*(t-1.d0)*d2
      return
end subroutine quint

!----------------------------------------------------------------------
! splinc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splinc.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Natural cubic spline coefficient generator, same algorithm as
! splinj but with x/y/y2/u dimensioned to the json=5000 maximum
! rather than to n, exactly as in the original file.
subroutine splinc(x, y, y2, n)
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: x(json), y(json)
      double precision, intent(out) :: y2(json)

      double precision :: u(json)
      integer :: i, k
      double precision :: sig, p, qn, un
      save

! natural spline
      y2(1) = 0.0d0
      u(1) = 0.0d0
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.0d0
         y2(i) = (sig-1.0d0)/p
         u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      qn = 0.0d0
      un = 0.0d0
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k = n-1, 1, -1
         y2(k) = y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine splinc

!----------------------------------------------------------------------
! splinj
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splinj.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Natural cubic spline coefficient generator, sized exactly to n
! (unlike splinc, which is dimensioned to the json=5000 maximum).
subroutine splinj(x, y, y2, n)
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: x(n), y(n)
      double precision, intent(out) :: y2(n)

      double precision :: u(n)
      integer :: i, k
      double precision :: sig, p, qn, un

! natural spline
      y2(1) = 0.0d0
      u(1) = 0.0d0
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.0d0
         y2(i) = (sig-1.0d0)/p
         u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      qn = 0.0d0
      un = 0.0d0
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k = n-1, 1, -1
         y2(k) = y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine splinj

!----------------------------------------------------------------------
! splnr
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splnr.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Single-precision (real, not double precision) natural/clamped cubic
! spline coefficient generator, the classic Numerical Recipes SPLINE
! routine kept in its original real precision (unlike cspline/splinj/
! splinc, which are real*8 ports of the same algorithm).
subroutine splnr(x, y, n, yp1, ypn, y2)
      implicit none
      integer, parameter :: nmax = 500

      integer, intent(in) :: n
      real, intent(in) :: x(n), y(n), yp1, ypn
      real, intent(out) :: y2(n)

      integer :: i, k
      real :: p, qn, sig, un, u(nmax)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
             /(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine splnr

!----------------------------------------------------------------------
! tridiag_gs
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original tridiag_gs.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0
! regression suite (examples/run_standard_solar_model).
!
! Thomas-algorithm tridiagonal solve (see also tridia.f90, the same
! algorithm operating through a shared common block instead of
! explicit dummy arguments): solves the tridiagonal system with
! sub-diagonal a, diagonal b, super-diagonal c, and right-hand side
! ex_prime, returning the solution in ex.
subroutine tridiag_gs(a, b, c, ex_prime, npt, ex)
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: a(json), b(json), c(json), &
           ex_prime(json)
      integer, intent(in) :: npt
      double precision, intent(out) :: ex(json)

      double precision :: gama(json)
      integer :: j
      double precision :: bet
      save

      bet = b(1)
      ex(1) = ex_prime(1)/bet
      do j = 2, npt
         gama(j) = c(j-1)/bet
         bet = b(j) - a(j)*gama(j)
         if (bet.eq.0) stop '#TRIDIA:SINGULAR MATRIX'
         ex(j) = (ex_prime(j) - a(j)*ex(j-1))/bet
      end do
      do j = npt-1, 1, -1
         ex(j) = ex(j) - gama(j+1)*ex(j+1)
      end do
!  911  format(1p6e10.2)

      return
end subroutine tridiag_gs

!----------------------------------------------------------------------
! ysplin
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ysplin.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! YCK 3/91. Find the coefficients for the natural cubic spline.
! On entry c(1,i) holds the tabulated function values at xi(i); on
! exit c(2,i)/c(3,i)/c(4,i) hold the first/second/third-order
! coefficients of the cubic on each sub-interval.
subroutine ysplin(xi, c, n)
      implicit none
      integer, parameter :: np = 100

      double precision, intent(in) :: xi(n)
      double precision, intent(inout) :: c(4,np)
      integer, intent(in) :: n

      double precision :: f(np), h(np), d(np), g, d3
      double precision :: const1, const2, const3, const4
      integer :: i
      save

! set the divided difference at each subinterval.
      do i = 2, n
         h(i) = xi(i)-xi(i-1)
         f(i) = (c(1,i)-c(1,i-1))/h(i)
      end do
! set the tridiagonally dominant matrix and the matrix equation
! for the natural spline.
!!!!! find the coefficients for the second order terms !!!!!
      i = 2
      const1 = 3.0d0
      const2 = 1.5d0
      const3 = 2.0d0
      const4 = 1.5d0
      c(2,i) = const1*h(i)*f(i+1)+const2*h(i+1)*f(i)
      d(i) = const3*h(i)+const4*h(i+1)
      do i = 3, n-2
         const1 = 3.0d0
         const2 = 3.0d0
         const3 = 2.0d0
         const4 = 2.0d0
         c(2,i) = const1*h(i)*f(i+1)+const2*h(i+1)*f(i)
         d(i) = const3*h(i)+const4*h(i+1)
      end do
      i = n-1
      const1 = 1.5d0
      const2 = 3.0d0
      const3 = 1.5d0
      const4 = 2.0d0
      c(2,i) = const1*h(i)*f(i+1)+const2*h(i+1)*f(i)
      d(i) = const3*h(i)+const4*h(i+1)
! solve the matrix equation with gauss method
! elimination of the sub-diagonal
      do i = 2, n-2
         g = h(i+2)/d(i)
         d(i+1) = d(i+1)-g*h(i)
         c(2,i+1) = c(2,i+1)-g*c(2,i)
      end do
! elimination of the super-diagonal
      c(2,n-1) = c(2,n-1)/d(n-1)
      do i = n-2, 2, -1
         c(2,i) = (c(2,i)-h(i)*c(2,i+1))/d(i)
      end do
! treatment for the first and last row.
      c(2,1) = 1.5d0*f(2)-c(2,2)/2.0d0
      c(2,n) = 1.5d0*f(n)-c(2,n-1)/2.0d0
!!! now, we have the coefficients for the first order terms
!
! find the coefficients for the second, and the third order terms
      do i = 1, n-1
         d3 = c(2,i)+c(2,i+1)-2.0d0*f(i+1)
         c(3,i) = (f(i+1)-c(2,i)-d3)/h(i+1)
         c(4,i) = d3/(h(i+1)*h(i+1))
      end do
!!!   now, we have complete set of coefficients for each sub-interval.

      return
end subroutine ysplin

!----------------------------------------------------------------------
! ctridi
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ctridi.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Solves a tridiagonal matrix system for the fractional abundance of a
! species. This routine is from Numerical Recipes, p.40.
!
! Analogue of tridia.f90 (the Henyey structure-equation tridiagonal
! solver): same Thomas-algorithm layout, used here for compositional
! (species-abundance) diffusion solves instead. Was originally shared
! with tridia.f90 via common/tridi/ (positional storage); converted
! (2026, GUIDELINES.md) to explicit arguments since this is real
! per-call data flow (matrix in, solution out), not global
! configuration -- see GUIDELINES.md's module-vs-argument distinction.
subroutine ctridi(n, sub_diag, diag, super_diag, rhs, solution)

      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution (originally U, the new species abundance) is
! solved for here.
      double precision, intent(in) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: solution(json)

! gamma_elim is solver-internal scratch, never read by any caller
! (unlike tridia.f90's gamma_elim(n), which doubled as a genuine
! cross-call input -- see tridia.f90's header note). SAVE preserved
! from the original common-block version though nothing here actually
! depends on it persisting across calls.
      double precision :: gamma_elim(json)
      save

      double precision :: bet
      integer :: j

      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2,n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) stop '#TRIDIA:SINGULAR MATRIX'
         solution(j) = (rhs(j) - sub_diag(j)*solution(j-1))/bet
      end do
      do j = n-1,1,-1
         solution(j) = solution(j) - gamma_elim(j+1)*solution(j+1)
      end do

      return
end subroutine ctridi

!----------------------------------------------------------------------
! tridia
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original tridia.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Thomas-algorithm tridiagonal solve (see also tridiag_gs.f90, the
! same algorithm operating on explicit dummy arguments instead of a
! shared common block, and ctridi.f90, the same algorithm for
! compositional diffusion). The caller fills sub_diag/diag/super_diag/
! rhs before calling; tridia solves for solution(:) in place and also
! returns dj(:), the per-zone contribution of the change in solution
! weighted by ei(:), and sumdj, their sum.
!
! Was originally shared with ctridi.f90 via common/tridi/ (positional
! storage); converted (2026, GUIDELINES.md) to explicit arguments
! since this is real per-call data flow, not global configuration.
! dj_n_seed replaces what was previously smuggled in via
! gamma_elim(n): the caller (dcoeft.f90, via seculr.f90) computes a
! surface wind-angular-momentum-loss term and used to stash it in the
! shared common block's gamma_elim(num_eq_points) slot specifically so
! this routine's very first statement (before gamma_elim is
! overwritten as pure solver scratch below) could pick it up as dj(n)'s
! initial value. That's now an explicit input instead of a COMMON
! side-channel; see dcoeft.f90's matching surface_wind_loss_term
! output argument.
! KC 2025-05-31 removed the unused ej dummy argument (see the
! commented-out original signature below).
!       SUBROUTINE TRIDIA(N,EI,EJ,DJ,SUMDJ)  ! KC 2025-05-31
subroutine tridia(n, ei, dj, sumdj, sub_diag, diag, super_diag, rhs, &
     solution, dj_n_seed)
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: ei(json)
      double precision, intent(out) :: dj(json)
      double precision, intent(out) :: sumdj
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution is solved for here. dj_n_seed is dj(n)'s initial
! value (see header note above).
      double precision, intent(in) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: solution(json)
      double precision, intent(in) :: dj_n_seed

! gamma_elim is solver-internal scratch (SAVE preserved from the
! original common-block version though nothing here actually depends
! on it persisting across calls, now that dj_n_seed carries the one
! value that used to be read from it across calls).
      double precision :: gamma_elim(json)

      double precision :: rhs_orig(json)
      integer :: i, j
      double precision :: bet, fj
      save

      dj(n) = dj_n_seed
      do i = 1, n
         rhs_orig(i) = rhs(i)
      end do
      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2, n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) stop '#TRIDIA:SINGULAR MATRIX'
         solution(j) = (rhs(j) - sub_diag(j)*solution(j-1))/bet
      end do
      dj(n) = dj(n)+(solution(n)-rhs_orig(n))*ei(n)
      sumdj = dj(n)
      fj = 1.0d0+(solution(n)-rhs_orig(n))/rhs_orig(n)
      solution(n) = rhs_orig(n)*fj
      do j = n-1, 1, -1
         solution(j) = solution(j) - gamma_elim(j+1)*solution(j+1)
         dj(j) = (solution(j)-rhs_orig(j))*ei(j)
         sumdj = sumdj + dj(j)
      end do

      return
end subroutine tridia

!----------------------------------------------------------------------
! bsstep
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original bsstep.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Bulirsch-Stoer step-size-control driver: repeatedly calls mmid at
! increasing substep counts (substep_sequence), extrapolates the
! results to zero step size via ratext, and either accepts the step
! (returning x0, h_did, h_next) or shrinks h and retries. This SR is
! the classic Numerical-Recipes BSSTEP algorithm.
! luminosity_linear/pressure_rotation_factor/.../saha_state are opaque
! pass-through arguments forwarded unchanged to mmid/deriv -- named to
! match the actual arguments used at the bsstep call sites in
! atm/atm_lib.f90.
subroutine bsstep(y, dydx, num_eqs, indep_var, h_step, tolerance, y_scale, &
     h_did, h_next, deriv, luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, in_atmosphere, &
     want_derivatives, conductive_opacity_flag, print_flag, log10_radius, &
     log10_teff, hydrogen_fraction, metal_fraction, call_count, saha_state, &
     step_err)
      use intpar_lib
      implicit none

      double precision, parameter :: one = 1.0d0, shrink_factor = 0.95d0, &
           grow_factor = 1.2d0

      double precision, intent(inout) :: y(3)
      double precision, intent(in) :: dydx(3)
      integer, intent(in) :: num_eqs
      double precision, intent(inout) :: indep_var
      double precision, intent(in) :: h_step, tolerance, y_scale(3)
      double precision, intent(out) :: h_did, h_next
      external deriv
      double precision, intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(inout) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag, print_flag
      double precision, intent(inout) :: log10_radius, log10_teff, &
           hydrogen_fraction, metal_fraction
      integer, intent(inout) :: call_count, saha_state
      double precision, intent(out) :: step_err(3)

      double precision :: y_err(3), y_sav(3), dy_sav(3), y_seq(3)
      integer :: substep_sequence(11)
      double precision :: h, x_sav, x_est, err_max
      integer :: i, j
      save
      data substep_sequence /2,4,6,8,12,16,24,32,48,64,96/

      h = h_step
      x_sav = indep_var
      do i = 1,num_eqs
       y_sav(i) = y(i)
       dy_sav(i) = dydx(i)
      end do
   20 do i = 1,max_stage_index
       call mmid(y_sav, dy_sav, num_eqs, x_sav, h, substep_sequence(i), &
            y_seq, deriv, luminosity_linear, pressure_rotation_factor, &
            temperature_rotation_factor, log10_gravity, in_atmosphere, &
            want_derivatives, conductive_opacity_flag, print_flag, &
            log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
            call_count, saha_state)
       x_est = (h/substep_sequence(i))**2
       call ratext(i, x_est, y_seq, y, y_err, num_eqs, extrap_order)
       err_max = 0.0d0
       do j = 1,num_eqs
          err_max = dmax1(err_max, dabs(y_err(j)/y_scale(j)))
          step_err(j) = dabs(y_err(j)/y_scale(j))
       end do
       err_max = err_max/tolerance
       if(err_max.lt.one) then
          indep_var = indep_var + h
          h_did = h
          if(i.eq.extrap_order) then
             h_next = h*shrink_factor
          else if (i.eq.extrap_order-1) then
             h_next = h*grow_factor
          else
             h_next = h*dfloat(substep_sequence(extrap_order-1))/ &
                  dfloat(substep_sequence(i))
          endif
          return
       endif
      end do
      h = 0.25d0*h/2.0d0**int((max_stage_index-extrap_order)/2)
!      H = 0.25D0*H/2**((IMAX-NUSE)/2)
      if(hydrogen_fraction+h.eq.hydrogen_fraction) then
         write(*,*) 'ERROR IN BSSTEP'
       stop
      end if
      goto 20

end subroutine bsstep

!----------------------------------------------------------------------
! intpol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intpol.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! ***** the interpolation routine *****
! this routine contains two interpolation methods; hermite,
! and spline.  each of these has own merit and demerit.
! if the interpolant is smooth, both of these will give good
! results.  generally, spline gives more smooth interpolation.
! when the interpolant contains abrupt variation in gradient,
! however, spline get worse at that part, while hermit still
! gives reasonable result.  unfortunately, there is no criterion
! for selection between these two methods.  therefore, this
! routine gives the right for selection to user.
! YCK 3/91
!
! x_grid  ; array of abscissa points
! y_grid  ; array of ordinate points
! n_grid  ; size of the arrays
! x_eval  ; a x-value at which we want to find y-value
! k_lo    ; the grid point smaller than and closest to x_eval
! y_eval  ; the value we want
! dy_eval ; the derivative value at x_eval
subroutine intpol(x_grid, y_grid, n_grid, x_eval, y_eval, dy_eval)
      use luout_lib
      implicit none
      integer, parameter :: np=100

      double precision, intent(in) :: x_grid(n_grid), y_grid(n_grid)
      integer, intent(in) :: n_grid
      double precision, intent(in) :: x_eval
      double precision, intent(out) :: y_eval, dy_eval

      double precision :: spline_coeff(4,np), dx_local, x_eval_copy
      double precision :: interp_value, interp_deriv
      integer :: i, k_lo, k_hi, k_mid
      data spline_coeff/400*0.0d0/
      save

! the coefficients for the zero-th order term
      do i=1,n_grid
         spline_coeff(1,i)=y_grid(i)
      end do
! find the coefficients for spline interploation
      call ysplin(x_grid, spline_coeff, n_grid)
      x_eval_copy=x_eval
! find the grid point,  k_lo, such that x_grid(k_lo)<=x_eval_copy, and
! abs(x_grid(k_lo)-x_eval_copy)<1.
      if(x_grid(1).gt.x_eval)then
         k_lo=1
         k_hi=2
         go to 522
      endif
      if(x_grid(n_grid).lt.x_eval)then
         k_lo=n_grid-1
         k_hi=n_grid
         go to 522
      endif
      k_lo=1
      k_hi=n_grid
    2 if((k_hi-k_lo).gt.1)then
         k_mid = (k_hi+k_lo)/2
         if(x_grid(k_mid).gt.x_eval_copy)then
            k_hi=k_mid
         else
            k_lo=k_mid
         endif
         go to 2
      endif
      if((k_hi-k_lo).le.0)then
         write(iowr, *) 'ERROR COX OP: INTERPOLATION'
         write(short_file_unit, *) 'ERROR COX OP: INTERPOLATION'
         stop
      endif
  522 continue
! now, (k_lo,k_hi) is sub-range of x_grid which contains x_eval_copy.
      dx_local=x_eval_copy-x_grid(k_lo)
! go on to the spline interpolation routine.
! evaluates the interpolation value in the sub-range we determined.
      interp_value=((spline_coeff(4,k_lo)*dx_local+spline_coeff(3,k_lo))*dx_local &
           +spline_coeff(2,k_lo))*dx_local+spline_coeff(1,k_lo)
      interp_deriv=(3.0d0*spline_coeff(4,k_lo)*dx_local+2.0d0*spline_coeff(3,k_lo)) &
           *dx_local+spline_coeff(2,k_lo)
! return the results from the spline routine
      y_eval=interp_value
      dy_eval=interp_deriv

      return
end subroutine intpol

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
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: xa(json), ya(json), y2a(json), x
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
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0
      return
end subroutine splintd2

!----------------------------------------------------------------------
! trapzd
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original trapzd.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Numerical Recipes-style trapzd: on the n=1 call, evaluates the
! integrand q at the single point r0=b2 and returns the crude
! 2-point trapezoidal estimate s over [b1,b2]. On subsequent calls
! (n>1) it refines s by adding the contribution of the newly
! inserted midpoints, with the interpolated intermediate rho/sm/w2/
! eta22 values linearly interpolated between the (rhop,smp,w2p,
! eta22p) state at b1 and the (rho,sm,w2,eta22) state at b2. it is
! the number of subintervals from the previous call, held across
! calls via SAVE (the Numerical Recipes call-count doubling scheme;
! it used to be passed in explicitly as a dummy argument i -- see
! the commented-out original signature below).
!       SUBROUTINE TRAPZD(B1,B2,S,N,RHO,RHOP,SM,SMP,W2,W2P,ETA22,
!      *ETA22P,Q,QP,I)  ! KC 2025-05-31
subroutine trapzd(b1, b2, s, n, rho, rhop, sm, smp, w2, w2p, eta22, &
     eta22p, q, qp)

      use const_lib
      implicit none



      double precision, intent(in) :: b1, b2
      double precision, intent(inout) :: s
      integer, intent(in) :: n
      double precision, intent(in) :: rho, rhop, sm, smp, w2, w2p, &
           eta22, eta22p
      double precision, intent(inout) :: q
      double precision, intent(in) :: qp

      integer :: it
      save

      double precision :: r0, r03, tnm, dr, del, y, sum, drho, dm, &
           deta2, dw2, r03t, rhot, smt, w2t, eta22t, q0
      integer :: j

      r0 = b2
      r03 = r0**3
      if (n.eq.1) then
!  aint = int(0=>r0) (rho/m)*r0'**7*omega**2*(5+eta2)/(2+eta2) dr0'
!  q is the integrand (ro'**7,etc.) evaluated at r0(i)
!  aint and its derivatives w/r/to r and theta are needed to find <g>
       q = (rho*w2*r03*(3.0d0+eta22)/(sm*eta22))*r03*r0
!        q(i) = dexp(cln*(hd(i)-hs(i)))*omega(i)**2*r0(i)**6
!    *   *(5.0d0+eta2(i))/(2.0d0+eta2(i))
       s = 0.5d0*(b2-b1)*(qp+q)
       it = 1
      else
       tnm = dfloat(it)
       dr = b2 - b1
       del = dr/tnm
       y = b1 + 0.5d0*del
       sum = 0.0d0
       drho = (rho - rhop)/dr
       dm = (sm - smp)/(b2**2 - b1**2)
       deta2 = (eta22 - eta22p)/dr
       dw2 = (w2 - w2p)/dr
       do j = 1, it
          r03t = y**3
! interpolate rho,m,omega,eta2+2 between shell i and shell i-1
          rhot = rhop+drho*del
          smt = smp+dm*(y**2 - b1**2)
          w2t = w2p + dw2*del
          eta22t = eta22p + deta2*del
! calculate q between shells
          q0 = (rhot*w2t*r03t*(3.0d0+eta22t)/(smt*eta22t))*r03t*y
! q0 = rho*w2*r07t*(3.0d0+eta22)/(sm*eta22)
          sum = sum + q0
          y = y+del
       end do
       s = 0.5d0*(s+del*sum)
       it = it*2
      end if

      return
end subroutine trapzd

!----------------------------------------------------------------------
! qgauss
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original qgauss.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 10-point (5-point-symmetric) Gauss-Legendre quadrature over [0,b]
! of three related integrands evaluated by the external routine
! func: g0g = integral of g*s, ginvg = integral of s/g, and sphig =
! integral of s, where (g,s) and (g2,s2) are func's outputs at the
! two symmetric abscissas about the interval midpoint. r0, hs, aint,
! q, w2, a, and i are passed through unchanged to func on each call.
subroutine qgauss(g0g, ginvg, sphig, b, r0, hs, aint, q, w2, a, i)

      use const_lib
      implicit none
      integer, parameter :: json = 5000



      double precision, intent(out) :: g0g, ginvg, sphig
      double precision, intent(in) :: b
      double precision, intent(in) :: r0(json), hs(json)
      double precision, intent(in) :: aint, q, w2, a
      integer, intent(in) :: i

      double precision :: x(5), w(5)
      data x/.14887433898163d0,.43339539412925d0,.67940956829902d0, &
           .86506336668899d0,.97390652851717d0/
      data w/.29552422471475d0,.26926671931d0,.21908636251598d0, &
           .14945134915058d0,.06667134430869d0/

      double precision :: xm, xr, dx, g, s, g2, s2
      integer :: j
      save

      xm = 0.5d0*b
      xr = xm
      g0g = 0.0d0
      ginvg = 0.0d0
      sphig = 0.0d0
      do j = 1, 5
       dx = xr*x(j)
       call func(xm+dx, g, s, r0, hs, aint, q, w2, a, i)
       call func(xm-dx, g2, s2, r0, hs, aint, q, w2, a, i)
       g0g = g0g+w(j)*(g*s+g2*s2)
       ginvg = ginvg+w(j)*(s/g+s2/g2)
       sphig = sphig+w(j)*(s+s2)
      end do
      g0g = g0g*xr
      ginvg = ginvg*xr
      sphig = sphig*xr
      return
end subroutine qgauss

!------------------------    GROUP: SR_P     -------------------------------
!----------------------------------------------------------------------
! intpt
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intpt.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). NOT exercised by that suite
! (only called from mhdpx2.f90, the MHD equation-of-state path, which
! the suite's 4 test cases do not select) -- verified by build +
! code review only.
!
! Bicubic-style interpolation in a (log10_pressure, log10_temperature)
! table_data(table_dim_t,table_dim_r,num_vars) grid, via two passes of
! the external 4-point Lagrangian interpolator lir: first at fixed
! t_indices columns to interpolate in pressure (work1 -> work2), then
! across the 4 selected t_indices to interpolate in temperature
! (work2 -> interp_vars). Dummy-argument names match the actual
! arguments used at the intpt call sites in mhdpx2.f90.
subroutine intpt(log10_pressure, log10_temperature, table_data, &
     table_dim_t, table_dim_r, num_vars, table_log10t, num_t, num_r, &
     work1, work2, y_work, interp_vars)
      implicit none
      integer, intent(in) :: table_dim_t, table_dim_r, num_vars
      double precision, intent(in) :: log10_pressure, log10_temperature
      double precision, intent(in) :: table_data(table_dim_t,table_dim_r,num_vars)
      double precision, intent(in) :: table_log10t(table_dim_t)
      integer, intent(in) :: num_t, num_r
      double precision, intent(inout) :: work1(num_vars,4), work2(num_vars,4)
      double precision, intent(inout) :: y_work(num_vars)
      double precision, intent(out) :: interp_vars(num_vars)

      integer :: r_lo_guess(4), r_indices(4,4), t_indices(4)
      double precision :: x_nodes(4)
! lir_order is INTEGER*4 in the original (overriding this file's
! IMPLICIT LOGICAL*4(L) for the single name L), holding a flag passed
! to the external routine lir; its exact meaning there is not
! established from this file alone.
      integer :: lir_order
      save

      integer :: n, i, m, j, t_col, iv, t_idx, r_idx, t_idx_max, r_idx_max
      double precision :: p_min, p_max
      integer :: lir_num_vars, lir_leading_dim, lir_num_points, lir_interp_mode

      do 100 n=1,num_t
         if(table_log10t(n).ge.log10_temperature) goto 101
         t_indices(1)=n
 100  continue
 101  if(t_indices(1).ge.2) t_indices(1)=t_indices(1)-1
      t_idx_max=num_t-3
      if(t_indices(1).gt.t_idx_max) t_indices(1)=t_idx_max
      do i=2,4
         t_indices(i)=t_indices(1)+i-1
      end do
      do i=1,4
         r_lo_guess(i)=1
         t_idx=t_indices(i)
         p_min=table_data(t_idx, 1,2)
         p_max=table_data(t_idx,num_r,2)
         if(log10_pressure.gt.p_max) then
            return
         end if
         do 200 m=1,num_r
            if(table_data(t_idx,m,2).ge.log10_pressure) goto 201
            r_lo_guess(i)=m
 200     continue
 201     if(r_lo_guess(i).ge.2) r_lo_guess(i)=r_lo_guess(i)-1
         r_idx_max=num_r-3
         if(r_lo_guess(i).gt.r_idx_max) r_lo_guess(i)=r_idx_max
      end do
      do i=1,4
         do j=1,4
            r_indices(j,i)=r_lo_guess(i)+j-1
         end do
      end do
      do t_col=1,4
         t_idx=t_indices(t_col)
         do i=1,4
            r_idx=r_indices(i,t_col)
            x_nodes(i)=table_data(t_idx,r_idx,2)
         end do
         do i=1,4
            r_idx=r_indices(i,t_col)
            do iv=1,num_vars
               work1(iv,i)=table_data(t_idx,r_idx,iv)
            end do
         end do

         lir_num_vars=num_vars
         lir_leading_dim=num_vars
         lir_num_points=4
         lir_order=1
         lir_interp_mode=1
         call lir(log10_pressure, x_nodes, y_work, work1, lir_num_vars, &
              lir_leading_dim, lir_num_points, lir_order, lir_interp_mode)
         do iv=1,num_vars
            work2(iv,t_col) = y_work(iv)
         end do

      end do
      do i=1,4
         t_idx=t_indices(i)
         x_nodes(i) = table_log10t(t_idx)
      end do

      lir_num_vars=num_vars
      lir_leading_dim=num_vars
      lir_num_points=4
      lir_order=1
      lir_interp_mode=1

      call lir(log10_temperature, x_nodes, interp_vars, work2, lir_num_vars, &
           lir_leading_dim, lir_num_points, lir_order, lir_interp_mode)

      return
end subroutine intpt

!----------------------------------------------------------------------
! lir
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lir.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Moved here (2026) from nuclear_lib.f90: a generic table-lookup
! interpolation/extrapolation routine with no nuclear-physics content
! -- called by this module's own intpt above, and by eos/mhd/mhdpx1.f90;
! it was mis-homed in nuclear_lib purely because nuclear/liburn.f90
! also happens to use it (see ratext below, which shares that
! history). Naming/module placement follows GUIDELINES.md's rule that
! folder/module placement should track function, not caller.
!
! Generic table-lookup interpolation/extrapolation routine, entered
! either as LIR (cubic interpolation/extrapolation, unless num_points
! < 4) or as LIR1 (always linear interpolation/extrapolation).
!
! FOR A SUCH THAT target_z=table_z(A),  SETS result_y(I)=table_y(I,A), I=1,num_y
! table_z(N),table_y(I,N) MUST BE SUPPLIED FOR N=1,num_points AND I=1,num_y
! y_stride IS FIRST DIMENSION OF table_y
! interp_flag IS SET TO 1 FOR INTERPOLATION AND 0 FOR EXTRAPOLATION
! IF continue_search.LE.1, SCAN TO FIND THE table_z(N) WHICH IMMEDIATELY
!   BOUND target_z, STARTING AT N=1
! IF continue_search.GT.1, SCAN STARTS FROM VALUE OF N FROM PREVIOUS
!   CALL OF LIR
! NOTE
! MOST OF THE COMPUTATION IS PERFORMED IN SINGLE PRECISION
subroutine lir(target_z,table_z,result_y,table_y,num_y,y_stride, &
     num_points,continue_search,interp_flag)

      implicit none

      double precision, intent(in) :: target_z
      double precision, intent(in) :: table_z(*)
      double precision, intent(out) :: result_y(*)
      double precision, intent(in) :: table_y(*)
      integer, intent(in) :: num_y, y_stride, num_points, continue_search
      integer, intent(out) :: interp_flag

      double precision :: weight(4)
      integer :: search_idx
      data search_idx/-1/
      save

      integer :: linear_mode
      integer :: stride, stride_m1, y_strided, num_y_strided, table_end
      double precision :: diff
      integer :: pivot, closest
      double precision :: y1, y2, y3, y4, z1, z2, z3, z4, z12, z34
      integer :: y_idx, j, k, base_idx
      double precision :: yy

      linear_mode=0
      go to 1
      entry lir1(target_z,table_z,result_y,table_y,num_y,y_stride, &
           num_points,continue_search,interp_flag)
      linear_mode=1
    1 continue
      stride=1
! CHECK NT AND RESET IL IF NECESSARY
      if(num_points.lt.2) go to 101
      if(num_points.lt.4) linear_mode=1
! ADDRESSING CONSTANTS
      interp_flag=1
      stride_m1=stride-1
      y_strided=stride*y_stride
      num_y_strided=(num_y-1)*stride+1
      table_end=(num_points-1)*stride+1
      diff=table_z(table_end)-table_z(1)
! SET INDEX FOR START OF SEARCH
      search_idx=(search_idx-2)*stride+1
      if(continue_search.le.1.or.search_idx.lt.1) search_idx=1
! DETERMINE POSITION OF target_z WITHIN table_z
    2 if(search_idx.gt.table_end) go to 8
! KC 2025-05-30 fixed "Arithmetic IF statement"
!       IF(DIFF) 4,102,3
!     3 IF(ZI(N)-Z) 5,6,9
!     4 IF(ZI(N)-Z) 9,6,5
      if (diff .lt. 0.0) then
         goto 4
      else if (diff .eq. 0.0) then
         goto 102
      else
         goto 3
      end if
    3 if (table_z(search_idx) .lt. target_z) then
         goto 5
      else if (table_z(search_idx) .eq. target_z) then
         goto 6
      else
         goto 9
      end if
    4 if (table_z(search_idx) .lt. target_z) then
         goto 9
      else if (table_z(search_idx) .eq. target_z) then
         goto 6
      else
         goto 5
      end if
    5 search_idx=search_idx+stride
      go to 2
! SET Y WHEN Z LIES ON A MESH POINT
    6 base_idx=(search_idx-1)*y_stride
      do 7 y_idx=1,num_y_strided
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       Y(I)=YI(I+J)
!     7 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         result_y(y_idx)=table_y(y_idx+base_idx)
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
    7 continue
      go to 30
! CONTROL WHEN Z DOES NOT LIE ON A MESH POINT
    8 interp_flag=0
    9 if(search_idx.le.1) interp_flag=0
      if(linear_mode.eq.1) go to 20
! CUBIC INTERPOLATION/EXTRAPOLATION
! PIVOTAL POINT (M) AND POINT (K) CLOSEST TO Z
!    10 M=N
      pivot=search_idx
      closest=3
      if(search_idx.gt.1+stride) go to 11
      pivot=1+stride+stride
      closest=search_idx
   11 if(search_idx.lt.table_end) go to 12
      pivot=table_end-stride
      closest=4
! WEIGHTING FACTORS
   12 y1=table_z(pivot-stride*2)
      y2=table_z(pivot-stride)
      y3=table_z(pivot)
      y4=table_z(pivot+stride)
      z1=target_z-y1
      z2=target_z-y2
      z3=target_z-y3
      z4=target_z-y4
!    13 Z12=Z1*Z2
      z12=z1*z2
      z34=z3*z4
!    14 A(1)=Z2*Z34/((Y1-Y2)*(Y1-Y3)*(Y1-Y4))
      weight(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      weight(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      weight(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      weight(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
! CORRECT A(K)
!    15 DIFF=A(1)+A(2)+A(3)+A(4)
      diff=weight(1)+weight(2)+weight(3)+weight(4)
      weight(closest)=(1.d0+weight(closest))-diff
! COMPUTE Y
!    16 M=(M-1)/IR-3
      pivot=(pivot-1)/stride-3
      pivot=pivot*y_strided
      do 18 y_idx=1,num_y_strided
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       K=I+M
!       YY=0.D0
!       DO 17 J=1,4
!       K=K+IRD
!       DIFF=YI(K)
!    17 YY=YY+A(J)*DIFF
!       Y(I)=YY
!    18 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         k=y_idx+pivot
         yy=0.d0
         do 17 j=1,4
            k=k+y_strided
            diff=table_y(k)
            yy=yy+weight(j)*diff
   17    continue
         result_y(y_idx)=yy
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
   18 continue
      go to 30
! LINEAR INTERPOLATION/EXTRAPOLATION
   20 if(search_idx.eq.1) search_idx=1+stride
      if(search_idx.gt.table_end) search_idx=table_end
      z1=table_z(search_idx)
      y1=(z1-target_z)/(z1-table_z(search_idx-stride))
      y2=1.0d0-y1
      base_idx=(search_idx-1)*y_stride
      pivot=base_idx-y_strided
      do 21 y_idx=1,num_y_strided,stride
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       Y(I)=Y1*YI(I+M)+Y2*YI(I+J)
!    21 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         result_y(y_idx)=y1*table_y(y_idx+pivot)+y2*table_y(y_idx+base_idx)
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
   21 continue
! RESET N
   30 search_idx=(search_idx+stride-1)/stride
      return
! DIAGNOSTICS
  101 continue
      return
  102 continue
      return
!  1001 FORMAT(/1X,10('*'),5X,'THERE ARE FEWER THAN TWO DATA POINTS IN',
!      *      ' LIR     NT =',I4,5X,10('*')/)
!  1002 FORMAT(/1X,10('*'),5X,'EXTREME VALUES OF INDEPENDENT VARIABLE',
!      *      ' EQUAL IN LIR',5X,10('*')/16X,'ZI(1) =',1PE13.5,',   ',
!      *       'ZI(',I4,') =',1PE13.5/)
end subroutine lir

!----------------------------------------------------------------------
! ratext
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ratext.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Moved here (2026) from nuclear_lib.f90: a diagonal rational-function
! extrapolator, the textbook companion to this module's own bsstep/
! mmid Bulirsch-Stoer stepper above (same algorithm as subroutine
! RZEXTR in Numerical Recipes, p.566) -- bsstep already called it via
! an errant `use nuclear_lib`. Also called directly by
! nuclear/liburn.f90 to extrapolate a sequence of sub-stepped
! lithium/beryllium burning-rate estimates (indexed by decreasing step
! size) to the zero-step-size limit as the sub-stepping is refined.
subroutine ratext(est_index, x_est, y_est, y_extrap, y_err, num_vars, &
     max_use)
      implicit none
      integer, parameter :: imax=11, nmax=15, ncol=7

      integer, intent(in) :: est_index, num_vars, max_use
      double precision, intent(in) :: x_est
      double precision, intent(in) :: y_est(num_vars)
      double precision, intent(out) :: y_extrap(num_vars), y_err(num_vars)

      double precision :: x_hist(imax), tableau(nmax,ncol), fx(ncol)
      save

! SAME AS SR RZEXTR FROM NUMERICAL RECIPES, P.566.
!
! SAVE CURRENT INDEPENDENT VARIABLE.
      integer :: var_idx, k, k2, num_use
      double precision :: yy, v, c, b1, b, ddy

      x_hist(est_index) = x_est
      if (est_index.eq.1) then
         do var_idx = 1, num_vars
            y_extrap(var_idx) = y_est(var_idx)
            tableau(var_idx,1) = y_est(var_idx)
            y_err(var_idx) = y_est(var_idx)
         end do
      else
!        USE AT MOST max_use PREVIOUS MEMBERS.
         num_use = min(est_index,max_use)
         do k = 1, num_use-1
            fx(k+1) = x_hist(est_index-k)/x_est
         end do
!        EVALUATE NEXT DIAGONAL IN TABLEAU.
         do var_idx=1,num_vars
            yy = y_est(var_idx)
            v = tableau(var_idx,1)
            c = yy
            tableau(var_idx,1) = yy
            do k2 = 2, num_use
               b1 = fx(k2)*v
               b = b1 - c
!              CARE NEEDED TO AVOID DIVISION BY ZERO.
               if (b.ne.0d0) then
                  b = (c - v)/b
                  ddy = c*b
                  c = b1*b
               else
                  ddy = v
               end if
               v = tableau(var_idx,k2)
               tableau(var_idx,k2) = ddy
               yy = yy + ddy
            end do
            y_err(var_idx) = ddy
            y_extrap(var_idx) = yy
         end do
      end if
      return
end subroutine ratext

end module numerics_lib
