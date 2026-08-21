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
! arguments used at the bsstep call sites in envint.f90.
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

end module numerics_lib
