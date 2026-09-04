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
! (2026) from net_lib.f90, where they'd originally been filed
! (along with the rest of nuclear/) purely because net_lib.f90's
! own liburn is one of their callers -- both are generic numerics with no
! nuclear-physics content, and this module's own bsstep/intpt already
! called them via an errant `use net_lib` (ratext is in fact the
! textbook Numerical-Recipes companion to bsstep/mmid above). See
! GUIDELINES.md's rule that folder/module placement should track
! function, not caller.
module numerics_lib
! 2026 numerics-gate opt-in: callers that pass the optional ierr of
! the gated procedures (bsstep/ksplint/splint/splintd2/intpol) can
! distinguish a NUMERICS TERMINATION -- the historical "solution
! diverged, stop" mode, e.g. a bsstep envelope-integration failure
! that is the normal pinned ending of some configurations -- from a
! hard configuration/table error. The driver maps this negative code
! to a clean process exit (the legacy stop exited 0), while
! yrec_capi surfaces it to pyyrec as a distinct status.
      integer, parameter, public :: numerics_termination = -2
! 2026 fold-in: the last standalone files in numerics/ (meval and its
! McAllister-Roulier quadratic-spline helpers, slopes/safedivide,
! simeqc) moved in as module procedures, giving them explicit
! interfaces (compile-time argument checking). Only simeqc has a
! caller outside this module (mixing/solve_composition); the rest
! are private. splinnr was deleted outright: zero callers tree-wide.
      private :: cases, choose, search, spline, meval, slopes, safedivide
contains

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
      use star_info_lib, only: json
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: x(n), y(n), yp1, ypn
      double precision, intent(out) :: y2(n)

      double precision :: u(json)
      integer :: i, k
      double precision :: sig, p, qn, un

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

! find the 'index'
      if(index.lt.1.or.index.gt.n_grid)index=1
      found_index=index
      if(x_eval.lt.grid_x(found_index))then
         do j=found_index-1,1,-1
            if(grid_x(j).le.x_eval)then
               found_index=j
               index=found_index
               
               return
            endif
         end do
         found_index=-1
      else
         do j=found_index,n_grid-1
            if(grid_x(j+1).gt.x_eval)then
               found_index=j
               index=found_index
               
               return
            endif
         end do
         found_index = -n_grid
      endif
      index=found_index

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
! lagrange4
!----------------------------------------------------------------------
! Added 2026 (readability sweep, R3). Applies the 4-point Lagrange
! weights returned by interp/intrp2 to four consecutive table values:
!     w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4)
! summed left to right, exactly as the former inline expressions in
! io/model_to_equal.f90 and io/equal_to_model.f90 did (byte-pinned).
! Callers pass the stencil as an array section, y(k0:k0+3).
pure function lagrange4(w, y) result(s)
      implicit none
      double precision, intent(in) :: w(4), y(4)
      double precision :: s

      s = w(1)*y(1)+w(2)*y(2)+w(3)*y(3)+w(4)*y(4)
end function lagrange4

!----------------------------------------------------------------------
! kspline
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original kspline.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Fixed-size (4-point) natural cubic spline coefficient generator:
! cspline hardwired to n=4 points and the natural-spline boundary
! condition.
subroutine kspline(x, y, y2)
      implicit none
      integer, parameter :: nm = 4
      double precision, parameter :: natural_bc = 1.0d30

      double precision, intent(in) :: x(nm), y(nm)
      double precision, intent(out) :: y2(nm)

! 2026 readability sweep (R3): the former private copy of the
! natural-spline body was token-identical to cspline's natural
! branch (yp1/ypn > 0.99d30 sets y2(1)=u(1)=0 and qn=un=0), so this
! is now a thin wrapper; byte-pinned.
      call cspline(x, y, nm, natural_bc, natural_bc, y2)
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
subroutine ksplint(xa, ya, y2a, x, y, ierr)
      implicit none
      integer, parameter :: nm = 4

      double precision, intent(in) :: xa(nm), ya(nm), y2a(nm), x
      double precision, intent(out) :: y

      double precision :: h, a, b
      integer :: klo, khi, k

      integer, intent(out) :: ierr

      ierr = 0
      klo = 1
      khi = nm
    do while (khi-klo .gt. 1)
         k = (khi+klo)/2
         if (xa(k) .gt. x) then
            khi = k
         else
            klo = k
         end if
    end do
!      write(*,*) khi, klo, xa(khi), xa(klo), x
      h = xa(khi) - xa(klo)
      if (h .eq. 0d0) then
            print*, 'Ksplint failure'
            ! 2026 (ROADMAP.md stage 3): the historical stop became an ierr
            ! return (numerics has no facade -- each public procedure carries
            ! its own gate).
            ierr = 1
            return
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
! slack allowed when x falls just outside the table before the
! bracket is clamped to the end intervals
      double precision, parameter :: below_table_slack = 0.99d0, &
           above_table_slack = 1.01d0

      jl = 0
      ju = n+1
   do while (ju-jl.gt.1)
         jm = (ju+jl)/2
         if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
   end do
      j = jl
      if ((j .eq. 0) .and. (x .gt. below_table_slack*xx(1))) then
         j = 1
      end if
      if ((j .eq. n) .and. (x .lt. above_table_slack*xx(n))) then
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
subroutine ludcmp(a, n, np, indx, d, ierr)
      implicit none
      integer, parameter :: nmax = 100
      double precision, parameter :: tiny = 1.0d-20

      integer, intent(in) :: n, np
      double precision, intent(inout) :: a(np,np)
      integer, intent(out) :: indx(n)
      double precision, intent(out) :: d
      integer, intent(out) :: ierr

      integer :: i, imax, j, k
      double precision :: vv(nmax)
      double precision :: aamax, dum, sum

      ierr = 0

      d = 1.d0
      do i = 1, n
          aamax = 0.d0
          do j = 1, n
              if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
          end do
          if (aamax.eq.0d0) then
! 2026 audit (section 7): library stop -> ierr (numerics must not
! kill the process; the drivers decide).
             write(*,*) 'ludcmp: singular matrix'
             ierr = 1
             return
          end if
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
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     call_count, saha_state, ierr)
      implicit none

      double precision, intent(in) :: y(3), dydx(3)
      integer, intent(in) :: n_var
      double precision, intent(in) :: x_start, h_total
      integer, intent(in) :: n_step
      integer, intent(out) :: ierr
      double precision, intent(out) :: y_out(3)
      external deriv
      double precision, intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(inout) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag
      double precision, intent(inout) :: log10_radius, log10_teff, &
           hydrogen_fraction, metal_fraction
      integer, intent(inout) :: call_count, saha_state

! h_sub is the size of each small step.
      double precision :: y_mid(3), y_new(3)
      double precision :: h_sub, h_sub2, x_current, y_swap
      integer :: i, step_index

      ierr = 0
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
           conductive_opacity_flag, log10_radius, log10_teff, &
           hydrogen_fraction, metal_fraction, call_count, saha_state, ierr)
      if (ierr /= 0) return
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
            conductive_opacity_flag, log10_radius, log10_teff, &
            hydrogen_fraction, metal_fraction, call_count, saha_state, ierr)
       if (ierr /= 0) return
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
      use star_info_lib, only: json
      implicit none

      double precision, intent(in) :: xval(json)
      double precision, intent(out) :: yval(json)
      double precision, intent(in) :: xtab(json), ytab(json)
      integer, intent(in) :: n, k

      double precision :: first_derivs(json)
! error tolerance used in subroutine 'choose' (via meval)
      double precision, parameter :: eps = 1.d-04
      integer :: err

! calculate the slopes at each data point.
      call slopes(xtab, ytab, first_derivs, n)

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
subroutine polint(xa, ya, n, x, y, dy, ierr)
      implicit none

      double precision, intent(in) :: xa(20), ya(20)
      integer, intent(in) :: n
      double precision, intent(in) :: x
      double precision, intent(out) :: y, dy
      integer, intent(out) :: ierr

      double precision :: c(20), d(20)
      integer :: ns, i, j
      double precision :: dif, dift, ho, hp, w, den

      ierr = 0
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
          if (dabs(den).lt.1.0d-20) then
! 2026 audit (section 7): library stop -> ierr (numerics must not
! kill the process; the drivers decide).
             write(*,*) 'polint: repeated abscissa (den ~ 0)'
             ierr = 1
             return
          end if
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
! Natural cubic spline coefficient generator: cspline with the
! natural boundary condition and x/y/y2 dimensioned to the json=5000
! maximum rather than to n, as in the original file.
subroutine splinc(x, y, y2, n)
      use star_info_lib, only: json
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: x(json), y(json)
      double precision, intent(out) :: y2(json)
      double precision, parameter :: natural_bc = 1.0d30

! 2026 readability sweep (R3): the former private copy of the
! natural-spline body was token-identical to cspline's natural
! branch, so this is now a thin wrapper (the json-shaped dummies are
! kept for the callers); byte-pinned.
      call cspline(x, y, n, natural_bc, natural_bc, y2)
      return
end subroutine splinc

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
subroutine tridiag_gs(a, b, c, ex_prime, npt, ex, ierr)
      use star_info_lib, only: json
      implicit none

      double precision, intent(in) :: a(json), b(json), c(json), &
           ex_prime(json)
      integer, intent(in) :: npt
      double precision, intent(out) :: ex(json)
      integer, intent(out) :: ierr

      double precision :: gama(json)
      integer :: j
      double precision :: bet

      ierr = 0
      bet = b(1)
      ex(1) = ex_prime(1)/bet
      do j = 2, npt
         gama(j) = c(j-1)/bet
         bet = b(j) - a(j)*gama(j)
         if (bet.eq.0) then
! 2026 audit (section 7): library stop -> ierr (numerics must not
! kill the process; the drivers decide).
            write(*,*) 'tridia: singular matrix'
            ierr = 1
            return
         end if
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
subroutine ctridi(n, sub_diag, diag, super_diag, rhs, solution, ierr)
      use star_info_lib, only: json

      implicit none

      integer, intent(in) :: n
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution (originally U, the new species abundance) is
! solved for here.
      double precision, intent(in) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: solution(json)
      integer, intent(out) :: ierr

! gamma_elim is solver-internal scratch, never read by any caller
! (unlike tridia's gamma_elim(n), which doubled as a genuine
! cross-call input -- see tridia's header note). The blanket SAVE
! below is kept from the original common-block version; whether any
! result depends on it has not been tested against the pins, so it
! stays (see audit/readability-sweep-2026-09-03/io.md).
      double precision :: gamma_elim(json)
      save   ! INTENTIONAL: tridiagonal solver carry (untested whether load-bearing); byte-pinned by Stage-0

      double precision :: bet
      integer :: j

      ierr = 0
      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2,n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) then
! 2026 audit (section 7): library stop -> ierr (numerics must not
! kill the process; the drivers decide).
            write(*,*) 'tridia: singular matrix'
            ierr = 1
            return
         end if
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
! gamma_elim(n): the caller (am_diffusion_coeffs.f90, via secular_transport.f90) computes a
! surface wind-angular-momentum-loss term and used to stash it in the
! shared common block's gamma_elim(num_eq_points) slot specifically so
! this routine's very first statement (before gamma_elim is
! overwritten as pure solver scratch below) could pick it up as dj(n)'s
! initial value. That's now an explicit input instead of a COMMON
! side-channel; see am_diffusion_coeffs.f90's matching surface_wind_loss_term
! output argument.
! KC 2025-05-31 removed the unused ej dummy argument (see the
! commented-out original signature below).
!       SUBROUTINE TRIDIA(N,EI,EJ,DJ,SUMDJ)  ! KC 2025-05-31
subroutine tridia(n, ei, dj, sumdj, sub_diag, diag, super_diag, rhs, &
     solution, dj_n_seed, ierr)
      use star_info_lib, only: json
      implicit none

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
      integer, intent(out) :: ierr

! gamma_elim is solver-internal scratch (no SAVE: nothing here depends
! on it persisting across calls, now that dj_n_seed carries the one
! value that used to be read from it across calls).
      double precision :: gamma_elim(json)

      double precision :: rhs_orig(json)
      integer :: i, j
      double precision :: bet, fj

      ierr = 0
      dj(n) = dj_n_seed
      do i = 1, n
         rhs_orig(i) = rhs(i)
      end do
      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2, n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) then
! 2026 audit (section 7): library stop -> ierr (numerics must not
! kill the process; the drivers decide).
            write(*,*) 'tridia: singular matrix'
            ierr = 1
            return
         end if
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
     want_derivatives, conductive_opacity_flag, log10_radius, &
     log10_teff, hydrogen_fraction, metal_fraction, call_count, saha_state, &
     step_err, ierr)
      use intpar_lib
      use math_lib
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
           conductive_opacity_flag
      double precision, intent(inout) :: log10_radius, log10_teff, &
           hydrogen_fraction, metal_fraction
      integer, intent(inout) :: call_count, saha_state
      double precision, intent(out) :: step_err(3)

      double precision :: y_err(3), y_sav(3), dy_sav(3), y_seq(3)
      integer :: substep_sequence(11)
      double precision :: h, x_sav, x_est, err_max
      integer :: i, j
      save   ! INTENTIONAL: blanket SAVE kept from the original (x_sav/h etc. persist between calls); byte-pinned by Stage-0
      data substep_sequence /2,4,6,8,12,16,24,32,48,64,96/

      integer, intent(out) :: ierr
      integer :: jerr_integrand

      ierr = 0
      h = h_step
      x_sav = indep_var
      do i = 1,num_eqs
       y_sav(i) = y(i)
       dy_sav(i) = dydx(i)
      end do
      do    ! step-halving retry loop (was label 20)
      do i = 1,max_stage_index
       call mmid(y_sav, dy_sav, num_eqs, x_sav, h, substep_sequence(i), &
            y_seq, deriv, luminosity_linear, pressure_rotation_factor, &
            temperature_rotation_factor, log10_gravity, in_atmosphere, &
            want_derivatives, conductive_opacity_flag, &
            log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
            call_count, saha_state, jerr_integrand)
       ! integrand (eos/kap/gradient) failure inside the midpoint
       ! substeps: same treatment as a diverged step
       if (jerr_integrand /= 0) then
          ierr = jerr_integrand
          return
       end if
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
      h = 0.25d0*h/pow(2.0d0, int((max_stage_index-extrap_order)/2))
!      H = 0.25D0*H/2**((IMAX-NUSE)/2)
      if(hydrogen_fraction+h.eq.hydrogen_fraction) then
         write(*,*) 'ERROR IN BSSTEP'
       ! 2026 (ROADMAP.md stage 3): the historical stop became an ierr
       ! return (numerics has no facade -- each public procedure carries
       ! its own gate).
       ierr = 1
       return
      end if
      end do

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
subroutine intpol(x_grid, y_grid, n_grid, x_eval, y_eval, dy_eval, ierr)
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

! the coefficients for the zero-th order term
      integer, intent(out) :: ierr

      ierr = 0
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
      else if(x_grid(n_grid).lt.x_eval)then
         k_lo=n_grid-1
         k_hi=n_grid
      else
      k_lo=1
      k_hi=n_grid
    do while ((k_hi-k_lo).gt.1)
         k_mid = (k_hi+k_lo)/2
         if(x_grid(k_mid).gt.x_eval_copy)then
            k_hi=k_mid
         else
            k_lo=k_mid
         endif
    end do
      if((k_hi-k_lo).le.0)then
         write(terminal_unit, *) 'ERROR COX OP: INTERPOLATION'
         write(run_log_unit, *) 'ERROR COX OP: INTERPOLATION'
         ! 2026 (ROADMAP.md stage 3): the historical stop became an ierr
         ! return (numerics has no facade -- each public procedure carries
         ! its own gate).
         ierr = 1
         return
      endif
      end if
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
subroutine splint(xa, ya, n, y2a, x, y, klo, khi, ierr)
      use luout_lib
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: xa(n), ya(n), y2a(n), x
      double precision, intent(out) :: y
      integer, intent(out) :: klo, khi


      integer :: k
      double precision :: h, a, b

      integer, intent(out) :: ierr

      ierr = 0
      klo = 1
      khi = n
    do while (khi-klo .gt. 1)
         k = (khi+klo)/2
         if (xa(k) .gt. x) then
            khi = k
         else
            klo = k
         end if
    end do
      h = xa(khi) - xa(klo)
      if (h .eq. 0d0) then
           write(run_log_unit,*) 'ERROR IN SPLINT ROUTINE.'
         ! 2026 (ROADMAP.md stage 3): the historical stop became an ierr
         ! return (numerics has no facade -- each public procedure carries
         ! its own gate).
         ierr = 1
         return
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
! maximum rather than to n, as in the original file; the work is
! done by splint.
subroutine splintd2(xa, ya, n, y2a, x, y, klo, khi, ierr)
      use star_info_lib, only: json
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: xa(json), ya(json), y2a(json), x
      double precision, intent(out) :: y
      integer, intent(out) :: klo, khi
      integer, intent(out) :: ierr

! 2026 readability sweep (R3): the body was token-identical to
! splint's (6d0 vs 6.0d0 is the same constant), so this is now a
! thin wrapper keeping the json-shaped dummies; byte-pinned.
      call splint(xa, ya, n, y2a, x, y, klo, khi, ierr)
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

      implicit none



      double precision, intent(in) :: b1, b2
      double precision, intent(inout) :: s
      integer, intent(in) :: n
      double precision, intent(in) :: rho, rhop, sm, smp, w2, w2p, &
           eta22, eta22p
      double precision, intent(inout) :: q
      double precision, intent(in) :: qp

      integer :: it
      save   ! INTENTIONAL: NR refinement accumulator carried between successive calls; byte-pinned by Stage-0

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
!        q(i) = exp(cln*(hd(i)-hs(i)))*omega(i)**2*r0(i)**6
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
! 2026 (bugsweep Batch 2): the linear interpolants must be evaluated
! at the offset (y - b1) of the current abscissa, as smt already is;
! the inherited `*del` used the constant sub-interval width instead.
          rhot = rhop+drho*(y - b1)
          smt = smp+dm*(y**2 - b1**2)
          w2t = w2p + dw2*(y - b1)
          eta22t = eta22p + deta2*(y - b1)
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
subroutine qgauss(integrand, g0g, ginvg, sphig, b, r0, hs, aint, q, w2, a, i)
      use star_info_lib, only: json

      implicit none

! 2026 (phase four, step 2 -- ROADMAP.md): the integrand used to be a
! hard-coded call to rotation's shape-integrand `func`, the one
! backwards dependency that kept numerics from being a pure leaf. It
! is now a procedure dummy; rotation/shape/rotation_shape_factors.f90 passes `func` at
! the call site. Same argument protocol as before (assumed-size for
! the two model-shaped arrays, matching the historical implicit
! interface).
      interface
         subroutine integrand(colatitude, local_gravity, area_element, &
              r0, log_mass, aint, q, w2, a, i)
            double precision, intent(in) :: colatitude
            double precision, intent(out) :: local_gravity, area_element
            double precision, intent(in) :: r0(*), log_mass(*)
            double precision, intent(in) :: aint, q, w2, a
            integer, intent(in) :: i
         end subroutine integrand
      end interface

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
      save   ! INTENTIONAL: quadrature state (empirically load-bearing); byte-pinned by Stage-0

      xm = 0.5d0*b
      xr = xm
      g0g = 0.0d0
      ginvg = 0.0d0
      sphig = 0.0d0
      do j = 1, 5
       dx = xr*x(j)
       call integrand(xm+dx, g, s, r0, hs, aint, q, w2, a, i)
       call integrand(xm-dx, g2, s2, r0, hs, aint, q, w2, a, i)
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
! lir_order is lir's continue_search argument (1 = fresh scan).
      integer :: lir_order

      integer :: n, i, m, j, t_col, iv, t_idx, r_idx, t_idx_max, r_idx_max
      double precision :: p_min, p_max
      integer :: lir_num_vars, lir_leading_dim, lir_num_points, lir_interp_mode

! Bracket the temperature: the original's GOTO 101 out of this scan
! is an EXIT, not a RETURN (a RETURN left every output unset for any
! in-range temperature).  t_indices(1) starts at 1 so a temperature
! below the first table row clamps to the bottom stencil instead of
! reading a stale/zero index.
      t_indices(1)=1
      do n=1,num_t
         if (table_log10t(n).ge.log10_temperature) exit
         t_indices(1)=n
      end do
      if(t_indices(1).ge.2) t_indices(1)=t_indices(1)-1
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
         do m=1,num_r
            if(table_data(t_idx,m,2).ge.log10_pressure) exit
            r_lo_guess(i)=m
         end do
         if(r_lo_guess(i).ge.2) r_lo_guess(i)=r_lo_guess(i)-1
         r_idx_max=num_r-3
         if(r_lo_guess(i).gt.r_idx_max) r_lo_guess(i)=r_idx_max
      end do
      do i=1,4
         do j=1,4
            r_indices(j,i)=r_lo_guess(i)+j-1
         end do
      end do
! lir call settings, shared by both interpolation passes
      lir_num_vars=num_vars
      lir_leading_dim=num_vars
      lir_num_points=4
      lir_order=1
      lir_interp_mode=1
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
! Moved here (2026) from net_lib.f90: a generic table-lookup
! interpolation/extrapolation routine with no nuclear-physics content
! -- called by this module's own intpt above, and by eos/mhd/mhdpx1.f90;
! it was mis-homed in net_lib purely because its own liburn
! also happens to use it (see ratext below, which shares that
! history). Naming/module placement follows GUIDELINES.md's rule that
! folder/module placement should track function, not caller.
!
! Generic table-lookup interpolation/extrapolation routine: cubic
! interpolation/extrapolation, falling back to linear when num_points
! < 4. (The historical always-linear ENTRY LIR1 had no callers and was
! deleted in the 2026 readability sweep; lir_impl keeps its
! linear-mode switch, which lir always passes as 0.)
!
! FOR A SUCH THAT target_z=table_z(A),  SETS result_y(I)=table_y(I,A), I=1,num_y
! table_z(N),table_y(I,N) MUST BE SUPPLIED FOR N=1,num_points AND I=1,num_y
! y_stride IS FIRST DIMENSION OF table_y
! interp_flag IS SET TO 1 FOR INTERPOLATION AND 0 FOR EXTRAPOLATION
! IF continue_search.LE.1, SCAN TO FIND THE table_z(N) WHICH IMMEDIATELY
!   BOUND target_z, STARTING AT N=1
! IF continue_search.GT.1, SCAN STARTS FROM VALUE OF N FROM PREVIOUS
!   CALL OF LIR
subroutine lir(target_z,table_z,result_y,table_y,num_y,y_stride, &
     num_points,continue_search,interp_flag)
      implicit none
      double precision, intent(in) :: target_z
      double precision, intent(in) :: table_z(*)
      double precision, intent(out) :: result_y(*)
      double precision, intent(in) :: table_y(*)
      integer, intent(in) :: num_y, y_stride, num_points, continue_search
      integer, intent(out) :: interp_flag
      call lir_impl(0, target_z, table_z, result_y, table_y, num_y, &
           y_stride, num_points, continue_search, interp_flag)
end subroutine lir

! Implementation behind lir (restructured 2026: the ENTRY statement
! became a wrapper plus this routine, and the label 2-9 search web
! became a structured loop; comparisons and arithmetic are unchanged).
! search_idx is SAVEd on purpose: it is the last-index interpolation
! memory that continue_search > 1 resumes from.
subroutine lir_impl(linear_mode_in, target_z,table_z,result_y,table_y, &
     num_y,y_stride,num_points,continue_search,interp_flag)

      implicit none

      integer, intent(in) :: linear_mode_in
      double precision, intent(in) :: target_z
      double precision, intent(in) :: table_z(*)
      double precision, intent(out) :: result_y(*)
      double precision, intent(in) :: table_y(*)
      integer, intent(in) :: num_y, y_stride, num_points, continue_search
      integer, intent(out) :: interp_flag

      double precision :: weight(4)
      integer :: search_idx
      data search_idx/-1/
      save   ! INTENTIONAL: last-index interpolation memory; byte-pinned by Stage-0

      integer :: linear_mode
      integer :: stride, stride_m1, y_strided, num_y_strided, table_end
      double precision :: diff
      integer :: pivot, closest
      double precision :: y1, y2, y3, y4, z1, z2, z3, z4, z12, z34
      integer :: y_idx, j, k, base_idx
      double precision :: yy

      linear_mode=linear_mode_in
      stride=1
! CHECK NT AND RESET IL IF NECESSARY
      if (num_points.lt.2) then
         return
      end if
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
      do
         if(search_idx.gt.table_end) then
            interp_flag=0
            exit
         end if
         if (diff .eq. 0.0) then
            return
         end if
         if (diff .lt. 0.0) then
            if (table_z(search_idx) .lt. target_z) exit
         else
            if (table_z(search_idx) .gt. target_z) exit
         end if
         if (table_z(search_idx) .eq. target_z) then
! SET Y WHEN Z LIES ON A MESH POINT
            base_idx=(search_idx-1)*y_stride
            do y_idx=1,num_y_strided
               result_y(y_idx)=table_y(y_idx+base_idx)
               if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
            end do
            search_idx=(search_idx+stride-1)/stride
            return
         end if
         search_idx=search_idx+stride
      end do
! CONTROL WHEN Z DOES NOT LIE ON A MESH POINT
      if(search_idx.le.1) interp_flag=0
      if (linear_mode.ne.1) then
! CUBIC INTERPOLATION/EXTRAPOLATION
! PIVOTAL POINT (M) AND POINT (K) CLOSEST TO Z
      pivot=search_idx
      closest=3
      if (search_idx.le.1+stride) then
      pivot=1+stride+stride
      closest=search_idx
      end if
      if (search_idx.ge.table_end) then
      pivot=table_end-stride
      closest=4
      end if
! WEIGHTING FACTORS
      y1=table_z(pivot-stride*2)
      y2=table_z(pivot-stride)
      y3=table_z(pivot)
      y4=table_z(pivot+stride)
      z1=target_z-y1
      z2=target_z-y2
      z3=target_z-y3
      z4=target_z-y4
      z12=z1*z2
      z34=z3*z4
      weight(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      weight(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      weight(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      weight(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
! CORRECT A(K)
      diff=weight(1)+weight(2)+weight(3)+weight(4)
      weight(closest)=(1.d0+weight(closest))-diff
! COMPUTE Y
      pivot=(pivot-1)/stride-3
      pivot=pivot*y_strided
      do y_idx=1,num_y_strided
         k=y_idx+pivot
         yy=0.d0
         do j=1,4
            k=k+y_strided
            diff=table_y(k)
            yy=yy+weight(j)*diff
         end do
         result_y(y_idx)=yy
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
      end do
      search_idx=(search_idx+stride-1)/stride
      return
! LINEAR INTERPOLATION/EXTRAPOLATION
      end if
      if(search_idx.eq.1) search_idx=1+stride
      if(search_idx.gt.table_end) search_idx=table_end
      z1=table_z(search_idx)
      y1=(z1-target_z)/(z1-table_z(search_idx-stride))
      y2=1.0d0-y1
      base_idx=(search_idx-1)*y_stride
      pivot=base_idx-y_strided
      do y_idx=1,num_y_strided,stride
         result_y(y_idx)=y1*table_y(y_idx+pivot)+y2*table_y(y_idx+base_idx)
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
      end do
! RESET N
      search_idx=(search_idx+stride-1)/stride
      return
end subroutine lir_impl

!----------------------------------------------------------------------
! ratext
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ratext.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Moved here (2026) from net_lib.f90: a diagonal rational-function
! extrapolator, the textbook companion to this module's own bsstep/
! mmid Bulirsch-Stoer stepper above (same algorithm as subroutine
! RZEXTR in Numerical Recipes, p.566) -- bsstep already called it via
! an errant `use net_lib`. Also called directly by
! net_lib.f90's own liburn to extrapolate a sequence of sub-stepped
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
      save   ! INTENTIONAL: rational-extrapolation state; byte-pinned by Stage-0

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

!----------------------------------------------------------------------
! simeqc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original simeqc.f; only variable names, source form, and comment
! style were updated.
!
! Gauss-Jordan elimination with partial pivoting, operating on a
! system matrix stored as a flat array: num_unknowns equations (rows),
! num_cols columns (num_cols > num_unknowns for one or more augmented
! right-hand-side columns), stored column-major in system_matrix(56).
subroutine simeqc(system_matrix, num_cols, num_unknowns, ierr)

      implicit none

      double precision, intent(inout) :: system_matrix(56)
      integer, intent(in) :: num_cols, num_unknowns
      integer :: jj, j, jy, it, i, ij, imax, ia, ib, ic, id, ix, jx, ny, &
           n1, ig, ih
      double precision :: biga, swap_val

      integer, intent(out) :: ierr

      ierr = 0

      jj=-num_unknowns
      do j=1,num_unknowns
      jy=j+1
      jj=jj+num_unknowns+1
      biga=0.0d0
      it=jj-j
      do i=j,num_unknowns
      ij=it+i
      if(dabs(biga).ge.dabs(system_matrix(ij))) cycle
      biga=system_matrix(ij)
      imax=i
      end do
      if (dabs(biga).eq.0.0d0) then
      write (5,1011)
 1011 format (1x,'STOPPED AT 1010')
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
      ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
      ! preserve the historical stop on a nonzero return.
      ierr = 1
      return
      end if
      ia=j+num_unknowns*(j-2)
      it=imax-j
      do i=j,num_cols
      ia=ia+num_unknowns
      ib=ia+it
      swap_val=system_matrix(ia)
      system_matrix(ia)=system_matrix(ib)
      system_matrix(ib)=swap_val
      system_matrix(ia)=system_matrix(ia)/biga
      end do
      if(j.eq.num_unknowns) exit
      ia=num_unknowns*(j-1)
      do ix=jy,num_unknowns
      ib=ia+ix
      it=j-ix
      do jx=jy,num_cols
      ic=num_unknowns*(jx-1)+ix
      id=ic+it
      system_matrix(ic)=system_matrix(ic)-system_matrix(ib)*system_matrix(id)
      end do
      end do
      end do
      ny=num_unknowns-1
      it=num_unknowns*num_unknowns
      do j=1,ny
      ia=it-j
      ic=num_unknowns*num_cols
      ib=ic-j
      do i=1,j
      system_matrix(ib)=system_matrix(ib)-system_matrix(ia)*system_matrix(ic)
      n1=num_cols-1
      ig=ib
      ih=ic
      do
      if(n1.le.num_unknowns) exit
      ig=ig-num_unknowns
      ih=ih-num_unknowns
      system_matrix(ig)=system_matrix(ig)-system_matrix(ia)*system_matrix(ih)
      n1=n1-1
      end do
      ia=ia-num_unknowns
      ic=ic-1
      end do
      end do
      return
end subroutine simeqc

!----------------------------------------------------------------------
! slopes
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original slopes.f; only variable names, source form, and comment
! style were updated.
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C.STATE UNIVERSITY
!
! SLOPES CALCULATES THE DERIVATIVE AT EACH OF THE DATA POINTS. THE
! SLOPES PROVIDED WILL INSURE THAT AN OSCULATORY QUADRATIC SPLINE WILL
! HAVE ONE ADDITIONAL KNOT BETWEEN TWO ADJACENT POINTS OF INTERPOLATION.
! CONVEXITY AND MONOTONICITY ARE PRESERVED WHEREVER THESE CONDITIONS
! ARE COMPATIBLE WITH THE DATA.
!
! ON INPUT--
!
!   table_x CONTAINS THE ABSCISSAS OF THE DATA POINTS.
!
!   table_y CONTAINS THE ORDINATES OF THE DATA POINTS.
!
!   num_points IS THE NUMBER OF DATA POINTS (DIMENSION OF table_x,
!   table_y).
!
!
! ON OUTPUT--
!
!   first_derivs CONTAINS THE VALUE OF THE FIRST DERIVATIVE AT EACH
!   DATA POINT.
!
! AND
!
!   SLOPES DOES NOT ALTER table_x,table_y,num_points.
!
! NOTE: like search.f, the original slopes.f has no blanket IMPLICIT
! REAL*8(A-H,O-Z) statement -- it relies on default Fortran implicit
! typing (I-N integer, else real) except where table_x/table_y/
! first_derivs and the M1/M2/... scratch variables are explicitly
! declared REAL*8. Types below match that original default typing.
!
!----------------------------------------------------------------------
!
subroutine slopes(table_x, table_y, first_derivs, num_points)
      use star_info_lib, only: json

      implicit none

      double precision, intent(in) :: table_x(json), table_y(json)
      double precision, intent(out) :: first_derivs(json)
      integer, intent(in) :: num_points
      integer :: num_points_m1, prev_idx, idx, next_idx
      double precision :: slope1, slope2, x_bar, x_hat, y_diff1, y_diff2, &
           y_x_mid, x_mid, slope1_saved, slope2_saved

      num_points_m1= num_points - 1
      prev_idx=1
      idx=2
      next_idx=3
!
! CALCULATE THE SLOPES OF THE TWO LINES JOINING THE FIRST THREE DATA
! POINTS
      y_diff1=table_y(2) - table_y(1)
      y_diff2=table_y(3) - table_y(2)
      slope1=y_diff1/(table_x(2) - table_x(1))
      slope1_saved=slope1
      slope2=y_diff2/(table_x(3)-table_x(2))
      slope2_saved=slope2
!
! MAIN LOOP OVER THE INTERIOR POINTS. (Restructured 2026 from the
! original goto flow at labels 10-50; arithmetic is unchanged.)
!
! IF ONE OF THE PRECEDING SLOPES IS ZERO OR IF THEY HAVE OPPOSITE SIGN,
! ASSIGN THE VALUE ZERO TO THE DERIVATIVE AT THE MIDDLE POINT.
      do
      if (slope1.eq.0.d0 .or. slope2.eq.0.d0 .or. (slope1*slope2).le.0.d0) then
         first_derivs(idx)= 0.d0
      else if (abs(slope1) .gt. abs(slope2)) then
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M1.
         x_bar=(y_diff2/slope1) + table_x(idx)
         x_hat= (x_bar + table_x(next_idx))/2.d0
         first_derivs(idx)=y_diff2/(x_hat - table_x(idx))
      else
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M2.
         x_bar=(-y_diff1/slope2) + table_x(idx)
         x_hat=(table_x(prev_idx) + x_bar)/2.d0
         first_derivs(idx)=y_diff1/(table_x(idx) - x_hat)
      end if
!
! INCREMENT COUNTERS
      prev_idx=idx
      idx=next_idx
      next_idx=next_idx+1
      if (idx .gt. num_points_m1) exit
!
! CALCULATE THE SLOPES OF THE TWO LINES JOINING THREE CONSECUTIVE DATA
! POINTS.
      y_diff1=table_y(idx) - table_y(prev_idx)
      y_diff2=table_y(next_idx) - table_y(idx)
      slope1=y_diff1/(table_x(idx) - table_x(prev_idx))
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!       M2=YDIF2/(XTAB(I1) - XTAB(I))
      call safedivide(y_diff2, (table_x(next_idx) - table_x(idx)), slope2)
      end do
!
! CALCULATE THE SLOPE AT THE LAST POINT, XTAB(NUM).
      if ((slope1*slope2) .lt. 0.d0) then
         first_derivs(num_points)=2.d0*slope2
      else
         x_mid= (table_x(num_points_m1)+table_x(num_points))/2.d0
         y_x_mid=first_derivs(num_points_m1)*(x_mid - table_x(num_points_m1)) + &
              table_y(num_points_m1)
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!       MTAB(NUM)=(YTAB(NUM)-YXMID)/(XTAB(NUM)-XMID)
         call safedivide((table_y(num_points)-y_x_mid), &
              (table_x(num_points)-x_mid), first_derivs(num_points))
         if ((first_derivs(num_points)*slope2) .lt. 0.d0) then
            first_derivs(num_points)=0.d0
         end if
      end if
!
! CALCULATE THE SLOPE AT THE FIRST POINT, XTAB(1).
      if ((slope1_saved*slope2_saved) .lt. 0.d0) then
         first_derivs(1)=2.d0*slope1_saved
         return
      end if
      x_mid=(table_x(1) + table_x(2))/2.d0
      y_x_mid=first_derivs(2)*(x_mid - table_x(2)) + table_y(2)
      first_derivs(1)=(y_x_mid - table_y(1))/(x_mid - table_x(1))
      if ((first_derivs(1) * slope1_saved) .lt. 0.d0) then
         first_derivs(1)=0.d0
      end if
      return
!
end subroutine slopes

!----------------------------------------------------------------------
! KC 2025-05-31 SAFEDIVIDE
!----------------------------------------------------------------------
subroutine safedivide(numerator, denominator, quotient)
      implicit none
      double precision, intent(in) :: numerator, denominator
      double precision, intent(out) :: quotient

      quotient = 0.d0
      if (numerator .ne. 0.d0 .and. denominator .ne. 0.d0) then
         quotient = numerator / denominator
      end if

      return
end subroutine safedivide

!----------------------------------------------------------------------
! cases
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original cases.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
! Computes the knots and other parameters of the spline on the
! interval (x_left,x_right).
!
! ON INPUT--
!
!   (x_left,y_left) AND (x_right,y_right) ARE THE COORDINATES OF THE
!   POINTS OF INTERPOLATION.
!
!   slope_left IS THE SLOPE AT (x_left,y_left).
!
!   slope_right IS THE SLOPE AT (x_right,y_right)
!
!   spline_case CONTROLS THE NUMBER AND LOCATION OF THE KNOTS.
!
!
! ON OUTPUT--
!
!   (knot_v_x,knot_v_y),(knot_w_x,knot_w_y),(knot_z_x,knot_z_y), AND
!   (knot_e_x,knot_e_y) ARE THE COORDINATES OF THE KNOTS AND OTHER
!   PARAMETERS OF THE SPLINE ON (x_left,x_right). (knot_e_x,knot_e_y)
!   AND (knot_y_x,knot_y_y) ARE USED ONLY IF spline_case=4.
!
! AND
!
!   CASES DOES NOT ALTER x_left,y_left,slope_left,slope_right,x_right,
!   y_right.
!
!----------------------------------------------------------------------
subroutine cases(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, knot_e_x, knot_e_y, knot_v_x, knot_v_y, knot_w_x, &
     knot_w_y, knot_z_x, knot_z_y, knot_y_x, knot_y_y, spline_case)

      double precision :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right, knot_v_x, knot_v_y, knot_z_x, knot_z_y, &
           knot_w_x, knot_w_y, knot_e_x, knot_e_y, &
           mbar1, mbar2, mbar3, c1, d1, h1, j1, knot_y_x, knot_y_y, &
           k1, ztwo
      integer :: spline_case
      if (spline_case .ne. 3 .and. spline_case .ne. 4) then
      if (spline_case .ne. 2) then

! CALCULATE THE PARAMETERS FOR CASE 1.
      knot_z_x=(y_left-y_right+slope_right*x_right-slope_left*x_left)/ &
           (slope_right-slope_left)
      ztwo=y_left+slope_left*(knot_z_x-x_left)
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=(y_left+ztwo)/2.d0
      knot_w_x=(knot_z_x+x_right)/2.d0
      knot_w_y=(ztwo+y_right)/2.d0
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
           (knot_z_x-knot_v_x)
      return

! CALCULATE THE PARAMETERS FOR CASE 2.
      end if
      knot_z_x=(x_left+x_right)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(knot_z_x+x_right)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=(knot_v_y+knot_w_y)/2.d0
      return

! CALCULATE THE PARAMETERS USED IN BOTH CASES 3 AND 4.
      end if
      c1=x_left+(y_right-y_left)/slope_left
      d1=x_right+(y_left-y_right)/slope_right
      h1=2.d0*c1-x_left
      j1=2.d0*d1-x_right
      mbar1=(y_right-y_left)/(h1-x_left)
      mbar2=(y_left-y_right)/(j1-x_right)

      if (spline_case .ne. 4) then

! CALCULATE THE PARAMETERS FOR CASE 3.
      k1=(y_left-y_right+x_right*mbar2-x_left*mbar1)/(mbar2-mbar1)
      if (abs(slope_left) .le. abs(slope_right)) then
      knot_z_x=(k1+x_right)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
      (knot_z_x-knot_v_x)
      return
      end if
      knot_z_x=(k1+x_left)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
           (knot_z_x-knot_v_x)
      return

! CALCULATE THE PARAMETERS FOR CASE 4.
      end if
      knot_y_x=(x_left+c1)/2.d0
      knot_v_x=(x_left+knot_y_x)/2.d0
      knot_v_y=slope_left*(knot_v_x-x_left) + y_left
      knot_z_x=(d1+x_right)/2.d0
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=slope_right*(knot_w_x-x_right) + y_right
      mbar3=(knot_w_y-knot_v_y)/(knot_w_x-knot_v_x)
      knot_y_y=mbar3*(knot_y_x-knot_v_x) + knot_v_y
      knot_z_y=mbar3*(knot_z_x-knot_v_x) + knot_v_y
      knot_e_x=(knot_y_x+knot_z_x)/2.d0
      knot_e_y=mbar3*(knot_e_x-knot_v_x) + knot_v_y
      return

end subroutine cases

!----------------------------------------------------------------------
! choose
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original choose.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A. ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!
! CHOOSE DETERMINES THE CASE NEEDED FOR THE COMPUTATION OF THE PARAME-
! TERS OF THE QUADRATIC SPLINE AND RETURNS THE VALUE IN THE VARIABLE
! spline_case.
!
! ON INPUT--
!
!   (x_left,y_left) GIVES THE COORDINATES OF ONE OF THE POINTS OF
!   INTERPOLATION.
!
!   slope_left SPECIFIES THE DERIVATIVE CONDITION AT (x_left,y_left).
!
!   (x_right,y_right) GIVES THE COORDINATES OF ONE OF THE POINTS OF
!   INTERPOLATION.
!
!   slope_right SPECIFIES THE DERIVATIVE CONDITION AT (x_right,y_right).
!
!   eps_tol IS AN ERROR TOLERANCE USED TO DISTINGUISH CASES WHEN
!   slope_left OR slope_right IS RELATIVELY CLOSE TO THE SLOPE OR TWICE
!   THE SLOPE OF THE LINE SEGMENT JOINING (x_left,y_left) AND
!   (x_right,y_right). IF eps_tol IS NOT EQUAL TO ZERO, THEN eps_tol
!   SHOULD BE GREATER THAN OR EQUAL TO MACHINE EPSILON.
!
!
! ON OUTPUT--
!
!   spline_case CONTAINS THE VALUE WHICH CONTROLS HOW THE PARAMETERS OF
!   THE QUADRATIC SPLINE ARE EVALUATED.
!
! AND
!
!   CHOOSE DOES NOT ALTER x_left,y_left,x_right,y_right,slope_left,
!   slope_right,eps_tol.
!
!----------------------------------------------------------------------
subroutine choose(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, eps_tol, spline_case)

      double precision :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right, mref, mref1, mref2, spq, prod, &
           prod1, prod2, eps_tol
      integer :: spline_case

! CALCULATE THE SLOPE spq OF THE LINE JOINING (x_left,y_left),(x_right,y_right).
      spq=(y_right-y_left)/(x_right-x_left)

! CHECK WHETHER OR NOT spq IS 0.
! ******MODIFICATION BY MARC PINSONNEAULT TO AVOID DIVISION BY ZERO
! ******IN SR CASES
      if (slope_left.eq.0.0d0 .or. slope_right.eq.0.0d0) then
         spline_case=2
         return
      end if
! ******
      if (spq .eq. 0.d0) then
         if ((slope_left*slope_right) .ge. 0.d0) then
            spline_case=2
         else
            spline_case=1
         end if
         return
      end if

      prod1=spq*slope_left
      prod2=spq*slope_right

! FIND THE ABSOLUTE VALUES OF THE SLOPES spq,slope_left,AND slope_right.
      mref=abs(spq)
      mref1=abs(slope_left)
      mref2=abs(slope_right)

! THE SIGN OF AT LEAST ONE OF THE SLOPES slope_left,slope_right DOES NOT
! AGREE WITH THE SIGN OF THE SLOPE spq.
      if ((prod1 .lt. 0.d0) .or. (prod2 .lt. 0.d0)) then
         if ((prod1 .lt. 0.d0) .and. (prod2 .lt. 0.d0)) then
            spline_case=2
            return
         end if
         if (prod1 .ge. 0.d0) then
            if (mref1 .gt. ((1.d0+eps_tol)*mref)) then
               spline_case=1
            else
               spline_case=2
            end if
            return
         end if
         if (mref2 .gt. ((1.d0+eps_tol)*mref)) then
            spline_case=1
         else
            spline_case=2
         end if
         return
      end if

! IF THE RELATIVE DEVIATION OF slope_left OR slope_right FROM spq IS LESS THAN
! eps_tol, THEN CHOOSE CASE 2 OR CASE 3.
      if (abs(spq-slope_left).gt.eps_tol*mref .and. abs(spq-slope_right).gt.eps_tol*mref) then
         prod=(mref-mref1)*(mref-mref2)
         if (prod .lt. 0.d0) then

! L1, THE LINE THROUGH (x_left,y_left) WITH SLOPE slope_left, AND L2, THE LINE
! THROUGH (x_right,y_right) WITH SLOPE slope_right, INTERSECT AT A POINT WHOSE
! ABSCISSA IS BETWEEN x_left AND x_right. THE ABSCISSA BECOMES A KNOT OF THE
! SPLINE.
            spline_case=1
            return
         end if
      end if

! IN CASES 3 AND 4, SIGN(slope_left)=SIGN(slope_right)=SIGN(spq).
! CHOOSE CASE 4 IF THE OTHER SLOPE IS GREATER THAN (2.-eps_tol)*mref
! (NEITHER L1 NOR L2 CROSSES THE MIDLINE: TWO KNOTS); OTHERWISE CASE 3
! (EITHER L1 OR L2 CROSSES THE MIDLINE, BUT NOT BOTH).
      if (mref1 .gt. (2.d0*mref)) then
         if (mref2 .gt. (2.d0-eps_tol)*mref) then
            spline_case=4
         else
            spline_case=3
         end if
         return
      end if
      if (mref2 .gt. (2.d0*mref)) then
         if (mref1 .gt. (2.d0-eps_tol)*mref) then
            spline_case=4
         else
            spline_case=3
         end if
         return
      end if

! BOTH L1 AND L2 CROSS THE LINE THROUGH (x_left+x_right/2.,y_left) AND
! (x_left+x_right/2.,y_right), WHICH IS THE MIDLINE OF THE RECTANGLE FORMED
! BY (x_left,y_left),(x_right,y_left),(x_right,y_right), AND (x_left,y_right),
! OR BOTH slope_left AND slope_right HAVE SIGNS DIFFERENT THAN THE SIGN OF
! spq, OR ONE OF slope_left AND slope_right HAS OPPOSITE SIGN FROM spq AND L1
! AND L2 INTERSECT TO THE LEFT OF x_left OR TO THE RIGHT OF x_right. THE
! POINT (x_left+x_right)/2. IS A KNOT OF THE SPLINE.
! (Restructured 2026 from the original goto decision tree; two
! unreachable branches at old labels 60/110 were dropped.)
      spline_case=2
      return

end subroutine choose

!----------------------------------------------------------------------
! search
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original search.f; only variable names, source form, and comment
! style were updated.
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
! SEARCH CONDUCTS A BINARY SEARCH FOR eval_point. SEARCH IS CALLED ONLY
! IF eval_point IS BETWEEN table_x(1) AND table_x(num_table_points).
!
! ON INPUT--
!
!   table_x CONTAINS THE ABSCISSAS OF THE DATA POINTS OF INTERPOLATION.
!
!   num_table_points IS THE DIMENSION OF table_x
!
!   eval_point IS THE VALUE WHOSE RELATIVE POSITION IN table_x IS
!   LOCATED BY SEARCH.
!
!
! ON OUTPUT--
!
!   found_flag IS SET EQUAL TO 1 IF eval_point IS FOUND IN table_x AND
!   IS SET EQUAL TO 0 OTHERWISE.
!
!   table_idx IS THE INDEX OF THE LARGEST VALUE IN table_x FOR WHICH
!   table_x(I) .LT. eval_point.
!
! AND
!
!   SEARCH DOES NOT ALTER table_x,num_table_points,eval_point.
!
! NOTE: unlike most YREC files, the original search.f has no blanket
! IMPLICIT REAL*8(A-H,O-Z) statement -- it relies on default Fortran
! implicit typing (I-N integer, else real) except where XTAB/S and
! FND/FIRST are explicitly declared. Types below are chosen to match
! that original default typing exactly (NUM/LCN/MIDDLE/LAST all fall
! in the I-N default-integer range).
!
subroutine search(table_x, num_table_points, eval_point, table_idx, &
     found_flag)
      use star_info_lib, only: json

      implicit none

      double precision, intent(in) :: table_x(json), eval_point
      integer, intent(in) :: num_table_points
      integer, intent(out) :: table_idx, found_flag
      integer :: first_idx, last_idx, middle_idx

      first_idx=1
      last_idx=num_table_points
      found_flag=0
!
! (Restructured 2026 from the original goto binary search at labels
! 10-60; comparisons are unchanged.)
      if (table_x(1) .eq. eval_point) then
         table_idx=1
         found_flag=1
         return
      end if
      if (table_x(num_table_points) .eq. eval_point) then
         table_idx=num_table_points
         found_flag=1
         return
      end if
!
      do
!
! IF (LAST-FIRST) .EQ. 1, S IS NOT IN XTAB.  SET POSITION EQUAL TO
! FIRST.
      if ((last_idx-first_idx) .eq. 1) then
         table_idx=first_idx
         return
      end if
!
      middle_idx=(first_idx+last_idx)/2
!
! CHECK IF S .EQ. XTAB(MIDDLE). IF NOT, CONTINUE THE SEARCH IN THE
! APPROPRIATE HALF OF THE VECTOR XTAB.
      if (table_x(middle_idx) .lt. eval_point) then
         first_idx=middle_idx
      else if (table_x(middle_idx) .eq. eval_point) then
         table_idx=middle_idx
         found_flag =1
         return
      else
         last_idx=middle_idx
      end if
      end do
end subroutine search

!----------------------------------------------------------------------
! spline
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original spline.f; only variable names, source form, and comment
! style were updated.
!
!                                 SHAPE PRESERVING QUDRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
! SPLINE FINDS THE IMAGE OF A POINT IN eval_point.
!
! ON INPUT--
!
!   eval_point CONTAINS THE VALUE AT WHICH THE SPLINE IS EVALUATED.
!
!   (x_left,y_left) ARE THE COORDINATES OF THE LEFT-HAND DATA POINT
!   USED IN THE EVALUATION OF eval_point.
!
!   (x_right,y_right) ARE THE COORDINATES OF THE RIGHT-HAND DATA POINT
!   USED IN THE EVALUATION OF eval_point.
!
!   z1,z2,y1,y2,e2,w2,v2 ARE THE PARAMETERS OF THE SPLINE.
!
!   spline_case CONTROLS THE EVALUATION OF THE SPLINE BY INDICATING
!   WHETHER ONE OR TWO KNOTS WERE PLACED IN THE INTERVAL
!   (x_left,x_right).
!
!
! ON OUTPUT--
!
!   SPLINE IS THE IMAGE OF eval_point.
!
! AND
!
!   SPLINE DOES NOT ALTER ANY OF THE INPUT PARAMETERS.
!
!----------------------------------------------------------------------
!
!  *****MODIFICATION DUE TO MARC PINSONNEAULT 6/87*****
!  IF DIVISION BY ZERO WOULD BE CAUSED,LINEAR INTERPOLATION IS USED
!  INSTEAD OF THE SPLINE.
! IF NCASE .EQ. 4, MORE THAN ONE KNOT WAS PLACED IN THE INTERVAL.
function spline(eval_point, z1, z2, x_left, y_left, x_right, y_right, &
     y1, y2, e2, w2, v2, spline_case)

      implicit none
      double precision :: spline
      double precision, intent(in) :: eval_point, z1, z2, x_left, y_left, &
           x_right, y_right, y1, y2, e2, w2, v2
      integer, intent(in) :: spline_case
      double precision :: linear_interp_frac

! (Restructured 2026 from the original arithmetic-IF goto fans at
! labels 10-100; arithmetic is unchanged.)
      if (spline_case .ne. 4) then
!
! CASES 1,2, OR 3.
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE KNOT.
      if (z1 .lt. eval_point) then
         if(x_right.ne.z1)then
         spline=(z2*(x_right-eval_point)**2+w2*2.d0*(eval_point-z1)*(x_right-eval_point) &
                 +y_right*(eval_point-z1)**2)/(x_right-z1)**2
         else
          linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
          spline = y_left + linear_interp_frac*(y_right - y_left)
         end if
      else if (z1 .eq. eval_point) then
         spline=z2
      else
         if(z1.ne.x_left)then
         spline=(y_left*(z1-eval_point)**2+v2*2.d0*(eval_point-x_left)*(z1-eval_point)+ &
                 z2*(eval_point-x_left)**2)/(z1-x_left)**2
         else
          linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
          spline = y_left + linear_interp_frac*(y_right - y_left)
         end if
      end if
      return
      end if
!
! CASE 4.
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE FIRST KNOT.
      if (y1 .lt. eval_point) then
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE SECOND KNOT.
         if (z1 .lt. eval_point) then
            if(x_right.ne.z1)then
            spline=(z2*(x_right-eval_point)**2+w2*2.d0*(eval_point-z1)*(x_right-eval_point) &
                    +y_right*(eval_point-z1)**2)/(x_right-z1)**2
            else
             linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
             spline = y_left + linear_interp_frac*(y_right - y_left)
            end if
         else if (z1 .eq. eval_point) then
            spline=z2
         else
            if(z1.ne.y1)then
            spline=(y2*(z1-eval_point)**2+e2*2.d0*(eval_point-y1)*(z1-eval_point)+z2*(eval_point &
                    -y1)**2)/(z1-y1)**2
            else
             linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
             spline = y_left + linear_interp_frac*(y_right - y_left)
            end if
         end if
      else if (y1 .eq. eval_point) then
         spline=y2
      else
         if(y1.ne.x_left)then
         spline=(y_left*(y1-eval_point)**2+v2*2.d0*(eval_point-x_left)*(y1-eval_point)+ &
                 y2*(eval_point-x_left)**2)/(y1-x_left)**2
         else
          linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
          spline = y_left + linear_interp_frac*(y_right - y_left)
         end if
      end if
      return
end function spline

!----------------------------------------------------------------------
! meval
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original meval.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Controls the evaluation of an osculatory (shape-preserving)
! quadratic spline, by D.F. McAllister & J.A. Roulier, coded by
! S.L. Dodd & M. Roulier, N.C. State University. The caller may
! provide slopes at the points of interpolation directly, or use
! subroutine SLOPES to compute slopes consistent with the shape of
! the data.
!
! On input --
!   eval_x must be a nondecreasing vector of points at which the
!   spline will be evaluated.
!   table_x contains the abscissas of the data points to be
!   interpolated. table_x must be increasing.
!   table_y contains the ordinates of the data points to be
!   interpolated.
!   table_slope contains the slope of the spline at each point of
!   interpolation.
!   num_table_points is the number of data points (dimension of
!   table_x and table_y).
!   num_eval_points is the number of points of evaluation (dimension
!   of eval_x and eval_y).
!   eps_tol is a relative error tolerance used in subroutine CHOOSE to
!   distinguish the situation table_slope(i) or table_slope(i+1) is
!   relatively close to the slope or twice the slope of the linear
!   segment between table_x(i) and table_x(i+1). If this situation
!   occurs, roundoff may cause a change in convexity or monotonicity
!   of the resulting spline and a change in the case number provided
!   by CHOOSE. If eps_tol is not equal to zero, then eps_tol should be
!   greater than or equal to machine epsilon.
!
! On output --
!   eval_y contains the images of the points in eval_x.
!   err_code is an error code --
!   err_code=0 - meval ran normally.
!   err_code=1 - eval_x(i) is less than table_x(1) for at least one i,
!                or eval_x(i) is greater than table_x(num_table_points)
!                for at least one i. meval will extrapolate to provide
!                function values for these abscissas.
!   err_code=2 - eval_x(i+1) .lt. eval_x(i) for some i.
!
! meval does not alter eval_x, table_x, table_y, table_slope,
! num_table_points, num_eval_points.
!
! meval calls the following subroutines or functions:
!    search
!    cases
!    choose
!    spline
!----------------------------------------------------------------------
subroutine meval(eval_x, eval_y, table_x, table_y, table_slope, &
     num_table_points, num_eval_points, eps_tol, err_code)

      implicit none

      integer, parameter :: max_points = 5000

      double precision, intent(in) :: eval_x(max_points), &
           table_x(max_points), table_y(max_points), &
           table_slope(max_points), eps_tol
      double precision, intent(out) :: eval_y(max_points)
      integer, intent(in) :: num_table_points, num_eval_points
      integer, intent(out) :: err_code

      double precision :: spline_v1, spline_v2, spline_w1, spline_w2, &
           spline_z1, spline_z2, spline_y1, spline_y2, spline_e1, spline_e2
      integer :: start_idx, start_idx1, end_idx, end_idx1, found_flag
      integer :: i, ind, loop_bound, table_idx, table_idx1, spline_case, &
           num_minus1
      logical :: recompute_tail_params
      start_idx = 1
      end_idx = num_eval_points
      err_code = 0
      if (num_eval_points .ne. 1) then

! Determine if eval_x is nondecreasing.
      loop_bound = num_eval_points - 1
      do i = 1, loop_bound
        if (eval_x(i+1) .ge. eval_x(i)) cycle
        err_code = 2
        return
      end do

! If eval_x(i) .lt. table_x(1), then eval_x(i)=table_y(1).
! If eval_x(i) .gt. table_x(num_table_points), then
! eval_x(i)=table_y(num_table_points).
!
! Determine if any of the points in eval_x are less than the abscissa
! of the first data point.
      end if
      do i = 1, num_eval_points
        if (eval_x(i) .ge. table_x(1)) exit
        start_idx = i + 1
      end do

      loop_bound = num_eval_points + 1

! Determine if any of the points in eval_x are greater than the
! abscissa of the last data point.
      do i = 1, num_eval_points
        ind = loop_bound - i
        if (eval_x(ind) .le. table_x(num_table_points)) exit
        end_idx = ind - 1
      end do

! Calculate the images of points of evaluation whose abscissas are
! less than the abscissa of the first data point.
      if (start_idx .ne. 1) then
! Set the error parameter to indicate that extrapolation has occurred.
      err_code = 1
      call choose(table_x(1), table_y(1), table_slope(1), table_slope(2), &
           table_x(2), table_y(2), eps_tol, spline_case)
      call cases(table_x(1), table_y(1), table_slope(1), table_slope(2), &
           table_x(2), table_y(2), spline_e1, spline_e2, spline_v1, &
           spline_v2, spline_w1, spline_w2, spline_z1, spline_z2, &
           spline_y1, spline_y2, spline_case)
      start_idx1 = start_idx - 1
      do i = 1, start_idx1
       eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, table_x(1), &
            table_y(1), table_x(2), table_y(2), spline_y1, spline_y2, &
            spline_e2, spline_w2, spline_v2, spline_case)
      end do
      if (num_eval_points .eq. 1) then
         return
      end if
      end if

! search locates the interval in which the first in-range point of
! evaluation lies. The single-out-of-range-high case skips straight
! to the high-extrapolation tail.
      recompute_tail_params = .true.
      if (num_eval_points .ne. 1 .or. end_idx .eq. num_eval_points) then
      call search(table_x, num_table_points, eval_x(start_idx), &
           table_idx, found_flag)

      table_idx1 = table_idx + 1

! If the first in-range point of evaluation is equal to one of the
! data points, assign the appropriate value from table_y. Continue
! until a point of evaluation is found which is not equal to a data
! point.
      if (found_flag .ne. 0) then
      do
         eval_y(start_idx) = table_y(table_idx)
      start_idx1 = start_idx
      start_idx = start_idx + 1
      if (start_idx .gt. num_eval_points) then
         return
      end if
      if (eval_x(start_idx1) .ne. eval_x(start_idx)) exit
      end do

! Advance the table pointers until the next point of evaluation lies
! inside the current table interval, assigning images to any points
! of evaluation that coincide with data points along the way.
! (Restructured 2026 from the original arithmetic-IF web at labels
! 100/110/120; comparisons and order of assignment are unchanged.)
      do
         if (eval_x(start_idx) .lt. table_x(table_idx1)) exit
         if (eval_x(start_idx) .eq. table_x(table_idx1)) then
            do
               eval_y(start_idx) = table_y(table_idx1)
               start_idx1 = start_idx
               start_idx = start_idx + 1
               if (start_idx .gt. num_eval_points) then
                  return
               end if
               if (eval_x(start_idx) .ne. eval_x(start_idx1)) exit
            end do
         end if
         table_idx = table_idx1
         table_idx1 = table_idx1 + 1
      end do

! Calculate the images of all the points which lie within range of
! the data.
      end if
      if (table_idx .ne. 1 .or. err_code .ne. 1) then
      call choose(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), eps_tol, spline_case)
      call cases(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), spline_e1, spline_e2, &
           spline_v1, spline_v2, spline_w1, spline_w2, spline_z1, &
           spline_z2, spline_y1, spline_y2, spline_case)
      end if

      do i = start_idx, end_idx

! If eval_x(i) is beyond the current table interval, advance the
! pointers (and recompute the spline parameters) first; if it is a
! data point, its image is known; otherwise evaluate the spline.
! (Restructured 2026 from the original arithmetic-IF web at labels
! 150/160/170/180; comparisons and evaluation order are unchanged.)
      if (eval_x(i) .gt. table_x(table_idx1)) then
         do
            table_idx = table_idx1
            table_idx1 = table_idx + 1
            if (eval_x(i) .le. table_x(table_idx1)) exit
         end do
         if (eval_x(i) .lt. table_x(table_idx1)) then
! Call choose to determine the appropriate case and then call cases
! to compute the parameters of the spline.
            call choose(table_x(table_idx), table_y(table_idx), &
                 table_slope(table_idx), table_slope(table_idx1), &
                 table_x(table_idx1), table_y(table_idx1), eps_tol, spline_case)
            call cases(table_x(table_idx), table_y(table_idx), &
                 table_slope(table_idx), table_slope(table_idx1), &
                 table_x(table_idx1), table_y(table_idx1), spline_e1, spline_e2, &
                 spline_v1, spline_v2, spline_w1, spline_w2, spline_z1, &
                 spline_z2, spline_y1, spline_y2, spline_case)
         end if
      end if
      if (eval_x(i) .eq. table_x(table_idx1)) then
! If eval_x(i) is a data point, its image is known.
         eval_y(i) = table_y(table_idx1)
         cycle
      end if
      eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
           table_x(table_idx), table_y(table_idx), table_x(table_idx1), &
           table_y(table_idx1), spline_y1, spline_y2, spline_e2, &
           spline_w2, spline_v2, spline_case)
      end do

! Calculate the images of the points of evaluation whose abscissas
! are greater than the abscissa of the last data point.
      if (end_idx .eq. num_eval_points) then
         return
      end if
      if ((table_idx1 .eq. num_table_points) .and. (eval_x(end_idx) .ne. table_x(num_table_points))) then
         recompute_tail_params = .false.
      end if
      end if

! Previously, when we arrived at 200 or 210, NUM1 could be improperly
! set. The NUM1= lines below protect from that. llp 8/19/08

! Set the error parameter to indicate that extrapolation has occurred.
      if (recompute_tail_params) then
      err_code = 1
      num_minus1 = max(num_table_points - 1, 1)
      call choose(table_x(num_minus1), table_y(num_minus1), &
           table_slope(num_minus1), table_slope(num_table_points), &
           table_x(num_table_points), table_y(num_table_points), &
           eps_tol, spline_case)
      call cases(table_x(num_minus1), table_y(num_minus1), &
           table_slope(num_minus1), table_slope(num_table_points), &
           table_x(num_table_points), table_y(num_table_points), &
           spline_e1, spline_e2, spline_v1, spline_v2, spline_w1, &
           spline_w2, spline_z1, spline_z2, spline_y1, spline_y2, &
           spline_case)
      end if
      end_idx1 = end_idx + 1
      num_minus1 = max(num_table_points - 1, 1)
      do i = end_idx1, num_eval_points
       eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
            table_x(num_minus1), table_y(num_minus1), &
            table_x(num_table_points), table_y(num_table_points), &
            spline_y1, spline_y2, spline_e2, spline_w2, spline_v2, &
            spline_case)
      end do

      return
end subroutine meval

end module numerics_lib
