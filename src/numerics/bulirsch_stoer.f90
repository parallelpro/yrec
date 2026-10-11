!----------------------------------------------------------------------
! bulirsch_stoer
!----------------------------------------------------------------------
! The Bulirsch-Stoer integrator (Numerical Recipes' bsstep, mmid and
! rzextr; YREC's numerics_lib, split 2026-10). bsstep tries one step of
! h with the modified midpoint rule (mmid) at increasing substep
! counts, extrapolates the results to zero step size (ratext), and
! accepts the step once the error is within tolerance, or retries a
! shorter one. The equations are an ode_system: the caller extends it
! with whatever its derivatives need and binds them as derivatives.
! (The envelope integration's twelve values used to pass through
! bsstep and mmid as arguments.)
module bulirsch_stoer
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      use math_lib
      implicit none
      private
      public :: ode_system, bsstep

! A system of first-order equations y' = f(x, y)
      type, abstract :: ode_system
      contains
         procedure(ode_derivatives), deferred :: derivatives
      end type ode_system

      abstract interface
! dydx = f(x, y); ierr /= 0 when they cannot be evaluated
         subroutine ode_derivatives(self, x, y, dydx, ierr)
            import :: ode_system, dp
            class(ode_system), intent(inout) :: self
            real(dp), intent(in) :: x, y(:)
            real(dp), intent(out) :: dydx(:)
            integer, intent(out) :: ierr
         end subroutine ode_derivatives
      end interface

! ratext's tableau, carried from one stage of a step to the next
! (formerly SAVEd in ratext): the step sizes squared and the diagonals
      integer, parameter :: imax = 11, nmax = 15, ncol = 7
      type :: rational_tableau
         real(dp) :: x_hist(imax) = 0.0d0
         real(dp) :: tableau(nmax,ncol) = 0.0d0
      end type rational_tableau

contains

!----------------------------------------------------------------------
! bsstep
!----------------------------------------------------------------------
! One Bulirsch-Stoer step of the num_eqs equations of system from
! indep_var, with y and its derivatives dydx there, trying h_step: on
! success indep_var advances by h_did, y holds the new values and
! h_next is the suggested next step; step_err is the error estimate
! scaled by y_scale. Stages i = 1..max_stage_index use
! substep_sequence(i) midpoint substeps, extrapolated over up to
! extrap_order stages; when none meets tolerance, h shrinks by
! 4*2**((max_stage_index - extrap_order)/2) and the step is retried.
! ierr: the derivatives' failure, or 1 when h underflows, i.e. when
! underflow_reference + h == underflow_reference. (The envelope passes
! its hydrogen fraction there: the F77 BSSTEP tested X+H.EQ.X with X
! its hydrogen-fraction argument, not the independent variable X0.)
subroutine bsstep(system, y, dydx, num_eqs, indep_var, h_step, tolerance, y_scale, &
     h_did, h_next, step_err, max_stage_index, extrap_order, &
     underflow_reference, ierr)
      real(dp), parameter :: one = 1.0d0, shrink_factor = 0.95d0, &
           grow_factor = 1.2d0

      class(ode_system), intent(inout) :: system
      integer, intent(in) :: num_eqs
      real(dp), intent(inout) :: y(num_eqs)
      real(dp), intent(in) :: dydx(num_eqs)
      real(dp), intent(inout) :: indep_var
      real(dp), intent(in) :: h_step, tolerance, y_scale(num_eqs)
      real(dp), intent(out) :: h_did, h_next
      real(dp), intent(out) :: step_err(num_eqs)
      integer, intent(in) :: max_stage_index, extrap_order
      real(dp), intent(in) :: underflow_reference
      integer, intent(out) :: ierr

! the number of midpoint substeps at each stage
      integer, parameter :: substep_sequence(11) = [2,4,6,8,12,16,24,32,48,64,96]
      real(dp) :: y_err(num_eqs), y_sav(num_eqs), dy_sav(num_eqs), y_seq(num_eqs)
      real(dp) :: h, x_sav, x_est, err_max
      integer :: i, j
      integer :: jerr_integrand
      type(rational_tableau) :: tableau

      ierr = 0
      h = h_step
      x_sav = indep_var
      do i = 1,num_eqs
       y_sav(i) = y(i)
       dy_sav(i) = dydx(i)
      end do
      do    ! step-halving retry loop (was label 20)
      do i = 1,max_stage_index
       call mmid(system, y_sav, dy_sav, num_eqs, x_sav, h, substep_sequence(i), &
            y_seq, jerr_integrand)
       if (jerr_integrand /= 0) then
          ierr = jerr_integrand
          return
       end if
       x_est = (h/substep_sequence(i))**2
       call ratext(tableau, i, x_est, y_seq, y, y_err, num_eqs, extrap_order)
       do j = 1,num_eqs
          step_err(j) = dabs(y_err(j)/y_scale(j))
       end do
! The largest scaled error among the components that are not NaN, and
! NaN only when every component is (then no stage is accepted and the
! step ends on underflow). This is how the compiled loop
! err_max = dmax1(err_max, ...) from err_max = 0 behaved: it folded the
! first max to the bare value and used a NaN-skipping max after it. The
! clean stop of an atmosphere integration (one equation) that starts
! from NaN, the pinned 0.03 Msun case, depends on it.
       err_max = step_err(1)
       do j = 2,num_eqs
          if (step_err(j) > err_max .or. ieee_is_nan(err_max)) err_max = step_err(j)
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
      if(underflow_reference+h.eq.underflow_reference) then
         write(*,*) 'ERROR IN BSSTEP'
       ierr = 1
       return
      end if
      end do

end subroutine bsstep

!----------------------------------------------------------------------
! mmid
!----------------------------------------------------------------------
! The modified midpoint rule (Numerical Recipes, p. 562): y at
! x_start + h_total from y and dydx at x_start in n_step substeps, into
! y_out, with system's derivatives.
subroutine mmid(system, y, dydx, n_var, x_start, h_total, n_step, y_out, ierr)
      class(ode_system), intent(inout) :: system
      integer, intent(in) :: n_var
      real(dp), intent(in) :: y(n_var), dydx(n_var)
      real(dp), intent(in) :: x_start, h_total
      integer, intent(in) :: n_step
      real(dp), intent(out) :: y_out(n_var)
      integer, intent(out) :: ierr

! h_sub is the size of each small step.
      real(dp) :: y_mid(n_var), y_new(n_var)
      real(dp) :: h_sub, h_sub2, x_current, y_swap
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
      call system%derivatives(x_current, y_new, y_out, ierr)
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
       call system%derivatives(x_current, y_new, y_out, ierr)
       if (ierr /= 0) return
      end do
! last step.
      do i = 1,n_var
       y_out(i) = 0.5d0*(y_mid(i) + y_new(i) + h_sub*y_out(i))
      end do
      return
end subroutine mmid

!----------------------------------------------------------------------
! ratext
!----------------------------------------------------------------------
! Diagonal rational-function extrapolation to zero step size
! (Numerical Recipes' rzextr, p. 566): adds the estimate y_est at
! x_est = (h/n)**2 as stage est_index of the tableau tab, using at most
! max_use stages, and returns the extrapolated y_extrap and its error
! y_err. Stage 1 starts a new tableau; the caller keeps tab between the
! stages of one step.
subroutine ratext(tab, est_index, x_est, y_est, y_extrap, y_err, num_vars, &
     max_use)
      type(rational_tableau), intent(inout) :: tab
      integer, intent(in) :: est_index, num_vars, max_use
      real(dp), intent(in) :: x_est
      real(dp), intent(in) :: y_est(num_vars)
      real(dp), intent(out) :: y_extrap(num_vars), y_err(num_vars)
      real(dp) :: fx(ncol)
      integer :: var_idx, k, k2, num_use
      real(dp) :: yy, v, c, b1, b, ddy

! SAVE CURRENT INDEPENDENT VARIABLE.
      tab%x_hist(est_index) = x_est
      if (est_index.eq.1) then
         do var_idx = 1, num_vars
            y_extrap(var_idx) = y_est(var_idx)
            tab%tableau(var_idx,1) = y_est(var_idx)
            y_err(var_idx) = y_est(var_idx)
         end do
      else
!        USE AT MOST max_use PREVIOUS MEMBERS.
         num_use = min(est_index,max_use)
         do k = 1, num_use-1
            fx(k+1) = tab%x_hist(est_index-k)/x_est
         end do
!        EVALUATE NEXT DIAGONAL IN TABLEAU.
         do var_idx=1,num_vars
            yy = y_est(var_idx)
            v = tab%tableau(var_idx,1)
            c = yy
            tab%tableau(var_idx,1) = yy
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
               v = tab%tableau(var_idx,k2)
               tab%tableau(var_idx,k2) = ddy
               yy = yy + ddy
            end do
            y_err(var_idx) = ddy
            y_extrap(var_idx) = yy
         end do
      end if
      return
end subroutine ratext

end module bulirsch_stoer
