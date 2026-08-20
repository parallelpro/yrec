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
