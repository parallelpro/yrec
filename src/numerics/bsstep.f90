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
! envint.f90.
subroutine bsstep(y, dydx, num_eqs, indep_var, h_step, tolerance, y_scale, &
     h_did, h_next, deriv, luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, in_atmosphere, &
     want_derivatives, conductive_opacity_flag, print_flag, log10_radius, &
     log10_teff, hydrogen_fraction, metal_fraction, call_count, saha_state, &
     step_err)
      use numerics_lib
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
