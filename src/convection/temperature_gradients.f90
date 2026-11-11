!----------------------------------------------------------------------
! temperature_gradients
!----------------------------------------------------------------------
! The temperature gradient del = dln T/dln P of a shell (formerly
! tpgrad): compute_temperature_gradients. Its steps:
!   1. the radiative gradient (radiative_gradient_of),
!        del_rad = 10^cdelrl kappa L P/(G m T^4) f_T/f_P,
!      L = (L/Lsun) Lsun, cdelrl = log10[3/(64 pi sigma)] and f_P, f_T
!      the rotation factors; the shell is radiative unless del_rad
!      exceeds del_ad by superadiabatic_tol: del = del_rad (del_ad in the
!      overshoot region below the envelope with adiabatic_overshoot);
!   2. above log10_t_cut (deep interiors) a convective shell takes
!      del = del_ad without mixing-length theory (take_adiabatic);
!   3. starspots (spotted_excess, Somers 2014): with a filling factor
!      f, the excess Delta = del_rad - del_ad becomes
!        Delta = del_rad/S - del_ad,  S = f x^4 + 1 - f,
!      with the spots' temperature contrast x, constant or, with
!      spot_depth_varies, x(T) = 1 - (1 - x_s) T_eff,u/T: magnetic
!      equipartition at the surface (P_mag = P/2) and P_mag constant with
!      depth, T_eff,u the unspotted surface's T_eff (from log_teff);
!   4. mixing-length theory (Bohm-Vitense 1958; mixing_length_coefficients
!      and solve_mlt_cubic): with g = G m/r^2, H_P = P/(rho g), the
!      eddy's optical thickness phi = alpha kappa rho H_P,
!      delta = -dln rho/dln T and c_P,
!        A1 = 16 sqrt(2) sigma T^3 kappa/[c_P (1 + phi^2/3) sqrt(g delta Delta/H_P)],
!        A3 = (3/4) phi^2/[(1 + phi^2/3) A1],
!      v in (0, 1] solves A3 v^3 + v^2 + A1 v - 1 = 0, and
!        del = del_ad + Delta v (v + A1) = del_ad + Delta (1 - A3 v^3),
!      so del lies between del_ad and del_ad + Delta. Dense shells
!      convect efficiently (A1 ~ rho^-1/2 -> 0, A3 ~ rho^1/2 -> infinity,
!      del - del_ad ~ Delta A3^-2/3 ~ rho^-1/3); tenuous ones do not
!      (A1 -> infinity, v -> 1/A1, del_rad - del ~ Delta A3/A1^3 ~ rho^2).
!      A negative radicand returns convection_unstable_radicand (with a
!      warning); a cubic that does not converge takes del_ad (logged);
!   5. the convective velocity v_c = v sqrt(g delta alpha^2 H_P Delta/8);
!   6. with want_derivatives, del's log derivatives for the Henyey
!      Jacobian (mlt_derivatives).
! The derivatives, X each of ln T, ln P and ln r (the others fixed):
!   dln del/dX = [ddel_ad/dX + (v^2 + A1 v) dDelta/dX
!                 + Delta (A1 + 2v) dv/dX + Delta A1 v dln A1/dX]/del,
!   dv/dX = -(A1 v dln A1/dX + A3 v^3 dln A3/dX)/(A1 + 2v + 3 A3 v^2)
!           (the cubic differentiated, at the last iterate's slope),
!   dln A1/dX = -w dln phi/dX + dln(T^3 kappa/c_P)/dX
!               - (1/2) dln(g delta Delta/H_P)/dX,
!   dln A3/dX = (2 - w) dln phi/dX - dln A1/dX,
! with w = 2 (phi^2/3)/(1 + phi^2/3) (so 2 - w = 2/(1 + phi^2/3)) and
! phi = alpha kappa P/g. The opacity's derivatives at constant P and T
! are k_T = dln kappa/dln T + (dln kappa/dln rho)(dln rho/dln T) and
! k_P = (dln kappa/dln rho)(dln rho/dln P), so dln phi/dln T = k_T,
! dln phi/dln P = 1 + k_P, dln phi/dln r = 2 (g ~ r^-2), and
!   dln Delta/dln T = [(k_T - 4) del_rad - ddel_ad/dln T]/Delta,
!   dln Delta/dln P = [(1 + k_P) del_rad - ddel_ad/dln P]/Delta,
! dDelta/dln r = 0 and dln A1/dln r = 2 - 2w.
! Above log10_t_cut and after an unconverged cubic the derivatives are
! del_ad's. The spot factor's own T derivative and the spotted del_rad
! in dDelta are not in them (M242, the B list); O1 lists the Jacobian's
! other gaps.
!
! In: log10 T, P, r and m of the shell (cgs), L/Lsun (not a log: the
! callers pass star%luminosity_lsun), the EOS and opacity results for
! the shell (eos_get, kap_get), the rotation factors f_P and f_T, log10
! Teff (the spots' depth-varying contrast) and in_overshoot_region (the
! shell lies between star%envelope_overshoot_bottom_index and
! _top_index, for adiabatic_overshoot; .false. for envelope and
! atmosphere points). Out: del (actual_gradient), del_rad, the three log
! derivatives (set only for a convective shell with want_derivatives),
! v_c (cm/s), is_convective and ierr. The state is all arguments: this
! routine is on the star-blind envelope integration's path.
module temperature_gradients
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use convection_configuration, only: convection_settings, convection_unstable_radicand
      use convection_messages, only: warn_unstable_radicand, log_cubic_nonconvergence
      use eos_lib, only: eos_result
      use kap_lib, only: i_kap, i_dlnkap_dlnrho, i_dlnkap_dlnt
      use physical_constants, only: ln10, cdelrl, cmixl2, cmixl3
      use math_lib
      implicit none
      private
      public :: compute_temperature_gradients

! A shell is radiative unless del_rad exceeds del_ad by more than this.
      real(dp), parameter :: superadiabatic_tol = 1.0d-6
! The cubic: at most max_cubic_iterations steps, converged when a step
! in v (dimensionless, at most 1) is below cubic_tolerance. It starts
! at v = 1/A1, or, when A3 exceeds cubic_large_a3 (the cubic term
! dominates), at A3^cubic_start_exponent, the F77's -0.333333333 for
! -1/3 (the results depend on the literal, so it stays).
      integer, parameter :: max_cubic_iterations = 25
      real(dp), parameter :: cubic_tolerance = 1.0d-10, cubic_large_a3 = 1.0d3, &
           cubic_start_exponent = -0.333333333d0
! The convective velocity (cm/s) of a convective shell whose
! mixing-length solution is not computed (above log10_t_cut, or a cubic
! that did not converge): a token nonzero value.
      real(dp), parameter :: unsolved_convective_velocity = 1.0d-11

! The mixing-length quantities of one shell (steps 4-6).
      type :: mixing_length_state
         real(dp) :: gravity, scale_height, phi_squared, phi_factor, radicand, &
              a1, a3, v, cubic_slope
      end type mixing_length_state

contains

!----------------------------------------------------------------------
! compute_temperature_gradients (formerly tpgrad): the steps above.
subroutine compute_temperature_gradients(log_temperature, log_pressure, &
     eos_res, kap_res, log_radius, log_mass, luminosity_lsun, &
     actual_gradient, radiative_gradient, dgrad_dt_component, &
     dgrad_dp_component, dgrad_dr_component, convective_velocity, &
     want_derivatives, is_convective, pressure_rotation_factor, &
     temperature_rotation_factor, log_teff, in_overshoot_region, ierr)
      real(dp), intent(in) :: log_temperature, log_pressure, log_radius, log_mass, &
           luminosity_lsun
      type(eos_result), intent(in) :: eos_res
      real(dp), intent(in) :: kap_res(:)
      real(dp), intent(out) :: actual_gradient, radiative_gradient
      real(dp), intent(out) :: dgrad_dt_component, dgrad_dp_component, dgrad_dr_component
      real(dp), intent(out) :: convective_velocity
      logical, intent(in) :: want_derivatives
      logical, intent(out) :: is_convective
      real(dp), intent(in) :: pressure_rotation_factor, temperature_rotation_factor, log_teff
      logical, intent(in) :: in_overshoot_region
      integer, intent(out) :: ierr
      type(mixing_length_state) :: mlt
      real(dp) :: excess, gradient_excess
      logical :: converged

      ierr = 0
! 1. The radiative gradient; a radiative shell.
      radiative_gradient = radiative_gradient_of(kap_res(i_kap), luminosity_lsun, &
           log_pressure, log_mass, log_temperature, temperature_rotation_factor, &
           pressure_rotation_factor)
      excess = radiative_gradient - eos_res%adiabatic_gradient
      if (excess <= superadiabatic_tol) then
         is_convective = .false.
         actual_gradient = radiative_gradient
         convective_velocity = 0.0d0
         if (convection_settings%adiabatic_overshoot .and. in_overshoot_region) &
              actual_gradient = eos_res%adiabatic_gradient
         return
      end if
! 2. A convective shell; adiabatic above log10_t_cut.
      is_convective = .true.
      actual_gradient = eos_res%adiabatic_gradient
      if (log_temperature > convection_settings%log10_t_cut) then
         call take_adiabatic(eos_res, want_derivatives, convective_velocity, &
              dgrad_dt_component, dgrad_dp_component, dgrad_dr_component)
         return
      end if
! 3. Starspots.
      if (convection_settings%spot_filling_factor /= 0.0d0) &
           excess = spotted_excess(radiative_gradient, eos_res%adiabatic_gradient, &
           log_teff, log_temperature)
! 4. Mixing-length theory: the coefficients, then the cubic.
      call mixing_length_coefficients(eos_res, kap_res(i_kap), log_radius, log_mass, &
           excess, mlt)
      if (.not. (mlt%radicand > 0.0d0)) then
         call warn_unstable_radicand(excess, eos_res%dlnrho_dlnt, mlt%radicand, &
              kap_res(i_kap), luminosity_lsun, log_pressure, log_mass, log_temperature, &
              temperature_rotation_factor, pressure_rotation_factor, radiative_gradient, &
              eos_res%adiabatic_gradient, log_radius)
         ierr = convection_unstable_radicand
         return
      end if
      mlt%a1 = cmixl3*mlt%phi_factor*eos_res%temperature**3*kap_res(i_kap)/ &
           (eos_res%specific_heat_cp*sqrt(mlt%radicand))
      call solve_mlt_cubic(mlt, converged)
      if (.not. converged) then
         call log_cubic_nonconvergence(log_pressure, log_temperature, kap_res(i_kap), &
              eos_res%specific_heat_cp, eos_res%dlnrho_dlnt)
         call take_adiabatic(eos_res, want_derivatives, convective_velocity, &
              dgrad_dt_component, dgrad_dp_component, dgrad_dr_component)
         return
      end if
      gradient_excess = excess*mlt%v*(mlt%v + mlt%a1)
      actual_gradient = eos_res%adiabatic_gradient + gradient_excess
! 5. The convective velocity. Its radicand has the sign of the one
! tested in step 4.
      convective_velocity = mlt%v*sqrt(mlt%gravity*(-eos_res%dlnrho_dlnt)* &
           (convection_settings%mixing_length_alpha**2)*mlt%scale_height*excess/8.0d0)
! 6. The derivatives.
      if (want_derivatives) call mlt_derivatives(eos_res, kap_res, radiative_gradient, &
           excess, gradient_excess, actual_gradient, mlt, dgrad_dt_component, &
           dgrad_dp_component, dgrad_dr_component)
end subroutine compute_temperature_gradients

!----------------------------------------------------------------------
! radiative_gradient_of: del_rad = 10^cdelrl kappa L P/(G m T^4) f_T/f_P,
! from the opacity, L/Lsun, log10 P, m and T, and the rotation factors.
pure real(dp) function radiative_gradient_of(opacity, luminosity_lsun, log_pressure, &
     log_mass, log_temperature, temperature_rotation_factor, pressure_rotation_factor)
      real(dp), intent(in) :: opacity, luminosity_lsun, log_pressure, log_mass, &
           log_temperature, temperature_rotation_factor, pressure_rotation_factor

      radiative_gradient_of = opacity*luminosity_lsun*exp(ln10*(log_pressure - &
           log_mass - 4d0*log_temperature + convection_settings%log10_solar_luminosity - &
           convection_settings%log10_gravitational_constant + cdelrl))* &
           temperature_rotation_factor/pressure_rotation_factor
end function radiative_gradient_of

!----------------------------------------------------------------------
! take_adiabatic: a convective shell whose mixing-length solution is not
! computed keeps del = del_ad (set by the caller), the token velocity,
! and, with want_derivatives, del_ad's T and P derivatives (none in r).
pure subroutine take_adiabatic(eos_res, want_derivatives, convective_velocity, &
     dgrad_dt_component, dgrad_dp_component, dgrad_dr_component)
      type(eos_result), intent(in) :: eos_res
      logical, intent(in) :: want_derivatives
      real(dp), intent(out) :: convective_velocity
      real(dp), intent(inout) :: dgrad_dt_component, dgrad_dp_component, dgrad_dr_component

      convective_velocity = unsolved_convective_velocity
      if (.not. want_derivatives) return
      dgrad_dt_component = eos_res%adiabatic_gradient_dt
      dgrad_dp_component = eos_res%adiabatic_gradient_dp
      dgrad_dr_component = 0.0d0
end subroutine take_adiabatic

!----------------------------------------------------------------------
! spotted_excess: Delta = del_rad/S - del_ad with S = f x^4 + 1 - f
! (step 3 above); x is the constant contrast, or, with
! spot_depth_varies, x(T) = 1 - (1 - x_s) T_eff,u/T with log10 T_eff,u =
! log10 Teff - (1/4) log10(f x_s^4 + 1 - f).
pure real(dp) function spotted_excess(radiative_gradient, adiabatic_gradient, log_teff, &
     log_temperature)
      real(dp), intent(in) :: radiative_gradient, adiabatic_gradient, log_teff, &
           log_temperature
      real(dp) :: filling_factor, contrast, log_teff_unspotted

      filling_factor = convection_settings%spot_filling_factor
      if (convection_settings%spot_depth_varies) then
         log_teff_unspotted = log_teff - 0.25d0*log10(filling_factor* &
              pow(convection_settings%spot_temp_contrast, 4.0d0) + 1.0d0 - filling_factor)
         contrast = 1.0d0 - (1.0d0 - convection_settings%spot_temp_contrast)* &
              (exp10(log_teff_unspotted))/(exp10(log_temperature))
      else
         contrast = convection_settings%spot_temp_contrast
      end if
      spotted_excess = radiative_gradient/(filling_factor*pow(contrast, 4.0d0) + &
           1.0d0 - filling_factor) - adiabatic_gradient
end function spotted_excess

!----------------------------------------------------------------------
! mixing_length_coefficients: g = G m/r^2, H_P = P/(rho g), phi^2 with
! phi = alpha rho kappa H_P, phi_factor = 1/(1 + phi^2/3) and the
! radicand g delta Delta/H_P of A1 (step 4 above), from the EOS result,
! the opacity, log10 r and m, and the excess Delta.
pure subroutine mixing_length_coefficients(eos_res, opacity, log_radius, log_mass, &
     excess, mlt)
      type(eos_result), intent(in) :: eos_res
      real(dp), intent(in) :: opacity, log_radius, log_mass, excess
      type(mixing_length_state), intent(out) :: mlt
      real(dp) :: phi

      mlt%gravity = exp(ln10*(convection_settings%log10_gravitational_constant + log_mass - &
           log_radius - log_radius))
      mlt%scale_height = eos_res%pressure/(eos_res%density*mlt%gravity)
      phi = convection_settings%mixing_length_alpha*eos_res%density*opacity*mlt%scale_height
      mlt%phi_squared = phi*phi
      mlt%phi_factor = 1.0d0/(1.0d0 + cmixl2*mlt%phi_squared)
      mlt%radicand = excess*mlt%gravity*(-eos_res%dlnrho_dlnt)/mlt%scale_height
      mlt%a1 = 0.0d0
      mlt%a3 = 0.0d0
      mlt%v = 0.0d0
      mlt%cubic_slope = 0.0d0
end subroutine mixing_length_coefficients

!----------------------------------------------------------------------
! solve_mlt_cubic: A3 = (3/4) phi^2 phi_factor/A1, then v from
! f(v) = A3 v^3 + v^2 + A1 v - 1 = 0 by third-order (Chebyshev) steps,
!   dv = [f/f'] [1 + (f/f') f''/(2 f')],  f' = A1 + 2v + 3 A3 v^2,
! with v capped at 1 before each step. cubic_slope is the f' of the last
! step (at the iterate before it), which the derivatives use.
pure subroutine solve_mlt_cubic(mlt, converged)
      type(mixing_length_state), intent(inout) :: mlt
      logical, intent(out) :: converged
      real(dp) :: three_a3, step
      integer :: iteration

      mlt%v = 1.0d0/mlt%a1
      mlt%a3 = 0.75d0*mlt%phi_squared*mlt%phi_factor/mlt%a1
      three_a3 = 3.0d0*mlt%a3
      if (mlt%a3 > cubic_large_a3) mlt%v = pow(mlt%a3, cubic_start_exponent)
      do iteration = 1, max_cubic_iterations
         mlt%v = min(mlt%v, 1.0d0)
         mlt%cubic_slope = mlt%a1 + mlt%v*(2.0d0 + mlt%v*three_a3)
         step = (-1.0d0 + mlt%v*(mlt%a1 + mlt%v*(1.0d0 + mlt%v*mlt%a3)))/mlt%cubic_slope
         step = step*(1.0d0 + step*(1.0d0 + mlt%v*three_a3)/mlt%cubic_slope)
         mlt%v = mlt%v - step
         if (abs(step) < cubic_tolerance) exit
      end do
      converged = iteration <= max_cubic_iterations
end subroutine solve_mlt_cubic

!----------------------------------------------------------------------
! mlt_derivatives: dln del/dln T, dln P and dln r (the header's
! equations), from the EOS's and the opacity's derivatives, the
! unspotted del_rad, the excess Delta, del - del_ad (gradient_excess),
! del and the mixing-length state. In the code: qa1* = dln A1/dX,
! qa3* = dln A3/dX, qv* = dv/dX, qddel* = dln Delta/dX.
pure subroutine mlt_derivatives(eos_res, kap_res, radiative_gradient, excess, &
     gradient_excess, actual_gradient, mlt, dgrad_dt_component, dgrad_dp_component, &
     dgrad_dr_component)
      type(eos_result), intent(in) :: eos_res
      real(dp), intent(in) :: kap_res(:), radiative_gradient, excess, gradient_excess, &
           actual_gradient
      type(mixing_length_state), intent(in) :: mlt
      real(dp), intent(out) :: dgrad_dt_component, dgrad_dp_component, dgrad_dr_component
      real(dp) :: dadiabatic_dlnt, dadiabatic_dlnp, dlnkap_dlnt_p, dlnkap_dlnp_t, &
           qddelt, qddelp, w, qa1t, qa1p, qa1r, qa3t, qa3p, qa3r, a1_v, a3_v3, &
           inverse_slope, qvt, qvp, qvr, dexcess_dv, excess_a1_v, inverse_gradient

      dadiabatic_dlnt = eos_res%adiabatic_gradient_dt*eos_res%adiabatic_gradient
      dadiabatic_dlnp = eos_res%adiabatic_gradient_dp*eos_res%adiabatic_gradient
      dlnkap_dlnt_p = kap_res(i_dlnkap_dlnrho)*eos_res%dlnrho_dlnt + kap_res(i_dlnkap_dlnt)
      dlnkap_dlnp_t = kap_res(i_dlnkap_dlnrho)*eos_res%dlnrho_dlnp
      qddelt = ((-4.0d0 + dlnkap_dlnt_p)*radiative_gradient - dadiabatic_dlnt)/excess
      qddelp = ((+1.0d0 + dlnkap_dlnp_t)*radiative_gradient - dadiabatic_dlnp)/excess
      w = 2.0d0*cmixl2*mlt%phi_squared*mlt%phi_factor
      qa1t = dlnkap_dlnt_p*(1d0 - w) - eos_res%specific_heat_cp_dt - 0.5d0* &
           (eos_res%dlnrho_dlnt + qddelt + eos_res%dlnrho_dlnt_dt/eos_res%dlnrho_dlnt) + 3d0
      qa1p = dlnkap_dlnp_t*(1d0 - w) - w - eos_res%specific_heat_cp_dp - 0.5d0* &
           (eos_res%dlnrho_dlnp + qddelp + eos_res%dlnrho_dlnt_dp/eos_res%dlnrho_dlnt) + 0.5d0
      qa1r = -2.0d0*w + 2.0d0
      qa3t = 2.0d0*mlt%phi_factor*dlnkap_dlnt_p - qa1t
      qa3p = 2.0d0*mlt%phi_factor*(1.0d0 + dlnkap_dlnp_t) - qa1p
      qa3r = 4.0d0*mlt%phi_factor - qa1r
      a1_v = mlt%a1*mlt%v
      a3_v3 = mlt%a3*mlt%v*mlt%v*mlt%v
      inverse_slope = 1.0d0/mlt%cubic_slope
      qvt = -(a1_v*qa1t + a3_v3*qa3t)*inverse_slope
      qvp = -(a1_v*qa1p + a3_v3*qa3p)*inverse_slope
      qvr = -(a1_v*qa1r + a3_v3*qa3r)*inverse_slope
      dexcess_dv = excess*(mlt%a1 + mlt%v + mlt%v)
      excess_a1_v = excess*mlt%a1*mlt%v
      inverse_gradient = 1.0d0/actual_gradient
      dgrad_dt_component = (dadiabatic_dlnt + qddelt*gradient_excess + dexcess_dv*qvt + &
           excess_a1_v*qa1t)*inverse_gradient
      dgrad_dp_component = (dadiabatic_dlnp + qddelp*gradient_excess + dexcess_dv*qvp + &
           excess_a1_v*qa1p)*inverse_gradient
      dgrad_dr_component = (dexcess_dv*qvr + excess_a1_v*qa1r)*inverse_gradient
end subroutine mlt_derivatives

! Sources: TPGRAD of the Yale code; the starspots are G. Somers's (9/14).
! In 2026 its 33 positional arguments became the EOS's and the
! opacity's result arrays, its stop on a negative radicand an ierr (the
! drivers stop on it), and the overshoot zone's index the logical
! in_overshoot_region.

end module temperature_gradients
