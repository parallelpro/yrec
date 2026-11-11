!----------------------------------------------------------------------
! t_tau_relations
!----------------------------------------------------------------------
! The analytic T-tau relations. Everything an atmosphere choice
! contributes is here: the relation T(tau; Teff), the starting
! temperature of the atmosphere integration (tau -> 0) and its
! photospheric limit, so adding a relation is one block in each of the
! three functions.
!
! Choices (atm_choice; 3, 4 and 5 are the tables, handled by
! atm_get_surface_pt, for which the three functions return NaN):
!   0  Eddington gray
!   1  Krishna-Swamy (ApJ 145, 176, 1966; solar fit)
!   2  Harvard-Smithsonian reference atmosphere (a polynomial fit;
!      hsra_log_t_offset is the run's HSRA offset, hsra_t_tau_offset())
!
! The starting guess for Krishna-Swamy uses the literal 0.550d0, not
! the relation at tau = 0 (1.39d0-0.815d0-0.025d0 need not equal 0.550d0
! bitwise).
module t_tau_relations
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: cc23
      use atm_choices, only: atm_eddington_gray, atm_krishna_swamy, atm_hsra
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use math_lib
      implicit none
      private
      public :: ttau_log10_temperature, ttau_start_log10_temperature, &
           ttau_photosphere_x_limit, hsra_t_tau_offset

contains

! The value of the three functions below for an atm_choice that is not a
! T-tau relation (a table, 3-5, or an unknown choice, which atm_init
! refuses): NaN. The envelope integrates a T-tau relation only when no
! table gave the surface, so with a choice of 0, 1 or 2.
pure function not_a_relation()
      real(dp) :: not_a_relation
      not_a_relation = ieee_value(1.0_dp, ieee_quiet_nan)
end function not_a_relation

! ---------------------------------------------------------------
! log10 T at optical depth tau for the chosen relation.
pure function ttau_log10_temperature(tau, log10_teff, atm_choice, hsra_log_t_offset) &
     result(log10_t)
      real(dp), intent(in) :: tau, log10_teff, hsra_log_t_offset
      integer, intent(in) :: atm_choice
      real(dp) :: log10_t

      select case (atm_choice)
      case (atm_eddington_gray)
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(tau + cc23)
      case (atm_krishna_swamy)
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(tau + &
              1.39d0 - 0.815d0*exp(-2.54d0*tau) - 0.025d0*exp(-30.0d0*tau))
      case (atm_hsra)
         log10_t = log10_teff + harvard_t_tau(tau) - hsra_log_t_offset
      case default
         log10_t = not_a_relation()
      end select
end function ttau_log10_temperature

! ---------------------------------------------------------------
! Starting log10 T for the atmosphere integration (tau near zero;
! the HSRA fit is evaluated at tau = 2/3, its historical anchor).
pure function ttau_start_log10_temperature(log10_teff, atm_choice, hsra_log_t_offset) &
     result(log10_t)
      real(dp), intent(in) :: log10_teff, hsra_log_t_offset
      integer, intent(in) :: atm_choice
      real(dp) :: log10_t

      select case (atm_choice)
      case (atm_eddington_gray)
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(cc23)
      case (atm_krishna_swamy)
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(0.550d0)
      case (atm_hsra)
         log10_t = log10_teff + harvard_t_tau(cc23) - hsra_log_t_offset
      case default
         log10_t = not_a_relation()
      end select
end function ttau_start_log10_temperature

! ---------------------------------------------------------------
! End of the atmosphere integration in x = log10 tau (+ opacity
! terms): tau = 2/3 at Teff for Eddington/HSRA, 0.312 for
! Krishna-Swamy.
pure function ttau_photosphere_x_limit(atm_choice) result(x_limit)
      integer, intent(in) :: atm_choice
      real(dp) :: x_limit

      select case (atm_choice)
      case (atm_krishna_swamy)
         x_limit = -0.505627854d0
      case (atm_eddington_gray, atm_hsra)
         x_limit = -0.176091259d0
      case default
         x_limit = not_a_relation()
      end select
end function ttau_photosphere_x_limit


! ---------------------------------------------------------------
! The HSRA anchor: the fit evaluated at tau = 2/3, which the star layer
! stores as the run's HSRA offset (star%hsra_log_t_offset).
pure function hsra_t_tau_offset() result(offset)
      real(dp) :: offset
      offset = harvard_t_tau(cc23)
end function hsra_t_tau_offset

! ---------------------------------------------------------------
! Polynomial fit to the Harvard-Smithsonian reference atmosphere
! T-tau relation.
pure function harvard_t_tau(optical_depth) result(fit_value)
      real(dp), intent(in) :: optical_depth
      real(dp) :: log10_tau, log10_tau2, log10_tau3, log10_tau4, &
           log10_tau5, log10_tau6, log10_tau7, log10_tau8, log10_tau9
      real(dp) :: fit_value
      log10_tau = log10(optical_depth)
      fit_value = 3.81152046471d0
      fit_value = fit_value + 0.146133736471d0*log10_tau
      log10_tau2 = log10_tau*log10_tau
      fit_value = fit_value + 0.0267719174279d0*log10_tau2
      log10_tau3 = log10_tau2*log10_tau
      fit_value = fit_value - 0.029280655317d0*log10_tau3
      log10_tau4 = log10_tau3*log10_tau
      fit_value = fit_value - 0.0123814456666d0*log10_tau4
      log10_tau5 = log10_tau4*log10_tau
      fit_value = fit_value + 0.00285734990893d0*log10_tau5
      log10_tau6 = log10_tau5*log10_tau
      fit_value = fit_value + 0.0024575213331d0*log10_tau6
      log10_tau7 = log10_tau6*log10_tau
      fit_value = fit_value + 0.000521560431455d0*log10_tau7
      log10_tau8 = log10_tau7*log10_tau
      fit_value = fit_value + 4.71770176883d-5*log10_tau8
      log10_tau9 = log10_tau8*log10_tau
      fit_value = fit_value + 1.58685112637d-6*log10_tau9
end function harvard_t_tau

end module t_tau_relations
