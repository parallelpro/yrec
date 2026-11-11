!----------------------------------------------------------------------
! t_tau_relations
!----------------------------------------------------------------------
! The analytic T-tau relations: everything a relation contributes to the
! atmosphere integration is here (a star-blind kernel):
!     ttau_log10_temperature        log10 T at optical depth tau
!     ttau_start_log10_temperature  log10 T where the integration starts
!                                   (tau -> 0)
!     ttau_photosphere_x_limit      log10 tau where it ends, at T = Teff
!     hsra_t_tau_offset             the HSRA fit at tau = 2/3, the run's
!                                   HSRA offset
! The relations (atm_choice; atm_choices names them):
!   0  Eddington gray: T^4 = (3/4) Teff^4 (tau + 2/3);
!   1  Krishna-Swamy (1966; a solar fit): T^4 = (3/4) Teff^4 (tau + q(tau)),
!      q(tau) = 1.39 - 0.815 exp(-2.54 tau) - 0.025 exp(-30 tau);
!   2  the Harvard-Smithsonian reference atmosphere: a polynomial fit in
!      log10 tau, shifted so that T = Teff at tau = 2/3
!      (hsra_log_t_offset, from hsra_t_tau_offset()).
! In log10, log T = log Teff + (1/4) log10(3/4) + (1/4) log10(tau + q).
! A table's choice (3-5), or an unknown one, gives NaN: the envelope
! integrates a relation only when no table gave the surface, and
! atm_init refuses an unknown choice. Adding a relation is a case in each
! of the three functions and a row in atm_choices.
module t_tau_relations
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use physical_constants, only: cc23
      use atm_choices, only: atm_eddington_gray, atm_krishna_swamy, atm_hsra
      use math_lib
      implicit none
      private
      public :: ttau_log10_temperature, ttau_start_log10_temperature, &
           ttau_photosphere_x_limit, hsra_t_tau_offset

! (1/4) log10(3/4) = -0.0312347, kept at the F77's five decimals: T =
! Teff at the photosphere to 3.2e-7 in log10 T (computing it is B)
      real(dp), parameter :: quarter_log10_three_quarters = -0.031235d0
! log10 tau at the photosphere, T = Teff: tau = 2/3 for Eddington and
! HSRA; for Krishna-Swamy the root of tau + q(tau) = 4/3, tau = 0.312156
      real(dp), parameter :: log10_tau_photosphere_eddington = -0.176091259d0
      real(dp), parameter :: log10_tau_photosphere_krishna_swamy = -0.505627854d0
! Krishna-Swamy's q(tau) = ks_q_far - ks_a1 exp(-ks_k1 tau) - ks_a2
! exp(-ks_k2 tau), and its start value q(0), which the F77 wrote as the
! literal 0.550 (1.39 - 0.815 - 0.025 need not equal it bitwise)
      real(dp), parameter :: ks_q_far = 1.39d0, ks_a1 = 0.815d0, ks_k1 = 2.54d0, &
           ks_a2 = 0.025d0, ks_k2 = 30.0d0
      real(dp), parameter :: ks_q_start = 0.550d0
! The HSRA fit's coefficients, of (log10 tau)^0 .. ^9
      real(dp), parameter :: hsra_coefficients(0:9) = [3.81152046471d0, 0.146133736471d0, &
           0.0267719174279d0, -0.029280655317d0, -0.0123814456666d0, 0.00285734990893d0, &
           0.0024575213331d0, 0.000521560431455d0, 4.71770176883d-5, 1.58685112637d-6]

contains

! the value of the three functions for a choice that is not a relation
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
         log10_t = log10_teff + quarter_log10_three_quarters + 0.25d0*log10(tau + cc23)
      case (atm_krishna_swamy)
         log10_t = log10_teff + quarter_log10_three_quarters + 0.25d0*log10(tau + &
              ks_q_far - ks_a1*exp(-ks_k1*tau) - ks_a2*exp(-ks_k2*tau))
      case (atm_hsra)
         log10_t = log10_teff + harvard_t_tau(tau) - hsra_log_t_offset
      case default
         log10_t = not_a_relation()
      end select
end function ttau_log10_temperature

! ---------------------------------------------------------------
! Starting log10 T for the atmosphere integration (tau near zero;
! the HSRA fit is evaluated at tau = 2/3, its anchor).
pure function ttau_start_log10_temperature(log10_teff, atm_choice, hsra_log_t_offset) &
     result(log10_t)
      real(dp), intent(in) :: log10_teff, hsra_log_t_offset
      integer, intent(in) :: atm_choice
      real(dp) :: log10_t

      select case (atm_choice)
      case (atm_eddington_gray)
         log10_t = log10_teff + quarter_log10_three_quarters + 0.25d0*log10(cc23)
      case (atm_krishna_swamy)
         log10_t = log10_teff + quarter_log10_three_quarters + 0.25d0*log10(ks_q_start)
      case (atm_hsra)
         log10_t = log10_teff + harvard_t_tau(cc23) - hsra_log_t_offset
      case default
         log10_t = not_a_relation()
      end select
end function ttau_start_log10_temperature

! ---------------------------------------------------------------
! End of the atmosphere integration, log10 tau at the photosphere (the
! integration's x = log10 tau): tau = 2/3 for Eddington and HSRA, 0.312
! for Krishna-Swamy.
pure function ttau_photosphere_x_limit(atm_choice) result(x_limit)
      integer, intent(in) :: atm_choice
      real(dp) :: x_limit

      select case (atm_choice)
      case (atm_krishna_swamy)
         x_limit = log10_tau_photosphere_krishna_swamy
      case (atm_eddington_gray, atm_hsra)
         x_limit = log10_tau_photosphere_eddington
      case default
         x_limit = not_a_relation()
      end select
end function ttau_photosphere_x_limit

! ---------------------------------------------------------------
! The HSRA anchor: the fit evaluated at tau = 2/3, which the star layer
! stores as the run's HSRA offset (envelope_integration's
! hsra_log_t_offset).
pure function hsra_t_tau_offset() result(offset)
      real(dp) :: offset
      offset = harvard_t_tau(cc23)
end function hsra_t_tau_offset

! ---------------------------------------------------------------
! The polynomial fit to the Harvard-Smithsonian reference atmosphere's
! log10 T in log10 tau, summed term by term in rising powers.
pure function harvard_t_tau(optical_depth) result(fit_value)
      real(dp), intent(in) :: optical_depth
      real(dp) :: fit_value
      real(dp) :: log10_tau, power
      integer :: k

      log10_tau = log10(optical_depth)
      fit_value = hsra_coefficients(0)
      power = 1.0d0
      do k = 1, ubound(hsra_coefficients, 1)
         power = power*log10_tau
         fit_value = fit_value + hsra_coefficients(k)*power
      end do
end function harvard_t_tau

!======================================================================
! Sources: the F77's QATM and ENVINT (the relations, KTTAU 0-2); S. Krishna
! Swamy, ApJ 145, 176 (1966); the Harvard-Smithsonian reference
! atmosphere (Gingerich et al. 1971), as the F77's fit. The 2026 port
! gathered the relations here (ttau_lib, 1ccc13d1).

end module t_tau_relations
