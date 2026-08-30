!----------------------------------------------------------------------
! ttau_lib
!----------------------------------------------------------------------
! The analytic T-tau relations, in ONE place (2026, follow-up to the
! envint purity split). Everything an atmosphere choice contributes
! lives here: the relation T(tau; Teff), the tau -> 0 starting
! temperature for the integration, and the photospheric integration
! limit -- so adding a relation is one self-contained block instead
! of three coordinated branches (atmosphere_derivs' dispatch, the
! kernel's starting guess, the kernel's x_limit).
!
! Choices (the atm_choice control; 3/4/5 are the TABULATED surface
! boundaries handled by prepare_surface_boundary, never seen here):
!   0  Eddington gray
!   1  Krishna-Swamy (Ap.J. 145, 176, 1966; solar fit)
!   2  Harvard-Smithsonian reference atmosphere (polynomial fit in
!      setup/harvard_t_tau.f90; atm_hras is the run's HSRA offset)
!
! Star-blind (enforced by check_boundaries' STAR_BLIND_FILES):
! atm_hras arrives as an argument. Expressions are verbatim from the
! former atmosphere_derivs statement functions and envint kernel
! branches -- including the start guess's 0.550d0 literal for
! Krishna-Swamy, which is NOT rewritten as the relation at tau = 0
! (1.39d0-0.815d0-0.025d0 need not equal 0.550d0 bitwise).
module ttau_lib
      use math_lib
      use phys_const_lib, only: cc23
      implicit none
      private
      public :: ttau_log10_temperature, ttau_start_log10_temperature, &
           ttau_photosphere_x_limit

      double precision :: harvard_t_tau
      external harvard_t_tau

contains

! ---------------------------------------------------------------
! log10 T at optical depth tau for the chosen relation.
function ttau_log10_temperature(tau, log10_teff, atm_choice, atm_hras) &
     result(log10_t)
      double precision, intent(in) :: tau, log10_teff, atm_hras
      integer, intent(in) :: atm_choice
      double precision :: log10_t

      if (atm_choice == 0) then
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(tau + cc23)
      else if (atm_choice == 1) then
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(tau + &
              1.39d0 - 0.815d0*exp(-2.54d0*tau) - 0.025d0*exp(-30.0d0*tau))
      else
         log10_t = log10_teff + harvard_t_tau(tau) - atm_hras
      end if
end function ttau_log10_temperature

! ---------------------------------------------------------------
! Starting log10 T for the atmosphere integration (tau near zero;
! the HSRA fit is evaluated at tau = 2/3, its historical anchor).
function ttau_start_log10_temperature(log10_teff, atm_choice, atm_hras) &
     result(log10_t)
      double precision, intent(in) :: log10_teff, atm_hras
      integer, intent(in) :: atm_choice
      double precision :: log10_t

      if (atm_choice == 0) then
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(cc23)
      else if (atm_choice == 1) then
         log10_t = log10_teff - 0.031235d0 + 0.25d0*log10(0.550d0)
      else
         log10_t = log10_teff + harvard_t_tau(cc23) - atm_hras
      end if
end function ttau_start_log10_temperature

! ---------------------------------------------------------------
! End of the atmosphere integration in x = log10 tau (+ opacity
! terms): tau = 2/3 at Teff for Eddington/HSRA, 0.312 for
! Krishna-Swamy.
function ttau_photosphere_x_limit(atm_choice) result(x_limit)
      integer, intent(in) :: atm_choice
      double precision :: x_limit

      if (atm_choice == 1) then
         x_limit = -0.505627854d0
      else
         x_limit = -0.176091259d0
      end if
end function ttau_photosphere_x_limit

end module ttau_lib
