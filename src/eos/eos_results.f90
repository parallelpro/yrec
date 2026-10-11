!----------------------------------------------------------------------
! eos_results
!----------------------------------------------------------------------
! eos_result: what the equation of state returns for one (log T, log P,
! X, Z), one named field per quantity (cgs units). eos_get fills it;
! the EOS families take it as their result argument r.
!
! The record is intent(inout) everywhere, because two of its fields are
! also inputs:
! - density and log10_density are the starting guess of the fully
!   ionized density iteration (the result depends on them within that
!   iteration's tolerance);
! - the six second derivatives are set only with want_derivatives, and
!   otherwise keep their values.
! Every other field is overwritten on each call.
module eos_results
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private

      type, public :: eos_result
! T (K) and P (dyn/cm^2), from the log T and log P asked for
         real(dp) :: temperature = 0.0d0, pressure = 0.0d0
! density (g/cm^3) and its log10
         real(dp) :: log10_density = 0.0d0, density = 0.0d0
! beta = P_gas/P, 1/beta, and beta14 (4(1-beta)/beta on the Yale path,
! 1-beta after OPAL or MHD)
         real(dp) :: beta = 0.0d0, beta_inverse = 0.0d0, beta14 = 0.0d0
! the ionization fractions of H, He+ and He++
         real(dp) :: ion_fraction(3) = 0.0d0
! R (1/mu_ion + 1/mu_e) in erg/(g K), 1/mu_ion, 1/mu_e
         real(dp) :: specific_gas_constant = 0.0d0, ion_mean_weight_inverse = 0.0d0, &
              electron_mean_weight_inverse = 0.0d0
! the electron degeneracy parameter eta (fully ionized gas; 0 otherwise)
         real(dp) :: electron_degeneracy_parameter = 0.0d0
! dln rho/dln T at constant P, dln rho/dln P at constant T, cp in
! erg/(g K), and the adiabatic gradient
         real(dp) :: dlnrho_dlnt = 0.0d0, dlnrho_dlnp = 0.0d0, &
              specific_heat_cp = 0.0d0, adiabatic_gradient = 0.0d0
! the second derivatives (set with want_derivatives):
! dlnrho_dlnt_dt = d(dln rho/dln T)/dln T, dlnrho_dlnt_dp =
! d(dln rho/dln T)/dln P, adiabatic_gradient_dt/_dp = dln grad_ad/dln T,
! dln grad_ad/dln P, specific_heat_cp_dt/_dp = dln cp/dln T, dln cp/dln P
         real(dp) :: dlnrho_dlnt_dt = 0.0d0, dlnrho_dlnt_dp = 0.0d0, &
              adiabatic_gradient_dt = 0.0d0, adiabatic_gradient_dp = 0.0d0, &
              specific_heat_cp_dt = 0.0d0, specific_heat_cp_dp = 0.0d0
      end type eos_result

end module eos_results
