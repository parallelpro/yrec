!----------------------------------------------------------------------
! yale_eos
!----------------------------------------------------------------------
! The Yale equation of state of Prather (1976): analytic gas physics,
! no tables but the Fermi-Dirac integrals.
!     compute_saha_eos           partial ionization: the Saha equations
!                                for H, He and 11 metals, solved by
!                                Newton iteration on the free-electron
!                                number
!     compute_fully_ionized_eos  full ionization, with a degenerate,
!                                relativistic electron gas and the
!                                optional Debye-Hueckel correction
! The Fermi-Dirac table and the Debye-Hueckel state live in
! yale_eos_tables (yale).
module yale_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use yale_eos_tables, only: debye_huckel_coefficient, debye_huckel_nu, &
           debye_huckel_x, debye_huckel_y, debye_huckel_z, &
           debye_huckel_z_total, yale
      use physical_constants, only: clni, cmkh, gas_constant, ln10
      use io_units, only: run_log_unit
      use eos_mixture, only: ix_he, n_mix_species
      use eos_configuration, only: eos_settings
      use math_lib
      implicit none
      private
      public :: compute_saha_eos, compute_fully_ionized_eos, read_fermi_table

! Single-precision literals, as in eqsaha.f: widened to double, they
! differ from the correctly rounded doubles in the last bits. The
! three saved values are the Saha solver's starting guesses (E, x1, x2),
! the previous call's solution.
      real(dp), parameter :: ionization_temp(11) = [real(dp) :: 59630., &
           69450., 88710., 91330., 94570., 130630., 157800., 158000., 168630., &
           183080., 250200.]
      real(dp), parameter :: helium_ionization_temp_1 = 285270.
      real(dp), parameter :: helium_ionization_temp_2 = 631370.
      real(dp), parameter :: saha_weight_term(13) = [real(dp) :: -.01, -.47, &
           .60, .49, .08, .11, .0, -.04, .64, .99, 1.03, .6, .0]
      integer, parameter :: nz = 11
      integer, parameter :: max_saha_iterations = 20
      real(dp), parameter :: saha_exponent_tol = 4.0D1
      real(dp), save :: saha_e_guess = 0.1D0, saha_x1_guess = 0.0D0, saha_x2_guess = 0.0D0
      real(dp), parameter :: saha_convergence_tol = 1.0D-09

! The Saha system of one compute_saha_eos call: its constants, the
! solution (E and the He fractions) and the reduced Newton matrix, which
! the derivative steps reuse.
      type :: saha_system
         real(dp) :: helium_mass_fraction = 0.0d0
         real(dp) :: max_electrons_per_ion = 0.0d0
         real(dp) :: temperature_inverse = 0.0d0
         real(dp) :: ionization_temp_over_t(11) = 0.0d0
         real(dp) :: saha_ratio(11) = 0.0d0
         integer :: nz0 = 0
         integer :: nz1 = 0
         logical :: skip_helium_i = .false.
         logical :: skip_helium_ii = .false.
         real(dp) :: helium_saha_ratio_1 = 0.0d0
         real(dp) :: helium_saha_ratio_2 = 0.0d0
         real(dp) :: mean_electrons_per_ion = 0.0d0
         real(dp) :: helium_ion_fraction_1 = 0.0d0
         real(dp) :: helium_ion_fraction_2 = 0.0d0
         real(dp) :: ep1 = 0.0d0
         real(dp) :: ep1e = 0.0d0
         real(dp) :: ep12 = 0.0d0
         real(dp) :: c11 = 0.0d0
         real(dp) :: c12 = 0.0d0
         real(dp) :: c13 = 0.0d0
         real(dp) :: c22 = 0.0d0
         real(dp) :: c23 = 0.0d0
         real(dp) :: c33 = 0.0d0
         real(dp) :: cr12 = 0.0d0
         real(dp) :: cr13 = 0.0d0
         real(dp) :: cr23 = 0.0d0
         real(dp) :: species_ion_fraction(11) = 0.0d0
         real(dp) :: species_weighted_term(11) = 0.0d0
      end type saha_system

! The first derivatives of the Saha solution that the second-derivative
! step reuses (sqet = dE/dln T, sqep = dE/dln P, sx1qt = dx1/dln T, ...).
      type :: saha_sensitivities
         real(dp) :: sqet = 0.0d0
         real(dp) :: sqep = 0.0d0
         real(dp) :: sx0qt = 0.0d0
         real(dp) :: sx0qp = 0.0d0
         real(dp) :: sx1qt = 0.0d0
         real(dp) :: sx1qp = 0.0d0
         real(dp) :: sx2qt = 0.0d0
         real(dp) :: sx2qp = 0.0d0
         real(dp) :: sk0qt = 0.0d0
         real(dp) :: sk1qt = 0.0d0
         real(dp) :: helium_neutral_fraction = 0.0d0
         real(dp) :: species_temp_deriv_term(11) = 0.0d0
         real(dp) :: qcpi = 0.0d0
         real(dp) :: pdt = 0.0d0
         real(dp) :: pdtq = 0.0d0
         real(dp) :: ramu = 0.0d0
         real(dp) :: betaut = 0.0d0
      end type saha_sensitivities


! The state of one fully ionized solution that its steps share: the
! Fermi-table interpolation (rows, factors, interpolated values), the
! pressure terms and the pressure derivatives.
      type :: fermi_solution
         real(dp) :: beta1 = 0.0d0
         real(dp) :: cl1 = 0.0d0
         real(dp) :: cl2 = 0.0d0
         real(dp) :: cl3 = 0.0d0
         real(dp) :: cmfdh0 = 0.0d0
         real(dp) :: cmfdh1 = 0.0d0
         real(dp) :: d8 = 0.0d0
         real(dp) :: dcl1 = 0.0d0
         real(dp) :: dcl2 = 0.0d0
         real(dp) :: dcl3 = 0.0d0
         real(dp) :: df1(5,3) = 0.0d0
         real(dp) :: df2(5,3) = 0.0d0
         real(dp) :: dl8 = 0.0d0
         real(dp) :: dml = 0.0d0
         real(dp) :: emul = 0.0d0
         real(dp) :: emul8 = 0.0d0
         real(dp) :: ff(5,3) = 0.0d0
         real(dp) :: ffe(5) = 0.0d0
         integer :: id1 = 0
         integer :: id1p = 0
         integer :: id2 = 0
         integer :: id3 = 0
         integer :: jt1 = 0
         integer :: jt2 = 0
         integer :: jt3 = 0
         integer :: nden = 0
         real(dp) :: pa = 0.0d0
         real(dp) :: pdh = 0.0d0
         real(dp) :: pdhp = 0.0d0
         real(dp) :: pe = 0.0d0
         real(dp) :: pinv = 0.0d0
         real(dp) :: pl8 = 0.0d0
         real(dp) :: pr = 0.0d0
         real(dp) :: qpd = 0.0d0
         real(dp) :: qpt = 0.0d0
         real(dp) :: qqpdd = 0.0d0
         real(dp) :: qqpdt = 0.0d0
         real(dp) :: qqptt = 0.0d0
         real(dp) :: ramut = 0.0d0
         real(dp) :: tl8 = 0.0d0
         real(dp) :: ttcu = 0.0d0
         real(dp) :: ttdcu = 0.0d0
         real(dp) :: ttud = 0.0d0
      end type fermi_solution

contains

!----------------------------------------------------------------------
! compute_saha_eos (formerly eqsaha)
!----------------------------------------------------------------------
! Prather partial-ionization equation of state: solves the coupled
! Saha ionization equations (11 metal/light-element species plus a
! 2-stage helium treatment) by Newton iteration on the mean number of
! free electrons per ion, then (if want_derivatives) differentiates
! the result with respect to log T and log P. Called from blend_eos.f90
! whenever the temperature is below the ionization cutoff, either
! directly (Prather EOS) or ahead of the SCV table lookup (to get
! derivatives SCV itself does not supply).
!
! Notation (Prather 1976): E (mean_electrons_per_ion) is the mean
! number of free electrons per ion; x1 and x2 (helium_ion_fraction_1/2)
! the He+ and He++ fractions; K_i (saha_ratio) the Saha constants
! n+ n_e/n0, scaled by 1/beta, for the 11 species and He. Each Newton
! step solves C (dx2, dx1, dE) = R for the residuals of the He++ and He+
! balances and of charge (r1, r2, r3); c11..c33 are C after reduction to
! triangular form, with multipliers cr12, cr13, cr23. The derivative
! steps differentiate the converged system: sqet = dE/dln T,
! sqep = dE/dln P, sx1qt = dx1/dln T, and so on.
subroutine compute_saha_eos(saha_mass_fractions, log10_temperature, &
     log10_pressure, r, want_derivatives, in_atmosphere, saha_state)

      real(dp), intent(in) :: saha_mass_fractions(:)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(eos_result), intent(inout) :: r
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
      type(saha_system) :: s
      type(saha_sensitivities) :: d

      s%helium_mass_fraction = saha_mass_fractions(ix_he)
      s%max_electrons_per_ion = 1.0d0 + s%helium_mass_fraction
      s%temperature_inverse = 1.0d0/r%temperature
      s%mean_electrons_per_ion = saha_e_guess
      s%helium_ion_fraction_1 = saha_x1_guess
      s%helium_ion_fraction_2 = saha_x2_guess
      call set_saha_constants(log10_temperature, log10_pressure, r, s)
      call solve_saha_system(saha_mass_fractions, log10_temperature, log10_pressure, s, saha_state)
      saha_e_guess = s%mean_electrons_per_ion
      saha_x1_guess = s%helium_ion_fraction_1
      saha_x2_guess = s%helium_ion_fraction_2
! the gas: electron weight, gas constant, density, ionization fractions
      r%electron_mean_weight_inverse = s%mean_electrons_per_ion*r%ion_mean_weight_inverse
      r%specific_gas_constant = gas_constant*(r%electron_mean_weight_inverse + &
           r%ion_mean_weight_inverse)
      r%density = r%beta*r%pressure/(r%specific_gas_constant*r%temperature)
      r%log10_density = log10(r%density)
      r%ion_fraction(1) = 0.0d0
      if(s%nz0 >= 7) r%ion_fraction(1) = s%species_ion_fraction(7)
      r%ion_fraction(2) = s%helium_ion_fraction_1
      r%ion_fraction(3) = s%helium_ion_fraction_2
      if (in_atmosphere) return
      call differentiate_saha_solution(s, r, d)
      if (.not. want_derivatives) return
      call differentiate_saha_solution_twice(s, d, r)
end subroutine compute_saha_eos

!----------------------------------------------------------------------
! set_saha_constants
!----------------------------------------------------------------------
! The Saha constants K_i of the 11 species and of He I/He II at (T, P):
! a species whose exponent is below -saha_exponent_tol is neutral and
! ends the list (nz0 species ionize); one above +saha_exponent_tol is
! fully ionized (K = 1e16; nz1 is the first species that is not). He is
! skipped when absent or neutral, and then its fractions are zero.
subroutine set_saha_constants(log10_temperature, log10_pressure, r, s)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(eos_result), intent(in) :: r
      type(saha_system), intent(inout) :: s
      integer :: i
      real(dp) :: stemp

! compute Saha K's
      stemp = 2.50d0*log10_temperature - log10_pressure + cmkh
      s%nz1 = 1
      do i=1,nz
         s%ionization_temp_over_t(i) = ionization_temp(i)*s%temperature_inverse
         s%saha_ratio(i) = ln10*(saha_weight_term(i) + stemp) - &
              s%ionization_temp_over_t(i)
         if(s%saha_ratio(i) < -saha_exponent_tol) exit
         if (s%saha_ratio(i) <= +saha_exponent_tol) then
         s%saha_ratio(i) = r%beta_inverse*exp(s%saha_ratio(i))
         cycle
         end if
         s%saha_ratio(i) = 1.0d16
         s%nz1 = i + 1
      end do
      if (i > nz) then
      i = nz + 1
      end if
      s%nz0 = i - 1
      s%nz1 = min0(nz,s%nz1)
      s%skip_helium_i = .true.
      s%skip_helium_ii= .true.
      if (s%helium_mass_fraction >= 1.0d-10) then
      s%helium_saha_ratio_1 = ln10*(saha_weight_term(12) + stemp) - &
           helium_ionization_temp_1/r%temperature
      if (s%helium_saha_ratio_1 >= -saha_exponent_tol) then
      s%skip_helium_i = .false.
      s%helium_saha_ratio_1 = r%beta_inverse*exp(s%helium_saha_ratio_1)
      s%helium_saha_ratio_2 = ln10*(saha_weight_term(13) + stemp) - &
           helium_ionization_temp_2/r%temperature
      if(s%helium_saha_ratio_2 < -saha_exponent_tol) then
      s%helium_ion_fraction_2 = 0.0d0
      else
      s%skip_helium_ii= .false.
      s%helium_saha_ratio_2 = r%beta_inverse*exp(s%helium_saha_ratio_2)
      end if
      else
      s%helium_ion_fraction_1 = 0.0d0
      s%helium_ion_fraction_2 = 0.0d0
      end if
      else
      s%helium_ion_fraction_1 = 0.0d0
      s%helium_ion_fraction_2 = 0.0d0
      end if
end subroutine set_saha_constants

!----------------------------------------------------------------------
! solve_saha_system
!----------------------------------------------------------------------
! Newton iteration on E and the He fractions (at most
! max_saha_iterations; saha_state counts the iterations), from the
! previous call's solution. The products in the elimination are floored
! at 1e-38 in magnitude, as on the DEC-20 the code was written for.
subroutine solve_saha_system(saha_mass_fractions, log10_temperature, log10_pressure, &
     s, saha_state)
      real(dp), intent(in) :: saha_mass_fractions(:)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(saha_system), intent(inout) :: s
      integer, intent(inout) :: saha_state
      real(dp) :: c12l
      real(dp) :: c13l
      real(dp) :: c21
      real(dp) :: c31
      real(dp) :: c32
      logical :: converged
      real(dp) :: cr12l
      real(dp) :: cr13l
      real(dp) :: delta_electrons_per_ion
      real(dp) :: deltel
      real(dp) :: deltx1
      real(dp) :: deltx2
      real(dp) :: delx1l
      real(dp) :: div
      real(dp) :: eep1
      real(dp) :: fact
      real(dp) :: fact1
      real(dp) :: fact2
      integer :: i
      integer :: isaha
      real(dp) :: r1
      real(dp) :: r1l
      real(dp) :: r2
      real(dp) :: r3
      real(dp) :: temp8

      do isaha=1,max_saha_iterations
         saha_state = saha_state + 1
         s%ep1 = 1.0d0/(1.0d0 + s%mean_electrons_per_ion)
         s%ep1e = s%ep1/s%mean_electrons_per_ion
         s%ep12 = s%ep1*s%ep1
         eep1 = s%mean_electrons_per_ion*s%ep1
         s%c33 = -s%mean_electrons_per_ion*(1.0d0 + s%mean_electrons_per_ion)
         r3 = s%mean_electrons_per_ion
         if (s%nz0 > 0) then
         do i=1,s%nz0
            div = 1.0d0/(s%saha_ratio(i) + eep1)
            temp8 = s%saha_ratio(i)*div
            s%species_ion_fraction(i) = temp8
            temp8 = temp8*saha_mass_fractions(i)
            r3 = r3 - temp8
            temp8 = temp8*eep1*div
            s%c33 = s%c33 - temp8
            s%species_weighted_term(i) = temp8
         end do
         end if
         s%c33 = s%c33*s%ep1e
         r3 = r3 - s%helium_mass_fraction*(s%helium_ion_fraction_1 + &
              s%helium_ion_fraction_2+s%helium_ion_fraction_2)
         converged = abs(r3) < saha_convergence_tol
         if (.not. s%skip_helium_i) then
         s%c22 = s%helium_saha_ratio_1 + eep1
         s%c23 = s%helium_ion_fraction_1*s%ep12
         c32 = s%helium_mass_fraction
         r2 = s%helium_saha_ratio_1*(1.0d0-s%helium_ion_fraction_1- &
              s%helium_ion_fraction_2) - s%helium_ion_fraction_1*eep1
         if (.not. s%skip_helium_ii) then
         s%c11 = eep1
         s%c12 = -s%helium_saha_ratio_2
         s%c13 = s%helium_ion_fraction_2*s%ep12
         c21 = s%helium_saha_ratio_1
         c31 = 2.0d0*s%helium_mass_fraction
         r1 = s%helium_ion_fraction_1*s%helium_saha_ratio_2 - &
              s%helium_ion_fraction_2*eep1
! Reduction of Saha matrix to diagonal form
         s%cr12 = - c21/s%c11
         s%cr13 = - c31/s%c11
         s%c22 = s%c22 + s%cr12*s%c12
         c32 = c32 + s%cr13*s%c12
         s%c23 = s%c23 + s%cr12*s%c13
         s%c33 = s%c33 + s%cr13*s%c13
         if (r1 /= 0.0d0) then
         r1l=log10(abs(r1))
         if (s%cr12 /= 0.0d0) then
         cr12l=log10(abs(s%cr12))
         fact=cr12l+r1l
         if (fact >= -38.0d0) then
         r2 = r2 + s%cr12*r1
         end if
         end if
         if (s%cr13 /= 0.0d0) then
         cr13l=log10(abs(s%cr13))
         fact=cr13l+r1l
         if (fact >= -38.0d0) then
         r3 = r3 + s%cr13*r1
         end if
         end if
         end if
! Entry for no fully ionized helium (SAHEX2 = 0.0)
         end if
         s%cr23 = - c32/s%c22
         s%c33 = s%c33 + s%cr23*s%c23
         r3 = r3 + s%cr23*r2
! Entry for neutral helium
         end if
         delta_electrons_per_ion = r3/s%c33
         if (.not. s%skip_helium_i) then
         deltx1 = (r2 - s%c23*delta_electrons_per_ion)/s%c22
         s%helium_ion_fraction_1 = s%helium_ion_fraction_1 + deltx1
         s%helium_ion_fraction_1 = max(0.0d0,min(1.0d0,s%helium_ion_fraction_1))
         if (.not. s%skip_helium_ii) then
! c12*deltx1, floored at 1e-38 in magnitude (written for the DEC-20's
! single-precision range)
         if (s%c12 /= 0.d0 .and. deltx1 /= 0.d0) then
         c12l=log10(abs(s%c12))
         delx1l=log10(abs(deltx1))
         fact1=c12l+delx1l
         if (fact1 < -38.0d0) then
         fact1=-38.0d0
         fact1=exp(ln10*fact1)*sign(1.0d0,s%c12)*sign(1.0d0,deltx1)
         else
         fact1=s%c12*deltx1
         end if
         else
         fact1=s%c12*deltx1
         end if
         if (s%c13 /= 0.d0 .and. delta_electrons_per_ion /= 0.d0) then
         c13l=log10(abs(s%c13))
         deltel=log10(abs(delta_electrons_per_ion))
         fact2=c13l+deltel
         if (fact2 < -38.0d0) then
         fact2=-38.0d0
         fact2=exp(ln10*fact2)*sign(1.0d0,s%c13)*sign(1.0d0,delta_electrons_per_ion)
         else
         fact2=s%c13*delta_electrons_per_ion
         end if
         else
         fact2=s%c13*delta_electrons_per_ion
         end if
         deltx2=(r1-fact1-fact2)/s%c11
         s%helium_ion_fraction_2 = s%helium_ion_fraction_2 + deltx2
         s%helium_ion_fraction_2 = max(0.0d0,min(1.0d0,s%helium_ion_fraction_2))
         end if
         end if
         s%mean_electrons_per_ion = s%mean_electrons_per_ion + delta_electrons_per_ion
         s%mean_electrons_per_ion = max(1.0d-11,min(s%max_electrons_per_ion, &
              s%mean_electrons_per_ion))
         converged = converged .and. abs(delta_electrons_per_ion) < saha_convergence_tol
         if(converged) exit
      end do
      if (isaha > max_saha_iterations) then
      write(run_log_unit,99) log10_temperature,log10_pressure, &
           s%mean_electrons_per_ion,delta_electrons_per_ion
 99   format(' -----SAHA FAILURE (TL,PL)=',2F10.6,'  (E,DE)=',2F20.12)
      end if
end subroutine solve_saha_system

!----------------------------------------------------------------------
! differentiate_saha_solution
!----------------------------------------------------------------------
! dE/dln T and dE/dln P from the converged system, then dln rho/dln T,
! dln rho/dln P, the internal energy derivative, cp and grad_ad.
subroutine differentiate_saha_solution(s, r, d)
      type(saha_system), intent(in) :: s
      type(eos_result), intent(inout) :: r
      type(saha_sensitivities), intent(inout) :: d
      integer :: i
      real(dp) :: qut
      real(dp) :: r1p
      real(dp) :: r1t
      real(dp) :: r2p
      real(dp) :: r2t
      real(dp) :: r3p
      real(dp) :: r3t
      real(dp) :: stemp
      real(dp) :: temp
      real(dp) :: usum

      r3t = 0.0d0
      r3p = 0.0d0
      if (s%nz0 > 0) then
      do i=s%nz1,s%nz0
         r3p = r3p - s%species_weighted_term(i)
         r3t = r3t -s%species_weighted_term(i)*(2.5d0 + r%beta14 + &
              s%ionization_temp_over_t(i))
      end do
      r3p = -r%beta_inverse*r3p
      end if
      if (.not. s%skip_helium_i) then
      d%helium_neutral_fraction = 1.0d0 - s%helium_ion_fraction_1 - s%helium_ion_fraction_2
      d%sk0qt = 2.5d0 + r%beta14 + helium_ionization_temp_1*s%temperature_inverse
      temp = s%helium_saha_ratio_1*d%helium_neutral_fraction
      r2t = temp*d%sk0qt
      r2p = -temp*r%beta_inverse
      if (.not. s%skip_helium_ii) then
      d%sk1qt = 2.5d0 + r%beta14 + helium_ionization_temp_2*s%temperature_inverse
      temp = s%helium_saha_ratio_2*s%helium_ion_fraction_1
      r1t = temp*d%sk1qt
      r1p = -temp*r%beta_inverse
      r2t = r2t + s%cr12*r1t
      r2p = r2p + s%cr12*r1p
      r3t = r3t + s%cr13*r1t
      r3p = r3p + s%cr13*r1p
      end if
      r3t = r3t + s%cr23*r2t
      r3p = r3p + s%cr23*r2p
      end if
      d%sqet = r3t/s%c33
      d%sqep = r3p/s%c33
      r%dlnrho_dlnt = -1.0d0 - r%beta14 - s%ep1*d%sqet
      r%dlnrho_dlnp = r%beta_inverse - s%ep1*d%sqep
! Compute internal energy temperature derivative(qut)
      usum = 0.0d0
      if (.not. s%skip_helium_i) then
      d%sx1qt = (r2t - s%c23*d%sqet)/s%c22
      d%sx1qp = (r2p - s%c23*d%sqep)/s%c22
      d%sx2qt = 0.0d0
      d%sx2qp = 0.0d0
      if (.not. s%skip_helium_ii) then
      d%sx2qt = (r1t - s%c12*d%sx1qt - s%c13*d%sqet)/s%c11
      d%sx2qp = (r1p - s%c12*d%sx1qp - s%c13*d%sqep)/s%c11
      end if
      d%sx0qt = -(d%sx1qt + d%sx2qt)
      d%sx0qp = -(d%sx1qp + d%sx2qp)
      usum = s%helium_mass_fraction*(d%sx1qt*helium_ionization_temp_1 + &
           d%sx2qt*(helium_ionization_temp_1 + helium_ionization_temp_2))
      end if
      if (s%nz0 > 0) then
      stemp = 2.50d0 + r%beta14 - s%ep1e*d%sqet
      do i=s%nz1,s%nz0
         d%species_temp_deriv_term(i) = s%species_weighted_term(i)*(stemp + &
              s%ionization_temp_over_t(i))
         usum = usum + d%species_temp_deriv_term(i)*ionization_temp(i)
      end do
      end if
      d%betaut = 0.75d0*(2.0d0 + r%beta14)
      d%ramu = gas_constant*r%ion_mean_weight_inverse
      qut =d%ramu*((1d0+s%mean_electrons_per_ion)*(d%betaut+3d0*r%beta_inverse*r%beta14) &
           +d%betaut*d%sqet+usum*s%temperature_inverse)
! Compute specific heat(qcp) and adiabatic gradient(dela)
      d%pdt = r%pressure/(r%density*r%temperature)
      d%pdtq = d%pdt*r%dlnrho_dlnt
      r%specific_heat_cp = qut - d%pdtq
      d%qcpi = 1.0d0/r%specific_heat_cp
      r%adiabatic_gradient = -d%pdtq*d%qcpi
end subroutine differentiate_saha_solution

!----------------------------------------------------------------------
! differentiate_saha_solution_twice
!----------------------------------------------------------------------
! The second derivatives: d(dln rho/dln T)/dln T and /dln P, and the log
! derivatives of cp and grad_ad in T and P.
subroutine differentiate_saha_solution_twice(s, d, r)
      type(saha_system), intent(in) :: s
      type(saha_sensitivities), intent(in) :: d
      type(eos_result), intent(inout) :: r
      real(dp) :: beta16
      real(dp) :: btemp
      real(dp) :: etemp
      integer :: i
      real(dp) :: qqutp
      real(dp) :: qqutt
      real(dp) :: r1p
      real(dp) :: r1t
      real(dp) :: r2p
      real(dp) :: r2t
      real(dp) :: r3p
      real(dp) :: r3t
      real(dp) :: sqetp
      real(dp) :: sqett
      real(dp) :: stemp
      real(dp) :: stemp1
      real(dp) :: stemp2
      real(dp) :: stemp3
      real(dp) :: stemp4
      real(dp) :: stemp5
      real(dp) :: stemp6
      real(dp) :: stemp7
      real(dp) :: sx1qtp
      real(dp) :: sx1qtt
      real(dp) :: sx2qtp
      real(dp) :: sx2qtt
      real(dp) :: upsum
      real(dp) :: utsum

      beta16 = 4.0d0*r%beta_inverse*r%beta14
      r3t = 0.0d0
      r3p = 0.0d0
      utsum = 0.0d0
      upsum = 0.0d0
      if (s%nz0 > 0) then
      etemp = (1.0d0+s%mean_electrons_per_ion+s%mean_electrons_per_ion)/ &
           (s%mean_electrons_per_ion*(1.0d0+s%mean_electrons_per_ion))**2
      stemp1 = 2.5d0 + r%beta14 - d%sqet*s%ep1e
      stemp2 = -r%beta_inverse - d%sqep*s%ep1e
      stemp3 = etemp*d%sqep*d%sqet - r%beta_inverse*r%beta14
      stemp4 = etemp*d%sqet**2 + beta16
      do i=s%nz1,s%nz0
         stemp5 = (stemp1 + s%ionization_temp_over_t(i))* &
              (1.0d0 - 2.0d0*s%species_ion_fraction(i))
         stemp6 = s%species_weighted_term(i)*(stemp4 - s%ionization_temp_over_t(i)) &
              + d%species_temp_deriv_term(i)*stemp5
         r3t = r3t - stemp6
         utsum = utsum + ionization_temp(i)*(stemp6 - d%species_temp_deriv_term(i))
         stemp7 = s%species_weighted_term(i)*(stemp3 + stemp2*stemp5)
         r3p = r3p - stemp7
         upsum = upsum + ionization_temp(i)*stemp7
      end do
      end if
      if (.not. s%skip_helium_i) then
      stemp = 2.0d0*s%helium_ion_fraction_1*s%ep1*d%sqet
      r2t = s%helium_saha_ratio_1*(d%helium_neutral_fraction*(beta16- &
           helium_ionization_temp_1*s%temperature_inverse+d%sk0qt**2)+2d0*d%sx0qt*d%sk0qt) &
           +s%ep12*(stemp*d%sqet - 2.0d0*d%sx1qt*d%sqet)
      r2p = s%helium_saha_ratio_1*(d%helium_neutral_fraction*r%beta_inverse* &
           (-r%beta14-d%sk0qt)+d%sx0qp*d%sk0qt-r%beta_inverse* &
           d%sx0qt) + s%ep12*(stemp*d%sqep - d%sx1qp*d%sqet - d%sx1qt*d%sqep)
      if (.not. s%skip_helium_ii) then
      stemp = 2.0d0*s%helium_ion_fraction_2*s%ep1*d%sqet
      r1t = s%helium_saha_ratio_2*(s%helium_ion_fraction_1*(beta16- &
           helium_ionization_temp_2*s%temperature_inverse+d%sk1qt**2)+2d0*d%sx1qt*d%sk1qt) &
           +s%ep12*(stemp*d%sqet - 2.0d0*d%sx2qt*d%sqet)
      r1p = s%helium_saha_ratio_2*(s%helium_ion_fraction_1*r%beta_inverse* &
           (-r%beta14-d%sk1qt)+d%sx1qp*d%sk1qt-r%beta_inverse* &
           d%sx1qt) + s%ep12*(stemp*d%sqep - d%sx2qp*d%sqet - d%sx2qt*d%sqep)
      r2t = r2t + s%cr12*r1t
      r2p = r2p + s%cr12*r1p
      r3p = r3p + s%cr13*r1p
      r3t = r3t + s%cr13*r1t
      end if
      r3p = r3p + s%cr23*r2p
      r3t = r3t + s%cr23*r2t
      end if
      sqett = r3t/s%c33
      sqetp = r3p/s%c33
      r%dlnrho_dlnt_dt = -beta16 + (s%ep1*d%sqet)**2 - s%ep1*sqett
      r%dlnrho_dlnt_dp = r%beta_inverse*r%beta14 + s%ep12*d%sqet*d%sqep - s%ep1*sqetp
      if (s%nz0 > 0) then
      stemp = s%ep1e*sqett
      stemp1 = s%ep1e*sqetp
      do i=s%nz1,s%nz0
         utsum = utsum - s%species_weighted_term(i)*ionization_temp(i)*stemp
         upsum = upsum - s%species_weighted_term(i)*ionization_temp(i)*stemp1
      end do
      end if
      if (.not. s%skip_helium_i) then
      sx1qtt = (r2t - s%c23*sqett)/s%c22
      sx1qtp = (r2p - s%c23*sqetp)/s%c22
      if (.not. s%skip_helium_ii) then
      sx2qtt = (r1t - s%c12*sx1qtt - s%c13*sqett)/s%c11
      sx2qtp = (r1p - s%c12*sx1qtp - s%c13*sqetp)/s%c11
      utsum = utsum + s%helium_mass_fraction*(helium_ionization_temp_1+ &
           helium_ionization_temp_2)*(sx2qtt - d%sx2qt)
      upsum = upsum + s%helium_mass_fraction*(helium_ionization_temp_1+ &
           helium_ionization_temp_2)*sx2qtp
      end if
      utsum = utsum + s%helium_mass_fraction*helium_ionization_temp_1*(sx1qtt - d%sx1qt)
      upsum = upsum + s%helium_mass_fraction*helium_ionization_temp_1*sx1qtp
      end if
      btemp = (1.0d0+s%mean_electrons_per_ion)*r%beta14*r%beta_inverse* &
           (24.0d0*r%beta_inverse - 9.0d0)
      qqutt = d%ramu*(btemp+0.75d0*(2.0d0+r%beta14*(8.0d0*r%beta_inverse+1d0))*d%sqet+ &
           d%betaut*sqett + utsum*s%temperature_inverse)
      qqutp = d%ramu*(-0.25d0*btemp+0.75d0*((2.0d0+r%beta14*(4.0d0*r%beta_inverse &
           +1.0d0))*d%sqep-r%beta14*r%beta_inverse*d%sqet) +d%betaut*sqetp+upsum*s%temperature_inverse)
      r%specific_heat_cp_dt = d%qcpi*(qqutt + d%pdtq*(1.0d0+r%dlnrho_dlnt) - d%pdt*r%dlnrho_dlnt_dt)
      r%specific_heat_cp_dp = d%qcpi*(qqutp - d%pdtq*(1.0d0-r%dlnrho_dlnp) - d%pdt*r%dlnrho_dlnt_dp)
      r%adiabatic_gradient_dt = -1.0d0 - r%dlnrho_dlnt - r%specific_heat_cp_dt + &
           r%dlnrho_dlnt_dt/r%dlnrho_dlnt
      r%adiabatic_gradient_dp = 1.0d0 - r%dlnrho_dlnp - r%specific_heat_cp_dp + &
           r%dlnrho_dlnt_dp/r%dlnrho_dlnt

end subroutine differentiate_saha_solution_twice

!----------------------------------------------------------------------
! compute_fully_ionized_eos (formerly eqrelv)
!----------------------------------------------------------------------
! Fully-ionized-gas equation of state ("relativistic" electron gas via
! the tabulated Fermi-Dirac integrals in yale_eos_tables's
! yale%fermi_table_data, following Prather's thesis). Given log10(T), log10(P), and the mean molecular
! weights, iterates on density until the electron+ion+radiation
! pressure sum matches the given P, then returns density and (if
! in_atmosphere is false, via blend_eos.f90's call) thermodynamic
! derivatives. Called from blend_eos.f90 both for the T >= saha cutoff
! branch (fully ionized gas assumed) and to interpolate against the
! compute_saha_eos/compute_scv_envelope_eos result near the ionization cutoff.
!
!          Qdtt,qdtp = nat-log derivatives of qdt=d(log d )/d(log T )
!          qat,qap = nat-log derivatives of the adiabatic derivative
! DBG 7/92 modified to include Debye-Hueckel correction.
!     Correction is applied to pressure and all thermodynamic variables.
!     Treatment is also valid for electron degeneracy.  A ramp function
!     between the nondegenerate and degernate electron gas cases of the
!     D-H correction is used (this is not correct but, hopefully is
!     close enough.)  The constant term is debye_huckel_coefficient.
subroutine compute_fully_ionized_eos(log10_temperature, log10_pressure, r)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(eos_result), intent(inout) :: r
      type(fermi_solution) :: f
      logical :: converged

      call start_fully_ionized_solution(log10_temperature, log10_pressure, r, f)
      call iterate_fully_ionized_density(log10_temperature, log10_pressure, r, f, converged)
      if (.not. converged) return
      r%log10_density = f%dl8
      r%density = f%d8
      call interpolate_fermi_table(f)
      call differentiate_pressure(f)
      call differentiate_density(f, r)
      call compute_internal_energy_and_cp(f, r)
      r%electron_degeneracy_parameter = f%cl1*yale%fermi_table_eta(f%id1) + &
           f%cl2*yale%fermi_table_eta(f%id2) + f%cl3*yale%fermi_table_eta(f%id3)

end subroutine compute_fully_ionized_eos

!----------------------------------------------------------------------
! start_fully_ionized_solution
!----------------------------------------------------------------------
! The pressures that do not depend on density (radiation, and the
! Debye-Hueckel mass terms), the log T row of the Fermi table and its
! interpolation factors, and the starting log density.
subroutine start_fully_ionized_solution(log10_temperature, log10_pressure, r, f)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(eos_result), intent(in) :: r
      type(fermi_solution), intent(inout) :: f
      real(dp) :: almix
      integer :: ii
      real(dp) :: ramu

      f%pinv = 1.0d0/r%pressure
      f%beta1 = 1.0d0 - r%beta
      f%tl8 = log10_temperature
      f%pl8 = log10_pressure
      f%pr = f%beta1*r%pressure
      ramu = gas_constant*r%ion_mean_weight_inverse
      f%ramut = ramu*r%temperature
      f%emul = log10(r%electron_mean_weight_inverse)
      f%emul8 = f%emul
!  Pe is tabulated as a function of X(=log d - log mue - 3/2log T) and
!  Y(= log T). Find index(JT1,JT2,JT3) for 3 pt interpolation in Y.
! MHP 10/02 corrected mixed type
      f%jt1 = int(5.0d0*(log10_temperature - 5.0d0))
      f%jt1 = max0(1,min0(f%jt1,18))  ! Force jt1 to be in the table
      f%jt2 = f%jt1 + 1
      f%jt3 = f%jt1 + 2
      f%ttud = 5.0d0*(f%tl8 - (5.0d0 + 0.20d0*real(f%jt1, dp)))
      f%ttcu = 0.50d0*f%ttud*(f%ttud - 1.0d0)
      f%ttdcu = f%ttud - 0.50d0
      f%dml = r%log10_density + f%emul
      f%nden = 1
      f%id1p = -100
! DBG 7/92 calculate D-H terms: Alpha-mixture, mass-FRAC0, mass-FRAC1
      if (eos_settings%use_debye_huckel_correction) then
         almix = 0.0d0
         do ii=1,18
            almix = almix + debye_huckel_nu(ii)*debye_huckel_z(ii)
         end do
         f%cmfdh0 = 2.0d0*debye_huckel_x+1.5d0*debye_huckel_y+0.5d0*(debye_huckel_z_total+almix)
         f%cmfdh1 = debye_huckel_x+debye_huckel_y+0.5d0*almix
      end if
end subroutine start_fully_ionized_solution

!----------------------------------------------------------------------
! iterate_fully_ionized_density
!----------------------------------------------------------------------
! Newton iteration on log rho until the electron (Fermi table), ion and
! radiation pressures, plus the Debye-Hueckel term, add up to P; steps
! are capped at 4 dex, and 20 steps without convergence is a failure
! (converged false, a message in the run log). For details consult
! Prather's thesis. The electron pressure pe is tabulated against
! X = log rho - log mu_e - 1.5 log T (3-point Lagrange, unequal spacing)
! and Y = log T (3-point Newton, equal spacing); pa is the ion pressure,
! pr the radiation pressure, pdh the Debye-Hueckel term.
subroutine iterate_fully_ionized_density(log10_temperature, log10_pressure, r, f, converged)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(eos_result), intent(inout) :: r
      type(fermi_solution), intent(inout) :: f
      logical, intent(out) :: converged
      real(dp) :: cden1
      real(dp) :: cden2
      real(dp) :: cden3
      real(dp) :: cmfdh
      real(dp) :: corr
      real(dp) :: dpel1
      real(dp) :: dpel2
      real(dp) :: dx1
      real(dp) :: dx2
      real(dp) :: dx3
      real(dp) :: fact
      real(dp) :: factl
      real(dp) :: factor
      integer :: id
      integer :: kk
      integer :: nn
      real(dp) :: pel
      real(dp) :: pr6
      real(dp) :: pr7
      real(dp) :: pr8
      real(dp) :: pr8l
      real(dp) :: pr9
      real(dp) :: pr9l
      real(dp) :: pt
      real(dp) :: ptl
      real(dp) :: slope
      real(dp) :: t0
      real(dp) :: t1
      real(dp) :: xx

      converged = .true.
      density_iter: do
! Begin iteration loop for correct density (was label 10)
!  find index (ID1,ID2,ID3) for 3-pt interpolation in X
      xx = f%dml - 1.50d0*f%tl8
      pr6= 20.0d0*(xx - yale%fermi_table_x_grid(1)) + 1.0d0
      kk = int(pr6)
      kk = min0(261,max0(1,kk))
      f%id1 = yale%fermi_table_x_lookup(kk)
!  If index unchanged from previous loop,skip this section
      if (f%id1 /= f%id1p) then
      f%id1p = f%id1
      f%id2 = f%id1 + 1
      f%id3 = f%id1 + 2
      pr9 = 1.0d0/( yale%fermi_table_x_grid(f%id3) - yale%fermi_table_x_grid(f%id2) )
      pr8 = 1.0d0/( yale%fermi_table_x_grid(f%id2) - yale%fermi_table_x_grid(f%id1) )
      pr7 = 1.0d0/( yale%fermi_table_x_grid(f%id1) - yale%fermi_table_x_grid(f%id3) )
      cden1 = -pr7*pr8
      cden2 = -pr8*pr9
      cden3 = -pr9*pr7
!  Interpolation in Y is Newtonian(equal spacing in Y)
      do nn = 1,3
       id = f%id1 + nn - 1
       f%df1(1,nn) = yale%fermi_table_data(1,id,f%jt2) - yale%fermi_table_data(1,id,f%jt1)
       f%df2(1,nn) = yale%fermi_table_data(1,id,f%jt3)-yale%fermi_table_data(1,id,f%jt2)-f%df1(1,nn)
       f%ff(1,nn) = yale%fermi_table_data(1,id,f%jt1) + f%ttud*f%df1(1,nn) + f%ttcu*f%df2(1,nn)
      end do
      dpel2 = 2.0d0*(f%ff(1,1)*cden1 + f%ff(1,2)*cden2 + f%ff(1,3)*cden3)
      end if
      dx1 = xx - yale%fermi_table_x_grid(f%id1)
      dx2 = xx - yale%fermi_table_x_grid(f%id2)
      dx3 = xx - yale%fermi_table_x_grid(f%id3)
!  Interpolation in X is Lagrangian(unequal table spacing in X)
      f%cl1 = dx2*dx3*cden1
      f%cl2 = dx3*dx1*cden2
      f%cl3 = dx1*dx2*cden3
      f%dcl1 = (dx2 + dx3)*cden1
      f%dcl2 = (dx3 + dx1)*cden2
      f%dcl3 = (dx1 + dx2)*cden3
      if (eos_settings%use_debye_huckel_correction) then
         r%electron_degeneracy_parameter = f%cl1*yale%fermi_table_eta(f%id1) + &
              f%cl2*yale%fermi_table_eta(f%id2) + f%cl3*yale%fermi_table_eta(f%id3)
         if (r%electron_degeneracy_parameter < eos_settings%debye_huckel_eta_min) then
            cmfdh = f%cmfdh0*sqrt(f%cmfdh0)
         else if (r%electron_degeneracy_parameter > eos_settings%debye_huckel_eta_max) then
            cmfdh = f%cmfdh0*sqrt(f%cmfdh1)
         else
            t0 = f%cmfdh0*sqrt(f%cmfdh0)
            t1 = f%cmfdh0*sqrt(f%cmfdh1)
            slope = (t1-t0)/(eos_settings%debye_huckel_eta_max-eos_settings%debye_huckel_eta_min)
            cmfdh = t0+(r%electron_degeneracy_parameter-eos_settings%debye_huckel_eta_min)*slope
         end if
      end if
!  Now find pe and its 1st and 2nd derivatives w.r.t. density
      pel =  f%cl1*f%ff(1,1) +  f%cl2*f%ff(1,2) +  f%cl3*f%ff(1,3)
      dpel1 = f%dcl1*f%ff(1,1) + f%dcl2*f%ff(1,2) + f%dcl3*f%ff(1,3)
      f%dl8 = f%dml - f%emul8
      f%d8 = exp(ln10*f%dl8)
      f%pa = f%ramut*f%d8
      f%pe = exp(ln10*pel)
      if (eos_settings%use_debye_huckel_correction) then
         f%pdh = debye_huckel_coefficient*cmfdh*sqrt(f%d8/r%temperature)*f%d8
      end if
      pt = f%pa + f%pe + f%pr
      if (eos_settings%use_debye_huckel_correction) then
         pt = pt + f%pdh
         f%pdhp = f%pdh/pt
      end if
      ptl = log10(pt)
      pr7 = (f%pa + f%pe*dpel1)/pt
      if (eos_settings%use_debye_huckel_correction) then
         pr7 = pr7 + 1.5d0*f%pdhp
      end if
      pr8 = (ptl - f%pl8)/pr7
      pr9 = ln10*( (f%pa + f%pe*(dpel1**2+dpel2/ln10))/(pt*pr7) - pr7 )
      if (eos_settings%use_debye_huckel_correction) then
         pr9 = pr9 + 2.25d0*f%pdhp/pr7
      end if
! Factor is the difference between the P predicted from the given rho
! and the true P; if too large increment rho and try again.
      factor=abs(pr8)
      if(factor < 1d-15) factor = 1d-15
      pr8l=log10(factor)
      pr9l=log10(abs(pr9))
      factl=log10(0.5d0)+pr9l+2.d0*pr8l
      fact=exp(ln10*factl)
      fact=sign(fact,pr9)
      corr=-pr8-fact
!  limit range of changes in density to a factor of 10**4   llp  1/31/07
      if(corr > 4d0) corr =  4d0
      if(corr < -4d0) corr = -4d0
      if(abs(corr) >= 1.0d-08) then
       f%dml = f%dml + corr
       f%nden = f%nden + 1
       if(f%nden <= 20) cycle density_iter
       write(run_log_unit,40) log10_temperature,log10_pressure,ptl,f%dml,corr
   40    format('EQRELV: Did not Converge: T,P,Pcalc,Dcalc,CORR', &
                4F10.6,F20.12)
       converged = .false.
       return
      end if
      exit density_iter
      end do density_iter
end subroutine iterate_fully_ionized_density

!----------------------------------------------------------------------
! interpolate_fermi_table
!----------------------------------------------------------------------
! The other Fermi-table quantities (rows 2-5) at the converged point.
subroutine interpolate_fermi_table(f)
      type(fermi_solution), intent(inout) :: f
      integer :: id
      integer :: kk
      integer :: nn

      do kk = 2,5
       do nn = 1,3
          id = f%id1 + nn - 1
          f%df1(kk,nn) = yale%fermi_table_data(kk,id,f%jt2) - yale%fermi_table_data(kk,id,f%jt1)
          f%df2(kk,nn)=(yale%fermi_table_data(kk,id,f%jt3)-yale%fermi_table_data(kk,id,f%jt2))-f%df1(kk,nn)
          f%ff(kk,nn)=yale%fermi_table_data(kk,id,f%jt1)+f%ttud*f%df1(kk,nn)+f%ttcu*f%df2(kk,nn)
       end do
       f%ffe(kk) = f%cl1*f%ff(kk,1) + f%cl2*f%ff(kk,2) + f%cl3*f%ff(kk,3)
      end do
end subroutine interpolate_fermi_table

!----------------------------------------------------------------------
! differentiate_pressure
!----------------------------------------------------------------------
! The log derivatives of P in T and rho, first and second: qpt =
! dln P/dln T, qpd = dln P/dln rho (the reciprocal of dlnrho_dlnp, kept
! as its own quantity as in the original), qqptt, qqpdt, qqpdd.
subroutine differentiate_pressure(f)
      type(fermi_solution), intent(inout) :: f
      real(dp) :: peq
      real(dp) :: pr1
      real(dp) :: pr2
      real(dp) :: pr3
      real(dp) :: pr4
      real(dp) :: pr5
      real(dp) :: qped
      real(dp) :: qpeqd
      real(dp) :: qpeqt
      real(dp) :: qqpedd

! Derivatives of P
      peq = exp(ln10*f%ffe(3))
      qped = f%ffe(4)
      qqpedd =(f%dcl1*f%ff(4,1) + f%dcl2*f%ff(4,2) + f%dcl3*f%ff(4,3))*clni
      pr1 = f%df1(3,1)+f%ttdcu*f%df2(3,1)
      pr2 = f%df1(3,2)+f%ttdcu*f%df2(3,2)
      pr3 = f%df1(3,3)+f%ttdcu*f%df2(3,3)
      pr4 = f%dcl1*f%ff(3,1) + f%dcl2*f%ff(3,2) + f%dcl3*f%ff(3,3)
      pr5 =   5.0d0*(f%cl1*pr1 + f%cl2*pr2 + f%cl3*pr3)
      qpeqd = peq*pr4
      qpeqt = peq*(pr5 - 1.50d0*pr4)
      f%qpt = (f%pa + peq)*f%pinv + 4.0d0*f%beta1
      if (eos_settings%use_debye_huckel_correction) then
         f%qpt = f%qpt - 0.5d0*f%pdhp
      end if
      f%qpd = (f%pa + f%pe*qped)*f%pinv
      if (eos_settings%use_debye_huckel_correction) then
         f%qpd = f%qpd + 1.5d0*f%pdhp
      end if
      f%qqptt = (f%pa + qpeqt)*f%pinv + 16.0d0*f%beta1 - f%qpt**2
      if (eos_settings%use_debye_huckel_correction) then
         f%qqptt = f%qqptt+(0.5d0*f%qpt+0.25d0)*f%pdhp
      end if
      f%qqpdt = (f%pa + qpeqd)*f%pinv - f%qpt*f%qpd
      if (eos_settings%use_debye_huckel_correction) then
         f%qqpdt = f%qqpdt+(0.5d0*f%qpt-0.75d0)*f%pdhp
      end if
      f%qqpdd = (f%pa + f%pe*(qped**2 + qqpedd))*f%pinv - f%qpd**2
      if (eos_settings%use_debye_huckel_correction) then
        f%qqpdd = f%qqpdd + (-1.5d0*f%qpd+2.25d0)*f%pdhp
      end if
end subroutine differentiate_pressure

!----------------------------------------------------------------------
! differentiate_density
!----------------------------------------------------------------------
! dln rho/dln T, dln rho/dln P and their T and P derivatives, from the
! pressure derivatives.
subroutine differentiate_density(f, r)
      type(fermi_solution), intent(in) :: f
      type(eos_result), intent(inout) :: r
      real(dp) :: pr1
      real(dp) :: qdpp

! Derivatives of density
      r%dlnrho_dlnp = 1.0d0/f%qpd
      r%dlnrho_dlnt = -f%qpt*r%dlnrho_dlnp
      pr1 = r%dlnrho_dlnp**2
      qdpp = -f%qqpdd*r%dlnrho_dlnp*pr1
      r%dlnrho_dlnt_dp = -f%qpt*qdpp - f%qqpdt*pr1
      r%dlnrho_dlnt_dt = -f%qpt*r%dlnrho_dlnt_dp - r%dlnrho_dlnp*(f%qqptt + r%dlnrho_dlnt*f%qqpdt)
end subroutine differentiate_density

!----------------------------------------------------------------------
! compute_internal_energy_and_cp
!----------------------------------------------------------------------
! The internal energy and its derivatives, cp and grad_ad, and their
! log derivatives in T and P.
subroutine compute_internal_energy_and_cp(f, r)
      type(fermi_solution), intent(in) :: f
      type(eos_result), intent(inout) :: r
      real(dp) :: acp
      real(dp) :: bcp
      real(dp) :: pd
      real(dp) :: pdt
      real(dp) :: pdu
      real(dp) :: pr1
      real(dp) :: pr2
      real(dp) :: pr3
      real(dp) :: pr4
      real(dp) :: pr5
      real(dp) :: qcpi
      real(dp) :: qqudt
      real(dp) :: qqutt
      real(dp) :: qudt
      real(dp) :: queqt
      real(dp) :: qutd
      real(dp) :: qutp
      real(dp) :: u
      real(dp) :: ua
      real(dp) :: udh
      real(dp) :: udhu
      real(dp) :: ue
      real(dp) :: ueq
      real(dp) :: ui
      real(dp) :: ur

! Internal energy term(u) and its derivatives
      ue  = exp(ln10*(f%ffe(2) + f%emul))
      ueq = exp(ln10*(f%ffe(5) + f%emul))
      pr1 = f%df1(5,1)+f%ttdcu*f%df2(5,1)
      pr2 = f%df1(5,2)+f%ttdcu*f%df2(5,2)
      pr3 = f%df1(5,3)+f%ttdcu*f%df2(5,3)
      pr4 = f%dcl1*f%ff(5,1) + f%dcl2*f%ff(5,2) + f%dcl3*f%ff(5,3)
      pr5 =   5.0d0*(f%cl1*pr1 + f%cl2*pr2 + f%cl3*pr3)
      queqt = ueq*(pr5 - 1.50d0*pr4)
      pd = r%pressure/r%density
      pdt = pd/r%temperature
      ua = 1.50d0*f%ramut
      ur = 3.0d0*f%beta1*pd
      if (eos_settings%use_debye_huckel_correction) then
         udh = 3.0d0*f%pdh/f%d8
      end if
      u = ua + ue + ur
      if (eos_settings%use_debye_huckel_correction) then
         u = u + udh
         udhu = udh/u
      end if
      ui = 1.0d0/u
      pdu = pd*ui
      qutd = (ua + ueq + 4.0d0*ur)*ui
      if(eos_settings%use_debye_huckel_correction) then
         qutd = qutd -0.5d0*udhu
      end if
      qudt = pdu*(1.0d0 - f%qpt)
      if (eos_settings%use_debye_huckel_correction) then
         qudt = qudt + 0.5d0*udhu
      end if
      qutp = qutd + r%dlnrho_dlnt*qudt
      qqutt = (ua + queqt + 16.0d0*ur)*ui - qutd**2
      if (eos_settings%use_debye_huckel_correction) then
         qqutt = qqutt + (0.5d0*qutd+0.25d0)*udhu
      end if
      qqudt = qudt*(f%qpt - qutd) - pdu*f%qqptt
      if (eos_settings%use_debye_huckel_correction) then
         qqudt = qqudt + (0.5d0*qudt-0.25d0)*udhu
      end if
! Specific heat(qcp) and its derivatives
      pr1 =-pdt*r%dlnrho_dlnt
      pr2 = u/r%temperature
      pr3 = -(1.0d0 + r%dlnrho_dlnt) +  r%dlnrho_dlnt_dt/r%dlnrho_dlnt
      pr4 = (1.0d0 - r%dlnrho_dlnp) +  r%dlnrho_dlnt_dp/r%dlnrho_dlnt
      acp = pr1*f%qpt
      bcp = pr2*qutd
      r%specific_heat_cp = acp + bcp
      qcpi = 1.0d0/r%specific_heat_cp
      r%specific_heat_cp_dt = acp*pr3 + pr1*(f%qqptt + r%dlnrho_dlnt*f%qqpdt) + &
              bcp*(qutp - 1.0d0) + pr2*(qqutt + r%dlnrho_dlnt*qqudt)
      r%specific_heat_cp_dp = acp*pr4 + r%dlnrho_dlnp*(pr1*f%qqpdt + bcp*qudt + pr2*qqudt)
      r%specific_heat_cp_dt = r%specific_heat_cp_dt*qcpi
      r%specific_heat_cp_dp = r%specific_heat_cp_dp*qcpi
! Adiabatic gradient(dela) and derivatives ( all derivs are ln(dela)/ln(
      r%adiabatic_gradient = pr1*qcpi
      r%adiabatic_gradient_dt = pr3 - r%specific_heat_cp_dt
      r%adiabatic_gradient_dp = pr4 - r%specific_heat_cp_dp
end subroutine compute_internal_energy_and_cp

!----------------------------------------------------------------------
! read_fermi_table
!----------------------------------------------------------------------
! Reads the Fermi-Dirac integral table of the degenerate electron gas
! (FERMI.TAB: 43 abscissas, 43 degeneracies, 860 cards of 5 values)
! into yale, and bins the abscissas into fermi_table_x_lookup, the
! interval index compute_fully_ionized_eos looks up. A non-increasing
! abscissa is a corrupt table: ierr = 1.
subroutine read_fermi_table(fermi_table_path, ierr)
      character(len=*), intent(in) :: fermi_table_path
      integer, intent(out) :: ierr
      integer :: fermi_unit, grid_idx, iden_idx, col_idx, card_idx
      integer :: bin_start, bin_width, bin_end

      ierr = 0
      open(newunit=fermi_unit, file=fermi_table_path, form='FORMATTED', status='OLD', &
           iostat=ierr)
      if (ierr /= 0) then
         write(run_log_unit,*) 'eos_init: cannot open the Fermi-Dirac table ', &
              trim(fermi_table_path)
         ierr = 1
         return
      end if
      read(fermi_unit,150)  (yale%fermi_table_x_grid(grid_idx), grid_idx = 1,43)
      read(fermi_unit,150)  (yale%fermi_table_eta(grid_idx), grid_idx = 1,43)
      read(fermi_unit,160) (grid_idx,iden_idx,(yale%fermi_table_data(col_idx,grid_idx,iden_idx), &
           col_idx=1,5),card_idx=1,860)
  150 format(8x,9f8.4)
  160 format(2i5,5f12.7)
      bin_start = 1
      do grid_idx=1,41
! a non-increasing abscissa is a corrupt table
       if (yale%fermi_table_x_grid(grid_idx+1) <= yale%fermi_table_x_grid(grid_idx)) then
          write(run_log_unit,1000) grid_idx
 1000       format(1x,39('>'),40('<')/1x,'read_fermi_table: abscissa',i4, &
           ' does not increase')
          close(fermi_unit)
          ierr = 1
          return
       endif
       bin_width = int((yale%fermi_table_x_grid(grid_idx+1) - yale%fermi_table_x_grid(grid_idx))*20.0d0 + 0.10d0)
       bin_end = bin_start + bin_width - 1
       if (grid_idx == 41)  bin_end = 261
       do iden_idx=bin_start,bin_end
          yale%fermi_table_x_lookup(iden_idx) = grid_idx
       end do
       bin_start = bin_end + 1
      end do

      close(fermi_unit)
end subroutine read_fermi_table

end module yale_eos
