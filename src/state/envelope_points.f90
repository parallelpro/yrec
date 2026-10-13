!----------------------------------------------------------------------
! envelope_points
!----------------------------------------------------------------------
! The records of the outer layers that the envelope integration
! (core/envelope/envelope_integration_kernel.f90) fills, and the star
! holds (star_info_lib: s%env_points, s%atm_points, s%atm_point):
!   envelope_structure_state   the envelope, point by point
!     (1..num_env_points): log P, log T, log m, log rho, log r, X, Z, the
!     convective flag, the gradients, the convective velocity, beta,
!     gamma1, cp, the ionization fractions, the opacity and its
!     derivatives (for the pulse writers' kap_kap_T and kap_kap_rho
!     columns), the luminosity and d ln rho/d ln T. read_starting_model,
!     envelope_refit, rebuild_envelope and the stitched model read it.
!   atmosphere_structure_state the atmosphere, point by point
!     (1..num_atm_points): log P, log T, log rho, the depth step, the
!     gradients, beta, gamma1, d ln rho/d ln T, the ionization
!     fractions, the opacity and its derivatives, and cp. The stitched
!     model reads it.
!   atmosphere_point_state     the atmosphere integration's current point
!     (log tau, log P, log T, log rho, the opacity, the ionization
!     fractions): the integrand writes it, the kernel stores it, and the
!     mass change falls back on its last values when the Allard lookup
!     fails.
! Types only; star-blind code (the kernel, the integrands) takes them
! from here rather than from star_info_lib.
module envelope_points
      use composition_layout, only: json
      implicit none
      private

      type, public :: envelope_structure_state
            double precision :: env_log10_pressure(json), &
                 env_log10_temperature(json)
            double precision :: env_log10_mass(json), env_log10_density(json), &
                 env_log10_radius(json)
            double precision :: env_hydrogen_fraction(json), &
                 env_metal_fraction(json)
            logical :: env_convective_flag(json)
            double precision :: env_gradients(3,json), &
                 env_convective_velocity(json), env_beta(json)
            double precision :: env_gamma1(json), env_specific_heat_cp(json), &
                 env_ion_fraction(3,json)
            double precision :: env_opacity(json), env_luminosity(json), &
                 env_dlnrho_dlnt(json)
! the opacity derivatives the envelope
! integration computes at every point (envelope_derivs: pt%qqod,
! pt%qqot), kept so that the pulse writers can fill the GYRE/GSM
! kap_kap_T and kap_kap_rho columns of the envelope points (they were
! written as zero, which switches off the kappa mechanism of the
! hydrogen ionisation zone in any nonadiabatic calculation on the file).
            double precision :: env_dlnkap_dlnrho(json), env_dlnkap_dlnt(json)
            integer :: num_env_points
      end type envelope_structure_state

      type, public :: atmosphere_structure_state
            double precision :: atmo_log10_pressure(json), &
                 atmo_log10_temperature(json)
            double precision :: atmo_log10_density(json), &
                 atmo_delta_depth(json)
            double precision :: atmo_gradients(3,json), atmo_beta(json)
            double precision :: atmo_gamma1(json), atmo_dlnrho_dlnt(json), &
                 atmo_ion_fraction(3,json)
            double precision :: atmo_opacity(json), atmo_specific_heat_cp(json)
! opacity derivatives at every atmosphere
! point (atmosphere_derivs: pt%qqod, pt%qqot), for the pulse
! writers' kap_kap_T and kap_kap_rho columns (see envelope_structure_state).
            double precision :: atmo_dlnkap_dlnrho(json), atmo_dlnkap_dlnt(json)
            integer :: num_atm_points
      end type atmosphere_structure_state

      type, public :: atmosphere_point_state
            double precision :: log10_tau, log10_pressure, log10_temperature, &
                 log10_density, opacity
            double precision :: ion_fraction(3)
      end type atmosphere_point_state

end module envelope_points
