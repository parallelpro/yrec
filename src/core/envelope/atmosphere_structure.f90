!----------------------------------------------------------------------
! atmosphere_structure
!----------------------------------------------------------------------
! The atmosphere of the current model, point by point (1..num_atm_points):
! log P, log T, log rho, the depth step, the gradients, beta, gamma1,
! d ln rho/d ln T, the ionization fractions, the opacity and its
! derivatives, and cp. The atmosphere integration
! (envelope_integration_kernel) fills it; the stitched model (and its
! pulse files) reads it. Also here: atm_point, the current point of that
! integration.
module atmosphere_structure
      use star_info_lib, only: json
      implicit none
      private

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
! point (atmosphere_derivs: pt_scr%qqod, pt_scr%qqot), for the pulse
! writers' kap_kap_T and kap_kap_rho columns (see envelope_structure).
            double precision :: atmo_dlnkap_dlnrho(json), atmo_dlnkap_dlnt(json)
            integer :: num_atm_points
      end type atmosphere_structure_state

      type(atmosphere_structure_state), public :: atmo_struct

! The current point of the atmosphere integration: log10 tau, log P,
! log T, log rho, the opacity and the ionization fractions. The surface
! table lookup fills log P and log T; atmosphere_derivs fills the rest at
! each step; the integration and the envelope start read it.
      type, public :: atmosphere_point_state
            double precision :: log10_tau, log10_pressure, log10_temperature, &
                 log10_density, opacity
            double precision :: ion_fraction(3)
      end type atmosphere_point_state

      type(atmosphere_point_state), public :: atm_point

end module atmosphere_structure
