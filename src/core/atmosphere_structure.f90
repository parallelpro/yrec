!----------------------------------------------------------------------
! atmosphere_structure
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! docs/GUIDELINES.md). Replaces common/atmstruct/: the per-point atmosphere
! structure profile (pressure/temperature/density/depth/gradients/
! ionization/opacity/specific heat in log space, indexed
! 1..num_atm_points), the atmosphere-region counterpart of
! envelope_structure's envelope profile. Computed by the atmosphere
! integrator (core/envelope_integration_kernel.f90, under core/envelope_integration.f90) and
! read by the stitched-model/profile writers (core/stitched_model.f90
! and its callers) needing the current atmosphere structure.
!
! Per docs/GUIDELINES.md's module-vs-argument test this is case 1b, same as
! envelope_structure/scrtch_lib/mdphy_lib: genuinely evolving per-model
! state, read from distant, unrelated points in the call graph. Field
! names are unchanged from the original COMMON member names, matching
! that same precedent.
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
! 2026 (pulse derivatives): opacity derivatives at every atmosphere
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
