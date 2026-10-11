!----------------------------------------------------------------------
! envelope_structure
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! docs/GUIDELINES.md). Replaces common/envstruct/: the per-point envelope
! structure profile (pressure/temperature/mass/density/radius/
! composition/convective-flag/gradients/velocity/opacity/luminosity in
! log space, indexed 1..num_env_points), computed by the envelope
! integrator (core/envint_kernel.f90, under core/envint_lib.f90) and
! read by core/read_starting_model.f90 and other callers that need the
! just-integrated envelope profile.
!
! Per docs/GUIDELINES.md's module-vs-argument test this is case 1b, same as
! scrtch_lib/mdphy_lib/temp_lib: genuinely evolving per-model state,
! read from distant, unrelated points in the call graph. Field names
! are unchanged from the original COMMON member names, matching that
! same precedent.
module envelope_structure
      use star_info_lib, only: json
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
! 2026 (pulse derivatives): the opacity derivatives the envelope
! integration computes at every point (envelope_derivs: pt_scr%qqod,
! pt_scr%qqot), kept so that the pulse writers can fill the GYRE/GSM
! kap_kap_T and kap_kap_rho columns of the envelope points (they were
! written as zero, which switches off the kappa mechanism of the
! hydrogen ionisation zone in any nonadiabatic calculation on the file).
            double precision :: env_dlnkap_dlnrho(json), env_dlnkap_dlnt(json)
            integer :: num_env_points
      end type envelope_structure_state

      type(envelope_structure_state), public :: env_struct

end module envelope_structure
