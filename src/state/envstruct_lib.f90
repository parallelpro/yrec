!----------------------------------------------------------------------
! envstruct_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/envstruct/: the per-point envelope
! structure profile (pressure/temperature/mass/density/radius/
! composition/convective-flag/gradients/velocity/opacity/luminosity in
! log space, indexed 1..num_env_points), computed by the envelope
! integrator (atm/envint.f90) and read by core/starin.f90 and other
! callers that need the just-integrated envelope profile.
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! scrtch_lib/mdphy_lib/temp_lib: genuinely evolving per-model state,
! read from distant, unrelated points in the call graph. Field names
! are unchanged from the original COMMON member names, matching that
! same precedent.
module envstruct_lib
      implicit none
      private
      integer, parameter :: json = 5000

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
            integer :: num_env_points
      end type envelope_structure_state

      type(envelope_structure_state), public :: env_struct

end module envstruct_lib
