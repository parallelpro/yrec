!----------------------------------------------------------------------
! atmstruct_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/atmstruct/: the per-point atmosphere
! structure profile (pressure/temperature/density/depth/gradients/
! ionization/opacity/specific heat in log space, indexed
! 1..num_atm_points), the atmosphere-region counterpart of
! envstruct_lib's envelope profile. Computed by the atmosphere
! integrator (atm/atm_lib.f90) and read by output writers
! (misc/stitch.f90, io/putstore.f90, io/wrtmod.f90, io/wrtout.f90) and
! rotation/ files needing the current atmosphere structure.
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! envstruct_lib/scrtch_lib/mdphy_lib: genuinely evolving per-model
! state, read from distant, unrelated points in the call graph. Field
! names are unchanged from the original COMMON member names, matching
! that same precedent.
module atmstruct_lib
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
            integer :: num_atm_points
      end type atmosphere_structure_state

      type(atmosphere_structure_state), public :: atmo_struct

end module atmstruct_lib
