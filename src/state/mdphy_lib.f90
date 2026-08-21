!----------------------------------------------------------------------
! mdphy_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/mdphy/: per-shell auxiliary physics
! quantities for the current (interpolated) intermediate structure
! between the beginning and end of a model timestep, fully set every
! call by setup/midmod.f90 (called once per rotational-diffusion
! sub-timestep from rotation/getw.f90) and read broadly across
! mixing/, rotation/, and nuclear/liburn.f90/liburn2.f90. Explicitly
! documented in-file as "the rotating-model counterpart of
! common/scrtch/'s shell_diag%del_grad".
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! scrtch_lib/fluxes_lib: genuinely evolving per-model state, read from
! distant, unrelated points in the call graph. Field names are
! unchanged from the original COMMON member names, matching that same
! precedent.
module mdphy_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: mdphy_state
            double precision :: amum(json), cpm(json), delm(json)
            double precision :: del_adiabatic_mix(json), &
                 del_radiative_mix(json)
            double precision :: esumm(json), om(json), qdtm(json)
            double precision :: thdifm(json), velm(json), viscm(json)
            double precision :: epsm(json)
      end type mdphy_state

      type(mdphy_state), public :: mix_phys

end module mdphy_lib
