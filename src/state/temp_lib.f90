!----------------------------------------------------------------------
! temp_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/temp/: per-shell specific heat, mean
! molecular weight, and diffusion coefficients, computed by
! misc/physic.f90 (cp/mean_molecular_weight/qdt) and rotation/viscos.f90
! (visc/thdif) every model, and read broadly (setup/hpoint.f90 copies
! them into oldmod_lib's previous-model arrays; rotation/getw.f90
! reads mean_molecular_weight during rotational mixing).
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! scrtch_lib/mdphy_lib: genuinely evolving per-model state, read from
! distant, unrelated points in the call graph. Field names are
! unchanged from the original COMMON member names, matching that same
! precedent.
module temp_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: shell_temp_state
            double precision :: cp(json), mean_molecular_weight(json)
            double precision :: qdt(json), thdif(json), visc(json)
      end type shell_temp_state

      type(shell_temp_state), public :: shell_temp

end module temp_lib
