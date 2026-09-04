!----------------------------------------------------------------------
! intpar_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/intpar/: Bulirsch-Stoer integrator
! control parameters (tolerance_fraction/max_stage_index/extrap_order),
! set once from NAMELIST /physics/ (io/read_controls.f90, defaulted there
! via a DATA statement) and read by numerics/bsstep.f90 and
! atm/atm_lib.f90. Global configuration, not per-call data -- each
! caller reads a different subset (atm_lib.f90 only tolerance_fraction,
! bsstep.f90 only max_stage_index/extrap_order), confirming these
! aren't actually data flowing between the two, just shared constants.
module intpar_lib
      implicit none

      double precision :: tolerance_fraction
      integer :: max_stage_index, extrap_order

end module intpar_lib
