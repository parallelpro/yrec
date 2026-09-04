!----------------------------------------------------------------------
! point_scratch_lib
!----------------------------------------------------------------------
! The shell-physics POINT SCRATCH (2026, envint purity split phase A):
! the per-point side channel through which the integrand callbacks
! (atmosphere_derivs / envelope_derivs, via their eos/kap/gradient
! calls) and henyey_coefficients publish the physics of the point
! they just evaluated, and through which envint_kernel's recorder reads
! it into atmo_struct / env_struct.
!
! Formerly star%pulse (the pulsation_diagnostics_state q* scratch,
! ex common/pulse1/pulse2) and the flat star%current_* members --
! solver working state that never belonged in star_info (same
! taxonomy call as rotation_scratch_lib: rot_scr/mix_scr/circ_scr).
! yrec_reset snapshots pt_scr alongside those for in-process
! re-entrancy.
module point_scratch_lib
      implicit none

      type, public :: point_physics_scratch
! the q* scratch (ex common/pulse1/pulse2; qqed/qqet/qfs deleted
! earlier as write-only)
            double precision :: qqdp, qqod, qqot, qdel, qdela, &
                 qqcp
            double precision :: qrmu, qtl, qpl, qdl, qo, qol, qt, qp
            double precision :: qqdt, qemu, qd
! the "current point" slots (ex flat star%current_*)
            double precision :: current_log10_pressure, current_log10_temperature, &
                 current_log10_radius, current_log10_mass, current_log10_density, &
                 current_opacity, current_beta, current_gradients(3), &
                 current_ion_fraction(3), current_velocity
      end type point_physics_scratch

      type(point_physics_scratch), save :: pt_scr

end module point_scratch_lib
