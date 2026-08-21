!----------------------------------------------------------------------
! light_burn_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/newrat/, common/oldrat/,
! common/liov/, common/prevcz/, and common/deuter/: per-model state for
! lithium-6/lithium-7/beryllium-9 and deuterium burning, read/written
! across nuclear/ (liburn.f90/liburn2.f90/lirate88.f90/dburn.f90/
! dburnm.f90/deutrate.f90), mixing/ (mix.f90/mixcz.f90/bursmix.f90),
! and several more distant files that only touch one or two of these
! members (e.g. setup/hpoint.f90's deuterium_burning_rate_start,
! wind/mdot.f90's accreted_mass_fraction).
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! oldmod_lib/scrtch_lib/turnover_lib: genuinely evolving per-model
! state, read from many distant, unrelated points in the call graph.
! Bundled into one type despite the five different original COMMON
! block names because they're all facets of the same light-element-
! burning subsystem and are already handed around together by every
! file that uses more than one of them. Field names are unchanged from
! the original COMMON member names, matching that same precedent.
module light_burn_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: light_element_burn_state
! former common/newrat/: Li6/Li7/Be9 burning rates at the end of the
! timestep, at the current (possibly overshoot-adjusted) depth.
            double precision :: rate_li6(json), rate_li7(json), &
                 rate_be9(json)
! former common/oldrat/: same, at the start of the timestep.
            double precision :: rate_li6_start(json), rate_li7_start(json), &
                 rate_be9_start(json)
! former common/liov/: pressure scale heights used to search downward
! from the CZ base for the true (overshoot-corrected) base location.
            double precision :: pressure_scale_height_start, &
                 pressure_scale_height_end
! former common/prevcz/: previous end-of-timestep values, used as the
! new beginning-of-timestep values.
            double precision :: cz_base_radius_prev, log_rate_li6_prev, &
                 log_rate_li7_prev, log_rate_be9_prev
            integer :: envelope_cz_base_zone_prev
! former common/deuter/: deuterium burning rate (current/start of
! timestep), accreted mass fraction, and the convection-zone base zone
! index used by the deuterium-burning routines.
            double precision :: deuterium_burning_rate(json), &
                 deuterium_burning_rate_start(json)
            double precision :: accreted_mass_fraction
            integer :: jcz
      end type light_element_burn_state

      type(light_element_burn_state), public :: light_burn

end module light_burn_lib
