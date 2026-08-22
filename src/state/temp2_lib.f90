!----------------------------------------------------------------------
! temp2_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/temp2/: per-shell rotational
! circulation/shear velocities (Eddington-Sweet, GSF-instability,
! secular-shear, mu-gradient) and their previous-iteration
! counterparts, computed by rotation/seculr/vcirc.f90 (see that file's own
! header for the physics) and read broadly across rotation/ and the
! output writers.
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! mdphy_lib/temp_lib: genuinely evolving per-model state, read from
! distant, unrelated points in the call graph. Field names are
! unchanged from the original COMMON member names, matching that same
! precedent.
module temp2_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: circulation_velocity_state
            double precision :: es_circulation_velocity(json), &
                 es_circulation_velocity_prev(json)
            double precision :: secular_shear_velocity(json), &
                 secular_shear_velocity_prev(json)
            double precision :: hle(json)
            double precision :: gsf_circulation_velocity(json), &
                 gsf_circulation_velocity_prev(json)
            double precision :: mu_gradient_velocity(json)
      end type circulation_velocity_state
! 2026 (phase six, step 1 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module temp2_lib
