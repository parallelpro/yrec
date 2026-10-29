!----------------------------------------------------------------------
! numerics_status
!----------------------------------------------------------------------
! The status codes every numerics routine returns in ierr: 0 on
! success, a positive code naming the failure otherwise. The routines
! print nothing; a caller that reports a failure writes
! '<routine>: '//numerics_status_text(ierr), a short description of the
! code, through its domain's messages.
module numerics_status
      implicit none
      private

      integer, parameter, public :: numerics_ok = 0
! a singular matrix: for ludcmp and simeqc singular to working
! precision (a zero row, or no pivot candidate above 1e-14 of its row's
! size); for solve_tridiagonal, which does not pivot, a zero pivot
      integer, parameter, public :: numerics_singular = 1
! two abscissae equal where they must differ (splint's bracketing
! points, polint's, lir's end points, osplin's end interval when a point
! lies beyond it)
      integer, parameter, public :: numerics_coincident_points = 2
! sizes or limits out of range (bsstep's, lir's stride), too few points
! (osplin, lir), a NaN pivot (ludcmp, simeqc)
      integer, parameter, public :: numerics_bad_input = 3
! bsstep: the step vanished against its underflow reference
      integer, parameter, public :: numerics_step_underflow = 4
! osplin: the evaluation points decrease somewhere
      integer, parameter, public :: numerics_unsorted_points = 5
! bsstep: the ode_system's derivatives failed (their own code comes back
! in bsstep's derivatives_ierr)
      integer, parameter, public :: numerics_derivatives_failed = 6

      public :: numerics_status_text

contains

!----------------------------------------------------------------------
! numerics_status_text
!----------------------------------------------------------------------
! A short description of a status code, for a caller's message.
pure function numerics_status_text(status) result(text)
      integer, intent(in) :: status
      character(len=:), allocatable :: text

      select case (status)
      case (numerics_ok)
         text = 'success'
      case (numerics_singular)
         text = 'singular matrix'
      case (numerics_coincident_points)
         text = 'coincident abscissae'
      case (numerics_bad_input)
         text = 'bad input (a size or limit out of range, too few points, or NaN)'
      case (numerics_step_underflow)
         text = 'step underflow'
      case (numerics_unsorted_points)
         text = 'evaluation points not sorted'
      case (numerics_derivatives_failed)
         text = 'derivatives failed'
      case default
         text = 'unknown numerics status'
      end select
end function numerics_status_text

end module numerics_status
