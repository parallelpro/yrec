!----------------------------------------------------------------------
! star_setup
!----------------------------------------------------------------------
! Added 2026 (phase five): the setup step of the former program main,
! as its own star-layer routine: constants and table loads via setups
! (which threads ierr into the kap/eos/atm facade inits), then the
! Monte-Carlo run-range and its per-run parameter read. Errors return
! via ierr.
subroutine star_setup(ierr)

      use monte_carlo_lib, only: setup_monte_carlo_runs
      implicit none

      integer, intent(out) :: ierr

      ierr = 0
! set up constants and read in tabular data
      call setups(ierr)
      if (ierr /= 0) return
! Monte-Carlo run range and per-run sampled-parameter read (2026:
! core/monte_carlo.f90 -- the standalone MC home).
      call setup_monte_carlo_runs

      return
end subroutine star_setup
