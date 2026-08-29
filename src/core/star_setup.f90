!----------------------------------------------------------------------
! star_setup
!----------------------------------------------------------------------
! Added 2026 (phase five): the setup step of the former program main,
! as its own star-layer routine: constants and table loads via setups
! (which threads ierr into the kap/eos/atm facade inits), then the
! Monte-Carlo run-range and its per-run parameter read. Errors return
! via ierr.
subroutine star_setup(ierr)

      use star_info_lib, only: star
      use monte_carlo_lib, only: setup_monte_carlo_runs
      implicit none

      integer, intent(out) :: ierr

      ierr = 0
! set up constants and read in tabular data
! MHP 8/25 directly pass file names instead of using common blocks
      call setups(star%job%mixture_weights,star%job%alex06_table_path,star%job%allard_table_path,star%job%atm_table_path,star%job%fermi_table_path,star%job%kurucz_table_path,star%job%kurucz_table2_path, &
           star%job%laol_table_path,star%job%laol_table2_path,star%job%opal95_table_path,star%job%opal92_table_path,star%job%zams_a_table_path,star%job%zams_b_table_path,star%job%zams_c_table_path,star%job%centre1_table_path,star%job%centre2_table_path,star%job%centre3_table_path, &
           star%job%centre4_table_path,star%job%centre5_table_path,star%job%opal92_table2_path,star%job%pure_z_table_path,star%job%scv_h_table_path,star%job%scv_he_table_path,star%job%scv_z_table_path,star%job%alex95_table_paths, ierr)
      if (ierr /= 0) return
! Monte-Carlo run range and per-run sampled-parameter read (2026:
! core/monte_carlo.f90 -- the standalone MC home).
      call setup_monte_carlo_runs

      return
end subroutine star_setup
