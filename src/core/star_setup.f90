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
      use luout_lib
      use const_lib
      implicit none

      integer, intent(out) :: ierr

      integer :: i

      ierr = 0
! set up constants and read in tabular data
! MHP 8/25 directly pass file names instead of using common blocks
      call setups(star%job%mixture_weights,star%job%alex06_table_path,star%job%allard_table_path,star%job%atm_table_path,star%job%fermi_table_path,star%job%kurucz_table_path,star%job%kurucz_table2_path, &
           star%job%laol_table_path,star%job%laol_table2_path,star%job%opal95_table_path,star%job%opal92_table_path,star%job%zams_a_table_path,star%job%zams_b_table_path,star%job%zams_c_table_path,star%job%centre1_table_path,star%job%centre2_table_path,star%job%centre3_table_path, &
           star%job%centre4_table_path,star%job%centre5_table_path,star%job%opal92_table2_path,star%job%pure_z_table_path,star%job%scv_h_table_path,star%job%scv_he_table_path,star%job%scv_z_table_path,star%job%alex95_table_paths, ierr)
      if (ierr /= 0) return
! MHP 3/96 changed I/O to read in only up to max run needed.
      if (lmonte) then
!c MHP 8/25 moved file open to parmin
!     OPEN(UNIT=IDYN,FILE=FDYN,FORM='FORMATTED',STATUS='OLD')
         star%job%mc_run_start = imbeg
         imend = min(imend,1000)
         star%job%mc_run_end = imend
! read in monte carlo data
         do i = 1,imend
            read(dynamics_unit,1511)star%run%s11_rate(i),star%run%s33_rate(i),star%run%s34_rate(i), &
                 star%run%s17_rate(i),star%run%metal_to_h_ratio(i),star%run%helium_fraction_param(i), &
                 star%run%luminosity_target(i),star%run%age_target(i)
 1511       format(7X,1P7E10.3/E9.3)
            write(iowr,*)i,star%run%s11_rate(i),star%run%s33_rate(i),star%run%s34_rate(i),star%run%s17_rate(i), &
                 star%run%metal_to_h_ratio(i),star%run%helium_fraction_param(i), &
                 star%run%luminosity_target(i),star%run%age_target(i)
            star%run%diffusion_factor(i) = star%run%helium_fraction_param(i)
         end do
      else
         star%job%mc_run_start = 1
         star%job%mc_run_end = 1
      endif

      return
end subroutine star_setup
