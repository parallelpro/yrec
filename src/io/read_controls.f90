!----------------------------------------------------------------------
! read_controls
!----------------------------------------------------------------------
! Added 2026 (phase five): the controls-read step of the former
! program main, as its own star-layer routine. parmin reads the
! control and physics namelists into const_lib's targets and fills
! the job-level file paths (star_job_lib). Errors return via ierr.
subroutine read_controls(ierr)

      use star_info_lib, only: star
      implicit none

      integer, intent(out) :: ierr

      ierr = 0
! read in user parameters
      call parmin(star%job%alex06_table_path,star%job%allard_table_path,star%job%atm_table_path,star%job%fermi_table_path,star%job%kurucz_table_path,star%job%kurucz_table2_path,star%job%laol_table_path, &
           star%job%laol_table2_path,star%job%opal95_table_path,star%job%opal92_table_path,star%job%zams_a_table_path,star%job%zams_b_table_path,star%job%zams_c_table_path,star%job%centre1_table_path,star%job%centre2_table_path,star%job%centre3_table_path,star%job%centre4_table_path, &
           star%job%centre5_table_path,star%job%opal92_table2_path,star%job%pulse_atm_path,star%job%pulse_env_path,star%job%pulse_mod_path,star%job%pure_z_table_path,star%job%scv_h_table_path,star%job%scv_he_table_path,star%job%scv_z_table_path,star%job%alex95_table_paths, ierr)
      if (ierr /= 0) return

      return
end subroutine read_controls
