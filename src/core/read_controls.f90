!----------------------------------------------------------------------
! read_controls
!----------------------------------------------------------------------
! Added 2026 (phase five): the controls-read step of the former
! program main, as its own star-layer routine. parmin reads the
! control and physics namelists into const_lib's targets and fills
! the job-level file paths (star_job_lib). Errors return via ierr.
subroutine read_controls(ierr)

      use star_job_lib, only: job
      implicit none

      integer, intent(out) :: ierr

      ierr = 0
! read in user parameters
      call parmin(job%alex06_table_path,job%allard_table_path,job%atm_table_path,job%fermi_table_path,job%kurucz_table_path,job%kurucz_table2_path,job%laol_table_path, &
           job%laol_table2_path,job%opal95_table_path,job%opal92_table_path,job%zams_a_table_path,job%zams_b_table_path,job%zams_c_table_path,job%centre1_table_path,job%centre2_table_path,job%centre3_table_path,job%centre4_table_path, &
           job%centre5_table_path,job%opal92_table2_path,job%pulse_atm_path,job%pulse_env_path,job%pulse_mod_path,job%pure_z_table_path,job%scv_h_table_path,job%scv_he_table_path,job%scv_z_table_path,job%alex95_table_paths, ierr)
      if (ierr /= 0) return

      return
end subroutine read_controls
