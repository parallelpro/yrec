!----------------------------------------------------------------------
! read_controls
!----------------------------------------------------------------------
! Added 2026 (phase five): the controls-read step of the former
! program main, as its own star-layer routine. parmin reads the
! control and physics namelists into the controls_lib BUFFER and
! fills the job-level file paths (star%job). Errors return via ierr.
!
! 2026 controls->star% campaign, phase B: star%ctrl is the
! authoritative post-read home. The four-step sequence below
! guarantees parmin always reads over pristine defaults (a namelist
! read only overwrites what the input file provides) -- it replaces
! the retired controls_reset_lib snapshot machinery structurally:
!   1. star%ctrl = controls_state()   pristine defaults (component
!                                     initializers, generated from
!                                     the buffer's declarations)
!   2. seed_controls_buffer           buffer <- star%ctrl
!   3. parmin                         namelist reads over the buffer
!   4. store_controls_to_star         star%ctrl <- buffer
! Until phase C migrates them, consumers still read the buffer's
! bare names; star%ctrl is written here but not yet read broadly.
subroutine read_controls(ierr)

      use star_info_lib, only: star, controls_state
      use controls_sync_lib, only: seed_controls_buffer, &
           store_controls_to_star
      implicit none

      integer, intent(out) :: ierr

      ierr = 0
      star%ctrl = controls_state()
      call seed_controls_buffer
! read in user parameters
      call parmin(star%job%alex06_table_path,star%job%allard_table_path,star%job%atm_table_path,star%job%fermi_table_path,star%job%kurucz_table_path,star%job%kurucz_table2_path,star%job%laol_table_path, &
           star%job%laol_table2_path,star%job%opal95_table_path,star%job%opal92_table_path,star%job%zams_a_table_path,star%job%zams_b_table_path,star%job%zams_c_table_path,star%job%centre1_table_path,star%job%centre2_table_path,star%job%centre3_table_path,star%job%centre4_table_path, &
           star%job%centre5_table_path,star%job%opal92_table2_path,star%job%pulse_atm_path,star%job%pulse_env_path,star%job%pulse_mod_path,star%job%pure_z_table_path,star%job%scv_h_table_path,star%job%scv_he_table_path,star%job%scv_z_table_path,star%job%alex95_table_paths, ierr)
      if (ierr /= 0) return
      call store_controls_to_star

      return
end subroutine read_controls
