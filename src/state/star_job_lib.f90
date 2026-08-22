!----------------------------------------------------------------------
! star_job_lib
!----------------------------------------------------------------------
! Added 2026 (phase five -- the embeddable engine; ROADMAP.md). The
! job-level configuration that defines a YREC run, as distinct from
! the model (star_info) and from the physics controls (const_lib's
! namelist targets, which stay where the namelists put them):
!
!  * the table/file paths parmin reads from the control namelist and
!    setups consumes (formerly ~26 locals of run_yrec threaded
!    through both calls),
!  * the LAOL mixture-weight work array setups fills and starin
!    consumes,
!  * the Monte-Carlo run-range bounds.
!
! One module-level instance (`job`), same single-instance pattern as
! star_info: one process, one job. read_controls fills it,
! star_setup consumes it, the run loop reads it.
module star_job_lib
      implicit none
      private

      type, public :: star_job
            character(len=256) :: alex06_table_path, allard_table_path, &
                 atm_table_path, fermi_table_path, kurucz_table_path, &
                 kurucz_table2_path, laol_table_path, laol_table2_path, &
                 opal95_table_path, opal92_table_path
            character(len=256) :: zams_a_table_path, zams_b_table_path, &
                 zams_c_table_path, centre1_table_path, centre2_table_path, &
                 centre3_table_path, centre4_table_path, centre5_table_path
            character(len=256) :: opal92_table2_path, pure_z_table_path, &
                 scv_h_table_path, scv_he_table_path, scv_z_table_path
            character(len=256) :: alex95_table_paths(7)
            character(len=256) :: pulse_atm_path, pulse_env_path, &
                 pulse_mod_path
            double precision :: mixture_weights(12)
            integer :: mc_run_start, mc_run_end
      end type star_job

      type(star_job), public, save :: job

      ! 2026 (libyrec): namelist-path overrides for embedded use. The C
      ! API (core/yrec_capi.f90) sets these before calling run_yrec;
      ! blank means the historical CLI behavior (parmin's getarg).
      ! Deliberately module-level rather than star_job members: they
      ! configure the NEXT run from outside the engine, so they must
      ! not be captured/restored by yrec_reset's job snapshot.
      character(len=256), public, save :: control_nml_override = ' '
      character(len=256), public, save :: physics_nml_override = ' '

end module star_job_lib
