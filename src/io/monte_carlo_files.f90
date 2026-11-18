!----------------------------------------------------------------------
! monte_carlo_files
!----------------------------------------------------------------------
! The files of the Monte-Carlo runs, read and written in their fixed
! layouts. run_files opens them with the run's other files: the sample
! file on dynamics_unit, the summary on mc_summary_unit and the models
! on mc_model_unit (written by model_files' write_mod_model).
!   read_monte_carlo_sample    one run's samples, a two-line record of
!                              the sample file
!   write_monte_carlo_summary  one run's summary, seven records, flushed
!                              and checked on disk (output_files)
!   rewind_run_cycle_files     the final model, the starting model and
!                              the run log rewound, so that the next
!                              calibration cycle or Monte-Carlo run
!                              rewrites them from the start
! Nothing here reads the star: core/run/monte_carlo.f90 gathers the
! values, calls these and reports a failed read; a failed write is
! reported by flush_output_file ("ERROR: cannot write <path>: <reason>",
! ierr = 1).
module monte_carlo_files
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use run_messages, only: run_log_unit
      use run_files, only: last_model_unit, starting_model_unit, &
           dynamics_unit, mc_summary_unit
      use output_files, only: flush_output_file
      implicit none
      private
      public :: read_monte_carlo_sample, write_monte_carlo_summary, &
           rewind_run_cycle_files

! The values of one run's record of the sample file, in the file's
! order: the S11, S33, S34 and S17 factors, Z/X, the metal settling
! factor, and the luminosity and age factors.
      integer, parameter, public :: num_mc_sample_values = 8

! One run's summary, as write_monte_carlo_summary writes it.
      type, public :: monte_carlo_summary
! record 1: the run number and its samples (num_mc_sample_values, in
! the sample file's order)
         integer :: run_number = 0
         real(dp) :: samples(num_mc_sample_values) = 0.0_dp
! record 2: the number of solar-calibration cycles and the
! calibration's partial derivatives d log L/dX, d log R/dX,
! d log L/d alpha, d log R/d alpha
         integer :: calibration_cycles = 0
         real(dp) :: calibration_derivatives(4) = 0.0_dp
! record 3: the Cl37 and Ga71 capture rates (SNU) and the neutrino
! fluxes pp, pep, hep, Be7, B8, N13, O15 and F17 (1e10 /cm^2/s)
         real(dp) :: cl37_snu = 0.0_dp, ga71_snu = 0.0_dp
         real(dp) :: neutrino_fluxes(8) = 0.0_dp
! record 4: the central T (1e6 K), density (g/cm^3), P (1e17 dyn/cm^2),
! X and Z
         real(dp) :: central_t6 = 0.0_dp, central_density = 0.0_dp, &
              central_p17 = 0.0_dp, central_x = 0.0_dp, central_z = 0.0_dp
! record 5: the number of shells, alpha, the initial Y and Z, and the
! final log(L/Lsun) and log(R/Rsun)
         integer :: num_shells = 0
         real(dp) :: alpha = 0.0_dp, initial_y = 0.0_dp, initial_z = 0.0_dp, &
              log_l_lsun = 0.0_dp, log_r_rsun = 0.0_dp
! record 6: the radius (Rsun) of the base of the surface convection
! zone and the mass above it (Msun), and the surface Y, Z and Z/X
         real(dp) :: cz_base_radius_rsun = 0.0_dp, cz_mass_msun = 0.0_dp, &
              surface_y = 0.0_dp, surface_z = 0.0_dp, surface_z_over_x = 0.0_dp
! record 7: the luminosities (Lsun) of the pp I, pp II, pp III and CNO
! channels and of gravitational energy release
         real(dp) :: luminosities_lsun(5) = 0.0_dp
      end type monte_carlo_summary

contains

! ---------------------------------------------------------------
! One run's samples: the next two-line record of the sample file
! (dynamics_unit), into samples in the file's order. ios is the read's
! status, nonzero with the runtime's message in msg when the file ends
! before the record or the record does not read; samples is then not to
! be used.
subroutine read_monte_carlo_sample(samples, ios, msg)
      real(dp), intent(out) :: samples(:)
      integer, intent(out) :: ios
      character(len=*), intent(inout) :: msg

      read(dynamics_unit,1511,iostat=ios,iomsg=msg) samples(1:num_mc_sample_values)
 1511 format(7X,1P7E10.3/E9.3)
end subroutine read_monte_carlo_sample

! ---------------------------------------------------------------
! One run's summary appended to the summary file (mc_summary_unit; path
! is its name, for the error message): seven records, each written
! unless an earlier one failed, then flushed and checked on disk
! (flush_output_file: the first failed write, or a file cut short). The
! metal settling factor is written twice in record 1, as WRTMONTE wrote
! FHE and its copy FZDIF. ierr = 1 when a write failed.
subroutine write_monte_carlo_summary(summary, path, ierr)
      type(monte_carlo_summary), intent(in) :: summary
      character(len=*), intent(in) :: path
      integer, intent(inout) :: ierr
      integer :: ios
      character(len=256) :: msg

      ios = 0
      if (ios == 0) write(mc_summary_unit,10,iostat=ios,iomsg=msg) summary%run_number, &
           summary%samples(1:6), summary%samples(6), summary%samples(7:8)
   10 format(I7,1P9E10.3)
      if (ios == 0) write(mc_summary_unit,20,iostat=ios,iomsg=msg) summary%calibration_cycles, &
           summary%calibration_derivatives
   20 format(1X,I2,1P4E11.4)
      if (ios == 0) write(mc_summary_unit,30,iostat=ios,iomsg=msg) summary%cl37_snu, &
           summary%ga71_snu, summary%neutrino_fluxes
   30 format(1X,2F8.3,1P8E10.3)
      if (ios == 0) write(mc_summary_unit,40,iostat=ios,iomsg=msg) summary%central_t6, &
           summary%central_density, summary%central_p17, summary%central_x, summary%central_z
   40 format(1X,F7.3,F7.2,F6.3,2F8.5)
      if (ios == 0) write(mc_summary_unit,50,iostat=ios,iomsg=msg) summary%num_shells, &
           summary%alpha, summary%initial_y, summary%initial_z, summary%log_l_lsun, &
           summary%log_r_rsun
   50 format(I5,F7.4,2F8.5,1P2E10.3)
      if (ios == 0) write(mc_summary_unit,60,iostat=ios,iomsg=msg) summary%cz_base_radius_rsun, &
           summary%cz_mass_msun, summary%surface_y, summary%surface_z, summary%surface_z_over_x
   60 format(F8.5,F9.6,2F8.5,F9.6)
      if (ios == 0) write(mc_summary_unit,70,iostat=ios,iomsg=msg) summary%luminosities_lsun
   70 format(1P5E10.3)
      call flush_output_file(mc_summary_unit, path, ios, msg, ierr)
end subroutine write_monte_carlo_summary

! ---------------------------------------------------------------
! The final model, the starting model and the run log rewound, so that
! the next calibration cycle or Monte-Carlo run rewrites the first and
! the last from the start and reads the second again.
subroutine rewind_run_cycle_files()
      rewind(last_model_unit)
      rewind(starting_model_unit)
      rewind(run_log_unit)
end subroutine rewind_run_cycle_files

end module monte_carlo_files
