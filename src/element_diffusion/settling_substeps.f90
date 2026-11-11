!----------------------------------------------------------------------
! settling_substeps
!----------------------------------------------------------------------
! How a step's settling is cut into sub-steps (MHP 6/90): none longer
! than a fraction (the control settling_timestep_fraction) of helium's
! settling timescale at the region's outer boundary, after Michaud et
! al. (1984, ApJ 282, 206),
!     tau = 4.348e21 yr q / T_b^1.5,  q = (M - m_b)/M,
! the mass fraction above the boundary and T_b its temperature (K). The
! sub-step count is max(1, ceiling(timestep/(fraction tau))): the
! quotient's integer part, plus one when the timestep is not an exact
! multiple of the longest sub-step, and at least one. (The F77 took
! mod(longest, timestep), arguments reversed, which added a sub-step at
! an exact multiple; docs/bug_sweep.md O16.)
! Both callers of diffuse_elements use it: core's composition step and
! rotation's sub-step, each with its own outer boundary.
module settling_substeps
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: ln10, seconds_per_year
      use math_lib
      implicit none
      private
      public :: settling_substep_plan, plan_settling_substeps

! Michaud et al.'s coefficient, in s: 4.348e21 yr
      real(dp), parameter :: michaud_timescale_coefficient = 4.348d21*seconds_per_year

! The sub-steps of one step and the quantities they come from (the
! rotation sub-step's solver diagnostics print them).
      type :: settling_substep_plan
         real(dp) :: mass_fraction_above = 0.0d0  ! q above the boundary
         real(dp) :: timescale = 0.0d0             ! tau, s
         real(dp) :: longest_substep = 0.0d0       ! fraction tau, s
         integer :: num_substeps = 1
         real(dp) :: substep = 0.0d0               ! timestep/num_substeps, s
      end type settling_substep_plan

contains

! The sub-steps of a settling step of timestep seconds: the boundary
! zone has enclosed mass boundary_mass (g) in a star of total_mass (g)
! and log10 temperature boundary_log_temperature; no sub-step exceeds
! timestep_fraction times the settling timescale there.
pure function plan_settling_substeps(timestep, timestep_fraction, total_mass, &
     boundary_mass, boundary_log_temperature) result(plan)
      real(dp), intent(in) :: timestep, timestep_fraction, total_mass, &
           boundary_mass, boundary_log_temperature
      type(settling_substep_plan) :: plan

      plan%mass_fraction_above = (total_mass - boundary_mass)/total_mass
      plan%timescale = michaud_timescale_coefficient*plan%mass_fraction_above/ &
           exp(ln10*1.5d0*boundary_log_temperature)
      plan%longest_substep = timestep_fraction*plan%timescale
      plan%num_substeps = int(timestep/plan%longest_substep)
      if (mod(timestep, plan%longest_substep) /= 0.0d0 .or. plan%num_substeps == 0) &
           plan%num_substeps = plan%num_substeps + 1
      plan%substep = timestep/real(plan%num_substeps, dp)
end function plan_settling_substeps

end module settling_substeps
