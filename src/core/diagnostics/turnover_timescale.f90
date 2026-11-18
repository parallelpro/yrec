!----------------------------------------------------------------------
! turnover_timescale
!----------------------------------------------------------------------
! The convective turnover timescale near the base of the surface
! convection zone, from the stitched model: compute_turnover_timescale
! gathers the points and stores the result; the walk itself is the
! star-blind kernel time_surface_convection_zone (turnover_walk.f90).
module turnover_timescale
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: compute_turnover_timescale

! The walk takes the envelope's points with the interior's when the
! surface convection zone's base lies less than this in log P below the
! interior's outer point [dex].
      real(dp), parameter :: shallow_zone_dex = 1.0d0

contains

!----------------------------------------------------------------------
! compute_turnover_timescale (formerly gettau; G. Somers 3/17)
!----------------------------------------------------------------------
! The convective turnover timescale star%convective_turnover_timescale
! (s), one pressure scale height above the base of the surface
! convection zone, and star%envelope_radius_rsun, the radius of that
! base (Rsun; 0 without a surface convection zone). It reads the stitched model
! (core/diagnostics/stitched_model.f90), the converged interior followed
! by the fixed-step envelope integration, built once per step in
! evolve_step before compute_observables, and walks it as one structure
! with time_surface_convection_zone. It integrates nothing itself, and
! star%photosphere_logP belongs to the stitched build (its atm_get).
! During the next step the wind reads the value between the copy the
! step stored at its end (the _prev member) and the current member, and
! the deuterium limiter reads the current member.
! Called by read_starting_model and by evolve_step's
! refresh_turnover_timescale (core/step/step_driver.f90). A spline of
! the walk that fails is reported here, and its status returned in ierr.
!   star reads: convective_flag, logP, nz, solar_mass_cgs,
!     solar_radius_cgs, stitched
!   star writes: convective_turnover_timescale, envelope_radius_rsun
subroutine compute_turnover_timescale(s, ierr)
      use stitched_model, only: &
           ip_mass, ip_logR, ip_logRho, ip_logP, ip_mixing_type, ip_gradr, &
           ip_grada, ip_conv_vel
      use star_info_lib, only: json
      use physical_constants, only: cgl, ln10
      use turnover_walk, only: turnover_result, time_surface_convection_zone, &
           velocity_spline_failed
      use numerics_lib, only: numerics_status_text
      use core_messages, only: core_report, message_length, unwritten_line
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      integer, intent(out) :: ierr

      real(dp) :: combined_radius(json), combined_pressure(json), &
           combined_density(json), combined_mass(json), &
           combined_gravity(json), combined_velocity(json)
      real(dp) :: combined_grad1(json), combined_grad2(json)
      logical :: combined_convective_flag(json)
      integer :: j, num_points
      type(turnover_result) :: turnover
      character(len=message_length) :: lines(1)

      ierr = 0
! Start from zero; photosphere_logP is left alone (the stitched build owns it).
      s%convective_turnover_timescale = 0.0

! The span. In the interior model, zones j+1..nz-1 are the convective
! ones below the surface (j is the outermost radiative zone below nz).
! When the base of that zone lies less than shallow_zone_dex in pressure
! below the interior's surface, log P(j+1) - log P(nz) < 1 (or the outer
! zones are radiative), the walk takes the stitched interior and
! envelope (regions 1-2 of the stitched grid, never the atmosphere),
! clamped to the json-sized arrays the walker takes. A deeper zone, or a
! fully convective star, is walked on the nz interior points alone.
      do j = s%nz-1, 1, -1
         if (.not. s%convective_flag(j)) exit
      end do
      if (s%logP(j+1) - s%logP(s%nz) < shallow_zone_dex) then
         num_points = min(s%stitched%n_ie, json)
      else
         num_points = s%nz
      end if
! The walker's inputs: m (g), log r, log P, log rho, g = G m/r^2, the
! convective velocity, gradr, grada and the convective flag.
      do j = 1, num_points
         combined_mass(j) = s%stitched%prof(ip_mass,j)*s%solar_mass_cgs
         combined_radius(j) = s%stitched%prof(ip_logR,j)
         combined_pressure(j) = s%stitched%prof(ip_logP,j)
         combined_density(j) = s%stitched%prof(ip_logRho,j)
         combined_gravity(j) = combined_mass(j)* &
              exp(ln10*(cgl - 2.0d0*combined_radius(j)))
         combined_velocity(j) = s%stitched%prof(ip_conv_vel,j)
         combined_grad1(j) = s%stitched%prof(ip_gradr,j)
         combined_grad2(j) = s%stitched%prof(ip_grada,j)
         combined_convective_flag(j) = s%stitched%prof(ip_mixing_type,j) /= 0.0d0
      end do

      call time_surface_convection_zone(combined_mass, combined_convective_flag, &
           combined_radius, combined_pressure, combined_density, combined_gravity, &
           combined_velocity, combined_grad1, combined_grad2, num_points, s%nz, &
           turnover, ierr)
      if (ierr /= 0) then
         lines = unwritten_line
         if (turnover%failed_spline == velocity_spline_failed) then
            write(lines,'(a)') ' ERROR: turnover time: the spline of the convective' &
                 //' velocity at the convection-zone base failed: ' &
                 //trim(numerics_status_text(ierr))
         else
            write(lines,'(a)') ' ERROR: turnover time: the spline of H_P at the' &
                 //' convection-zone base failed: '//trim(numerics_status_text(ierr))
         end if
         call core_report(lines)
         return
      end if
      s%convective_turnover_timescale = turnover%turnover_time
      s%envelope_radius_rsun = turnover%base_radius/s%solar_radius_cgs
end subroutine compute_turnover_timescale

end module turnover_timescale
