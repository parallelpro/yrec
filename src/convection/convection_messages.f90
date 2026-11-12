!----------------------------------------------------------------------
! convection_messages
!----------------------------------------------------------------------
! Convection's messages, in one place (rotation_messages is the model).
! The kernels print nothing: each message below is formatted here from
! the values the kernel passes, in sentence case after the public entry
! point that met it ("find_convection_zones: ...",
! "extend_by_semiconvection: ...", "compute_temperature_gradients:
! ..."), a message's values on indented lines after its first:
!   to the run log:
!     log_too_many_zones        more convection zones than the zone list
!                               has rows (find_convection_zones)
!     log_mixed_regions         the mixed regions without and with
!                               overshoot, four to a line
!                               (extend_by_overshoot)
!     log_zones_merged          two mixed regions merged by overshoot or
!                               by semiconvection
!     log_semiconvective_edge   the core's edge moved by semiconvection
!     log_cubic_nonconvergence  the mixing-length cubic did not converge
!                               (the gradient is taken adiabatic)
!   to the terminal and the run log:
!     warn_unstable_radicand    the mixing-length square root of a
!                               negative number (the call returns ierr)
! The run log is run_messages' run_log_unit, written directly, as the
! kernels did (a test points it at standard output); the terminal is
! standard output. Each message is written into an array of lines with
! an internal write (one record per element) after every element is set
! to unwritten_line, and the lines the write reached are printed, so the
! output has the records the format makes, as a write to a unit would.
! Shell indices take five digits (json, 5000 shells).
module convection_messages
      use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
      use run_messages, only: run_log_unit
      implicit none
      private
      public :: log_too_many_zones, log_mixed_regions, log_zones_merged, &
           log_semiconvective_edge, log_cubic_nonconvergence, warn_unstable_radicand

! What merged two mixed regions (log_zones_merged's heading).
      integer, parameter, public :: merged_by_overshoot = 1, merged_by_semiconvection = 2

      integer, parameter :: message_length = 200
      character(len=*), parameter :: unwritten_line = achar(0)

contains

! The zone list has num_rows rows, all of them used.
subroutine log_too_many_zones(num_rows)
      integer, intent(in) :: num_rows
      character(len=message_length) :: lines(1)
      lines = unwritten_line
      write(lines, '("find_convection_zones: more convection zones than the zone list''s ", &
           &i0," rows")') num_rows
      call write_convection_lines(run_log_unit, lines)
end subroutine log_too_many_zones

! The first num_zones regions of each list (first and last shell), four
! to a line; the format stops after the last region.
subroutine log_mixed_regions(bounds_no_overshoot, bounds, num_zones)
      integer, intent(in) :: bounds_no_overshoot(:,:), bounds(:,:), num_zones
      character(len=message_length) :: lines(max(1, (num_zones + 3)/4))
      integer :: k, j
      lines = unwritten_line
      write(lines, '("find_convection_zones: mixed regions without overshoot", &
           &4(:," [",i5,"-",i5,"]"))') ((bounds_no_overshoot(k,j), j = 1, 2), k = 1, num_zones)
      call write_convection_lines(run_log_unit, lines)
      lines = unwritten_line
      write(lines, '("find_convection_zones: mixed regions with overshoot   ", &
           &4(:," [",i5,"-",i5,"]"))') ((bounds(k,j), j = 1, 2), k = 1, num_zones)
      call write_convection_lines(run_log_unit, lines)
end subroutine log_mixed_regions

! The two regions before the merge (first and last shell of each, from
! the centre out) and the merged one.
subroutine log_zones_merged(cause, lower, upper, merged)
      integer, intent(in) :: cause, lower(:), upper(:), merged(:)
      character(len=message_length) :: lines(1)
      lines = unwritten_line
      if (cause == merged_by_overshoot) then
         write(lines, '("find_convection_zones: mixed regions",2(" [",i5,"-",i5,"]"), &
              &" merged by overshoot into [",i5,"-",i5,"]")') lower, upper, merged
      else
         write(lines, '("extend_by_semiconvection: mixed regions",2(" [",i5,"-",i5,"]"), &
              &" merged by semiconvection into [",i5,"-",i5,"]")') lower, upper, merged
      end if
      call write_convection_lines(run_log_unit, lines)
end subroutine log_zones_merged

! The core's upper edge moved from shell old_edge to new_edge;
! reached_limit when the distance the instability propagates in the
! step, or the surface, stopped it. The radiative gradients at the first
! shell above the old edge, with the core's composition and with the
! shell's own.
subroutine log_semiconvective_edge(old_edge, new_edge, reached_limit, &
     radiative_gradient_zone_mix, radiative_gradient_local)
      integer, intent(in) :: old_edge, new_edge
      logical, intent(in) :: reached_limit
      real(dp), intent(in) :: radiative_gradient_zone_mix, radiative_gradient_local
      character(len=message_length) :: lines(1)
      character(len=16) :: limit_note
      limit_note = ""
      if (reached_limit) limit_note = ", at its limit"
      lines = unwritten_line
      write(lines, '("extend_by_semiconvection: the core''s edge moved from shell",i5, &
           &" to",i5,a,"; del_rad above the old edge with the core''s composition",f8.4, &
           &", with its own",f8.4)') old_edge, new_edge, trim(limit_note), &
           radiative_gradient_zone_mix, radiative_gradient_local
      call write_convection_lines(run_log_unit, lines)
end subroutine log_semiconvective_edge

subroutine log_cubic_nonconvergence(log_pressure, log_temperature, opacity, &
     specific_heat_cp, dlnrho_dlnt)
      real(dp), intent(in) :: log_pressure, log_temperature, opacity, &
           specific_heat_cp, dlnrho_dlnt
      character(len=message_length) :: lines(2)
      lines = unwritten_line
      write(lines, '("compute_temperature_gradients: the mixing-length cubic did not ", &
           &"converge; del = del_ad"/"  log P",f10.6,", log T",f10.6,", kappa",1pe13.5, &
           &", c_P",e13.5,", dln rho/dln T",e13.5)') log_pressure, log_temperature, &
           opacity, specific_heat_cp, dlnrho_dlnt
      call write_convection_lines(run_log_unit, lines)
end subroutine log_cubic_nonconvergence

! The superadiabatic excess, dln rho/dln T and the radicand, then the
! point's opacity, L/Lsun, log P, log m, log T, the two rotation
! factors, the two gradients and log r.
subroutine warn_unstable_radicand(superadiabatic_excess, dlnrho_dlnt, radicand, &
     opacity, luminosity_lsun, log_pressure, log_mass, log_temperature, &
     temperature_rotation_factor, pressure_rotation_factor, radiative_gradient, &
     adiabatic_gradient, log_radius)
      real(dp), intent(in) :: superadiabatic_excess, dlnrho_dlnt, radicand, &
           opacity, luminosity_lsun, log_pressure, log_mass, log_temperature, &
           temperature_rotation_factor, pressure_rotation_factor, radiative_gradient, &
           adiabatic_gradient, log_radius
      character(len=message_length) :: lines(4)
      lines = unwritten_line
      write(lines, '("compute_temperature_gradients: the mixing-length square root of ", &
           &"a negative number; ierr returned"/ &
           &"  del_rad - del_ad",1pe12.3,", dln rho/dln T",e12.3,", radicand",e12.3/ &
           &"  kappa",e12.3,", L/Lsun",e12.3,", log P",e12.3,", log m",e12.3,", log T",e12.3/ &
           &"  f_T",e12.3,", f_P",e12.3,", del_rad",e12.3,", del_ad",e12.3,", log r",e12.3)') &
           superadiabatic_excess, dlnrho_dlnt, radicand, opacity, luminosity_lsun, &
           log_pressure, log_mass, log_temperature, temperature_rotation_factor, &
           pressure_rotation_factor, radiative_gradient, adiabatic_gradient, log_radius
      call write_convection_lines(output_unit, lines)
      call write_convection_lines(run_log_unit, lines)
end subroutine warn_unstable_radicand

subroutine write_convection_lines(unit, lines)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: lines(:)
      integer :: k
      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         write(unit, '(a)') trim(lines(k))
      end do
end subroutine write_convection_lines

end module convection_messages
