!----------------------------------------------------------------------
! convection_messages
!----------------------------------------------------------------------
! Convection's messages, in one place (rotation_messages is the model).
! The kernels print nothing: each message below is formatted here from
! the values the kernel passes, with the text and the format it has
! always had, and written to the destination it has always had:
!   to the run log:
!     log_too_many_zones        more convection zones than the bounds
!                               array holds (find_convection_zones)
!     log_mixed_regions         the mixed regions without and with
!                               overshoot (extend_by_overshoot)
!     log_zones_merged          two mixed regions merged by overshoot or
!                               by semiconvection
!     log_semiconvective_edge   an edge moved by semiconvection
!     log_cubic_nonconvergence  the mixing-length cubic did not converge
!                               (the gradient is taken adiabatic)
!   to the terminal:
!     warn_unstable_radicand    the mixing-length square root of a
!                               negative number (the call returns ierr)
! The run log is run_messages' run_log_unit, written directly, as the
! kernels did (a test points it at standard output); the terminal is
! standard output. Each message is written into an array of lines with
! an internal write (one record per element) after every element is set
! to unwritten_line, and the lines the write reached are printed, so the
! output has the records the format makes, as a write to a unit would.
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

subroutine log_too_many_zones()
      character(len=message_length) :: lines(1)
      lines = unwritten_line
      write(lines, '(" -----TOO MANY MIX ZONES")')
      call write_convection_lines(run_log_unit, lines)
end subroutine log_too_many_zones

! The first num_zones regions of each list, four to a line.
subroutine log_mixed_regions(bounds_no_overshoot, bounds, num_zones)
      integer, intent(in) :: bounds_no_overshoot(:,:), bounds(:,:), num_zones
      character(len=message_length) :: lines(max(1, (num_zones + 3)/4))
      integer :: k, j
      lines = unwritten_line
      write(lines, '(1x,"MIXED REGIONS WITHOUT OVERSHOOT",4("[",i4,"-",i4," ]"))') &
           ((bounds_no_overshoot(k,j), j = 1, 2), k = 1, num_zones)
      call write_convection_lines(run_log_unit, lines)
      lines = unwritten_line
      write(lines, '(1x,"MIXED REGIONS WITH OVERSHOOT   ",4("[",i4,"-",i4," ]"))') &
           ((bounds(k,j), j = 1, 2), k = 1, num_zones)
      call write_convection_lines(run_log_unit, lines)
end subroutine log_mixed_regions

! The two regions before the merge (first and last shell of each, from
! the centre out) and the merged one.
subroutine log_zones_merged(cause, lower, upper)
      integer, intent(in) :: cause, lower(:), upper(:)
      character(len=message_length) :: lines(2)
      lines = unwritten_line
      if (cause == merged_by_overshoot) then
         write(lines, '(2x,"CONVECTION ZONES MERGED DUE TO OVERSHOOT"/2x, &
              &"OLD",2("[",i3,"-",i3,"]")," NEW","[",i3,"-",i3,"]")') &
              lower, upper, lower(1), upper(2)
      else
         write(lines, '(2x,"MIXED ZONES MERGED DUE TO SEMICONVECTION"/2x, &
              &"OLD",2("[",i3,"-",i3,"]")," NEW","[",i3,"-",i3,"]")') &
              lower, upper, lower(1), upper(2)
      end if
      call write_convection_lines(run_log_unit, lines)
end subroutine log_zones_merged

! An edge of a mixed region moved from old_edge to new_edge;
! reached_limit when the distance limit stopped it. The radiative
! gradients at the first shell outside the old edge, with the
! convection zone's composition and with the shell's own.
subroutine log_semiconvective_edge(old_edge, new_edge, reached_limit, &
     radiative_gradient_zone_mix, radiative_gradient_local)
      integer, intent(in) :: old_edge, new_edge
      logical, intent(in) :: reached_limit
      real(dp), intent(in) :: radiative_gradient_zone_mix, radiative_gradient_local
      character(len=message_length) :: lines(1)
      lines = unwritten_line
      write(lines, '(1x,"CZ OLD EDGE ",i3," EXTENDED TO-",i3," LIMIT=",l1, &
           &" RAD.GRADS-IN/OUT",2f8.4)') old_edge, new_edge, reached_limit, &
           radiative_gradient_zone_mix, radiative_gradient_local
      call write_convection_lines(run_log_unit, lines)
end subroutine log_semiconvective_edge

subroutine log_cubic_nonconvergence(log_pressure, log_temperature, opacity, &
     specific_heat_cp, dlnrho_dlnt)
      real(dp), intent(in) :: log_pressure, log_temperature, opacity, &
           specific_heat_cp, dlnrho_dlnt
      character(len=message_length) :: lines(2)
      lines = unwritten_line
      write(lines, '(" -----CUBIC NON-CONVERGENCE(PL,TL,CAPPA,CP,QDT)=",2f10.6/1p3e12.6)') &
           log_pressure, log_temperature, opacity, specific_heat_cp, dlnrho_dlnt
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
      write(lines, '("***WARNING***"/"TPGRAD TRIED TO TAKE SQUARE ROOT", &
           &" OF NEGATIVE NUMBER - ASSUMED ADIABATIC CONVECTION"/ &
           &" DELDEL ",1pe12.3," QDT ",e12.3," DISC. ",e12.3)') &
           superadiabatic_excess, dlnrho_dlnt, radicand
      write(lines(4), '(1p10e12.3)') opacity, luminosity_lsun, log_pressure, log_mass, &
           log_temperature, temperature_rotation_factor, pressure_rotation_factor, &
           radiative_gradient, adiabatic_gradient, log_radius
      call write_convection_lines(output_unit, lines)
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
