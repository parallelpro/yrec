!----------------------------------------------------------------------
! dburn
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original dburn.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Compute the abundance changes resulting from deuterium burning, by
! sub-stepping the burning rate across the timestep (capped so no
! single sub-step burns more than 0.5 in natural log, keeping the
! discretization error under ~1%) rather than solving an ODE for
! deuterium directly. If mass accretion is active, accreted and
! pre-existing deuterium are burned separately (accreted matter is
! only exposed for ~half the timestep on average) and mass-weighted
! back together.
subroutine dburn(zone_begin, zone_end, num_zones, shell_mass, &
     composition, timestep)

      use light_burn_lib
      use oldmod_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: timestep





      save

      double precision :: hydrogen_fraction, deuterium_fraction, &
           helium3_fraction, rate_start, rate_end
      double precision :: total_shell_mass, rate_start_sum, rate_end_sum
      integer :: zone_idx, substep_idx
      double precision :: deuterium_fraction_test
      double precision :: burning_rate, num_substeps_real, substep_dt, &
           rate_accum, rate_increment, deuterium_fraction_burned
      integer :: num_substeps
      double precision :: accreted_deuterium_fraction, &
           accreted_deuterium_burned, newly_accreted_deuterium
      double precision :: deuterium_change_original, &
           deuterium_change_accreted, deuterium_fraction_new, &
           hydrogen_fraction_new, helium3_fraction_new, deuterium_change

      if (zone_begin.eq.zone_end) then
         hydrogen_fraction = prev_model%old_composition(1,zone_begin)
         deuterium_fraction = prev_model%old_composition(12,zone_begin)
         helium3_fraction = prev_model%old_composition(4,zone_begin)
         rate_start = light_burn%deuterium_burning_rate_start(zone_begin)
         rate_end = light_burn%deuterium_burning_rate(zone_begin)
      else
         total_shell_mass = 0.0d0
         rate_start_sum = 0.0d0
         rate_end_sum = 0.0d0
         hydrogen_fraction = 0.0d0
         deuterium_fraction = 0.0d0
         helium3_fraction = 0.0d0
         do zone_idx = zone_begin, zone_end
            total_shell_mass = total_shell_mass + shell_mass(zone_idx)
            rate_start_sum = rate_start_sum + &
                 shell_mass(zone_idx)*light_burn%deuterium_burning_rate_start(zone_idx)
            rate_end_sum = rate_end_sum + &
                 shell_mass(zone_idx)*light_burn%deuterium_burning_rate(zone_idx)
            hydrogen_fraction = hydrogen_fraction + &
                 prev_model%old_composition(1,zone_idx)*shell_mass(zone_idx)
            deuterium_fraction = deuterium_fraction + &
                 prev_model%old_composition(12,zone_idx)*shell_mass(zone_idx)
            helium3_fraction = helium3_fraction + &
                 prev_model%old_composition(4,zone_idx)*shell_mass(zone_idx)
         end do
         rate_start = rate_start_sum/total_shell_mass
         rate_end = rate_end_sum/total_shell_mass
         hydrogen_fraction = hydrogen_fraction/total_shell_mass
         deuterium_fraction = deuterium_fraction/total_shell_mass
         helium3_fraction = helium3_fraction/total_shell_mass
      end if
      if (use_mass_accretion .and. zone_end.eq.num_zones .and. &
           mass_accretion_rate.gt.0.0d0) then
         deuterium_fraction_test = (deuterium_fraction*total_shell_mass + &
              accreted_composition(12)*light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + light_burn%accreted_mass_fraction)
      else
         deuterium_fraction_test = deuterium_fraction
      end if
      if (deuterium_fraction_test.lt.1.0d-11) then
         do zone_idx = zone_begin, zone_end
            composition(12,zone_idx) = 0.0d0
         end do
         goto 9999
      end if
!     BURN DEUTERIUM IN A SERIES OF SMALL STEPS,
!     INCREMENTING THE BURNING RATE FROM THE STARTING
!     TO ENDING VALUES.
!     RESTRICT THE MAXIMUM BURNING TO 0.5 IN THE
!     NATURAL LOG (<1% ERROR FROM A SIMPLE SUM
!     OF DISCRETE STEPS OF THIS SIZE).
!     MHP 10/02 INITIALIZE RATE
      burning_rate = max(rate_start, rate_end)
      num_substeps = int(burning_rate*timestep/0.1d0) + 1
      num_substeps_real = dfloat(num_substeps)
      substep_dt = timestep/num_substeps_real
      rate_accum = rate_start
      rate_increment = (rate_end - rate_start)/num_substeps_real
      deuterium_fraction_burned = deuterium_fraction
      do substep_idx = 1, num_substeps
         burning_rate = rate_accum + 0.5d0*rate_increment
         deuterium_fraction_burned = deuterium_fraction_burned* &
              exp(-2.0d0*burning_rate*substep_dt)
         rate_accum = rate_accum + rate_increment
      end do
!     INCLUDE MASS ACCRETION FROM DEUTERIUM BURNING
      if (use_mass_accretion .and. zone_end.eq.num_zones) then
!        ACCRETED MATTER IS EXPOSED TO BURNING FOR, ON
!        AVERAGE. 1/2 OF THE TIMESTEP.  light_burn%accreted_mass_fraction IS
!        DEFINED AS DMDT*DT/ORIGINAL CZ MASS.
!        BURN BOTH THE ACCRETED D AND THE ORIGINAL D
!        SEPARATELY AND FIND THE NEW MASS-WEIGHTED
!        AVERAGE D.  NOTE THAT ACCRETION OF ELEMENTS
!        1-11 HAS ALREADY BEEN TREATED.
         rate_accum = rate_start
         accreted_deuterium_fraction = accreted_composition(12)
         accreted_deuterium_burned = 0.0d0
         do substep_idx = 1, num_substeps
            burning_rate = rate_accum + 0.5d0*rate_increment
!           FULLY BURN PRIOR ACCRETED MATTER
            accreted_deuterium_burned = accreted_deuterium_burned* &
                 exp(-2.0d0*burning_rate*substep_dt)
!           EXPOSE NEWLY ACCRETED MATTER TO 1/2 BURNING
            newly_accreted_deuterium = accreted_deuterium_fraction/ &
                 num_substeps_real*exp(-burning_rate*substep_dt)
            accreted_deuterium_burned = accreted_deuterium_burned + &
                 newly_accreted_deuterium
            rate_accum = rate_accum + rate_increment
         end do
!        NOW PERFORM A WEIGHTED SUM OF THE ORIGINAL
!        AND NEWLY ACCRETED MATTER
         deuterium_change_original = deuterium_fraction_burned - &
              deuterium_fraction
         deuterium_change_accreted = accreted_deuterium_burned - &
              accreted_composition(12)
         deuterium_fraction_new = (deuterium_fraction_burned*total_shell_mass &
              + accreted_deuterium_burned*light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + light_burn%accreted_mass_fraction)
         hydrogen_fraction_new = hydrogen_fraction + 0.5d0* &
              (deuterium_change_original*total_shell_mass + &
              deuterium_change_accreted*light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + light_burn%accreted_mass_fraction)
         helium3_fraction_new = helium3_fraction - 1.5d0* &
              (deuterium_change_original*total_shell_mass + &
              deuterium_change_accreted*light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + light_burn%accreted_mass_fraction)
      else
!        INCREMENT H,D,HE3 WITHOUT MASS ACCRETION
         deuterium_change = deuterium_fraction_burned - deuterium_fraction
         deuterium_fraction_new = deuterium_fraction_burned
         hydrogen_fraction_new = hydrogen_fraction + 0.5d0*deuterium_change
         helium3_fraction_new = helium3_fraction - 1.5d0*deuterium_change
      end if
      do zone_idx = zone_begin, zone_end
         composition(1,zone_idx) = hydrogen_fraction_new
         composition(4,zone_idx) = helium3_fraction_new
         composition(12,zone_idx) = deuterium_fraction_new
      end do
 9999 continue
      return
end subroutine dburn
