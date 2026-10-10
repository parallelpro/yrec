!----------------------------------------------------------------------
! light_element_lib
!----------------------------------------------------------------------
! Deuterium burning, operator-split from the implicit network
! (2026-10: moved here from core/burn_lib with the model state as
! arguments): deuterium_rate, the per-shell rate with the
! convective-overturn cap; dburn, the exponential sub-stepped burn over
! a step in Gyr (mix); dburnm, the same over a rotation sub-step in
! seconds (mid_timestep_model). They read configuration (accretion
! settings) from star%job/star%ctrl and nothing else from star.
module light_element_lib
      use net_lib, only: dpg_rate, avogadro_dpg, deuteron_mass_amu
      implicit none
      private
      public :: deuterium_rate, dburn, dburnm

! dburn/dburnm interface constants (2026 W3). The two deuterium
! drivers are deliberately kept as separate bodies (see the dburn
! header table); the literals that distinguish their interfaces are
! named here so the split is visible in one place.
! dburn_min_deuterium: dburn's "nothing to burn" threshold on the
! accretion-mixed deuterium mass fraction.
      double precision, parameter :: dburn_min_deuterium = 1.0d-11
! dburnm_min_deuterium: dburnm's threshold on the plain CZ-averaged
! deuterium mass fraction.
      double precision, parameter :: dburnm_min_deuterium = 1.0d-14
! gyr_per_year: dburnm's seconds -> Gyr timestep conversion factor,
! applied as timestep*gyr_per_year/seconds_per_year.
      double precision, parameter :: gyr_per_year = 1.0d-9

contains

!----------------------------------------------------------------------
! deuterium_rate (formerly deutrate)
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original deutrate.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! The rate of nonequilibrium deuterium burning at one shell, per Gyr
! (X * N_A<sigma v> rho / m_d, excluding the deuterium abundance), for
! dburn/dburnm. In the surface convection zone (in_surface_cz) it is
! capped so that deuterium cannot burn faster than the convective
! overturn (turnover_timescale, s, when > 1). 2026-10: deutrate
! stored the result in star%deuterium_burning_rate(_start) itself; the
! caller (mix) stores it now.
double precision function deuterium_rate(log_density, log_temperature, &
     hydrogen_fraction, in_surface_cz, turnover_timescale)
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log_density, log_temperature, &
           hydrogen_fraction
      logical, intent(in) :: in_surface_cz
      double precision, intent(in) :: turnover_timescale

! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
!  POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
!  OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
!  CGS UNITS.
      double precision :: rho, t9, t9p13, t9p23, t9m13, t9m23, t9m1
      double precision :: rdeut, rdeutmax, rdeut2

      rho=exp(ln10*log_density)
      t9 = exp(ln10*(log_temperature - 9.0d0))
      t9p13 = pow(t9, cc13)
      t9p23 = t9p13**2
      t9m13=1.0d0/t9p13
      t9m23=t9m13**2
      t9m1=1.0d0/t9
! MHP 5/02 ADD DEUTERIUM BURNING TERM TO THE CODE
! IF DEUTERIUM IS ABOVE A MINIMUM THRESHOLD VALUE.
! RDEUT IS THE RATE (EXCLUDING FACTORS OF THE
! ABUNDANCES) AND QRTDEUT IS THE DERIVATIVE W/R/T T.
! NOTE THAT SCREENING IS EXCLUDED - REASONABLE GIVEN
! THE LOW TEMPERATURES INVOLVED.
! 2026 (bugsweep Batch 2): T9**(-2/3), as in CF88's d(p,g)3He fit;
! was T9**(+2/3) (deutrate.f). 2026-10: the rate is dpg_rate, shared
! with engeb; this copy used to divide N_A by A = 2 (3.0115D23) where
! engeb divides by the deuteron mass, 0.68% apart.
      rdeut = dpg_rate(rho, t9, t9p13, t9p23, t9m13, t9m23)
! NOW LIMIT DEUTERIUM BURNING IN A SURFACE CZ TO BE ON A TIME SCALE
! NO SHORTER THAN THE CONVECTIVE OVERTURN TIMESCALE.
      if(in_surface_cz .and. turnover_timescale.gt.1.0d0)then
         rdeutmax = avogadro_dpg/deuteron_mass_amu/turnover_timescale
         rdeut2 = rdeut*hydrogen_fraction
         if(rdeut2.gt.rdeutmax)then
            if(hydrogen_fraction.gt.1.0d-6)then
               rdeut = rdeutmax/hydrogen_fraction
            endif
         endif
      endif
      deuterium_rate = hydrogen_fraction*rdeut*gyr_amu_per_sec_gram
end function deuterium_rate

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
!
! dburn and dburnm are near-twins with the same sub-stepping body.
! They are kept as two routines on purpose: unifying them would change
! the numbers (units, thresholds and weighting differ). The four
! interface differences, side by side:
!
!   difference        dburn                      dburnm
!   ----------------  -------------------------  ---------------------------
!   timestep units    used as passed (Gyr)       seconds, converted to Gyr
!                                                (timestep_gyr =
!                                                timestep*gyr_per_year/
!                                                seconds_per_year)
!   "nothing to       accretion-mixed            plain CZ average
!   burn" threshold   deuterium_fraction_test    deuterium_fraction
!                     .lt. dburn_min_deuterium   .lt. dburnm_min_deuterium
!                     (1.0d-11)                  (1.0d-14)
!   accretion         absolute masses:           normalised to the CZ mass:
!   weighting         (D*M_cz + D_acc*f)/        (D + D_acc*f_sub)/
!                     (M_cz + f) with            (1 + f_sub) with f_sub =
!                     M_cz = total_shell_mass,   step_fraction*
!                     f = star%accreted_         accreted_mass
!                     mass_fraction
!   rate source       star%deuterium_burning_    dummy arguments
!                     rate_start / star%         deuterium_rate_start /
!                     deuterium_burning_rate     deuterium_rate_end
!
! The named constants live at the top of burn_lib.
subroutine dburn(zone_begin, zone_end, num_zones, shell_mass, &
     composition, timestep, xa_start, rate_start_zone, rate_end_zone, &
     accreted_mass)

      use star_info_lib, only: star, json, i_h2
      use math_lib
      implicit none

      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: timestep
! model state, passed in since the move to net/ (2026-10): the
! start-of-step composition, the deuterium rates per shell at the start
! and end of the step (deuterium_rate), and the mass accreted this step
! (g; accreted_mass)
      double precision, intent(in) :: xa_start(15,json), &
           rate_start_zone(json), rate_end_zone(json), accreted_mass
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
         hydrogen_fraction = xa_start(1,zone_begin)
         deuterium_fraction = xa_start(12,zone_begin)
         helium3_fraction = xa_start(4,zone_begin)
         rate_start = rate_start_zone(zone_begin)
         rate_end = rate_end_zone(zone_begin)
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
                 shell_mass(zone_idx)*rate_start_zone(zone_idx)
            rate_end_sum = rate_end_sum + &
                 shell_mass(zone_idx)*rate_end_zone(zone_idx)
            hydrogen_fraction = hydrogen_fraction + &
                 xa_start(1,zone_idx)*shell_mass(zone_idx)
            deuterium_fraction = deuterium_fraction + &
                 xa_start(12,zone_idx)*shell_mass(zone_idx)
            helium3_fraction = helium3_fraction + &
                 xa_start(4,zone_idx)*shell_mass(zone_idx)
         end do
         rate_start = rate_start_sum/total_shell_mass
         rate_end = rate_end_sum/total_shell_mass
         hydrogen_fraction = hydrogen_fraction/total_shell_mass
         deuterium_fraction = deuterium_fraction/total_shell_mass
         helium3_fraction = helium3_fraction/total_shell_mass
      end if
      if (star%job%use_mass_accretion .and. zone_end.eq.num_zones .and. &
           star%ctrl%mass_accretion_rate.gt.0.0d0) then
         deuterium_fraction_test = (deuterium_fraction*total_shell_mass + &
              star%ctrl%accreted_composition(i_h2)*accreted_mass)/ &
              (total_shell_mass + accreted_mass)
      else
         deuterium_fraction_test = deuterium_fraction
      end if
      if (deuterium_fraction_test.lt.dburn_min_deuterium) then
         do zone_idx = zone_begin, zone_end
            composition(i_h2,zone_idx) = 0.0d0
         end do
         return
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
      if (star%job%use_mass_accretion .and. zone_end.eq.num_zones) then
!        ACCRETED MATTER IS EXPOSED TO BURNING FOR, ON
!        AVERAGE. 1/2 OF THE TIMESTEP.  accreted_mass IS
!        DEFINED AS DMDT*DT/ORIGINAL CZ MASS.
!        BURN BOTH THE ACCRETED D AND THE ORIGINAL D
!        SEPARATELY AND FIND THE NEW MASS-WEIGHTED
!        AVERAGE D.  NOTE THAT ACCRETION OF ELEMENTS
!        1-11 HAS ALREADY BEEN TREATED.
         rate_accum = rate_start
         accreted_deuterium_fraction = star%ctrl%accreted_composition(i_h2)
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
              star%ctrl%accreted_composition(i_h2)
         deuterium_fraction_new = (deuterium_fraction_burned*total_shell_mass &
              + accreted_deuterium_burned*accreted_mass)/ &
              (total_shell_mass + accreted_mass)
         hydrogen_fraction_new = hydrogen_fraction + 0.5d0* &
              (deuterium_change_original*total_shell_mass + &
              deuterium_change_accreted*accreted_mass)/ &
              (total_shell_mass + accreted_mass)
         helium3_fraction_new = helium3_fraction - 1.5d0* &
              (deuterium_change_original*total_shell_mass + &
              deuterium_change_accreted*accreted_mass)/ &
              (total_shell_mass + accreted_mass)
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
         composition(i_h2,zone_idx) = deuterium_fraction_new
      end do
      return
end subroutine dburn

!----------------------------------------------------------------------
! dburnm
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original dburnm.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Near-duplicate of dburn (compute the abundance changes resulting
! from deuterium burning by sub-stepping the burning rate across the
! timestep) used from the sub-stepped call site: here the start/end
! burning rates are passed in explicitly as dummy arguments
! (deuterium_rate_start/deuterium_rate_end) rather than read from
! star%deuterium_burning_rate(_start), and the accreted-mass fraction
! used for this call is only a fraction (step_fraction) of the
! full-timestep accreted_mass. See the dburn header for
! the full list of differences and for the sub-stepping algorithm.
subroutine dburnm(zone_begin, zone_end, num_zones, shell_mass, &
     composition, timestep, deuterium_rate_end, deuterium_rate_start, &
     step_fraction, xa_start, accreted_mass)
      use star_info_lib, only: star, json, i_h2
      use phys_const_lib
      use math_lib
      implicit none

      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: timestep
      double precision, intent(in) :: deuterium_rate_end(json), &
           deuterium_rate_start(json)
      double precision, intent(in) :: step_fraction
! the start-of-step composition and the mass accreted this step (g),
! passed in since the move to net/ (2026-10)
      double precision, intent(in) :: xa_start(15,json), accreted_mass
      double precision :: hydrogen_fraction, deuterium_fraction, &
           helium3_fraction, rate_start, rate_end
      double precision :: total_shell_mass, rate_start_sum, rate_end_sum
      integer :: zone_idx, substep_idx
      double precision :: timestep_gyr
      double precision :: burning_rate, num_substeps_real, substep_dt, &
           rate_accum, rate_increment, deuterium_fraction_burned
      integer :: num_substeps
      double precision :: accreted_deuterium_fraction, &
           accreted_deuterium_burned, newly_accreted_deuterium
      double precision :: accreted_mass_fraction_substep
      double precision :: deuterium_change_original, &
           deuterium_change_accreted, deuterium_fraction_new, &
           hydrogen_fraction_new, helium3_fraction_new, deuterium_change

! timestep_gyr converts timestep (assumed in seconds) to Gyr, matching
! the units used elsewhere for the burning rates (dburn receives Gyr
! directly; see the table in the dburn header).
      timestep_gyr = timestep*gyr_per_year/seconds_per_year
      if(zone_begin.eq.zone_end)then
         hydrogen_fraction = xa_start(1,zone_begin)
         deuterium_fraction = xa_start(12,zone_begin)
         helium3_fraction = xa_start(4,zone_begin)
         rate_start = deuterium_rate_start(zone_begin)
         rate_end = deuterium_rate_end(zone_begin)
      else
         total_shell_mass = 0.0d0
         rate_start_sum = 0.0d0
         rate_end_sum = 0.0d0
         hydrogen_fraction = 0.0d0
         deuterium_fraction = 0.0d0
         helium3_fraction = 0.0d0
         do zone_idx = zone_begin,zone_end
            total_shell_mass = total_shell_mass + shell_mass(zone_idx)
            rate_start_sum = rate_start_sum + &
                 shell_mass(zone_idx)*deuterium_rate_start(zone_idx)
            rate_end_sum = rate_end_sum + &
                 shell_mass(zone_idx)*deuterium_rate_end(zone_idx)
            hydrogen_fraction = hydrogen_fraction + &
                 xa_start(1,zone_idx)*shell_mass(zone_idx)
            deuterium_fraction = deuterium_fraction + &
                 xa_start(12,zone_idx)*shell_mass(zone_idx)
            helium3_fraction = helium3_fraction + &
                 xa_start(4,zone_idx)*shell_mass(zone_idx)
         end do
         rate_start = rate_start_sum/total_shell_mass
         rate_end = rate_end_sum/total_shell_mass
         hydrogen_fraction = hydrogen_fraction/total_shell_mass
         deuterium_fraction = deuterium_fraction/total_shell_mass
         helium3_fraction = helium3_fraction/total_shell_mass
      endif
      if(deuterium_fraction.lt.dburnm_min_deuterium)then
         do zone_idx = zone_begin,zone_end
            composition(i_h2,zone_idx) = 0.0d0
         end do
         return
      endif
! BURN DEUTERIUM IN A SERIES OF SMALL STEPS,
! INCREMENTING THE BURNING RATE FROM THE STARTING
! TO ENDING VALUES.
! RESTRICT THE MAXIMUM BURNING TO 0.5 IN THE
! NATURAL LOG (<1% ERROR FROM A SIMPLE SUM
! OF DISCRETE STEPS OF THIS SIZE).
! MHP 10/02 INITIALIZE RATE
      burning_rate = max(rate_start,rate_end)
      num_substeps = int(burning_rate*timestep_gyr/0.1d0)+1
      num_substeps_real = dfloat(num_substeps)
      substep_dt = timestep_gyr/num_substeps_real
      rate_accum = rate_start
      rate_increment = (rate_end-rate_start)/num_substeps_real
      deuterium_fraction_burned = deuterium_fraction
      do substep_idx = 1,num_substeps
         burning_rate = rate_accum + 0.5d0*rate_increment
         deuterium_fraction_burned = deuterium_fraction_burned* &
              exp(-2.0d0*burning_rate*substep_dt)
         rate_accum = rate_accum + rate_increment
      end do
! INCLUDE MASS ACCRETION FROM DEUTERIUM BURNING
      if(star%job%use_mass_accretion .and. zone_end.eq.num_zones)then
! ACCRETED MATTER IS EXPOSED TO BURNING FOR, ON
! AVERAGE. 1/2 OF THE TIMESTEP.  accreted_mass IS
! DEFINED AS DMDT*DT/ORIGINAL CZ MASS.
! BURN BOTH THE ACCRETED D AND THE ORIGINAL D
! SEPARATELY AND FIND THE NEW MASS-WEIGHTED
! AVERAGE D.  NOTE THAT ACCRETION OF ELEMENTS
! 1-11 HAS ALREADY BEEN TREATED.
         rate_accum = rate_start
         accreted_deuterium_fraction = star%ctrl%accreted_composition(i_h2)
         accreted_deuterium_burned = 0.0d0
         do substep_idx = 1,num_substeps
            burning_rate = rate_accum + 0.5d0*rate_increment
! FULLY BURN PRIOR ACCRETED MATTER
            accreted_deuterium_burned = accreted_deuterium_burned* &
                 exp(-2.0d0*burning_rate*substep_dt)
! EXPOSE NEWLY ACCRETED MATTER TO 1/2 BURNING
            newly_accreted_deuterium = accreted_deuterium_fraction/ &
                 num_substeps_real*exp(-burning_rate*substep_dt)
            accreted_deuterium_burned = accreted_deuterium_burned + &
                 newly_accreted_deuterium
            rate_accum = rate_accum + rate_increment
         end do
! NOW PERFORM A WEIGHTED SUM OF THE ORIGINAL
! AND NEWLY ACCRETED MATTER
! UNLIKE THE ROUTINE CALLED BY MIX, THE TIMESTEP
! HERE IS ONLY A PORTION OF THE TOTAL MODEL DT.
! DO ONLY THE SMALLER TIMESTEP PORTION OF THE
! TOTAL MASS ACCRETION.
         accreted_mass_fraction_substep = step_fraction*accreted_mass
         deuterium_change_original = deuterium_fraction_burned - &
              deuterium_fraction
         deuterium_change_accreted = accreted_deuterium_burned - &
              star%ctrl%accreted_composition(i_h2)
         deuterium_fraction_new = (deuterium_fraction_burned+ &
              accreted_deuterium_burned*accreted_mass_fraction_substep)/ &
              (1.0d0+accreted_mass_fraction_substep)
         hydrogen_fraction_new = hydrogen_fraction + 0.5d0* &
              (deuterium_change_original+ &
              deuterium_change_accreted*accreted_mass_fraction_substep)/ &
              (1.0d0+accreted_mass_fraction_substep)
         helium3_fraction_new = helium3_fraction - 1.5d0* &
              (deuterium_change_original+ &
              deuterium_change_accreted*accreted_mass_fraction_substep)/ &
              (1.0d0+accreted_mass_fraction_substep)
      else
! INCREMENT H,D,HE3 WITHOUT MASS ACCRETION
         deuterium_change = deuterium_fraction_burned - deuterium_fraction
         deuterium_fraction_new = deuterium_fraction_burned
         hydrogen_fraction_new = hydrogen_fraction + 0.5d0*deuterium_change
         helium3_fraction_new = helium3_fraction - 1.5d0*deuterium_change
      endif
      do zone_idx = zone_begin,zone_end
         composition(1,zone_idx) = hydrogen_fraction_new
         composition(4,zone_idx) = helium3_fraction_new
         composition(i_h2,zone_idx) = deuterium_fraction_new
      end do
      return
end subroutine dburnm

end module light_element_lib
