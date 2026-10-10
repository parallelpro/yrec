!----------------------------------------------------------------------
! light_element_lib
!----------------------------------------------------------------------
! The light elements, operator-split from the implicit network
! (2026-10: moved here from core/burn_lib with the model state as
! arguments). Deuterium: deuterium_rate, the per-shell rate with the
! convective-overturn cap; dburn, the exponential sub-stepped burn over
! a step (burn_deuterium_region, for the main step and the rotation
! sub-steps). Li6/Li7/Be9: li_be_rates, the
! proton-capture rates (the star-layer store_li_be_rates/burn_li_be store and
! integrate them). They read configuration from star%job/star%ctrl and
! nothing else from star.
module light_element_lib
      use net_lib, only: deuterium_proton_capture_rate, avogadro_dpg, deuteron_mass_amu
      implicit none
      private
      public :: deuterium_rate, burn_deuterium_region, li_be_rates, &
           min_deuterium_step, min_deuterium_substep

! burn_deuterium_region's thresholds as its two callers have always used
! them: mix's step (min_deuterium_step) and the rotation sub-steps
! (min_deuterium_substep), deuterium mass fraction.
      double precision, parameter :: min_deuterium_step = 1.0d-11
      double precision, parameter :: min_deuterium_substep = 1.0d-14

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
! was T9**(+2/3) (deutrate.f). 2026-10: the rate is deuterium_proton_capture_rate, shared
! with compute_zone_energy; this copy used to divide N_A by A = 2 (3.0115D23) where
! compute_zone_energy divides by the deuteron mass, 0.68% apart.
      rdeut = deuterium_proton_capture_rate(rho, t9, t9p13, t9p23, t9m13, t9m23)
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
! burn_deuterium_region (formerly dburn and dburnm)
!----------------------------------------------------------------------
! Deuterium burning over a step dt_gyr in one radiative shell
! (zone_begin = zone_end) or one convection zone (homogenized), by
! exponential sub-steps with the rate interpolated linearly between its
! start- and end-of-step values (rate_start_zone, rate_end_zone, per
! Gyr, from deuterium_rate). D + p -> He3: X changes by D/2 and He3 by
! -3D/2 of the deuterium burned. Under mass accretion the outermost
! region blends in accreted_mass_step grams of gas of composition
! ctrl%accreted_composition, mass-weighted. Below min_deuterium (after
! that blend) the region's deuterium is set to zero.
!
! 2026-10: one routine for the main step (mix's burn_deuterium, dt in
! Gyr, min_deuterium = 1e-11) and the rotation sub-steps
! (mid_timestep_model, dt converted from s, the sub-step's share of the
! accreted mass, min_deuterium = 1e-14) -- formerly dburn and dburnm,
! two copies that differed in units, in dburnm's grams-plus-one
! accretion blend and in its accretion-blind threshold test.
subroutine burn_deuterium_region(zone_begin, zone_end, num_zones, shell_mass, &
     composition, timestep, xa_start, rate_start_zone, rate_end_zone, &
     accreted_mass, min_deuterium)

      use star_info_lib, only: star, json, i_h2
      use math_lib
      implicit none

      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: timestep
! timestep in Gyr; the start-of-step composition; the deuterium rates
! per shell at the start and end of the step; the mass accreted over
! this step, g; the threshold below which deuterium is zeroed
      double precision, intent(in) :: xa_start(15,json), &
           rate_start_zone(json), rate_end_zone(json), accreted_mass, &
           min_deuterium
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
! the region mass weights the accreted material below (2026-10: it was
! left unset here, i.e. zero, so a lone accreting surface shell took
! the accreted deuterium and dropped its own)
         total_shell_mass = shell_mass(zone_begin)
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
      if (deuterium_fraction_test.lt.min_deuterium) then
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
end subroutine burn_deuterium_region


!----------------------------------------------------------------------
! li_be_rates
!----------------------------------------------------------------------
! The Li6, Li7 and Be9 proton-capture rates at T9, per unit rho*X
! (Caughlan & Fowler 1988; no screening, no energy release), scaled by
! the deck's ctrl%li6/li7/be9_*_rate_scale: Li6(p,He3)He4, Li7(p,a)He4,
! and the sum of Be9(p,g)B10, Be9(p,d)2He4 and Be9(p,a)Li6. Moved from
! burn_lib's store_li_be_rates (2026-10), formulas unchanged.
subroutine li_be_rates(t9, fli6, fli7, fsbe9)
      use star_info_lib, only: star
      use phys_const_lib
      use math_lib
      double precision, intent(in) :: t9
      double precision, intent(out) :: fli6, fli7, fsbe9
      double precision :: t913, t923, t943, t953, t932, t934
      double precision :: fbe91, fbe92, fbe93, fx, ex, t9a, c56

         t913=pow(t9, cc13)
         t923=t913*t913
         t943=t923*t923
! MHP 10/91 ADDED DEFINITION
         t953 = t943*t913
         t932=pow(t9, 1.5d0)
         t934=pow(t9, 7.5d-1)
!
! LI6(P,HE3)ALPHA
!
!           FLI6=3.73D10/T923*EXP(-8.413/T913-(T9/5.50)**2)*
!    1      (1.+%TYPO%0.50*T913-0.061*T923-0.021*T9+0.006*T943+0.006*T953)
!    2      +1.33D10/T932*EXP(-17.763/T9)+1.29D9/T9*EXP(-21.820/T9)
! ***TYPO IN OLD RATE***
         fli6=3.73d10/t923*exp(-8.413d0/t913-(t9/5.50d0)**2)* &
         (1.d0+0.050d0*t913-0.061d0*t923-0.021d0*t9+0.006d0*t943 &
          +0.005d0*t953)+1.33d10/t932*exp(-17.763d0/t9) &
          +1.29d9/t9*exp(-21.820d0/t9)
!
! LI7(P,ALFA)HE4
!
!           FLI7=8.04D08/T923*EXP(-8.471/T913-(T9/30.068)**2)*
!    1      (1.+0.049*T913+0.230*T923+0.079*T9-0.027*T943-0.023*T953)
!    2      +1.54D06/T932*EXP(-4.479/T9)+1.07D10/T932*EXP(-30.443/T9)
         t9a = t9/(1.0d0+0.759d0*t9)
         c56 = 1.25d0*cc23
         fli7=1.096d9/t923*exp(-8.472d0/t913)-4.830d8*pow(t9a, c56)/ &
         t932*exp(-8.472d0/pow(t9a, cc13))+1.06d10/t932*exp(-30.442d0/t9)
!
! BE9(P,GAMMA)B10
!
         fbe91=1.33d+07/t923*exp(-10.359d0/t913-(t9/0.846d0)**2) &
         *(1.0d0+4.0d-2*t913+1.52d0*t923+4.28d-1*t9+2.15d0*t943+ &
         1.54d0*t953) + 9.64d+4/t932*exp(-3.445d0/t9) + &
         2.72d+6/t932*exp(-10.62d0/t9)
!    3      2.76D+6/T932*EXP(-10.62D0/T9) %OLD LAST LINE%
!
! BE9(P,D)2HE4 - UNCHANGED.
!
         fx=2.11d+11/t923*exp(-10.359d0/t913-(t9/5.2d-1)**2) &
         *(1.0d0+4.0d-2*t913+1.09d0*t923+3.07d-1*t9+3.21d0*t943 &
         +2.3d0*t953)
         ex = exp(-3.046d0/t9)/t9
         fbe92=fx+5.79d+8*ex+8.5d+8/t934*exp(-5.8d0/t9)
!
! BE9(P,ALPHA)LI6
!
         fbe93=fx+4.51d+8*ex+6.7d+8/t934*exp(-5.16d0/t9)
! G Somers 6/14, SCALE BY THE NEW CROSS SECTIONS
         fli6=star%ctrl%li6_rate_scale*fli6
         fli7=star%ctrl%li7_rate_scale*fli7
         fbe91=star%ctrl%be9_pg_rate_scale*fbe91
         fbe92=star%ctrl%be9_pd_rate_scale*fbe92
         fbe93=star%ctrl%be9_palpha_rate_scale*fbe93
! G Somers END
! SUM RATES
         fsbe9 = fbe91 + fbe92 + fbe93
end subroutine li_be_rates

end module light_element_lib
