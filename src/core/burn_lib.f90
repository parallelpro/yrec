!----------------------------------------------------------------------
! burn_lib
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass -- ROADMAP.md "Decoupling the
! physics domains from the model"). The BURN DRIVERS, split out of
! net_lib: the routines that advance the model's composition and
! diagnostics over a timestep -- engeb (the per-zone energy-generation
! driver, with its eqburn equilibrium helper), the deuterium drivers
! dburn/dburnm and rate bookkeeping deutrate, and the light-element
! drivers liburn/liburn2 with lirate88. All of them read/write
! star_info (the previous model, the light-burn rate store, the flux
! and energy diagnostics), which is exactly why they are STAR-LAYER
! code: in MESA terms these are struct_burn_mix, not net. net_lib
! keeps the pure kernels (rates, sneut, nulosses, neutrino, azbar,
! safedivexp, the Fermi inverses), which are functions of
! (logRho, logT, composition, controls) only -- the surface
! test_nuclear pins. Procedure bodies are unchanged; only the module
! boundary moved.
module burn_lib
      use net_lib
      implicit none
contains



!----------------------------------------------------------------------
! eqburn
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqburn.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Computes explicit hydrogen- and helium-burning rates (dX/dt, dY/dt,
! dC/dt, dO/dt) assuming equilibrium He3 and CN-cycle abundances --
! the classic approximation that lets YREC track only a handful of
! "slow" abundances instead of every CNO isotope individually (see
! dburn.f90 for the analogous deuterium-burning treatment). Used as
! the initial guess for the burning rates at the start of a timestep;
! may be supplemented elsewhere by a fully implicit non-equilibrium
! calculation.
subroutine eqburn(rate_pp, rate_he3_he3, rate_he3_he4, rate_c12_p, &
     rate_c13_p, rate_n14_p, rate_o16_p, rate_c12_alpha, &
     rate_triple_alpha, shell_mass, shell_temperature, zone_begin, &
     zone_end, dc_dt, do_dt, dx_dt, dy_dt, equilibrium_xc12, &
     equilibrium_xo16, hydrogen_fraction, metal_fraction)

      use star_info_lib, only: star, json
      use const_lib
      implicit none

      double precision, intent(in) :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c12_alpha(json), &
           rate_triple_alpha(json)
      double precision, intent(in) :: shell_mass(json), &
           shell_temperature(json)
      integer, intent(in) :: zone_begin, zone_end
      double precision, intent(out) :: dc_dt, do_dt, dx_dt, dy_dt, &
           equilibrium_xc12, equilibrium_xo16, hydrogen_fraction, &
           metal_fraction



      double precision :: zone_avg_abundance(11)
      double precision :: total_shell_mass
      integer :: species_idx, zone_idx
      double precision :: pp_rate, he3_he3_rate, he3_he4_rate, c12_p_rate, &
           c13_p_rate, n14_p_rate, o16_p_rate, c12_alpha_rate, &
           triple_alpha_rate
      double precision :: local_xhe3, helium_fraction, local_xc13, &
           local_xn14
      double precision :: pp_reaction_term, he3_he4_reaction_term, &
           equilibrium_xhe3
      double precision :: cno_sum, cn_ratio_c12, cn_ratio_c13
      double precision :: triple_alpha_term, c12_alpha_term

!     SKIP BURNING CALCULATIONS IF STARTING SHELL BELOW T CUTOFF FOR REACTIONS.
      dx_dt = 0.0d0
      dy_dt = 0.0d0
      dc_dt = 0.0d0
      do_dt = 0.0d0
      if (shell_temperature(zone_begin).lt.star%ctrl%tcut(1)) return
!     COMPUTE EXPLICIT HYDROGEN AND HELIUM BURNING RATES ASSUMING EQUILIBRIUM
!     HE3 AND CN CYCLE ABUNDANCES.
!     THIS IS USED FOR AN INITIAL GUESS AT THE BURNING RATES AT THE START
!     OF THE TIME STEP, AND CAN BE SUPPLEMENTED BY A FULLY IMPLICIT
!     NON-EQUILIBRIUM CALCULATION IN THE 3RD AND 4TH LEVELS OF ITERATION.
      if (zone_begin.ne.zone_end) then
!        HOMOGENIZE CONVECTION ZONES.
!        zone_avg_abundance IS THE MASS-WEIGHTED AVERAGE ABUNDANCE FOR
!        THE CZ. total_shell_mass IS THE TOTAL MASS OF THE CZ.
!        INITIALIZE SUMS.
         total_shell_mass = 0.0d0
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = 0.0d0
         end do
         do zone_idx = zone_begin, zone_end
            total_shell_mass = total_shell_mass + shell_mass(zone_idx)
            do species_idx = 1, 11
               zone_avg_abundance(species_idx) = &
                    zone_avg_abundance(species_idx) + &
                    star%prev%xa_start(species_idx,zone_idx)*shell_mass(zone_idx)
            end do
         end do
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = &
                 zone_avg_abundance(species_idx)/total_shell_mass
         end do
      else
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = &
                 star%prev%xa_start(species_idx,zone_begin)
         end do
      end if
      if (zone_begin.eq.zone_end) then
!        PP
         pp_rate = rate_pp(zone_begin)
!        HE3,HE3
         he3_he3_rate = rate_he3_he3(zone_begin)
!        HE3,HE4
         he3_he4_rate = rate_he3_he4(zone_begin)
!        C12,P
         c12_p_rate = rate_c12_p(zone_begin)
!        C13,P
         c13_p_rate = rate_c13_p(zone_begin)
!        N14,P + N15,P
         n14_p_rate = rate_n14_p(zone_begin)
!        O16,P + O17,P.
         o16_p_rate = rate_o16_p(zone_begin)
!        C13,ALPHA -- not used (rate_c13_alpha not passed in)
!        O16,ALPHA (NOT USED)
!        C12,ALPHA
         c12_alpha_rate = rate_c12_alpha(zone_begin)
!        N14,ALPHA -- not used (rate_n14_alpha not passed in)
!        TRIPLE ALPHA
         triple_alpha_rate = rate_triple_alpha(zone_begin)
!        C12,C12 (NOT USED)
!        BRANCHING RATIO FOR N15,P :
!        f = FRACTION GOING TO C12+ALPHA, 1-f = FRACTION GOING TO O16
!        -- not used (not passed in)
      else
!        USE THE MASS-WEIGHTED AVERAGE RATES FOR THE CZ
         pp_rate = 0.0d0
         he3_he3_rate = 0.0d0
         he3_he4_rate = 0.0d0
         c12_p_rate = 0.0d0
         c13_p_rate = 0.0d0
         n14_p_rate = 0.0d0
         o16_p_rate = 0.0d0
         c12_alpha_rate = 0.0d0
         triple_alpha_rate = 0.0d0
         do zone_idx = zone_begin, zone_end
            pp_rate = pp_rate + shell_mass(zone_idx)*rate_pp(zone_idx)
            he3_he3_rate = he3_he3_rate + &
                 shell_mass(zone_idx)*rate_he3_he3(zone_idx)
            he3_he4_rate = he3_he4_rate + &
                 shell_mass(zone_idx)*rate_he3_he4(zone_idx)
            c12_p_rate = c12_p_rate + shell_mass(zone_idx)*rate_c12_p(zone_idx)
            c13_p_rate = c13_p_rate + shell_mass(zone_idx)*rate_c13_p(zone_idx)
            n14_p_rate = n14_p_rate + shell_mass(zone_idx)*rate_n14_p(zone_idx)
            o16_p_rate = o16_p_rate + shell_mass(zone_idx)*rate_o16_p(zone_idx)
            c12_alpha_rate = c12_alpha_rate + &
                 shell_mass(zone_idx)*rate_c12_alpha(zone_idx)
            triple_alpha_rate = triple_alpha_rate + &
                 shell_mass(zone_idx)*rate_triple_alpha(zone_idx)
         end do
         pp_rate = pp_rate/total_shell_mass
         he3_he3_rate = he3_he3_rate/total_shell_mass
         he3_he4_rate = he3_he4_rate/total_shell_mass
         c12_p_rate = c12_p_rate/total_shell_mass
         c13_p_rate = c13_p_rate/total_shell_mass
         n14_p_rate = n14_p_rate/total_shell_mass
         o16_p_rate = o16_p_rate/total_shell_mass
         c12_alpha_rate = c12_alpha_rate/total_shell_mass
         triple_alpha_rate = triple_alpha_rate/total_shell_mass
      end if

!     INITIAL ABUNDANCES OF SPECIES.

      hydrogen_fraction = zone_avg_abundance(1)
      local_xhe3 = zone_avg_abundance(4)
      helium_fraction = zone_avg_abundance(2)
      metal_fraction = zone_avg_abundance(3)
      equilibrium_xc12 = zone_avg_abundance(5)
      local_xc13 = zone_avg_abundance(6)
      local_xn14 = zone_avg_abundance(7)
      equilibrium_xo16 = zone_avg_abundance(9)
      if (hydrogen_fraction.gt.1.0d-10 .and. &
           shell_temperature(zone_begin).gt.star%ctrl%tcut(2)) then
!        FIND EQUILIBRIUM HELIUM-3 ABUNDANCE, USING THE QUADRATIC FORMULA.
!        THE EQUATION IS -2*R(3,3)*XHE3**2 - R(3,4)*XHE3*XHE4 + R(1,1)X**2 = 0.
         pp_reaction_term = pp_rate*hydrogen_fraction**2
         he3_he4_reaction_term = he3_he4_rate*helium_fraction
         equilibrium_xhe3 = 0.25d0*(sqrt(he3_he4_reaction_term**2 + &
              8.0d0*he3_he3_rate*pp_reaction_term) - he3_he4_reaction_term)/ &
              he3_he3_rate
!        USE THE MINIMUM OF THE CURRENT AND EQUILIBRIUM HE3 ABUNDANCE -
!        NEEDED TO AVOID SPURIOUSLY HIGH BURNING RATES IN OUTER LAYERS,
!        WHERE THE TIMESCALE FOR HELIUM-3 TO REACH EQUILIBRIUM IS LONG.
         local_xhe3 = min(local_xhe3, equilibrium_xhe3)
!        PP BURNING.
         dx_dt = -3.0d0*pp_reaction_term + 2.0d0*he3_he3_rate*local_xhe3**2 &
              - he3_he4_reaction_term*local_xhe3
! (Restructured 2026: the CN block below runs only on the hydrogen-
! burning branch; the old `goto 100` skipped it from the else.)
      if (shell_temperature(zone_begin).gt.star%ctrl%tcut(3)) then
!        FIND EQUILIBRIUM C12,C13,N14 ABUNDANCES TREATING CN PROCESSING AS
!        A CLOSED LOOP.
         cno_sum = equilibrium_xc12/1.2d1 + local_xc13/1.3d1 + &
              local_xn14/1.4d1
!        FIND EXPLICIT RATES FOR P+C12,P+C13,P+N14(==R121,R131,R141)
         cn_ratio_c12 = n14_p_rate/c12_p_rate
         cn_ratio_c13 = n14_p_rate/c13_p_rate
!        IN EQUILIBRIUM WE HAVE
!        XC12 = XN14{R(1,14)/R(1,12)}
!        XC13 = XN14{R(1,14)/R(1,13)}
!        AND FROM CONSERVATION OF SPECIES XC12/12 +XC13/13 + XN14/14 IS
!        UNCHANGED.
!        EXPRESSING C12 AND C13 AS FUNCTIONS OF N14, SOLVE FOR THE NEW N14
!        EQUILIBRIUM ABUNDANCE.
         local_xn14 = cno_sum/(cn_ratio_c12/1.2d1 + cn_ratio_c13/1.3d1 + &
              1.0d0/1.4d1)
         equilibrium_xc12 = local_xn14*cn_ratio_c12
         local_xc13 = local_xn14*cn_ratio_c13
         dx_dt = dx_dt - c12_p_rate*hydrogen_fraction*equilibrium_xc12 - &
              c13_p_rate*hydrogen_fraction*local_xc13 - &
              2.0d0*n14_p_rate*hydrogen_fraction*local_xn14 - &
              2.0d0*o16_p_rate*hydrogen_fraction*equilibrium_xo16
      end if
      else
         dx_dt = 0.0d0
      end if
!     HELIUM BURNING REACTIONS.
      if (helium_fraction.gt.1.0d-10 .and. &
           shell_temperature(zone_begin).gt.star%ctrl%tcut(4)) then
         triple_alpha_term = triple_alpha_rate*helium_fraction**3
         c12_alpha_term = c12_alpha_rate*helium_fraction*equilibrium_xc12
         dy_dt = -4.0d0*(c12_alpha_term + 3.0d0*triple_alpha_term)
         dc_dt = 1.2d1*(triple_alpha_term - c12_alpha_term)
         do_dt = 1.6d1*c12_alpha_term
      else
         dy_dt = 0.0d0
         dc_dt = 0.0d0
         do_dt = 0.0d0
      end if

      return
end subroutine eqburn


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

      use star_info_lib, only: star, json
      use const_lib
      implicit none

      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: timestep
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
         hydrogen_fraction = star%prev%xa_start(1,zone_begin)
         deuterium_fraction = star%prev%xa_start(12,zone_begin)
         helium3_fraction = star%prev%xa_start(4,zone_begin)
         rate_start = star%light_burn%deuterium_burning_rate_start(zone_begin)
         rate_end = star%light_burn%deuterium_burning_rate(zone_begin)
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
                 shell_mass(zone_idx)*star%light_burn%deuterium_burning_rate_start(zone_idx)
            rate_end_sum = rate_end_sum + &
                 shell_mass(zone_idx)*star%light_burn%deuterium_burning_rate(zone_idx)
            hydrogen_fraction = hydrogen_fraction + &
                 star%prev%xa_start(1,zone_idx)*shell_mass(zone_idx)
            deuterium_fraction = deuterium_fraction + &
                 star%prev%xa_start(12,zone_idx)*shell_mass(zone_idx)
            helium3_fraction = helium3_fraction + &
                 star%prev%xa_start(4,zone_idx)*shell_mass(zone_idx)
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
              star%ctrl%accreted_composition(12)*star%light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + star%light_burn%accreted_mass_fraction)
      else
         deuterium_fraction_test = deuterium_fraction
      end if
      if (deuterium_fraction_test.lt.1.0d-11) then
         do zone_idx = zone_begin, zone_end
            composition(12,zone_idx) = 0.0d0
         end do
         continue
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
!        AVERAGE. 1/2 OF THE TIMESTEP.  star%light_burn%accreted_mass_fraction IS
!        DEFINED AS DMDT*DT/ORIGINAL CZ MASS.
!        BURN BOTH THE ACCRETED D AND THE ORIGINAL D
!        SEPARATELY AND FIND THE NEW MASS-WEIGHTED
!        AVERAGE D.  NOTE THAT ACCRETION OF ELEMENTS
!        1-11 HAS ALREADY BEEN TREATED.
         rate_accum = rate_start
         accreted_deuterium_fraction = star%ctrl%accreted_composition(12)
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
              star%ctrl%accreted_composition(12)
         deuterium_fraction_new = (deuterium_fraction_burned*total_shell_mass &
              + accreted_deuterium_burned*star%light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + star%light_burn%accreted_mass_fraction)
         hydrogen_fraction_new = hydrogen_fraction + 0.5d0* &
              (deuterium_change_original*total_shell_mass + &
              deuterium_change_accreted*star%light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + star%light_burn%accreted_mass_fraction)
         helium3_fraction_new = helium3_fraction - 1.5d0* &
              (deuterium_change_original*total_shell_mass + &
              deuterium_change_accreted*star%light_burn%accreted_mass_fraction)/ &
              (total_shell_mass + star%light_burn%accreted_mass_fraction)
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
! Near-duplicate of dburn.f90 (compute the abundance changes resulting
! from deuterium burning by sub-stepping the burning rate across the
! timestep) used from a different call site: here the start/end
! burning rates are passed in explicitly as dummy arguments
! (deuterium_rate_start/deuterium_rate_end) rather than read from
! common/deuter/, and the accreted-mass fraction used for this call is
! only a fraction (step_fraction) of the full-timestep value in
! common/masschg/. See dburn.f90 for the fuller commentary on the
! sub-stepping algorithm.
subroutine dburnm(zone_begin, zone_end, num_zones, shell_mass, &
     composition, timestep, deuterium_rate_end, deuterium_rate_start, &
     step_fraction)
      use star_info_lib, only: star, json
      use const_lib
      implicit none

      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: timestep
      double precision, intent(in) :: deuterium_rate_end(json), &
           deuterium_rate_start(json)
      double precision, intent(in) :: step_fraction
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
! the units used elsewhere for the burning rates.
      timestep_gyr = timestep*1.0d-9/seconds_per_year
      if(zone_begin.eq.zone_end)then
         hydrogen_fraction = star%prev%xa_start(1,zone_begin)
         deuterium_fraction = star%prev%xa_start(12,zone_begin)
         helium3_fraction = star%prev%xa_start(4,zone_begin)
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
                 star%prev%xa_start(1,zone_idx)*shell_mass(zone_idx)
            deuterium_fraction = deuterium_fraction + &
                 star%prev%xa_start(12,zone_idx)*shell_mass(zone_idx)
            helium3_fraction = helium3_fraction + &
                 star%prev%xa_start(4,zone_idx)*shell_mass(zone_idx)
         end do
         rate_start = rate_start_sum/total_shell_mass
         rate_end = rate_end_sum/total_shell_mass
         hydrogen_fraction = hydrogen_fraction/total_shell_mass
         deuterium_fraction = deuterium_fraction/total_shell_mass
         helium3_fraction = helium3_fraction/total_shell_mass
      endif
      if(deuterium_fraction.lt.1.0d-14)then
         do zone_idx = zone_begin,zone_end
            composition(12,zone_idx) = 0.0d0
         end do
         continue
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
! AVERAGE. 1/2 OF THE TIMESTEP.  star%light_burn%accreted_mass_fraction IS
! DEFINED AS DMDT*DT/ORIGINAL CZ MASS.
! BURN BOTH THE ACCRETED D AND THE ORIGINAL D
! SEPARATELY AND FIND THE NEW MASS-WEIGHTED
! AVERAGE D.  NOTE THAT ACCRETION OF ELEMENTS
! 1-11 HAS ALREADY BEEN TREATED.
         rate_accum = rate_start
         accreted_deuterium_fraction = star%ctrl%accreted_composition(12)
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
         accreted_mass_fraction_substep = step_fraction*star%light_burn%accreted_mass_fraction
         deuterium_change_original = deuterium_fraction_burned - &
              deuterium_fraction
         deuterium_change_accreted = accreted_deuterium_burned - &
              star%ctrl%accreted_composition(12)
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
         composition(12,zone_idx) = deuterium_fraction_new
      end do
      return
end subroutine dburnm


!----------------------------------------------------------------------
! deutrate
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original deutrate.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Compute the rate of nonequilibrium deuterium burning (excluding the
! abundance factor) at shell i, storing it in common/deuter/ for later
! use by dburn/dburnm. If the shell lies within (or below) the surface
! convection zone, the rate is capped so that deuterium burning cannot
! proceed faster than the local convective overturn timescale.
subroutine deutrate(dl,tl,x,i,itlvl)
      use star_info_lib, only: star, json
      use const_lib
      implicit none

      double precision, intent(in) :: dl, tl, x
      integer, intent(in) :: i, itlvl





      double precision :: c21
      data c21/5.240358E-8/
! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
!  POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
!  OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
!  CGS UNITS.
      double precision :: rho, t9, t9p13, t9p23, t9m13, t9m23, t9m1
      double precision :: z, tfacdeut, rdeut, rdeutmax, rdeut2

      rho=exp(ln10*dl)
      t9 = exp(ln10*(tl - 9.0d0))
      t9p13 = t9**cc13
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
      z = -3.72d0*t9m13
      tfacdeut = 1.0d0+0.112d0*t9p13+3.38d0*t9p23+2.65d0*t9
! FACTOR OF 3.0115D23 REFLECTS AVAGADROS NUMBER DIVIDED BY THE
! MASS OF THE DEUTERON IN AMU
      rdeut = rho*2.240d3*t9p23*exp(z)*tfacdeut*3.0115d23
! NOW LIMIT DEUTERIUM BURNING IN A SURFACE CZ TO BE ON A TIME SCALE
! NO SHORTER THAN THE CONVECTIVE OVERTURN TIMESCALE.
      if(i.ge.star%light_burn%jcz .and. star%turnover%convective_turnover_timescale.gt.1.0d0)then
         rdeutmax = 3.0115d23/star%turnover%convective_turnover_timescale
         rdeut2 = rdeut*x
         if(i.eq.star%light_burn%jcz)then
! JVS 0712 Commented out write command
!            WRITE(*,911)RDEUT2,RDEUTMAX,TAUCZ
!  911        FORMAT(1P3E15.8)
         endif
         if(rdeut2.gt.rdeutmax)then
! JVS 0712 Commented out write command
!            WRITE(*,*)RDEUT2,RDEUTMAX
            if(x.gt.1.0d-6)then
               rdeut = rdeutmax/x
            endif
         endif
      endif
      star%light_burn%deuterium_burning_rate(i) = x*rdeut*c21
      if(itlvl.eq.1)then
         star%light_burn%deuterium_burning_rate_start(i) = star%light_burn%deuterium_burning_rate(i)
      endif
      return
end subroutine deutrate


!----------------------------------------------------------------------
! engeb
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original engeb.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 3/92 DBG Added new neutrino loss calculation routines
! NOVEMBER 6, 1990 (JNB)
! THE FOLLOWING SUBROUTINE CALCULATES THE NUCLEAR REACTION RATES WITH
!   SPECIAL ATTENTION TO DETAIL REQUIRED FOR CALCULATING SOLAR NEUTRINO
!   FLUXES. THE NEUTRINO FLUXES ARE EVALUATED ALSO IN THIS SUBROUTINE.
! ALL PHYSICAL QUANTITIES ARE CGS. THE REACTION RATES ARE IN
!  GM^{-1}S^{-1}.
!  MHP 10/97
!
!   On 10/13/97, JNB converted the nuclear masses from neutral nuclear
!   masses to bare nuclear masses by subtracting Z(I)*(m_e)*c^2 from
!   the neutral nuclear masses. This caused changes in a number of
!   places: in ANUC(I), in Q1(I)-Q7(I), in the calculation of the Be7electron
!   Be7proton rates, and in the calculation of the N15p branching ratio.
!
!   JNB made some purely cosmetic changes on 1/20/96. Revised all of the
!   input Q's on 9/23-25/97 to agree with submitted version of Solar
!   Fusion Workshop paper. The SStandard are fixed to agree with
!   the Workshop paper. JNB recalculated all of the EG(I) to determine
!   the best values for the energy generation for all the reactions,
!   taking account of my improved calculations of neutrino energy loss.
!   The calculations are documented in Vol. 19, 132-141, 1997 of my notes.
!   The neutral atom mass differences are taken from Table of Isotopes,
!   8th Ed, 1996 and the neutrino energy losses from Bahcall, Gallium
!   solar neutrino experiments, Phys Rev C, in press, 1997.
!
!  All numbers in this subroutine have been calculated by John Bahcall
!  so that they agree with the modern numbers in Neutrino
!  Astrophysics (1989) or much more recent results, as indicated in
!  comments lines.
!
!  PREVIOUS VERSIONS OF THE YALE CODE
!  USED THE NUCLEAR REACTION RATES PRIMARILY FROM FOWLER, CAUGHLAN, AND
!  ZIMMERMAN IN ANN. REV. ASTR. AND ASTROPHYS. (1975) OR, FOR THE BE7 AND
!  B8 RATES, THE VALUES FROM BAHCALL AND SEARS (1972), ALSO. ANN. REV. ASTR.
!  AND ASTROPHYS. SOME VALUES WERE FROM HARRIS ET AL (1983, ANN. REV. ASTR.
!  IN A FEW CASES, NUMERICAL ERRORS (NOT DUE TO REVISIONS IN THE NUCLEAR DATA)
!  HAVE BEEN CORRECTED.  IN A SMALL NUMBER OF OTHER PLACES, THE PREVIOUS
!  PROGRAM CONTAINED PHYSICALLY INCORRECT STATEMENTS, WHICH HAVE BEEN CORRECTED
!  IN THE PRESENT VERSION. IN THE COURSE OF THE CHECKING AND REVISIONS, I
!  HAVE ADDED MANY COMMENT STATEMENTS IN ORDER TO MAKE THE PROGRAM MORE
!  TRANSPARENT AND SIMPLER TO REVISE.  SOME STATEMENTS HAVE BEEN ADDED
!  TO SIMPLIFY THE PROCESS OF REVISING THE SUBROUTINE AS IMPROVED DATA
!  BECOME AVAILABLE.  FOR EXAMPLE, I HAVE INSERTED EXPLICITLY
!  THE STANDARD CROSS SECTION FACTORS IN A WAY THAT IS EASY TO REVISE;
!  SSTANDARD(I) = 1.0 WHEN THE ITH CROSS SECTION HAS THE STANDARD VALUE
!  MEASUREMENTS USED IN NEUTRINO ASTROPHYSICS.
!  MHP 10/97
!
!  The PREVIOUS standard values of SStandard(I) used were those that
!   were given in Table 1 of Bahcall and Pinsonneault (1992).  This work
!   was  published in the Reviews of Modern Physics, 64,885, 1992. The
!   current version of SStandard(I) refers to Bahcall and
!   Pinsonneault(1995) , Rev. Mod. Phys. 67, 781 (1995), Table 1, if
!   changes have occured.  Otherwise, they are the same as in the
!   Bahcall-Pinsonnealt (1992) Rev. Mod. Phys. article.
!
!  Weakscreening is a parameter set in this subroutine. To obtain
!  the Graboske et al. and Salpeter standard results, use:
!  weakscreening = 0.03.  To investigate the effect of always using
!  weak screening, use a value for weakscreening greater than unity,
!  e. g., 30. The transition region between weak and strong screening
!  is defined by Graboske et al. as $\Lambda_{12} = 0.1$, see page
!  ApJ 181, 465 (1973) [5/15/97].
!
!  The value of SStandard(17) for hep is taken from Carlson et al (1991)
!  Phys. Rev. C 44, 619. It corresponds to an S sub 0 =
!  1.3 E-20 keV-barns, a factor of 0.1625  smaller than indicated by the
!  older measurements and calculations used in Neutrino Astrophysics.
!
! THE NUCLEAR ENERGY RELEASE TO THE STAR FROM EACH REACTION
!  IS TAKEN FROM BAHCALL AND ULRICH, RMP 60, 297 (1988) AND TAKES ACCOUNT
!  ACCURATELY OF NEUTRINO ENERGY LOSS.
!
! I HAVE ADDED A NEW SECTION AT THE END OF THE SUBROUTINE THAT
!  CALCULATES THE SOLAR NEUTRINO FLUXES AT THE EARTH.  THESE FLUXES
!  ARE IN THE UNITS OF CM^-2 SEC^-2 PER GM. TO GET THE FLUX
!  FROM A SHELL, MULTIPLY BY THE MASS OF THE SHELL IN UNITS OF
!  GRAMS.  THE FLUXES ARE IN A COMMON BLOCK, FLUXES.  I ALSO CALCULATE
!  THE FICTIONAL NEUTRINO FLUXES ASSOCIATED WITH THE HE3 + HE3 AND WITH
!  THE HE3 + HE4 RECTIONS; THESE FICTIONAL FLUXES ARE USEFUL DIAGNOSTICS
!  OF THE SOLAR MODEL.
!
! *****************************************
! CHANGING NUCLEAR REACTION CROSS SECTIONS.
! *****************************************
! NUCLEAR CROSS SECTIONS CAN BE CHANGED SIMPLY BY INSERTING NEW NUMBERS
!  FOR THE DATA VALUES OF SSTANDARD(I), WHICH ARE THE RATIOS OF THE
!  DESIRED CROSS SECTION FACTORS TO THE VALUES GIVEN IN TABLES 3.2 AND
!  3.4 OF NEUTRINO ASTROPHYSICS (1989). IF THE VALUE GIVEN IN NEUTRINO
!  ASTROPHYSICS IS USED, THEN SSTANDARD(I) = 1.0 .  TO INCREASE THE CROSS
!  SECTION FOR REACTION K BY A FACTOR OF TWO COMPARED TO THE STANDARD
!  VALUE, SET SSTANDARD(K) = 2.0 .  TO DETERMINE WHICH VALUE OF I GOES
!  WITH WHICH REACTION, SEE THE SECTION JUST BELOW.
!  THE ENERGY DERIVATIVES ENTER IN A FORM IN WHICH THEY ARE DIVIDED BY
!  THE ABSOLUTE VALUES OF THE CROSS SECTIONS AT ZERO ENERGY.  THUS IF
!  THE SHAPE OF THE CROSS SECTION EXTRAPOLATION IS UNCHANGED AND ONLY
!  THE INTERCEPT OF S(E) AT ZERO ENERGY IS CHANGED, THEN NO CORRECTION
!  NEED BE MADE FOR THE DERIVATIVES.  THEY ARE AUTOMATICALLY SCALED
!  CORRECTLY.  THE EXACT WAY THAT THE DERIVATIVES ENTER THE
!  CALCULATIONS IS DESCRIBED IN THE SECTION LABELED ``DEFINING THE
!  Q(I)'' THAT IS PRESENTED BELOW.
! ************************************
! IDENTIFYING THE REACTIONS.
! ************************************
!  THE VALUE OF J DENOTES WHICH OF THE REACTIONS THE COEFFICIENTS
!  REFER TO:
!  J = 1, PP; J = 2, HE3+HE3; J = 3, HE3+ HE4; J =4, P + C12;  J = 5, P+C13;
!  J = 6. P + N14; J = 7, P + O16.
!  REACTIONS J = 8, 13 ARE NOT RELEVANT FOR THE SOLAR INTERIOR; THEY ARE
!   HOLDOVERS FROM THE EARLIER YALE CODE.
!  REACTION 14 IS PEP; REACTION 15 IS BE7 ELECTRON CAPTURE; REACTION 16 IS
!   BE7 PROTON CAPTURE; REACTION 17 IS THE HEP REACTION.
!  do not change sstandard(14) unless you want to change the ratio of
!  pep to pp.
!   REACTIONS 14-17 WERE NOT EXPLICITLY INCLUDED IN THE YALE
!    PREVIOUS VERSION OF THE CODE, BUT THEY ARE IMPORTANT FOR
!    SOLAR NEUTRINO CALCULATIONS.
!   THE BRANCHING OF THE N15 + P REACTIONS IS TREATED IN A SERIES OF
!    SEPARATE STATEMENTS FOLLOWING THE CALCULATION OF THE BE7 + P
!    REACTION. SEE THE DEFINITIONS OF F3 AND F4.  IF THE CROSS-SECTION
!    FACTORS OF THE N15 + P REACTIONS ARE REVISED, THEN THE NUMERICAL
!    COEFFICIENTS MUST BE CHANGED IN THE DEFINITION OF C12ALPHA AND
!    O16GAMMA.
! FOR Q1(I), ...,Q(5(I), I = 8 CORRESPONDS TO THE BE7 +  P REACTION.
!  THIS ASSIGNMENT FOR I = 8 IS ONLY VALID FOR THE LISTED Q'S AND NOT
!  FOR OTHER ARRAYS IN THE PROGRAM.
!  IU IS THE SHELL NUMBER.
!
! Dummy-argument renaming used throughout this modernized version:
!   EPP1 -> pp_chain_energy_gen        EPP2 -> he3he4_be7_electron_energy_gen
!   EPP3 -> he3he4_be7_proton_energy_gen   ECN -> cno_cycle_energy_gen
!   E3AL -> triple_alpha_energy_gen    PEP -> dlnepsilon_dlnrho
!   PET -> dlnepsilon_dlnt             SUM1 -> total_energy_gen_rate
!   DL -> log_density                  TL -> log_temperature
!   X -> hydrogen_fraction             Y -> helium_fraction
!   XHE3 -> he3_fraction               XC12/XC13 -> c12_fraction/c13_fraction
!   XN14 -> n14_fraction               XO16/XO18 -> o16_fraction/o18_fraction
!   XH2 -> deuterium_fraction          IU -> shell_index
!   HR1..HR13 -> reaction_rate_1..13 (yr^-1, amu^-1, output for kemcom)
!   HF1 -> n15_alpha_branch_fraction   HF2 -> be7_electron_capture_fraction
subroutine engeb(pp_chain_energy_gen, he3he4_be7_electron_energy_gen, &
     he3he4_be7_proton_energy_gen, cno_cycle_energy_gen, &
     triple_alpha_energy_gen, dlnepsilon_dlnrho, dlnepsilon_dlnt, &
     total_energy_gen_rate, log_density, &
     log_temperature, hydrogen_fraction, helium_fraction, he3_fraction, &
     c12_fraction, c13_fraction, n14_fraction, o16_fraction, &
     o18_fraction, deuterium_fraction, shell_index, reaction_rate_1, &
     reaction_rate_2, reaction_rate_3, reaction_rate_4, reaction_rate_5, &
     reaction_rate_6, reaction_rate_7, reaction_rate_8, reaction_rate_9, &
     reaction_rate_10, reaction_rate_11, reaction_rate_12, &
     reaction_rate_13, n15_alpha_branch_fraction, &
     be7_electron_capture_fraction)

      use star_info_lib, only: star, i_nu_b8, i_nu_be7, i_nu_f17, i_nu_hep, i_nu_n13, i_nu_o15, i_nu_pep, i_nu_pp, json
      use luout_lib
      use const_lib
      implicit none

      double precision, intent(out) :: pp_chain_energy_gen, &
           he3he4_be7_electron_energy_gen, he3he4_be7_proton_energy_gen, &
           cno_cycle_energy_gen, triple_alpha_energy_gen, &
           dlnepsilon_dlnrho, dlnepsilon_dlnt
! total_energy_gen_rate (originally SUM1) is intent(inout), not
! intent(out): in the log_temperature.le.tcut(1) early-return branch
! below (preserved verbatim from the original), it is never assigned,
! so the caller's incoming value is left untouched on that path -- a
! real property of the original F77 code, not a bug introduced here.
      double precision, intent(inout) :: total_energy_gen_rate
      double precision, intent(in) :: log_density, log_temperature, &
           hydrogen_fraction, helium_fraction, he3_fraction, c12_fraction, &
           c13_fraction, n14_fraction, o16_fraction, o18_fraction, &
           deuterium_fraction
      integer, intent(in) :: shell_index
      double precision, intent(out) :: reaction_rate_1(json), &
           reaction_rate_2(json), reaction_rate_3(json), &
           reaction_rate_4(json), reaction_rate_5(json), &
           reaction_rate_6(json), reaction_rate_7(json), &
           reaction_rate_8(json), reaction_rate_9(json), &
           reaction_rate_10(json), reaction_rate_11(json), &
           reaction_rate_12(json), reaction_rate_13(json), &
           n15_alpha_branch_fraction(json), &
           be7_electron_capture_fraction(json)










! 9/06 GN --- New neutrino loss common block
! KC 2025-05-30 reordered common block elements
!       COMMON/NULOSS/LNULOS1,DSNUDT,DSNUDD
! former common/nuloss/: use_itoh_neutrino_loss (switch selecting the
! Itoh 1996 neutrino-loss routines below) is real shared configuration
! -- now use-associated from const_lib. neutrino_dlnq_dlnt/
! neutrino_dlnq_dlnd (the log-derivatives of the resulting loss rate
! w.r.t. T/rho) are set and consumed entirely within this file --
! core/parmin.f90, which declared the same common block, never
! actually touches them -- so they're genuinely local, not shared
! state, and become plain locals here rather than moving to a module.
      double precision :: neutrino_dlnq_dlnt, neutrino_dlnq_dlnd





      double precision :: mass_fraction(13), reaction_rate(13), &
           dlnrate_dlnrho(13), dlnrate_dlnt(13), screening_factor(13), &
           dscreen_dlnrho(13), dscreen_dlnt(13), reaction_energy_gen(13), &
           charge_product(13), z53(13), z43(13), z23(13), &
           z86(13), q1(8), q2(8), q3(8), q4(8), q5(8), q6(7), q7(7), &
           q8(7), eg(50)
      double precision :: atomic_mass_amu(13), atomic_number(13)
      double precision :: v1(7), v2(7), v3(7)
      double precision :: years_per_sec_over_amu

      integer :: num_isotopes, nrxns

! ***MHP 5/91 STATEMENTS ADDED FOR EVOLVED STAR NEUTRINO LOSSES.
! NEUTRINO COEFFICIENTS

      data v1,v2,v3/ &
           6.002d19,2.084d20,1.872d21,9.383d-1,-4.141d-1,5.829d-2,5.5924, &
           4.886d10,7.580d10,6.023d10,6.290d-3,7.483d-3,3.061d-4,1.5654, &
           2.320d-7,8.449d-8,1.787d-8,2.581d-2,1.734d-2,6.990d-4,0.56457/
! ***************************
! ANUC ARE ATOMIC MASS UNITS.
! ***************************
!  On 10/13/97, JNB converted ANUC(I) from neutral atomic masses to
!  bare nuclear masses.
!
!  The scale is the mass of C12 divided by 12 or 931.49432 MeV,
!  which is 1.6605402 times 10^{-24} gm. The values are obtained
!  by dividing the mass excess (expressed in MeV) by 931.49432 MeV and
!  adding to this the atomic mass number, A.  The value for Be7, which
!  is used implicitly in this subroutine was, 7.016930, until 13/13/97.
!  We now use the bare nuclear mass of 7.014735 .
!
      data atomic_mass_amu/1.008665,1.007276,2.013553,3.015501,3.014933, &
           4.001506,11.996709,13.000064,13.999233,15.990526,17.994772, &
           19.986954,23.978458/, &
           atomic_number/0.,1.,1.,1.,2.,2.,6.,6.,7.,8.,8.,10.,12./, &
           num_isotopes/13/

! THE ISOTOPES ARE NEUTRON, H1, D, H3, HE3, HE4, C12, C13, N14, O16, O18,
!  NE20, MG24, RESPECTIVELY. ALL OF THESE NUMBERS WERE CHECKED.
! NELEM IS THE NUMBER OF ISOTOPES INCLUDED.
! **************************************************************************
! THE QUANTITIES Q1(J), Q2(J), ...,Q5(J) ARE THE TERMS IN EQUATION 3.14 OF
!  NEUTRINO ASTROPHYSICS AND IN EQUATION 53 OF FOWLER, CAUGHLAN, AND
!  ZIMMERMAN (1967) EQ. 53, IN BOTH CASES MULTIPLIED BY T SUB 9 ^(-2/3).
!  THE REACTIONS CORRESPONDING TO EACH J ARE LISTED ABOVE, UNDER:
!   IDENTIFYING THE RECTIONS.
!  FOR THIS SET OF PARAMETERS, AND ONLY FOR THIS SET OF PARAMETERS,
!  J = 8 CORRESPONDS TO THE BE7 + P REACTION.
! **************************************************************************
! GENERAL EXPRESSION FOR Q'S:
! T_9^(-2/3)[S_EFF/S(0)] =
!   [T_9^(-2/3) + Q1(I)T_9^(-1/3) + Q2(I) + Q3(I)T_9^(1/3) +
!    Q4(I)T_9^(2/3) + Q(5)T_9 ]
!
! BY COMPARISON WITH EQUATION 3.14 OF NEUTRINO ASTROPHYSICS, WE SEE THAT
!
!
! EACH OF THE Q'S IS INDEPENDENT OF TEMPERATURE (T), AS CAN BE SEEN FROM
!  EQUATIONS 3.10 AND 3.11 .
!
! MHP 10/97
! All of the values of the Q1, ...., Q5 have been recalculated, using
!  where needed nuclear cross sections given in Tables 3.2 and 3.4 of
!  Neutrino Astrophysics.  They have been updated on 9/16/97 to reflect
!  the derivatives and cross section factors in Adelberger et al. 1998,
!  the paper on the Solar Fusion Workshop.
!
! ******************************************************************
! Q6 IS THE COEFFICIENT OF THE TEMPERATURE TERM IN THE DEFINITION OF
!  TAU, EQUATION 3.10 OF NEUTRINO ASTROPHYSICS.
!  TAU = Q6*(T SUB 9 TO THE (-1/3) POWER ).
! ******************************************************************
! SLIGHT CHANGES HAVE BEEN MADE IN THE PREVIOUS VALUES OF Q6 TO MAKE
!  THE DATA MORE ACCURATE.
! NOTE THAT Q6 IS NEGATIVE.
! ********************************************************************
! Q7 IS THE CONSTANT IN FRONT OF THE REACTION RATE. THE LOGARITHMS ARE
!  RELATED TO S FACTORS BY A LOGARITHM AND VARIOUS NUMERICAL FACTORS.
! ********************************************************************
! THE GENERAL RELATION IS:
!  Q7 = 70.62860 + ((LN(Z0*Z1/A))/3) -LN(A0*A1) - LN(S SUB 0)
!       -LN(1 + DELTA_01)
!  HERE S SUB 0 IS THE CROSS SECTION FACTOR IN UNITS OF KEV-BARNS.
!  THE NUMERICAL VALUES USED HERE ARE TAKEN FROM TABLES 3.4 AND 3.2
!   OF NUCLEAR ASTROPHYSICS.
!  THE QUANTITY DELTA_01 IS NON-ZERO (EQUAL TO UNITY) ONLY WHEN THE
!   TWO REACTING NUCLEI, 0 AND 1, ARE IDENTICAL.
! Q8 REFLECTS A TERM IN THE EXPONENTIAL THAT OCCURS IN THE RATES, THE
!  TERM BEING PROPORTIONAL TO E^( CONSTANT*T_9^2). SEE HARRIS, ET AL.
!  1983, ANNUAL REV. ASTRON. ASTROPHYS.21, 165 (1983) FOR THE MEANING
!  OF THIS OBSCURE TERM.
! THE VALUES OF Q1(I),...Q7(I) GIVEN BELOW WERE OBTAINED USING THE DATA
!  IN NEUTRINO ASTROPHYSICS WITH THE AID OF AUXILARY COMPUTER CODES
!  THAT GENERATED ACCURATE EVALUATIONS. THE NUMBERS GIVEN HERE ARE
!  IN MANY CASES COMPLETELY DIFFERENT FROM THE VALUES IN THE ORIGINAL
!  YALE CODE.
! NRXNS IS THE NUMBER OF REACTIONS BEING TRACKED.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! 10/13/97. Changed Q(I) so that are now calculated for the bare nuclear
! masses. The calculates were made cues.f . Last date on the calculations
! in my notes is 10/10/97.  Note Q6 is -tau(T_9 = 1); in cues.f, I calculate
! tau(T_6 = 1). The connection is:  Q6 = 0.1*tau(cues.f).
!
!
! BP00 values
!      DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
!     1Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, -0.2095/,
!     2Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/,
!     3Q4/0.,0.,0.,4.79,7.456,0.0,.69703,.16762/,
!     4Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,.12114/,
!     5Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/,
!     6Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
!     7Q8/0.,0.,0.,0.0,0.0,0.0,0./,
!     8NRXNS/13/
! Revised values 7/21/03
!      DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
!     1Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, 0.0/,
!     2Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,0.0/,
!     3Q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.0/,
!     4Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.0/,
!     5Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/,
!     6Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
!     7Q8/0.,0.,0.,0.0,0.0,0.0,0./,
!     8NRXNS/13/
! MHP 9/14 RESTORED BE7+P DERIVATIVES THAT WERE ZEROED OUT
      data q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/, &
           q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224,-0.2095/, &
           q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/, &
           q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.16762/, &
           q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.12114/, &
           q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/, &
           q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/, &
           q8/0.,0.,0.,0.0,0.0,0.0,0./, &
           nrxns/13/
! ***NOTE THAT SSTANDARD IS AN INPUT PARAMETER SET IN THE NAMELIST;
! PREVIOUS PUBLISHED SETS OF SSTANDARD ARE INDICATED BELOW.
! Changed slightly 3He-3He on 9/25/97 to take account of the S'.
!
!****************FD Feb 09 ************************************
! Original data from Neutrino Astrophysics 1989 Bahcall
! Table 3.2 and 3.4
! Table 3.2 pp
!         Q(MeV)  S_0(KeV barns)
! H-p     1.44    4.07E-022
! He3-He3 12.86   5.15E+003
! He3-He4 1.59    0.54
! 7Be-p   0.14    0.02
! He-p    19.8    8.00E-020
!
!
! Table 3.4 cno   Q(MeV)  S_0(MeV barns)
!
! C12-p   1.94    1.45E-003
! C13-p   7.55    5.50E-003
! N14-p   7.3     3.32E-003
! N15-p O16       12.13   6.40E-002
! N15-p C12       4.97    78
! O16-p   0.6     9.40E-003
!
!****************************************************************************
! FD-MP Fev 2009 values changed to included the results presented at Gran Sasso
! by A. Formicola
! Cross sections are:
!  s11      pp      = 3.92 +- 0.08 D-25 MeVb      Physun talk Formicola
!  s13      Hep     = 8.6  +- 1.3  D-20 KeVb
!  s33      He3-He3 = 5.32 +- 0.08 MeVb
!  s34      He3_He4 = 0.568 +- 0.014 KeVb         mean LUNA 2007 - Brown 2007 - Singh 2004
!  s17      Be7+p   = 22.1 +-0.6 +-0.6 eVb        A.Junghans PRC70(2004)045501
!  s114     N14-p   = 1.57 +- 0.13 KeVb           Marta et al. Luna coll. PRC 78 (2008) 022802)
!
! This cross section factor gives the following SSTANDARD
!
!  SSTANDARD=0.9631,1.0330,1.0519,0.9241,1.3818,0.47290,1.0,
! 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.84770,1.0750
!**********************************************************************
! BP03 values, changed 7/21/03
!     SStandard/0.9681,1.0485,0.9815,0.9241,1.3818,1.0542,1.0,
!    $          1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.8807,1.075/
!  Previously (6/16/97) used S at Gamow Peak. Agrees with Workshop paper.
!
!     SStandard/0.9828,1.0485,0.9815,0.9241,1.3818,1.0542,1.0,
!    $          1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.7819,0.2875/
!  This is the revised version created from the Seattle Workshop paper
!  on Solar Fusion Reactions. Corrections inserted 9/16/1997. See pg.
!  132, Volume 19, of notes.
!     9  SStandard/0.9828,1.0291,0.9815,0.9241,1.3818,1.0542,1.0,
!     $ 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.018,0.7819,0.2875/
! VALUES BASED ON PARKER REVIEW ARTICLE.
!    9  SSTANDARD/1.0049,0.9709,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,1.0,1.913/
! BP 92 VALUES
!     9  SSTANDARD/0.9828,0.9709,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!     $  1.0,1.0,1.0,1.0,0.9218,0.1625/
! BP 95 VALUES
!     9  SSTANDARD/0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0,1.0,1.0,
!     $  1.0,1.0,1.0,1.0,1.0,0.92088,0.1625/
! LATEST VALUES, INCLUDING NEW PP RATE
!    9  SSTANDARD/0.9558,0.9709,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,0.9218,0.1625/
! Dar & Shaviv Values
!    9  SSTANDARD/0.9828,0.9709,0.8333,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,0.6996,0.1625/
! Schramm & Shi Values
!    9  SSTANDARD/0.9828,1.0874,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,0.8313,0.1625/
! *********************************************************************
! THE VALUES OF SSTANDARD(I) ARE TO BE CHANGED FROM UNITY IF THE CROSS
!  SECTION FACTORS ARE NOT THE ONES GIVEN IN NEUTRINO ASTROPHYSICS,
!  TABLE 3.2 AND TABLE 3.4 .
! THE CURRENT VALUES ARE CHANGED TO BE THE PREFERRED VALUES LISTED IN
!  THE LAST COLUMN OF TABLE 1 OF BAHCALL AND PINSONNEAULT (1991),
!   6/14/91.
! *********************************************************************
! ZPRD IS USED IN SCREENING CALCULATIONS. IT IS THE PRODUCT OF THE
!  CHARGES OF THE INTERACTING IONS. ZPRD WAS CHECKED. Z86 IS USED
!  IN CALCULATING INTERMEDIATE SCREENING AND IS DEFINED BY GRABOSKE ET
!  AL, AP. J. 181, PAGE 465 (1973), IN TABLE 4.  Z86 WAS CHECKED AND
!  SOME NUMERICAL VALUES WERE MADE SLIGHTLY MORE ACCURATE.  Z53, Z43,
!  AND Z23 ARE ALSO DEFINED IN TABLE 4 (SEE ABOVE).  SINCE THEY ARE
!  ONLY USED IN STRONG SCREENING, THE VALUES OF Z53, Z43, AND Z23
!  WERE NOT CHECKED.
      data charge_product/1.,4.,4.,6.,6.,7.,8.,12.,16.,12.,14.,12.,36./, &
           z53/1.175, &
           3.73,3.73,4.804,4.804,5.385,5.941,9.014,11.24,9.014,10.15, &
           9.104,23.28/,z43/0.52,1.31,1.31,1.488,1.488,1.61,1.721,2.577, &
           3.025,2.577,2.81,2.577,5.668/,z23/-0.413,-0.655,-0.655,-0.643, &
           -0.643,-0.659,-0.673,-0.889,-0.946,-0.889,-0.92,-0.889,-1.36/, &
           z86/1.630,5.917,5.917,8.302,8.302,9.520,10.716,16.192,20.978, &
           16.192,18.606,16.192,45.6635/, &
           years_per_sec_over_amu/5.240358d-8/
      double precision :: term, ion_mean_weight_inverse, &
           electron_mean_weight_inverse, xtr, zeta_sum, &
           electron_number_density_na, dd, density, t9, t9_p13, t9_p23, &
           t9_m13, t9_m23, t9_m1, t9_m2, t9_m12, t9_m32
      double precision :: dgdeut, qrtdeut, en, rdeut, qdeut, zz, &
           tfacdeut, tfacdeut2, rdeutmax, rdeut2
      double precision :: pfmc2, efmkt, fprf, degd
      double precision :: xxl, xxl6, xxl8, zcurl, zbar, z58, z28, z33, tm1
      double precision :: uwk, uint, ustr
      double precision :: r1, r2, a1, a2, a3, a4, a5, dr1, da1
      double precision :: be7electron, be7proton, temp3, qrbe7, &
           zprdbe7p, z86be7p, utotbe7p, camube7
      double precision :: f1, f2, f3, f4, o16gamma, c12alpha
      double precision :: convert, egdeut
      double precision :: sum2, sum3
      double precision :: fourpiau2, q6hep, zprdhe3p, z86he3p, utothe3p
      double precision :: flux_value
      double precision :: carbon_fraction_total, oxygen_fraction_total, &
           neutrino_temp, neutrino_density, neutrino_loss_snu
      double precision :: el, eli, ez, ez3, emue, ex1, ex2, ex3, el2
      double precision :: polx10, polx11, polx12, polx21, polx22, polx31, &
           polx32, qedn, qetn, qetnx, qednx
      integer :: i, k, nz

      call setup_abundances_and_composition
! The reaction-rate / screening / energy-release block below stays
! inline: extracting it into internal subroutines perturbs the
! compiler's floating-point instruction scheduling enough to shift
! the last bit of the pp-chain rates, breaking the byte-identical
! regression pin (measured, 2026). Sections that do split cleanly
! (setup, energy totals, neutrino fluxes) are internal subroutines
! under contains.
! SET RATES EQUAL TO ZERO FOR THE LOG_10(T) < 6.0.
! REPLACED FIXED 1 MILLION K THRESHOLD WITH TCUT(1).
!      IF(TL.LE.6.0) THEN
      if (log_temperature.le.star%ctrl%tcut(1)) then
! MHP 5/02 DEUTERIUM BURNING
         dgdeut = 0.0d0
         qrtdeut = 0.0d0
         en = -20.
         dlnepsilon_dlnrho = 0.
         dlnepsilon_dlnt = 0.
         do i = 1,nrxns
            eg(i) = 0.
            reaction_rate(i) = 0.
            reaction_energy_gen(i) = 0.
         end do
         continue
         return
      end if
! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
!  POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
!  OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
!  CGS UNITS.
      density=exp(ln10*dd)
      t9 = exp(ln10*(log_temperature - 9.0d0))
      t9_p13 = t9**cc13
      t9_p23 = t9_p13**2
      t9_m13=1./t9_p13
      t9_m23=t9_m13**2
      t9_m1=1./t9
      t9_m2=t9_m1**2
      t9_m12=1./dsqrt(t9)
      t9_m32=t9_m1*t9_m12
! MHP 5/02 ADD DEUTERIUM BURNING TERM TO THE CODE
! IF DEUTERIUM IS ABOVE A MINIMUM THRESHOLD VALUE.
! RDEUT IS THE RATE (EXCLUDING FACTORS OF THE
! ABUNDANCES) AND QRTDEUT IS THE DERIVATIVE W/R/T T.
! NOTE THAT SCREENING IS EXCLUDED - REASONABLE GIVEN
! THE LOW TEMPERATURES INVOLVED.
      if (deuterium_fraction.le.1.0d-11) then
         rdeut = 0.0d0
         qrtdeut = 0.0d0
      else
! ENERGY YIELD FOR DEUTERIUM BURNING
        qdeut = 5.494d0
        zz = -3.72d0*t9_m13
        tfacdeut = 1.0d0+0.112d0*t9_p13+3.38d0*t9_p23+2.65d0*t9
! FACTOR OF 6.023D23/ REFLECTS AVAGADROS NUMBER DIVIDED BY THE
! MASS OF THE DEUTERON IN AMU
        rdeut = density*2.240d3*t9_p23*exp(zz)*tfacdeut*6.023d23/ &
             atomic_mass_amu(3)
        tfacdeut2 = 0.112d0*t9_p13+6.76d0*t9_p23+7.95d0*t9
        qrtdeut = cc13*((tfacdeut2/tfacdeut) -2.0d0 - zz)
! NOW LIMIT DEUTERIUM BURNING IN A SURFACE CZ TO BE ON A TIME SCALE
! NO SHORTER THAN THE CONVECTIVE OVERTURN TIMESCALE.
        if (shell_index.ge.star%light_burn%jcz .and. star%turnover%convective_turnover_timescale.gt.1.0d0) then
           rdeutmax = 6.023d23/atomic_mass_amu(3)/star%turnover%convective_turnover_timescale
           rdeut2 = rdeut*hydrogen_fraction
           if (rdeut2.gt.rdeutmax) then
! JVS 0712 Commented out write command
!              WRITE(*,*)RDEUT2,RDEUTMAX
              if (hydrogen_fraction.gt.1.0d-6) then
                 rdeut = rdeutmax/hydrogen_fraction
              end if
           end if
         end if
      end if
! ***********************************
! F PRIME/F
! ***********************************
! THE NEXT PIECE OF CODE COMPUTES FIRST THE FERMI ENERGY DIVIDED BY KT, WHERE
!  PFMC2 IS THE FERMI MOMENTUM DIVIDED BY MASS OF ELECTRON TIMES C, ALL
!  SQUARED AND EFMKT IS THE FERMI ENERGY DIVIDED BY KT.  THE QUANTITY
!  FPRF IS THE RATIO OF F PRIME TO F IN SALPETER'S SCREENING CORRECTION.
!  THE VALUE OF FPRF IS DETERMINED IN THE INTERMEDIATE CASE BY AN
!  INTERPOLATION FORMULA DEPENDING UPON THE DEGREE OF DEGENERACY, AS
!  MEASURED BY DEGD = LOG10(E_F/KT). THE NUMERICAL VALUES FOR THE FIT
!  WERE TAKEN FROM SALPETER'S ORIGINAL PAPER, FIGURE 1. THE ONLY CHANGES
!  IN THIS PART OF THE SUBROUTINE WERE THE CORRECTION OF THE ERROR IN
!  THE DEFINITION OF RWE (SEE ABOVE) AND REFINEMENTS OF THE COEFFICIENTS
!  IN THE EXPRESSIONS FOR PFMC2 AND EFMKT.
      pfmc2=1.017677e-4*electron_number_density_na**0.6666667
      efmkt=5.92986*t9_m1*(dsqrt(1.+pfmc2)-1.)
      if (efmkt.le.1.e-2) then
         fprf=1.0
      else
         degd=dlog10(efmkt)
         if (degd.ge.1.5) then
            fprf=0.0
         else
            fprf=0.75793-0.54621*degd-0.30964*degd**2+0.12535*degd**3+ &
            0.1203*degd**4-0.012857*degd**5-0.014768*degd**6
         end if
      end if
! END OF CALCULATION OF FERMI ENERGY DIVIDED BY KT AND THE INTERPOLATION
!  FORMULA FOR (F PRIME/F), WHICH APPEARS IN SALPETER'S SCREENING
!  FORMULA.
! NOW WE GET TO THE COMPUTATION OF WEAK SCREENING (SEE ALSO PAGE 61 OF
!  NEUTRINO ASTROPHYSICS, WHICH GIVES ONLY A SIMPLIFIED FORMULA) AND
!  ALSO THE MORE COMPLICATED INTERMEDIATE AND STRONG SCREENING CASES
!  (SEE REFERENCES TO AP. J. 181 ABOVE FOR INTERMEDIATE
!  AND STRONG SCREENING), ESPECIALLY PAGE 465.
! XXL IS USED IN ALL THE SCREENING FORMULAE.  THIS QUANTITY IS THE
!  IS THE FUNCTION CALLED LAMBDA SUB ZERO BY GRABOSKE ET AL. THE
!  QUANTITY XXL**0.86 = XXL8 IS USED IN CALCULATING INTERMEDIATE
!  SCREENING.
! ZCURL IS THE QUANTITY DEFINED BY EQUATION (4) OF DEWITT ET AL; IT
!  IS THEIR Z WITH A CURLY SYMBOL ON ITS TOP. ZCURL IS USED IN WEAK
!  AND IN INTERMEDIATE SCREENING AND WAS FIRST DEFINED BY SALPETER.
!  ZCURL IS THE SAME AS SALPETER'S ZETA EXCEPT FOR THE FACTOR OF 1/AMU.
! (Z SUB 1 TIMES Z SUB 2)*XXL*ZCURL GIVES THE WEAK SCREENING FACTOR,
!  THE SAME AS SALPETER OR AS EQUATION (19) OF GRABOSKE ET AL.
! Z BAR IS THE SAME AS Z BAR OF DE WITT ET AL AND IS THE AVERAGE CHARGE
!  OF THE IONS.  IT IS EQUAL TO EMU/AMU.  THE QUANTITY (Z BAR)**0.28 =
!  Z28 OCCURS IN THE COMPUTATION OF INTERMEDIATE SCREENING.
! XXL6 IS USED FOR COMPUTING STRONG STRONG SCREENING.
! THE NOTATION USED HERE IS EXPLAINED IN LARGE PART BY
!  EQUATION (4) OF DEWITT ET AL., AP. J. 181, PAGE 439.
! THE FINAL EXPRESSION FOR WEAK SCREENING IS EXACTLY EQUAL TO SALPETER'S
!  FORMULA, WHICH INCLUDES A DEGENERACY CORRECTION. THE MORE GENERAL
!  EXPRESSIONS ARE GIVEN IN TABLE 4 AND EQUATION (19) OF GRABOSKE ET AL.
      xxl=5.9426e-6*t9_m32*dsqrt(density*ion_mean_weight_inverse)
      xxl6=xxl**0.666667
      xxl8=xxl**0.86
      zcurl=dsqrt((zeta_sum+fprf*electron_mean_weight_inverse)/ &
           ion_mean_weight_inverse)
      zbar=electron_mean_weight_inverse/ion_mean_weight_inverse
      z58=zcurl**0.58
      z28=zbar**0.28
      z33=zbar**cc13
      tm1=xxl*zcurl
! COMPUTE SCREENING FOR EACH OF THE REACTIONS.
      do i=1,nrxns
         uwk=tm1*charge_product(i)
         if (uwk.le.star%ctrl%weak_screening_threshold) then
! WEAKSCREENING IS A NUMERICAL PARAMETER PASSED IN THE FLUX COMMON
!  BLOCK. TO OBTAIN THE GRABOSKE ET AL. AND SALPETER STANDARD RESULTS,
!  USE: WEAKSCREENING = 0.03.  FOR THE STANDARD SOLAR MODEL, THIS IS THE
!  VALUE THAT SHOULD BE ADOPTED. TO INVESTIGATE THE EFFECT OF ALWAYS USING
!  WEAK SCREENING, USE A LARGE VALUE FOR WEAKSCREENING, E. G., 30.  AS
!  LONG AS WEAKSCREENING IS ASSUMED TO BE BIGGER THAN ONE, THE PROGRAM
!  WILL ALWAYS CALCULATE FOR THE SUN WITH THE WEAK SCREENING
!  APPROXIMATION.
! UTOT IS THE FINAL SCREENING CORRECTION WHICH APPEARS IN THE
!  RATE EXPRESSION AS: EXP(UTOT) .
! DSCR IS THE LOGARITHMIC DERIVATIVE WITH RESPECT TO DENSITY OF THE
!  SCREENING CORRECTION, D LOG (E^(U_TOT)) /D LOG RHO .
! DSCT IS THE LOGARITHMIC DERIVATIVE OF THE SCREENING WITH RESPECT TO T,
! D LOG (E^(U_TOT)) /D LOG T .  THE LIMITING FORMULAE GIVEN IN THE FIRST
! OPTION ARE OBVIOUS SINCE IN THE WEAK LIMIT SALPETER'S FORMULA SHOWS
! THAT U IS PROPORTIONAL TO THE SQUARE ROOT OF RHO TIMES T^(-3/2).
            screening_factor(i)=uwk
            dscreen_dlnrho(i)=0.5*uwk
            dscreen_dlnt(i)=-1.5*uwk
         else
            uint=0.38*xxl8*xtr*z86(i)/(ion_mean_weight_inverse*z58*z28)
            if (uwk.le.2.) then
               screening_factor(i)=uint
               dscreen_dlnrho(i)=0.43*uint
               dscreen_dlnt(i)=-1.29*uint
            else
               ustr=0.624*z33*xxl6*(z53(i)+0.316*z33*z43(i)+0.737* &
               z23(i)/(zbar*xxl6))
               if (ustr.lt.uint.or.uwk.ge.5.) then
                  screening_factor(i)=ustr
                  dscreen_dlnrho(i)=0.208*z33*(z53(i)+0.316*z33*z43(i))*xxl6
                  dscreen_dlnt(i)=-3.*dscreen_dlnrho(i)
               else
                  screening_factor(i)=uint
                  dscreen_dlnrho(i)=0.43*uint
                  dscreen_dlnt(i)=-1.29*uint
               end if
            end if
         end if
      end do
! ****************************************************************
! END OF SCREENING CALCULATION. WEAK AND INTERMEDIATE SCREENING FORMS
!  ARE GIVEN CORRECTLY.  STRONG SCREENING WAS NOT CHECKED BECAUSE IT IS
!  NOT RELEVANT FOR THE SUN.
! ****************************************************************
      nz=1
      if (hydrogen_fraction.eq.0.0) then
         f1=0.
         f2=0.
         f3=0.
         f4=0.
      else
      nz=8
! **************************************************************
!  CALCULATE REACTION RATES FOR THE THREE PRINCIPAL RECTIONS OF
!   THE PP CHAIN: PP, HE3+HE3, HE3 +HE4, AND THE FOUR PROTON
!   BURNING REACTIONS ON C12: C13, N14, AND O16.
! **************************************************************
! R1 IS (T SUB 9)^(-3/2) TIMES (S SUB EFF)/(S SUB 0). THE CORRECT
!  EXPRESSION FOR (S SUB EFF)/(S SUB 0) IS GIVEN IN EQUATION 3.14 OF
!  NEUTRINO ASTROPHYSICS.  THE NUMERICAL FORM THAT IS USED IS EQUATION
!  52 OF FOWLER, CAUHLAN, AND ZIMMERMAN, VOL. 5, 1967.
! RATE(I) IS THE RATE OF THE DIFFERENT REACTIONS PER SECOND PER GRAM,
!  EXCEPT THAT THE MASS FRACTIONS ARE OMITTED AT THIS POINT AND PUT IN
!  LATER.
! DRATT IS LOGARITHMIC DERIVATIVE OF RATE WITH RESPECT TO TEMPERATURE,
!  D LOG RATE DIVIDED BY D LOG T, LOG TO BASE 10.
! DRATRO IS THE LOGARITHMIC DERIVATIVE OF THE RATE WITH RESPECT TO
!  DENSITY, D LOG RATE/D LOG RHO, LOG TO BASE 10.
! THE PREVIOUS YALE VERSION HAD THE RATE FOR THE O16 + P REACTION
!  MULTIPLIED BY T9**(-1/7). THIS FACTOR IS INCORRECT AND HAS BEEN
!  REMOVED; IT APPEARED BEFORE AS AN IF STATEMENT REFERRING ONLY TO
!  RATE(7).
      do i=1,7
!         R1=T9M23+Q1(I)*T9M13+Q2(I)+Q3(I)*T9P13+Q4(I)*T9P23+Q5(I)*T9
! MHP 8/14 RATES CORRECTED TO PERMIT USER MODIFICATION OF REACTION
! RATE DERIVATIVES
         r1=t9_m23+q1(i)*t9_m13+star%qs0e_scale(i)*(q2(i)+q3(i)*t9_p13)+ &
              star%qqs0ee_scale(i)*(q4(i)*t9_p23+q5(i)*t9)
         reaction_rate(i)=density*r1*exp(q6(i)*t9_m13+q7(i)+(q8(i)*t9)**2+ &
              screening_factor(i))
         reaction_rate(i) = reaction_rate(i)*star%cross_section_scale(i)
         if (reaction_rate(i).lt.1.e-30) then
            reaction_rate(i)=0.
            dlnrate_dlnt(i)=0.
         else
            dlnrate_dlnrho(i)=1.+dscreen_dlnrho(i)
            dlnrate_dlnt(i)=dscreen_dlnt(i)-(q6(i)*t9_m13+(2.*t9_m23+ &
            q1(i)*t9_m13-q3(i)*t9_p13-2.*q4(i)*t9_p23-3.*q5(i)*t9)/r1)/3.+ &
            2.*(q8(i)*t9)**2
         end if
      end do
! ***************************************************************
! END OF CALCULATION OF REACTION RATES FOR FIRST 7 REACTIONS.
! ***************************************************************
! ***********************************************
! BE7 BURNING
! ***********************************************
! CALCULATE THE BURNING OF BE7 BY PROTONS (WHICH PRODUCES THE MOST
!  EXPERIMENTALLY ACCESSIBLE SOLAR NEUTRINOS) AND THE BURNING OF BE7
!  BY ELECTRON CAPTURE (THE DOMINANT PROCESS).
!  THE ELECTRON CAPTURE RATE IN SEC^(-1) IS GIVEN BY EQUATION (3.18) OF
!  NEUTRINO ASTROPHYSICS.  CAN OMIT THE FACTOR OF RHO (IN CGS UNITS)
!  WHICH ALSO APPEARS IN THE BE7PROTON RATE. THE BE7 PROTON CAPTURE
!  RATE IS GIVEN BY TABLE 3.2 AND EQUATION (3.12).  NOTE THAT 1 OVER MU
!  SUB E IS EQUAL TO EMU.
! HERE WE USE THE NOTATION BE7 + PROTON = BE7PROTON AND
!  BE7 + E = BE7ELECTRON.
! F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY ELECTRON CAPTURE.
! F2 IS THE FRACTION OF THE BE7 THAT IS BURNED BY PROTON CAPTURE.
! F3 IS THE FRACTION OF THE N14 THAT IS BURNED BY P, ALPHA REACTION.
! F4 IS THE FRACTION OF THE N15 THAT IS BURNED BY P, GAMMA REACTION.
! SEE TABLE 21 OF BAHCALL AND ULRICH (1988), REV. MOD. PHYS. 60.
! 10/13/97. I changed Temp3 (i.e., tau) by 5/10^5 as a result of
! using bare nuclear masses. Previously coefficient was -10.26202.
! I did not change Be7electron for the slightly different S0 for
! electron capture, since that is done in SStandard. Previous coefficient
! in Be7electron expression was (3.126571E+5). 10/14/97.
!
         be7electron = (1.752e-10)*t9_m12*(1.0 + 0.004*(1000.*t9 - 16.))
         be7electron = be7electron*electron_mean_weight_inverse*star%cross_section_scale(15)
         temp3 = (-10.2625*t9_m13)
         be7proton = (3.128813e+5)*hydrogen_fraction*star%cross_section_scale(16)*exp(temp3)
! INCLUDE FOR BE7PROTON THE T9M23 FACTOR AND ALL CORRECTIONS PROPORTIONAL TO
!  Q1,...,Q5 FROM EQUATION 3.14 OF NEUTRINO ASTROPHYSICS. THESE
!  CORRECTIONS ARE DEFINED EARLIER IN THIS SUBROUTINE.
!         QRBE7 = T9M23 + Q1(8)*T9M13 + Q2(8)+ Q3(8)*T9P13
!     $          + Q4(8)*T9P23 + Q5(8)*T9
! MHP 9/14 ADDED THE ABILITY TO ALTER DERIVATIVES INDEPENDENTLY
         qrbe7 = t9_m23 + q1(8)*t9_m13 + star%qs0e_scale(8)*(q2(8)+ q3(8)*t9_p13) &
              + star%qqs0ee_scale(8)*(q4(8)*t9_p23 + q5(8)*t9)
         be7proton = be7proton*qrbe7
! CALCULATE THE SCREENING CORRECTION FOR BE7 + P REACTION.  USE WEAK AND
!  INTERMEDIATE SCREENING FORMULAE.
         zprdbe7p = 4.0
         z86be7p = 5.7790
         uwk = tm1*zprdbe7p
         if (uwk.le.star%ctrl%weak_screening_threshold) then
            utotbe7p = uwk
         else
            uint = 0.38*xxl8*xtr*z86be7p/(ion_mean_weight_inverse*z58*z28)
            utotbe7p = uint
          end if
         be7proton = be7proton*exp(utotbe7p)
! END OF CALCULATION OF SCREENING CORRECTION FOR BE7 + P REACTION.
! MULTIPLY RATES BY FACTOR OF RHO/[(ATOMIC MASS UNIT)*A(BE7)] TO GET
!  IN UNITS OF GM^{-1}.  CALL FACTOR CAMUBE7. Corrected for Be7 bare
! mass number. Previously used 8.582295E+22 in the expression below.
! These corrections were made on 10/14/97.
!
         camube7 = density*8.584981e+22
         be7proton = camube7*be7proton
         be7electron = camube7*be7electron
! END OF MULTIPLICATION INSERTED NOVEMBER 6, 1990.
         f1 = be7electron/(be7electron + be7proton)
         f2 = be7proton/(be7electron + be7proton)
! *********************************************************************
! END OF CALCULATION OF CRUCIAL BE7 ELECTRON CAPTURE AND PROTON CAPTURE
!  RATES AND THEIR RATIO.
! *********************************************************************
! ***************************
! N15 + P BRANCHING.
! ***************************
! THE FOLLOWING STATEMENTS COMPUTE THE BRANCHING OF N15(P,ALPHA)C12
!  AND N15(P,GAMMA)O16. THESE STATEMENTS REPLACE OUTDATED STATEMENTS
!  IN THE YALE CODE. THE CNO CROSS-SECTION FACTORS ARE FROM TABLE 3.4
!  OF NEUTRINO ASTROPHYSICS.  THE RATIO OF THE REACTIONS DEPENDS ONLY
!  UPON THE EFFECTIVE ZERO ENERGY S-FACTOR, WHICH IS S0(ZERO ENERGY)
!  TIMES THE COMBINATION OF TEMPERATURE AND S-FACTOR DERIVATIVES/S0
!  THAT WAS USED PREVIOUSLY AS R1 IN THE RATE CALCULATIONS. THE
!  NUMERICAL COEFFICIENTS THAT APPEAR IN THE RATE WERE REPRESENTED
!  BY THE Q1(J),...Q5(J) FOR THE OTHER REACTIONS.
!  THE QVALUES FOR THE N15 REACTIONS HAVE BEEN COMPUTED SEPARATELY.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! To agree with Solar Fusion Workshop paper, the value of S0 in keV-b
! has been changed from 78000. to 675000. on 9/25/97. On 10/14/97,
! JNB changed cues.f so as to compute the Q-coefficients for the N15 + p
! reactions.  Also, checked that the coefficients are the same as the
! ones given earlier in energy.f when we use the older CNO data.
!
      o16gamma = t9_m23 + 0.0273016*t9_m13 + 0.14374 + 0.027490*t9_p13 &
                + 6.14685*t9_p23 + 2.98940*t9
! MULTIPLY BY THE VALUE OF S0 IN KEV-B.
!      O16GAMMA = O16GAMMA*64.
! MHP 9/14 ADDED THE OPTION TO MODIFY THE RELATIVE CROSS SECTIONS
! FOR N15+P -> C12+ALPHA AND O16+GAMMA
      o16gamma = o16gamma*64*star%o16_gamma_scale
!
      c12alpha = t9_m23 + 0.0273016*t9_m13 + 2.01186 + 0.384763*t9_p13 &
                + 17.0579*t9_p23 + 8.29580*t9
!      C12ALPHA = C12ALPHA*67500
      c12alpha = c12alpha*67500*star%c12_alpha_scale
      f3 = c12alpha/(c12alpha + o16gamma)
      f4 = 1.0d0 - f3
! END OF NEW ROUTINE FOR THE BRANCHING OF N15 + P .
      end if
      do i=nz,nrxns
         reaction_rate(i)=0.
         dlnrate_dlnrho(i)=0.
         dlnrate_dlnt(i)=0.
   end do
! ***MHP 3/91 ALPHA CAPTURE REACTIONS UPDATED TO CAUGHLAN AND FOWLER(1988)
!    RATES.  THE RATES ARE EXPRESSED IN THE SAME TERMS USED BY CZ, WITH
!    THE CONVERSION FACTOR IN THE FRONT OBTAINED FROM VANDENBERG'S
!    NOTES ON THE REACTION RATES.
!  RATE(8)  HE4+C13
!  RATE(10) HE4+C12=>O16
!  RATE(11) HE4+N14=>O18
!  RATE(12) TRIPLE ALPHA
      if (log_temperature.ge.star%ctrl%tcut(4)) then
! C13(ALPHA,N) O16
      r1=t9_m23+0.0129d0*t9_m13+2.04d0+0.184d0*t9_p13
      a1 = 6.77d15*exp(-32.329d0*t9_m13-(t9/1.284d0)**2)
      a2 = 3.82d5*exp(-9.373*t9_m1)
      a3 = 1.41d6*exp(-11.873*t9_m1)
      a4 = 2.00d9*exp(-20.409*t9_m1)
      a5 = 2.92d9*exp(-29.283*t9_m1)
      reaction_rate(8) = 1.157126d22*density*exp(screening_factor(8))* &
           (a1*r1+t9_m32*(a2+a3+a4+a5))
      dlnrate_dlnrho(8)=1.0d0+dscreen_dlnrho(8)
      dr1 = cc13*(-2.0d0*t9_m23-0.0129d0*t9_m13+0.184d0*t9_p13)
      da1 = a1*(cc13*32.329d0*t9_m13 - 2.0d0*(t9/1.284d0)**2)
      dlnrate_dlnt(8) = dscreen_dlnt(8)+density/reaction_rate(8)*(dr1*a1 + &
           r1*da1 + &
           a2*(9.373*t9_m1-1.5d0)+a3*(11.873*t9_m1-1.5d0)+ &
           a4*(20.409*t9_m1-1.5d0)+a5*(29.283*t9_m1-1.5d0))
! C12(ALPHA,GAMMA)O16
      r1 = 1.0d0/(1.0d0+0.0489d0*t9_m23)
      r2 = 1.0d0/(1.0d0+0.2654d0*t9_m23)
      a1 = t9_m2*exp(-32.120*t9_m13)
      a2 = 1.04d8*r1**2*exp(-(t9/3.496)**2)
      a3 = 1.76d8*r2**2
      a4 = 1.25d3*t9_m32*exp(-27.499*t9_m1)
      a5 = 1.43d-2*t9**5*exp(-15.541*t9_m1)
      reaction_rate(10) = 1.25388d22*density*exp(screening_factor(10))* &
           (a1*(a2+a3)+a4+a5)
      dlnrate_dlnrho(10) = 1.0d0+dscreen_dlnrho(10)
      dlnrate_dlnt(10) = dscreen_dlnt(10)+density/reaction_rate(10)* &
           (a1*((cc13*32.120*t9_m13-2.0d0)* &
           (a2+a3)+a2*(r1*cc13*0.1956-2.0d0*(t9/3.496)**2)+a3* &
           (cc13*1.0616d0*r2))+a4*(27.499*t9_m1-1.5d0)+a5* &
           (5.0d0+15.541*t9_m1))
! N14(ALPHA,GAMMA)F18 + F18=>O18+EPLUS+NU
      r1 = t9_m23+0.012d0*t9_m13+1.45d0+0.177d0*t9_p13+1.97d0*t9_p23 &
           +0.406d0*t9
      a1 = 7.78d9*exp(-36.031d0*t9_m13-(t9/0.881d0)**2)
      a2 = t9_m32*2.36d-10*exp(-2.798d0*t9_m1)
      a3 = t9_m32*2.03d0*exp(-5.054d0*t9_m1)
      a4 = t9_m23*1.15d4*exp(-12.310*t9_m1)
      reaction_rate(11)= 1.07452d22*density*exp(screening_factor(11))* &
           (a1*r1+a2+a3+a4)
      dlnrate_dlnrho(11)=1.+dscreen_dlnrho(11)
      dr1 = cc13*(-2.0d0*t9_m23-0.012d0*t9_m13+0.177d0*t9_p13+ &
            3.94d0*t9_p23)+0.406d0*t9
      da1 = a1*(cc13*36.031d0*t9_m13-2.0d0*(t9/0.881d0)**2)
      dlnrate_dlnt(11) = dscreen_dlnt(11)+density/reaction_rate(11)*(dr1*a1+ &
              r1*da1+a2* &
              (2.798d0*t9_m1-1.5d0)+a3*(5.054d0*t9_m1-1.5d0)+ &
              a4*(12.310d0*t9_m1-cc23))
! TRIPLE ALPHA
      reaction_rate(12) = 1.565315d21*density**2*t9_m1*t9_m2*2.79e-8* &
                 exp(-4.4027*t9_m1+screening_factor(12))
      dlnrate_dlnrho(12) = 2.0d0+dscreen_dlnrho(12)
      dlnrate_dlnt(12) = -3.0d0+dscreen_dlnt(12)+4.4027d0*t9_m1
! *******************
! EG(I)
! *******************
! MULTIPLY THE RATES PER GRAM, RATE(I), BY THE ABUNDANCES OF THE
!  REACTING SPECIES BY MASS, TO GET THE TOTAL RATES PER GRAM, EG.
      end if
      eg(1)=reaction_rate(1)*hydrogen_fraction*hydrogen_fraction
! MHP 5/02 ADD DEUTERIUM BURNING IF RELEVANT
      if (deuterium_fraction.gt.1.0d-11) then
         egdeut = rdeut*hydrogen_fraction*deuterium_fraction
      else
         egdeut = 0.0d0
      end if
      eg(2)=reaction_rate(2)*he3_fraction*he3_fraction
      eg(3)=reaction_rate(3)*he3_fraction*helium_fraction
      eg(4)=reaction_rate(4)*hydrogen_fraction*c12_fraction
      eg(5)=reaction_rate(5)*hydrogen_fraction*c13_fraction
      eg(6)=reaction_rate(6)*hydrogen_fraction*n14_fraction
      eg(7)=reaction_rate(7)*hydrogen_fraction*o16_fraction
      eg(8)=reaction_rate(8)*helium_fraction*c13_fraction
!     EG(9)=RATE(9)*Y*X016
      eg(10)=reaction_rate(10)*helium_fraction*c12_fraction
      eg(11)=reaction_rate(11)*helium_fraction*n14_fraction
      eg(12)=reaction_rate(12)*helium_fraction**3
!     EG(13)=RATE(13)*XC12*XC12
! ******************************************************************
! ****************************************
      call compute_energy_generation
      call compute_neutrino_emission
      return

contains


! ---------------------------------------------------------------
! Zero the neutrino/alpha-capture yields, form the fractional
! abundances of the burning species, the ion and electron mean
! molecular weights, and the screening precursors (xtr, zet).
subroutine setup_abundances_and_composition
! ZERO OUT THE ENERGY YIELDS FROM NEUTRINOS(ENU) AND ALPHA CAPTURE
! REACTIONS (EALPCA).
      star%engeb%neutrino_loss_rate = 0.0d0
      star%engeb%alpha_capture_energy = 0.0d0
! DEFINE NEXT THE FRACTIONAL ABUNDANCES BY MASS OF THE IMPORTANT
!  ISOTOPES.
! X, Y, Z, XHE3,..., XBE9 ARE THE MASS FRACTIONS OF THE ISOTOPES.
!  THE ABUNDANCES OF NEUTRONS, H2, H3, NE20,AND MG24, WHICH ARE,
!  RESPECTIVELY, XFRAC(I) FOR I = 1,3,4,12,13, ARE NO LONGER USED.
      mass_fraction(1) = 0.0
      mass_fraction(2) = hydrogen_fraction
! MHP 5/02 ADDED DEUTERIUM
!      XFRAC(3) = 0.0
      mass_fraction(3) = deuterium_fraction
      mass_fraction(4) = 0.0
      mass_fraction(5) = he3_fraction
      mass_fraction(6) = helium_fraction
      mass_fraction(7) = c12_fraction
      mass_fraction(8) = c13_fraction
      mass_fraction(9) = n14_fraction
      mass_fraction(10) = o16_fraction
      mass_fraction(11) = o18_fraction
      mass_fraction(12) = 0.0
      mass_fraction(13) = 0.0
! *******************************************************************
! BEGIN CALCULATION OF SCREENING CORRECTION.
! *******************************************************************
!  THE BASIC REFERENCES ARE SALPETER, AUSTRALIAN JOURNAL OF PHYSICS,
!  VOL. 7, 373 (1954). THE FORMULA FOR WEAK SCREENING THAT IS BEING
!  PROGRAMMED IS EQUATION (25) OF THIS PAPER. THE OTHER IMPORTANT
!  REFERENCES ARE: DEWITT, GRABOSKE, AND COOPER, AP. J. 181, 439 (1973)
!  AND GRABOSKE ET AL., AP. J. 181, 457 (1973). THE VALUES OF EMU AND
!  ZET ARE ESSENTIAL FOR COMPUTING WEAK SCREENING; THE VALUE OF AMU IS
!  USED IN AN NON-ESSENTIAL WAY IN THIS COMPUTATION. XTR IS USED IN
!  COMPUTING INTERMEDIATE SCREENING.
! AMU IS ONE OVER THE MEAN MOLECULAR WEIGHT OF THE IONS, MU SUB I .
! EMU IS ONE OVER THE ELECTRON MEAN MOLECULAR WEIGHT, MU SUB E.
!  EMU IS USED HERE AS THE NAME FOR THE SECOND PART OF THE ZETA FUNCTION
!  IN THE SALPETER EXPRESSION FOR WEAK SCREENING.
! XTR IS USED LATER IN THE INTERMEDIATE SCREENING CALCULATION.  THE AVERAGE
!  OF THE QUANTITY Z**(3B -1) IS EQUAL TO XTR/AMU.
! ZET IS THE FIRST PART OF THE SALPETER SCREENING ZETA VARIABLE.
! MU = SUM OVER I OF [X(I)/A(I)].
! MU SUB E = SUM OVER I [ Z(I)*X(I)/A(I)].
      ion_mean_weight_inverse = 0.
      electron_mean_weight_inverse = 0.
      xtr = 0.
      zeta_sum = 0.
      do i = 1,num_isotopes
         term = mass_fraction(i)/atomic_mass_amu(i)
         ion_mean_weight_inverse = ion_mean_weight_inverse+term
         electron_mean_weight_inverse = electron_mean_weight_inverse+ &
              term*atomic_number(i)
         xtr = xtr+term*atomic_number(i)**1.58
         zeta_sum = zeta_sum+term*atomic_number(i)**2
      end do
! DL AND DT ARE THE THE LOG10 OF THE DENSITY AND TEMPERATURE.
!  THE UNIT OF TEMPERATURE IS 10^9 K AND THE UNIT OF DENSITY IS
!  GM PER CM^3 .
! PDT AND PDP ARE THE DERIVATIVES OF THE DENSITY WITH RESPECT TO
!  TEMPERATURE AND DENSITY.
! DD = LOG RHO TO THE BASE 10.
! CLN = LN10.  CLN IS CONVERSION BETWEEN LOG10 AND LN.
! CONVERT DENSITY TO UNLOGGED FORM.
! RWE = RHO/(MU SUB E), I. E., THE NUMBER OF ELECTRONS DIVIDED BY
!  AVOGADRO'S NUMBER.
      electron_number_density_na = ( exp(ln10*log_density) )* &
           electron_mean_weight_inverse
! THE EXPRESSION FOR RWE WAS INCORRECT IN THE ORIGINAL YALE SUBROUTINE.
!  THE ORIGINAL VERSION HAD ( EXP(CLN*DL) ) DIVIDED BY EMU INSTEAD OF
!  MULTIPLIED BY EMU.  RWE IS USED LATER IN COMPUTING THE SCREENING
!  CORRECTION.
      dd = log_density
end subroutine setup_abundances_and_composition

! ---------------------------------------------------------------
! Energy release per reaction (MeV -> erg), reactions 9 and 13
! zeroed as in the original Yale code, the total energy generation
! with d ln eps / d ln rho and d ln eps / d ln T, and the global
! outputs returned to the caller.
subroutine compute_energy_generation
! ENERGY GENERATION.
! ****************************************
! CALCULATE ENERGY GENERATION BY MULTIPLYING RATES PER GRAM PER SEC BY
!  THE ENERGY RELEASE.  THE ENERGIES ARE TAKEN FROM TABLE 21 OF BAHCALL AND
!  ULRICH (1988), REV. MOD. PHYS. 60, 297. THIS TABLE IS BASED UPON A CAREFUL
!  CALCULATION OF THE AVERAGE AMOUNT OF ENERGY LOSS BY NEUTRINOS FOR
!  EACH REACTION. THE NUMBERS FOR THE C12 + P REACTION SEQUENCE AND THE
!  C13 + P REACTION ARE BROKEN DOWN SEPARATELY FOR THIS
!  SUBROUTINE.
! THE FINAL NUMBERS ARE IN ERG PER GM PER SECOND.
! DEFINE THE CONSTANT TO CONVERT MEV'S TO ERGS. THE NUMBERS THAT APPEAR
!  ARE IN MEV SO THEY CAN BE EASILY IDENTIFIED.
      convert = 1.602177e-6
! THE MULTIPLYING CONSTANTS BELOW ARE IN MEV.
! JNB changed the pp energy release by 0.002 MeV because of a better
! estimate of the neutrino energy loss on 9/25/97. See pg. 139 of
! Vol. 19 of my notes.  On pgs. 139-141, I document other small changes
! to this energy generation. No large changes; all of order keV changes
! except for the rare 8B reaction.  9/28/97.
!
      reaction_energy_gen(1)=eg(1)*6.664*convert
! MHP 5/02 ADD DEUTERIUM BURNING
      if (deuterium_fraction.gt.1.0d-11) then
         dgdeut = egdeut*qdeut*convert
      else
         dgdeut = 0.0d0
      end if
      reaction_energy_gen(2)=eg(2)*12.860*convert
      reaction_energy_gen(3)=eg(3)*(1.586+f1*17.394+f2*11.499)*convert
      reaction_energy_gen(4)=eg(4)*3.457372*convert
      reaction_energy_gen(5)=eg(5)*7.550628*convert
      reaction_energy_gen(6)=eg(6)*(9.054+f3*4.966+f4*12.128)*convert
      reaction_energy_gen(7)=eg(7)*3.553*convert
      reaction_energy_gen(8)=eg(8)*2.216*convert
      reaction_energy_gen(10)=eg(10)*7.162*convert
      reaction_energy_gen(11)=eg(11)*5.815*convert
      reaction_energy_gen(12)=eg(12)*7.275*convert
! JVS 10/11 Need to grab He3 energy generation
      star%engeb%he3_he3_energy_rate = reaction_energy_gen(2)
      star%engeb%he3_burning_energy_rate = reaction_energy_gen(2)+reaction_energy_gen(3)
! JVS end

! *******************************************************************
! END OF CALCULATION OF ENERGY RELEASE.
! *******************************************************************
! SET TO ZERO O16+ALPHA AND C12+C12 RATES.
      dlnrate_dlnrho(9) = 0.0
      dlnrate_dlnt(9) = 0.0
      reaction_energy_gen(9) = 0.0
      dlnrate_dlnrho(13) = 0.0
      dlnrate_dlnt(13) = 0.0
      reaction_energy_gen(13) = 0.0
! END OF XEROING OUT OF REACTIONS 9 AND 13.
      total_energy_gen_rate=0.0
      sum2=0.0
      sum3=0.0
      do i=1,nrxns
! *******************************************************************
! SUM OF THE TOTAL ENERGY GENERATION IN ERGS PER GRM PER SECOND WITH
! DERIVATIVES WITH RESPECT TO DENSITY AND TO TEMPERATURE.
! *******************************************************************
! SUM1 = SUM OF ALL ENERGY GENERATION RATES. NOTE THAT THE BURNING OF
!  BE7 IS INCLUDED IN DG(3) ABOVE.
! SUM2 = SUM OVER I OF DG(I)* [D LOG RATE(I) / D LOG RHO ].
! SUM3 = SUM OVER I OF DG(I)* [D LOG RATE(I) / D LOG T ].
         total_energy_gen_rate=total_energy_gen_rate+reaction_energy_gen(i)
         sum2=sum2+reaction_energy_gen(i)*dlnrate_dlnrho(i)
         sum3=sum3+reaction_energy_gen(i)*dlnrate_dlnt(i)
      end do
! MHP 5/02 ADD DEUTERIUM BURNING
      total_energy_gen_rate = total_energy_gen_rate + dgdeut
      sum2 = sum2 + dgdeut
      sum3 = sum3 + dgdeut*qrtdeut
      if (total_energy_gen_rate.le.1.e-12) then
         en=-20.
         dlnepsilon_dlnrho=0.
         dlnepsilon_dlnt=0.
         do i=1,nrxns
            eg(i)=0.
         end do
      else
! ******************************************************
! GLOBAL QUANTITIES THAT ARE RETURNED BY THE SUBROUTINE.
! ******************************************************
! PEP AND PET ARE THE DERIVATIVES OF THE TOTAL ENERGY GENERATION RATE
!  WITH RESPECT TO DENSITY AND TEMPERATURE.
! MHP 5/90 CHANGE DERIVATIVES TO BE D LN EPS/D LN RHO AND D LN EPS/D LN T
! TO PUT THEM IN THE SAME FORM AS PRATHER DERIVATIVES.
         dlnepsilon_dlnrho = sum2
         dlnepsilon_dlnt = sum3
      end if
! PDP = D LOG RHO/ D LOG P; PDT = D LOG RHO/ D LOG T.
! *****************************************************
! END OF COMPUTATION OF THE GLOBAL QUANTITIES.
! *****************************************************
      do i=1,nrxns
         if (reaction_rate(i).le.1.e-5) reaction_rate(i) = 0.0
      end do
! ******************************************************
end subroutine compute_energy_generation

! ---------------------------------------------------------------
! Rates per 1e9 yr per amu (hrk), the pp/CNO energy split, and --
! above tcut(5) -- the eight solar neutrino fluxes (pp, pep, hep,
! Be7, B8, N13, O15, F17) with hep screening. The tcut(5) early
! RETURN is equivalent in or out of the section: nothing follows
! this call in engeb.
subroutine compute_neutrino_emission
! RATES PER 10^9 YEARS PER ATOMIC MASS UNIT: HRK(IU)
! ******************************************************
! HR1, ..., HR13 ARE THE RATES OF THE INDIVIDUAL REACTIONS.
!  THE INTERPRETATION OF WHICH REACTION GOES WITH WHICH SYMBOL CAN BE
!  MADE EASILY BY LOOKING AT THE DEFINITIONS OF THE EG(I)'S.
!  THE ABUNDANCES ARE UPDATED IN SUBROUTINE KEMCOM USING THESE MATRICES.
! C21 IS THE PRODUCT OF (10^9 YEARS/1 SECOND)*(1 ATOMIC MASS UNIT/1
!  GRAM). I HAVE USED HERE SIDEREAL YEAR IN CONVERTING TO SECONDS.
      reaction_rate_1(shell_index)=reaction_rate(1)*years_per_sec_over_amu
      reaction_rate_2(shell_index)=reaction_rate(2)*years_per_sec_over_amu
      reaction_rate_3(shell_index)=reaction_rate(3)*years_per_sec_over_amu
      reaction_rate_4(shell_index)=reaction_rate(4)*years_per_sec_over_amu
      reaction_rate_5(shell_index)=reaction_rate(5)*years_per_sec_over_amu
      reaction_rate_6(shell_index)=reaction_rate(6)*years_per_sec_over_amu
      reaction_rate_7(shell_index)=reaction_rate(7)*years_per_sec_over_amu
      reaction_rate_8(shell_index)=reaction_rate(8)*years_per_sec_over_amu
      reaction_rate_9(shell_index)=reaction_rate(9)*years_per_sec_over_amu
      reaction_rate_10(shell_index)=reaction_rate(10)*years_per_sec_over_amu
      reaction_rate_11(shell_index)=reaction_rate(11)*years_per_sec_over_amu
      reaction_rate_12(shell_index)=reaction_rate(12)*years_per_sec_over_amu
      reaction_rate_13(shell_index)=reaction_rate(13)*years_per_sec_over_amu
      n15_alpha_branch_fraction(shell_index)=f3
      be7_electron_capture_fraction(shell_index)=f1
! ****************************************
! END OF COMPUTATION OF HRK(IU).
! ****************************************
! ****************************************
! CALCULATING THE TOTAL ENERGY GENERATION.
! ****************************************
! THE SUMMATION OF THE ENERGIES IS GIVEN IN TABLE 21 OF NEUTRINO
!  ASTROPHYSICS.
! THE ORIGINAL YALE SUBROUTINE CONTAINED SERIOUS ERRORS.  THE
!  CALCULATION OF THE RATE OF ENERGY THROUGH EPP2 AND EPP3 (SEE
!  BELOW) CONTAINED TWO FACTORS OF BRANCHING RATIOS, RATHER THAN
!  THE SINGLE FACTOR THAT SHOULD BE PRESENT. THIS HAD THE EFFECT
!  OF REDUCING ARTIFICIALLY THE ENERGY CALCULATED FROM THESE
!  REACTIONS.
! EPP1 INCLUDES THE ENERGY GENERATED BY THE PP REACTION, BY THE H2 + P
!  REACTION, AND BY THE HE3 + HE3 REACTION.  SEE TABLE 21 OF NEUTRINO
!  ASTROPHYSICS.
      pp_chain_energy_gen = reaction_energy_gen(1)+reaction_energy_gen(2)+dgdeut
!      EPP1 = DG(1)+DG(2)
! EPP3 INCLUDES THE ENERGY GENERATED BY THE HE3 + HE4 REACTION AND BY
!  THE BURNING OF BE7 THROUGH PROTON CAPTURE.
!      EPP3 = EG(3)*(1.586 + F2*11.499)*CONVERT
      he3he4_be7_proton_energy_gen = eg(3)*f2*(1.586 + 11.499)*convert
! EPP2 INCLUDES THE ENERGY GENERATED BY THE HE3 + HE4 REACTION AND BY
!  THE BURNING OF BE7 THROUGH ELECTRON CAPTURE.
!      EPP2 = EG(3)*(1.586 + F1*17.394)*CONVERT
      he3he4_be7_electron_energy_gen = reaction_energy_gen(3) - &
           he3he4_be7_proton_energy_gen
! ECN IS THE ENERGY GENERATED THROUGH THE CNO CYCLE.
      cno_cycle_energy_gen=reaction_energy_gen(4)+reaction_energy_gen(5)+ &
           reaction_energy_gen(6)+reaction_energy_gen(7)
! E3AL IS THE ENERGY GENERATED THROUGH THE TRIPLE-ALPHA REACTION AND
!  IS NEGLIGIBLE FOR THE SUN.
      triple_alpha_energy_gen = reaction_energy_gen(12)


! ENERGY FROM ALPHA CAPTURE REACTIONS.
      star%engeb%alpha_capture_energy=reaction_energy_gen(8)+reaction_energy_gen(10)+ &
           reaction_energy_gen(11)
      if (star%ctrl%lsnu) then
! MHP 9/91 CHANGE TO TURN OFF NEUTRINO CALC FOR HYDROGEN-EXHAUSTED CORE.
         if (hydrogen_fraction.le.1.0d-6) then
            do i=1,10
               star%flux%neutrino_flux(i)=0.0d0
            end do
         else
! ****************************************************************
! CALCULATION OF NEUTRINO FLUXES
! ****************************************************************
! THIS PART OF THE SUBROUTINE CALCULATES THE NEUTRINOS FLUXES IN
!  NUMBER PER GRAM PER SQUARE CENTIMETER PER SECOND AT THE EARTH'S SURFACE
!  (ASSUMING NOTHING HAPPENS TO THE NEUTRINOS AFTER THEY ARE CREATED).
! SEE TABLES 3.1 AND 3.3 OF NEUTRINOS ASTROPHYSICS OR EQUATIONS 6.1-6.8
!  FOR THE REACTIONS. THE ORDER OF THE REACTIONS IS THE SAME AS IN
!  EQUATIONS 6.1-6.8 .
! DEFINE 4*PI*(AU)**2 .
         fourpiau2 = 2.812295e+27
! FLUX OF PP NEUTRINOS.
         star%flux%neutrino_flux(i_nu_pp) = eg(1)/fourpiau2
! FLUX OF PEP NEUTRINOS. USE EQUATION 3.17 OF NEUTRINO ASTROPHYSICS.
! Note that should not change SStandard(14) unless the ratio of pep to pp
!  is changed.  Pep rate is explicitly scaled here with respect to the pp
!  rate.
         star%flux%neutrino_flux(i_nu_pep) = (3.4848e-6)*electron_number_density_na*t9_m12* &
              (1.0 + 20.*t9)*eg(1)
         star%flux%neutrino_flux(i_nu_pep) = star%flux%neutrino_flux(i_nu_pep)*star%cross_section_scale(14)/fourpiau2
! FLUX OF HEP NEUTRINOS.  USE EQUATION 3.12 DIRECTLY.
         q6hep = -6.1399
! Q6 IS THE NEGATIVE OF THE COEFFICIENT OF T9M13 IN TAU, EQUATION 3.10.
         star%flux%neutrino_flux(i_nu_hep) = (1.71724e+11)*density*t9_m23*exp(q6hep*t9_m13)
! THE DERIVATIVES OF THE CROSS SECTION FACTOR ARE NOT KNOWN AND ARE
!  TAKEN TO BE ZERO.  THE ONLY TERM FROM EQUATION 3.14 THAT SURVIVES
!  IS 5/(12*TAU).
         star%flux%neutrino_flux(i_nu_hep) = (1.0 + 0.067862*t9_p13)*star%cross_section_scale(17)* &
              star%flux%neutrino_flux(i_nu_hep)
! CALCULATE WEAK OR INTERMEDIATE SCREENING FOR HEP NEUTRINOS.
         zprdhe3p = 2.0
         z86he3p = 3.08687
         uwk = tm1*zprdhe3p
         if (uwk.le.star%ctrl%weak_screening_threshold) then
            utothe3p = uwk
         else
            uint = 0.38*xxl8*xtr*z86he3p/(ion_mean_weight_inverse*z58*z28)
            utothe3p = uint
         end if
! END OF CALCULATION OF SCREENING CORRECTION FOR HE3 + P REACTION.
         star%flux%neutrino_flux(i_nu_hep) = star%flux%neutrino_flux(i_nu_hep)*exp(utothe3p)
         star%flux%neutrino_flux(i_nu_hep) = star%flux%neutrino_flux(i_nu_hep)*hydrogen_fraction*he3_fraction/ &
              fourpiau2
! COMPUTE BE7MASSFRACTION. THIS IS NOT REQUIRED FOR THE NEUTRINO
!  FLUXES SINCE BE7 IS ALWAYS IN EQUILIBRIUM WITH THE SLOWER PRODUCTION
!  RATE OF HE3 + HE4.  HOWEVER, IT IS OF INTEREST IN SOME APPLICATIONS
!  TO KNOW THE BE7 MASS FRACTION, SO I COMPUTE IT HERE AND IT CAN BE
!  EXTRACTED WITH A COMMON STATEMENT IF DESIRED.
         star%engeb%be7_mass_fraction = eg(3)/(be7proton + be7electron)
! END OF NOVEMBER 6, 1990  ADDITION.
! FLUX OF BE7 NEUTRINOS.
         star%flux%neutrino_flux(i_nu_be7) = eg(3)*f1/fourpiau2
! FLUX OF B8 NEUTRINOS.
         star%flux%neutrino_flux(i_nu_b8) = eg(3)*f2/fourpiau2
! FLUX OF N13 NEUTRINOS.
         star%flux%neutrino_flux(i_nu_n13) = eg(4)/fourpiau2
! FLUX OF O15 NEUTRINOS.
         star%flux%neutrino_flux(i_nu_o15) = eg(6)/fourpiau2
! FLUX OF F17 NEUTRINOS.
         star%flux%neutrino_flux(i_nu_f17) = eg(7)/fourpiau2
! FLUX OF FICTIONAL HE3 + HE3 NEUTRINOS.
         star%flux%neutrino_flux(9) = eg(2)/fourpiau2
! FLUX OF FICTIONAL HE3 + HE4 NEUTRINOS.
         star%flux%neutrino_flux(10) = eg(3)/fourpiau2
! SET UNITS OF NEUTRINO FLUXES TO BE 10**10 PER CM^2 PER SEC PER GM AT THE
!  EARTH. MULTIPLY BY 10**-10.
!  IF THE VALUE FOR THIS SHELL IS NEGLIGIBLY SMALL, SET EQUAL TO ZERO.
         do k = 1,10
            star%flux%neutrino_flux(k) = (1.0e-10)*star%flux%neutrino_flux(k)
            flux_value = star%flux%neutrino_flux(k)
! KC 2025-05-30 CHANGED 1.E-50 TO 0.0 TO AVOID UNDERFLOW
            if (flux_value.le.0.0) then
              star%flux%neutrino_flux(k) = 0.0
            end if
         end do
! MHP 9/91 ENDIF INSERTED HERE.
         end if
      end if
! END OF NEUTRINO FLUX ROUTINE.
!
! ***MHP 5/91
!C CALCULATE NEUTRINO LOSSES FOR NEUTRINO-COOLED CORES OF EVOLVED STARS.
! 3/92 DBG Added option to use new (more sophisticated) neutrino loss
! routines.  See subroutine NEUTR for complete description.


      if (log_temperature.le.star%ctrl%tcut(5)) return


          carbon_fraction_total = c12_fraction+c13_fraction
          oxygen_fraction_total = o16_fraction+o18_fraction
          neutrino_temp=10.0**log_temperature
          neutrino_density=10.0**log_density


!**** Itoh 1996 Neutrino loss routines - Grant Newsham 9/06 *****


      if (star%ctrl%use_itoh_neutrino_loss) then



          call neutrino(neutrino_temp,neutrino_density,hydrogen_fraction, &
               helium_fraction,carbon_fraction_total,oxygen_fraction_total, &
               neutrino_loss_snu,neutrino_dlnq_dlnt,neutrino_dlnq_dlnd)


          star%engeb%neutrino_loss_rate = -neutrino_loss_snu
          neutrino_dlnq_dlnt = -neutrino_dlnq_dlnt*neutrino_temp/star%engeb%neutrino_loss_rate
          neutrino_dlnq_dlnd = -neutrino_dlnq_dlnd*neutrino_density/star%engeb%neutrino_loss_rate


          total_energy_gen_rate = total_energy_gen_rate + star%engeb%neutrino_loss_rate


          dlnepsilon_dlnrho = dlnepsilon_dlnrho + neutrino_dlnq_dlnd
          dlnepsilon_dlnt = dlnepsilon_dlnrho + neutrino_dlnq_dlnt


!****************************************************************


       else



!     THESE ARE OLD NEUTRINO LOSS ROUTINES


         el = t9/5.9302
         eli = 1.0/el
         ez = exp(cc13*ln10*(dd-9.0))*eli*(0.7937+0.2063*hydrogen_fraction)
         ez3 = ez**3
         emue = 0.5*(1.0+hydrogen_fraction)
         ex1 = 0.0
         if (t9.ge.0.2) then
!C PAIR NEUTRINOS
            el2 = el*el
            polx10=(1.+el2*(-13.04+el2*(133.5+el2*(1534.+el2*918.6))))
            polx11 = v1(1) + ez*(v1(2) + ez*v1(3))
            polx12 = ez3 + eli*(v1(4) + eli*(v1(5) + eli*v1(6)))
            ex1 = dexp(-ez*v1(7)-eli-eli-ln10*dd)*polx10*polx11/polx12
         end if
!C PHOTO NEUTRINOS
         polx21 = v2(1) + ez*(v2(2) + ez*v2(3))
         polx22 = ez3 + eli*(v2(4) + eli*(v2(5) + eli*v2(6)))
         ex2 = emue*el**5*dexp(-ez*v2(7))*polx21/polx22
!C PLASMA NEUTRINOS
         polx31 = v3(1) + ez*(v3(2) + ez*v3(3))
         polx32 = ez3 + eli*(v3(4) + eli*(v3(5) + eli*v3(6)))
         ex3 = emue**3*dexp(-ez*v3(7)+ln10*(dd+dd))*polx31/polx32
         star%engeb%neutrino_loss_rate = -(ex1 + ex2 + ex3)
         total_energy_gen_rate = total_energy_gen_rate + star%engeb%neutrino_loss_rate
         qetnx = 0.0
         qednx = 0.0
         if (t9.ge.0.2) then
! MHP 10/02 fixed column 72 problem
            qedn= ez*(v1(2)+2.*ez*v1(3))/polx11 - ez*v1(7) &
                  - 3.0d0*ez3/polx12
            qetn= -qedn+ eli*(v1(4)+eli*(2.*v1(5)+3.*eli*v1(6)))/polx12
            qetnx= el2*(-26.08+el2*(534.+el2*(9204.+el2*7348.8)))/polx10
            qetnx = (qetnx + qetn + eli+eli)*ex1
            qednx = (-1.0 +cc13*qedn)*ex1
         end if
         qedn= ez*(v2(2)+2.*ez*v2(3))/polx21 - ez*v2(7) - 3.*ez3/polx22
         qetn= -qedn+ eli*(v2(4)+eli*(2.*v2(5)+3.*eli*v2(6)))/polx22
         qetnx = qetnx + (5.0 + qetn)*ex2
         qednx = qednx +cc13*qedn*ex2
         qedn= ez*(v3(2)+2.*ez*v3(3))/polx31 - ez*v3(7) - 3.*ez3/polx32
         qetn=-qedn + eli*(v3(4)+eli*(2.*v3(5)+3.*eli*v3(6)))/polx32
         qetnx = qetnx + qetn*ex3
         qednx = qednx + (2.0 +cc13*qedn)*ex3



         dlnepsilon_dlnt = dlnepsilon_dlnt - qetnx
         dlnepsilon_dlnrho = dlnepsilon_dlnrho - qednx



      end if


end subroutine compute_neutrino_emission

end subroutine engeb


!----------------------------------------------------------------------
! liburn
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original liburn.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Determine lithium-6, lithium-7, and beryllium-9 burning in a
! non-rotating model (see liburn2.f90 for the simpler single-pass
! rotating-model variant). The burning rates depend on the local T and
! rho, the abundance, and (in convection zones) the mixing.
!
! In radiative regions burning is done implicitly, correcting the
! rates for the change in abundance as a function of time, by
! repeatedly halving the sub-step size and using a rational-function
! extrapolation (ratext) to converge on the zero-step-size limit.
!
! In the convective region the program iterates between burning and
! mixing in the same way, to correct the rates for the change in
! abundance as a function of time.
!
! 11/91 HR added to call.
!
! shell_mass genuinely needs intent(inout), not intent(in): it's
! temporarily perturbed for the CZ-base mass adjustment
! (shell_mass_save = shell_mass(cz_base_zone) ... shell_mass(cz_base_zone)
! = shell_mass_save) and the perturbed value is read by the
! immediately-following rate/abundance sums before being restored, so
! the write is real, not dead code -- unlike eos/eqstat.f90's
! metal_fraction fix, this can't just be narrowed away. Externally the
! net effect on the caller's array is still zero (every path that
! writes also restores before return), which is exactly why
! rotation/getw.f90 -- one of this routine's two callers -- had
! declared its own shell_mass intent(in) and needed widening to
! intent(inout) to match, once this routine gained an explicit
! interface (moved into net_lib.f90).
subroutine liburn(timestep, composition, radius, mass_coordinate, &
     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
      use star_info_lib, only: star, i_grad_ad, i_grad_rad, json
      use luout_lib
      use const_lib
      use numerics_lib
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: radius(json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones














      double precision :: li6_substep_depletion(json), &
           li7_substep_depletion(json), be9_substep_depletion(json)
      double precision :: light_element_save(3,json)
      integer :: substep_counts(11)
      double precision :: extrap_tol(3), extrap_y(3), extrap_err(3), &
           extrap_result(3)
      data substep_counts/2,4,6,8,12,16,24,32,48,64,96/
      data extrap_tol/1.0d-6,1.0d-6,1.0d-6/
      integer :: zone_idx, refine_idx, substep_idx, species_idx
      integer :: cz_base_zone, cz_base_zone_old, min_zone, max_zone
      double precision :: del_diff, del_diff_below, cz_base_frac
      double precision :: search_radius, shell_radius, delta_radius, &
           cz_base_radius
      double precision :: substep_frac, substep_dt
      logical :: converged
      double precision :: log_rate_li6, log_rate_li7, log_rate_be9
      double precision :: extrap_x
      double precision :: li6_cz_start, li7_cz_start, be9_cz_start, &
           cz_mass_start
      double precision :: log_rate_li6_cz_start, log_rate_li7_cz_start, &
           log_rate_be9_cz_start
      double precision :: shell_mass_save
      double precision :: li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end
      double precision :: log_rate_li6_cz_end, log_rate_li7_cz_end, &
           log_rate_be9_cz_end
      logical :: accretion_active
      double precision :: li6_accreted, li7_accreted, be9_accreted
      double precision :: cz_log_rate_li6, cz_log_rate_li7, cz_log_rate_be9
      double precision :: li6_depletion, li7_depletion, be9_depletion
      double precision :: li6_added, li7_added, be9_added
      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
           radiative_frac

! SAVE ORIGINAL ABUNDANCES.
      do zone_idx = 1,num_zones
         light_element_save(1,zone_idx) = composition(13,zone_idx)
         light_element_save(2,zone_idx) = composition(14,zone_idx)
         light_element_save(3,zone_idx) = composition(15,zone_idx)
      end do
! THE DEGREE OF LITHIUM BURNING IN A SURFACE CZ DEPENDS SENSITIVELY
! ON THE TEMPERATURE AT ITS BASE - SO ACCURATELY LOCATING IS IMPORTANT.
! DETERMINE THE TRUE LOCATION (FX) OF THE BASE OF THE CZ AT THE END
! OF THE TIMESTEP, AND THE LOCATION OF THE EDGE OF OVERSHOOT REGIONS
! IF APPLICABLE.
      if(env_cz_zone.gt.1.and.env_cz_zone.lt.num_zones)then
         if(star%job%rotation_active.and.star%job%instability_transport_active)then
            del_diff = star%mix_phys%del_adiabatic_mix(env_cz_zone) - &
                 star%mix_phys%del_radiative_mix(env_cz_zone)
            del_diff_below = star%mix_phys%del_adiabatic_mix(env_cz_zone-1) - &
                 star%mix_phys%del_radiative_mix(env_cz_zone-1)
         else
! EVALUATE DEL(AD) - DEL(RAD) AT THE LAST CONVECTIVE POINT AND THE ONE
! BELOW IT.
            del_diff = star%diag%del_grad(i_grad_ad,env_cz_zone)-star%diag%del_grad(i_grad_rad,env_cz_zone)
            del_diff_below = star%diag%del_grad(i_grad_ad,env_cz_zone-1)-star%diag%del_grad(i_grad_rad,env_cz_zone-1)
         endif
! USE LINEAR INTERPOLATION TO FIND THE DISTANCE OF THE TRUE LOCATION
! OF THE BASE FROM THE ZONE MIDPOINT. IF FX IS NEGATIVE,THEN THE TRUE
! BASE IS HIGHER; IF IT IS POSITIVE, THE TRUE BASE IS LOWER.
         cz_base_frac = max(-0.5d0,0.5d0-del_diff_below/ &
              (del_diff_below-del_diff))
         cz_base_frac = min(0.5d0,cz_base_frac)
         if(.not.star%job%envelope_overshoot_active)then
            cz_base_zone = env_cz_zone
            cz_base_zone_old = env_cz_zone_old
         else
! STARTING CZ DEPTH
            if(star%light_burn%cz_base_radius_prev.eq.0.0d0)then
               star%light_burn%cz_base_radius_prev = 0.5d0*(exp(ln10*star%prev%logR_start(env_cz_zone_old)) &
                        +exp(ln10*star%prev%logR_start(env_cz_zone_old-1)))
               search_radius = star%light_burn%cz_base_radius_prev - star%light_burn%pressure_scale_height_start
               do zone_idx = env_cz_zone_old-1,1,-1
                  shell_radius = exp(ln10*star%prev%logR_start(zone_idx))
                  if(shell_radius.lt.search_radius)then
                     cz_base_zone_old = zone_idx + 1
                     exit
                  endif
               end do
               if (zone_idx .lt. 1) cz_base_zone_old = 1
            else
               cz_base_zone_old = star%light_burn%envelope_cz_base_zone_prev
            endif
! ENDING CZ DEPTH : DETERMINE OVERSHOOT FROM TRUE CZ BASE.
            delta_radius = exp(ln10*radius(env_cz_zone))-exp(ln10*radius(env_cz_zone-1))
            cz_base_radius = 0.5d0*(exp(ln10*radius(env_cz_zone))+exp(ln10*radius(env_cz_zone-1))) &
                    -cz_base_frac*delta_radius
            star%light_burn%cz_base_radius_prev = cz_base_radius
            search_radius = cz_base_radius - star%light_burn%pressure_scale_height_end
            do zone_idx = env_cz_zone-1,1,-1
               shell_radius = exp(ln10*radius(zone_idx))
               if(shell_radius.lt.search_radius)then
                  cz_base_zone = zone_idx + 1
                  delta_radius = exp(ln10*star%prev%logR_start(zone_idx+1))-shell_radius
                  cz_base_frac = 0.5d0-((search_radius-shell_radius)/delta_radius)
                  cz_base_frac = max(-0.5d0,cz_base_frac)
                  cz_base_frac = min(0.5d0,cz_base_frac)
                  exit
               endif
            end do
            if (zone_idx .lt. 1) cz_base_zone = 1
         endif
      else
         cz_base_zone = env_cz_zone
         cz_base_zone_old = env_cz_zone_old
      endif
      star%light_burn%envelope_cz_base_zone_prev = cz_base_zone
! RADIATIVE INTERIOR.
      min_zone = min(cz_base_zone,cz_base_zone_old)
      max_zone = max(cz_base_zone,cz_base_zone_old)
      do zone_idx = 1,min_zone-1
         if(star%light_burn%rate_be9(zone_idx).le.1.0d-32 .or. star%light_burn%rate_be9_start(zone_idx).le.1.0d-32)exit
         if(composition(13,zone_idx).lt.1.0d-24.and.composition(14,zone_idx).lt.1.0d-24 .and.composition(15,zone_idx).lt.1.0d-24)cycle
         if(log_temperature(zone_idx).gt.7.0d0)then
            composition(13,zone_idx) = 0.0d0
            composition(14,zone_idx) = 0.0d0
            composition(15,zone_idx) = 0.0d0
            cycle
         endif
         do refine_idx = 1,11
! PERFORM BURNING IN substep_counts SMALLER TIMESTEPS.
!  DEFINE FRAC - THE FRACTION OF THE TOTAL TIMESTEP IN EACH SMALL ONE -
!  AND DELT - THE FRACTIONAL TIMESTEP IN YEARS.
            substep_frac = 1.0d0/dfloat(substep_counts(refine_idx))
            substep_dt = substep_frac*timestep
            converged = .false.
! RESTORE INITIAL ABUNDANCES.
            composition(13,zone_idx) = light_element_save(1,zone_idx)
            composition(14,zone_idx) = light_element_save(2,zone_idx)
            composition(15,zone_idx) = light_element_save(3,zone_idx)
! STORE STARTING REACTION RATES.
            log_rate_li6 = log(star%light_burn%rate_li6_start(zone_idx))-0.5d0*substep_frac* &
                   (log(star%light_burn%rate_li6(zone_idx))-log(star%light_burn%rate_li6_start(zone_idx)))
            log_rate_li7 = log(star%light_burn%rate_li7_start(zone_idx))-0.5d0*substep_frac* &
                   (log(star%light_burn%rate_li7(zone_idx))-log(star%light_burn%rate_li7_start(zone_idx)))
            log_rate_be9 = log(star%light_burn%rate_be9_start(zone_idx))-0.5d0*substep_frac* &
                   (log(star%light_burn%rate_be9(zone_idx))-log(star%light_burn%rate_be9_start(zone_idx)))
            do substep_idx = 1,substep_counts(refine_idx)
! INCREMENT THE REACTION RATES.
               log_rate_li6 = log_rate_li6+substep_frac* &
                    (log(star%light_burn%rate_li6(zone_idx))-log(star%light_burn%rate_li6_start(zone_idx)))
               log_rate_li7 = log_rate_li7+substep_frac* &
                    (log(star%light_burn%rate_li7(zone_idx))-log(star%light_burn%rate_li7_start(zone_idx)))
               log_rate_be9 = log_rate_be9+substep_frac* &
                    (log(star%light_burn%rate_be9(zone_idx))-log(star%light_burn%rate_be9_start(zone_idx)))
               li6_substep_depletion(zone_idx) = substep_dt*exp(log_rate_li6)
               li7_substep_depletion(zone_idx) = substep_dt*exp(log_rate_li7)
               be9_substep_depletion(zone_idx) = substep_dt*exp(log_rate_be9)
! THE REACTION RATES ARE OF THE FORM
!    DLI/DT = F(RHO,T,X)*LI =>DLN LI = DT*F(RHO,T,X)
! SOLVE FOR D LN (SPECIES)/DT AND ZERO OUT IF THE DEPLETION IS TOO HIGH.
               if(li6_substep_depletion(zone_idx).lt.3.0d1)then
                  composition(13,zone_idx) = composition(13,zone_idx)/ &
                       exp(li6_substep_depletion(zone_idx))
               else
                  composition(13,zone_idx) = 0.0d0
               endif
               if(li7_substep_depletion(zone_idx).lt.3.0d1)then
                  composition(14,zone_idx) = composition(14,zone_idx)/ &
                       exp(li7_substep_depletion(zone_idx))
               else
                  composition(14,zone_idx) = 0.0d0
               endif
               if(be9_substep_depletion(zone_idx).lt.3.0d1)then
                  composition(15,zone_idx) = composition(15,zone_idx)/ &
                       exp(be9_substep_depletion(zone_idx))
               else
                  composition(15,zone_idx) = 0.0d0
               endif
            end do
! STORE ABUNDANCES AS A FUNCTION OF TIMESTEP IN VECTOR YEST.
            extrap_y(1) = composition(13,zone_idx)
            extrap_y(2) = composition(14,zone_idx)
            extrap_y(3) = composition(15,zone_idx)
            extrap_x = substep_frac**2
! USE A RATIONAL POLYNOMIAL EXTRAPOLATOR TO EXTRAPOLATE TO THE SOLUTION FOR
! ZERO TIMESTEP.
            call ratext(refine_idx,extrap_x,extrap_y,extrap_result,extrap_err,3,7)
! CHECK IF THE SOLUTION HAS CONVERGED.
            if(refine_idx.gt.1)then
               converged = .true.
               do species_idx=1,3
                  if(extrap_y(species_idx).lt.1.0d-24)cycle
                  extrap_err(species_idx) = abs(extrap_err(species_idx)/extrap_y(species_idx))
                  if(extrap_err(species_idx).gt.extrap_tol(species_idx))converged=.false.
               end do
               if(converged)then
!                 WRITE(ISHORT,912)I,J,(YEXT(K2)K2=1,3)
!   912             FORMAT(1X,'CONVERGED',I5,' LEVEL ',I2,' LI6 ',1P,
!      *            E12.5,' LI7 ',E12.5,' BE9 ',E12.5)
                  exit
               endif
            endif
         end do
         if (refine_idx > 11) then
! IF THE PROGRAM GETS HERE THEN IT FAILED TO CONVERGE TO WITHIN
! THE SPECIFIED TOLERANCE IN THE MAXIMUM NUMBER OF ITERATIONS.
         write(short_file_unit,911)zone_idx,(extrap_err(species_idx),species_idx=1,3)
  911    format(1x,'***LIBURN CONVERGENCE FAILURE IN SHELL ',i4, &
         'ERRORS '/1p3e10.3)
         end if
! WRITE NEW ABUNDANCES AND EXIT.
         composition(13,zone_idx)=extrap_result(1)
         if(composition(13,zone_idx).lt.1.0d-24)composition(13,zone_idx)=0.0d0
         composition(14,zone_idx)=extrap_result(2)
         if(composition(14,zone_idx).lt.1.0d-24)composition(14,zone_idx)=0.0d0
         composition(15,zone_idx)=extrap_result(3)
         if(composition(15,zone_idx).lt.1.0d-24)composition(15,zone_idx)=0.0d0
      end do
! CONVECTION ZONE.
!
! SKIP IF WHOLE CZ IS BELOW THE BURNING THRESHOLD.
      if (star%light_burn%rate_be9_start(cz_base_zone_old).le.1.0d-32.or.star%light_burn%rate_be9(cz_base_zone).le.1.0d-32) then
         continue
         return
      end if
! FIND RATES AT THE BEGINNING OF THE TIMESTEP (USING THE DEPTH AT THE START).
      li6_cz_start = 0.0d0
      li7_cz_start = 0.0d0
      be9_cz_start = 0.0d0
      cz_mass_start = 0.0d0
      do zone_idx = cz_base_zone_old,num_zones
         li6_cz_start = li6_cz_start+composition(13,zone_idx)*shell_mass(zone_idx)
         li7_cz_start = li7_cz_start+composition(14,zone_idx)*shell_mass(zone_idx)
         be9_cz_start = be9_cz_start+composition(15,zone_idx)*shell_mass(zone_idx)
         cz_mass_start = cz_mass_start + shell_mass(zone_idx)
      end do
!    67 CONTINUE
      li6_cz_start = li6_cz_start/cz_mass_start
      li7_cz_start = li7_cz_start/cz_mass_start
      be9_cz_start = be9_cz_start/cz_mass_start
      if(star%light_burn%log_rate_li6_prev.le.0.0d0)then
! COMPUTE MASS-WEIGHTED AVERAGE RATES AT THE START OF THE STEP.
         log_rate_li6_cz_start = 0.0d0
         log_rate_li7_cz_start = 0.0d0
         log_rate_be9_cz_start = 0.0d0
         do zone_idx = cz_base_zone_old,num_zones
            log_rate_li6_cz_start = log_rate_li6_cz_start + star%light_burn%rate_li6_start(zone_idx)*shell_mass(zone_idx)
            log_rate_li7_cz_start = log_rate_li7_cz_start + star%light_burn%rate_li7_start(zone_idx)*shell_mass(zone_idx)
            log_rate_be9_cz_start = log_rate_be9_cz_start + star%light_burn%rate_be9_start(zone_idx)*shell_mass(zone_idx)
         end do
         log_rate_li6_cz_start = log(log_rate_li6_cz_start/cz_mass_start)
         log_rate_li7_cz_start = log(log_rate_li7_cz_start/cz_mass_start)
         log_rate_be9_cz_start = log(log_rate_be9_cz_start/cz_mass_start)
      else
! USE THE RATE FROM THE END OF THE PREVIOUS TIMESTEP.
         log_rate_li6_cz_start = star%light_burn%log_rate_li6_prev
         log_rate_li7_cz_start = star%light_burn%log_rate_li7_prev
         log_rate_be9_cz_start = star%light_burn%log_rate_be9_prev
      endif
! USE THE LOCATION OF THE TRUE EDGE OF THE CONVECTION ZONE (FX, FOUND AT
! BEGINNING OF SR) TO ADJUST THE BURNING RATE AND MASS OF THE BOTTOM POINT
! SUCH THAT IT INCLUDES THE ENTIRE C.Z.
      shell_mass_save = shell_mass(cz_base_zone)
      if(cz_base_zone.gt.1.and.cz_base_zone.lt.num_zones)then
         shell_mass(cz_base_zone) = shell_mass(cz_base_zone)+cz_base_frac* &
              (mass_coordinate(cz_base_zone)-mass_coordinate(cz_base_zone-1))
         star%light_burn%rate_li6(cz_base_zone) = star%light_burn%rate_li6(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%light_burn%rate_li6(cz_base_zone-1)-star%light_burn%rate_li6(cz_base_zone))
         star%light_burn%rate_li7(cz_base_zone) = star%light_burn%rate_li7(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%light_burn%rate_li7(cz_base_zone-1)-star%light_burn%rate_li7(cz_base_zone))
         star%light_burn%rate_be9(cz_base_zone) = star%light_burn%rate_be9(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%light_burn%rate_be9(cz_base_zone-1)-star%light_burn%rate_be9(cz_base_zone))
      endif
! FIND RATES AT THE END OF THE TIMESTEP (USING THE DEPTH AT THE END).
! ALSO STORE INITIAL ABUNDANCES(FLI60,FLI70,FBE90).
! FM IS THE TOTAL MASS IN THE CZ.
      li6_cz_end = 0.0d0
      li7_cz_end = 0.0d0
      be9_cz_end = 0.0d0
      cz_mass_end = 0.0d0
      log_rate_li6_cz_end = 0.0d0
      log_rate_li7_cz_end = 0.0d0
      log_rate_be9_cz_end = 0.0d0
      do zone_idx = cz_base_zone,num_zones
         log_rate_li6_cz_end = log_rate_li6_cz_end + star%light_burn%rate_li6(zone_idx)*shell_mass(zone_idx)
         log_rate_li7_cz_end = log_rate_li7_cz_end + star%light_burn%rate_li7(zone_idx)*shell_mass(zone_idx)
         log_rate_be9_cz_end = log_rate_be9_cz_end + star%light_burn%rate_be9(zone_idx)*shell_mass(zone_idx)
         li6_cz_end = li6_cz_end+composition(13,zone_idx)*shell_mass(zone_idx)
         li7_cz_end = li7_cz_end+composition(14,zone_idx)*shell_mass(zone_idx)
         be9_cz_end = be9_cz_end+composition(15,zone_idx)*shell_mass(zone_idx)
         cz_mass_end = cz_mass_end + shell_mass(zone_idx)
      end do
!    75 CONTINUE
      li6_cz_end = li6_cz_end/cz_mass_end
      li7_cz_end = li7_cz_end/cz_mass_end
      be9_cz_end = be9_cz_end/cz_mass_end
! MASS-WEIGHTED AVERAGE RATES AT THE END OF THE STEP.
      log_rate_li6_cz_end = log(log_rate_li6_cz_end/cz_mass_end)
      log_rate_li7_cz_end = log(log_rate_li7_cz_end/cz_mass_end)
      log_rate_be9_cz_end = log(log_rate_be9_cz_end/cz_mass_end)
! RESTORE ORIGINAL MASS TO BASE OF C.Z.
      shell_mass(cz_base_zone) = shell_mass_save
! FOR THE STARTING ABUNDANCE, USE THE MIXED ABUNDANCE FOR THE
! DEEPER CZ.
      if(cz_base_zone.lt.cz_base_zone_old)then
         li6_cz_start = li6_cz_end
         li7_cz_start = li7_cz_end
         be9_cz_start = be9_cz_end
      endif
      do refine_idx = 1,11
! PERFORM BURNING IN substep_counts SMALLER TIMESTEPS.
!  DEFINE FRAC - THE FRACTION OF THE TOTAL TIMESTEP IN EACH SMALL ONE -
!  AND DELT - THE FRACTIONAL TIMESTEP IN YEARS.
         substep_frac = 1.0d0/dfloat(substep_counts(refine_idx))
         substep_dt = substep_frac*timestep
         converged = .false.
! INITIALIZE ABUNDANCES.
         li6_cz_end = li6_cz_start
         li7_cz_end = li7_cz_start
         be9_cz_end = be9_cz_start
! INITIALIZE RATES.
         cz_log_rate_li6 = log_rate_li6_cz_start-0.5d0*substep_frac* &
              (log_rate_li6_cz_end-log_rate_li6_cz_start)
         cz_log_rate_li7 = log_rate_li7_cz_start-0.5d0*substep_frac* &
              (log_rate_li7_cz_end-log_rate_li7_cz_start)
         cz_log_rate_be9 = log_rate_be9_cz_start-0.5d0*substep_frac* &
              (log_rate_be9_cz_end-log_rate_be9_cz_start)
! BURN IN substep_counts SMALL TIMESTEPS.
! MHP 05/02 ADDED THE IMPACT OF MASS ACCRETION ON LIGHT
! ELEMENT ABUNDANCES.  INITIALIZE ARRAYS; SPECIES ARE
! ASSUMED TO EXPERIENCE BURNING FOR ON AVERAGE 1/2 OF
! THE STEP WHERE THEY ARE INITIALLY ACCRETED AND ARE
! THEN FULLY BURNED.
         if(star%job%use_mass_accretion.and.star%ctrl%mass_accretion_rate.gt.0.0d0)then
            accretion_active = .true.
            li6_accreted = 0.0d0
            li7_accreted = 0.0d0
            be9_accreted = 0.0d0
         else
            accretion_active = .false.
         endif
         do substep_idx = 1,substep_counts(refine_idx)
! INCREMENT RATES BY LINEAR INTERPOLATION IN THE LOG.
            cz_log_rate_li6 = cz_log_rate_li6 + substep_frac* &
                 (log_rate_li6_cz_end-log_rate_li6_cz_start)
            cz_log_rate_li7 = cz_log_rate_li7 + substep_frac* &
                 (log_rate_li7_cz_end-log_rate_li7_cz_start)
            cz_log_rate_be9 = cz_log_rate_be9 + substep_frac* &
                 (log_rate_be9_cz_end-log_rate_be9_cz_start)
            li6_depletion = substep_dt*exp(cz_log_rate_li6)
            li7_depletion = substep_dt*exp(cz_log_rate_li7)
            be9_depletion = substep_dt*exp(cz_log_rate_be9)
            if(li6_depletion.lt.3.0d1 .and. li6_cz_end.gt.1.0d-24)then
               li6_cz_end = li6_cz_end/exp(li6_depletion)
            else
               li6_cz_end = 0.0d0
            endif
            if(li6_cz_end.lt.1.0d-24)li6_cz_end = 0.0d0
            if(li7_depletion.lt.3.0d1 .and. li7_cz_end.gt.1.0d-24)then
               li7_cz_end = li7_cz_end/exp(li7_depletion)
            else
               li7_cz_end = 0.0d0
            endif
            if(li7_cz_end.lt.1.0d-24)li7_cz_end = 0.0d0
            if(be9_depletion.lt.3.0d1 .and. be9_cz_end.gt.1.0d-24)then
               be9_cz_end = be9_cz_end/exp(be9_depletion)
            else
               be9_cz_end = 0.0d0
            endif
            if(be9_cz_end.lt.1.0d-24)be9_cz_end = 0.0d0
! MHP 05/02 MASS ACCRETION.  THE NET LITHIUM ADDED IN
! A GIVEN STEP IS XXXADD WHILE THE LITHIUM ADDED IN PRIOR
! STEPS IS XXXA; THE TWO ARE ADDED AFTER BURNING.
! N.B. SINCE THIS IS A HALF-LIFE PROBLEM IT IS OK TO BURN
! THE COMPONENTS SEPARATELY.
            if(accretion_active)then
               li6_added = star%ctrl%accreted_composition(13)*substep_frac/exp(0.5d0*li6_depletion)
               li6_accreted = li6_accreted/exp(li6_depletion) + li6_added
               li7_added = star%ctrl%accreted_composition(14)*substep_frac/exp(0.5d0*li7_depletion)
               li7_accreted = li7_accreted/exp(li7_depletion) + li7_added
               be9_added = star%ctrl%accreted_composition(15)*substep_frac/exp(0.5d0*be9_depletion)
               be9_accreted = be9_accreted/exp(be9_depletion) + be9_added
            endif
         end do
! DO A MASS-WEIGHTED AVERAGE OF THE ORIGINAL LIGHT
! ELEMENT CONTENT AND THE NET AMOUNT ADDED.
! FMASSACC = DMDT*DT/ORIGINAL CZ MASS
! NOTE: THIS FORMULATION ASSUMES THAT ALL ACCRETED MATTER
! LANDS IN A CZ WITH THE END-OF-TIMESTEP DEPTH.
! NOT STRICTLY TRUE, BUT NOT A BAD APPROXIMATION EITHER.
         if(accretion_active)then
            write(*,913)li6_cz_start,li6_cz_end,li6_accreted,li7_cz_start, &
                 li7_cz_end,li7_accreted,star%light_burn%accreted_mass_fraction
 913        format(1p7e12.3)
            li6_cz_end = (li6_cz_end*cz_mass_end+li6_accreted*star%light_burn%accreted_mass_fraction)/ &
                 (cz_mass_end+star%light_burn%accreted_mass_fraction)
            li7_cz_end = (li7_cz_end*cz_mass_end+li7_accreted*star%light_burn%accreted_mass_fraction)/ &
                 (cz_mass_end+star%light_burn%accreted_mass_fraction)
            be9_cz_end = (be9_cz_end*cz_mass_end+be9_accreted*star%light_burn%accreted_mass_fraction)/ &
                 (cz_mass_end+star%light_burn%accreted_mass_fraction)
         endif
! STORE THE ABUNDANCE AS A FUNCTION OF THE TIME STEP.
! STORE ABUNDANCES AS A FUNCTION OF TIMESTEP IN VECTOR YEST.
         extrap_y(1) = li6_cz_end
         extrap_y(2) = li7_cz_end
         extrap_y(3) = be9_cz_end
         extrap_x = substep_frac**2
! USE A RATIONAL POLYNOMIAL EXTRAPOLATOR TO EXTRAPOLATE TO THE SOLUTION FOR
! ZERO TIMESTEP.
         call ratext(refine_idx,extrap_x,extrap_y,extrap_result,extrap_err,3,7)
! CHECK IF THE SOLUTION HAS CONVERGED.
         if(refine_idx.gt.1)then
            converged = .true.
            do species_idx=1,3
               if(extrap_y(species_idx).lt.1.0d-24)cycle
               extrap_err(species_idx) = abs(extrap_err(species_idx)/extrap_y(species_idx))
               if(extrap_err(species_idx).gt.extrap_tol(species_idx))converged=.false.
            end do
            if(converged)then
!              WRITE(ISHORT,912)JENV,J,(YEXT(K2),YERR(K2),K2=1,3)
               exit
            endif
         endif
      end do
      if (refine_idx .gt. 11) then
! IF THE PROGRAM GETS HERE THEN IT FAILED TO CONVERGE TO WITHIN
! THE SPECIFIED TOLERANCE IN THE MAXIMUM NUMBER OF ITERATIONS.
      write(short_file_unit,911)cz_base_zone,(extrap_err(species_idx),species_idx=1,3)
      end if
! WRITE NEW ABUNDANCES AND EXIT.
      li6_cz_end = extrap_result(1)
      if(li6_cz_end.lt.1.0d-24)li6_cz_end=0.0d0
      li7_cz_end = extrap_result(2)
      if(li7_cz_end.lt.1.0d-24)li7_cz_end=0.0d0
      be9_cz_end = extrap_result(3)
      if(be9_cz_end.lt.1.0d-24)be9_cz_end=0.0d0
      do zone_idx = max_zone,num_zones
         composition(13,zone_idx) = li6_cz_end
         composition(14,zone_idx) = li7_cz_end
         composition(15,zone_idx) = be9_cz_end
      end do
! STORE ENDING RATE FOR USE AT THE BEGINNING OF THE NEXT STEP.
      star%light_burn%log_rate_li6_prev = log_rate_li6_cz_end
      star%light_burn%log_rate_li7_prev = log_rate_li7_cz_end
      star%light_burn%log_rate_be9_prev = log_rate_be9_cz_end
! NOW SOLVE FOR ABUNDANCES IN THE REGION WHICH BEGAN CONVECTIVE AND
! ENDED RADIATIVE.
      if (cz_base_zone.le.cz_base_zone_old) then
         continue
         return
      end if
! FIND STARTING AND ENDING LOCATION IN MASS OF THE CZ BASE, AND
! ASSUME THAT THE FRACTION OF TIME SPENT IN THE CZ IS PROPORTIONAL TO
! THE LOCATION IN MASS (I.E. A POINT 1/3 FROM THE OLD TO THE NEW CZ
! BASE IN MASS SPENT 1/3 OF THE TIME IN THE CZ AND 2/3 OUT OF IT).
      if(cz_base_zone_old.gt.1)then
         mass_coord_beg = 0.5d0*(mass_coordinate(cz_base_zone_old)+mass_coordinate(cz_base_zone_old-1))
      else
         mass_coord_beg = 0.0d0
      endif
      mass_coord_end = 0.5d0*(mass_coordinate(cz_base_zone)+mass_coordinate(cz_base_zone-1))
      delta_mass = mass_coord_beg - mass_coord_end
      do zone_idx = cz_base_zone_old,cz_base_zone-1
! MHP 9/91 CHANGE TO AVOID DIVISION BY ZERO.
! SKIP IF SHELL TEMPERATURE DROPS BELOW BURNING THRESHOLD.
         if(star%light_burn%rate_be9(zone_idx).le.1.0d-32)exit
         radiative_frac = (mass_coordinate(zone_idx)-mass_coord_beg)/delta_mass
! USE FRAD*RADIATIVE RATE AND (1-FRAD)*CONVECTIVE RATE.
         li6_depletion = timestep*exp(radiative_frac*log(star%light_burn%rate_li6(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li6_cz_start)
         li7_depletion = timestep*exp(radiative_frac*log(star%light_burn%rate_li7(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li7_cz_start)
         be9_depletion = timestep*exp(radiative_frac*log(star%light_burn%rate_be9(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_be9_cz_start)
!***REMEMBER TO ADD FAILSAFES FOR LARGE DEPLETION***
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!          HCOMP(13,I) = HCOMP(13,I)/EXP(DLI6)
         call safedivexp(composition(13,zone_idx),li6_depletion)
         if(composition(13,zone_idx).lt.1.0d-24)composition(13,zone_idx)=0.0d0
!          HCOMP(14,I) = HCOMP(14,I)/EXP(DLI7)
         call safedivexp(composition(14,zone_idx),li7_depletion)
         if(composition(14,zone_idx).lt.1.0d-24)composition(14,zone_idx)=0.0d0
         composition(15,zone_idx) = composition(15,zone_idx)/exp(be9_depletion)
         if(composition(15,zone_idx).lt.1.0d-24)composition(15,zone_idx)=0.0d0
      end do
      return
end subroutine liburn


!----------------------------------------------------------------------
! liburn2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original liburn2.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Determine lithium-6, lithium-7, and beryllium-9 burning. The burning
! rates depend on the local T and rho, the abundance, and (in
! convection zones) the mixing.
!
! In radiative regions burning is done implicitly, using a single
! trapezoidal (midpoint-in-the-log) estimate of the burning rate over
! the timestep -- unlike liburn.f90, this variant does not iterate
! sub-step refinement with a rational-function extrapolation.
!
! In the convective region the program likewise uses a single
! time-averaged (midpoint-in-the-log) rate over the timestep, again
! without the iterative sub-step refinement used by liburn.f90.
!
! The degree of lithium burning in a surface CZ depends sensitively
! on the temperature at its base -- so accurately locating it is
! important. Determine the true location (cz_base_frac) of the base
! of the CZ at the end of the timestep, and the location of the edge
! of overshoot regions if applicable.
!
! 11/91 HR added to call.
!
! shell_mass genuinely needs intent(inout) -- same reasoning as
! liburn's identical fix above: the write (a save-then-restore pair
! around the CZ-base mass adjustment) is read by the following
! rate/abundance sums before being restored, so it's real, not dead
! code, even though the net effect on the caller's array is zero.
subroutine liburn2(timestep, composition, radius, mass_coordinate, &
     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
      use star_info_lib, only: star, i_grad_ad, i_grad_rad, json
      use luout_lib
      use const_lib
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: radius(json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones












      double precision :: li6_substep_depletion(json), &
           li7_substep_depletion(json), be9_substep_depletion(json)
      integer :: zone_idx, min_zone, max_zone
      double precision :: del_diff, del_diff_below, cz_base_frac
      double precision :: search_radius, shell_radius, delta_radius, &
           cz_base_radius
      integer :: cz_base_zone, cz_base_zone_old
      double precision :: log_rate_li6, log_rate_li7, log_rate_be9
      double precision :: li6_cz_start, li7_cz_start, be9_cz_start, &
           cz_mass_start
      double precision :: log_rate_li6_cz_start, log_rate_li7_cz_start, &
           log_rate_be9_cz_start
      double precision :: shell_mass_save
      double precision :: li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end
      double precision :: log_rate_li6_cz_end, log_rate_li7_cz_end, &
           log_rate_be9_cz_end
      double precision :: cz_log_rate_li6, cz_log_rate_li7, cz_log_rate_be9
      double precision :: li6_depletion, li7_depletion, be9_depletion
      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
           radiative_frac

! THE DEGREE OF LITHIUM BURNING IN A SURFACE CZ DEPENDS SENSITIVELY
! ON THE TEMPERATURE AT ITS BASE - SO ACCURATELY LOCATING IS IMPORTANT.
! DETERMINE THE TRUE LOCATION (FX) OF THE BASE OF THE CZ AT THE END
! OF THE TIMESTEP, AND THE LOCATION OF THE EDGE OF OVERSHOOT REGIONS
! IF APPLICABLE.
      if(env_cz_zone.gt.1.and.env_cz_zone.lt.num_zones)then
         if(star%job%rotation_active.and.star%job%instability_transport_active)then
            del_diff = star%mix_phys%del_adiabatic_mix(env_cz_zone) - &
                 star%mix_phys%del_radiative_mix(env_cz_zone)
            del_diff_below = star%mix_phys%del_adiabatic_mix(env_cz_zone-1) - &
                 star%mix_phys%del_radiative_mix(env_cz_zone-1)
         else
! EVALUATE DEL(AD) - DEL(RAD) AT THE LAST CONVECTIVE POINT AND THE ONE
! BELOW IT.
            del_diff = star%diag%del_grad(i_grad_ad,env_cz_zone)-star%diag%del_grad(i_grad_rad,env_cz_zone)
            del_diff_below = star%diag%del_grad(i_grad_ad,env_cz_zone-1)-star%diag%del_grad(i_grad_rad,env_cz_zone-1)
         endif
! USE LINEAR INTERPOLATION TO FIND THE DISTANCE OF THE TRUE LOCATION
! OF THE BASE FROM THE ZONE MIDPOINT. IF FX IS NEGATIVE,THEN THE TRUE
! BASE IS HIGHER; IF IT IS POSITIVE, THE TRUE BASE IS LOWER.
         cz_base_frac = max(-0.5d0,0.5d0-del_diff_below/ &
              (del_diff_below-del_diff))
         cz_base_frac = min(0.5d0,cz_base_frac)
         if(.not.star%job%envelope_overshoot_active)then
            cz_base_zone = env_cz_zone
            cz_base_zone_old = env_cz_zone_old
         else
! STARTING CZ DEPTH
            if(star%light_burn%cz_base_radius_prev.eq.0.0d0)then
               star%light_burn%cz_base_radius_prev = 0.5d0*(exp(ln10*star%prev%logR_start(env_cz_zone_old)) &
                        +exp(ln10*star%prev%logR_start(env_cz_zone_old-1)))
               search_radius = star%light_burn%cz_base_radius_prev - star%light_burn%pressure_scale_height_start
               do zone_idx = env_cz_zone_old-1,1,-1
                  shell_radius = exp(ln10*star%prev%logR_start(zone_idx))
                  if(shell_radius.lt.search_radius)then
                     cz_base_zone_old = zone_idx + 1
                     exit
                  endif
               end do
               if (zone_idx .lt. 1) cz_base_zone_old = 1
            else
               cz_base_zone_old = star%light_burn%envelope_cz_base_zone_prev
            endif
! ENDING CZ DEPTH : DETERMINE OVERSHOOT FROM TRUE CZ BASE.
            delta_radius = exp(ln10*radius(env_cz_zone))-exp(ln10*radius(env_cz_zone-1))
            cz_base_radius = 0.5d0*(exp(ln10*radius(env_cz_zone))+exp(ln10*radius(env_cz_zone-1))) &
                    -cz_base_frac*delta_radius
            star%light_burn%cz_base_radius_prev = cz_base_radius
            search_radius = cz_base_radius - star%light_burn%pressure_scale_height_end
            do zone_idx = env_cz_zone-1,1,-1
               shell_radius = exp(ln10*radius(zone_idx))
               if(shell_radius.lt.search_radius)then
                  cz_base_zone = zone_idx + 1
                  delta_radius = exp(ln10*star%prev%logR_start(zone_idx+1))-shell_radius
                  cz_base_frac = 0.5d0-((search_radius-shell_radius)/delta_radius)
                  cz_base_frac = max(-0.5d0,cz_base_frac)
                  cz_base_frac = min(0.5d0,cz_base_frac)
                  exit
               endif
            end do
            if (zone_idx .lt. 1) cz_base_zone = 1
         endif
      else
         cz_base_zone = env_cz_zone
         cz_base_zone_old = env_cz_zone_old
      endif
      star%light_burn%envelope_cz_base_zone_prev = cz_base_zone
! RADIATIVE INTERIOR.
      min_zone = min(cz_base_zone,cz_base_zone_old)
      max_zone = max(cz_base_zone,cz_base_zone_old)
      do zone_idx = 1,min_zone-1
         if(star%light_burn%rate_be9(zone_idx).le.1.0d-32 .or. star%light_burn%rate_be9_start(zone_idx).le.1.0d-32)cycle
         if(composition(13,zone_idx).lt.1.0d-24.and.composition(14,zone_idx).lt.1.0d-24 &
         .and.composition(15,zone_idx).lt.1.0d-24)cycle
         if(log_temperature(zone_idx).gt.7.0d0)then
            composition(13,zone_idx) = 0.0d0
            composition(14,zone_idx) = 0.0d0
            composition(15,zone_idx) = 0.0d0
            cycle
         endif
         log_rate_li6 = 0.5d0*(log(star%light_burn%rate_li6(zone_idx)) + log(star%light_burn%rate_li6_start(zone_idx)))
         log_rate_li7 = 0.5d0*(log(star%light_burn%rate_li7(zone_idx)) + log(star%light_burn%rate_li7_start(zone_idx)))
         log_rate_be9 = 0.5d0*(log(star%light_burn%rate_be9(zone_idx)) + log(star%light_burn%rate_be9_start(zone_idx)))
         li6_substep_depletion(zone_idx) = timestep*exp(log_rate_li6)
         li7_substep_depletion(zone_idx) = timestep*exp(log_rate_li7)
         be9_substep_depletion(zone_idx) = timestep*exp(log_rate_be9)
! THE REACTION RATES ARE OF THE FORM
!    DLI/DT = F(RHO,T,X)*LI =>DLN LI = DT*F(RHO,T,X)
! SOLVE FOR D LN (SPECIES)/DT AND ZERO OUT IF THE DEPLETION IS TOO HIGH.
         if(li6_substep_depletion(zone_idx).lt.3.0d1)then
            composition(13,zone_idx) = composition(13,zone_idx)/ &
                 exp(li6_substep_depletion(zone_idx))
         else
            composition(13,zone_idx) = 0.0d0
         endif
         if(li7_substep_depletion(zone_idx).lt.3.0d1)then
            composition(14,zone_idx) = composition(14,zone_idx)/ &
                 exp(li7_substep_depletion(zone_idx))
         else
            composition(14,zone_idx) = 0.0d0
         endif
         if(be9_substep_depletion(zone_idx).lt.3.0d1)then
            composition(15,zone_idx) = composition(15,zone_idx)/ &
                 exp(be9_substep_depletion(zone_idx))
         else
            composition(15,zone_idx) = 0.0d0
         endif
! WRITE NEW ABUNDANCES AND EXIT.
         if(composition(13,zone_idx).lt.1.0d-24)composition(13,zone_idx)=0.0d0
         if(composition(14,zone_idx).lt.1.0d-24)composition(14,zone_idx)=0.0d0
         if(composition(15,zone_idx).lt.1.0d-24)composition(15,zone_idx)=0.0d0
      end do
! CONVECTION ZONE.
!
! SKIP IF WHOLE CZ IS BELOW THE BURNING THRESHOLD.
      if (star%light_burn%rate_be9_start(cz_base_zone_old).le.1.0d-32.or.star%light_burn%rate_be9(cz_base_zone).le.1.0d-32) then
         continue
         return
      end if
! FIND RATES AT THE BEGINNING OF THE TIMESTEP (USING THE DEPTH AT THE START).
      li6_cz_start = 0.0d0
      li7_cz_start = 0.0d0
      be9_cz_start = 0.0d0
      cz_mass_start = 0.0d0
      do zone_idx = cz_base_zone_old,num_zones
         li6_cz_start = li6_cz_start+composition(13,zone_idx)*shell_mass(zone_idx)
         li7_cz_start = li7_cz_start+composition(14,zone_idx)*shell_mass(zone_idx)
         be9_cz_start = be9_cz_start+composition(15,zone_idx)*shell_mass(zone_idx)
         cz_mass_start = cz_mass_start + shell_mass(zone_idx)
      end do
!    67 CONTINUE
      li6_cz_start = li6_cz_start/cz_mass_start
      li7_cz_start = li7_cz_start/cz_mass_start
      be9_cz_start = be9_cz_start/cz_mass_start
      if(star%light_burn%log_rate_li6_prev.le.0.0d0)then
! COMPUTE MASS-WEIGHTED AVERAGE RATES AT THE START OF THE STEP.
         log_rate_li6_cz_start = 0.0d0
         log_rate_li7_cz_start = 0.0d0
         log_rate_be9_cz_start = 0.0d0
         do zone_idx = cz_base_zone_old,num_zones
            log_rate_li6_cz_start = log_rate_li6_cz_start + star%light_burn%rate_li6_start(zone_idx)*shell_mass(zone_idx)
            log_rate_li7_cz_start = log_rate_li7_cz_start + star%light_burn%rate_li7_start(zone_idx)*shell_mass(zone_idx)
            log_rate_be9_cz_start = log_rate_be9_cz_start + star%light_burn%rate_be9_start(zone_idx)*shell_mass(zone_idx)
         end do
         log_rate_li6_cz_start = log(log_rate_li6_cz_start/cz_mass_start)
         log_rate_li7_cz_start = log(log_rate_li7_cz_start/cz_mass_start)
         log_rate_be9_cz_start = log(log_rate_be9_cz_start/cz_mass_start)
      else
! USE THE RATE FROM THE END OF THE PREVIOUS TIMESTEP.
         log_rate_li6_cz_start = star%light_burn%log_rate_li6_prev
         log_rate_li7_cz_start = star%light_burn%log_rate_li7_prev
         log_rate_be9_cz_start = star%light_burn%log_rate_be9_prev
      endif
! USE THE LOCATION OF THE TRUE EDGE OF THE CONVECTION ZONE (FX, FOUND AT
! BEGINNING OF SR) TO ADJUST THE BURNING RATE AND MASS OF THE BOTTOM POINT
! SUCH THAT IT INCLUDES THE ENTIRE C.Z.
      shell_mass_save = shell_mass(cz_base_zone)
      if(cz_base_zone.gt.1.and.cz_base_zone.lt.num_zones)then
         shell_mass(cz_base_zone) = shell_mass(cz_base_zone)+cz_base_frac* &
              (mass_coordinate(cz_base_zone)-mass_coordinate(cz_base_zone-1))
         star%light_burn%rate_li6(cz_base_zone) = star%light_burn%rate_li6(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%light_burn%rate_li6(cz_base_zone-1)-star%light_burn%rate_li6(cz_base_zone))
         star%light_burn%rate_li7(cz_base_zone) = star%light_burn%rate_li7(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%light_burn%rate_li7(cz_base_zone-1)-star%light_burn%rate_li7(cz_base_zone))
         star%light_burn%rate_be9(cz_base_zone) = star%light_burn%rate_be9(cz_base_zone)+0.5d0*cz_base_frac* &
              (star%light_burn%rate_be9(cz_base_zone-1)-star%light_burn%rate_be9(cz_base_zone))
      endif
! FIND RATES AT THE END OF THE TIMESTEP (USING THE DEPTH AT THE END).
! ALSO STORE INITIAL ABUNDANCES(FLI60,FLI70,FBE90).
! FM IS THE TOTAL MASS IN THE CZ.
      li6_cz_end = 0.0d0
      li7_cz_end = 0.0d0
      be9_cz_end = 0.0d0
      cz_mass_end = 0.0d0
      log_rate_li6_cz_end = 0.0d0
      log_rate_li7_cz_end = 0.0d0
      log_rate_be9_cz_end = 0.0d0
      do zone_idx = cz_base_zone,num_zones
         log_rate_li6_cz_end = log_rate_li6_cz_end + star%light_burn%rate_li6(zone_idx)*shell_mass(zone_idx)
         log_rate_li7_cz_end = log_rate_li7_cz_end + star%light_burn%rate_li7(zone_idx)*shell_mass(zone_idx)
         log_rate_be9_cz_end = log_rate_be9_cz_end + star%light_burn%rate_be9(zone_idx)*shell_mass(zone_idx)
         li6_cz_end = li6_cz_end+composition(13,zone_idx)*shell_mass(zone_idx)
         li7_cz_end = li7_cz_end+composition(14,zone_idx)*shell_mass(zone_idx)
         be9_cz_end = be9_cz_end+composition(15,zone_idx)*shell_mass(zone_idx)
         cz_mass_end = cz_mass_end + shell_mass(zone_idx)
      end do
!    75 CONTINUE
      li6_cz_end = li6_cz_end/cz_mass_end
      li7_cz_end = li7_cz_end/cz_mass_end
      be9_cz_end = be9_cz_end/cz_mass_end
! MASS-WEIGHTED AVERAGE RATES AT THE END OF THE STEP.
      log_rate_li6_cz_end = log(log_rate_li6_cz_end/cz_mass_end)
      log_rate_li7_cz_end = log(log_rate_li7_cz_end/cz_mass_end)
      log_rate_be9_cz_end = log(log_rate_be9_cz_end/cz_mass_end)
! RESTORE ORIGINAL MASS TO BASE OF C.Z.
      shell_mass(cz_base_zone) = shell_mass_save
! FOR THE STARTING ABUNDANCE, USE THE MIXED ABUNDANCE FOR THE
! DEEPER CZ.
      if(cz_base_zone.lt.cz_base_zone_old)then
         li6_cz_start = li6_cz_end
         li7_cz_start = li7_cz_end
         be9_cz_start = be9_cz_end
      endif
! INITIALIZE ABUNDANCES.
         li6_cz_end = li6_cz_start
         li7_cz_end = li7_cz_start
         be9_cz_end = be9_cz_start
! TIME-AVERAGED RATES USING LINEAR INTERPOLATION IN THE LOG.
         cz_log_rate_li6 = 0.5d0*(log_rate_li6_cz_end+log_rate_li6_cz_start)
         cz_log_rate_li7 = 0.5d0*(log_rate_li7_cz_end+log_rate_li7_cz_start)
         cz_log_rate_be9 = 0.5d0*(log_rate_be9_cz_end+log_rate_be9_cz_start)
         li6_depletion = timestep*exp(cz_log_rate_li6)
         li7_depletion = timestep*exp(cz_log_rate_li7)
         be9_depletion = timestep*exp(cz_log_rate_be9)
         if(li6_depletion.lt.3.0d1 .and. li6_cz_end.gt.1.0d-24)then
            li6_cz_end = li6_cz_end/exp(li6_depletion)
         else
            li6_cz_end = 0.0d0
         endif
         if(li6_cz_end.lt.1.0d-24)li6_cz_end = 0.0d0
         if(li7_depletion.lt.3.0d1 .and. li7_cz_end.gt.1.0d-24)then
            li7_cz_end = li7_cz_end/exp(li7_depletion)
         else
            li7_cz_end = 0.0d0
         endif
         if(li7_cz_end.lt.1.0d-24)li7_cz_end = 0.0d0
         if(be9_depletion.lt.3.0d1 .and. be9_cz_end.gt.1.0d-24)then
            be9_cz_end = be9_cz_end/exp(be9_depletion)
         else
            be9_cz_end = 0.0d0
         endif
         if(be9_cz_end.lt.1.0d-24)be9_cz_end = 0.0d0
      do zone_idx = max_zone,num_zones
         composition(13,zone_idx) = li6_cz_end
         composition(14,zone_idx) = li7_cz_end
         composition(15,zone_idx) = be9_cz_end
      end do
! STORE ENDING RATE FOR USE AT THE BEGINNING OF THE NEXT STEP.
      star%light_burn%log_rate_li6_prev = log_rate_li6_cz_end
      star%light_burn%log_rate_li7_prev = log_rate_li7_cz_end
      star%light_burn%log_rate_be9_prev = log_rate_be9_cz_end
! NOW SOLVE FOR ABUNDANCES IN THE REGION WHICH BEGAN CONVECTIVE AND
! ENDED RADIATIVE.
      if (cz_base_zone.le.cz_base_zone_old) then
         continue
         return
      end if
! FIND STARTING AND ENDING LOCATION IN MASS OF THE CZ BASE, AND
! ASSUME THAT THE FRACTION OF TIME SPENT IN THE CZ IS PROPORTIONAL TO
! THE LOCATION IN MASS (I.E. A POINT 1/3 FROM THE OLD TO THE NEW CZ
! BASE IN MASS SPENT 1/3 OF THE TIME IN THE CZ AND 2/3 OUT OF IT).
      if(cz_base_zone_old.gt.1)then
         mass_coord_beg = 0.5d0*(mass_coordinate(cz_base_zone_old)+mass_coordinate(cz_base_zone_old-1))
      else
         mass_coord_beg = 0.0d0
      endif
      mass_coord_end = 0.5d0*(mass_coordinate(cz_base_zone)+mass_coordinate(cz_base_zone-1))
      delta_mass = mass_coord_beg - mass_coord_end
      do zone_idx = cz_base_zone_old,cz_base_zone-1
! MHP 9/91 CHANGE TO AVOID DIVISION BY ZERO.
! SKIP IF SHELL TEMPERATURE DROPS BELOW BURNING THRESHOLD.
         if(star%light_burn%rate_be9(zone_idx).le.1.0d-32)exit
         radiative_frac = (mass_coordinate(zone_idx)-mass_coord_beg)/delta_mass
! USE FRAD*RADIATIVE RATE AND (1-FRAD)*CONVECTIVE RATE.
         li6_depletion = timestep*exp(radiative_frac*log(star%light_burn%rate_li6(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li6_cz_start)
         li7_depletion = timestep*exp(radiative_frac*log(star%light_burn%rate_li7(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li7_cz_start)
         be9_depletion = timestep*exp(radiative_frac*log(star%light_burn%rate_be9(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_be9_cz_start)
!***REMEMBER TO ADD FAILSAFES FOR LARGE DEPLETION***
         composition(13,zone_idx) = composition(13,zone_idx)/exp(li6_depletion)
         if(composition(13,zone_idx).lt.1.0d-24)composition(13,zone_idx)=0.0d0
         composition(14,zone_idx) = composition(14,zone_idx)/exp(li7_depletion)
         if(composition(14,zone_idx).lt.1.0d-24)composition(14,zone_idx)=0.0d0
         composition(15,zone_idx) = composition(15,zone_idx)/exp(be9_depletion)
         if(composition(15,zone_idx).lt.1.0d-24)composition(15,zone_idx)=0.0d0
      end do
      return
end subroutine liburn2


!----------------------------------------------------------------------
! lirate88
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lirate88.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine computes the burning rates for Li6, Li7, Be9 at the
! beginning and end of a timestep and stores them in common blocks
! newrat and oldrat.
!
! Burning rates from Caughlin and Fowler (1988).
subroutine lirate88(composition, log_density, log_temperature, num_zones, &
     use_current_model)
      use star_info_lib, only: star, json
      use const_lib
      implicit none

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_density(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: num_zones, use_current_model





! G Somers END

      double precision :: tlim
      data tlim/6.0d0/
      integer :: zone_idx, tail_idx
      double precision :: rhox, t9, t913, t923, t943, t953, t932, t934
      double precision :: fli6, fli7, fbe91, fbe92, fbe93, fx, ex, fsbe9
      double precision :: t9a, c56

      do zone_idx = 1,num_zones
         if(log_temperature(zone_idx).lt.tlim.and.star%prev%logT_start(zone_idx).lt.tlim)exit
         if(use_current_model.eq.1)then
            rhox = exp(ln10*log_density(zone_idx))*composition(1,zone_idx)
            t9=exp(ln10*(log_temperature(zone_idx)-9.0d0))
         else
            rhox = exp(ln10*star%prev%logRho_start(zone_idx))*star%prev%xa_start(1,zone_idx)
            t9=exp(ln10*(star%prev%logT_start(zone_idx)-9.0d0))
         endif
         t913=t9**cc13
         t923=t913*t913
         t943=t923*t923
! MHP 10/91 ADDED DEFINITION
         t953 = t943*t913
         t932=t9**1.5d0
         t934=t9**7.5d-1
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
         fli7=1.096d9/t923*exp(-8.472d0/t913)-4.830d8*t9a**c56/ &
         t932*exp(-8.472d0/t9a**cc13)+1.06d10/t932*exp(-30.442d0/t9)
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
         if(use_current_model.eq.1)then
            star%light_burn%rate_li6(zone_idx) = fli6*rhox
            star%light_burn%rate_li7(zone_idx) = fli7*rhox
            star%light_burn%rate_be9(zone_idx) = fsbe9*rhox
         else
            star%light_burn%rate_li6_start(zone_idx) = fli6*rhox
            star%light_burn%rate_li7_start(zone_idx) = fli7*rhox
            star%light_burn%rate_be9_start(zone_idx) = fsbe9*rhox
         endif
      end do
      do tail_idx = zone_idx,num_zones
         star%light_burn%rate_li6(tail_idx)=0.0d0
         star%light_burn%rate_li7(tail_idx)=0.0d0
         star%light_burn%rate_be9(tail_idx)=0.0d0
         star%light_burn%rate_li6_start(tail_idx)=0.0d0
         star%light_burn%rate_li7_start(tail_idx)=0.0d0
         star%light_burn%rate_be9_start(tail_idx)=0.0d0
      end do
      return
end subroutine lirate88


end module burn_lib
