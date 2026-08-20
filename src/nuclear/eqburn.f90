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

      implicit none
      integer, parameter :: json = 5000

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

! common/ctlim/: only tcut is used here. Naming matches eqstat2.f90.
      double precision :: atime(14), tcut(5), saha_log10t_cutoff, &
           tenv0, tenv1, tenv, tgcut
      common/ctlim/ atime, tcut, saha_log10t_cutoff, tenv0, tenv1, tenv, tgcut

! common/oldmod/: previous-timestep model snapshot; only
! old_composition is used here. Naming matches dburn.f90.
      double precision :: old_pressure(json), old_temperature(json), &
           old_radius(json), old_luminosity(json), old_density(json), &
           old_composition(15,json), old_shell_mass(json), old_teff
      logical :: old_convective_flag(json), old_cz_flag(json)
      integer :: old_num_zones
      common/oldmod/ old_pressure, old_temperature, old_radius, &
           old_luminosity, old_density, old_composition, old_shell_mass, &
           old_convective_flag, old_cz_flag, old_teff, old_num_zones

      double precision :: zone_avg_abundance(11)
      save

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
      if (shell_temperature(zone_begin).lt.tcut(1)) return
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
                    old_composition(species_idx,zone_idx)*shell_mass(zone_idx)
            end do
         end do
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = &
                 zone_avg_abundance(species_idx)/total_shell_mass
         end do
      else
         do species_idx = 1, 11
            zone_avg_abundance(species_idx) = &
                 old_composition(species_idx,zone_begin)
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
           shell_temperature(zone_begin).gt.tcut(2)) then
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
      else
         dx_dt = 0.0d0
         goto 100
      end if
      if (shell_temperature(zone_begin).gt.tcut(3)) then
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
  100 continue
!     HELIUM BURNING REACTIONS.
      if (helium_fraction.gt.1.0d-10 .and. &
           shell_temperature(zone_begin).gt.tcut(4)) then
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
