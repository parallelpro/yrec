!----------------------------------------------------------------------
! test_net
!----------------------------------------------------------------------
! Standalone net-domain test program (2026), companion to
! eos/kap/atm's (see eos/test/test_eos.f90 for the initialization
! conventions). net/ is a function of (log T, log rho, composition)
! once its settings are installed, so each piece is tested on its own.
! It prints values for byte-comparison against expected_test_net.out,
! and checks that print ok or FAIL (each confirmed to fail under a
! mutation of the code it covers):
!   * set_net_settings installs net_settings_from_star's record (the
!     table self-check passes); the cross-section scales of
!     map_user_inputs (Solar Fusion II inputs);
!   * compute_network_rates on a (log T, log rho) grid, and
!     compute_reaction_rates' analytic derivatives against centred
!     differences (5e-3: the screening derivatives hold the plasma fixed);
!   * sneut (Itoh et al. 1996) and azbar; deuterium_rate;
!   * the reaction table: baryon balance, and each row's Q against the
!     atomic mass excesses (the difference is the row's neutrinos);
!   * burn_implicitly against backward Euler's closed form for one
!     reaction at a time, its mass conservation with all of them, and its
!     ierr when it runs out of iterations;
!   * the trace integrator (surviving_fraction, accreted_survival)
!     against closed forms and an independent quadrature;
!   * compute_nuclear_energy at three reference points, and its
!     derivatives against centred differences; the thermal neutrino
!     losses alone (C/O, nothing burning), both fits, with their
!     derivatives; deuterium burning outside and inside the surface
!     convection zone (the cap binding), with their derivatives; nothing
!     below the burning cutoff;
!   * strong screening: the regime, its effect, the rates' derivatives;
!   * compute_abundance_rates (the network's dX/dt) on a hydrogen- and
!     a helium-burning zone, burn_deuterium_region (shell, region,
!     accreting surface region, below the threshold), compute_li_be_rates,
!     the deuterium-rate cap (exactly N_A/(m_d tau)), and the solar
!     neutrino fluxes;
!   * the configuration: the table's self-check, the legacy sstandard
!     map, and the scale reported as without effect.
program test_net
      use run_setup, only: set_up_constants_and_tables
      use user_inputs, only: map_user_inputs
      use controls_buffer
      use io_units
      use star_info_lib, only: star_info, i_h1, i_h2, i_he3, i_he4, i_metals, &
           i_c12, i_c13, i_n14, i_o16, i_o18
      use net_lib, only: compute_network_rates, deuterium_rate, &
           num_rate_rows, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, &
           r_po16, r_pp, burn_implicitly, n_network, s_h1, s_he3, s_he4, &
           surviving_fraction, accreted_survival, linear_in_time, &
           linear_in_log, compute_nuclear_energy, nuclear_energy_result, &
           net_settings_t, set_net_settings, compute_abundance_rates, &
           burn_deuterium_region, surface_accretion, compute_li_be_rates, &
           network_species_slot, chain_pp1, chain_cno, chain_triple_alpha
! the domain's own test: it also reaches net/'s internal modules
      use screening, only: num_isotopes, iso_h1, iso_h2, iso_he3, iso_he4, &
           iso_c12, iso_c13, iso_n14, iso_o16, iso_o18
      use reaction_rates, only: nuclear_rates, zone_conditions, &
           compute_zone_conditions
      use net_configuration, only: num_reactions, num_extra_scales, r_3alpha, &
           r_ac12, r_ac13, r_an14, rv_frac_be7_ecap, net_settings, num_sstandard, &
           cross_section_scales_from_sstandard
      use network, only: s_c12, s_o16, write_network_table, compute_reaction_rates, &
           network_table_errors, inert_scale_count
      use thermal_neutrinos, only: sneut, azbar
      use net_settings_setup, only: net_settings_from_star
      implicit none
! the test's star (star_info_lib declares no instance); static, as it
! is too large for the stack
      type(star_info), save :: star

      character(len=256) :: yrec_input
      character(len=256) :: fermi_path, dummy_path
      character(len=256) :: dummy_paths7(7)
      double precision :: laol_work(12)
      integer :: setups_ierr, ipt, i

      integer, parameter :: npts = 4
      double precision :: grid_logt(npts), grid_logd(npts)
      data grid_logt /6.90d0, 7.10d0, 7.30d0, 7.60d0/
      data grid_logd /1.0d0,  1.8d0,  2.2d0,  3.0d0/

! rate_vec: the one-zone rate vector of rates (slots r_*, rv_*)
      double precision :: rate_vec(num_rate_rows)
! the composition of the rates check, star%xa layout
      double precision :: test_xa(15)
! derivative check: analytic vs centred differences of the rates
      type(nuclear_rates) :: nr0, nr_tp, nr_tm, nr_dp, nr_dm
      double precision :: fd_mass(num_isotopes), fd_err_t(num_reactions), &
           fd_err_d(num_reactions), fd_num
      integer, parameter :: nfd = 3
      double precision :: fd_logt(nfd), fd_logd(nfd)
      data fd_logt /7.10d0, 7.30d0, 8.00d0/
      data fd_logd /1.80d0, 2.20d0, 3.00d0/
      double precision, parameter :: fd_h = 1.0d-4
      integer :: k
      double precision :: snu, dsnudt, dsnudd, dsnuda, dsnudz
      double precision :: xmass(3), aion(3), zion(3), ymass(3), &
           abar, zbar
! the eps of the last point energy_derivative_errors evaluated
      double precision :: last_eps

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_net: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      fermi_path = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      dummy_path = ""
      dummy_paths7 = ""
      laol_work = 0.0d0

! the run log's unit, as io/run_files.f90 assigns it
      run_log_unit = 20
      star%ctrl%fermi_unit = 15
      open(run_log_unit, file="test_net.short", status="replace")

! everything gated off; only the Fermi table is a hard requirement
      star%ctrl%use_mhd_eos = .false.
      star%ctrl%use_scv_eos = .false.
      star%ctrl%use_opal95_eos = .false.
      star%ctrl%use_opal2001_eos = .false.
      star%ctrl%use_opal2006_eos = .false.
      star%ctrl%use_opal95_tables = .false.
      star%ctrl%use_opal92_tables = .false.
      star%ctrl%use_laol89_tables = .false.
      star%ctrl%use_alex06_tables = .false.
      star%ctrl%use_alex95_tables = .false.
      star%ctrl%use_kurucz90_tables = .false.
      star%use_two_z_tables = .false.
      star%ctrl%use_pure_z_table = .false.
      star%ctrl%use_conductive_opacity = .false.
      star%job%atm_choice = 0

      star%job%mixture_weights = laol_work
      star%job%alex06_table_path = dummy_path
      star%job%allard_table_path = dummy_path
      star%job%atm_table_path = dummy_path
      star%job%fermi_table_path = fermi_path
      star%job%kurucz_table_path = dummy_path
      star%job%kurucz_table2_path = dummy_path
      star%job%laol_table_path = dummy_path
      star%job%laol_table2_path = dummy_path
      star%job%opal95_table_path = dummy_path
      star%job%opal92_table_path = dummy_path
      star%job%zams_a_table_path = dummy_path
      star%job%zams_b_table_path = dummy_path
      star%job%zams_c_table_path = dummy_path
      star%job%centre1_table_path = dummy_path
      star%job%centre2_table_path = dummy_path
      star%job%centre3_table_path = dummy_path
      star%job%centre4_table_path = dummy_path
      star%job%centre5_table_path = dummy_path
      star%job%opal92_table2_path = dummy_path
      star%job%pure_z_table_path = dummy_path
      star%job%scv_h_table_path = dummy_path
      star%job%scv_he_table_path = dummy_path
      star%job%scv_z_table_path = dummy_path
      star%job%alex95_table_paths = dummy_paths7
      call set_up_constants_and_tables(star, setups_ierr)
      if (setups_ierr /= 0) then
         write(*,'(a)') "test_net: FAIL (setups error)"
         stop 1
      end if

! The nuclear cross sections: a run's controls read fills the S-factor
! controls with the Solar Fusion II values (Adelberger et al. 2011;
! defaults/physics.defaults) before map_user_inputs turns them into the
! star's scales. This test bypasses the read, so it sets the same values
! in the controls buffer, with the new-rates path on.
      use_new_nuclear_rates = .true.
      s0_pp = 4.01d-22
      s0_he3he3 = 5.21d3
      s0_he3he4 = 5.6d-1
      s0_p_c12 = 1.34d0
      s0_p_c13 = 7.6d0
      s0_p_n14 = 1.66d0
      s0_p_o16 = 1.06d1
      s0_be7_electron = 1.7709d-10
      s0_be7_p = 0.0208d0
      s0_n15_p_c12_branch = 7.3d4
      s0_n15_p_o16_branch = 3.6d1
      s0p_pp = 4.49d-24
      s0p_he3he3 = -4.9d0
      s0p_he3he4 = -3.6d-4
      s0p_p_c12 = 2.6d-3
      s0p_p_c13 = -7.83d-3
      s0p_p_n14 = -3.3d-3
      s0p_p_o16 = -5.4d-2
      s0p_be7_p = -3.12d-5
      s0pp_p_c12 = 8.3d-5
      s0pp_p_c13 = 7.29d-4
      s0pp_p_o16 = 0.0d0
      s0pp_be7_p = -2.288d-7

! cross-section scales via the real remap, then the network's settings
      call map_user_inputs(star)
      call set_net_settings(net_settings_from_star(star), setups_ierr)
      call report("set_net_settings: the table checks out, the settings installed", &
           setups_ierr == 0)
      write(*,'(a)') "# test_net: cross_section_scale (reaction slots), then " // &
           "extra_cross_section_scale (pep, Be7+e, Be7+p, hep), from map_user_inputs " // &
           "(newcross defaults)"
      write(*,'(4(1pe20.12))') (star%cross_section_scale(i), i = 1, num_reactions)
      write(*,'(4(1pe20.12))') (star%extra_cross_section_scale(i), i = 1, num_extra_scales)

! reaction rates over a (logT, logRho) grid at a solar-ish mixture
      write(*,'(a)') "# test_net: rates, X=0.70 Y3=3e-5 C12=3.5e-3 " // &
           "C13=4e-5 N14=1e-3 O16=8e-3 O18=2e-5"
      test_xa = 0.0d0
      test_xa(1) = 0.70d0
      test_xa(2) = 0.28d0
      test_xa(4) = 3.0d-5
      test_xa(5) = 3.5d-3
      test_xa(6) = 4.0d-5
      test_xa(7) = 1.0d-3
      test_xa(9) = 8.0d-3
      test_xa(11) = 2.0d-5
      do ipt = 1, npts
         call compute_network_rates(grid_logd(ipt), grid_logt(ipt), test_xa, rate_vec)
         write(*,'(a,i2)') "rates ", ipt
         write(*,'(4(1pe20.12))') rate_vec(r_pp), rate_vec(r_he3he3), &
              rate_vec(r_he3he4), rate_vec(r_pc12)
         write(*,'(4(1pe20.12))') rate_vec(r_pc13), rate_vec(r_pn14), &
              rate_vec(r_po16), rate_vec(r_ac13)
         write(*,'(4(1pe20.12))') rate_vec(r_ac12), rate_vec(r_an14), &
              rate_vec(r_3alpha), rate_vec(rv_frac_be7_ecap)
      end do

! Itoh neutrino losses
      write(*,'(a)') "# test_net: sneut (Itoh et al. 1996)"
      do ipt = 1, npts
         call sneut(10.0d0**grid_logt(ipt), 10.0d0**grid_logd(ipt), &
              1.30d0, 1.10d0, snu, dsnudt, dsnudd, dsnuda, dsnudz)
         write(*,'(a,i2,3(1pe20.12))') "snu ", ipt, snu, dsnudt, dsnudd
      end do

! azbar bookkeeping on a 3-species H/He4/C12 mixture
      xmass = (/0.70d0, 0.28d0, 0.02d0/)
      aion  = (/1.0d0, 4.0d0, 12.0d0/)
      zion  = (/1.0d0, 2.0d0, 6.0d0/)
      call azbar(xmass, aion, zion, 3, ymass, abar, zbar)
      write(*,'(a)') "# test_net: azbar H/He4/C12"
      write(*,'(5(1pe20.12))') abar, zbar, ymass(1), ymass(2), ymass(3)

! the deuterium burning rate of a shell (per s per amu) outside the
! surface convection zone, so uncapped
      write(*,'(a)') "# test_net: deuterium_rate, log T 6.2, log rho 1.5, outside the surface CZ"
      write(*,'(1pe20.12)') deuterium_rate(1.5d0, 6.2d0, 0.70d0, .false., 0.0d0)

! analytic dln(rate)/dlnT and dln(rate)/dln(rho) of compute_reaction_rates
! against centred differences (step 1e-4 dex), worst relative error per
! reaction slot over three (log T, log rho) points; slots that are zero
! everywhere print 0. Guards the hand-derived Jacobian terms; the
! screening derivatives are approximate (their plasma terms are held
! fixed), so the check allows 2e-4 in d/d ln T and 2e-3 in d/d ln rho.
      write(*,'(a)') "# test_net: rate derivatives, max |analytic/numerical - 1| (dlnT, dlnrho)"
      fd_mass = 0.0d0
      fd_mass(iso_h1) = 0.70d0
      fd_mass(iso_he3) = 3.0d-5
      fd_mass(iso_he4) = 0.28d0
      fd_mass(iso_c12) = 3.5d-3
      fd_mass(iso_c13) = 4.0d-5
      fd_mass(iso_n14) = 1.0d-3
      fd_mass(iso_o16) = 8.0d-3
      fd_mass(iso_o18) = 2.0d-5
      fd_err_t = 0.0d0
      fd_err_d = 0.0d0
      do ipt = 1, nfd
         call compute_reaction_rates(fd_logd(ipt), fd_logt(ipt), fd_mass, nr0)
         call compute_reaction_rates(fd_logd(ipt), fd_logt(ipt)+fd_h, fd_mass, nr_tp)
         call compute_reaction_rates(fd_logd(ipt), fd_logt(ipt)-fd_h, fd_mass, nr_tm)
         call compute_reaction_rates(fd_logd(ipt)+fd_h, fd_logt(ipt), fd_mass, nr_dp)
         call compute_reaction_rates(fd_logd(ipt)-fd_h, fd_logt(ipt), fd_mass, nr_dm)
         do k = 1, num_reactions
            if (nr0%rate(k) <= 0.0d0 .or. nr_tp%rate(k) <= 0.0d0 .or. nr_tm%rate(k) <= 0.0d0) cycle
            fd_num = (log(nr_tp%rate(k)) - log(nr_tm%rate(k)))/(2.0d0*fd_h*log(10.0d0))
            fd_err_t(k) = max(fd_err_t(k), abs(nr0%dlnrate_dlnt(k)/fd_num - 1.0d0))
            fd_num = (log(nr_dp%rate(k)) - log(nr_dm%rate(k)))/(2.0d0*fd_h*log(10.0d0))
            fd_err_d(k) = max(fd_err_d(k), abs(nr0%dlnrate_dlnrho(k)/fd_num - 1.0d0))
         end do
      end do
      do k = 1, num_reactions
         write(*,'(i3, 2es10.2)') k, fd_err_t(k), fd_err_d(k)
      end do
      call report("rate derivatives vs centred differences (2e-4 in T, 2e-3 in rho)", &
           maxval(fd_err_t) < 2.0d-4 .and. maxval(fd_err_d) < 2.0d-3)

! the burning network's reaction table and its baryon-balance check
      write(*,'(a)') "# test_net: network table"
      call write_network_table(6)

      call check_burn_implicitly()
      call check_burn_nonconvergence()
      call check_trace_integrator()
      call check_nuclear_energy()
      call check_thermal_neutrino_losses()
      call check_deuterium_energy()
      call check_below_cutoff()
      call check_strong_screening()
      call check_abundance_rates()
      call check_light_elements()
      call check_solar_neutrino_fluxes()
      call check_configuration()

      close(run_log_unit)
      write(*,'(a)') "test_net: done"

contains

      subroutine report(label, passed)
         character(len=*), intent(in) :: label
         logical, intent(in) :: passed
         write(*,'(a,a,a)') label, ': ', trim(merge('ok  ', 'FAIL', passed))
      end subroutine report

! burn_implicitly. With one reaction on, backward Euler has a closed
! form; tight tolerances make the comparison test the step's equation
! rather than its stopping rule. pp alone: dX/dt = -3 g X^2, so
! X1 = 2 X0/(1 + sqrt(1 + 12 g dt X0)) and He3 gains X0 - X1.
! C12(a,g) alone: with u = g dt Y1 C1, Y1 = Y0 - 4u, C1 = C0 - 12u,
! O1 = O0 + 16u, and 48 g dt u^2 - (1 + g dt (12 Y0 + 4 C0)) u
! + g dt Y0 C0 = 0. Then every rate of a real zone over a long step, at
! the production tolerances: the burned abundances, and the change in
! the sum of the mass fractions (no reaction emits mass there).
      subroutine check_burn_implicitly()
         double precision, parameter :: dt = 1.0d0, tight_abs = 1.0d-14, &
              tight_rel = 1.0d-13, min_abund = 1.0d-20
         double precision :: g(num_rate_rows), x0(n_network), x(n_network)
         double precision :: exact, a, b, c, u, gdt
         integer :: ierr

         write(*,'(a)') "# test_net: burn_implicitly"
         x0 = [0.70d0, 3.0d-5, 0.28d0, 3.5d-3, 4.0d-5, 1.0d-3, 8.0d-3, 2.0d-5]
         g = 0.0d0
         g(r_pp) = 10.0d0/(12.0d0*x0(s_h1))
         x = x0
         call burn_implicitly(x0, g, dt, min_abund, tight_abs, tight_rel, 50, x, ierr)
         exact = 2.0d0*x0(s_h1)/(1.0d0 + sqrt(1.0d0 + 12.0d0*g(r_pp)*dt*x0(s_h1)))
         call report("pp alone vs closed form (1e-12)", ierr == 0 .and. &
              abs(x(s_h1)/exact - 1.0d0) < 1.0d-12 .and. &
              abs(x(s_he3)/(x0(s_he3) + x0(s_h1) - exact) - 1.0d0) < 1.0d-12 .and. &
              all(x(s_he4:) == x0(s_he4:)))

         x0 = [0.0d0, 0.0d0, 0.5d0, 0.3d0, 0.0d0, 0.0d0, 0.2d0, 0.0d0]
         g = 0.0d0
         g(r_ac12) = 2.0d0
         gdt = g(r_ac12)*dt
         x = x0
         call burn_implicitly(x0, g, dt, min_abund, tight_abs, tight_rel, 50, x, ierr)
         a = 48.0d0*gdt
         b = 1.0d0 + gdt*(12.0d0*x0(s_he4) + 4.0d0*x0(s_c12))
         c = gdt*x0(s_he4)*x0(s_c12)
         u = 2.0d0*c/(b + sqrt(b*b - 4.0d0*a*c))
         call report("C12(a,g) alone vs closed form (1e-12)", ierr == 0 .and. &
              abs(x(s_he4)/(x0(s_he4) - 4.0d0*u) - 1.0d0) < 1.0d-12 .and. &
              abs(x(s_c12)/(x0(s_c12) - 12.0d0*u) - 1.0d0) < 1.0d-12 .and. &
              abs(x(s_o16)/(x0(s_o16) + 16.0d0*u) - 1.0d0) < 1.0d-12)

         x0 = [0.70d0, 3.0d-5, 0.28d0, 3.5d-3, 4.0d-5, 1.0d-3, 8.0d-3, 2.0d-5]
         call compute_network_rates(2.2d0, 7.3d0, test_xa, g)
         x = x0
         call burn_implicitly(x0, g, 1.0d17, 1.0d-20, 1.0d-5, 1.0d-4, 50, x, ierr)
         write(*,'(a,i0,a)') "all rates, log T 7.3, log rho 2.2, dt 1e17 s (ierr ", ierr, &
              "): H He3 He4 C12 C13 N14 O16 O18"
         write(*,'(4(1pe20.12))') x
         call report("sum X conserved (1e-14)", abs(sum(x) - sum(x0)) < 1.0d-14)
      end subroutine check_burn_implicitly

! The trace integrator: the surviving fraction of gas burning at a rate
! linear in time (exposure (r0 + r1) dt/2) or linear in log
! (exposure (r1 - r0) dt/ln(r1/r0); r0 dt when r1 = r0), and the mean
! surviving fraction of gas accreted uniformly over the step, against
! (1 - exp(-r dt))/(r dt) for a constant rate and against a fine
! Simpson rule (200000 intervals) otherwise.
      subroutine check_trace_integrator()
         integer, parameter :: ncase = 4
         double precision, parameter :: r0(ncase) = [1.0d-9, 5.0d-9, 2.0d-8, 3.0d-9], &
              r1(ncase) = [3.0d-9, 5.0d-9, 1.0d-10, 3.0d-11], &
              step(ncase) = [1.0d9, 1.0d9, 5.0d8, 2.0d9]
         double precision :: exposure, exact, err_survive, err_accreted, ref
         integer :: k, law

         write(*,'(a)') "# test_net: trace integrator"
         err_survive = 0.0d0
         err_accreted = 0.0d0
         do law = linear_in_time, linear_in_log
            do k = 1, ncase
               exposure = lambda(law, r0(k), r1(k), step(k), 1.0d0)
               exact = exp(-exposure)
               err_survive = max(err_survive, abs(surviving_fraction(r0(k), r1(k), &
                    step(k), law)/exact - 1.0d0))
               if (r0(k) == r1(k)) then
                  ref = (1.0d0 - exp(-r0(k)*step(k)))/(r0(k)*step(k))
               else
                  ref = accreted_reference(law, r0(k), r1(k), step(k))
               end if
               err_accreted = max(err_accreted, abs(accreted_survival(r0(k), r1(k), &
                    step(k), law)/ref - 1.0d0))
            end do
         end do
         call report("surviving_fraction vs closed form (1e-13)", err_survive < 1.0d-13)
         call report("accreted_survival vs reference (1e-6)", err_accreted < 1.0d-6)
         write(*,'(a,es8.1)') "accreted_survival max relative error:", err_accreted
      end subroutine check_trace_integrator

! the exposure from the start of the step to the fraction f of it
      double precision function lambda(law, r0, r1, dt, f)
         integer, intent(in) :: law
         double precision, intent(in) :: r0, r1, dt, f
         if (law == linear_in_log .and. r0 /= r1) then
            lambda = dt*r0*(exp(f*log(r1/r0)) - 1.0d0)/log(r1/r0)
         else
            lambda = dt*(r0*f + 0.5d0*(r1 - r0)*f*f)
         end if
      end function lambda

      double precision function accreted_reference(law, r0, r1, dt)
         integer, intent(in) :: law
         double precision, intent(in) :: r0, r1, dt
         integer, parameter :: n = 200000
         double precision :: total, w
         integer :: i
         total = lambda(law, r0, r1, dt, 1.0d0)
         accreted_reference = 0.0d0
         do i = 0, n
            w = merge(1.0d0, merge(4.0d0, 2.0d0, mod(i,2) == 1), i == 0 .or. i == n)
            accreted_reference = accreted_reference + &
                 w*exp(-(total - lambda(law, r0, r1, dt, dble(i)/dble(n))))
         end do
         accreted_reference = accreted_reference/(3.0d0*dble(n))
      end function accreted_reference

! compute_nuclear_energy at the production settings (cutoffs, weak
! screening threshold; neutrino fluxes and Itoh losses on) at three
! points: the solar centre, a CNO-burning core, a helium-burning core.
! Each prints the total, the pp, CNO and triple-alpha parts and the
! neutrino loss (erg/g/s), then |analytic/numerical - 1| of d eps/d ln T
! and d eps/d ln rho, against centred differences of 1e-4 dex. The
! derivatives hold the screening plasma and the branchings fixed, so the
! check allows 5e-3.
      subroutine check_nuclear_energy()
         integer, parameter :: npoint = 3
         double precision, parameter :: h = 1.0d-4
         double precision, parameter :: logt(npoint) = [7.196d0, 7.35d0, 8.10d0], &
              logd(npoint) = [2.176d0, 1.75d0, 4.0d0]
         type(nuclear_energy_result) :: e, ep, em
         double precision :: xa(15, npoint), fd
         double precision :: err_t, err_d, worst
         integer :: k

         call install_production_settings()
         xa = 0.0d0
         xa([i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, i_n14, i_o16, i_o18], 1) = &
              [0.34d0, 0.64d0, 0.02d0, 8.0d-6, 2.5d-5, 1.5d-6, 5.5d-3, 8.0d-3, 1.5d-5]
         xa([i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, i_n14, i_o16, i_o18], 2) = &
              [0.70d0, 0.28d0, 0.02d0, 1.0d-5, 3.0d-3, 4.0d-5, 1.0d-3, 8.0d-3, 2.0d-5]
         xa([i_h1, i_he4, i_metals, i_c12, i_n14, i_o16], 3) = &
              [0.0d0, 0.97d0, 0.03d0, 1.0d-2, 1.0d-2, 1.0d-2]

         write(*,'(a)') "# test_net: compute_nuclear_energy: eps, pp, CNO, " // &
              "triple alpha, neutrino loss; derivative errors (dlnT, dlnrho)"
         worst = 0.0d0
         do k = 1, npoint
            call compute_nuclear_energy(logd(k), logt(k), xa(:,k), .false., 0.0d0, e)
            write(*,'(a,i2,f7.3,f7.3)') "point", k, logt(k), logd(k)
            write(*,'(5(1pe20.12))') e%total_energy_gen_rate, e%eps_by_chain(chain_pp1), &
                 e%eps_by_chain(chain_cno), e%eps_by_chain(chain_triple_alpha), e%eps_neutrino_loss
            call compute_nuclear_energy(logd(k), logt(k)+h, xa(:,k), .false., 0.0d0, ep)
            call compute_nuclear_energy(logd(k), logt(k)-h, xa(:,k), .false., 0.0d0, em)
            fd = (ep%total_energy_gen_rate - em%total_energy_gen_rate)/(2.0d0*h*log(10.0d0))
            err_t = abs(e%deps_dlnt/fd - 1.0d0)
            call compute_nuclear_energy(logd(k)+h, logt(k), xa(:,k), .false., 0.0d0, ep)
            call compute_nuclear_energy(logd(k)-h, logt(k), xa(:,k), .false., 0.0d0, em)
            fd = (ep%total_energy_gen_rate - em%total_energy_gen_rate)/(2.0d0*h*log(10.0d0))
            err_d = abs(e%deps_dlnrho/fd - 1.0d0)
            write(*,'(2es10.2)') err_t, err_d
            worst = max(worst, err_t, err_d)
         end do
         call report("compute_nuclear_energy derivatives vs centred differences (5e-3)", &
              worst < 5.0d-3)
      end subroutine check_nuclear_energy

! The settings of a solar run: the default cutoffs and weak-screening
! threshold, with the solar neutrino fluxes and the Itoh et al. losses
! on (the defaults have both off), light-element rate scales 1, keeping
! the cross-section scales.
      subroutine install_production_settings()
         type(net_settings_t) :: s
         integer :: ierr
         s = net_settings
         s%logT_cutoffs = [6.5d0, 6.5d0, 6.82d0, 7.7d0, 7.5d0]
         s%weak_screening_threshold = 0.03d0
         s%calc_neutrinos = .true.
         s%use_itoh_neutrino_loss = .true.
         s%li6_rate_scale = 1.0d0
         s%li7_rate_scale = 1.0d0
         s%be9_pg_rate_scale = 1.0d0
         s%be9_pd_rate_scale = 1.0d0
         s%be9_palpha_rate_scale = 1.0d0
         call set_net_settings(s, ierr)
      end subroutine install_production_settings

! compute_abundance_rates, dX_i/dt / A_i of the network species (s_*
! order), on a hydrogen-burning zone (log T 7.3, the rates test's
! mixture) and a helium-burning one (log T 8.1, He4 0.97). The hydrogen
! shell timestep limiter reads the first entry.
      subroutine check_abundance_rates()
         double precision :: rates(num_rate_rows), xa(15), dxdt(n_network)

         call install_production_settings()
         write(*,'(a)') "# test_net: compute_abundance_rates, dX/dt / A: H He3 He4 C12 C13 N14 O16 O18"
         call compute_network_rates(2.2d0, 7.3d0, test_xa, rates)
         call compute_abundance_rates(rates, test_xa(network_species_slot), dxdt)
         write(*,'(4(1pe20.12))') dxdt
         xa = 0.0d0
         xa([i_he4, i_metals, i_c12, i_n14, i_o16]) = [0.97d0, 0.03d0, 1.0d-2, 1.0d-2, 1.0d-2]
         call compute_network_rates(4.0d0, 8.1d0, xa, rates)
         call compute_abundance_rates(rates, xa(network_species_slot), dxdt)
         write(*,'(4(1pe20.12))') dxdt
      end subroutine check_abundance_rates

! The light elements. burn_deuterium_region over 3e11 s on three shells
! with rates rising through the step: shell 1 alone, the three as one
! region, the three as the accreting surface region, and a region below
! the threshold; each prints the region's new X, He3 and D.
! compute_li_be_rates at three temperatures; deuterium_rate inside and
! outside the surface convection zone where the turnover cap binds (the
! capped rate is N_A/(m_d tau) per gram, times amu_cgs = 1/N_A: one
! deuteron mass per turnover time, 1/(m_d tau) per amu).
      subroutine check_light_elements()
         integer, parameter :: nz = 3
         double precision, parameter :: t9(3) = [3.0d-3, 5.0d-3, 1.0d-2]
         double precision :: xa_start(15, nz), mass(nz), rate0(nz), &
              rate1(nz), li6, li7, be9
         type(surface_accretion) :: none, accreting
         integer :: k

         call install_production_settings()
         write(*,'(a)') "# test_net: burn_deuterium_region: X, He3, D"
         mass = [1.0d30, 2.0d30, 3.0d30]
         rate0 = [1.0d-12, 2.0d-12, 4.0d-12]
         rate1 = [2.0d-12, 3.0d-12, 5.0d-12]
         xa_start = 0.0d0
         do k = 1, nz
            xa_start([i_h1, i_he4, i_he3, i_h2], k) = [0.70d0, 0.28d0, 3.0d-5, 1.0d-5*(k+1)]
         end do
         accreting%on = .true.
         accreting%inflow = .true.
         accreting%deuterium_fraction = 3.0d-5
         accreting%mass = 1.0d28
         call deuterium_case("shell 1", 1, 1, none, 1.0d-11, xa_start, mass, rate0, rate1)
         call deuterium_case("region 1-3", 1, 3, none, 1.0d-11, xa_start, mass, rate0, rate1)
         call deuterium_case("accreting surface region 1-3", 1, 3, accreting, 1.0d-11, xa_start, mass, rate0, rate1)
         call deuterium_case("below threshold", 1, 3, none, 1.0d-4, xa_start, mass, rate0, rate1)


         write(*,'(a)') "# test_net: compute_li_be_rates (T9 0.003, 0.005, 0.01): Li6, Li7, Be9"
         do k = 1, 3
            call compute_li_be_rates(t9(k), li6, li7, be9)
            write(*,'(3(1pe20.12))') li6, li7, be9
         end do

         write(*,'(a)') "# test_net: deuterium_rate, log T 6.5, log rho 1.5: outside, inside the surface CZ (turnover 1e6 s)"
         write(*,'(2(1pe20.12))') deuterium_rate(1.5d0, 6.5d0, 0.70d0, .false., 1.0d6), &
              deuterium_rate(1.5d0, 6.5d0, 0.70d0, .true., 1.0d6)
         call report("deuterium_rate capped at N_A/(m_d tau) in the surface CZ (1e-12)", &
              abs(deuterium_rate(1.5d0, 6.5d0, 0.70d0, .true., 1.0d6)/ &
              (1.0d0/(2.013553d0*1.0d6)) - 1.0d0) < 1.0d-12 .and. &
              deuterium_rate(1.5d0, 6.5d0, 0.70d0, .false., 1.0d6) > &
              deuterium_rate(1.5d0, 6.5d0, 0.70d0, .true., 1.0d6))
      end subroutine check_light_elements

      subroutine deuterium_case(label, first, last, accretion, threshold, &
           xa_start, mass, rate0, rate1)
         character(len=*), intent(in) :: label
         integer, intent(in) :: first, last
         type(surface_accretion), intent(in) :: accretion
         double precision, intent(in) :: threshold, xa_start(:,:), mass(:), &
              rate0(:), rate1(:)
         double precision :: comp(size(xa_start,1), size(xa_start,2))
         comp = xa_start
         call burn_deuterium_region(first, last, size(mass), mass, comp, 3.0d11, &
              xa_start, rate0, rate1, accretion, threshold)
         write(*,'(a)') label
         write(*,'(3(1pe20.12))') comp(i_h1, first), comp(i_he3, first), comp(i_h2, first)
      end subroutine deuterium_case

! The solar neutrino fluxes and the Be7 abundance at the solar centre
! (compute_nuclear_energy, neutrino fluxes on).
      subroutine check_solar_neutrino_fluxes()
         type(nuclear_energy_result) :: e
         double precision :: xa(15)

         call install_production_settings()
         xa = 0.0d0
         xa([i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, i_n14, i_o16, i_o18]) = &
              [0.34d0, 0.64d0, 0.02d0, 8.0d-6, 2.5d-5, 1.5d-6, 5.5d-3, 8.0d-3, 1.5d-5]
         call compute_nuclear_energy(2.176d0, 7.196d0, xa, .false., 0.0d0, e)
         write(*,'(a)') "# test_net: solar neutrino fluxes at the solar centre, then the Be7 mass fraction"
         write(*,'(5(1pe20.12))') e%neutrino_flux
         write(*,'(1pe20.12)') e%be7_mass_fraction
      end subroutine check_solar_neutrino_fluxes
! The thermal neutrino losses of both fits, on a carbon-oxygen mixture
! where nothing burns (no H, no He4), so eps is the loss alone: its
! value and |analytic/numerical - 1| of d eps/d ln T and d eps/d ln rho
! (centred differences of 1e-4 dex), at log T 8.1 (photo and plasma
! neutrinos) and 8.6 (pair too, T9 >= 0.2), log rho 4 and 5: Beaudet,
! Petrosian & Salpeter (the default) and Itoh et al. (to 1e-5: sneut's
! bremsstrahlung derivative is good to about 2e-6 of the total there).
! The loss is negative.
      subroutine check_thermal_neutrino_losses()
         integer, parameter :: npoint = 2
         double precision, parameter :: logt(npoint) = [8.1d0, 8.6d0], &
              logd(npoint) = [4.0d0, 5.0d0]
         type(net_settings_t) :: s
         double precision :: xa(15), err_t, err_d, worst(2)
         logical :: negative(2)
         integer :: k, fit, ierr

         xa = 0.0d0
         xa([i_metals, i_c12, i_o16]) = [1.0d0, 0.5d0, 0.5d0]
         write(*,'(a)') "# test_net: thermal neutrino losses on C/O (BPS, then Itoh): " // &
              "eps, derivative errors (dlnT, dlnrho)"
         do fit = 1, 2
            call install_production_settings()
            s = net_settings
            s%use_itoh_neutrino_loss = fit == 2
            call set_net_settings(s, ierr)
            worst(fit) = 0.0d0
            negative(fit) = .true.
            do k = 1, npoint
               call energy_derivative_errors(logd(k), logt(k), xa, .false., 0.0d0, &
                    err_t, err_d)
               worst(fit) = max(worst(fit), err_t, err_d)
               negative(fit) = negative(fit) .and. last_eps < 0.0d0
               write(*,'(i2,f6.2,f6.2,1pe20.12,2es10.2)') fit, logt(k), logd(k), &
                    last_eps, err_t, err_d
            end do
         end do
         call report("BPS neutrino loss negative, derivatives vs centred differences (1e-6)", &
              negative(1) .and. worst(1) < 1.0d-6)
         call report("Itoh neutrino loss negative, derivatives vs centred differences (1e-5)", &
              negative(2) .and. worst(2) < 1.0d-5)
      end subroutine check_thermal_neutrino_losses

! |analytic/numerical - 1| of compute_nuclear_energy's d eps/d ln T and
! d eps/d ln rho at one point, against centred differences of 1e-4 dex;
! last_eps keeps the point's eps.
      subroutine energy_derivative_errors(logd, logt, xa, in_surface_cz, &
           turnover, err_t, err_d)
         double precision, intent(in) :: logd, logt, xa(:), turnover
         logical, intent(in) :: in_surface_cz
         double precision, intent(out) :: err_t, err_d
         double precision, parameter :: h = 1.0d-4
         type(nuclear_energy_result) :: e, ep, em
         double precision :: fd

         call compute_nuclear_energy(logd, logt, xa, in_surface_cz, turnover, e)
         last_eps = e%total_energy_gen_rate
         call compute_nuclear_energy(logd, logt+h, xa, in_surface_cz, turnover, ep)
         call compute_nuclear_energy(logd, logt-h, xa, in_surface_cz, turnover, em)
         fd = (ep%total_energy_gen_rate - em%total_energy_gen_rate)/(2.0d0*h*log(10.0d0))
         err_t = abs(e%deps_dlnt/fd - 1.0d0)
         call compute_nuclear_energy(logd+h, logt, xa, in_surface_cz, turnover, ep)
         call compute_nuclear_energy(logd-h, logt, xa, in_surface_cz, turnover, em)
         fd = (ep%total_energy_gen_rate - em%total_energy_gen_rate)/(2.0d0*h*log(10.0d0))
         err_d = abs(e%deps_dlnrho/fd - 1.0d0)
      end subroutine energy_derivative_errors

! Deuterium burning in a pre-main-sequence surface convection zone
! (log T 6.6, log rho 1, D 2e-5): outside the zone, and inside it with a
! turnover time of 1e7 s, where the cap binds (the rate per unit X there
! exceeds N_A/(m_d tau) about 500 times). Capped, eps_D = (N_A/m_d/tau)
! X_D Q depends on neither rho nor T, so d eps/d ln rho and d eps/d ln T
! are pp's alone; both points are checked against centred differences
! (the capped point's to 1e-3: there the differences of an eps of 5e6
! resolve pp's derivatives of order 0.1 only to about 1e-4).
      subroutine check_deuterium_energy()
         double precision :: xa(15), eps_free, err_t, err_d, worst_free, worst_capped
         double precision, parameter :: turnover = 1.0d7
         double precision :: capped_eps

         call install_production_settings()
         xa = 0.0d0
         xa([i_h1, i_he4, i_metals, i_he3, i_h2]) = [0.70d0, 0.28d0, 0.02d0, 3.0d-5, 2.0d-5]
         write(*,'(a)') "# test_net: deuterium energy, log T 6.6, log rho 1: eps and " // &
              "derivative errors (dlnT, dlnrho), outside and inside the surface CZ (capped)"
         call energy_derivative_errors(1.0d0, 6.6d0, xa, .false., turnover, err_t, err_d)
         eps_free = last_eps
         worst_free = max(err_t, err_d)
         write(*,'(1pe20.12,2es10.2)') last_eps, err_t, err_d
         call energy_derivative_errors(1.0d0, 6.6d0, xa, .true., turnover, err_t, err_d)
         capped_eps = last_eps
         worst_capped = max(err_t, err_d)
         write(*,'(1pe20.12,2es10.2)') last_eps, err_t, err_d
         call report("deuterium energy derivatives vs centred differences (1e-6)", &
              worst_free < 1.0d-6)
         call report("capped deuterium: eps below the free rate's, derivatives vs " // &
              "centred differences (1e-3)", capped_eps < eps_free .and. worst_capped < 1.0d-3)
      end subroutine check_deuterium_energy

! Below the burning cutoff (log T <= logT_cutoffs(1), 6.5) nothing is
! computed: burning is off and every value is zero, at the cutoff too.
      subroutine check_below_cutoff()
         type(nuclear_energy_result) :: e1, e2
         double precision :: xa(15)

         call install_production_settings()
         xa = 0.0d0
         xa([i_h1, i_he4, i_metals, i_he3, i_h2]) = [0.70d0, 0.28d0, 0.02d0, 3.0d-5, 2.0d-5]
         call compute_nuclear_energy(1.0d0, 6.4d0, xa, .false., 0.0d0, e1)
         call compute_nuclear_energy(1.0d0, 6.5d0, xa, .false., 0.0d0, e2)
         call report("below the burning cutoff: nothing burns, everything zero", &
              .not. e1%burning .and. .not. e2%burning .and. &
              e1%total_energy_gen_rate == 0.0d0 .and. e2%total_energy_gen_rate == 0.0d0 .and. &
              all(e1%eps_by_chain == 0.0d0) .and. e1%deps_dlnt == 0.0d0 .and. &
              e1%deps_dlnrho == 0.0d0 .and. e1%eps_neutrino_loss == 0.0d0)
      end subroutine check_below_cutoff

! Strong screening (Graboske et al. 1973): at log T 7, log rho 6 the
! weak-screening exponent Lambda_0 z_tilde Z1 Z2 is at least 5 for every
! reaction, so each is strongly screened. Its rates' d ln rate/d ln T and
! d ln rate/d ln rho against centred differences (1e-4 dex), and the
! screening itself: ln(pp rate) less its value with weak screening
! (forced with a threshold of 1e3) is the strong exponent less the weak
! one, the strong one from Graboske et al.'s formula with pp's charges
! (Z1 = Z2 = 1; z_p = (Z1+Z2)^p - Z1^p - Z2^p, here at run time), and it
! is negative; the same for triple alpha at log T 7.8, log rho 8 (Z = 2,
! 6, which give z53 = 9.014).
      subroutine check_strong_screening()
         double precision, parameter :: h = 1.0d-4, logt = 7.0d0, logd = 6.0d0
         type(nuclear_rates) :: nr0, nrp, nrm, nr_weak
         type(zone_conditions) :: zone
         type(net_settings_t) :: s
         double precision :: fd, worst, weak_min, strong, ln_ratio, &
              strong_3a, weak_3a, ln_ratio_3a
         integer :: k, ierr

         call install_production_settings()
         call compute_zone_conditions(logd, logt, fd_mass, zone)
         weak_min = zone%plasma%weak_unit
         associate (p => zone%plasma)
         strong = 0.624d0*p%mean_charge_13*p%lambda0_23*(zp(5.0d0/3.0d0, 1, 1) + &
              0.316d0*p%mean_charge_13*zp(4.0d0/3.0d0, 1, 1) + &
              0.737d0*zp(2.0d0/3.0d0, 1, 1)/(p%mean_charge*p%lambda0_23))
         end associate
         call compute_reaction_rates(logd, logt, fd_mass, nr0)
         worst = 0.0d0
         do k = 1, num_reactions
            if (nr0%rate(k) <= 0.0d0) cycle
            call compute_reaction_rates(logd, logt+h, fd_mass, nrp)
            call compute_reaction_rates(logd, logt-h, fd_mass, nrm)
            fd = (log(nrp%rate(k)) - log(nrm%rate(k)))/(2.0d0*h*log(10.0d0))
            worst = max(worst, abs(nr0%dlnrate_dlnt(k)/fd - 1.0d0))
            call compute_reaction_rates(logd+h, logt, fd_mass, nrp)
            call compute_reaction_rates(logd-h, logt, fd_mass, nrm)
            fd = (log(nrp%rate(k)) - log(nrm%rate(k)))/(2.0d0*h*log(10.0d0))
            worst = max(worst, abs(nr0%dlnrate_dlnrho(k)/fd - 1.0d0))
         end do
         s = net_settings
         s%weak_screening_threshold = 1.0d3
         call set_net_settings(s, ierr)
         call compute_reaction_rates(logd, logt, fd_mass, nr_weak)
         ln_ratio = log(nr0%rate(r_pp)) - log(nr_weak%rate(r_pp))
         call compute_reaction_rates(8.0d0, 7.8d0, fd_mass, nrm)
         call install_production_settings()
         call compute_reaction_rates(8.0d0, 7.8d0, fd_mass, nrp)
         call compute_zone_conditions(8.0d0, 7.8d0, fd_mass, zone)
         weak_3a = 12.0d0*zone%plasma%weak_unit
         associate (p => zone%plasma)
         strong_3a = 0.624d0*p%mean_charge_13*p%lambda0_23*(zp(5.0d0/3.0d0, 2, 6) + &
              0.316d0*p%mean_charge_13*zp(4.0d0/3.0d0, 2, 6) + &
              0.737d0*zp(2.0d0/3.0d0, 2, 6)/(p%mean_charge*p%lambda0_23))
         end associate
         ln_ratio_3a = log(nrp%rate(r_3alpha)) - log(nrm%rate(r_3alpha))
         write(*,'(a)') "# test_net: strong screening, log T 7, log rho 6: " // &
              "Lambda_0 z_tilde, pp rate (strong, weak), worst derivative error"
         write(*,'(1pe20.12,1pe20.12,1pe20.12,es10.2)') weak_min, nr0%rate(r_pp), &
              nr_weak%rate(r_pp), worst
         call report("strong screening: in the strong regime, below weak, derivatives " // &
              "vs centred differences (1e-6)", weak_min >= 5.0d0 .and. &
              nr0%rate(r_pp) < nr_weak%rate(r_pp) .and. worst < 1.0d-6)
         call report("strong screening: the exponent is Graboske et al.'s, pp and " // &
              "triple alpha (1e-10)", abs(ln_ratio - (strong - weak_min)) < 1.0d-10 .and. &
              weak_3a >= 5.0d0 .and. abs(ln_ratio_3a - (strong_3a - weak_3a)) < 1.0d-10)
      end subroutine check_strong_screening

! Graboske et al.'s z_p = (Z1+Z2)^p - Z1^p - Z2^p, the test's own copy
      double precision function zp(p, z1, z2)
         double precision, intent(in) :: p
         integer, intent(in) :: z1, z2
         zp = dble(z1 + z2)**p - dble(z1)**p - dble(z2)**p
      end function zp

! burn_implicitly's error contract: a step that needs more Newton
! iterations than it is given returns ierr = 1 (here the every-rate
! step of check_burn_implicitly with a single iteration); with enough,
! ierr = 0.
      subroutine check_burn_nonconvergence()
         double precision :: g(num_rate_rows), x0(n_network), x(n_network)
         integer :: ierr1, ierr50

         x0 = [0.70d0, 3.0d-5, 0.28d0, 3.5d-3, 4.0d-5, 1.0d-3, 8.0d-3, 2.0d-5]
         call compute_network_rates(2.2d0, 7.3d0, test_xa, g)
         x = x0
         call burn_implicitly(x0, g, 1.0d17, 1.0d-20, 1.0d-5, 1.0d-4, 1, x, ierr1)
         x = x0
         call burn_implicitly(x0, g, 1.0d17, 1.0d-20, 1.0d-5, 1.0d-4, 50, x, ierr50)
         call report("burn_implicitly: ierr 1 when out of iterations, 0 when converged", &
              ierr1 == 1 .and. ierr50 == 0)
      end subroutine check_burn_nonconvergence

! The configuration's contracts: the reaction table passes its
! self-check (network_table_errors: every row consistent, every slot
! computed by exactly one row); the legacy sstandard vector lands on the
! reaction and extra scales by its old numbering (entry i = 0.1 i;
! entries 9 and 13, which belong to no reaction, give a warning); and a
! scale set on a reaction whose rate law takes none (C12(a,g)) is the
! one set_net_settings reports as having no effect.
      subroutine check_configuration()
         double precision :: sstandard(num_sstandard), xs(num_reactions), &
              extra(num_extra_scales)
         type(net_settings_t) :: s
         integer :: i

         write(*,'(a)') "# test_net: configuration"
         call report("reaction table self-check (no errors)", network_table_errors() == 0)
         sstandard = [(0.1d0*i, i = 1, num_sstandard)]
         call cross_section_scales_from_sstandard(sstandard, xs, extra)
         call report("sstandard onto the reaction and extra scales", &
              all(xs(r_pp:r_po16) == sstandard(1:7)) .and. xs(r_ac13) == sstandard(8) .and. &
              xs(r_ac12) == sstandard(10) .and. xs(r_an14) == sstandard(11) .and. &
              xs(r_3alpha) == sstandard(12) .and. all(extra == sstandard(14:17)))
         s = net_settings
         s%cross_section_scale(r_ac12) = 2.0d0
         call report("a scale on C12(a,g) is reported as without effect", &
              inert_scale_count(s) == 1 .and. inert_scale_count(net_settings) == 0)
      end subroutine check_configuration
end program test_net
