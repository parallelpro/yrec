!----------------------------------------------------------------------
! test_eos
!----------------------------------------------------------------------
! Standalone eos-domain test program (2026, phase three -- docs/ROADMAP.md
! stage 2), following MESA's per-module test convention: boot only
! what the eos domain needs, evaluate the public facade entries over
! a fixed grid of points, and print the results in a stable format
! for byte-comparison against the checked-in expected output
! (expected_test_eos.out, regenerated deliberately and reviewed like
! any other diff).
!
! Initialization notes -- everything here replicates what the full
! code does at startup, with the source of each piece cited:
!  * unit numbers: core/read_input.f90's assignments (run_log_unit=20,
!    fermi_unit=15, scv units 72-74, iopale=49).
!  * flags/cutoffs: explicit assignments below (TSCUT=6.0 matching
!    the reference solar-model namelists; all kap table flags false
!    so setups' kap_init is a no-op; atm_choice=0 so atm_init is a
!    no-op).
!  * constants + eos table loads: the real core/run/run_setup.f90 +
!    eos_lib's eos_init, with real Fermi/SCV table paths and dummy
!    paths for everything gated off.
!  * env_comp (the envelope-composition state blend_eos consumes):
!    core/model/starting_model.f90's exact mixture algorithm (its statements at the
!    "COMPUTE SURFACE MIX VALUES" block), seeded from const_lib's
!    mixture_weights_seed defaults (the G&N93 12-species mixture) and starin's
!    atomic_weight data, for X=0.70, Z=0.016492 (matching the shipped
!    EOSOPAL06Z0.016492 table).
!  * the OPAL-2006 EOS file is opened on iopale exactly as
!    core/read_input.f90 does; the table itself is lazily read by esac06
!    on first use, as in production.
!
! The MHD path is exercised only if YREC_MHD_TABLES is set to a
! directory holding the 8 MHD table files -- the repository ships no
! MHD tables, so by default that section reports SKIPPED. This is
! the infrastructure half of closing docs/ROADMAP.md's LMHD coverage gap;
! the data half needs tables from outside the repo.
!
! Then (2026-10, check_configurations) every EOS configuration a deck
! can choose, over the same grid: the OPAL 1995 and 2001 tables over
! SCV, SCV alone, the Yale solvers alone, numerical derivatives,
! use_diffusion_z, and the Debye-Hueckel correction.
!
! And the MHD table interpolation on a synthetic table
! (check_mhd_table_lookup), which runs without the MHD tables.
program test_eos
      use io_units
      use eos_lib
      use eos_configuration, only: eos_settings
      use opal_eos_tables, only: opal06
      use opal_eos, only: esac06
      use mhd_eos, only: interpolate_mhd_table
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      type(eos_settings_t) :: settings
      integer :: gerr

      character(len=256) :: yrec_input, mhd_dir
      character(len=256) :: fermi_path, scvh_path, scvhe_path, scvz_path, &
           opal06_path, dummy_path
! core/model/starting_model.f90's atomic weights and the default of the
! mixture_weights_seed control (the G&N93 12-species mixture), in the
! ix_* order of the envelope mixture
      double precision, parameter :: atomic_weight(12) = [23.0d0, 26.99d0, &
           24.32d0, 55.86d0, 28.1d0, 12.015d0, 1.008d0, 16.0d0, 14.01d0, &
           39.96d0, 20.19d0, 4.004d0]
      double precision, parameter :: mixture_weights_seed(12) = [0.001999d0, &
           0.003238d0, 0.037573d0, 0.071794d0, 0.040520d0, 0.173285d0, 0.0d0, &
           0.482273d0, 0.053152d0, 0.005379d0, 0.098668d0, 0.0d0]
! the envelope mixture handed to eos_set_mixture
      double precision :: envelope_hydrogen_fraction, envelope_metal_fraction, &
           mean_atomic_mass, species_fractions(12)

      double precision :: w(12), wsum, scale
      integer :: i, ipt, eos_ierr

! the EOS variables the sections print (see eval_eos)
      double precision :: logt, t, logp, p, logd, d
      double precision :: x_frac, z_frac
      double precision :: beta, betai, beta14, fxion(3), rmu, amu, emu, eta
      double precision :: qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp
      logical :: lderiv, latmo, in_table
      integer :: ksaha

! eos_get_gamma1 argument set
      double precision :: t6, rho, pres, gamma1, grad_ad

! test point grids
      integer, parameter :: npts = 8
      double precision :: grid_logt(npts), grid_logp(npts)
      data grid_logt /3.70d0, 3.90d0, 4.50d0, 5.50d0, 6.30d0, 6.90d0, &
           7.10d0, 7.161d0/
      data grid_logp /4.5d0, 6.0d0, 8.0d0, 11.5d0, 14.0d0, 16.0d0, &
           16.9d0, 17.35d0/

      integer, parameter :: ng1 = 4
      double precision :: g1_t6(ng1), g1_rho(ng1), g1_p(ng1)
! the synthetic MHD tables' variable count and X step (check_mhd_reader)
      integer, parameter :: mhd_ivarc_test = 20
      real(dp), parameter :: dx_centre_test = 0.01d0
      data g1_t6  /0.006d0, 0.10d0,  2.0d0,   12.0d0/
      data g1_rho /1.0d-7,  1.0d-4,  1.0d-1,  80.0d0/
      data g1_p   /3.16d4,  1.0d8,   3.16d13, 1.6d17/

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_eos: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      fermi_path  = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      scvh_path   = trim(yrec_input)//"/eos/scv/h_tab_i.dat"
      scvhe_path  = trim(yrec_input)//"/eos/scv/he_tab_i.dat"
      scvz_path   = trim(yrec_input)//"/eos/scv/z_tab_i.dat"
      opal06_path = trim(yrec_input)//"/eos/opal2006/EOSOPAL06Z0.016492"
      dummy_path = ""

      run_log_unit = 20
      open(run_log_unit, file="test_eos.short", status="replace")

! the envelope mixture, by core/model/starting_model.f90's algorithm
! (its "COMPUTE SURFACE MIX VALUES" block), for X = 0.70, Z = 0.016492
! (the Z of the shipped EOSOPAL06Z0.016492 table)
      envelope_hydrogen_fraction = 0.70d0
      envelope_metal_fraction = 0.016492d0
      w = mixture_weights_seed
      wsum = w(1)+w(2)+w(3)+w(4)+w(5)+w(6)+w(8)+w(9)+w(10)+w(11)
      scale = envelope_metal_fraction/wsum
      w(7) = envelope_hydrogen_fraction/scale
      w(12) = (1.0d0 - envelope_hydrogen_fraction - envelope_metal_fraction)/scale
      wsum = 0.0d0
      do i = 1, 12
         w(i) = scale*w(i)/atomic_weight(i)
         wsum = wsum + w(i)
      end do
      mean_atomic_mass = wsum
      scale = 1.0d0/mean_atomic_mass
      do i = 1, 12
         species_fractions(i) = w(i)*scale
      end do
      call eos_set_mixture(envelope_hydrogen_fraction, envelope_metal_fraction, &
           mean_atomic_mass, species_fractions)

! the production EOS (SCV and OPAL 2006, the Saha cutoff at log T = 6,
! as in the reference solar decks), installed with the tables it reads;
! the OPAL table is read by esac06 on first use, as in a run
      settings%use_scv_eos = .true.
      settings%use_opal2006_eos = .true.
      settings%saha_log10t_cutoff = 6.0d0
      call eos_init(settings, fermi_path, scvh_path, scvhe_path, scvz_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, opal06_path, gerr)
      if (gerr /= 0) then
         write(*,'(a)') "test_eos: FAIL (eos_init error)"
         stop 1
      end if
! the SCV envelope-mixture table, which a run builds when it reads the
! starting model (eos_init reads the SCV tables but does not build it)
      call eos_update_envelope_tables(envelope_hydrogen_fraction, envelope_metal_fraction)

      x_frac = 0.70d0
      z_frac = 0.016492d0

      write(*,'(a)') "# test_eos: eos_eval over (logT, logP) grid, " // &
           "X=0.70 Z=0.016492, SCV+OPAL2006, no MHD, no DH"
      do ipt = 1, npts
         logt = grid_logt(ipt)
         logp = grid_logp(ipt)
         logd = 0.0d0
         d = 0.0d0
         ksaha = 0
         lderiv = .true.
         latmo = .false.
         fxion = 0.0d0
         call eval_eos(x_frac, z_frac, lderiv, latmo, gerr)
      if (gerr /= 0) stop 1
         write(*,'(a,i2)') "point ", ipt
         write(*,'(4(1pe24.15))') logt, logp, logd, beta
         write(*,'(4(1pe24.15))') qdt, qdp, qcp, dela
         write(*,'(4(1pe24.15))') rmu, amu, emu, eta
         write(*,'(3(1pe24.15))') fxion(1), fxion(2), fxion(3)
         write(*,'(4(1pe24.15))') qdtt, qdtp, qat, qap
         write(*,'(2(1pe24.15),i4)') qcpt, qcpp, ksaha
      end do

      write(*,'(a)') "# test_eos: eos_get_gamma1, OPAL2006 branch"
      do ipt = 1, ng1
         t6 = g1_t6(ipt)
         rho = g1_rho(ipt)
         pres = g1_p(ipt)
         ksaha = 0
         call eos_get_gamma1(x_frac, z_frac, t6, rho, pres, gamma1, &
              grad_ad, ksaha, ierr=gerr)
      if (gerr /= 0) stop 1
         write(*,'(a,i2,2(1pe24.15))') "g1 ", ipt, gamma1, grad_ad
      end do

! Only the two hotter points here: at the two low-T points the
! Yale/SCV branch's Cox & Guili gamma1 reconstruction consumes
! quantities blend_eos leaves undefined under the branch's historical
! lderiv=.false. setting (inherited verbatim from calcad's original
! formulation; NaN/uninitialized in practice) -- a pre-existing
! property of that code path, unreachable in production configs
! (which set LOPALE06), and not meaningfully pinnable.
      write(*,'(a)') "# test_eos: eos_get_gamma1, Yale/SCV branch " // &
           "(defined regime only)"
      settings = eos_settings
      settings%use_opal2006_eos = .false.
      call set_eos_settings(settings)
      do ipt = 3, ng1
         t6 = g1_t6(ipt)
         rho = g1_rho(ipt)
         pres = g1_p(ipt)
         ksaha = 0
         call eos_get_gamma1(x_frac, z_frac, t6, rho, pres, gamma1, &
              grad_ad, ksaha, ierr=gerr)
      if (gerr /= 0) stop 1
         write(*,'(a,i2,2(1pe24.15))') "g1 ", ipt, gamma1, grad_ad
      end do
      settings%use_opal2006_eos = .true.
      call set_eos_settings(settings)

! MHD path: infrastructure for closing the LMHD coverage gap; runs
! only when tables are supplied from outside the repository.
      call get_environment_variable("YREC_MHD_TABLES", mhd_dir)
      if (len_trim(mhd_dir) == 0) then
         write(*,'(a)') "# test_eos: MHD path SKIPPED " // &
              "(YREC_MHD_TABLES not set; no MHD tables ship with YREC)"
      else
         settings = eos_settings
         settings%use_mhd_eos = .true.
         call eos_init(settings, fermi_path, scvh_path, scvhe_path, scvz_path, &
              trim(mhd_dir)//"/zams_a", trim(mhd_dir)//"/zams_b", &
              trim(mhd_dir)//"/zams_c", trim(mhd_dir)//"/centre1", &
              trim(mhd_dir)//"/centre2", trim(mhd_dir)//"/centre3", &
              trim(mhd_dir)//"/centre4", trim(mhd_dir)//"/centre5", opal06_path, gerr)
      if (gerr /= 0) stop 1
         write(*,'(a)') "# test_eos: eos_eval over the same grid, MHD"
         do ipt = 1, npts
            logt = grid_logt(ipt)
            logp = grid_logp(ipt)
            logd = 0.0d0
            ksaha = 0
            lderiv = .true.
            latmo = .false.
            fxion = 0.0d0
            call eval_eos(x_frac, z_frac, lderiv, latmo, gerr)
      if (gerr /= 0) stop 1
            write(*,'(a,i2)') "mhd point ", ipt
            write(*,'(4(1pe24.15))') logt, logp, logd, beta
            write(*,'(4(1pe24.15))') qdt, qdp, qcp, dela
         end do
      end if

! Error paths (2026, docs/ROADMAP.md stage 3): in the eos domain the error
! conditions behind the converted stops are corrupt-table /
! misconfiguration cases (out-of-table lookups stay *recoverable* by
! design here -- in_table flags and sentinels -- unlike kap), so
! with healthy tables the facade cannot be driven to ierr /= 0. The
! first check asserts the success path threads ierr = 0 end to end
! through compute_blended_eos/blend_eos/compute_opal2006_eos/esac06. The second drives esac06 (a
! domain internal -- this test is in-domain, white-box by design)
! with an invalid rad_flag, the same class of internal error the
! converted stops guarded, and asserts ierr = 1 with no crash; its
! diagnostic goes to run_log_unit as always.
      write(*,'(a)') "# test_eos: error paths via ierr"
      logt = 6.30d0
      logp = 14.0d0
      logd = 0.0d0
      d = 0.0d0
      ksaha = 0
      lderiv = .true.
      latmo = .false.
      fxion = 0.0d0
      call eval_eos(x_frac, z_frac, lderiv, latmo, eos_ierr)
      write(*,'(a,i4)') "err facade-success ierr = ", eos_ierr
      call esac06(opal06, 0.70d0, 2.0d0, 0.1d0, 9, 2, eos_ierr, in_table)
      write(*,'(a,i4)') "err bad-rad-flag   ierr = ", eos_ierr

      call check_configurations()
      call check_derivatives()
      call check_mhd_table_lookup()
      call check_mhd_reader()

      close(run_log_unit)
      write(*,'(a)') "test_eos: done"
contains

! Every EOS configuration a deck can choose, each evaluated with eos_get
! over the same (log T, log P) grid as the first section:
!   opal95, opal2001   the OPAL 1995 and 2001 tables, over SCV
!   scv_only           SCV at low T, no OPAL table
!   yale_only          the Yale Saha and fully ionized solvers alone
!   numerical          OPAL 2006 + SCV with numerical derivatives
!   diffusion_z        OPAL 2006 + SCV with use_diffusion_z
!   debye_huckel       the Yale solvers with the Debye-Hueckel correction,
!                      the composition passed per zone
! Each point prints log rho, beta, dln rho/dln T, dln rho/dln P, cp,
! grad_ad, eta and the Saha state.
      subroutine check_configurations()
         character(len=256) :: legacy
         legacy = trim(yrec_input)//"/eos/legacy/"
         call configure_eos("opal95", opal95=.true., scv=.true., file=trim(legacy)//"opal/OPALeos.0188")
         call evaluate_grid(0.70d0, 0.0188d0)
         call configure_eos("opal2001", opal2001=.true., scv=.true., file=trim(legacy)//"opal2001/O01EOS0188z")
         call evaluate_grid(0.70d0, 0.0188d0)
         call configure_eos("scv_only", scv=.true.)
         call evaluate_grid(0.70d0, 0.016492d0)
         call configure_eos("yale_only")
         call evaluate_grid(0.70d0, 0.016492d0)
         call configure_eos("numerical", opal2006=.true., scv=.true., numerical=.true., file=opal06_path)
         call evaluate_grid(0.70d0, 0.016492d0)
         call configure_eos("diffusion_z", opal2006=.true., scv=.true., diffusion_z=.true., file=opal06_path)
         call evaluate_grid(0.70d0, 0.016492d0)
         call configure_eos("debye_huckel", debye_huckel=.true.)
         call evaluate_grid(0.70d0, 0.016492d0, with_composition=.true.)
      end subroutine check_configurations

! Set one configuration (the controls a deck would set), opening an OPAL
! file on the unit read_controls would use when one is named.
      subroutine configure_eos(label, opal95, opal2001, opal2006, scv, numerical, &
           diffusion_z, debye_huckel, file)
         character(len=*), intent(in) :: label
         logical, intent(in), optional :: opal95, opal2001, opal2006, scv, &
              numerical, diffusion_z, debye_huckel
         character(len=*), intent(in), optional :: file
         type(eos_settings_t) :: s
         character(len=256) :: opal_path
         integer :: ierr
         write(*,'(a,a)') "# test_eos: configuration ", label
         s = eos_settings
         s%use_opal95_eos = present(opal95)
         s%use_opal2001_eos = present(opal2001)
         s%use_opal2006_eos = present(opal2006)
         s%use_scv_eos = present(scv)
         s%use_numerical_derivatives = present(numerical)
         s%use_diffusion_z = present(diffusion_z)
         s%use_debye_huckel_correction = present(debye_huckel)
         s%debye_huckel_eta_min = -1.0d0
         s%debye_huckel_eta_max = 1.0d0
         opal_path = ""
         if (present(file)) opal_path = file
         call eos_init(s, fermi_path, scvh_path, scvhe_path, scvz_path, &
              dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
              dummy_path, dummy_path, dummy_path, opal_path, ierr)
         if (ierr /= 0) write(*,'(a,i0)') "eos_init ierr = ", ierr
      end subroutine configure_eos

      subroutine evaluate_grid(x, z, with_composition)
         double precision, intent(in) :: x, z
         logical, intent(in), optional :: with_composition
         double precision :: comp(15)
         integer :: k, ierr
         comp = 0.0d0
         comp(1) = x
         comp(2) = 1.0d0 - x - z
         comp(3) = z
         comp(5) = 0.17d0*z
         comp(7) = 0.05d0*z
         comp(9) = 0.48d0*z
         do k = 1, npts
            logt = grid_logt(k)
            logp = grid_logp(k)
            logd = 0.0d0
            d = 0.0d0
            ksaha = 0
            fxion = 0.0d0
            if (present(with_composition)) then
               call eval_eos(x, z, .true., .false., ierr, comp)
            else
               call eval_eos(x, z, .true., .false., ierr)
            end if
            if (ierr /= 0) then
               write(*,'(a,2f7.3,a,i0)') "eos logT logP", logt, logp, " ierr ", ierr
            else
               write(*,'(a,2f7.3,4(1pe23.15))') "eos", logt, logp, logd, beta, qdt, qdp
               write(*,'(a,14x,3(1pe23.15),i3)') "   ", qcp, dela, eta, ksaha
            end if
         end do
      end subroutine evaluate_grid

! eos_get on the test's variables: they go into an eos_result and come
! back out (log T and log P are inputs only)
      subroutine eval_eos(x, z, want, atm, ierr, comp)
         double precision, intent(in) :: x, z
         logical, intent(in) :: want, atm
         integer, intent(out) :: ierr
         double precision, intent(in), optional :: comp(:)
         type(eos_result) :: res
         res = eos_result(temperature=t, pressure=p, log10_density=logd, density=d, &
              beta=beta, beta_inverse=betai, beta14=beta14, ion_fraction=fxion, &
              specific_gas_constant=rmu, ion_mean_weight_inverse=amu, &
              electron_mean_weight_inverse=emu, electron_degeneracy_parameter=eta, &
              dlnrho_dlnt=qdt, dlnrho_dlnp=qdp, specific_heat_cp=qcp, &
              adiabatic_gradient=dela, dlnrho_dlnt_dt=qdtt, dlnrho_dlnt_dp=qdtp, &
              adiabatic_gradient_dt=qat, adiabatic_gradient_dp=qap, &
              specific_heat_cp_dt=qcpt, specific_heat_cp_dp=qcpp)
         call eos_get(logt, logp, x, z, res, want, atm, ksaha, &
              composition_at_zone=comp, ierr=ierr)
         t = res%temperature
         p = res%pressure
         logd = res%log10_density
         d = res%density
         beta = res%beta
         betai = res%beta_inverse
         beta14 = res%beta14
         fxion = res%ion_fraction
         rmu = res%specific_gas_constant
         amu = res%ion_mean_weight_inverse
         emu = res%electron_mean_weight_inverse
         eta = res%electron_degeneracy_parameter
         qdt = res%dlnrho_dlnt
         qdp = res%dlnrho_dlnp
         qcp = res%specific_heat_cp
         dela = res%adiabatic_gradient
         qdtt = res%dlnrho_dlnt_dt
         qdtp = res%dlnrho_dlnt_dp
         qat = res%adiabatic_gradient_dt
         qap = res%adiabatic_gradient_dp
         qcpt = res%specific_heat_cp_dt
         qcpp = res%specific_heat_cp_dp
      end subroutine eval_eos

! The derivative check: at points in each regime, every configuration's
! analytic derivatives against centred finite differences of its own
! values (a step of 1e-4 in log T and in log P). Each row prints the
! analytic value, the finite difference and their relative difference.
! The second-derivative outputs are checked as
!   dlnrho_dlnt_dt   d(dln rho/dln T)/dln T
!   dlnrho_dlnt_dp   d(dln rho/dln T)/dln P
!   adiabatic_gradient_dt, _dp   dln grad_ad/dln T, dln grad_ad/dln P
!   specific_heat_cp_dt, _dp     dln cp/dln T, dln cp/dln P
! Table interpolation adds noise to finite differences, so agreement is
! expected only where the values are analytic; the numbers are pinned
! either way. log T = 5.8 lies in the ramp between the Saha and the
! fully ionized solutions (saha_log10t_cutoff - 0.5 to the cutoff).
      subroutine check_derivatives()
         character(len=256) :: legacy
         legacy = trim(yrec_input)//"/eos/legacy/"
         write(*,'(a)') "# test_eos: derivative check (analytic, finite difference, relative difference)"
         call configure_eos("yale_only")
         call derivative_points("yale_only")
         call configure_eos("scv_only", scv=.true.)
         call derivative_points("scv_only")
         call configure_eos("production", opal2006=.true., scv=.true., file=opal06_path)
         call derivative_points("production")
         call configure_eos("opal95", opal95=.true., scv=.true., file=trim(legacy)//"opal/OPALeos.0188")
         call derivative_points("opal95")
         call configure_eos("opal2001", opal2001=.true., scv=.true., file=trim(legacy)//"opal2001/O01EOS0188z")
         call derivative_points("opal2001")
         call configure_eos("numerical", opal2006=.true., scv=.true., numerical=.true., file=opal06_path)
         call derivative_points("numerical")
         call configure_eos("debye_huckel", debye_huckel=.true.)
         call derivative_points("debye_huckel")
      end subroutine check_derivatives

         subroutine derivative_points(label)
            character(len=*), intent(in) :: label
            integer, parameter :: nd = 6
            double precision, parameter :: dlogt(nd) = [3.90d0, 4.50d0, 5.50d0, 5.80d0, 6.90d0, 7.161d0]
            double precision, parameter :: dlogp(nd) = [6.0d0, 8.0d0, 11.5d0, 12.6d0, 16.0d0, 17.35d0]
            double precision, parameter :: h = 1.0d-4, ln10 = 2.302585092994045684d0
            double precision :: c(11), tp(11), tm(11), pp(11), pm(11)
            integer :: k
            do k = 1, nd
               call eval_point(dlogt(k), dlogp(k), .true., c)
               call eval_point(dlogt(k) + h, dlogp(k), .false., tp)
               call eval_point(dlogt(k) - h, dlogp(k), .false., tm)
               call eval_point(dlogt(k), dlogp(k) + h, .false., pp)
               call eval_point(dlogt(k), dlogp(k) - h, .false., pm)
               write(*,'(a,1x,a,2f7.3)') "deriv", trim(label), dlogt(k), dlogp(k)
! slots of c: 1 log rho, 2 dlnrho_dlnt, 3 dlnrho_dlnp, 4 cp, 5 grad_ad,
! 6-11 the six second derivatives
               call row("dlnrho_dlnt", c(2), (tp(1) - tm(1))/(2.0d0*h))
               call row("dlnrho_dlnp", c(3), (pp(1) - pm(1))/(2.0d0*h))
               call row("dlnrho_dlnt_dt", c(6), (tp(2) - tm(2))/(2.0d0*h*ln10))
               call row("dlnrho_dlnt_dp", c(7), (pp(2) - pm(2))/(2.0d0*h*ln10))
               call row("grad_ad_dt", c(8), (log10(tp(5)) - log10(tm(5)))/(2.0d0*h))
               call row("grad_ad_dp", c(9), (log10(pp(5)) - log10(pm(5)))/(2.0d0*h))
               call row("cp_dt", c(10), (log10(tp(4)) - log10(tm(4)))/(2.0d0*h))
               call row("cp_dp", c(11), (log10(pp(4)) - log10(pm(4)))/(2.0d0*h))
            end do
         end subroutine derivative_points

         subroutine eval_point(lt, lp, want, v)
            double precision, intent(in) :: lt, lp
            logical, intent(in) :: want
            double precision, intent(out) :: v(11)
            integer :: ierr
            logt = lt
            logp = lp
            logd = 0.0d0
            d = 0.0d0
            ksaha = 0
            fxion = 0.0d0
            call eval_eos(0.70d0, 0.016492d0, want, .false., ierr)
            if (ierr /= 0) write(*,'(a,2f8.4,a,i0)') "eos logT logP", lt, lp, " ierr ", ierr
            v = [logd, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp]
         end subroutine eval_point

         subroutine row(name, analytic, fd)
            character(len=*), intent(in) :: name
            double precision, intent(in) :: analytic, fd
            double precision :: rel
            rel = (analytic - fd)/max(abs(fd), 1.0d-30)
            write(*,'(2x,a16,3(1pe13.4))') name, analytic, fd, rel
         end subroutine row
! The MHD tables' (log P, log T) interpolation on a 7 (T) x 6 (rho) x 3
! table whose column 2 is log P: inside, at the table corners (the
! stencils clamp), and above the highest pressure (it returns with the
! output unchanged).
subroutine check_mhd_table_lookup()
      integer, parameter :: nt = 7, nr = 6, nv = 3
      double precision :: table(nt,nr,nv), tlog(nt), w1(nv,4), w2(nv,4), yw(nv), &
           out(nv), log_p(4), log_t(4)
      logical :: in_range(4)
      integer :: it, ir, i
      write(*,'(a)') "# MHD table lookup (synthetic table)"
      do it = 1, nt
         tlog(it) = 3.5d0 + 0.25d0*it
         do ir = 1, nr
            table(it,ir,1) = -8.0d0 + ir
            table(it,ir,2) = 2.0d0*tlog(it) + 0.8d0*ir
            table(it,ir,3) = tlog(it)*ir - 0.1d0*ir*ir
         end do
      end do
      log_p = [11.3d0, 9.0d0, 13.7d0, 30.0d0]
      log_t = [4.6d0, 3.6d0, 5.3d0, 4.6d0]
      do i = 1, 4
         out = -9.0d0
         w1 = 0.0d0
         w2 = 0.0d0
         yw = 0.0d0
         call interpolate_mhd_table(log_p(i), log_t(i), table, nt, nr, nv, tlog, nt, nr, w1, w2, yw, &
              out, in_range(i))
         write(*,'(a,2f6.2,a,3(1pe24.15))') "lookup at", log_p(i), log_t(i), " =", out
      end do
      call check("mhd lookup: in range inside and at the corners, not above the table", &
           all(in_range(1:3)) .and. .not. in_range(4))
end subroutine check_mhd_table_lookup

! The MHD reader on a synthetic table set (the real tables are not
! shipped): eight unformatted files in the layout mhdst reads, written
! into the working directory -- three ZAMS-type compositions (X = 0.6,
! 0.7, 0.8; lower rows log T = 3.5-4.75, upper 4.5-5.75) and five
! centre-type ones (X = 0.1-0.9, log T = 5.5-6.75), each with its X -+ dX
! companions. Every row tabulates an ideal gas, log P = log rho + log T
! + c(X) with c = 7.9 + 0.3 X, so the interpolation in log P, log T and
! X reproduces log rho = log P - log T - c(X) exactly. Checks: eos_init
! reads the set, eos_get returns that density in each family (and at
! the ZAMS/centre boundary), and points outside the tables, a missing
! file and a truncated file return ierr.
subroutine check_mhd_reader()
      real(dp), parameter :: x_zams(3) = [0.6d0, 0.7d0, 0.8d0]
      real(dp), parameter :: x_centre(5) = [0.1d0, 0.3d0, 0.5d0, 0.7d0, 0.9d0]
      real(dp), parameter :: z_mhd = 0.02d0
      character(len=16), parameter :: zams_names(3) = [character(len=16) :: &
           "mhd_zams_a", "mhd_zams_b", "mhd_zams_c"]
      character(len=16), parameter :: centre_names(5) = [character(len=16) :: &
           "mhd_centre1", "mhd_centre2", "mhd_centre3", "mhd_centre4", "mhd_centre5"]
      real(dp), parameter :: test_log_t(4) = [4.1d0, 5.1d0, 6.1d0, 5.5d0]
      real(dp), parameter :: x_test = 0.65d0, log_rho_test = -6.0d0
      type(eos_settings_t) :: s
      type(eos_result) :: res
      integer :: i, u, ierr, ks
      real(dp) :: lp
      logical :: ok_read, ok_rho, ok_grad, ok_out_t, ok_out_p

      write(*,'(a)') "# MHD reader (synthetic tables)"
      do i = 1, 3
         call write_mhd_table(trim(zams_names(i)), 0, x_zams(i))
      end do
      do i = 1, 5
         call write_mhd_table(trim(centre_names(i)), 1, x_centre(i))
      end do
      s = eos_settings
      s%use_mhd_eos = .true.
      call eos_init(s, fermi_path, scvh_path, scvhe_path, scvz_path, &
           zams_names(1), zams_names(2), zams_names(3), centre_names(1), &
           centre_names(2), centre_names(3), centre_names(4), centre_names(5), &
           opal06_path, ierr)
      ok_read = ierr == 0
      call check("mhd reader: eos_init reads the eight tables", ok_read)
      ok_rho = ok_read
      ok_grad = ok_read
      do i = 1, size(test_log_t)
         if (.not. ok_read) exit
         lp = log_rho_test + test_log_t(i) + mhd_c(x_test)
         ks = 0
         res = eos_result()
         call eos_get(test_log_t(i), lp, x_test, z_mhd, res, .true., .false., ks, ierr=ierr)
         ok_rho = ok_rho .and. ierr == 0 .and. abs(res%log10_density - log_rho_test) < 1.0d-10
         ok_grad = ok_grad .and. ierr == 0 .and. abs(res%adiabatic_gradient - 0.4d0) < 1.0d-12 &
              .and. abs(res%dlnrho_dlnt + 1.0d0) < 1.0d-12 .and. abs(res%dlnrho_dlnp - 1.0d0) < 1.0d-12
      end do
      call check("mhd reader: log rho = log P - log T - c(X) in each family and at the boundary", ok_rho)
      call check("mhd reader: grad_ad, dln rho/dln T and dln rho/dln P as tabulated", ok_grad)
! Each point outside the tables follows one inside at the same log P,
! so that a stale result would pass the mean-weight check as success.
      ks = 0
      res = eos_result()
      call eos_get(6.9d0, lp, x_test, z_mhd, res, .true., .false., ks, ierr=ierr)
      ok_out_t = ierr /= 0
      call check("mhd reader: log T above the tables returns ierr", ok_out_t)
      lp = 5.6d0 + mhd_c(x_test)
      ks = 0
      res = eos_result()
      call eos_get(6.1d0, lp, x_test, z_mhd, res, .true., .false., ks, ierr=ierr)
      ok_out_p = ierr == 0 .and. abs(res%log10_density + 0.5d0) < 1.0d-10
      ks = 0
      res = eos_result()
      call eos_get(5.1d0, lp, x_test, z_mhd, res, .true., .false., ks, ierr=ierr)
      ok_out_p = ok_out_p .and. ierr /= 0
      call check("mhd reader: log P above the tables returns ierr", ok_out_p)
      call eos_init(s, fermi_path, scvh_path, scvhe_path, scvz_path, &
           zams_names(1), "mhd_missing", zams_names(3), centre_names(1), &
           centre_names(2), centre_names(3), centre_names(4), centre_names(5), &
           opal06_path, ierr)
      call check("mhd reader: a missing table returns ierr = 1", ierr == 1)
      open(newunit=u, file="mhd_truncated", form="unformatted", status="replace")
      write(u) mhd_ivarc_test, 0, 0, 0.0d0
      close(u)
      call eos_init(s, fermi_path, scvh_path, scvhe_path, scvz_path, &
           zams_names(1), "mhd_truncated", zams_names(3), centre_names(1), &
           centre_names(2), centre_names(3), centre_names(4), centre_names(5), &
           opal06_path, ierr)
      call check("mhd reader: a truncated table returns ierr = 1", ierr == 1)
      s%use_mhd_eos = .false.
      call set_eos_settings(s)
end subroutine check_mhd_reader

! the synthetic tables' log P - log rho - log T at hydrogen fraction x
pure real(dp) function mhd_c(x)
      real(dp), intent(in) :: x
      mhd_c = 7.9d0 + 0.3d0*x
end function mhd_c

! One synthetic MHD table file: kind 0 (ZAMS: one composition, a lower
! and an upper part) or kind 1 (centre: the composition and its X -+ dX
! companions, each one part).
subroutine write_mhd_table(path, kind, x)
      character(len=*), intent(in) :: path
      integer, intent(in) :: kind
      real(dp), intent(in) :: x
      integer, parameter :: n_rows = 6
      integer :: u, pass, n_pass
      real(dp) :: xp, yp, dx
      open(newunit=u, file=path, form="unformatted", status="replace")
      n_pass = merge(1, 3, kind == 0)
      dx = merge(0.0d0, dx_centre_test, kind == 0)
      do pass = 1, n_pass
         xp = x
         if (pass == 2) xp = x - dx
         if (pass == 3) xp = x + dx
         yp = 1.0d0 - xp - 0.02d0
         write(u) mhd_ivarc_test, kind, merge(0, merge(-1, 1, pass == 2), pass == 1), dx
         write(u) 3, 1.008d0, xp, xp, 4.003d0, yp, yp, 16.0d0, 0.02d0, 0.02d0, 0.6d0
         if (kind == 0) then
            write(u) n_rows, n_rows, 0.0d0, 0.0d0
            call write_mhd_part(u, 3.5d0, xp)
            call write_mhd_part(u, 4.5d0, xp)
         else
            write(u) 0, n_rows, 0.0d0, 0.0d0
            call write_mhd_part(u, 5.5d0, xp)
         end if
      end do
      close(u)
end subroutine write_mhd_table

! n_rows T rows from log T = log_t0 in steps of 0.25, each with six
! densities (log rho = -8.5 + 1.5 j) of an ideal gas at hydrogen
! fraction x: the 20 tabulated variables, of which compute_mhd_eos reads
! log rho, log P, chi_rho, chi_T, grad_ad, log cp, the ionization
! fractions, eta and log P_gas
subroutine write_mhd_part(u, log_t0, x)
      integer, intent(in) :: u
      real(dp), intent(in) :: log_t0, x
      integer, parameter :: n_rows = 6, n_rho = 6
      real(dp) :: v(mhd_ivarc_test), lt, lr
      integer :: it, j
      do it = 1, n_rows
         lt = log_t0 + 0.25d0*(it - 1)
         write(u) n_rho, lt
         do j = 1, n_rho
            lr = -8.5d0 + 1.5d0*j
            v = 0.0d0
            v(1) = lr
            v(2) = lr + lt + mhd_c(x)
            v(3) = 12.0d0 + 0.1d0*lt
            v(4) = 1.0d0
            v(5) = 1.0d0
            v(8) = 0.4d0
            v(9) = 8.5d0
            v(14) = 1.0d0
            v(15) = 0.5d0
            v(16) = 0.25d0
            v(18) = -3.0d0
            v(20) = v(2)
            write(u) v
         end do
      end do
end subroutine write_mhd_part

subroutine check(label, ok)
      character(len=*), intent(in) :: label
      logical, intent(in) :: ok
      write(*,'(a,1x,a)') merge("ok  ", "FAIL", ok), label
end subroutine check

end program test_eos
