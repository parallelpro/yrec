!----------------------------------------------------------------------
! test_eos
!----------------------------------------------------------------------
! The equation of state, through eos_lib, without a star: the program
! installs the settings and the tables itself (eos_init), sets the
! envelope mixture a run would (eos_set_mixture: core's mixture
! algorithm on the default mixture seed, X = 0.70, Z = 0.016492, the Z
! of the shipped OPAL 2006 table) and builds the SCV envelope table
! (eos_update_envelope_tables). Its output is byte-compared with
! expected_test_eos.out (regenerated only on purpose). In order:
!   the production EOS (SCV and OPAL 2006) over a (log T, log P) grid;
!   eos_get_gamma1 on its two branches; the MHD EOS when YREC_MHD_TABLES
!   names a directory with the real tables (none ship with YREC: SKIPPED
!   otherwise); the error paths; every configuration a deck can choose
!   (check_configurations); the derivatives against finite differences
!   (check_derivatives); the MHD table lookup on a synthetic table; then
!   the ok/FAIL checks:
!   check_mhd_reader  the MHD reader on eight synthetic tables written
!                     to the working directory (densities in each family,
!                     the out-of-range errors, a missing or truncated file)
!   check_physics     the ideal-gas and radiation limits, continuity
!                     across the Saha ramp and the OPAL edge ramp,
!                     in_atmosphere, the envelope table after a change
!                     of composition, the Debye-Hueckel correction
! Each ok/FAIL check was confirmed to fail under a mutation of the code
! it covers; a failed check makes the program stop with status 1. The
! run log (test_eos.short) gets eos's messages.
program test_eos
      use io_units
      use eos_lib
      use eos_configuration, only: eos_settings
      use opal_eos_tables, only: opal06
      use opal_eos, only: esac06, compute_opal_eos, iv_opal06
      use physical_constants, only: gas_constant
      use mhd_eos, only: interpolate_mhd_table
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      type(eos_settings_t) :: settings
! the tables eos_init reads: the shipped Fermi, SCV and OPAL 2006 tables
      type(eos_table_paths_t) :: table_paths
      integer :: gerr
! the ok/FAIL checks that failed (the program then stops with status 1)
      integer :: failures = 0

      character(len=256) :: yrec_input, mhd_dir
      type(eos_table_paths_t) :: mhd_paths
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
      table_paths%fermi = fermi_path
      table_paths%scv_h = scvh_path
      table_paths%scv_he = scvhe_path
      table_paths%scv_z = scvz_path
      table_paths%opal = opal06_path
      call eos_init(settings, table_paths, gerr)
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

! The Yale/SCV branch at the same four points (it once asked blend_eos
! to skip the Saha solver's derivatives, so its two cool points read
! unset values).
      write(*,'(a)') "# test_eos: eos_get_gamma1, Yale/SCV branch"
      settings = eos_settings
      settings%use_opal2006_eos = .false.
      call set_eos_settings(settings)
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
         mhd_paths = table_paths
         mhd_paths%mhd_zams = [character(len=256) :: trim(mhd_dir)//"/zams_a", &
              trim(mhd_dir)//"/zams_b", trim(mhd_dir)//"/zams_c"]
         mhd_paths%mhd_centre = [character(len=256) :: trim(mhd_dir)//"/centre1", &
              trim(mhd_dir)//"/centre2", trim(mhd_dir)//"/centre3", &
              trim(mhd_dir)//"/centre4", trim(mhd_dir)//"/centre5"]
         call eos_init(settings, mhd_paths, gerr)
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
! through compute_blended_eos/blend_eos/compute_opal_eos/esac06. The second drives esac06 (a
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
      call check_physics()

      close(run_log_unit)
      if (failures > 0) then
         write(*,'(i0,a)') failures, " check(s) failed"
         stop 1
      end if
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
         type(eos_table_paths_t) :: paths
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
         paths = table_paths
         paths%opal = ""
         if (present(file)) paths%opal = file
         call eos_init(s, paths, ierr)
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
      type(eos_table_paths_t) :: paths
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
      paths = table_paths
      paths%mhd_zams = zams_names
      paths%mhd_centre = centre_names
      call eos_init(s, paths, ierr)
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
      paths%mhd_zams(2) = "mhd_missing"
      call eos_init(s, paths, ierr)
      call check("mhd reader: a missing table returns ierr = 1", ierr == 1)
      open(newunit=u, file="mhd_truncated", form="unformatted", status="replace")
      write(u) mhd_ivarc_test, 0, 0, 0.0d0
      close(u)
      paths%mhd_zams(2) = "mhd_truncated"
      call eos_init(s, paths, ierr)
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

! ---------------------------------------------------------------
! The EOS against the physics it must reproduce (ok/FAIL), each check
! confirmed to fail under a mutation of the code it covers:
!   ideal_gas_checks       the limits of the Yale solvers: a neutral
!                          gas and a weakly degenerate fully ionized one
!                          (rho = mu P_gas/(R T), with the first-order
!                          Fermi-Dirac correction), grad_ad and cp
!                          against the ideal gas plus radiation, and
!                          grad_ad -> 1/4 when radiation dominates
!   continuity_checks      rho and grad_ad across the edges of the
!                          Saha-to-fully-ionized ramp (cutoff - 0.5 and
!                          the cutoff) and across the OPAL 2006 table's
!                          edge ramp
!   atmosphere_checks      eos_get with in_atmosphere = .true.
!   envelope_table_checks  eos_update_envelope_tables after a change of
!                          the envelope composition
!   debye_huckel_checks    the Debye-Hueckel correction and its metal
!                          mixture (eos_set_debye_huckel_z)
!   gamma1_checks          eos_get_gamma1 outside the OPAL table
!   diffusion_z_checks     the diffusion-of-Z mixture against the default
!                          one where Z is the envelope's
subroutine check_physics()
      write(*,'(a)') "# test_eos: physics checks"
      call ideal_gas_checks()
      call continuity_checks()
      call atmosphere_checks()
      call envelope_table_checks()
      call debye_huckel_checks()
      call gamma1_checks()
      call diffusion_z_checks()
end subroutine check_physics

! One eos_get at (log T, log P) and (x, z), from a fresh record.
subroutine eos_at(lt, lp, x, z, res, ierr, atmosphere, comp)
      real(dp), intent(in) :: lt, lp, x, z
      type(eos_result), intent(out) :: res
      integer, intent(out) :: ierr
      logical, intent(in), optional :: atmosphere
      real(dp), intent(in), optional :: comp(:)
      integer :: iterations
      logical :: atm
      atm = .false.
      if (present(atmosphere)) atm = atmosphere
      res = eos_result()
      iterations = 0
      call eos_get(lt, lp, x, z, res, .true., atm, iterations, &
           composition_at_zone=comp, ierr=ierr)
end subroutine eos_at

! cp and grad_ad of an ideal gas plus radiation with gas-pressure
! fraction b and R/mu = r_mu: with delta = (4 - 3b)/b and alpha = 1/b,
! cp = (R/mu)(3/2 + 12(1 - b)/b + delta^2) and grad_ad = (R/mu) delta/(b cp)
! (Kippenhahn & Weigert sec. 13.2)
! With use_diffusion_z the metals scale with Z/Z_env; where Z is the
! envelope's and only X differs (X = 0.65 against 0.70), the zone's
! 1/mu_ion is the default mixture's (the same H and He terms), on the
! Saha path (log T = 4.5) and the fully ionized one (log T = 6.5).
subroutine diffusion_z_checks()
      real(dp), parameter :: x = 0.65d0, lt(2) = [4.5d0, 6.5d0], lp(2) = [9.0d0, 15.0d0]
      type(eos_result) :: plain, diffusion
      integer :: k, ierr1, ierr2
      logical :: ok
      ok = .true.
      do k = 1, 2
         call configure_eos("yale_only")
         call eos_at(lt(k), lp(k), x, envelope_metal_fraction, plain, ierr1)
         call configure_eos("yale_diffusion_z", diffusion_z=.true.)
         call eos_at(lt(k), lp(k), x, envelope_metal_fraction, diffusion, ierr2)
         ok = ok .and. ierr1 == 0 .and. ierr2 == 0 .and. &
              abs(diffusion%ion_mean_weight_inverse/plain%ion_mean_weight_inverse - 1.0d0) < 1.0d-14 &
              .and. abs(diffusion%log10_density - plain%log10_density) < 1.0d-12
      end do
      call check("diffusion of Z: at the envelope's Z, the default mixture's 1/mu_ion and rho", ok)
end subroutine diffusion_z_checks

! Below the OPAL 2006 table (T = 1500 K) eos_get_gamma1 takes the
! Yale/SCV value (to 1e-6: the Saha solver starts from the previous
! call's solution), with OPAL selected as without it.
subroutine gamma1_checks()
      type(eos_settings_t) :: s
      real(dp) :: g1_yale, ga_yale, g1, ga
      integer :: ierr1, ierr2, iterations
      call configure_eos("production", opal2006=.true., scv=.true., file=opal06_path)
      s = eos_settings
      s%use_opal2006_eos = .false.
      call set_eos_settings(s)
      iterations = 0
      call eos_get_gamma1(0.70d0, 0.016492d0, 0.0015d0, 1.0d-9, 1.0d3, g1_yale, ga_yale, &
           iterations, ierr=ierr1)
      s%use_opal2006_eos = .true.
      call set_eos_settings(s)
      iterations = 0
      call eos_get_gamma1(0.70d0, 0.016492d0, 0.0015d0, 1.0d-9, 1.0d3, g1, ga, &
           iterations, ierr=ierr2)
      call check("eos_get_gamma1: below the OPAL table, the Yale/SCV value", &
           ierr1 == 0 .and. ierr2 == 0 .and. abs(g1/g1_yale - 1.0d0) < 1.0d-6 .and. &
           abs(ga/ga_yale - 1.0d0) < 1.0d-6)
end subroutine gamma1_checks

subroutine ideal_radiation_gas(b, r_mu, cp, grad_ad)
      real(dp), intent(in) :: b, r_mu
      real(dp), intent(out) :: cp, grad_ad
      real(dp) :: delta
      delta = (4.0d0 - 3.0d0*b)/b
      cp = r_mu*(1.5d0 + 12.0d0*(1.0d0 - b)/b + delta*delta)
      grad_ad = r_mu*delta/(b*cp)
end subroutine ideal_radiation_gas

subroutine ideal_gas_checks()
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0
! weakly degenerate, fully ionized points (eta from -4.9 to -2.6)
      real(dp), parameter :: fi_lt(4) = [6.0d0, 6.3d0, 6.5d0, 7.0d0]
      real(dp), parameter :: fi_lp(4) = [13.0d0, 14.3d0, 15.0d0, 16.5d0]
      type(eos_result) :: r
      real(dp) :: mu_e_inv, mu_inv, rho_ideal, correction, cp, grad, worst
      integer :: k, ierr
      logical :: ok
      call configure_eos("yale_only")
! a neutral gas (log T = 3.5): rho = mu P_gas/(R T) with 1/mu = 1/mu_ion +
! 1/mu_e, the ions' from the envelope mixture and few free electrons;
! grad_ad and cp of an ideal monatomic gas plus its little radiation,
! but for the metals' partial ionization (1e-3)
      call eos_at(3.5d0, 3.0d0, x, z, r, ierr)
      rho_ideal = r%beta*r%pressure/(gas_constant*(mean_atomic_mass + &
           r%electron_mean_weight_inverse)*r%temperature)
      call ideal_radiation_gas(r%beta, r%specific_gas_constant, cp, grad)
      call check("ideal gas, neutral: rho = mu P_gas/(R T)", ierr == 0 .and. &
           abs(r%density/rho_ideal - 1.0d0) < 1.0d-12 .and. &
           r%electron_mean_weight_inverse < 1.0d-4*mean_atomic_mass)
      call check("ideal gas, neutral: grad_ad -> 0.4, cp -> 5R/(2 mu) (+ radiation)", ierr == 0 .and. &
           abs(r%adiabatic_gradient - 0.4d0) < 1.0d-3 .and. &
           abs(r%adiabatic_gradient/grad - 1.0d0) < 2.0d-3 .and. &
           abs(r%specific_heat_cp/cp - 1.0d0) < 2.0d-3)
! fully ionized: rho = mu beta P/(R T) (1 - 2^(-5/2) e^eta P_e/P_gas),
! the weak-degeneracy expansion of the electron pressure, with 1/mu_e
! from the composition (X/1.008 + (1 - X)/2.002)
      mu_e_inv = x/1.008d0 + (1.0d0 - x)/2.002d0
      mu_inv = mean_atomic_mass + mu_e_inv
      worst = 0.0d0
      ok = .true.
      do k = 1, size(fi_lt)
         call eos_at(fi_lt(k), fi_lp(k), x, z, r, ierr)
         ok = ok .and. ierr == 0
         correction = 2.0d0**(-2.5d0)*exp(r%electron_degeneracy_parameter)*mu_e_inv/mu_inv
         rho_ideal = r%beta*r%pressure/(gas_constant*mu_inv*r%temperature)*(1.0d0 - correction)
         worst = max(worst, abs(r%density/rho_ideal - 1.0d0))
      end do
      call check("ideal gas, fully ionized: rho = mu P_gas/(R T) with the Fermi-Dirac correction", &
           ok .and. worst < 2.0d-4)
! grad_ad and cp against the ideal gas plus radiation: nearly ideal
! (beta = 0.9997, eta = -4.9) and radiation-dominated (beta = 0.1)
      call eos_at(6.0d0, 13.0d0, x, z, r, ierr)
      call ideal_radiation_gas(r%beta, r%specific_gas_constant, cp, grad)
      call check("ideal gas, fully ionized: grad_ad -> 0.4, cp -> 5R/(2 mu) (+ radiation)", &
           ierr == 0 .and. abs(r%adiabatic_gradient/grad - 1.0d0) < 5.0d-4 .and. &
           abs(r%specific_heat_cp/cp - 1.0d0) < 2.0d-3 .and. abs(r%adiabatic_gradient - 0.4d0) < 1.0d-3)
      call eos_at(7.5d0, 15.45d0, x, z, r, ierr)
      call ideal_radiation_gas(r%beta, r%specific_gas_constant, cp, grad)
      call check("radiation-dominated: grad_ad and cp of the gas plus radiation", &
           ierr == 0 .and. r%beta < 0.2d0 .and. abs(r%adiabatic_gradient/grad - 1.0d0) < 1.0d-4 .and. &
           abs(r%specific_heat_cp/cp - 1.0d0) < 1.0d-4)
      call eos_at(7.5d0, 15.41d0, x, z, r, ierr)
      call check("radiation-dominated: grad_ad -> 1/4", ierr == 0 .and. r%beta < 0.03d0 .and. &
           abs(r%adiabatic_gradient - 0.25d0) < 1.0d-3)
! radiation pressure above P (log P_rad = 15.402 at log T = 7.5): no gas
! pressure, an error
      call eos_at(7.5d0, 15.40d0, x, z, r, ierr)
      call check("radiation pressure above P: beta <= 0 returns ierr", ierr == 1)
end subroutine ideal_gas_checks

! rho, grad_ad, the ionization fractions and d(dln rho/dln T)/dln T and
! /dln P on both sides of a boundary in log T (fixed log P): the jump
! over 2e-7 in log T, against 1e-6 (1e-4 for the second derivatives,
! which carry the fully ionized iteration's tolerance; a step at the
! boundary would be the gap between the two solutions, 1e-3 to 1e-1)
subroutine continuity_in_t(label, lt0, lp)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: lt0, lp
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0, h = 1.0d-7
      type(eos_result) :: below, above
      integer :: ierr1, ierr2
      call eos_at(lt0 - h, lp, x, z, below, ierr1)
      call eos_at(lt0 + h, lp, x, z, above, ierr2)
      call check(label, ierr1 == 0 .and. ierr2 == 0 .and. &
           abs(above%log10_density - below%log10_density) < 1.0d-6 .and. &
           abs(above%adiabatic_gradient - below%adiabatic_gradient) < 1.0d-6 .and. &
           all(abs(above%ion_fraction - below%ion_fraction) < 1.0d-6) .and. &
           abs(above%dlnrho_dlnt_dt - below%dlnrho_dlnt_dt) < 1.0d-4 .and. &
           abs(above%dlnrho_dlnt_dp - below%dlnrho_dlnt_dp) < 1.0d-4)
end subroutine continuity_in_t

subroutine continuity_checks()
      call configure_eos("yale_only")
      call continuity_in_t("continuity: the blended quantities at the Saha ramp's foot (cutoff - 0.5)", &
           5.5d0, 13.0d0)
      call continuity_in_t("continuity: the blended quantities at the Saha cutoff", 6.0d0, 13.0d0)
      call continuity_in_t("continuity: the blended quantities at the Saha cutoff, denser", 6.0d0, 15.0d0)
! the OPAL 2006 table's edge at log T = 5, in log P: the blend is
! OPAL's own result below log P = 12.87, the Yale/SCV result above
! 13.01, and ramps between. Each end of the ramp is bisected (to 1e-11
! in log P) and the two sides compared: rho to 1e-8 and grad_ad to 1e-7
! in the log (the gap between the two results there is 1e-3 to 1e-2).
      call configure_eos("production", opal2006=.true., scv=.true., file=opal06_path)
      call check("continuity: rho and grad_ad where the OPAL 2006 edge ramp reaches OPAL", &
           ramp_end_is_continuous(12.80d0, 12.94d0, .true.))
      call check("continuity: rho and grad_ad at the OPAL 2006 table's edge (the Yale/SCV result beyond)", &
           ramp_end_is_continuous(13.05d0, 12.94d0, .false.))
end subroutine continuity_checks

! At log T = 5 between log P = lp_equal (where the blend equals the
! reference) and lp_differ (inside the ramp, where it does not), bisect
! for the end of the OPAL 2006 edge ramp and compare rho and grad_ad on
! its two sides. The reference is OPAL's own result (compute_opal_eos)
! at the ramp's inner end (inner_end), and the Yale/SCV result (OPAL
! switched off) at the table's edge.
logical function ramp_end_is_continuous(lp_equal, lp_differ, inner_end)
      real(dp), intent(in) :: lp_equal, lp_differ
      logical, intent(in) :: inner_end
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0, lt = 5.0d0
      real(dp) :: a, b, m
      type(eos_result) :: below, above
      integer :: k, ierr1, ierr2
      logical :: ok
      a = lp_equal
      b = lp_differ
      ok = blend_equals_reference(a, inner_end) .and. .not. blend_equals_reference(b, inner_end)
      do k = 1, 40
         m = 0.5d0*(a + b)
         if (blend_equals_reference(m, inner_end)) then
            a = m
         else
            b = m
         end if
      end do
      call eos_at(lt, a, x, z, below, ierr1)
      call eos_at(lt, b, x, z, above, ierr2)
      ramp_end_is_continuous = ok .and. ierr1 == 0 .and. ierr2 == 0 .and. &
           abs(above%log10_density - below%log10_density) < 1.0d-8 .and. &
           abs(above%adiabatic_gradient - below%adiabatic_gradient) < 1.0d-7
end function ramp_end_is_continuous

! Whether the production blend at (log T = 5, lp) equals OPAL's own
! result (opal_reference) or the Yale/SCV result, to 1e-12 in log rho.
logical function blend_equals_reference(lp, opal_reference)
      real(dp), intent(in) :: lp
      logical, intent(in) :: opal_reference
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0, lt = 5.0d0
      type(eos_result) :: blend, reference
      type(eos_settings_t) :: s
      integer :: ierr
      logical :: covers
      call eos_at(lt, lp, x, z, blend, ierr)
      if (opal_reference) then
         reference = eos_result()
         call compute_opal_eos(iv_opal06, lt, lp, x, z, reference, ierr, covers)
      else
         s = eos_settings
         s%use_opal2006_eos = .false.
         call set_eos_settings(s)
         call eos_at(lt, lp, x, z, reference, ierr)
         s%use_opal2006_eos = .true.
         call set_eos_settings(s)
      end if
      blend_equals_reference = abs(blend%log10_density - reference%log10_density) < 1.0d-12
end function blend_equals_reference

! in_atmosphere = .true. at Saha points, in the ramp and below it, with
! and without SCV: the same density, beta, ionization fractions, first
! derivatives, cp and grad_ad as a full evaluation (to 1e-10, and the
! derivatives to 1e-8: the Saha solver starts from the previous call's
! solution)
subroutine atmosphere_checks()
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0
      real(dp), parameter :: lt(3) = [3.7d0, 4.6d0, 5.8d0], lp(3) = [4.5d0, 9.0d0, 12.6d0]
      type(eos_result) :: full, atm
      integer :: k, c, ierr1, ierr2
      logical :: ok
      ok = .true.
      do c = 1, 2
         if (c == 1) call configure_eos("yale_only")
         if (c == 2) call configure_eos("scv_only", scv=.true.)
         do k = 1, size(lt)
            call eos_at(lt(k), lp(k), x, z, full, ierr1)
            call eos_at(lt(k), lp(k), x, z, atm, ierr2, atmosphere=.true.)
            ok = ok .and. ierr1 == 0 .and. ierr2 == 0 .and. &
                 abs(atm%log10_density - full%log10_density) < 1.0d-10 .and. &
                 atm%beta == full%beta .and. &
                 all(abs(atm%ion_fraction - full%ion_fraction) < 1.0d-10) .and. &
                 abs(atm%dlnrho_dlnt/full%dlnrho_dlnt - 1.0d0) < 1.0d-8 .and. &
                 abs(atm%dlnrho_dlnp/full%dlnrho_dlnp - 1.0d0) < 1.0d-8 .and. &
                 abs(atm%specific_heat_cp/full%specific_heat_cp - 1.0d0) < 1.0d-8 .and. &
                 abs(atm%adiabatic_gradient/full%adiabatic_gradient - 1.0d0) < 1.0d-8
         end do
      end do
      call check("in_atmosphere: the density, ionization and first derivatives of a full evaluation", ok)
end subroutine atmosphere_checks

! After the envelope composition changes (eos_set_mixture) and
! eos_update_envelope_tables rebuilds the SCV envelope table, the table
! at the new composition agrees with SCV's general-composition path
! there (to 1e-4 in log rho; a table left at the old composition is off
! by 0.05)
subroutine envelope_table_checks()
      real(dp), parameter :: x_new = 0.60d0, z_new = 0.03d0
      real(dp), parameter :: lt(3) = [3.7d0, 3.9d0, 4.1d0], lp(3) = [5.0d0, 6.5d0, 8.0d0]
      type(eos_result) :: table, general
      integer :: k, ierr1, ierr2
      logical :: ok
      call configure_eos("scv_only", scv=.true.)
      call set_envelope_mixture(x_new, z_new)
      call eos_update_envelope_tables(x_new, z_new)
      ok = .true.
      do k = 1, size(lt)
         call eos_at(lt(k), lp(k), x_new, z_new, table, ierr1)
         call eos_at(lt(k), lp(k), x_new + 2.0d-5, z_new, general, ierr2)
         ok = ok .and. ierr1 == 0 .and. ierr2 == 0 .and. &
              abs(table%log10_density - general%log10_density) < 1.0d-4
      end do
      call check("envelope tables: rebuilt at a new composition, the SCV table agrees with SCV's general path", ok)
      call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
      call eos_update_envelope_tables(envelope_hydrogen_fraction, envelope_metal_fraction)
end subroutine envelope_table_checks

! The envelope mixture of (x, z) by the main program's algorithm, from
! the default seed or another, handed to eos_set_mixture.
subroutine set_envelope_mixture(x, z, seed)
      real(dp), intent(in) :: x, z
      real(dp), intent(in), optional :: seed(:)
      real(dp) :: wm(12), ws, sc, mam, fractions(12)
      integer :: j
      wm = mixture_weights_seed
      if (present(seed)) wm = seed
      ws = wm(1)+wm(2)+wm(3)+wm(4)+wm(5)+wm(6)+wm(8)+wm(9)+wm(10)+wm(11)
      sc = z/ws
      wm(7) = x/sc
      wm(12) = (1.0d0 - x - z)/sc
      ws = 0.0d0
      do j = 1, 12
         wm(j) = sc*wm(j)/atomic_weight(j)
         ws = ws + wm(j)
      end do
      mam = ws
      fractions = wm/mam
      call eos_set_mixture(x, z, mam, fractions)
end subroutine set_envelope_mixture

! The Debye-Hueckel correction lowers the pressure at a given density,
! so at a given pressure it raises the density, the more the higher
! the metals' charges: its rise of rho is larger when the envelope
! mixture's metals beyond C, N and O are all Fe (eos takes those metals
! from the envelope mixture, eos_set_mixture).
subroutine debye_huckel_checks()
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0, lt = 6.8d0, lp = 17.0d0
      real(dp) :: comp(15), seed(12), raise(2)
      type(eos_result) :: none, with_dh
      integer :: ierr1, ierr2, k
      logical :: ok
      comp = 0.0d0
      comp(1) = x
      comp(2) = 1.0d0 - x - z
      comp(3) = z
      comp(5) = 0.17d0*z
      comp(7) = 0.05d0*z
      comp(9) = 0.48d0*z
      seed = mixture_weights_seed
      seed([1, 2, 3, 5, 10, 11]) = 0.0d0
      seed(4) = sum(mixture_weights_seed([1, 2, 3, 4, 5, 10, 11]))
      ok = .true.
      do k = 1, 2
         if (k == 1) call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
         if (k == 2) call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction, seed)
         call configure_eos("yale_only")
         call eos_at(lt, lp, x, z, none, ierr1)
         call configure_eos("debye_huckel", debye_huckel=.true.)
         call eos_at(lt, lp, x, z, with_dh, ierr2, comp=comp)
         ok = ok .and. ierr1 == 0 .and. ierr2 == 0
         raise(k) = with_dh%density/none%density - 1.0d0
      end do
      call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
      call check("Debye-Hueckel: it raises rho, the more with the envelope's heavier metals as Fe", &
           ok .and. raise(1) > 0.0d0 .and. raise(2) > raise(1))
end subroutine debye_huckel_checks

subroutine check(label, ok)
      character(len=*), intent(in) :: label
      logical, intent(in) :: ok
      write(*,'(a,1x,a)') merge("ok  ", "FAIL", ok), label
      if (.not. ok) failures = failures + 1
end subroutine check

end program test_eos
