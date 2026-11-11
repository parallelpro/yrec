!----------------------------------------------------------------------
! test_eos
!----------------------------------------------------------------------
! The equation of state, through eos_lib, without a star: the program
! installs the settings and the tables itself (eos_init), sets the
! envelope mixture a run would (eos_set_mixture: core's mixture
! algorithm on the default mixture seed, X = 0.70, Z = 0.016492, the Z
! of the shipped OPAL 2006 table), which builds the SCV envelope table
! once the SCV tables are read. Its output is byte-compared with
! expected_test_eos.out (regenerated only on purpose). In order:
!   the production EOS (SCV and OPAL 2006) over a (log T, log P) grid;
!   gamma1_at_density (gamma1 and grad_ad at T, rho) on its two branches; the MHD EOS when YREC_MHD_TABLES
!   names a directory with the real tables (none ship with YREC: SKIPPED
!   otherwise); the error paths; every configuration a deck can choose
!   (check_configurations); the derivatives against finite differences
!   (check_derivatives); the MHD table lookup on a synthetic table; then
!   the ok/FAIL checks (the first, check_scv_general_path, before the
!   others: SCV's general path before and after the envelope table is
!   built, which needs the program's first evaluations):
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
      use run_messages, only: run_log_unit
      use eos_lib
      use eos_configuration, only: eos_settings, set_eos_settings
      use blended_eos, only: gamma1_at_density
      use eos_mixture, only: debye_huckel_composition_of
      use yale_eos_tables, only: debye_huckel_composition, idh_c, idh_n, idh_o, idh_ne, &
           idh_na, idh_mg, idh_al, idh_si, idh_ar, idh_fe
      use opal_eos_tables, only: opal_tables, opal95, opal_shared, i_opal_cv, i_opal01_cv, &
           i_opal_p, i_opal01_p
      use opal_eos, only: esac06, esac01, compute_opal_eos, iv_opal01, iv_opal06
      use opal95_eos, only: esac
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
           inverse_atomic_weight, number_fractions(12)

      double precision :: w(12), wsum, scale
      integer :: i, ipt, eos_ierr

! the point the sections evaluate, and its result (eval_eos)
      double precision :: logt, logp
      double precision :: x_frac, z_frac
      type(eos_result) :: ev
      logical :: lderiv, latmo, in_table
      integer :: saha_iterations

! SCV's general path at a cool point before the envelope table is
! built, and the Saha result there (scv_general_path_before_table)
      type(eos_result) :: general_before_table(2)

! the derivative check's points, and the largest relative difference of
! the first and of the second derivatives at each (derivative_points)
      integer, parameter :: nd = 6
      double precision :: worst_first(nd), worst_second(nd)

! gamma1_at_density argument set
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

! a table source's reference result at (log T, log P, X, Z), for
! ramp_end_is_continuous
      abstract interface
         subroutine reference_eos(lt, lp, x, z, r)
            import :: dp, eos_result
            real(dp), intent(in) :: lt, lp, x, z
            type(eos_result), intent(out) :: r
         end subroutine reference_eos
      end interface

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
      inverse_atomic_weight = wsum
      scale = 1.0d0/inverse_atomic_weight
      do i = 1, 12
         number_fractions(i) = w(i)*scale
      end do
      call eos_set_mixture(envelope_hydrogen_fraction, envelope_metal_fraction, &
           inverse_atomic_weight, number_fractions)

! the production EOS (SCV and OPAL 2006, the Saha cutoff at log T = 6,
! as in the reference solar decks), installed with the tables it reads;
! the OPAL table is read by esac06 on first use, as in a run
      settings%use_scv_eos = .true.
      settings%opal_vintage = opal_vintage_2006
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
! SCV's general path before the envelope table exists (a run's first
! card reaches it so), then the mixture again, which now builds the SCV
! envelope table (the first eos_set_mixture came before eos_init read
! the SCV tables), as a run's starting model does
      call scv_general_path_before_table(general_before_table)
      call eos_set_mixture(envelope_hydrogen_fraction, envelope_metal_fraction, &
           inverse_atomic_weight, number_fractions)
      call check_scv_general_path(general_before_table)

      x_frac = 0.70d0
      z_frac = 0.016492d0

      write(*,'(a)') "# test_eos: eos_eval over (logT, logP) grid, " // &
           "X=0.70 Z=0.016492, SCV+OPAL2006, no MHD, no DH"
      do ipt = 1, npts
         logt = grid_logt(ipt)
         logp = grid_logp(ipt)
         ev = eos_result()
         saha_iterations = 0
         lderiv = .true.
         latmo = .false.
         call eval_eos(x_frac, z_frac, lderiv, latmo, gerr)
      if (gerr /= 0) stop 1
         write(*,'(a,i2)') "point ", ipt
         write(*,'(4(1pe24.15))') logt, logp, ev%log10_density, ev%beta
         write(*,'(4(1pe24.15))') ev%dlnrho_dlnt, ev%dlnrho_dlnp, ev%specific_heat_cp, &
              ev%adiabatic_gradient
         write(*,'(4(1pe24.15))') ev%specific_gas_constant, ev%ion_mean_weight_inverse, &
              ev%electron_mean_weight_inverse, ev%electron_degeneracy_parameter
         write(*,'(3(1pe24.15))') ev%ion_fraction(1), ev%ion_fraction(2), ev%ion_fraction(3)
         write(*,'(4(1pe24.15))') ev%dlnrho_dlnt_dt, ev%dlnrho_dlnt_dp, ev%adiabatic_gradient_dt, &
              ev%adiabatic_gradient_dp
         write(*,'(2(1pe24.15),i4)') ev%specific_heat_cp_dt, ev%specific_heat_cp_dp, saha_iterations
      end do

      write(*,'(a)') "# test_eos: gamma1_at_density, OPAL2006 branch"
      do ipt = 1, ng1
         t6 = g1_t6(ipt)
         rho = g1_rho(ipt)
         pres = g1_p(ipt)
         saha_iterations = 0
         call gamma1_at_density(x_frac, z_frac, t6, rho, pres, gamma1, &
              grad_ad, saha_iterations, ierr=gerr)
      if (gerr /= 0) stop 1
         write(*,'(a,i2,2(1pe24.15))') "g1 ", ipt, gamma1, grad_ad
      end do

! The Yale/SCV branch at the same four points (it once asked blend_eos
! to skip the Saha solver's derivatives, so its two cool points read
! unset values).
      write(*,'(a)') "# test_eos: gamma1_at_density, Yale/SCV branch"
      settings = eos_settings
      settings%opal_vintage = opal_vintage_none
      call set_eos_settings(settings)
      do ipt = 1, ng1
         t6 = g1_t6(ipt)
         rho = g1_rho(ipt)
         pres = g1_p(ipt)
         saha_iterations = 0
         call gamma1_at_density(x_frac, z_frac, t6, rho, pres, gamma1, &
              grad_ad, saha_iterations, ierr=gerr)
      if (gerr /= 0) stop 1
         write(*,'(a,i2,2(1pe24.15))') "g1 ", ipt, gamma1, grad_ad
      end do
      settings%opal_vintage = opal_vintage_2006
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
            ev = eos_result()
            saha_iterations = 0
            lderiv = .true.
            latmo = .false.
            call eval_eos(x_frac, z_frac, lderiv, latmo, gerr)
      if (gerr /= 0) stop 1
            write(*,'(a,i2)') "mhd point ", ipt
            write(*,'(4(1pe24.15))') logt, logp, ev%log10_density, ev%beta
            write(*,'(4(1pe24.15))') ev%dlnrho_dlnt, ev%dlnrho_dlnp, ev%specific_heat_cp, &
                 ev%adiabatic_gradient
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
      ev = eos_result()
      saha_iterations = 0
      lderiv = .true.
      latmo = .false.
      call eval_eos(x_frac, z_frac, lderiv, latmo, eos_ierr)
      write(*,'(a,i4)') "err facade-success ierr = ", eos_ierr
      call esac06(opal_tables(iv_opal06), 0.70d0, 2.0d0, 0.1d0, 9, 2, eos_ierr, in_table)
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
! grad_ad, eta and the Saha iteration count.
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

! SCV's general path (a composition off the envelope table's) at a cool
! point SCV covers, before the envelope table is built: kept in
! results(1), with the Saha result there (SCV off) in results(2).
      subroutine scv_general_path_before_table(results)
         type(eos_result), intent(out) :: results(2)
         type(eos_settings_t) :: s, before
         integer :: ierr
         before = eos_settings
         s = eos_settings
         s%opal_vintage = opal_vintage_none
         call set_eos_settings(s)
         call eos_at(3.7d0, 5.0d0, 0.70d0 + 2.0d-5, 0.016492d0, results(1), ierr)
         s%use_scv_eos = .false.
         call set_eos_settings(s)
         call eos_at(3.7d0, 5.0d0, 0.70d0 + 2.0d-5, 0.016492d0, results(2), ierr)
         call set_eos_settings(before)
      end subroutine scv_general_path_before_table

! SCV's general path reads the three tables, not the envelope table: before
! that table is built it covers the point (it is not the Saha result,
! 9e-4 dex away there) and gives what it gives after the build.
      subroutine check_scv_general_path(before_table)
         type(eos_result), intent(in) :: before_table(2)
         type(eos_result) :: after_table(2)
         call scv_general_path_before_table(after_table)
         call check("SCV general path: before the envelope table is built, SCV's result, as after", &
              abs(before_table(1)%log10_density - before_table(2)%log10_density) > 1.0d-4 .and. &
              abs(before_table(1)%log10_density - after_table(1)%log10_density) < 1.0d-12)
      end subroutine check_scv_general_path

! Set one configuration (the controls a deck would set) with eos_init,
! the OPAL table of the vintage when one is named; with SCV, the
! program's envelope mixture again, which builds the envelope table, as
! a run's starting model does after eos_init (which reads the SCV tables
! anew).
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
         s%opal_vintage = opal_vintage_none
         if (present(opal95)) s%opal_vintage = opal_vintage_1995
         if (present(opal2001)) s%opal_vintage = opal_vintage_2001
         if (present(opal2006)) s%opal_vintage = opal_vintage_2006
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
         if (s%use_scv_eos) call set_envelope_mixture(envelope_hydrogen_fraction, &
              envelope_metal_fraction)
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
            ev = eos_result()
            saha_iterations = 0
            if (present(with_composition)) then
               call eval_eos(x, z, .true., .false., ierr, comp)
            else
               call eval_eos(x, z, .true., .false., ierr)
            end if
            if (ierr /= 0) then
               write(*,'(a,2f7.3,a,i0)') "eos logT logP", logt, logp, " ierr ", ierr
            else
               write(*,'(a,2f7.3,4(1pe23.15))') "eos", logt, logp, ev%log10_density, ev%beta, &
                    ev%dlnrho_dlnt, ev%dlnrho_dlnp
               write(*,'(a,14x,3(1pe23.15),i3)') "   ", ev%specific_heat_cp, ev%adiabatic_gradient, &
                    ev%electron_degeneracy_parameter, saha_iterations
            end if
         end do
      end subroutine evaluate_grid

! eos_get at the program's (log T, log P) and (x, z), into its record ev
! (which the caller sets, eos_result() for a fresh one)
      subroutine eval_eos(x, z, want, atm, ierr, comp)
         double precision, intent(in) :: x, z
         logical, intent(in) :: want, atm
         integer, intent(out) :: ierr
         double precision, intent(in), optional :: comp(:)
         call eos_get(logt, logp, x, z, ev, want, atm, saha_iterations, &
              composition_at_zone=comp, ierr=ierr)
      end subroutine eval_eos

! The derivative check: at points in each regime, every configuration's
! analytic derivatives against centred finite differences of its own
! values (a step of 1e-4 in log T and in log P). Each row prints the
! analytic value and the finite difference, to five digits. The
! second-derivative outputs are checked as
!   dlnrho_dlnt_dt   d(dln rho/dln T)/dln T
!   dlnrho_dlnt_dp   d(dln rho/dln T)/dln P
!   adiabatic_gradient_dt, _dp   dln grad_ad/dln T, dln grad_ad/dln P
!   specific_heat_cp_dt, _dp     dln cp/dln T, dln cp/dln P
! Their relative difference, (analytic - fd)/|fd|, is judged, not
! printed. Where the derivatives are analytic it is the finite
! difference's own error, truncation and round-off amplified by 1/h
! (5e-12 to 7e-7 at the Saha points), and its leading digits differ
! between platforms (grad_ad_dp at log T 3.9: -1.4137e-8 on macOS,
! -1.4133e-8 on Linux). It is checked, ok or FAIL, where the two must
! agree: the Saha solution's (log T 3.9 to 5.5), also with the
! Debye-Hueckel correction on (which acts on the fully ionized gas);
! the fully ionized gas's (6.9 and 7.161, through the Fermi table's
! interpolation); the OPAL tables' first derivatives (through their
! interpolation; the numerical configuration's are OPAL 2006's); and
! the numerical second derivatives. Elsewhere no agreement is expected,
! and the two printed columns pin what there is: table interpolation
! adds noise to finite differences, the second derivatives inside a
! table are the Yale solution's rather than the table's, the
! Debye-Hueckel gas's eta ramp is not differentiated, and log T = 5.8
! lies in the ramp between the Saha and the fully ionized solutions
! (saha_log10t_cutoff - 0.5 to the cutoff), which blends derivatives
! rather than differentiating the blend.
      subroutine check_derivatives()
         character(len=256) :: legacy
         legacy = trim(yrec_input)//"/eos/legacy/"
         write(*,'(a)') "# test_eos: derivative check (analytic, finite difference; ok/FAIL on their relative difference)"
         call configure_eos("yale_only")
         call derivative_points("yale_only")
! the Saha solver's derivatives are analytic (log T 3.9 to 5.5: 2e-8 in
! the first, 7e-7 in the second now), the fully ionized gas's carry the
! Fermi table's interpolation (6.9 and 7.161: 6e-4 and 0.04)
         call check("derivatives: the Saha solution's first and second, against finite differences", &
              all(worst_first(1:3) < 1.0d-6) .and. all(worst_second(1:3) < 1.0d-4))
         call check("derivatives: the fully ionized gas's first (1e-3) and second (0.1)", &
              all(worst_first(5:6) < 1.0d-3) .and. all(worst_second(5:6) < 0.1d0))
         call configure_eos("scv_only", scv=.true.)
         call derivative_points("scv_only")
! the OPAL tables' first derivatives agree to 1e-3 (2006, 2001) and
! 3e-3 (1995) now; a derivative bug (a sign, a factor ln 10, the wrong
! variable) is of order one
         call configure_eos("production", opal2006=.true., scv=.true., file=opal06_path)
         call derivative_points("production")
         call check("derivatives: OPAL 2006's first, against finite differences (3e-3)", &
              all(worst_first < 3.0d-3))
         call configure_eos("opal95", opal95=.true., scv=.true., file=trim(legacy)//"opal/OPALeos.0188")
         call derivative_points("opal95")
         call check("derivatives: OPAL 1995's first, against finite differences (1e-2)", &
              all(worst_first < 1.0d-2))
         call configure_eos("opal2001", opal2001=.true., scv=.true., file=trim(legacy)//"opal2001/O01EOS0188z")
         call derivative_points("opal2001")
         call check("derivatives: OPAL 2001's first, against finite differences (3e-3)", &
              all(worst_first < 3.0d-3))
         call configure_eos("numerical", opal2006=.true., scv=.true., numerical=.true., file=opal06_path)
         call derivative_points("numerical")
         call check("derivatives: the numerical second derivatives, against finite differences (1e-4)", &
              all(worst_second < 1.0d-4))
         call configure_eos("debye_huckel", debye_huckel=.true.)
         call derivative_points("debye_huckel")
         call check("derivatives: the Saha solution's with the Debye-Hueckel correction on (1e-6, 1e-4)", &
              all(worst_first(1:3) < 1.0d-6) .and. all(worst_second(1:3) < 1.0d-4))
      end subroutine check_derivatives

! Each point's analytic derivatives against centred differences, one
! row each (the two values, not their relative difference: see
! check_derivatives); the largest relative difference of the first and
! of the second derivatives at each point stays in
! worst_first/worst_second for the checks.
         subroutine derivative_points(label)
            character(len=*), intent(in) :: label
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
               worst_first(k) = 0.0d0
               worst_second(k) = 0.0d0
               call row("dlnrho_dlnt", c(2), (tp(1) - tm(1))/(2.0d0*h), worst_first(k))
               call row("dlnrho_dlnp", c(3), (pp(1) - pm(1))/(2.0d0*h), worst_first(k))
               call row("dlnrho_dlnt_dt", c(6), (tp(2) - tm(2))/(2.0d0*h*ln10), worst_second(k))
               call row("dlnrho_dlnt_dp", c(7), (pp(2) - pm(2))/(2.0d0*h*ln10), worst_second(k))
               call row("grad_ad_dt", c(8), (log10(tp(5)) - log10(tm(5)))/(2.0d0*h), worst_second(k))
               call row("grad_ad_dp", c(9), (log10(pp(5)) - log10(pm(5)))/(2.0d0*h), worst_second(k))
               call row("cp_dt", c(10), (log10(tp(4)) - log10(tm(4)))/(2.0d0*h), worst_second(k))
               call row("cp_dp", c(11), (log10(pp(4)) - log10(pm(4)))/(2.0d0*h), worst_second(k))
            end do
         end subroutine derivative_points

         subroutine eval_point(lt, lp, want, v)
            double precision, intent(in) :: lt, lp
            logical, intent(in) :: want
            double precision, intent(out) :: v(11)
            integer :: ierr
            logt = lt
            logp = lp
            ev = eos_result()
            saha_iterations = 0
            call eval_eos(0.70d0, 0.016492d0, want, .false., ierr)
            if (ierr /= 0) write(*,'(a,2f8.4,a,i0)') "eos logT logP", lt, lp, " ierr ", ierr
            v = [ev%log10_density, ev%dlnrho_dlnt, ev%dlnrho_dlnp, ev%specific_heat_cp, &
                 ev%adiabatic_gradient, ev%dlnrho_dlnt_dt, ev%dlnrho_dlnt_dp, &
                 ev%adiabatic_gradient_dt, ev%adiabatic_gradient_dp, ev%specific_heat_cp_dt, &
                 ev%specific_heat_cp_dp]
         end subroutine eval_point

         subroutine row(name, analytic, fd, worst)
            character(len=*), intent(in) :: name
            double precision, intent(in) :: analytic, fd
            double precision, intent(inout) :: worst
            double precision :: rel
            rel = (analytic - fd)/max(abs(fd), 1.0d-30)
            write(*,'(2x,a16,2(1pe13.4))') name, analytic, fd
            worst = max(worst, abs(rel))
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
      integer :: it, ir, i, lir_status(4)
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
              out, in_range(i), lir_status(i))
         write(*,'(a,2f6.2,a,3(1pe24.15))') "lookup at", log_p(i), log_t(i), " =", out
      end do
      call check("mhd lookup: in range inside and at the corners, not above the table", &
           all(in_range(1:3)) .and. .not. in_range(4) .and. all(lir_status(1:3) == 0))
! a table whose T rows all coincide: lir's status (coincident nodes)
! comes back (it used to be dropped)
      tlog = 4.0d0
      call interpolate_mhd_table(log_p(1), log_t(1), table, nt, nr, nv, tlog, nt, nr, w1, w2, yw, &
           out, in_range(1), lir_status(1))
      call check("mhd lookup: coincident table rows return lir's status", lir_status(1) /= 0)
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
!   continuity_checks      every blended field across the ends of the
!                          Saha-to-fully-ionized ramp (cutoff - 0.5 and
!                          the cutoff; Yale alone and with SCV) and of
!                          the OPAL 2006 table's edge ramp
!   atmosphere_checks      eos_get with in_atmosphere = .true.
!   envelope_table_checks  the SCV envelope table eos_set_mixture builds
!                          for a new envelope composition
!   debye_huckel_checks    the Debye-Hueckel correction and its metal
!                          mixture (the envelope's, eos_set_mixture; in
!                          a zone, the metals beyond O summing to its Z
!                          less its C, N and O)
!   gamma1_checks          gamma1_at_density outside the OPAL table
!   scv_mean_weight_checks R/mu and 1/mu_e where SCV gives the density
!   diffusion_z_checks     the diffusion-of-Z mixture against the default
!                          one where Z is the envelope's
!   settings_checks        eos_init's check of the OPAL vintage
!   opal_trial_call_checks the 1995 and 2001 density searches' trial calls
!   scv_reader_checks      the SCV reader: a re-read drops the envelope
!                          table; a missing or truncated table is ierr
!   table_reader_checks    the Fermi-Dirac and OPAL 2006 readers: a
!                          missing or truncated table is ierr
subroutine check_physics()
      write(*,'(a)') "# test_eos: physics checks"
      call ideal_gas_checks()
      call continuity_checks()
      call atmosphere_checks()
      call envelope_table_checks()
      call debye_huckel_checks()
      call gamma1_checks()
      call scv_mean_weight_checks()
      call diffusion_z_checks()
      call settings_checks()
      call opal_trial_call_checks()
      call scv_reader_checks()
      call table_reader_checks()
end subroutine check_physics

! The Fermi-Dirac reader (at eos_init) and the OPAL 2006 reader (on
! first use): a missing table and one cut short return ierr = 1, the
! cut table's unit closed. The OPAL table read before stays loaded
! until eos_init names another file, which its reader then reads (the
! cut copy here, then the real table back).
subroutine table_reader_checks()
      type(eos_settings_t) :: s
      type(eos_table_paths_t) :: paths
      type(eos_result) :: r
      integer :: ierr, ierr_missing, ierr_cut, ierr_card, ierr_wide, it
      logical :: cut_open
      s = eos_settings
      s%opal_vintage = opal_vintage_none
      s%use_scv_eos = .false.
      paths = table_paths
      paths%fermi = "fermi_missing"
      call eos_init(s, paths, ierr_missing)
      call copy_head(fermi_path, "fermi_truncated", 40)
      paths%fermi = "fermi_truncated"
      call eos_init(s, paths, ierr_cut)
      inquire(file="fermi_truncated", opened=cut_open)
      call check("Fermi table: a missing or truncated table returns ierr = 1, its unit closed", &
           ierr_missing == 1 .and. ierr_cut == 1 .and. .not. cut_open)
! a card whose abscissa is out of the table (card 100's 14 made 44), and
! abscissae spanning more than the X lookup's 261 bins (the first, -16,
! made -30): the reader stores neither
      call copy_editing(fermi_path, "fermi_bad_card", 110, 1, "   44")
      paths%fermi = "fermi_bad_card"
      call eos_init(s, paths, ierr_card)
      call copy_editing(fermi_path, "fermi_wide_x", 1, 9, "-30.0000")
      paths%fermi = "fermi_wide_x"
      call eos_init(s, paths, ierr_wide)
      call check("Fermi table: a card out of its place, or abscissae wider than the lookup, returns ierr = 1", &
           ierr_card == 1 .and. ierr_wide == 1)
      call copy_head(opal06_path, "opal06_truncated", 3000)
      call configure_eos("production", opal2006=.true., scv=.true., file="opal06_truncated")
      it = 0
      r = eos_result()
      call eos_get(5.0d0, 12.8d0, 0.70d0, 0.016492d0, r, .true., .false., it, ierr=ierr)
      inquire(file="opal06_truncated", opened=cut_open)
      call check("OPAL 2006 table: a truncated table returns ierr = 1 on first use, its unit closed", &
           ierr == 1 .and. .not. cut_open)
      paths = table_paths
      paths%opal = "opal06_missing"
      s = eos_settings
      call eos_init(s, paths, ierr_missing)
      call check("OPAL 2006 table: a missing table returns ierr = 1 at eos_init", ierr_missing == 1)
      call configure_eos("production", opal2006=.true., scv=.true., file=opal06_path)
      r = eos_result()
      call eos_get(5.0d0, 12.8d0, 0.70d0, 0.016492d0, r, .true., .false., it, ierr=ierr)
      if (ierr /= 0) write(*,'(a)') "table_reader_checks: the OPAL table did not read back"
end subroutine table_reader_checks

! the file src into dst, with line line_no's columns from col on
! replaced by text
subroutine copy_editing(src, dst, line_no, col, text)
      character(len=*), intent(in) :: src, dst, text
      integer, intent(in) :: line_no, col
      character(len=400) :: line
      integer :: u_in, u_out, k, ios
      open(newunit=u_in, file=src, status="old", action="read")
      open(newunit=u_out, file=dst, status="replace")
      k = 0
      do
         read(u_in, '(a)', iostat=ios) line
         if (ios /= 0) exit
         k = k + 1
         if (k == line_no) line(col:col+len(text)-1) = text
         write(u_out, '(a)') trim(line)
      end do
      close(u_in)
      close(u_out)
end subroutine copy_editing

! the first n_lines lines of the file src, into dst
subroutine copy_head(src, dst, n_lines)
      character(len=*), intent(in) :: src, dst
      integer, intent(in) :: n_lines
      character(len=400) :: line
      integer :: u_in, u_out, k, ios
      open(newunit=u_in, file=src, status="old", action="read")
      open(newunit=u_out, file=dst, status="replace")
      do k = 1, n_lines
         read(u_in, '(a)', iostat=ios) line
         if (ios /= 0) exit
         write(u_out, '(a)') trim(line)
      end do
      close(u_in)
      close(u_out)
end subroutine copy_head

! The SCV reader. After eos_init reads the tables again, the envelope
! table, built from the earlier ones, is not used: a lookup at the
! envelope composition takes SCV's general path, the value it has when
! the table is built for another composition (and not the envelope
! table's, 1e-5 dex away). A missing He table, and an H table cut short,
! return ierr = 1, and the tables opened before are closed.
subroutine scv_reader_checks()
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0, lt = 3.73d0, lp = 5.1d0
      type(eos_result) :: table, after_read, general
      type(eos_settings_t) :: s
      type(eos_table_paths_t) :: paths
      integer :: ierr, ierr1, ierr2, ierr3, u_in, u_out, k, ios
      logical :: h_open
      character(len=200) :: line
      call configure_eos("scv_only", scv=.true.)
      call eos_at(lt, lp, x, z, table, ierr1)
      s = eos_settings
      call eos_init(s, table_paths, ierr)
      call eos_at(lt, lp, x, z, after_read, ierr2)
      call set_envelope_mixture(0.60d0, 0.03d0)
      call eos_at(lt, lp, x, z, general, ierr3)
      call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
      call check("SCV reader: after a re-read the envelope table waits for its rebuild", &
           ierr == 0 .and. ierr1 == 0 .and. ierr2 == 0 .and. ierr3 == 0 .and. &
           after_read%log10_density == general%log10_density .and. &
           after_read%log10_density /= table%log10_density)
      paths = table_paths
      paths%scv_he = "scv_missing"
      call eos_init(s, paths, ierr)
      inquire(file=scvh_path, opened=h_open)
      call check("SCV reader: a missing table returns ierr = 1 and leaves no table open", &
           ierr == 1 .and. .not. h_open)
      open(newunit=u_in, file=scvh_path, status="old", action="read")
      open(newunit=u_out, file="scv_truncated", status="replace")
      do k = 1, 300
         read(u_in, '(a)', iostat=ios) line
         if (ios /= 0) exit
         write(u_out, '(a)') trim(line)
      end do
      close(u_in)
      close(u_out)
      paths = table_paths
      paths%scv_h = "scv_truncated"
      call eos_init(s, paths, ierr)
      call check("SCV reader: a truncated table returns ierr = 1", ierr == 1)
      call configure_eos("scv_only", scv=.true.)
end subroutine scv_reader_checks

! The density searches of the 1995 and 2001 tables call the interpolator
! for the pressure alone (deriv_order = 1) with the radiation correction;
! the other slots then hold the last full call's values, and the trial
! calls leave them so (the radiation correction once rewrote cv from
! them on every trial call, growing it a hundredfold each time). The
! pressure is the full call's.
subroutine opal_trial_call_checks()
      real(dp), parameter :: x = 0.70d0, t6 = 2.0d0, rho = 0.1d0
      real(dp) :: cv_full, p_full
      integer :: k, ierr
      logical :: covers, ok
      character(len=256) :: legacy
      legacy = trim(yrec_input)//"/eos/legacy/"
      call configure_eos("opal95", opal95=.true., scv=.true., file=trim(legacy)//"opal/OPALeos.0188")
      call esac(opal95, x, t6, rho, 10, 1, ierr, covers)
      ok = ierr == 0 .and. covers
      cv_full = opal95%eos_output(i_opal_cv)
      p_full = opal95%eos_output(i_opal_p)
      do k = 1, 20
         call esac(opal95, x, t6, rho, 1, 1, ierr, covers)
         ok = ok .and. ierr == 0
      end do
      call check("OPAL 1995: trial calls keep the other slots and give the full call's pressure", &
           ok .and. opal95%eos_output(i_opal_cv) == cv_full .and. opal95%eos_output(i_opal_p) == p_full)
      call configure_eos("opal2001", opal2001=.true., scv=.true., file=trim(legacy)//"opal2001/O01EOS0188z")
      call esac01(opal_tables(iv_opal01), x, t6, rho, 9, 1, ierr, covers)
      ok = ierr == 0 .and. covers
      cv_full = opal_tables(iv_opal01)%eos_output(i_opal01_cv)
      p_full = opal_tables(iv_opal01)%eos_output(i_opal01_p)
      do k = 1, 20
         call esac01(opal_tables(iv_opal01), x, t6, rho, 1, 1, ierr, covers)
         ok = ok .and. ierr == 0
      end do
      call check("OPAL 2001: trial calls keep the other slots and give the full call's pressure", &
           ok .and. opal_tables(iv_opal01)%eos_output(i_opal01_cv) == cv_full .and. opal_tables(iv_opal01)%eos_output(i_opal01_p) == p_full)
end subroutine opal_trial_call_checks

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
! an envelope without metals (Z_env = 0): nothing to scale, so the zone's
! metals are counted as carbon, as without diffusion of Z (once a
! division by Z_env = 0)
      call set_metal_free_envelope(envelope_hydrogen_fraction)
      ok = .true.
      do k = 1, 2
         call configure_eos("yale_only")
         call eos_at(lt(k), lp(k), x, 0.01d0, plain, ierr1)
         call configure_eos("yale_diffusion_z", diffusion_z=.true.)
         call eos_at(lt(k), lp(k), x, 0.01d0, diffusion, ierr2)
         ok = ok .and. ierr1 == 0 .and. ierr2 == 0 .and. &
              diffusion%ion_mean_weight_inverse == plain%ion_mean_weight_inverse .and. &
              diffusion%log10_density == plain%log10_density
      end do
      call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
      call check("diffusion of Z: an envelope without metals takes the default mixture", ok)
end subroutine diffusion_z_checks

! an envelope of hydrogen (x) and helium alone, as eos_set_mixture takes it
subroutine set_metal_free_envelope(x)
      real(dp), intent(in) :: x
      real(dp) :: fractions(12), inverse_weight
      fractions = 0.0d0
      inverse_weight = x/atomic_weight(7) + (1.0d0 - x)/atomic_weight(12)
      fractions(7) = x/atomic_weight(7)/inverse_weight
      fractions(12) = (1.0d0 - x)/atomic_weight(12)/inverse_weight
      call eos_set_mixture(x, 0.0d0, inverse_weight, fractions)
end subroutine set_metal_free_envelope

! eos_init refuses settings that select more than one OPAL vintage (the
! builder's opal_vintage_conflict) and installs nothing.
subroutine settings_checks()
      type(eos_settings_t) :: s, before
      integer :: ierr
      before = eos_settings
      s = eos_settings
      s%opal_vintage = opal_vintage_conflict
      call eos_init(s, table_paths, ierr)
      call check("settings: two OPAL vintages selected, eos_init returns ierr = 1", &
           ierr == 1 .and. eos_settings%opal_vintage == before%opal_vintage)
      s = eos_settings
      s%use_debye_huckel_correction = .true.
      s%debye_huckel_eta_min = 1.0d0
      s%debye_huckel_eta_max = 1.0d0
      call eos_init(s, table_paths, ierr)
      call check("settings: a Debye-Hueckel range with eta_max <= eta_min, eos_init returns ierr = 1", &
           ierr == 1 .and. (eos_settings%use_debye_huckel_correction .eqv. before%use_debye_huckel_correction))
! the defaults are the controls' (defaults/controls_registry.tsv)
      s = eos_settings_t()
      call check("settings: the record's defaults are the controls' (cutoff 6, eta -1 to 1)", &
           s%saha_log10t_cutoff == 6.0d0 .and. s%debye_huckel_eta_min == -1.0d0 .and. &
           s%debye_huckel_eta_max == 1.0d0 .and. s%opal_vintage == opal_vintage_none .and. &
           .not. (s%use_scv_eos .or. s%use_mhd_eos .or. s%use_debye_huckel_correction))
end subroutine settings_checks

! Below the OPAL 2006 table (T = 1500 K) gamma1_at_density takes the
! Yale/SCV value (to 1e-6), with OPAL selected as without it.
subroutine gamma1_checks()
      type(eos_settings_t) :: s
      real(dp) :: g1_yale, ga_yale, g1, ga
      integer :: ierr1, ierr2, iterations
      call configure_eos("production", opal2006=.true., scv=.true., file=opal06_path)
      s = eos_settings
      s%opal_vintage = opal_vintage_none
      call set_eos_settings(s)
      iterations = 0
      call gamma1_at_density(0.70d0, 0.016492d0, 0.0015d0, 1.0d-9, 1.0d3, g1_yale, ga_yale, &
           iterations, ierr=ierr1)
      s%opal_vintage = opal_vintage_2006
      call set_eos_settings(s)
      iterations = 0
      call gamma1_at_density(0.70d0, 0.016492d0, 0.0015d0, 1.0d-9, 1.0d3, g1, ga, &
           iterations, ierr=ierr2)
      call check("gamma1_at_density: below the OPAL table, the Yale/SCV value", &
           ierr1 == 0 .and. ierr2 == 0 .and. abs(g1/g1_yale - 1.0d0) < 1.0d-6 .and. &
           abs(ga/ga_yale - 1.0d0) < 1.0d-6)
end subroutine gamma1_checks

! Where SCV gives the density (below the Saha ramp, inside its table),
! R/mu is the gas constant that density implies, beta P/(rho T), and
! 1/mu_e = (R/mu)/R - 1/mu_ion, both on the envelope table's path and on
! the general one (the Saha solution's own R/mu would differ by percents
! where SCV has molecules), and R/mu is continuous at the ramp's foot.
subroutine scv_mean_weight_checks()
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0
      real(dp), parameter :: lt(3) = [3.6d0, 4.0d0, 5.2d0], lp(3) = [5.0d0, 7.0d0, 12.0d0]
      real(dp), parameter :: h = 1.0d-7
      type(eos_result) :: r, below, above
      integer :: k, c, ierr, ierr1, ierr2
      real(dp) :: x_c
      logical :: ok
      call configure_eos("scv_only", scv=.true.)
      ok = .true.
      do c = 1, 2
         x_c = x
         if (c == 2) x_c = x + 2.0d-5
         do k = 1, size(lt)
            call eos_at(lt(k), lp(k), x_c, z, r, ierr)
            ok = ok .and. ierr == 0 .and. &
                 abs(r%specific_gas_constant/(r%beta*r%pressure/(r%density*r%temperature)) - 1.0d0) < 1.0d-12 &
                 .and. abs(r%electron_mean_weight_inverse - (r%specific_gas_constant/gas_constant - &
                 r%ion_mean_weight_inverse)) < 1.0d-12
         end do
      end do
      call check("SCV: R/mu = beta P/(rho T) and 1/mu_e from it, on both SCV paths", ok)
      call eos_at(5.5d0 - h, 13.0d0, x, z, below, ierr1)
      call eos_at(5.5d0 + h, 13.0d0, x, z, above, ierr2)
      call check("SCV: R/mu and 1/mu_e continuous at the Saha ramp's foot", ierr1 == 0 .and. ierr2 == 0 .and. &
           abs(above%specific_gas_constant/below%specific_gas_constant - 1.0d0) < 1.0d-6 .and. &
           abs(above%electron_mean_weight_inverse - below%electron_mean_weight_inverse) < 1.0d-6)
end subroutine scv_mean_weight_checks

! cp and grad_ad of an ideal gas plus radiation with gas-pressure
! fraction b and R/mu = r_mu: with delta = (4 - 3b)/b,
! cp = (R/mu)(3/2 + 12(1 - b)/b + delta^2) and grad_ad = (R/mu) delta/(b cp)
! (Kippenhahn & Weigert sec. 13.2)
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
      rho_ideal = r%beta*r%pressure/(gas_constant*(inverse_atomic_weight + &
           r%electron_mean_weight_inverse)*r%temperature)
      call ideal_radiation_gas(r%beta, r%specific_gas_constant, cp, grad)
      call check("ideal gas, neutral: rho = mu P_gas/(R T)", ierr == 0 .and. &
           abs(r%density/rho_ideal - 1.0d0) < 1.0d-12 .and. &
           r%electron_mean_weight_inverse < 1.0d-4*inverse_atomic_weight)
      call check("ideal gas, neutral: grad_ad -> 0.4, cp -> 5R/(2 mu) (+ radiation)", ierr == 0 .and. &
           abs(r%adiabatic_gradient - 0.4d0) < 1.0d-3 .and. &
           abs(r%adiabatic_gradient/grad - 1.0d0) < 2.0d-3 .and. &
           abs(r%specific_heat_cp/cp - 1.0d0) < 2.0d-3)
! fully ionized: rho = mu beta P/(R T) (1 - 2^(-5/2) e^eta P_e/P_gas),
! the weak-degeneracy expansion of the electron pressure, with 1/mu_e
! from the composition (X/1.008 + (1 - X)/2.002)
      mu_e_inv = x/1.008d0 + (1.0d0 - x)/2.002d0
      mu_inv = inverse_atomic_weight + mu_e_inv
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

! Every field a blend sets, on both sides of a boundary in log T (fixed
! log P), 2e-7 in log T apart: log rho, rho, R/mu and 1/mu_e
! (mean_weights false leaves these two out), the ionization fractions,
! the four first derivatives and the six second derivatives. The
! tolerances: 1e-6 (relative for rho, R/mu and cp), 1e-4 for the
! second derivatives, which carry the fully ionized iteration's
! tolerance; a field left out of a blend jumps by the gap between the
! two solutions there, 1e-3 to 1e-1. The fields that jump are named.
subroutine continuity_in_t(label, lt0, lp, mean_weights)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: lt0, lp
      logical, intent(in) :: mean_weights
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0, h = 1.0d-7
      type(eos_result) :: below, above
      integer :: ierr1, ierr2
      call eos_at(lt0 - h, lp, x, z, below, ierr1)
      call eos_at(lt0 + h, lp, x, z, above, ierr2)
      call check(label, ierr1 == 0 .and. ierr2 == 0 .and. &
           results_agree(below, above, 1.0d-6, 1.0d-4, mean_weights, .true.))
end subroutine continuity_in_t

! Whether two results agree field by field (a to tol, the second
! derivatives to tol2 when second is set; rho, R/mu and cp relative),
! naming in the output the fields that do not.
logical function results_agree(a, b, tol, tol2, mean_weights, second)
      type(eos_result), intent(in) :: a, b
      real(dp), intent(in) :: tol, tol2
      logical, intent(in) :: mean_weights, second
      character(len=200) :: off
      off = ""
      call note_gap(off, "log10_density", abs(a%log10_density - b%log10_density), tol)
      call note_gap(off, "density", abs(a%density/b%density - 1.0d0), tol)
      call note_gap(off, "beta", abs(a%beta - b%beta), tol)
      if (mean_weights) then
         call note_gap(off, "R/mu", abs(a%specific_gas_constant/b%specific_gas_constant - 1.0d0), tol)
         call note_gap(off, "1/mu_e", abs(a%electron_mean_weight_inverse - b%electron_mean_weight_inverse), tol)
      end if
      call note_gap(off, "ion_fraction", maxval(abs(a%ion_fraction - b%ion_fraction)), tol)
      call note_gap(off, "dlnrho_dlnt", abs(a%dlnrho_dlnt - b%dlnrho_dlnt), tol)
      call note_gap(off, "dlnrho_dlnp", abs(a%dlnrho_dlnp - b%dlnrho_dlnp), tol)
      call note_gap(off, "cp", abs(a%specific_heat_cp/b%specific_heat_cp - 1.0d0), tol)
      call note_gap(off, "grad_ad", abs(a%adiabatic_gradient - b%adiabatic_gradient), tol)
      if (second) then
         call note_gap(off, "dlnrho_dlnt_dt", abs(a%dlnrho_dlnt_dt - b%dlnrho_dlnt_dt), tol2)
         call note_gap(off, "dlnrho_dlnt_dp", abs(a%dlnrho_dlnt_dp - b%dlnrho_dlnt_dp), tol2)
         call note_gap(off, "grad_ad_dt", abs(a%adiabatic_gradient_dt - b%adiabatic_gradient_dt), tol2)
         call note_gap(off, "grad_ad_dp", abs(a%adiabatic_gradient_dp - b%adiabatic_gradient_dp), tol2)
         call note_gap(off, "cp_dt", abs(a%specific_heat_cp_dt - b%specific_heat_cp_dt), tol2)
         call note_gap(off, "cp_dp", abs(a%specific_heat_cp_dp - b%specific_heat_cp_dp), tol2)
      end if
      results_agree = len_trim(off) == 0
      if (.not. results_agree) write(*,'(a,a)') "     jump in:", trim(off)
end function results_agree

! add name to the list off when gap is not below limit
subroutine note_gap(off, name, gap, limit)
      character(len=*), intent(inout) :: off
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: gap, limit
      if (.not. gap < limit) off = trim(off)//" "//name
end subroutine note_gap

! Continuity at both ends of the Saha ramp (its foot, cutoff - 0.5, and
! the cutoff), with the Yale solvers alone and with SCV under the ramp's
! foot; and at both ends of the OPAL 2006 table's edge ramp. At the
! cutoff R/mu and 1/mu_e are left out: below it they are beta P/(rho T),
! above it the composition's, 7e-4 apart (an open question to the
! authors, docs/bug_sweep.md).
subroutine continuity_checks()
      call configure_eos("yale_only")
      call continuity_in_t("continuity: the blended quantities at the Saha ramp's foot (cutoff - 0.5)", &
           5.5d0, 13.0d0, .true.)
      call continuity_in_t("continuity: the blended quantities at the Saha cutoff", 6.0d0, 13.0d0, .false.)
      call continuity_in_t("continuity: the blended quantities at the Saha cutoff, denser", &
           6.0d0, 15.0d0, .false.)
      call configure_eos("scv_only", scv=.true.)
      call continuity_in_t("continuity: the blended quantities at the Saha ramp's foot, SCV under it", &
           5.5d0, 13.0d0, .true.)
! the OPAL 2006 table's edge at log T = 5, in log P: the blend is
! OPAL's own result below log P = 12.87, the Yale/SCV result above
! 13.01, and ramps between
      call configure_eos("production", opal2006=.true., scv=.true., file=opal06_path)
      call check("continuity: where the OPAL 2006 edge ramp reaches OPAL", &
           ramp_end_is_continuous(5.0d0, 12.80d0, 12.94d0, opal2006_alone))
      call check("continuity: at the OPAL 2006 table's edge (the Yale/SCV result beyond)", &
           ramp_end_is_continuous(5.0d0, 13.05d0, 12.94d0, without_opal))
end subroutine continuity_checks

! A table source's edge ramp at log T = lt, in log P, between lp_equal
! (where the blend is the reference: the source alone at the ramp's inner
! end, the blend without the source at its outer end) and lp_differ
! (inside the ramp): the end is bisected to 1e-11 in log P and every
! field the table blend sets (log rho, rho, beta, the mean weights, the
! first derivatives) compared on its two sides, to 1e-7 (the gap between
! the two results there is 1e-3 to 1e-2; beta is OPAL's and Yale's alike
! there, so its omission would not show). Any source passes its own two
! references.
logical function ramp_end_is_continuous(lt, lp_equal, lp_differ, reference)
      real(dp), intent(in) :: lt, lp_equal, lp_differ
      procedure(reference_eos) :: reference
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0
      real(dp) :: a, b, m
      type(eos_result) :: below, above
      integer :: k, ierr1, ierr2
      logical :: ok
      a = lp_equal
      b = lp_differ
      ok = blend_equals_reference(lt, a, reference) .and. .not. blend_equals_reference(lt, b, reference)
      do k = 1, 40
         m = 0.5d0*(a + b)
         if (blend_equals_reference(lt, m, reference)) then
            a = m
         else
            b = m
         end if
      end do
      call eos_at(lt, a, x, z, below, ierr1)
      call eos_at(lt, b, x, z, above, ierr2)
      ramp_end_is_continuous = ok .and. ierr1 == 0 .and. ierr2 == 0 .and. &
           results_agree(below, above, 1.0d-7, 0.0d0, .true., .false.)
end function ramp_end_is_continuous

! Whether the blend at (lt, lp) is the reference's result: the same log
! rho, bit for bit (adopting a table's result copies it; a blend with a
! weight below 1 does not), so that the bisection finds where the blend
! starts, not where it has moved rho by a rounding error.
logical function blend_equals_reference(lt, lp, reference)
      real(dp), intent(in) :: lt, lp
      procedure(reference_eos) :: reference
      real(dp), parameter :: x = 0.70d0, z = 0.016492d0
      type(eos_result) :: blend, ref
      integer :: ierr
      call eos_at(lt, lp, x, z, blend, ierr)
      call reference(lt, lp, x, z, ref)
      blend_equals_reference = blend%log10_density == ref%log10_density
end function blend_equals_reference

! The references of the OPAL 2006 edge ramp: OPAL's own result, and the
! blend with OPAL switched off.
subroutine opal2006_alone(lt, lp, x, z, r)
      real(dp), intent(in) :: lt, lp, x, z
      type(eos_result), intent(out) :: r
      integer :: ierr
      logical :: covers
      r = eos_result()
      call compute_opal_eos(iv_opal06, lt, lp, x, z, r, ierr, covers)
end subroutine opal2006_alone

subroutine without_opal(lt, lp, x, z, r)
      real(dp), intent(in) :: lt, lp, x, z
      type(eos_result), intent(out) :: r
      type(eos_settings_t) :: s, before
      integer :: ierr
      before = eos_settings
      s = eos_settings
      s%opal_vintage = opal_vintage_none
      call set_eos_settings(s)
      call eos_at(lt, lp, x, z, r, ierr)
      call set_eos_settings(before)
end subroutine without_opal

! in_atmosphere = .true. at Saha points, in the ramp and below it, with
! and without SCV: the same density, beta, ionization fractions, first
! derivatives, cp and grad_ad as a full evaluation (to 1e-10, and the
! derivatives to 1e-8)
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

! When the envelope composition changes, eos_set_mixture builds the SCV
! envelope table for it: a lookup at the new composition is the table's
! (not the general path's value) and agrees with SCV's general path to
! 5e-4 in log rho (3e-4 at log T = 4.13; a table left at the old
! composition is off by 0.05);
! and a zone at the composition of the table before takes the general
! path, the value it had before that table was built. Off the table's
! nodes in log T, so that the general path's cell does not depend on
! where its search starts.
subroutine envelope_table_checks()
      real(dp), parameter :: x_new = 0.60d0, z_new = 0.03d0
      real(dp), parameter :: lt(3) = [3.73d0, 3.93d0, 4.13d0], lp(3) = [5.0d0, 6.5d0, 8.0d0]
      type(eos_result) :: table, general, general_before(3)
      integer :: k, ierr1, ierr2
      logical :: ok
      call configure_eos("scv_only", scv=.true.)
      ok = .true.
      do k = 1, size(lt)
         call eos_at(lt(k), lp(k), x_new, z_new, general_before(k), ierr1)
         ok = ok .and. ierr1 == 0
      end do
      call set_envelope_mixture(x_new, z_new)
      do k = 1, size(lt)
         call eos_at(lt(k), lp(k), x_new, z_new, table, ierr1)
         call eos_at(lt(k), lp(k), x_new + 2.0d-5, z_new, general, ierr2)
         ok = ok .and. ierr1 == 0 .and. ierr2 == 0 .and. &
              abs(table%log10_density - general%log10_density) < 5.0d-4 .and. &
              table%log10_density /= general_before(k)%log10_density
      end do
      call check("envelope tables: eos_set_mixture builds the SCV table for a new composition", ok)
      call set_envelope_mixture(0.65d0, 0.025d0)
      ok = .true.
      do k = 1, size(lt)
         call eos_at(lt(k), lp(k), x_new, z_new, general, ierr1)
         ok = ok .and. ierr1 == 0 .and. &
              general%log10_density == general_before(k)%log10_density
      end do
      call check("envelope tables: a table built for another composition is not used", ok)
      call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
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
! the heavier metals follow the zone's Z: a zone at twice the envelope's
! Z (its C, N and O doubled too) gets the density rise of an envelope at
! that Z, to 3e-3 (the gas under the correction differs a little: the
! envelope's 1/mu_ion; 3e-4 here); with the heavier metals left at the
! envelope's Z the rise is 3% smaller
      ok = .true.
      do k = 1, 2
         if (k == 1) call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
         if (k == 2) call set_envelope_mixture(envelope_hydrogen_fraction, 2.0d0*envelope_metal_fraction)
         call configure_eos("yale_only")
         call eos_at(lt, lp, x, 2.0d0*z, none, ierr1)
         call configure_eos("debye_huckel", debye_huckel=.true.)
         call eos_at(lt, lp, x, 2.0d0*z, with_dh, ierr2, comp=doubled_z(comp))
         ok = ok .and. ierr1 == 0 .and. ierr2 == 0
         raise(k) = with_dh%density/none%density - 1.0d0
      end do
      call set_envelope_mixture(envelope_hydrogen_fraction, envelope_metal_fraction)
      call check("Debye-Hueckel: the heavier metals scale with the zone's Z", &
           ok .and. abs(raise(1)/raise(2) - 1.0d0) < 3.0d-3)
      call record_metals_check(x, z, comp)
end subroutine debye_huckel_checks

! The record debye_huckel_composition_of builds from a zone's
! composition: its metals add up to the zone's Z, those beyond O in the
! envelope mixture's proportions. A zone whose Z exceeds another's by
! carbon alone, as where a core burns helium, has the other's Ne..Fe
! (scaling them with the zone's Z would multiply them by 31 here); a
! zone with every metal doubled has them doubled.
subroutine record_metals_check(x, z, comp)
      real(dp), intent(in) :: x, z, comp(:)
      integer, parameter :: beyond_o(7) = [idh_ne, idh_na, idh_mg, idh_al, idh_si, &
           idh_ar, idh_fe]
      real(dp), parameter :: burnt_carbon = 0.5d0
      real(dp) :: carbon_rich(size(comp))
      type(debye_huckel_composition) :: zone, burnt, doubled
      logical :: ok

      carbon_rich = comp
      carbon_rich(2) = comp(2) - burnt_carbon
      carbon_rich(3) = comp(3) + burnt_carbon
      carbon_rich(5) = comp(5) + burnt_carbon
      zone = debye_huckel_composition_of(x, z, comp)
      burnt = debye_huckel_composition_of(x, z + burnt_carbon, carbon_rich)
      doubled = debye_huckel_composition_of(x, 2.0d0*z, doubled_z(comp))
      ok = abs(sum(zone%metals) - z) < 1.0d-12*z .and. &
           abs(sum(burnt%metals) - (z + burnt_carbon)) < 1.0d-12 .and. &
           all(abs(burnt%metals(beyond_o) - zone%metals(beyond_o)) <= &
           1.0d-10*zone%metals(beyond_o)) .and. &
           all(abs(doubled%metals(beyond_o) - 2.0d0*zone%metals(beyond_o)) <= &
           1.0d-10*zone%metals(beyond_o)) .and. &
           burnt%metals(idh_c) == carbon_rich(5) .and. zone%metals(idh_n) == comp(7) .and. &
           zone%metals(idh_o) == comp(9) .and. all(zone%metals(beyond_o) > 0.0d0)
      call check("Debye-Hueckel: a zone's metals add up to its Z; carbon a core makes leaves Ne..Fe alone", ok)
end subroutine record_metals_check

! a composition (star%xa layout) with its metals doubled, at the expense
! of He-4
pure function doubled_z(c) result(d)
      real(dp), intent(in) :: c(:)
      real(dp) :: d(size(c))
      d = c
      d(2) = c(2) - c(3)
      d(3:15) = 2.0d0*c(3:15)
end function doubled_z

subroutine check(label, ok)
      character(len=*), intent(in) :: label
      logical, intent(in) :: ok
      write(*,'(a,1x,a)') merge("ok  ", "FAIL", ok), label
      if (.not. ok) failures = failures + 1
end subroutine check

end program test_eos
