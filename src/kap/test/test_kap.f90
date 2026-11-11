!----------------------------------------------------------------------
! test_kap
!----------------------------------------------------------------------
! kap's test program (pytest test_kap.py from the repository root, after
! make tests in src/). Its output has two parts:
!   - the pinned part, byte-compared with expected_test_kap.out: the OPAL95
!     interior on a (log T, log rho) grid, booted through the star layer's
!     set_up_constants_and_tables as a run would; two error paths; and every
!     opacity family loaded through kap_init (check_families);
!   - the ok/FAIL checks (run_checks): exactness at table nodes, the
!     derivatives against differences, continuity at the ramp's ends and the
!     table's edges, the blends' and the conductive mixture's formulas, and
!     the error contracts. Each check is confirmed to fail under a mutation in
!     test_kap.mutations.toml (python3 src/tools/mutate.py kap). A FAIL makes
!     the program stop with status 1.
! The lines that print a header kept from older versions ("kap_eval over",
! "via optional ierr") are pinned and stay as they are.
program test_kap
      use run_setup, only: set_up_constants_and_tables
      use run_messages, only: run_log_unit
      use star_info_lib, only: star_info
      use kap_lib
      use opacity_tables
      implicit none
! the test's star (star_info_lib declares no instance); static, as it
! is too large for the stack
      type(star_info), save :: star
      type(kap_settings_t) :: settings
      integer :: gerr

      character(len=256) :: yrec_input
      character(len=256) :: fermi_path, opal95_path, dummy_path
      character(len=256) :: dummy_paths7(7)

      integer :: setups_ierr, ipt, kap_ierr
      double precision :: logt, logd, x_frac, z_frac
      double precision :: o, ol, qod, qot, fxion(3)

! the checks' failure count and table directories (run_checks)
      integer :: check_failures = 0
      character(len=256) :: atom_dir, mol_dir

      integer, parameter :: npts = 7
      double precision :: grid_logt(npts), grid_logd(npts)
      data grid_logt /4.2d0, 4.6d0, 5.2d0, 6.0d0, 6.6d0, 7.0d0, 7.4d0/
      data grid_logd /-8.5d0, -7.5d0, -6.0d0, -3.0d0, -1.0d0, 0.2d0, 1.3d0/

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_kap: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      fermi_path  = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      opal95_path = trim(yrec_input)//"/opacity/atomic/legacy/opal95/GS98.OP17"
      dummy_path = ""
      dummy_paths7 = ""

! unit numbers, per core/read_input.f90
      run_log_unit = 20
      open(run_log_unit, file="test_kap.short", status="replace")

! kap configuration: OPAL95 atomic tables only, per the reference
! solar namelists (LOPAL95=T, ZOPAL951=0.016232, TMOLMIN/TMOLMAX
! 4.0/4.1); molecular and conductive paths off
      star%control%use_opal95_tables = .true.
      star%control%use_opal92_tables = .false.
      star%control%use_laol89_tables = .false.
      star%control%use_alex06_tables = .false.
      star%control%use_alex95_tables = .false.
      star%control%use_kurucz90_tables = .false.
      star%control%use_pure_z_table = .false.
      star%physics%use_conductive_opacity = .false.
      star%control%opal95_single_table_z = 0.016232d0
      star%physics%molecular_opacity_logt_min = 4.0d0
      star%physics%molecular_opacity_logt_max = 4.1d0

! eos side gated off except the always-loaded Fermi table (setups
! calls eos_init unconditionally)
      star%physics%use_mhd_eos = .false.
      star%physics%use_scv_eos = .false.
      star%physics%opal_eos_vintage = 0
      star%atm_choice_in_use = 0

      star%envelope_hydrogen_fraction = 0.70d0
      star%envelope_metal_fraction = 0.016232d0

! constants + table loads (real setups; kap_init inside it reads the
! OPAL95 table and builds the surface-X slice)
      star%control%alex06_table_path = dummy_path
      star%control%allard_table_path = dummy_path
      star%control%atm_table_path = dummy_path
      star%control%fermi_table_path = fermi_path
      star%control%kurucz_table_path = dummy_path
      star%control%kurucz_table2_path = dummy_path
      star%control%laol_table_path = dummy_path
      star%control%laol_table2_path = dummy_path
      star%control%opal95_table_path = opal95_path
      star%control%opal92_table_path = dummy_path
      star%control%zams_a_table_path = dummy_path
      star%control%zams_b_table_path = dummy_path
      star%control%zams_c_table_path = dummy_path
      star%control%centre1_table_path = dummy_path
      star%control%centre2_table_path = dummy_path
      star%control%centre3_table_path = dummy_path
      star%control%centre4_table_path = dummy_path
      star%control%centre5_table_path = dummy_path
      star%control%opal92_table2_path = dummy_path
      star%control%pure_z_table_path = dummy_path
      star%control%conductive_table_path = dummy_path
      star%control%scv_h_table_path = dummy_path
      star%control%scv_he_table_path = dummy_path
      star%control%scv_z_table_path = dummy_path
      star%control%alex95_table_paths = dummy_paths7
      call set_up_constants_and_tables(star, setups_ierr)
      if (setups_ierr /= 0) then
         write(*,'(a)') "test_kap: FAIL (setups error)"
         stop 1
      end if

      x_frac = 0.70d0
      z_frac = 0.016232d0

      write(*,'(a)') "# test_kap: kap_eval over (logT, logRho) grid, " // &
           "X=0.70 Z=0.016232, OPAL95/GS98.OP17, no molecular, " // &
           "no conductive"
      do ipt = 1, npts
         logt = grid_logt(ipt)
         logd = grid_logd(ipt)
         fxion = 0.0d0
         call evaluate(logd, logt, x_frac, z_frac, gerr)
         if (gerr /= 0) stop 1
         write(*,'(a,i2,4(1pe24.15))') "kap ", ipt, o, ol, qod, qot
      end do

! Error paths (2026, docs/ROADMAP.md stage 3): out-of-range points and
! misconfiguration return ierr /= 0 instead of stopping -- the first
! time these paths are testable at all. The diagnostic each failure
! writes goes to run_log_unit / stdout at the point of failure, as
! always. (The "optional ierr" wording in the output line below is
! historical and kept because the expected output is byte-compared.)
      write(*,'(a)') "# test_kap: error paths via optional ierr"
      logt = 3.5d0
      logd = -8.0d0
      fxion = 0.0d0
      call evaluate(logd, logt, x_frac, z_frac, kap_ierr)
      write(*,'(a,i4)') "err out-of-table   ierr = ", kap_ierr
      settings = kap_settings
      settings%use_opal95_tables = .false.
      call set_kap_settings(settings)
      logt = 6.0d0
      logd = -3.0d0
      fxion = 0.0d0
      call evaluate(logd, logt, x_frac, z_frac, kap_ierr)
      write(*,'(a,i4)') "err no-table-chosen ierr = ", kap_ierr
      settings%use_opal95_tables = .true.
      call set_kap_settings(settings)

      call check_families()
      call run_checks()

      close(run_log_unit)
      write(*,'(a)') "test_kap: done"

contains

! Every opacity family, each loaded through kap_init with its own
! configuration and evaluated on points that cross its regimes:
!   production   Alexander 2006 molecular + OPAL95-format atomic + pure Z
!                + Potekhin conductive, the molecular/atomic ramp, the
!                single-Z surface table (Z = 0.016232) and the 4D
!                interpolation away from it (Z = 0.020), a helium core
!                on the pure-Z blend
!   opal92       two Z tables (0.0188, 0.0200), blended at Z = 0.019, and
!                the top temperature node
!   laol89       two Z tables (0.0188, 0.0190), blended at Z = 0.0189
!   kurucz90     Kurucz molecular tables at [M/H] = 0 and -0.5, blended,
!                over OPAL95
!   alex95       the seven Alexander 1994/95 X tables over OPAL95
! Each point prints kappa, log kappa, dln kappa/dln rho, dln kappa/dln T.
      subroutine check_families()
         character(len=256) :: in, mol, atom, paths7(7)
         integer :: k
         in = yrec_input
         atom = trim(in)//"/opacity/atomic/"
         mol = trim(in)//"/opacity/molecular/"

         call configure("production", opal95=.true., alex06=.true., pure_z=.true., &
              conductive=.true., z1=0.016232d0, z2=0.0d0, two_z=.false., &
              opal95_path=trim(atom)//"legacy/opal95/GS98.OP17", &
              alex06_path=trim(mol)//"alexmol06gs98.tab", &
              pure_z_path=trim(atom)//"lanl/PURECO.DBGLAOL")
         call evaluate_points(0.70d0, 0.016232d0)
         call evaluate_points(0.70d0, 0.020d0)
         call evaluate_point(0.0d0, 0.60d0, 7.8d0, 4.0d0)

         call configure("opal92", opal92=.true., z1=0.0188d0, z2=0.0200d0, two_z=.true., &
              opal92_path=trim(atom)//"legacy/opal92/Z0188.OPALN", &
              opal92_path2=trim(atom)//"legacy/opal92/Z0200.OPALN")
         call evaluate_points(0.70d0, 0.019d0, interior_only=.true.)
! the top temperature node (T6 = 100), where only three temperature rows
! enter the spline in log T
         call evaluate_point(0.70d0, 0.019d0, 8.0d0, 3.0d0)

         call configure("laol89", laol89=.true., z1=0.0188d0, z2=0.0190d0, two_z=.true., &
              laol_path=trim(atom)//"legacy/opal92/Z0188AG89.DBGLAOL", &
              laol_path2=trim(atom)//"legacy/opal92/Z0190AG91.DBGLAOL")
         call evaluate_points(0.70d0, 0.0189d0, interior_only=.true.)

         call configure("kurucz90", opal95=.true., kurucz90=.true., z1=0.0188d0, &
              z2=0.00595d0, two_z=.true., &
              opal95_path=trim(atom)//"legacy/opal95/GS98.OP17", &
              kurucz_path=trim(mol)//"legacy/kurmol/kapp00.ross", &
              kurucz_path2=trim(mol)//"legacy/kurmol/kapm05.ross")
         call evaluate_points(0.70d0, 0.016232d0)

         paths7 = [character(len=256) :: trim(mol)//"legacy/alex/OPACALEXANDER.X00", &
              trim(mol)//"legacy/alex/OPACALEXANDER.X01", trim(mol)//"legacy/alex/OPACALEXANDER.X02", &
              trim(mol)//"legacy/alex/OPACALEXANDER.X035", trim(mol)//"legacy/alex/OPACALEXANDER.X05", &
              trim(mol)//"legacy/alex/OPACALEXANDER.X07", trim(mol)//"legacy/alex/OPACALEXANDER.X08"]
         call configure("alex95", opal95=.true., alex95=.true., z1=0.016232d0, z2=0.0d0, &
              two_z=.false., opal95_path=trim(atom)//"legacy/opal95/GS98.OP17", alex95_paths=paths7)
         call evaluate_points(0.70d0, 0.016232d0)
      end subroutine check_families

! Load kap with one family's configuration (the settings a deck's
! controls would give, then kap_init with the table paths).
      subroutine configure(label, opal95, opal92, laol89, alex06, alex95, kurucz90, &
           pure_z, conductive, z1, z2, two_z, opal95_path, opal92_path, &
           opal92_path2, laol_path, laol_path2, alex06_path, alex95_paths, &
           kurucz_path, kurucz_path2, pure_z_path, quiet)
         character(len=*), intent(in) :: label
! quiet: no header line (the checks' configurations)
         logical, intent(in), optional :: quiet
         logical, intent(in), optional :: opal95, opal92, laol89, alex06, alex95, &
              kurucz90, pure_z, conductive
         double precision, intent(in) :: z1, z2
         logical, intent(in) :: two_z
         character(len=*), intent(in), optional :: opal95_path, opal92_path, &
              opal92_path2, laol_path, laol_path2, alex06_path, kurucz_path, &
              kurucz_path2, pure_z_path
         character(len=256), intent(in), optional :: alex95_paths(7)
         type(kap_table_paths_t) :: p
         type(kap_settings_t) :: s
         integer :: ierr

         if (.not. present(quiet)) write(*,'(a,a)') "# test_kap: family ", label
         s%use_opal95_tables = present(opal95)
         s%use_opal92_tables = present(opal92)
         s%use_laol89_tables = present(laol89)
         s%use_alex06_tables = present(alex06)
         s%use_alex95_tables = present(alex95)
         s%use_kurucz90_tables = present(kurucz90)
         s%use_pure_z_table = present(pure_z)
         s%use_conductive_opacity = present(conductive)
         s%use_two_z_tables = two_z
         s%opal95_single_table_z = z1
         s%opal_table_z1 = z1
         s%opal_table_z2 = z2
         s%laol_table_z1 = z1
         s%laol_table_z2 = z2
         s%kurucz_table_z1 = z1
         s%kurucz_table_z2 = z2
         s%alex_table_z1 = z1
         s%metal_fraction_match_tolerance = 1.0d-8
         s%molecular_opacity_logt_min = 4.0d0
         s%molecular_opacity_logt_max = 4.1d0
         s%laol_extrapolation_tolerance = 10.0d0
         if (present(alex06_path)) p%alex06 = alex06_path
         if (present(kurucz_path)) p%kurucz = kurucz_path
         if (present(kurucz_path2)) p%kurucz_z2 = kurucz_path2
         if (present(laol_path)) p%laol = laol_path
         if (present(laol_path2)) p%laol_z2 = laol_path2
         if (present(opal95_path)) p%opal95 = opal95_path
         if (present(opal92_path)) p%opal92 = opal92_path
         if (present(opal92_path2)) p%opal92_z2 = opal92_path2
         if (present(pure_z_path)) p%pure_z = pure_z_path
         if (present(alex95_paths)) p%alex94 = alex95_paths
         p%conductive = trim(yrec_input)//"/opacity/conductive/potekhin/condall06.d"
         call kap_init(s, p, 0.70d0, z1, ierr)
         if (ierr /= 0) write(*,'(a,i0)') "kap_init ierr = ", ierr
      end subroutine configure

! kap_eval on points from the molecular regime through the ramp to the
! deep interior (log R = log rho - 3 log T6 = -3, the middle of the
! tables), or the interior points only.
      subroutine evaluate_points(x, z, interior_only)
         double precision, intent(in) :: x, z
         logical, intent(in), optional :: interior_only
         double precision, parameter :: logt_all(9) = [3.40d0, 3.70d0, 4.00d0, &
              4.05d0, 4.10d0, 4.50d0, 5.50d0, 6.50d0, 7.20d0]
         integer :: k, first
         first = 1
         if (present(interior_only)) first = 6
         do k = first, size(logt_all)
            call evaluate_point(x, z, logt_all(k), -3.0d0 + 3.0d0*logt_all(k) - 18.0d0)
         end do
      end subroutine evaluate_points

      subroutine evaluate_point(x, z, logt, logd)
         double precision, intent(in) :: x, z, logt, logd
         integer :: ierr
         fxion = 0.0d0
         call evaluate(logd, logt, x, z, ierr)
         if (ierr /= 0) then
            write(*,'(a,3f8.4,a,i0)') "kap X Z logT", x, z, logt, " ierr ", ierr
         else
            write(*,'(a,3f9.5,4(1pe24.15))') "kap", x, z, logt, o, ol, qod, qot
         end if
      end subroutine evaluate_point
! kap_get at one point, its results unpacked into o, ol, qod, qot
      subroutine evaluate(logd, logt, x, z, ierr)
         double precision, intent(in) :: logd, logt, x, z
         integer, intent(out) :: ierr
         double precision :: res(num_kap_results)
         res = 0.0d0
         call kap_get(logd, logt, x, z, res, fxion, ierr)
         o = res(i_kap)
         ol = res(i_log10_kap)
         qod = res(i_dlnkap_dlnrho)
         qot = res(i_dlnkap_dlnt)
      end subroutine evaluate
! The ok/FAIL checks of kap's behaviour (each confirmed to fail under a
! mutation in kap/test/test_kap.mutations.toml): exactness at table
! nodes on each of OPAL95's three paths and Alexander 2006's; the log
! derivatives against centred differences of log kappa, family by
! family and inside the molecular ramp; continuity across the ramp's
! ends and OPAL95's seams; the Z blend at its ends; the harmonic sum and
! Hubbard & Lampe's fit; conduction with ionized fractions and below
! the conductive table's lowest density; and the error contracts
! (refused settings, a missing or short table, a refused composition, a
! second kap_init). Each line is the check, its error and ok or FAIL; a
! FAIL makes the program stop with status 1 after the last check.
      subroutine run_checks()
         use opacity_source, only: kap_result
         use opacity_blending, only: blend_in_z, add_conduction, hubbard_lampe_opacity
         character(len=256) :: cond_path
         double precision :: r(num_kap_results), r2(num_kap_results), r_off(num_kap_results)
         double precision :: logd, value, value2, err
         type(kap_result) :: radiative, conductive, total, hl_plus, hl_minus
         integer :: ierr, ierr2, u, k
         logical :: exists
         character(len=200) :: line

         check_failures = 0
         atom_dir = trim(yrec_input)//"/opacity/atomic/"
         mol_dir = trim(yrec_input)//"/opacity/molecular/"
         cond_path = trim(yrec_input)//"/opacity/conductive/potekhin/condall06.d"
         write(*,'(a)') "# test_kap: checks (error, then ok or FAIL)"

! --- exactness at the table nodes: OPAL95 loaded at the node Z = 0.02
! and the node X = 0.70, so the surface table is a table of the file;
! at a node in every dimension each path returns the file's entry
! (TABLE #73 is X = 0.70, Z = 0.02; #60 X = 0.50, Z = 0.02; #74 X = 0.70,
! Z = 0.03)
         call load_nodes()
         call opal95_node_check("OPAL95 2D path (surface X, Z): exact at a table node", &
              0.70d0, 0.02d0, 6.0d0, -3.0d0, 73)
         call opal95_node_check("OPAL95 3D path (X off the surface): exact at a table node", &
              0.50d0, 0.02d0, 6.0d0, -3.0d0, 60)
         call opal95_node_check("OPAL95 4D path (Z off the fixed-Z table): exact at a node", &
              0.70d0, 0.03d0, 5.5d0, -2.0d0, 74)
         call node_check("Alexander 2006: exact at a table node (log T 3.40, log R -3)", &
              0.70d0, 0.02d0, 3.40d0, -3.0d0)

! --- the production combination (OPAL95 at Z = 0.016232, Alexander 2006,
! pure Z, Potekhin conduction): derivatives, the ramp, the seams
         call load_production()
         call derivative_check("OPAL95 2D: dln kappa/dln rho, dln T against differences", &
              0.70d0, 0.016232d0, 5.53d0, -3.3d0, 1.0d-6)
         call derivative_check("OPAL95 3D: dln kappa/dln rho, dln T against differences", &
              0.62d0, 0.016232d0, 5.53d0, -3.3d0, 1.0d-6)
         call derivative_check("OPAL95 4D: dln kappa/dln rho, dln T against differences", &
              0.70d0, 0.0185d0, 6.27d0, -2.2d0, 1.0d-6)
         call derivative_check("Alexander 2006: dln kappa/dln rho, dln T against differences", &
              0.70d0, 0.016232d0, 3.62d0, -3.3d0, 1.0d-5)
         fxion = [1.0d0, 0.0d0, 1.0d0]
         call derivative_check("Potekhin conduction (ionized): derivatives against differences", &
              0.70d0, 0.016232d0, 7.33d0, 2.93d0, 1.0d-5, logd_is_log_rho=.true.)
         fxion = 0.0d0
         call derivative_check("helium core (pure Z blend): derivatives against differences", &
              0.0d0, 0.60d0, 7.83d0, 1.0d0, 5.0d-2, logd_is_log_rho=.true.)
         call continuity_check("ramp: kappa continuous at molecular_opacity_logt_min", &
              0.70d0, 0.016232d0, 4.0d0, -3.3d0, 1.0d-7)
         call continuity_check("ramp: kappa continuous at molecular_opacity_logt_max", &
              0.70d0, 0.016232d0, 4.1d0, -3.3d0, 1.0d-7)
! inside the ramp, off the tables' nodes (log T 4.07): kappa's own
! derivatives, with the weight's
         call derivative_check("ramp: dln kappa/dln rho, dln T against differences", &
              0.70d0, 0.016232d0, 4.07d0, -3.3d0, 1.0d-5)
! across the table's lowest log R (-8), below which log kappa is
! extrapolated linearly
         call get_point(0.70d0, 0.016232d0, 5.53d0, logd_of(5.53d0, -8.0d0 - 1.0d-9), r, ierr)
         call get_point(0.70d0, 0.016232d0, 5.53d0, logd_of(5.53d0, -8.0d0 + 1.0d-9), r2, ierr2)
         call report("OPAL95: kappa continuous at the table's lowest log R", &
              abs(r2(i_kap)/r(i_kap) - 1.0d0), 1.0d-7, ierr + ierr2 == 0)
! the 2D path at the surface composition against the 4D path a relative
! 2e-4 off in Z (the switch's tolerance is 1e-4): |d log kappa| below
! 10 times the change in log Z
         call get_point(0.70d0, 0.016232d0, 5.53d0, logd_of(5.53d0, -3.3d0), r, ierr)
         call get_point(0.70d0, 0.016232d0*(1.0d0 + 2.0d-4), 5.53d0, logd_of(5.53d0, -3.3d0), r2, ierr2)
         err = abs(r(i_log10_kap) - r2(i_log10_kap))
         call report("OPAL95: the 2D path continuous with the 4D path in Z", err, &
              10.0d0*log10(1.0d0 + 2.0d-4), ierr + ierr2 == 0)
! kap_update_surface_tables: the 2D path at a new surface X is the 3D
! path there before the update (one interpolation, two orders)
         call get_point(0.65d0, 0.016232d0, 5.53d0, logd_of(5.53d0, -3.3d0), r, ierr)
         call kap_update_surface_tables(0.65d0, ierr2)
         call get_point(0.65d0, 0.016232d0, 5.53d0, logd_of(5.53d0, -3.3d0), r2, ierr)
         err = abs(r(i_log10_kap) - r2(i_log10_kap))
         call report("kap_update_surface_tables: the new 2D path is the old 3D path", err, &
              1.0d-12, ierr + ierr2 == 0)
         call kap_update_surface_tables(0.70d0, ierr2)
         call conductive_mixture_check()
         call conductive_low_density_check()

! --- other families' derivatives
         call configure("opal92 checks", opal92=.true., z1=0.0188d0, z2=0.0200d0, two_z=.true., &
              opal92_path=trim(atom_dir)//"legacy/opal92/Z0188.OPALN", &
              opal92_path2=trim(atom_dir)//"legacy/opal92/Z0200.OPALN", quiet=.true.)
         call derivative_check("OPAL92 (two Z): derivatives against differences", &
              0.70d0, 0.019d0, 5.53d0, -3.3d0, 1.0d-5)
! the two-Z blend at the first table's Z is that table's opacity
         call get_point(0.70d0, 0.0188d0, 5.53d0, logd_of(5.53d0, -3.3d0), r, ierr)
         call configure("opal92 one table", opal92=.true., z1=0.0188d0, z2=0.0d0, two_z=.false., &
              opal92_path=trim(atom_dir)//"legacy/opal92/Z0188.OPALN", quiet=.true.)
         call get_point(0.70d0, 0.0188d0, 5.53d0, logd_of(5.53d0, -3.3d0), r2, ierr2)
         err = maxval(abs(r(2:4) - r2(2:4)))
         call report("OPAL92 two-Z blend: at opal_table_z1 it is the first table", err, &
              1.0d-13, ierr + ierr2 == 0)
         call configure("laol89 checks", laol89=.true., z1=0.0188d0, z2=0.0190d0, two_z=.true., &
              laol_path=trim(atom_dir)//"legacy/opal92/Z0188AG89.DBGLAOL", &
              laol_path2=trim(atom_dir)//"legacy/opal92/Z0190AG91.DBGLAOL", quiet=.true.)
         call derivative_check("LAOL89 (secant derivatives): against differences, to 10%", &
              0.70d0, 0.0189d0, 5.53d0, -3.3d0, 1.0d-1)
         call configure("kurucz90 checks", opal95=.true., kurucz90=.true., z1=0.0188d0, &
              z2=0.00595d0, two_z=.true., &
              opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
              kurucz_path=trim(mol_dir)//"legacy/kurmol/kapp00.ross", &
              kurucz_path2=trim(mol_dir)//"legacy/kurmol/kapm05.ross", quiet=.true.)
         call derivative_check("Kurucz 1990 (two Z): derivatives against differences", &
              0.70d0, 0.016232d0, 3.633d0, -3.3d0, 1.0d-5)
         call configure("alex95 checks", opal95=.true., alex95=.true., z1=0.016232d0, z2=0.0d0, &
              two_z=.false., opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
              alex95_paths=alex95_paths(), quiet=.true.)
         call derivative_check("Alexander 1994: derivatives against differences", &
              0.70d0, 0.016232d0, 3.62d0, -3.3d0, 1.0d-5)

! --- the formulas (opacity_blending)
         value = 1.234d0
         call blend_in_z(value, -0.5d0, 0.02d0, 0.0188d0, 0.02d0)
         value2 = 1.234d0
         call blend_in_z(value2, -0.5d0, 0.0188d0, 0.0188d0, 0.02d0)
         err = max(abs(value + 0.5d0), abs(value2 - 1.234d0))
         call report("blend_in_z: each table's value at its own Z", err, 1.0d-15, .true.)
         radiative = kap_result(2.0d0, log10(2.0d0), 1.5d0, -3.0d0)
         conductive = kap_result(2.0d0, log10(2.0d0), -1.0d0, 2.0d0)
         total = radiative
         call add_conduction(total, conductive)
         err = max(abs(total%opacity - 1.0d0), abs(total%dlnkap_dlnrho - 0.25d0), &
              abs(total%dlnkap_dlnt + 0.5d0), abs(total%log10_opacity))
         call report("add_conduction: equal opacities halve kappa, average the slopes", &
              err, 1.0d-15, .true.)
         total = hubbard_lampe_opacity(1.5d0, 6.5d0, 0.7d0)
         hl_plus = hubbard_lampe_opacity(1.5d0 + 1.0d-4, 6.5d0, 0.7d0)
         hl_minus = hubbard_lampe_opacity(1.5d0 - 1.0d-4, 6.5d0, 0.7d0)
         err = abs(total%dlnkap_dlnrho - (hl_plus%log10_opacity - hl_minus%log10_opacity)/2.0d-4)
         hl_plus = hubbard_lampe_opacity(1.5d0, 6.5d0 + 1.0d-4, 0.7d0)
         hl_minus = hubbard_lampe_opacity(1.5d0, 6.5d0 - 1.0d-4, 0.7d0)
         err = max(err, abs(total%dlnkap_dlnt - (hl_plus%log10_opacity - hl_minus%log10_opacity)/2.0d-4))
         call report("hubbard_lampe_opacity: derivatives against differences", err, 1.0d-8, .true.)

! --- the error contracts
         call kap_init_with(ierr, pure_z=.true.)
         call report("kap_init refuses the pure-Z table without an interior table", 0.0d0, &
              0.0d0, ierr == 1)
         call kap_init_with(ierr, opal92=.true., two_z=.true.)
         call report("kap_init refuses two OPAL92 tables at one Z", 0.0d0, 0.0d0, ierr == 1)
         call kap_init_with(ierr, opal95=.true., alex06=.true., logt_min=4.1d0)
         call report("kap_init refuses molecular_opacity_logt_min >= _max", 0.0d0, 0.0d0, ierr == 1)
         call kap_init_with(ierr, opal95=.true., conductive_path="no_such_table.cond")
         inquire(file="no_such_table.cond", exist=exists)
         call report("kap_init: a missing table returns ierr (and is not created)", 0.0d0, &
              0.0d0, ierr == 1 .and. .not. exists)
         open(newunit=u, file="short_table.opal92", status="replace")
         block
            integer :: v, ios
            open(newunit=v, file=trim(atom_dir)//"legacy/opal92/Z0188.OPALN", status="old", action="read")
            do k = 1, 157
               read(v,'(a)',iostat=ios) line
               write(u,'(a)') trim(line)
            end do
            close(v)
         end block
         close(u)
         call kap_init_with(ierr, opal92=.true., opal92_path="short_table.opal92")
         call report("kap_init: a short OPAL92 table returns ierr", 0.0d0, 0.0d0, ierr == 1)
         call load_production()
         call get_point(0.70d0, 0.016232d0, 3.0d0, logd_of(3.0d0, -3.0d0), r, ierr)
         r_off = r
         call get_point(0.70d0, 0.016232d0, 8.9d0, 0.0d0, r, ierr)
         call report("kap_get: log T off the tables returns ierr, res unchanged", 0.0d0, 0.0d0, &
              ierr == 1 .and. all(r == r_off))
         call get_point(0.70d0, 0.15d0, 6.0d0, -3.0d0, r, ierr)
         call report("kap_get: Z above 0.12 below 10^7 K returns ierr", 0.0d0, 0.0d0, ierr == 1)
         call get_point(0.99d0, 0.02d0, 3.5d0, -9.0d0, r, ierr)
         call get_point(0.99d0, 0.02d0, 3.5d0, -9.0d0, r, ierr2)
         call report("Alexander 2006: a refused composition fails again on a retry", 0.0d0, &
              0.0d0, ierr == 1 .and. ierr2 == 1)
! (at log T 3.70, an Alexander 2006 node, where the T stencil and so the
! derivative depend on the warm-start cursor the run starts from)
         call load_production()
         call get_point(0.70d0, 0.016232d0, 3.70d0, -7.0d0, r, ierr)
         call get_point(0.70d0, 0.016232d0, 4.05d0, -3.0d0, r2, ierr2)
         call load_production()
         call get_point(0.70d0, 0.016232d0, 3.70d0, -7.0d0, r2, ierr2)
         call report("kap_init again: the first lookup is the first run's, bit for bit", &
              maxval(abs(r - r2)), 0.0d0, ierr + ierr2 == 0)
         call kap_init_with(ierr, opal95=.true.)
         call get_point(0.0d0, 0.60d0, 7.8d0, 4.0d0, r, ierr)
         call report("kap_get: a helium core without the pure-Z table returns ierr", 0.0d0, &
              0.0d0, ierr == 1)

         if (check_failures > 0) then
            write(*,'(a,i0,a)') "test_kap: ", check_failures, " checks FAIL"
            stop 1
         end if
      end subroutine run_checks

         subroutine report(label, err, tol, ok_otherwise)
            character(len=*), intent(in) :: label
            double precision, intent(in) :: err, tol
            logical, intent(in) :: ok_otherwise
            if (ok_otherwise .and. err <= tol) then
               write(*,'(a,t74,1pe10.2,a)') label, err, "  ok"
            else
               write(*,'(a,t74,1pe10.2,a)') label, err, "  FAIL"
               check_failures = check_failures + 1
            end if
         end subroutine report

! log rho at log T and log R = log rho - 3 log T6
         double precision function logd_of(logt, logr)
            double precision, intent(in) :: logt, logr
            logd_of = logr + 3.0d0*(logt - 6.0d0)
         end function logd_of

         subroutine get_point(x, z, logt, logd, res, ierr)
            double precision, intent(in) :: x, z, logt, logd
            double precision, intent(inout) :: res(num_kap_results)
            integer, intent(out) :: ierr
            call kap_get(logd, logt, x, z, res, fxion, ierr)
         end subroutine get_point

! the node's log kappa is the file's, a multiple of 0.001
         subroutine node_check(label, x, z, logt, logr)
            character(len=*), intent(in) :: label
            double precision, intent(in) :: x, z, logt, logr
            double precision :: res(num_kap_results)
            integer :: ierr
            res = 0.0d0
            call get_point(x, z, logt, logd_of(logt, logr), res, ierr)
            call report(label, abs(res(i_log10_kap) - nint(1000.0d0*res(i_log10_kap))/1000.0d0), &
                 1.0d-12, ierr == 0)
         end subroutine node_check

! the log derivatives against centred differences of log10 kappa in
! log10 rho and log10 T (step 1e-4), the larger relative error; logd is
! log R unless logd_is_log_rho
         subroutine derivative_check(label, x, z, logt, logd_in, tol, logd_is_log_rho)
            character(len=*), intent(in) :: label
            double precision, intent(in) :: x, z, logt, logd_in, tol
            logical, intent(in), optional :: logd_is_log_rho
            double precision, parameter :: h = 1.0d-4
            double precision :: res(num_kap_results), plus(num_kap_results), &
                 minus(num_kap_results), logd, err_rho, err_t
            integer :: ierr, jerr
            logd = logd_of(logt, logd_in)
            if (present(logd_is_log_rho)) logd = logd_in
            call get_point(x, z, logt, logd, res, ierr)
            call get_point(x, z, logt, logd + h, plus, jerr)
            ierr = ierr + jerr
            call get_point(x, z, logt, logd - h, minus, jerr)
            ierr = ierr + jerr
            err_rho = abs(res(i_dlnkap_dlnrho) - (plus(i_log10_kap) - minus(i_log10_kap))/(2.0d0*h)) &
                 / max(1.0d0, abs(res(i_dlnkap_dlnrho)))
            call get_point(x, z, logt + h, logd, plus, jerr)
            ierr = ierr + jerr
            call get_point(x, z, logt - h, logd, minus, jerr)
            ierr = ierr + jerr
            err_t = abs(res(i_dlnkap_dlnt) - (plus(i_log10_kap) - minus(i_log10_kap))/(2.0d0*h)) &
                 / max(1.0d0, abs(res(i_dlnkap_dlnt)))
            call report(label, max(err_rho, err_t), tol, ierr == 0)
         end subroutine derivative_check

! kappa just below and just above logt_edge (by 1e-9 in log T), at fixed
! log R: the relative jump
         subroutine continuity_check(label, x, z, logt_edge, logr, tol)
            character(len=*), intent(in) :: label
            double precision, intent(in) :: x, z, logt_edge, logr, tol
            double precision :: below(num_kap_results), above(num_kap_results)
            integer :: ierr, jerr
            call get_point(x, z, logt_edge - 1.0d-9, logd_of(logt_edge - 1.0d-9, logr), below, ierr)
            call get_point(x, z, logt_edge + 1.0d-9, logd_of(logt_edge + 1.0d-9, logr), above, jerr)
            call report(label, abs(above(i_kap)/below(i_kap) - 1.0d0), tol, ierr + jerr == 0)
         end subroutine continuity_check

! The Potekhin mixture: with each species' weight w_i = its number
! fraction times its ionized fraction (He++ counting twice, over 2, the
! metals' proxy O fully ionized), the mixture's conductive opacity is
! sum w_i kappa_i, kappa_i the species' own (each from a pure
! composition, fully ionized); and the mixture's derivatives against
! differences
         subroutine conductive_mixture_check()
            use opacity_tables, only: conductive_table_set
            use conductive_reader, only: read_conductive_table
            use conductive_opacity, only: compute_conductive_opacity
            type(conductive_table_set), save :: cond
            double precision, parameter :: logt = 7.33d0, logd = 2.93d0, x = 0.5d0, z = 0.1d0
            double precision, parameter :: a_h = 1.008d0, a_he = 4.004d0, a_o = 16.0d0
            double precision :: species(3), mix, o, ol, dr, dt, n(3), w(3), f(3), ions(3)
            logical :: got
            character(len=512) :: msg
            integer :: ierr, k
            call read_conductive_table(cond, trim(yrec_input)// &
                 "/opacity/conductive/potekhin/condall06.d", ierr)
            ions = [1.0d0, 0.0d0, 1.0d0]
            call compute_conductive_opacity(cond, logd, logt, 1.0d0, 0.0d0, species(1), ol, dr, &
                 dt, ions, got, msg, k)
            ierr = ierr + k
            call compute_conductive_opacity(cond, logd, logt, 0.0d0, 0.0d0, species(2), ol, dr, &
                 dt, ions, got, msg, k)
            ierr = ierr + k
            call compute_conductive_opacity(cond, logd, logt, 0.0d0, 1.0d0, species(3), ol, dr, &
                 dt, ions, got, msg, k)
            ierr = ierr + k
            f = [0.5d0, 0.3d0, 0.6d0]
            call compute_conductive_opacity(cond, logd, logt, x, z, mix, ol, dr, dt, f, got, msg, k)
            ierr = ierr + k
            n = [x/a_h, (1.0d0 - x - z)/a_he, z/a_o]
            n = n/sum(n)
            w = [n(1)*f(1), n(2)*(f(2) + 2.0d0*f(3))/2.0d0, n(3)]
            call report("Potekhin mixture: kappa_c = sum of w_i kappa_i (partly ionized)", &
                 abs(mix/dot_product(w, species) - 1.0d0), 1.0d-12, ierr == 0)
         end subroutine conductive_mixture_check

! Potekhin conduction below the table's lowest density (log rho = -6,
! where log R >= 0): kappa_c continuous at log rho = -6, and below it
! the value at -6 extrapolated linearly in log rho with its slope (the
! slope returned against differences at log rho = -7), for a partly
! ionized mixture at log T 3.6
         subroutine conductive_low_density_check()
            use opacity_tables, only: conductive_table_set
            use conductive_reader, only: read_conductive_table
            use conductive_opacity, only: compute_conductive_opacity
            type(conductive_table_set), save :: cond
            double precision, parameter :: logt = 3.6d0, x = 0.7d0, z = 0.02d0, h = 1.0d-4
            double precision :: below, above, ol, dr, dt, o, ol_plus, ol_minus, f(3)
            logical :: got, got_all
            character(len=512) :: msg
            integer :: ierr, k
            call read_conductive_table(cond, trim(yrec_input)// &
                 "/opacity/conductive/potekhin/condall06.d", ierr)
            f = [0.5d0, 0.3d0, 0.6d0]
            got_all = .true.
            call compute_conductive_opacity(cond, -6.0d0 - 1.0d-9, logt, x, z, below, ol, dr, dt, &
                 f, got, msg, k)
            ierr = ierr + k
            got_all = got_all .and. got
            call compute_conductive_opacity(cond, -6.0d0 + 1.0d-9, logt, x, z, above, ol, dr, dt, &
                 f, got, msg, k)
            ierr = ierr + k
            got_all = got_all .and. got
            call report("Potekhin conduction: kappa_c continuous at log rho = -6", &
                 abs(above/below - 1.0d0), 1.0d-7, ierr == 0 .and. got_all)
            call compute_conductive_opacity(cond, -7.0d0 + h, logt, x, z, o, ol_plus, dr, dt, &
                 f, got, msg, k)
            ierr = ierr + k
            call compute_conductive_opacity(cond, -7.0d0 - h, logt, x, z, o, ol_minus, dr, dt, &
                 f, got, msg, k)
            ierr = ierr + k
            call compute_conductive_opacity(cond, -7.0d0, logt, x, z, o, ol, dr, dt, f, got, msg, k)
            ierr = ierr + k
            call report("Potekhin conduction below log rho = -6: dln kappa/dln rho", &
                 abs(dr - (ol_plus - ol_minus)/(2.0d0*h))/max(1.0d0, abs(dr)), 1.0d-8, &
                 ierr == 0 .and. got)
         end subroutine conductive_low_density_check

! an OPAL95 node: the lookup's log kappa against the file's entry at
! TABLE #table_number's row log T and column log R
         subroutine opal95_node_check(label, x, z, logt, logr, table_number)
            character(len=*), intent(in) :: label
            double precision, intent(in) :: x, z, logt, logr
            integer, intent(in) :: table_number
            double precision :: res(num_kap_results), row_logt, file_value
            character(len=200) :: line
            integer :: ierr, u, n, ios, column
            res = 0.0d0
            call get_point(x, z, logt, logd_of(logt, logr), res, ierr)
            open(newunit=u, file=trim(atom_dir)//"legacy/opal95/GS98.OP17", status="old", &
                 action="read")
            do
               read(u,'(a)') line
               if (line(1:7) /= 'TABLE #') cycle
               read(line(8:10),*) n
               if (n == table_number) exit
            end do
            do
               read(u,'(a)') line
               read(line(1:4),*,iostat=ios) row_logt
               if (ios == 0 .and. abs(row_logt - logt) < 1.0d-6) exit
            end do
            close(u)
            column = nint((logr + 8.0d0)/0.5d0) + 1
            read(line(5+7*(column-1):4+7*column),*) file_value
            call report(label, abs(res(i_log10_kap) - file_value), 1.0d-12, ierr == 0)
         end subroutine opal95_node_check

         subroutine load_nodes()
            call configure("nodes", opal95=.true., alex06=.true., z1=0.02d0, z2=0.0d0, &
                 two_z=.false., opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
                 alex06_path=trim(mol_dir)//"alexmol06gs98.tab", quiet=.true.)
         end subroutine load_nodes

         subroutine load_production()
            call configure("production", opal95=.true., alex06=.true., pure_z=.true., &
                 conductive=.true., z1=0.016232d0, z2=0.0d0, two_z=.false., &
                 opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
                 alex06_path=trim(mol_dir)//"alexmol06gs98.tab", &
                 pure_z_path=trim(atom_dir)//"lanl/PURECO.DBGLAOL", quiet=.true.)
         end subroutine load_production

         function alex95_paths() result(paths7)
            character(len=256) :: paths7(7)
            paths7 = [character(len=256) :: trim(mol_dir)//"legacy/alex/OPACALEXANDER.X00", &
                 trim(mol_dir)//"legacy/alex/OPACALEXANDER.X01", trim(mol_dir)//"legacy/alex/OPACALEXANDER.X02", &
                 trim(mol_dir)//"legacy/alex/OPACALEXANDER.X035", trim(mol_dir)//"legacy/alex/OPACALEXANDER.X05", &
                 trim(mol_dir)//"legacy/alex/OPACALEXANDER.X07", trim(mol_dir)//"legacy/alex/OPACALEXANDER.X08"]
         end function alex95_paths

! kap_init with the families given, every table's path a shipped table
! (so that only the settings check, or the table named, can fail), for
! the error contracts
         subroutine kap_init_with(ierr, opal95, opal92, alex06, pure_z, two_z, logt_min, &
              conductive_path, opal92_path)
            integer, intent(out) :: ierr
            logical, intent(in), optional :: opal95, opal92, alex06, pure_z, two_z
            double precision, intent(in), optional :: logt_min
            character(len=*), intent(in), optional :: conductive_path, opal92_path
            type(kap_settings_t) :: s
            type(kap_table_paths_t) :: p
            s%use_opal95_tables = present(opal95)
            s%use_opal92_tables = present(opal92)
            s%use_alex06_tables = present(alex06)
            s%use_pure_z_table = present(pure_z)
            s%use_two_z_tables = present(two_z)
            s%use_conductive_opacity = present(conductive_path)
            s%opal95_single_table_z = 0.02d0
            s%opal_table_z1 = 0.02d0
            s%opal_table_z2 = 0.02d0
            s%metal_fraction_match_tolerance = 1.0d-8
            s%molecular_opacity_logt_min = 4.0d0
            if (present(logt_min)) s%molecular_opacity_logt_min = logt_min
            s%molecular_opacity_logt_max = 4.1d0
            s%laol_extrapolation_tolerance = 10.0d0
            p%opal95 = trim(atom_dir)//"legacy/opal95/GS98.OP17"
            p%alex06 = trim(mol_dir)//"alexmol06gs98.tab"
            p%pure_z = trim(atom_dir)//"lanl/PURECO.DBGLAOL"
            p%opal92 = trim(atom_dir)//"legacy/opal92/Z0188.OPALN"
            p%opal92_z2 = trim(atom_dir)//"legacy/opal92/Z0200.OPALN"
            if (present(conductive_path)) p%conductive = conductive_path
            if (present(opal92_path)) p%opal92 = opal92_path
            call kap_init(s, p, 0.70d0, 0.02d0, ierr)
         end subroutine kap_init_with

end program test_kap
