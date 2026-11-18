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
      use run_messages, only: run_log_unit, flush_pending_log
      use star_info_lib, only: star_info
      use kap_settings_builder, only: kap_settings_from_star
      use kap_lib, only: kap_init, kap_get, kap_update_surface_tables, &
           kap_set_envelope_metal_fraction, set_kap_settings, kap_settings_t, kap_table_paths_t, &
           num_kap_results, i_kap, i_log10_kap, i_dlnkap_dlnrho, i_dlnkap_dlnt
      use, intrinsic :: iso_fortran_env, only: dp => real64
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
      real(dp) :: logt, logd, x_frac, z_frac
      real(dp) :: o, ol, qod, qot, fxion(3)

! the checks' failure count and table directories (run_checks)
      integer :: check_failures = 0
      character(len=256) :: atom_dir, mol_dir

      integer, parameter :: npts = 7
      real(dp) :: grid_logt(npts), grid_logd(npts)
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

! the run log's unit and file (a run's controls read assigns the unit,
! run_files' assign_unit_numbers)
      run_log_unit = 20
      open(run_log_unit, file="test_kap.short", status="replace")

! kap configuration: OPAL95 atomic tables only, as the reference solar
! decks (use_opal95_tables, opal95_single_table_z = 0.016232,
! molecular_opacity_logt_min and _max = 4.0 and 4.1); no molecular
! family, no conduction
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

! Two error paths: a point off the table, and no interior family, return
! ierr with the message on standard output and in the run log. (The
! header line's "via optional ierr" is pinned text.)
      write(*,'(a)') "# test_kap: error paths via optional ierr"
      logt = 3.5d0
      logd = -8.0d0
      fxion = 0.0d0
      call evaluate(logd, logt, x_frac, z_frac, kap_ierr)
      write(*,'(a,i4)') "err out-of-table   ierr = ", kap_ierr
      settings = kap_settings_from_star(star)
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
         real(dp), intent(in) :: z1, z2
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
         real(dp), intent(in) :: x, z
         logical, intent(in), optional :: interior_only
         real(dp), parameter :: logt_all(9) = [3.40d0, 3.70d0, 4.00d0, &
              4.05d0, 4.10d0, 4.50d0, 5.50d0, 6.50d0, 7.20d0]
         integer :: k, first
         first = 1
         if (present(interior_only)) first = 6
         do k = first, size(logt_all)
            call evaluate_point(x, z, logt_all(k), -3.0d0 + 3.0d0*logt_all(k) - 18.0d0)
         end do
      end subroutine evaluate_points

      subroutine evaluate_point(x, z, logt, logd)
         real(dp), intent(in) :: x, z, logt, logd
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
         real(dp), intent(in) :: logd, logt, x, z
         integer, intent(out) :: ierr
         real(dp) :: res(num_kap_results)
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
         real(dp) :: r(num_kap_results), r2(num_kap_results), r_off(num_kap_results)
         real(dp) :: value, value2, err
         type(kap_result) :: radiative, conductive, total, hl_plus, hl_minus
         integer :: ierr, ierr2, u, log_unit
         logical :: exists

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
! the tenth X table is at X = 1 - Z (TABLE #122 is X = 0.98, Z = 0.02)
         call opal95_node_check("OPAL95 3D path: exact at a node of the X = 1 - Z table", &
              0.98d0, 0.02d0, 6.0d0, -3.0d0, 122)

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
! the same above X = 0.8, where the surface table takes the X tables that
! exist at its Z, the tenth at X = 1 - Z: at X = 0.92 (tables 7 to 10)
! and at Z = 0.06 (no X = 0.95 table: 6, 7, 8 and 10); and at Z = 0.04,
! the first Z whose fixed-Z table has no X = 0.95 slot, against the 4D
! path a relative 2e-4 off in Z, which reads the full set
         call surface_table_check("OPAL95 surface table at X = 0.92: the 2D path is the 3D path", &
              0.016232d0, 0.92d0)
         call surface_table_check("OPAL95 surface table at Z = 0.06, X = 0.85: 2D is 3D", 0.06d0, 0.85d0)
         call configure("surface table", opal95=.true., z1=0.04d0, z2=0.0d0, two_z=.false., &
              opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", quiet=.true.)
         call kap_update_surface_tables(0.85d0, ierr)
         call get_point(0.85d0, 0.04d0, 5.53d0, logd_of(5.53d0, -3.3d0), r, ierr)
         call get_point(0.85d0, 0.04d0*(1.0d0 + 2.0d-4), 5.53d0, logd_of(5.53d0, -3.3d0), r2, ierr2)
         call report("OPAL95 at Z = 0.04, X = 0.85: 2D continuous with 4D in Z", &
              abs(r(i_log10_kap) - r2(i_log10_kap)), 10.0d0*log10(1.0d0 + 2.0d-4), ierr + ierr2 == 0)
         call high_r_shift_check()
         call high_r_edge_checks()
         call precedence_checks()
! without a molecular family the interior opacity is used as it is below
! molecular_opacity_logt_max too (no ramp toward a table that is not
! there): OPAL95 alone, after a lookup elsewhere, at the node log T 4.05
         call configure("opal95 alone", opal95=.true., z1=0.02d0, z2=0.0d0, two_z=.false., &
              opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", quiet=.true.)
         call get_point(0.70d0, 0.02d0, 6.0d0, logd_of(6.0d0, -3.0d0), r, ierr)
         call opal95_node_check("no molecular family: OPAL95 exact at a node below log T 4.1", &
              0.70d0, 0.02d0, 4.05d0, -3.0d0, 73)
         call load_production()
         call conductive_mixture_check()
         call conductive_low_density_check()
         call conductive_no_electrons_check()

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
         call more_family_checks()

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
! (by its message: the first table's file is not at that Z either, which
! its reader would refuse next)
         call report("kap_init refuses two OPAL92 tables at one Z", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'need two Zs') > 0)
         call kap_init_with(ierr, opal95=.true., alex06=.true., logt_min=4.1d0)
         call report("kap_init refuses molecular_opacity_logt_min >= _max", 0.0d0, 0.0d0, ierr == 1)
         call kap_init_with(ierr, opal95=.true., conductive_path="no_such_table.cond")
         inquire(file="no_such_table.cond", exist=exists)
         call report("kap_init: a missing table returns ierr (and is not created)", 0.0d0, &
              0.0d0, ierr == 1 .and. .not. exists)
! the readers' shared open (open_table): a missing OPAL 1992 path
         call kap_init_with(ierr, opal92=.true., opal92_path="no_such_table.opal92")
         inquire(file="no_such_table.opal92", exist=exists)
         call report("kap_init: a missing OPAL92 table returns ierr (open_table)", 0.0d0, &
              0.0d0, ierr == 1 .and. .not. exists .and. &
              index(last_log_line(), 'cannot open the OPAL 1992 table') > 0)
         call copy_table(trim(atom_dir)//"legacy/opal92/Z0188.OPALN", "short_table.opal92", 157)
         call kap_init_with(ierr, opal92=.true., opal92_path="short_table.opal92", table_z=0.0188d0)
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
! what a family checks at its load: the Z its table file states against
! its control's (the first lines of Z0188.OPALN and Z0188AG89.DBGLAOL
! with a table Z of 0.02), the LAOL tolerance, a Kurucz 1990 file of
! header lines only (the copies are short: the message, the run log's
! last line, tells the refusal from a short table's)
         call kap_init_with(ierr, opal92=.true., opal92_path="short_table.opal92")
         call report("kap_init refuses an OPAL92 table that is not at opal_table_z1", 0.0d0, &
              0.0d0, ierr == 1 .and. index(last_log_line(), 'not at opal_table_z1') > 0)
         call copy_table(trim(atom_dir)//"legacy/opal92/Z0188AG89.DBGLAOL", "head_of_table.laol", 20)
         call kap_init_with(ierr, laol_path="head_of_table.laol")
         call report("kap_init refuses a LAOL89 table that is not at laol_table_z1", 0.0d0, &
              0.0d0, ierr == 1 .and. index(last_log_line(), 'not at laol_table_z1') > 0)
         call kap_init_with(ierr, opal95=.true., pure_z=.true., laol_tolerance=0.0d0)
         call report("kap_init refuses a laol_z_tolerance that is not positive", 0.0d0, 0.0d0, &
              ierr == 1)
         open(newunit=u, file="no_rows.kurucz", status="replace")
         write(u,'(a)') "KURUCZ ROSSELAND OPACITIES", " LOG T  LOG P"
         close(u)
         call kap_init_with(ierr, opal95=.true., kurucz_path="no_rows.kurucz")
         call report("kap_init refuses a Kurucz 1990 table without rows", 0.0d0, 0.0d0, ierr == 1)
! one doctored copy per reader, each refused with the reader's own
! message (the run log's last line; the copies are short, so a reader
! that missed its check would end on "ends early"): an OPAL95 table whose
! first table states X = 0.1; an Alexander 2006 file whose first table
! states X = 0.1; an Alexander 1994 X = 0 file whose first table states
! Z = 0.08; a LAOL89 file that states 12 X values; a Kurucz 1990 file
! with a row that is not numbers
         call copy_table(trim(atom_dir)//"legacy/opal95/GS98.OP17", "doctored_header.op17", 250, &
              edit_line=241, edit_column=37, edit_text="0.1000")
         call kap_init_with(ierr, opal95=.true., opal95_path="doctored_header.op17")
         call report("kap_init refuses an OPAL95 table at another X than its slot's", 0.0d0, &
              0.0d0, ierr == 1 .and. index(last_log_line(), 'is expected') > 0)
         call copy_table(trim(mol_dir)//"alexmol06gs98.tab", "doctored_header.alex06", 10, &
              edit_line=1, edit_column=31, edit_text=" 0.100000")
         call kap_init_with(ierr, opal95=.true., alex06=.true., alex06_path="doctored_header.alex06")
         call report("kap_init refuses an Alexander 2006 table at another X", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'expects X = 0.0000 and has 0.1000') > 0)
         call copy_table(trim(mol_dir)//"legacy/alex/OPACALEXANDER.X00", "doctored_header.alex94", 5, &
              edit_line=1, edit_column=27, edit_text="0.08")
         call kap_init_with(ierr, opal95=.true., alex94_first_path="doctored_header.alex94")
         call report("kap_init refuses an Alexander 1994 table at another Z", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'expects Z = 0.1000 and has 0.0800') > 0)
         call copy_table(trim(atom_dir)//"legacy/opal92/Z0188AG89.DBGLAOL", "doctored_size.laol", 20, &
              edit_line=2, edit_column=19, edit_text="12")
         call kap_init_with(ierr, laol_path="doctored_size.laol")
         call report("kap_init refuses a LAOL89 table with more X than fit", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'has 12 X') > 0)
         call copy_table(trim(mol_dir)//"legacy/kurmol/kapp00.ross", "doctored_row.kurucz", 60, &
              edit_line=30, edit_column=1, edit_text=" not a row of numbers")
         call kap_init_with(ierr, opal95=.true., kurucz_path="doctored_row.kurucz")
         call report("kap_init refuses a Kurucz 1990 table with an unreadable row", 0.0d0, &
              0.0d0, ierr == 1 .and. index(last_log_line(), 'unreadable row') > 0)
! more of the readers' refusals, each by its message: an OPAL95 file cut
! inside its first table (300 lines), cut after it (316 lines), and with
! a row that is not numbers; an Alexander 2006 table whose first log R
! column is -8.5 and one whose first row is at log T 4.4; an Alexander
! 1994 row that states 16 values; a pure-Z file that states 2 X values
         call copy_table(trim(atom_dir)//"legacy/opal95/GS98.OP17", "cut_inside.op17", 300)
         call kap_init_with(ierr, opal95=.true., opal95_path="cut_inside.op17")
         call report("kap_init: an OPAL95 file that ends inside a table says so", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'ends inside table 1') > 0)
         call copy_table(trim(atom_dir)//"legacy/opal95/GS98.OP17", "cut_after.op17", 316)
         call kap_init_with(ierr, opal95=.true., opal95_path="cut_after.op17")
         call report("kap_init: an OPAL95 file that ends after a whole table says so", 0.0d0, &
              0.0d0, ierr == 1 .and. index(last_log_line(), 'ends after 1 of its 126 tables') > 0)
         call copy_table(trim(atom_dir)//"legacy/opal95/GS98.OP17", "doctored_row.op17", 320, &
              edit_line=250, edit_column=5, edit_text=" not a number")
         call kap_init_with(ierr, opal95=.true., opal95_path="doctored_row.op17")
         call report("kap_init: an OPAL95 row that is not numbers is refused", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'cannot read a row of table 1') > 0)
         call copy_table(trim(mol_dir)//"alexmol06gs98.tab", "doctored_logr.alex06", 10, &
              edit_line=4, edit_column=7, edit_text=" -8.500")
         call kap_init_with(ierr, opal95=.true., alex06=.true., alex06_path="doctored_logr.alex06")
         call report("kap_init refuses an Alexander 2006 table with other log R columns", 0.0d0, &
              0.0d0, ierr == 1 .and. index(last_log_line(), 'in column 1 and has') > 0)
         call copy_table(trim(mol_dir)//"alexmol06gs98.tab", "doctored_logt.alex06", 10, &
              edit_line=5, edit_column=1, edit_text="4.400")
         call kap_init_with(ierr, opal95=.true., alex06=.true., alex06_path="doctored_logt.alex06")
         call report("kap_init refuses an Alexander 2006 row at another log T", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'in row 85 and has') > 0)
         call copy_table(trim(mol_dir)//"legacy/alex/OPACALEXANDER.X00", "doctored_row.alex94", 5, &
              edit_line=2, edit_column=3, edit_text=" 16")
         call kap_init_with(ierr, opal95=.true., alex94_first_path="doctored_row.alex94")
         call report("kap_init refuses an Alexander 1994 row of another length", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'and has 16 at log T') > 0)
         call copy_table(trim(atom_dir)//"lanl/PURECO.DBGLAOL", "doctored_size.purez", 20, &
              edit_line=2, edit_column=19, edit_text=" 2")
         call kap_init_with(ierr, opal95=.true., pure_z=.true., pure_z_path="doctored_size.purez")
         call report("kap_init refuses a pure-Z table with more than one X", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'has 2 X') > 0)
! where the Kurucz 1990 table has no value (above log rho = -3) the
! interior family's opacity is used alone: OPAL95's (loaded at the same
! Z), bit for bit
         call configure("opal95", opal95=.true., z1=0.0188d0, z2=0.0d0, two_z=.false., &
              opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", quiet=.true.)
         call get_point(0.70d0, 0.016232d0, 4.05d0, -2.5d0, r, ierr)
         call configure("kurucz90", opal95=.true., kurucz90=.true., z1=0.0188d0, &
              z2=0.00595d0, two_z=.true., &
              opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
              kurucz_path=trim(mol_dir)//"legacy/kurmol/kapp00.ross", &
              kurucz_path2=trim(mol_dir)//"legacy/kurmol/kapm05.ross", quiet=.true.)
         call get_point(0.70d0, 0.016232d0, 4.05d0, -2.5d0, r2, ierr2)
         call report("Kurucz 1990 without a value: the interior opacity alone", &
              maxval(abs(r - r2)), 0.0d0, ierr + ierr2 == 0)
! without an interior family the reason is that, whatever the zone's Z
! (here off the envelope's 0.02, with no table to be at it)
         call kap_init_with(ierr, alex06=.true.)
         call get_point(0.70d0, 0.03d0, 6.0d0, -3.0d0, r, ierr)
         call report("kap_get without an interior family says so, at any Z", 0.0d0, 0.0d0, &
              ierr == 1 .and. index(last_log_line(), 'no interior opacity table chosen') > 0)
! a message while the run log is not open (its unit one that nothing has
! opened) creates no file for that unit, and reaches the log once it is
! open (flush_pending_log, which open_run_log calls)
         log_unit = run_log_unit
         run_log_unit = 77
         call kap_init_with(ierr, pure_z=.true.)
         inquire(file="fort.77", exist=exists)
         run_log_unit = log_unit
         call flush_pending_log()
         call report("kap_report: a message before the run log opens waits for it", 0.0d0, &
              0.0d0, ierr == 1 .and. .not. exists .and. &
              index(last_log_line(), 'use_pure_z_table needs an interior table') > 0)
! a single-Z interior table serves the zones at the envelope's Z, which
! follows the model: after the envelope moves to Z = 0.02 a zone there
! has an opacity (the one table's), and one at the old Z has none
         call configure("opal92 one table", opal92=.true., z1=0.0188d0, z2=0.0d0, two_z=.false., &
              opal92_path=trim(atom_dir)//"legacy/opal92/Z0188.OPALN", quiet=.true.)
         call kap_set_envelope_metal_fraction(0.02d0)
         call get_point(0.70d0, 0.02d0, 5.53d0, logd_of(5.53d0, -3.3d0), r, ierr)
         call get_point(0.70d0, 0.0188d0, 5.53d0, logd_of(5.53d0, -3.3d0), r2, ierr2)
         call report("kap_set_envelope_metal_fraction: one Z table serves the new Z only", 0.0d0, &
              0.0d0, ierr == 0 .and. ierr2 == 1)

         if (check_failures > 0) then
            write(*,'(a,i0,a)') "test_kap: ", check_failures, " checks FAIL"
            stop 1
         end if
      end subroutine run_checks

         subroutine report(label, err, tol, ok_otherwise)
            character(len=*), intent(in) :: label
            real(dp), intent(in) :: err, tol
            logical, intent(in) :: ok_otherwise
            if (ok_otherwise .and. err <= tol) then
               write(*,'(a,t74,1pe10.2,a)') label, err, "  ok"
            else
               write(*,'(a,t74,1pe10.2,a)') label, err, "  FAIL"
               check_failures = check_failures + 1
            end if
         end subroutine report

! log rho at log T and log R = log rho - 3 log T6
         real(dp) function logd_of(logt, logr)
            real(dp), intent(in) :: logt, logr
            logd_of = logr + 3.0d0*(logt - 6.0d0)
         end function logd_of

         subroutine get_point(x, z, logt, logd, res, ierr)
            real(dp), intent(in) :: x, z, logt, logd
            real(dp), intent(inout) :: res(num_kap_results)
            integer, intent(out) :: ierr
            call kap_get(logd, logt, x, z, res, fxion, ierr)
         end subroutine get_point

! the node's log kappa is the file's, a multiple of 0.001
         subroutine node_check(label, x, z, logt, logr)
            character(len=*), intent(in) :: label
            real(dp), intent(in) :: x, z, logt, logr
            real(dp) :: res(num_kap_results)
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
            real(dp), intent(in) :: x, z, logt, logd_in, tol
            logical, intent(in), optional :: logd_is_log_rho
            real(dp), parameter :: h = 1.0d-4
            real(dp) :: res(num_kap_results), plus(num_kap_results), &
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
            real(dp), intent(in) :: x, z, logt_edge, logr, tol
            real(dp) :: below(num_kap_results), above(num_kap_results)
            integer :: ierr, jerr
            call get_point(x, z, logt_edge - 1.0d-9, logd_of(logt_edge - 1.0d-9, logr), below, ierr)
            call get_point(x, z, logt_edge + 1.0d-9, logd_of(logt_edge + 1.0d-9, logr), above, jerr)
            call report(label, abs(above(i_kap)/below(i_kap) - 1.0d0), tol, ierr + jerr == 0)
         end subroutine continuity_check

! The Potekhin mixture (Potekhin et al. 1999, section 3.3): 1/K =
! sum x_i/K_i(rho_i), with x_i = n_i Z_i/sum n_k Z_k each species' share
! of the electrons (Z_i the mean charge: H+, He+ + 2 He++, 8 for the
! metals' proxy O) and K_i read at rho_i = rho (A_i/Z_i) Y_e, where the
! pure, fully ionized species has the mixture's electron density Y_e =
! sum X_k Z_k/A_k (A/Z = 1, 2, 2 for the table's H, He and O). Each K_i
! comes from a pure composition, fully ionized, at the density rho'_i
! = rho_i/s_i that its own call shifts to rho_i (s_i = (A_i/Z_i) Y_e of
! the pure species): kappa'_i = 16 sigma T^3/(3 rho'_i K_i(rho_i)), so
! the mixture's kappa_c = sum x_i (rho'_i/rho) kappa'_i.
         subroutine conductive_mixture_check()
            use opacity_tables, only: conductive_table_set
            use conductive_reader, only: read_conductive_table
            use conductive_opacity, only: compute_conductive_opacity
            type(conductive_table_set), save :: cond
            real(dp), parameter :: logt = 7.33d0, logd = 2.93d0, x = 0.5d0, z = 0.1d0
            real(dp), parameter :: a_h = 1.008d0, a_he = 4.004d0, a_o = 16.0d0
            real(dp), parameter :: a_over_z(3) = [1.0d0, 2.0d0, 2.0d0]
            real(dp) :: species(3), mix, ol, dr, dt, n(3), w(3), f(3), ions(3), &
                 charge(3), y_e, pure_y_e(3), logd_pure(3)
            logical :: got
            character(len=512) :: msg
            integer :: ierr, k
            call read_conductive_table(cond, trim(yrec_input)// &
                 "/opacity/conductive/potekhin/condall06.d", ierr)
            f = [0.5d0, 0.3d0, 0.6d0]
            charge = [f(1), f(2) + 2.0d0*f(3), 8.0d0]
            n = [x/a_h, (1.0d0 - x - z)/a_he, z/a_o]
            w = n*charge/sum(n*charge)
            y_e = x*charge(1)/a_h + (1.0d0 - x - z)*charge(2)/a_he + z*charge(3)/a_o
            pure_y_e = [1.0d0/a_h, 2.0d0/a_he, 8.0d0/a_o]
            logd_pure = logd + log10(a_over_z*y_e) - log10(a_over_z*pure_y_e)
            ions = [1.0d0, 0.0d0, 1.0d0]
            call compute_conductive_opacity(cond, logd_pure(1), logt, 1.0d0, 0.0d0, species(1), &
                 ol, dr, dt, ions, got, msg, k)
            ierr = ierr + k
            call compute_conductive_opacity(cond, logd_pure(2), logt, 0.0d0, 0.0d0, species(2), &
                 ol, dr, dt, ions, got, msg, k)
            ierr = ierr + k
            call compute_conductive_opacity(cond, logd_pure(3), logt, 0.0d0, 1.0d0, species(3), &
                 ol, dr, dt, ions, got, msg, k)
            ierr = ierr + k
            call compute_conductive_opacity(cond, logd, logt, x, z, mix, ol, dr, dt, f, got, msg, k)
            ierr = ierr + k
            call report("Potekhin mixture: 1/K = sum x_i/K_i(rho_i), electron fractions", &
                 abs(mix/sum(w*10.0d0**(logd_pure - logd)*species) - 1.0d0), 1.0d-12, ierr == 0)
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
            real(dp), parameter :: logt = 3.6d0, x = 0.7d0, z = 0.02d0, h = 1.0d-4
            real(dp) :: below, above, ol, dr, dt, o, ol_plus, ol_minus, f(3)
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

! A gas without free electrons (X = 1, Z = 0, hydrogen neutral) does
! not conduct: no conductive opacity (got false), ierr 0 and a finite
! value returned, inside the table's reach (log rho -3, log T 3.7)
         subroutine conductive_no_electrons_check()
            use opacity_tables, only: conductive_table_set
            use conductive_reader, only: read_conductive_table
            use conductive_opacity, only: compute_conductive_opacity
            type(conductive_table_set), save :: cond
            real(dp) :: o, ol, dr, dt, f(3)
            logical :: got
            character(len=512) :: msg
            integer :: ierr, k
            call read_conductive_table(cond, trim(yrec_input)// &
                 "/opacity/conductive/potekhin/condall06.d", ierr)
            f = 0.0d0
            call compute_conductive_opacity(cond, -3.0d0, 3.7d0, 1.0d0, 0.0d0, o, ol, dr, dt, &
                 f, got, msg, k)
            call report("Potekhin conduction: none for a gas without free electrons", 0.0d0, &
                 0.0d0, ierr + k == 0 .and. .not. got .and. o == o .and. ol == ol)
         end subroutine conductive_no_electrons_check

! OPAL95 loaded with its fixed-Z table at z (surface X = 0.70): the 3D
! path at (x, z), then the 2D path there once kap_update_surface_tables
! has rebuilt the surface table at x (one interpolation, two orders)
         subroutine surface_table_check(label, z, x)
            character(len=*), intent(in) :: label
            real(dp), intent(in) :: z, x
            real(dp) :: r3d(num_kap_results), r2d(num_kap_results)
            integer :: ierr, ierr2, ierr3
            call configure("surface table", opal95=.true., z1=z, z2=0.0d0, two_z=.false., &
                 opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", quiet=.true.)
            call get_point(x, z, 5.53d0, logd_of(5.53d0, -3.3d0), r3d, ierr)
            call kap_update_surface_tables(x, ierr2)
            call get_point(x, z, 5.53d0, logd_of(5.53d0, -3.3d0), r2d, ierr3)
            call report(label, abs(r3d(i_log10_kap) - r2d(i_log10_kap)), 1.0d-12, &
                 ierr + ierr2 + ierr3 == 0)
         end subroutine surface_table_check

! The high-R shift where the tables' empty cells differ by X: a copy of
! GS98.OP17 whose X = 0.70 tables have no value at log T 6.30 in the last
! log R column. At X = 0.75, log T 6.15, log R 0.3 the stencil's fourth T
! row is that one and its last column the empty one, in one of its four
! X tables (0.5, 0.7, 0.8, 0.9). The lookup must shift that row's log R
! stencil whatever X the previous lookup had (after X = 0.85, whose
! stencil starts at the X = 0.70 table, and after X = 0.20), so that no
! empty cell is read: log kappa within 0.002 of the unedited table's (it
! is 0.013 off with the empty cell in the stencil), and
! the 2D path (the surface table rebuilt at X = 0.75) the 3D path's.
         subroutine high_r_shift_check()
            real(dp), parameter :: x = 0.75d0, z = 0.02d0, logt = 6.15d0, logr = 0.3d0
            real(dp) :: original(num_kap_results), after_high_x(num_kap_results), &
                 after_low_x(num_kap_results), on_surface(num_kap_results), other(num_kap_results)
            integer :: ierr, k
            call doctor_opal95_table("doctored.op17", "X=0.7000", "6.30", 19)
            call configure("original", opal95=.true., z1=z, z2=0.0d0, two_z=.false., &
                 opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", quiet=.true.)
            call get_point(x, z, logt, logd_of(logt, logr), original, ierr)
            call configure("doctored", opal95=.true., z1=z, z2=0.0d0, two_z=.false., &
                 opal95_path="doctored.op17", quiet=.true.)
            call get_point(0.85d0, z, 5.53d0, logd_of(5.53d0, -3.3d0), other, k)
            ierr = ierr + k
            call get_point(x, z, logt, logd_of(logt, logr), after_high_x, k)
            ierr = ierr + k
            call configure("doctored", opal95=.true., z1=z, z2=0.0d0, two_z=.false., &
                 opal95_path="doctored.op17", quiet=.true.)
            call get_point(0.20d0, z, 5.53d0, logd_of(5.53d0, -3.3d0), other, k)
            ierr = ierr + k
            call get_point(x, z, logt, logd_of(logt, logr), after_low_x, k)
            ierr = ierr + k
            call report("OPAL95 high-R shift: independent of the previous lookup's X", &
                 maxval(abs(after_high_x - after_low_x)), 0.0d0, ierr == 0)
            call kap_update_surface_tables(x, k)
            ierr = ierr + k
            call get_point(x, z, logt, logd_of(logt, logr), on_surface, k)
            ierr = ierr + k
            call report("OPAL95 high-R shift: no empty cell of the stencil's X tables read", &
                 abs(after_low_x(i_log10_kap) - original(i_log10_kap)), 2.0d-3, ierr == 0 .and. &
                 abs(on_surface(i_log10_kap) - after_low_x(i_log10_kap)) <= 1.0d-12)
         end subroutine high_r_shift_check

! OPAL95's empty high-R corner, at log T 7.45 (X = 0.70, Z = 0.016232,
! OPAL95 alone): the stencil's rows at log T 7.3 and 7.4 end at log R 0,
! those at 7.5 and 7.6 at -0.5.
!   - log R -0.7: the stencil's top column is empty on the last two rows,
!     so their four columns shift down by one (no extrapolation): the
!     derivatives, and log kappa against its four neighbours' bilinear
!     mean (the nodes log T 7.4 and 7.5, log R -1 and -0.5; the cubic
!     is 0.02 from it), which an empty cell in a stencil would move by
!     tenths of a dex;
!   - log R 0.2: past the last two rows' last column by more than a
!     column, they are extrapolated linearly from it: the derivatives.
         subroutine high_r_edge_checks()
            real(dp), parameter :: x = 0.70d0, z = 0.016232d0, logt = 7.45d0
            real(dp) :: r(num_kap_results), node(num_kap_results)
            real(dp) :: neighbours
            integer :: ierr, k, it, ir
            call configure("opal95 alone", opal95=.true., z1=z, z2=0.0d0, two_z=.false., &
                 opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", quiet=.true.)
            call derivative_check("OPAL95 high-R shift (rows end at different log R): derivatives", &
                 x, z, logt, -0.7d0, 1.0d-6)
            call get_point(x, z, logt, logd_of(logt, -0.7d0), r, ierr)
            neighbours = 0.0d0
            do it = 0, 1
               do ir = 0, 1
                  call get_point(x, z, 7.4d0 + 0.1d0*it, logd_of(7.4d0 + 0.1d0*it, -1.0d0 + 0.5d0*ir), &
                       node, k)
                  ierr = ierr + k
                  neighbours = neighbours + 0.5d0*merge(0.6d0, 0.4d0, ir == 1)*node(i_log10_kap)
               end do
            end do
            call report("OPAL95 high-R shift: log kappa within 0.05 of its neighbours'", &
                 abs(r(i_log10_kap) - neighbours), 5.0d-2, ierr == 0)
            call derivative_check("OPAL95 past a row's last log R (extrapolated): derivatives", &
                 x, z, logt, 0.2d0, 1.0d-6)
         end subroutine high_r_edge_checks

! The precedence lists (opacity_families): with two families of a role
! selected, the lookup is the first one's alone, bit for bit (the second
! is neither loaded nor used).
         subroutine precedence_checks()
            real(dp) :: first(num_kap_results), both(num_kap_results)
            integer :: ierr, k
            call configure("opal95", opal95=.true., z1=0.016232d0, z2=0.0d0, two_z=.false., &
                 opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", quiet=.true.)
            call get_point(0.70d0, 0.016232d0, 5.53d0, logd_of(5.53d0, -3.3d0), first, ierr)
            call configure("opal95 and opal92", opal95=.true., opal92=.true., z1=0.016232d0, &
                 z2=0.0d0, two_z=.false., opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
                 opal92_path=trim(atom_dir)//"legacy/opal92/Z0188.OPALN", quiet=.true.)
            call get_point(0.70d0, 0.016232d0, 5.53d0, logd_of(5.53d0, -3.3d0), both, k)
            call report("precedence: OPAL95 with OPAL92 also selected is OPAL95 alone", &
                 maxval(abs(both - first)), 0.0d0, ierr + k == 0)
            call configure("alex06", opal95=.true., alex06=.true., z1=0.016232d0, z2=0.0d0, &
                 two_z=.false., opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
                 alex06_path=trim(mol_dir)//"alexmol06gs98.tab", quiet=.true.)
            call get_point(0.70d0, 0.016232d0, 3.62d0, logd_of(3.62d0, -3.3d0), first, ierr)
            call configure("alex06 and alex95", opal95=.true., alex06=.true., alex95=.true., &
                 z1=0.016232d0, z2=0.0d0, two_z=.false., &
                 opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
                 alex06_path=trim(mol_dir)//"alexmol06gs98.tab", alex95_paths=alex95_paths(), &
                 quiet=.true.)
            call get_point(0.70d0, 0.016232d0, 3.62d0, logd_of(3.62d0, -3.3d0), both, k)
            call report("precedence: Alexander 2006 with 1994 also selected is 2006's", &
                 maxval(abs(both - first)), 0.0d0, ierr + k == 0)
         end subroutine precedence_checks

! The families the decks do not use, off the pinned points. Called with
! Alexander 1994 over OPAL95 loaded (alex_table_z1 = 0.016232, surface
! X = 0.70).
         subroutine more_family_checks()
            use opacity_tables, only: pure_z_table_set
            use laol_reader, only: read_pure_z_table
            use laol_opacity, only: build_pure_z_splines, compute_pure_z_opacity
            type(pure_z_table_set), save :: pure_z
            real(dp) :: r(num_kap_results), r2(num_kap_results)
            real(dp) :: pz, pz_log, pz_drho, pz_dt
            character(len=256) :: pure_z_path
            character(len=512) :: msg
            integer :: ierr, ierr2, ierr3
! --- Alexander 1994: the path off the surface X; the tables rebuilt at
! another Z (at the nodes Z = 0.02, X = 0.70, log T 3.40, log R -3 the
! file's 3-decimal entry); the surface table refreshed at X = 0.62; and
! above the top log R (1.0), at the surface X and a T node, log kappa
! linear from the last two columns: kappa(1.5) = 2 kappa(1.0) - kappa(0.5)
            call derivative_check("Alexander 1994 off the surface X: derivatives against differences", &
                 0.62d0, 0.016232d0, 3.62d0, -3.3d0, 1.0d-5)
            call node_check("Alexander 1994 at another Z: its tables rebuilt, exact at a node", &
                 0.70d0, 0.02d0, 3.40d0, -3.0d0)
! (a lookup at the surface X first: the lookup at another Z above left
! the tables at that Z, and the first lookup back at this Z rebuilds
! them with the surface table at its own X)
            call get_point(0.70d0, 0.016232d0, 3.62d0, logd_of(3.62d0, -3.3d0), r, ierr)
            call get_point(0.62d0, 0.016232d0, 3.62d0, logd_of(3.62d0, -3.3d0), r, ierr)
            call kap_update_surface_tables(0.62d0, ierr2)
            call get_point(0.62d0, 0.016232d0, 3.62d0, logd_of(3.62d0, -3.3d0), r2, ierr3)
            call report("Alexander 1994 surface table at X = 0.62: 2D is 3D", &
                 abs(r(i_log10_kap) - r2(i_log10_kap)), 1.0d-12, ierr + ierr2 + ierr3 == 0)
            call kap_update_surface_tables(0.70d0, ierr2)
            call top_logr_check("Alexander 1994 above its top log R: linear from the last two")
! --- Alexander 2006: the same above its top log R (a 0.03 Msun deck's
! photosphere is there)
            call configure("alex06", opal95=.true., alex06=.true., z1=0.016232d0, z2=0.0d0, &
                 two_z=.false., opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
                 alex06_path=trim(mol_dir)//"alexmol06gs98.tab", quiet=.true.)
            call top_logr_check("Alexander 2006 above its top log R: linear from the last two")
! --- Kurucz 1990: no row of the table has the density
            call configure("kurucz90", opal95=.true., kurucz90=.true., z1=0.0188d0, z2=0.0d0, &
                 two_z=.false., opal95_path=trim(atom_dir)//"legacy/opal95/GS98.OP17", &
                 kurucz_path=trim(mol_dir)//"legacy/kurmol/kapp00.ross", quiet=.true.)
            call get_point(0.70d0, 0.016232d0, 3.6d0, -25.0d0, r, ierr)
            call report("Kurucz 1990: a density no row has returns ierr", 0.0d0, 0.0d0, &
                 ierr == 1 .and. index(last_log_line(), 'has no value') > 0)
! --- OPAL 1992 (one table at Z = 0.0188, surface X = 0.70; its X tables
! are 0, 0.35 and 0.70): between two X tables; at the grid's X = 0.35,
! one table, against the linear X interpolation 2e-5 beside it; log T
! off the table; and a helium core, where the opacity is the Z blend of
! the table's own (at opal_table_z1) and the pure-Z table's (at Z = 1)
            pure_z_path = trim(atom_dir)//"lanl/PURECO.DBGLAOL"
            call read_pure_z_table(pure_z, pure_z_path, 10.0d0, ierr)
            call build_pure_z_splines(pure_z)
            call compute_pure_z_opacity(pure_z, 2.0d0, 7.8d0, pz, pz_log, pz_drho, pz_dt, msg, ierr2)
            call configure("opal92 and pure Z", opal92=.true., pure_z=.true., z1=0.0188d0, z2=0.0d0, &
                 two_z=.false., opal92_path=trim(atom_dir)//"legacy/opal92/Z0188.OPALN", &
                 pure_z_path=pure_z_path, quiet=.true.)
            call derivative_check("OPAL92 between two X tables: derivatives against differences", &
                 0.50d0, 0.0188d0, 5.53d0, -3.3d0, 1.0d-5)
            call get_point(0.35d0, 0.0188d0, 5.53d0, logd_of(5.53d0, -3.3d0), r, ierr3)
            call get_point(0.35d0 + 2.0d-5, 0.0188d0, 5.53d0, logd_of(5.53d0, -3.3d0), r2, ierr3)
            call report("OPAL92 at an X of its grid: the X interpolation beside it", &
                 abs(r(i_log10_kap) - r2(i_log10_kap)), 1.0d-4, ierr3 == 0)
            call get_point(0.70d0, 0.0188d0, 3.0d0, logd_of(3.0d0, -3.3d0), r, ierr3)
            call report("OPAL92: log T off the table returns ierr", 0.0d0, 0.0d0, &
                 ierr3 == 1 .and. index(last_log_line(), 'is outside the OPAL 1992 tables') > 0)
            call helium_core_blend_check("helium core (OPAL92): the Z blend of its table and the pure-Z", &
                 pz_log, pz_drho, pz_dt, ierr2)
! --- LAOL 1989 (one table at Z = 0.0188): the same helium core
            call configure("laol89 and pure Z", laol89=.true., pure_z=.true., z1=0.0188d0, z2=0.0d0, &
                 two_z=.false., laol_path=trim(atom_dir)//"legacy/opal92/Z0188AG89.DBGLAOL", &
                 pure_z_path=pure_z_path, quiet=.true.)
            call helium_core_blend_check("helium core (LAOL89): the Z blend of its table and the pure-Z", &
                 pz_log, pz_drho, pz_dt, ierr2)
            call laol_column_checks()
         end subroutine more_family_checks

! at the surface X and the T node log T 3.60, above the molecular
! table's top log R (1.0): kappa(1.5) against 2 kappa(1.0) - kappa(0.5)
         subroutine top_logr_check(label)
            character(len=*), intent(in) :: label
            real(dp) :: at_05(num_kap_results), at_10(num_kap_results), at_15(num_kap_results)
            integer :: ierr, k
            call get_point(0.70d0, 0.016232d0, 3.60d0, logd_of(3.60d0, 0.5d0), at_05, ierr)
            call get_point(0.70d0, 0.016232d0, 3.60d0, logd_of(3.60d0, 1.0d0), at_10, k)
            ierr = ierr + k
            call get_point(0.70d0, 0.016232d0, 3.60d0, logd_of(3.60d0, 1.5d0), at_15, k)
            call report(label, abs(at_15(i_log10_kap) - (2.0d0*at_10(i_log10_kap) - at_05(i_log10_kap))), &
                 1.0d-12, ierr + k == 0)
         end subroutine top_logr_check

! X = 0, Z = 0.6 at log T 7.8, log rho 2 (no conduction there: below
! Hubbard & Lampe's density): the interior table's opacity
! there (a lookup at its own Z, 0.0188) blended linearly in Z with the
! pure-Z table's (pz_log, pz_drho, pz_dt, at Z = 1, read with status
! pz_ierr): log kappa and both derivatives
         subroutine helium_core_blend_check(label, pz_log, pz_drho, pz_dt, pz_ierr)
            use opacity_blending, only: blend_in_z
            character(len=*), intent(in) :: label
            real(dp), intent(in) :: pz_log, pz_drho, pz_dt
            integer, intent(in) :: pz_ierr
            real(dp) :: interior(num_kap_results), core(num_kap_results), expected(3)
            integer :: ierr, k
            call get_point(0.0d0, 0.0188d0, 7.8d0, 2.0d0, interior, ierr)
            call get_point(0.0d0, 0.60d0, 7.8d0, 2.0d0, core, k)
            expected = interior(2:4)
            call blend_in_z(expected(1), pz_log, 0.60d0, 0.0188d0, 1.0d0)
            call blend_in_z(expected(2), pz_drho, 0.60d0, 0.0188d0, 1.0d0)
            call blend_in_z(expected(3), pz_dt, 0.60d0, 0.0188d0, 1.0d0)
            call report(label, maxval(abs(core(2:4) - expected)), 1.0d-12, ierr + k + pz_ierr == 0)
         end subroutine helium_core_blend_check

! The LAOL column kernel (the LAOL 1989 and pure-Z lookups' own) on a
! table the test makes: log kappa = 1 + 0.4 log rho - 2.5 (log T - 6) at
! eight densities (log rho -3 to 4) and six temperatures (log T 6 to 7).
! A spline through linear data is that line, and so is the linear
! extrapolation past a row's ends, so inside the table and 1 dex past
! either end log kappa and both derivatives are the line's, exactly;
! beyond the tolerance (ln 10 = 2.3 dex) the lookup returns ierr; a row
! left with three points is refused; and log kappa above 35 is capped.
         subroutine laol_column_checks()
            use opacity_tables, only: pure_z_table_set, laol_opacity_cap
            use laol_opacity, only: build_pure_z_splines, compute_pure_z_opacity
            type(pure_z_table_set), save :: tbl
            real(dp), parameter :: a = 1.0d0, b = 0.4d0, c = -2.5d0
            real(dp) :: o, ol, dr, dt, err, logrho(3)
            character(len=512) :: msg
            integer :: ierr, k, it
            call fill_linear_table(tbl, a, b, c)
            logrho = [0.3d0, -4.0d0, 5.0d0]
            err = 0.0d0
            ierr = 0
            do k = 1, 3
               call compute_pure_z_opacity(tbl, logrho(k), 6.5d0, o, ol, dr, dt, msg, it)
               ierr = ierr + it
               err = max(err, abs(ol - (a + b*logrho(k) + c*0.5d0)), abs(dr - b), abs(dt - c))
            end do
            call report("LAOL column: exact on a linear table, 1 dex past a row's ends", err, 1.0d-11, ierr == 0)
            call compute_pure_z_opacity(tbl, -3.0d0 - 2.4d0, 6.5d0, o, ol, dr, dt, msg, ierr)
            call compute_pure_z_opacity(tbl, -3.0d0 - 2.2d0, 6.5d0, o, ol, dr, dt, msg, k)
            call report("LAOL column: beyond the tolerance returns ierr", 0.0d0, 0.0d0, &
                 ierr == 1 .and. k == 0)
            call fill_linear_table(tbl, a + 40.0d0, b, c)
            call compute_pure_z_opacity(tbl, 0.3d0, 6.5d0, o, ol, dr, dt, msg, ierr)
            call report("LAOL column: log kappa above 35 is capped", 0.0d0, 0.0d0, &
                 ierr == 0 .and. o == laol_opacity_cap .and. ol == 35.0d0)
            call fill_linear_table(tbl, a, b, c, short_row=3)
            call compute_pure_z_opacity(tbl, 0.3d0, 6.5d0, o, ol, dr, dt, msg, ierr)
            call report("LAOL column: a row of fewer than four points returns ierr", 0.0d0, 0.0d0, &
                 ierr == 1 .and. index(msg, 'in log rho') > 0)
         end subroutine laol_column_checks

! tbl made the table with log kappa = offset + b log rho + c (log T - 6)
! at log rho -3 to 4 and log T 6 to 7, splined; with short_row, its T
! row of that index left with three densities
         subroutine fill_linear_table(tbl, offset, b, c, short_row)
            use opacity_tables, only: pure_z_table_set
            use laol_opacity, only: build_pure_z_splines
            type(pure_z_table_set), intent(inout) :: tbl
            real(dp), intent(in) :: offset, b, c
            integer, intent(in), optional :: short_row
            integer :: ir, it
            tbl%num_rho = 8
            tbl%num_t = 6
            tbl%log_extrapolation_tolerance = log(10.0d0)
            tbl%opacity = 0.0d0
            do it = 1, tbl%num_t
               tbl%logt_grid(it) = 10.0d0**(6.0d0 + 0.2d0*(it - 1))
               do ir = 1, tbl%num_rho
                  tbl%rho_grid(ir) = 10.0d0**(-3.0d0 + (ir - 1))
                  tbl%opacity(ir,it) = 10.0d0**(offset + b*(-3.0d0 + (ir - 1)) + c*0.2d0*(it - 1))
               end do
            end do
            if (present(short_row)) tbl%opacity(4:8,short_row) = 0.0d0
            call build_pure_z_splines(tbl)
         end subroutine fill_linear_table

! A copy of GS98.OP17 at target in which, in every table whose header
! holds x_text, the row that starts with row_text has no value (9.999)
! from column first_column to the last
         subroutine doctor_opal95_table(target, x_text, row_text, first_column)
            character(len=*), intent(in) :: target, x_text, row_text
            integer, intent(in) :: first_column
            character(len=200) :: line
            integer :: u, v, ios, column
            logical :: in_table
            in_table = .false.
            open(newunit=u, file=target, status="replace")
            open(newunit=v, file=trim(atom_dir)//"legacy/opal95/GS98.OP17", status="old", &
                 action="read")
            do
               read(v,'(a)',iostat=ios) line
               if (ios /= 0) exit
               if (line(1:7) == 'TABLE #') in_table = index(line, x_text) > 0
               if (in_table .and. line(1:4) == row_text) then
                  do column = first_column, 19
                     line(5+7*(column-1):4+7*column) = '  9.999'
                  end do
               end if
               write(u,'(a)') trim(line)
            end do
            close(v)
            close(u)
         end subroutine doctor_opal95_table

! an OPAL95 node: the lookup's log kappa against the file's entry at
! TABLE #table_number's row log T and column log R
         subroutine opal95_node_check(label, x, z, logt, logr, table_number)
            character(len=*), intent(in) :: label
            real(dp), intent(in) :: x, z, logt, logr
            integer, intent(in) :: table_number
            real(dp) :: res(num_kap_results), row_logt, file_value
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

! The first num_lines lines of the file at source, written to target;
! with edit_line, that line's text from column edit_column on is
! overwritten with edit_text (a doctored header or row).
         subroutine copy_table(source, target, num_lines, edit_line, edit_column, edit_text)
            character(len=*), intent(in) :: source, target
            integer, intent(in) :: num_lines
            integer, intent(in), optional :: edit_line, edit_column
            character(len=*), intent(in), optional :: edit_text
            character(len=200) :: line
            integer :: u, v, k, ios
            open(newunit=u, file=target, status="replace")
            open(newunit=v, file=source, status="old", action="read")
            do k = 1, num_lines
               read(v,'(a)',iostat=ios) line
               if (ios /= 0) exit
               if (present(edit_line)) then
                  if (k == edit_line) line(edit_column:edit_column+len(edit_text)-1) = edit_text
               end if
               write(u,'(a)') trim(line)
            end do
            close(v)
            close(u)
         end subroutine copy_table

! the last line of the run log (kap_report writes each message there
! too): what the last message said
         function last_log_line() result(last)
            character(len=512) :: last, line
            integer :: v, ios
            last = ''
            flush(run_log_unit)
            open(newunit=v, file="test_kap.short", status="old", action="read", iostat=ios)
            if (ios /= 0) return
            do
               read(v,'(a)',iostat=ios) line
               if (ios /= 0) exit
               last = line
            end do
            close(v)
         end function last_log_line

! kap_init with the families given, every table's path a shipped table
! (so that only the settings check, or the table named, can fail), for
! the error contracts. The table Zs are 0.02, or table_z.
         subroutine kap_init_with(ierr, opal95, opal92, alex06, pure_z, two_z, logt_min, &
              conductive_path, opal92_path, laol_path, kurucz_path, laol_tolerance, table_z, &
              opal95_path, alex06_path, alex94_first_path, pure_z_path)
            integer, intent(out) :: ierr
            logical, intent(in), optional :: opal95, opal92, alex06, pure_z, two_z
            real(dp), intent(in), optional :: logt_min, laol_tolerance, table_z
            character(len=*), intent(in), optional :: conductive_path, opal92_path, laol_path, &
                 kurucz_path, opal95_path, alex06_path, alex94_first_path, pure_z_path
            type(kap_settings_t) :: s
            type(kap_table_paths_t) :: p
            s%use_opal95_tables = present(opal95)
            s%use_opal92_tables = present(opal92)
            s%use_laol89_tables = present(laol_path)
            s%use_alex95_tables = present(alex94_first_path)
            s%use_kurucz90_tables = present(kurucz_path)
            s%use_alex06_tables = present(alex06)
            s%use_pure_z_table = present(pure_z)
            s%use_two_z_tables = present(two_z)
            s%use_conductive_opacity = present(conductive_path)
            s%opal95_single_table_z = 0.02d0
            s%opal_table_z1 = 0.02d0
            if (present(table_z)) s%opal_table_z1 = table_z
            s%opal_table_z2 = 0.02d0
            s%laol_table_z1 = 0.02d0
            s%kurucz_table_z1 = 0.02d0
            s%alex_table_z1 = 0.02d0
            s%metal_fraction_match_tolerance = 1.0d-8
            s%molecular_opacity_logt_min = 4.0d0
            if (present(logt_min)) s%molecular_opacity_logt_min = logt_min
            s%molecular_opacity_logt_max = 4.1d0
            s%laol_extrapolation_tolerance = 10.0d0
            if (present(laol_tolerance)) s%laol_extrapolation_tolerance = laol_tolerance
            p%opal95 = trim(atom_dir)//"legacy/opal95/GS98.OP17"
            p%alex06 = trim(mol_dir)//"alexmol06gs98.tab"
            p%pure_z = trim(atom_dir)//"lanl/PURECO.DBGLAOL"
            p%opal92 = trim(atom_dir)//"legacy/opal92/Z0188.OPALN"
            p%opal92_z2 = trim(atom_dir)//"legacy/opal92/Z0200.OPALN"
            if (present(laol_path)) p%laol = laol_path
            if (present(opal95_path)) p%opal95 = opal95_path
            if (present(alex06_path)) p%alex06 = alex06_path
            if (present(pure_z_path)) p%pure_z = pure_z_path
            if (present(alex94_first_path)) then
               p%alex94 = alex95_paths()
               p%alex94(1) = alex94_first_path
            end if
            if (present(conductive_path)) p%conductive = conductive_path
            if (present(opal92_path)) p%opal92 = opal92_path
            if (present(kurucz_path)) p%kurucz = kurucz_path
            call kap_init(s, p, 0.70d0, 0.02d0, ierr)
         end subroutine kap_init_with

end program test_kap
