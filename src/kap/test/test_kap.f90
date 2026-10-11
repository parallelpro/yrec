!----------------------------------------------------------------------
! test_kap
!----------------------------------------------------------------------
! Standalone kap-domain test program (2026, phase three -- docs/ROADMAP.md
! stage 2), companion to eos/test/test_eos.f90 (see there for the
! initialization conventions; each piece below cites its source the
! same way). Boots the OPAL95 atomic-opacity configuration used by
! the reference solar models (LOPAL95 with the GS98.OP17 table,
! ZOPAL951=0.016232), evaluates kap_eval over a fixed (logT, logRho)
! grid entirely above the molecular ramp (logT >= TMOLMAX=4.1, so no
! molecular tables are needed), with the conductive correction off,
! and prints the results for byte-comparison against the checked-in
! expected_test_kap.out. Then (2026-10) every opacity family, each
! loaded through kap_init: the production combination (Alexander 2006,
! OPAL95-format, pure Z, Potekhin conductive), OPAL92 and LAOL89 with
! two Z tables, and the Kurucz 1990 and Alexander 1994/95 molecular
! tables (check_families).
program test_kap
      use io_units
      use star_info_lib, only: star
      use kap_lib
      use opacity_tables
      implicit none
      type(kap_settings_t) :: settings
      integer :: gerr

      character(len=256) :: yrec_input
      character(len=256) :: fermi_path, opal95_path, dummy_path
      character(len=256) :: dummy_paths7(7)
      double precision :: laol_work(12)

      integer :: setups_ierr, ipt, kap_ierr
      double precision :: logt, logd, x_frac, z_frac
      double precision :: o, ol, qod, qot, fxion(3)

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
      laol_work = 0.0d0

! unit numbers, per core/read_input.f90
      run_log_unit = 20
      star%ctrl%fermi_unit = 15
      star%ctrl%opal95_table_unit = 48
      open(run_log_unit, file="test_kap.short", status="replace")

! kap configuration: OPAL95 atomic tables only, per the reference
! solar namelists (LOPAL95=T, ZOPAL951=0.016232, TMOLMIN/TMOLMAX
! 4.0/4.1); molecular and conductive paths off
      star%ctrl%use_opal95_tables = .true.
      star%ctrl%use_opal92_tables = .false.
      star%ctrl%use_laol89_tables = .false.
      star%ctrl%use_alex06_tables = .false.
      star%ctrl%use_alex95_tables = .false.
      star%ctrl%use_kurucz90_tables = .false.
      star%use_two_z_tables = .false.
      star%ctrl%use_pure_z_table = .false.
      star%ctrl%use_conductive_opacity = .false.
      star%ctrl%opal95_single_table_z = 0.016232d0
      star%ctrl%molecular_opacity_logt_min = 4.0d0
      star%ctrl%molecular_opacity_logt_max = 4.1d0

! eos side gated off except the always-loaded Fermi table (setups
! calls eos_init unconditionally)
      star%ctrl%use_mhd_eos = .false.
      star%ctrl%use_scv_eos = .false.
      star%ctrl%use_opal95_eos = .false.
      star%ctrl%use_opal2001_eos = .false.
      star%ctrl%use_opal2006_eos = .false.
      star%job%atm_choice = 0

      star%envelope_hydrogen_fraction = 0.70d0
      star%envelope_metal_fraction = 0.016232d0

! constants + table loads (real setups; kap_init inside it reads the
! OPAL95 table and builds the surface-X slice)
      star%job%mixture_weights = laol_work
      star%job%alex06_table_path = dummy_path
      star%job%allard_table_path = dummy_path
      star%job%atm_table_path = dummy_path
      star%job%fermi_table_path = fermi_path
      star%job%kurucz_table_path = dummy_path
      star%job%kurucz_table2_path = dummy_path
      star%job%laol_table_path = dummy_path
      star%job%laol_table2_path = dummy_path
      star%job%opal95_table_path = opal95_path
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
      star%job%conductive_table_path = dummy_path
      star%job%scv_h_table_path = dummy_path
      star%job%scv_he_table_path = dummy_path
      star%job%scv_z_table_path = dummy_path
      star%job%alex95_table_paths = dummy_paths7
      call setups(setups_ierr)
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
         call kap_eval(logd, logt, x_frac, z_frac, o, ol, qod, qot, fxion, gerr)
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
      call kap_eval(logd, logt, x_frac, z_frac, o, ol, qod, qot, fxion, &
           kap_ierr)
      write(*,'(a,i4)') "err out-of-table   ierr = ", kap_ierr
      settings = kap_settings
      settings%use_opal95_tables = .false.
      call set_kap_settings(settings)
      logt = 6.0d0
      logd = -3.0d0
      fxion = 0.0d0
      call kap_eval(logd, logt, x_frac, z_frac, o, ol, qod, qot, fxion, &
           kap_ierr)
      write(*,'(a,i4)') "err no-table-chosen ierr = ", kap_ierr
      settings%use_opal95_tables = .true.
      call set_kap_settings(settings)

      call check_families()

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
!   opal92       two Z tables (0.0188, 0.0200), blended at Z = 0.019
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
           kurucz_path, kurucz_path2, pure_z_path)
         character(len=*), intent(in) :: label
         logical, intent(in), optional :: opal95, opal92, laol89, alex06, alex95, &
              kurucz90, pure_z, conductive
         double precision, intent(in) :: z1, z2
         logical, intent(in) :: two_z
         character(len=*), intent(in), optional :: opal95_path, opal92_path, &
              opal92_path2, laol_path, laol_path2, alex06_path, kurucz_path, &
              kurucz_path2, pure_z_path
         character(len=256), intent(in), optional :: alex95_paths(7)
         character(len=256) :: p(9), p7(7)
         type(kap_settings_t) :: s
         integer :: ierr

         write(*,'(a,a)') "# test_kap: family ", label
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
         p = ""
         p7 = ""
         if (present(alex06_path)) p(1) = alex06_path
         if (present(kurucz_path)) p(2) = kurucz_path
         if (present(kurucz_path2)) p(3) = kurucz_path2
         if (present(laol_path)) p(4) = laol_path
         if (present(laol_path2)) p(5) = laol_path2
         if (present(opal95_path)) p(6) = opal95_path
         if (present(opal92_path)) p(7) = opal92_path
         if (present(opal92_path2)) p(8) = opal92_path2
         if (present(pure_z_path)) p(9) = pure_z_path
         if (present(alex95_paths)) p7 = alex95_paths
         laol_work = 0.0d0
         call kap_init(s, 0.70d0, z1, laol_work, p(1), p(2), p(3), p(4), p(5), &
              p(6), p(7), p(8), p(9), p7, &
              trim(yrec_input)//"/opacity/conductive/potekhin/condall06.d", ierr)
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
         call kap_eval(logd, logt, x, z, o, ol, qod, qot, fxion, ierr)
         if (ierr /= 0) then
            write(*,'(a,3f8.4,a,i0)') "kap X Z logT", x, z, logt, " ierr ", ierr
         else
            write(*,'(a,3f9.5,4(1pe24.15))') "kap", x, z, logt, o, ol, qod, qot
         end if
      end subroutine evaluate_point
end program test_kap
