!----------------------------------------------------------------------
! test_kap
!----------------------------------------------------------------------
! Standalone kap-domain test program (2026, phase three -- ROADMAP.md
! stage 2), companion to eos/test/test_eos.f90 (see there for the
! initialization conventions; each piece below cites its source the
! same way). Boots the OPAL95 atomic-opacity configuration used by
! the reference solar models (LOPAL95 with the GS98.OP17 table,
! ZOPAL951=0.016232), evaluates kap_get over a fixed (logT, logRho)
! grid entirely above the molecular ramp (logT >= TMOLMAX=4.1, so no
! molecular tables are needed), with the conductive correction off,
! and prints the results for byte-comparison against the checked-in
! expected_test_kap.out.
program test_kap
      use const_lib
      use luout_lib
      use envelope_comp_lib
      use kap_lib
      implicit none

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

! unit numbers, per core/parmin.f90
      short_file_unit = 20
      fermi_unit = 15
      opal95_table_unit = 48
      open(short_file_unit, file="test_kap.short", status="replace")

! kap configuration: OPAL95 atomic tables only, per the reference
! solar namelists (LOPAL95=T, ZOPAL951=0.016232, TMOLMIN/TMOLMAX
! 4.0/4.1); molecular and conductive paths off
      use_opal95_tables = .true.
      use_opal92_tables = .false.
      use_laol89_tables = .false.
      use_alex06_tables = .false.
      use_alex95_tables = .false.
      use_kurucz90_tables = .false.
      use_two_z_tables = .false.
      use_pure_z_table = .false.
      use_conductive_opacity = .false.
      opal95_single_table_z = 0.016232d0
      molecular_opacity_logt_min = 4.0d0
      molecular_opacity_logt_max = 4.1d0

! eos side gated off except the always-loaded Fermi table (setups
! calls eos_init unconditionally)
      use_mhd_eos = .false.
      use_scv_eos = .false.
      use_opal95_eos = .false.
      use_opal2001_eos = .false.
      use_opal2006_eos = .false.
      atm_choice = 0

      env_comp%envelope_hydrogen_fraction = 0.70d0
      env_comp%envelope_metal_fraction = 0.016232d0

! constants + table loads (real setups; kap_init inside it reads the
! OPAL95 table and builds the surface-X slice)
      call setups(laol_work, dummy_path, dummy_path, dummy_path, &
           fermi_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           opal95_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_paths7, setups_ierr)
      if (setups_ierr /= 0) then
         write(*,'(a)') "test_kap: FAIL (setups error)"
         stop 1
      end if

      x_frac = 0.70d0
      z_frac = 0.016232d0

      write(*,'(a)') "# test_kap: kap_get over (logT, logRho) grid, " // &
           "X=0.70 Z=0.016232, OPAL95/GS98.OP17, no molecular, " // &
           "no conductive"
      do ipt = 1, npts
         logt = grid_logt(ipt)
         logd = grid_logd(ipt)
         fxion = 0.0d0
         call kap_get(logd, logt, x_frac, z_frac, o, ol, qod, qot, fxion)
         write(*,'(a,i2,4(1pe24.15))') "kap ", ipt, o, ol, qod, qot
      end do

! Error paths (2026, ROADMAP.md stage 3): with the optional ierr
! passed, out-of-range points and misconfiguration return ierr /= 0
! instead of stopping -- the first time these paths are testable at
! all. The diagnostic each failure writes goes to short_file_unit /
! stdout at the point of failure, as always.
      write(*,'(a)') "# test_kap: error paths via optional ierr"
      logt = 3.5d0
      logd = -8.0d0
      fxion = 0.0d0
      call kap_get(logd, logt, x_frac, z_frac, o, ol, qod, qot, fxion, &
           kap_ierr)
      write(*,'(a,i4)') "err out-of-table   ierr = ", kap_ierr
      use_opal95_tables = .false.
      logt = 6.0d0
      logd = -3.0d0
      fxion = 0.0d0
      call kap_get(logd, logt, x_frac, z_frac, o, ol, qod, qot, fxion, &
           kap_ierr)
      write(*,'(a,i4)') "err no-table-chosen ierr = ", kap_ierr
      use_opal95_tables = .true.

      close(short_file_unit)
      write(*,'(a)') "test_kap: done"
end program test_kap
