!----------------------------------------------------------------------
! test_atm
!----------------------------------------------------------------------
! Standalone atm-domain test program (2026, phase three -- ROADMAP.md
! stage 2), companion to eos/test/test_eos.f90 (see there for the
! initialization conventions). Exercises the surface-boundary-
! condition table machinery for all three tabulated atmosphere
! options the YREC paper names, loading each through atm_init:
!   atm_choice=3 -- Kurucz (1993) solar table, looked up via surfp;
!   atm_choice=5 -- Castelli & Kurucz (2003) solar table, via kcsurfp
!     (kcsurfp reads the choice-5 table loaded by atm_init);
!   atm_choice=4 -- Allard NextGen solar table (loaded together with
!     the Kurucz table, as in production), looked up through the
!     public facade entry atm_get_surface_pt.
! surfp/kcsurfp are atm-internal, called directly here because this
! test lives inside the atm domain -- the same intra-domain access
! the domain's own code has. The full atm_get envelope integration is
! NOT exercised here: it needs the entire eos+kap+solver state booted
! (a full-model concern, covered by Stage-0); noted in ROADMAP.md.
! Results print for byte-comparison against expected_test_atm.out.
program test_atm
      use const_lib
      use luout_lib
      use atm_table_lib
      use atm_lib
      use opacity_table_lib
      use scv_eos_lib
      implicit none

      character(len=256) :: yrec_input
      character(len=256) :: fermi_path, kurucz_atm_path, castelli_path, &
           allard_path, dummy_path
      character(len=256) :: dummy_paths7(7)
      double precision :: laol_work(12)

      integer :: setups_ierr, ipt, atm_ierr
      double precision :: teffl, gl
      logical :: failed

      integer, parameter :: npts = 4
      double precision :: grid_teffl(npts), grid_gl(npts)
      data grid_teffl /3.60d0, 3.70d0, 3.762d0, 3.90d0/
      data grid_gl    /4.9d0,  4.6d0,  4.438d0, 4.2d0/

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_atm: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      fermi_path      = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      kurucz_atm_path = trim(yrec_input)//"/atmos/kurucz/atmk1990p00.tab"
      castelli_path   = trim(yrec_input)//"/atmos/CastelliKurucz/atmk2004p00.tab"
      allard_path     = trim(yrec_input)//"/atmos/allard/Nextgen2.all"
      dummy_path = ""
      dummy_paths7 = ""
      laol_work = 0.0d0

! unit numbers, per core/parmin.f90
      short_file_unit = 20
      fermi_unit = 15
      atm_table_file_unit = 38
      allard_table_unit = 66
      open(short_file_unit, file="test_atm.short", status="replace")
! surfp/kcsurfp write their out-of-table diagnostic to iowr (the main
! output unit in production, per parmin) as well as short_file_unit;
! point both at the scratch .short file so the byte-compared stdout
! stays deterministic when the error paths below fire.
      iowr = short_file_unit

! everything else gated off (setups calls eos_init/kap_init/atm_init
! unconditionally; only the Fermi table is a hard requirement)
      use_mhd_eos = .false.
      use_scv_eos = .false.
      use_opal95_eos = .false.
      use_opal2001_eos = .false.
      use_opal2006_eos = .false.
      use_opal95_tables = .false.
      use_opal92_tables = .false.
      use_laol89_tables = .false.
      use_alex06_tables = .false.
      use_alex95_tables = .false.
      use_kurucz90_tables = .false.
      use_two_z_tables = .false.
      use_pure_z_table = .false.
      use_conductive_opacity = .false.

! constants (real setups; its atm_init call is a no-op at
! atm_choice=0 -- the per-option loads happen explicitly below)
      atm_choice = 0
      call setups(laol_work, dummy_path, dummy_path, dummy_path, &
           fermi_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_paths7, setups_ierr)
      if (setups_ierr /= 0) then
         write(*,'(a)') "test_atm: FAIL (setups error)"
         stop 1
      end if

! --- Kurucz (1993), atm_choice=3, looked up via surfp ---
      atm_choice = 3
      call atm_init(kurucz_atm_path, dummy_path)
      write(*,'(a)') "# test_atm: Kurucz (choice 3), surfp lookups"
      do ipt = 1, npts
         teffl = grid_teffl(ipt)
         gl = grid_gl(ipt)
         call surfp(teffl, gl, .false., atm_ierr)
         write(*,'(a,i2,2(1pe24.15))') "krz ", ipt, &
              atm_table%atm_log10_pressure, atm_table%atm_log10_temperature
      end do

! --- Castelli & Kurucz (2003), atm_choice=5, via kcsurfp ---
      atm_choice = 5
      call atm_init(castelli_path, dummy_path)
      write(*,'(a)') "# test_atm: Castelli/Kurucz (choice 5), " // &
           "kcsurfp lookups"
      do ipt = 1, npts
         teffl = grid_teffl(ipt)
         gl = grid_gl(ipt)
         call kcsurfp(teffl, gl, .false., atm_ierr)
         write(*,'(a,i2,2(1pe24.15))') "ck  ", ipt, &
              atm_table%atm_log10_pressure, atm_table%atm_log10_temperature
      end do

! --- Allard NextGen, atm_choice=4, through the facade entry ---
      atm_choice = 4
      call atm_init(kurucz_atm_path, allard_path)
      write(*,'(a)') "# test_atm: Allard (choice 4), " // &
           "atm_get_surface_pt lookups"
      do ipt = 1, npts
         teffl = grid_teffl(ipt)
         gl = grid_gl(ipt)
         call atm_get_surface_pt(teffl, gl, .false., failed)
         write(*,'(a,i2,l2,2(1pe24.15))') "al  ", ipt, failed, &
              atm_table%atm_log10_pressure, atm_table%atm_log10_temperature
      end do

! Error paths (2026, ROADMAP.md stage 3): surfp's out-of-table check
! (logTeff < 3.5 or logG < -0.5) used to stop the run; with the new
! required ierr it returns 1 -- diagnostics go to iowr/short_file_unit
! as always (both pointed at the scratch .short file here). The
! facade check asserts the optional-ierr success path threads
! ierr = 0 end to end (atm_get_surface_pt -> alsurfp).
      write(*,'(a)') "# test_atm: error paths via ierr"
      teffl = 3.40d0
      gl = 4.5d0
      call surfp(teffl, gl, .false., atm_ierr)
      write(*,'(a,i4)') "err out-of-table   ierr = ", atm_ierr
      teffl = grid_teffl(1)
      gl = grid_gl(1)
      call atm_get_surface_pt(teffl, gl, .false., failed, ierr=atm_ierr)
      write(*,'(a,i4)') "err facade-success ierr = ", atm_ierr

      close(short_file_unit)
      write(*,'(a)') "test_atm: done"
end program test_atm
