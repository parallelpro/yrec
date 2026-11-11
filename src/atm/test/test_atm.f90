!----------------------------------------------------------------------
! test_atm
!----------------------------------------------------------------------
! Standalone test of atm/. Two parts:
! The pins, through atm_lib (their output is compared byte for byte):
! each tabulated atmosphere loaded with atm_init and looked up with
! atm_get_surface_pt:
!   atm_choice=3 -- the Kurucz (1993) solar table;
!   atm_choice=5 -- the Castelli & Kurucz (2003) solar table;
!   atm_choice=4 -- the Allard NextGen solar tables: NextGen 2, NextGen 2
!     at tau = 100, the old NextGen 1 format, and a lookup outside the
!     tables;
!   atm_choice=0/1/2 -- the analytic T-tau relations (Eddington,
!     Krishna-Swamy, HSRA) on a grid of optical depths;
! and the Kurucz error path.
! The checks, each printed with its error (or 0/1 for a contract) and ok
! or FAIL against its tolerance, at the kernels (the readers'
! tables, synthetic tables and files) and through atm_lib:
!   the T-tau relations: T = Teff at each photosphere, the start as
!     tau -> 0, the gray law, T rising with tau, the HSRA fit's value,
!     NaN for a choice that is not a relation;
!   the Kurucz and Castelli/Kurucz lookups: exact at the nodes,
!     continuous across the rows, exact on a plane in every branch, and
!     the refusals of M231, M232 and M234;
!   the Allard lookups: exact at the nodes, continuous across the rows,
!     T(tau = 100) > Teff, the outcome's flags, exact on a polynomial
!     (a synthetic file), the low-gravity extension (M233), and the
!     reader's errors (M221, too few rows);
!   atm_lib's contracts: an unknown atm_choice, a table not read, a
!     choice that is not a table, and the choice descriptor and its
!     switch to gray.
! Each check fails under a mutation of the code it covers
! (test_atm.mutations.toml; python3 src/tools/mutate.py atm); the program
! stops with status 1 on any FAIL. The envelope integration that uses
! these lives in core and is covered by the full-model tests.
program test_atm
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      use run_messages, only: run_log_unit, terminal_unit
      use physical_constants, only: cc23
      use math_lib
      use atm_lib
      use atm_tables, only: surface_pressure_table, allard_table, kurucz_num_teff, &
           castelli_num_teff, surface_pressure_num_g
      use kurucz_table_reader, only: read_surface_pressure_table
      use kurucz_atmosphere, only: interpolate_surface_pressure, surface_pressure_outside
      use allard_table_reader, only: read_allard_tables, allard_composition_missing
      use allard_atmosphere, only: interpolate_allard_surface, allard_lookup_result
      implicit none
      integer :: gerr
      character(len=256) :: yrec_input
      character(len=256) :: kurucz_atm_path, castelli_path, allard_path, &
           allard_old_path
      integer :: ipt, atm_ierr
      real(dp) :: teffl, gl, hras, tau, log10_p, log10_t
      integer :: choice, itau
      real(dp), parameter :: grid_tau(6) = [1.0d-4, 1.0d-2, 0.312d0, &
           2.0d0/3.0d0, 1.0d0, 10.0d0]
      logical :: failed
      type(atm_settings_t) :: settings
      integer :: failures = 0

      integer, parameter :: npts = 4
      real(dp), parameter :: grid_teffl(npts) = [3.60d0, 3.70d0, 3.762d0, 3.90d0]
      real(dp), parameter :: grid_gl(npts) = [4.9d0, 4.6d0, 4.438d0, 4.2d0]

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_atm: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      kurucz_atm_path = trim(yrec_input)//"/atmos/kurucz/atmk1990p00.tab"
      castelli_path   = trim(yrec_input)//"/atmos/CastelliKurucz/atmk2004p00.tab"
      allard_path     = trim(yrec_input)//"/atmos/allard/Nextgen2.all"
      allard_old_path = trim(yrec_input)//"/atmos/allard/Nextgen.p00"

! The diagnostics (the table readers' notes, the out-of-table messages)
! go to the run log and terminal_unit; both point at a scratch file so
! that stdout stays deterministic.
      open(newunit=run_log_unit, file="test_atm.short", status="replace")
      terminal_unit = run_log_unit

! --- before any table is read: a lookup in a table atm_init did not
! read is an error, not a read of an empty table
      write(*,'(a)') "# test_atm: before any table is read"
      call atm_get_surface_pt(atm_kurucz, 3.7d0, 4.4d0, log10_p, log10_t, failed, atm_ierr)
      call atm_get_surface_pt(atm_castelli_kurucz, 3.7d0, 4.4d0, log10_p, log10_t, failed, gerr)
      call report_true("atm_get_surface_pt: a table atm_init did not read returns ierr 1", &
           atm_ierr == 1 .and. gerr == 1)

! --- Kurucz (1993), atm_choice=3 ---
      settings%atm_choice = 3
      call atm_init(settings, atm_table_paths_t(kurucz_atm_path, ""), gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Kurucz (choice 3), surfp lookups"
      do ipt = 1, npts
         call lookup(3, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,2(1pe24.15))') "krz ", ipt, log10_p, log10_t
      end do

! --- Castelli & Kurucz (2003), atm_choice=5 ---
      settings%atm_choice = 5
      call atm_init(settings, atm_table_paths_t(castelli_path, ""), gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Castelli/Kurucz (choice 5), " // &
           "kcsurfp lookups"
      do ipt = 1, npts
         call lookup(5, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,2(1pe24.15))') "ck  ", ipt, log10_p, log10_t
      end do

! --- Allard NextGen 2, atm_choice=4 (the Kurucz path is not read) ---
      settings%atm_choice = 4
      call atm_init(settings, atm_table_paths_t(kurucz_atm_path, allard_path), gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Allard (choice 4), " // &
           "atm_get_surface_pt lookups"
      do ipt = 1, npts
         call lookup(4, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,l2,2(1pe24.15))') "al  ", ipt, failed, log10_p, log10_t
      end do

! --- Allard NextGen 2 with P and T at tau = 100 (allard_use_tau100) ---
      settings%allard_use_tau100 = .true.
      call atm_init(settings, atm_table_paths_t(kurucz_atm_path, allard_path), gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Allard at tau = 100, atm_get_surface_pt lookups"
      do ipt = 1, npts
         call lookup(4, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,l2,2(1pe24.15))') "a100", ipt, failed, log10_p, log10_t
      end do
      settings%allard_use_tau100 = .false.

! --- Allard NextGen 1, the old table format ---
      call atm_init(settings, atm_table_paths_t(kurucz_atm_path, allard_old_path), gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Allard NextGen 1 (old format), atm_get_surface_pt lookups"
      do ipt = 1, npts
         call lookup(4, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,l2,2(1pe24.15))') "al1 ", ipt, failed, log10_p, log10_t
      end do
! a hot star is outside the Allard tables: the lookup fails (the caller
! falls back to a gray atmosphere)
      call atm_get_surface_pt(4, 4.5d0, 4.0d0, log10_p, log10_t, failed, gerr)
      write(*,'(a,l2,i4)') "al1 out of range: failed, ierr =", failed, gerr

! --- the analytic T-tau relations (atm_choice 0 Eddington, 1
! Krishna-Swamy, 2 HSRA): log T at a grid of optical depths, the
! starting temperature and the photospheric limit, at log Teff 3.762
      write(*,'(a)') "# test_atm: T-tau relations at log Teff = 3.762"
      hras = hsra_t_tau_offset()
      write(*,'(a,1pe24.15)') "hsra offset ", hras
      do choice = 0, 2
         write(*,'(a,i2,2(1pe24.15))') "ttau start, x_limit", choice, &
              ttau_start_log10_temperature(3.762d0, choice, hras), &
              ttau_photosphere_x_limit(choice)
         do itau = 1, size(grid_tau)
            tau = grid_tau(itau)
            write(*,'(a,i2,2(1pe24.15))') "ttau", choice, tau, &
                 ttau_log10_temperature(tau, 3.762d0, choice, hras)
         end do
      end do

! --- error paths: a Kurucz lookup outside the table (log Teff < 3.5 or
! log g < -0.5) returns ierr = 1 with its diagnostic in the run log; a
! lookup inside it returns 0
      write(*,'(a)') "# test_atm: error paths via ierr"
      settings%atm_choice = 3
      call atm_init(settings, atm_table_paths_t(kurucz_atm_path, ""), gerr)
      if (gerr /= 0) stop 1
      teffl = 3.40d0
      gl = 4.5d0
      call atm_get_surface_pt(3, teffl, gl, log10_p, log10_t, failed, atm_ierr)
      write(*,'(a,i4)') "err out-of-table   ierr = ", atm_ierr
      call atm_get_surface_pt(3, grid_teffl(1), grid_gl(1), log10_p, log10_t, &
           failed, atm_ierr)
      write(*,'(a,i4)') "err facade-success ierr = ", atm_ierr

! --- the checks
      write(*,'(a)') "# test_atm: checks"
      call ttau_checks()
      call surface_pressure_checks()
      call allard_checks()
      call facade_checks()
      call allard_facade_checks()

      close(run_log_unit)
      write(*,'(a)') "test_atm: done"
      if (failures > 0) stop 1
contains

subroutine lookup(atm_choice, log_teff, log_g)
      integer, intent(in) :: atm_choice
      real(dp), intent(in) :: log_teff, log_g
      call atm_get_surface_pt(atm_choice, log_teff, log_g, log10_p, log10_t, &
           failed, gerr)
      if (gerr /= 0) stop 1
end subroutine lookup

! one check: its error against its tolerance
subroutine report(label, err, tolerance)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: err, tolerance
      if (err <= tolerance) then
         write(*,'(a,t84,1pe10.2,a)') label, err, "  ok"
      else
         write(*,'(a,t84,1pe10.2,a)') label, err, "  FAIL"
         failures = failures + 1
      end if
end subroutine report

! a contract check: holds or not
subroutine report_true(label, holds)
      character(len=*), intent(in) :: label
      logical, intent(in) :: holds
      call report(label, merge(0.0d0, 1.0d0, holds), 0.0d0)
end subroutine report_true

! ======================================================================
! The checks (ok/FAIL): each prints its error (or 0/1 for a contract)
! and ok or FAIL against its tolerance; each has a mutation in
! test_atm.mutations.toml that makes it print FAIL.
! ======================================================================

subroutine ttau_checks()
      real(dp), parameter :: log_teff = 3.762d0
      real(dp) :: hsra_offset, err, tau, previous, value
      integer :: choice, k
      logical :: rising

      hsra_offset = hsra_t_tau_offset()
! the photosphere: T = Teff where each relation ends its integration
      call report("ttau eddington: T = Teff at tau = 2/3", &
           abs(ttau_log10_temperature(cc23, log_teff, atm_eddington_gray, 0.0d0) - log_teff), 1.0d-6)
      tau = exp10(ttau_photosphere_x_limit(atm_krishna_swamy))
      call report("ttau krishna-swamy: T = Teff at its photospheric tau", &
           abs(ttau_log10_temperature(tau, log_teff, atm_krishna_swamy, 0.0d0) - log_teff), 1.0d-6)
      call report("ttau hsra: T = Teff at tau = 2/3", &
           abs(ttau_log10_temperature(cc23, log_teff, atm_hsra, hsra_offset) - log_teff), 1.0d-14)
      call report("ttau: the photosphere of Eddington and HSRA is log10 tau = log10(2/3)", &
           max(abs(ttau_photosphere_x_limit(atm_eddington_gray) - log10(cc23)), &
           abs(ttau_photosphere_x_limit(atm_hsra) - log10(cc23))), 1.0d-9)
! the start of the integration: the relation as tau -> 0
      err = 0.0d0
      do choice = atm_eddington_gray, atm_krishna_swamy
         err = max(err, abs(ttau_start_log10_temperature(log_teff, choice, 0.0d0) - &
              ttau_log10_temperature(1.0d-14, log_teff, choice, 0.0d0)))
      end do
      call report("ttau: Eddington's and Krishna-Swamy's start is their T at tau -> 0", err, 1.0d-12)
! the gray law: T^4 linear in tau + 2/3
      err = abs((ttau_log10_temperature(10.0d0, log_teff, atm_eddington_gray, 0.0d0) - &
           ttau_log10_temperature(1.0d0, log_teff, atm_eddington_gray, 0.0d0)) - &
           0.25d0*log10((10.0d0 + cc23)/(1.0d0 + cc23)))
      call report("ttau eddington: T^4 is linear in tau + 2/3 (the gray law)", err, 1.0d-14)
! T rises with depth, for each relation (the HSRA fit has the solar
! temperature minimum near tau = 3e-4, so from tau = 1e-3)
      rising = .true.
      do choice = atm_eddington_gray, atm_hsra
         previous = -huge(1.0d0)
         do k = 0, 40
            tau = exp10(-3.0d0 + 0.1d0*k)
            value = ttau_log10_temperature(tau, log_teff, choice, hsra_offset)
            rising = rising .and. value > previous
            previous = value
         end do
      end do
      call report_true("ttau: T rises with tau from 1e-3 to 10 in each relation", rising)
! the HSRA fit itself, against its pinned value
      call report("ttau hsra: the fit at tau = 1 (log Teff 3.762) is 3.786755164513735", &
           abs(ttau_log10_temperature(1.0d0, log_teff, atm_hsra, hsra_offset) - &
           3.786755164513735d0), 2.0d-15)
! a choice that is not a relation
      call report_true("ttau: a table's choice gives NaN, not a relation", &
           ieee_is_nan(ttau_log10_temperature(1.0d0, log_teff, atm_kurucz, 0.0d0)) .and. &
           ieee_is_nan(ttau_start_log10_temperature(log_teff, atm_allard, 0.0d0)) .and. &
           ieee_is_nan(ttau_photosphere_x_limit(atm_castelli_kurucz)))
end subroutine ttau_checks

! The Kurucz and Castelli/Kurucz lookups, at the kernel: the tables read
! with their reader, and a synthetic table whose log P is a plane.
subroutine surface_pressure_checks()
      type(surface_pressure_table) :: tbl, plane
      real(dp) :: log_p, err, below, above, eps
      integer :: ierr, spline_status, row, column, k, n
      logical :: holds

      call read_surface_pressure_table(kurucz_atm_path, kurucz_num_teff, tbl, ierr)
      if (ierr /= 0) stop 1
! exact at every node of the rows below log Teff 4.5 (cubic splines
! through the nodes)
      err = 0.0d0
      do row = 1, tbl%num_teff
         if (tbl%log_teff(row) > 4.5d0) exit
         do column = tbl%first_column(row), tbl%last_column(row)
            call interpolate_surface_pressure(tbl, tbl%log_teff(row), tbl%log_g(column), &
                 log_p, ierr, spline_status)
            if (ierr /= 0) then
               err = huge(1.0d0)
            else
               err = max(err, abs(log_p - tbl%log10_pressure(row,column)))
            end if
         end do
      end do
      call report("kurucz: exact at every node below log Teff 4.5", err, 1.0d-13)
! continuous across the log Teff nodes, where the stencil moves
      eps = 1.0d-9
      err = 0.0d0
      do row = 3, tbl%num_teff - 3
         call interpolate_surface_pressure(tbl, tbl%log_teff(row) - eps, 4.37d0, below, ierr, spline_status)
         call interpolate_surface_pressure(tbl, tbl%log_teff(row) + eps, 4.37d0, above, ierr, spline_status)
         err = max(err, abs(above - below))
      end do
      call report("kurucz: continuous across the log Teff rows (log g 4.37)", err, 1.0d-6)
! M232: more than one step above the last row or column is refused
      call interpolate_surface_pressure(tbl, 4.64d0, 4.4d0, log_p, ierr, spline_status)
      holds = ierr == surface_pressure_outside
      call interpolate_surface_pressure(tbl, 3.7d0, 5.51d0, log_p, ierr, spline_status)
      holds = holds .and. ierr == surface_pressure_outside
      call interpolate_surface_pressure(tbl, 4.62d0, 5.49d0, log_p, ierr, spline_status)
      holds = holds .and. ierr == 0
      call report_true("kurucz: a point more than a step above the table is refused (M232)", holds)
! M234: a row ending below the top column is exact at its last valid
! node; an all -999 row is refused
      call read_surface_pressure_table(trim(yrec_input)//"/atmos/kurucz/atmk1990m50.tab", &
           kurucz_num_teff, tbl, ierr)
      if (ierr /= 0) stop 1
      row = 1
      column = tbl%last_column(row)
      call interpolate_surface_pressure(tbl, tbl%log_teff(row), tbl%log_g(column), log_p, &
           ierr, spline_status)
      call report("kurucz m50: exact at the last valid node of a row ending at log g 4 (M234)", &
           abs(log_p - tbl%log10_pressure(row,column)), 1.0d-13)
      call read_surface_pressure_table(trim(yrec_input)//"/atmos/kurucz/atmk1990p01.tab", &
           kurucz_num_teff, tbl, ierr)
      if (ierr /= 0) stop 1
      call interpolate_surface_pressure(tbl, 4.04d0, 4.0d0, log_p, ierr, spline_status)
      call report_true("kurucz p01: a stencil on its all -999 row is refused (M234)", &
           ierr == surface_pressure_outside)
! M231: a Castelli/Kurucz stencil reaching the last row (log g 5 only)
      call read_surface_pressure_table(castelli_path, castelli_num_teff, tbl, ierr)
      if (ierr /= 0) stop 1
      call interpolate_surface_pressure(tbl, 4.69d0, 4.9d0, log_p, ierr, spline_status)
      holds = ierr == surface_pressure_outside
      call interpolate_surface_pressure(tbl, 4.68d0, 4.9d0, log_p, ierr, spline_status)
      holds = holds .and. ierr == 0
      call report_true("castelli: a stencil reaching its log-g-5-only last row is refused (M231)", holds)
! exact on a plane: log P = 1 + 2 log Teff - 0.7 log g on the Castelli
! grid, every branch (splines; 3-point Lagrange and linear above log
! Teff 4.5)
      plane = tbl
      do row = 1, plane%num_teff
         plane%log10_pressure(row,:) = 1.0d0 + 2.0d0*plane%log_teff(row) - 0.7d0*plane%log_g
      end do
      plane%first_column = 1
      plane%last_column = surface_pressure_num_g
      err = 0.0d0
      do k = 0, 40
         do n = 0, 10
            call interpolate_surface_pressure(plane, 3.55d0 + 0.0275d0*k, 0.05d0 + 0.49d0*n, &
                 log_p, ierr, spline_status)
            if (ierr /= 0) then
               err = huge(1.0d0)
            else
               err = max(err, abs(log_p - (1.0d0 + 2.0d0*(3.55d0 + 0.0275d0*k) - &
                    0.7d0*(0.05d0 + 0.49d0*n))))
            end if
         end do
      end do
      call report("castelli grid: exact on a plane in log Teff and log g, every branch", err, 1.0d-12)
! a row of the 3-point branch (log Teff 4.5-4.55) without its top three
! columns refuses the point
      do row = 1, plane%num_teff
         if (plane%log_teff(row) > 4.5d0) exit
      end do
      plane%first_column(row) = surface_pressure_num_g - 1
      call interpolate_surface_pressure(plane, plane%log_teff(row), 4.9d0, log_p, ierr, &
           spline_status)
      call report_true("castelli grid: a 3-point row without its top three columns is refused", &
           ierr == surface_pressure_outside)
! the reader's errors: a missing file, a truncated one
      call read_surface_pressure_table("no_such_table.tab", kurucz_num_teff, tbl, ierr)
      holds = ierr == 1 .and. .not. tbl%loaded
      call write_truncated_copy(kurucz_atm_path, "truncated.tab", 20)
      call read_surface_pressure_table("truncated.tab", kurucz_num_teff, tbl, ierr)
      call report_true("kurucz reader: a missing file and a truncated one return ierr 1", &
           holds .and. ierr == 1 .and. .not. tbl%loaded)
end subroutine surface_pressure_checks

! The Allard lookups, at the kernel: Nextgen2.all read by the reader, and
! synthetic files in the Allard layout written here.
subroutine allard_checks()
! (allocatable: the tables are too large for the stack)
      type(allard_table), allocatable :: tbl
      type(allard_lookup_result) :: lookup
      real(dp) :: err, below, above, eps, log_teff, log_g
      integer :: ierr, row, column, k, n, num_inside
      logical :: holds

      allocate(tbl)
      call read_allard_tables(allard_path, 0.0d0, 0.0d0, tbl, ierr)
      if (ierr /= 0) stop 1
! exact at the nodes (4-point Lagrange through them), where the stencil's
! rows all hold the node's log g
      err = 0.0d0
      num_inside = 0
      do row = 1, tbl%num_teff
         do column = tbl%first_column(row), tbl%last_column(row)
            call interpolate_allard_surface(tbl, .false., tbl%log_teff(row), tbl%log_g(column), lookup)
            if (lookup%teff_above .or. lookup%gravity_outside .or. lookup%row_outside > 0) cycle
            if (lookup%teff_below .or. lookup%polint_status /= 0) then
               err = huge(1.0d0)
               cycle
            end if
            num_inside = num_inside + 1
            err = max(err, abs(lookup%log10_pressure - tbl%log10_pressure(row,column)))
         end do
      end do
      if (num_inside < 100) err = huge(1.0d0)
      call report("allard: exact at the nodes of Nextgen2 (log P at T = Teff)", err, 1.0d-12)
! continuous across the log Teff rows
      eps = 1.0d-9
      err = 0.0d0
      do row = 3, tbl%num_teff - 3
         call interpolate_allard_surface(tbl, .false., tbl%log_teff(row) - eps, 4.5d0, lookup)
         below = lookup%log10_pressure
         call interpolate_allard_surface(tbl, .false., tbl%log_teff(row) + eps, 4.5d0, lookup)
         above = lookup%log10_pressure
         err = max(err, abs(above - below))
      end do
      call report("allard: continuous across the log Teff rows (log g 4.5)", err, 1.0d-6)
! tau = 100 lies below the photosphere: T(tau = 100) > Teff
      holds = .true.
      do k = 0, 20
         log_teff = 3.4d0 + 0.03d0*k
         call interpolate_allard_surface(tbl, .true., log_teff, 4.5d0, lookup)
         holds = holds .and. lookup%polint_status == 0 .and. lookup%log10_temperature > log_teff
      end do
      call report_true("allard: T at tau = 100 is above Teff (log Teff 3.4-4.0, log g 4.5)", holds)
! the outcome's contract: above, outside and below the tables
      call interpolate_allard_surface(tbl, .false., tbl%log_teff_max, 4.5d0, lookup)
      holds = lookup%teff_above .and. .not. lookup%gravity_outside
      call interpolate_allard_surface(tbl, .false., 3.6d0, tbl%log_g_max, lookup)
      holds = holds .and. lookup%gravity_outside .and. .not. lookup%teff_above
      call interpolate_allard_surface(tbl, .false., tbl%log_teff_min, 4.5d0, lookup)
      holds = holds .and. lookup%teff_below
      call report_true("allard: above, outside and below the tables are told apart", holds)

! a synthetic file: log P cubic in log g and quadratic in log Teff, which
! 4-point Lagrange reproduces
      call write_synthetic_allard("synthetic_allard.all", 9)
      call read_allard_tables("synthetic_allard.all", 0.0d0, 0.0d0, tbl, ierr)
      if (ierr /= 0) stop 1
      err = 0.0d0
      do k = 0, 20
         do n = 0, 20
            log_teff = log10(3000.0d0) + (log10(5000.0d0) - log10(3000.0d0))*k/20.0d0
            log_g = 3.0d0 + 0.125d0*n
            call interpolate_allard_surface(tbl, .false., log_teff, log_g, lookup)
            if (lookup%teff_above .or. lookup%gravity_outside .or. lookup%teff_below .or. &
                 lookup%row_outside > 0 .or. lookup%polint_status /= 0) then
               err = huge(1.0d0)
            else
               err = max(err, abs(lookup%log10_pressure - synthetic_log_p(log_teff, log_g)))
            end if
         end do
      end do
      call report("allard synthetic: exact on a polynomial (to the file's 9 digits)", err, 1.0d-6)
! the extension: four column widths below a row's first log g, one above
! its last (M233)
      call interpolate_allard_surface(tbl, .false., 3.6d0, 1.01d0, lookup)
      holds = .not. (lookup%gravity_outside .or. lookup%row_outside > 0)
      call interpolate_allard_surface(tbl, .false., 3.6d0, 0.99d0, lookup)
      holds = holds .and. lookup%gravity_outside
      call interpolate_allard_surface(tbl, .false., 3.6d0, 5.99d0, lookup)
      holds = holds .and. .not. (lookup%gravity_outside .or. lookup%row_outside > 0)
      call interpolate_allard_surface(tbl, .false., 3.6d0, 6.01d0, lookup)
      holds = holds .and. lookup%gravity_outside
      call report_true("allard synthetic: lookups reach 4 widths below log g 3, 1 above 5.5 (M233)", holds)
! M221: a composition the file lacks; a table of 3 rows
      call read_allard_tables(allard_path, -0.3d0, 0.0d0, tbl, ierr)
      holds = ierr == allard_composition_missing
      call read_allard_tables(allard_path, 0.0d0, 0.4d0, tbl, ierr)
      call report_true("allard: an [Fe/H] or [alpha/Fe] the file lacks is its own error (M221)", &
           holds .and. ierr == allard_composition_missing)
      call write_synthetic_allard("synthetic_allard_3.all", 3)
      call read_allard_tables("synthetic_allard_3.all", 0.0d0, 0.0d0, tbl, ierr)
      call report_true("allard: a table of 3 rows returns ierr 1", ierr == 1)
! a hole inside a row (one record missing), and a row of three valid
! columns, are refused
      call write_synthetic_allard("synthetic_allard_hole.all", 9, skip_row=5, skip_columns=[3])
      call read_allard_tables("synthetic_allard_hole.all", 0.0d0, 0.0d0, tbl, ierr)
      call report_true("allard: a -999 inside a row is refused", ierr == 1 .and. .not. tbl%loaded)
      call write_synthetic_allard("synthetic_allard_short.all", 9, skip_row=5, &
           skip_columns=[4, 5, 6])
      call read_allard_tables("synthetic_allard_short.all", 0.0d0, 0.0d0, tbl, ierr)
      call report_true("allard: a row of three valid columns is refused", ierr == 1 .and. .not. tbl%loaded)
! the reader's errors: a missing file, an empty one, one in neither
! layout, a record whose log P cannot be read
      call read_allard_tables("no_such_table.all", 0.0d0, 0.0d0, tbl, ierr)
      holds = ierr == 1
      call write_text_file("empty.all", [character(len=1) ::])
      call read_allard_tables("empty.all", 0.0d0, 0.0d0, tbl, ierr)
      holds = holds .and. ierr == 1
      call write_text_file("neither.all", [character(len=40) :: 'not an atmosphere table'])
      call read_allard_tables("neither.all", 0.0d0, 0.0d0, tbl, ierr)
      call report_true("allard reader: a missing file, an empty one, one in neither layout: ierr 1", &
           holds .and. ierr == 1)
      call write_synthetic_allard("synthetic_allard_bad.all", 9, bad_record=.true.)
      call read_allard_tables("synthetic_allard_bad.all", 0.0d0, 0.0d0, tbl, ierr)
      call report_true("allard reader: a record whose log P cannot be read returns ierr 1", &
           ierr == 1 .and. .not. tbl%loaded)
end subroutine allard_checks

! log P of the synthetic Allard file: quadratic in log Teff, cubic in
! log g
pure function synthetic_log_p(log_teff, log_g) result(log_p)
      real(dp), intent(in) :: log_teff, log_g
      real(dp) :: log_p
      log_p = 2.0d0 + 3.0d0*log_teff - 0.5d0*log_g + 0.3d0*log_teff*log_g + &
           0.2d0*log_teff**2 + 0.1d0*log_g**2 - 0.05d0*log_g**3
end function synthetic_log_p

! A file in the Allard layout: num_rows Teff from 3000 K by 250 K, log g
! 3.0-5.5 by 0.5, at [Fe/H] = [alpha/Fe] = 0, and one row at [Fe/H] =
! -1 that the reader must skip; log P, and log P + 1 and log Teff + 0.1
! at tau = 100.
! skip_row, skip_columns: records left out (an entry the file does not
! give); bad_record: one record whose log P is not a number.
subroutine write_synthetic_allard(path, num_rows, skip_row, skip_columns, bad_record)
      character(len=*), intent(in) :: path
      integer, intent(in) :: num_rows
      integer, intent(in), optional :: skip_row, skip_columns(:)
      logical, intent(in), optional :: bad_record
      integer :: unit, row, column
      real(dp) :: teff, log_g
      open(newunit=unit, file=path, status='replace')
      write(unit,'(a)') ' 2 ALLARD # synthetic test table'
      write(unit,'(a)') 'Teff  logg  [M/H] [a/M] log P            log P(tau=100)   log T(tau=100)'
      if (present(bad_record)) write(unit,'(f6.0,3f6.2,a)') 3500.0d0, 4.0d0, 0.0d0, 0.0d0, &
           '  not a number   not a number   not a number'
      do row = num_rows, 1, -1
         teff = 3000.0d0 + 250.0d0*(row - 1)
         do column = 1, 6
            if (present(skip_row) .and. present(skip_columns)) then
               if (row == skip_row .and. any(skip_columns == column)) cycle
            end if
            log_g = 2.5d0 + 0.5d0*column
            write(unit,'(f6.0,3f6.2,1p3d16.8)') teff, log_g, 0.0d0, 0.0d0, &
                 synthetic_log_p(log10(teff), log_g), synthetic_log_p(log10(teff), log_g) + 1.0d0, &
                 log10(teff) + 0.1d0
         end do
      end do
      write(unit,'(f6.0,3f6.2,1p3d16.8)') 4000.0d0, 4.0d0, -1.0d0, 0.0d0, 9.0d0, 9.0d0, 9.0d0
      close(unit)
end subroutine write_synthetic_allard

! a file of the given lines
subroutine write_text_file(path, lines)
      character(len=*), intent(in) :: path, lines(:)
      integer :: unit, k
      open(newunit=unit, file=path, status='replace')
      do k = 1, size(lines)
         write(unit,'(a)') trim(lines(k))
      end do
      close(unit)
end subroutine write_text_file

! the first num_lines lines of the file at source
subroutine write_truncated_copy(source, path, num_lines)
      character(len=*), intent(in) :: source, path
      integer, intent(in) :: num_lines
      character(len=256) :: line
      integer :: in, out, k
      open(newunit=in, file=source, status='old', action='read')
      open(newunit=out, file=path, status='replace')
      do k = 1, num_lines
         read(in,'(a)') line
         write(out,'(a)') trim(line)
      end do
      close(in)
      close(out)
end subroutine write_truncated_copy

! atm_lib's contracts: the settings check, lookups it cannot make, and
! the choice descriptor.
subroutine facade_checks()
      type(atm_settings_t) :: bad
      real(dp) :: log_p, log_t
      logical :: failed, holds
      integer :: ierr, choice

      bad%atm_choice = atm_last_choice + 1
      call atm_init(bad, atm_table_paths_t("", ""), ierr)
      holds = ierr == 1
      bad%atm_choice = -1
      call atm_init(bad, atm_table_paths_t("", ""), ierr)
      call report_true("atm_init: an atm_choice that is none of the choices returns ierr 1 (M236)", &
           holds .and. ierr == 1)
      log_p = 7.0d0
      call atm_get_surface_pt(atm_hsra, 3.7d0, 4.4d0, log_p, log_t, failed, ierr)
      call report_true("atm_get_surface_pt: a choice that is not a table returns ierr 1", &
           ierr == 1 .and. .not. failed .and. log_p == 7.0d0)
! the descriptor
      holds = .true.
      do choice = -2, atm_last_choice + 2
         holds = holds .and. (atm_choice_is_known(choice) .eqv. (choice >= 0 .and. choice <= atm_last_choice))
         holds = holds .and. (atm_choice_is_tabulated(choice) .eqv. &
              any(choice == [atm_kurucz, atm_allard, atm_castelli_kurucz]))
         holds = holds .and. ((len_trim(atm_choice_name(choice)) > 0) .eqv. atm_choice_is_known(choice))
      end do
      call report_true("atm_choices: the known choices, the tables (3-5), a name for each known one", holds)
! the switch to gray and back: Kurucz at log Teff 3.95, Allard at its
! last log Teff, never for the others
      holds = atm_leaves_table(atm_kurucz, 3.95d0) .and. .not. atm_leaves_table(atm_kurucz, 3.9499d0)
      holds = holds .and. atm_reenters_table(atm_kurucz, 3.9499d0) .and. &
           .not. atm_reenters_table(atm_kurucz, 3.95d0)
      holds = holds .and. atm_leaves_table(atm_allard, atm_allard_teff_max()) .and. &
           atm_reenters_table(atm_allard, atm_allard_teff_max() - 1.0d-6)
      do choice = 0, atm_last_choice
         if (choice == atm_kurucz .or. choice == atm_allard) cycle
         holds = holds .and. .not. atm_leaves_table(choice, 9.0d0) .and. &
              .not. atm_reenters_table(choice, 0.0d0)
      end do
      call report_true("atm_choices: Kurucz and Allard switch to gray at their limit, the rest never", holds)
! the terminal's note of each switch (pinned)
      call atm_note_choice_switch(atm_kurucz, atm_eddington_gray, 3.951d0)
      call atm_note_choice_switch(atm_allard, atm_eddington_gray, 4.004d0)
      call atm_note_choice_switch(atm_eddington_gray, atm_kurucz, 3.949d0)
      call atm_note_choice_switch(atm_eddington_gray, atm_allard, 3.999d0)
end subroutine facade_checks

! atm_get_surface_pt on the Allard tables (Nextgen2): a point below them
! is an error; a point inside their log g range but outside a stencil
! row's is lookup_failed, not an error, and leaves the outputs alone
subroutine allard_facade_checks()
      type(allard_table), allocatable :: tbl
      type(allard_lookup_result) :: lookup
      type(atm_settings_t) :: allard_settings
      real(dp) :: log_p, log_t, log_g
      logical :: failed, holds
      integer :: ierr, k

      allocate(tbl)
      call read_allard_tables(allard_path, 0.0d0, 0.0d0, tbl, ierr)
      allard_settings%atm_choice = atm_allard
      call atm_init(allard_settings, atm_table_paths_t("", allard_path), ierr)
      if (ierr /= 0) stop 1
      call atm_get_surface_pt(atm_allard, tbl%log_teff_min - 0.01d0, 4.5d0, log_p, log_t, failed, ierr)
      holds = ierr == 1 .and. .not. failed
! the first log g at log Teff 3.9 that a stencil row's range leaves out
      log_g = 0.0d0
      do k = 0, 100
         call interpolate_allard_surface(tbl, .false., 3.9d0, tbl%log_g_min + 0.05d0*(k + 1), lookup)
         if (lookup%row_outside > 0 .and. .not. lookup%gravity_outside) then
            log_g = tbl%log_g_min + 0.05d0*(k + 1)
            exit
         end if
      end do
      log_p = 7.0d0
      call atm_get_surface_pt(atm_allard, 3.9d0, log_g, log_p, log_t, failed, ierr)
      call report_true("atm_get_surface_pt: below the Allard tables ierr 1, off a row's log g failed", &
           holds .and. log_g > 0.0d0 .and. ierr == 0 .and. failed .and. log_p == 7.0d0)
end subroutine allard_facade_checks

end program test_atm
