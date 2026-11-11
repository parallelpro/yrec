!----------------------------------------------------------------------
! opacity_tables
!----------------------------------------------------------------------
! The types of kap/'s state and configuration, and what the families
! share:
!   - one table set per family (opal95_table_set, opal92_table_set,
!     laol_table_set, pure_z_table_set, alex06_table_set,
!     alex94_table_set, kurucz_table_set, conductive_table_set): the
!     tables as read, what is built from them, and the lookup's
!     warm-start cursors and composition caches; opacity_table_state
!     collects them, and kap_lib holds its one instance, which only
!     kap_lib writes (the family routines get their sets as arguments);
!   - the table dimensions, empty-cell markers and composition-cache
!     tolerances;
!   - the settings record (kap_settings_t, and kap_settings, installed by
!     kap_init) and the tables' paths (kap_table_paths_t);
!   - the four-point stencil locators the families share.
! The grids, cursors and caches start from the F77 code's DATA values,
! as component initializers.
module opacity_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
! what the family modules and kap_lib use (each imports its names)
      public :: alex_composition_tol, alex94_surface_slot, &
           kurucz_max_num_densities, kurucz_max_num_temps, kurucz_table_set, &
           laol_log10_opacity_cap, laol_opacity_cap, laol_table_set, n_alex06_d, n_alex06_t, &
           n_alex06_x, n_alex06_z, n_alex94_d, n_alex94_t, n_alex94_x, &
           n_alex94_z, n_laol_rho, n_laol_t, n_laol_x, n_opal92_d, &
           n_opal92_t, n_opal92_x, n_opal95_d, n_opal95_t, n_opal95_tables, &
           n_opal95_x, n_opal95_z, opal92_missing_opacity, opal92_table_set, &
           opal92_x_match_tol, opal95_composition_tol, &
           opal95_missing_opacity, opal95_missing_test, stencil4_locate, &
           stencil4_locate_opal95
! OPAL92 (Livermore) table dimensions
      integer, parameter :: n_opal92_t = 50, n_opal92_d = 17, &
           n_opal92_x = 3, n_opal92_xt = 150, n_opal92_4d = 68
! Alexander 1994/95 table dimensions
      integer, parameter :: n_alex94_x = 7, n_alex94_z = 15, &
           n_alex94_t = 23, n_alex94_d = 17, n_alex94_xt = 8, &
           n_alex94_xz = 105
! the slot of alex94_table_set%opacity that holds the surface-X table
! (after the seven fixed-Z tables, one per X)
      integer, parameter :: alex94_surface_slot = n_alex94_xt
! Alexander 2006 table dimensions
      integer, parameter :: n_alex06_x = 9, n_alex06_z = 16, &
           n_alex06_t = 85, n_alex06_d = 19, n_alex06_xz = 143
! Kurucz table dimensions
      integer, parameter :: kurucz_max_num_temps = 60, &
           kurucz_max_num_densities = 50, kurucz_num_x_tables = 1, &
           kurucz_num_x_temp_entries = kurucz_max_num_temps*kurucz_num_x_tables, &
           kurucz_num_spline_coeffs = 4*kurucz_max_num_densities
! OPAL95 opacity table dimensions. n_opal95_xz is the slot count of the
! unpacked layout (13 Z rows of 10 X slots, slot = table_start_index(iz)
! + ix), not the number of tables in the file (n_opal95_tables = 126).
      integer, parameter :: n_opal95_t = 70, n_opal95_d = 19, &
           n_opal95_x = 10, n_opal95_z = 13, n_opal95_xz = 130
! Number of tables in an OPAL95 table file (read_opal95_tables reads exactly this
! many; 8 X columns at all 13 Z, X=0.95 at 10 Z, X=1-Z at 12 Z).
      integer, parameter :: n_opal95_tables = 126
! Empty-cell markers. OPAL95 cells with no data are stored as
! opal95_missing_opacity and detected with ".ge./.gt. opal95_missing_test".
! OPAL92 tables use opal92_missing_opacity ("<= -9.999" tests).
      real(dp), parameter :: opal95_missing_opacity = 9.999d0
      real(dp), parameter :: opal95_missing_test = 9.9d0
      real(dp), parameter :: opal92_missing_opacity = -9.999d0
! LAOL89 table dimensions (array extents; the file-size checks in
! read_laol_tables/read_pure_z_table accept at most 11 X columns).
      integer, parameter :: n_laol_x = 12, n_laol_rho = 104, n_laol_t = 52
! The LAOL89 lookups' cap: where the interpolated log10 kappa exceeds
! laol_log10_opacity_cap, kappa is laol_opacity_cap (cm^2/g)
      real(dp), parameter :: laol_opacity_cap = 1.0d35, laol_log10_opacity_cap = 35.0d0
! The Potekhin conductive table's grid: log T, log rho, and the ion
! charges (log Z)
      integer, parameter, public :: cond_n_temp = 19, cond_n_rho = 64, &
           cond_n_z = 15
! Composition-cache tolerances: a lookup whose X (and Z) is within this
! of the cached value reuses the previously built fixed-composition
! table. Each family keeps the tolerance its original code used.
      real(dp), parameter :: alex_composition_tol = 1.0d-8
      real(dp), parameter :: opal92_x_match_tol = 1.0d-5
      real(dp), parameter :: opal95_composition_tol = 1.0d-4

! kurucz_table_set: one Kurucz 1990 table (the first or the second Z's):
! its T grid, each row's densities and opacities, the rows' splines in
! log rho, and the lookup's cursors (ix_t, ix_rho). check_range: the
! lookup's early return above log rho = -3 or log T = 4.1, the first
! table's only, as in the F77 code (read_kurucz_tables clears it for the
! second).
      type, public :: kurucz_table_set
           real(dp) :: grid_logt(kurucz_max_num_temps)
! opacity: kappa (cm^2/g) as read (the files' log10, at no
! microturbulence, taken to kappa); log10_rho: the row's densities
           real(dp) :: opacity(kurucz_num_x_temp_entries,kurucz_max_num_densities), &
                log10_rho(kurucz_num_x_temp_entries,kurucz_max_num_densities)
           integer :: num_temps
           integer :: ix_t = 1, ix_rho = 1
           real(dp) :: spline_coeffs(kurucz_num_x_temp_entries,kurucz_num_spline_coeffs)
           integer :: density_start_index(kurucz_num_x_temp_entries), &
                density_count(kurucz_num_x_temp_entries)
           logical :: check_range = .true.
      end type kurucz_table_set

! laol_table_set: one LAOL89 table (the first or the second Z's): the
! table as read (opacity, grid_x, grid_rho in g/cm^3, grid_t, which
! build_laol_splines converts to log10 T in place, and the extents), the
! rows' splines in log rho (slaol_*), and table_number (2 for the second
! Z's, whose run-log lines say "#2" and "SECOND").
      type, public :: laol_table_set
           real(dp) :: opacity(n_laol_x,n_laol_rho,n_laol_t), &
                grid_x(n_laol_x), grid_t(n_laol_t), grid_rho(n_laol_rho)
           integer :: num_x, num_rho, num_t
           real(dp) :: slaol_opacity(n_laol_x,n_laol_rho,n_laol_t), &
                slaol_log_rho(n_laol_x,n_laol_rho,n_laol_t), &
                slaol_d2opacity(n_laol_x,n_laol_rho,n_laol_t)
           integer :: slaol_num_points(n_laol_x,n_laol_t)
           integer :: table_number = 1
! ln of laol_extrapolation_tolerance: how far in log rho the lookup
! extrapolates past a row's ends
           real(dp) :: log_extrapolation_tolerance
      end type laol_table_set

! opal92_table_set: one OPAL92 (Livermore) table (the first or the
! second Z's): its grids, log kappa in (X, T) rows by log R, the rows'
! splines, the surface table at the envelope X (surface_*), and the
! lookup's warm-start cursors index_x, index_t and index_rho (findex
! resets a cursor outside 1..n to 1, so they start at 1).
      type, public :: opal92_table_set
           real(dp) :: grid_logt(n_opal92_t), &
                grid_x(n_opal92_x), grid_logr(n_opal92_d)
           real(dp) :: log10_opacity(n_opal92_xt,n_opal92_d)
           integer :: num_x, num_temps
           real(dp) :: spline_coeffs(n_opal92_xt,n_opal92_4d)
           integer :: density_start_index(n_opal92_xt), &
                density_count(n_opal92_xt)
! former common/llot4/ (the surface-X table built by
! build_opal92_surface_table)
           real(dp) :: surface_x, surface_z
           real(dp) :: surface_spline_coeffs(n_opal92_t,n_opal92_4d)
           integer :: surface_x_index
! former common/kipmll/ (cached grid indices)
           integer :: index_x = 1, index_t = 1, index_rho = 1
      end type opal92_table_set

! opal95_table_set: the OPAL 1995-format tables (former common
! /llot95a/, /llot95/, /llot95e/ and /op95indx/): the full set in X and
! Z, the fixed-Z table built from it at opal95_single_table_z, the
! fixed-X, fixed-Z surface table at the envelope composition, and the
! lookup's warm-start cursors. The grids in X and Z and the slot layout
! are the files' (read_opal95_tables checks each table's X and Z).
      type, public :: opal95_table_set
           real(dp) :: grid_logt(n_opal95_t)
           real(dp) :: grid_x(n_opal95_x) = &
                [0.0d0, 0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0,0.9d0, &
                 0.95d0,1.0d0]
           real(dp) :: grid_logr(n_opal95_d)
           real(dp) :: grid_z(n_opal95_z) = &
                [0.0d0, 0.0001d0, 0.0003d0, 0.001d0, 0.002d0, &
                 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
                 0.04d0, 0.06d0, 0.08d0, 0.10d0]
           real(dp) :: full_opacity(n_opal95_xz,n_opal95_t,n_opal95_d)
           integer :: num_x_at_z(n_opal95_z) = &
                [10,10,10,10,10,10,10,10,10,10,9,9,8]
! the slot of (X index ix, Z index iz) is table_start_index(iz) + ix:
! ten slots per Z row, so the rows with fewer X tables leave slots 109,
! 119, 129 and 130 unwritten (no stencil addresses them)
           integer :: table_start_index(n_opal95_z) = &
                [0,10,20,30,40,50,60,70,80,90,100,110,120]
           real(dp) :: fixed_z_opacity(n_opal95_x,n_opal95_t,n_opal95_d)
           real(dp) :: fixed_z
           real(dp) :: surface_opacity(n_opal95_t,n_opal95_d)
           real(dp) :: surface_x
! the warm-start cursors: the first index of the Z stencil, of each Z
! row's X stencil, of the T stencil and of each T row's log R stencil
           integer :: index_z = 1
           integer :: index_x(4,4) = 1
           integer :: index_t = 1
           integer :: index_rho(4) = 1
      end type opal95_table_set

! alex94_table_set: the seven Alexander & Ferguson (1994) files (former
! common /galot/, /alot/, /alotall/): every table in X and Z
! (full_opacity, slot = X index + 7 (Z index - 1)), the fixed-Z tables
! per X built from them (opacity(1:7,:,:)) and the surface-X table
! (opacity(surface_slot,:,:)); the compositions those were built at, and
! the lookup's cursors.
      type, public :: alex94_table_set
           real(dp) :: grid_logt(n_alex94_t) = &
                [3.00d0,3.05d0,3.10d0,3.15d0,3.20d0,3.25d0,3.30d0, &
                 3.35d0,3.40d0,3.45d0,3.50d0,3.55d0,3.60d0,3.65d0, &
                 3.70d0,3.75d0,3.80d0,3.85d0,3.90d0,3.95d0,4.00d0, &
                 4.05d0,4.10d0]
           real(dp) :: grid_x(n_alex94_x) = &
                [0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0]
           real(dp) :: grid_logr(n_alex94_d) = &
                [-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0,-4.5d0,-4.0d0, &
                 -3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0, &
                 0.0d0, 0.5d0, 1.0d0]
           real(dp) :: grid_z(n_alex94_z) = &
                [0.0d0, 0.00001d0, 0.00003d0, 0.0001d0, 0.0003d0, &
                 0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
                 0.04d0, 0.06d0, 0.08d0, 0.10d0]
           real(dp) :: opacity(n_alex94_xt,n_alex94_t,n_alex94_d)
           real(dp) :: cached_x = 0.0d0, cached_z = 0.0d0
           integer :: index_x = 4, index_t = 12, index_r = 9
           real(dp) :: full_opacity(n_alex94_xz,n_alex94_t,n_alex94_d)
      end type alex94_table_set

! alex06_table_set: the Alexander 2006 (Ferguson et al. 2005) file
! (former common /galot06/, /alot06/, /alot06all/): every table in X and
! Z (full_opacity, slot = n_alex06_z (X index - 1) + Z index), the table
! at one composition built from them (opacity), that composition (the
! cache), and the lookup's cursors.
      type, public :: alex06_table_set
           real(dp) :: grid_logt(n_alex06_t) = &
                [2.70d0,2.75d0,2.80d0,2.85d0,2.90d0,2.91d0,2.92d0, &
                 2.93d0,2.94d0,2.95d0,2.96d0,2.97d0,2.98d0,2.99d0,3.00d0, &
                 3.01d0,3.02d0,3.03d0,3.04d0,3.05d0,3.06d0,3.07d0,3.08d0, &
                 3.09d0,3.10d0,3.11d0,3.12d0,3.13d0,3.14d0,3.15d0,3.16d0, &
                 3.17d0,3.18d0,3.19d0,3.20d0,3.21d0,3.22d0,3.23d0,3.24d0, &
                 3.25d0,3.26d0,3.27d0,3.28d0,3.29d0,3.30d0,3.31d0,3.32d0, &
                 3.33d0,3.34d0,3.35d0,3.36d0,3.37d0,3.38d0,3.39d0,3.40d0, &
                 3.41d0,3.42d0,3.43d0,3.44d0,3.45d0,3.46d0,3.47d0,3.48d0, &
                 3.49d0,3.50d0,3.55d0,3.60d0,3.65d0,3.70d0,3.75d0,3.80d0, &
                 3.85d0,3.90d0,3.95d0,4.00d0,4.05d0,4.10d0,4.15d0,4.20d0, &
                 4.25d0,4.30d0,4.35d0,4.40d0,4.45d0,4.50d0]
           real(dp) :: grid_x(n_alex06_x) = &
                [0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0,0.9d0,1.0d0]
           real(dp) :: grid_logr(n_alex06_d) = &
                [-8.0d0,-7.5d0,-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0, &
                 -4.5d0,-4.0d0,-3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0, &
                 -0.5d0, 0.0d0, 0.5d0, 1.0d0]
           real(dp) :: grid_z(n_alex06_z) = &
                [0.0d0, 0.00001d0, 0.00003d0, 0.0001d0, 0.0003d0, &
                 0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
                 0.04d0, 0.05d0,0.06d0, 0.08d0, 0.10d0]
           real(dp) :: opacity(n_alex06_t,n_alex06_d)
           real(dp) :: cached_x = 0.0d0, cached_z = 0.0d0
           integer :: index_x = 4, index_t = 43, index_r = 10
           real(dp) :: full_opacity(n_alex06_xz,n_alex06_t,n_alex06_d)
      end type alex06_table_set

! pure_z_table_set: the LAOL pure-Z table (Z = 1, pure C and O; former
! common /zlaol/, /zslaol/): the table as read (opacity in density and T,
! the density grid in g/cm^3 and the T grid, converted to log10 T by
! build_pure_z_splines) and its splines in log rho, one per T row.
      type, public :: pure_z_table_set
           real(dp) :: opacity(n_laol_rho,n_laol_t), logt_grid(n_laol_t), &
                rho_grid(n_laol_rho)
           integer :: num_rho, num_t
           real(dp) :: spline_opacity(n_laol_rho,n_laol_t), spline_log_rho(n_laol_rho,n_laol_t), &
                spline_d2opacity(n_laol_rho,n_laol_t)
           integer :: spline_num_points(n_laol_t)
! ln of laol_extrapolation_tolerance: how far in log rho the lookup
! extrapolates past a row's ends
           real(dp) :: log_extrapolation_tolerance
      end type pure_z_table_set

! conductive_table_set: the Potekhin table (condall06.d): log10 of the
! electron thermal conductivity on a grid of log T, log rho and log Z
! (the ion charge).
      type, public :: conductive_table_set
           real(dp) :: temp_grid(cond_n_temp)
           real(dp) :: rho_grid(cond_n_rho)
           real(dp) :: z_grid(cond_n_z)
           real(dp) :: log10_kappa_table(cond_n_temp, cond_n_rho, cond_n_z)
           logical :: loaded = .false.
      end type conductive_table_set

! opacity_table_state: every family's tables. kap_lib holds the one
! instance; the family modules get their sets as arguments.
      type, public :: opacity_table_state
! (the components run from the smallest to the largest: on arm64 Mach-O
! a static instance's components must start within 8 MB of it, and the
! whole is about 8.7 MB)
! the envelope metal fraction the tables were loaded for (kap_init's),
! which a single-Z interior table (OPAL92, LAOL89 without a second Z)
! must match
           real(dp) :: envelope_metal_fraction = 0.0d0
           type(conductive_table_set) :: conductive
! OPAL92 (Livermore) tables: (1) at opal_table_z1, (2) at
! opal_table_z2 (read only when kap_settings%use_two_z_tables). See
! opal92_table_set.
           type(opal92_table_set) :: opal92(2)
! Kurucz90 tables: (1) at kurucz_table_z1, (2) at kurucz_table_z2
! (read only when kap_settings%use_two_z_tables). See kurucz_table_set.
           type(kurucz_table_set) :: kurucz(2)
           type(pure_z_table_set) :: pure_z
           type(alex94_table_set) :: alex94
           type(opal95_table_set) :: opal95
           type(alex06_table_set) :: alex06
! LAOL89 tables: (1) at laol_table_z1, (2) at laol_table_z2 (read
! only when kap_settings%use_two_z_tables). See laol_table_set.
           type(laol_table_set) :: laol(2)
      end type opacity_table_state

! ---- the opacity settings ----
! Everything kap needs from the run's configuration, in one record,
! installed by kap_init (from the star layer's kap_settings_from_star,
! setup/kap_settings_builder.f90) or kap_lib's set_kap_settings, through
! install_kap_settings; nothing else writes it.
      type, public :: kap_settings_t
! the interior tables (OPAL95-format, OPAL92, LAOL89) and the pure-Z
! table for helium-burning cores; the molecular tables (Alexander 2006
! or 1994/95, or Kurucz 1990); the Potekhin conductive opacities
         logical :: use_opal95_tables = .false., use_opal92_tables = .false., &
              use_laol89_tables = .false., use_pure_z_table = .false.
         logical :: use_alex06_tables = .false., use_alex95_tables = .false., &
              use_kurucz90_tables = .false.
         logical :: use_conductive_opacity = .false.
! read a second table at another Z and blend linearly in Z
         logical :: use_two_z_tables = .false.
! the tables' metal fractions: the OPAL95 surface table's, the two
! OPAL92, LAOL89 and Kurucz tables', the Alexander 1994/95 tables'
         real(dp) :: opal95_single_table_z = 0.0d0, opal_table_z1 = 0.0d0, &
              opal_table_z2 = 0.0d0, laol_table_z1 = 0.0d0, &
              laol_table_z2 = 0.0d0, kurucz_table_z1 = 0.0d0, &
              kurucz_table_z2 = 0.0d0, alex_table_z1 = 0.0d0
! log T below which only the molecular tables are used, and above which
! only the interior ones; the opacity ramps between them
         real(dp) :: molecular_opacity_logt_min = 0.0d0, &
              molecular_opacity_logt_max = 0.0d0
! how far Z may stray from the envelope's when only the surface table
! is loaded; the LAOL tables' extrapolation limit, whose natural log is the
! distance in log rho (laol_z_tolerance, legacy tollaol)
         real(dp) :: metal_fraction_match_tolerance = 0.0d0, &
              laol_extrapolation_tolerance = 0.0d0
      end type kap_settings_t
      type(kap_settings_t), public, protected :: kap_settings
      public :: install_kap_settings

! The tables kap_init reads, by path (blank: not given): the interior
! tables (OPAL 1995-format; OPAL 1992 and LAOL 1989, each with a second-Z
! table for use_two_z_tables; the pure-Z table), the molecular tables
! (Alexander 2006; the seven Alexander & Ferguson 1994 tables; Kurucz
! 1990 and its second Z) and the Potekhin conductive table. The star
! layer fills it from the controls (kap_table_paths_from_star,
! setup/kap_settings_builder.f90).
      type, public :: kap_table_paths_t
         character(len=256) :: opal95 = '', opal92 = '', opal92_z2 = '', &
              laol = '', laol_z2 = '', pure_z = ''
         character(len=256) :: alex06 = '', alex94(7) = '', kurucz = '', &
              kurucz_z2 = ''
         character(len=256) :: conductive = ''
      end type kap_table_paths_t

contains

! install_kap_settings: install the run's opacity settings (kap_lib's
! set_kap_settings, which also chooses their families).
subroutine install_kap_settings(settings)
      type(kap_settings_t), intent(in) :: settings
      kap_settings = settings
end subroutine install_kap_settings


! stencil4_locate
! Warm-start search for the first index idx of the 4-point stencil
! grid(idx:idx+3) around x, starting from the previous idx: walk down
! while x is below grid(idx+2), otherwise walk up, clamping to
! 1..n-3 (the Alexander lookups' T and log R searches); callers clamp idx
! to 1..n-2 first.
pure subroutine stencil4_locate(grid, n, x, idx)
      integer, intent(in) :: n
      real(dp), intent(in) :: grid(:), x
      integer, intent(inout) :: idx
      integer :: i

      if (x < grid(idx+2)) then
         do i = idx+1,2,-1
            if (x > grid(i)) then
               idx = i - 1
               exit
            end if
         end do
         if (i < (2)) then
         idx = 1
         end if
      else
         do i = idx+3,n
            if (x < grid(i)) then
               idx = i - 2
               idx = min(n-3,idx)
               exit
            end if
         end do
         if (i > n) then
         idx = n - 3
         end if
      end if
end subroutine stencil4_locate

! stencil4_locate_opal95
! The OPAL95 flavour of stencil4_locate (compute_opal95_opacity's T and
! log R searches): it walks up when x >= grid(idx+2), the upward scan
! stops at n-1 and there is no min() clamp on the way up. The two are
! kept apart so that each family keeps its F77 tie-breaking at nodes.
pure subroutine stencil4_locate_opal95(grid, n, x, idx)
      integer, intent(in) :: n
      real(dp), intent(in) :: grid(:), x
      integer, intent(inout) :: idx
      integer :: i

      if (x >= grid(idx+2)) then
         do i = idx+3,n-1
            if (x < grid(i)) then
               idx = i - 2
               exit
            end if
         end do
         if (i > (n-1)) then
         idx = n - 3
         end if
      else
         do i = idx+1,2,-1
            if (x > grid(i)) then
               idx = i - 1
               exit
            end if
         end do
         if (i < (2)) then
         idx = 1
         end if
      end if
end subroutine stencil4_locate_opal95


end module opacity_tables
