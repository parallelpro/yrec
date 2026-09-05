!----------------------------------------------------------------------
! opacity_table_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles the low/high-temperature opacity
! lookup tables -- former common/gllot/, llot/, lintpl/, gllot2/,
! llot2/, lintpl2/ (OPAL92/Livermore), galot/, alot/, alotall/
! (Alexander 1994/95), galot06/, alot06/, alot06all/ (Alexander 2006),
! gkrz/, krz/, kipm/, intpl2/, gkrz2/, krz2/, kipm2/, intpl22/
! (Kurucz), slaol/, slaol2/, nwlaol2/, zlaol/, zslaol/ (LAOL89) -- into
! one derived type. These are Category 1 (loaded once at startup, read
! broadly thereafter), kept in their own module (like mhd_eos_lib and
! atm_table_lib) rather than folded into const_lib to keep large table
! arrays out of the globally-`use`d const_lib.
!
! DATA-statement defaults (kipmll/kipm/kipm2's cached-index scalars,
! and the Alexander 94/06 table-grid arrays) are preserved here as
! declaration-time initializers, transcribed verbatim from the
! original DATA statements in kap/alex94/read_alex94_tables.f90,
! kap/alex06/readalex06.f90, kap/kurucz90/kurucz.f90 (and its former
! kurucz2.f90 clone). kap/opal92/opal92_interp3d_z2.f90's mirror of
! kap/opal92/opal92_interp3d.f90's common/kipmll/ (its own common/kipmll2/) never had a
! DATA statement in the original, so its members are left without an
! initializer here too, preserving that asymmetry.
module opacity_table_lib
      implicit none
! OPAL92 (Livermore) table dimensions
      integer, parameter :: n_opal92_t = 50, n_opal92_d = 17, &
           n_opal92_x = 3, n_opal92_xt = 150, n_opal92_4d = 68
! Alexander 1994/95 table dimensions
      integer, parameter :: n_alex94_x = 7, n_alex94_z = 15, &
           n_alex94_t = 23, n_alex94_d = 17, n_alex94_xt = 8, &
           n_alex94_xz = 105
! Alexander 2006 table dimensions
      integer, parameter :: n_alex06_x = 9, n_alex06_z = 16, &
           n_alex06_t = 85, n_alex06_d = 19, n_alex06_xz = 143
! Kurucz table dimensions
      integer, parameter :: kurucz_max_num_temps = 60, &
           kurucz_max_num_densities = 50, kurucz_num_x_tables = 1, &
           kurucz_num_x_temp_entries = kurucz_max_num_temps*kurucz_num_x_tables, &
           kurucz_num_spline_coeffs = 4*kurucz_max_num_densities
! OPAL95 opacity table dimensions. n_opal95_xz is the SLOT count of
! the unpacked (13 Z rows x 10 X columns) layout used by ll95tbl and
! opal95_fixed_z_table (slot = start_index(iz) + ix), not the number
! of tables in the file (n_opal95_tables = 126). 2026 (bugsweep
! Batch 2): was 126 with a packed start_index -- see below.
      integer, parameter :: n_opal95_t = 70, n_opal95_d = 19, &
           n_opal95_x = 10, n_opal95_z = 13, n_opal95_xz = 130
! Number of tables in an OPAL95 table file (ll95tbl reads exactly this
! many; 8 X columns at all 13 Z, X=0.95 at 10 Z, X=1-Z at 12 Z).
      integer, parameter :: n_opal95_tables = 126
! Empty-cell markers. OPAL95 cells with no data are stored as
! opal95_missing_opacity and detected with ".ge./.gt. opal95_missing_test".
! OPAL92 tables use opal92_missing_opacity ("<= -9.999" tests).
      double precision, parameter :: opal95_missing_opacity = 9.999d0
      double precision, parameter :: opal95_missing_test = 9.9d0
      double precision, parameter :: opal92_missing_opacity = -9.999d0
! LAOL89 table dimensions (array extents; the file-size checks in
! rdlaol/rdzlaol accept at most 11 X columns).
      integer, parameter :: n_laol_x = 12, n_laol_rho = 104, n_laol_t = 52
! LAOL89 opacity cap returned by gtlaol/gtlaol2/gtpurz when the
! interpolated log10(opacity) overflows.
      double precision, parameter :: laol_opacity_cap = 1.0d35
! Composition-cache tolerances: a lookup whose X (and Z) is within this
! of the cached value reuses the previously built fixed-composition
! table. Each family keeps the tolerance its original code used.
      double precision, parameter :: alex_composition_tol = 1.0d-8
      double precision, parameter :: opal92_x_match_tol = 1.0d-5
      double precision, parameter :: opal95_composition_tol = 1.0d-4

! kurucz_table_set: one Kurucz90 molecular-opacity table (former
! common/gkrz/, krz/, kipm/, intpl2/ and their gkrz2/, krz2/, kipm2/,
! intpl22/ second-Z mirrors). 2026 wave 3 (R5): the two member sets
! kurucz_*/kurucz2_* became opacity_table%kurucz(1)/(2) so that one
! reader/spliner/interpolator serves both. Component names are the
! former kurucz_* names with the prefix dropped; shapes unchanged.
! check_range: kurucz.f90's original early-out (skip the table above
! log rho = -3 or log T = 4.1) was never in kurucz2.f; the merged
! kurucz() applies it only when this is .true. (read_kurucz_tables
! clears it for set 2).
      type, public :: kurucz_table_set
           double precision :: grid_logt(kurucz_max_num_temps)
           double precision :: log10_opacity(kurucz_num_x_temp_entries,kurucz_max_num_densities), &
                log10_rho(kurucz_num_x_temp_entries,kurucz_max_num_densities)
           integer :: num_temps
           integer :: ix_t = 1, ix_rho = 1
           double precision :: spline_coeffs(kurucz_num_x_temp_entries,kurucz_num_spline_coeffs)
           integer :: density_start_index(kurucz_num_x_temp_entries), &
                density_count(kurucz_num_x_temp_entries)
           logical :: check_range = .true.
      end type kurucz_table_set

      type, public :: opacity_table_state
! former common/gllot/, llot/, lintpl/ (OPAL92, first Z table)
           double precision :: opal92_grid_logt(n_opal92_t), &
                opal92_grid_x(n_opal92_x), opal92_grid_logr(n_opal92_d)
           double precision :: opal92_log10_opacity(n_opal92_xt,n_opal92_d)
           integer :: opal92_num_x, opal92_num_temps
           double precision :: opal92_spline_coeffs(n_opal92_xt,n_opal92_4d)
           integer :: opal92_density_start_index(n_opal92_xt), &
                opal92_density_count(n_opal92_xt)
! former common/gllot2/, llot2/, lintpl2/ (OPAL92, second Z table)
           double precision :: opal92_grid_logt_z2(n_opal92_t), &
                opal92_grid_x_z2(n_opal92_x), opal92_grid_logr_z2(n_opal92_d)
           double precision :: opal92_log10_opacity_z2(n_opal92_xt,n_opal92_d)
           integer :: opal92_num_x_z2, opal92_num_temps_z2
           double precision :: opal92_spline_coeffs_z2(n_opal92_xt,n_opal92_4d)
           integer :: opal92_density_start_index_z2(n_opal92_xt), &
                opal92_density_count_z2(n_opal92_xt)
! former common/llot4/ (OPAL92 surface table)
           double precision :: opal92_surface_x, opal92_surface_z
           double precision :: opal92_surface_spline_coeffs(n_opal92_t,n_opal92_4d)
           integer :: opal92_surface_x_index
! former common/llot42/ (OPAL92 surface table, second Z)
           double precision :: opal92_surface_x_z2, opal92_surface_z_z2
           double precision :: opal92_surface_spline_coeffs_z2(n_opal92_t,n_opal92_4d)
           integer :: opal92_surface_x_index_z2
! former common/kipmll/ (cached grid indices, DATA-initialized in
! kap/opal92/opal92_interp3d.f90) and common/kipmll2/ (its second-Z mirror in
! kap/opal92/opal92_interp3d_z2.f90, never DATA-initialized in the original -- no
! initializer here either, preserving that asymmetry).
           integer :: opal92_index_x = 1, opal92_index_t = 1, opal92_index_rho = 1
           integer :: opal92_index_x_z2, opal92_index_t_z2, opal92_index_rho_z2
! former common/galot/, alot/, alotall/ (Alexander 1994/95)
           double precision :: alex94_grid_logt(n_alex94_t) = &
                [3.00d0,3.05d0,3.10d0,3.15d0,3.20d0,3.25d0,3.30d0, &
                 3.35d0,3.40d0,3.45d0,3.50d0,3.55d0,3.60d0,3.65d0, &
                 3.70d0,3.75d0,3.80d0,3.85d0,3.90d0,3.95d0,4.00d0, &
                 4.05d0,4.10d0]
           double precision :: alex94_grid_x(n_alex94_x) = &
                [0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0]
           double precision :: alex94_grid_logr(n_alex94_d) = &
                [-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0,-4.5d0,-4.0d0, &
                 -3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0, &
                 0.0d0, 0.5d0, 1.0d0]
           double precision :: alex94_grid_z(n_alex94_z) = &
                [0.0d0, 0.00001d0, 0.00003d0, 0.0001d0, 0.0003d0, &
                 0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
                 0.04d0, 0.06d0, 0.08d0, 0.10d0]
           double precision :: alex94_opacity(n_alex94_xt,n_alex94_t,n_alex94_d)
           double precision :: alex94_cached_x = 0.0d0, alex94_cached_z = 0.0d0
           integer :: alex94_index_x = 4, alex94_index_t = 12, alex94_index_r = 9
           double precision :: alex94_full_opacity(n_alex94_xz,n_alex94_t,n_alex94_d)
! former common/galot06/, alot06/, alot06all/ (Alexander 2006)
           double precision :: alex06_grid_logt(n_alex06_t) = &
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
           double precision :: alex06_grid_x(n_alex06_x) = &
                [0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0,0.9d0,1.0d0]
           double precision :: alex06_grid_logr(n_alex06_d) = &
                [-8.0d0,-7.5d0,-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0, &
                 -4.5d0,-4.0d0,-3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0, &
                 -0.5d0, 0.0d0, 0.5d0, 1.0d0]
           double precision :: alex06_grid_z(n_alex06_z) = &
                [0.0d0, 0.00001d0, 0.00003d0, 0.0001d0, 0.0003d0, &
                 0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
                 0.04d0, 0.05d0,0.06d0, 0.08d0, 0.10d0]
           double precision :: alex06_opacity(n_alex06_t,n_alex06_d)
           double precision :: alex06_cached_x = 0.0d0, alex06_cached_z = 0.0d0
           integer :: alex06_index_x = 4, alex06_index_t = 43, alex06_index_r = 10
           double precision :: alex06_full_opacity(n_alex06_xz,n_alex06_t,n_alex06_d)
! Kurucz90 tables: (1) at kurucz_table_z1, (2) at kurucz_table_z2
! (read only when star%use_two_z_tables). See kurucz_table_set.
           type(kurucz_table_set) :: kurucz(2)
! former common/slaol/, slaol2/, nwlaol2/, zlaol/, zslaol/ (LAOL89/
! SLAOL: the spline-prepared first-Z, second-Z and pure-Z tables)
           double precision :: slaol_opacity(n_laol_x,n_laol_rho,n_laol_t), &
                slaol_log_rho(n_laol_x,n_laol_rho,n_laol_t), &
                slaol_d2opacity(n_laol_x,n_laol_rho,n_laol_t)
           integer :: slaol_num_points(n_laol_x,n_laol_t)
           double precision :: slaol2_opacity(n_laol_x,n_laol_rho,n_laol_t), &
                slaol2_log_rho(n_laol_x,n_laol_rho,n_laol_t), &
                slaol2_d2opacity(n_laol_x,n_laol_rho,n_laol_t)
           integer :: slaol2_num_points(n_laol_x,n_laol_t)
! former common/nwlaol/ and nwlaol2/: the LAOL89 tables as read by
! rdlaol (opacity, X grid, rho grid, T grid -- sulaol converts the T
! grid to log10 in place) and their extents.
           double precision :: laol_opacity(n_laol_x,n_laol_rho,n_laol_t), &
                laol_grid_x(n_laol_x), laol_grid_t(n_laol_t), laol_grid_rho(n_laol_rho)
           integer :: laol_num_x, laol_num_rho, laol_num_t
           double precision :: laol2_opacity(n_laol_x,n_laol_rho,n_laol_t), &
                laol2_grid_x(n_laol_x), laol2_grid_t(n_laol_t), laol2_grid_rho(n_laol_rho)
           integer :: laol2_num_x, laol2_num_rho, laol2_num_t
           double precision :: zlaol_opacity(n_laol_rho,n_laol_t), zlaol_logt_grid(n_laol_t), &
                zlaol_logrho_grid(n_laol_rho)
           integer :: zlaol_num_rho, zlaol_num_t
           double precision :: zslaol_opacity(n_laol_rho,n_laol_t), zslaol_log_rho(n_laol_rho,n_laol_t), &
                zslaol_d2opacity(n_laol_rho,n_laol_t)
           integer :: zslaol_num_points(n_laol_t)
! former common/llot95a/: the OPAL95 opacity table grid and full
! (all-Z) opacity array. opal95_grid_x/opal95_grid_z/
! opal95_table_start_index/opal95_num_x_at_z DATA-initialized in
! kap/opal95/ll95tbl.f90 (values transcribed verbatim below);
! opal95_grid_logt/opal95_full_opacity filled by file I/O there.
           double precision :: opal95_grid_logt(n_opal95_t)
           double precision :: opal95_grid_x(n_opal95_x) = &
                [0.0d0, 0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0,0.9d0, &
                 0.95d0,1.0d0]
           double precision :: opal95_grid_logr(n_opal95_d)
           double precision :: opal95_grid_z(n_opal95_z) = &
                [0.0d0, 0.0001d0, 0.0003d0, 0.001d0, 0.002d0, &
                 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
                 0.04d0, 0.06d0, 0.08d0, 0.10d0]
           double precision :: opal95_full_opacity(n_opal95_xz,n_opal95_t,n_opal95_d)
           integer :: opal95_num_x_at_z(n_opal95_z) = &
                [10,10,10,10,10,10,10,10,10,10,9,9,8]
! 2026 (bugsweep Batch 2): the inherited packed offsets
! [...,100,109,118] assumed rows 11-13 hold 9/9/8 tables, but the
! reader/consumers add the UNPACKED ix (1..10, with ix = 9/10 absent
! at high Z), so the X = 1-Z table of Z = 0.06 (slot 110) landed on
! (X = 0, Z = 0.08) and that of Z = 0.08 (slot 119) on (X = 0,
! Z = 0.10), overwriting them. Unpacked offsets of 10 per Z row
! remove the collision; slots 109, 119, 129, 130 stay unwritten and
! are never addressed by the stencils.
           integer :: opal95_table_start_index(n_opal95_z) = &
                [0,10,20,30,40,50,60,70,80,90,100,110,120]
! former common/llot95/: the single-Z OPAL95 opacity table, sliced at
! the model's actual Z. atm/turnover/acoustic_depths.f90 declared a mismatched single-
! scalar layout for this block (never read/set there) -- the majority
! (7-file) array+scalar layout below is used; acoustic_depths.f90's reference
! now maps onto (an unused corner of) this same type.
           double precision :: opal95_fixed_z_opacity(n_opal95_x,n_opal95_t,n_opal95_d)
           double precision :: opal95_fixed_z
! former common/llot95e/: the OPAL95 surface-X opacity table.
           double precision :: opal95_surface_opacity(n_opal95_t,n_opal95_d)
           double precision :: opal95_surface_x
! former common/op95indx/: cached Z/X/T/rho interpolation indices, all
! 22 scalar values DATA-initialized to 1 in kap/opal95/ll95tbl.f90.
           integer :: opal95_index_z = 1
           integer :: opal95_index_x(4,4) = 1
           integer :: opal95_index_t = 1
           integer :: opal95_index_rho(4) = 1
! former common/op95fact/: Z/X/T/rho interpolation weights and
! T/rho derivative weights, recomputed fresh each call -- no DATA.
           double precision :: opal95_weight_z(4), opal95_weight_x(4,4), &
                opal95_weight_t(4), opal95_dweight_t(4), &
                opal95_weight_rho(4,4), opal95_dweight_rho(4,4)
! former common/op95ext/: log(R) extrapolation-edge state, set fresh
! each call -- no DATA.
           double precision :: opal95_logr, opal95_logr_lo_edge, &
                opal95_logr_hi_edge(4)
           logical :: opal95_extrap_lo, opal95_extrap_hi, opal95_extrap_hi_row(4)
      end type opacity_table_state

      type(opacity_table_state), public :: opacity_table


! 2026 (phase six, step 3 -- ROADMAP.md): evicted here from
! const_lib, where this table/working data had landed during the
! phase-one COMMON conversion; it belongs with this domain's state.
! former common/nwlaol/'s one remaining member (the table arrays are
! opacity_table%laol_*; the namelist controls tollaol/use_pure_z_table
! are star%ctrl members since 2026 wave 2; the table unit numbers are
! newunit locals of kap/laol89/rdlaol.f90/rdzlaol.f90). llaol is not a
! namelist value and nothing in src/ assigns it: it stays at its
! compile-time .false., so core/read_starting_model.f90's
! `if (.not.llaol)` guard (use the LAOL-table metal mix instead of
! mixture_weights_seed) is never skipped.
      logical :: llaol = .false.


contains

! stencil4_locate
! Warm-start search for the first index idx of the 4-point stencil
! grid(idx:idx+3) around x, starting from the previous idx: walk down
! while x is below grid(idx+2), otherwise walk up, clamping to
! 1..n-3. 2026 readability: extracted token-for-token from the four
! identical copies in alex94_interp3d (T, R) and getalex06 (T, R);
! callers must pre-clamp idx to 1..n-2 as they did before.
subroutine stencil4_locate(grid, n, x, idx)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: grid(n), x
      integer, intent(inout) :: idx
      integer :: i

      if (x.lt.grid(idx+2)) then
         do i = idx+1,2,-1
            if (x.gt.grid(i)) then
               idx = i - 1
               exit
            endif
         end do
         if (i < (2)) then
         idx = 1
         end if
      else
         do i = idx+3,n
            if (x.lt.grid(i)) then
               idx = i - 2
               idx = min(n-3,idx)
               exit
            endif
         end do
         if (i > n) then
         idx = n - 3
         end if
      endif
end subroutine stencil4_locate

! stencil4_locate_opal95
! The OPAL95 flavour of stencil4_locate (getopal95's T and rho
! searches): the up/down decision is ".ge. grid(idx+2)", the upward
! scan stops at n-1 and there is no min() clamp on the way up.
! Extracted token-for-token; the two flavours are deliberately not
! merged so that each family keeps its own tie-breaking.
subroutine stencil4_locate_opal95(grid, n, x, idx)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: grid(n), x
      integer, intent(inout) :: idx
      integer :: i

      if (x.ge.grid(idx+2)) then
         do i = idx+3,n-1
            if (x.lt.grid(i)) then
               idx = i - 2
               exit
            endif
         end do
         if (i > (n-1)) then
         idx = n - 3
         end if
      else
         do i = idx+1,2,-1
            if (x.gt.grid(i)) then
               idx = i - 1
               exit
            endif
         end do
         if (i < (2)) then
         idx = 1
         end if
      endif
end subroutine stencil4_locate_opal95

! alex_x_stencil_start
! First of the four Alexander-1994 X tables bracketing hydrogen_fraction
! (the 7-column X grid is walked as a fixed ladder: 1, 2, 3 or 4).
! Extracted token-for-token from alex94_interp3d and alex94_surface_table.
integer function alex_x_stencil_start(hydrogen_fraction)
      implicit none
      double precision, intent(in) :: hydrogen_fraction

      if (hydrogen_fraction.lt.opacity_table%alex94_grid_x(4)) then
         if (hydrogen_fraction.gt.opacity_table%alex94_grid_x(3)) then
            alex_x_stencil_start = 2
         else
            alex_x_stencil_start = 1
         endif
      else
         if (hydrogen_fraction.gt.opacity_table%alex94_grid_x(5)) then
            alex_x_stencil_start = 4
         else
            alex_x_stencil_start = 3
         endif
      endif
end function alex_x_stencil_start

end module opacity_table_lib
