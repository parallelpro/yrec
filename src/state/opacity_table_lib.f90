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
! original DATA statements in kap/alxtbl.f90, kap/readalex06.f90,
! kap/kurucz.f90, kap/kurucz2.f90. kap/yllo3d2.f90's mirror of
! kap/yllo3d.f90's common/kipmll/ (its own common/kipmll2/) never had a
! DATA statement in the original, so its members are left without an
! initializer here too, preserving that asymmetry.
module opacity_table_lib
      implicit none
! OPAL92 (Livermore) table dimensions
      integer, parameter :: n_opal92_t = 50, n_opal92_d = 17, &
           n_opal92_x = 3, n_opal92_xt = 150, n_opal92_4d = 68
! Alexander 1994/95 table dimensions
      integer, parameter :: n_alex95_x = 7, n_alex95_z = 15, &
           n_alex95_t = 23, n_alex95_d = 17, n_alex95_xt = 8, &
           n_alex95_xz = 105
! Alexander 2006 table dimensions
      integer, parameter :: n_alex06_x = 9, n_alex06_z = 16, &
           n_alex06_t = 85, n_alex06_d = 19, n_alex06_xz = 143
! Kurucz table dimensions
      integer, parameter :: kurucz_max_num_temps = 60, &
           kurucz_max_num_densities = 50, kurucz_num_x_tables = 1, &
           kurucz_num_x_temp_entries = kurucz_max_num_temps*kurucz_num_x_tables, &
           kurucz_num_spline_coeffs = 4*kurucz_max_num_densities

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
! kap/yllo3d.f90) and common/kipmll2/ (its second-Z mirror in
! kap/yllo3d2.f90, never DATA-initialized in the original -- no
! initializer here either, preserving that asymmetry).
           integer :: abund_index = 1, temp_index = 1, dens_index = 1
           integer :: abund_index_z2, temp_index_z2, dens_index_z2
! former common/galot/, alot/, alotall/ (Alexander 1994/95)
           double precision :: alex95_grid_logt(n_alex95_t) = &
                [3.00d0,3.05d0,3.10d0,3.15d0,3.20d0,3.25d0,3.30d0, &
                 3.35d0,3.40d0,3.45d0,3.50d0,3.55d0,3.60d0,3.65d0, &
                 3.70d0,3.75d0,3.80d0,3.85d0,3.90d0,3.95d0,4.00d0, &
                 4.05d0,4.10d0]
           double precision :: alex95_grid_x(n_alex95_x) = &
                [0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.7d0,0.8d0]
           double precision :: alex95_grid_logr(n_alex95_d) = &
                [-7.0d0,-6.5d0,-6.0d0,-5.5d0,-5.0d0,-4.5d0,-4.0d0, &
                 -3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0, &
                 0.0d0, 0.5d0, 1.0d0]
           double precision :: alex95_grid_z(n_alex95_z) = &
                [0.0d0, 0.00001d0, 0.00003d0, 0.0001d0, 0.0003d0, &
                 0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
                 0.04d0, 0.06d0, 0.08d0, 0.10d0]
           double precision :: alex95_opacity(n_alex95_xt,n_alex95_t,n_alex95_d)
           double precision :: alex95_cached_x = 0.0d0, alex95_cached_z = 0.0d0
           integer :: alex95_index_x = 4, alex95_index_t = 12, alex95_index_r = 9
           double precision :: alex95_full_opacity(n_alex95_xz,n_alex95_t,n_alex95_d)
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
! former common/gkrz/, krz/, kipm/, intpl2/ (Kurucz, first table)
           double precision :: kurucz_grid_logt(kurucz_max_num_temps)
           double precision :: kurucz_log10_opacity(kurucz_num_x_temp_entries,kurucz_max_num_densities), &
                kurucz_log10_rho(kurucz_num_x_temp_entries,kurucz_max_num_densities)
           integer :: kurucz_num_temps
           integer :: kurucz_ix_x = 1, kurucz_ix_t = 1, kurucz_ix_rho = 1
           double precision :: kurucz_spline_coeffs(kurucz_num_x_temp_entries,kurucz_num_spline_coeffs)
           integer :: kurucz_density_start_index(kurucz_num_x_temp_entries), &
                kurucz_density_count(kurucz_num_x_temp_entries)
! former common/gkrz2/, krz2/, kipm2/, intpl22/ (Kurucz, second table)
           double precision :: kurucz2_grid_logt(kurucz_max_num_temps)
           double precision :: kurucz2_log10_opacity(kurucz_num_x_temp_entries,kurucz_max_num_densities), &
                kurucz2_log10_rho(kurucz_num_x_temp_entries,kurucz_max_num_densities)
           integer :: kurucz2_num_temps
           integer :: kurucz2_ix_x = 1, kurucz2_ix_t = 1, kurucz2_ix_rho = 1
           double precision :: kurucz2_spline_coeffs(kurucz_num_x_temp_entries,kurucz_num_spline_coeffs)
           integer :: kurucz2_density_start_index(kurucz_num_x_temp_entries), &
                kurucz2_density_count(kurucz_num_x_temp_entries)
! former common/slaol/, slaol2/, nwlaol2/, zlaol/, zslaol/ (LAOL89/
! SLAOL, dimensions are raw literals in every declaring file, not
! named PARAMETERs -- preserved as literals here too)
           double precision :: slaol_opacity(12,104,52), slaol_log_rho(12,104,52), &
                slaol_d2opacity(12,104,52)
           integer :: slaol_num_points(12,52)
           double precision :: slaol2_opacity(12,104,52), slaol2_log_rho(12,104,52), &
                slaol2_d2opacity(12,104,52)
           integer :: slaol2_num_points(12,52)
           double precision :: olaol2(12,104,52), oxa2(12), ot2(52), orho2(104)
           integer :: nxyz2, nrho2, nt2
           double precision :: zlaol_opacity(104,52), zlaol_logt_grid(52), &
                zlaol_logrho_grid(104)
           integer :: zlaol_num_rho, zlaol_num_t
           double precision :: zslaol_opacity(104,52), zslaol_log_rho(104,52), &
                zslaol_d2opacity(104,52)
           integer :: zslaol_num_points(52)
      end type opacity_table_state

      type(opacity_table_state), public :: opacity_table

end module opacity_table_lib
