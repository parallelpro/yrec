!----------------------------------------------------------------------
! atm_tables
!----------------------------------------------------------------------
! The types of atm/: its settings record, the paths of the tables it
! reads, and the tables themselves, with their sizes. atm_lib keeps the
! run's instances (atm_init fills them); a reader fills the table it is
! given, and an interpolator reads the table it is given.
!   atm_settings_t          the settings (setup/atm_settings_builder.f90)
!   atm_table_paths_t       the table paths (the same builder)
!   surface_pressure_table  a Kurucz or Castelli/Kurucz table: log P at
!                           (log Teff, log g)
!   allard_table            the Allard tables: log P, and log P and log T
!                           at tau = 100, at (log Teff, log g)
module atm_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: atm_settings_t, atm_table_paths_t, surface_pressure_table, &
           allard_table
      public :: kurucz_num_teff, castelli_num_teff, surface_pressure_num_g, &
           allard_max_teff, allard_max_g

! The settings record of atm/: which tables atm_init loads (atm_choice
! 3 the Kurucz table, 4 the Allard tables, 5 the Castelli/Kurucz
! table), whether the Allard lookups return P and T at tau = 100, and
! the [Fe/H] and [alpha/Fe] of the Allard table to read.
      type :: atm_settings_t
         integer :: atm_choice = 0
         logical :: allard_use_tau100 = .false.
         real(dp) :: allard_target_feh = 0.0d0, allard_target_alpha = 0.0d0
      end type atm_settings_t

! The tables atm_init reads, by path: the Kurucz or Castelli/Kurucz
! surface-pressure table (atm_choice 3, or 5) and the Allard tables
! (4). The star layer fills it from the controls
! (atm_table_paths_from_star, setup/atm_settings_builder.f90).
      type :: atm_table_paths_t
         character(len=256) :: surface_pressure = '', allard = ''
      end type atm_table_paths_t

! The surface-pressure tables' rows in log Teff (Kurucz 1993; Castelli &
! Kurucz 2003) and their columns in log g (both).
      integer, parameter :: kurucz_num_teff = 57, castelli_num_teff = 76
      integer, parameter :: surface_pressure_num_g = 11

! A surface-pressure table: log P at (log Teff row, log g column), with
! each row's first and last valid log g column (a row may hold -999 at
! either end). The reader allocates the arrays (the two tables have
! different row counts), zero-filled.
      type :: surface_pressure_table
         logical :: loaded = .false.
         integer :: num_teff = 0
         real(dp), allocatable :: log_teff(:), log_g(:)
         real(dp), allocatable :: log10_pressure(:,:)
         integer, allocatable :: first_column(:), last_column(:)
      end type surface_pressure_table

! The Allard tables' largest extents (rows in log Teff, columns in log g).
      integer, parameter :: allard_max_teff = 250, allard_max_g = 25

! The Allard tables at one [Fe/H] and [alpha/Fe]: the log Teff and log g
! grids, log P at T = Teff, and log P and log T at tau = 100; -999 marks
! an entry the file does not give. Each row's first and last valid log g
! column and the log g range a lookup may use (the valid range extended
! below and above), and the table's ranges. Everything is 0 until the
! reader runs (atm_lib refuses a lookup in tables it did not read).
      type :: allard_table
         logical :: loaded = .false.
! the old NextGen 1 layout (one composition, no tau = 100 columns)
         logical :: is_old_nextgen = .false.
         integer :: num_teff = 0, num_g = 0, num_feh = 0, num_alpha = 0
         real(dp) :: log_teff(allard_max_teff) = 0.0d0
         real(dp) :: log_g(allard_max_g) = 0.0d0
! the [Fe/H] and [alpha/Fe] values counted (sized as the log g grid,
! which bounds their count; the F77's)
         real(dp) :: feh(allard_max_g) = 0.0d0, alpha(allard_max_g) = 0.0d0
         real(dp) :: log10_pressure(allard_max_teff,allard_max_g) = 0.0d0
         real(dp) :: log10_pressure_tau100(allard_max_teff,allard_max_g) = 0.0d0
         real(dp) :: log10_temperature_tau100(allard_max_teff,allard_max_g) = 0.0d0
         integer :: first_column(allard_max_teff) = 0, last_column(allard_max_teff) = 0
         real(dp) :: row_log_g_min(allard_max_teff) = 0.0d0, &
              row_log_g_max(allard_max_teff) = 0.0d0
! the range a lookup may use: one row below the first and above the last
! log Teff, and the widest of the rows' log g ranges
         real(dp) :: log_teff_min = 0.0d0, log_teff_max = 0.0d0
         real(dp) :: log_g_min = 0.0d0, log_g_max = 0.0d0
! the first and last tabulated log Teff
         real(dp) :: log_teff_first = 0.0d0, log_teff_last = 0.0d0
      end type allard_table

!======================================================================
! Sources: the surface-pressure tables of Kurucz (1993) and Castelli &
! Kurucz (2003), and the model atmospheres of Allard & Hauschildt. These
! types replace the COMMON blocks fac, atmprt, alatm01, alatm02, alatm05,
! atmos2 and atmos2c.

end module atm_tables
