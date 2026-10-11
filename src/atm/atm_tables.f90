!----------------------------------------------------------------------
! atm_tables
!----------------------------------------------------------------------
! The settings record of atm/ and its tables (atm_table): the Kurucz and
! Castelli/Kurucz surface-pressure tables (each a surface_p_table), the
! Allard tables with their per-row bounds, and the result of the last
! table lookup.
module atm_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: atm_settings, atm_settings_t, atm_table, atm_table_ng, &
           atm_table_nga, atm_table_nt, atm_table_nta, atm_table_ntc, &
           set_atm_settings, surface_p_table
! The settings record of atm/: which tables atm_init loads (atm_choice
! 3 and 4 the Kurucz table, 4 also the Allard tables, 5 the
! Castelli/Kurucz table), whether the Allard lookups return P and T at
! tau = 100, and the [Fe/H] and [alpha/Fe] of the Allard table to read.
      type :: atm_settings_t
         integer :: atm_choice = 0
         logical :: allard_use_tau100 = .false.
         real(dp) :: allard_target_feh = 0.0d0, allard_target_alpha = 0.0d0
      end type atm_settings_t
      type(atm_settings_t), protected :: atm_settings
! Kurucz (nt) and Kurucz-Castelli (ntc) atmosphere tables: rows in
! log Teff, atm_table_ng columns in log g.
      integer, parameter :: atm_table_nt = 57, atm_table_ntc = 76
      integer, parameter :: atm_table_ng = 11
! Allard table maximum extents (rows in log Teff, columns in log g).
      integer, parameter :: atm_table_nta = 250, atm_table_nga = 25

! surface_p_table: one surface-pressure table -- log P at (log Teff row,
! log g column) -- with its log Teff rows, the atm_table_ng log g
! columns and each row's lowest and highest valid log g column. The
! Kurucz and Castelli/Kurucz tables have different row counts, so the
! arrays are allocated by the reader, zero-filled.
      type :: surface_p_table
           integer :: num_teff = 0
           real(dp), allocatable :: teff(:), logg(:)
           real(dp), allocatable :: log10_pressure(:,:)
           integer, allocatable :: gmin_index(:), gmax_index(:)
      end type surface_p_table

      type :: atm_table_state
! the Kurucz (atm_table_nt rows) and Castelli/Kurucz (atm_table_ntc
! rows) surface-pressure tables
           type(surface_p_table) :: kurucz_surface_p, castelli_surface_p
! the result of the last surface-table lookup: log P, and log T (log
! Teff, or T at tau = 100 from the Allard tables)
           real(dp) :: atm_log10_pressure, atm_log10_temperature
! the Allard surface-pressure table: the log Teff, log g, [Fe/H] and
! [alpha/Fe] grids, log P, and log P and log T at tau = 100
           real(dp) :: allard_teffl_grid(atm_table_nta), &
                allard_gl_grid(atm_table_nga), allard_feh_grid(atm_table_nga), &
                allard_alpha_grid(atm_table_nga)
           real(dp) :: allard_log10_pressure(atm_table_nta,atm_table_nga), &
                allard_log10_pressure_tau100(atm_table_nta,atm_table_nga), &
                allard_log10_temp_tau100(atm_table_nta,atm_table_nga)
           logical :: allard_is_old_nextgen
           integer :: allard_num_teff, allard_num_gl, allard_num_feh, &
                allard_num_alpha
! each log Teff row's valid log g range, and the table's ranges
           real(dp) :: allard_gl_row_min(atm_table_nta), &
                allard_gl_row_max(atm_table_nta)
           integer :: allard_gl_index_min(atm_table_nta), &
                allard_gl_index_max(atm_table_nta)
           real(dp) :: allard_teffl_min, allard_teffl_max, &
                allard_gl_min, allard_gl_max
! the log Teff range of the table read
           real(dp) :: allard_al_teffl_min, allard_al_teffl_max
      end type atm_table_state

      type(atm_table_state) :: atm_table

contains

subroutine set_atm_settings(settings)
      type(atm_settings_t), intent(in) :: settings
      atm_settings = settings
end subroutine set_atm_settings

!======================================================================
! Sources: the surface-pressure tables of Kurucz (1993) and Castelli &
! Kurucz (2003), and the model atmospheres of Allard & Hauschildt. This
! state was the COMMON blocks fac, atmprt, alatm01, alatm02, alatm05,
! atmos2 and atmos2c.

end module atm_tables
