!----------------------------------------------------------------------
! atm_lib
!----------------------------------------------------------------------
! The public interface of atm/, the surface boundary condition:
! everything code outside atm/ may use, here and only here
! (check_boundaries enforces it).
!     atm_init             check and install the settings (atm_settings_t,
!                          built by the star layer's atm_settings_from_star)
!                          and read the table they select
!     atm_get_surface_pt   log P and log T at (log Teff, log g) from the
!                          Kurucz, Allard or Castelli/Kurucz table
!     atm_allard_teff_max  the Allard tables' highest log Teff
!     atm_eddington_gray ... atm_castelli_kurucz, atm_choice_is_known,
!     atm_choice_is_tabulated, atm_choice_name
!                          what each atm_choice is (atm_choices)
!     atm_leaves_table, atm_reenters_table, atm_note_choice_switch
!                          when the star layer replaces a table by the
!                          gray relation and returns to it, and the note
!     ttau_log10_temperature, ttau_start_log10_temperature,
!     ttau_photosphere_x_limit, hsra_t_tau_offset
!                          the analytic T-tau relations (t_tau_relations)
! It keeps the run's settings and tables, and it reports what the table
! lookups return (atm_messages); the lookups and the T-tau relations are
! star-blind kernels. The atmosphere integration that uses the T-tau
! relations, and its results, belong to the star layer
! (core/envelope/envelope_integration_kernel.f90). atm reads no star.
module atm_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use atm_tables, only: atm_settings_t, atm_table_paths_t, surface_pressure_table, &
           allard_table, kurucz_num_teff, castelli_num_teff
      use atm_choices, only: atm_eddington_gray, atm_krishna_swamy, atm_hsra, atm_kurucz, &
           atm_allard, atm_castelli_kurucz, atm_choice_is_known, atm_choice_is_tabulated, &
           atm_choice_name, has_gray_switch, kurucz_gray_switch_log_teff
      use kurucz_table_reader, only: read_surface_pressure_table
      use allard_table_reader, only: read_allard_tables
      use kurucz_atmosphere, only: interpolate_surface_pressure
      use allard_atmosphere, only: interpolate_allard_surface, allard_lookup_result
      use t_tau_relations, only: ttau_log10_temperature, &
           ttau_start_log10_temperature, ttau_photosphere_x_limit, hsra_t_tau_offset
      use atm_messages, only: report_atm_error, report_surface_pressure_outside, &
           log_spline_failure, report_allard_teff_above, report_allard_gravity_outside, &
           report_allard_teff_below, report_allard_gravity_outside_row, &
           report_polint_failure, log_allard_lookup, note_choice_switch
      implicit none
      private
      public :: atm_init, atm_get_surface_pt, atm_allard_teff_max
      public :: atm_settings_t, atm_table_paths_t
      public :: atm_eddington_gray, atm_krishna_swamy, atm_hsra, atm_kurucz, atm_allard, &
           atm_castelli_kurucz, atm_choice_is_known, atm_choice_is_tabulated, atm_choice_name
      public :: atm_leaves_table, atm_reenters_table, atm_note_choice_switch
      public :: ttau_log10_temperature, ttau_start_log10_temperature, &
           ttau_photosphere_x_limit, hsra_t_tau_offset

! The run's settings and tables: atm_init installs them, and nothing else
! writes them. A table atm_init did not read keeps its type's defaults
! (not loaded; the Allard tables' limits 0).
      type(atm_settings_t) :: installed_settings
      type(surface_pressure_table) :: kurucz_table, castelli_table
      type(allard_table) :: allard_tables

contains

!----------------------------------------------------------------------
! atm_init
!----------------------------------------------------------------------
! Checks the settings, installs them, and reads the table they select:
! the Kurucz table (atm_choice 3), the Allard tables (4) or the
! Castelli/Kurucz table (5). An atm_choice outside 0-5, or a table that
! cannot be read, returns ierr = 1 after its message.
subroutine atm_init(settings, paths, ierr)
      type(atm_settings_t), intent(in) :: settings
      type(atm_table_paths_t), intent(in) :: paths
      integer, intent(out) :: ierr

      ierr = 0
      if (.not. atm_choice_is_known(settings%atm_choice)) then
         call report_unknown_choice(settings%atm_choice)
         ierr = 1
         return
      end if
      installed_settings = settings

      select case (settings%atm_choice)
      case (atm_kurucz)
         call read_surface_pressure_table(paths%surface_pressure, kurucz_num_teff, &
              kurucz_table, ierr)
      case (atm_allard)
         call read_allard_tables(paths%allard, settings%allard_target_feh, &
              settings%allard_target_alpha, allard_tables, ierr)
      case (atm_castelli_kurucz)
         call read_surface_pressure_table(paths%surface_pressure, castelli_num_teff, &
              castelli_table, ierr)
      end select
end subroutine atm_init

! atm_choice is none of the six: an error (M236: it used to load no table
! and run the HSRA relation).
subroutine report_unknown_choice(atm_choice)
      integer, intent(in) :: atm_choice
      character(len=200) :: text
      write(text,'(a,i0,a)') 'atm_choice = ', atm_choice, ' is not an atmosphere choice' // &
           ' (0 Eddington gray, 1 Krishna-Swamy, 2 HSRA, 3 Kurucz, 4 Allard, 5 Castelli/Kurucz)'
      call report_atm_error(trim(text))
end subroutine report_unknown_choice

!----------------------------------------------------------------------
! atm_get_surface_pt
!----------------------------------------------------------------------
! The tabulated surface boundary condition at (log Teff, log g) for
! atm_choice 3 (Kurucz), 4 (Allard) or 5 (Castelli/Kurucz): log P, and
! log T (log Teff, or T at tau = 100 from the Allard tables with
! allard_use_tau100). Both are set only when the lookup succeeds.
! lookup_failed: the point is outside the Allard tables (the caller
! falls back to a gray atmosphere), as is every point when atm_init did
! not read them (their limits are 0; O19). ierr /= 0: the point is
! outside a Kurucz table or below the Allard tables (1), an interpolation
! failed (its numerics_status), the Kurucz or Castelli/Kurucz table was
! not read, or atm_choice is not a table (1). The messages go through
! atm_messages. atm_choice is the caller's current choice, which the
! star layer switches between tabulated and gray as Teff changes.
subroutine atm_get_surface_pt(atm_choice, log_teff, log_g, log10_pressure, &
     log10_temperature, lookup_failed, ierr)
      integer, intent(in) :: atm_choice
      real(dp), intent(in) :: log_teff, log_g
      real(dp), intent(inout) :: log10_pressure, log10_temperature
      logical, intent(out) :: lookup_failed
      integer, intent(out) :: ierr

      ierr = 0
      lookup_failed = .false.
      select case (atm_choice)
      case (atm_kurucz)
         call look_up_surface_pressure(kurucz_table)
      case (atm_castelli_kurucz)
         call look_up_surface_pressure(castelli_table)
      case (atm_allard)
         call look_up_allard()
      case default
         call report_not_a_table()
         ierr = 1
      end select

contains

! a Kurucz or Castelli/Kurucz table: log P, with log T = log Teff
subroutine look_up_surface_pressure(tbl)
      type(surface_pressure_table), intent(in) :: tbl
      integer :: spline_status

      if (.not. tbl%loaded) then
         call report_not_loaded()
         ierr = 1
         return
      end if
      call interpolate_surface_pressure(tbl, log_teff, log_g, log10_pressure, ierr, &
           spline_status)
      if (spline_status /= 0) then
         call log_spline_failure(spline_status)
      else if (ierr /= 0) then
         call report_surface_pressure_outside(log_teff, log_g)
      end if
      if (ierr /= 0) return
      log10_temperature = log_teff
end subroutine look_up_surface_pressure

! the Allard tables: each reason the point is outside them, in the order
! the lookup tests them; then the result, noted in the run log (M235)
subroutine look_up_allard()
      type(allard_lookup_result) :: lookup
      integer :: row

      call interpolate_allard_surface(allard_tables, installed_settings%allard_use_tau100, &
           log_teff, log_g, lookup)
      if (lookup%teff_above) call report_allard_teff_above(allard_tables%log_teff_min, &
           log_teff, allard_tables%log_teff_max)
      if (lookup%gravity_outside) call report_allard_gravity_outside(allard_tables%log_g_min, &
           log_g, allard_tables%log_g_max)
      if (lookup%teff_below) call report_allard_teff_below(allard_tables%log_teff_min, &
           log_teff, allard_tables%log_teff_max)
      row = lookup%row_outside
      if (row > 0) call report_allard_gravity_outside_row(allard_tables%row_log_g_min(row), &
           log_g, allard_tables%row_log_g_max(row))
      if (lookup%polint_status /= 0) call report_polint_failure(lookup%polint_status)

      lookup_failed = lookup%teff_above .or. lookup%gravity_outside .or. row > 0
      if (lookup%teff_below) ierr = 1
      if (lookup%polint_status /= 0) ierr = lookup%polint_status
      if (lookup_failed .or. ierr /= 0) return
      call log_allard_lookup(installed_settings%allard_use_tau100, log_teff, log_g, &
           lookup%log10_temperature, lookup%log10_pressure)
      log10_pressure = lookup%log10_pressure
      log10_temperature = lookup%log10_temperature
end subroutine look_up_allard

subroutine report_not_loaded()
      character(len=200) :: text
      write(text,'(a,i0,a)') 'atm_get_surface_pt: the table of atm_choice ', atm_choice, &
           ' was not read (atm_init reads only the deck''s choice''s table)'
      call report_atm_error(trim(text))
end subroutine report_not_loaded

subroutine report_not_a_table()
      character(len=200) :: text
      write(text,'(a,i0,a)') 'atm_get_surface_pt: atm_choice ', atm_choice, &
           ' is not a table (3 Kurucz, 4 Allard, 5 Castelli/Kurucz)'
      call report_atm_error(trim(text))
end subroutine report_not_a_table

end subroutine atm_get_surface_pt

!----------------------------------------------------------------------
! atm_allard_teff_max
!----------------------------------------------------------------------
! The highest log Teff of the Allard tables read (the star layer
! switches to a gray atmosphere above it).
pure function atm_allard_teff_max() result(log_teff_max)
      real(dp) :: log_teff_max
      log_teff_max = allard_tables%log_teff_last
end function atm_allard_teff_max

!----------------------------------------------------------------------
! atm_leaves_table, atm_reenters_table, atm_note_choice_switch
!----------------------------------------------------------------------
! The star layer's switch of its working choice (star%atm_choice_in_use,
! core/structure/surface_boundary.f90): a table whose choice has a gray
! switch (has_gray_switch: Kurucz, Allard) hands the surface over to the
! Eddington gray relation at and above its switch's log Teff, and takes
! it back below: kurucz_gray_switch_log_teff for the Kurucz table, the
! last log Teff read for the Allard tables. Neither holds for any other
! choice.
pure logical function atm_leaves_table(atm_choice, log10_teff)
      integer, intent(in) :: atm_choice
      real(dp), intent(in) :: log10_teff
      atm_leaves_table = .false.
      if (has_gray_switch(atm_choice)) atm_leaves_table = log10_teff >= gray_switch_log_teff(atm_choice)
end function atm_leaves_table

pure logical function atm_reenters_table(atm_choice, log10_teff)
      integer, intent(in) :: atm_choice
      real(dp), intent(in) :: log10_teff
      atm_reenters_table = .false.
      if (has_gray_switch(atm_choice)) atm_reenters_table = log10_teff < gray_switch_log_teff(atm_choice)
end function atm_reenters_table

! the log Teff of a choice's gray switch
pure real(dp) function gray_switch_log_teff(atm_choice)
      integer, intent(in) :: atm_choice
      if (atm_choice == atm_allard) then
         gray_switch_log_teff = allard_tables%log_teff_last
      else
         gray_switch_log_teff = kurucz_gray_switch_log_teff
      end if
end function gray_switch_log_teff

! the terminal's note of a switch from from_choice to to_choice at
! log10_teff (atm_messages' note_choice_switch)
subroutine atm_note_choice_switch(from_choice, to_choice, log10_teff)
      integer, intent(in) :: from_choice, to_choice
      real(dp), intent(in) :: log10_teff
      call note_choice_switch(from_choice, to_choice, log10_teff, allard_tables%log_teff_last)
end subroutine atm_note_choice_switch

!======================================================================
! Sources: the F77's SETUPS (the table loads, by KTTAU), ENVINT (the
! lookup by KTTAU, the gray fallback of ALSURFP's LAlFail) and SURFBC
! (the switch to gray and back, MHP 9/01); the 2026 port
! gathered them behind this facade.

end module atm_lib
