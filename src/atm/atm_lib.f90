!----------------------------------------------------------------------
! atm_lib
!----------------------------------------------------------------------
! The public interface of atm/, the surface boundary condition:
! everything code outside atm/ may use, here and only here
! (check_boundaries enforces it).
!     atm_init             install the settings (atm_settings_t, built
!                          by the star layer's atm_settings_from_star) and
!                          read the tables they select
!     atm_get_surface_pt   log P and log T at (log Teff, log g) from the
!                          Kurucz, Allard or Castelli/Kurucz table
!     atm_allard_teff_max  the Allard tables' highest log Teff
!     set_atm_settings     install settings without reading tables
!     ttau_log10_temperature, ttau_start_log10_temperature,
!     ttau_photosphere_x_limit, hsra_t_tau_offset
!                          the analytic T-tau relations (t_tau_relations)
! The atmosphere integration that uses the T-tau relations, and its
! results, belong to the star layer (core/envelope_integration_kernel,
! core/atmosphere_derivs). atm reads no star.
module atm_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use kurucz_atmosphere, only: read_surface_pressure_table, interpolate_surface_pressure
      use allard_atmosphere, only: read_allard_tables, interpolate_allard_surface
      use t_tau_relations, only: ttau_log10_temperature, &
           ttau_start_log10_temperature, ttau_photosphere_x_limit, hsra_t_tau_offset
      use atm_tables, only: atm_settings_t, atm_table, atm_table_nt, &
           atm_table_ntc, set_atm_settings, atm_table_paths_t
      use run_messages, only: report_error
      implicit none
      private
      public :: atm_init, atm_get_surface_pt, atm_allard_teff_max
      public :: atm_settings_t, set_atm_settings, atm_table_paths_t
      public :: ttau_log10_temperature, ttau_start_log10_temperature, &
           ttau_photosphere_x_limit, hsra_t_tau_offset
contains

!----------------------------------------------------------------------
! atm_init
!----------------------------------------------------------------------
! Installs the settings and reads the tables they select: the Kurucz
! table (atm_choice 3), the Allard tables (4) or the Castelli/Kurucz
! table (5). An atm_choice outside 0-5, or a table that
! cannot be read, returns ierr = 1 with the diagnostic printed.
subroutine atm_init(settings, paths, ierr)
      type(atm_settings_t), intent(in) :: settings
      type(atm_table_paths_t), intent(in) :: paths
      integer, intent(out) :: ierr

      ierr = 0
      if (settings%atm_choice < 0 .or. settings%atm_choice > 5) then
         call report_unknown_choice(settings%atm_choice)
         ierr = 1
         return
      end if
      call set_atm_settings(settings)

      if (settings%atm_choice == 3) then
         call read_surface_pressure_table(paths%surface_pressure, atm_table%kurucz_surface_p, &
              atm_table_nt, ierr)
         if (ierr /= 0) return
! MHP 6/97: the Allard model atmospheres (M239: the Kurucz table, which
! no choice-4 lookup reads, is no longer loaded with them)
      else if (settings%atm_choice == 4) then
         call read_allard_tables(paths%allard, ierr)
         if (ierr /= 0) return
! JNT 6/14: the Kurucz/Castelli atmospheres
      else if (settings%atm_choice == 5) then
         call read_surface_pressure_table(paths%surface_pressure, atm_table%castelli_surface_p, &
              atm_table_ntc, ierr)
         if (ierr /= 0) return
      end if
end subroutine atm_init

! atm_choice is none of the six: an error, on the terminal and in the
! run log (M236: it used to load no table and run the HSRA relation).
subroutine report_unknown_choice(atm_choice)
      integer, intent(in) :: atm_choice
      character(len=200) :: text
      write(text,'(a,i0,a)') 'atm_choice = ', atm_choice, ' is not an atmosphere choice' // &
           ' (0 Eddington gray, 1 Krishna-Swamy, 2 HSRA, 3 Kurucz, 4 Allard, 5 Castelli/Kurucz)'
      call report_error(trim(text))
end subroutine report_unknown_choice

!----------------------------------------------------------------------
! atm_get_surface_pt
!----------------------------------------------------------------------
! The tabulated surface boundary condition at (log Teff, log g) for
! atm_choice 3 (Kurucz), 4 (Allard) or 5 (Castelli/Kurucz): log P, and
! log T (log Teff, or T at tau = 100 from the Allard tables with
! allard_use_tau100). Both are set only when the lookup succeeds.
! lookup_failed: the point is outside the Allard tables (the caller
! falls back to a gray atmosphere). ierr /= 0: a table interpolation
! failed. atm_choice is the caller's current choice, which the star
! layer switches between tabulated and gray as Teff changes.
subroutine atm_get_surface_pt(atm_choice, log_teff, log_g, log10_pressure, &
     log10_temperature, lookup_failed, ierr)
      integer, intent(in) :: atm_choice
      real(dp), intent(in) :: log_teff, log_g
      real(dp), intent(inout) :: log10_pressure, log10_temperature
      logical, intent(out) :: lookup_failed
      integer, intent(out) :: ierr

      ierr = 0
      lookup_failed = .false.
      if (atm_choice == 3) then
         call interpolate_surface_pressure(atm_table%kurucz_surface_p, log_teff, &
              log_g, .false., ierr)
      else if (atm_choice == 5) then
         call interpolate_surface_pressure(atm_table%castelli_surface_p, log_teff, &
              log_g, .false., ierr)
      else if (atm_choice == 4) then
         call interpolate_allard_surface(log_teff, log_g, .false., lookup_failed, ierr)
      end if
      if (ierr /= 0 .or. lookup_failed) return
      log10_pressure = atm_table%atm_log10_pressure
      log10_temperature = atm_table%atm_log10_temperature
end subroutine atm_get_surface_pt

!----------------------------------------------------------------------
! atm_allard_teff_max
!----------------------------------------------------------------------
! The highest log Teff of the Allard tables read (the star layer
! switches to a gray atmosphere above it).
pure function atm_allard_teff_max() result(log_teff_max)
      real(dp) :: log_teff_max
      log_teff_max = atm_table%allard_al_teffl_max
end function atm_allard_teff_max

end module atm_lib
