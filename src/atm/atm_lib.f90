!----------------------------------------------------------------------
! atm_lib
!----------------------------------------------------------------------
! Added 2026 as part of the YREC readability refactor's phase two
! (disentangling the solver from the physics domains -- see
! docs/GUIDELINES.md's "Physics domains still entangled with the solver").
! Public facade of the atm domain (same shape as eos_lib/kap_lib):
! atm_init loads the tabulated-atmosphere surface-pressure tables at
! startup, and atm_get_surface_pt is the Allard-table surface lookup
! used from outside the domain. The envelope/atmosphere integrator
! itself (formerly envint, "atm_get") lives in core/envelope_integration.f90;
! the table interpolators are in atm/tables/ (surfp, kcsurfp,
! interpolate_allard_surface) and the T-tau relations in t_tau_relations.f90.
module atm_lib
      use atm_tables, only: atm_settings_t, set_atm_settings, atm_table, &
           atm_table_nt, atm_table_ntc
      use kurucz_atmosphere, only: read_surface_pressure_table, interpolate_surface_pressure
      use allard_atmosphere, only: read_allard_tables, interpolate_allard_surface
      use t_tau_relations, only: ttau_log10_temperature, &
           ttau_start_log10_temperature, ttau_photosphere_x_limit, hsra_t_tau_offset
      implicit none
      private
      public :: atm_init, atm_get_surface_pt, atm_allard_teff_max
      public :: atm_settings_t, set_atm_settings
      public :: ttau_log10_temperature, ttau_start_log10_temperature, &
           ttau_photosphere_x_limit, hsra_t_tau_offset
contains

!----------------------------------------------------------------------
! atm_init
!----------------------------------------------------------------------
! Added 2026 (phase three, docs/ROADMAP.md stage 1): the atm domain's
! startup-time table-load lifecycle entry, following MESA's
! <mod>_init convention. Covers, per the atm_choice option in force:
! the Kurucz surface-pressure table (atm_choice 3 and 4, into
! const_lib's kurucz_* tables plus atm_table's gmin/gmax edge
! indices), the Allard NextGen atmosphere tables (atm_choice 4, via
! tables/read_allard_tables.f90), and the Kurucz/Castelli surface-pressure
! table (atm_choice 5). All of this previously lived as an inline
! open/read block in setup/setups.f90 -- moved here verbatim,
! including the G Somers 5/15 invalid-pressure-edge catches and the
! table sizes atm_table_nt/ntc x atm_table_ng from atm_tables.
!
! 2026 wave 3 (R5): the Kurucz and Kurucz/Castelli tables are the two
! surface_p_table instances atm_table%kurucz_surface_p (atm_table_nt
! rows) and castelli_surface_p (atm_table_ntc rows). Each is
! allocated on first use with exactly the extents the fixed arrays
! had and zero-filled as they were (allocate_surface_p); the READ
! statements are unchanged apart from the names; the identical
! G Somers 5/15 edge scans became the internal find_gravity_edges,
! its working scalars host-associated.
subroutine atm_init(settings, atm_table_path, allard_table_path, ierr)
      type(atm_settings_t), intent(in) :: settings
      character(len=*), intent(in) :: atm_table_path
      character(len=*), intent(in) :: allard_table_path
! 2026 (docs/ROADMAP.md stage 3): ierr-not-stop; table-load failures print
! their diagnostic at the point of failure and come back as ierr /= 0.
      integer, intent(out) :: ierr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT PRESSURE TABLE FOR SURFACE BOUNDARY CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = 0
      call set_atm_settings(settings)

      if (settings%atm_choice == 3 .or. settings%atm_choice == 4) then
         call read_surface_pressure_table(atm_table_path, atm_table%kurucz_surface_p, &
              atm_table_nt, ierr)
         if (ierr /= 0) return
! MHP 6/97: the Allard model atmospheres, used instead of Kurucz below
! Teff = 10,000 K
         if (settings%atm_choice == 4) then
            call read_allard_tables(allard_table_path, ierr)
            if (ierr /= 0) return
         end if
! JNT 6/14: the Kurucz/Castelli atmospheres
      else if (settings%atm_choice == 5) then
         call read_surface_pressure_table(atm_table_path, atm_table%castelli_surface_p, &
              atm_table_ntc, ierr)
         if (ierr /= 0) return
      end if
end subroutine atm_init

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
      double precision, intent(in) :: log_teff, log_g
      double precision, intent(inout) :: log10_pressure, log10_temperature
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
function atm_allard_teff_max() result(log_teff_max)
      double precision :: log_teff_max
      log_teff_max = atm_table%allard_al_teffl_max
end function atm_allard_teff_max

end module atm_lib
