!----------------------------------------------------------------------
! eos_lib
!----------------------------------------------------------------------
! The public interface of eos/, the equation of state: everything code
! outside eos/ may use, here and only here (check_boundaries enforces
! it).
!     eos_init        install the settings (eos_settings_t, built by the
!                     star layer's eos_settings_from_star) and load
!                     the tables they select
!     eos_get         density, the thermodynamic quantities and their
!                     derivatives at (log T, log P, X, Z), in an
!                     eos_result record (eos_results says what each
!                     field is)
!     eos_set_mixture the envelope mixture (and, with SCV, the envelope
!                     table built from it)
! eos_get takes a selected whole-EOS table source's result alone (MHD,
! when use_mhd_eos is set), else the blended EOS (blended_eos): the Yale
! EOS, with the table sources (eos_source_list: the SCV table at low
! temperature, an OPAL table blended in where it covers the point)
! applied over it. Each EOS family is one module; eos reads no star.
module eos_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use eos_configuration, only: eos_settings, eos_settings_t, &
           set_eos_settings, eos_table_paths_t, opal_vintage_none, opal_vintage_1995, &
           opal_vintage_2001, opal_vintage_2006, opal_vintage_conflict, saha_ramp_width
      use blended_eos, only: compute_blended_eos
      use eos_table_source, only: eos_table_source_t, stage_whole_eos
      use eos_source_list, only: table_sources, n_table_sources
      use scv_eos, only: build_scv_envelope_table
      use fermi_table_reader, only: read_fermi_table
      use eos_mixture, only: eos_mix, debye_huckel_composition_of
      use opal_eos_tables, only: n_opal_vintages
      use yale_eos_tables, only: debye_huckel_composition, fermi_log10t_min
      use eos_messages, only: eos_report, reset_eos_warnings, message_length, unwritten_line
      implicit none
      private
      public :: eos_init, eos_get, eos_result
      public :: eos_set_mixture
      public :: eos_settings_t, eos_table_paths_t
      public :: opal_vintage_none, opal_vintage_1995, opal_vintage_2001, &
           opal_vintage_2006, opal_vintage_conflict

contains

!----------------------------------------------------------------------
! eos_get
!----------------------------------------------------------------------
! The EOS at (log T, log P, X, Z), in res (see eos_results for what each
! field is, and which fields are also inputs). want_derivatives asks for
! the six second derivatives; in_atmosphere skips the Saha solver's
! second derivatives (the first ones are always computed);
! saha_iterations is a counter the Saha solver adds its Newton iterations
! to (eos reads it nowhere; the solver's starting point is res's saha_*
! fields, the record's previous solution).
! composition_at_zone (a zone's composition, star%xa layout) gives the
! Debye-Hueckel correction its composition (debye_huckel_composition_of;
! every star-layer call passes it, and without it, as in the unit
! tests, the envelope's mixture at the zone's X and Z stands in). ierr
! /= 0 reports a corrupt table, a bad configuration or a point no
! source covers; the message is already out (eos_report).
subroutine eos_get(log10_temperature, log10_pressure, &
     hydrogen_fraction, metal_fraction, res, want_derivatives, &
     in_atmosphere, saha_iterations, composition_at_zone, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: res
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_iterations
      real(dp), intent(in), optional :: composition_at_zone(:)
      integer, intent(out) :: ierr

      type(debye_huckel_composition) :: dh
      type(eos_table_source_t) :: sources(n_table_sources)
      logical :: covers
      integer :: k

      ierr = 0
      sources = table_sources()
      do k = 1, n_table_sources
         if (sources(k)%stage /= stage_whole_eos) cycle
         if (.not. sources(k)%selected()) cycle
         call sources(k)%lookup(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, res, covers, ierr)
         return
      end do
      if (eos_settings%use_debye_huckel_correction) dh = &
           debye_huckel_composition_of(hydrogen_fraction, metal_fraction, composition_at_zone)
      call compute_blended_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, dh, res, want_derivatives, in_atmosphere, saha_iterations, ierr)
end subroutine eos_get

!----------------------------------------------------------------------
! eos_init
!----------------------------------------------------------------------
! Installs the settings and loads the tables they select, from paths:
! reads the Fermi-Dirac table (always), then each selected table
! source's tables in eos_source_list's order (the SCV tables; the OPAL
! table, which is opened here and read on first use by its vintage's
! reader; the MHD tables). Each file is opened on a new unit; a missing
! or corrupt one returns ierr = 1, the later ones not read, and so do
! settings that select more than one OPAL
! vintage, a Debye-Hueckel degeneracy range with eta_max <= eta_min, or
! a Saha cutoff whose ramp would reach below the Fermi-Dirac table's
! first row (log T 5.2: the cutoff must be at least 5.7), so that the
! fully ionized gas never extrapolates the table below it (nothing is
! installed then).
subroutine eos_init(settings, paths, ierr)
      type(eos_settings_t), intent(in) :: settings
      type(eos_table_paths_t), intent(in) :: paths
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)
      type(eos_table_source_t) :: sources(n_table_sources)
      integer :: k

      ierr = 0
      if (settings%opal_vintage < opal_vintage_none .or. &
           settings%opal_vintage > n_opal_vintages) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'eos_init: more than one OPAL EOS vintage selected'
         call eos_report(eos_lines)
         ierr = 1
         return
      end if
      if (settings%use_debye_huckel_correction .and. &
           .not. settings%debye_huckel_eta_max > settings%debye_huckel_eta_min) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'eos_init: the Debye-Hueckel degeneracy range needs eta_max > eta_min:', &
              settings%debye_huckel_eta_min, settings%debye_huckel_eta_max
         call eos_report(eos_lines)
         ierr = 1
         return
      end if
      if (.not. settings%saha_log10t_cutoff - saha_ramp_width >= fermi_log10t_min) then
         eos_lines = unwritten_line
         write(eos_lines,'(a,f4.2,a,f4.2,a,f4.2,a,f0.3)') 'eos_init: saha_log10t_cutoff must be at least ', &
              fermi_log10t_min + saha_ramp_width, ' (the Fermi-Dirac table''s first row, log T ', &
              fermi_log10t_min, ', plus the Saha ramp''s ', saha_ramp_width, ' dex): ', &
              settings%saha_log10t_cutoff
         call eos_report(eos_lines)
         ierr = 1
         return
      end if
      call set_eos_settings(settings)
      call reset_eos_warnings()
      call read_fermi_table(trim(paths%fermi), ierr)
      if (ierr /= 0) return
      sources = table_sources()
      do k = 1, n_table_sources
         if (.not. sources(k)%selected()) cycle
         call sources(k)%read_tables(paths, ierr)
         if (ierr /= 0) return
      end do
end subroutine eos_init



!----------------------------------------------------------------------
! eos_set_mixture
!----------------------------------------------------------------------
! The envelope mixture the Saha and SCV paths read through eos_mix: X,
! Z, the inverse mean atomic weight sum_i X_i/A_i and the species'
! number fractions (X_i/A_i over that sum). The star layer
! sets it whenever it recomputes the envelope composition, and with SCV
! on the SCV envelope table is built for its X and Z here (once the
! tables are read; a table already built for them stays). The
! Debye-Hueckel correction takes its metals heavier than O from it, in
! its proportions (debye_huckel_composition_of).
subroutine eos_set_mixture(envelope_hydrogen_fraction, &
     envelope_metal_fraction, inverse_atomic_weight, number_fractions)
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, inverse_atomic_weight
      real(dp), intent(in) :: number_fractions(:)

      eos_mix%envelope_hydrogen_fraction = envelope_hydrogen_fraction
      eos_mix%envelope_metal_fraction = envelope_metal_fraction
      eos_mix%inverse_atomic_weight = inverse_atomic_weight
      eos_mix%number_fractions = number_fractions
      if (eos_settings%use_scv_eos) call build_scv_envelope_table( &
           envelope_hydrogen_fraction, envelope_metal_fraction)
end subroutine eos_set_mixture



end module eos_lib
