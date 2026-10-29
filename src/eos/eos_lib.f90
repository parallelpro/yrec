!----------------------------------------------------------------------
! eos_lib
!----------------------------------------------------------------------
! The public interface of eos/, the equation of state: everything code
! outside eos/ may use, here and only here (check_boundaries enforces
! it).
!     eos_init        install the settings (eos_settings_t, built by the
!                     star layer's equation_of_state_settings) and load
!                     the tables they select
!     eos_get         density, the thermodynamic quantities and their
!                     derivatives at (log T, log P, X, Z), in an
!                     eos_result record (eos_results says what each
!                     field is)
!     eos_set_mixture the envelope mixture, and what eos derives from
!                     it: the heavier metals of the Debye-Hueckel
!                     correction and, with SCV, the envelope table
!     eos_get_gamma1  gamma1 and the adiabatic gradient at (T, rho)
!                     (for test_eos)
!     set_eos_settings   install settings without reloading tables
! eos_get takes the MHD EOS (mhd_eos) when use_mhd_eos is set, else
! the blended EOS (blended_eos): the Yale EOS, with the SCV table at
! low temperature and an OPAL table blended in where it covers the
! point. Each EOS family is one module; eos reads no star.
module eos_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use eos_configuration, only: eos_settings, eos_settings_t, &
           set_eos_settings, eos_table_paths_t, opal_vintage_none, opal_vintage_1995, &
           opal_vintage_2001, opal_vintage_2006, opal_vintage_conflict
      use blended_eos, only: compute_blended_eos, blend_eos
      use mhd_eos, only: compute_mhd_eos, read_mhd_tables
      use opal_eos, only: esac06
      use scv_eos, only: build_scv_envelope_table, read_scv_tables
      use yale_eos, only: read_fermi_table
      use eos_mixture, only: eos_mix
      use opal_eos_tables, only: i_opal_gamma1, i_opal_gamma2_ratio, opal06, &
           open_opal_eos_table, n_opal_vintages
      use yale_eos_tables, only: debye_huckel_x, debye_huckel_y, debye_huckel_z, &
           debye_huckel_z_total, idh_c, idh_n, idh_o, idh_ne, idh_na, idh_mg, idh_al, &
           idh_si, idh_ar, idh_fe
      use eos_mixture, only: ix_ne, ix_na, ix_mg, ix_al, ix_si, ix_ar, ix_fe, &
           mixture_atomic_weights
      use composition_layout, only: i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, &
           i_n14, i_n15, i_o16, i_o17, i_o18
      use physical_constants, only: ln10
      use eos_messages, only: eos_report, message_length, unwritten_line
      use math_lib
      implicit none
      private
      public :: eos_init, eos_get, eos_get_gamma1, eos_result
      public :: eos_set_mixture
      public :: eos_settings_t, set_eos_settings, eos_table_paths_t
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
! to (eos reads it nowhere; the solver's starting guess is its own
! module state, the previous call's solution).
! composition_at_zone (a zone's composition, star%xa layout) sets the
! Debye-Hueckel composition; every star-layer call passes it (omitted, as
! by the unit tests, the previous call's composition stays). ierr /= 0
! reports a corrupt table or a bad configuration; the diagnostic is
! already in the run log.
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

      ierr = 0
      if (eos_settings%use_mhd_eos) then
         call compute_mhd_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, res, ierr)
      else
         if (present(composition_at_zone)) &
              call eos_set_debye_huckel_composition(composition_at_zone)
         call compute_blended_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, res, want_derivatives, in_atmosphere, saha_iterations, ierr)
      end if
end subroutine eos_get

!----------------------------------------------------------------------
! eos_init
!----------------------------------------------------------------------
! Installs the settings and loads the tables they select, from paths:
! opens the OPAL table (its vintage's reader reads it on first use) and
! reads the MHD tables, the Fermi-Dirac table (always) and the SCV
! tables. Each file is opened on a new unit; a missing or corrupt one
! returns ierr = 1, and so do settings that select more than one OPAL
! vintage or a Debye-Hueckel degeneracy range with eta_max <= eta_min
! (nothing is installed then).
subroutine eos_init(settings, paths, ierr)
      type(eos_settings_t), intent(in) :: settings
      type(eos_table_paths_t), intent(in) :: paths
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)

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
      call set_eos_settings(settings)
      if (settings%opal_vintage /= opal_vintage_none) then
         call open_opal_eos_table(trim(paths%opal), ierr)
         if (ierr /= 0) return
      end if
      if (settings%use_mhd_eos) then
         call read_mhd_tables(paths%mhd_zams, paths%mhd_centre, ierr)
         if (ierr /= 0) return
      end if
      call read_fermi_table(trim(paths%fermi), ierr)
      if (ierr /= 0) return
      if (settings%use_scv_eos) call read_scv_tables(trim(paths%scv_h), &
           trim(paths%scv_he), trim(paths%scv_z), ierr)
end subroutine eos_init

!----------------------------------------------------------------------
! eos_get_gamma1
!----------------------------------------------------------------------
! gamma1 and the adiabatic gradient at (T6, rho, P) and (X, Z): from
! the OPAL 2006 table (esac06, with the radiation correction) when
! the 2006 vintage is selected and the table covers the point, else from
! blend_eos, with gamma1 built from chi_rho, chi_T and cv (Cox & Giuli
! ch. 9). test_eos is its only caller.
subroutine eos_get_gamma1(hydrogen_fraction, metal_fraction, &
     temperature_1e6k, density, pressure, gamma1, adiabatic_gradient, &
     saha_iterations, ierr)

      real(dp), intent(in) :: hydrogen_fraction, metal_fraction, &
           temperature_1e6k, density, pressure
      real(dp), intent(out) :: gamma1, adiabatic_gradient
      integer, intent(inout) :: saha_iterations
! ierr as for eos_get; an out-of-table esac06 lookup is not an error.
      integer, intent(out) :: ierr
! esac06's full evaluation (its nine variables) with the radiation
! correction
      integer, parameter :: opal06_variables = 9, with_radiation = 1

      integer :: jerr
      real(dp) :: log10_temperature, log10_pressure
      type(eos_result) :: g
      real(dp) :: chi_rho, chi_t, specific_heat_cv
      logical :: in_table

      ierr = 0
      jerr = 0
      in_table = .false.

      if (eos_settings%opal_vintage == opal_vintage_2006) then
         call esac06(opal06, hydrogen_fraction, temperature_1e6k, density, &
              opal06_variables, with_radiation, jerr, in_table)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      end if
      if (eos_settings%opal_vintage == opal_vintage_2006 .and. in_table) then
         gamma1 = opal06%eos_output(i_opal_gamma1)
         adiabatic_gradient = 1.0d0/opal06%eos_output(i_opal_gamma2_ratio)
      else
! (rho is the density guess; P and T are set from their logs)
         g%density = density
         g%log10_density = log10(density)
         log10_temperature = log10(temperature_1e6k*1.0d6)
         log10_pressure = log10(pressure)
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, g, .false., .false., saha_iterations, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         chi_rho = 1.0d0/g%dlnrho_dlnp
         chi_t = -chi_rho*g%dlnrho_dlnt
         specific_heat_cv = g%specific_heat_cp - exp(ln10*(log10_pressure - &
              g%log10_density - log10_temperature))*chi_t**2/chi_rho
         gamma1 = chi_rho*g%specific_heat_cp/specific_heat_cv
         adiabatic_gradient = g%adiabatic_gradient
      end if

end subroutine eos_get_gamma1


!----------------------------------------------------------------------
! eos_set_mixture
!----------------------------------------------------------------------
! The envelope mixture the Saha and SCV paths read through eos_mix: X,
! Z, the inverse mean atomic weight sum_i X_i/A_i and the species'
! number fractions (X_i/A_i over that sum). The star layer
! sets it whenever it recomputes the envelope composition, and with SCV
! on the SCV envelope table is built for its X and Z here (once the
! tables are read; a table already built for them stays). The
! Debye-Hueckel correction takes its metals heavier than O from it, as
! mass fractions X_i = fx_i A_i sum_j(X_j/A_j): Ne, Na, Mg, Al, Si, Ar
! and Fe (its P, S, Cl, Ca, Ti, Cr, Mn and Ni are not in the mixture and
! stay zero; C, N and O come from each zone's composition); here the
! mixture's own, and in each zone in the mixture's proportions with the
! zone's Z less its C, N and O for their sum
! (eos_set_debye_huckel_composition).
subroutine eos_set_mixture(envelope_hydrogen_fraction, &
     envelope_metal_fraction, inverse_atomic_weight, number_fractions)
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, inverse_atomic_weight
      real(dp), intent(in) :: number_fractions(:)

      eos_mix%envelope_hydrogen_fraction = envelope_hydrogen_fraction
      eos_mix%envelope_metal_fraction = envelope_metal_fraction
      eos_mix%inverse_atomic_weight = inverse_atomic_weight
      eos_mix%number_fractions = number_fractions
      call set_debye_huckel_heavy_metals()
      if (eos_settings%use_scv_eos) call build_scv_envelope_table( &
           envelope_hydrogen_fraction, envelope_metal_fraction)
end subroutine eos_set_mixture

! The Debye-Hueckel correction's metals heavier than O: the envelope
! mixture's mass fractions; given beyond_o, a zone's Z less its C, N
! and O, in the mixture's proportions with beyond_o for their sum, so
! that the correction's metals add up to the zone's Z (where a core
! burns helium its Z grows by the C and O it makes, not by Ne..Fe; none
! when the mixture has none).
subroutine set_debye_huckel_heavy_metals(beyond_o)
      real(dp), intent(in), optional :: beyond_o
      real(dp) :: metal_scale, envelope_beyond_o

      metal_scale = 1.0d0
      if (present(beyond_o)) then
         envelope_beyond_o = envelope_mass_fraction(ix_ne) + envelope_mass_fraction(ix_na) + &
              envelope_mass_fraction(ix_mg) + envelope_mass_fraction(ix_al) + &
              envelope_mass_fraction(ix_si) + envelope_mass_fraction(ix_ar) + &
              envelope_mass_fraction(ix_fe)
         metal_scale = 0.0d0
         if (envelope_beyond_o > 0.0d0) metal_scale = max(beyond_o, 0.0d0)/envelope_beyond_o
      end if
      debye_huckel_z(idh_ne) = metal_scale*envelope_mass_fraction(ix_ne)
      debye_huckel_z(idh_na) = metal_scale*envelope_mass_fraction(ix_na)
      debye_huckel_z(idh_mg) = metal_scale*envelope_mass_fraction(ix_mg)
      debye_huckel_z(idh_al) = metal_scale*envelope_mass_fraction(ix_al)
      debye_huckel_z(idh_si) = metal_scale*envelope_mass_fraction(ix_si)
      debye_huckel_z(idh_ar) = metal_scale*envelope_mass_fraction(ix_ar)
      debye_huckel_z(idh_fe) = metal_scale*envelope_mass_fraction(ix_fe)

contains

real(dp) function envelope_mass_fraction(i)
      integer, intent(in) :: i
      envelope_mass_fraction = eos_mix%number_fractions(i)*eos_mix%inverse_atomic_weight*mixture_atomic_weights(i)
end function envelope_mass_fraction

end subroutine set_debye_huckel_heavy_metals

!----------------------------------------------------------------------
! eos_set_debye_huckel_composition
!----------------------------------------------------------------------
! The composition the Debye-Hueckel correction uses, from a zone's
! composition (star%xa layout), set by eos_get from its
! composition_at_zone at every call (private since 2026-10: the star
! layer used to set the surface zone's before an envelope integration,
! whose calls then passed none, so the integration used whatever the
! last setter left): X, Y, Z, C, N and O the zone's, the heavier metals
! the envelope mixture's in its proportions, summing to the zone's Z less
! its C, N and O (set_debye_huckel_heavy_metals). Does nothing when the
! correction is off.
subroutine eos_set_debye_huckel_composition(composition)
      real(dp), intent(in) :: composition(:)

      if (.not. eos_settings%use_debye_huckel_correction) return
      debye_huckel_x = composition(i_h1)
      debye_huckel_y = composition(i_he4)+composition(i_he3)
      debye_huckel_z_total = composition(i_metals)
      debye_huckel_z(idh_c) = composition(i_c12)+composition(i_c13)
      debye_huckel_z(idh_n) = composition(i_n14)+composition(i_n15)
      debye_huckel_z(idh_o) = composition(i_o16)+composition(i_o17)+composition(i_o18)
      call set_debye_huckel_heavy_metals(composition(i_metals) - debye_huckel_z(idh_c) - &
           debye_huckel_z(idh_n) - debye_huckel_z(idh_o))
end subroutine eos_set_debye_huckel_composition


end module eos_lib
