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
!     eos_set_mixture, eos_set_debye_huckel_z   the envelope mixture and
!                     the metal mixture of the Debye-Hueckel correction
!     eos_set_debye_huckel_composition   the zone composition the
!                     Debye-Hueckel correction uses
!     eos_update_envelope_tables   rebuild the SCV envelope table when
!                     the envelope composition changes
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
           set_eos_settings
      use blended_eos, only: compute_blended_eos, blend_eos
      use mhd_eos, only: compute_mhd_eos, read_mhd_tables
      use opal2006_eos, only: esac06
      use scv_eos, only: build_scv_envelope_table, read_scv_tables
      use yale_eos, only: read_fermi_table
      use eos_mixture, only: eos_mix
      use opal_eos_tables, only: i_opal_gamma1, i_opal_gamma2_ratio, opal06, &
           open_opal_eos_table
      use yale_eos_tables, only: debye_huckel_x, debye_huckel_y, debye_huckel_z, &
           debye_huckel_z_total
      use physical_constants, only: ln10
      use math_lib
      implicit none
      private
      public :: eos_init, eos_get, eos_get_gamma1, eos_result
      public :: eos_set_mixture, eos_set_debye_huckel_z, &
           eos_set_debye_huckel_composition, eos_update_envelope_tables
      public :: eos_settings_t, set_eos_settings

contains

!----------------------------------------------------------------------
! eos_get
!----------------------------------------------------------------------
! The EOS at (log T, log P, X, Z), in res (see eos_results for what each
! field is, and which fields are also inputs). want_derivatives asks for
! the six second derivatives; in_atmosphere skips the Saha solver's
! derivatives; saha_state carries that solver's state between calls.
! composition_at_zone (a zone's composition, star%xa layout) sets the
! Debye-Hueckel composition; callers with only X and Z omit it. ierr /= 0
! reports a corrupt table or a bad configuration; the diagnostic is
! already in the run log.
subroutine eos_get(log10_temperature, log10_pressure, &
     hydrogen_fraction, metal_fraction, res, want_derivatives, &
     in_atmosphere, saha_state, composition_at_zone, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: res
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
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
              metal_fraction, res, want_derivatives, in_atmosphere, saha_state, ierr)
      end if
end subroutine eos_get

!----------------------------------------------------------------------
! eos_init
!----------------------------------------------------------------------
! Installs the settings and loads the tables they select: opens the
! OPAL table (its vintage's reader reads it on first use) and reads the
! MHD tables, the Fermi-Dirac table (always) and the SCV tables. Each
! file is opened on a new unit; a missing one returns ierr = 1.
subroutine eos_init(settings, fermi_table_path, scv_h_table_path, &
     scv_he_table_path, scv_z_table_path, zams_a_table_path, &
     zams_b_table_path, zams_c_table_path, centre1_table_path, &
     centre2_table_path, centre3_table_path, centre4_table_path, &
     centre5_table_path, opal_eos_table_path, ierr)

      character(len=*), intent(in) :: fermi_table_path
      character(len=*), intent(in) :: scv_h_table_path, &
           scv_he_table_path, scv_z_table_path
      character(len=*), intent(in) :: zams_a_table_path, &
           zams_b_table_path, zams_c_table_path, centre1_table_path, &
           centre2_table_path, centre3_table_path, centre4_table_path, &
           centre5_table_path
! the run's EOS settings (installed here); the OPAL table the settings
! select (opened here, read lazily by its vintage's reader)
      type(eos_settings_t), intent(in) :: settings
      character(len=*), intent(in) :: opal_eos_table_path
      integer, intent(out) :: ierr

      ierr = 0
      call set_eos_settings(settings)
      if (settings%use_opal95_eos .or. settings%use_opal2001_eos .or. &
           settings%use_opal2006_eos) then
         call open_opal_eos_table(opal_eos_table_path, ierr)
         if (ierr /= 0) return
      end if
      if (settings%use_mhd_eos) then
         call read_mhd_tables(zams_a_table_path, zams_b_table_path, &
              zams_c_table_path, centre1_table_path, centre2_table_path, &
              centre3_table_path, centre4_table_path, centre5_table_path, &
              ierr)
         if (ierr /= 0) return
      end if
      call read_fermi_table(fermi_table_path, ierr)
      if (ierr /= 0) return
      if (settings%use_scv_eos) call read_scv_tables(scv_h_table_path, &
           scv_he_table_path, scv_z_table_path, ierr)
end subroutine eos_init

!----------------------------------------------------------------------
! eos_get_gamma1
!----------------------------------------------------------------------
! gamma1 and the adiabatic gradient at (T6, rho, P) and (X, Z): from
! the OPAL 2006 table (esac06, with the radiation correction) when
! use_opal2006_eos is set, else from blend_eos, with gamma1 built from
! chi_rho, chi_T and cv (Cox & Giuli ch. 9). Outside the OPAL table it
! returns what esac06 left from the previous point. test_eos is its
! only caller.
subroutine eos_get_gamma1(hydrogen_fraction, metal_fraction, &
     temperature_1e6k, density, pressure, gamma1, adiabatic_gradient, &
     saha_state, ierr)

      real(dp), intent(in) :: hydrogen_fraction, metal_fraction, &
           temperature_1e6k, density, pressure
      real(dp), intent(out) :: gamma1, adiabatic_gradient
      integer, intent(inout) :: saha_state
! ierr as for eos_eval; an out-of-table esac06 lookup is not an error.
      integer, intent(out) :: ierr

      integer :: jerr

! per-call copies: esac06 and blend_eos write into several of their
! input arguments
      integer :: eos_interp_order, eos_rad_flag
      real(dp) :: x_local, t6_local, d_local
      real(dp) :: log10_temperature, log10_pressure
      type(eos_result) :: g
      real(dp) :: chi_rho, chi_t, specific_heat_cv
      logical :: in_table

      ierr = 0
      jerr = 0

      if (eos_settings%use_opal2006_eos) then
         eos_interp_order = 9
         eos_rad_flag = 1
         x_local = hydrogen_fraction
         t6_local = temperature_1e6k
         d_local = density
         call esac06(opal06, x_local, t6_local, d_local, eos_interp_order, &
              eos_rad_flag, jerr, in_table)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         gamma1 = opal06%eos_output(i_opal_gamma1)
         adiabatic_gradient = 1.0d0/opal06%eos_output(i_opal_gamma2_ratio)
      else
! (rho is the density guess; P and T are set from their logs)
         g%density = density
         g%log10_density = log10(density)
         log10_temperature = log10(temperature_1e6k*1.0d6)
         log10_pressure = log10(pressure)
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, g, .false., .true., saha_state, jerr)
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
! Z, the mean atomic mass and the species fractions. The star layer
! sets it whenever it recomputes the envelope composition.
subroutine eos_set_mixture(envelope_hydrogen_fraction, &
     envelope_metal_fraction, mean_atomic_mass, species_fractions)
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, mean_atomic_mass
      real(dp), intent(in) :: species_fractions(:)

      eos_mix%envelope_hydrogen_fraction = envelope_hydrogen_fraction
      eos_mix%envelope_metal_fraction = envelope_metal_fraction
      eos_mix%amuenv = mean_atomic_mass
      eos_mix%fxenv = species_fractions
end subroutine eos_set_mixture

!----------------------------------------------------------------------
! eos_set_debye_huckel_z
!----------------------------------------------------------------------
! The 18-element metal mixture (C, N, O, Ne, Na, Mg, Al, Si, P, S, Cl,
! Ar, Ca, Ti, Cr, Mn, Fe, Ni, summing to Z) that the Debye-Hueckel
! correction weights with debye_huckel_nu. kap_lib passes the LAOL89
! table's mixture after reading it; eos_set_debye_huckel_composition
! then resets elements 1-3 for each zone.
subroutine eos_set_debye_huckel_z(metal_mixture)
      real(dp), intent(in) :: metal_mixture(:)

      debye_huckel_z = metal_mixture
end subroutine eos_set_debye_huckel_z

!----------------------------------------------------------------------
! eos_set_debye_huckel_composition
!----------------------------------------------------------------------
! The composition the Debye-Hueckel correction uses, from a zone's
! composition (star%xa layout). eos_eval sets it from
! composition_at_zone; the star layer sets the surface zone's before
! the envelope integration, whose calls pass none. Does nothing when
! the correction is off.
subroutine eos_set_debye_huckel_composition(composition)
      real(dp), intent(in) :: composition(:)

      if (.not. eos_settings%use_debye_huckel_correction) return
      debye_huckel_x = composition(1)
      debye_huckel_y = composition(2)+composition(4)
      debye_huckel_z_total = composition(3)
      debye_huckel_z(1) = composition(5)+composition(6)
      debye_huckel_z(2) = composition(7)+composition(8)
      debye_huckel_z(3) = composition(9)+composition(10)+composition(11)
end subroutine eos_set_debye_huckel_composition

!----------------------------------------------------------------------
! eos_update_envelope_tables
!----------------------------------------------------------------------
! Rebuild the tables tied to the envelope composition: the SCV
! envelope-mixture table, when SCV is on.
subroutine eos_update_envelope_tables(envelope_hydrogen_fraction, &
     envelope_metal_fraction)
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction

      if (eos_settings%use_scv_eos) call build_scv_envelope_table( &
           envelope_hydrogen_fraction, envelope_metal_fraction)
end subroutine eos_update_envelope_tables

end module eos_lib
