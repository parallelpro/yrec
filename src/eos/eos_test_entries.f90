!----------------------------------------------------------------------
! eos_test_entries
!----------------------------------------------------------------------
! What eos's test program (eos/test/test_eos.f90) reaches inside eos
! beyond eos_lib, in one place, so that no production module carries a
! procedure for the test alone. Its own procedures:
!     gamma1_at_density      gamma1 and grad_ad at (T6, rho): the OPAL
!                            2006 table, else the blended EOS
!     esac95, esac01, esac06 the OPAL interpolators' entry points, by
!                            OPAL's names
! and, re-exported from the modules that hold them, the internals the
! checks call or inspect: the settings' installer, the OPAL lookup and
! its tables' instances, the MHD table lookup, SCV's search start, its
! envelope table (whose never-built entries a check fills) and its H
! table (whose particle fractions a check reads), the warnings'
! flags, the table sources' list and the Debye-Hueckel record.
module eos_test_entries
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use physical_constants, only: ln10
      use eos_configuration, only: eos_settings, set_eos_settings, opal_vintage_2006
      use blended_eos, only: compute_blended_eos
      use opal_eos, only: compute_opal_eos, interpolate_opal_eos
      use opal_eos_tables, only: opal_tables, opal_eos_vintage, iv_opal95, &
           iv_opal01, iv_opal06, i_opal_cv, i_opal01_cv, i_opal_p, i_opal01_p, i_opal_gamma1, &
           i_opal_gamma2_ratio
      use mhd_eos, only: interpolate_mhd_table
      use scv_eos_tables, only: scv_cell, last_scv_cell, tablenv, nptsx, tablex, tlogx, &
           iscv_log10_p, iscv_frac_neutral, iscv_frac_atom_ion1
      use eos_messages, only: eos_warned, once_fermi_extrapolation
      use eos_table_source, only: eos_table_source_t, stage_in_place_of_saha, stage_over_the_gas, &
           stage_whole_eos
      use eos_source_list, only: table_sources, n_table_sources
      use eos_mixture, only: debye_huckel_composition_of
      use yale_eos_tables, only: debye_huckel_composition, idh_c, idh_n, idh_o, idh_ne, &
           idh_na, idh_mg, idh_al, idh_si, idh_ar, idh_fe
      use math_lib
      implicit none
      private
      public :: gamma1_at_density, esac95, esac01, esac06
      public :: eos_settings, set_eos_settings, compute_opal_eos, opal_tables, &
           iv_opal95, iv_opal01, iv_opal06, i_opal_cv, i_opal01_cv, i_opal_p, i_opal01_p, &
           interpolate_mhd_table, scv_cell, last_scv_cell, tablenv, nptsx, tablex, tlogx, &
           iscv_log10_p, iscv_frac_neutral, iscv_frac_atom_ion1, eos_warned, &
           once_fermi_extrapolation, &
           eos_table_source_t, stage_in_place_of_saha, stage_over_the_gas, stage_whole_eos, &
           table_sources, n_table_sources, debye_huckel_composition_of, debye_huckel_composition, &
           idh_c, idh_n, idh_o, idh_ne, idh_na, idh_mg, idh_al, idh_si, idh_ar, idh_fe

contains

!----------------------------------------------------------------------
! gamma1_at_density (formerly eos_lib's eos_get_gamma1)
!----------------------------------------------------------------------
! gamma1 and the adiabatic gradient at (T6, rho, P) and (X, Z): from
! the OPAL 2006 table (esac06, with the radiation correction) when
! the 2006 vintage is selected and the table covers the point, else from
! the blended EOS (compute_blended_eos, no second derivatives), with
! gamma1 built from chi_rho, chi_T and cv (Cox & Giuli
! ch. 9). test_eos is its only caller.
subroutine gamma1_at_density(hydrogen_fraction, metal_fraction, &
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
      type(debye_huckel_composition) :: dh
      real(dp) :: chi_rho, chi_t, specific_heat_cv
      logical :: in_table

      ierr = 0
      jerr = 0
      in_table = .false.

      if (eos_settings%opal_vintage == opal_vintage_2006) then
         call esac06(opal_tables(iv_opal06), hydrogen_fraction, temperature_1e6k, density, &
              opal06_variables, with_radiation, jerr, in_table)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      end if
      if (eos_settings%opal_vintage == opal_vintage_2006 .and. in_table) then
         gamma1 = opal_tables(iv_opal06)%eos_output(i_opal_gamma1)
         adiabatic_gradient = 1.0d0/opal_tables(iv_opal06)%eos_output(i_opal_gamma2_ratio)
      else
! (rho is the density guess; P and T are set from their logs)
         g%density = density
         g%log10_density = log10(density)
         log10_temperature = log10(temperature_1e6k*1.0d6)
         log10_pressure = log10(pressure)
         if (eos_settings%use_debye_huckel_correction) &
              dh = debye_huckel_composition_of(hydrogen_fraction, metal_fraction)
         call compute_blended_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, dh, g, .false., .false., saha_iterations, jerr)
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

end subroutine gamma1_at_density

!----------------------------------------------------------------------
! esac95, esac01, esac06
!----------------------------------------------------------------------
! OPAL's entry points for the 1995, 2001 and 2006 tables
! (opal_eos's interpolate_opal_eos), as test_eos calls them.
subroutine esac95(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)
      type(opal_eos_vintage), intent(inout) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, density
      integer, intent(in) :: deriv_order, rad_flag
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
      call interpolate_opal_eos(v, iv_opal95, hydrogen_fraction, t6_temperature, &
           density, deriv_order, rad_flag, ierr, in_table)
end subroutine esac95

subroutine esac01(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)
      type(opal_eos_vintage), intent(inout) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, density
      integer, intent(in) :: deriv_order, rad_flag
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
      call interpolate_opal_eos(v, iv_opal01, hydrogen_fraction, t6_temperature, &
           density, deriv_order, rad_flag, ierr, in_table)
end subroutine esac01

subroutine esac06(v, hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, in_table)
      type(opal_eos_vintage), intent(inout) :: v
      real(dp), intent(in) :: hydrogen_fraction, t6_temperature, density
      integer, intent(in) :: deriv_order, rad_flag
      integer, intent(out) :: ierr
      logical, intent(out) :: in_table
      call interpolate_opal_eos(v, iv_opal06, hydrogen_fraction, t6_temperature, &
           density, deriv_order, rad_flag, ierr, in_table)
end subroutine esac06

end module eos_test_entries
