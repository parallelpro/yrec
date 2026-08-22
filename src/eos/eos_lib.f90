!----------------------------------------------------------------------
! eos_lib
!----------------------------------------------------------------------
! Added 2026 as part of the YREC readability refactor's phase two
! (disentangling the solver from the physics domains -- see
! GUIDELINES.md's "Physics domains still entangled with the solver").
!
! eos_get is the first facade of this phase: a single explicit-
! interface entry point that replaces the `if (use_mhd_eos) call meqos
! else call eqstat` dispatch that used to be duplicated at every call
! site, and additionally centralizes the Debye-Huckel composition
! setup (debye_huckel_x/y/z_total/z(3), read only by eos/yale/eqrelv.f90
! on the non-MHD path) that 6 of those call sites also used to
! duplicate verbatim. eos/eqstat.f90 (which also hosts eqstat2, its
! co-located pair), eos/yale/eqrelv.f90, and eos/mhd/meqos.f90 are all
! unchanged; this is a pure dispatch/setup wrapper.
!
! composition_at_zone is OPTIONAL: callers that carry full per-species
! composition (misc/coefft.f90, io/wrtout.f90, core/starin.f90,
! misc/physic.f90, mixing/hsubp.f90, mixing/sconvec.f90) pass
! composition(:,idx) and get the same Debye-Huckel setup they used to
! compute themselves. Callers that only ever carried bulk
! hydrogen_fraction/metal_fraction (atm/atm_lib.f90, atm/qatm.f90,
! atm/qenv.f90, wind/massloss.f90) omit it -- exactly matching their
! prior behavior, since none of them set these fields before either.
!
! core/starin.f90 previously had a bug here (confirmed against the
! original F77 source via git history): a missing ELSE meant it called
! *both* meqos and eqstat when MHD was on, and *neither* when MHD was
! off. eos_get's if/else has the structurally-correct form, so
! migrating that call site fixes the bug by construction.
!
! mixing/hsubp.f90, mixing/sconvec.f90, and wind/massloss.f90
! previously called eqstat unconditionally, with no LMHD check at all
! -- confirmed authentic original YREC behavior (unchanged since the
! very first commit, not a modernization artifact). Migrating them to
! eos_get extends real MHD support to these three secondary/diagnostic
! calculations for the first time, per explicit user sign-off; this is
! an acknowledged numerics change for use_mhd_eos=.true. runs, which
! the Stage-0 regression suite cannot verify (no test case sets LMHD).
!
! atm/turnover/calcad.f90 is NOT migrated: it calls esac06 directly (bypassing
! eqstat2's boundary-ramping) under its own use_opal2006_eos check for
! a self-contained acoustic-depth diagnostic, and never checks
! use_mhd_eos at all -- confirmed to match the original F77 exactly,
! a deliberate design choice, not part of this dispatch pattern.
module eos_lib
      implicit none
contains

subroutine eos_get(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
     in_atmosphere, saha_state, composition_at_zone)

      use const_lib
      implicit none

      double precision, intent(inout) :: log10_temperature
      double precision, intent(out) :: temperature
      double precision, intent(inout) :: log10_pressure
      double precision, intent(out) :: pressure
      double precision, intent(inout) :: log10_density
      double precision, intent(out) :: density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(inout) :: beta
      double precision, intent(out) :: beta_inverse, beta14
      double precision, intent(inout) :: ion_fraction(3)
      double precision, intent(out) :: specific_gas_constant
      double precision, intent(inout) :: ion_mean_weight_inverse
      double precision, intent(out) :: electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision, intent(inout) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      double precision, intent(out) :: dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
      double precision, intent(in), optional :: composition_at_zone(15)

      if (use_mhd_eos) then
         call meqos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp)
      else
         if (present(composition_at_zone) .and. use_debye_huckel_correction) then
            debye_huckel_x = composition_at_zone(1)
            debye_huckel_y = composition_at_zone(2)+composition_at_zone(4)
            debye_huckel_z_total = composition_at_zone(3)
            debye_huckel_z(1) = composition_at_zone(5)+composition_at_zone(6)
            debye_huckel_z(2) = composition_at_zone(7)+composition_at_zone(8)
            debye_huckel_z(3) = composition_at_zone(9)+composition_at_zone(10)+ &
                 composition_at_zone(11)
         end if
         call eqstat(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
              in_atmosphere, saha_state)
      end if

      return
end subroutine eos_get

!----------------------------------------------------------------------
! eos_init
!----------------------------------------------------------------------
! Added 2026 (phase three, ROADMAP.md stage 1): the eos domain's
! startup-time table-load lifecycle entry, following MESA's
! <mod>_init convention. Currently covers the MHD EOS tables (the
! only eos tables loaded from named files at startup; the OPAL
! 1995/2001/2006 EOS tables are lazily read on first use inside
! esac/esac01/esac06 via their lreadco guards). The use_mhd_eos gate
! previously lived at the call site (setup/setups.f90) -- moved
! inside so the caller needs no knowledge of which eos configuration
! requires which tables. Behavior is identical: when use_mhd_eos is
! false this is a no-op, exactly as the gated call was.
subroutine eos_init(zams_a_table_path, zams_b_table_path, &
     zams_c_table_path, centre1_table_path, centre2_table_path, &
     centre3_table_path, centre4_table_path, centre5_table_path)

      use const_lib
      implicit none

      character(len=256), intent(in) :: zams_a_table_path, &
           zams_b_table_path, zams_c_table_path, centre1_table_path, &
           centre2_table_path, centre3_table_path, centre4_table_path, &
           centre5_table_path

      if (use_mhd_eos) then
         call mhdtbl(zams_a_table_path, zams_b_table_path, &
              zams_c_table_path, centre1_table_path, centre2_table_path, &
              centre3_table_path, centre4_table_path, centre5_table_path)
      end if

      return
end subroutine eos_init

!----------------------------------------------------------------------
! eos_get_gamma1
!----------------------------------------------------------------------
! Added 2026 (phase three, ROADMAP.md stage 1): a component-accessor
! entry in the spirit of MESA's eosDT_get_component, created so
! atm/turnover/calcad.f90 (the acoustic-depth diagnostic, previously
! the last file calling eos internals -- esac06/eqstat2 -- directly)
! can go through the facade. Given composition and an (unlogged)
! T/rho/P point, returns gamma1 and the adiabatic gradient,
! dispatching on use_opal2006_eos exactly as calcad historically did:
! the OPAL 2006 table (via esac06, density-basis, radiation
! correction on, 9 variables; a failed table lookup falls through
! with whatever the previous point left behind, preserved verbatim
! from calcad's original alternate-return handling) when enabled,
! otherwise the Yale/SCV path (via eqstat2, pressure-basis) with
! gamma1 built from chi_rho/chi_T/cv per Cox & Giuli ch. 9 --
! statement-for-statement the arithmetic calcad performed at its own
! call sites, so the migration is byte-identical.
subroutine eos_get_gamma1(hydrogen_fraction, metal_fraction, &
     temperature_1e6k, density, pressure, gamma1, adiabatic_gradient, &
     saha_state)

      use atm_table_lib
      use const_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, metal_fraction, &
           temperature_1e6k, density, pressure
      double precision, intent(out) :: gamma1, adiabatic_gradient
      integer, intent(inout) :: saha_state

! --- locals (per-call copies: esac06/eqstat2 write into several of
!     the positions the inputs occupy) ---
      integer :: eos_interp_order, eos_rad_flag
      double precision :: x_local, t6_local, d_local
      double precision :: temperature_local, log10_temperature, &
           log10_pressure, pressure_local, log10_density
      double precision :: beta_local, beta_inverse_local, beta14_local, &
           ion_fraction_local(3), specific_gas_constant_local, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision :: qdt_eos, qdp_eos, qcp_eos, dela_eos, qdtt_eos, &
           qdtp_eos, qat_eos, qap_eos, qcpt_eos, qcpp_eos
      double precision :: chi_rho, chi_t, specific_heat_cv
      logical :: eos_deriv_flag, eos_atmosphere_flag

      if (use_opal2006_eos) then
         eos_interp_order = 9
         eos_rad_flag = 1
         x_local = hydrogen_fraction
         t6_local = temperature_1e6k
         d_local = density
         call esac06(x_local, t6_local, d_local, eos_interp_order, &
              eos_rad_flag, *100)
  100    continue
         gamma1 = atm_table%eos_output(8)
         adiabatic_gradient = 1.0d0/atm_table%eos_output(9)
      else
         x_local = hydrogen_fraction
         temperature_local = temperature_1e6k*1.0d6
         d_local = density
         log10_density = dlog10(density)
         log10_temperature = dlog10(temperature_local)
         pressure_local = pressure
         log10_pressure = dlog10(pressure_local)
         eos_deriv_flag = .false.
         eos_atmosphere_flag = .true.
         call eqstat2(log10_temperature, temperature_local, &
              log10_pressure, pressure_local, log10_density, d_local, &
              x_local, metal_fraction, beta_local, beta_inverse_local, &
              beta14_local, ion_fraction_local, &
              specific_gas_constant_local, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              qdt_eos, qdp_eos, qcp_eos, dela_eos, qdtt_eos, qdtp_eos, &
              qat_eos, qap_eos, qcpt_eos, qcpp_eos, eos_deriv_flag, &
              eos_atmosphere_flag, saha_state)
         chi_rho = 1.0d0/qdp_eos
         chi_t = -chi_rho*qdt_eos
         specific_heat_cv = qcp_eos - exp(ln10*(log10_pressure - &
              log10_density - log10_temperature))*chi_t**2/chi_rho
         gamma1 = chi_rho*qcp_eos/specific_heat_cv
         adiabatic_gradient = dela_eos
      end if

      return
end subroutine eos_get_gamma1

end module eos_lib
