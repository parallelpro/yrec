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
! atm/turnover/calcad.f90 (the acoustic-depth diagnostic) was NOT
! migrated to eos_get during phase two: it deliberately bypasses
! eqstat2's boundary-ramping, calling esac06 directly under its own
! use_opal2006_eos check, confirmed to match the original F77 -- not
! part of this dispatch pattern. As of phase three (ROADMAP.md stage
! 1) it goes through eos_get_gamma1 below instead, which preserves
! that same deliberate dispatch inside the facade boundary.
module eos_lib
      use scv_eos_lib
      use yale_eos_lib
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
     in_atmosphere, saha_state, composition_at_zone, ierr)

      use const_lib
      use scv_eos_lib
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
! 2026 (ROADMAP.md stage 3): OPTIONAL ierr, the transitional form of
! MESA's ierr-not-stop discipline (same shape as kap_lib's kap_get).
! When the caller passes ierr, any error condition that used to stop
! deep inside the eos internals is returned here instead, ierr /= 0,
! and no stop occurs. When the caller omits ierr, behavior is exactly
! historical: the same diagnostic messages, then a stop -- now located
! in this facade's funnel rather than at the point of failure.
      integer, intent(out), optional :: ierr

      integer :: jerr

      if (present(ierr)) ierr = 0
      jerr = 0

      if (use_mhd_eos) then
         call meqos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, jerr)
         if (jerr /= 0) go to 900
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
              in_atmosphere, saha_state, jerr)
         if (jerr /= 0) go to 900
      end if

      return

! error funnel: reached only when a callee reported jerr /= 0. With
! ierr present the caller takes responsibility; without it, preserve
! the historical stop (the diagnostic already printed at the point of
! failure).
  900 continue
      if (present(ierr)) then
         ierr = jerr
         return
      end if
      stop
end subroutine eos_get

!----------------------------------------------------------------------
! eos_init
!----------------------------------------------------------------------
! Added 2026 (phase three, ROADMAP.md stage 1): the eos domain's
! startup-time table-load lifecycle entry, following MESA's
! <mod>_init convention. Covers the eos tables loaded from named
! files at startup: the MHD EOS tables (gated on use_mhd_eos, moved
! inside from the call site so the caller needs no knowledge of
! which eos configuration requires which tables -- a no-op when MHD
! is off, exactly as the gated call was) and the Yale EOS's
! Fermi-Dirac degenerate-electron table (read unconditionally, into
! yale_eos_lib; this block previously lived inline in
! setup/setups.f90, writing into atm_table_lib -- moved here
! verbatim, including its bin-construction loop, its glitch check,
! and that check's `stop`, which ROADMAP.md stage 3 will convert to
! ierr). The OPAL 1995/2001/2006 EOS tables are not loaded here:
! they are lazily read on first use inside esac/esac01/esac06 via
! their lreadco guards. The SCV (Saumon-Chabrier-Van Horn) EOS tables
! are also read here (gated on use_scv_eos, into const_lib's former
! common/scveos/ state consumed by eos/scv/eqscve+eqscvg and the io
! writers) -- another inline setups.f90 block moved in verbatim; its
! read is provably independent of the atm tables read later in
! startup (separate file units, disjoint state), so absorbing it here
! is order-safe.
subroutine eos_init(fermi_table_path, scv_h_table_path, &
     scv_he_table_path, scv_z_table_path, zams_a_table_path, &
     zams_b_table_path, zams_c_table_path, centre1_table_path, &
     centre2_table_path, centre3_table_path, centre4_table_path, &
     centre5_table_path, ierr)

      use yale_eos_lib
      use luout_lib
      use const_lib
      use scv_eos_lib
      implicit none
      integer, parameter :: nts = 63

      character(len=256), intent(in) :: fermi_table_path
      character(len=256), intent(in) :: scv_h_table_path, &
           scv_he_table_path, scv_z_table_path
      character(len=256), intent(in) :: zams_a_table_path, &
           zams_b_table_path, zams_c_table_path, centre1_table_path, &
           centre2_table_path, centre3_table_path, centre4_table_path, &
           centre5_table_path
! 2026 (ROADMAP.md stage 3): OPTIONAL ierr, same contract as eos_get's.
      integer, intent(out), optional :: ierr

      integer :: jerr
      integer :: grid_idx, iden_idx, col_idx, card_idx
      integer :: bin_start, bin_width, bin_end
      integer :: t_idx, p_idx
      double precision :: scvhe_dummy_val, scvz_dummy_val
      integer :: scvhe_dummy_npts, scvz_dummy_npts

      if (present(ierr)) ierr = 0
      jerr = 0

      if (use_mhd_eos) then
         call mhdtbl(zams_a_table_path, zams_b_table_path, &
              zams_c_table_path, centre1_table_path, centre2_table_path, &
              centre3_table_path, centre4_table_path, centre5_table_path, &
              jerr)
         if (jerr /= 0) go to 900
      end if

!  INPUT F-TABLES FOR DEGENERATE EQUATION OF STATE
!  OPEN DATA FILES
      open(unit=fermi_unit,file=fermi_table_path,form='FORMATTED',status='OLD')
      read(fermi_unit,150)  (yale_eos%fermi_table_x_grid(grid_idx), grid_idx = 1,43)
      read(fermi_unit,150)  (yale_eos%fermi_table_eta(grid_idx), grid_idx = 1,43)
      read(fermi_unit,160) (grid_idx,iden_idx,(yale_eos%fermi_table_data(col_idx,grid_idx,iden_idx), &
           col_idx=1,5),card_idx=1,860)
  150 format(8x,9f8.4)
  160 format(2i5,5f12.7)
      rewind fermi_unit
      bin_start = 1
      do grid_idx=1,41
! MHP 10/02 SHOULD BE INT(DVAL,ETC.)
       if (yale_eos%fermi_table_x_grid(grid_idx+1).le.yale_eos%fermi_table_x_grid(grid_idx)) then
!       IF (INT(DVAL(I+1)).LE.INT(DVAL(I))) THEN
          write(short_file_unit,1000) grid_idx
 1000       format(1x,39('>'),40('<')/1x,'ERROR IN SUBROUTINE SETUPS'/ &
           1x,'GLITCH IN FERMI TABLE ELEMENT',i4/1x,'RUN STOPPED')
! 2026 (ROADMAP.md stage 3): stop converted to the ierr funnel below.
          jerr = 1
          go to 900
       endif
       bin_width = int((yale_eos%fermi_table_x_grid(grid_idx+1) - yale_eos%fermi_table_x_grid(grid_idx))*20.0d0 + 0.10d0)
       bin_end = bin_start + bin_width - 1
       if (grid_idx.eq.41)  bin_end = 261
       do iden_idx=bin_start,bin_end
          yale_eos%fermi_table_x_lookup(iden_idx) = grid_idx
  170    continue
       end do
       bin_start = bin_end + 1
  180 continue
      end do

!  CLOSE EQUATION OF STATE FILE.
      close(fermi_unit)

! MHP 5/97 ADDED OPTION FOR NEW SCV EQUATION OF STATE TABLES.
      if(use_scv_eos)then
         open(unit=scv_h_unit,file=scv_h_table_path,status='OLD')
         open(unit=scv_he_unit,file=scv_he_table_path,status='OLD')
         open(unit=scv_z_unit,file=scv_z_table_path,status='OLD')
!  READ IN EQUATION OF STATE TABLES FOR HYDROGEN AND HELIUM
         do t_idx = 1, nts
            read(scv_h_unit,1) tlogx(t_idx),nptsx(t_idx)
            read(scv_he_unit,1) scvhe_dummy_val,scvhe_dummy_npts
            read(scv_z_unit,1) scvz_dummy_val,scvz_dummy_npts
    1       format(f5.2,i4)
! TABLE GRID POINTS IN T, P(T) ARE THE SAME - NPTSY AND TLOGX
! READ IN TO RETAIN PARALLEL COMMON BLOCK STRUCTURE.
            do p_idx = 1, nptsx(t_idx)
               read(scv_h_unit,2) (tablex(t_idx,p_idx,col_idx),col_idx=1,11)
               read(scv_he_unit,2) (tabley(t_idx,p_idx,col_idx),col_idx=1,11)
    2          format(f6.2,1p2e13.5,0p,8f9.4)
            end do
!  READ IN METAL EQUATION OF STATE TABLE; COMPUTED USING THE PRATHER
! EQUATION OF STATE IN THE OLD YALE CODE.
            do p_idx = nptsx(t_idx),1,-1
               read(scv_z_unit,3)(tablez(t_idx,p_idx,col_idx),col_idx=1,13)
 3             format(f6.2,12f9.4)
            end do
         end do
      endif

      return

! error funnel: same contract as eos_get's.
  900 continue
      if (present(ierr)) then
         ierr = jerr
         return
      end if
      stop
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
     saha_state, ierr)

      use opal_eos_lib
      use const_lib
      use scv_eos_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, metal_fraction, &
           temperature_1e6k, density, pressure
      double precision, intent(out) :: gamma1, adiabatic_gradient
      integer, intent(inout) :: saha_state
! 2026 (ROADMAP.md stage 3): OPTIONAL ierr, same contract as eos_get's.
! Note the out-of-table alternate return of esac06 is NOT an error
! here: it falls through with the previous point's results, preserved
! verbatim from calcad's original handling.
      integer, intent(out), optional :: ierr

      integer :: jerr

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

      if (present(ierr)) ierr = 0
      jerr = 0

      if (use_opal2006_eos) then
         eos_interp_order = 9
         eos_rad_flag = 1
         x_local = hydrogen_fraction
         t6_local = temperature_1e6k
         d_local = density
         call esac06(x_local, t6_local, d_local, eos_interp_order, &
              eos_rad_flag, jerr, *100)
         if (jerr /= 0) go to 900
  100    continue
         gamma1 = opal_eos%eos_output_06(8)
         adiabatic_gradient = 1.0d0/opal_eos%eos_output_06(9)
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
              eos_atmosphere_flag, saha_state, jerr)
         if (jerr /= 0) go to 900
         chi_rho = 1.0d0/qdp_eos
         chi_t = -chi_rho*qdt_eos
         specific_heat_cv = qcp_eos - exp(ln10*(log10_pressure - &
              log10_density - log10_temperature))*chi_t**2/chi_rho
         gamma1 = chi_rho*qcp_eos/specific_heat_cv
         adiabatic_gradient = dela_eos
      end if

      return

! error funnel: same contract as eos_get's.
  900 continue
      if (present(ierr)) then
         ierr = jerr
         return
      end if
      stop
end subroutine eos_get_gamma1


!----------------------------------------------------------------------
! eos_set_mixture
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass): the star layer pushes the envelope
! mixture here whenever it recomputes the envelope composition
! (starin's mixture blocks, rscale's rescale); the eos internals read
! only eos_mixture_lib's eos_mix. This removed the eos domain's last
! reads of star_info.
subroutine eos_set_mixture(envelope_hydrogen_fraction, &
     envelope_metal_fraction, mean_atomic_mass, species_fractions)

      use eos_mixture_lib, only: eos_mix
      implicit none

      double precision, intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, mean_atomic_mass
      double precision, intent(in) :: species_fractions(12)

      eos_mix%envelope_hydrogen_fraction = envelope_hydrogen_fraction
      eos_mix%envelope_metal_fraction = envelope_metal_fraction
      eos_mix%amuenv = mean_atomic_mass
      eos_mix%fxenv = species_fractions

      return
end subroutine eos_set_mixture

end module eos_lib
