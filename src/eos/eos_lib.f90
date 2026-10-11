!----------------------------------------------------------------------
! eos_lib
!----------------------------------------------------------------------
! Added 2026 as part of the YREC readability refactor's phase two
! (disentangling the solver from the physics domains -- see
! docs/GUIDELINES.md's "Physics domains still entangled with the solver").
!
! eos_eval is the first facade of this phase: a single explicit-
! interface entry point that replaces the `if (use_mhd_eos) call compute_mhd_eos
! else call compute_blended_eos` dispatch that used to be duplicated at every call
! site, and additionally centralizes the Debye-Huckel composition
! setup (debye_huckel_x/y/z_total/z(3), read only by compute_fully_ionized_eos
! on the non-MHD path) that 6 of those call sites also used to
! duplicate verbatim. compute_blended_eos (which also hosts blend_eos, its
! co-located pair), compute_fully_ionized_eos, and compute_mhd_eos are all
! unchanged; this is a pure dispatch/setup wrapper.
!
! composition_at_zone is OPTIONAL: callers that carry full per-species
! composition (misc/henyey_coefficients.f90, io/write_legacy_output.f90, core/read_starting_model.f90,
! misc/shell_physics.f90, mixing/compute_scale_height.f90, mixing/semiconvection.f90) pass
! composition(:,idx) and get the same Debye-Huckel setup they used to
! compute themselves. Callers that only ever carried bulk
! hydrogen_fraction/metal_fraction (atm/atm_lib.f90, atm/atmosphere_derivs.f90,
! atm/envelope_derivs.f90, wind/massloss.f90) omit it -- exactly matching their
! prior behavior, since none of them set these fields before either.
!
! core/read_starting_model.f90 previously had a bug here (confirmed against the
! original F77 source via git history): a missing ELSE meant it called
! *both* compute_mhd_eos and compute_blended_eos when MHD was on, and *neither* when MHD was
! off. eos_eval's if/else has the structurally-correct form, so
! migrating that call site fixes the bug by construction.
!
! mixing/compute_scale_height.f90, mixing/semiconvection.f90, and wind/massloss.f90
! previously called compute_blended_eos unconditionally, with no LMHD check at all
! -- confirmed authentic original YREC behavior (unchanged since the
! very first commit, not a modernization artifact). Migrating them to
! eos_eval extends real MHD support to these three secondary/diagnostic
! calculations for the first time, per explicit user sign-off; this is
! an acknowledged numerics change for use_mhd_eos=.true. runs, which
! the Stage-0 regression suite cannot verify (no test case sets LMHD).
!
! eos_get_gamma1 below was created for the acoustic-depth diagnostic
! (calcad), which has since been retired (commit 55c372b); its only
! remaining caller is eos/test/test_eos.f90.
module eos_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opal_eos_tables, only: i_opal_gamma1, i_opal_gamma2_ratio, opal06, &
           opal_eos
      use scv_eos_tables, only: nptsx, scv_nt, tablex, tabley, tablez, tlogx
      use yale_eos_tables, only: debye_huckel_coefficient, debye_huckel_nu, &
           debye_huckel_x, debye_huckel_y, debye_huckel_z, &
           debye_huckel_z_total, yale
      use physical_constants, only: cpi, ln10
      use io_units, only: run_log_unit
      use eos_mixture, only: eos_mix, n_mix_species
      use eos_configuration, only: eos_settings, eos_settings_t, &
           set_eos_settings
      use math_lib
! 2026 named-index result array (ROADMAP "Named-index result arrays"):
! eos_get packs eos_eval's 24 thermodynamic outputs into one
! res(num_eos_results) array indexed by the constants below (MESA's
! eosDT_get shape). The independent variables (log10 T, log10 P),
! composition, mode flags, and the saha state stay explicit
! arguments. res is intent(inout): the historically-inout slots
! (i_log10_density, i_beta, the ionization fractions, i_mu_ion_inv,
! i_dlnrho_dlnt/p, i_cp, i_grada) carry their previous-call values in
! exactly as the old per-caller locals did; eos_eval itself is called
! unchanged, so results are byte-identical.
      use blended_eos, only: compute_blended_eos, blend_eos
      use mhd_eos, only: compute_mhd_eos, read_mhd_tables
      use opal2006_eos, only: esac06
      use scv_eos, only: build_scv_envelope_table
      implicit none
! Everything is private unless exported below (2026: eos_get is the
! single public query -- the named-index result-array form; the
! long-argument engine underneath it, eos_eval, is exported ONLY for
! the domain's own standalone test and stays outside the
! check_boundaries.py cross-domain allowlist). The i_*/num_eos_results
! index parameters carry their own public attributes.
! eos_get_gamma1 is exported for eos/test/test_eos.f90 only (no
! production caller; see the header note above).
      private
      public :: eos_get, eos_get_gamma1, eos_init, eos_set_mixture, eos_set_debye_huckel_z
      public :: eos_eval
! the settings record (eos_configuration) and its setter; the per-zone
! Debye-Hueckel composition; the envelope-mixture SCV table
      public :: eos_settings_t, set_eos_settings, eos_set_debye_huckel_composition, &
           eos_update_envelope_tables
! result-array slots for eos_get
      integer, parameter, public :: &
           i_temperature = 1, i_pressure = 2, i_log10_density = 3, &
           i_density = 4, i_beta = 5, i_beta_inverse = 6, i_beta14 = 7, &
           i_fxion = 8, &          ! three slots, 8-10: ion_fraction(1:3)
           i_gas_constant = 11, i_mu_ion_inv = 12, i_mu_e_inv = 13, &
           i_eta = 14, i_dlnrho_dlnt = 15, i_dlnrho_dlnp = 16, &
           i_cp = 17, i_grada = 18, i_dlnrho_dlnt_dt = 19, &
           i_dlnrho_dlnp_dt = 20, i_grada_dt = 21, i_grada_dp = 22, &
           i_cp_dt = 23, i_cp_dp = 24
      integer, parameter, public :: num_eos_results = 24

contains

! ---------------------------------------------------------------
! Named-index result-array form of eos_eval (see module header).
subroutine eos_get(log10_temperature, log10_pressure, &
     hydrogen_fraction, metal_fraction, res, want_derivatives, &
     in_atmosphere, saha_state, composition_at_zone, ierr)
      real(dp), intent(inout) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(inout) :: res(:)
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
      real(dp), intent(in), optional :: composition_at_zone(:)
      integer, intent(out) :: ierr
      real(dp) :: fxion_local(3)

! The three ionization fractions occupy res(i_fxion:i_fxion+2); they
! go through a contiguous local because eos_eval takes them as one
! 3-element inout array. An absent composition_at_zone is passed on as
! absent (an optional actual may be an absent optional dummy).
      fxion_local = res(i_fxion:i_fxion+2)
      call eos_eval(log10_temperature, res(i_temperature), &
           log10_pressure, res(i_pressure), res(i_log10_density), &
           res(i_density), hydrogen_fraction, metal_fraction, &
           res(i_beta), res(i_beta_inverse), res(i_beta14), &
           fxion_local, res(i_gas_constant), res(i_mu_ion_inv), &
           res(i_mu_e_inv), res(i_eta), res(i_dlnrho_dlnt), &
           res(i_dlnrho_dlnp), res(i_cp), res(i_grada), &
           res(i_dlnrho_dlnt_dt), res(i_dlnrho_dlnp_dt), &
           res(i_grada_dt), res(i_grada_dp), res(i_cp_dt), &
           res(i_cp_dp), want_derivatives, in_atmosphere, saha_state, &
           composition_at_zone=composition_at_zone, ierr=ierr)
      res(i_fxion:i_fxion+2) = fxion_local
end subroutine eos_get


subroutine eos_eval(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
     in_atmosphere, saha_state, composition_at_zone, ierr)

      real(dp), intent(inout) :: log10_temperature
      real(dp), intent(out) :: temperature
      real(dp), intent(inout) :: log10_pressure
      real(dp), intent(out) :: pressure
      real(dp), intent(inout) :: log10_density
      real(dp), intent(out) :: density
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(inout) :: beta
      real(dp), intent(out) :: beta_inverse, beta14
      real(dp), intent(inout) :: ion_fraction(:)
      real(dp), intent(out) :: specific_gas_constant
      real(dp), intent(inout) :: ion_mean_weight_inverse
      real(dp), intent(out) :: electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      real(dp), intent(inout) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      real(dp), intent(out) :: dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
      real(dp), intent(in), optional :: composition_at_zone(:)
! 2026 (docs/ROADMAP.md stage 3): ierr follows MESA's ierr-not-stop
! discipline (same shape as kap_lib's kap_get): any error condition
! that used to stop deep inside the eos internals is returned here
! instead, ierr /= 0, with the diagnostic already printed at the point
! of failure.
      integer, intent(out) :: ierr

      integer :: jerr

      ierr = 0
      jerr = 0

      if (eos_settings%use_mhd_eos) then
         call compute_mhd_eos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      else
         if (present(composition_at_zone) .and. eos_settings%use_debye_huckel_correction) then
            debye_huckel_x = composition_at_zone(1)
            debye_huckel_y = composition_at_zone(2)+composition_at_zone(4)
            debye_huckel_z_total = composition_at_zone(3)
            debye_huckel_z(1) = composition_at_zone(5)+composition_at_zone(6)
            debye_huckel_z(2) = composition_at_zone(7)+composition_at_zone(8)
            debye_huckel_z(3) = composition_at_zone(9)+composition_at_zone(10)+ &
                 composition_at_zone(11)
         end if
         call compute_blended_eos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
              in_atmosphere, saha_state, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      end if

end subroutine eos_eval

!----------------------------------------------------------------------
! eos_init
!----------------------------------------------------------------------
! Added 2026 (phase three, docs/ROADMAP.md stage 1): the eos domain's
! startup-time table-load lifecycle entry, following MESA's
! <mod>_init convention. Covers the eos tables loaded from named
! files at startup: the MHD EOS tables (gated on use_mhd_eos, moved
! inside from the call site so the caller needs no knowledge of
! which eos configuration requires which tables -- a no-op when MHD
! is off, exactly as the gated call was) and the Yale EOS's
! Fermi-Dirac degenerate-electron table (read unconditionally, into
! yale_eos_tables; this block previously lived inline in
! setup/setups.f90, writing into atm_tables -- moved here
! verbatim, including its bin-construction loop, its glitch check,
! and that check's `stop`, which docs/ROADMAP.md stage 3 will convert to
! ierr). The OPAL 1995/2001/2006 EOS tables are not loaded here:
! they are lazily read on first use inside esac/esac01/esac06 via
! their table_loaded_flag guards. The SCV (Saumon-Chabrier-Van Horn)
! EOS tables are also read here (gated on use_scv_eos, into
! scv_eos_tables's state consumed by eos/scv/compute_scv_envelope_eos+compute_scv_eos and the io
! writers) -- another inline setups.f90 block moved in verbatim; its
! read is provably independent of the atm tables read later in
! startup (separate file units, disjoint state), so absorbing it here
! is order-safe.
subroutine eos_init(settings, fermi_table_path, scv_h_table_path, &
     scv_he_table_path, scv_z_table_path, zams_a_table_path, &
     zams_b_table_path, zams_c_table_path, centre1_table_path, &
     centre2_table_path, centre3_table_path, centre4_table_path, &
     centre5_table_path, opal_eos_table_path, ierr)

      integer, parameter :: nts = scv_nt

      character(len=256), intent(in) :: fermi_table_path
      character(len=256), intent(in) :: scv_h_table_path, &
           scv_he_table_path, scv_z_table_path
      character(len=256), intent(in) :: zams_a_table_path, &
           zams_b_table_path, zams_c_table_path, centre1_table_path, &
           centre2_table_path, centre3_table_path, centre4_table_path, &
           centre5_table_path
! the run's EOS settings (installed here); the OPAL table the settings
! select (opened here, read lazily by its vintage's reader)
      type(eos_settings_t), intent(in) :: settings
      character(len=*), intent(in) :: opal_eos_table_path
! 2026 (docs/ROADMAP.md stage 3): ierr, same contract as eos_eval's.
      integer, intent(out) :: ierr

      integer :: jerr
      integer :: grid_idx, iden_idx, col_idx, card_idx
      integer :: bin_start, bin_width, bin_end
      integer :: t_idx, p_idx
      real(dp) :: scvhe_dummy_val, scvz_dummy_val
      integer :: scvhe_dummy_npts, scvz_dummy_npts
      integer :: fermi_unit, scv_h_unit, scv_he_unit, scv_z_unit
      logical :: is_open
! the Debye-Hueckel coefficient's constants (cgs; formerly computed in
! setups)
      real(dp), parameter :: boltzmann_constant = 1.380658d-16, &
           electron_charge_esu = 4.802d-10
      real(dp) :: hydrogen_atom_mass

      ierr = 0
      jerr = 0
      call set_eos_settings(settings)
! the Debye-Hueckel coefficient and the charges of the 18 metals it
! weights (C, N, O, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, Ca, Ti, Cr, Mn, Fe, Ni)
      hydrogen_atom_mass = 1.0d0/6.0222137d23
      debye_huckel_coefficient = -sqrt(cpi/(boltzmann_constant*hydrogen_atom_mass* &
           hydrogen_atom_mass*hydrogen_atom_mass))*electron_charge_esu* &
           electron_charge_esu*electron_charge_esu/3.0d0
      debye_huckel_nu = [6.0d0, 7.0d0, 8.0d0, 10.0d0, 11.0d0, 12.0d0, 13.0d0, &
           14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0, 20.0d0, 22.0d0, 24.0d0, &
           25.0d0, 26.0d0, 28.0d0]
! the OPAL table, opened for its reader
      if (settings%use_opal95_eos .or. settings%use_opal2001_eos .or. &
           settings%use_opal2006_eos) then
         if (opal_eos%table_unit /= 0) then
            inquire(unit=opal_eos%table_unit, opened=is_open)
            if (is_open) close(opal_eos%table_unit)
         end if
         open(newunit=opal_eos%table_unit, file=opal_eos_table_path, status='OLD', &
              iostat=jerr)
         if (jerr /= 0) then
            write(run_log_unit,*) 'eos_init: cannot open the OPAL EOS table ', &
                 trim(opal_eos_table_path)
            ierr = 1
            return
         end if
      end if
      if (settings%use_mhd_eos) then
         call read_mhd_tables(zams_a_table_path, zams_b_table_path, &
              zams_c_table_path, centre1_table_path, centre2_table_path, &
              centre3_table_path, centre4_table_path, centre5_table_path, &
              jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      end if

!  INPUT F-TABLES FOR DEGENERATE EQUATION OF STATE
!  OPEN DATA FILES
      open(newunit=fermi_unit,file=fermi_table_path,form='FORMATTED',status='OLD')
      read(fermi_unit,150)  (yale%fermi_table_x_grid(grid_idx), grid_idx = 1,43)
      read(fermi_unit,150)  (yale%fermi_table_eta(grid_idx), grid_idx = 1,43)
      read(fermi_unit,160) (grid_idx,iden_idx,(yale%fermi_table_data(col_idx,grid_idx,iden_idx), &
           col_idx=1,5),card_idx=1,860)
  150 format(8x,9f8.4)
  160 format(2i5,5f12.7)
      rewind fermi_unit
      bin_start = 1
      do grid_idx=1,41
! MHP 10/02 SHOULD BE INT(DVAL,ETC.)
       if (yale%fermi_table_x_grid(grid_idx+1) <= yale%fermi_table_x_grid(grid_idx)) then
!       IF (INT(DVAL(I+1)).LE.INT(DVAL(I))) THEN
          write(run_log_unit,1000) grid_idx
 1000       format(1x,39('>'),40('<')/1x,'ERROR IN SUBROUTINE SETUPS'/ &
           1x,'GLITCH IN FERMI TABLE ELEMENT',i4/1x,'RUN STOPPED')
! 2026 (docs/ROADMAP.md stage 3): stop converted to the ierr funnel below.
          jerr = 1
          ierr = jerr
          return
       endif
       bin_width = int((yale%fermi_table_x_grid(grid_idx+1) - yale%fermi_table_x_grid(grid_idx))*20.0d0 + 0.10d0)
       bin_end = bin_start + bin_width - 1
       if (grid_idx == 41)  bin_end = 261
       do iden_idx=bin_start,bin_end
          yale%fermi_table_x_lookup(iden_idx) = grid_idx
       end do
       bin_start = bin_end + 1
      end do

!  CLOSE EQUATION OF STATE FILE.
      close(fermi_unit)

! MHP 5/97 ADDED OPTION FOR NEW SCV EQUATION OF STATE TABLES.
      if(settings%use_scv_eos)then
         open(newunit=scv_h_unit,file=scv_h_table_path,status='OLD')
         open(newunit=scv_he_unit,file=scv_he_table_path,status='OLD')
         open(newunit=scv_z_unit,file=scv_z_table_path,status='OLD')
!  READ IN EQUATION OF STATE TABLES FOR HYDROGEN AND HELIUM
         do t_idx = 1, nts
            read(scv_h_unit,1) tlogx(t_idx),nptsx(t_idx)
            read(scv_he_unit,1) scvhe_dummy_val,scvhe_dummy_npts
            read(scv_z_unit,1) scvz_dummy_val,scvz_dummy_npts
    1       format(f5.2,i4)
! The He and Z tables share the H table's T grid and P(T) point
! counts; their copies are read into dummies.
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

end subroutine eos_init

!----------------------------------------------------------------------
! eos_get_gamma1
!----------------------------------------------------------------------
! Added 2026 (phase three, docs/ROADMAP.md stage 1): a component-accessor
! entry in the spirit of MESA's eosDT_get_component, created for the
! acoustic-depth diagnostic (calcad). That diagnostic was retired in
! commit 55c372b, so this entry now has NO production caller: it is
! exercised only by eos/test/test_eos.f90 and is kept as a pinned
! test surface. Given composition and an (unlogged) T/rho/P point,
! returns gamma1 and the adiabatic gradient, dispatching on
! use_opal2006_eos exactly as calcad historically did: the OPAL 2006
! table (via esac06, density-basis, radiation correction on, 9
! variables; a failed table lookup falls through with whatever the
! previous point left behind, as calcad did) when enabled, otherwise the
! Yale/SCV path (via blend_eos, pressure-basis) with gamma1 built from
! chi_rho/chi_T/cv per Cox & Giuli ch. 9.
subroutine eos_get_gamma1(hydrogen_fraction, metal_fraction, &
     temperature_1e6k, density, pressure, gamma1, adiabatic_gradient, &
     saha_state, ierr)

      real(dp), intent(in) :: hydrogen_fraction, metal_fraction, &
           temperature_1e6k, density, pressure
      real(dp), intent(out) :: gamma1, adiabatic_gradient
      integer, intent(inout) :: saha_state
! 2026 (docs/ROADMAP.md stage 3): ierr, same contract as eos_eval's.
! Note an out-of-table esac06 lookup (in_table false) is NOT an
! error here: it falls through with the previous point's results, preserved
! verbatim from calcad's original handling.
      integer, intent(out) :: ierr

      integer :: jerr

! --- locals (per-call copies: esac06/blend_eos write into several of
!     the positions the inputs occupy) ---
      integer :: eos_interp_order, eos_rad_flag
      real(dp) :: x_local, t6_local, d_local
      real(dp) :: temperature_local, log10_temperature, &
           log10_pressure, pressure_local, log10_density
      real(dp) :: beta_local, beta_inverse_local, beta14_local, &
           ion_fraction_local(3), specific_gas_constant_local, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      real(dp) :: qdt_eos, qdp_eos, qcp_eos, dela_eos, qdtt_eos, &
           qdtp_eos, qat_eos, qap_eos, qcpt_eos, qcpp_eos
      real(dp) :: chi_rho, chi_t, specific_heat_cv
      logical :: eos_deriv_flag, eos_atmosphere_flag, in_table

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
         x_local = hydrogen_fraction
         temperature_local = temperature_1e6k*1.0d6
         d_local = density
         log10_density = log10(density)
         log10_temperature = log10(temperature_local)
         pressure_local = pressure
         log10_pressure = log10(pressure_local)
         eos_deriv_flag = .false.
         eos_atmosphere_flag = .true.
         call blend_eos(log10_temperature, temperature_local, &
              log10_pressure, pressure_local, log10_density, d_local, &
              x_local, metal_fraction, beta_local, beta_inverse_local, &
              beta14_local, ion_fraction_local, &
              specific_gas_constant_local, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              qdt_eos, qdp_eos, qcp_eos, dela_eos, qdtt_eos, qdtp_eos, &
              qat_eos, qap_eos, qcpt_eos, qcpp_eos, eos_deriv_flag, &
              eos_atmosphere_flag, saha_state, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         chi_rho = 1.0d0/qdp_eos
         chi_t = -chi_rho*qdt_eos
         specific_heat_cv = qcp_eos - exp(ln10*(log10_pressure - &
              log10_density - log10_temperature))*chi_t**2/chi_rho
         gamma1 = chi_rho*qcp_eos/specific_heat_cv
         adiabatic_gradient = dela_eos
      end if

end subroutine eos_get_gamma1


!----------------------------------------------------------------------
! eos_set_mixture
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass): the star layer pushes the envelope
! mixture here whenever it recomputes the envelope composition
! (starin's mixture blocks, rscale's rescale); the eos internals read
! only eos_mixture's eos_mix. This removed the eos domain's last
! reads of star_info.
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
! Added 2026 (wave 2): the LAOL89 opacity table header carries the
! 18-element metal mixture (C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,Ca,Ti,Cr,
! Mn,Fe,Ni, scaled to sum to ZHIT) that the Debye-Huckel correction
! weights with debye_huckel_nu; kap/kap_lib.f90 pushes it here
! right after read_laol_tables reads it (formerly read_laol_tables wrote yale_eos_tables's
! debye_huckel_z directly). Elements 1-3 are recomputed per shell by
! the star layer (core/henyey_iteration.f90 etc.); 4-18 come only from
! here.
subroutine eos_set_debye_huckel_z(metal_mixture)
      real(dp), intent(in) :: metal_mixture(:)

      debye_huckel_z = metal_mixture

end subroutine eos_set_debye_huckel_z

!----------------------------------------------------------------------
! eos_set_debye_huckel_composition
!----------------------------------------------------------------------
! The composition the Debye-Hueckel correction uses on calls that pass
! none (composition_at_zone absent), from a zone's composition (star%xa
! layout): the star layer sets the surface zone's before the envelope
! integration. Does nothing when the correction is off. (2026-10: the
! three callers wrote yale_eos_tables' variables themselves.)
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
