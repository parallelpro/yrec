!----------------------------------------------------------------------
! parmin
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original parmin.f; only source form was updated -- see the naming
! note below for why variable names are treated differently here.
!
! Reads the CONTROL and PHYSICS namelists (yrec8.nml1/yrec8.nml2, or
! the two files named on the command line) that configure an entire
! YREC run: envelope/atmosphere physics switches, opacity/EOS table
! selection, diffusion and mixing options, rotation and wind-law
! parameters, nuclear cross sections, and the sequence of "kind" run
! cards (evolve / rescale / rescale+evolve) executed by the driver.
! Also opens every logical unit used for the run's input/output files
! and writes a human-readable summary of the active run parameters to
! the .short log.
!
! NAMING NOTE: unlike other files in this modernization project, the
! members of every COMMON block below, and every variable listed in
! NAMELIST /control/ or /physics/, keep their EXACT original spelling
! (only lowercased). Fortran's namelist reader matches input-file
! entries (e.g. "LCALS = .TRUE." in a *.nml1 file) to variables by
! name, against every *.nml1 and *.nml2 file in the repository -- so
! renaming any of these would silently break every existing input
! file that sets it. This also means every dummy argument of this
! particular subroutine happens to be a NAMELIST /control/ member
! too (the namelist read sets these arguments' values directly, which
! are then forwarded to other routines), so none of them could be
! renamed either. Only true local variables (never referenced by
! COMMON, NAMELIST, or the argument list) were renamed for
! readability; see also the two SUBROUTINE EXPAND_VALUE locals below,
! a fully-local scope with no such restriction.
subroutine parmin(falex06, fallard, fatm, ffermi, fkur, fkur2, flaol, &
     flaol2, fliv95, flldat, fmhd1, fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, &
     fmhd7, fmhd8, fopal2, fpatm, fpenv, fpmod, fpurez, fscvh, fscvhe, &
     fscvz, opecalex)

      use envelope_comp_lib
      use const_lib
      use luout_lib
      use intpar_lib
      implicit none

! PARAMETERS for Allard model surface pressures (n_allard_teff/
! n_allard_logg) and the shared array length used by the variable-FC
! and acoustic-depth diagnostics (max_diag_pts). n_atm_teff/n_atm_logg/
! n_katm_teff/n_katm_logg/n_scv_teff/n_scv_press (former parameters for
! the tabulated Kurucz/Castelli surface pressures and SCV EOS tables,
! now handled via const_lib's own array dimensions) are no longer used
! in this file.
      integer, parameter :: n_allard_teff = 54, n_allard_logg = 5
      integer, parameter :: max_diag_pts = 5000

! --- CONTROL/PHYSICS namelist variables (including this routine's
!     dummy arguments) and their explicit type overrides (CHARACTER
!     length / array DIMENSION); kept at original spelling -- see
!     NAMING NOTE above. Grouped as in the original file; a few
!     members declared here (OLAOL/OXA/OT/ORHO/TOLLAOL, YRECVER,
!     GITHASH, DESCRIP, FMONTE1/FMONTE2, ICLCD/IACAT/IJLAST/IJVS/
!     IJENT/IJDEL, AWIND) are also COMMON members and get their type
!     declared once, below, alongside their COMMON block instead. ---
      character(len=256) :: version_fmt
      character(len=256) :: flaol, fpurez
      character(len=256) :: flaol2, fopal2, fkur2
      character(len=3) :: anewcp, atmp
      character(len=8) :: amix, aiso
      character(len=256) :: fiso
      character(len=256) :: fatm
      character(len=256) :: fstch
      character(len=256) :: fallard, fscvh, fscvhe, fscvz
      character(len=256) :: flast, ffirst, ffermi, &
           fdebug, ftrack, fshort, fmilne, fmodpt, &
           fstor, fpmod, fpenv, fpatm, fdyn, &
           flldat, fsnu, fscomp, fkur, &
           fmhd1, fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, fmhd7, fmhd8
      integer :: kindrn(50)
      double precision :: rsclm(50), rsclx(50), rsclz(50), rsclcm(50)
      character(len=256) :: fliv95
      character(len=256) :: fopale, fopale01, fcondopacp, fopale06
      character(len=256) :: opecalex(7)
      character(len=256) :: falex06

! --- NAMELIST-only variables that had no explicit declaration in the
!     original (relying on IMPLICIT typing); kept at original
!     spelling, typed per the implicit rule their first letter implies
!     (A-H,O-Z => double precision; L => logical) ---
      double precision :: alfa, fk, wmax_sun, pmmwmax, zalex2, zopal952
      logical :: lpmm, lsolwind

! --- true local variables (never referenced by COMMON, NAMELIST, or
!     the argument list): freely renamed for readability ---
      character(len=3) :: element_id(12)
      character(len=8) :: mixture_id_table(4)
      double precision :: zx_mix_table(4), frac_c_table(4), &
           frac_n_table(4), frac_o_table(4)
      character(len=256) :: fcalcad
      character(len=256) :: control_nml_file, physics_nml_file
      character(len=256) :: shell_cmd
      integer :: i, j, last_slash_idx
      integer :: short_prefix_len
      double precision :: one_third, two_thirds
! parmin_ln10: this file's own private ln(10) (never read after being
! set -- see the assignment below), distinct from const_lib's ln10
! (which this file now also uses via const3's `use const_lib`, hence
! the rename needed here to avoid a name collision).
      double precision :: parmin_ln10
      double precision :: sum_frac
      integer :: nkind
      integer :: first_model_binary_lu, last_model_binary_lu, &
           stored_models_binary_lu

! former common/vnewcb/: vnew is a NAMELIST /physics/ value spelled
! identically to its const_lib canonical name -- use-associated
! directly rather than locally declared.

! former common/lunum/: none of these 12 members are namelist values
! (they're hardcoded unit numbers this file assigns unconditionally
! further down), so they're simply renamed in place to their canonical
! const_lib names (first_unit/run_unit/standard_unit/fermi_unit/
! opal_model_unit/opal_envelope_unit/opal_atm_unit/dynamics_unit/
! laol_table_unit/neutrino_unit/composition_unit/kurucz_table_unit),
! now use-associated rather than locally declared.

! fmonte1/fmonte2: NAMELIST /physics/ members, each with a different
! canonical const_lib spelling (monte_carlo_file1_path/
! monte_carlo_file2_path), so kept local under their NAMELIST spelling
! here and copy-assigned after the namelist read below. imonte1/imonte2
! (former common/iomonte/'s remaining members) are not namelist values
! -- this file assigns them fixed unit numbers unconditionally further
! down -- so they're simply renamed in place to their canonical
! const_lib names (monte_carlo_unit1/monte_carlo_unit2), now
! use-associated rather than locally declared.
      character(len=256) :: fmonte1, fmonte2

! descrip is declared only in this file (no other file shares
! common/desc/), so it's simply a plain local rather than needing any
! module treatment.
      character(len=256) :: descrip(2)

! former common/ccout/: lstore/lstatm/lstenv/lstmod/lstphys/lstrot/
! lscrib/lstch are NAMELIST /physics/ values, spelled identically to
! their const_lib canonical names, so they're use-associated directly
! (no separate local copy or rename needed -- the NAMELIST statement
! below can bind straight to the use-associated module variable, since
! its spelling matches the .nml2 files exactly). lphhd (former
! common/ccout/'s remaining member) is not a namelist value -- this
! file sets it directly further down -- also now use-associated.

! npenv/nprtmod/npoint: NAMELIST /physics/ members, spelled identically
! to their const_lib canonical names -- use-associated directly, same
! reasoning as lstore etc above. nprtpt is also a NAMELIST member but
! needs a different canonical spelling (print_point_interval), so it
! keeps its local NAMELIST-spelled name here and is copy-assigned
! after the namelist read below, per this file's usual pattern.
      integer :: nprtpt

! former common/pulsegyre/: new (2026) dedicated block for the
! GYRE-format periodic pulsation-structure output feature
! (io/write_gyre_pulse.f90, triggered from io/wrtout.f90 every
! pulse_gyre_interval converged models). pulse_gyre_interval is a
! NAMELIST /control/ value spelled identically to its const_lib
! canonical name -- use-associated directly.

! former common/ccout2/: ldebug/lcorr/lmilne/ltrack/lstpch are
! NAMELIST /physics/ values, spelled identically to their const_lib
! canonical names -- use-associated directly, same reasoning as
! lstore etc above.

! lnew0: NAMELIST /physics/ member spelled identically to its
! const_lib canonical name -- use-associated directly. tridt/tridl are
! also NAMELIST members but need different canonical spellings
! (tri_delta_teffl/tri_delta_logl), so they keep their local
! NAMELIST-spelled names here and are copy-assigned after the
! namelist read below. senv0/lsenv0 (former common/cenv/'s remaining
! members, renamed to requested_envelope_mass/
! change_envelope_mass_flag) are not namelist values and unused in
! this file -- core/starin.f90 computes them -- so they're dropped
! entirely.
      double precision :: tridt, tridl

! rescal/iresca are not namelist values -- computed elsewhere in this
! file (rescal from rsclm/rsclx/rsclz/rsclcm, iresca from kindrn), so
! renamed in place to their canonical const_lib names (rescale_params/
! rescale_kind), now use-associated. nmodls/lfirst/numrun (former
! common/ckind/'s remaining members) are NAMELIST values with
! different canonical spellings (num_models/first_call_flag/
! num_runs), kept local under their NAMELIST spelling here and
! copy-assigned after the namelist read below.
      integer :: nmodls(50), numrun
      logical :: lfirst(50)


! clsun/crsun: NAMELIST /physics/ members (must keep this exact
! spelling); copied into const_lib's solar_luminosity_cgs/
! solar_radius_cgs after the namelist read below. The other six former
! common/const/ members (clsunl/clnsun/cmsun/cmsunl/crsunl/cmbol) are
! unused in this file -- setup/setups.f90 computes them from
! solar_luminosity_cgs/solar_radius_cgs at startup -- so they're
! dropped entirely rather than carried forward.
      double precision :: clsun, crsun



! dtwind: NAMELIST /physics/ member with a different canonical
! spelling (max_domega_global), so kept local under its NAMELIST
! spelling here and copy-assigned after the namelist read below.
      double precision :: dtwind

! lptime: NAMELIST /physics/ member with a different canonical
! spelling (use_structure_dt_limits), so kept local under its
! NAMELIST spelling here and copy-assigned after the namelist read
! below.
      logical :: lptime

! htoler/fcorr0/fcorri/niter1/niter2/niter3: NAMELIST /physics/ members
! spelled identically to their const_lib canonical names -- use-
! associated directly. hpttol is also a NAMELIST member (kept local
! purely so a .nml2 file could still set HPTTOL directly by that name)
! but is NOT copy-assigned into const_lib's chi_grid_scale the usual
! way: when lnewvars is set (as in every current example), setup/
! remap.f90 recomputes chi_grid_scale from the more user-friendly
! tol_dm_min/tol_dm_max/etc namelist inputs when `call remap` runs
! below, overwriting whatever hpttol/chi_grid_scale started as -- so
! copying hpttol into chi_grid_scale here would just get clobbered
! (this file's own diagnostic writes read chi_grid_scale directly,
! after `call remap`, not hpttol, to reflect that). fcorr (former
! common/ctol/'s remaining member) is unused in this file --
! core/crrect.f90/core/main.f90 compute it at runtime -- so it's
! dropped entirely.
      double precision :: hpttol(12)

! dtdif/itdif1/itdif2: NAMELIST /physics/ members spelled identically
! to their const_lib canonical names -- use-associated directly. djok
! is also a NAMELIST member but needs a different canonical spelling
! (convergence_tolerance), so it keeps its local NAMELIST-spelled name
! here and is copy-assigned after the namelist read below.
      double precision :: djok

! lovste: NAMELIST /physics/ member (must keep this exact spelling).
! Former common /dpmix/; every other member already matches
! const_lib's canonical spelling, so only lovste needs to stay local
! (copied into const_lib's envelope_overshoot_active right after the
! namelist read below).
      logical :: lovste

! atmstp/envstp/lenvg: NAMELIST /physics/ members, each with a
! different canonical const_lib spelling (atm_step_size/
! envelope_step_size/envelope_generation_flag), so kept local under
! their NAMELIST spelling here and copy-assigned after the namelist
! read below.
      double precision :: atmstp, envstp
      logical :: lenvg

! lexcom: NAMELIST /physics/ member (must keep this exact spelling);
! copied into const_lib's use_extended_composition after the namelist
! read below.
      logical :: lexcom

! lkuthe: NAMELIST /physics/ member with a different canonical
! spelling (helium_flash_active), so kept local under its NAMELIST
! spelling here and copy-assigned after the namelist read below.
      logical :: lkuthe

! atmerr/atmd0/atmbeg/atmmin/atmmax: NAMELIST /physics/ members, each
! with a different canonical const_lib spelling (atm_error_tol/
! atm_step_initial/atm_step_begin/atm_step_min/atm_step_max), so kept
! local under their NAMELIST spelling here and copy-assigned after the
! namelist read below.
      double precision :: atmerr, atmd0, atmbeg, atmmin, atmmax

! enverr/envbeg/envmin/envmax: NAMELIST /physics/ members, each with a
! different canonical const_lib spelling (env_error_tol/env_step_begin/
! env_step_min/env_step_max), so kept local under their NAMELIST
! spelling here and copy-assigned after the namelist read below.
      double precision :: enverr, envbeg, envmin, envmax

! stolr0/imax/nuse: NAMELIST /physics/ members (must keep this exact
! spelling, see this file's naming note at the top). Former common
! /intpar/; copied into intpar_lib's canonically-named variables after
! the namelist read below, since intpar_lib's names differ from these
! namelist-fixed ones.
      double precision :: stolr0
      integer :: imax, nuse

! tscut: NAMELIST /physics/ member (must keep this exact spelling).
! Former common/ctlim/ member alongside atime/tcut/tenv0/tenv1/tgcut;
! those five kept their spelling when ctlim moved to const_lib, but
! tscut did not (const_lib calls it saha_log10t_cutoff), so it stays
! local here and is copied into that canonical name after the
! namelist read below, same treatment as stolr0/imax/nuse above.
      double precision :: tscut

! former common/label/: xenv0/zenv0 are not namelist values, and
! genuinely used in this file (initial envelope H/Z, read in the
! run-parameters dump below) -- renamed in place to their canonical
! const_lib names (initial_envelope_x/initial_envelope_z), now
! use-associated rather than locally declared.

! xnewcp/lnewcp: NAMELIST /physics/ members, each with a different
! canonical const_lib spelling (new_species_value/
! rescale_species_active), so kept local under their NAMELIST spelling
! here and copy-assigned after the namelist read below (lnewcp is
! also overridden further down in this file, in the CNO-rescaling
! validation block, so it gets a second copy-assignment there too).
! inewcp/lrel (former common/newcmp/'s remaining members) are not
! namelist values and are genuinely used in this file, so renamed in
! place to their canonical const_lib names (new_species_index/
! value_relative_to_h), now use-associated rather than locally
! declared.
      double precision :: xnewcp
      logical :: lnewcp

! isetmix/isetiso/r12_13/r14_15/r16_17/r16_18/xh2_ini/xhe3_ini/
! xli6_ini/xli7_ini/xbe9_ini/xb10_ini/xb11_ini: NAMELIST /physics/
! members, each with a different canonical const_lib spelling, kept
! local under their NAMELIST spelling here and copy-assigned after the
! namelist read below. frac_c/frac_n/frac_o/zxmix are likewise
! NAMELIST members with different canonical spellings, but are also
! overridden further down in this file (the CNO-mixture validation
! block), so they get a second copy-assignment there too. lmixture/
! lisotope (former common/newmx/'s remaining members) are not
! namelist values -- set by this file's own CNO/isotope validation
! logic -- so renamed in place to their canonical const_lib names
! (change_cno_mixture_active/change_isotope_ratios_active), now
! use-associated rather than locally declared.
      integer :: isetmix, isetiso
      double precision :: frac_c, frac_n, frac_o, r12_13, r14_15, r16_17, r16_18, zxmix, xh2_ini, &
           xhe3_ini, xli6_ini, xli7_ini, xbe9_ini, xb10_ini, xb11_ini

! optol: NAMELIST /physics/ member (must keep this exact spelling);
! copied into const_lib's metal_fraction_match_tolerance after the
! namelist read below. zsi (former common/optab/'s remaining
! non-dead member) is not a namelist value -- its default moved to
! const_lib.f90 (see that file's header note) since it's now
! use-associated rather than locally declared. idt/idd are unused in
! this file, so dropped entirely.
      double precision :: optol

! lrot/linstb: NAMELIST /physics/ members (must keep this exact
! spelling). Former common /rot/; every other member already matches
! const_lib's canonical spelling, so only these two need to stay local
! (copied into const_lib's rotation_active/instability_transport_active
! right after the namelist read below).
      logical :: lrot, linstb

! endage/setdt/end_dcen/end_xcen/end_ycen: NAMELIST /physics/ members
! (must keep this exact spelling); copied into const_lib's
! target_end_age/timestep_override/central_deuterium_stop/
! central_hydrogen_stop/central_helium_stop after the namelist read
! below. lendag/lsetdt (former common/sett/'s remaining two members)
! are not namelist values -- they're derived from the above further
! down in this file -- so they're simply renamed in place to their
! canonical const_lib names (end_age_stop_active/
! timestep_override_active), now use-associated from const_lib rather
! than locally declared.
      double precision :: endage(50), setdt(50), end_dcen(50), end_xcen(50), end_ycen(50)

! fo: NAMELIST /physics/ member spelled identically to its const_lib
! canonical name -- use-associated directly. fw/fc/fes/fgsf/fmu/fss/
! rcrit (former common/vmult/'s remaining members) are also NAMELIST
! members but need different canonical spellings (difad_velocity_scale/
! mixing_velocity_scale/es_velocity_scale/gsf_velocity_scale/
! mu_gradient_scale/secular_shear_velocity_scale/critical_reynolds),
! so they keep their local NAMELIST-spelled names here and are
! copy-assigned after the namelist read below.
      double precision :: fw, fc, fes, fgsf, fmu, fss, rcrit

! etadh0/etadh1/ldh: NAMELIST /physics/ members, each with a different
! canonical const_lib spelling (debye_huckel_eta_min/
! debye_huckel_eta_max/use_debye_huckel_correction -- no readable
! rename had ever been established for this block, so these names are
! new (2026) rather than picked up from an existing convention), so
! kept local under their NAMELIST spelling here and copy-assigned
! after the namelist read below. cdh/zdh/xxdh/yydh/zzdh/dhnue (former
! common/debhu/'s remaining members) are unused in this file --
! setup/setups.f90 computes cdh/dhnue at startup, mixing/hsubp.f90
! computes zdh/xxdh/yydh/zzdh per-shell -- so they're dropped entirely
! rather than carried forward.
      double precision :: etadh0, etadh1
      logical :: ldh

! ies/imu: NAMELIST /physics/ members spelled identically to their
! const_lib canonical names -- use-associated directly. fesc/fssc/
! fgsfc/igsf (former common/vmult2/'s remaining members) are also
! NAMELIST members but need different canonical spellings
! (es_mixing_scale/secular_shear_mixing_scale/gsf_mixing_scale/
! gsf_inhibition_mode -- the last chosen over parmin's own igsf since
! it's the majority spelling, used directly by 2 of igsf's 3 peer
! files), so they keep their local NAMELIST-spelled names here and are
! copy-assigned after the namelist read below.
      double precision :: fesc, fssc, fgsfc
      integer :: igsf

! grtol/ilambda/niter_gs/ldify: NAMELIST /physics/ members, each with
! a different canonical const_lib spelling (settling_tolerance/
! coulomb_log_choice/settling_num_iterations/diffuse_helium_active),
! so kept local under their NAMELIST spelling here and copy-assigned
! after the namelist read below.
      double precision :: grtol
      integer :: ilambda, niter_gs
      logical :: ldify


! cmin/abstol/reltol/kemmax: NAMELIST /physics/ members, each with a
! different canonical const_lib spelling (min_abundance/
! absolute_tolerance/relative_tolerance/max_burn_iterations), so kept
! local under their NAMELIST spelling here and copy-assigned after the
! namelist read below.
      double precision :: cmin, abstol, reltol
      integer :: kemmax

! former common/lopal95/: iliv95 is not a namelist value and genuinely
! used in this file -- renamed in place to its canonical const_lib
! name (opal95_table_unit), now use-associated rather than locally
! declared.

! dt_gs/xmin/ymin/lthoulfit: NAMELIST /physics/ members, each with a
! different canonical const_lib spelling (settling_timestep_fraction/
! hydrogen_diffusion_floor/helium_diffusion_min/use_thoul_fit), so
! kept local under their NAMELIST spelling here and copy-assigned
! after the namelist read below.
      double precision :: dt_gs, xmin, ymin
      logical :: lthoulfit

! fgry/fgrz/lthoul: NAMELIST /physics/ members spelled identically to
! their const_lib canonical names -- use-associated directly. ldifz is
! also a NAMELIST member but needs a different canonical spelling
! (use_diffusion_z), so it keeps its local NAMELIST-spelled name here
! and is copy-assigned after the namelist read below.
      logical :: ldifz

! lnewdif: NAMELIST /physics/ member with a different canonical
! spelling (use_new_diffusion_routines), so kept local under its
! NAMELIST spelling here and copy-assigned after the namelist read
! below. ldifli is a NAMELIST member spelled identically to its
! const_lib canonical name -- use-associated directly.
      logical :: lnewdif

! lpulse/ipver: NAMELIST /physics/ members, each with a different
! canonical const_lib spelling (pulsation_output_active/
! pulsation_file_version), so kept local under their NAMELIST spelling
! here and copy-assigned after the namelist read below. xmsol (former
! common/pulse/'s remaining member) is unused in this file, so it's
! dropped entirely.
      logical :: lpulse
      integer :: ipver

! po: NAMELIST /physics/ members, all renamed in const_lib (poa/pob/poc/
! pomax/lpout -> po_weight_l/po_weight_teff/po_weight_age/po_max_len_sq/
! po_output_enabled), so kept local under their NAMELIST spelling and
! copy-assigned after the namelist read below.
      double precision :: poa, pob, poc, pomax
      logical :: lpout

! track: NAMELIST /physics/ member, renamed in const_lib (itrver ->
! track_file_version), kept local and copy-assigned below.
      integer :: itrver

! kttau: NAMELIST /physics/ member with a different canonical const_lib
! spelling (atm_choice), so kept local under its NAMELIST spelling
! here and copy-assigned after the namelist read below. kttau0/lttau
! (former common/atmos/'s remaining members) are not namelist values
! -- both computed from kttau right after the namelist read (kttau0 =
! kttau; lttau = .false.) -- so they're simply renamed in place to
! their canonical const_lib names (atm_choice_initial/
! use_ttau_relation), now use-associated rather than locally declared.
! hras is unused in this file, so it's dropped entirely.
      integer :: kttau

! lmhd: NAMELIST /physics/ member (must keep this exact spelling);
! copied into const_lib's use_mhd_eos after the namelist read below.
! iomhd1-8 (former common/mhd/'s remaining members, hardcoded unit
! numbers assigned further down in this file) are not namelist values,
! so they're simply renamed in place to their canonical const_lib
! names (unit_zams_a/b/c/unit_centre1-5), now use-associated from
! const_lib rather than locally declared.
      logical :: lmhd

! core: NAMELIST /physics/ members, all renamed in const_lib (lcore/
! mcore/fcore -> extend_core_inward/num_core_shells_added/
! core_mass_reduction_factor), kept local and copy-assigned below.
      logical :: lcore
      integer :: mcore
      double precision :: fcore

! nwlaol: olaol/oxa/ot/orho/tollaol/iolaol/numofxyz/numrho/numt/llaol/
! iopurez are spelled identically to their const_lib canonical names --
! use-associated directly. lpurez is a NAMELIST /physics/ member with a
! different canonical spelling (use_pure_z_table), so kept local under
! its NAMELIST spelling here and copy-assigned after the namelist read
! below.
      logical :: lpurez

! lrwsh/liso: NAMELIST /physics/ members, each with a different
! canonical const_lib spelling (lrwsh_placeholder/
! isochrone_output_active), so kept local under their NAMELIST spelling
! here and copy-assigned after the namelist read below. iiso (former
! common/chrone/'s remaining member) is not a namelist value and
! genuinely used in this file -- renamed in place to its canonical
! const_lib name (isochrone_file_unit), now use-associated rather than
! locally declared.
      logical :: lrwsh, liso

! newxym: NAMELIST /physics/ members, all renamed in const_lib (xenv0a/
! zenv0a/cmixla/lsenv0a/senv0a -> initial_x_array/initial_z_array/
! mixing_length_array/has_senv0_array/senv0_array), kept local and
! copy-assigned below.
      double precision :: xenv0a(50), zenv0a(50), cmixla(50), senv0a(50)
      logical :: lsenv0a(50)

! former common/atmos2/: the Kurucz surface-pressure table (atmpl/
! atmtl/atmgl/atmz) is unused in this file, renamed to its canonical
! const_lib names (kurucz_log10_pressure_table/kurucz_teff_table/
! kurucz_logg_table/kurucz_table_z) and dropped entirely. ioatm is not
! a namelist value and genuinely used in this file -- renamed in place
! to its canonical const_lib name (atm_table_file_unit), now
! use-associated rather than locally declared.

! former common/atmos2c/: the Kurucz/Castelli surface-pressure table
! (atmplc/atmtlc/atmglc) is unused in this file, renamed to its
! canonical const_lib names (kurucz_castelli_log10_pressure_table/
! kurucz_castelli_teff_table/kurucz_castelli_logg_table) and dropped
! entirely.

! lnulos1: NAMELIST /physics/ member (must keep this exact spelling);
! copied into const_lib's use_itoh_neutrino_loss after the namelist
! read below. dsnudt/dsnudd (former common/nuloss/'s other two
! members) are unused in this file and confirmed dead everywhere else
! too (nuclear/engeb.f90's neutrino_dlnq_dlnt/neutrino_dlnq_dlnd,
! which shared this block only by position, are genuinely local to
! that file and were made plain locals there), so they're dropped
! entirely rather than carried forward.
      logical :: lnulos1

! cals2: NAMELIST /physics/ members, all renamed in const_lib (toll/
! tolr/tolz/lcals/lcalsolzx/calsolzx/calsolage -> luminosity_tolerance/
! radius_tolerance/zx_tolerance/calibrate_solar_model/
! calibrate_solar_zx/target_solar_zx/target_solar_age), kept local and
! copy-assigned below.
      double precision :: toll, tolr, tolz, calsolzx, calsolage
      logical :: lcals, lcalsolzx

! zramp: rsclzc/rsclzm1/rsclzm2/iolaol2/ioopal2/nk are spelled
! identically to their const_lib canonical names -- use-associated
! directly. lzramp is a NAMELIST /physics/ member with a different
! canonical spelling (use_z_ramp), so kept local under its NAMELIST
! spelling here and copy-assigned after the namelist read below.
      logical :: lzramp

! calstar: xls/xlstol/steff/sr/lteff/lcalst are NAMELIST /physics/
! members with different canonical const_lib spellings
! (target_luminosity_lsun/target_star_luminosity_tolerance/
! target_teff/target_radius_rsun/specify_teff_flag/
! calibrate_star_flag -- xlstol's canonical name is disambiguated from
! common/cals2/'s own luminosity_tolerance member, which wins the
! shorter name), so kept local under their NAMELIST spelling here and
! copy-assigned after the namelist read below. bli/alri/ager/blr/blrp/
! agei/lstar/lpassr (former common/calstar/'s remaining members) are
! unused in this file, so they're dropped entirely.
      double precision :: xls, xlstol, steff, sr
      logical :: lteff, lcalst

! opaleos: lopale/lopale01/lopale06/lnumderiv are NAMELIST /physics/
! members with different canonical spellings (use_opal95_eos/
! use_opal2001_eos/use_opal2006_eos/use_numerical_derivatives), kept
! local under their NAMELIST spelling and copy-assigned below (also
! re-synced after the "disable older OPAL EOS" validation block may
! override lopale/lopale01). iopale is spelled identically to its
! canonical name -- use-associated directly.
      logical :: lopale, lopale01, lopale06, lnumderiv

! newopac: zlaol1/zlaol2/zopal1/zopal2/zopal951/zalex1/zkur1/zkur2/
! tmolmin/tmolmax/lalex06/llaol89/lopal92/lopal95/lkur90/lalex95 are
! NAMELIST /physics/ members with different canonical spellings, kept
! local under their NAMELIST spelling and copy-assigned below
! (lalex95/lkur90 also re-synced after the "disable older Alexander
! opacities" validation block may override them). l2z (former
! common/newopac/'s remaining member) is unused in this file, so it's
! dropped entirely.
      double precision :: zlaol1, zlaol2, zopal1, zopal2, zopal951, zalex1, zkur1, zkur2, tmolmin, &
           tmolmax
      logical :: lalex06, llaol89, lopal92, lopal95, lkur90, lalex95

! miscopac: ikur2/icondopacp are spelled identically to their canonical
! names -- use-associated directly. lcondopacp is a NAMELIST /physics/
! member with a different canonical spelling
! (use_conductive_opacity), so kept local under its NAMELIST spelling
! here and copy-assigned after the namelist read below.
      logical :: lcondopacp

! former common/alexo/: ialxo is not a namelist value and genuinely
! used in this file -- renamed in place to its canonical const_lib
! name (alex95_table_unit), now use-associated rather than locally
! declared.

! former common/alex06/: ialex06 is not a namelist value and genuinely
! used in this file -- renamed in place to its canonical const_lib
! name (alex06_table_unit), now use-associated rather than locally
! declared.

! former common/alexmix/: xalex/zalex are not namelist values and
! genuinely used in this file -- renamed in place to their canonical
! const_lib names (alex_mixture_x/alex_mixture_z), now use-associated
! rather than locally declared. Their DATA defaults moved to
! const_lib.f90 since DATA can no longer target them here.

! varfc: vfc/lvfc are spelled identically to their canonical names --
! use-associated directly; lvfc's DATA default moved to const_lib.f90
! since DATA can no longer target it here. ldifad is a NAMELIST
! /physics/ member with a different canonical spelling
! (use_diffusion_advection_transport), so kept local under its
! NAMELIST spelling here and copy-assigned after the namelist read
! below.
      logical :: ldifad

! notran: lnoj is a NAMELIST /physics/ member with a different
! canonical spelling (no_am_transport_in_core), so kept local under
! its NAMELIST spelling here and copy-assigned after the namelist read
! below.
      logical :: lnoj

! sstandard/lnewnuc: NAMELIST /physics/ members (must keep this exact
! spelling, see this file's naming note at the top). Former common
! /cross/; only these two are actually read anywhere (sstandard's
! namelist value is never referenced outside its own declaration --
! setup/remap.f90 fully recomputes const_lib's cross_section_scale
! from other inputs regardless -- so it's dropped here rather than
! copied). lnewnuc is copied into const_lib's use_new_nuclear_rates
! right after the namelist read below, since remap.f90 needs it.
      double precision :: sstandard(17)
      logical :: lnewnuc

! newcross: s0_1_1/s0_3_3/s0_3_4/s0_1_12/s0_1_13/s0_1_14/s0_1_16/
! s0_1_be7e/s0_1_be7p/s0_1_15_c12alp/s0_1_15_o16/s0p_1_1/s0p_3_3/
! s0p_3_4/s0p_1_12/s0p_1_13/s0p_1_14/s0p_1_16/s0pp_1_12/s0pp_1_13/
! s0pp_1_16/s0p_1_be7p/s0pp_1_be7p are NAMELIST /physics/ members, each
! with a different canonical const_lib spelling, kept local under
! their NAMELIST spelling here and copy-assigned below. s0_pep/s0_hep
! are NAMELIST /physics/ members spelled identically to their const_lib
! canonical names -- use-associated directly; their DATA defaults
! moved to const_lib.f90.
      double precision :: s0_1_1, s0_3_3, s0_3_4, s0_1_12, s0_1_13, s0_1_14, s0_1_16, &
           s0_1_be7e, s0_1_be7p, s0_1_15_c12alp, s0_1_15_o16, s0p_1_1, s0p_3_3, s0p_3_4, &
           s0p_1_12, s0p_1_13, s0p_1_14, s0p_1_16, s0pp_1_12, s0pp_1_13, s0pp_1_16, s0p_1_be7p, &
           s0pp_1_be7p

! newparam: all 29 NAMELIST /physics/ members are spelled identically
! to their const_lib canonical names -- use-associated directly.

! former common/monte/: lmonte/imbeg/imend are NAMELIST /physics/
! values spelled identically to their const_lib canonical names --
! use-associated directly rather than locally declared.

! scveos: tlogx/tablex/tabley/smix/tablez/tablenv/nptsx/idtt/idp are
! spelled identically to their const_lib canonical names --
! use-associated directly. lscv is a NAMELIST /physics/ member with a
! different canonical spelling (use_scv_eos), so kept local under its
! NAMELIST spelling here and copy-assigned after the namelist read
! below.
      logical :: lscv

! former common/scv2/: iscvh/iscvhe/iscvz are not namelist values and
! genuinely used in this file -- renamed in place to their canonical
! const_lib names (scv_h_unit/scv_he_unit/scv_z_unit), now
! use-associated rather than locally declared.

! alatm_feh/alatm_alpha/laltptau100: NAMELIST /physics/ members, each
! with a different canonical const_lib spelling (allard_target_feh/
! allard_target_alpha/allard_use_tau100), so kept local under their
! NAMELIST spelling here and copy-assigned after the namelist read
! below. ioatma (former common/alatm03/'s remaining member) is not a
! namelist value and genuinely used in this file -- renamed in place
! to its canonical const_lib name (allard_table_unit), now
! use-associated rather than locally declared.
      double precision :: alatm_feh, alatm_alpha
      logical :: laltptau100

! former common/alatm04/: dummy1-4 are unused in this file, so they're
! dropped entirely.

! tdisk/pdisk/ldisk: NAMELIST /physics/ members (must keep this exact
! spelling); copied into const_lib's disk_temperature/disk_pressure/
! disk_locking_active after the namelist read below. sage (former
! common/disk/'s remaining member) is unused in this file, so it's
! dropped entirely rather than carried forward.
      double precision :: tdisk, pdisk
      logical :: ldisk

! weakscreening: NAMELIST /physics/ member (must keep this exact
! spelling). Former common /weak/; copied into const_lib's
! weak_screening_threshold right after the namelist read below.
      double precision :: weakscreening

! sbrot: NAMELIST /physics/ members, both renamed in const_lib (lsolid/
! impjmod -> force_solid_body_rotation/solid_body_mode_flag), kept
! local and copy-assigned below.
      logical :: lsolid
      integer :: impjmod

! dmdt0/compacc/lmdot: NAMELIST /physics/ members (must keep this
! exact spelling). Former common /masschg/; fczdmdt/ftotdmdt/creim/
! lreimer already match const_lib's canonical spelling, so only these
! three need to stay local (copied into const_lib's
! mass_accretion_rate/accreted_composition/use_mass_accretion right
! after the namelist read below).
      double precision :: dmdt0, compacc(15)
      logical :: lmdot

! cmixing: cstmixing/cstdiffmix are NAMELIST /physics/ members spelled
! identically to their const_lib canonical names -- use-associated
! directly.

! former common/acdpth/: ageout/lclcd/iclcd/ljlast/ljwrt/lacout are not
! namelist values and genuinely used in this file -- renamed in place
! to their canonical const_lib names (ageout_placeholder/
! lclcd_placeholder/iclcd_placeholder/ljlast_placeholder/
! ljwrt_placeholder/acoustic_depth_output), now use-associated rather
! than locally declared. Their DATA defaults moved to const_lib.f90
! since DATA can no longer target them here. tauczn/deladj/tauhe/
! tnorm/tcz/whe/acatmr/acatmd/acatmp/acatmt/tatmos/iacat/ijlast/ladon/
! laoly/ijvs/ijent/ijdel (former common/acdpth/'s remaining members)
! are unused in this file, so they're dropped entirely.

! govs: ltrist is a NAMELIST /physics/ member with a different
! canonical const_lib spelling (use_envelope_triangle_dt), so kept
! local under its NAMELIST spelling here and copy-assigned after the
! namelist read below.
      logical :: ltrist

! pmmwind: NAMELIST /physics/ members, all renamed in const_lib (pmma/
! pmmb/pmmc/pmmd/pmmm/pmmjd/pmmmd/pmmsolp/pmmsolw/pmmsoltau/lmwind/
! lrossby/lbscale/awind -> pmm_exponent_a/pmm_exponent_b/
! pmm_exponent_c/pmm_exponent_d/pmm_exponent_m/pmm_norm_jdot/
! pmm_norm_mdot/pmm_solar_pressure/pmm_solar_omega/
! pmm_solar_turnover_timescale/use_pmm_wind_law/scale_by_rossby_number/
! scale_by_b_field/wind_law_name), kept local and copy-assigned below.
      double precision :: pmma, pmmb, pmmc, pmmd, pmmm, pmmjd, pmmmd, pmmsolp, pmmsolw, pmmsoltau
      logical :: lmwind, lrossby, lbscale
      character(len=3) :: awind

! wmax: NAMELIST /physics/ member with a different canonical const_lib
! spelling (wind_saturation_omega), so kept local under its NAMELIST
! spelling here and copy-assigned after the namelist read below.
! exmd/extau/exr/exm/exl/expr/constfactor/excen/c_2/ljdot0 (former
! common/cwind/'s other NAMELIST-visible or tautological members) are
! spelled identically to their const_lib canonical names --
! use-associated directly; ljdot0's DATA default moved to const_lib.f90
! since DATA can no longer target it here. exw is not a namelist value
! and genuinely used in this file -- renamed in place to its canonical
! const_lib name (wind_law_omega_exponent), now use-associated rather
! than locally declared. structfactor (former common/cwind/'s
! remaining member) is unused in this file, so it's dropped entirely.
      double precision :: wmax

! lnewtcz/lcalcenv: NAMELIST /physics/ members (must keep this exact
! spelling); copied into const_lib's use_new_turnover_timescale/
! calc_envelope_flag after the namelist read below. taucz/taucz0/
! pphot/pphot0/fracstep (former common/ovrtrn/'s other five members)
! are unused in this file -- they're genuinely evolving per-model
! state read/written by many distant files, now state/turnover_lib.f90
! -- so they're dropped from this file's own declarations entirely.
      logical :: lnewtcz, lcalcenv

! mag: NAMELIST /physics/ members, both renamed in const_lib (codm/
! lcodm -> constant_background_diffusion_coeff/
! use_constant_background_diffusion), kept local and copy-assigned
! below.
      double precision :: codm
      logical :: lcodm

! former common/xsect/: xsli6/xsli7/xsbe91/xsbe92/xsbe93/lxli6/lxli7/
! lxbe91/lxbe92/lxbe93 are NAMELIST /physics/ members used only within
! this file (as raw inputs immediately transformed into sli6/sli7/
! sbe91/sbe92/sbe93 below) -- not shared with any other file, so they
! stay ordinary local variables rather than becoming const_lib members.
      double precision :: xsli6, xsli7, xsbe91, xsbe92, xsbe93
      logical :: lxli6, lxli7, lxbe91, lxbe92, lxbe93

! sli6/sli7/sbe91/sbe92/sbe93 are themselves NAMELIST /physics/ members
! (see the "G Somers 6/14" list below) so must keep their exact
! spelling (this file's naming note at the top); the canonical
! const_lib names (li6_rate_scale etc, former common/burnscs/) are set
! via copy-assignment once these are computed below.
      double precision :: sli6, sli7, sbe91, sbe92, sbe93

! spotf/spotx/lsdepth: NAMELIST /physics/ members (must keep this
! exact spelling); copied into const_lib's spot_filling_factor/
! spot_temp_contrast/spot_depth_varies after the namelist read below.
      double precision :: spotf, spotx
      logical :: lsdepth

! former common/version/: yrecver/githash are not namelist values and
! genuinely used in this file -- renamed in place to their canonical
! const_lib names (yrec_version_string/git_hash_string), now
! use-associated rather than locally declared.

      save
!
!
! SPLIT NAMELIST INTO TWO: CONTROL and PHYSICS
      namelist /control/ &
           &    cmixla, calsolage, calsolzx, &
           &    descrip, &
           &    endage, &
           &    flaol, fpurez,flaol2, fopal2, &
           &    flast, ffirst, ffermi, fdebug, ftrack, fshort, fstch, &
           &    fmilne, fmodpt, fstor, fpmod, fpatm, fpenv, &
           &    fdyn, flldat, fsnu, fscomp, fkur, fmhd1, &
           &    fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, fmhd7, fmhd8, fiso, fatm, &
           &    fkur2, fallard, fscvh, fscvhe, fscvz, fopale, fliv95, &
           &    fmonte1,fmonte2, &
           &    ipver, itrver, &
           &    kindrn, &
           &    ldebug, lcorr, lmilne, ltrack, lstore, lfirst, &
           &    lstpch, lscrib, lstch, &
! G Somers 11/14
           &    lstatm,lstenv,lstmod,lstphys,lstrot, &
! G Somers END
           &    lpulse, lzramp, lteff, lcalst, lpurez, &
! MHP 9/24 add LCALSOLZX to namelist
           &    liso, lrwsh, lsenv0a, lpout,lcals,lcalsolzx, &
           &    llaol89,lopal92,lopal95,lkur90,lalex95, &
           &    npoint, &
           &    npenv, nprtmod, nprtpt, numrun, nmodls, &
           &    opecalex, &
           &    poa, pob, poc, pomax, &
           &    rsclm, rsclx, rsclz, rsclcm, rsclzc, rsclzm1, rsclzm2, &
           &    setdt, senv0a,steff,sr, &
           &    tolr, toll,tolz, &
           &    xenv0a, xls, xlstol, &
           &    zenv0a, &
           &    zlaol1,zlaol2,zopal1,zopal2, zopal951, &
           &    zopal952, zalex1, zalex2, zkur1, zkur2, &
           &      fopale01,fcondopacp,fopale06,falex06,lalex06, &
! MHP 10/24 ADDED END_DCEN,END_XCEN,END_YCEN VECTORS TO NML1, USED IN MAIN
! MHP 10/24 ADDED HEAVY ELEMENT MIXTURE CONTROLS TO NML1,USED IN STARIN
           &  end_dcen,end_xcen,end_ycen,isetmix,isetiso, &
           &  amix,aiso,frac_c,frac_n,frac_o,r12_13,r14_15,r16_17,r16_18,zxmix, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini, &
! new (2026): GYRE-format periodic pulsation output interval, additive
! only -- absent from every existing *.nml1 file, so it simply keeps
! its default (off) there.
           &  pulse_gyre_interval
!
      namelist /physics/ &
           &    atmmin, atmbeg, atmerr, atmmax, atmd0, anewcp, atmp, acfpft, &
           &    atime, alphac, alphae, alfa, alpham, atmstp, abstol, betac, &
           &    cmin, clsun, crsun, &
           &    dpenv, dtdif, dtwind, djok, dt_gs, &
           &    enverr, envmax, envmin, envbeg, envstp,etadh0,etadh1, &
           &    fcorr0, fcorri, fk,  fw, fc, fo, fmu, fes, &
           &    fcore, fgsf, fss, fesc, fssc, fgsfc, fgry, fgrz, &
           &    grtol, &
           &    htoler, hpttol, &
           &    itfp1, itfp2, imax, itdif1, itdif2, ies, igsf, imu, ilambda, &
           &    kttau, kemmax, &
           &    lvfc, ldifad, lnoj, lnewdif, ldify, ldifz, ldifli, lsnu, ldh, &
           &    lnewcp, lkuthe, lovstc, lovste, lovstm, lovmax, &
           &    lexcom, lrot, lnew0, linstb, lwnew, ljdot0, lptime,ladov,ltrist, &
           &    lenvg, lnulos1, lthoul, lthoulfit, &
           &    lopale, lmhd, lcore, lsemic, lnews, &
           &    mcore, &
           &    niter1, niter2, niter3, niter4, nuse, niter_gs, &
           &    optol, &
           &    rcrit, reltol, &
           &    stolr0, &
           &    tcut, tscut, tenv0, tenv1, tgcut, tridt, tridl, &
           &    tollaol, &
           &    vnew, &
           &    walpcz, wnew, weakscreening, &
           &    xnewcp, xmin, &
           &    ymin, tmolmin,tmolmax, &
           &    lmonte,imbeg,imend,sstandard,lscv, &
           &    ldisk,tdisk,pdisk,wmax,lsolid,impjmod,  &  !JNT 09/2025 FOR 05/15
           &    dmdt0,fczdmdt,ftotdmdt,compacc,creim,lreimer,lmdot, &
           &    lopale01,lcondopacp,lopale06,lnumderiv, &
           &    alatm_feh,alatm_alpha,laltptau100, &  ! For new Allard Atmospheres
           &    cstmixing, cstdiffmix,       &  !CFD oct2009 To mimic mixing(reduce settling)
           &    lsolwind,lmwind,lrossby,lpmm,lbscale, &
           &    awind,pmma,pmmb,pmmc,pmmd,pmmm,pmmjd,pmmmd,pmmwmax, &
! MHP 8/17 ADDED WMAX_SUN
           &    pmmsolp,pmmsolw,pmmsoltau,lcodm,codm,wmax_sun, &
! G Somers 6/14
           &    xsli6,xsli7,xsbe91,xsbe92,xsbe93, &
           &    lxli6,lxli7,lxbe91,lxbe92,lxbe93, &
           &    sli6,sli7,sbe91,sbe92,sbe93,lnewnuc, &
           &    spotf, spotx, lsdepth, &
! G Somers END
! MHP 09/14 ADDED CROSS SECTIONS
           &    s0_1_1,s0_3_3,s0_3_4,s0_1_12,s0_1_13,s0_1_14,s0_1_16, &
           &    s0_pep,s0_1_be7e,s0_1_be7p,s0_hep,s0_1_15_c12alp,s0_1_15_o16, &
           &    s0p_1_1,s0p_3_3,s0p_3_4,s0p_1_12,s0p_1_13,s0p_1_14, &
           &    s0p_1_16,s0pp_1_12,s0pp_1_13,s0pp_1_16,s0p_1_be7p,s0pp_1_be7p, &
           &    flag_dx,flag_dw,flag_dz,lstruct_time, &
           &    time_core_min,time_dl,time_dp,time_dr,time_dt,time_dw_global, &
           &    time_dw_mix,time_dx_core_frac,time_dx_core_tot,time_dx_shell, &
           &    time_dx_total,time_dy_core_frac,time_dy_core_tot,time_dy_shell, &
           &    time_dy_total,tol_czbase_fine_width,tol_dl_max,tol_dm_max, &
           &    tol_dm_min,tol_dp_core_max,tol_dp_czbase_max,tol_dp_env_max, &
           &    tol_dx_max,tol_dz_max,time_max_dt_frac,lnewvars, &
! G Somers 3/17 USE NEW OVERTURN TIMESCALE CALC?
           &    lnewtcz, lcalcenv
!
! DBG DATA CARDS FOR THE RUN PARAMETERS
! MHP DATA FOR MONTE CARLO OPTION, ETC
! lmonte/imbeg/imend defaults moved to const_lib.f90 (former
! common/monte/).
! Changed slightly 3He-3He on 9/25/97 to take account of the S'.
!  Previously (6/16/97) used S at Gamow Peak. Agrees with Workshop paper.
!
      data weakscreening/0.03d0/
      data sstandard/0.9828,1.0485,0.9815,0.9241,1.3818,1.0542,1.0, &
           &  1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.7819,0.2875/
! MHP 7/93 VARIABLE FC OPTION
! MHP 9/94 COMBINED DIFFUSION/ADVECTION OPTION
! lvfc's default moved to const_lib.f90 (former common/varfc/): DATA
! can no longer target it here now that it's use-associated.
      data ldifad/.false./
! MHP 9/93
      data lnoj/.false./
      data tdisk,pdisk,ldisk/0.0d0,7.2722d-6,.false./
! alex_mixture_x/alex_mixture_z defaults moved to const_lib.f90 (former
! common/alexmix/): DATA can no longer target them here now that
! they're use-associated.
      data lsenv0a, senv0a /50*.false.,50*1.26d-4/
      ! xenv0/zenv0 defaults moved to const_lib.f90 (former common/label/).
! ldebug/lcorr/npoint/lmilne/ltrack/lstore/lstpch/lscrib/lstch/nprtmod
! defaults moved to const_lib.f90 (see its own header note): DATA can
! no longer target them here now that they're use-associated from
! const_lib.
      data lenvg, atmstp, envstp/.false.,0.5,0.5/
      data nprtpt/1/
! pulse_gyre_interval default moved to const_lib.f90 (former
! common/pulsegyre/).
      data numrun, kindrn, lfirst, nmodls &
           &      /1,50*1,50*.true.,50*0/
! MHP 10/24 ADDED NEW DEFAULTS FOR END CONDITIONS ON CENTRAL D,X,Y
      data endage, setdt, rsclm, rsclx, rsclz, rsclcm &
           &      /50*0.0,50*0.0,50*0.0,50*0.0,50*0.0,50*0.0/
      data end_dcen,end_xcen, end_ycen &
           &      /50*0.0,50*0.0,50*0.0/
      data element_id/'HE3','C12','C13','N14','N15','O16','O17','O18','H2 ', &
           & 'LI6','LI7','BE9'/
      data optol/1.0d-8/
! zsi's default (0.0d0) moved to const_lib.f90 -- DATA can no longer
! target it here now that it's use-associated.
! tcut/saha_log10t_cutoff/tenv0/tenv1/tgcut defaults moved to
! const_lib.f90 (see its own header note): DATA can no longer target
! them here now that they're use-associated from const_lib.
      data atmerr,atmd0,atmbeg,atmmin,atmmax/3.0d-4,1.0d-10,1.0d-1, &
           &      1.0d-1,5.0d-1/
      data enverr,envbeg,envmin,envmax/3.0d-4,1.0d-1,1.0d-1,5.0d-1/
      data stolr0,imax,nuse/1.0d-3,11,7/
! dtdif/itdif1/itdif2/htoler defaults moved to const_lib.f90 (former
! common/difus/, common/ctol/): DATA can no longer target them here
! now that they're use-associated.
      data djok/1.0d-4/
      data hpttol/1.0d-8,8.0d-2,5.0d-2,5.0d-2,1.0d0,1.0d0,0.0d0,5.0d-2, &
           &             2.0d-2,5.0d-2,5.0d-2,1.0d-1/
      data lnewcp,anewcp,xnewcp/.false.,'   ',1.3d1/
! MHP 10/24 ADDED NEW MIXTURE CONTROL ISETISO CONTROLS CNO ISOTOPE RATIOS AND
! LIGHT ELEMENT ABUNDANCES D,HE3,LI6,LI7,BE9,B10,B11 (1=USED)
! ISETMIX CONTROLS C+N+O MASS FRACTIONS (1=USED)
! AMIX AND AISO ARE STRINGS IDENTIFYING EITHER A PRESET MIXTURE OR A CUSTOM ONE ('CUS')
! SUPPORTED AMIX AT PRESENT ARE 'GS98','AAG21',M22M','M22P'. SUPPORTED AISO IS 'L21'.
! FOR A CUSTOM MIXTURE YOU CAN ENTER INDIVIDUAL VALUES.
! TO BE ADDED - AMIX FROM ATOMIC OPACITY TABLES (INEWMIX=2) AND OTHER MIXTURES/ISOTOPES
      data isetiso,isetmix,amix,aiso/0,0,'GS98','L21'/
!     L21 DEFAULT ISOTOPE DATA LODDERS ET AL. 2021 SSRV 217,44
      data r12_13,r14_15,r16_17,r16_18,xh2_ini,xhe3_ini,xli6_ini, &
           &      xli7_ini,xbe9_ini,xb10_ini,xb11_ini/88.27d0,411.9d0,471.4d0, &
           &      2693.0d0,2.781d-5,3.461d-5,7.187d-10,1.025d-8,1.595d-10, &
           &      1.002d-9,4.405d-9/
!     GS98 DEFAULT CNO FRACTIONS OF METALS GREVESSE&sAUVAL 1998 SSRV 85,161
      data zxmix,frac_c,frac_n,frac_o/0.02292d0,0.172148d0,0.050426d0, &
           &      0.468195/
! MHP 10/24 DATA FOR CNO FRACTIONS AND Z/X OF DIFFERNT SOLAR MIXTURES.
!     AMIXT IS THE LIST OF IDS,EACH OF WHICH HAS A ZX AND CNO FRACS
!     ENTRY 1 =GS98(IN PARMIN),2=ASPLUND ET AL. 2021 A&A 653,141
!     3,4=MAGG ET AL. 2021 (MET,PHOT) A&A 661,140
      data mixture_id_table,zx_mix_table,frac_c_table,frac_n_table,frac_o_table/'GS98','AAG21', &
           &      'M22P','M22M',0.02292d0,0.0187d0,0.0225d0,0.0226d0, &
           &  0.172148d0,0.184156d0,0.19239d0,0.191425d0, &
           &  0.050426d0,0.050344d0,0.059012d0,0.058716d0, &
           &  0.468195d0,0.416592d0,0.415631d0,0.413545d0/
! acfpft/itfp1/itfp2 defaults moved to const_lib.f90 (former
! common/rot/).
      data tridt,tridl/1.0d-2,8.0d-2/
! niter1/niter2/niter3/fcorr0/fcorri defaults moved to const_lib.f90
! (former common/ctol/).
! atime's default moved to const_lib.f90 -- see the tcut/etc. note
! above; ATIME(13) was orginally = 1.5.
      data dtwind /1.0d1/
      data lptime/.true./
! JVS 04/14
      data ltrist/.false./
! cmixl's default (1.4d0) moved to const_lib.f90's own declaration: a
! DATA statement can no longer target it here now that it's
! use-associated from const_lib rather than locally declared/COMMON.
      data lkuthe/.false./
!       DATA DPENV,LNSTDMX,LOVSTC,ALPHAC,LOVSTE,ALPHAE, LOVSTM, ALPHAM
!      */1.0D0,.FALSE.,.FALSE., 0.0D0, .FALSE.,0.0D0, .FALSE., 0.0/
! dpenv/lovstc/alphac/alphae/lovstm/alpham defaults moved to
! const_lib.f90 (former common/dpmix/); lovste stays local (NAMELIST
! spelling).
      data lovste/.false./
! ladov/lovmax/betac defaults moved to const_lib.f90 (former
! common/dpmix/).
! JVS 07/13
! END JVS
! lnew0 default moved to const_lib.f90 (former common/cenv/).
      data lexcom/.false./
! walpcz/lwnew/wnew defaults moved to const_lib.f90 (former
! common/rot/); lrot/linstb stay local (NAMELIST spelling).
      data lrot,linstb/.false.,.false./
! ljdot0's default moved to const_lib.f90 (former common/cwind/): DATA
! can no longer target it here now that it's use-associated.
      data alfa,fk/1.5d0,1.0d0/
! fo's default moved to const_lib.f90 (former common/vmult/): DATA can
! no longer target it here now that it's use-associated.
      data fw,fc,fes,fgsf,fmu,fss,rcrit/1.0d0, &
           &      1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d3/
! MHP 8/17 INITIALIZED WMAX_SUN
      data wmax,wmax_sun/3.0d-4,1000.0/
! DBG PULSE DATA CARD FOR PULSATION
      data lpulse/.false./
      data ipver/1/
      data itrver/1/
      data kttau/0/
      data clsun,crsun/3.8515d33,6.9598d10/
! YC  If LMHD is TRUE use MHD equation of state tables.  LU numbers
!     are stored in IOMHDi.
! DBG If LCORE is TRUE then calculate shells interior to start up
!     model's inner most shell.
      data lcore,mcore,fcore/.false.,0,1.0/
      data lmhd/.false./
! MHP 5/90 NEW DATA STATEMENTS FOR NEW PARAMETERS
      data grtol,ilambda,niter_gs,ldify/1.0d-8,1,10, &
           &      .false./
! niter4/lnews/lsnu defaults moved to const_lib.f90 (former
! common/neweng/).
      data dt_gs,xmin,ymin,lthoulfit/0.1d0,1.0d-3,1.0d-3,.false./
! lthoul/fgrz/fgry/ldifli defaults moved to const_lib.f90 (former
! common/gravs3/, common/gravs4/): DATA can no longer target them here
! now that they're use-associated.
      data ldifz/.false./
      data lnewdif/.false./
      data cmin,abstol,reltol,kemmax/1.0d-20,1.0d-5,1.0d-4,50/
      data etadh0, etadh1, ldh/-1.0d0, 1.0d0, .false./
      data fesc,fssc,fgsfc/1.0d0,1.0d0,1.0d0/
! ies/imu defaults moved to const_lib.f90 (former common/vmult2/):
! DATA can no longer target them here now that they're use-associated.
      data igsf/1/
! lsemic's default moved to const_lib.f90 (former common/dpmix/).
! tollaol/llaol defaults moved to const_lib.f90 (former common/nwlaol/):
! DATA can no longer target them here now that they're use-associated.
! DBGLAOL
      data lpurez/.false./
! DBG 11/11/91
      data lrwsh,liso/.false.,.false./
! 3/92 DBG
      data lnulos1/.false./
! DBG PULSE OUT 7/92
      data pomax,poa,pob,poc,lpout/0.1d0,1.0d0,10.0d0,0.0d0,.false./
! MHP 06/13 ADDED FLAG TO CALIBRATE TO SOLAR Z/X, SOLAR Z/X, SOLAR AGE
      data toll,tolr,tolz,lcals,lcalsolzx,calsolage,calsolzx/1.0d-5, &
           &      1.0d-4,1.0d-3,.false.,.false.,4.57d9,0.02292d0/
!      DATA TOLL,TOLR,LCALS/1.0D-5,1.0D-4,.FALSE./
! DBG 4/94 ZRAMP STUFF
! rsclzc/rsclzm1/rsclzm2 defaults moved to const_lib.f90 (former
! common/zramp/): DATA can no longer target them here now that they're
! use-associated.
      data lzramp/.false./
! DBG 12/94 CALIBRATED STELLAR MODEL STUFF
      data lcalst, lteff/.false., .false./
! YCK >>>  2/95 OPAL eos
! LLP >>> OPAL 2001 EOS, Potekhin Conductive Opacities,
!         OPAL 2006 EOS, Use Numerical Derivitives switches
      data lopale, lopale01,lcondopacp,lopale06,lnumderiv &
           &      /.false.,.false.,.false.,.false.,.false./
! MHP 8/25 Removed hard coded defaults
!     Alex low T opacity
!      DATA OPECALEX/'OPACALEXANDER.X00',
!     +              'OPACALEXANDER.X01',
!     +              'OPACALEXANDER.X02',
!     +              'OPACALEXANDER.X035',
!     +              'OPACALEXANDER.X05',
!     +              'OPACALEXANDER.X07',
!     +              'OPACALEXANDER.X08'/

! DBG 1/96 THE ARRAY V, READ IN VIA RDLAOL, CONTAINED THE MASS FRACTIONS
! OF THE ENVELOPE ELEMENTS. IT WAS USED IN STARIN TO DEFINE FXENV,
! WHICH ARE THE NUMBER FRACTION OF THE ENVELOPE ELEMENTS. FXENV WAS
! THEN UPDATED IN EQSTAT AND EQSAHA. HERE WE DEFINE VNEW PASSED
! IN A COMMON BLOCK VNEWCB. IT IS IDENTICAL TO V EXCEPT THAT THE NUMBERS
! ARE DEFINED HERE EXPLICITY FOR A G&N93 SOLAR MIXTURE. YOU CAN
! CHANGE THEM VIA THE PHYSICS NAMELIST. V IS SET EQUAL TO VNEW IN
!  STARIN EXCEPT WHEN LLAOL=T (TO MAINTAIN BACKWARD COMPATIBILITY.
!            Na          Al          Mg          Fe
!            Si          C           H           O
!            N           Ar          Ne          He
! vnew default moved to const_lib.f90 (former common/vnewcb/).
! MHP 5/97 OPTION FOR SAUMON, CHABRIER, AND VAN HORN EOS ADDED
      data lscv/.false./
! MHP 3/99 OPTION FOR SB ROTATION ENFORCED IN THE ENTIRE STAR AT
! ALL TIMES
! JNT 09/2025 FOR 05/15 IMPJMOD default set to 0
      data lsolid, impjmod/.false.,0/   !JNT 09/2025
! fczdmdt/ftotdmdt/creim/lreimer defaults moved to const_lib.f90
! (former common/masschg/); dmdt0/compacc/lmdot stay local (NAMELIST
! spelling).
      data dmdt0,compacc,lmdot &
           &      /-1.0d-14, &
           &      0.71668d0,0.265721d0,0.01757d0,2.9d-5, &
           &      3.013d-3,3.385d-5,9.346d-4,0.0d0,8.462d-3,0.0d0, &
           &      1.696d-5,0.0d0,2.0d-9,2.0d-9,3.0d-11, &
           &      .false./
! mhp 8/10 added scaled solar wind mass loss option
!      DATA LSOLWIND,DMSUN,DMWSUN,DMWMAX/.FALSE.,-2.0D-14,2.863E-6,9.054E-5/
      data lsolwind/.false./
! 3/09 Alexander 2006 opacity table options and opacity ramp options
      data tmolmin,tmolmax,lalex06/4.0d0,4.1d0,.false./
!FD 10/09 Mimic mixing options - acting on setling and differential settling
! cstmixing/cstdiffmix's defaults moved to const_lib.f90 (former
! common/cmixing/): DATA can no longer target them here now that
! they're use-associated.
! JVS 02/11 Initialize acoustic depth common block values appropriately
! ageout_placeholder/lclcd_placeholder/ljlast_placeholder/
! ljwrt_placeholder/acoustic_depth_output defaults moved to
! const_lib.f90 (former common/acdpth/): DATA can no longer target
! them here now that they're use-associated.
! JVS end
! MHP 02/12 NEW PARAMETERIZATION OF ANGULAR MOMENTUM AND MASS LOSS
! FROM MAGNETIZED SOLAR-LIKE WINDS
!      DATA LPMM,PMMA,PMMB,PMMC,PMMM,PMMJD,PMMMD,PMMWMAX,
!     *     PMMSOLP,PMMSOLW,PMMSOLTAU
!     *   /.FALSE.,2.0D0,1.0D0,0.0D0,0.22D0,1.32E30,1.27E12,
!     *     2.836E-5,4.9304D0,2.836E-6,1.065E6/
! G Somers 6/16 NEW PARAMETERIZATION OF ANGULAR MOMENTUM AND MASS LOSS FROM
! SOLAR WINDS. FOLLOW MATT ET AL (2012) FORMULATION, BUT DEFAULT TO KAWALER
! TYPE LAW.
      data lmwind,lrossby,lbscale/.false.,.true.,.false./
! MHP 8/17 CHANGED DEFAULT FOR PMMA TO 2 FROM 0
      data awind,pmma,pmmb,pmmc,pmmd,pmmm/'K97',2.0,1.0,2.0,0.0,0.5/
      data pmmjd,pmmmd,pmmsolp,pmmsolw,pmmsoltau/1.32e30,1.27e12, &
           &      4.9304d0,2.836e-6,9.4805d5/
! MHP 02/12 PERMIT CONSTANT DIFFUSION COEFFICIENT
      data lcodm,codm/.false.,2.5d4/
!
! G Somers 06/14 ALLOW NEW LI DESTRUCTION CROSS SECTIONS
!           NEW VALUES SHOULD BE IN UNITS OF keV b.
!           DEFAULT LI6 = 5,500 keV b FROM FOWLER ET AL. 1967
!           DEFAULT LI7 = 52 keV b FROM ROLFS & KAVANAGH 1986
!           DEFAULT BE9(P,G)B10 (BE91) = 1.1 keV b INFERRED FROM FOWLER 1988
!               THEY DON'T SHOW THIS REACTION, SO THIS IS APPROXIMATE.
!           DEFAULT BE9(P,D)2HE4 (BE92) = 15,000 keV b FROM FOWLER ET AL. 1967
!           DEFAULT BE9(P,A)LI6 (BE93) = 15,000 keV b FROM FOWLER ET AL. 1967
      data xsli6, xsli7, xsbe91, xsbe92, xsbe93 &
           &      /5.5d3, 5.2d1, 1.1d0, 1.5d4, 1.5d4/
      data lxli6, lxli7, lxbe91, lxbe92, lxbe93, &
           &      sli6, sli7, sbe91, sbe92, sbe93 &
           &      /.false.,.false.,.false.,.false.,.false., &
           &      1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/
      data spotf, spotx, lsdepth/0.00, 1.00, .false./
! G Somers END
! MHP 8/14 DEFAULT CROSS-SECTIONS ARE TAKEN FROM THE SOLAR FUSION II PAPER
! REFERENCE ADELBERGER ET AL. 2011. UNITS ARE KeV b
      data s0_1_1,s0_3_3,s0_3_4,s0_1_12,s0_1_13/ &
           &      4.01d-22,5.21d3,5.6d-1,1.34d0,7.6d0/
! s0_pep/s0_hep's defaults moved to const_lib.f90 (former
! common/newcross/): DATA can no longer target them here now that
! they're use-associated.
      data s0_1_14,s0_1_16,s0_1_be7e/ &
           &      1.66d0,1.06d1,1.7709d-10/
! NOTE: PEP IS THE PROPORTIONALITY CONSTANT RELATIVE TO PP
! NOTE: BE7+E- IS THE PROPORTIONALITY CONSTANT IN THE LINEAR TERM
! THE CODE USES T9, NOT T6, SO ANY EXPRESSION IN TERMS OF T/10^6 K
! NEEDS TO BE DIVIDED BY 1000^0.5 (FOR BOTH PEP AND BE7+E-)
      data s0_1_be7p,s0_1_15_c12alp,s0_1_15_o16/ &
           &      0.0208d0,7.3d4,3.6d1/
! REFERENCE FIRST DERRIVATIVES OF CROSS-SECTIONS (ADELBERGER ET AL. 2011)
! UNITS ARE b
      data s0p_1_1,s0p_3_3,s0p_3_4,s0p_1_12/ &
           &      4.49d-24,-4.9d0,-3.6d-4,2.6d-3/
      data s0p_1_13,s0p_1_14,s0p_1_16,s0p_1_be7p/ &
           &      -7.83d-3,-3.3d-3,-5.4d-2,-3.12d-5/
! REFERENCE SECOND DERIVATIVES OF CROSS SECTIONS (ADELBERGER ET AL. 2011)
      data s0pp_1_12,s0pp_1_13,s0pp_1_16,s0pp_1_be7p/ &
           &      8.3d-5,7.29d-4,0.0d0,-2.288d-7/
      data lnewnuc /.false./
! All 30 former common/newparam/ members' DATA defaults moved to
! const_lib.f90: DATA can no longer target them here now that they're
! use-associated.
! DEFAULT TO YES FOR NEW TAUCZ CALCULATION
      data lnewtcz,lcalcenv /.true.,.true./
!
! THIS SUBROUTINE READS ALL USER DEFINED QUANTITIES FROM THE
! FILES yrec8.nml1 and yrec8.nml2
! VALUES FOR LOGICAL UNITS USED IN READS AND WRITES SET IN DATA
! STATEMENT; THEY CAN BE CHANGED IN THE NAMELIST IF NEEDED.
! LOGICAL UNIT 5 = READ FROM SCREEN
! LOGICAL UNIT 6 = WRITE TO SCREEN: FOR BATCH USE STATUS FILE INSTEAD
! SPECIFY ALL LOGICAL UNIT NUMBERS HERE:
! OUTPUT: STATUS FILE
      iowr = 6
! OUTPUT: LAST MODEL (TEXT)
      ilast = 11
! INPUT: FIRST MODEL (TEXT)
      first_unit = 12
! INPUT: PHYSICS NAMELIST
      run_unit = 13
! INPUT: CONTROL NAMELIST
      standard_unit = 14
! INPUT: FERMI TABLES
      fermi_unit = 15
! OUTPUT: RESERVED for DEBUGGING
      idebug = 18
! OUTPUT: TRACK
      itrack = 19
! OUTPUT: ALL DIAGNOSTIC INFO
      short_file_unit = 20
! OUTPUT: MILNE INVARIANT VARIABLES
      imilne = 21
! OUTPUT: SHELL BY SHELL INFO ON MODELS
      imodpt = 22
! OUTPUT: SAVED MODELS, CAN BE USED AS STARTING MODEL
      istor = 23
! OUTPUT: FOR PULSATION CODE, INTERIOR
      opal_model_unit = 24
! OUTPUT: FOR PULSATION CODE, ENVELOPE
      opal_envelope_unit = 25
! OUTPUT: FOR PULSATION CODE, ATMOSPHERE
      opal_atm_unit = 26
! OUTPUT: BINARY OUTPUT OF LAST MODEL
      last_model_binary_lu = 27
! OUTPUT: BINARY OUTPUT OF STORED MODELS
      stored_models_binary_lu = 28
! INPUT: BINARY STARTING MODEL
      first_model_binary_lu = 29
! OUTPUT: INFO RELAVENT TO DYNAMO
      dynamics_unit = 30
! YCK INPUT: OPAL92 OPACITY TABLES
      laol_table_unit = 32
! YCK INPUT: OPAL95 OPACITY TABLE
      opal95_table_unit = 48
! OUTPUT: SNU FLUXES
      neutrino_unit = 33
! OUTPUT: EXTENDED COMPOSITION INFO
      composition_unit = 34
! YCK INPUT: KURUCZ LOW T OPACITIES
      kurucz_table_unit = 36
! OUPUT: ISOCHRONE INFORMATION
      isochrone_file_unit = 37
! INPUT: KURUCZ ATMOSPHER TABLE
      atm_table_file_unit = 38
! YCK INPUT: Alex LOW T OPACITIES
      alex95_table_unit = 39
! INPUT: MHD EQU. OF STATE TABLES
      unit_zams_a = 40
      unit_zams_b = 41
      unit_zams_c = 42
      unit_centre1 = 43
      unit_centre2 = 44
      unit_centre3 = 45
      unit_centre4 = 46
      unit_centre5 = 47
! INPUT: OPAL EQUATION OF STATE
      iopale = 49
! INPUT LAOL OPACITIES IN DENSE GRID FORMAT
      iolaol = 61
! INPUT: LAOL OPACITIES FOR PURE CN IN DENSE GRID FORMAT
      iopurez = 62
! DBG 4/94
!     INPUT:
! DBG 8/95 SECOND OPOACITY TABLES FOR ZRAMP AND Z DIFFUSION
      iolaol2 = 63
      ioopal2 = 64
      ikur2 = 65
! MHP 6/97 ADDED OPTION FOR ALLARD MODEL ATMOSPHERES
      allard_table_unit = 66
! MHP 6/98 MONTE CARLO FOR SNUs
      monte_carlo_unit1 = 70
      monte_carlo_unit2 = 71
! INPUT FILES FOR THE SCV EOS
      scv_h_unit=72
      scv_he_unit=73
      scv_z_unit=74
! Input files for Potekhin conductive opacities. LLP 7/8/06
      icondopacp = 75
!      FcondOpacP = 'condall.d'
! Default FeH and Alpha for new Allard Atmospheres
      alatm_feh = 0d0
      alatm_alpha = 0d0
      laltptau100 = .false.
! 3/09 Input file for 2006 Alexander opacities
      alex06_table_unit = 90

      print *,''
      print *,'Yale Rotating Evolution Code - YREC, v',yrec_version_string(1:len_trim(yrec_version_string)),' (',git_hash_string(1:len_trim(git_hash_string)),')'

! JVS 02/11 Altered the yrec8 input format so that files can be entered
! on the command line, with the *.nml1 as the first argument, and *.nml2 as
! the second. Defaults to yrec8.nml1 and yrec8.nml2 if none are provided.
!      OPEN(UNIT=ISTAND, FILE='yrec8.nml1', STATUS='OLD')
!      OPEN(UNIT=IRUN, FILE='yrec8.nml2', STATUS='OLD')
!      READ(UNIT=ISTAND, NML=CONTROL)
!      READ(UNIT=IRUN, NML=PHYSICS)
!      CLOSE(ISTAND)
!      CLOSE(IRUN)

! Dynamically create format string so version info is nicely spaced
      write(version_fmt, 315) len_trim(yrec_version_string), len_trim(git_hash_string)
      315 format('(''# YREC v'', A', i2.2, ', '' ('', A', i2.2, &
           &        ', '')'')')

      call getarg(1, control_nml_file)
      if (control_nml_file(1:2) .eq. ' ') control_nml_file = 'yrec8.nml1'
      print *, ' '
      write(*,*) 'CONTROL namelist :  ',control_nml_file(1:len_trim(control_nml_file))

      call getarg(2, physics_nml_file)
      if (physics_nml_file(1:2) .eq. ' ') physics_nml_file = 'yrec8.nml2'
      write(*,*) 'PHYSICS namelist :  ',physics_nml_file(1:len_trim(physics_nml_file))

      open(unit=standard_unit, file=control_nml_file, status='OLD')
      open(unit=run_unit, file=physics_nml_file, status='OLD')
      read(unit=standard_unit, nml=control)
      read(unit=run_unit, nml=physics)
      close(standard_unit)
      close(run_unit)
! stolr0/imax/nuse must keep their exact NAMELIST /physics/ spelling
! (see this file's naming note at the top), so intpar_lib's
! canonically-named variables are set by copying from them here,
! after the namelist read, rather than by renaming in place.
      tolerance_fraction = stolr0
      max_stage_index = imax
      extrap_order = nuse
! tscut must likewise keep its NAMELIST spelling; copy into const_lib's
! saha_log10t_cutoff here.
      saha_log10t_cutoff = tscut
! use_new_nuclear_rates/weak_screening_threshold: same reasoning,
! copied from their NAMELIST-spelled locals before remap runs, since
! remap.f90 reads use_new_nuclear_rates to decide how to compute
! const_lib's cross-section-scale members.
      use_new_nuclear_rates = lnewnuc
      weak_screening_threshold = weakscreening
! envelope_overshoot_active/rotation_active/instability_transport_
! active/mass_accretion_rate/accreted_composition/use_mass_accretion:
! same reasoning, copied from their NAMELIST-spelled locals. Order
! relative to `call remap` doesn't matter for these -- remap.f90
! doesn't read any of them.
      envelope_overshoot_active = lovste
      rotation_active = lrot
      instability_transport_active = linstb
      mass_accretion_rate = dmdt0
      accreted_composition = compacc
      use_mass_accretion = lmdot
! use_itoh_neutrino_loss/use_new_turnover_timescale/calc_envelope_flag:
! same reasoning, copied from their NAMELIST-spelled locals.
      use_itoh_neutrino_loss = lnulos1
      use_new_turnover_timescale = lnewtcz
      calc_envelope_flag = lcalcenv
! clsun/crsun must likewise keep their NAMELIST spelling; copy into
! const_lib's solar_luminosity_cgs/solar_radius_cgs here (before
! setup/setups.f90 computes the rest of former common/const/ from
! these two).
      solar_luminosity_cgs = clsun
      solar_radius_cgs = crsun
! lexcom must likewise keep its NAMELIST spelling; copy into
! const_lib's use_extended_composition here.
      use_extended_composition = lexcom
! spotf/spotx/lsdepth, tdisk/pdisk/ldisk, endage/setdt/end_dcen/
! end_xcen/end_ycen, lmhd, and optol must likewise keep their NAMELIST
! spelling; copy into const_lib's canonical names here.
      spot_filling_factor = spotf
      spot_temp_contrast = spotx
      spot_depth_varies = lsdepth
      disk_temperature = tdisk
      disk_pressure = pdisk
      disk_locking_active = ldisk
      target_end_age = endage
      timestep_override = setdt
      central_deuterium_stop = end_dcen
      central_hydrogen_stop = end_xcen
      central_helium_stop = end_ycen
      use_mhd_eos = lmhd
      metal_fraction_match_tolerance = optol
      print_point_interval = nprtpt
      convergence_tolerance = djok
      settling_tolerance = grtol
      coulomb_log_choice = ilambda
      settling_num_iterations = niter_gs
      diffuse_helium_active = ldify
      settling_timestep_fraction = dt_gs
      hydrogen_diffusion_floor = xmin
      helium_diffusion_min = ymin
      use_thoul_fit = lthoulfit
      use_diffusion_z = ldifz
      use_new_diffusion_routines = lnewdif
      atm_error_tol = atmerr
      atm_step_initial = atmd0
      atm_step_begin = atmbeg
      atm_step_min = atmmin
      atm_step_max = atmmax
      env_error_tol = enverr
      env_step_begin = envbeg
      env_step_min = envmin
      env_step_max = envmax
      pulsation_output_active = lpulse
      pulsation_file_version = ipver
      atm_choice = kttau
      debye_huckel_eta_min = etadh0
      debye_huckel_eta_max = etadh1
      use_debye_huckel_correction = ldh
      monte_carlo_file1_path = fmonte1
      monte_carlo_file2_path = fmonte2
      num_models = nmodls
      first_call_flag = lfirst
      num_runs = numrun
      max_domega_global = dtwind
      use_structure_dt_limits = lptime
      tri_delta_teffl = tridt
      tri_delta_logl = tridl
      atm_step_size = atmstp
      envelope_step_size = envstp
      envelope_generation_flag = lenvg
      helium_flash_active = lkuthe
      new_species_value = xnewcp
      rescale_species_active = lnewcp
      mixture_change_mode = isetmix
      isotope_change_mode = isetiso
      target_carbon_cno_fraction = frac_c
      target_nitrogen_cno_fraction = frac_n
      target_oxygen_cno_fraction = frac_o
      c12_to_c13_ratio = r12_13
      n14_to_n15_ratio = r14_15
      o16_to_o17_ratio = r16_17
      o16_to_o18_ratio = r16_18
      target_metal_fraction = zxmix
      initial_h2_fraction = xh2_ini
      initial_he3_fraction = xhe3_ini
      initial_li6_fraction = xli6_ini
      initial_li7_fraction = xli7_ini
      initial_be9_fraction = xbe9_ini
      initial_b10_fraction = xb10_ini
      initial_b11_fraction = xb11_ini
      difad_velocity_scale = fw
      mixing_velocity_scale = fc
      es_velocity_scale = fes
      gsf_velocity_scale = fgsf
      mu_gradient_scale = fmu
      secular_shear_velocity_scale = fss
      critical_reynolds = rcrit
      es_mixing_scale = fesc
      secular_shear_mixing_scale = fssc
      gsf_mixing_scale = fgsfc
      gsf_inhibition_mode = igsf
      min_abundance = cmin
      absolute_tolerance = abstol
      relative_tolerance = reltol
      max_burn_iterations = kemmax
      po_weight_l = poa
      po_weight_teff = pob
      po_weight_age = poc
      po_max_len_sq = pomax
      po_output_enabled = lpout
      track_file_version = itrver
      extend_core_inward = lcore
      num_core_shells_added = mcore
      core_mass_reduction_factor = fcore
      use_pure_z_table = lpurez
      lrwsh_placeholder = lrwsh
      isochrone_output_active = liso
      initial_x_array = xenv0a
      initial_z_array = zenv0a
      mixing_length_array = cmixla
      has_senv0_array = lsenv0a
      senv0_array = senv0a
      luminosity_tolerance = toll
      radius_tolerance = tolr
      zx_tolerance = tolz
      calibrate_solar_model = lcals
      calibrate_solar_zx = lcalsolzx
      target_solar_zx = calsolzx
      target_solar_age = calsolage
      use_z_ramp = lzramp
      target_luminosity_lsun = xls
      target_star_luminosity_tolerance = xlstol
      target_teff = steff
      target_radius_rsun = sr
      specify_teff_flag = lteff
      calibrate_star_flag = lcalst
      use_opal95_eos = lopale
      use_opal2001_eos = lopale01
      use_opal2006_eos = lopale06
      use_numerical_derivatives = lnumderiv
      laol_table_z1 = zlaol1
      laol_table_z2 = zlaol2
      opal_table_z1 = zopal1
      opal_table_z2 = zopal2
      opal95_single_table_z = zopal951
      alex_table_z1 = zalex1
      kurucz_table_z1 = zkur1
      kurucz_table_z2 = zkur2
      molecular_opacity_logt_min = tmolmin
      molecular_opacity_logt_max = tmolmax
      use_alex06_tables = lalex06
      use_laol89_tables = llaol89
      use_opal92_tables = lopal92
      use_opal95_tables = lopal95
      use_kurucz90_tables = lkur90
      use_alex95_tables = lalex95
      use_conductive_opacity = lcondopacp
      use_diffusion_advection_transport = ldifad
      no_am_transport_in_core = lnoj
      allard_target_feh = alatm_feh
      allard_target_alpha = alatm_alpha
      allard_use_tau100 = laltptau100
      force_solid_body_rotation = lsolid
      solid_body_mode_flag = impjmod
      use_envelope_triangle_dt = ltrist
      pmm_exponent_a = pmma
      pmm_exponent_b = pmmb
      pmm_exponent_c = pmmc
      pmm_exponent_d = pmmd
      pmm_exponent_m = pmmm
      pmm_norm_jdot = pmmjd
      pmm_norm_mdot = pmmmd
      pmm_solar_pressure = pmmsolp
      pmm_solar_omega = pmmsolw
      pmm_solar_turnover_timescale = pmmsoltau
      use_pmm_wind_law = lmwind
      scale_by_rossby_number = lrossby
      scale_by_b_field = lbscale
      wind_law_name = awind
      wind_saturation_omega = wmax
      constant_background_diffusion_coeff = codm
      use_constant_background_diffusion = lcodm
      s0_pp = s0_1_1
      s0_he3he3 = s0_3_3
      s0_he3he4 = s0_3_4
      s0_p_c12 = s0_1_12
      s0_p_c13 = s0_1_13
      s0_p_n14 = s0_1_14
      s0_p_o16 = s0_1_16
      s0_be7_electron = s0_1_be7e
      s0_be7_p = s0_1_be7p
      s0_n15_p_c12_branch = s0_1_15_c12alp
      s0_n15_p_o16_branch = s0_1_15_o16
      s0p_pp = s0p_1_1
      s0p_he3he3 = s0p_3_3
      s0p_he3he4 = s0p_3_4
      s0p_p_c12 = s0p_1_12
      s0p_p_c13 = s0p_1_13
      s0p_p_n14 = s0p_1_14
      s0p_p_o16 = s0p_1_16
      s0pp_p_c12 = s0pp_1_12
      s0pp_p_c13 = s0pp_1_13
      s0pp_p_o16 = s0pp_1_16
      s0p_be7_p = s0p_1_be7p
      s0pp_be7_p = s0pp_1_be7p
      use_scv_eos = lscv
! MHP 8/14 SUBROUTINE TO CONVERT MORE USER-FRIENDLY INPUT VARIABLES
! INTO THE VECTORS USED IN THE CODE (SUPERCEDES OLDER INPUTS)
      call remap
! MHP 06/13 Added memory of whether the choice of atmospheres has
! been changed during the run, and what the original setting was
      atm_choice_initial = kttau
      use_ttau_relation = .false.
! DBG WRITE OUT ENTIRE NAMELIST TO ISHORT
      write(short_file_unit,nml=physics)
      write(short_file_unit,nml=control)

! Post-process all CONTROL namelist vars that hold path values.
! Expand any placeholders found in the string with the value taken from a
! corresponding environment variable, if one is defined.
      call expand_value(falex06)
      call expand_value(fallard)
      call expand_value(fatm)
      call expand_value(fcondopacp)
      call expand_value(ffermi)
      call expand_value(ffirst)
      call expand_value(flast)
      call expand_value(fliv95)
      call expand_value(fmodpt)
      call expand_value(fopale06)
      call expand_value(fpatm)
      call expand_value(fpenv)
      call expand_value(fpmod)
      call expand_value(fpurez)
      call expand_value(fscomp)
      call expand_value(fscvh)
      call expand_value(fscvhe)
      call expand_value(fscvz)
      call expand_value(fshort)
      call expand_value(fsnu)
      call expand_value(fstor)
      call expand_value(ftrack)

! Create output directory as specified in the FTRACK value of CONTROL
! namelist if it doesn't already exist.
      shell_cmd = 'mkdir -p '
        ! find index of last '/' char. Use that to snip out the directory name.
      do i = len_trim(ftrack), 1, -1
          if (ftrack(i:i) .eq. '/') then
              last_slash_idx = i
              goto 1250
          endif
      end do
      1250 continue
      shell_cmd(len_trim(shell_cmd)+2:len_trim(ftrack(1:last_slash_idx))+len_trim(shell_cmd)) = ftrack(1:last_slash_idx)
      print *,"OUTPUT placed in :  ",ftrack(1:last_slash_idx)
      print *, ''
      call system(shell_cmd)


! JVS 02/11 Acoustic depth/ Asteroseismic glitch output. Puts output
! in the same directory as all other output, and names it with the
! same conventions
      if (acoustic_depth_output) then
            iclcd_placeholder = 91
            short_prefix_len=index(fshort,'short')
            fcalcad=fshort(1:short_prefix_len-1)//'calcad'
            open(unit=iclcd_placeholder, file=fcalcad, status='UNKNOWN')
            write(iclcd_placeholder,*) 'Acoustic depth calculation output file'
            write(iclcd_placeholder,*) 'age (Gyr),radius(cm),1/sound speed(s/cm),radius (CZ), 1/cs (CZ)                                         &
&       delad,gamma1,P, T, X'


!      IACAT = 92
!      FACAT=FSHORT(1:MRK-1)//'acatm'
!      OPEN(UNIT=IACAT, FILE=FACAT, STATUS='UNKNOWN')
!      WRITE(IACAT,*) 'Acoustic depth calculation output file: atmosphere integration'
!      WRITE(IACAT,*) 'age (Gyr),radius(cm),1/sound speed(s/cm),delad,gamma1,
!     * P, T, X'

!            IJLAST = 93
!            FJLAST=FSHORT(1:MRK-1)//'jlast'
!            OPEN(UNIT=IJLAST, FILE=FJLAST, STATUS='UNKNOWN')

!      IJVS = 94
!      FJVS=FSHORT(1:MRK-1)//'jvs'
!      OPEN(UNIT=IJVS, FILE=FJVS, STATUS='UNKNOWN')
!      WRITE(IJVS,*) 'The dels'

!      IJENT = 95
!      FJENT=FSHORT(1:MRK-1)//'ent'
!      OPEN(UNIT=IJENT, FILE=FJENT, STATUS='UNKNOWN')
!      WRITE(IJENT,*) 'Profiles: Age(Gyr), log(R), log(L), Log(LHe3),
!     * Conv flag, log(T), log(P), log(D)'

!      IJDEL = 96
!      FJDEL=FSHORT(1:MRK-1)//'del'
!      OPEN(UNIT=IJDEL, FILE=FJDEL, STATUS='UNKNOWN')
!      WRITE(IJDEL,*) 'Profiles: Age(Gyr), depth, del, del, del'
      endif


! JVS END
!
! G Somers 6/14, DEFINE SCALING COEFFICIENT FOR LI/BE CROSS SECTIONS
!          DEFAULT LI6 = 5.5 MeV b FROM FOWLER ET AL. 1967
!          DEFAULT LI7 = 52 keV b FROM ROLFS & KAVANAGH 1986
!           DEFAULT BE9(P,G)B10 (BE91) = 1.1 keV b INFERRED FROM FOWLER 1988
!           DEFAULT BE9(P,D)2HE4 (BE92) = 15,000 keV b FROM FOWLER ET AL. 1967
!           DEFAULT BE9(P,A)LI6 (BE93) = 15,000 keV b FROM FOWLER ET AL. 1967
      if (lxli6) then
          sli6 = xsli6/5.5d3
      endif
      if (lxli7) then
          sli7 = xsli7/5.2d1
      endif
      if (lxbe91) then
          sbe91 = xsbe91/1.1d0
      endif
      if (lxbe92) then
          sbe92 = xsbe92/1.5d4
      endif
      if (lxbe93) then
          sbe93 = xsbe93/1.5d4
      endif
! sli6/sli7/sbe91/sbe92/sbe93 must keep their NAMELIST spelling (see
! declaration above), so copy into const_lib's canonical names here.
      li6_rate_scale = sli6
      li7_rate_scale = sli7
      be9_pg_rate_scale = sbe91
      be9_pd_rate_scale = sbe92
      be9_palpha_rate_scale = sbe93
! G Somers END
! MHP 8/25 open relevant table as well as ensuring that only one is selected
!  Disable Older OPAL EOS's if a newer one is specified
      if (lopale06) then
          lopale01 = .false.
          lopale = .false.
         open(iopale, file=fopale06,status='OLD')
      endif
      if (lopale01) then
         open(iopale, file=fopale01,status='OLD')
         lopale = .false.
      else if(lopale) then
         open(iopale, file=fopale,status='OLD')
      endif
! 3/09 Disable older Alexander opacities if a newer one is specified
      if(lalex06)then
         lalex95 = .false.
         lkur90 = .false.
      endif
      if(lalex95)then
         lkur90 = .false.
      endif
! lopale/lopale01/lalex95/lkur90 can all be overridden above, after
! they were already copy-assigned into their const_lib canonical names
! (use_opal95_eos/use_opal2001_eos/use_alex95_tables/
! use_kurucz90_tables) earlier in this subroutine -- re-sync now so the
! const_lib values reflect any override from this block, not just the
! raw namelist read.
      use_opal95_eos = lopale
      use_opal2001_eos = lopale01
      use_alex95_tables = lalex95
      use_kurucz90_tables = lkur90

      open(istor,file=fstor,form='FORMATTED',status='UNKNOWN')
      rewind(istor)
! G Somers 11/14 write the new header for the .store file, if LSTORE = TRUE.
      if(lstore)then
! JvS 08/25 Added stitched interior and envelope option      
         if(lstch)then
            lphhd = .true.
         else
            write(istor,version_fmt) yrec_version_string, git_hash_string
            write(istor,1012)
         endif   
      endif
      1012 format('# Header Key',/,'# ModType    ModNum    #Shells    ', &
           &  'M/Msun    log(Teff)    log(L/Lsun)    log(M/gram)    Age/Gyr', &
           &  '    Timestep/yr    log(M_inner/gram)    log(M_outer/gram)',/, &
           &  '# JCORE  JENV  CMIXL  EOS  ATM  ALOK HIK  LPUREZ  COMPMIX', &
           &  '  LEXCOM  LDIFY  LDIFZ  LSEMIC  LOVSTC  LOVSTE  LOVSTM', &
           &  '  LROT  LINSTB  LJDOT0  LDISK  TDISK  PDISK  WMAX  LSTORE', &
           &  '  LSTATM  LSTENV  LSTMOD  LSTPHYS  LSTROT',/,'# PPI_lum', &
           &  '    PPII_lum    PPIII_lum    CNO_lum    3He_lum    ', &
           &  'OldNeu_lum    Grav_lum',/,'# ENV1-3   log(Teff)', &
           &  '    log(L/Lsun)    P_base    T_base    R_base    MatrixElements' &
           &  ,/)
! 1013 FORMAT('# JCORE  JENV  CMIXL  EOS  ATM  ALOK HIK  LPUREZ  COMPMIX',
!     1 '  LEXCOM  LDIFY  LDIFZ  LSEMIC  LOVSTC  LOVSTE  LOVSTM',
!     1 '  LROT  LINSTB  LJDOT0  LDISK  TDISK  PDISK  WMAX  LSTORE',
!     1 '  LSTATM  LSTENV  LSTMOD  LSTPHYS  LSTROT'
!     1 ,/)
! 1014     FORMAT(
!     1'MODEL SHELL MASS RADIUS LUMINOSITY PRESSURE TEMPERATURE DENSITY OMEGA ',
!     1'CONVECTIVE INTERIOR_PT H1 He4 METALS He3 C12 C13 N14 N15 O16 O17 O18 H2 Li6 Li7 ',
!     1'Be9 OPACITY GRAV DELR DEL DELAD V_CONV GAM1 HII HEII HEIII BETA ',
!     1'ETA PPI PPII PPIII CNO TRIPLE_ALPHA E_NUC E_NEU E_GRAV CP DLNRHODLNT A RP/RE FP ',
!     1'FT J/M MOMENT DEL_KE V_ES V_GSF V_SS VTOT ')
! G Somers END

      open(short_file_unit,file=fshort,form='FORMATTED',status='UNKNOWN')
      rewind(short_file_unit)
      if (ltrack) then
          open(unit=itrack,file=ftrack, form='FORMATTED', &
           &         status='UNKNOWN')
          rewind(itrack)
      endif

! SNU OUTPUT
      if(lsnu) then
          open(unit=neutrino_unit,file=fsnu, form='FORMATTED', &
           &         status='UNKNOWN')
          rewind(neutrino_unit)
      endif
      if(ldebug) then
            open(idebug,file=fdebug,form='FORMATTED', &
           &          status='UNKNOWN')
      end if
!     MHP Moved opening of FMILNE here, no need to do elsewhere
      if(lmilne)then
         open(unit=imilne,file=fmilne,form='FORMATTED', &
           &      status='UNKNOWN',access='APPEND')
      endif
!     MHP 10/02 LBNIN never set, ignore loop
!      IF (.NOT.LBNIN) THEN
         open(unit=first_unit,file=ffirst,form='FORMATTED',status='OLD')
!      END IF
      open(unit=ilast,file=flast,form='FORMATTED',status='UNKNOWN')
      open(unit=imodpt,file=fmodpt,form='FORMATTED',status='UNKNOWN')
!     OPEN ALL PULSE FILES
      if(lpulse) then
      open(opal_model_unit, file=fpmod,status='UNKNOWN',form='FORMATTED')
      open(opal_envelope_unit, file=fpenv,status='UNKNOWN',form='FORMATTED')
      open(opal_atm_unit, file=fpatm,status='UNKNOWN',form='FORMATTED')
      end if
! MHP 6/98
! MHP 8/25 Moved call from main to here for opening dynamics_unit
      if(lmonte)then
         open(unit=dynamics_unit,file=fdyn,form='FORMATTED',status='OLD')
         open(monte_carlo_unit1, file=fmonte1,status='UNKNOWN',form='FORMATTED')
         open(monte_carlo_unit2, file=fmonte2,status='UNKNOWN',form='FORMATTED')
      endif
!     MHP 8/25 Moved opening of conductive opacity and EoS tables here, to avoid complicated passages of declared variables.
      if(lcondopacp)then
         open(icondopacp,file=fcondopacp,status='OLD')
      endif
      if(liso) then
         open(isochrone_file_unit, file=fiso,status='UNKNOWN', form='FORMATTED')
      endif
      if(lsemic)then
         if(lovstc.or.lovste.or.lovstm)then
            write(short_file_unit,2)lsemic,lovste,lovstc,lovstm
      2       format(1x,'ERROR IN SUBROUTINE PARMIN'/'SEMI-CONVECTION', &
           &  ' AND OVERSHOOT FLAGS BOTH TURNED ON'/'FLAGS LSEMIC',l2, &
           &  ' OVERSHOOT - CORE,ENVELOPE,INTERMEDIATE-',3l2/'RUN STOPPED')
            stop
         endif
      endif
      write(short_file_unit,1)(chi_grid_scale(i),i=1,12),alphae,alphac,linstb,ljdot0, &
           &               alfa,fk,fw,fc,fo,fmu,rcrit
      1 format(1x,'PT TOL',12f6.3/1x,'O.S.ENV',f6.3,' O.S.CORE',f6.3, &
           &        ' LINSTB ',l1,' LJDOT ',l1,' WIND IND.',f6.3,' FK', &
           &        1pe8.2/1x,' FV',0pf5.2,' FC',f5.2,' COUPLING' &
           &        ,f6.3, ' F MU',f5.2,' RCRIT',f9.1)
      tenv = 0.5d0*(tenv0 + tenv1)
      if(lrot) then
         lnew0 = .true.
         if(walpcz.lt.-2.0d0) walpcz = -2.0d0
         if(walpcz.gt.0.0d0) walpcz = 0.0d0
! JNT 09/2025 LSOLID OVERWRITES IMPJMOD
         if(lsolid) impjmod = 1
!CCCCC OLD OR NEW WINDLAW.
         if(.not.lmwind)then
!CCCCC INSTRUCTIONS FOR THE OLD WINDLAW
!CCCCC
!CCCCC SET UP COEFFICIENTS FOR ANGULAR MOMENTUM LOSS VIA WINDS.
!CCCCC GIVEN THE INDEX ALFA, THE FORMULA FOR JDOT IS
!CCCCC JDOT = FK*2.036D33*1.452D9**ALFA*(MDOT/1.0D-14)**(1-2ALFA/3)
!CCCCC *OMEGA**(1+4ALFA/3)*(R/RSUN)**(2-ALFA)*(M/MSUN)**-ALFA/3
!CCCCC EXMD = EXPONENT OF MDOT TERM; EXW SAME FOR OMEGA;EXR FOR R;EXM FOR M.
!CCCCC FK IS A FUDGE FACTOR,SET TO 1 TO REPRODUCE THE OBSERVED SOLAR ANGULAR
!CCCCC MOMENTUM LOSS.
            one_third = 1.0d0/3.0d0
            two_thirds = 2.0d0/3.0d0
            constfactor = fk*2.036d33*1.452d9**alfa
            exmd = 1.0d0 - two_thirds*alfa
            wind_law_omega_exponent = 1.0d0 + 2.0d0*two_thirds*alfa
            exr = 2.0d0 - alfa
            exm = -one_third*alfa
!
! INSTRUCTIONS FOR THE NEW WINDLAW
!
! SET UP INDICES FOR LOSS LAW IN TERMS OF PMM A, B, C, M.
! INCLUDE A PMMD AS WELL, FOR TURNING TAUCZ and CENTRIFUG DEPENDENCE ON/OFF.
! EXCEN = EXPONENT FOR CENTRIFUGAL TERM, {K2^2/(K2^2+0.5*W^2 R^3/GM)}^PMMM
!
! AWIND = 'K97' ENFORCES KRISHNAMURTHI (1997) WIND LAW.
! AWIND = 'V13' ENFORCES VAN SADERS + PINSONNEAULT (2013) WIND LAW, AS ADAPTED
!         FROM MATT ET AL. (2008,2012).
! AWIND = 'CUS' ADOPTS THE M-D VALUES GIVEN IN THE NAMELIST.
!
!
! THEN CONVERT THE PMM VALUES INTO EXPONENTS FOR THE TORQUE CALCULATION.
! BASIC EQUATIONS, WHERE EVERY TERM IS SCALED TO SOLAR, ARE:
!
! Jdot = w * Bmag^4m * Mdot^1-2m * R^5m+2 * M^-m * Fcen^md
! Bmag = Pphot^0.5 * (w * R^c)^b * Tcz^d
! GM/R * Mdot = (Lx/Lbol * Lbol)^a -> Mdot = R * M^-1 * Lbol * w^a
! ...where Lx/Lbol is assumed to scale with w to the power of a.
!
! THIS GIVES:
! JDOT ~ w ^ 1+a-2ma+4mb
!      ~ Tcz ^ 4md
!      ~ R ^ 3+3m+4mbc
!      ~ M ^ m-1
!      ~ Lbol ^ 1-2m
!      ~ Pphot ^ 2m
!      ~ Fcen ^ md
!
! For moment, ignore GM/R term in Mdot. This makes:
! JDOT ~ R ^ 2+5m+4mbc
! JDOT ~ M ^ -m
!
         else
            if(awind.eq.'K97')then
               pmmm = alfa/3.0
! MHP 8/17 CORRECT DEFAULT FOR A = 2, NOT 0
               pmma = 2.0
               pmmb = 1.0
               pmmc = 2.0
               pmmd = 0.0
               lbscale = .false.
            elseif(awind.eq.'V13')then
               pmmm = 0.22
               pmma = 2.0
               pmmb = 1.0
               pmmc = 0.0
               pmmd = 1.0
               lbscale = .true.
            endif
            wind_law_omega_exponent   = 1.0d0 + pmma - 2.0d0*pmma*pmmm + 4.0d0*pmmm*pmmb
! G Somers 8/17 ZERO'D OUT EXTAU. TAUCZ TERM NOW COMPUTED IN
! MWIND/MCOWIND, NOT IN AMCALC.
!            EXTAU = 4.0D0*PMMM
            extau = 0.0
! JvS 09/25 REMOVED TYPO IN EXR = 2.0D0+5.0D0*PMMM-4.0D0*PMMM*PMMB*PMMC
            exr   = 2.0d0+5.0d0*pmmm-4.0d0*pmmm*pmmc
            exm   = -pmmm
            exl   = 1.0d0 - 2.0d0*pmmm
! G Somers 11/17, ADDED PMMD TO SWTICH OFF IN K97 FORM.
            expr  = 2.0d0*pmmm*pmmd
            excen = pmmm*pmmd
! INITIALIZE CONSTANT FACTOR FOR CENTRIFUGAL TERM
            c_2 = 0.0506
! SET THE CONSTANT FACTOR
            constfactor = fk*pmmjd/pmmsolw**wind_law_omega_exponent
! IF RELEVANT RESET THE SATURATION THRESHOLD IN
! TERMS OF THE SOLAR ROTATION RATE.  WMAX_SUN<1000
! INDICATES SATURATION (AT THE SUN), SO
! WSAT = WMAX_SUN*PMMSOLW
            if(wmax_sun.lt.1.0e3)then
               wmax = wmax_sun*pmmsolw
            endif
         endif
      endif
!     WINDLAW END
! impjmod/pmma/pmmb/pmmc/pmmd/pmmm/lbscale/wmax can all be overridden
! above (the LSOLID-overwrites-IMPJMOD line and the K97/V13 PMM-windlaw
! branch), after they were already copy-assigned into their const_lib
! canonical names (solid_body_mode_flag/pmm_exponent_a/pmm_exponent_b/
! pmm_exponent_c/pmm_exponent_d/pmm_exponent_m/scale_by_b_field/
! wind_saturation_omega) earlier in this subroutine -- re-sync now so
! the const_lib values reflect any override from this block, not just
! the raw namelist read.
      solid_body_mode_flag = impjmod
      pmm_exponent_a = pmma
      pmm_exponent_b = pmmb
      pmm_exponent_c = pmmc
      pmm_exponent_d = pmmd
      pmm_exponent_m = pmmm
      scale_by_b_field = lbscale
      wind_saturation_omega = wmax
!
      parmin_ln10 = dlog(10.0d0)
      if(lnewcp) then
       value_relative_to_h = .true.
       if(atmp.eq.'ABS') value_relative_to_h = .false.
! DECIDE WHICH ELEMENT IN ARRAY HCOMP TO BE RESCALED
! USING CHARACTER ARRAY AID AND INPUT CHARACTER VARIABLE ANEWCP
       do 10 i = 1,12
          if(anewcp.eq.element_id(i)) then
! INEWCP IS THE INDEX OF THE ELEMENT BEING ALTERED
            new_species_index = i + 3
            goto 30
          endif
      10    continue
! ANEWCP NOT A RECOGNIZED ELEMENT
       lnewcp = .false.
       write(short_file_unit,20) anewcp
      20    format(1x,'VARIABLE',a4,1x,'NOT A RECOGNIZED ELEMENT'/1x, &
           &    'RESCALING NOT PERFORMED')
      30    continue
      endif
      change_cno_mixture_active = .false.
      change_isotope_ratios_active = .false.
      if(isetmix.eq.1)then
! IF DEFAULT MIX (GS98) CNOFRACS ARE ALREADY SET.
!         IF(AMIX.EQ.'GS98')THEN
!            LMIXTURE = .TRUE.
!            GOTO 606
!         ENDIF
! FOR A CUSTOM MIX,DISABLE IF THE SUM OF CNO MASS FRACTIONS EXCEEDS ONE
! OR IF ANY MASS FRACTION IS NEGATIVE
         if(amix.eq.'CUS')then
            if(frac_c.lt.0.0d0.or.frac_n.lt.0.0d0.or.frac_o.lt.0.0d0)then
               write(*,591)frac_c,frac_n,frac_o
               write(short_file_unit,591)frac_c,frac_n,frac_o
      591          format('NEGATIVE INPUT CNO FRACTION ',3e12.4, &
           &              ' MIX NOT MODIFIED')
               goto 602
            endif
            sum_frac=frac_c+frac_n+frac_o
            if(sum_frac.ge.1.0d0)then
               write(*,598)frac_c,frac_n,frac_o
               write(short_file_unit,598)frac_c,frac_n,frac_o
      598          format('INPUT CNO FRACTION ',3e12.4, &
           &              ' EXCEEDS 1. MIX NOT MODIFIED')
               goto 602
            endif
! VALID MIXTURE, USE CUSTOM ENTRIES FROM .NML1
            change_cno_mixture_active = .true.
            goto 606
         endif
! SEARCH THROUGH OTHER VALID MIXTURE ENTRIES;IF FOUND,ASSIGN
!         DO I = 2,4
         do i = 1,4
            if(amix.eq.mixture_id_table(i))then
               zxmix = zx_mix_table(i)
               frac_c = frac_c_table(i)
               frac_n = frac_n_table(i)
               frac_o = frac_o_table(i)
               change_cno_mixture_active = .true.
               goto 606
            endif
         end do
!     NO VALID MIX SPECIFIED
         write(*,589)amix
         write(short_file_unit,589)amix
      589    format('DESIRED CNO MIXTURE ',a8,' NOT FOUND. MIX NOT ALTERED.')
      endif
      606 if(change_cno_mixture_active)then
         write(*,604)amix,frac_c,frac_n,frac_o
         write(short_file_unit,604)amix,frac_c,frac_n,frac_o
      604    format('CNO MIXTURE ',a8,' C ',e12.4,' N ',e12.4,' O ', &
           &         e12.4,' APPLIED TO STARTING MODEL.')
      endif
!     CHECK IF ISOTOPE RATIOS NEED TO BE ALTERED
      602 if(isetiso.eq.1)then
! FOR A CUSTOM MIX,DISABLE IF THE SUM OF CNO MASS FRACTIONS EXCEEDS ONE
! OR IF ANY MASS FRACTION IS NEGATIVE
         if(aiso.eq.'CUS')then
            if(r12_13.lt.0.0d0 .or. r14_15.lt.0.0d0 .or. r16_17.lt.0.0d0 &
           &  .or. r16_18.lt.0.0d0 .or. xh2_ini.lt.0.0d0 .or. xhe3_ini.lt.0.0d0 &
           &  .or. xli6_ini.lt.0.0d0 .or. xli7_ini.lt.0.0d0 .or. &
           &  xbe9_ini.lt.0.0d0.or.xb10_ini.lt.0.0d0 .or.xb11_ini.lt.0.0d0)then
               write(*,596)r12_13,r14_15,r16_17,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini
               write(short_file_unit,596)r12_13,r14_15,r16_17,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini
      596          format('NEGATIVE INPUT ISOTOPE RATIO OR LIGHT ELEMENT' &
           &   ' MASS FRACTION ',11e12.4, &
           &              ' MIX NOT MODIFIED')
               goto 603
            endif
            sum_frac= xh2_ini+xhe3_ini+xli6_ini+xli7_ini+xbe9_ini+ &
           &                 xb10_ini+xb11_ini
            if(sum_frac.ge.1.0d0)then
               write(*,595)xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini, &
           &  xb10_ini,xb11_ini
               write(short_file_unit,595)xh2_ini,xhe3_ini,xli6_ini,xli7_ini, &
           &  xbe9_ini,xb10_ini,xb11_ini
      595          format('SUM OF LIGHT ELEMENT MASS FRACTIONS EXCEEDS 1', &
           &  11e12.4,' MIX NOT MODIFIED')
               goto 603
            endif
!     CURRENTLY THERE ARE ONLY 2 VALID OPTIONS - THE DEFAULT (L21) OR
! A CUSTOM MIXTURE (CUS) - IF NEITHER IS TRUE, EXIT
         else if(aiso.ne.'L21')then
            goto 603
         endif
!     PASSED ALL CHECKS - EITHER THE DEFAULT OR THE CUSTOM SETTINGS WILL BE APPLIED
         change_isotope_ratios_active = .true.
      endif
      if(change_isotope_ratios_active)then
         write(*,605)aiso,r12_13,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini
         write(short_file_unit,605)aiso,r12_13,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini
      605    format('ISOTOPE AND LIGHT ELEMENT MIXTURE ',a8,' C12/C13 ', &
           &    e12.4,' O16/O18 ',e12.4,' H2 ',e12.4,' HE3 ',e12.4,' LI6 ', &
           &    e12.4,' LI7 ',e12.4,' BE9 ',e12.4,' APPLIED TO STARTING MODEL.')
      endif
! DBG 12/95 ENSURE CORRECT PARAMETERS FOR Z DIFFUSION
      603 if (ldifz) then
           ldify=.true.
         lthoul=.true.
      end if

! G SOMERS 04/15; ENSURE CORRECT PARAMETERS FOR LIGHT ELEMENT DIFFUSION.
      if(ldifli)then
          lnewdif=.true.
          ldify=.true.
          ldifz=.true.
          lthoul=.true.
          lthoulfit=.false.
          ilambda=4
      endif
! ldify/ldifz/lnewdif/lthoulfit/ilambda can all be overridden above,
! after they were already copy-assigned into their const_lib canonical
! names (diffuse_helium_active/use_diffusion_z/
! use_new_diffusion_routines/use_thoul_fit/coulomb_log_choice) earlier
! in this subroutine -- re-sync now so the const_lib values reflect
! any override from this block, not just the raw namelist read.
      diffuse_helium_active = ldify
      use_diffusion_z = ldifz
      use_new_diffusion_routines = lnewdif
      use_thoul_fit = lthoulfit
      coulomb_log_choice = ilambda
! lnewcp/frac_c/frac_n/frac_o/zxmix can likewise be overridden above
! (the ANEWCP-rescaling and CNO-mixture validation blocks), after
! already being copy-assigned into their const_lib canonical names
! earlier in this subroutine -- re-sync for the same reason as above.
      rescale_species_active = lnewcp
      target_carbon_cno_fraction = frac_c
      target_nitrogen_cno_fraction = frac_n
      target_oxygen_cno_fraction = frac_o
      target_metal_fraction = zxmix

!     WRITE OUT RUN PARAMETERS.

      write(short_file_unit,50)initial_envelope_x,initial_envelope_z
      50 format(30x,'RUN DATA VALUES'/3x, &
           &        'LINE  1     XENV0     ZENV0       ZSI'/2x, &
           &        'STANDARD       N/A       N/A  0.00E+00'/3x, &
           &        'CURRENT',1p2e10.0,2x,e8.2)
      if(npoint.le.0) npoint = 9999
      write(short_file_unit,70) ldebug,lcorr,npoint,lmilne,ltrack,lstore,lstpch
      70 format(3x,'LINE  2    LDEBUG     LCORR    NPOINT    LMILNE    LTRA                                                            &
&CK    LSTORE    LSTPCH'/2x,'STANDARD',2(9x,'T'),6x, &
           & '9999',9x,'F',2(9x,'T'),9x,'T'/3x,'CURRENT',2(9x,l1),6x, &
           & i4,3(9x,l1),9x,l1)
      if(npenv.le.0) npenv = 9999
      write(short_file_unit,90) lscrib,lstatm,lstenv,lstmod,lstphys,lstrot
      90 format(3x,'LINE  3    LSCRIB    LSTATM    LSTENV    LSTMOD',5x, &
           & 'LSTPHYS    LSTROT'/2x,'STANDARD',9x,'T',2(9x,'F'), &
           & 9x,'T',2(9x,'F'),9x,'2'/3x,'CURRENT',5(9x,l1),9x,l1)
      write(short_file_unit,110)lenvg,atmstp,envstp
      110 format(3x,'LINE 4     LENVG    ATMSTP    ENVSTP'/2x,'STANDARD',7x, &
           &        'N/A',6x,'0.50',6x,'0.50'/3x,'CURRENT',9x,l1,2(4x,f6.3))
      if(nprtmod.le.0) nprtmod = 9999
      if(nprtpt.le.0) nprtpt = 9999
      if(pulse_gyre_interval.lt.0) pulse_gyre_interval = 0
      write(short_file_unit,130) nprtmod,nprtpt
      130 format(3x,'LINE  4  NPRTMOD    NPRTPT'/2x,'STANDARD', &
           & 1(7x,'N/A'),9x,'5'/3x,'CURRENT',2(6x,i4))

!     SPIT OUT NAMELIST VARIABLES TO ISHORT

      write(short_file_unit,25) (tcut(j),j=1,5),tscut,tenv0,tenv1,tgcut
      25 format(3x,'LINE  2    TCUT-  E  TCUT- PP  TCUT-CNO  TCUT-                                                                     &
&3A  TCUT- NU TCUT-SAHA     TENV0     TENV1     TGCUT'/2x, 'STANDAR                                                            &
&D',9x,'6.50',6x,'6.50',6x,'6.82',6x,'7.70',6x,'7.50',6x, &
           & '6.00',6x,'3.00',6x,'9.00',6x,'6.90'/3x,'CURRENT', &
           & 9(5x,f5.2))
      write(short_file_unit,35) atmerr,atmmax,atmd0,enverr,envmax,envmin
! MHP 10/02 obsolete variables removed
!      WRITE(ISHORT,35) NIATM,ATMERR,ATMMAX,ATMD0,NIENV,ENVERR,ENVMAX,
!     *ENVMIN
      35 format(3x,'LINE  3     NIATM    ATMERR    ATMMAX     ATMD0',5x, &
           & 'NIENV    ENVERR    ENVMAX    ENVMIN'/2x,'STANDARD',9x,'8',2x, &
           & '3.00E-04  5.00E-01  1.00E-10',8x,'10  3.00E-04  5.00E-01  2.50E-0                                                            &
&1'/3x,'CURRENT',10x,3(1pe10.2),10x,3(1pe10.2))
!     *1'/3X,'CURRENT',7X,I3,3(1PE10.2),7X,I3,3(1PE10.2))
      write(short_file_unit,45) tridt,tridl/tridt/tridl
      45 format(3x,'LINE  4     TRIDT     TRIDL    LSENV0',5x, &
           & 'SENV0'/2x,'STANDARD',6x,'0.01',6x,'0.08',9x,'F',2x, &
           & '1.00E-07'/3x,'CURRENT',2f10.4)
      write(short_file_unit,55)(htoler(5,j),j=1,2),((htoler(i,j),i=1,4),j=1,2)
      55 format(3x,'LINE  5 TOL.RHS-P TOL.RHS-T MIN.COR-P MIN.COR-T MIN.COR                                                            &
&-R MIN.COR-L MAX.COR-P MAX.COR-T MAX.COR-R MAX.COR-L'/2x,'STANDARD                                                            &
&  3.00E-05  2.50E-06  6.00E-05  4.50E-05  3.00E-05  9.00E-05  9.00                                                            &
&E-01  5.00E-01  5.00E-01  2.00E+00'/3x,'CURRENT',10(2x,1pe8.2))
      write(short_file_unit,65)(chi_grid_scale(j),j=1,8)
      65 format(3x,'LINE  6  D(S)-MIN  D(S)-MAX   FLAG-DX   FLAG-DZ    MAX                                                             &
&DP MAX DL/LT    MAX DX    MAX DZ'/2x,'STANDARD  1.00E-08  8.00E-02                                                            &
&  5.00E-02  1.00E+00  5.00E-02  2.00E-02  1.00E+00  1.00E+00'/3x, &
           & 'CURRENT',8(2x,1pe8.2))
      write(short_file_unit,75)lnewcp,anewcp,value_relative_to_h,xnewcp
      75 format(3x,'LINE  7    LNEWCP    ANEWCP      LREL    XNEWCP'/2x, &
           & 'STANDARD',9x,'F',7x,'N/A',9x,'T',7x,'N/A'/3x,'CURRENT',9x,l1,7x, &
           & a3,9x,l1,1pe10.2)
      write(short_file_unit,85)acfpft,itfp1,itfp2
      85 format(3x,'LINE  8    ACFPFT     ITFP1     ITFP2'/2x,'STANDARD', &
           & 1x,'1.000E-20',9x,'2',8x,'20'/3x,'CURRENT',1pe10.3,6x,i4,6x,i4)
      write(short_file_unit,105) niter1,niter2,fcorr0,fcorri
      105 format(32x,'RUN DATA VALUES'/3x, &
           & 'LINE  1    NITER1    NITER2    FCORR0    FCORRI'/2x, &
           & 'STANDARD',9x,'2',8x,'20',6x,'0.80',6x,'0.10'/3x, &
           & 'CURRENT',2(6x,i4),2(5x,f5.2))
      write(short_file_unit,145)(atime(i),i=1,3),atime(7)
      145 format(3x,'LINE  5 XCORE MIN  DEL.XCORE FRAC.XCORE DEL.XSHELL'/2x, &
           & 'STANDARD     0.001     0.020     0.500     0.100'/3x,'CURRENT', &
           & 4(4x,f6.3))
      write(short_file_unit,155)(atime(i),i=4,6)
      155 format(3x,'LINE  7 DEL.YCORE  FRAC.YCORE DEL.YSHELL'/2x,'STANDARD                                                             &
&    0.020     0.300    0.0015'/3x,'CURRENT',2(5x,f5.3),4x,f6.4)
      write(short_file_unit,165) lkuthe
      165 format(3x,'LINE  8    LKUTHE'/,2x, &
           & 'STANDARD',9x,'F'/3x,'CURRENT',9x,l1)
      write(short_file_unit,175) cmixl,dpenv,lovstc,alphac,lovste,alphae
      175 format(3x,'LINE  9   CMIXL     DPENV    LOVSTC    ALPHAC    LOVSTE                                                            &
&    ALPHAE'/2x,'STANDARD',7x,'N/A',6x,'1.00',2(9x,'F',6x,'0.00')/ &
           & 3x,'CURRENT',2(5x,f5.2),2(9x,l1,5x,f5.2))
      write(short_file_unit,185) lnew0,lexcom
      185 format(3x,'LINE 10   LNEW0    LEXCOM'/2x, &
           & 'STANDARD',2(9x,'F')/3x,'CURRENT',2(9x,l1))
      write(short_file_unit,195) lrot,walpcz,linstb
      195 format(3x,'LINE 11      LROT    WALPCZ    LINSTB'/2x,'STANDARD', &
           & 7x,'N/A',6x,'0.00',7x,'N/A'/3x,'CURRENT',9x,l1,5x,f5.2,9x,l1)
      if(kttau .eq. 0) then
           write(short_file_unit, 197)
      else if (kttau .eq. 1) then
           write(short_file_unit, 198)
      else if (kttau .eq. 2) then
           write(short_file_unit, 1999)
      else if (kttau .eq. 3) then
           write(short_file_unit, 1888)
      else if (kttau .eq. 4) then
           write(short_file_unit, 1889)
! JNT 6/14 ADD FOR NEW KURUCZ/CASTELLI ATMOSPHERE TABLES
      else if (kttau .eq. 5) then
           write(short_file_unit, 1887)
      end if
      197 format(' USING EDDINGTON T-TAU RELATION.')
      198 format(' USING KRISHNA-SWAMY T-TAU RELATION.')
      1999 format(' USING HARVARD-SMITHSONIAN REFERENCE ATMOSPHERE')
      1888 format(' USING KURUCZ ATMOSPHERE TABLE')
      1889 format(' USING ALLARD ATMOSPHERE TABLE')
      1887 format(' USING KURUCZ/CASTELLI ATMOSPHERE TABLE')

! DBG PULSE
      if (lpulse) then
          write(short_file_unit,196)
      196     format(/,' CALCULATE PULSATION OUTPUT ON LAST MODEL')
      end if
      if (lpurez) then
          write(short_file_unit,*) ' USING PURE C AND N OPACITY TABLES'
      end if

      write(imodpt,310) descrip(1),  descrip(2)

      write(short_file_unit,314)
      write(short_file_unit,version_fmt) yrec_version_string, git_hash_string
      write(short_file_unit,310) descrip(1),  descrip(2)
      310    format('# DESCRIPTION OF RUN:',a80,/, '#',9x,'  ',8x,': ', &
           &           a80,/,'#', 100('='))

      if(ltrack) then

         write(itrack,314)
      314    format('#',/,'#',100('='))
         write(itrack,version_fmt) yrec_version_string, git_hash_string
         write(itrack,320) descrip(1), descrip(2)
      320    format('# DESCRIPTION OF RUN:',a80,/, '#',9x,'  ',8x,': ', &
           &           a80,/,'#', 100('='))

      endif
      write(short_file_unit,323)
      323 format(' USING OSCILATORY SPLINE INTERPOLATION IN HPOINT')

!     INTERPRET RUN FROM SEQUENCE OF "KIND" CARDS

      write(short_file_unit,200)
      200 format(/35x,'RUN CARDS'/)

      lfirst(1) = .true.

!     RUN LOOP
      do 1000 nkind=1, numrun
! READ IN NMODLS AND MODEL SOURCE(MEMORY OR FIRST MODEL)-SAME FOR ALL
       rescale_kind(nkind) = kindrn(nkind)
       if(kindrn(nkind).eq.1) then
! EVOLVE CARD
! MHP 10/24 GENERALIZE STOP CONDITIONS
!          LENDAG(NKIND) = ENDAGE(NKIND).GT.0D0
            if(endage(nkind).gt.0.0d0 .or. end_dcen(nkind).gt.0.0d0 &
           &  .or. end_xcen(nkind).gt.0.0d0 .or. end_ycen(nkind).gt.0.0d0)then
               end_age_stop_active(nkind)=.true.
            else
               end_age_stop_active(nkind)=.false.
            endif
          timestep_override_active(nkind) = setdt(nkind).gt.0d0
            if (nmodls(nkind).gt.0) then
          if (lfirst(nkind)) then
             write(iowr,350) nkind,nmodls(nkind)
             write(short_file_unit,350) nkind,nmodls(nkind)
      350          format(/1x,'RUN #',i3,'   EVOLVE ',i5, &
           &          ' MODELS, STARTING', &
           &          ' WITH THE INPUT "FIRST MODEL".')
          else
             write(iowr,351) nkind, nmodls(nkind)
             write(short_file_unit,351) nkind, nmodls(nkind)
      351          format(/1x,'RUN #',i3,'   EVOLVE ',i5, &
           &          ' MODELS, STARTING', &
           &          ' WITH THE PREVIOUS RUN''S LAST MODEL.')
          end if
! GENERALIZE STOP CONDITIONS
          if(end_age_stop_active(nkind).or.timestep_override_active(nkind)) then
             write(iowr,370)end_age_stop_active(nkind),timestep_override_active(nkind), &
           &               endage(nkind), setdt(nkind),end_dcen(nkind), &
           &          end_xcen(nkind),end_ycen(nkind)
             write(short_file_unit,370)end_age_stop_active(nkind),timestep_override_active(nkind), &
           &          endage(nkind), setdt(nkind),end_dcen(nkind), &
           &          end_xcen(nkind),end_ycen(nkind)
      370          format(1x,'EVOLVE TO AGE ',l1,' SET DELT ', &
           &               l1,' FINAL AGE ', e9.2,' FIXED TSTEP ',e9.2, &
           &   ' CENTRAL D ',e10.4,' CENTRAL X ',e12.4,' CENTRAL Y ',e12.4)
          endif
            end if
       else if(kindrn(nkind).eq.2) then
! RESCALE CARD:  RESCALE STARTING MODEL
! QUANTITIES TO BE RESCALED STORED IN ARRAY RESCALE(4,50)
! WHERE THE ELEMENTS MASS,X,Z,CORE MASS ARE STORED IN ORDER
          rescale_params(1,nkind) = rsclm(nkind)
          rescale_params(2,nkind) = rsclx(nkind)
          rescale_params(3,nkind) = rsclz(nkind)
          rescale_params(4,nkind) = rsclcm(nkind)
            if (nmodls(nkind) .gt. 0) then
          if (lfirst(nkind)) then
             write(iowr,450) nkind
             write(short_file_unit,450) nkind
      450          format(/1x,'RUN #',i3, &
           &          '   RESCALE THE INPUT MODEL: "FIRST MODEL".')
          else
             write(iowr,451) nkind
             write(short_file_unit,451) nkind
      451          format(/1x,'RUN #',i3, &
           &          '   RESCALE THE PREVIOUS RUN''S LAST MODEL.')
          end if
          write(iowr,452) nmodls(nkind),(rescale_params(i,nkind),i = 1,4)
          write(short_file_unit,452) nmodls(nkind), &
           &       (rescale_params(i,nkind),i = 1,4)
      452       format(1x,'RELAX RESCALED MODEL',i3, &
           &       ' TIMES. RESCALE THE FOLLOW', &
           &       'ING(0=USE CURRENT VALUE):'/1x,'MASS ', &
           &       f9.6,3x,'X',f9.6,3x,'Z', &
           &       f9.6,3x,'CORE MASS',f9.6)
            end if
         else if(kindrn(nkind).eq.3) then
! RESCALE AND EVOLVE CARD:  RESCALE STARTING MODEL
! QUANTITIES TO BE RESCALED STORED IN ARRAY RESCALE(4,50)
! WHERE THE ELEMENTS MASS,X,Z,CORE MASS ARE STORED IN ORDER
            rescale_params(1,nkind) = rsclm(nkind)
            rescale_params(2,nkind) = rsclx(nkind)
            rescale_params(3,nkind) = rsclz(nkind)
            rescale_params(4,nkind) = rsclcm(nkind)
            if (lfirst(nkind)) then
               write(iowr,550) nkind
               write(short_file_unit,550) nkind
      550          format(/1x,'RUN #',i3, &
           &          '   RESCALE & EVOLVE THE INPUT MODEL: "FIRST MODEL".')
            else
               write(iowr,451) nkind
               write(short_file_unit,451) nkind
!   551          FORMAT(/1X,'RUN #',I3,
!      1         '   RESCALE & EVOLVE THE PREVIOUS RUN''S LAST MODEL.')
            end if
            write(iowr,452) nmodls(nkind),(rescale_params(i,nkind),i = 1,4)
            write(short_file_unit,452) nmodls(nkind), &
           &       (rescale_params(i,nkind),i = 1,4)
       end if
         if(rescale_params(3,nkind).ge.0.0d0)  env_comp%envelope_metal_fraction=rescale_params(3,nkind)
      1000 continue
      return

end subroutine parmin


!----------------------------------------------------------------------
! expand_value
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original parmin.f (this subroutine lived in the same source file);
! only variable names, source form, and comment style were updated.
!
! Replaces any defined "{YREC_XXX}" placeholder string in the passed
! namelist path variable with the value of the corresponding
! environment variable (falling back to a built-in default path if
! that environment variable is not set).
subroutine expand_value(path_value)

      implicit none

      integer, parameter :: n_env_vars = 3  ! Number of possible env. vars.
      character(len=256) :: path_value
      character(len=256) :: temp_value
      character(len=256) :: placeholder_names(n_env_vars)
      character(len=256) :: default_values(n_env_vars)
      character(len=256) :: placeholder
      integer :: env_val_len, orig_val_len, placeholder_len, default_len
      integer :: i
! Each placeholder (env var name enclosed in curly braces) to be
! supported for expansion in namelists, along with default value in the
! case where the var is not defined in the execution environment.
! The number of assignment pairs here must match the value of the
! parameter n_env_vars, above.
! Each placeholder (env var name enclosed in curly braces) to be
! supported for expansion in namelists, along with default value in the
! case where the var is not defined in the execution environment.
! The number of assignment pairs here must match the value of the
! parameter NUM_ENVVARS, above.
      placeholder_names(1) = "{YREC_INPUT}"
      default_values(1) = "../../input"

      placeholder_names(2) = "{YREC_START}"
      default_values(2) = "../../startmodels"

      placeholder_names(3) = "{YREC_OUTPUT}"
      default_values(3) = "output"

      do 5000 i=1, n_env_vars
        placeholder_len = len_trim(placeholder_names(i))
        default_len = len_trim(default_values(i))
        placeholder = placeholder_names(i)(1:placeholder_len)

        if (path_value(1:placeholder_len) .eq. placeholder) then    ! If placeholder
            call getenv(placeholder(2:placeholder_len-1), temp_value)
            env_val_len = len_trim(temp_value)
            orig_val_len = len_trim(path_value)
            if (env_val_len .ne. 0) then       ! If env var set
                temp_value(env_val_len+1:env_val_len+orig_val_len-placeholder_len) = path_value(placeholder_len+1:orig_val_len)
                  !print *,"Override: ",TEMP(1:LEN_TRIM(TEMP))
            else                          ! No env var. Use default path.
                temp_value(1:default_len) = default_values(i)(1:len_trim(default_values(i)))
                temp_value(default_len+1:default_len+orig_val_len) = path_value(placeholder_len+1:orig_val_len)
                  !print *,"default: ",TEMP(1:LEN_TRIM(TEMP))
            end if
            path_value = temp_value
        end if
      5000 continue

end subroutine expand_value
