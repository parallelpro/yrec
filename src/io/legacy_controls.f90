!----------------------------------------------------------------------
! legacy_controls
!----------------------------------------------------------------------
! The controls reader's working state: the variables the legacy
! CONTROL and PHYSICS namelists bind to, under their original
! spellings (a namelist matches input entries by variable name, so
! renaming one would break every deck that sets it), the two namelist
! groups, the declarations of the new-style &star_job/&controls
! variables (inlist_new_decl.inc, generated), and the reader's other
! working variables. read_input copies the values into their
! canonical homes (controls_buffer, star%job) after the read;
! reset_legacy_controls clears everything at the start of each read.
module legacy_controls
      use controls_buffer
      use integration_parameters
      implicit none


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
           fshort, &
           fdyn, &
           flldat, fscomp, fkur, &
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
      character(len=256) :: control_nml_file, physics_nml_file
      character(len=256) :: shell_cmd
      integer :: i, last_slash_idx
      double precision :: one_third, two_thirds
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
! opal92_table_unit/neutrino_unit/composition_unit/kurucz_table_unit),
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

! former common/pulsegyre/: new (2026) dedicated block for the
! GYRE-format periodic pulsation-structure output feature
! (io/write_gyre_pulse.f90, triggered from io/yrec_output.f90 every
! pulse_gyre_interval converged models). pulse_gyre_interval is a
! NAMELIST /control/ value spelled identically to its const_lib
! canonical name -- use-associated directly.


! lnew0: NAMELIST /physics/ member spelled identically to its
! const_lib canonical name -- use-associated directly. tridt/tridl are
! also NAMELIST members but need different canonical spellings
! (tri_delta_teffl/tri_delta_logl), so they keep their local
! NAMELIST-spelled names here and are copy-assigned after the
! namelist read below. senv0/lsenv0 (former common/cenv/'s remaining
! members, renamed to requested_envelope_mass/
! change_envelope_mass_flag) are not namelist values and unused in
! this file -- core/read_starting_model.f90 computes them -- so they're dropped
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
! map_user_inputs.f90 recomputes chi_grid_scale from the more user-friendly
! tol_dm_min/tol_dm_max/etc namelist inputs when `call remap` runs
! below, overwriting whatever hpttol/chi_grid_scale started as -- so
! copying hpttol into chi_grid_scale here would just get clobbered
! (this file's own diagnostic writes read chi_grid_scale directly,
! after `call remap`, not hpttol, to reflect that). fcorr (former
! common/ctol/'s remaining member) is unused in this file --
! core/henyey_iteration.f90/core/main.f90 compute it at runtime -- so it's
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
! /intpar/; copied into integration_parameters's canonically-named variables after
! the namelist read below, since integration_parameters's names differ from these
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
! setup/setups.f90 computes cdh/dhnue at startup, convection's compute_scale_height
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

! lpulse: NAMELIST /physics/ member with a different
! canonical const_lib spelling (pulsation_output_active/
! pulsation_file_version), so kept local under their NAMELIST spelling
! here and copy-assigned after the namelist read below. xmsol (former
! common/pulse/'s remaining member) is unused in this file, so it's
! dropped entirely.


! track: NAMELIST /physics/ member, renamed in const_lib (itrver ->
! track_file_version), kept local and copy-assigned below.

! kttau: NAMELIST /physics/ member with a different canonical const_lib
! spelling (atm_choice), so kept local under its NAMELIST spelling
! here and copy-assigned after the namelist read below. kttau0/lttau
! (former common/atmos/'s remaining members) are not namelist values
! -- both computed from kttau right after the namelist read (kttau0 =
! kttau; lttau = .false.) -- computed working state, so they live on
! the star structure (star%atm_choice_initial/star%use_ttau_relation
! since the 2026 phase-A controls eviction).
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

! nwlaol: tollaol is a NAMELIST /physics/ member spelled identically
! to its canonical controls_buffer name -- use-associated directly (the
! table arrays are opacity_table%laol_* now; the unit numbers are
! newunit locals of the readers). lpurez is a NAMELIST /physics/ member with a
! different canonical spelling (use_pure_z_table), so kept local under
! its NAMELIST spelling here and copy-assigned after the namelist read
! below.
      logical :: lpurez

! lrwsh/liso: NAMELIST /physics/ members, each with a different
! canonical const_lib spelling (rewind_short_file/
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
! too (net_lib.f90's compute_zone_energy's neutrino_dlnq_dlnt/neutrino_dlnq_dlnd,
! which shared this block only by position, are genuinely local to
! that routine and were made plain locals there), so they're dropped
! entirely rather than carried forward.
      logical :: lnulos1

! cals2: NAMELIST /physics/ members, all renamed in const_lib (toll/
! tolr/tolz/lcals/lcalsolzx/calsolzx/calsolage -> luminosity_tolerance/
! radius_tolerance/zx_tolerance/calibrate_solar_model/
! calibrate_solar_zx/target_solar_zx/target_solar_age), kept local and
! copy-assigned below.
      double precision :: toll, tolr, tolz, calsolzx, calsolage
      logical :: lcals, lcalsolzx

! zramp: rsclzc/rsclzm1/rsclzm2/nk are spelled
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

! miscopac: kurucz_table2_unit (originally ikur2, the second Kurucz
! table's unit) and icondopacp are use-associated directly (neither is
! a namelist value). lcondopacp is a NAMELIST /physics/
! member with a different canonical spelling
! (use_conductive_opacity), so kept local under its NAMELIST spelling
! here and copy-assigned after the namelist read below.
      logical :: lcondopacp

! former common/alexo/: ialxo is not a namelist value and genuinely
! used in this file -- renamed in place to its canonical const_lib
! name (alex94_table_unit), now use-associated rather than locally
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
! /cross/. With lnewnuc = .false. sstandard becomes
! star%cross_section_scale (adopt_canonical_names, before
! map_user_inputs); with .true. map_user_inputs recomputes the scales
! from the s0_* inputs. lnewnuc is copied into use_new_nuclear_rates
! right after the namelist read below, since map_user_inputs needs it.
      double precision :: sstandard(17)
      logical :: lnewnuc

! 2026 phase-C rename campaign: these NAMELIST-spelled locals carry
! the legacy /control//physics spellings for controls whose
! controls_buffer buffers now bear their canonical names (copied into
! them in adopt_canonical_names below). Their defaults duplicate the
! buffer defaults exactly so both read paths agree, and are assigned
! at entry (see the reset block at the top of the executable part)
! rather than by initializer, so nothing is SAVEd across in-process
! re-reads.
      double precision :: alphac, alphae
      double precision :: dtdif
      double precision :: wnew, creim
      double precision :: cstmixing, cstdiffmix
      double precision :: tcut(5)
      double precision :: vnew(12)
      integer :: itdif1, itdif2
      integer :: niter1, niter2, niter3, niter4
      logical :: ldifli, lsnu, lovstc
      logical :: lwnew, ljdot0, lthoul
      logical :: lsemic, lnews, lmonte
      logical :: lreimer
      integer :: imbeg, imend

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

! former common/acdpth/ (calcad acoustic-depth machinery): retired
! 2026 with the LACOUT/.calcad output mode.

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

! lcalcenv: NAMELIST /physics/ member (must keep this exact
! spelling); copied into calc_envelope_flag after the namelist read
! below. (2026: lnewtcz retired -- the legacy taucal turnover mode
! is gone; the newer calculation is the only one.) taucz/taucz0/
! pphot/pphot0/fracstep (former common/ovrtrn/'s other five members)
! are unused in this file -- they're genuinely evolving per-model
! state read/written by many distant files, now state/turnover_lib.f90
! -- so they're dropped from this file's own declarations entirely.
      logical :: lcalcenv

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
!
!
! SPLIT NAMELIST INTO TWO: CONTROL and PHYSICS

      namelist /control/ &
           &    cmixla, calsolage, calsolzx, &
           &    descrip, &
           &    endage, &
           &    flaol, fpurez,flaol2, fopal2, &
           &    flast, ffirst, ffermi, fshort, fstch, &
           &    fdyn, flldat, fscomp, fkur, fmhd1, &
           &    fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, fmhd7, fmhd8, fiso, fatm, &
           &    fkur2, fallard, fscvh, fscvhe, fscvz, fopale, fliv95, &
           &    fmonte1,fmonte2, &
           &    kindrn, &
           &    lfirst, set_initial_model_number, &
           &    terminal_interval, report_solver_diagnostics, &
           &    log_L_upper_limit, log_L_lower_limit, &
           &    Teff_upper_limit, Teff_lower_limit, &
           &    log_g_upper_limit, log_g_lower_limit, &
           &    nu_max_upper_limit, nu_max_lower_limit, &
           &    nu_max_sun, delta_nu_sun, Teff_sun, &
           &    inlist_used_file, profile_data_prefix, &
           &    lzramp, lteff, lcalst, lpurez, &
! MHP 9/24 add LCALSOLZX to namelist
           &    liso, lrwsh, lsenv0a,lcals,lcalsolzx, &
           &    llaol89,lopal92,lopal95,lkur90,lalex95, &
           &    npenv, numrun, nmodls, &
           &    opecalex, &
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
           &    cmin, clsun, crsun, solar_mass_cgs, G_cgs, &
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
           &    lcalcenv

! 2026 inlist revamp: new-style (&star_job/&controls) machinery,
! generated from defaults/controls_registry.tsv.
      include 'inlist_new_decl.inc'


contains

! ---------------------------------------------------------------
! Every variable declared above, cleared at the start of each read: they
! were read_input's locals, which -finit-local-zero zeroed on entry
! (character ones to NUL bytes), and an in-process rerun (pyyrec,
! test_reentry) must start from the same state.
subroutine reset_legacy_controls
      version_fmt = repeat(achar(0), len(version_fmt))
      flaol = repeat(achar(0), len(flaol))
      fpurez = repeat(achar(0), len(fpurez))
      flaol2 = repeat(achar(0), len(flaol2))
      fopal2 = repeat(achar(0), len(fopal2))
      fkur2 = repeat(achar(0), len(fkur2))
      anewcp = repeat(achar(0), len(anewcp))
      atmp = repeat(achar(0), len(atmp))
      amix = repeat(achar(0), len(amix))
      aiso = repeat(achar(0), len(aiso))
      fiso = repeat(achar(0), len(fiso))
      fatm = repeat(achar(0), len(fatm))
      fstch = repeat(achar(0), len(fstch))
      fallard = repeat(achar(0), len(fallard))
      fscvh = repeat(achar(0), len(fscvh))
      fscvhe = repeat(achar(0), len(fscvhe))
      fscvz = repeat(achar(0), len(fscvz))
      flast = repeat(achar(0), len(flast))
      ffirst = repeat(achar(0), len(ffirst))
      ffermi = repeat(achar(0), len(ffermi))
      fshort = repeat(achar(0), len(fshort))
      fdyn = repeat(achar(0), len(fdyn))
      flldat = repeat(achar(0), len(flldat))
      fscomp = repeat(achar(0), len(fscomp))
      fkur = repeat(achar(0), len(fkur))
      fmhd1 = repeat(achar(0), len(fmhd1))
      fmhd2 = repeat(achar(0), len(fmhd2))
      fmhd3 = repeat(achar(0), len(fmhd3))
      fmhd4 = repeat(achar(0), len(fmhd4))
      fmhd5 = repeat(achar(0), len(fmhd5))
      fmhd6 = repeat(achar(0), len(fmhd6))
      fmhd7 = repeat(achar(0), len(fmhd7))
      fmhd8 = repeat(achar(0), len(fmhd8))
      kindrn = 0
      rsclm = 0
      rsclx = 0
      rsclz = 0
      rsclcm = 0
      fliv95 = repeat(achar(0), len(fliv95))
      fopale = repeat(achar(0), len(fopale))
      fopale01 = repeat(achar(0), len(fopale01))
      fcondopacp = repeat(achar(0), len(fcondopacp))
      fopale06 = repeat(achar(0), len(fopale06))
      opecalex = repeat(achar(0), len(opecalex))
      falex06 = repeat(achar(0), len(falex06))
      alfa = 0
      fk = 0
      wmax_sun = 0
      pmmwmax = 0
      zalex2 = 0
      zopal952 = 0
      lpmm = .false.
      lsolwind = .false.
      element_id = repeat(achar(0), len(element_id))
      mixture_id_table = repeat(achar(0), len(mixture_id_table))
      zx_mix_table = 0
      frac_c_table = 0
      frac_n_table = 0
      frac_o_table = 0
      control_nml_file = repeat(achar(0), len(control_nml_file))
      physics_nml_file = repeat(achar(0), len(physics_nml_file))
      shell_cmd = repeat(achar(0), len(shell_cmd))
      i = 0
      last_slash_idx = 0
      one_third = 0
      two_thirds = 0
      sum_frac = 0
      nkind = 0
      first_model_binary_lu = 0
      last_model_binary_lu = 0
      stored_models_binary_lu = 0
      fmonte1 = repeat(achar(0), len(fmonte1))
      fmonte2 = repeat(achar(0), len(fmonte2))
      descrip = repeat(achar(0), len(descrip))
      tridt = 0
      tridl = 0
      nmodls = 0
      numrun = 0
      lfirst = .false.
      clsun = 0
      crsun = 0
      dtwind = 0
      lptime = .false.
      hpttol = 0
      djok = 0
      lovste = .false.
      atmstp = 0
      envstp = 0
      lenvg = .false.
      lexcom = .false.
      lkuthe = .false.
      atmerr = 0
      atmd0 = 0
      atmbeg = 0
      atmmin = 0
      atmmax = 0
      enverr = 0
      envbeg = 0
      envmin = 0
      envmax = 0
      stolr0 = 0
      imax = 0
      nuse = 0
      tscut = 0
      xnewcp = 0
      lnewcp = .false.
      isetmix = 0
      isetiso = 0
      frac_c = 0
      frac_n = 0
      frac_o = 0
      r12_13 = 0
      r14_15 = 0
      r16_17 = 0
      r16_18 = 0
      zxmix = 0
      xh2_ini = 0
      xhe3_ini = 0
      xli6_ini = 0
      xli7_ini = 0
      xbe9_ini = 0
      xb10_ini = 0
      xb11_ini = 0
      optol = 0
      lrot = .false.
      linstb = .false.
      endage = 0
      setdt = 0
      end_dcen = 0
      end_xcen = 0
      end_ycen = 0
      fw = 0
      fc = 0
      fes = 0
      fgsf = 0
      fmu = 0
      fss = 0
      rcrit = 0
      etadh0 = 0
      etadh1 = 0
      ldh = .false.
      fesc = 0
      fssc = 0
      fgsfc = 0
      igsf = 0
      grtol = 0
      ilambda = 0
      niter_gs = 0
      ldify = .false.
      cmin = 0
      abstol = 0
      reltol = 0
      kemmax = 0
      dt_gs = 0
      xmin = 0
      ymin = 0
      lthoulfit = .false.
      ldifz = .false.
      lnewdif = .false.
      kttau = 0
      lmhd = .false.
      lcore = .false.
      mcore = 0
      fcore = 0
      lpurez = .false.
      lrwsh = .false.
      liso = .false.
      xenv0a = 0
      zenv0a = 0
      cmixla = 0
      senv0a = 0
      lsenv0a = .false.
      lnulos1 = .false.
      toll = 0
      tolr = 0
      tolz = 0
      calsolzx = 0
      calsolage = 0
      lcals = .false.
      lcalsolzx = .false.
      lzramp = .false.
      xls = 0
      xlstol = 0
      steff = 0
      sr = 0
      lteff = .false.
      lcalst = .false.
      lopale = .false.
      lopale01 = .false.
      lopale06 = .false.
      lnumderiv = .false.
      zlaol1 = 0
      zlaol2 = 0
      zopal1 = 0
      zopal2 = 0
      zopal951 = 0
      zalex1 = 0
      zkur1 = 0
      zkur2 = 0
      tmolmin = 0
      tmolmax = 0
      lalex06 = .false.
      llaol89 = .false.
      lopal92 = .false.
      lopal95 = .false.
      lkur90 = .false.
      lalex95 = .false.
      lcondopacp = .false.
      ldifad = .false.
      lnoj = .false.
      sstandard = 0
      lnewnuc = .false.
      alphac = 0
      alphae = 0
      dtdif = 0
      wnew = 0
      creim = 0
      cstmixing = 0
      cstdiffmix = 0
      tcut = 0
      vnew = 0
      itdif1 = 0
      itdif2 = 0
      niter1 = 0
      niter2 = 0
      niter3 = 0
      niter4 = 0
      ldifli = .false.
      lsnu = .false.
      lovstc = .false.
      lwnew = .false.
      ljdot0 = .false.
      lthoul = .false.
      lsemic = .false.
      lnews = .false.
      lmonte = .false.
      lreimer = .false.
      imbeg = 0
      imend = 0
      s0_1_1 = 0
      s0_3_3 = 0
      s0_3_4 = 0
      s0_1_12 = 0
      s0_1_13 = 0
      s0_1_14 = 0
      s0_1_16 = 0
      s0_1_be7e = 0
      s0_1_be7p = 0
      s0_1_15_c12alp = 0
      s0_1_15_o16 = 0
      s0p_1_1 = 0
      s0p_3_3 = 0
      s0p_3_4 = 0
      s0p_1_12 = 0
      s0p_1_13 = 0
      s0p_1_14 = 0
      s0p_1_16 = 0
      s0pp_1_12 = 0
      s0pp_1_13 = 0
      s0pp_1_16 = 0
      s0p_1_be7p = 0
      s0pp_1_be7p = 0
      lscv = .false.
      alatm_feh = 0
      alatm_alpha = 0
      laltptau100 = .false.
      tdisk = 0
      pdisk = 0
      ldisk = .false.
      weakscreening = 0
      lsolid = .false.
      impjmod = 0
      dmdt0 = 0
      compacc = 0
      lmdot = .false.
      ltrist = .false.
      pmma = 0
      pmmb = 0
      pmmc = 0
      pmmd = 0
      pmmm = 0
      pmmjd = 0
      pmmmd = 0
      pmmsolp = 0
      pmmsolw = 0
      pmmsoltau = 0
      lmwind = .false.
      lrossby = .false.
      lbscale = .false.
      awind = repeat(achar(0), len(awind))
      wmax = 0
      lcalcenv = .false.
      codm = 0
      lcodm = .false.
      xsli6 = 0
      xsli7 = 0
      xsbe91 = 0
      xsbe92 = 0
      xsbe93 = 0
      lxli6 = .false.
      lxli7 = .false.
      lxbe91 = .false.
      lxbe92 = .false.
      lxbe93 = .false.
      sli6 = 0
      sli7 = 0
      sbe91 = 0
      sbe92 = 0
      sbe93 = 0
      spotf = 0
      spotx = 0
      lsdepth = .false.
end subroutine reset_legacy_controls

end module legacy_controls
