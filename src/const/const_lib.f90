!----------------------------------------------------------------------
! const_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/const1/, common/const2/, and
! common/const3/: physical/mixing-length constants that are set once
! (by setup/setups.f90, at run startup) and read broadly across the
! codebase, never varying per call. This is global configuration, not
! per-call data, so it becomes a module of plain (non-parameter)
! module-level variables rather than subroutine arguments -- matching
! MESA's own const_def/chem_def convention. Every file that used to
! declare any of these three COMMON blocks now does `use const_lib`
! instead; setup/setups.f90's existing assignment statements
! (unchanged) now set these module variables directly instead of the
! old COMMON slots.
!
! const1 members: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
! const2 members: gas_constant, radiation_constant_over_3, ca3l, csig,
!                 csigl, cgl, cmkh, cmkhn
! const3 members: cdelrl, cmixl, cmixl2, cmixl3, clndp,
!                 seconds_per_year -- cmixl (the mixing length) is the
!                 one member here that isn't a pure physical constant:
!                 it's copied from the per-kind-card namelist array
!                 cmixla(nk) each time a new kind card starts, rather
!                 than computed once at the very start of the run like
!                 the rest of const1-3 -- still "occasional
!                 configuration read broadly," not per-call data, so
!                 the same module treatment applies.
! ctlim members: atime, tcut, saha_log10t_cutoff, tenv0, tenv1, tenv,
!                tgcut -- NAMELIST /physics/ values (core/parmin.f90),
!                same "set once, read broadly" character as const1-3;
!                tenv is the one computed member (0.5*(tenv0+tenv1)),
!                still only once, not per call.
module const_lib
      implicit none

! former common/const1/
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! former common/const2/
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn

! former common/const3/. cmixl's default (1.4d0) was previously set by
! a DATA statement in core/parmin.f90 (data lkuthe,cmixl/.false.,
! 1.4d0/) -- moved here since parmin.f90 can no longer target a
! use-associated variable with DATA; main.f90 overwrites it with the
! per-kind-card mixing length (cmixla(nk)) before it's ever read for a
! real model, same as before.
      double precision :: cdelrl, cmixl2, cmixl3, clndp, seconds_per_year
      double precision :: cmixl = 1.4d0

! former common/ctlim/. Defaults (previously two DATA statements in
! core/parmin.f90, now illegal there since these are use-associated
! rather than locally declared) moved here as declaration-time
! initializers. tenv is the one member computed at runtime (see
! core/parmin.f90: tenv = 0.5d0*(tenv0+tenv1)), so it has no default.
      double precision :: atime(14) = (/1.0d-3,2.0d-2,5.0d-1,2.0d-2, &
           3.0d-1,1.5d-3,1.0d-1,2.0d-2,4.0d-2,2.0d-2,2.0d-2,0.25d0, &
           1.5d0,0.25d0/)
      double precision :: tcut(5) = (/6.5d0,6.5d0,6.82d0,7.7d0,7.5d0/)
      double precision :: saha_log10t_cutoff = 6.0d0
      double precision :: tenv0 = 3.0d0, tenv1 = 9.0d0, tgcut = 6.9d0
      double precision :: tenv

! former common/cross/, common/weak/. All 6 cross members and
! weak_screening_threshold get their real, final values from
! setup/remap.f90 (called from core/parmin.f90 right after the
! NAMELIST /physics/ read) -- cross_section_scale/qs0e_scale/
! qqs0ee_scale/o16_gamma_scale/c12_alpha_scale are fully recomputed
! there from other inputs regardless of any namelist value; only
! use_new_nuclear_rates and weak_screening_threshold are read
! (unchanged) from the namelist itself (lnewnuc/weakscreening in
! core/parmin.f90, which can't be renamed -- see the naming note at
! that file's top -- so parmin.f90 copies them into these canonical
! names right after the namelist read, before remap runs).
      double precision :: cross_section_scale(17), qs0e_scale(8), &
           qqs0ee_scale(8), o16_gamma_scale, c12_alpha_scale
      logical :: use_new_nuclear_rates
      double precision :: weak_screening_threshold

! former common/dpmix/: overshoot/semiconvection mixing-length
! parameters and on/off flags, all NAMELIST /physics/ values. iov1/
! iov2/iovim have no default here, matching before (no DATA statement
! for them in core/parmin.f90 either -- COMMON left them at whatever
! the loader zero-filled).
      double precision :: dpenv = 1.0d0
      logical :: lovstc = .false.
      double precision :: alphac = 0.0d0
      logical :: envelope_overshoot_active = .false.
      double precision :: alphae = 0.0d0
      logical :: lovstm = .false.
      double precision :: alpham = 0.0d0
      logical :: ladov = .false.
      logical :: lovmax = .false.
      double precision :: betac = 0.15d0
      logical :: lsemic = .false.
      integer :: iov1, iov2, iovim

! former common/rot/: rotation control parameters/flags, all
! NAMELIST /physics/ values. walpcz is clamped once at startup
! (core/parmin.f90: if(walpcz.lt.-2.0d0) walpcz=-2.0d0, etc.) --
! still configuration, not per-call data.
      double precision :: acfpft = 1.0d-36
      integer :: itfp1 = 5, itfp2 = 20
      logical :: rotation_active = .false.
      double precision :: walpcz = 0.0d0
      logical :: instability_transport_active = .false., lwnew = .false.
      double precision :: wnew = 0.0d0

! former common/masschg/: mass-accretion/Reimers-wind parameters, all
! NAMELIST /physics/ values. use_mass_accretion can be flipped off
! mid-run by wind/mdot.f90 (if a circumstellar disk is exhausted) --
! still a shared global flag touched from one distant point, same
! character as walpcz's startup clamp, just possibly more than once.
      double precision :: mass_accretion_rate = -1.0d-14
      double precision :: fczdmdt = 0.1d0, ftotdmdt = 1.0d-2
      double precision :: accreted_composition(15) = (/0.71668d0, &
           0.265721d0,0.01757d0,2.9d-5,3.013d-3,3.385d-5,9.346d-4, &
           0.0d0,8.462d-3,0.0d0,1.696d-5,0.0d0,2.0d-9,2.0d-9,3.0d-11/)
      double precision :: creim = -4.0d-13
      logical :: lreimer = .false., use_mass_accretion = .false.

! former common/neweng/: all 3 members are NAMELIST /physics/ values,
! already canonically spelled in core/parmin.f90 (no rename needed
! there).
      integer :: niter4 = 0
      logical :: lnews = .false., lsnu = .false.

! former common/burnscs/: light-element cross-section scale factors.
! core/parmin.f90's own local names (sli6 etc) are themselves NAMELIST
! /physics/ members (see that file's naming note at the top) and can't
! be renamed, so parmin.f90 keeps them local, computes them from the
! real namelist inputs (xsli6 etc, held in its own internal
! common/xsect/) via e.g. sli6 = xsli6/5.5d3, and then copy-assigns
! into these canonical names.
      double precision :: li6_rate_scale, li7_rate_scale, &
           be9_pg_rate_scale, be9_pd_rate_scale, be9_palpha_rate_scale

! former common/nuloss/'s one config member: switch selecting the Itoh
! 1996 neutrino-loss routines (nuclear/engeb.f90), a NAMELIST /physics/
! value. core/parmin.f90's local name for it (lnulos1) is itself
! namelist-visible and can't be renamed, so it stays local there and
! is copied into this canonical name after the namelist read -- same
! treatment as li6_rate_scale etc above. common/nuloss/'s other two
! members (dsnudt/dsnudd in core/parmin.f90, neutrino_dlnq_dlnt/
! neutrino_dlnq_dlnd in nuclear/engeb.f90) were dead/purely-local and
! dropped/delocalized rather than moved here -- see those files.
      logical :: use_itoh_neutrino_loss = .false.

! former common/ovrtrn/'s two config members: NAMELIST /physics/
! values selecting the newer convective-turnover-timescale calculation
! and whether to run the full envelope integration. Same
! namelist-can't-rename treatment as use_itoh_neutrino_loss above --
! core/parmin.f90 keeps its local lnewtcz/lcalcenv and copies into
! these canonical names. The other five former common/ovrtrn/ members
! are genuinely evolving per-model state, not configuration -- see
! state/turnover_lib.f90.
      logical :: use_new_turnover_timescale = .true.
      logical :: calc_envelope_flag = .true.

! former common/const/: solar physical constants (luminosity, mass,
! radius in cgs, plus their log10/ln/bolometric-magnitude
! derivatives), all computed once at startup in setup/setups.f90 --
! same "set once, read broadly" character as const1-3. Two of the
! eight (solar_luminosity_cgs, solar_radius_cgs) start from
! core/parmin.f90's NAMELIST /physics/ values (clsun/crsun, which
! can't be renamed -- see that file's naming note -- so parmin.f90
! copies them into these canonical names right after the namelist
! read, before setups.f90 runs); the rest are pure derived quantities
! computed by setups.f90 from those two. No declaration-time
! initializers here since every member is set at runtime before first
! use.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude

! former common/flag/: single NAMELIST /physics/ member selecting the
! extended (15-species) vs. default composition array. Same
! namelist-can't-rename treatment as use_itoh_neutrino_loss etc above
! -- core/parmin.f90 keeps its local lexcom and copies into this
! canonical name.
      logical :: use_extended_composition = .false.

! former common/spots/: all 3 members are NAMELIST /physics/ values,
! same namelist-can't-rename treatment as above -- core/parmin.f90
! keeps its local spotf/spotx/lsdepth and copies into these canonical
! names.
      double precision :: spot_filling_factor = 0.00d0
      double precision :: spot_temp_contrast = 1.00d0
      logical :: spot_depth_varies = .false.

! former common/disk/: disk_temperature/disk_pressure/
! disk_locking_active are NAMELIST /physics/ values (core/parmin.f90's
! tdisk/pdisk/ldisk, kept local there and copy-assigned). disk_lifetime
! (former common/disk/'s remaining member, originally "sage") is not a
! namelist value -- core/main.f90 sets it at runtime -- so it has no
! declaration-time default here.
      double precision :: disk_temperature = 0.0d0
      double precision :: disk_pressure = 7.2722d-6
      logical :: disk_locking_active = .false.
      double precision :: disk_lifetime

! former common/sett/: target_end_age/timestep_override/
! central_deuterium_stop/central_hydrogen_stop/central_helium_stop are
! NAMELIST /physics/ values (core/parmin.f90's endage/setdt/end_dcen/
! end_xcen/end_ycen, kept local there and copy-assigned).
! end_age_stop_active/timestep_override_active (former common/sett/'s
! remaining two members, originally lendag/lsetdt) are not namelist
! values -- core/parmin.f90 computes them from the above -- so they
! have no declaration-time default here.
      double precision :: target_end_age(50) = 0.0d0
      double precision :: timestep_override(50) = 0.0d0
      double precision :: central_deuterium_stop(50) = 0.0d0
      double precision :: central_hydrogen_stop(50) = 0.0d0
      double precision :: central_helium_stop(50) = 0.0d0
      logical :: end_age_stop_active(50), timestep_override_active(50)

! former common/mhd/: use_mhd_eos is a NAMELIST /physics/ value
! (core/parmin.f90's lmhd, kept local there and copy-assigned).
! unit_zams_a/b/c/unit_centre1-5 (former common/mhd/'s remaining
! members, originally iomhd1-8) are not namelist values --
! core/parmin.f90 unconditionally assigns them fixed unit numbers at
! runtime -- so they have no declaration-time default here.
      logical :: use_mhd_eos = .false.
      integer :: unit_zams_a, unit_zams_b, unit_zams_c
      integer :: unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

! former common/optab/: metal_fraction_match_tolerance is a NAMELIST
! /physics/ value (core/parmin.f90's optol, kept local there and
! copy-assigned). zsi/idt/idd (former common/optab/'s remaining
! members) are not namelist values -- dead in core/parmin.f90 (dropped
! there) but genuinely set-and-consumed-locally in several other
! files (misc/coefft.f90, misc/physic.f90, atm/envint.f90,
! core/starin.f90), each independently assigning the same constants
! (idt=15, idd(:)=5) -- kept here rather than deleted since removing
! the assignment would be a logic change, not a mechanical conversion.
      double precision :: metal_fraction_match_tolerance
      double precision :: zsi = 0.0d0
      integer :: idt, idd(4)

! former common/ccout/: lstore/lstatm/lstenv/lstmod/lstphys/lstrot/
! lscrib/lstch/lphhd are all NAMELIST /physics/ values in
! core/parmin.f90 (lphhd is set directly there, not namelist-read, but
! still lives in this same former block) -- all 9 keep their original
! COMMON member spelling as their canonical const_lib name (no
! separate readable rename ever established for this block), so
! core/parmin.f90 use-associates them directly rather than keeping
! separate locals.
      logical :: lstore = .false., lstatm, lstenv, lstmod, lstphys
      logical :: lstrot, lscrib = .true., lstch = .false., lphhd

! former common/ccout1/: npenv/nprtmod/npoint are NAMELIST /physics/
! values, spelled identically to their canonical names (use-associated
! directly in core/parmin.f90). print_point_interval (originally
! nprtpt) is also a NAMELIST value but needed a different, more
! readable name established elsewhere, so core/parmin.f90 keeps a
! local nprtpt and copy-assigns.
      integer :: npenv, nprtmod = 1, npoint = 1
      integer :: print_point_interval

! former common/ccout2/: ldebug/lcorr/lmilne/ltrack/lstpch are all
! NAMELIST /physics/ values, spelled identically to their canonical
! names -- same treatment as common/ccout/ above.
      logical :: ldebug = .false., lcorr = .true., lmilne = .false.
      logical :: ltrack = .true., lstpch = .false.

! former common/lunum/: logical unit numbers for the various input/
! output files, none of them NAMELIST values -- core/parmin.f90
! assigns them all unconditionally at startup (fixed unit numbers),
! so no declaration-time defaults are needed here.
      integer :: first_unit, run_unit, standard_unit, fermi_unit
      integer :: opal_model_unit, opal_envelope_unit, opal_atm_unit
      integer :: dynamics_unit, laol_table_unit, neutrino_unit
      integer :: composition_unit, kurucz_table_unit

! former common/monte/: lmonte/imbeg/imend are all NAMELIST /physics/
! values spelled identically to their canonical names -- same
! treatment as common/ccout/ above.
      logical :: lmonte = .false.
      integer :: imbeg = 1, imend = 1

! former common/ctol/: htoler/fcorr0/fcorri/niter1/niter2/niter3 are
! NAMELIST /physics/ values spelled identically to their canonical
! names -- use-associated directly in core/parmin.f90. chi_grid_scale
! (originally hpttol) is also a NAMELIST value but needed a different,
! more readable name established elsewhere, so core/parmin.f90 keeps
! a local hpttol and copy-assigns (no declaration-time default here,
! since it's always set that way before first use). fcorr (former
! common/ctol/'s remaining member) is not a namelist value --
! core/crrect.f90/core/main.f90 compute it at runtime -- so it has no
! declaration-time default either.
      double precision :: htoler(5,2) = reshape((/6.0d-5,4.5d-5,3.0d-5, &
           9.0d-5,3.0d-5,9.0d-1,5.0d-1,5.0d-1,2.0d0,2.5d-6/), (/5,2/))
      double precision :: fcorr0 = 0.8d0, fcorri = 0.1d0, fcorr
      double precision :: chi_grid_scale(12)
      integer :: niter1 = 2, niter2 = 20, niter3 = 2

! former common/difus/: dtdif/itdif1/itdif2 are NAMELIST /physics/
! values spelled identically to their canonical names -- use-
! associated directly. convergence_tolerance (originally djok) is also
! a NAMELIST value but needed a different, more readable name
! established elsewhere, so core/parmin.f90 keeps a local djok and
! copy-assigns.
      double precision :: dtdif = 1.0d-2
      double precision :: convergence_tolerance
      integer :: itdif1 = 1, itdif2 = 1

! former common/gravst/: all 4 members are NAMELIST /physics/ values,
! each with a different canonical spelling than core/parmin.f90's
! terse names (grtol/ilambda/niter_gs/ldify), so core/parmin.f90 keeps
! those local and copy-assigns.
      double precision :: settling_tolerance
      integer :: coulomb_log_choice, settling_num_iterations
      logical :: diffuse_helium_active

! former common/gravs2/: all 4 members are NAMELIST /physics/ values,
! each with a different canonical spelling than core/parmin.f90's
! terse names (dt_gs/xmin/ymin/lthoulfit), so core/parmin.f90 keeps
! those local and copy-assigns.
      double precision :: settling_timestep_fraction, &
           hydrogen_diffusion_floor, helium_diffusion_min
      logical :: use_thoul_fit

! former common/gravs3/: fgry/fgrz/lthoul are NAMELIST /physics/
! values spelled identically to their canonical names -- use-
! associated directly. use_diffusion_z (originally ldifz) is also a
! NAMELIST value but needed a different, more readable name, so
! core/parmin.f90 keeps a local ldifz and copy-assigns.
      double precision :: fgry = 1.0d0, fgrz = 1.0d0
      logical :: lthoul = .false.
      logical :: use_diffusion_z

! former common/gravs4/: use_new_diffusion_routines (originally
! lnewdif) is a NAMELIST value with a different canonical spelling, so
! core/parmin.f90 keeps a local lnewdif and copy-assigns. ldifli is a
! NAMELIST value spelled identically to its canonical name -- use-
! associated directly.
      logical :: use_new_diffusion_routines
      logical :: ldifli = .false.

! former common/intatm/: all 5 members are NAMELIST /physics/ values
! with different canonical spellings than core/parmin.f90's terse
! names (atmerr/atmd0/atmbeg/atmmin/atmmax -> atm_error_tol/
! atm_step_initial/atm_step_begin/atm_step_min/atm_step_max), all kept
! local there and copy-assigned.
      double precision :: atm_error_tol, atm_step_initial, atm_step_begin
      double precision :: atm_step_min, atm_step_max

! former common/intenv/: all 4 members are NAMELIST /physics/ values
! with different canonical spellings than core/parmin.f90's terse
! names (enverr/envbeg/envmin/envmax -> env_error_tol/env_step_begin/
! env_step_min/env_step_max), all kept local there and copy-assigned.
      double precision :: env_error_tol, env_step_begin, env_step_min
      double precision :: env_step_max

! former common/pulse/: pulsation_output_active/pulsation_file_version
! (originally lpulse/ipver) are NAMELIST values with different
! canonical spellings, kept local in core/parmin.f90 and copy-assigned.
! pulsation_mass_msun (originally xmsol) is not a namelist value --
! core/main.f90 sets it at runtime (from total_mass_msun) -- so it has
! no declaration-time default; it was unused in core/parmin.f90 itself
! and dropped there.
      logical :: pulsation_output_active
      integer :: pulsation_file_version
      double precision :: pulsation_mass_msun

! former common/atmos/: atm_choice (originally kttau) is a NAMELIST
! value with a different canonical spelling, kept local in
! core/parmin.f90 and copy-assigned. atm_choice_initial/
! use_ttau_relation (originally kttau0/lttau) are not namelist values
! -- core/parmin.f90 computes them from kttau right after the
! namelist read (atm_choice_initial = kttau; use_ttau_relation =
! .false.) -- so no declaration-time defaults are needed for them.
! atm_hras (former common/atmos/'s remaining member, originally hras)
! is not a namelist value -- setup/setups.f90 computes it at runtime
! -- so it has no declaration-time default; it was unused in
! core/parmin.f90 itself and dropped there.
      integer :: atm_choice
      integer :: atm_choice_initial
      logical :: use_ttau_relation
      double precision :: atm_hras

! former common/debhu/: Debye-Huckel EOS correction data. No readable
! rename had ever been established for this block (every file used
! the original cryptic COMMON member spelling, unlike most other
! blocks), so these canonical names are new (2026), not picked up from
! an existing convention. use_debye_huckel_correction/
! debye_huckel_eta_min/debye_huckel_eta_max (originally ldh/etadh0/
! etadh1) are NAMELIST /physics/ values, kept local in
! core/parmin.f90 and copy-assigned. debye_huckel_coefficient/
! debye_huckel_nu (originally cdh/dhnue) are not namelist values --
! setup/setups.f90 computes them once at startup -- so they have no
! declaration-time default. debye_huckel_x/debye_huckel_y/
! debye_huckel_z_total/debye_huckel_z (originally xxdh(or xxdy)/yydh/
! zzdh/zdh) are likewise not namelist values -- mixing/hsubp.f90
! recomputes them per shell from the local composition -- so they too
! have no declaration-time default.
      double precision :: debye_huckel_coefficient
      double precision :: debye_huckel_eta_min, debye_huckel_eta_max
      double precision :: debye_huckel_z(18)
      double precision :: debye_huckel_x, debye_huckel_y, &
           debye_huckel_z_total
      double precision :: debye_huckel_nu(18)
      logical :: use_debye_huckel_correction

! former common/vnewcb/: vnew (initial abundances for a 12-species
! set, Na/Al/Mg/Fe/Si/C/H/O/N/Ar/Ne/He for a G&N93 solar mixture) is a
! NAMELIST /physics/ value spelled identically to its canonical name.
      double precision :: vnew(12) = (/0.001999d0, 0.003238d0, &
           0.037573d0, 0.071794d0, 0.040520d0, 0.173285d0, 0.000000d0, &
           0.482273d0, 0.053152d0, 0.005379d0, 0.098668d0, 0.000000d0/)

! former common/iomonte/: monte_carlo_file1_path/monte_carlo_file2_path
! (originally fmonte1/fmonte2) are NAMELIST values with different
! canonical spellings, kept local in core/parmin.f90 and copy-assigned
! -- no declaration-time default (character strings, no DATA default
! in the original either). monte_carlo_unit1/monte_carlo_unit2
! (originally imonte1/imonte2) are not namelist values -- core/parmin.f90
! assigns them fixed unit numbers (70/71) unconditionally at startup
! -- so no declaration-time default is needed for them either.
      character(len=256) :: monte_carlo_file1_path, monte_carlo_file2_path
      integer :: monte_carlo_unit1, monte_carlo_unit2

! former common/pulsegyre/: pulse_gyre_interval is a NAMELIST
! /control/ value spelled identically to its canonical name.
      integer :: pulse_gyre_interval = 0

! former common/cenv/: lnew0 is a NAMELIST /physics/ value spelled
! identically to its canonical name. tri_delta_teffl/tri_delta_logl
! (originally tridt/tridl) are also NAMELIST values but need different
! canonical spellings, so core/parmin.f90 keeps tridt/tridl local
! (with their own DATA defaults there) and copy-assigns -- no
! declaration-time default needed here. requested_envelope_mass/
! change_envelope_mass_flag (originally senv0/lsenv0) are not namelist
! values and unused in core/parmin.f90 -- core/starin.f90 computes
! them -- so no declaration-time default here either.
      double precision :: tri_delta_teffl, tri_delta_logl
      logical :: lnew0 = .false.
      double precision :: requested_envelope_mass
      logical :: change_envelope_mass_flag

! former common/ckind/: rescale_params/rescale_kind (originally rescal/
! iresca) are not namelist values -- core/parmin.f90 computes them at
! runtime (rescale_params from rsclm/rsclx/rsclz/rsclcm, rescale_kind
! from kindrn) -- so no declaration-time default. num_models/
! first_call_flag/num_runs (originally nmodls/lfirst/numrun) are
! NAMELIST values with different canonical spellings, kept local in
! core/parmin.f90 (with their own DATA defaults there) and
! copy-assigned, so likewise no default needed here.
      double precision :: rescale_params(4,50)
      integer :: rescale_kind(50)
      integer :: num_models(50), num_runs
      logical :: first_call_flag(50)

! former common/ct2/: max_domega_global (originally dtwind) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/parmin.f90 and copy-assigned.
      double precision :: max_domega_global

! former common/ct3/: use_structure_dt_limits (originally lptime) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/parmin.f90 and copy-assigned.
      logical :: use_structure_dt_limits

! former common/envgen/: all 3 members are NAMELIST values with
! different canonical spellings, kept local in core/parmin.f90 and
! copy-assigned.
      double precision :: atm_step_size, envelope_step_size
      logical :: envelope_generation_flag

! former common/heflsh/: helium_flash_active (originally lkuthe) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/parmin.f90 and copy-assigned.
      logical :: helium_flash_active

! former common/label/: initial_envelope_x/initial_envelope_z
! (originally xenv0/zenv0) are not namelist values -- genuinely used
! in core/parmin.f90, renamed in place there -- so their DATA defaults
! moved here as declaration-time initializers.
      double precision :: initial_envelope_x = 0.7d0
      double precision :: initial_envelope_z = 0.02d0

! former common/newcmp/: new_species_value/rescale_species_active
! (originally xnewcp/lnewcp) are NAMELIST values with different
! canonical spellings, kept local in core/parmin.f90 and copy-assigned
! (rescale_species_active also re-synced after the ANEWCP-rescaling
! validation block may override it). new_species_index/
! value_relative_to_h (originally inewcp/lrel) are not namelist values
! -- genuinely used in core/parmin.f90, renamed in place there -- no
! declaration-time default (none in the original either).
      double precision :: new_species_value
      integer :: new_species_index
      logical :: rescale_species_active
      logical :: value_relative_to_h

! former common/newmx/: mixture_change_mode/isotope_change_mode/
! target_carbon_cno_fraction/target_nitrogen_cno_fraction/
! target_oxygen_cno_fraction/c12_to_c13_ratio/n14_to_n15_ratio/
! o16_to_o17_ratio/o16_to_o18_ratio/target_metal_fraction/
! initial_h2_fraction/initial_he3_fraction/initial_li6_fraction/
! initial_li7_fraction/initial_be9_fraction/initial_b10_fraction/
! initial_b11_fraction (originally isetmix/isetiso/frac_c/frac_n/
! frac_o/r12_13/r14_15/r16_17/r16_18/zxmix/xh2_ini/xhe3_ini/xli6_ini/
! xli7_ini/xbe9_ini/xb10_ini/xb11_ini) are NAMELIST values with
! different canonical spellings, kept local in core/parmin.f90 and
! copy-assigned (the CNO-fraction/metal-fraction four are also
! re-synced after the CNO-mixture validation block may override
! them). change_cno_mixture_active/change_isotope_ratios_active
! (originally lmixture/lisotope) are not namelist values -- set by
! core/parmin.f90's own CNO/isotope validation logic, renamed in place
! there -- so no declaration-time default for those two.
      integer :: mixture_change_mode, isotope_change_mode
      logical :: change_cno_mixture_active, change_isotope_ratios_active
      double precision :: target_carbon_cno_fraction, &
           target_nitrogen_cno_fraction, target_oxygen_cno_fraction
      double precision :: c12_to_c13_ratio, n14_to_n15_ratio, &
           o16_to_o17_ratio, o16_to_o18_ratio
      double precision :: target_metal_fraction
      double precision :: initial_h2_fraction, initial_he3_fraction, &
           initial_li6_fraction, initial_li7_fraction, &
           initial_be9_fraction, initial_b10_fraction, &
           initial_b11_fraction

! former common/vmult/: difad_velocity_scale/mixing_velocity_scale/
! es_velocity_scale/gsf_velocity_scale/mu_gradient_scale/
! secular_shear_velocity_scale/critical_reynolds (originally fw/fc/
! fes/fgsf/fmu/fss/rcrit) are NAMELIST values with different canonical
! spellings, kept local in core/parmin.f90 and copy-assigned. fo is a
! NAMELIST value spelled identically to its canonical name --
! use-associated directly; its DATA default moved here since DATA can
! no longer target a use-associated entity.
      double precision :: fo = 1.0d0
      double precision :: difad_velocity_scale, mixing_velocity_scale, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds

! former common/vmult2/: es_mixing_scale/secular_shear_mixing_scale/
! gsf_mixing_scale/gsf_inhibition_mode (originally fesc/fssc/fgsfc/
! igsf) are NAMELIST values with different canonical spellings, kept
! local in core/parmin.f90 and copy-assigned. ies/imu are NAMELIST
! values spelled identically to their canonical names --
! use-associated directly; their DATA defaults moved here since DATA
! can no longer target use-associated entities.
      integer :: ies = 1
      integer :: imu = 1
      double precision :: es_mixing_scale, secular_shear_mixing_scale, &
           gsf_mixing_scale
      integer :: gsf_inhibition_mode

! former common/burtol/: min_abundance/absolute_tolerance/
! relative_tolerance/max_burn_iterations (originally cmin/abstol/
! reltol/kemmax) are NAMELIST values with different canonical
! spellings, kept local in core/parmin.f90 and copy-assigned.
      double precision :: min_abundance, absolute_tolerance, &
           relative_tolerance
      integer :: max_burn_iterations

! former common/lopal95/: opal95_table_unit (originally iliv95) is not
! a namelist value -- genuinely used in core/parmin.f90, renamed in
! place there.
      integer :: opal95_table_unit

! former common/po/: po_weight_l/po_weight_teff/po_weight_age/
! po_max_len_sq/po_output_enabled (originally poa/pob/poc/pomax/
! lpout) are NAMELIST values with different canonical spellings, kept
! local in core/parmin.f90 and copy-assigned.
      double precision :: po_weight_l, po_weight_teff, po_weight_age, &
           po_max_len_sq
      logical :: po_output_enabled

! former common/track/: track_file_version (originally itrver) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/parmin.f90 and copy-assigned.
      integer :: track_file_version

! former common/core/: extend_core_inward/num_core_shells_added/
! core_mass_reduction_factor (originally lcore/mcore/fcore) are
! NAMELIST values with different canonical spellings, kept local in
! core/parmin.f90 and copy-assigned.
      logical :: extend_core_inward
      integer :: num_core_shells_added
      double precision :: core_mass_reduction_factor

! former common/nwlaol/: the LAOL pure-Z opacity table (olaol/oxa/ot/
! orho/tollaol/iolaol/numofxyz/numrho/numt/llaol/iopurez) is spelled
! identically to its canonical name everywhere -- use-associated
! directly. tollaol/llaol's DATA defaults moved here from
! core/parmin.f90 since DATA can no longer target use-associated
! entities. use_pure_z_table (originally lpurez) is a NAMELIST value
! with a different canonical spelling, kept local in core/parmin.f90
! and copy-assigned.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104)
      double precision :: tollaol = 10.0d0
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol = .false.
      logical :: use_pure_z_table

! former common/chrone/: lrwsh_placeholder/isochrone_output_active
! (originally lrwsh/liso) are NAMELIST values with different canonical
! spellings, kept local in core/parmin.f90 and copy-assigned.
! isochrone_file_unit (originally iiso) is not a namelist value --
! genuinely used in core/parmin.f90, renamed in place there.
      logical :: lrwsh_placeholder, isochrone_output_active
      integer :: isochrone_file_unit

! former common/newxym/: initial_x_array/initial_z_array/
! mixing_length_array/has_senv0_array/senv0_array (originally xenv0a/
! zenv0a/cmixla/lsenv0a/senv0a) are NAMELIST values with different
! canonical spellings, kept local in core/parmin.f90 and copy-assigned.
      double precision :: initial_x_array(50), initial_z_array(50), &
           mixing_length_array(50), senv0_array(50)
      logical :: has_senv0_array(50)

! former common/atmos2/: the Kurucz surface-pressure table
! (kurucz_log10_pressure_table/kurucz_teff_table/kurucz_logg_table/
! kurucz_table_z, originally atmpl/atmtl/atmgl/atmz) and
! atm_table_file_unit (originally ioatm) are spelled identically to
! their canonical names -- use-associated directly.
      double precision :: kurucz_log10_pressure_table(57,11), &
           kurucz_teff_table(57), kurucz_logg_table(11), kurucz_table_z
      integer :: atm_table_file_unit

! former common/atmos2c/: the Kurucz/Castelli surface-pressure table
! (kurucz_castelli_log10_pressure_table/kurucz_castelli_teff_table/
! kurucz_castelli_logg_table, originally atmplc/atmtlc/atmglc) is
! spelled identically to its canonical name everywhere -- use-associated
! directly. Unused in core/parmin.f90.
      double precision :: kurucz_castelli_log10_pressure_table(76,11), &
           kurucz_castelli_teff_table(76), kurucz_castelli_logg_table(11)

! former common/cals2/: luminosity_tolerance/radius_tolerance/
! zx_tolerance/calibrate_solar_model/calibrate_solar_zx/
! target_solar_zx/target_solar_age (originally toll/tolr/tolz/lcals/
! lcalsolzx/calsolzx/calsolage) are NAMELIST values with different
! canonical spellings, kept local in core/parmin.f90 and copy-assigned.
! luminosity_tolerance is the winner of a name collision with
! common/calstar/'s xlstol (disambiguated there as
! target_star_luminosity_tolerance).
      double precision :: luminosity_tolerance, radius_tolerance, &
           zx_tolerance, target_solar_zx, target_solar_age
      logical :: calibrate_solar_model, calibrate_solar_zx

! former common/zramp/: rsclzc/rsclzm1/rsclzm2/iolaol2/ioopal2/nk are
! spelled identically to their canonical names everywhere --
! use-associated directly; their DATA defaults moved here from
! core/parmin.f90 since DATA can no longer target use-associated
! entities. use_z_ramp (originally lzramp) is a NAMELIST value with a
! different canonical spelling, kept local in core/parmin.f90 and
! copy-assigned.
      double precision :: rsclzc(50) = -1.0d0, rsclzm1(50) = -1.0d0, &
           rsclzm2(50) = -1.0d0
      integer :: iolaol2, ioopal2, nk
      logical :: use_z_ramp

! former common/calstar/: target_luminosity_lsun/
! target_star_luminosity_tolerance/target_teff/target_radius_rsun/
! specify_teff_flag/calibrate_star_flag (originally xls/xlstol/steff/
! sr/lteff/lcalst) are NAMELIST values with different canonical
! spellings, kept local in core/parmin.f90 and copy-assigned.
! target_star_luminosity_tolerance is disambiguated from
! common/cals2/'s own luminosity_tolerance member (see that block's
! note above). log_l_prev_model/log_r_prev_model/age_at_target_radius/
! log_l_at_target_radius/log_l_at_target_radius_prev_run/
! age_prev_model/star_found_flag/just_passed_target_radius_flag
! (former common/calstar/'s remaining members) are unused in
! core/parmin.f90 -- genuinely used in misc/chkscal.f90/
! setup/setscal.f90/core/main.f90, so still declared here.
      double precision :: target_luminosity_lsun, &
           target_star_luminosity_tolerance, target_teff, &
           target_radius_rsun, log_l_prev_model, log_r_prev_model, &
           age_at_target_radius, log_l_at_target_radius, &
           log_l_at_target_radius_prev_run, age_prev_model
      logical :: star_found_flag, specify_teff_flag, &
           just_passed_target_radius_flag, calibrate_star_flag

end module const_lib
