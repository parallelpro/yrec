!----------------------------------------------------------------------
! controls_lib
!----------------------------------------------------------------------
! Added 2026 (phase six, step 3 -- ROADMAP.md). The run controls split
! out of const_lib: essentially every NAMELIST /control/ and /physics/
! target (core/read_input.f90 reads into these), grouped by their former
! COMMON blocks with the original per-block commentary retained.
! 2026 phase D: RELOCATED to io/ and PARMIN-PRIVATE. The only
! legitimate users are the read path -- io/read_input.f90 (namelist
! targets), setup/map_user_inputs.f90 (reads+writes the buffer), the generated
! state/controls_sync_lib.f90 (seed/store), and net/test/test_net.f90
! (sets remap inputs). Everything else reads star%ctrl / star%job.
! The const_lib umbrella is gone.
!
! Contract: read-only after read_controls/parmin. Phase A of the
! 2026 controls->star% campaign evicted every non-namelist straggler
! the original COMMON blocks mixed in (ctlim's tenv, atmos's
! kttau0/lttau/hras, the cross/weak scale arrays, the solar octet,
! cmixl, nk -- flat star% members now; iolaol2/ioopal2 to luout_lib).
!
! Phase B: this module is now formally the namelist read BUFFER.
! read_controls re-seeds it from pristine star%ctrl defaults before
! every parmin read and stores it back into star%ctrl (the
! authoritative home) afterwards -- see io/read_controls.f90 and the
! generated state/controls_state_def.inc / controls_sync_lib.f90
! (tools/gen_controls_state.py; regenerate on any member change
! here). Consumers migrate to star%ctrl%... in phase C.

module controls_lib
      implicit none

! 2026 MESA-style output: .true. -> the historical per-model streams
! (.track/.store); .false. -> MESA-format output (history.data-layout
! .history file; profiles to follow). Compile default .true. keeps
! every legacy deck byte-pinned; parmin's new-style inlist path
! flips the default to .false. before the read (either format can set
! it explicitly).
      logical, public :: use_legacy_output = .true.
! MESA-style output controls (new-format inlists only):
      character(len=256), public :: star_history_name = 'history.data'
      character(len=256), public :: history_columns_file = ' '
      character(len=256), public :: profile_columns_file = ' '
      integer, public :: profile_interval = 50
      character(len=8), public :: pulse_format = 'GYRE'
! Independent toggles: profiles only, pulse only, both, or neither
! (cadence stays with profile_interval / pulse_gyre_interval):
      logical, public :: write_profile_flag = .false.
      logical, public :: write_pulse_flag = .false.
! diagnostic: per-zone solar neutrino production table (engeb per
! shell) written to the short/log stream at the start of each run
      logical, public :: compute_neutrino_fluxes = .false.

! former common/ctlim/. Defaults (previously two DATA statements in
! core/read_input.f90, now illegal there since these are use-associated
! rather than locally declared) moved here as declaration-time
! initializers. tenv, the one member computed at runtime
! (0.5*(tenv0+tenv1)), was evicted to star%tenv (2026 phase A).
      double precision :: atime(14) = (/1.0d-3,2.0d-2,5.0d-1,2.0d-2, &
           3.0d-1,1.5d-3,1.0d-1,2.0d-2,4.0d-2,2.0d-2,2.0d-2,0.25d0, &
           1.5d0,0.25d0/)
      double precision :: tcut(5) = (/6.5d0,6.5d0,6.82d0,7.7d0,7.5d0/)
      double precision :: saha_log10t_cutoff = 6.0d0
      double precision :: tenv0 = 3.0d0, tenv1 = 9.0d0, tgcut = 6.9d0

! former common/cross/, common/weak/. Only use_new_nuclear_rates and
! weak_screening_threshold are read (unchanged) from the NAMELIST
! /physics/ itself (lnewnuc/weakscreening in core/read_input.f90, which
! can't be renamed -- see the naming note at that file's top -- so
! read_input.f90 copies them into these canonical names right after the
! namelist read, before remap runs). The five scale members
! (cross_section_scale/qs0e_scale/qqs0ee_scale/o16_gamma_scale/
! c12_alpha_scale) were fully recomputed by setup/map_user_inputs.f90 from
! other inputs regardless of any namelist value -- working state,
! evicted to flat star% members (2026 phase A).
      logical :: use_new_nuclear_rates
      double precision :: weak_screening_threshold

! former common/dpmix/: overshoot/semiconvection mixing-length
! parameters and on/off flags, all NAMELIST /physics/ values. iov1/
! iov2/iovim have no default here, matching before (no DATA statement
! for them in core/read_input.f90 either -- COMMON left them at whatever
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

! former common/rot/: rotation control parameters/flags, all
! NAMELIST /physics/ values. walpcz is clamped once at startup
! (core/read_input.f90: if(walpcz.lt.-2.0d0) walpcz=-2.0d0, etc.) --
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
! already canonically spelled in core/read_input.f90 (no rename needed
! there).
      integer :: niter4 = 0
      logical :: lnews = .false., lsnu = .false.

! former common/burnscs/: light-element cross-section scale factors.
! core/read_input.f90's own local names (sli6 etc) are themselves NAMELIST
! /physics/ members (see that file's naming note at the top) and can't
! be renamed, so read_input.f90 keeps them local, computes them from the
! real namelist inputs (xsli6 etc, held in its own internal
! common/xsect/) via e.g. sli6 = xsli6/5.5d3, and then copy-assigns
! into these canonical names.
      double precision :: li6_rate_scale, li7_rate_scale, &
           be9_pg_rate_scale, be9_pd_rate_scale, be9_palpha_rate_scale

! former common/nuloss/'s one config member: switch selecting the Itoh
! 1996 neutrino-loss routines (net_lib.f90's engeb), a NAMELIST
! /physics/ value. core/read_input.f90's local name for it (lnulos1) is
! itself namelist-visible and can't be renamed, so it stays local there
! and is copied into this canonical name after the namelist read --
! same treatment as li6_rate_scale etc above. common/nuloss/'s other
! two members (dsnudt/dsnudd in core/read_input.f90, neutrino_dlnq_dlnt/
! neutrino_dlnq_dlnd in net_lib.f90's engeb) were dead/purely-local
! and dropped/delocalized rather than moved here -- see those files.
      logical :: use_itoh_neutrino_loss = .false.

! former common/ovrtrn/'s two config members: NAMELIST /physics/
! values selecting the newer convective-turnover-timescale calculation
! and whether to run the full envelope integration. Same
! namelist-can't-rename treatment as use_itoh_neutrino_loss above --
! core/read_input.f90 keeps its local lcalcenv and copies into the
! canonical name. (2026: use_new_turnover_timescale/LNEWTCZ retired
! along with the legacy taucal turnover mode.) The other five former
! common/ovrtrn/ members are genuinely evolving per-model state, not
! configuration -- see state/turnover_lib.f90.
      logical :: calc_envelope_flag = .true.

! former common/const/: solar physical constants (luminosity, mass,
! radius in cgs, plus their log10/ln/bolometric-magnitude
! derivatives): computed at startup by setup/setups.f90 from the
! NAMELIST /physics/ clsun/crsun values (parmin copies those in) and
! overwritten per Monte-Carlo run -- computed state, evicted to flat
! star% members (2026 phase A batch 5).

! former common/flag/: single NAMELIST /physics/ member selecting the
! extended (15-species) vs. default composition array. Same
! namelist-can't-rename treatment as use_itoh_neutrino_loss etc above
! -- core/read_input.f90 keeps its local lexcom and copies into this
! canonical name.
      logical :: use_extended_composition = .false.

! former common/spots/: all 3 members are NAMELIST /physics/ values,
! same namelist-can't-rename treatment as above -- core/read_input.f90
! keeps its local spotf/spotx/lsdepth and copies into these canonical
! names.
      double precision :: spot_filling_factor = 0.00d0
      double precision :: spot_temp_contrast = 1.00d0
      logical :: spot_depth_varies = .false.

! former common/disk/: disk_locking_age_gyr/disk_omega_rad_s/
! disk_locking_active are NAMELIST /physics/ values (core/read_input.f90's
! tdisk/pdisk/ldisk, kept local there and copy-assigned). disk_lifetime
! (former common/disk/'s remaining member, originally "sage") is not a
! namelist value -- core/main.f90 sets it at runtime -- so it has no
! declaration-time default here.
      double precision :: disk_locking_age_gyr = 0.0d0
      double precision :: disk_omega_rad_s = 7.2722d-6
      logical :: disk_locking_active = .false.

! former common/sett/: target_end_age/timestep_override/
! central_deuterium_stop/central_hydrogen_stop/central_helium_stop are
! NAMELIST /physics/ values (core/read_input.f90's endage/setdt/end_dcen/
! end_xcen/end_ycen, kept local there and copy-assigned).
! end_age_stop_active/timestep_override_active (former common/sett/'s
! remaining two members, originally lendag/lsetdt) are not namelist
! values -- core/read_input.f90 computes them from the above -- so they
! have no declaration-time default here.
      double precision :: target_end_age(50) = 0.0d0
      double precision :: timestep_override(50) = 0.0d0
      double precision :: central_deuterium_stop(50) = 0.0d0
      double precision :: central_hydrogen_stop(50) = 0.0d0
      double precision :: central_helium_stop(50) = 0.0d0
      logical :: end_age_stop_active(50), timestep_override_active(50)

! former common/mhd/: use_mhd_eos is a NAMELIST /physics/ value
! (core/read_input.f90's lmhd, kept local there and copy-assigned).
! unit_zams_a/b/c/unit_centre1-5 (former common/mhd/'s remaining
! members, originally iomhd1-8) are not namelist values --
! core/read_input.f90 unconditionally assigns them fixed unit numbers at
! runtime -- so they have no declaration-time default here.
      logical :: use_mhd_eos = .false.
      integer :: unit_zams_a, unit_zams_b, unit_zams_c
      integer :: unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

! former common/optab/: metal_fraction_match_tolerance is a NAMELIST
! /physics/ value (core/read_input.f90's optol, kept local there and
! copy-assigned). zsi/idt/idd (former common/optab/'s remaining
! members) are not namelist values -- dead in core/read_input.f90 (dropped
! there) but genuinely set-and-consumed-locally in several other
! files. DELETED 2026 phase D (user-approved dead cleanup): every
! site only ASSIGNED the same constants (idt=15, idd(:)=5); nothing
! anywhere read them -- dead stores since the COMMON era.
      double precision :: metal_fraction_match_tolerance
      double precision :: zsi = 0.0d0

! former common/ccout/: lstore/lstatm/lstenv/lstmod/lstphys/lstrot/
! lscrib/lstch/lphhd are all NAMELIST /physics/ values in
! core/read_input.f90 (lphhd is set directly there, not namelist-read, but
! still lives in this same former block) -- all 9 keep their original
! COMMON member spelling as their canonical const_lib name (no
! separate readable rename ever established for this block), so
! core/read_input.f90 use-associates them directly rather than keeping
! separate locals.
      logical :: lstore = .false., lstatm, lstenv, lstmod, lstphys
      logical :: lstrot, lscrib = .true., lstch = .false., lphhd

! former common/ccout1/: npenv/nprtmod/npoint are NAMELIST /physics/
! values, spelled identically to their canonical names (use-associated
! directly in core/read_input.f90). print_point_interval (originally
! nprtpt) is also a NAMELIST value but needed a different, more
! readable name established elsewhere, so core/read_input.f90 keeps a
! local nprtpt and copy-assigns.
      integer :: npenv, nprtmod = 1, npoint = 1
      integer :: print_point_interval

! former common/ccout2/: ldebug/lcorr/lmilne/ltrack/lstpch are all
! NAMELIST /physics/ values, spelled identically to their canonical
! names -- same treatment as common/ccout/ above.
      logical :: ldebug = .false., lcorr = .true., lmilne = .false.
      logical :: ltrack = .true., lstpch = .false.

! former common/lunum/: logical unit numbers for the various input/
! output files, none of them NAMELIST values -- core/read_input.f90
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
! names -- use-associated directly in core/read_input.f90. chi_grid_scale
! (originally hpttol) is also a NAMELIST value but needed a different,
! more readable name established elsewhere, so core/read_input.f90 keeps
! a local hpttol and copy-assigns (no declaration-time default here,
! since it's always set that way before first use). fcorr (former
! common/ctol/'s remaining member) is not a namelist value --
! core/henyey_iterate.f90/core/main.f90 compute it at runtime -- so it has no
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
! established elsewhere, so core/read_input.f90 keeps a local djok and
! copy-assigns.
      double precision :: dtdif = 1.0d-2
      double precision :: convergence_tolerance
      integer :: itdif1 = 1, itdif2 = 1

! former common/gravst/: all 4 members are NAMELIST /physics/ values,
! each with a different canonical spelling than core/read_input.f90's
! terse names (grtol/ilambda/niter_gs/ldify), so core/read_input.f90 keeps
! those local and copy-assigns.
      double precision :: settling_tolerance
      integer :: coulomb_log_choice, settling_num_iterations
      logical :: diffuse_helium_active

! former common/gravs2/: all 4 members are NAMELIST /physics/ values,
! each with a different canonical spelling than core/read_input.f90's
! terse names (dt_gs/xmin/ymin/lthoulfit), so core/read_input.f90 keeps
! those local and copy-assigns.
      double precision :: settling_timestep_fraction, &
           hydrogen_diffusion_floor, helium_diffusion_min
      logical :: use_thoul_fit

! former common/gravs3/: fgry/fgrz/lthoul are NAMELIST /physics/
! values spelled identically to their canonical names -- use-
! associated directly. use_diffusion_z (originally ldifz) is also a
! NAMELIST value but needed a different, more readable name, so
! core/read_input.f90 keeps a local ldifz and copy-assigns.
      double precision :: fgry = 1.0d0, fgrz = 1.0d0
      logical :: lthoul = .false.
      logical :: use_diffusion_z

! former common/gravs4/: use_new_diffusion_routines (originally
! lnewdif) is a NAMELIST value with a different canonical spelling, so
! core/read_input.f90 keeps a local lnewdif and copy-assigns. ldifli is a
! NAMELIST value spelled identically to its canonical name -- use-
! associated directly.
      logical :: use_new_diffusion_routines
      logical :: ldifli = .false.

! former common/intatm/: all 5 members are NAMELIST /physics/ values
! with different canonical spellings than core/read_input.f90's terse
! names (atmerr/atmd0/atmbeg/atmmin/atmmax -> atm_error_tol/
! atm_step_initial/atm_step_begin/atm_step_min/atm_step_max), all kept
! local there and copy-assigned.
      double precision :: atm_error_tol, atm_step_initial, atm_step_begin
      double precision :: atm_step_min, atm_step_max

! former common/intenv/: all 4 members are NAMELIST /physics/ values
! with different canonical spellings than core/read_input.f90's terse
! names (enverr/envbeg/envmin/envmax -> env_error_tol/env_step_begin/
! env_step_min/env_step_max), all kept local there and copy-assigned.
      double precision :: env_error_tol, env_step_begin, env_step_min
      double precision :: env_step_max

! former common/pulse/: pulsation_output_active/pulsation_file_version
! (originally lpulse/ipver) are NAMELIST values with different
! canonical spellings, kept local in core/read_input.f90 and copy-assigned.
! pulsation_mass_msun (originally xmsol) is not a namelist value --
! core/main.f90 sets it at runtime (from total_mass_msun) -- so it has
! no declaration-time default; it was unused in core/read_input.f90 itself
! and dropped there.
      logical :: pulsation_output_active
      integer :: pulsation_file_version

! former common/atmos/: atm_choice (originally kttau) is a NAMELIST
! value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned. NOTE: atm_choice is clobbered
! at runtime (surfbc's hot-edge fallback sets it from
! star%atm_choice_initial; run_yrec's calibration cycles restore it)
! -- when controls move to star%ctrl it will need a working copy on
! star%, like cmixl. The block's other members (kttau0/lttau/hras ->
! atm_choice_initial/use_ttau_relation/atm_hras) were computed or
! runtime-toggled, never namelist values -- evicted to flat star%
! members (2026 phase A).
      integer :: atm_choice

! former common/vnewcb/: vnew (initial abundances for a 12-species
! set, Na/Al/Mg/Fe/Si/C/H/O/N/Ar/Ne/He for a G&N93 solar mixture) is a
! NAMELIST /physics/ value spelled identically to its canonical name.
      double precision :: vnew(12) = (/0.001999d0, 0.003238d0, &
           0.037573d0, 0.071794d0, 0.040520d0, 0.173285d0, 0.000000d0, &
           0.482273d0, 0.053152d0, 0.005379d0, 0.098668d0, 0.000000d0/)

! former common/iomonte/: monte_carlo_file1_path/monte_carlo_file2_path
! (originally fmonte1/fmonte2) are NAMELIST values with different
! canonical spellings, kept local in core/read_input.f90 and copy-assigned
! -- no declaration-time default (character strings, no DATA default
! in the original either). monte_carlo_unit1/monte_carlo_unit2
! (originally imonte1/imonte2) are not namelist values -- core/read_input.f90
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
! canonical spellings, so core/read_input.f90 keeps tridt/tridl local
! (with their own DATA defaults there) and copy-assigns -- no
! declaration-time default needed here. requested_envelope_mass/
! change_envelope_mass_flag (originally senv0/lsenv0) are not namelist
! values and unused in core/read_input.f90 -- core/read_starting_model.f90 computes
! them -- so no declaration-time default here either.
      double precision :: tri_delta_teffl, tri_delta_logl
      logical :: lnew0 = .false.
      double precision :: requested_envelope_mass
      logical :: change_envelope_mass_flag

! former common/ckind/: rescale_params/rescale_kind (originally rescal/
! iresca) are not namelist values -- core/read_input.f90 computes them at
! runtime (rescale_params from rsclm/rsclx/rsclz/rsclcm, rescale_kind
! from kindrn) -- so no declaration-time default. num_models/
! first_call_flag/num_runs (originally nmodls/lfirst/numrun) are
! NAMELIST values with different canonical spellings, kept local in
! core/read_input.f90 (with their own DATA defaults there) and
! copy-assigned, so likewise no default needed here.
      double precision :: rescale_params(4,50)
      integer :: rescale_kind(50)
      integer :: num_models(50), num_runs
      logical :: first_call_flag(50)

! former common/ct2/: max_domega_global (originally dtwind) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned.
      double precision :: max_domega_global

! former common/ct3/: use_structure_dt_limits (originally lptime) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned.
      logical :: use_structure_dt_limits

! former common/envgen/: all 3 members are NAMELIST values with
! different canonical spellings, kept local in core/read_input.f90 and
! copy-assigned.
      double precision :: atm_step_size, envelope_step_size
      logical :: envelope_generation_flag

! former common/heflsh/: helium_flash_active (originally lkuthe) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned.
      logical :: helium_flash_active

! former common/label/: initial_envelope_x/initial_envelope_z
! (originally xenv0/zenv0) are not namelist values -- genuinely used
! in core/read_input.f90, renamed in place there -- so their DATA defaults
! moved here as declaration-time initializers.
      double precision :: initial_envelope_x = 0.7d0
      double precision :: initial_envelope_z = 0.02d0

! former common/newcmp/: new_species_value/rescale_species_active
! (originally xnewcp/lnewcp) are NAMELIST values with different
! canonical spellings, kept local in core/read_input.f90 and copy-assigned
! (rescale_species_active also re-synced after the ANEWCP-rescaling
! validation block may override it). new_species_index/
! value_relative_to_h (originally inewcp/lrel) are not namelist values
! -- genuinely used in core/read_input.f90, renamed in place there -- no
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
! different canonical spellings, kept local in core/read_input.f90 and
! copy-assigned (the CNO-fraction/metal-fraction four are also
! re-synced after the CNO-mixture validation block may override
! them). change_cno_mixture_active/change_isotope_ratios_active
! (originally lmixture/lisotope) are not namelist values -- set by
! core/read_input.f90's own CNO/isotope validation logic, renamed in place
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
! spellings, kept local in core/read_input.f90 and copy-assigned. fo is a
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
! local in core/read_input.f90 and copy-assigned. ies/imu are NAMELIST
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
! spellings, kept local in core/read_input.f90 and copy-assigned.
      double precision :: min_abundance, absolute_tolerance, &
           relative_tolerance
      integer :: max_burn_iterations

! former common/lopal95/: opal95_table_unit (originally iliv95) is not
! a namelist value -- genuinely used in core/read_input.f90, renamed in
! place there.
      integer :: opal95_table_unit

! former common/po/: po_weight_l/po_weight_teff/po_weight_age/
! po_max_len_sq/po_output_enabled (originally poa/pob/poc/pomax/
! lpout) are NAMELIST values with different canonical spellings, kept
! local in core/read_input.f90 and copy-assigned.
      double precision :: po_weight_l, po_weight_teff, po_weight_age, &
           po_max_len_sq
      logical :: po_output_enabled

! former common/track/: track_file_version (originally itrver) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned.
      integer :: track_file_version

! former common/core/: extend_core_inward/num_core_shells_added/
! core_mass_reduction_factor (originally lcore/mcore/fcore) are
! NAMELIST values with different canonical spellings, kept local in
! core/read_input.f90 and copy-assigned.
      logical :: extend_core_inward
      integer :: num_core_shells_added
      double precision :: core_mass_reduction_factor

! former common/chrone/: rewind_short_file/isochrone_output_active
! (originally lrwsh/liso) are NAMELIST values with different canonical
! spellings, kept local in core/read_input.f90 and copy-assigned.
! isochrone_file_unit (originally iiso) is not a namelist value --
! genuinely used in core/read_input.f90, renamed in place there.
      logical :: rewind_short_file, isochrone_output_active
      integer :: isochrone_file_unit

! former common/newxym/: initial_x_array/initial_z_array/
! mixing_length_array/has_senv0_array/senv0_array (originally xenv0a/
! zenv0a/cmixla/lsenv0a/senv0a) are NAMELIST values with different
! canonical spellings, kept local in core/read_input.f90 and copy-assigned.
      double precision :: initial_x_array(50), initial_z_array(50), &
           mixing_length_array(50), senv0_array(50)
      logical :: has_senv0_array(50)

! former common/cals2/: luminosity_tolerance/radius_tolerance/
! zx_tolerance/calibrate_solar_model/calibrate_solar_zx/
! target_solar_zx/target_solar_age (originally toll/tolr/tolz/lcals/
! lcalsolzx/calsolzx/calsolage) are NAMELIST values with different
! canonical spellings, kept local in core/read_input.f90 and copy-assigned.
! luminosity_tolerance is the winner of a name collision with
! common/calstar/'s xlstol (disambiguated there as
! target_star_luminosity_tolerance).
      double precision :: luminosity_tolerance, radius_tolerance, &
           zx_tolerance, target_solar_zx, target_solar_age
      logical :: calibrate_solar_model, calibrate_solar_zx

! former common/zramp/: rsclzc/rsclzm1/rsclzm2 are spelled
! identically to their canonical names everywhere -- use-associated
! directly; their DATA defaults moved here from core/read_input.f90
! since DATA can no longer target use-associated entities.
! use_z_ramp (originally lzramp) is a NAMELIST value with a
! different canonical spelling, kept local in core/read_input.f90 and
! copy-assigned. Evicted 2026 phase A: iolaol2/ioopal2 to luout_lib
! (parmin-assigned units), nk to star%job%nk (the run-list cursor).
      double precision :: rsclzc(50) = -1.0d0, rsclzm1(50) = -1.0d0, &
           rsclzm2(50) = -1.0d0
      logical :: use_z_ramp

! former common/calstar/: target_luminosity_lsun/
! target_star_luminosity_tolerance/target_teff/target_radius_rsun/
! specify_teff_flag/calibrate_star_flag (originally xls/xlstol/steff/
! sr/lteff/lcalst) are NAMELIST values with different canonical
! spellings, kept local in core/read_input.f90 and copy-assigned.
! target_star_luminosity_tolerance is disambiguated from
! common/cals2/'s own luminosity_tolerance member (see that block's
! note above). log_l_prev_model/log_r_prev_model/age_at_target_radius/
! log_l_at_target_radius/log_l_at_target_radius_prev_run/
! age_prev_model/star_found_flag/just_passed_target_radius_flag
! (former common/calstar/'s remaining members) are unused in
! core/read_input.f90 -- genuinely used in misc/check_star_calibration.f90/
! setup/setup_star_calibration.f90/core/main.f90, so still declared here.
      double precision :: target_luminosity_lsun, target_star_luminosity_tolerance, target_teff, target_radius_rsun
      logical :: specify_teff_flag, calibrate_star_flag

! former common/opaleos/: use_opal95_eos/use_opal2001_eos/
! use_opal2006_eos/use_numerical_derivatives (originally lopale/
! lopale01/lopale06/lnumderiv) are NAMELIST values with different
! canonical spellings, kept local in core/read_input.f90 and copy-assigned
! (use_opal95_eos/use_opal2001_eos also re-synced after the "disable
! older OPAL EOS" validation block may override them). iopale is
! spelled identically to its canonical name -- use-associated directly.
      logical :: use_opal95_eos, use_opal2001_eos, use_opal2006_eos, &
           use_numerical_derivatives
      integer :: iopale

! former common/newopac/: laol_table_z1/laol_table_z2/opal_table_z1/
! opal_table_z2/opal95_single_table_z/alex_table_z1/kurucz_table_z1/
! kurucz_table_z2/molecular_opacity_logt_min/
! molecular_opacity_logt_max/use_alex06_tables/use_laol89_tables/
! use_opal92_tables/use_opal95_tables/use_kurucz90_tables/
! use_alex95_tables (originally zlaol1/zlaol2/zopal1/zopal2/zopal951/
! zalex1/zkur1/zkur2/tmolmin/tmolmax/lalex06/llaol89/lopal92/lopal95/
! lkur90/lalex95) are NAMELIST values with different canonical
! spellings, kept local in core/read_input.f90 and copy-assigned
! (use_alex95_tables/use_kurucz90_tables also re-synced after the
! "disable older Alexander opacities" validation block may override
! them). use_two_z_tables (originally l2z) is unused in
! core/read_input.f90, so still declared here for its other users.
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, &
           molecular_opacity_logt_min, molecular_opacity_logt_max
      logical :: use_alex06_tables, use_laol89_tables, use_opal92_tables, use_opal95_tables, use_kurucz90_tables, use_alex95_tables

! former common/miscopac/: ikur2/icondopacp are spelled identically to
! their canonical names -- use-associated directly.
! use_conductive_opacity (originally lcondopacp) is a NAMELIST value
! with a different canonical spelling, kept local in core/read_input.f90
! and copy-assigned.
      integer :: ikur2, icondopacp
      logical :: use_conductive_opacity

! former common/alexo/: alex95_table_unit (originally ialxo) is not a
! namelist value -- genuinely used in core/read_input.f90, renamed in
! place there.
      integer :: alex95_table_unit

! former common/alex06/: alex06_table_unit (originally ialex06) is not
! a namelist value -- genuinely used in core/read_input.f90, renamed in
! place there.
      integer :: alex06_table_unit

! former common/alexmix/: alex_mixture_x/alex_mixture_z (originally
! xalex/zalex) are not namelist values -- genuinely used in
! core/read_input.f90, renamed in place there. Their DATA defaults moved
! here since DATA can no longer target use-associated entities.
      double precision :: alex_mixture_x = 0.7d0, alex_mixture_z = 0.02d0

! former common/varfc/: lvfc is spelled identically to its canonical
! name -- use-associated directly; its DATA default moved here from
! core/read_input.f90 since DATA can no longer target a use-associated
! entity. use_diffusion_advection_transport (originally ldifad) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned. vfc, the block's per-zone WORK
! ARRAY (written by getfc/seculr every step), was never a control --
! evicted to flat star%vfc (2026 phase C batch 2).
      logical :: lvfc = .false.
      logical :: use_diffusion_advection_transport

! former common/notran/: no_am_transport_in_core (originally lnoj) is
! a NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned.
      logical :: no_am_transport_in_core

! former common/scv2/: scv_h_unit/scv_he_unit/scv_z_unit (originally
! iscvh/iscvhe/iscvz) are not namelist values -- genuinely used in
! core/read_input.f90, renamed in place there.
      integer :: scv_h_unit, scv_he_unit, scv_z_unit

! former common/alatm03/: allard_target_feh/allard_target_alpha/
! allard_use_tau100 (originally alatm_feh/alatm_alpha/laltptau100) are
! NAMELIST values with different canonical spellings, kept local in
! core/read_input.f90 and copy-assigned. allard_table_unit (originally
! ioatma) is not a namelist value -- genuinely used in
! core/read_input.f90, renamed in place there.
      double precision :: allard_target_feh, allard_target_alpha
      logical :: allard_use_tau100
      integer :: allard_table_unit

! former common/sbrot/: force_solid_body_rotation/solid_body_mode_flag
! (originally lsolid/impjmod) are NAMELIST values with different
! canonical spellings, kept local in core/read_input.f90 and copy-assigned
! (solid_body_mode_flag also re-synced after the "LSOLID overwrites
! IMPJMOD" line may override it).
      logical :: force_solid_body_rotation
      integer :: solid_body_mode_flag

! former common/cmixing/: cstmixing/cstdiffmix are NAMELIST values
! spelled identically to their canonical names -- use-associated
! directly.
      double precision :: cstmixing = 1.0d0, cstdiffmix = 1.0d0

! former common/acdpth/: output_ages_gyr/calcad_ageout_output_active/
! calcad_file_unit/ageout_model_output_flag/ageout_bracket_armed/
! acoustic_depth_output (originally ageout/lclcd/iclcd/ljlast/ljwrt/
! lacout) are not namelist values -- genuinely used in
! core/read_input.f90, renamed in place there. Their DATA defaults moved
! here since DATA can no longer target use-associated entities.
! acoustic_depth_cz_fraction/eos_adiabatic_gradient/acoustic_depth_heii/
! acoustic_crossing_time_seconds/acoustic_depth_cz_seconds/heii_zone_acoustic_width/acatmr_placeholder/
! acatmd_placeholder/acatmp_placeholder/acatmt_placeholder/
! atmosphere_sound_travel_time/iacat_placeholder/ageout_model_unit/
! laoly_placeholder/ijvs_placeholder/ijent_placeholder/ijdel_placeholder
! (former common/acdpth/'s remaining members, originally tauczn/
! deladj/tauhe/tnorm/tcz/whe/acatmr/acatmd/acatmp/acatmt/tatmos/iacat/
! ijlast/laoly/ijvs/ijent/ijdel) are unused in core/read_input.f90 --
! genuinely used (or carried as placeholders) in atm/turnover/acoustic_depths.f90,
! io/write_legacy_output.f90, core/main.f90, atm/atm_lib.f90, kap/opal95/getopal95.f90, so
! still declared here, using the majority spelling across those five
! files (atm/turnover/acoustic_depths.f90 and io/write_legacy_output.f90 give several of these real
! semantic names instead -- e.g. normalized_acoustic_depth,
! adiabatic_gradient -- but the 3-file placeholder spelling wins
! per the majority-name convention). compute_acoustic_depth (LADON) is
! spelled identically to its canonical name everywhere -- used, not a
! placeholder.
      double precision :: acatmr_placeholder(5000), acatmd_placeholder(5000), acatmp_placeholder(5000), acatmt_placeholder(5000)
      integer :: iacat_placeholder, ageout_model_unit
      logical :: laoly_placeholder
      integer :: ijvs_placeholder, ijent_placeholder, ijdel_placeholder
      double precision :: output_ages_gyr(5) = &
           [0.5d0, 1.0d0, 5.0d0, 10.0d0, 20.0d0]
      logical :: calcad_ageout_output_active = .false., ageout_model_output_flag = .false., &
           ageout_bracket_armed = .false., acoustic_depth_output = .false.
      integer :: calcad_file_unit

! former common/govs/: use_envelope_triangle_dt (originally ltrist) is
! a NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned.
      logical :: use_envelope_triangle_dt

! former common/pmmwind/: pmm_exponent_a/pmm_exponent_b/
! pmm_exponent_c/pmm_exponent_d/pmm_exponent_m/pmm_norm_jdot/
! pmm_norm_mdot/pmm_solar_pressure/pmm_solar_omega/
! pmm_solar_turnover_timescale/use_pmm_wind_law/scale_by_rossby_number/
! scale_by_b_field/wind_law_name (originally pmma/pmmb/pmmc/pmmd/pmmm/
! pmmjd/pmmmd/pmmsolp/pmmsolw/pmmsoltau/lmwind/lrossby/lbscale/awind)
! are NAMELIST values with different canonical spellings, kept local
! in core/read_input.f90 and copy-assigned (pmm_exponent_a/b/c/d/m and
! scale_by_b_field also re-synced after the K97/V13 PMM-windlaw branch
! may override them).
      double precision :: pmm_exponent_a, pmm_exponent_b, pmm_exponent_c, &
           pmm_exponent_d, pmm_exponent_m, pmm_norm_jdot, pmm_norm_mdot, &
           pmm_solar_pressure, pmm_solar_omega, pmm_solar_turnover_timescale
      logical :: use_pmm_wind_law, scale_by_rossby_number, scale_by_b_field
      character(len=3) :: wind_law_name

! former common/cwind/: wind_saturation_omega (originally wmax) is a
! NAMELIST value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned (also re-synced after the K97/V13
! PMM-windlaw branch may override it). exmd/extau/exr/exm/exl/expr/
! constfactor/excen/c_2/ljdot0 are spelled identically to their
! canonical names -- use-associated directly; ljdot0's DATA default
! moved here since DATA can no longer target a use-associated entity.
! wind_law_omega_exponent (originally exw) is not a namelist value --
! genuinely used in core/read_input.f90, renamed in place there.
! structfactor (former common/cwind/'s remaining member) is unused in
! core/read_input.f90, so still declared here for its other users.
      double precision :: exmd, extau, exr, exm, exl, expr, constfactor, &
           structfactor, excen, c_2, wind_law_omega_exponent
      double precision :: wind_saturation_omega
      logical :: ljdot0 = .true.

! former common/mag/: constant_background_diffusion_coeff/
! use_constant_background_diffusion (originally codm/lcodm) are
! NAMELIST values with different canonical spellings, kept local in
! core/read_input.f90 and copy-assigned.
      double precision :: constant_background_diffusion_coeff
      logical :: use_constant_background_diffusion

! former common/newparam/: all 30 NAMELIST values (the "intuitively
! named" timestep/tolerance parameters replacing DTDIF/DTWIND/HPTTOL/
! ATIME) are spelled identically to their canonical names everywhere --
! use-associated directly. Their DATA defaults moved here from
! core/read_input.f90 since DATA can no longer target use-associated
! entities.
      double precision :: flag_dx = 0.05d0, flag_dw = 0.10d0, flag_dz = 0.05d0
      double precision :: time_core_min = 1.0d-3, time_dl = 2.0d-2, &
           time_dp = 4.0d-2, time_dr = 2.0d-2, time_dt = 2.0d-2, &
           time_dw_global = 8.0d-2, time_dw_mix = 8.0d-2, &
           time_dx_core_frac = 0.5d0, time_dx_core_tot = 2.0d-2, &
           time_dx_shell = 0.1d0, time_dx_total = 1.5d-3, &
           time_dy_core_frac = 0.5d0, time_dy_core_tot = 2.0d-2, &
           time_dy_shell = 0.1d0, time_dy_total = 1.5d-3, &
           time_max_dt_frac = 1.5d0
      double precision :: tol_czbase_fine_width = 0.0d0, tol_dl_max = 0.02d0, &
           tol_dm_max = 0.08d0, tol_dm_min = 1.0d-8, &
           tol_dp_core_max = 0.05d0, tol_dp_czbase_max = 0.05d0, &
           tol_dp_env_max = 0.05d0, tol_dx_max = 1.0d0, tol_dz_max = 1.0d0
      logical :: lstruct_time = .false., lnewvars = .false.

! former common/newcross/: user-supplied nuclear reaction S-factors
! (and first/second derivative ratios relative to the Adelberger et
! al. 1998 Solar Fusion I values). s0_pp/s0_p_c12/s0_p_c13/s0_p_n14/
! s0_p_o16/s0_be7_electron/s0_be7_p/s0_n15_p_c12_branch/
! s0_n15_p_o16_branch/s0p_pp/s0p_he3he3/s0p_he3he4/s0p_p_c12/
! s0p_p_c13/s0p_p_n14/s0p_p_o16/s0pp_p_c12/s0pp_p_c13/s0pp_p_o16/
! s0p_be7_p/s0pp_be7_p (originally s0_1_1/s0_1_12/s0_1_13/s0_1_14/
! s0_1_16/s0_1_be7e/s0_1_be7p/s0_1_15_c12alp/s0_1_15_o16/s0p_1_1/
! s0p_3_3/s0p_3_4/s0p_1_12/s0p_1_13/s0p_1_14/s0p_1_16/s0pp_1_12/
! s0pp_1_13/s0pp_1_16/s0p_1_be7p/s0pp_1_be7p) are NAMELIST values with
! different canonical spellings, kept local in core/read_input.f90 and
! copy-assigned. s0_he3he3 (originally s0_3_3, renamed) is likewise
! kept local and copy-assigned. s0_he3he4 (originally s0_3_4, renamed)
! is likewise kept local and copy-assigned. s0_pep/s0_hep are NAMELIST
! values spelled identically to their canonical names --
! use-associated directly; their DATA defaults moved here since DATA
! can no longer target use-associated entities.
      double precision :: s0_pp, s0_he3he3, s0_he3he4, s0_p_c12, s0_p_c13, &
           s0_p_n14, s0_p_o16, s0_be7_electron, s0_be7_p, &
           s0_n15_p_c12_branch, s0_n15_p_o16_branch
      double precision :: s0p_pp, s0p_he3he3, s0p_he3he4, s0p_p_c12, &
           s0p_p_c13, s0p_p_n14, s0p_p_o16, s0pp_p_c12, s0pp_p_c13, &
           s0pp_p_o16, s0p_be7_p, s0pp_be7_p
      double precision :: s0_pep = 3.5734d-6, s0_hep = 8.6d-20

! (former common/optab/). use_scv_eos (originally lscv) is a NAMELIST
! value with a different canonical spelling, kept local in
! core/read_input.f90 and copy-assigned.
! (the SCV table declarations that followed this comment moved to
! state/scv_eos_lib.f90 -- phase six, step 3)


end module controls_lib
