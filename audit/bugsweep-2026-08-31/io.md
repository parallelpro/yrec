# Bug sweep: /Applications/YREC/src/io/*.f90

Files reviewed in full: check_controls.f90, controls_lib.f90, equal_to_model.f90,
history_output.f90, luout_lib.f90, model_to_equal.f90, output_columns.f90,
print_allard_tables.f90, profile_output.f90, read_controls.f90, read_mod_model.f90,
read_model2.f90, read_yrec7.f90, run_log.f90, write_fgong_pulse.f90,
write_gsm_pulse.f90, write_gyre_pulse.f90, write_mod_model.f90,
write_output_headers.f90, yrec_output.f90.

Special checks requested and performed:
- read_controls.f90 legacy-name adoption copies: mechanically extracted all 224
  copy assignments in adopt_canonical_names (plus the re-syncs in
  derive_options_and_open_files and echo_settings) and cross-checked against the
  controls_lib comments, the namelist member lists, and the original parmin.f
  COMMON orders (DIFUS/ITDIF1-2 verified against consumers, s0_* table, vmult,
  gravs*, intatm/intenv, calstar/cals2, newopac, pmmwind). No transposed or
  duplicated pairs found (no duplicate LHS; only kttau appears twice as RHS,
  intentionally: atm_choice and star%atm_choice_initial). Every namelist-spelled
  local is either copied or use-associated; the only never-consumed locals
  (fstch, zalex2, zopal952) are inherited dead namelist slots. ONE completeness
  defect found: the lfirst(1) force (finding 2).
- write_fgong_pulse var()/glob() slots checked against the FGONG standard
  (JCD ivar=40 layout) and the species_slot table in core/stitched_model.f90:
  var(1-15,17,18,19,21-25,29-31,33-35) all land on the correct FGONG columns
  (21 He3, 22 C12, 23 C13, 24 N14, 25 O16, 29 H2, 30 He4, 31 Li7, 33 N15,
  34 O17, 35 O18); var(2) is correctly ln(m/M); var(15)=N^2 r/g; glob(13-15)
  age/Teff/G correct. No slot bugs found.
- write_mod_model / read_mod_model symmetry: same 15-entry slots() table, same
  columns header string, same per-row order (conv flag before omega), same
  globals; symmetric. read_yrec7 and read_model2 were line-checked against
  getyrec7.f / getmodel2.f at 6cd5673 (formats, CFENV indexing, the odd-pair
  extended-composition read, the omega 58.9 sentinel, the TLUMX max over
  |neu| — all faithful; i_lum_* index order verified = 1..7).

## model_to_equal.f90:125,233,255,316,329 -- N15 (slot 8) used as the metal mass fraction in metal diffusion
- class: physical
- severity: high
- confidence: high
- provenance: inherited (git show 6cd5673:src/model_to_equal.f: "C MASS FRACTION OF METALS / EZ_H(I)=FACINTERP(1)*HCOMP(8,K0)+...", same in the zone-center block and first-point block)
- detail: model_to_equal builds the equal-grid metal-abundance run (EZ /
  star%metal_abundance_change[_mid]) by interpolating composition(8,...) under
  the comment "MASS FRACTION OF METALS". The composition array is the standard
  15-slot HCOMP (caller src/rotation/microdiff/gravitational_settling.f90
  documents composition(1)=X, (2)=Y, (3)=Z; star_info_lib pins i_metals=3,
  i_n15=8), so slot 8 is the N15 mass fraction (~1e-8..1e-6), not Z (~0.02).
  The settling solver then evolves this "Z" profile and equal_to_model applies
  the change EZ-EZ_ORIG to composition(3,...) (real Z) — an asymmetry between
  the two mirrored transforms. Net effect: whenever use_diffusion_z is on, the
  computed metal-settling change is based on the wrong (vastly smaller and
  differently-shaped) abundance profile, so metal diffusion is effectively
  inert/incorrect. Present verbatim in the F77 original, so it is a long-
  standing physics defect, not a modernization slip — but it directly negates
  the feature the LDIFZ flag claims to provide.

## read_controls.f90:2317 -- forced lfirst(1)=.true. no longer reaches first_call_flag
- class: logical
- severity: medium
- confidence: high
- provenance: modernization (6cd5673:src/parmin.f:1282 "LFIRST(1) = .TRUE." acted on the single COMMON/CKIND/ LFIRST that all consumers shared)
- detail: interpret_kind_cards sets the local lfirst(1)=.true. as a safety
  against decks that (nonsensically) mark card 1 as a continuation run. But
  the canonical copy first_call_flag = lfirst happens earlier, at line 1613 in
  adopt_canonical_names, and is never re-synced after interpret_kind_cards;
  the final store_controls_to_star then persists the un-forced value into
  star%job%first_call_flag, which is what read_starting_model/run_yrec/the
  calibration setups actually consult. A deck with LFIRST(1)=.FALSE. (or
  LFIRST=50*.FALSE.) now skips loading the starting model on card 1 instead
  of being silently corrected as the F77 code guaranteed. The forced local is
  still used for interpret_kind_cards' own run-plan printout, so the log would
  even claim "from the starting model" while the run does otherwise. Fix: move
  the force before line 1613, or re-copy after interpret_kind_cards.

## output_columns.f90:97-98 -- append_column writes past sel(max_cols) for long columns files
- class: logical
- severity: medium
- confidence: high
- provenance: new-code (2026 output machinery, no F77 ancestor)
- detail: append_column does nsel = nsel + 1; sel(nsel) = j with no check
  against max_cols (=128). parse_columns feeds it one call per non-comment
  line of a user columns file; names may repeat, so a history/profile columns
  file with more than 128 selecting lines overflows the fixed-size sel arrays
  in history_output/profile_output (module state) and corrupts adjacent module
  variables. Purely config-triggered, but it is user input with no guard and
  no diagnostic.

## read_controls.f90:1824-1830 -- last_slash_idx used uninitialized when the log path has no '/'
- class: logical
- severity: low
- confidence: medium
- provenance: unclear (the slash-scan pattern predates the modernization but postdates 6cd5673; the current fshort-based form is 2026)
- detail: derive_options_and_open_files scans fshort backwards for '/' and
  only assigns last_slash_idx inside the if; if the expanded log path contains
  no slash (e.g. a deck sets FSHORT='run.log', overriding the
  '{YREC_OUTPUT}/run.log' default), last_slash_idx is read undefined at
  fshort(1:last_slash_idx) and in the mkdir command construction — undefined
  behavior, potentially a bogus mkdir or a crash. yrec_output's own
  output_init_mesa handles the islash=0 case explicitly (line 85-90), so the
  two path-splitters are inconsistent.

## read_controls.f90:1792-1807 -- expand_value applied asymmetrically: fopale/fopale01 (and others) never expanded
- class: logical
- severity: low
- confidence: high (asymmetry), medium (real-world impact)
- provenance: inherited (feature commit e5c35e5 already had exactly this list, modulo retired paths)
- detail: Placeholder expansion ({YREC_INPUT}/...) is applied to fopale06 and
  fcondopacp but not to their sibling EOS/opacity paths fopale and fopale01,
  which are opened directly at lines 1878-1884; likewise fdyn, fiso, fmonte1/2,
  fkur, fkur2, flaol, flaol2, flldat, fopal2, fmhd1-8 and opecalex are opened
  or stored without expansion. Parallel branches that should mirror each other
  don't: a deck using {YREC_INPUT} for LOPALE01's table fails with a
  file-not-found on the literal placeholder string, while the LOPALE06 spelling
  works. Inherited from the expansion feature's introduction, worth closing.

## read_controls.f90:2327-2333 (interpret_kind_cards) -- stop conditions never armed for kind-3 (rescale+evolve) cards
- class: logical
- severity: low
- confidence: medium
- provenance: inherited (6cd5673:src/parmin.f computes LENDAG/LSETDT only inside the KINDRN.EQ.1 branch, same as now)
- detail: end_age_stop_active/timestep_override_active are computed only for
  kind-1 evolve cards. A kind-3 card evolves too, but its slots keep whatever
  they held (module storage without initializer: zero-fill on first run, stale
  values on in-process re-reads), so ENDAGE/SETDT/END_*CEN on a rescale+evolve
  card are silently ignored — or, worse, a re-run in the same process can
  inherit the previous deck's flags for that slot. Inherited behavior, but the
  module-variable re-read path makes the stale-value case new-ish.

## write_gyre_pulse.f90:86-89 -- brunt_n2(2) read uninitialized when num_shells==2
- class: numerical
- severity: low
- confidence: high (mechanism), low (reachability)
- provenance: new-code (2026 pulsation output)
- detail: the centered-difference loop runs i=2..num_shells-1, which is empty
  for num_shells==2, yet the boundary fill (guarded only by num_shells>=2)
  copies brunt_n2(2) into brunt_n2(1) — an uninitialized read. Unreachable for
  any physical model (nz is always >> 3), so cosmetic in practice.

## write_gyre_pulse.f90:24-50 -- dummy argument log_luminosity actually carries linear L/Lsun
- class: logical (name-vs-code mismatch)
- severity: low
- confidence: high
- provenance: modernization (naming introduced in the 2026 rename; the numeric treatment matches the caller)
- detail: yrec_output.f90:185 passes star%luminosity_lsun (linear, because
  shell luminosity can be negative), and the code multiplies it linearly by
  solar_luminosity_cgs (lines 50, 94) — numerically correct — but the dummy is
  named log_luminosity alongside genuinely-log log_radius/log_pressure/etc.
  that ARE exponentiated. No wrong output today; a classic trap for the next
  editor of a file the sweep flags as a known duplicated-physics site.

## read_model2.f90:94 -- TLUMX skip changed 12X -> 10X (silent fix, contradicts "unchanged" header)
- class: logical (documentation/provenance mismatch)
- severity: low
- confidence: high
- provenance: modernization (6cd5673:src/getmodel2.f format 40 is "FORMAT(12X,1P7E17.9)"; putmodel2.f writes "('TLUMX',5X,1P7E17.9)", i.e. a 10-char prefix)
- detail: the modern reader skips 10 columns, which exactly matches what
  putmodel2 writes; the F77 reader's 12X misaligned every field by 2 columns
  (benign for positive 2-digit-exponent values under blank-null rules, but it
  dropped the sign of a negative TLUMX(6)). So the modern code is BETTER than
  the original — but the file header claims "Logic and numerics are unchanged
  from the original getmodel2.f", which is false for this record and undermines
  the byte-pinning audit trail. Since TLUMX is only a cosmetic seed
  (recomputed on the first energy call), impact is nil; the header should be
  corrected rather than the code.

## read_controls.f90:2407 -- Z rescale guard .ge.0 contradicts the "0 = keep current" contract
- class: logical
- severity: low
- confidence: medium
- provenance: inherited (6cd5673:src/parmin.f: "IF(RESCAL(3,NKIND).GE.0.0D0) ZENV=RESCAL(3,NKIND)"; format 452 documents "0=USE CURRENT VALUE")
- detail: rescale_params(3,*) defaults to 0 (and is never set at all for
  kind-1 cards, whose slots are zero-filled module storage), yet the .ge.0
  guard overwrites star%envelope_metal_fraction with 0 on every card and then
  (2026 addition) calls eos_set_mixture with Z=0. Downstream model loading
  recomputes the envelope composition, which is why this has never visibly
  fired, but the 2026 eos_set_mixture call extends the window in which the EOS
  mixture is configured for Z=0. The printed contract and the guard disagree;
  .gt.0 matches the documentation.

## Weak/uncertain observations
- model_to_equal.f90:63 `mod(total_radius_span,min_radius_spacing).ne.0.0d0` — the known mod(dp,dp) ceiling idiom (float equality on computed doubles), same family as the already-reported rezone finding; inherited.
- equal_to_model.f90:57,115,138 metal_scale_ratio = metal_new/composition(3,...) divides by Z with no zero guard — Z=0 deck with LDIFZ on gives NaN; inherited (ZZ2=ZZ/HCOMP(3,I)).
- equal_to_model.f90:84 `if (k0 .eq. 0) k0=1` fires only after min(k0,num_equal_points-3) made k0 0 (NPT=3); NPT=3 still reads equal_radius(4) beyond the meaningful range — inherited "JVS fix for NPT=3?" acknowledged as dubious in the original comment.
- model_to_equal.f90:87 first-point search loops interp_search_index=2..num_equal_points over the MODEL radius array (bound should be a model-grid index) — inherited verbatim (DO 15 IU=2,NPT over HRU), works only because DR<=DRMIN makes the first hit come early.
- read_controls.f90:1936-1940 the LSEMIC error prints lovste before lovstc under the label "OVERSHOOT - CORE,ENVELOPE,INTERMEDIATE" (envelope value shown under CORE) — inherited cosmetic transposition in an abort message.
- yrec_output.f90:149 star%m(star%h_shell_zone_begin-1) indexes zone 0 if a flagged H-shell ever begins at zone 1; guarded in practice by how has_h_shell is computed.
- yrec_output.f90:183 gyre_suffix is char(len=5) with '(I5.5)' — model numbers >99999 print '*****' into the filename; cosmetic.
- write_fgong_pulse.f90:6 header says "5 values per line in 1P5E16.9" but the code writes the ivers-1300 wide E26.18E3 layout (correctly); stale comment.
- write_fgong_pulse.f90:49-50 glob(4)/glob(5) use the outermost point's CURRENT Z and X as the FGONG "initial" Z/X0; acceptable proxy, mildly wrong after heavy diffusion.
- write_fgong_pulse.f90:67,117 bare `X` edit descriptor (no count) is nonstandard Fortran (extension); compiles under gfortran.
- write_gsm_pulse.f90:54-57 on h5fcreate failure the routine warns and returns, but later h5lt calls' herr values are never checked — a mid-write HDF5 failure produces a silently truncated GSM file.
- print_allard_tables.f90:21-24 a "print" routine mutates atm_table%allard_teffl_min/max as a side effect — inherited verbatim from alprint.f (which recomputed TEFFLmin/max the same way).
- read_controls.f90:2429 nml_file_has_group matches '&controls' as a prefix, so a group named &controls_extra would also trigger single-inlist mode; no such group exists today.
- read_controls.f90 fstch/zalex2/zopal952 are namelist-visible but dead (never stored or opened) — inherited dead slots kept for input-file compatibility.
- history_output.f90:201 surface_z_div_x divides by surface X with no guard — only an issue for hydrogen-free surfaces.
- controls_lib.f90:707 comment claims cstmixing/cstdiffmix are "spelled identically to their canonical names" but the canonical names are constant_mixing_coeff/constant_settling_reduction (the copy in read_controls:1544-1545 is correct); stale comment only.

## Summary
- high severity: 1 (model_to_equal N15-as-Z metal diffusion, inherited)
- medium severity: 2 (lfirst(1) force lost in modernization; append_column overflow)
- low severity: 6 findings + 15 weak/uncertain one-liners
