# Readability wave 1 — io domain (worktree /Applications/YREC-wt/io, branch rs/io)

Assignment: /Applications/YREC-wt/tools/wave1_io.md (batches R1, R2, R3 of
audit/readability-sweep-2026-09-03 for the io domain).

## Commits

| batch | hash | verification |
|-------|------|--------------|
| R1 | e91d7d6 | gate1 + 37-case pin selection byte-identical |
| R2 | 554509d | gate1 + 37-case pin selection byte-identical |
| R3 | 6e5f55a | gate1 + 37-case pin selection byte-identical, aux battery 22 passed |

All builds with `-O3 -ffp-contract=off -finit-local-zero -Wall`, USE_HDF5=1.
No file outside the assigned list was edited; deps.mk did not need
regeneration (no source file, module or `use` statement was added or
removed; gate1 reports `deps.mk up to date` after every batch).

## R1 — dead code and stale comments (e91d7d6)

- src/numerics/numerics_lib.f90: deleted `boole`, `splinj`, `splnr`,
  `lir1` (no caller anywhere under src/, only a mention in GUIDELINES.md;
  grep of the whole src tree recorded in the commit message). Fixed the
  five "OPTIONAL ierr" comments (ierr is not optional; the historical
  `stop` became an ierr return). `ctridi` comment now says the SAVE is
  kept as-is and untested for being load-bearing; `tridia` comment
  states it has no SAVE; `bsstep` comment explains the blanket SAVE
  (x_sav/h etc. persist between calls). Header list of contents updated.
- src/numerics/intpar_lib.f90: stale path comment fixed
  (core/read_input.f90 -> io/read_controls.f90).
- src/io/read_controls.f90: removed unused `n_allard_teff`,
  `n_allard_logg`, `max_diag_pts`, `parmin_ln10`; replaced the 60-line
  DATA-statement archaeology block by a pointer to controls_lib plus the
  still-useful mixture / VNEW / cross-section notes; removed the
  commented-out LBNIN code, the orphaned `DO I = 2,4` comment and the
  duplicate expand_value comment; moved `shell_cmd = 'mkdir -p '` and
  its comment next to its only use in derive_options_and_open_files;
  header paths fixed.
- src/io/controls_lib.f90: comment-only tidy (generated
  controls_sync_lib.f90 unchanged because gen_controls_state.py strips
  comments — verified by regenerating and diffing).
- src/io/equal_to_model.f90: removed the debug `if (j .eq. 0) print*`;
  the `k0` clamp is now the plain `if (k0 .lt. 1) k0 = 1` (SUMMARY #6).
- src/io/model_to_equal.f90: removed the debug print, the dead CENTER
  DERIVATIVE block and its locals dr1/dr2/fac1/fac2/delr and
  deriv_factors; the two `intrp2` calls no longer pass the dead
  derivative-factor argument.
- src/io/write_gyre_pulse.f90: signature reduced to the arguments it
  reads (`num_shells, mass_coordinate, log_density, log_luminosity,
  log_pressure, log_radius, log_temperature, omega, pulse_path`) with
  `use star_info_lib, only: star, json`; src/io/yrec_output.f90 call
  site updated and two unused `use star_info_lib, only: star` removed.
- src/io/output_columns.f90: `if (ierr /= 0) then; close(u); return`
  split into a normal if-block (behaviour unchanged).
- src/io/print_allard_tables.f90: unused nta/nga removed.
- src/io/write_fgong_pulse.f90: header comment now states the actual
  ivers (1300) and E26.18E3 format.
- src/io/write_gsm_pulse.f90: stub comment says what it does ("returns
  without writing").
- src/io/read_model2.f90, read_yrec7.f90: stale comments fixed.
- src/state/phys_const_lib.f90: stale comment fixed.
- src/state/star_info_lib.f90: placeholder comments removed and the
  "NOT here" list trimmed to things that are actually elsewhere.

## R2 — named indices / constants (554509d)

- src/io/profile_output.f90: `integer, parameter :: col_zone = 1,
  col_mixing_type = 9` (the two integer-valued profile columns); used
  in the two `prof_sel(i) == ...` tests in `write_profile`. R2 also
  used them in `profile_column_names` (`names(col_zone) = 'zone'`),
  which broke the aux test
  `test_mesa_output.py::test_default_columns_lists_in_sync` (it
  harvests `names(<digits>) = '...'` with a regex; the 37 pins do not
  run it). The R3 commit restores the literal `names(1)`/`names(9)`
  with a comment saying why; the parameters remain.
- src/io/read_yrec7.f90: `double precision, parameter ::
  omega_log10_zero_sentinel = 58.9D0` with a comment (value of
  -log10(omega) marking a non-rotating shell in legacy model files);
  `envelope_fit_coeffs(i+i+i-3+j)` -> `envelope_fit_coeffs(3*i-3+j)`
  (integer arithmetic, identical value).
- src/numerics/numerics_lib.f90: `locate` — `below_table_slack =
  0.99d0`, `above_table_slack = 1.01d0` named parameters in the two
  clamp tests; `osplin` — `eps` is now a named parameter (was a local
  assigned on entry; `meval`'s `eps_tol` is intent(in) and was left).

Left alone in R2 (already named, or would not read better): `slots(15)`
(already `n_slots`-style), the `kttau` ladder in read_controls, the
`const1..const4` locals in `ysplin`, the intpt column indices.

## R3 — small shared helpers, token-identical (6e5f55a)

- src/numerics/numerics_lib.f90:
  - `kspline` is now a thin wrapper `call cspline(x, y, nm, natural_bc,
    natural_bc, y2)` (nm=4, natural_bc=1.0d30); its former body was
    token-identical to cspline's natural branch (yp1/ypn > 0.99d30 sets
    y2(1)=u(1)=0, qn=un=0).
  - `splinc` likewise wraps `cspline` (json-shaped dummies kept for the
    callers).
  - `splintd2` wraps `splint` (bodies identical; `6d0` vs `6.0d0` is the
    same constant).
  - new `pure function lagrange4(w, y)` = `w(1)*y(1)+w(2)*y(2)+w(3)*y(3)
    +w(4)*y(4)` (left-to-right sum, same tokens as the former inline
    expressions).
- src/io/model_to_equal.f90: the 16 inline 4-point Lagrange sums
  `interp_factors(1)*Y(k0)+...+interp_factors(4)*Y(k0+3)` replaced by
  `lagrange4(interp_factors, Y(k0:k0+3))` (also `composition(1,k0:k0+3)`
  and `composition(i_metals,k0:k0+3)`).
- src/io/equal_to_model.f90: the 2 sums for `delta_x` / `delta_z`
  replaced the same way.
- src/io/write_fgong_pulse.f90: `exp(ln10*cgl)` computed once into a
  local `grav_const_cgs` (used for `glob(15)` and for the per-shell
  `grav`); the product `G*m/(r*r)` keeps its operand order.
- src/io/write_gyre_pulse.f90: same hoist for the per-shell `grav`.
- src/io/profile_output.f90: literal `names(1)`/`names(9)` restored
  (see R2 above); comment added next to the parameters.

interp vs intrp2: weights are token-identical, but they are already the
same module and the remappers call the one they need; no merge done.

## Deferred (cross-domain)

- src/rotation/microdiff/microdiff_mte.f90 (lines ~157-167) and
  src/rotation/microdiff/microdiff_etm.f90 contain the same inline
  4-point Lagrange sums with weights `facinterp(1..4)`; after R3 they can
  become `lagrange4(facinterp, Y(k0:k0+3))` (numerics_lib is already
  used there). Owner: rotation domain; not touched here.
- None of my edits required changes outside the assigned files.

## Reverted (changed numbers)

- None. Every batch was byte-identical on the first pin run. The one
  revert (profile_output `names(col_zone)` -> `names(1)`) was for the
  aux test's source regex, not for numbers; the aux battery was
  1 failed / 21 passed with the R2 form and 22 passed after the
  revert.

## Skipped (disagree with reviewer / too risky)

- numerics_lib `tridia` no-op rescale: class B (changes numbers).
- numerics_lib `bsstep` underflow test: known issue, class B.
- Single-precision literal defaults in controls_lib: class B.
- Dead namelist members (pmmwmax, zalex2, zopal952, lpmm, lsolwind,
  sstandard): removing them changes accepted input (namelist would
  reject old inlists) — not a readability change.
- `'tridia:'` message text and the FGONG `X` descriptor: output/terminal
  text reaching files; left verbatim.
- `private` default in controls_lib: would touch the generated
  controls_sync_lib.f90 interface and every consumer; cross-domain.
- `ysplin` const1..const4 rewrite: would change parenthesization.
- `end=9999` -> `iostat`: control-flow change in an input path; not R1-R3.
- 'GSM' vs 'gsm' case-fold: changes accepted input.
- Fixed scratch sizes (#9), simeqc `write(5,..)`: not readability, and
  the write unit is output behaviour.
- read_controls loop-index idiom / `kttau` ladder: rewriting would not
  be token-identical.
- print_allard_tables side-effect move: changes call order relative to
  file opens; left for a class-B batch.

## Things the audit missed

- The aux test `test_mesa_output.py::test_default_columns_lists_in_sync`
  parses `history_output.f90` / `profile_output.f90` source text with
  `names\((\d+)\)`; any readability edit that replaces those literal
  indices by named constants fails it. Worth a note in GUIDELINES.md
  (not my file) or a more tolerant regex in the test (not my file).

- read_controls.f90: the `mkdir -p` command string was built in one
  internal procedure and used in another, with the explanatory comment
  dangling next to the wrong one (fixed in R1).
- read_controls.f90: `last_slash_idx` (a read_controls local, host-
  associated into derive_options_and_open_files, lines ~1805-1813) is
  only assigned when FSHORT contains a `/`; with a bare file name it is
  read uninitialised and the code relies on `-finit-local-zero` giving
  0 (and then `fshort(1:0)` is an empty directory). Left as is per the
  rules (adding an explicit init would be a semantic decision); flagged
  for the owner.
