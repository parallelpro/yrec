# eos/ -- fourth evaluation (2026-10-31, yrec-next e8c847ba, read-only)

Scope: every file of src/eos/ read whole (README, eos_lib, eos_configuration, eos_results, eos_messages, eos_mixture, blended_eos, yale/yale_eos + yale_eos_tables, scv/scv_eos + scv_eos_tables, opal/opal_eos + opal/opal95_eos + opal_eos_tables, mhd/mhd_eos + mhd_eos_tables, test/test_eos), against src/rotation/README.md, GUIDELINES "Naming", "The star object: six rules", "Declarations and interfaces", "Code style", and TARGET_ARCHITECTURE 4.5(c). Since 6e058523 the only change under src/eos/ is one README line (`star%ctrl` -> `star%physics`), so the code is the one the third evaluation scored. Scores below are my own calibration. Where one is half a point under the third evaluation's, that is because of the new findings and a stricter reading of the F77-shaped kernels. Nothing regressed.

Checked and found true: the facade's three entries (check_boundaries `FACADE`: eos_get, eos_init, eos_set_mixture); STAR_FREE_DIRS and the seven STAR_BLIND_FILES; the settings and paths builders and the install site (setup/microphysics_settings.f90, settings_install.f90, run_setup.f90); every star-layer eos_get call passes `composition_at_zone` (20 of 20); the defaults (physics.defaults: cutoff 6, eta -1/1, switches off); constants are physical_constants parameters; the step sizes 3e-4/1.5e-3; the ramp weight; the OPAL 2006 T floor (log T 3.27); 49 ok/FAIL lines in the pin; the recipes' file and routine names (io/README §5 "Add a table path", expand_paths, derive_controls' newest-wins chain).

Scratch evidence: SP/review4_eos/drv.f90 linked against the tree's existing eos/numerics/math objects (nothing in the repo built or touched), output in SP/review4_eos/drv.out.

## Scores (1-10)

| Criterion | 1st -> 2nd -> 3rd -> now | Justification |
|---|---|---|
| Readability | 6 -> 7 -> 8 -> **7.5** | The facade, blend_eos's four named steps, the mixture and the settings read at rotation's bar. But about 4,500 of the ~5,700 source lines still have the F77 shape. These are the Fermi-Dirac solver (a 40-member `fermi_solution`, unnamed table rows `ff(1..5)`, `pr1..pr9` reused for unrelated quantities), the Saha derivative step (`sqet`, `sx1qt`, ...), OPAL's 350-line `interpolate_opal_eos` and the `esac` transliteration, and MHD's 33-argument `mhdst1`. History notes and commented-out code also remain in procedure bodies. |
| Organization | 6 -> 7 -> 8 -> **8** | A clean facade, one settings record, a module per family with separate table modules, and the domain is star-free and star-blind. SCV is nested inside `gas_eos` and MHD is a branch of `eos_get` rather than sources. The test-only entry points (`gamma1_at_density`, `esac01/esac06`, `set_eos_settings`) are public in production modules. Writable public state remains (`eos_mix`, the SCV and OPAL search state), and OPAL 1995 is still a separate family (M94). |
| Maintainability | 4 -> 6 -> 7 -> **7** | There are 49 mutation-confirmed checks and a byte pin. Against that, nine latent defects are open (M92-M94, M105-M107; M119 and M87 are weak tests), the kernels write their own messages (GUIDELINES code style 7), OPAL coverage is 71%, and with diffusion of Z the zone's 1/mu_ion follows two different rules depending on the path (new finding 2). |
| Extensibility | 4 -> 6 -> 7 -> **6.5** | Recipes exist for a control, a source, an OPAL vintage and a blended region, and the OPAL vintages are rows of `opal_layout`. But a source still means edits in eight places, with no source descriptor (TARGET_ARCHITECTURE 4.5(c) is not started), and the source recipe as written gives NaN or zero mean weights (new finding 1) and drops ion fractions (M123). |
| Docs accuracy | 6 -> 7 -> 8 -> **7.5** | The README is thorough and right on the main paths. Against it: the recipe trap above, deck counts off by six, five stale comments that name the wrong module, file, size or unit, gaps where a reader would assume otherwise (which Z OPAL uses, which mixture rule each path uses, where the second derivatives and eta come from), and the open M119, M122 and M123. |

## Summary

eos/ is unchanged since the third evaluation and remains the second-best-documented domain after rotation. Its front half (eos_lib, blended_eos, eos_mixture, the settings and results records) meets the rotation standard. The family kernels (Yale fully ionized, SCV, OPAL, MHD) keep their F77 structure under modern names, and so do most of the open defects. The review found ten new items, none severe:
- **Recipe trap:** the source recipe tells a new source to compute its mean weights with `mean_weights_from_density` on a fresh record. That routine reads T, P and 1/mu_ion, which a fresh record holds at zero, so a newcomer gets NaN R/mu and a zero 1/mu_ion, which `adopt_table_result` then copies into the zone.
- **Two mixture rules:** with diffusion of Z, the OPAL and MHD paths count a zone's excess metals as carbon, while the Yale and SCV paths scale the mixture. The scratch run shows a 3.9e-5 step in 1/mu_ion at the OPAL edge; it is inherited from the F77 MU and affects only the 1/mu_e output.
- **SCV H+ fraction:** the denominator counts H2 once, overstating the H+ fraction by up to 33% where H2 and H+ coexist (F77).
- **Smaller items:** reader checks, the cutoff's range, deck counts, stale comments, doc gaps and style residue.

To reach 9 or more, eos needs the table-source descriptor, the latent defects fixed (most fixes are byte-identical), and the Fermi-Dirac, OPAL and MHD kernels split into named steps with named table rows.

## New findings (not in bug_sweep), ranked

1. **src/eos/README.md:122 and :128** (the source recipe, steps 1 and 4), with eos_mixture.f90:145-150 and blended_eos.f90:461-476. *Issue:* step 1 says a source fills "the mean weights (`mean_weights_from_density` gives R/mu and 1/mu_e from the density)", and step 4 says it is called "on a fresh record". `mean_weights_from_density` computes beta P/(rho T) and keeps r's 1/mu_ion. A fresh `eos_result` has temperature, pressure and ion_mean_weight_inverse all 0 (eos_results.f90:28, :37), so R/mu = 0/0 = NaN. `adopt_table_result` then copies the source's 1/mu_ion (0) and 1/mu_e into the zone's result. The recipe's field list names neither T, P nor 1/mu_ion. OPAL avoids the trap by setting T and P itself and calling `compute_mean_weights` from X and Z (opal_eos.f90:82-83, :157-160); the README's "traced" recipe evidently did the same without saying so. *Fix:* the recipe says the lookup sets T and P and takes 1/mu_ion from `compute_mean_weights` (or is handed a copy of r), and says which fields `adopt_table_result`/`blend_table_result` read and write (together with M123). *Class:* A. *Size:* S.

2. **eos_mixture.f90:112-136** (`compute_mean_weights`), called by opal_eos.f90:158-160 and mhd_eos.f90:105-107, against blended_eos.f90:329-360 (`zone_mixture`'s diffusion branch). *Issue:* `compute_mean_weights` ignores `use_diffusion_z`. With diffusion of Z (40 of the 310 testsuite and examples decks), a zone's 1/mu_ion on the OPAL and MHD paths counts its excess metals as carbon (1/12.015), while on the Yale and SCV paths all the metals scale with Z/Z_env. *Measured* (scratch driver: production with diffusion of Z, X = 0.70, Z = 1.1 Z_env, log T 5): 1/mu_ion is 0.7659619 at log P 12.80, where OPAL is adopted, and 0.7659232 at log P 13.05, on the Yale/SCV result. The step is 3.9e-5, which is dZ times (0.0832 minus the mixture's mean 1/A of 0.0598). R/mu is unaffected (beta P/(rho T)), and 1/mu_e, which is backed out of it, absorbs the step. Only 1/mu_ion and 1/mu_e move, and they feed only output (the profile column `mu_e_inv`, the pulse files); semiconvection reads their sum. The README (§2.1, steps 1 and 4) implies a single mixture rule. test_eos's OPAL edge-ramp check (`results_agree`, test_eos.f90:1192-1220) does not compare 1/mu_ion, and the test runs without diffusion of Z, so this is unseen. *In the F77?* Yes: `MU` always uses excess-as-carbon, and EQSTAT2 alone branches on LDIFZ. *Fix:* compute_opal_eos takes the 1/mu_ion that `zone_mixture` already put in r (one rule for every path), and `results_agree` compares 1/mu_ion, with one ramp run under diffusion of Z. Alternatively, document the two rules. *Class:* B (output only). *Size:* S.

3. **scv_eos.f90:643-651** (`scv_ionization_fractions`). *Issue:* `ion_fraction(1) = x_H+/(x_H2 + x_H + x_H+)`, but the H nuclei number 2x_H2 + x_H + x_H+. The routine's own comment calls the fractions "of all particles" (:624-625) and "of the H nuclei" (:647-648), and the envelope build counts H2 twice (n = 2(1-Y)rho/m_H/(1 + 3x_H2 + x_H), :715-717). The right value is (1 - x_H2 - x_H)/(1 + 3x_H2 + x_H); the code gives (1 - x_H2 - x_H)/(1 + x_H2 + x_H). Over the shipped h_tab_i.dat nodes with an H+ fraction above 1e-4 and log P_gas <= 10, the overstatement reaches 33% (log T 4.02, log P 10: 3.8e-4 against 2.9e-4), and it approaches a factor 2 at log P 12. OPAL's adopt step leaves `ion_fraction` alone (M123), so the SCV value is what production decks report wherever SCV covers a zone below log T 6. It reaches kap's conductive opacity weight (conductive_opacity.f90:138) and the envelope and atmosphere ionization outputs. *In the F77?* Yes (eqscve.f:444). *Fix:* 2 xtf_h2 in the denominator, after confirming that the tables' X_H2 is a per-particle fraction, as the comment and the build assume. *Class:* B (negligible in structure). *Size:* S.

4. **scv_eos.f90:909-911 and opal95_eos.f90:1007-1014.** *Issue:* `read_scv_rows` reads each row's log T and point count from the He and Z tables and discards them without comparing them to the H table's, so a He or Z table on another grid is read misaligned without a word. The OPAL 1995 reader exits its loop silently on a zero T6, which leaves `t6_grid` entries at 0 and 1/(0-0) spacings; the 2001/2006 reader rejects the same table with an error (opal_eos.f90:881-888). README §3 says "every reader checks its reads". (A 2026-08-31 audit note records the shared-grid assumption, but bug_sweep has no entry.) *Fix:* compare the copies and reject a mismatch, and make the 1995 zero T6 an error, both with the reader's message. *Class:* A. *Size:* S.

5. **yale_eos.f90:733-737, eos_configuration.f90:36-39 and eos_lib.f90:102-118.** *Issue:* `saha_log10t_cutoff` is never range-checked. The Fermi-Dirac table's T rows start at log T 5.2: `jt1 = int(5(log T - 5))` is clamped to 1, so the interpolation offset `ttud` goes negative below 5.2. The Saha ramp reaches cutoff - 0.5, so any cutoff below 5.7 makes the fully ionized solver extrapolate the table quadratically in log T, silently. In the scratch run with cutoff 5.0 (log T 5.0, log P 9) the density stays within 2.5e-5 of ideal, so at low density the extrapolation is benign; at degenerate points it is untested. No shipped deck sets the control, and README §4 gives no range. *Fix:* `eos_init` refuses (or warns on) cutoff - 0.5 < 5.2, with the table's T origin as a named parameter, and §4 states the range. *Class:* A. *Size:* S.

6. **src/eos/README.md:45 and :111.** *Issue:* "298 of the 304 shipped decks (testsuite and examples)" is wrong. testsuite/ and examples/ hold 310 decks (git ls-files): 304 on OPAL 2006 + SCV + Debye-Hueckel, 2 on SCV, 4 on Yale. The 298/304 count leaves out examples/class_exercises' six decks, all of them production configuration. *Fix:* "304 of the 310". *Class:* A. *Size:* S.

7. **Stale comments.** *Issue:*
   - yale_eos.f90:184: "Called from blend_eos.f90". No such file exists; the caller is `gas_eos` in blended_eos.f90.
   - yale_eos_tables.f90:32-33: "eos_lib's debye_huckel_composition_of". It is in eos_mixture.
   - opal_eos_tables.f90:344-346: "those (48 MB apiece)". The object's symbol table gives 15.5 MB per 2001/2006 instance and 2.0 MB for 1995.
   - mhd_eos.f90:45-46, :126-127, :159-160: errors "with a message in the run log". These are `eos_report` messages, which go to the terminal and the run log.
   - mhd_eos.f90:510: mhdpx's input "PGL = log10 (gas pressure)". The call passes total log P (:70, as the F77 MEQOS did), and slot 2 is documented as "log10 total pressure" (mhd_eos_tables.f90:52). One of the two is wrong, and if the header is the right one, the MHD density ignores radiation pressure.

   *Fix:* correct each. For the MHD pressure, settle it against the table documentation or by the beta = 10^(slot20 - slot2) identity. *Class:* A. *Size:* S.

8. **README §2.1 step 4 and §3 (:43, :63), and eos_results.f90:18-19, :39.** *Issue:* the documentation omits or misstates things a reader would assume otherwise:
   - OPAL is evaluated at its table's single Z with the zone's X; the zone's Z enters only the mean weights (`interpolate_opal_eos` takes no Z; opal_eos.f90:114-115, :691). This matters under diffusion of Z and Z rescales.
   - Inside OPAL "the second derivatives are the Yale (Saha) solution's" (README :63), but above the cutoff, which is most of OPAL's domain, they are the fully ionized gas's.
   - eos_results' header says the second derivatives are set only with `want_derivatives`, which the MHD path does not honour (the README says so; the record's header does not).
   - eta is not blended across the Saha ramp (it jumps to 0 at the ramp's foot). Under OPAL it is the Yale solve's value at the Yale density, while eos_results says "fully ionized gas; 0 otherwise".

   *Fix:* a short "which source fills which field" table in §3, plus one line each in §2.1 and in eos_results' header. *Class:* A. *Size:* S.

9. **test_eos.f90:1028-1033 against README.md:107.** *Issue:* the check "the record's defaults are the controls'" compares `eos_settings_t()` with literals (6, -1, 1, none, false), not with physics.defaults or the registry. Changing a default in physics.defaults would leave eos_settings_t stale and the check would still pass, so README §4's "Its defaults are the controls'" is kept true by hand. *Fix:* generate the defaults, or compare them with physics.defaults in test_controls or check_settings_fields. *Class:* A. *Size:* S.

10. **Style residue against GUIDELINES "Code style" rule 5.** *Issue:*
    - Commented-out code and open questions sit in radsub (opal95_eos.f90:826-836: "cvt=(eos(5)*molenak/...", "!      EN=EN*REVISE", "should this be ET=ET*REVISE???").
    - History sits in procedure bodies rather than Sources sections: scv_eos.f90:127-130 and :152 ("2026 (bugsweep Batch 2)..."), opal_eos.f90:245-247 and :758-761, opal95_eos.f90:630-633 and :739-743, test_eos.f90:274-284 ("Error paths (2026, docs/ROADMAP.md stage 3)...").

    *Fix:* move the history to each file's Sources section; turn the open questions into a bug_sweep "Decisions needed" entry or drop them. *Class:* A. *Size:* S.

## What would raise each score (concrete)

**Readability 7.5 -> 9+**
- Name the Fermi-Dirac table's five rows (`i_fd_*`). Split `fermi_solution` into the interpolation stencil, the pressure terms and the derivatives. Replace the reused `pr1..pr9`/`cden*` temporaries in `iterate_fully_ionized_density` and `compute_internal_energy_and_cp` with named quantities.
- Give the Saha sensitivities names that read (`de_dlnt`, `dx_he1_dlnt`, ...), or document them field by field in `saha_sensitivities`.
- Cut `interpolate_opal_eos` into named steps (the X stencil, coverage and orders, the X interpolation, T6-rho, scaling and radiation), as `blend_eos` is. Write `compute_scv_envelope_eos`'s four-way smoothing as one helper applied in T and then in P.
- Give `mhdst1` a per-table record instead of 33 arguments, with assumed-shape arrays (Declarations rule 4).
- Move the in-body history and commented-out code (finding 10) to Sources sections.

**Organization 8 -> 9+**
- Make SCV a step of `blend_eos` with its own edge weight; this also fixes M88. Make MHD a source.
- Introduce the table-source descriptor (name, window, reader, lookup, edge weight) that TARGET_ARCHITECTURE 4.5(c) asks for and ROADMAP defers "to a third source". With it, `blend_eos` loops over an ordered list.
- Move the test-only entry points (`gamma1_at_density`, `esac01/esac06`, `set_eos_settings`) into a test-support module or behind a test interface.
- Make `eos_mix` protected and settable only by `eos_set_mixture`.
- Carry the SCV and OPAL search state in the result record or a lookup record (M107), and put OPAL 1995 on the shared layout with clamped stencils (M94).

**Maintainability 7 -> 9**
- Fix the open latent defects. M93, M105, M106 and M107 are small, and the fixes for M105 and M106 are byte-identical. Fix M92 and M94 with their tests.
- Return the kernels' failures (OPAL density search, interpolators, Saha non-convergence) as codes, and report them in `opal_step`/`eos_get` (code style rule 7).
- Use one 1/mu_ion rule on every path (finding 2).
- Run test_eos under `make devel` (`-fcheck=all -ffpe-trap`) in CI, which would catch M92, M94 and M107 mechanically.
- Raise OPAL coverage (71%) with synthetic corrupt tables, and fix M119 and M87.

**Extensibility 6.5 -> 9**
- Correct the source recipe (finding 1 and M123): state the lookup's input contract (T, P, 1/mu_ion) and the fields the adopt and blend helpers carry. Better, let a source declare the fields it fills.
- With the descriptor, a new source becomes one module plus one list entry and its settings, instead of edits to eos_configuration, eos_lib, blended_eos, setup, io and the test. Re-trace the recipe from a fresh record.

**Docs accuracy 7.5 -> 9**
- Fix findings 1 and 6-9 and the documentation items M119, M122 and M123.
- Add a per-path table to §3: which source fills each `eos_result` field (density, mean weights, ion fractions, eta, second derivatives), and the Z and mixture rule each source uses.
- Date and re-measure the coverage figure in the same commit as the next pin change.
