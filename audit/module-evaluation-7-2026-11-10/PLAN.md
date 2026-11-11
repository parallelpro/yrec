# Plan: eos, kap, atm, convection, element_diffusion and core past 8 (2026-11-10)

From the seventh evaluation in this directory. Finding numbers are each report's "New findings":
E-n (eos.md), C-n (core.md), KA-n (kap_atm.md), CD-n (convection_diffusion.md). M, O and L numbers
are docs/bug_sweep.md's. **A** = byte-identical on every deck and pin; **B** = changes results or
pinned text, so it goes in its own measured commit (section 7).

## 1. Target

Each package's average above 8 at a fresh read-only evaluation with the seventh's rubric, and no
criterion below 8. Aim at 8.5: two reviewers of the same code differ by about half a point.

| Package | Now | R / O / M / E / D now | Aim, from the report's "what would raise" | What carries it |
|---|---|---|---|---|
| eos | 7.2 | 7.5 / 8 / 7 / 6.5 / 7 | 9 / 9 / 8.5 / 8.5 / 9 (~8.8) | the source descriptor; the latent fixes; kernels as named steps |
| convection | 6.4 | 6 / 7 / 6 / 6.5 / 6.5 | 8.5 / 8.5 / 8.5 / 8.5 / 8.5 (~8.5) | ok/FAIL tests; overshoot owned here; the prose |
| atm | 6.3 | 6.5 / 6.5 / 6 / 6 / 6.5 | 8.5 / 8.5 / 8.5 / 8 / 8 (~8.3) | ok/FAIL tests; lookups through arguments; one choice descriptor |
| kap | 6.2 | 6.5 / 7 / 5.5 / 6 / 6 | 8.5 / 8.5 / 8.5 / 8 / 8 (~8.3) | ok/FAIL tests; kap_messages; kap_eval as steps; kap_source |
| core | 6.7 | 7 / 7 / 6.5 / 6.5 / 6.5 | 8.5 / 8.5 / 8 / 8 / 8.5 (~8.3) | messages; envelope and Henyey tests; README; the F77 shape |
| element_diffusion | 5.9 | 5.5 / 6.5 / 5.5 / 6 / 6 | 8.5 / 8.5 / 8.5 / 8 / 8 (~8.3) | ok/FAIL tests; one pipeline (the original routines fixed or retired); the prose |

element_diffusion was added on 2026-11-10 (section 11, decision 4). Its report's list reaches about 8.1,
so its round adds the items marked "beyond the report's list".

core is the one at risk. Its report's own list reaches about 8.0 (maintainability and extensibility
"to 7.5-8"), so section 6 adds items beyond that list.

## 2. Why a groundwork round first

All five reports name the same missing middle layer:
- kernels that print (eos about 80 messages, kap 58, atm 115, convection 11, core 139);
- tests that pin printed numbers with no ok/FAIL check or mutation run (kap, atm, convection; core's two
  largest kernels have no test at all);
- READMEs without the flow, conventions, coverage figure and traced recipes;
- the F77 shape: `double precision`, `.lt.`, ALL-CAPS prose, history in the code.

Building the tools for these once makes each package round shorter. It also lets a checker hold the
gains, the way check_boundaries holds the architecture.

The deck-comparison harness also has to be rebuilt first. The macOS tmp cleaner has deleted the
scratchpad's `check.py`, its extra rotating decks and its baseline. The source survives in a session
transcript and the deck generator in another.

## 3. Round 0: groundwork (A only; one agent plus the integrator)

1. **The harness in the repository.** `src/tools/compare_decks.py` (today's check.py rules: build, run
   the 56 decks, byte-compare every output with a baseline, compare logs with the `standard/` pins).
   - The extra decks it needs go into `testsuite/`: the rotating variants `rv_*`, `cm_accretion`
     and `o24_rescale_rotating`.
   - The baseline goes to a git-ignored directory outside /private/tmp.
   - Newcomers get the same check the rounds use.
2. **bug_sweep §14.** Register the seventh evaluation's new findings for all 13 packages, so every later
   item has a number. Also register:
   - kap's nine and atm's five unregistered sweep defects (KA-4, KA-5);
   - the four original-settling defects (CD-3).
3. **`src/tools/check_style.py`, a ratchet.** It counts per directory:
   - `double precision`, `.lt.`-style operators, `endif`/`enddo`, `dmax1`-style intrinsics;
   - ALL-CAPS comment lines and commented-out code;
   - modules without `private`, `use` without `only:`, explicit-shape dummies;
   - I/O statements outside a `*_messages` module.

   A table holds each directory's allowed counts, which can only fall. It runs in `make check` and in
   CI, and reviewers get the counts without writing scripts.
4. **Kernels print nothing.** check_no_io extends to STAR_BLIND_FILES and to a per-package kernel list.
   The list grows as each package's messages module lands.
5. **`src/tools/mutate.py`.** It generalizes math's and numerics' mutate scripts, with one manifest per
   test program: file, text, replacement, and the check that must FAIL. Each round adds its manifest,
   so "confirmed by a mutation" can be re-run (eos maintainability).
6. **CI's devel job.** It builds the test programs with `-fcheck=all -ffpe-trap` and runs the unit tests
   (eos maintainability; it would catch M92, M94 and M107).
   - Run it locally first, and register what it finds.
   - A known failure is skipped by name until fixed, not hidden.
7. **The README standard as a checklist** in docs/GUIDELINES.md: rotation/README's sections, and the
   coverage command (`src/tools/coverage.py --dir <pkg>`) whose figure each README quotes with its date.
8. **The rubric into the repository.** Copy SP/notes/REVIEW7.md into audit/, so the eighth evaluation
   uses the same text.

## 4. Each package round: the same phases

One worktree and one agent per package, outside /private/tmp. The phases run in order, and each
commit is checked as in section 8:

- **a.** A-kind fixes (latent bugs, error paths).
- **b.** Structure: the messages module, state out of module variables, the descriptor, named steps.
- **c.** Style, down to the ratchet's zero for the package.
- **d.** Tests: ok/FAIL checks, the mutation manifest (every mutation caught), a coverage figure.
- **e.** The README to the standard, each recipe traced on a throwaway branch.
- **f.** The B list handed to the integrator for section 7.

Shared files are edited only by the integrator at merge; an agent lists what it needs there. The
shared files are bug_sweep, ROADMAP, the controls registry and generated files, the checker lists,
and other packages' root `test_*.py`.

## 5. Round 1: kap, atm, convection, element_diffusion and eos in parallel

The five directories are disjoint. atm's round also owns the three core call sites of its choice
descriptor (below). convection and element_diffusion share the composition step: element_diffusion's
round owns the settling sub-step policy's move out of core and rotation, convection's round owns the
overshoot's.

### kap (6.2 -> ~8.3)
- **a.**
  - KA-1: key the Alexander 2006 cache on the composition, and clear it on a failed build and at
    kap_init.
  - KA-2: check `iostat` on every open and return ierr; OPAL92 opens with `status='old'`. This also
    closes M170.
  - KA-6: validate the settings at kap_init.
  - The A items of KA-4 (all but 4a); M104; L10's gap.
- **b.**
  - kap_messages (KA-9), with the text unchanged.
  - OPAL95's scratch made local, and `opacity_table` protected (KA-12).
  - `kap_eval` as named steps returning a `kap_result`: molecular, interior, helium-core blend, ramp,
    conduction. It becomes private behind a test hook.
  - A `kap_source` descriptor (name, window, loader, lookup), so that family precedence is one ordered
    list instead of three if-chains (TARGET_ARCHITECTURE 4.5(c)).
- **c.**
  - One log R extrapolation helper for OPAL95's three interpolators, and one row/T kernel for LAOL89's
    two.
  - The misleading names: `log10_opacity`, `zlaol_logrho_grid`, alex95 against alex94.
  - Named thresholds and the conductive constants (-3.5194 = log10(16σ/3)).
  - Sources sections in place of the 104 history comments; KA-13, KA-14; retired/.
- **d.** ok/FAIL checks:
  - exactness at the table nodes;
  - continuity across the ramp and the family seams;
  - the derivatives against finite differences;
  - the Z blend at its ends;
  - the error contracts: a missing table, bad settings, two equal Zs.

  Plus the paths no test takes: OPAL95's 3-D path, `kap_update_surface_tables`, and conduction with
  nonzero ion fractions.
- **e.** README §3 (ierr and messages), a tests section, and KA-11's recipe traced to every file.
  With the descriptor, a family is one module, one list entry and its controls.
- **f (B).**
  - O18's ramp derivatives.
  - KA-4a, the conductive opacity's low-density sign (reached by every deck, negligible in κ).
  - KA-8's banner, "pure C and N" -> "pure C and O" (35 pinned logs move by one line).

### atm (6.3 -> ~8.3)
- **a.**
  - KA-3: an [Fe/H] or [α/Fe] the Allard file lacks returns ierr.
  - KA-6: `atm_choice` validated at init, so an unknown choice no longer runs HSRA silently.
  - KA-5's items: the Castelli/Kurucz sentinels above about 48 kK, the upper range check, the
    extension width, the gravity-search gap.
  - KA-7 (O19's wider half, no shipped deck reaches it); KA-12's unused Kurucz load for choice 4.
- **b.**
  - The lookups return through arguments, not through `atm_table`'s members (KA-12).
  - atm_messages (KA-9), with the text unchanged.
  - One home for what a choice is (tabulated or not, its fallback, its banner), asked by
    core/envelope/envelope_integration_kernel.f90:306, core/structure/surface_boundary.f90:136-174
    and io/derived_controls.f90:327-339.
- **c.**
  - allard_atmosphere's headers in today's names, and the typos.
  - `fail_9999` and the banners gone; the old layout's 54 and 5 named.
  - The T-tau constants named (¼ log10(3/4), log10 τ_ph); the dead print parameters removed.
  - Sources sections; retired/.
- **d.** ok/FAIL checks:
  - each T-tau relation gives T = Teff at its photospheric τ;
  - the tables are exact at their nodes and continuous across the fallback to gray;
  - the error contracts: an unknown choice, a missing [Fe/H], out of range.
- **e.** README §3, a tests section, and KA-11's recipes traced, including the core files.
- **f (B).**
  - O19, both halves.
  - The per-lookup notes dropped from the Allard run log (26% of it).

### convection (6.4 -> ~8.5)
- **a.**
  - CD-5: the twelfth convection zone returns ierr or is a stated limit, not silently dropped.
  - CD-11's formats.
  - The dead `only_check_core` branch and the unreachable fully convective test.
- **b.**
  - CD-14: the overshoot extent and H_P owned by convection/ and called from core and rotation. The
    identical copies are A; the differing ones go to f.
  - One zone scan and one merge.
  - Semiconvection's EOS-opacity-gradient triplet as one procedure.
  - convection_messages (CD-15), with the text unchanged.
- **c.**
  - `compute_temperature_gradients` as named steps (radiative gradient, spot correction, cubic,
    velocity, derivatives), with the derivative equations written once.
  - The 132 ALL-CAPS lines as sentences; history to Sources; the commented-out line removed.
  - `<`, `real(dp)`, `max`/`abs`/`sqrt`; the two misleading names.
  - Slots from composition_layout; named numbers (the 1e-11 velocity, the cubic's start, the 0.99e-3
    O16 floor, `vtol`).
- **d.** ok/FAIL checks:
  - ∇ between ∇_ad and ∇_rad;
  - the efficient and inefficient MLT limits;
  - the derivatives against finite differences, spots included;
  - homogenization conserves each species and leaves a uniform composition alone;
  - overshoot extents against α H_P;
  - a merge;
  - a semiconvective case that actually extends (CD-12: the cited deck is fully convective).
- **e.** A flow section (the composition step's passes, in order), the messages convention, the
  coverage figure, what each test checks, CD-13's recipe traced, and CD-12 and CD-16 corrected.
- **f (B).** CD-1 (semiconvection's core), CD-2 (the starspot Jacobian), CD-6, CD-14's differing
  copies, O10, M96, and CD-15's run-log text.

### element_diffusion (5.9 -> ~8.3)
- **a.**
  - CD-4: both pipelines search the model's zones only up to the equal grid's point count.
  - CD-8, CD-10 (microdiff's "iteration" repeats one identical solve), CD-9 (microdiff_grid's three
    hard-coded light elements, taken from species_table instead).
  - `light(3,json)` and settling_grid's literal slots 1-11 named.
- **b.**
  - One region finder, one Coulomb logarithm and one grid builder for both pipelines. The Bahcall &
    Loeb temperature unit and the year are defined once.
  - The settling sub-step policy (Michaud's timescale, the count) moved into the domain from core's
    settle and rotation's settle_helium.
  - element_diffusion_messages (CD-15), with the text unchanged.
  - **Beyond the report's list:** the original settling routines either fixed and tested or moved to
    retired/, per the owner's decision in its B list, so the domain has one pipeline that decks run.
- **c.**
  - The 353 ALL-CAPS lines as sentences; the legacy names in comments and the typos.
  - `private` in the four modules without it; `only:` everywhere; assumed shape for the 67
    explicit-shape lines; `real(dp)`.
  - The misleading scratch names (`metal_abundance_change` holds Z, ...); `data` statements as
    parameters; CD-16.
- **d.** ok/FAIL checks:
  - each species' total mass conserved, in a consistent cell convention;
  - X+Y+Z+He3 = 1;
  - the light elements sink;
  - the change is linear in a small timestep, and zero at a zero timestep;
  - Thoul's coefficients against the published fit at its sample points;
  - the Coulomb logarithm against Paquette's limits;
  - a case with a convective core (zone_begin > 1), the path CD-4 is on.

  M33 fixed, so that the pins are deterministic.
- **e.** A flow section (microdiff's stages), "which settings act where" (CD-17), the messages
  convention, the coverage figure, CD-13's recipe traced, CD-16 and CD-7's documentation.
- **f (B).** O21; M33 if it moves results; CD-7's change; CD-3's four defects if the original routines
  are kept; M181 (diffuse_lithium's overrides) with a recommendation.

### eos (7.2 -> ~8.8)
- **a.** M93, M105, M106, M107, M136, M137, E-2, E-3 (the range warning), E-4; M92, M94 and M128 with
  their tests; M87 and M144.
- **b.**
  - The table-source descriptor (name, window, reader, lookup, edge weight), with `blend_eos` looping
    over an ordered list. SCV becomes a source with an edge weight.
  - OPAL 1995 on the shared layout with clamped stencils (M94), deleting about 1,000 lines.
  - The test-only entries behind a test module; `eos_mix` protected.
  - The search state in a lookup record (M107).
  - The kernels' failures returned as codes and reported once, through eos_messages.
- **c.**
  - The OPAL, SCV and Fermi-Dirac kernels as named steps (the report's list).
  - Named thresholds: E-7's sentinel family; yale_eos's 1e-15, 4 dex and 1e-8; MHD's 5e-7 and 1e-4.
  - SCV's radiation correction written once.
  - E-8: `pure`, assumed shape, `mhdst1`'s record.
  - Sources sections; retired/.
- **d.** A mutation manifest for test_eos's 49 checks. OPAL coverage (71%) raised with synthetic
  corrupt tables. test_eos run under the devel job.
- **e.** M148, M156, M157, M122, M123, E-3, E-5, E-6; the per-path table in §3; §4's long paragraph
  split; the source and vintage recipes re-traced.
- **f (B).** E-1 where it is reached; M88 if the edge weight moves results; the capitalized messages;
  M129's 1/mu_ion rules once decided.

## 6. Round 2: core, after round 1 is merged

core goes after round 1 for two reasons. atm's choice descriptor and convection's overshoot change
core's call sites, and core's README recipes name them.

Two agents work on disjoint sub-directories, (envelope, structure, step) and (model, run,
diagnostics). The integrator handles main.f90 and the README.

- **a.**
  - C-1: bound the calibration's arrays at the 25th pair and stop with a verdict.
  - C-2: set the divergence flag on every path.
  - C-4: the C entry returns ierr for a path over 256 characters and closes only the units it opened.
  - C-7: the vestigial save.
- **b.**
  - core_messages (C-5's 139 writes). The two printing star-blind kernels go silent, and check_no_io
    covers them.
  - C-8: the star-member arguments removed.
  - Explicit-argument steps for the four largest hosts: `henyey_iterate`, `compute_timestep`, `rezone`
    and `read_starting_model`.
  - `private` in henyey_elimination, envelope_integration and its kernel.
  - The turnover time's refresh moved from observables into the step.
- **c.** The F77 shape, one commit per sub-directory:
  - `real(dp)` (459 lines), the operators (315), `end if`/`end do` (229), `only:` lists (44);
  - C-10's names, and the card kinds from io/run_cards;
  - C-12's names (`dlogl_dlogr`, `log_luminosity`, ...);
  - the history headers as Sources sections; C-9's commented-out code to retired/.
- **d.**
  - C-6: an ok/FAIL test of the envelope integration kernel (a case with a known solution, and its
    error contracts).
  - An ok/FAIL test of the Henyey elimination (a block system against a dense solve, and the
    singular-pivot contract).
  - Mutation manifests for these and for the four kernel tests that exist.
  - A coverage figure for core under the unit tests and the decks.
- **e.** The run's flow, the messages and negative-ierr conventions, every test with its mutation, the
  coverage figure, and four recipes, traced: an energy source, a stop condition, a timestep limit and
  a mesh criterion. C-11, C-12, and TARGET_ARCHITECTURE's Phase 2 status.
- **Beyond the report's list** (what lifts maintainability and extensibility past 8):
  - the timestep limits as kernels behind a settings record, so a new limit is one kernel, one line in
    `take_tightest_limit` and a test;
  - the species counts and composition dimensions named (the "literal 15");
  - the Stage B open items fixed in f.
- **f (B).**
  - C-3's warning text.
  - O4's out-of-range read, O6, O7 and M126, each measured.
  - M85 and M96.
  - M176, M190 and M191 once decided.

## 7. Kind B: measured commits, after each round's A work is merged

- They run one package at a time, so each reseed starts from a known base.
- Each item gets its own commit. The commit message gives the measured change on the decks and pins
  (the largest relative change per column). The pins are reseeded in the same commit, and the
  bug_sweep entry is closed there too.
- Message text changes (capitals, wrong words, dropped notes) go in one re-pin commit per package.
  Those commits change log text only; every other output file stays byte-identical.
- Each package's B list goes to the owner before any of it starts.

## 8. Checks for every commit and merge

- compare_decks: ALL IDENTICAL for A commits, and the 35 pins OK.
- `make tests`; pytest exact and with `YREC_TEST_TOLERANT=1`.
- check_boundaries, check_star_contracts, check_style (counts only fall), and the generators' `--check`.
- The package's mutation manifest: every mutation caught.
- Each round is merged into yrec-next with `--no-ff`. Pushes happen on the owner's word only.

## 9. The eighth evaluation

Fresh read-only reviewers score the six packages with the copied rubric and the same report format,
in a new audit/ directory. A package still under 8 gets its "Left" list as a short follow-up round.

## 10. Order and size

| Step | Work | Agents | Size |
|---|---|---|---|
| 0 | groundwork, the harness first | 1 + integrator | M |
| 1 | kap, atm, convection, element_diffusion, eos | 5 in parallel | L each; eos's descriptor and OPAL 1995 merge are the largest items |
| 1B | round 1's measured fixes and re-pins | integrator, serial | M |
| 2 | core | 2 in parallel | L |
| 2B | core's measured fixes | integrator, serial | M |
| 3 | eighth evaluation | 6 read-only reviewers | S |

## 11. Decisions for the owner

Answered on 2026-11-10: yes to all four.
- Kind B fixes are allowed. Each package's list still goes to the owner before its B work starts.
- Message rewording is allowed, with one re-pin commit per package.
- M190 and M191 are fixed in core's B work. M176 and M129 come to the owner with a recommendation
  in core's and eos's B lists.
- element_diffusion joins round 1 (section 5). Its original settling routines' fix-or-retire choice
  goes in its B list.

The questions as they were asked:

1. Kind B fixes in these rounds, as measured commits, from a list approved per package before the
   start? Recommended: yes. Without them kap's and core's maintainability stay under 8.
2. Message rewording with re-pins, one commit per package? Recommended: yes. The reviewers marked
   the capitalized, stale and wrong messages down in all five.
3. core: M176, M190 and M191 (recommended: fix M190 and M191). eos: M129's 1/mu_ion rules.
4. element_diffusion (5.9) and rotation (8.0) are outside this plan.
   - element_diffusion shares the composition step and the test approach with convection, so it would
     add one package to round 1.
   - Its original settling routines would then need a fix-or-retire decision.
