# Module evaluation 8, 2026-11-17 -- the six packages rounds 0 to 2B worked on

Six read-only reviews at yrec-next 1179afbf (pushed), out of 10 against the bar of the seventh
evaluation (RUBRIC.md here: 10 = at or above the standard rotation/ set with nothing material left;
9 = only minor items left; 8 = a few real gaps; 7 = several; below 7 = structural problems). They
cover the packages that rounds 0, 1, 1B, 2 and 2B worked on since the seventh evaluation (823326a7):
atm, convection, element_diffusion, eos, kap and core. core has two reports: core_a.md (structure/,
step/, envelope/) and core_b.md (model/, run/, diagnostics/, and core's README as a whole). The
other seven packages were not reviewed again and keep the seventh's scores.

The reviewers read only and ran nothing. Each of their findings was then read again in the code,
and against the F77 import, before it was registered (docs/bug_sweep.md section 19, M319-M355);
FINDING_MAP.md maps every finding to its number and says what that check corrected.

| Package | Fourth (10-31) | Sixth (11-02) | Seventh (11-10) | Now | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|---|---|
| setup | 8.2 | 9.2 | 9.5 | 9.5 (not re-scored) | 10 / 9 / 9.5 / 9.5 / 9.5 |
| io | 7.5 | 8.9 | 9.3 | 9.3 (not re-scored) | 9.5 / 9 / 9.5 / 9.5 / 9 |
| math | 9.3 | -- | 9.3 | 9.3 (not re-scored) | 10 / 10 / 8.5 / 10 / 8 |
| state | 7.8 | 8.9 | 9.1 | 9.1 (not re-scored) | 9.5 / 9 / 9 / 9 / 9 |
| atm | -- | -- | 6.3 | 8.6 | 9 / 8.5 / 8.5 / 8.5 / 8.5 |
| convection | -- | -- | 6.4 | 8.6 | 8.5 / 8.5 / 8.5 / 8.5 / 9 |
| numerics | 8.5 | -- | 8.5 | 8.5 (not re-scored) | 9 / 9 / 7.5 / 9 / 8 |
| eos | 7.3 | -- | 7.2 | 8.4 | 8 / 8.5 / 8.5 / 8.5 / 8.5 |
| net | 8.2 | -- | 8.1 | 8.1 (not re-scored) | 9 / 8.5 / 8.5 / 7.5 / 7 |
| rotation | -- | -- | 8.0 | 8.0 (not re-scored) | 7.5 / 8 / 8 / 8.5 / 8 |
| kap | -- | -- | 6.2 | 7.9 | 7.5 / 8.5 / 7.5 / 8 / 8 |
| element_diffusion | -- | -- | 5.9 | 7.9 | 8 / 8 / 7 / 8 / 8.5 |
| core | -- | -- | 6.7 | 7.7 | 7.6 / 7.7 / 7.4 / 7.8 / 8.2 |

core is the mean of its six directories, and its five criteria are the means of theirs:

| core directory | Seventh | Now | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|
| structure | 7.7 | 8.1 | 8 / 8 / 8 / 8 / 8.5 |
| run | 6.8 | 8.0 | 8 / 8 / 7.5 / 8 / 8.5 |
| step | 7.0 | 7.7 | 7.5 / 7.5 / 7.5 / 8 / 8 |
| model | 6.6 | 7.6 | 7.5 / 7.5 / 7 / 8 / 8 |
| diagnostics | 7.0 | 7.6 | 7.5 / 7.5 / 7 / 8 / 8 |
| envelope | 6.2 | 7.4 | 7 / 7.5 / 7.5 / 7 / 8 |

core's README and in-code documentation, the whole package: 6.5 -> 8.5 (the README alone 9).

## Reading the table

- Every package the rounds worked on rose: atm by 2.3, convection by 2.2, element_diffusion by 2.0,
  kap by 1.7, eos by 1.2 and core by 1.0.
- atm, convection and eos are at the standard: none of their criteria is under 8. What keeps each
  from 9 is small and listed in its report (for eos: four kernels that keep their F77 shape inside,
  and solvers that send their own messages).
- kap (7.9) has the standard's structure: one facade, a descriptor per table family in one
  precedence list, private state, checked settings and reads, a messages module, organization 8.5.
  Its readability and maintainability are at 7.5: the seven lookups under the facade are the F77
  routines in modern syntax, and the tests check each family's derivative once, not the tables'
  edges, the readers' header checks or the shared open helper.
- element_diffusion (7.9) rose as far as byte-identical work could take it. Its maintainability is
  at 7 because of its open defects, and those are a decision: the owner froze its measured fixes on
  2026-11-11. The score says what the code is, as the seventh's did.
- core (7.7): structure/ (8.1) and run/ (8.0) are at the bar. step/, model/, diagnostics/ and
  envelope/ are under it for the same few reasons: long hosts whose internal steps share their
  locals (`evolve_step`, `rescale_model`, `turnover_from_interior_new`, the 870-line envelope
  kernel), drivers without a test of their own (model/'s three drivers, all of diagnostics/, the
  Li/Be burn), the dead accretion-entropy path, and the defects below.
- The docs criterion is at 8 or above in every package reviewed. Each README has the standard's six
  parts and matched the code where the reviewers checked it.
- The seven packages not reviewed are where the seventh left them. rotation (8.0), which set the
  standard, net (8.1) and numerics (8.5) are now below atm and convection; rotation's and net's
  lists in the seventh's SUMMARY and PLAN stand.

## What the rounds' claims looked like when sampled

The rubric told the reviewers not to grade on the rounds' own claims, and to check a sample of what
docs/ROADMAP.md says was done against the code.

**What held.**
- core (structure/, step/, envelope/): C-2, C-7, the messages route, `private` in every module, the
  Henyey kernels and their test, C-8, `henyey_iterate`'s steps, the timestep kernels, the F77 shapes
  at 0, retired/, the species counts, O5, O7, O8, O9, O19's core side, M132 and M273's part B. Each
  of the 179 mutations' `old` strings occurs exactly once in its file.
- core (model/, run/, diagnostics/): sixteen ROADMAP claims held as stated (M34, M85/M96, M86, M97,
  M115, M126, M190, M191, M268, M270, M271, O4, O22, O6's part A, M31's narrow form, the turnover
  refresh in the step). `check_boundaries`, `check_star_contracts` (209 declared) and
  `mutate.py core --check` (179 of 179 apply) pass. The "Add a stop condition" recipe, traced again,
  names every hand-edited file.
- eos: the fixes of 24 bug_sweep entries are where their entries say; the manifest's 82 mutations
  apply and cover all 76 checks; OPAL 1995 runs on the shared code; the README's deck figures hold
  (310 decks: 304 with OPAL 2006, SCV and Debye-Hueckel, 2 SCV only, 4 Yale only).
- kap: the fixes of M219, M222-M226, M228-M230, M236, M277, M278 and O18's ramp are in the code as
  described; 42 mutations and 37 ok lines; I/O statements only in the readers, the open helper and
  the messages module; none of the 44 reads and 8 opens is unchecked.
- atm: the kernels take their table and return their result; the tables are private to atm_lib and
  reset at each `atm_init`; 39 checks, 45 mutations, each mutation's text found once and each check
  named by one.
- convection and element_diffusion: one zone scan and one merge, the gradient's derivatives written
  once, overshoot.f90's ownership of H_P and the envelope extent, M241, M242, M248, M249, M285, the
  dropped `diffuse_helium`, the single implicit solve, the table-sized light rows and the sub-step
  plan in both callers: each as described.

**Three claims that did not fully hold** (all core's; registered as M332 and M333):
1. **core/README.md:143, "Each check of the nine programs with a manifest was confirmed to fail
   under a mutation of the code it covers."** It holds for the programs whose manifest names every
   check (test_henyey_elimination, test_timestep_limits, test_regrid, test_model_kernels,
   test_run_contracts). It does not for four: test_triangle_fit has 15 mutations for 83 checks,
   test_mass_change 14 for 78 printed checks (12 of 21 distinct ones), test_structure 10 that name 8
   of its 15 partials (the radiative `qt` partials have none), test_envelope_integration 32 for 35.
   The manifests' own headers say "each family of checks", which is what was run.
2. **"so that the diagnostics feed no physics"** (ROADMAP's round 2; observables.f90's header;
   core/README.md:66 and :101). `renormalize_luminosity_breakdown`, inside `compute_observables`,
   rescales `star%luminosity_breakdown_lsun` after the step's `compute_timestep` has run. The model
   file stores the rescaled values, and the next card's first `compute_timestep` reads them (the
   hydrogen and helium burning limits, the sign of L_grav). The factor is 1 to within the solver's
   last correction, so the effect is small; the statement is still not true.
3. **Item 24's message form.** core/README.md:114 and core_messages.f90:17-19 state it ("ERROR: " or
   "WARNING: ", then the source, in sentence case; "; run stopped" on an error that stops the run).
   Seventeen messages of structure/, step/ and envelope/ and ten of model/, run/ and diagnostics/
   do not have it, and four events reach the terminal and not the run log (the solar calibration's
   Newton step among them).

**Smaller, of the same kind.**
- Five README coverage figures predate the rounds that followed them: core's (2026-11-13, before
  round 2B), and kap's, atm's, convection's and eos's (2026-11-10, before round 1B).
- Two recipes were traced before a later commit changed what they name: eos's "Add a table source"
  gives `edge_weight` the signature it had before M88, and atm's two recipes miss the hand-written
  lists of the choices that the rewording of 2026-11-11 added.
- kap/README.md:76 still describes the conductive mixture as it was before M278.
- docs/bug_sweep.md listed M308 as open; ebf4d507 fixed it (corrected with this registration).

## New findings worth the authors' attention first (details in each report and in bug_sweep §19)

Most are latent for the shipped decks, by the reviewers' deck counts and the ones repeated for the
registration. Two are not: M325 acts at the first step of every card, and M324 at every model that
takes its branch; what either moves in a deck was not measured. Whether a shipped deck meets M320's
or M321's trigger was not checked. "F77" = inherited. Each was confirmed in the code, and four by a
scratch run as well (M340, M341, M345, M354).

- **core/step and core/structure:**
  - an accreting, disk-locked model older than the disk-locking age (whose default is 0) takes one
    step of negative length (M319; F77);
  - a non-positive luminosity at the core edge gives a negative burning limit, which wins the
    minimum and becomes the next timestep by `abs` (M320; F77);
  - an Allard deck whose gravity leaves the tables while its Teff stays inside switches to gray and
    back at every model, with three message lines and a new triangle each time (M321; the switch
    F77's, the new triangle 2026's; also atm's finding 2);
  - the envelope kernel reads locals nothing set, defined only by gfortran's `-finit-local-zero`
    (M322; 2026);
  - accretion spends an atmosphere lookup or integration and two EOS calls per step on an entropy
    path whose results nothing reads, but whose failures stop the run (M323; F77).
- **core/diagnostics, core/run and core/model:**
  - the turnover time's innermost-point branch is overwritten by the bottom cell's crossing time: a
    2026 goto-elimination regression that the first bug sweep recorded in August and nobody
    registered; a comment now describes it (M324);
  - each card's first step runs with the turnover time zeroed after the starting model computed it:
    no deuterium cap, and a wind that interpolates toward zero (M325; F77);
  - the core extension divides 0 by 0 for one added point, which is the default factor's case, and
    never shifts the density profile (M326; F77; all 213 shipped decks that extend the core set the
    factor to 10);
  - the stitched model swallows an envelope-integration failure after its ERROR line (M327; 2026);
  - the star calibration reports a crossing at the second model of every pass when the model starts
    above R* (M328; F77);
  - two accesses beyond what was set, at the walk's last points and at a json-point model (M329,
    M330; F77), and a NULL path from a C host (M331).
- **kap:**
  - kap's copy of the envelope Z is the run plan's at `kap_init` and never follows the model, so an
    OPAL 1992 or LAOL 1989 deck with one Z table is refused at its first interior lookup unless a
    card rescales Z to the model's own (M338; 2026; no shipped deck uses either family);
  - OPAL 1995's surface table puts the X = 1 − Z table at X = 1, and at Z = 0.04 or 0.10 reads an
    X = 0.95 slot that was never built (M339; F77; the defect M225 fixed in Alexander 2006);
  - the conductive opacity is NaN with ierr 0 for a gas with no electrons, and its density
    derivative is not the value's where a species' shifted density is clamped: both from M278's
    rewrite in round 2B (M340, M341; run: the derivative is off by 0.06 to 0.11 in the clamped band);
  - the Alexander lookups extrapolate beyond their tables with ierr 0 (M342), and `kap_init` still
    accepts a table file whose Z is not the control's (M343).
- **atm:**
  - the Kurucz lookup has no per-row log g limit: a point below a row's first valid column is
    extrapolated with ierr 0 (M345; F77; run: log Teff 3.93, log g 0.5 returns a pressure from rows
    that start at log g 1.0 and 1.5);
  - four unguarded paths of the Allard reader (M346).
- **eos:**
  - the MHD reader reads a composition record past its arrays before it checks the count (M334;
    F77; no MHD tables ship);
  - two table-shape assumptions that no reader checks (M335);
  - SCV's two interpolations differ by up to 5e-4 dex in log rho at the composition handover, and no
    list says so (M336; F77).
- **convection:**
  - a write one element past the flag arrays for a model of exactly 5000 shells (M348; F77);
  - the overshoot merge is not the union when the upper region reaches below the lower one's first
    shell (M349; F77; M285 fixed the same at the top);
  - semiconvection's propagation limit never applies to the first shell above the core (M350; F77).
- **element_diffusion** (recorded only; the package is frozen): a 0/0 in the CNO tracers of a zone
  with Z = 0 (M352), the original routines' derivative factors (M353), a one-zone settling region
  returned as found (M354).

## What is left under 8

Per package, the reviewers' "Left" lists for the criteria under 8, as a short follow-up round's work
list. Kind A = no result changes (byte-identical on every deck and pin); B = results or pinned text
change, a measured commit of its own. The bug_sweep numbers are section 19's.

### kap (readability 7.5, maintainability 7.5)

Readability, all Kind A:
- `compute_opal95_opacity` as four named steps (locate T and R, the high-R corner test as its own
  procedure, the weights in R, the 2D/3D/4D choice);
- `compute_kurucz_opacity` with one header and the row search as a procedure;
  `compute_conductive_opacity` with one `species_conductivities(logrho)` for its two branches;
- opal95_reader on `read_failed`, without labels; alex06_reader's table body written once;
- today's names for the 23 lines of F77 locals in opal92, alex06, alex94 and kurucz;
  `opal92_surface_slot` and the other unnamed values, with unit and meaning;
- the internal procedures given their own locals; the interim "(formerly ...)" names dropped and the
  F77 ones moved to the Sources sections;
- `pure` on the lookups and builds that qualify; the three unused imports.

Maintainability:
- M338, M339, M340 and M343, each with its check and mutation (A for the shipped decks);
- the checks kap's finding 7 lists (A): `open_table` (with the conductive reader moved onto it),
  OPAL 1995's high-R edge, the precedence list, one doctored header per reader, the single-Z
  contract;
- M298 (A with the tables the decks use) and, as one measured commit, M299 with M341 (B, round-off
  in every deck's photosphere);
- M342 decided: a refusal after a count over the decks (A if none goes off the tables), or the
  README's statement (A);
- M112's kap part (A);
- the coverage re-measured, with the lookups above 85% (A).

### core/step (readability 7.5, organization 7.5, maintainability 7.5)

- M323, the dead accretion-entropy path moved to step/retired/ (B: the EOS's carried state may move
  cm_accretion at round-off);
- `apply_mass_change`'s and `advance_composition`'s steps with arguments; `evolve_step`'s stages
  over a small step record in place of 25 shared locals (A);
- burn_step's and step_driver's thresholds named; the 17 unused names; M333's three stale
  step_driver comments (A);
- the Li/Be burn as a kernel (the rates and the base as arguments), tested for the closed-form
  survival, the mixed average and the conservation of the accreted share (A);
- convection's `homogenize_region` reused; one atomic-weight table (A where the values agree, B
  where they do not); `burn_network_region` through core_messages (A; its text with M332, B);
- M319 and M320, each with its check and mutation (B where they fire, byte-identical otherwise);
- a check of the shell branch's dX/dt limit (A); M89 decided (B); M146's rest;
- the timestep recipe's four missing lines, and recipes for a new kernel with its test program and
  for a new stage of the step (A; M333).

### core/envelope (readability 7.0, organization 7.5, maintainability 7.5, extensibility 7.0)

- the kernel's host split into the three parts its header names, as module procedures over one
  integration-state record (the 70 locals), with one point store and one Gamma1 (A);
- `point_physics`' fields renamed; `atm_get`'s vertex group as one optional record and its inputs
  `intent(in)` (A);
- M322 (A); `cfg` as `intent(in)` with the choice in use as a local; `taucz_env_accum` removed; one
  named no-photosphere sentinel (A);
- the refit's interpolation and the rebuild's shell masses as kernels, with a unit test of
  `append_envelope_points`' kernel: continuity at the fitting mass, the two branches, the degenerate
  interval (A);
- checks of a tabulated surface and of the fallback to gray, on a synthetic table through atm's
  test entry (A), with M321 (B for an Allard deck that reaches the tables' edge);
- a cap on `atm_retry` with its event and error (A unless a deck reaches the cap); M309 (B: a
  diagnostics line); M310 (A, descriptions);
- a traced recipe, "record a quantity at each envelope point" (A);
- M332's envelope texts (B, message text).

### core/model (readability 7.5, organization 7.5, maintainability 7.0)

- `rescale_model`'s ten steps with explicit arguments and a state record, as `rezone` has (A);
- the step bodies re-indented and their F77 names replaced, one mechanical commit per file (A); one
  `add_points_for_jumps` for the X and Z copies (A); names for the species ranges, the zone-bound
  sizes and the thresholds (A); the eight lines over 132 columns (A);
- one shell-mass kernel (regridding's), so that model/ calls diagnostics/ only for the starting
  model's stitch and turnover time, and those two calls made by the caller (A); the mesh criteria's
  arithmetic on arrays (A);
- M326 (A for the one-point case and the control check; B for the density shift) and M330 (A);
- in test_model_kernels (A): the mesh pieces that need no star once they take arrays, the core
  extension for N = 1 and N = 13, the horizontal-branch mass and core rescales;
- a pinned deck that rescales a horizontal-branch model (A: a new pin);
- M202 (B); M316 and M317 (not traced); O6's and M176's parts B (the owner's schedule).

### core/diagnostics (readability 7.5, organization 7.5, maintainability 7.0)

- `turnover_from_interior_new`'s steps with explicit arguments and a result record, star-free, with
  a name without "new", moved beside the step that owns it (A);
- `build_stitched_model`'s locals named; observables' unnamed values, and a source for the
  capture-rate tables (A);
- the luminosity renormalization moved into the step before `compute_timestep` (B, at the solver's
  tolerance), or documented where it is (A; M333);
- M324 (B: every model that takes the branch; measure which decks do), M327 (A where no deck fails
  there, B for a run that now continues past the error), M329 (A);
- `test_turnover` on synthetic profiles and `test_stitched` (the index map, the acoustic radius of a
  constant sound speed, a missing column's error, M327's error path), each with its mutations (A);
- M139.

### core/run (average 8.0; maintainability 7.5)

- M325 (B: the first step of each card) and M328 (A), each with its check;
- checks of `check_solar_calibration`'s Newton step against its 2x2 system and of the abundance
  stops' disarm and trigger (A);
- M312 and M331 (A); the Monte-Carlo files' handling moved to io (A).

### element_diffusion (maintainability 7), if the owner lifts the freeze

Each a measured commit or a decision; nothing here is to be done while the freeze stands.
- M33, so that every pin is deterministic (B: the pin "old: without the Thoul et al. coefficients");
- O21 and M247, with a convective-core case that reaches M247's path and a 2- and a 3-point grid
  case (B: M247's fix moves a pinned deck);
- M243-M246, M288 and M353 fixed with a check each, or the original routines retired with their
  pins (B: test pins only; the retire-or-fix decision);
- M287 and M257's 2.2 (B: every settling deck, about 2e-8);
- M352 and M354 (A: guards, byte-identical on the shipped decks);
- M289 and M250 (decisions; B);
- checks of the microdiff coefficients with their mutations (A): D1 proportional to r² T^(5/2),
  `fgry` and `fgrz` linear, the hydrogen-gradient term's sign, a settling timescale against
  Michaud's estimate; the bundled checks split so that each assertion has its mutation (A);
- the messages reworded (B, text; B list item 6); M258's two descriptions and M355 (A).
