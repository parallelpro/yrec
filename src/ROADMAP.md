# YREC phase-three roadmap: adopting MESA's module architecture

Written 2026-08-21, immediately after the phase-two domain sweep
finished (see GUIDELINES.md for phases one and two). The reference
architecture studied is MESA 26.04.1 (local copy at $MESA_DIR),
whose every physics module follows one identical skeleton:
`public/` (exactly `<mod>_def.f90` + `<mod>_lib.f90`, the only
surface other modules may touch), `private/` (all implementation),
`defaults/` (namelist controls), `other/` (user hook templates),
`test/` (standalone per-module test program), and `make/` (each
module builds its own static library). MESA even ships a
`package_template/` directory for stamping out new modules.

Phase two gave YREC the *semantic* half of this: one public
`<domain>_lib` facade per domain, domain folders with
pipeline/component subfolders, and derived-type state modules. This
roadmap covers the *mechanical* half -- what actually stops a file
from bypassing the facade, what makes a domain testable on its own,
and what makes errors survivable. Each stage below is independently
shippable and ends with the standing verification (full clean build +
Stage-0 byte-identical regression), plus, once stage 2 lands, the new
per-module tests.

**Explicit non-goal: handles / multi-instance support.** MESA's
`alloc_eos_handle` machinery exists so one process can run multiple
differently-configured stars (binaries). YREC runs one star per
process; the phase-one decision to use single module-level instances
(GUIDELINES.md, `prev_model` discussion) stands. Do not add handles.

---

## Stage 1 -- finish the public/private boundary (fix the bypasses)

The facades exist but nothing routes around-the-side calls through
them. Current inventory of bypasses (verified by grep, 2026-08-21):

1. `atm/turnover/calcad.f90` calls `eqstat2` and `esac06` directly
   (documented exception since the eos sweep). MESA's answer is
   `eosDT_get_component` -- a public accessor for *one named
   component's* results, bypassing the blend but not the boundary.
   Fix: add an `eos_get_component`-style entry to `eos_lib`
   (log T, log P in; the OPAL-2006 or eqstat2 result set out) and
   migrate calcad to it.
2. `wind/massloss.f90` calls `atm/tables/alsurfp.f90` directly.
   Fix: add a small `atm_lib` public entry (surface P/T lookup at
   given Teff/g -- the thing `alsurfp` does) and migrate.
3. `setup/setups.f90` performs three domains' table loads itself
   (`setupopac` for kap, `mhdtbl` for eos/MHD, `alfilein` for
   atm/Allard). MESA's answer is the `<mod>_init` lifecycle entry.
   Fix: add `kap_init`, extend `eos_lib` with an init that loads MHD
   tables when enabled, add `atm_init`; `setups.f90` then calls three
   facade inits instead of three internals. (`setupopac.f90` and
   friends move from de-facto-public to genuinely internal.)
4. `core/starin.f90` and `setup/hpoint.f90` call `kap/surfopac.f90`
   (refresh cached surface-composition table slices). Fix: make this
   a public `kap_lib` entry (`kap_update_surface_tables` or similar)
   -- it is a legitimate lifecycle operation, just unnamed as such.
5. `setup/grsett.f90` / `setup_grsett.f90` call
   `rotation/microdiff/` internals (`lax_wendrof1/2`,
   `get_imp_diffco`, `thdiff`). This is the legacy
   backwards-compatibility diffusion path. RESOLVED (2026-08-21, user
   choice): a third option beat the two anticipated here -- the pair
   is settling physics that lived in `setup/` only by name, so it was
   relocated into `rotation/microdiff/` itself, dissolving the bypass
   (the calls became intra-folder, symmetric with `microdiff`'s own
   use of the same kernels).
6. `core/main.f90` does `use opacity_table_lib` directly -- check
   what it touches and either route through `kap_lib` or drop the
   use-statement if vestigial.
7. State-boundary oddity: the OPAL-2006 EOS files
   (`eos/opal/*06.f90`) store their working arrays in
   `atm_table_lib` (historical accident of the original COMMON
   layout). Move those members to `opal_eos_lib` so eos state lives
   in eos's state module. Pure member relocation, byte-identical
   verifiable.

Deliverable: zero cross-domain calls that do not go through a
`<domain>_lib` facade, verified by a grep-based checker script (see
stage 2 for making that enforcement automatic).

Found during stage-1 execution (2026-08-21): `setup/setups.f90` also
read the **Fermi-Dirac table** (the degenerate-electron EOS table,
former `common/ccr/`) inline, into members that lived in
`atm_table_lib` but are consumed only by `eos/yale/eqrelv.f90` --
the same two-part pattern as the OPAL-2006 oddity. RESOLVED same day:
the `fermi_table_*` members moved to a new `state/yale_eos_lib.f90`
(the Yale/Prather EOS's own state module) and the load into
`eos_init`; the Kurucz/Castelli/Allard surface-table loads (the
inline reads noted under item 3) moved into a full `atm_init` at the
same time; and the SCV EOS table reads (a third inline setups.f90
block, spotted last) moved into `eos_init` as well -- order-safe
since the SCV read shares no file units or state with the atm reads
it now precedes, and directly regression-verified since the Stage-0
OPALSCV cases exercise the SCV tables. **Stage 1 is complete**:
every inventoried bypass and all three execution-time discoveries
are resolved, each verified byte-identical. `setup/setups.f90` now
performs no domain table I/O at all -- it computes constants and
calls `kap_init`/`eos_init`/`atm_init`.

## Stage 2 -- per-module standalone builds and tests

MESA compiles each module into its own static library and ships a
standalone test program per module; the library boundary *is* the
enforcement (private objects simply are not exported).

- Makefile: give each domain a library target (`libeos.a`,
  `libkap.a`, `libatm.a`, `libnuclear.a`, ...), link `yrec` from the
  libraries. The existing recursive-find single-build keeps working
  during the transition; the library split can go domain by domain.
- Add `test/` programs per domain, starting with eos: a small main
  that loads tables, evaluates `eos_get` (and `eos_get_component`)
  over a checked-in grid of (rho, T, X, Z) points, and byte-compares
  against stored expected output. Same pattern for kap
  (`kap_get` over a T/rho/composition grid) and atm (`atm_get` for a
  few (Teff, L) vertices).
- **This is what finally closes the LMHD coverage gap** -- with one
  caveat found during execution (2026-08-21): **no MHD table files
  ship with the repository** (`find` turns up only source/objects
  under `eos/mhd/`), so the gap splits into an infrastructure half
  and a data half. `test_eos` implements the infrastructure half: its
  MHD section runs `eos_get` under `use_mhd_eos` when the
  `YREC_MHD_TABLES` environment variable points at the 8 table files,
  and self-reports SKIPPED otherwise. The pinned baseline covers the
  shipped-data configuration only; actually pinning the MHD path
  needs tables obtained from outside the repo. No Stage-0 case sets
  LMHD, so until then the MHD path of `eos_get` (including the
  starin.f90 bug fix and the hsubp/sconvec/massloss MHD extension
  from the eos sweep) remains verifiable only in structure, not in
  numbers.
- Add the boundary-checker script from stage 1 to the test suite so
  facade bypasses fail CI rather than accumulating again.

## Stage 3 -- the ierr-not-stop error discipline

Inventory (2026-08-21): 90 `stop` statements in library domains --
eos 39, kap 29, atm 7, rotation 5, numerics 5, wind 2, mixing 1,
misc 2, nuclear 0. In MESA no library routine ever stops; everything
returns `ierr` and the *application* decides.

- Add an `intent(out) :: ierr` argument to each facade
  (`eos_get`, `kap_get`, `atm_get`, and the new stage-1 entries),
  defaulting to the current behavior at the top-level caller: on
  nonzero ierr, `core/main.f90` prints the context and stops -- so
  the observable behavior of a failing run is unchanged, but the
  library itself becomes embeddable.
- Convert leaf `stop`s bottom-up per domain (kap's table-miss stops,
  eos's out-of-range stops, ...), threading ierr through the
  intermediate calls. Do this per domain, one verified commit each.
- Sequencing note: this stage deliberately comes *after* stage 2
  because Stage-0 regression runs never trigger any of these stops
  (a passing run proves nothing about the error paths). The
  standalone per-module tests can deliberately feed out-of-range
  inputs and assert `ierr /= 0` with no crash -- the first time
  YREC's error paths become testable at all.

## Stage 4 -- named-index result arrays

`eos_get` currently has 27 positional arguments; MESA's `eosDT_get`
returns one `res(:)` array indexed by named constants (`i_lnPgas`,
`i_Cp`, `i_gamma1`, ... in `eos_def`). Adopt the same:

- Add the index-constant block (`i_pressure`, `i_density`,
  `i_grad_ad`, `i_cp`, ... plus a `num_eos_results` count) to
  `eos_lib` (YREC's analogue of `eos_def` -- a separate `_def` file
  is optional; the constants matter, the file split does not).
- New facade signature: inputs stay explicit (T, P, composition),
  the ~20 outputs collapse into `res(num_eos_results)` (+
  `d_dlnT(:)`/`d_dlnP(:)` arrays for the derivative outputs).
  `eqstat`/`eqstat2`/`meqos` keep their historical signatures; the
  facade packs/unpacks -- same thin-wrapper principle as phase two.
- Migrate the ~10 caller files; each site's local variable soup
  (`beta, beta_inverse, beta14, ...` declared in every caller)
  collapses to one array + named indexing, which is the real
  readability win.
- Then the same for `kap_get` (smaller: 4 outputs) if it proves
  worth it there.
- This is the churniest stage (every call site rewritten), which is
  why it comes after stage 2's tests exist. It remains byte-identical
  verifiable throughout: packing/unpacking the same doubles changes
  no arithmetic.

## Stage 5 (optional) -- the other/ hook system

MESA's `other/` modules hold null-default procedure pointers
(`s% eos_rq % other_eos_component => my_eos`) so users override
physics without editing the library. For YREC this is the lowest-
value item: the user base edits the source directly, and without
handles the pointers would live as module-level state. Only build it
if a concrete use case appears (e.g. an experimental EOS someone
wants to A/B without forking eqstat2). If built: one procedure
pointer per facade (`other_eos_get`, `other_kap_get`), checked
before the standard dispatch, null by default -- behavior-identical
when unset.

---

## Sequencing rationale and cross-cutting notes

Order is 1 -> 2 -> 3 -> 4 (-> 5) because each stage makes the next
one safer: 1 gives every cross-domain interaction a named facade
entry, so 2's library boundaries have nothing to strand; 2's
standalone tests are what make 3's error paths and 4's signature
churn verifiable beyond Stage-0's happy path; 4 is the largest
mechanical churn and should not go first on that ground alone. 5 is
optional and independent.

The deferred blanket-`SAVE` cleanup (GUIDELINES.md, "Deferred:
unneeded blanket SAVE statements") is independent of all five stages
and can interleave anywhere; the natural moment for each file is
whenever a stage already has that file open for edits.

Verification discipline is unchanged throughout: full clean build +
Stage-0 byte-identical diff for every commit, per-module tests
additionally from stage 2 on, and the boundary-checker keeps stage
1's win locked in. As with phase two, investigate before assuming:
each stage's first commit should re-verify this roadmap's inventory
against the then-current source rather than trusting these counts.
