# compare_decks' extra decks

`src/tools/compare_decks.py` runs these 17 decks with the shipped ones (its
`quick`, `rot` and `all` sets). Each takes a path that no pinned deck takes,
so a refactor that changes what that path computes shows up as a byte
difference against the baseline. They have no pins: only the baseline
checks them.

They live here, not in `testsuite/` or `examples/`, because pytest must not
collect them: `test_all.py` reads only the directories `pytest.ini` lists
(non-recursively), but `test_controls.py` and `test_generated_files.py`
search `examples/` and `testsuite/` recursively, and `test_controls` would
fail on a deck that has no row in `src/io/test/expected_controls/hashes.txt`.

Each deck's first line names its base deck and the settings it changes. The
paths are the shipped decks' `{YREC_INPUT}`, `{YREC_START}` and
`{YREC_OUTPUT}`; compare_decks runs each deck from this directory with
`YREC_INPUT` and `YREC_START` set to the tree's `input/` and `startmodels/`,
so the output goes to `output/` here (git-ignored). By hand:

```
cd src/tools/compare_decks
YREC_INPUT=$PWD/../../../input YREC_START=$PWD/../../../startmodels ../../yrec rv_ref.inlist
```

## The rotation variants (rv_)

`rv_ref` is `testsuite/Test_solarGS98_difrotmixper18p4d1p0.inlist` (a solar
model with GS98, diffusion, differential rotation, rotational mixing, disk
locking and the Matt et al. 2012 wind) stopped at 3e8 yr in place of 4.568 Gyr.
Each other `rv_` deck is `rv_ref` with one rotation switch changed:

| Deck | Setting (rv_ref's) | Path it exercises |
|---|---|---|
| `rv_ref` | (none) | the reference the others differ from |
| `rv_zahn` | `velocity_form_correction = .TRUE.` (F) | composition mixing scaled by the Zahn (1992) coupling factor |
| `rv_gsf0` | `gsf_inhibition_mode = 0` (1) | the GSF and ABCD instabilities' other stability conditions |
| `rv_gsf2` | `gsf_inhibition_mode = 2` (1) | as rv_gsf0 |
| `rv_sbm1` | `solid_body_mode_flag = 1` (0) | solid-body rotation everywhere |
| `rv_sbm2` | `solid_body_mode_flag = 2`, `cz_rotation_law_exponent = -0.5D0` (0, 0.0) | solid-body core with a decoupled convection zone |
| `rv_sbm3` | `solid_body_mode_flag = 3`, `cz_rotation_law_exponent = -0.5D0` (0, 0.0) | solid-body core coupled at the convection zone's base |
| `rv_alpham05` | `cz_rotation_law_exponent = -0.5D0` (0.0) | a differentially rotating convection zone, omega ~ r^-0.5 |
| `rv_nocore` | `no_am_transport_in_core = .TRUE.` (F) | no angular momentum transport in the convective core |
| `rv_nobg` | `use_constant_background_diffusion = .FALSE.` (T) | transport without the constant background diffusion |
| `rv_basic` | `use_extended_composition = .FALSE.` (T) | the basic composition (no H2, Li6, Li7, Be9) under rotation |
| `rv_iter1` | `max_diffusion_iters = 1` (10) | one iteration per transport sub-step |
| `rv_nodisk` | `disk_locking_active = .FALSE.` (T) | the wind from the start, no disk locking |
| `rv_kawaler` | `use_pmm_wind_law = .FALSE.` (T) | the Kawaler wind law |
| `rv_notrans` | `instability_transport_active = .FALSE.`, `disk_locking_active = .FALSE.` (T, T) | rotation without instability transport: the wind-only spin-down |

`cz_rotation_law_exponent` is the inlist name of what `rv_ref` sets as
`new_zone_omega_weight` (an alias the inlist format still accepts; the
inlist name wins when a deck sets both).

## The other two

| Deck | Base deck and settings | Path it exercises |
|---|---|---|
| `cm_accretion` | `examples/run_config_matrix/cm_massloss_const.inlist` with `mass_change = 1.0D-9` (-1.0D-9) | constant mass gain; the configuration matrix has only loss |
| `o24_rescale_rotating` | `testsuite/Test_solarGS98_difrotmixsolidper18p4d1p0.inlist` with `max_age(2) = 1.0D8`, `instability_transport_active = .FALSE.` and two more cards: a rescale (`run_kind(3) = 2`, every target negative, so kept) of the rotating model, then 5 models | the un-aged models of a rescale card, with the wind but no instability transport (docs/bug_sweep.md O24) |

## Where they come from

They were made on 2026-10-28 for the rotation rounds, and remade on
2026-11-10 from the same base decks by the same recipe: one setting
changed on its existing line, its comment kept. The remade decks' `.log`
and `.mod` files are byte-identical (version line aside) to those the
2026-10-28 decks wrote on 2026-11-02. A `test_controls` dump of each deck
differs from its base deck's in exactly the settings above (for `rv_ref`
the age limit; for `o24_rescale_rotating` also the two added cards), and
every `rv_` deck's `.log` and `.mod` differ from `rv_ref`'s and from each
other's.
