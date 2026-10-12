# io/ — controls in, output and models out

This directory holds the run's input and output:
- the controls reader, which turns a deck into `star%ctrl` and `star%job`;
- the output: the run log, the history and profile files, and the pulsation files;
- the stellar model files.

The physics never reads a deck or writes a file itself. It reads its controls from `star%ctrl` and `star%job`, and the output reads what the step computed.

---

## 1. What lives where

**The controls reader**

| Module (file) | What it holds |
|---|---|
| `controls_reader` (`read_controls.f90`) | `read_controls`, the entry point, and `read_input`, the read in order (§2). Also the namelist read itself (`read_namelist_files`), the copy to the canonical names (`adopt_canonical_names`) and the table paths' load and store. |
| `legacy_controls` | The variables the namelists bind to, under their legacy spellings; the legacy `/control/` and `/physics/` groups; the new-format groups (generated `inlist_new_decl.inc`); and `set_legacy_defaults`, which gives each variable its default at the start of every read. |
| `controls_buffer` | Every control under its canonical name, with its default as the initializer: the buffer the deck is read into. Also the named constants `max_runs`, `itime_*` and `ichi_*`. |
| `run_files` | The fixed unit numbers, the `{YREC_INPUT}`/`{YREC_START}`/`{YREC_OUTPUT}` placeholders in paths, the output directory with its copy of the deck (`inlist_used`), and opening the run's files. |
| `derived_controls` | `derive_controls`: the controls computed from the ones read (rate scales, which EOS and opacity tables, the wind law's exponents, the mixture changes, the diffusion switches), the semiconvection-with-overshoot stop, and the run log's header. |
| `run_cards` | `interpret_kind_cards`: the run plan from the kind cards (evolve, rescale, rescale and evolve). |
| `controls_check` | Warnings about controls that are set but inert. run_yrec calls it once, after `star_setup`. |
| `io_units` | The three fixed units written from everywhere (the run log, the terminal and the final model), and `check_io`, the report for a failed open, read or write. |

**The output and the models**

| Module (file) | What it holds |
|---|---|
| `yrec_output` | The coordinator: `output_init_mesa` opens the job's files, `output_run_header` writes each card's header, and `output_write_model` decides each model what is due. |
| `history_output`, `profile_output` | The history and profile files in MESA's layouts: their column tables and writers. |
| `output_columns` | The column selection both use (`history_columns_file`, `profile_columns_file`). The defaults are in `defaults/{history,profile}_columns.list`. |
| `run_log` | The progress line, and the gate for the solver diagnostics. |
| `pulsation_output` | The FGONG, GSM (HDF5) and GYRE writers. Their columns are ported from MESA. |
| `model_files` | The `.mod` model writer and reader, and the readers of the older YREC7 and MODEL2 starting models. |

## 2. How a deck is read

A deck is either a single-file inlist (`&star_job` and `&controls`) or a legacy CONTROL/PHYSICS namelist pair (`.nml1`/`.nml2`). `read_controls` runs:
1. `star%ctrl = controls_state()`: the pristine defaults, generated from `controls_buffer`'s initializers.
2. `seed_controls_buffer`: buffer ← `star%ctrl`.
3. `read_input`:
   - `set_legacy_defaults`, `assign_unit_numbers`, and loading the table paths;
   - `read_namelist_files`: the namelist read, over the defaults;
   - `adopt_canonical_names`: each legacy-spelled variable copied to its canonical name;
   - `map_user_inputs` (`setup/`): the readable timestep and zoning inputs mapped onto the code's vectors, and the cross sections from the S-factors;
   - `store_controls_to_star` (the output setup reads `star%ctrl`);
   - `expand_paths`, `open_run_files`, then storing the table paths;
   - `derive_controls`;
   - `interpret_kind_cards`.
4. `store_controls_to_star`: `star%ctrl` and `star%job` ← buffer.

Every read starts from the same state, so an in-process rerun (pyyrec, `test_reentry`) reads exactly as a fresh run does.

**Errors.** Every open and namelist read in io/, and every output writer's open, checks its status. `check_io` (`io_units`) reports `cannot <action> <path>: <reason>` on the terminal, and in the run log once it is open, and returns ierr = 1, which stops the run with exit status 1. A misspelled control, for example, gives `cannot read &star_job from deck.inlist: Cannot match namelist object name ...`.

**A control's three names.** A namelist matches a deck's entries by variable name, so the legacy decks fix the legacy spellings (`KTTAU`). Inside the code a control has its canonical name (`atm_choice`, in `controls_buffer` and `star%ctrl`). Its new-format name is usually the canonical one. `defaults/controls_registry.tsv` maps all three and documents each control. Derived code works on the canonical names; only the reader sees the legacy ones.

**Generated files** (`test_generated_files.py` checks they are current; each generator takes `--check`):

| File | Generator | From |
|---|---|---|
| `state/controls_state_def.inc`, `state/job_controls_def.inc`, `state/controls_sync.f90` | `tools/gen_controls_state.py` | `controls_buffer.f90` |
| `io/inlist_new_decl.inc`, `io/inlist_new_read.inc`, `defaults/star_job.defaults`, `defaults/controls.defaults` | `tools/gen_inlist.py` | the registry, and the reader's defaults for the `.defaults` docs |
| `io/test/controls_dump.inc` | `tools/gen_controls_dump.py` | the `star%ctrl` and `star%job` declarations |
| `io/default_columns.inc` | `tools/gen_default_columns.py` | `defaults/*_columns.list` |

`tools/upgrade_inlist.py` converts a legacy pair to a single-file inlist.

## 3. Output

Each job writes to its output directory (the run log's directory):
- the history file: one row per converged model;
- `profile<N>.data` every `profile_interval` models, with a pulsation file in `pulse_format` beside it;
- the run log;
- the `.mod` model every model. It is the restart file, and the solver's recovery source when a step fails.

The interior-only GYRE file is written every `pulse_gyre_interval` models. The writers only read: core's `compute_observables` and `stitched_model` compute every quantity first.

## 4. Tests

- `test_controls` (`make tests` in `src/`, then `pytest test_controls.py`) reads every deck in `examples/` and `testsuite/` through `read_controls` and compares all of `star%ctrl` and `star%job` with `io/test/expected_controls/`. Any change to the reader shows up here. Regenerate the expected files with `YREC_REGEN_CONTROLS=1` after a deliberate change.
- `test_generated_files.py`: the generated files are current.
- `test_inlist_convert.py`: a converted legacy deck runs identically.
- `test_reentry.py`: in-process reruns.
- `test_mesa_output.py`: the history and profile contract, and GSM output with `YREC_TEST_HDF5=1` on an HDF5 build.

## 5. Extending

**Add a control:**
1. Declare it in `controls_buffer.f90`, with its default as the initializer.
2. Add its row to `defaults/controls_registry.tsv`.
3. Run `tools/gen_controls_state.py`, `tools/gen_inlist.py` and `tools/gen_controls_dump.py`.
4. Read it from `star%ctrl` (or `star%job`) where it is used.
5. Regenerate the expected controls (`YREC_REGEN_CONTROLS=1 pytest test_controls.py`) and check that only the new control changed.

A control that only the new format reads needs nothing in `legacy_controls`.

**Add an output column:** add its name and value to the table in `history_output` or `profile_output`, and list it in `defaults/*_columns.list` if it belongs in the default selection.
