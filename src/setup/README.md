# setup/ — the domains' settings records: their builders, and where they are installed

The physics domains do not read the star (`docs/GUIDELINES.md`, star rule 1): each takes its configuration as one settings record, which the star layer builds from the run's controls and installs at the start of the run, of each Monte-Carlo run or of each run card; the few values that follow the star's surface mixture are installed whenever it changes.

This directory holds the builders, pure functions from the star, read-only, to one domain's record (`<x>_settings_from_star`, in module `<x>_settings_builder`), and the one module that installs them, `settings_install`. (`map_user_inputs`, which maps the readable controls onto the code's vectors during the controls read, is the read's: `io/user_inputs.f90`.)

Every procedure here takes the star as its argument `s`, `intent(in)` (the one exception, `install_surface_x_settings`, takes the surface X instead); the callers pass their own `s` (`run_yrec`'s, or the test programs' `star`). Nothing here names a module-level star (there is none since 2026-10).

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `kap_settings_builder` | `kap_settings_from_star` and `kap_table_paths_from_star`: kap's record and the paths of the tables `kap_init` loads. |
| `eos_settings_builder` | `eos_settings_from_star` and `eos_table_paths_from_star`. |
| `atm_settings_builder` | `atm_settings_from_star` and `atm_table_paths_from_star`. |
| `convection_settings_builder` | `convection_settings_from_star`. |
| `diffusion_settings_builder` | `diffusion_settings_from_star` (element diffusion). |
| `rotation_settings_builder` | `rotation_settings_from_star`. |
| `net_settings_builder` | `net_settings_from_star`. |
| `settings_install` | Every install of a domain's record, and of the domains' state that follows the surface mixture, one procedure per lifetime (below): `install_run_settings`, `install_monte_carlo_settings`, `install_card_settings`, and per surface mixture `install_surface_mixture_settings`, `install_surface_z_settings`, `install_surface_x_settings`. `tools/check_boundaries.py` (`INSTALL_FILES`) fails on an installer the star layer calls anywhere else. |

| Builder | Builds (the type, where it is declared) | Installed by (`settings_install`), with | When |
|---|---|---|---|
| `kap_settings_from_star`, `kap_table_paths_from_star` | kap's `kap_settings_t` and its tables' paths `kap_table_paths_t` (`kap/opacity_tables.f90`, through `kap_lib`), from the controls; `use_two_z_tables` derived from two of them | `install_run_settings`, with `kap_init` (`kap_lib`), which also loads the tables | per run |
| `eos_settings_from_star`, `eos_table_paths_from_star` | eos's `eos_settings_t` and its tables' paths `eos_table_paths_t` (`eos/eos_configuration.f90`), from the controls | `install_run_settings`, with `eos_init` (`eos_lib`) | per run |
| `atm_settings_from_star`, `atm_table_paths_from_star` | atm's `atm_settings_t` and its tables' paths `atm_table_paths_t` (`atm/atm_tables.f90`), from the controls (`atm_choice` is the deck's, `physics%atm_choice`: it selects the tables `atm_init` loads) | `install_run_settings`, with `atm_init` (`atm_lib`) | per run |
| `net_settings_from_star` | net's `net_settings_t` (`net/net_configuration.f90`, through `net_lib`), from the controls and the star's rate scales | `install_monte_carlo_settings`, with `set_net_settings` (`net_lib`; it also checks the reaction table and returns `ierr`) | per Monte-Carlo run, after `begin_monte_carlo_run` has applied that run's rate factors |
| `convection_settings_from_star` | convection's `convection_settings_t` (`convection/convection_configuration.f90`), with the card's mixing length | `install_card_settings`, with `set_convection_settings` (`convection_lib`) | per card |
| `rotation_settings_from_star` | rotation's `rotation_settings_t` (`rotation/rotation_configuration.f90`), from the controls, the solar units and log10 G. For a rotating run it also derives what the controls only imply: the wind law's constant, exponents and saturation (`magnetic_braking`'s pure `derive_wind_law`), the convection zone's rotation law clamped to [-2, 0], and `solid_body_mode_flag = 1` under `force_solid_body_rotation` | `install_card_settings`, with `set_rotation_settings` (`rotation_configuration`) | per card |
| `diffusion_settings_from_star` | element diffusion's `diffusion_settings_t` (`element_diffusion/element_diffusion_configuration.f90`), with the Monte-Carlo run's `settling_z_factor_in_use` | `install_card_settings`, with `set_diffusion_settings` (`element_diffusion_lib`) | per card |

The state that follows the surface mixture is not a record but three setters, which `settings_install` calls from the star's surface mixture (the lifetime "per surface mixture"):

| Install (`settings_install`) | Sets | When |
|---|---|---|
| `install_surface_mixture_settings` | convection's `envelope_metals_no_cno` (`set_convection_envelope_metals`; `install_card_settings` installed 0, the value before the first card's model read), then the EOS mixture (`install_surface_z_settings`) | the first card's model read (`update_surface_mixture`, `core/model/starting_model.f90`) |
| `install_surface_z_settings` | the EOS mixture (`eos_set_mixture`: the envelope X, Z, 1/mu_ion and number fractions; with SCV on the EOS builds its envelope table for them) | the run plan's envelope Z before `kap_init` (`set_run_plan_envelope_z`, `core/run/run_setup.f90`), a Z rescale (`rescale_interior_z`, `core/model/model_rescaling.f90`), and with the whole mixture |
| `install_surface_x_settings` | kap's surface opacity tables for a surface X (`kap_update_surface_tables`; returns `ierr`) | every card's model read (`update_surface_mixture`, the envelope X), a rezone that has moved the surface X (`refresh_surface_opacity`, `core/model/rezoning.f90`) |

Besides the controls (`s%control`, `s%physics`), the builders read these star members (each builder's header lists exactly what it reads: `! star reads:`, `! star writes: none`):
- the run copy `settling_z_factor_in_use` (element diffusion: the Monte-Carlo run's metal settling factor);
- `mixing_length_alpha` (convection: the card's value) and `envelope_metals_no_cno` (convection: the first card's surface Z without C, N and O; still 0 when the first card's record is built, so `install_surface_mixture_settings` passes it again after that card's model read);
- the solar units: `solar_luminosity_cgs` and `log10_solar_luminosity` (convection, rotation), `solar_radius_cgs` and `solar_mass_cgs` (rotation);
- net's rate scales, which `map_user_inputs` sets during the read and each Monte-Carlo run partly rescales: `cross_section_scale`, `extra_cross_section_scale`, `s_prime_scale`, `s_double_prime_scale`, `be7p_s_prime_scale`, `be7p_s_double_prime_scale`, `n15_pg_scale`, `n15_pa_scale`.

| Not a builder (file) | What it does |
|---|---|
| `map_user_inputs` (`io/user_inputs.f90`, part of the controls read since 2026-10-29) | Inside the controls read: with `use_new_nuclear_rates`, turns the per-reaction S0 values and their derivatives (`s0_*`, `s0p_*`, `s0pp_*`) into the star's scale arrays by reaction slot (`s%cross_section_scale`, `s%s_prime_scale`, `s%s_double_prime_scale`), the extra scales by name (`s%extra_cross_section_scale(xs_pep)`, ...), Be7 + p's derivative scales (`s%be7p_s_prime_scale`, `s%be7p_s_double_prime_scale`) and the two N15 + p branch scales (the derivative and branch scales are reset to 1 first); with `use_new_rate_variables` (legacy `lnewvars`), writes the `tol_*`, `flag_*` and `time_*` controls into `chi_grid_scale`, `timestep_factors`, `use_structure_dt_limits`, `max_domega_global` and `diffusion_timestep_factor` in the controls buffer (without it, `chi_grid_scale` is the deck's `hpttol`, inlist name `zoning_tolerances`, which `adopt_canonical_names` copies). It may use the buffer (`io/controls_buffer.f90`) because it runs inside the read (`CONTROLS_BUFFER_IMPORTERS` in `tools/check_boundaries.py`). The slot names of `chi_grid_scale` and `timestep_factors` (`ichi_*`, `itime_*`) are `state/control_slots.f90`'s. |

core keeps records of its own for its star-blind kernels: the structure equations' (`set_structure_settings`, `core/structure/structure_equations.f90`) and the envelope integration's (`envelope_integration`'s `cfg`, `core/envelope/envelope_integration.f90`), both filled from the star at every call. They are core's, not a domain's, so they are neither built nor installed here.

## 2. When each runs

In a run (`run_yrec`, `core/run/run_driver.f90`), in order:
1. `read_controls` (`io/controls_reader.f90`): the namelists into the buffer, `adopt_canonical_names`, `map_user_inputs`, `check_controls` (`io/controls_check.f90`), `derive_controls` (`io/derived_controls.f90`), then the store into `s%control` and `s%physics`, and the run copies armed from them.
2. `prepare_run` (`core/run/run_setup.f90`) sets the run plan's envelope Z and the EOS mixture (`set_run_plan_envelope_z`, `install_surface_z_settings`), then calls `set_up_constants_and_tables`, which sets the constants and then calls `install_run_settings`: `kap_init(kap_settings_from_star(s), kap_table_paths_from_star(s), ...)`, `eos_init(eos_settings_from_star(s), eos_table_paths_from_star(s), ierr)`, `atm_init(atm_settings_from_star(s), atm_table_paths_from_star(s), ierr)`, each of which also loads its tables.
3. For each Monte-Carlo run: `begin_monte_carlo_run` (the run copies re-armed, the run's samples applied), then `install_monte_carlo_settings`: `set_net_settings(net_settings_from_star(s), ierr)`.
4. For each run card, `begin_kind_card`: the card's `mixing_length_alpha`, then `install_card_settings`: `set_convection_settings(convection_settings_from_star(s))`, `set_rotation_settings(rotation_settings_from_star(s))` and `set_diffusion_settings(diffusion_settings_from_star(s))`, before the card's starting model is read.
5. Per surface mixture (the table above): the first card's model read calls `install_surface_mixture_settings`, every card's model read `install_surface_x_settings`, a Z rescale `install_surface_z_settings`, a rezone that has moved the surface X `install_surface_x_settings`.

How the three kinds of setup code differ:
- `io/derived_controls.f90` (`derive_controls`) runs inside the controls read, after `map_user_inputs` and `controls_check`'s `check_controls`, before the final store. It derives controls from controls, on the buffer's names: the light-element rate scales, the choice of one OPAL EOS and one low-temperature opacity source (and `s%physics%opal_eos_table_path`), a fresh envelope triangle for rotating runs, the starting-model mixture changes, the switches that diffusion options imply, `s%envelope_store_log_temperature`, and the run log's header. Its results become part of `s%control` and `s%physics`. It changes the controls; it builds no settings record. `map_user_inputs` belongs with it.
- `setup/` turns the finished controls into per-domain settings records (the builders, pure) and installs them (`settings_install`, whose `install_run_settings` also loads the tables).
- `core/run/run_setup.f90` acts once per run: the run plan's envelope Z and EOS mixture, the constants (`physical_constants`' `cgl`, the derived solar units in the star), then `install_run_settings`, and the Monte-Carlo samples (`prepare_run`).

## 3. The rules

- **Copies, and the few derivations listed.** A builder copies controls and the listed star members into the record, one field each, under the domain's readable names (`shape_tolerance` for `fpft_average_weight`, legacy `acfpft`; `log10_t_cut` for `gradient_logt_cutoff`). Three builders derive as well, and their headers say so: `rotation_settings_from_star` derives the wind law (`magnetic_braking`'s pure `derive_wind_law`) and clamps the convection zone's rotation law, `kap_settings_from_star` derives `use_two_z_tables` (`control%use_z_ramp .or. physics%use_diffusion_z`), and convection's and rotation's builders take log10 G as `log10_newton_g(G_cgs)`. The OPAL vintage is a derived control (`opal_eos_vintage`, set by `derive_controls`) that `eos_settings_from_star` copies. Any other derived quantity is computed by `derive_controls` (from controls, once per read), by a pure function of the domain that the builder calls, or by the domain from its record.
- **Pure, and `intent(in)`.** A builder is a `pure function` that takes `type(star_info), intent(in) :: s` and returns the record, so the compiler rejects a side effect; its contract says `star writes: none`, and `tools/check_star_contracts.py` checks the lists (`test_boundaries.py`).
- **Every field set, no defaults invented.** A builder sets every field of its record (`tools/check_boundaries.py`'s `check_settings_fields` fails on a field it leaves unset, unless `SETTINGS_DEFAULT_ONLY` lists it with the reason; none is listed), and never supplies a value of its own: every field comes from a control or a star member, whose default is written once (a control's in `io/controls_buffer.f90`). The initializers of the domain's type are what a test that fills a record by hand starts from. Nothing checks yet that a field copies the right control (§6).
- **One reader.** The domain reads only its record. `tools/check_boundaries.py` enforces it: the directories in `STAR_FREE_DIRS` never use the star at all, rotation (`CONTROLS_FREE_DIRS`) keeps the star for its state but reads neither `s%control` nor `s%physics` (the settling latch it shares with core is the star's `settling_active`), no domain and no star-blind kernel reads `physical_constants`' run-time `cgl`, and outside users reach a domain in `FACADE_ONLY_DOMAINS` only through its `<domain>_lib`.
- **One install module.** Every installer call (`set_<domain>_settings`, `kap_init`, `eos_init`, `atm_init`, and the setters of the state that follows the surface mixture, `set_convection_envelope_metals`, `eos_set_mixture`, `kap_update_surface_tables`) is in `settings_install.f90`, in the procedure of its lifetime; `run_setup`, `run_driver`, `starting_model`, `model_rescaling` and `rezoning` call those procedures. `tools/check_boundaries.py` (`INSTALL_FILES`, `INSTALLER_CALL_RE`) enforces it in the star layer (a domain's calls inside its own directory, such as `eos_init`'s of its own setter, are its own; the tests install records themselves).
- **One name.** Module `<x>_settings_builder` (file `setup/<x>_settings_builder.f90`) holds `<x>_settings_from_star`, and `<x>_table_paths_from_star` for the paths of the tables a domain's init reads (kap, eos, atm); `<x>` is the prefix of the domain's settings type `<x>_settings_t` (`diffusion` for element diffusion). The name differs from the domains' installed records (`kap_settings`, `rotation_settings`, ...), which are the domains' module variables.
- Builders print nothing and cannot fail; the per-run installs return `ierr` for a table they cannot load, `set_net_settings` for a reaction table it rejects, and `install_surface_x_settings` for a surface table.

## 4. Tests

The builders are exercised by every run; the tests build records the same way the run does:
- `core/test/test_structure.f90` (`pytest test_structure.py`): `convection_settings_from_star`, `rotation_settings_from_star` and `net_settings_from_star`, after `set_up_constants_and_tables` (whose `install_run_settings` installs the kap, eos and atm records);
- `rotation/test/test_rotation.f90` (`test_rotation.py`): `rotation_settings_from_star`, then sets `rotation_settings` fields directly for the braking checks;
- `element_diffusion/test/test_element_diffusion.f90` (`test_element_diffusion.py`): `diffusion_settings_from_star`;
- `net/test/test_net.f90` (`test_net.py`): `map_user_inputs` and `net_settings_from_star`, against the legacy rate scales;
- `kap/test`: `set_up_constants_and_tables` on a test star, i.e. the per-run builders and `install_run_settings`; `eos/test` and `convection/test` boot no star, they fill an `eos_settings_t` and an `eos_table_paths_t` themselves and call `eos_init`;
- `test_boundaries.py`: the builders' and the installs' star contracts, and the boundary rules of §3 (the install sites among them).

`python3 src/tools/coverage.py --work DIR --dir setup --test test_net --test test_rotation --test test_element_diffusion --deck examples/run_config_matrix/cm_base.inlist --deck examples/run_config_matrix/cm_rot_base.inlist --deck examples/run_config_matrix/cm_settle_he.inlist` measures the line coverage of setup/: 100% of 226 lines (2026-10-31): every builder and every install runs. (`test_structure` needs a deck argument, which `coverage.py` does not give it, so it is left out of the command; `map_user_inputs`, counted here until 2026-10-29, is io/'s.)

## 5. How to extend

**Add a field to a domain's settings record** (traced on 2026-10-31 on a throwaway branch with a copy of `gradient_logt_cutoff` in convection's record, read in `temperature_gradients`):
1. Declare the field in the domain's settings type, with its unit and meaning and a default initializer: `kap/opacity_tables.f90` (`kap_settings_t`), `eos/eos_configuration.f90`, `atm/atm_tables.f90`, `convection/convection_configuration.f90`, `element_diffusion/element_diffusion_configuration.f90`, `rotation/rotation_configuration.f90`, or `net/net_configuration.f90` (`net_settings_t`). If the domain re-exports the type through its `<domain>_lib`, nothing changes there.
2. Copy it in the builder, `setup/<x>_settings_builder.f90`: `settings%<field> = s%physics%<control>` (or `s%control%...`, its group's record; or a star member such as a run copy). `tools/check_boundaries.py` fails until the builder sets it.
3. Run `python3 tools/check_star_contracts.py --fix` in `src/`, which adds the control to the builder's `star reads:` list; if the field copies a star member rather than a control, also `python3 tools/gen_star_docs.py`, as the builder is then one of that member's readers in `state/README.md`.
4. Read it in the domain as `<record>%<field>` (`convection_settings%...`, `rotation_settings%...`).
5. If the control is new, add it first ("Add a control", below).
6. Give it a check: a test that sets the field directly in the record and sees its effect (as `test_rotation` sets `rotation_settings`).
7. `make -j8`, the unit tests and the decks. A new field that only copies an existing control changes no result.

**Add a physics domain with its own settings record** (traced on 2026-10-31 on a throwaway branch with a domain `dummy/` installed per card: it built, linked, passed `check_boundaries.py`, the contracts and its own test):
1. A directory `<dom>/` with its facade `<dom>_lib.f90` (the public entry points; it re-exports the settings type and `set_<dom>_settings`) and `<dom>_configuration.f90`: the type `<dom>_settings_t`, the module-level record `type(<dom>_settings_t), save :: <dom>_settings` and the setter `set_<dom>_settings`, as `convection/` and `element_diffusion/` have them. The domain uses no `star_info_lib` (its slot constants come from `state/composition_layout` and the other layout modules), and its state comes in as arguments.
2. Its builder: a new file `setup/<dom>_settings_builder.f90`, module `<dom>_settings_builder`, holding `pure function <dom>_settings_from_star(s) result(settings)` with `type(star_info), intent(in) :: s` (`setup/convection_settings_builder.f90` is the pattern); its header ends with `!   star reads: none` and `!   star writes: none`, and `check_star_contracts.py --fix` fills the list in. If it reads star members, run `python3 tools/gen_star_docs.py` too (`state/README.md` lists the builder among their readers).
3. Install it: its call `call set_<dom>_settings(<dom>_settings_from_star(s))` in `setup/settings_install.f90`, in the procedure of its lifetime (`install_run_settings` per run, `install_monte_carlo_settings` per Monte-Carlo run, `install_card_settings` per card), with the `use <dom>_lib, only: set_<dom>_settings` and the builder's `use` in that procedure. Nothing else calls the setter (`INSTALL_FILES` in `tools/check_boundaries.py`). State of the domain that follows the star's surface mixture is not in the record: see "Add per-surface-mixture state", below.
4. `tools/check_boundaries.py`: the directory in `DOMAINS`, its entry points in `PUBLIC`, and the directory in `STAR_FREE_DIRS` and `FACADE_ONLY_DOMAINS` (or `CONTROLS_FREE_DIRS`, if it keeps the star for its state, as rotation does).
5. The build: the Makefile finds the sources by itself; add the directory to `DOMAIN_DIRS` and its library to `YREC_LIBS_ORDER` (`src/Makefile`), then `python3 tools/gen_makefile_deps.py`. Its test program: a `test_<dom>` target like `test_convection`'s, added to `tests` and to `clean`'s list of binaries, its binary in `.gitignore`, and a `test_<dom>.py` at the root that runs it (`test_convection.py` is the pattern); CI (`.github/workflows/ci.yml`) builds and runs the test programs it names. The constants of `state/physical_constants` are compile-time parameters, so a domain's unit test needs no setup for them; log10 G is the one run-time value (`cgl`, set by `run_setup` from the control `G_cgs`), and a domain takes it from its settings record, as convection's and rotation's `log10_gravitational_constant` (`tools/check_boundaries.py` fails on a domain or a star-blind kernel that reads `cgl`).
6. Its controls, by "Add a control" below, and a `README.md` for the directory (`rotation/README.md` is the standard).

**Add per-surface-mixture state** (a value of a domain that follows the star's surface mixture, as convection's `envelope_metals_no_cno`, the EOS mixture and kap's surface tables do; traced on 2026-10-31 on a throwaway branch with a convection value set from the envelope Z):
1. In the domain, the value and its setter: a module variable (or a field of the installed record) and a public `set_<dom>_<what>(...)` that takes plain values, re-exported by `<dom>_lib`, as `convection_lib`'s `set_convection_envelope_metals`.
2. Its call in `setup/settings_install.f90`, in the surface-mixture procedure of the change it follows: `install_surface_mixture_settings` (the first card's whole mixture), `install_surface_z_settings` (the envelope Z: the run plan's, a Z rescale) or `install_surface_x_settings` (the surface X: every card, a rezone), from the star's surface members (`envelope_*`; `state/README.md`, "Run and card bookkeeping"). Nothing outside `settings_install` calls the setter.
3. Add the setter's name to `INSTALLER_CALL_RE` in `tools/check_boundaries.py`, so that a call elsewhere fails, and to the "One install module" rule of §3.
4. `python3 tools/check_star_contracts.py --fix` (the install's `star reads:` list), `python3 tools/check_boundaries.py`, `make -j8`, the unit tests and the decks.

**Add a control:** the full recipe is `io/README.md` §5, "Add a control" (traced end to end on 2026-10-29). In short:
1. Declare it in `io/controls_buffer.f90` with its default as the initializer, the only place its default is written.
2. Its row in `defaults/controls_registry.tsv`, with `canonical` its name, `home` its `group` (`control` or `physics`; nothing outside the read assigns it: a value the run changes gets a run copy, `state/README.md` §7) and `used` `run`.
3. From `src/`: `python3 tools/gen_controls_state.py`, `python3 tools/gen_inlist.py`, `python3 tools/gen_controls_dump.py`.
4. The field in the domain's record and its copy in the builder here: "Add a field to a domain's settings record", above, from step 1.
5. Values the run cannot use are rejected in `io/controls_check.f90`'s `check_controls`.
6. `make`, `python3 tools/check_boundaries.py` (it checks the row's `home` and `used` against the code), and `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py` from the repository root: only the new member's lines change.

## 6. Open

- No check ties a field to the right control: `check_settings_fields` sees that every field is set, not from what. The plan (it needs the controls registry, the controls read's file): a registry column naming each control's settings field (`<x>%<field>`), filled once from the builders' copies; `check_settings_fields` then compares every builder copy `settings%f = s%<group>%c` with the row of `c`, and `gen_inlist.py` lists the field in `controls_registry.md`.
