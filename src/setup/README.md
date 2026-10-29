# setup/ — the domains' settings records, built from the controls

The physics domains do not read the star (`docs/GUIDELINES.md`, star rule 1): each takes its configuration as one settings record, which the star layer builds from the run's controls and installs at the start of the run or of each run card. This directory holds those builders: functions from the star, read-only, to one domain's record. It also holds `map_user_inputs`, which is part of the controls read rather than a builder (§1).

Every procedure here takes the star as its argument `s`, `intent(in)` for the builders; the callers pass their own `s` (`run_yrec`'s, or the test programs' `star`). Nothing here names a module-level star (there is none since 2026-10).

---

## 1. What lives where

| Builder (file) | Builds (the type, where it is declared) | Installed by, where | When |
|---|---|---|---|
| `opacity_settings` (`microphysics_settings.f90`) | kap's `kap_settings_t` (`kap/opacity_tables.f90`, through `kap_lib`) | `kap_init` (`kap_lib`), called from `set_up_constants_and_tables` in `core/run/run_setup.f90` | per run |
| `equation_of_state_settings` (`microphysics_settings.f90`) | eos's `eos_settings_t` (`eos/eos_configuration.f90`) | `eos_init` (`eos_lib`), from `set_up_constants_and_tables` | per run |
| `atmosphere_settings` (`microphysics_settings.f90`) | atm's `atm_settings_t` (`atm/atm_tables.f90`) | `atm_init` (`atm_lib`), from `set_up_constants_and_tables` | per run |
| `convection_settings_from_star` (`microphysics_settings.f90`) | convection's `convection_settings_t` (`convection/convection_configuration.f90`) | `set_convection_settings` (`convection_lib`), from `begin_kind_card` in `core/run/run_driver.f90` | per card |
| `diffusion_settings_from_star` (`microphysics_settings.f90`) | element diffusion's `diffusion_settings_t` (`element_diffusion/element_diffusion_configuration.f90`) | `set_diffusion_settings` (`element_diffusion_lib`), from `begin_kind_card` | per card |
| `rotation_settings_from_star` (`rotation_settings_setup.f90`) | rotation's `rotation_settings_t` (`rotation/rotation_configuration.f90`), from the controls and the solar units. For a rotating run it also derives what the controls only imply: the wind law's constant, exponents and saturation (`magnetic_braking`'s pure `derive_wind_law`), the convection zone's rotation law clamped to [-2, 0], and `solid_body_mode_flag = 1` under `force_solid_body_rotation` | `set_rotation_settings` (`rotation_configuration`), from `begin_kind_card` | per card |
| `net_settings_from_star` (`net_settings_setup.f90`) | net's `net_settings_t` (`net/net_configuration.f90`, through `net_lib`), from the controls and the star's cross-section scales | `set_net_settings` (`net_lib`; it also checks the reaction table and returns `ierr`), from `run_yrec` in `core/run/run_driver.f90` | per Monte-Carlo run |

Besides the controls (`s%ctrl`, `s%job`), the builders read a few star members: `use_two_z_tables` (kap), `mixing_length_alpha` (convection: the card's value), `zenvm` (convection: the first card's surface Z without C, N and O; still 0 when the first card's record is built, so it is passed again after that card's model read, below), and the solar units `solar_luminosity_cgs`, `log10_solar_luminosity`, `solar_radius_cgs`, `solar_mass_cgs` (convection, rotation). Each builder's header lists exactly what it reads (`! star reads:`, `! star writes: none`).

| Not a builder (file) | What it does |
|---|---|
| `map_user_inputs` (`user_inputs.f90`) | Inside the controls read: with `use_new_nuclear_rates`, turns the per-reaction S0 values and their derivatives (`s0_*`, `s0p_*`, `s0pp_*`) into the star's scale arrays by reaction slot (`s%cross_section_scale`, `s%qs0e_scale`, `s%qqs0ee_scale`), the extra scales by name (`s%extra_cross_section_scale(xs_pep)`, ...), Be7 + p's derivative scales (`s%be7p_qs0e_scale`, `s%be7p_qqs0ee_scale`) and the two N15 + p branch scales (the derivative and branch scales are reset to 1 first); with `lnewvars`, writes the `tol_*`, `flag_*` and `time_*` controls into `chi_grid_scale`, `atime`, `use_structure_dt_limits`, `max_domega_global` and `diffusion_timestep_factor` in the controls buffer (without it, `chi_grid_scale` is the deck's `hpttol`, inlist name `zoning_tolerances`, which `adopt_canonical_names` copies). It may use the buffer (`io/controls_buffer.f90`) because it runs inside the read (`CONTROLS_BUFFER_IMPORTERS` in `tools/check_boundaries.py`). |

Configuration set elsewhere:
- net's record is installed once per Monte-Carlo run rather than per card, after `apply_monte_carlo_parameters` has scaled that run's rate factors.
- Convection's `envelope_metals_no_cno` is known only after the first card's starting model: `update_surface_mixture` (`core/model/starting_model.f90`) computes `s%zenvm` and passes it with `set_convection_envelope_metals`, after `begin_kind_card` has installed the record (with the value `s%zenvm` had then: 0 on the first card).
- core keeps records of its own for its star-blind kernels: the structure equations' (`set_structure_settings`, `core/structure/structure_equations.f90`) and the envelope integration's (`envelope_integration`'s `cfg`, `core/envelope/envelope_integration.f90`), both filled from the star at every call. They are core's, not a domain's, so they are not built here.

## 2. When each runs

In a run (`run_yrec`, `core/run/run_driver.f90`), in order:
1. `read_controls` (`io/controls_reader.f90`): the namelists into the buffer, `adopt_canonical_names`, `map_user_inputs`, `check_controls` (`io/controls_check.f90`), `derive_controls` (`io/derived_controls.f90`), then the store into `s%ctrl` and `s%job`.
2. `run_yrec` sets `s%use_two_z_tables`; `prepare_run` (`core/run/run_setup.f90`) calls `set_up_constants_and_tables`, which installs the three per-run records: `kap_init(opacity_settings(s), ...)`, `eos_init(equation_of_state_settings(s), ...)`, `atm_init(atmosphere_settings(s), ...)`, each of which also loads its tables.
3. For each Monte-Carlo run: `apply_monte_carlo_parameters`, then `set_net_settings(net_settings_from_star(s), ierr)`.
4. For each run card, `begin_kind_card`: the card's `mixing_length_alpha`, then `set_convection_settings(convection_settings_from_star(s))`, `set_rotation_settings(rotation_settings_from_star(s))` and `set_diffusion_settings(diffusion_settings_from_star(s))`, before the card's starting model is read; on the first card the model read's `update_surface_mixture` adds `set_convection_envelope_metals`.

How the three kinds of setup code differ:
- `io/derived_controls.f90` (`derive_controls`) runs inside the controls read, after `map_user_inputs` and `controls_check`'s `check_controls`, before the final store. It derives controls from controls, on the buffer's names: the light-element rate scales, the choice of one OPAL EOS and one low-temperature opacity source (and `s%job%opal_eos_table_path`), a fresh envelope triangle for rotating runs, the starting-model mixture changes, the switches that diffusion options imply, `s%tenv`, and the run log's header. Its results become part of `s%ctrl` and `s%job`. It changes the controls; it builds no settings record. (Derivations that only a domain needs belong to that domain's builder here: the wind law is `rotation_settings_from_star`'s.) `map_user_inputs` belongs with it.
- `setup/` turns the finished controls into per-domain settings records, with no side effects.
- `core/run/run_setup.f90` acts once per run: the constants (`physical_constants`, the derived solar units in the star, `s%atm_hras`), the tables, with the records from `setup/` (`kap_init`, `eos_init`, `atm_init`), and the Monte-Carlo samples (`prepare_run`). It has side effects.

## 3. The rules

- **Pure copies.** A builder copies controls and the listed star members into the record, one field each, under the domain's readable names (`laol_extrapolation_tolerance` for `tollaol`, `shape_tolerance` for `acfpft`); it computes nothing itself. A derived quantity is computed by `derive_controls` (from controls, once per read), by a pure function of the domain that the builder calls (`rotation_settings_from_star` calls `magnetic_braking`'s `derive_wind_law`), or by the domain from its record.
- **`intent(in)`.** A builder takes `type(star_info), intent(in) :: s` and returns the record; its contract says `star writes: none`, and `tools/check_star_contracts.py` checks the lists (`test_boundaries.py`).
- **No defaults invented.** A builder never supplies a value of its own: every field it sets comes from a control or a star member. A field it does not set keeps the declared initializer of the domain's type, which is then the domain's default and is documented there.
- **One reader.** The domain reads only its record. `tools/check_boundaries.py` enforces it: the directories in `STAR_FREE_DIRS` never use the star at all, rotation (`CONTROLS_FREE_DIRS`) keeps the star for its state but reads no `s%ctrl` or `s%job` (the settling latch it shares with core is the star's `settling_active`), and outside users reach a domain in `FACADE_ONLY_DOMAINS` only through its `<domain>_lib`.
- **Two install points.** Per-run records are installed in `set_up_constants_and_tables` (`core/run/run_setup.f90`), per-card records in `begin_kind_card` (`core/run/run_driver.f90`); the network's per-Monte-Carlo-run record and convection's `envelope_metals_no_cno` are the two exceptions (§1).
- **Names.** Two naming schemes are in use: `opacity_settings`, `equation_of_state_settings`, `atmosphere_settings` (per-run) and `<domain>_settings_from_star` (per-card). New builders use `<domain>_settings_from_star`.
- Builders print nothing and cannot fail; a domain's `*_init` returns `ierr` for a table it cannot load.

## 4. Tests

The builders are exercised by every run; the tests build records the same way the run does:
- `core/test/test_structure.f90` (`pytest test_structure.py`): `convection_settings_from_star`, `rotation_settings_from_star` and `net_settings_from_star`, after `set_up_constants_and_tables` (which installs the kap, eos and atm records);
- `rotation/test/test_rotation.f90` (`test_rotation.py`): `rotation_settings_from_star`, then sets `rotation_settings` fields directly for the braking checks;
- `element_diffusion/test/test_element_diffusion.f90` (`test_element_diffusion.py`): `diffusion_settings_from_star`;
- `net/test/test_net.f90` (`test_net.py`): `map_user_inputs` and `net_settings_from_star`, against the legacy rate scales;
- `eos/test`, `kap/test`, `convection/test`: `set_up_constants_and_tables` on a test star, i.e. the per-run builders;
- `test_boundaries.py`: the builders' star contracts, and the boundary rules of §3.

`python3 src/tools/coverage.py --work DIR --dir setup --test test_structure --test test_net --test test_rotation --test test_element_diffusion --deck examples/run_config_matrix/cm_base.inlist --deck examples/run_config_matrix/cm_rot_base.inlist --deck examples/run_config_matrix/cm_settle_he.inlist` measures the line coverage of setup/: 95.2% of 252 lines (2026-10-29). The builders are covered completely; the 12 lines not run are `map_user_inputs`'s branches for a reaction whose input S0 is zero (its derivative scales are then set to 0).

## 5. How to extend

**Add a field to a domain's settings record:**
1. Declare the field in the domain's settings type, with its unit and meaning and a default initializer: `kap/opacity_tables.f90` (`kap_settings_t`), `eos/eos_configuration.f90`, `atm/atm_tables.f90`, `convection/convection_configuration.f90`, `element_diffusion/element_diffusion_configuration.f90`, `rotation/rotation_configuration.f90`, or `net/net_configuration.f90` (`net_settings_t`). If the domain re-exports the type through its `<domain>_lib`, nothing changes there.
2. Copy it in the builder: `settings%<field> = s%ctrl%<control>` (or `s%job%...`, or a star member) in `setup/microphysics_settings.f90`, `setup/rotation_settings_setup.f90` or `setup/net_settings_setup.f90`.
3. Run `python3 tools/check_star_contracts.py --fix` in `src/`, which adds the control to the builder's `star reads:` list.
4. Read it in the domain as `<record>%<field>` (`kap_settings%...`, `rotation_settings%...`).
5. If the control is new, add it first ("Add a control", below).
6. Give it a check: a test that sets the field directly in the record and sees its effect (as `test_rotation` sets `rotation_settings`).
7. `make -j8`, the unit tests and the decks. A new field that only copies an existing control changes no result.

**Add a physics domain with its own settings record:**
1. A directory `<dom>/` with its facade `<dom>_lib.f90` (the public entry points; it re-exports the settings type and `set_<dom>_settings`) and `<dom>_configuration.f90`: the type `<dom>_settings_t`, the module-level record `type(<dom>_settings_t), save :: <dom>_settings` and the setter `set_<dom>_settings`, as `convection/` and `element_diffusion/` have them. The domain uses no `star_info_lib` (its slot constants come from `state/composition_layout` and the other layout modules), and its state comes in as arguments.
2. Its builder in `setup/`: `function <dom>_settings_from_star(s) result(settings)` with `type(star_info), intent(in) :: s`, a new file `setup/<dom>_settings_setup.f90` (as `rotation_settings_setup.f90`) or an addition to `microphysics_settings.f90`; its header ends with `!   star reads: none` and `!   star writes: none`, and `check_star_contracts.py --fix` fills the list in.
3. Install it: per run in `set_up_constants_and_tables` (`core/run/run_setup.f90`), next to `kap_init`; or per card in `begin_kind_card` (`core/run/run_driver.f90`), next to `set_convection_settings`. Add the `use <dom>_lib, only: set_<dom>_settings` and the builder's `use` there.
4. `tools/check_boundaries.py`: the directory in `DOMAINS`, its entry points in `PUBLIC`, and the directory in `STAR_FREE_DIRS` and `FACADE_ONLY_DOMAINS` (or `CONTROLS_FREE_DIRS`, if it keeps the star for its state, as rotation does).
5. The build: the Makefile finds the sources by itself; add the directory to `DOMAIN_DIRS` and its library to `YREC_LIBS_ORDER` (`src/Makefile`), then `python3 tools/gen_makefile_deps.py`. Its test program: a `test_<dom>` target like `test_convection`'s, added to `tests`, its binary in `.gitignore`, and a `test_<dom>.py` at the root that runs it (`test_convection.py` is the pattern). A unit test of a domain that reads the run-time constants of `state/physical_constants` (`gas_constant`, `cgl`, `cdelrl`, ...) must call `set_up_constants_and_tables` on its star first, as `test_convection` does.
6. Its controls, by "Add a control" below, and a `README.md` for the directory (`rotation/README.md` is the standard).

**Add a control:** the full recipe is `io/README.md` §7, "Add a control" (traced end to end on 2026-10-29). In short:
1. Declare it in `io/controls_buffer.f90` with its default as the initializer, the only place its default is written.
2. Its row in `defaults/controls_registry.tsv`, with `canonical` its name, `home` `ctrl` (`job` only if code outside the read will assign it) and `used` `run`.
3. From `src/`: `python3 tools/gen_controls_state.py`, `python3 tools/gen_inlist.py`, `python3 tools/gen_controls_dump.py`.
4. The field in the domain's record and its copy in the builder here: "Add a field to a domain's settings record", above, from step 1.
5. Values the run cannot use are rejected in `io/controls_check.f90`'s `check_controls`.
6. `make`, `python3 tools/check_boundaries.py` (it checks the row's `home` and `used` against the code), and `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py` from the repository root: only the new member's lines change.
