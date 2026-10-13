# setup/ — the domains' settings records, built from the controls

The physics domains do not read `star` (`docs/GUIDELINES.md`, star rule 1): each takes its configuration as one settings record, which the star layer installs at the start of the run or of each run card. This directory builds those records from the run's controls (`star%ctrl`, `star%job`) and a few other `star` members. It also holds `map_user_inputs`, which, during the controls read, turns the descriptively named rate and tolerance controls into the arrays the code uses.

| Module (file) | What it holds |
|---|---|
| `microphysics_settings` | `opacity_settings` (kap's `kap_settings_t`), `equation_of_state_settings` (eos's `eos_settings_t`), `atmosphere_settings` (atm's `atm_settings_t`) and `convection_settings_from_star` (convection's `convection_settings_t`). Each is a function that copies controls into the record; besides `ctrl` and `job` they read `star%use_two_z_tables` (kap) and `star%mixing_length_alpha`, `star%solar_luminosity_cgs`, `star%log10_solar_luminosity` and `star%zenvm` (convection). |
| `rotation_settings_setup` | `rotation_settings_from_star`: rotation's `rotation_settings_t`, from controls only. |
| `user_inputs` | `map_user_inputs`: with `use_new_nuclear_rates`, the per-reaction S0 values and their derivatives (`s0_*`, `s0p_*`, `s0pp_*`) become `star%cross_section_scale`, `star%qs0e_scale`, `star%qqs0ee_scale` and the two N15 + p branch scales; with `lnewvars`, the `tol_*`, `flag_*` and `time_*` controls overwrite `chi_grid_scale`, `atime`, `use_structure_dt_limits`, `max_domega_global` and `diffusion_timestep_factor` in the controls buffer. The derivative and branch scales are reset to 1 first. |

The settings functions write nothing (their headers say `star writes: none`); `map_user_inputs` writes the five rate-scale members of `star` and the controls buffer (`io/controls_buffer.f90`), which it may import because it runs inside the read (`CONTROLS_BUFFER_IMPORTERS` in `tools/check_boundaries.py`).

## When each runs

| Procedure | Called from | When |
|---|---|---|
| `map_user_inputs` | `io/read_controls.f90` (`read_input`) | Once per controls read, i.e. once per `run_yrec` call: after `adopt_canonical_names`, before the buffer is first stored into `star%ctrl`/`star%job`. |
| `opacity_settings`, `equation_of_state_settings`, `atmosphere_settings` | `core/run/run_setup.f90` (`set_up_constants_and_tables`), as the first argument of `kap_init`, `eos_init` and `atm_init` | Once per run, after the controls read (`prepare_run`), before the first run card. |
| `convection_settings_from_star`, `rotation_settings_from_star` | `core/run/run_driver.f90` (`begin_kind_card`), through `set_convection_settings` and `set_rotation_settings` | At the start of every run card, after the card's `mixing_length_alpha` is set and before its starting model is read. |

The tests build records the same way: `core/test/test_structure.f90` (convection and rotation), `rotation/test/test_rotation.f90` (rotation) and `net/test/test_net.f90` (`map_user_inputs`).

Two pieces of domain configuration are set elsewhere. The network's record, `net_settings`, is built by `configure_network` in `core/step/burn_step.f90` (once per Monte-Carlo run, from the rate scales above). Convection's `envelope_metals_no_cno` is computed later than the record: `update_surface_mixture` (`core/model/starting_model.f90`) computes `star%zenvm` on the first card, after `begin_kind_card` has installed the record, and passes it with `set_convection_envelope_metals`. `element_diffusion/` still reads `star%ctrl` and `star%job` directly.

## How setup/ differs from the other "setup" code

Three places prepare a run:

1. **`io/derived_controls.f90`** (`derive_controls`) runs inside the controls read, after `map_user_inputs` and before the final store. It derives controls from controls, on the buffer's names: the light-element rate scales, the choice of one OPAL EOS and one low-temperature opacity source (and `star%job%opal_eos_table_path`), the wind law's coefficients and exponents, the starting-model mixture changes, the switches that diffusion options imply, `star%tenv`, and the run log's header; it also rejects semiconvection combined with overshoot. Its results become part of `star%ctrl` and `star%job`. It changes the controls; it builds no settings record.
2. **`setup/`** reads the finished controls and hands each domain its record. `map_user_inputs` is the exception: it runs inside the read, like `derive_controls`, and writes controls and rate scales rather than a settings record.
3. **`core/run/run_setup.f90`** acts on the controls once per run: it sets the constants (the variables of `physical_constants`, the derived solar units in `star`, `star%atm_hras`), loads the tables by calling `kap_init`, `eos_init` and `atm_init` with the records from `microphysics_settings`, and sets up the Monte-Carlo runs (`prepare_run`). It has side effects; the settings functions have none.

In short: `io/derived_controls` turns controls into more controls, `setup/` turns controls into per-domain settings records, and `core/run/run_setup` uses those records to initialize the domains.
