# YREC inlist controls registry (draft for review)

Proposed namelist revamp: `&control`/`&physics` -> `&star_job`/`&controls`
with readable names; where a genuine MESA control equivalent exists, the
MESA name is adopted for the inlist.

Status legend: **canonical** = name already established internally during
the modernization (parmin's post-read copy); **proposed** = hand-named in
this draft; **keep** = already readable, unchanged; **todo** = named
provisionally, semantics need verification before freezing.

Defaults marked `*` are implicit static zero/false/'' (no explicit
assignment before the namelist read; pinned by -finit-local-zero).

## &control -> &star_job (146 entries)

| legacy | type | default | proposed name | status | MESA | doc |
|---|---|---|---|---|---|---|
| `cmixla` | real(50) | `0.0*` | **`mixing_length_alpha`** | canonical | mesa: mixing_length_alpha |  |
| `calsolage` | real | `0.0*` | **`target_solar_age`** | canonical |  |  |
| `calsolzx` | real | `0.0*` | **`target_solar_zx`** | canonical |  |  |
| `descrip` | string(2) | `''*` | **`run_description`** | proposed |  | Per-run description strings echoed in output headers |
| `endage` | real(50) | `0.0*` | **`max_age`** | canonical | mesa: max_age | former common/sett/: target_end_age/timestep_override/ central_deuterium_stop/central_hydrogen_stop/central_he |
| `flaol` | string | `''*` | **`laol_opacity_table_file`** | proposed |  | LAOL 1989 opacity table |
| `fpurez` | string | `''*` | **`pure_z_opacity_table_file`** | proposed |  | Pure-Z opacity table (metal diffusion support) |
| `flaol2` | string | `''*` | **`laol_opacity_table2_file`** | proposed |  | LAOL 1989 opacity table, second Z |
| `fopal2` | string | `''*` | **`opal92_opacity_table2_file`** | proposed |  | OPAL 1992 opacity table, second Z |
| `flast` | string | `''*` | **`last_model_file`** | proposed |  | Path for the final model written at end of job |
| `ffirst` | string | `''*` | **`first_model_file`** | proposed |  | Starting model read when read_first_model is set for a run |
| `ffermi` | string | `''*` | **`fermi_table_file`** | proposed |  | Fermi-Dirac integral table (Yale/Prather degenerate-electron EOS) |
| `fdebug` | string | `''*` | **`debug_output_file`** | proposed |  | Debug output stream |
| `ftrack` | string | `''*` | **`track_output_file`** | proposed |  | Evolution track (.track) output path |
| `fshort` | string | `''*` | **`short_output_file`** | proposed |  | Model summary (.short) output path |
| `fstch` | string | `''*` | **`stored_change_file`** | todo |  | TODO verify: stored structure-change stream |
| `fmilne` | string | `''*` | **`milne_output_file`** | todo |  | TODO verify: Milne atmosphere output |
| `fmodpt` | string | `''*` | **`model_print_file`** | proposed |  | Per-model print stream (rotational-mixing composition changes etc.) |
| `fstor` | string | `''*` | **`store_output_file`** | proposed |  | Full stored-model (.store) output path |
| `fpmod` | string | `''*` | **`pulse_model_file`** | proposed |  | Pulsation (oscillation-code) model output |
| `fpatm` | string | `''*` | **`pulse_atm_file`** | proposed |  | Pulsation output: atmosphere part |
| `fpenv` | string | `''*` | **`pulse_env_file`** | proposed |  | Pulsation output: envelope part |
| `fdyn` | string | `''*` | **`dynamo_output_file`** | todo |  | TODO verify: rotation/dynamo diagnostic output |
| `flldat` | string | `''*` | **`ll_data_file`** | proposed |  | Livermore opacity auxiliary data file |
| `fsnu` | string | `''*` | **`neutrino_output_file`** | proposed |  | Solar neutrino flux (.snu) output path |
| `fscomp` | string | `''*` | **`surface_composition_file`** | todo |  | TODO verify: surface composition output |
| `fkur` | string | `''*` | **`kurucz_atm_table_file`** | proposed |  | Kurucz model-atmosphere boundary-condition table |
| `fmhd1` | string | `''*` | **`mhd_table_file_1`** | proposed |  | MHD EOS table 1 of 8 |
| `fmhd2` | string | `''*` | **`mhd_table_file_2`** | proposed |  | MHD EOS table 2 of 8 |
| `fmhd3` | string | `''*` | **`mhd_table_file_3`** | proposed |  | MHD EOS table 3 of 8 |
| `fmhd4` | string | `''*` | **`mhd_table_file_4`** | proposed |  | MHD EOS table 4 of 8 |
| `fmhd5` | string | `''*` | **`mhd_table_file_5`** | proposed |  | MHD EOS table 5 of 8 |
| `fmhd6` | string | `''*` | **`mhd_table_file_6`** | proposed |  | MHD EOS table 6 of 8 |
| `fmhd7` | string | `''*` | **`mhd_table_file_7`** | proposed |  | MHD EOS table 7 of 8 |
| `fmhd8` | string | `''*` | **`mhd_table_file_8`** | proposed |  | MHD EOS table 8 of 8 |
| `fiso` | string | `''*` | **`isotope_output_file`** | todo |  | TODO verify: isotope output stream |
| `fatm` | string | `''*` | **`atm_table_file`** | proposed |  | Atmosphere P(T) boundary-condition table |
| `fkur2` | string | `''*` | **`kurucz_castelli_atm_table_file`** | proposed |  | Castelli-Kurucz model-atmosphere table |
| `fallard` | string | `''*` | **`allard_atm_table_file`** | proposed |  | Allard NextGen model-atmosphere table |
| `fscvh` | string | `''*` | **`scv_h_table_file`** | proposed |  | Saumon-Chabrier-Van Horn EOS: hydrogen table |
| `fscvhe` | string | `''*` | **`scv_he_table_file`** | proposed |  | SCV EOS: helium table |
| `fscvz` | string | `''*` | **`scv_z_table_file`** | proposed |  | SCV EOS: metals table |
| `fopale` | string | `''*` | **`opal_eos_table_file`** | proposed |  | OPAL 1995/2001 EOS table |
| `fliv95` | string | `''*` | **`opal95_opacity_table_file`** | proposed |  | OPAL/Livermore 1995 opacity tables |
| `fmonte1` | string | `''*` | **`monte_carlo_file1_path`** | canonical |  | former common/iomonte/: monte_carlo_file1_path/monte_carlo_file2_path (originally fmonte1/fmonte2) are NAMELIS |
| `fmonte2` | string | `''*` | **`monte_carlo_file2_path`** | canonical |  | former common/iomonte/: monte_carlo_file1_path/monte_carlo_file2_path (originally fmonte1/fmonte2) are NAMELIS |
| `ipver` | integer | `1` | **`pulsation_file_version`** | canonical |  |  |
| `itrver` | integer | `1` | **`track_file_version`** | canonical |  | former common/track/: track_file_version (originally itrver) is a NAMELIST value with a different canonical sp |
| `kindrn` | integer(50) | `0*` | **`run_kind`** | todo |  | Per-run kind card: what run NK does (rescale-and-evolve, evolve, ...). TODO enumerate codes |
| `ldebug` | logical | `.false.*` | **`write_debug_output`** | proposed |  | Enable the debug output stream |
| `lcorr` | logical | `.false.*` | **`write_correction_output`** | todo |  | TODO verify: print solver correction diagnostics |
| `lmilne` | logical | `.false.*` | **`write_milne_output`** | todo |  | TODO verify |
| `ltrack` | logical | `.false.*` | **`write_track_output`** | proposed | ~ MESA history | Write the .track evolution file |
| `lstore` | logical | `.false.*` | **`write_store_output`** | proposed | ~ MESA profiles | Write full stored models (.store) |
| `lfirst` | logical(50) | `.false.*` | **`first_call_flag`** | canonical |  |  |
| `lstpch` | logical | `.false.*` | **`store_final_model_flag`** | proposed |  | Store the last model of each run (with write_store_output) |
| `lscrib` | logical | `.false.*` | **`scribe_flag`** | todo |  | TODO verify |
| `lstch` | logical | `.false.*` | **`store_changes_flag`** | todo |  | TODO verify |
| `lstatm` | logical | `.false.*` | **`store_atm_profile`** | proposed |  | Include the atmosphere layers in stored models |
| `lstenv` | logical | `.false.*` | **`store_env_profile`** | proposed |  | Include the envelope in stored models |
| `lstmod` | logical | `.false.*` | **`store_interior_profile`** | proposed |  | Include the interior in stored models |
| `lstphys` | logical | `.false.*` | **`store_physics_columns`** | proposed |  | Include physics (EOS/opacity) columns in stored models |
| `lstrot` | logical | `.false.*` | **`store_rotation_columns`** | proposed |  | Include rotation columns in stored models |
| `lpulse` | logical | `.false.` | **`pulsation_output_active`** | canonical |  | former common/pulse/: pulsation_output_active/pulsation_file_version (originally lpulse/ipver) are NAMELIST va |
| `lzramp` | logical | `.false.` | **`use_z_ramp`** | canonical |  |  |
| `lteff` | logical | `.false.*` | **`specify_teff_flag`** | canonical |  |  |
| `lcalst` | logical | `.false.*` | **`calibrate_star_flag`** | canonical |  |  |
| `lpurez` | logical | `.false.` | **`use_pure_z_table`** | canonical |  |  |
| `liso` | logical | `.false.*` | **`isochrone_output_active`** | canonical |  | former common/chrone/: lrwsh_placeholder/isochrone_output_active (originally lrwsh/liso) are NAMELIST values w |
| `lrwsh` | logical | `.false.*` | **`lrwsh_placeholder`** | canonical |  | former common/chrone/: lrwsh_placeholder/isochrone_output_active (originally lrwsh/liso) are NAMELIST values w |
| `lsenv0a` | logical(50) | `.false.*` | **`has_senv0_array`** | canonical |  |  |
| `lpout` | logical | `.false.*` | **`po_output_enabled`** | canonical |  |  |
| `lcals` | logical | `.false.*` | **`calibrate_solar_model`** | canonical |  |  |
| `lcalsolzx` | logical | `.false.*` | **`calibrate_solar_zx`** | canonical |  |  |
| `llaol89` | logical | `.false.*` | **`use_laol89_tables`** | canonical |  |  |
| `lopal92` | logical | `.false.*` | **`use_opal92_tables`** | canonical |  |  |
| `lopal95` | logical | `.false.*` | **`use_opal95_tables`** | canonical |  |  |
| `lkur90` | logical | `.false.*` | **`use_kurucz90_tables`** | canonical |  |  |
| `lalex95` | logical | `.false.*` | **`use_alex95_tables`** | canonical |  |  |
| `npoint` | integer | `0*` | **`profile_point_interval`** | todo |  | TODO verify: print every Nth mass point |
| `npenv` | integer | `0*` | **`envelope_print_points`** | todo |  | TODO verify: envelope points printed |
| `nprtmod` | integer | `0*` | **`profile_model_interval`** | proposed |  | Write a stored model every N models |
| `nprtpt` | integer | `1` | **`print_point_interval`** | canonical |  |  |
| `numrun` | integer | `0*` | **`num_runs`** | canonical |  |  |
| `nmodls` | integer(50) | `0*` | **`max_model_number`** | canonical | mesa: max_model_number |  |
| `opecalex` | string(7) | `'OPACALEXANDER.X00',` | **`alex_z_value`** | todo |  | TODO verify: Z of the Alexander opacity table in use |
| `poa` | real | `0.0*` | **`po_weight_l`** | canonical |  | former common/po/: po_weight_l/po_weight_teff/po_weight_age/ po_max_len_sq/po_output_enabled (originally poa/p |
| `pob` | real | `0.0*` | **`po_weight_teff`** | canonical |  | former common/po/: po_weight_l/po_weight_teff/po_weight_age/ po_max_len_sq/po_output_enabled (originally poa/p |
| `poc` | real | `0.0*` | **`po_weight_age`** | canonical |  | former common/po/: po_weight_l/po_weight_teff/po_weight_age/ po_max_len_sq/po_output_enabled (originally poa/p |
| `pomax` | real | `0.0*` | **`po_max_len_sq`** | canonical |  |  |
| `rsclm` | real(50) | `0.0*` | **`rescale_mass`** | proposed | mesa: new_mass | Per-run target mass (Msun) for rescale runs |
| `rsclx` | real(50) | `0.0*` | **`rescale_x`** | proposed |  | Per-run target hydrogen fraction for rescale runs |
| `rsclz` | real(50) | `0.0*` | **`rescale_z`** | proposed | mesa: new_Z | Per-run target metallicity for rescale runs |
| `rsclcm` | real(50) | `0.0*` | **`rescale_core_mass`** | proposed |  | Per-run rescale: core mass bound |
| `rsclzc` | real(50) | `0.0*` | **`rescale_core_z`** | todo |  | TODO verify: core Z rescale |
| `rsclzm1` | real(50) | `0.0*` | **`rescale_z_mix1`** | todo |  | TODO verify |
| `rsclzm2` | real(50) | `0.0*` | **`rescale_z_mix2`** | todo |  | TODO verify |
| `setdt` | real(50) | `0.0*` | **`timestep_override`** | canonical |  |  |
| `senv0a` | real(50) | `0.0*` | **`senv0_array`** | canonical |  |  |
| `steff` | real | `0.0*` | **`target_teff`** | canonical |  |  |
| `sr` | real | `0.0*` | **`target_radius_rsun`** | canonical |  |  |
| `tolr` | real | `0.0*` | **`radius_tolerance`** | canonical |  | former common/cals2/: luminosity_tolerance/radius_tolerance/ zx_tolerance/calibrate_solar_model/calibrate_sola |
| `toll` | real | `0.0*` | **`luminosity_tolerance`** | canonical |  | former common/cals2/: luminosity_tolerance/radius_tolerance/ zx_tolerance/calibrate_solar_model/calibrate_sola |
| `tolz` | real | `0.0*` | **`zx_tolerance`** | canonical |  |  |
| `xenv0a` | real(50) | `0.0*` | **`initial_x_array`** | canonical |  | former common/newxym/: initial_x_array/initial_z_array/ mixing_length_array/has_senv0_array/senv0_array (origi |
| `xls` | real | `0.0*` | **`target_luminosity_lsun`** | canonical |  | former common/calstar/: target_luminosity_lsun/ target_star_luminosity_tolerance/target_teff/target_radius_rsu |
| `xlstol` | real | `0.0*` | **`target_star_luminosity_tolerance`** | canonical |  |  |
| `zenv0a` | real(50) | `0.0*` | **`initial_z_array`** | canonical |  | former common/newxym/: initial_x_array/initial_z_array/ mixing_length_array/has_senv0_array/senv0_array (origi |
| `zlaol1` | real | `0.0*` | **`laol_table_z1`** | canonical |  | former common/newopac/: laol_table_z1/laol_table_z2/opal_table_z1/ opal_table_z2/opal95_single_table_z/alex_ta |
| `zlaol2` | real | `0.0*` | **`laol_table_z2`** | canonical |  | former common/newopac/: laol_table_z1/laol_table_z2/opal_table_z1/ opal_table_z2/opal95_single_table_z/alex_ta |
| `zopal1` | real | `0.0*` | **`opal_table_z1`** | canonical |  | former common/newopac/: laol_table_z1/laol_table_z2/opal_table_z1/ opal_table_z2/opal95_single_table_z/alex_ta |
| `zopal2` | real | `0.0*` | **`opal_table_z2`** | canonical |  |  |
| `zopal951` | real | `0.0*` | **`opal95_single_table_z`** | canonical |  |  |
| `zopal952` | real | `0.0*` | **`opal95_z2`** | proposed |  | Second Z of the OPAL95 opacity table pair |
| `zalex1` | real | `0.0*` | **`alex_table_z1`** | canonical |  |  |
| `zalex2` | real | `0.0*` | **`alex_z2`** | proposed |  | Second Z of the Alexander table pair |
| `zkur1` | real | `0.0*` | **`kurucz_table_z1`** | canonical |  |  |
| `zkur2` | real | `0.0*` | **`kurucz_table_z2`** | canonical |  |  |
| `fopale01` | string | `''*` | **`opal01_eos_table_file`** | proposed |  | OPAL 2001 EOS table |
| `fcondopacp` | string | `''*` | **`conductive_opacity_table_file`** | proposed |  | Potekhin conductive opacity table |
| `fopale06` | string | `''*` | **`opal06_eos_table_file`** | proposed |  | OPAL 2006 EOS table |
| `falex06` | string | `''*` | **`alex06_opacity_table_file`** | proposed |  | Ferguson/Alexander 2006 low-T opacity tables |
| `lalex06` | logical | `.false.*` | **`use_alex06_tables`** | canonical |  |  |
| `end_dcen` | real(50) | `0.0*` | **`central_deuterium_stop`** | canonical | mesa: log_center_density_limit |  |
| `end_xcen` | real(50) | `0.0*` | **`central_hydrogen_stop`** | canonical | mesa: xa_central_lower_limit |  |
| `end_ycen` | real(50) | `0.0*` | **`central_helium_stop`** | canonical | mesa: xa_central_lower_limit (He) |  |
| `isetmix` | integer | `0*` | **`mixture_change_mode`** | canonical |  | former common/newmx/: mixture_change_mode/isotope_change_mode/ target_carbon_cno_fraction/target_nitrogen_cno_ |
| `isetiso` | integer | `0*` | **`isotope_change_mode`** | canonical |  | former common/newmx/: mixture_change_mode/isotope_change_mode/ target_carbon_cno_fraction/target_nitrogen_cno_ |
| `amix` | string | `''*` | **`cno_mixture_name`** | proposed |  | Heavy-element mixture selector (GS98/CUS/...) applied to the start model |
| `aiso` | string | `''*` | **`isotope_mixture_name`** | proposed |  | Isotope-ratio set selector (L21 default or CUS custom) |
| `frac_c` | real | `0.0*` | **`target_carbon_cno_fraction`** | canonical |  |  |
| `frac_n` | real | `0.0*` | **`target_nitrogen_cno_fraction`** | canonical |  |  |
| `frac_o` | real | `0.0*` | **`target_oxygen_cno_fraction`** | canonical |  |  |
| `r12_13` | real | `0.0*` | **`c12_to_c13_ratio`** | canonical |  |  |
| `r14_15` | real | `0.0*` | **`n14_to_n15_ratio`** | canonical |  |  |
| `r16_17` | real | `0.0*` | **`o16_to_o17_ratio`** | canonical |  |  |
| `r16_18` | real | `0.0*` | **`o16_to_o18_ratio`** | canonical |  |  |
| `zxmix` | real | `0.0*` | **`target_metal_fraction`** | canonical |  |  |
| `xh2_ini` | real | `0.0*` | **`initial_h2_fraction`** | canonical |  |  |
| `xhe3_ini` | real | `0.0*` | **`initial_he3_fraction`** | canonical |  |  |
| `xli6_ini` | real | `0.0*` | **`initial_li6_fraction`** | canonical |  |  |
| `xli7_ini` | real | `0.0*` | **`initial_li7_fraction`** | canonical |  |  |
| `xbe9_ini` | real | `0.0*` | **`initial_be9_fraction`** | canonical |  |  |
| `xb10_ini` | real | `0.0*` | **`initial_b10_fraction`** | canonical |  |  |
| `xb11_ini` | real | `0.0*` | **`initial_b11_fraction`** | canonical |  |  |
| `pulse_gyre_interval` | integer | `0*` | **`pulse_gyre_interval`** | keep |  | former common/pulsegyre/: pulse_gyre_interval is a NAMELIST /control/ value spelled identically to its canonic |

## &physics -> &controls (245 entries)

| legacy | type | default | proposed name | status | MESA | doc |
|---|---|---|---|---|---|---|
| `atmmin` | real | `0.0*` | **`atm_step_min`** | canonical |  |  |
| `atmbeg` | real | `0.0*` | **`atm_step_begin`** | canonical |  | former common/intatm/: all 5 members are NAMELIST /physics/ values with different canonical spellings than cor |
| `atmerr` | real | `0.0*` | **`atm_error_tol`** | canonical |  | former common/intatm/: all 5 members are NAMELIST /physics/ values with different canonical spellings than cor |
| `atmmax` | real | `0.0*` | **`atm_step_max`** | canonical |  |  |
| `atmd0` | real | `0.0*` | **`atm_step_initial`** | canonical |  | former common/intatm/: all 5 members are NAMELIST /physics/ values with different canonical spellings than cor |
| `anewcp` | string | `''*` | **`new_cp_factor`** | todo |  | TODO verify (pairs with xnewcp/lnewcp) |
| `atmp` | string | `''*` | **`atm_print_flag`** | todo |  | TODO verify |
| `acfpft` | real | `0.0*` | **`fpft_average_weight`** | todo |  | TODO verify: averaging coefficient in fp/ft update |
| `atime` | real(14) | `0.0*` | **`age_time_factor`** | todo |  | TODO verify |
| `alphac` | real | `0.0*` | **`overshoot_alpha_core`** | proposed | ~ MESA overshoot f | Core overshoot extent in pressure scale heights |
| `alphae` | real | `0.0*` | **`overshoot_alpha_envelope`** | proposed |  | Envelope overshoot extent in pressure scale heights |
| `alfa` | real | `0.0*` | **`grid_alpha`** | todo |  | TODO verify (hpoint grid weighting) |
| `alpham` | real | `0.0*` | **`overshoot_alpha_min`** | todo |  | TODO verify |
| `atmstp` | real | `0.0*` | **`atm_step_size`** | canonical |  | former common/envgen/: all 3 members are NAMELIST values with different canonical spellings, kept local in cor |
| `abstol` | real | `0.0*` | **`absolute_tolerance`** | canonical |  | former common/burtol/: min_abundance/absolute_tolerance/ relative_tolerance/max_burn_iterations (originally cm |
| `betac` | real | `0.0*` | **`overshoot_beta_core`** | todo |  | TODO verify |
| `cmin` | real | `0.0*` | **`min_abundance`** | canonical |  | former common/burtol/: min_abundance/absolute_tolerance/ relative_tolerance/max_burn_iterations (originally cm |
| `clsun` | real | `0.0*` | **`solar_luminosity_cgs`** | canonical |  | former common/const/: solar physical constants (luminosity, mass, radius in cgs, plus their log10/ln/bolometri |
| `crsun` | real | `0.0*` | **`solar_radius_cgs`** | canonical |  |  |
| `dpenv` | real | `0.0*` | **`envelope_logP_step`** | todo |  | TODO verify: envelope integration pressure step |
| `dtdif` | real | `0.0*` | **`diffusion_timestep_factor`** | proposed |  | Sub-timestep factor for rotation diffusion |
| `dtwind` | real | `1.0d1` | **`max_domega_global`** | canonical |  | former common/ct2/: max_domega_global (originally dtwind) is a NAMELIST value with a different canonical spell |
| `djok` | real | `1.0d-4` | **`convergence_tolerance`** | canonical |  |  |
| `dt_gs` | real | `0.0*` | **`settling_timestep_fraction`** | canonical |  | former common/gravs2/: all 4 members are NAMELIST /physics/ values, each with a different canonical spelling t |
| `enverr` | real | `0.0*` | **`env_error_tol`** | canonical |  | former common/intenv/: all 4 members are NAMELIST /physics/ values with different canonical spellings than cor |
| `envmax` | real | `0.0*` | **`env_step_max`** | canonical |  |  |
| `envmin` | real | `0.0*` | **`env_step_min`** | canonical |  | former common/intenv/: all 4 members are NAMELIST /physics/ values with different canonical spellings than cor |
| `envbeg` | real | `0.0*` | **`env_step_begin`** | canonical |  | former common/intenv/: all 4 members are NAMELIST /physics/ values with different canonical spellings than cor |
| `envstp` | real | `0.0*` | **`envelope_step_size`** | canonical |  | former common/envgen/: all 3 members are NAMELIST values with different canonical spellings, kept local in cor |
| `etadh0` | real | `0.0*` | **`debye_huckel_eta_min`** | canonical |  |  |
| `etadh1` | real | `0.0*` | **`debye_huckel_eta_max`** | canonical |  |  |
| `fcorr0` | real | `0.0*` | **`solver_correction_factor0`** | todo |  | TODO verify: under-relaxation factor, first iterations |
| `fcorri` | real | `0.0*` | **`solver_correction_factor`** | todo |  | TODO verify: under-relaxation factor, later iterations |
| `fk` | real | `0.0*` | **`dynamical_shear_factor`** | todo |  | TODO verify: instability strength factor |
| `fw` | real | `0.0*` | **`difad_velocity_scale`** | canonical |  |  |
| `fc` | real | `0.0*` | **`mixing_velocity_scale`** | canonical |  |  |
| `fo` | real | `0.0*` | **`horizontal_shear_factor`** | todo |  | TODO verify |
| `fmu` | real | `0.0*` | **`mu_gradient_scale`** | canonical |  |  |
| `fes` | real | `0.0*` | **`es_velocity_scale`** | canonical |  |  |
| `fcore` | real | `0.0*` | **`core_mass_reduction_factor`** | canonical |  |  |
| `fgsf` | real | `0.0*` | **`gsf_velocity_scale`** | canonical |  |  |
| `fss` | real | `0.0*` | **`secular_shear_velocity_scale`** | canonical |  |  |
| `fesc` | real | `0.0*` | **`es_mixing_scale`** | canonical |  |  |
| `fssc` | real | `0.0*` | **`secular_shear_mixing_scale`** | canonical |  |  |
| `fgsfc` | real | `0.0*` | **`gsf_mixing_scale`** | canonical |  |  |
| `fgry` | real | `0.0*` | **`settling_y_factor`** | todo |  | TODO verify: gravitational settling factor (Y) |
| `fgrz` | real | `0.0*` | **`settling_z_factor`** | todo |  | TODO verify: gravitational settling factor (Z) |
| `grtol` | real | `0.0*` | **`settling_tolerance`** | canonical |  | former common/gravst/: all 4 members are NAMELIST /physics/ values, each with a different canonical spelling t |
| `htoler` | real(5,2) | `0.0*` | **`solver_tolerance`** | todo |  | TODO verify: Henyey convergence tolerance |
| `hpttol` | real(12) | `1.0d-8,8.0d-2,5.0d-2,5.0d-2,1.0d0,1.0d0,0.0d0,5.0d-2, &` | **`rezone_mass_spacing`** | proposed |  | Desired mass spacing for rezoning (hpoint) |
| `itfp1` | integer | `0*` | **`fpft_iter1`** | todo |  | TODO verify: fp/ft iteration bound |
| `itfp2` | integer | `0*` | **`fpft_iter2`** | todo |  | TODO verify |
| `imax` | integer | `0*` | **`max_stage_index`** | canonical |  |  |
| `itdif1` | integer | `0*` | **`num_rotation_structure_iters`** | proposed |  | Structure<->rotation iterations per step (itrot loop) |
| `itdif2` | integer | `0*` | **`max_diffusion_iters`** | proposed |  | Rotation diffusion-solve iteration limit (seculr/checkj) |
| `ies` | integer | `0*` | **`enable_eddington_sweet`** | proposed |  | Eddington-Sweet circulation switch (integer enable) |
| `igsf` | integer | `1` | **`gsf_inhibition_mode`** | canonical |  |  |
| `imu` | integer | `0*` | **`enable_mu_gradient`** | proposed |  | Mean-molecular-weight gradient effects switch |
| `ilambda` | integer | `0*` | **`coulomb_log_choice`** | canonical |  |  |
| `kttau` | integer | `0` | **`atm_choice`** | canonical | mesa: atm_option (int code here) | former common/atmos/: atm_choice (originally kttau) is a NAMELIST value with a different canonical spelling, k |
| `kemmax` | integer | `0*` | **`max_burn_iterations`** | canonical |  |  |
| `lvfc` | logical | `.false.*` | **`velocity_check_flag`** | todo |  | TODO verify |
| `ldifad` | logical | `.false.` | **`use_diffusion_advection_transport`** | canonical |  |  |
| `lnoj` | logical | `.false.` | **`no_am_transport_in_core`** | canonical |  | former common/notran/: no_am_transport_in_core (originally lnoj) is a NAMELIST value with a different canonica |
| `lnewdif` | logical | `.false.` | **`use_new_diffusion_routines`** | canonical |  | former common/gravs4/: use_new_diffusion_routines (originally lnewdif) is a NAMELIST value with a different ca |
| `ldify` | logical | `.false.*` | **`diffuse_helium_active`** | canonical |  |  |
| `ldifz` | logical | `.false.` | **`use_diffusion_z`** | canonical |  |  |
| `ldifli` | logical | `.false.*` | **`diffuse_lithium`** | proposed |  | Include Li/Be in gravitational settling |
| `lsnu` | logical | `.false.*` | **`calc_neutrinos`** | proposed |  | Compute solar neutrino fluxes |
| `ldh` | logical | `.false.*` | **`use_debye_huckel_correction`** | canonical |  |  |
| `lnewcp` | logical | `.false.*` | **`rescale_species_active`** | canonical |  |  |
| `lkuthe` | logical | `.false.` | **`helium_flash_active`** | canonical |  | former common/heflsh/: helium_flash_active (originally lkuthe) is a NAMELIST value with a different canonical  |
| `lovstc` | logical | `.false.*` | **`core_overshoot_active`** | proposed |  | Enable convective-core overshoot |
| `lovste` | logical | `.false.` | **`envelope_overshoot_active`** | canonical |  |  |
| `lovstm` | logical | `.false.*` | **`min_overshoot_active`** | todo |  | TODO verify: third overshoot switch (with lovstc/lovste) |
| `lovmax` | logical | `.false.*` | **`overshoot_max_extent_flag`** | todo |  | TODO verify |
| `lexcom` | logical | `.false.` | **`use_extended_composition`** | canonical |  | former common/flag/: single NAMELIST /physics/ member selecting the extended (15-species) vs. default composit |
| `lrot` | logical | `.false.*` | **`rotation_flag`** | canonical | mesa: rotation_flag |  |
| `lnew0` | logical | `.false.*` | **`new_triangle_init_flag`** | todo |  | TODO verify: envelope-fit triangle initialization control |
| `linstb` | logical | `.false.*` | **`instability_transport_active`** | canonical |  |  |
| `lwnew` | logical | `.false.*` | **`set_initial_omega`** | proposed |  | Initialize solid-body rotation at omega = initial_omega |
| `ljdot0` | logical | `.false.*` | **`use_wind_torque`** | proposed |  | Apply magnetic wind angular-momentum loss |
| `lptime` | logical | `.true.` | **`use_structure_dt_limits`** | canonical |  | former common/ct3/: use_structure_dt_limits (originally lptime) is a NAMELIST value with a different canonical |
| `ladov` | logical | `.false.*` | **`adiabatic_overshoot_active`** | todo |  | TODO verify: adiabatic extension of overshoot region |
| `ltrist` | logical | `.false.` | **`use_envelope_triangle_dt`** | canonical |  | former common/govs/: use_envelope_triangle_dt (originally ltrist) is a NAMELIST value with a different canonic |
| `lenvg` | logical | `.false.*` | **`envelope_generation_flag`** | canonical |  |  |
| `lnulos1` | logical | `.false.` | **`use_itoh_neutrino_loss`** | canonical |  | former common/nuloss/'s one config member: switch selecting the Itoh 1996 neutrino-loss routines (net_lib.f90' |
| `lthoul` | logical | `.false.*` | **`use_thoul_diffusion`** | proposed |  | Thoul et al. diffusion coefficients |
| `lthoulfit` | logical | `.false.*` | **`use_thoul_fit`** | canonical |  |  |
| `lopale` | logical | `.false.*` | **`use_opal95_eos`** | canonical |  | former common/opaleos/: use_opal95_eos/use_opal2001_eos/ use_opal2006_eos/use_numerical_derivatives (originall |
| `lmhd` | logical | `.false.` | **`use_mhd_eos`** | canonical |  | former common/mhd/: use_mhd_eos is a NAMELIST /physics/ value (core/parmin.f90's lmhd, kept local there and co |
| `lcore` | logical | `.false.*` | **`extend_core_inward`** | canonical |  | former common/core/: extend_core_inward/num_core_shells_added/ core_mass_reduction_factor (originally lcore/mc |
| `lsemic` | logical | `.false.*` | **`use_semiconvection`** | proposed |  | Semiconvection treatment (conflicts with core overshoot) |
| `lnews` | logical | `.false.*` | **`improved_first_guess_flag`** | proposed |  | Use stored structure changes for the next model's first guess |
| `mcore` | integer | `0*` | **`num_core_shells_added`** | canonical |  |  |
| `niter1` | integer | `0*` | **`max_iter_level1`** | proposed |  | Henyey iteration limit, level 1 (fresh surface BC) |
| `niter2` | integer | `0*` | **`max_iter_level2`** | proposed |  | Henyey iteration limit, level 2 |
| `niter3` | integer | `0*` | **`max_iter_level3`** | proposed |  | Henyey iteration limit, level 3 |
| `niter4` | integer | `0*` | **`max_iter_level4`** | proposed |  | Henyey iteration limit, level 4 (high precision) |
| `nuse` | integer | `0*` | **`extrap_order`** | canonical |  |  |
| `niter_gs` | integer | `0*` | **`settling_num_iterations`** | canonical |  |  |
| `optol` | real | `1.0d-8` | **`metal_fraction_match_tolerance`** | canonical |  | former common/optab/: metal_fraction_match_tolerance is a NAMELIST /physics/ value (core/parmin.f90's optol, k |
| `rcrit` | real | `0.0*` | **`critical_reynolds`** | canonical |  |  |
| `reltol` | real | `0.0*` | **`relative_tolerance`** | canonical |  |  |
| `stolr0` | real | `0.0*` | **`tolerance_fraction`** | canonical |  |  |
| `tcut` | real(5) | `0.0*` | **`nuclear_logT_cutoffs`** | proposed |  | log T cutoffs (5) below which burning groups are switched off |
| `tscut` | real | `0.0*` | **`saha_log10t_cutoff`** | canonical |  |  |
| `tenv0` | real | `0.0*` | **`envelope_logT0`** | todo |  | TODO verify: envelope temperature bound |
| `tenv1` | real | `0.0*` | **`envelope_logT1`** | todo |  | TODO verify |
| `tgcut` | real | `0.0*` | **`grav_energy_logT_cutoff`** | todo |  | TODO verify: T cutoff for gravitational-energy handling |
| `tridt` | real | `0.0*` | **`tri_delta_teffl`** | canonical |  | former common/cenv/: lnew0 is a NAMELIST /physics/ value spelled identically to its canonical name. tri_delta_ |
| `tridl` | real | `0.0*` | **`tri_delta_logl`** | canonical |  | former common/cenv/: lnew0 is a NAMELIST /physics/ value spelled identically to its canonical name. tri_delta_ |
| `tollaol` | ? | `?` | **`laol_z_tolerance`** | todo |  | TODO verify: LAOL table Z-match tolerance |
| `vnew` | real(12) | `0.0*` | **`mixture_weights_seed`** | proposed |  | Species mixture-weight vector seeding starin's mixture algorithm |
| `walpcz` | real | `0.0*` | **`cz_alpha_weight`** | todo |  | TODO verify |
| `wnew` | real | `0.0*` | **`initial_omega`** | proposed |  | Initial rotation rate when set_initial_omega |
| `weakscreening` | real | `0.03d0` | **`weak_screening_threshold`** | canonical |  |  |
| `xnewcp` | real | `0.0*` | **`new_species_value`** | canonical |  | former common/newcmp/: new_species_value/rescale_species_active (originally xnewcp/lnewcp) are NAMELIST values |
| `xmin` | real | `0.0*` | **`hydrogen_diffusion_floor`** | canonical |  |  |
| `ymin` | real | `0.0*` | **`helium_diffusion_min`** | canonical |  |  |
| `tmolmin` | real | `0.0*` | **`molecular_opacity_logt_min`** | canonical |  |  |
| `tmolmax` | real | `0.0*` | **`molecular_opacity_logt_max`** | canonical |  |  |
| `lmonte` | logical | `.false.*` | **`monte_carlo_active`** | proposed |  | Monte-Carlo multi-realization mode |
| `imbeg` | integer | `0*` | **`mc_run_start`** | proposed |  | First Monte-Carlo realization index |
| `imend` | integer | `0*` | **`mc_run_end`** | proposed |  | Last Monte-Carlo realization index |
| `sstandard` | real(17) | `0.9828,1.0485,0.9815,0.9241,1.3818,1.0542,1.0, &` | **`mc_standard_seed`** | todo |  | TODO verify |
| `lscv` | logical | `.false.` | **`use_scv_eos`** | canonical |  |  |
| `ldisk` | logical | `.false.*` | **`disk_locking_active`** | canonical |  |  |
| `tdisk` | real | `0.0*` | **`disk_temperature`** | canonical |  | former common/disk/: disk_temperature/disk_pressure/ disk_locking_active are NAMELIST /physics/ values (core/p |
| `pdisk` | real | `0.0*` | **`disk_pressure`** | canonical |  |  |
| `wmax` | real | `0.0*` | **`wind_saturation_omega`** | canonical |  |  |
| `lsolid` | logical | `.false.*` | **`force_solid_body_rotation`** | canonical |  | former common/sbrot/: force_solid_body_rotation/solid_body_mode_flag (originally lsolid/impjmod) are NAMELIST  |
| `impjmod` | integer | `0*` | **`solid_body_mode_flag`** | canonical |  |  |
| `dmdt0` | real | `0.0*` | **`mass_change`** | canonical | mesa: mass_change | former common/masschg/: mass-accretion/Reimers-wind parameters, all NAMELIST /physics/ values. use_mass_accret |
| `fczdmdt` | real | `0.0*` | **`accretion_cz_fraction`** | todo |  | TODO verify: accreted-mass fraction into the CZ |
| `ftotdmdt` | real | `0.0*` | **`accretion_total_fraction`** | todo |  | TODO verify |
| `compacc` | real(15) | `0.0*` | **`accreted_composition`** | canonical |  |  |
| `creim` | real | `0.0*` | **`reimers_scaling_factor`** | proposed | mesa: Reimers_scaling_factor | Reimers mass-loss efficiency eta |
| `lreimer` | logical | `.false.*` | **`use_reimers_wind`** | proposed |  | Enable Reimers red-giant mass loss |
| `lmdot` | logical | `.false.*` | **`use_mass_change`** | canonical | mesa: mass_change on/off |  |
| `lopale01` | logical | `.false.*` | **`use_opal2001_eos`** | canonical |  | former common/opaleos/: use_opal95_eos/use_opal2001_eos/ use_opal2006_eos/use_numerical_derivatives (originall |
| `lcondopacp` | logical | `.false.*` | **`use_conductive_opacity`** | canonical |  |  |
| `lopale06` | logical | `.false.*` | **`use_opal2006_eos`** | canonical |  | former common/opaleos/: use_opal95_eos/use_opal2001_eos/ use_opal2006_eos/use_numerical_derivatives (originall |
| `lnumderiv` | logical | `.false.*` | **`use_numerical_derivatives`** | canonical |  |  |
| `alatm_feh` | real | `0d0` | **`allard_target_feh`** | canonical |  | former common/alatm03/: allard_target_feh/allard_target_alpha/ allard_use_tau100 (originally alatm_feh/alatm_a |
| `alatm_alpha` | real | `0d0` | **`allard_target_alpha`** | canonical |  | former common/alatm03/: allard_target_feh/allard_target_alpha/ allard_use_tau100 (originally alatm_feh/alatm_a |
| `laltptau100` | logical | `.false.` | **`allard_use_tau100`** | canonical |  |  |
| `cstmixing` | real | `0.0*` | **`constant_mixing_coeff`** | proposed |  | Constant added diffusion coefficient (mixing) |
| `cstdiffmix` | real | `0.0*` | **`constant_settling_reduction`** | proposed |  | Constant coefficient mimicking mixing to reduce settling |
| `lsolwind` | logical | `.false.` | **`calibrate_solar_wind`** | todo |  | TODO verify: solar wind-loss calibration mode |
| `lmwind` | logical | `.false.*` | **`use_pmm_wind_law`** | canonical |  |  |
| `lrossby` | logical | `.false.*` | **`scale_by_rossby_number`** | canonical |  |  |
| `lpmm` | logical | `.false.*` | **`use_pmm_wind`** | proposed |  | Use the PMM wind law (Somers pinned Matt-style) |
| `lbscale` | logical | `.false.*` | **`scale_by_b_field`** | canonical |  |  |
| `awind` | string | `''*` | **`wind_law_name`** | canonical |  |  |
| `pmma` | real | `0.0*` | **`pmm_exponent_a`** | canonical |  | former common/pmmwind/: pmm_exponent_a/pmm_exponent_b/ pmm_exponent_c/pmm_exponent_d/pmm_exponent_m/pmm_norm_j |
| `pmmb` | real | `0.0*` | **`pmm_exponent_b`** | canonical |  | former common/pmmwind/: pmm_exponent_a/pmm_exponent_b/ pmm_exponent_c/pmm_exponent_d/pmm_exponent_m/pmm_norm_j |
| `pmmc` | real | `0.0*` | **`pmm_exponent_c`** | canonical |  | former common/pmmwind/: pmm_exponent_a/pmm_exponent_b/ pmm_exponent_c/pmm_exponent_d/pmm_exponent_m/pmm_norm_j |
| `pmmd` | real | `0.0*` | **`pmm_exponent_d`** | canonical |  |  |
| `pmmm` | real | `0.0*` | **`pmm_exponent_m`** | canonical |  |  |
| `pmmjd` | real | `0.0*` | **`pmm_norm_jdot`** | canonical |  |  |
| `pmmmd` | real | `0.0*` | **`pmm_norm_mdot`** | canonical |  |  |
| `pmmwmax` | real | `0.0*` | **`pmm_omega_sat`** | proposed |  | PMM wind: saturation omega |
| `pmmsolp` | real | `0.0*` | **`pmm_solar_pressure`** | canonical |  |  |
| `pmmsolw` | real | `0.0*` | **`pmm_solar_omega`** | canonical |  |  |
| `pmmsoltau` | real | `0.0*` | **`pmm_solar_turnover_timescale`** | canonical |  |  |
| `lcodm` | logical | `.false.*` | **`use_constant_background_diffusion`** | canonical |  |  |
| `codm` | real | `0.0*` | **`constant_background_diffusion_coeff`** | canonical |  | former common/mag/: constant_background_diffusion_coeff/ use_constant_background_diffusion (originally codm/lc |
| `wmax_sun` | real | `0.0*` | **`solar_omega_sat`** | proposed |  | Saturation omega for the solar-calibrated wind |
| `xsli6` | real | `0.0*` | **`li6_rate_scale`** | todo |  | TODO verify: Li6 burning rate scale |
| `xsli7` | real | `0.0*` | **`li7_rate_scale`** | todo |  | TODO verify |
| `xsbe91` | real | `0.0*` | **`be9_rate_scale1`** | todo |  | TODO verify |
| `xsbe92` | real | `0.0*` | **`be9_rate_scale2`** | todo |  | TODO verify |
| `xsbe93` | real | `0.0*` | **`be9_rate_scale3`** | todo |  | TODO verify |
| `lxli6` | logical | `.false.*` | **`scale_li6_rate`** | todo |  | TODO verify |
| `lxli7` | logical | `.false.*` | **`scale_li7_rate`** | todo |  | TODO verify |
| `lxbe91` | logical | `.false.*` | **`scale_be9_rate1`** | todo |  | TODO verify |
| `lxbe92` | logical | `.false.*` | **`scale_be9_rate2`** | todo |  | TODO verify |
| `lxbe93` | logical | `.false.*` | **`scale_be9_rate3`** | todo |  | TODO verify |
| `sli6` | real | `0.0*` | **`li6_rate_scale`** | canonical |  | former common/burnscs/: light-element cross-section scale factors. core/parmin.f90's own local names (sli6 etc |
| `sli7` | real | `0.0*` | **`li7_rate_scale`** | canonical |  | former common/burnscs/: light-element cross-section scale factors. core/parmin.f90's own local names (sli6 etc |
| `sbe91` | real | `0.0*` | **`be9_pg_rate_scale`** | canonical |  |  |
| `sbe92` | real | `0.0*` | **`be9_pd_rate_scale`** | canonical |  |  |
| `sbe93` | real | `0.0*` | **`be9_palpha_rate_scale`** | canonical |  |  |
| `lnewnuc` | logical | `.false.` | **`use_new_nuclear_rates`** | canonical |  |  |
| `spotf` | real | `0.0*` | **`spot_filling_factor`** | canonical |  | former common/spots/: all 3 members are NAMELIST /physics/ values, same namelist-can't-rename treatment as abo |
| `spotx` | real | `0.0*` | **`spot_temp_contrast`** | canonical |  |  |
| `lsdepth` | logical | `.false.*` | **`spot_depth_varies`** | canonical |  |  |
| `s0_1_1` | real | `0.0*` | **`s0_pp`** | canonical |  | former common/newcross/: user-supplied nuclear reaction S-factors (and first/second derivative ratios relative |
| `s0_3_3` | real | `0.0*` | **`s0_he3he3`** | canonical |  | former common/newcross/: user-supplied nuclear reaction S-factors (and first/second derivative ratios relative |
| `s0_3_4` | real | `0.0*` | **`s0_he3he4`** | canonical |  | former common/newcross/: user-supplied nuclear reaction S-factors (and first/second derivative ratios relative |
| `s0_1_12` | real | `0.0*` | **`s0_p_c12`** | canonical |  | former common/newcross/: user-supplied nuclear reaction S-factors (and first/second derivative ratios relative |
| `s0_1_13` | real | `0.0*` | **`s0_p_c13`** | canonical |  | former common/newcross/: user-supplied nuclear reaction S-factors (and first/second derivative ratios relative |
| `s0_1_14` | real | `0.0*` | **`s0_p_n14`** | canonical |  |  |
| `s0_1_16` | real | `0.0*` | **`s0_p_o16`** | canonical |  |  |
| `s0_pep` | real | `0.0*` | **`s0_pep`** | keep |  |  |
| `s0_1_be7e` | real | `0.0*` | **`s0_be7_electron`** | canonical |  |  |
| `s0_1_be7p` | real | `0.0*` | **`s0_be7_p`** | canonical |  |  |
| `s0_hep` | real | `0.0*` | **`s0_hep`** | keep |  |  |
| `s0_1_15_c12alp` | real | `0.0*` | **`s0_n15_p_c12_branch`** | canonical |  |  |
| `s0_1_15_o16` | real | `0.0*` | **`s0_n15_p_o16_branch`** | canonical |  |  |
| `s0p_1_1` | real | `0.0*` | **`s0p_pp`** | canonical |  |  |
| `s0p_3_3` | real | `0.0*` | **`s0p_he3he3`** | canonical |  |  |
| `s0p_3_4` | real | `0.0*` | **`s0p_he3he4`** | canonical |  |  |
| `s0p_1_12` | real | `0.0*` | **`s0p_p_c12`** | canonical |  |  |
| `s0p_1_13` | real | `0.0*` | **`s0p_p_c13`** | canonical |  |  |
| `s0p_1_14` | real | `0.0*` | **`s0p_p_n14`** | canonical |  |  |
| `s0p_1_16` | real | `0.0*` | **`s0p_p_o16`** | canonical |  |  |
| `s0pp_1_12` | real | `0.0*` | **`s0pp_p_c12`** | canonical |  |  |
| `s0pp_1_13` | real | `0.0*` | **`s0pp_p_c13`** | canonical |  |  |
| `s0pp_1_16` | real | `0.0*` | **`s0pp_p_o16`** | canonical |  |  |
| `s0p_1_be7p` | real | `0.0*` | **`s0p_be7_p`** | canonical |  |  |
| `s0pp_1_be7p` | real | `0.0*` | **`s0pp_be7_p`** | canonical |  |  |
| `flag_dx` | real | `0.0*` | **`flag_dx`** | keep |  | former common/newparam/: all 30 NAMELIST values (the "intuitively named" timestep/tolerance parameters replaci |
| `flag_dw` | real | `0.0*` | **`flag_dw`** | keep |  | former common/newparam/: all 30 NAMELIST values (the "intuitively named" timestep/tolerance parameters replaci |
| `flag_dz` | real | `0.0*` | **`flag_dz`** | keep |  | former common/newparam/: all 30 NAMELIST values (the "intuitively named" timestep/tolerance parameters replaci |
| `lstruct_time` | logical | `.false.*` | **`lstruct_time`** | keep |  |  |
| `time_core_min` | real | `0.0*` | **`time_core_min`** | keep |  |  |
| `time_dl` | real | `0.0*` | **`time_dl`** | keep |  |  |
| `time_dp` | real | `0.0*` | **`time_dp`** | keep |  |  |
| `time_dr` | real | `0.0*` | **`time_dr`** | keep |  |  |
| `time_dt` | real | `0.0*` | **`time_dt`** | keep |  |  |
| `time_dw_global` | real | `0.0*` | **`time_dw_global`** | keep |  |  |
| `time_dw_mix` | real | `0.0*` | **`time_dw_mix`** | keep |  |  |
| `time_dx_core_frac` | real | `0.0*` | **`time_dx_core_frac`** | keep |  |  |
| `time_dx_core_tot` | real | `0.0*` | **`time_dx_core_tot`** | keep |  |  |
| `time_dx_shell` | real | `0.0*` | **`time_dx_shell`** | keep |  |  |
| `time_dx_total` | real | `0.0*` | **`time_dx_total`** | keep |  |  |
| `time_dy_core_frac` | real | `0.0*` | **`time_dy_core_frac`** | keep |  |  |
| `time_dy_core_tot` | real | `0.0*` | **`time_dy_core_tot`** | keep |  |  |
| `time_dy_shell` | real | `0.0*` | **`time_dy_shell`** | keep |  |  |
| `time_dy_total` | real | `0.0*` | **`time_dy_total`** | keep |  |  |
| `tol_czbase_fine_width` | real | `0.0*` | **`tol_czbase_fine_width`** | keep |  |  |
| `tol_dl_max` | real | `0.0*` | **`tol_dl_max`** | keep |  |  |
| `tol_dm_max` | real | `0.0*` | **`tol_dm_max`** | keep |  |  |
| `tol_dm_min` | real | `0.0*` | **`tol_dm_min`** | keep |  |  |
| `tol_dp_core_max` | real | `0.0*` | **`tol_dp_core_max`** | keep |  |  |
| `tol_dp_czbase_max` | real | `0.0*` | **`tol_dp_czbase_max`** | keep |  |  |
| `tol_dp_env_max` | real | `0.0*` | **`tol_dp_env_max`** | keep |  |  |
| `tol_dx_max` | real | `0.0*` | **`tol_dx_max`** | keep |  |  |
| `tol_dz_max` | real | `0.0*` | **`tol_dz_max`** | keep |  |  |
| `time_max_dt_frac` | real | `0.0*` | **`time_max_dt_frac`** | keep |  |  |
| `lnewvars` | logical | `.false.*` | **`use_new_variables`** | todo |  | TODO verify |
| `lnewtcz` | logical | `.false.*` | **`use_new_turnover_timescale`** | canonical |  | former common/ovrtrn/'s two config members: NAMELIST /physics/ values selecting the newer convective-turnover- |
| `lcalcenv` | logical | `.false.*` | **`calc_envelope_flag`** | canonical |  |  |

