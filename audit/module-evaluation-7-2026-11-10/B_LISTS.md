# Round 1's measured fixes (Kind B), for the owner's approval (2026-11-10)

PLAN.md section 7: these run after round 1's byte-identical work is merged, one package at a time. Each is its own
commit with its measured change on the 56 decks and 35 pins, and the pins reseeded in it. Numbers are
docs/bug_sweep.md's. "Decision" marks the items that depend on what the code was meant to do; each carries a
recommendation. core's list comes with round 2.

**The owner's answer (2026-11-11):** "agreed on all except don't do anything to element_diffusion".
- kap, atm, convection and eos: every item is approved, with its recommendation as written.
  - An item marked "no recommendation yet" (kap item 8) is researched and brought back with one before any
    code changes.
  - convection item 5 (M96) waits for core's M85 in round 2.
  - The control descriptions in the last section go with their packages' B work.
- element_diffusion: nothing is done. Its section below stays as a record, and so do its controls'
  descriptions and its open bug_sweep entries.

**Status after round 1B (2026-11-12).** kap, atm, convection and eos were done one package at a time, each item
its own measured commit, merged into yrec-next in d66dd005, c5dd0ab3, 6fcc96cb and 1395075b (docs/ROADMAP.md,
round 1B; each bug_sweep entry gives its commit).
- Done:
  - kap: 1 O18 (0025d96f), 2 M222 (c797b801), 3 M237's banner (0f6f4978), 4 M227 (c92127bd, the code kept),
    5 the messages (1220787b), 6 KA-14 (2fec5fed; test_structure re-pinned in c162117f), 7 M277 (758819a0).
  - atm: 1 O19 (1f2b14fb), 2 M235 (aae0218b), 3 M233 (f037c9bd, the code kept), 4 the messages (12210cb9),
    5 M283 (67cd7815), 6 the gray constant (c9e6ec4e).
  - convection: 1 M241 (5a370511), 2 M242 (e6f70a34), 3 M249 (a21950d9), 4 CD-14 with O10's gate and width
    (4c1c730a, 86b54640, 4cf4e0ac), 6 M254 and the messages (76a19fc7), 7 M285 (f80adb7b).
  - eos: 1 M206 (95126040; its check, b8ecbef3), 2 M88 (a0039581), 3 M92 (346e1af2), 4 M94 (e6a317cf),
    5 M128 (b8c4b8cf), 6 M129 (15b8b199, the recommended rule), 7 the messages (e5a03a6f), 9 M290 (46bcab10),
    10 M291 (03be24fe), 11 M292 (4b9aecff); and, not on the list, the fully ionized iteration's overflow capped
    (1bd927eb, bug_sweep M305).
  - The descriptions of use_pure_z_table and laol_extrapolation_tolerance (c92127bd) and atm_table_file
    (4499f60a). element_diffusion's descriptions are not done.
- Waiting:
  - kap item 8, M278, on the owner: researched (bug_sweep M278: Potekhin et al.'s rule; B recommended, A the
    weights alone, or leave it); its ZAHB trial is to be measured again on today's base.
  - convection item 5, M96: in round 2, with core's M85.
  - eos item 8 (the Saha non-convergence), on the owner: counted (2 of the 56 decks, 2,217 times, all in the
    starting model's envelope refit, whose order is M86, proposed for core's round 2).
  - Also for the owner: keeping M285's fix, which is not latent (ZAHB_m1p0feh+0p0_TAHB's core helium burning
    9.8% longer).
- Two statements below proved wrong: convection item 1's TAHB sample moves the lower edge of an intermediate
  zone, not below the envelope's base, and item 7's merge is reached by ZAHB_m1p0feh+0p0_TAHB.

## kap
1. O18: the opacity ramp's derivatives made the derivatives of the blend. Moves every deck whose models cross the
   ramp, at the solver's level.
2. M222: the conductive opacity's low-density branch takes the mixture term's sign and extrapolates the right way.
   Every deck's photosphere reaches it; conduction is negligible there, so the change is at round-off in κ, but
   every pin moves.
3. M237: the banner says "pure C and O tables" for the C and O table (35 pinned logs move by one line).
4. Decision, M227: the LAOL89 extrapolation tolerance is applied in ln and described in log10. **Recommended:**
   keep the code and correct the description (nothing moves).
5. kap's messages reworded: sentence case, no "RUN STOPPED" on paths that return ierr (one re-pin commit).
   Every error goes to both the terminal and the run log, and the lookups return a reason instead of writing where
   they fail. test_kap's pinned "LOG T OF … OUT OF OPAL 95" line moves.
6. KA-14: the conductive opacity's −3.5194 computed from σ as log10(16σ/3). It differs at 1e-5, so conduction moves
   at round-off (every pin, as with M222; one measured commit with it).
7. New (round 1): the OPAL95 high-R density shift (opal95_opacity.f90:368) searches the Z = 0 slot of the full table
   set instead of the fixed-Z table. The F77 has it (getopal95.f). **Recommended:** fix it in its own measured
   commit, and decide from what it moves.
8. Decision, new (round 1): the conductive mixture (conductive_opacity.f90:110-111, 143) weights each species by
   number fraction times ionized fraction without renormalizing, so a less ionized mixture conducts more. At
   log T 7.3 and log ρ 4.5, κ is 0.137 fully ionized and 0.0013 with no ionization. The F77 has it (CONDOPACPINT).
   No recommendation yet: the fix should follow the published mixing rule the routine implements. I would check
   that rule first, then propose the change with its measured effect, next to item 2.

## atm
1. O19, both halves: the Allard reader tests its real flag and returns an error; the accretion surface keeps its
   lookup in locals, checks it, and takes the envelope integration's photosphere unless the Allard tables are
   loaded. No pinned deck accretes.
   - Round 1 measured it. The tau = 100 check (now in atm/retired/) moves no deck.
   - The accretion half moves cm_accretion, a compare deck whose log has 24 "ALSURFP" pairs.
   - atm's side is ready, with a named error for a lookup in an unread table. What accretion's surface should use
     is core's decision.
2. M235: the Allard kernel's per-lookup notes leave the run log (1,920 of the 7,383 lines of the Allard pinned log).
3. Decision, M233: the Allard extension is four column widths in the code and one in its headers.
   **Recommended:** keep the code and correct the headers (nothing moves).
4. atm's messages reworded ("Program Terminated" banners, capitals), one re-pin commit. Since round 1 every text is
   in atm_messages. The commit also sorts which messages go to standard output and which to the terminal unit. It
   moves test_atm's pinned output and cm_accretion's log.
5. New (round 1): atm_init keeps the tables a previous run in the same process read. It moves no shipped deck.
   **Recommended:** restart from a clean state at each atm_init, as kap's round did for kap_init.
6. Optional: the gray T-tau constant written as ¼ log10(3/4) instead of −0.031235. It moves every gray and
   Krishna-Swamy deck by 3.2e-7 dex. **Recommended:** measure it on its own, then commit it straight after kap's
   M222, which reseeds every pin anyway. It is not worth a reseed of its own.

## convection
1. Decision, M241: semiconvection's core-only switch takes the first mixed zone as the core whatever it is.
   **Recommended:** keep semiconvection for the core only, as the 2025 change meant, and act only when zone 1
   starts at the centre; no envelope semiconvection. Also return when there is no mixed zone. Moves decks with
   semiconvection and no convective core (no pinned deck).
   - It is not latent (round 1, measured at 916ee23a): the shipped
     `sample_cases/UpperMS/Sample_m9p0feh+0p0GS98_yaleeos_sc_TAMS` stops in model 110. From card 2 on there is no
     mixed zone, and zone 1's bounds are left over from an earlier model (an upper edge at shell 707 of 648).
   - Each call moves that edge one shell further and writes densities into the shells past the model's end.
     The run then exits 1 on a NaN gradient.
   - The TAHB sample moves its edge 16 times below the envelope's base.
   - The fix lets the TAMS sample run and moves the TAHB sample; neither is pinned.
2. M242: the starspot Jacobian takes the spotted gradient and its d ln S/d ln T. Moves iteration counts and
   round-off in the three starspot testsuite decks.
3. M249: the adiabatic overshoot region is cleared when there are no zones (no shipped deck reaches it).
4. One H_P and one overshoot extent, owned by convection/ (the differing copies of the plan's CD-14), with O10's
   gate and width. Measured; moves the decks with envelope overshoot where the copies differ.
5. M96: the O16 floor's metals recomputed at a later card's Z rescale (with core's M85 in round 2). Moves only
   decks with a later Z rescale, where the floor binds.
6. M254's `[` format (all 12 pinned testsuite logs) and convection's messages reworded, one re-pin commit.
   - The plan's CD-11 formats join it (`1p3e12.6` to `1p3e13.5`, `i3` to `i5`): each changes message text.
   - The square-root warning also goes to the run log; today it reaches only the terminal.
7. New (round 1): the overshoot merge can shrink a region. When two regions merge, the merged one ends at the
   upper region's top even if the lower region's overshoot reached higher; the fix takes the higher of the two.
   - Where: convection_zones.f90's merge_overlapping_zones.
   - The F77 has it too (convec.f).
   - No deck is known to reach it.

## element_diffusion
1. Decision, M243-M246: the original settling routines (no shipped deck selects them: 0 of 354 settling decks)
   carry four inherited defects. **Recommended:** retire them to element_diffusion/retired/ with a header, and make
   the control that selects them stop the run with an error that names the retired file. The alternative is to fix
   the four and test them, which keeps a second pipeline nothing runs. Moves no shipped deck either way.
   - Round 1 found a fifth defect, which the F77 has too (grsett.f): they never solve the metals' diffusion term.
     The metals' D2 is computed and not used, so the metals settle by drift alone.
2. M33: the metal coefficients zeroed outside the Thoul branch, as the F77 does, so that test_element_diffusion's
   pin is deterministic (the test's pin is reseeded; no deck moves).
3. O21: both settling grids get at least 4 points. Moves only the small regions that misread today.
4. Decision, M250: settling_z_factor also scales Li and Be. **Recommended:** keep it (Li and Be settle as metals)
   and say so in the description (nothing moves).
5. Decision, M181: `diffuse_lithium` turns metal settling and five other switches on. **Recommended:** option (c),
   a warning whenever the override changes a value: the behaviour stays, 9 pinned logs gain a line, and a reader
   of the log sees what ran. Option (b), touching only Li and Be, would move 169 decks' results (9 pinned).
6. element_diffusion's messages reworded, one re-pin commit. The numerics-failure message also goes to the run log
   (CD-15).
7. Moved from round 1's A list: M247 (CD-4). Its fix moves one pinned deck,
   Test_solar_GS_norot_gray_OPALSCV_SF2_GS98_OP_CF10, whose pre-main-sequence convective core reaches the path.
   - From model 250 on, log L and log Teff move by about 1e-5 dex; its .mod file is identical.
   - The fix is one line in diffusion_grid's build_equal_grid_map, with a convective-core check added to the test.
8. Decision, new (round 1): the settling grids get one point too few. They take ceil(span / narrowest zone) points,
   so the intervals are wider than the narrowest zone; the F77 comment ("DR <= DRMIN") wants one more.
   - The F77 has it (microdiff_mte.f, model_to_equal.f).
   - **Recommended:** fix it (one more point, as the comment says). It moves every settling deck.
9. The single-precision 2.2: microdiff divides by `2.2` and the original routines by `2.2d0`. Making it `2.2d0`
   moves every settling deck by about 2e-8 relative. **Recommended:** do it in the same measured commit as item 8.
10. `settling_num_iterations` < 1 skips the implicit solve and reports a non-convergence with change 0. The F77
    has it too. **Recommended:** the control check rejects values below 1. No shipped deck sets it.

## eos
1. M206: SCV's smoothing guarded against the unbuilt entries below a short row. Moves only lookups above the lower
   row's highest pressure (no shipped deck known).
2. M88: SCV's table edge becomes a weighted ramp, as OPAL's edges are (with the source descriptor). Moves models
   that cross the SCV edge.
3. M92: SCV's general path takes the stencil's mixing entropy at the zone's X before the envelope table is built.
   Moves the SCV decks' first-card refits and general-path lookups away from the envelope's X.
4. M94: OPAL 1995's loops bounded by its table (as OPAL 2001 and 2006 since M41). test_eos's OPAL 1995 pins may
   move; no deck uses OPAL 1995.
   - Since round 1, OPAL 1995 runs on the shared OPAL code. It keeps its old numbers by emulating the reads past
     its table (`stencil_entry`, `unclamped_coverage` and the T6-node line).
   - The fix sets `clamps_stencil` for 1995 and deletes those three.
   - Measured: test_eos's 1995 pin at log T 3.70, log P 4.5 moves by 2.4e-6 dex in log ρ.
5. M128: SCV's H+ fraction counts an H2 molecule as two nuclei. Moves the low-mass decks where SCV covers H2 and H+,
   through the conductive opacity.
6. Decision, M129: one 1/μ_ion rule on every path under diffusion of Z. **Recommended:** OPAL and MHD take the
   Yale and SCV scaling (zone_mixture's rule). That moves only the 1/μ_ion and 1/μ_e output columns of the 355
   diffusion-of-Z decks, not their structure. The alternative, the zone's composition everywhere, moves the Yale
   and SCV densities of those decks.
7. eos's capitalized messages and "WARNNING... OPAL TBL FAIL" reworded, one re-pin commit. Also in it:
   - the SCV read-failure message, which always names the H table;
   - OPAL 2001's ihi/ilo line, which prints no indices;
   - M93's error, which since round 1 goes to the run log only and should reach the terminal too.
8. Decision, new (round 1): the Saha iteration's non-convergence is a warning, and decks reach it. Making it an
   error would stop those runs. **Recommended:** keep it a warning, count how often each deck reaches it, and
   decide from that count.
9. New (round 1): the Yale fully ionized iteration diverges from a fresh record (log ρ 0) at log T 8.4 to 9.35
   and log P 21 to 23.5; from a nearby start it converges. Before M93 this returned success; now it is an error.
   No shipped deck reaches it. The F77 has it. **Recommended:** start from the previous zone's density or a
   bracketing guess, measured.
10. New (round 1): OPAL's table reader (readcoeos) keeps the last X block's T6 grid and never compares the grids
    across blocks; a check would refuse a table whose blocks differ. The shipped tables agree, so no deck moves.
    The F77 has it.
11. New (round 1): in OPAL's X interpolation, a huge entry that enters with a negative Lagrange weight passes all
    three range checks, and a negative value comes back. Real tables cannot reach this; the F77 has it. Recommended:
    a check on the result's sign, which moves no deck.

## Control descriptions (all packages, one commit)
The registry's descriptions are copied into the generated `defaults/controls_registry.md` and `control.defaults`,
so these go in one commit with the generators rerun, after the owner's decisions above:
- `use_pure_z_table` (io/controls_buffer.f90:323-324, M237): "a pure mixture of equal numbers of C and O, at Z = 1".
- `laol_extrapolation_tolerance` (:809-810, M227, if the code is kept): the distance in log10 ρ is the natural log
  of this value (2.3 dex for 10).
- `settling_tolerance` and `settling_num_iterations` (:1144, :1146, M253): only the original routines iterate;
  the new routines solve once, and these decide only the solver-diagnostics lines and the "failed to converge"
  message.
- `constant_mixing_coeff` and `constant_settling_reduction` (:1162, :1164, M258): read only by the original
  routines with Thoul et al.'s coefficients.
- `settling_z_factor` (:1141-1142, M250): also scales Li and Be in the new routines.
- `atm_table_file` (:405, M239): "Kurucz (atm_choice 3)", not "3 and 4" (and setup/settings_install.f90:55-56's
  comment, which moves nothing).
