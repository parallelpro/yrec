# Round 1's measured fixes (Kind B), for the owner's approval (2026-11-10)

PLAN.md section 7: these run after round 1's byte-identical work is merged, one package at a time. Each is its own
commit with its measured change on the 56 decks and 35 pins, and the pins reseeded in it. Numbers are
docs/bug_sweep.md's. "Decision" marks the items that depend on what the code was meant to do; each carries a
recommendation. core's list comes with round 2.

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

## atm
1. O19, both halves: the Allard reader tests its real flag and returns an error; the accretion surface keeps its
   lookup in locals, checks it, and takes the envelope integration's photosphere unless the Allard tables are
   loaded. No pinned deck accretes.
2. M235: the Allard kernel's per-lookup notes leave the run log (1,920 of the 7,383 lines of the Allard pinned log).
3. Decision, M233: the Allard extension is four column widths in the code and one in its headers.
   **Recommended:** keep the code and correct the headers (nothing moves).
4. atm's messages reworded ("Program Terminated" banners, capitals), one re-pin commit.

## convection
1. Decision, M241: semiconvection's core-only switch takes the first mixed zone as the core whatever it is.
   **Recommended:** keep semiconvection for the core only, as the 2025 change meant, and act only when zone 1
   starts at the centre; no envelope semiconvection. Moves decks with semiconvection and no convective core (no
   pinned deck).
2. M242: the starspot Jacobian takes the spotted gradient and its d ln S/d ln T. Moves iteration counts and
   round-off in the three starspot testsuite decks.
3. M249: the adiabatic overshoot region is cleared when there are no zones (no shipped deck reaches it).
4. One H_P and one overshoot extent, owned by convection/ (the differing copies of the plan's CD-14), with O10's
   gate and width. Measured; moves the decks with envelope overshoot where the copies differ.
5. M96: the O16 floor's metals recomputed at a later card's Z rescale (with core's M85 in round 2). Moves only
   decks with a later Z rescale, where the floor binds.
6. M254's `[` format (all 12 pinned testsuite logs) and convection's messages reworded, one re-pin commit.

## element_diffusion
1. Decision, M243-M246: the original settling routines (no shipped deck selects them: 0 of 354 settling decks)
   carry four inherited defects. **Recommended:** retire them to element_diffusion/retired/ with a header, and make
   the control that selects them stop the run with an error that names the retired file. The alternative is to fix
   the four and test them, which keeps a second pipeline nothing runs. Moves no shipped deck either way.
2. M33: the metal coefficients zeroed outside the Thoul branch, as the F77 does, so that test_element_diffusion's
   pin is deterministic (the test's pin is reseeded; no deck moves).
3. O21: both settling grids get at least 4 points. Moves only the small regions that misread today.
4. Decision, M250: settling_z_factor also scales Li and Be. **Recommended:** keep it (Li and Be settle as metals)
   and say so in the description (nothing moves).
5. Decision, M181: `diffuse_lithium` turns metal settling and five other switches on. **Recommended:** option (c),
   a warning whenever the override changes a value: the behaviour stays, 9 pinned logs gain a line, and a reader
   of the log sees what ran. Option (b), touching only Li and Be, would move 169 decks' results (9 pinned).
6. element_diffusion's messages reworded, one re-pin commit.

## eos
1. M206: SCV's smoothing guarded against the unbuilt entries below a short row. Moves only lookups above the lower
   row's highest pressure (no shipped deck known).
2. M88: SCV's table edge becomes a weighted ramp, as OPAL's edges are (with the source descriptor). Moves models
   that cross the SCV edge.
3. M92: SCV's general path takes the stencil's mixing entropy at the zone's X before the envelope table is built.
   Moves the SCV decks' first-card refits and general-path lookups away from the envelope's X.
4. M94: OPAL 1995's loops bounded by its table (as OPAL 2001 and 2006 since M41). test_eos's OPAL 1995 pins may
   move; no deck uses OPAL 1995.
5. M128: SCV's H+ fraction counts an H2 molecule as two nuclei. Moves the low-mass decks where SCV covers H2 and H+,
   through the conductive opacity.
6. Decision, M129: one 1/μ_ion rule on every path under diffusion of Z. **Recommended:** OPAL and MHD take the
   Yale and SCV scaling (zone_mixture's rule). That moves only the 1/μ_ion and 1/μ_e output columns of the 355
   diffusion-of-Z decks, not their structure. The alternative, the zone's composition everywhere, moves the Yale
   and SCV densities of those decks.
7. eos's capitalized messages and "WARNNING... OPAL TBL FAIL" reworded, one re-pin commit.
