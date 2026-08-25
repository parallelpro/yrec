!----------------------------------------------------------------------
! main
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original main.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). This is the LAST file converted
! in the project (256 of 257 files were already converted before this
! one); every callee below was cross-checked against its own current
! (already-converted) signature, argument by argument.
!
! This is the top-level PROGRAM unit (not a SUBROUTINE) -- the original
! main.f has no PROGRAM statement (a legal, unnamed F77 main program);
! it is given the name "main" here, matching the file name, since
! free-form Fortran expects "end program <name>" to pair with
! "program <name>".
!
! MAIN drives the whole simulation: it reads user parameters (parmin)
! and physics tables (setups), then for each Monte-Carlo realization
! (if any) and each "kind card" run (the DO 200 NK loop) reads or
! rescales the starting model (starin), evolves it for the requested
! number of models (the DO 100 MODELN loop) via the Newton-Raphson
! relaxation driver (crrect, 4 iteration levels), mass loss (massloss),
! mixing (mix/mixcz), light-element burning (lirate88/liburn), the
! rotation curve (getw/fpft), and re-meshing (hpoint), writes output
! after every converged model (wrtout, plus wrtlst/putstore for stored
! models), and (for automatic solar/stellar calibration) checks and
! iterates the calibration parameters (chkcal/chkscal/setcal/setscal).
!
! CROSS-CALLEE NAMING NOTE: several locals here are threaded into more
! than one already-converted callee, and those callees do not always
! agree on a name for the same physical data (this file, being the
! PROGRAM unit, is free to choose its own local names; the callees'
! names were fixed by earlier batches). Judgment calls made below, all
! verified against the actual physics/usage in this file and cross-
! checked against each callee's current signature:
!   - BL (log10 L/Lsun, surface value) is named star%log_L,
!     matching crrect.f90/starin.f90/chkcal.f90/chkscal.f90's own
!     slot names for this exact scalar; wrtlst.f90/massloss.f90/
!     getnewenv.f90/putstore.f90/wrtout.f90 instead call the same
!     scalar "log_luminosity_lsun".
!   - HL (linear L/Lsun, per-shell array) is named star%luminosity_lsun,
!     matching crrect.f90/starin.f90; findsh.f90/htimer.f90/mix.f90/
!     hpoint.f90/getnewenv.f90/wrtlst.f90/wrtout.f90/putstore.f90/
!     massloss.f90/getw.f90 instead call the same array
!     "log_luminosity" (or plain "luminosity"/"star%luminosity_lsun" in
!     findsh.f90's own case) -- a misnomer inherited from those files'
!     own earlier conversions in most cases; out of scope to fix here.
!   - M (number of mass zones) is named star%nz, the majority
!     spelling among callees (massloss/mix/hpoint/getw/convec/
!     lirate88/liburn/getnewenv/mixcz); starin.f90/wrtlst.f90/
!     wrtout.f90/putstore.f90 instead call it num_shells,
!     crrect.f90/fpft.f90/findsh.f90/htimer.f90 call it num_points,
!     wrtmonte.f90 keeps it cryptic as "m".
!   - JENV (envelope CZ base zone index) is named
!     star%envelope_cz_bottom_index, matching wrtlst.f90/wrtout.f90/
!     putstore.f90/wrtmonte.f90 (majority); findsh.f90 calls it
!     envelope_edge, convec.f90/mix.f90 call it envelope_cz_edge,
!     hpoint.f90 calls it convective_envelope_edge_zone,
!     massloss.f90 calls it envelope_boundary_zone, starin.f90 calls
!     it envelope_zone_index.
!   - JCORE (core CZ top zone index) is named star%core_cz_top_index,
!     matching wrtlst.f90/wrtout.f90/putstore.f90/wrtmonte.f90
!     (majority); findsh.f90 calls it core_edge, convec.f90/mix.f90
!     call it core_cz_edge, hpoint.f90 calls it
!     convective_core_edge_zone, htimer.f90 calls it
!     convective_core_edge_zone as well.
!   - LARGE (T if the correction routine has diverged / the model has
!     failed) is named model_diverged_flag; crrect.f90's own slot name
!     for this same local is corrections_too_large, starin.f90's is
!     model_failed_flag.
!   - LNEW (T if a fresh envelope-fit triangle must be generated) is
!     named recompute_envelope_triangle; crrect.f90's own slot name is
!     start_new_triangle, starin.f90's is envelope_recomputed_flag,
!     getnewenv.f90's is new_points_added_flag.
!   - LRESET is named reset_triangle, matching crrect.f90 (used in 4 of
!     its 5 call sites here); hpoint.f90's own slot name is
!     point_reset_flag.
!   - DELTSH is a dual-purpose scratch slot: starin.f90 treats it as
!     delta_time_abs (= |delta_time|, set by starin itself) while
!     htimer.f90 treats it as hydrogen_dt (the H-burning-shell
!     timestep governor, computed fresh by htimer). It is named
!     hydrogen_dt here since that is its dominant role (read back by
!     the next model's htimer.f90 call); starin's transient use of the
!     same storage as |delta_time| is noted here for the record.
!   - ISTORE is named istore_flag, matching starin.f90; hpoint.f90's
!     own slot name for the same scalar is envelope_store_index.
!   - SENVOLD (set by massloss, read by getnewenv) is named
!     target_envelope_mass, matching getnewenv.f90 (how it is actually
!     consumed); massloss.f90's own slot name for the same scalar is
!     old_log_envelope_mass_fraction.
!   - V(12) is a dual-purpose work array: passed to setups.f90 as
!     laol_work_array (a scratch/table-loading array) and to
!     starin.f90 as species_mix_weights (the VNEW mixture-weight
!     copy). Named mixture_weights here as the more physically
!     meaningful of the two uses; flagged here since the two callees'
!     own names disagree on its purpose.
!   - FPATM/FPENV/FPMOD are named pulse_atm_path/pulse_env_path/
!     pulse_mod_path, matching pdist.f90 (their only consumer besides
!     parmin.f90, which keeps its own dummy spelling fpatm/fpenv/fpmod
!     per that file's NAMELIST-driven naming constraint -- COMMON/
!     argument passing is positional, so this does not need to match).
!   - The remaining PARMIN/SETUPS file-path locals (FALEX06, FALLARD,
!     FATM, FFERMI, FKUR, FKUR2, FLAOL, FLAOL2, FLIV95, FLLDAT, FMHD1-
!     FMHD8, FOPAL2, FPUREZ, FSCVH, FSCVHE, FSCVZ, OPECALEX) are named
!     to match setups.f90's own (descriptive) dummy-argument spelling
!     for the same slot, since parmin.f90 keeps its own dummy names at
!     cryptic lowercase spelling only (per that file's NAMELIST
!     constraint, documented in parmin.f90's own header) -- again,
!     positional COMMON/argument passing means this file's own local
!     names are free to differ from parmin.f90's.
!
! COMMON BLOCK NOTE: every COMMON block below was grepped against all
! other already-converted .f90 files and, where established, reuses
! the exact member names found there (per the project's COMMON-block
! rule). Two conflicts were found and resolved, both documented again
! at the point of declaration below:
!   - common/difus/'s 4th member (originally ITDIF2) is established as
!     "max_iterations" in checkc.f90, but this file already uses
!     max_iterations as a local name (matching crrect.f90's own slot
!     name for the NITER argument, threaded through all 4 of this
!     file's crrect calls) -- reusing it for the (here-unused) common
!     member would alias two different entities to the same name, so
!     the common member is instead kept at its lowercased-cryptic
!     spelling, itdif2, here.
!   - common/cals2/'s 1st member and common/calstar/'s 2nd member are
!     both independently established as "luminosity_tolerance" (in
!     chkcal.f90 and chkscal.f90 respectively) -- distinct physical
!     quantities (solar-calibration L tolerance vs. target-star L
!     tolerance) that never previously co-occurred in the same file.
!     common/cals2/'s member keeps the established "luminosity_
!     tolerance"; common/calstar/'s member is disambiguated here as
!     "target_star_luminosity_tolerance".
!   - Several common blocks establish "_placeholder"-suffixed names
!     for members that are unused in the file(s) that first declared
!     them, but that this file actually reads/writes (common/acdpth/'s
!     output_ages_gyr/calcad_ageout_output_active/ageout_model_output_flag/
!     ageout_bracket_armed; common/rotprt/'s star%lprt0_placeholder;
!     common/chrone/'s lrwsh_placeholder; common/cenv/'s
!     lnew0). Per the precedent set by atm_lib.f90 (keeps
!     calcad_ageout_output_active despite noting its own active use) and
!     getw.f90 (keeps star%print_rotation_diagnostics despite noting its own active
!     use), these established names are reused verbatim here too
!     rather than renamed, even though this file actively assigns/
!     reads them; each active use is called out in a comment at its
!     point of use below.
program main
! 2026 (phase five, step A): the body of this program moved verbatim
! into core/run_yrec.f90 -- the first cut of the embeddable engine.
! This wrapper is the CLI: run the configured job, and stop nonzero
! if the engine ever reports an error (phase B will start doing so).
      implicit none
      integer :: ierr

      ierr = 0
      call run_yrec(ierr)
      if (ierr /= 0) stop 1
end program main
