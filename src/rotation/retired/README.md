# rotation/retired — commented-out code, kept for reference

Not compiled: the `.txt` suffix keeps it out of the Makefile's
`find . -name '*.f90'`. Both paths were commented out by the author's
decision, and were moved here from the live files on 2026-10-28 so that
those read as the code that runs.

| File | What |
|---|---|
| `advection_diffusion.f90.txt` | The advection-diffusion transport of angular momentum (`use_diffusion_advection_transport`, legacy LDIFAD; Maeder & Zahn 1998), commented out 2026-10-12: `am_advection_diffusion_coeffs` and `banded_solver`, then every fragment of it that was in another file, each under a header naming the file and the live line it preceded, and the scratch fields it used. It never worked in the modernized code (most of its coefficient arrays were never filled); a rotating run that sets the flag stops. |
| `burn_settle_mix.f90.txt`, `burn_mix_extrapolated.f90.txt`, `burn_li_be_midpoint.f90.txt` | The Bulirsch–Stoer burn-and-mix path (MHP 6/00), commented out; each header says how to revive it. |
| `bulirsch_stoer_calls.f90.txt` | That path's two call sites in `evolve_rotation.f90`. |

To revive either path, put its pieces back where their headers say,
uncomment them, and (for the advection-diffusion transport) fill the
arrays its banner lists as unset.
