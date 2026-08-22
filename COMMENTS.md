Where YREC's speed actually comes from — not laziness, but real specialization: a fixed ~15-species composition (no configurable nuclear network), a 2-variable (P,T) Henyey matrix with the outer envelope/atmosphere solved separately by analytic integration rather than folded into the same linear system, and narrower-purpose EOS/opacity tables. MESA pays for its generality (arbitrary isotopes, WDs, extreme metallicities, radiative levitation, etc.) with a much bigger state vector and matrix per zone. This tradeoff is legitimate engineering, not something to "fix."

Where the code shows its age, concretely:

Global COMMON-block state with no enforced interface. ~90 COMMON blocks, 500+ members, matched positionally across files with zero compiler-checked consistency. I hit this constantly while building the pulsation feature — the only way to know what sesum/so/del_grad actually held was to trace assignment sites by hand in coefft.f90, because nothing declares "here's what this subroutine reads and writes."
No derived types / no allocatable arrays. Every subroutine that touches shell data re-declares double precision :: x(json) with json=5000 hardcoded — duplicated 40+ times. Changing the max grid size means touching dozens of files.
Output is fixed FORTRAN WRITE format statements, not data-driven. MESA lets you add a history/profile column by editing a text list; in YREC, adding one column means editing a format string and a write-list in lockstep in the exact right file (I just did this for the GYRE writer — real friction, not hypothetical).
EOS/opacity/nuclear are not swappable modules — they're interleaved with the solver via the same COMMON blocks the Henyey coefficients are built from (coefft.f90 reads/writes scrtch/pulse1 directly). MESA's eos/kap/net are behind documented call interfaces precisely so they can be swapped.

Highest-leverage single improvement, given your stated goal of eventually swapping in MESA physics modules: don't rewrite the solver — it works and is fast. Instead, wrap the EOS/opacity/nuclear COMMON-block accesses behind a small number of explicit-interface subroutines (rho,T in → P,S,Γ1,χT,χρ,cp,κ,dlnκ/dlnT,dlnκ/dlnρ,εnuc,... out), called from coefft.f90 instead of reading globals directly. That's a bounded, mechanically verifiable refactor (same regression-diff discipline as this session), and it's the actual precondition for physics-module swapping — right now that goal is blocked not by the solver but by the fact that "the EOS" isn't a thing with a boundary, it's smeared across a dozen COMMON blocks.




