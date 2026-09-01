# eos/opal sweep findings (reconstructed from the agent's final report
# after the scratchpad purge; the headline finding independently
# verified in the main session)

## eos/opal/esac.f90:117-123 -- OPAL95 EOS first call always errors (VERIFIED)
- class: logical
- severity: high
- confidence: high (verified by main session against 6cd5673:src/esac.f)
- provenance: modernization
- detail: the original's `IF(Z+XH-1.D-6 .GT. 1.0D0) GO TO 61` became a
  single-statement logical IF guarding only the first write; the second
  diagnostic write, `ierr = 1`, and `return` are UNCONDITIONAL inside the
  first-call table-load block. Every use_opal95_eos run dies at its first
  OPAL-regime EOS call (eqstat.f90:614 propagates). esac01/esac06 converted
  the identical construct correctly (block IF), proving it accidental.
  Invisible to the battery: no pinned case uses the 1995 vintage.

## eos/opal/esac.f90:219 -- X-dimension over-read for X > 0.6
- class: logical
- severity: medium
- confidence: high
- provenance: inherited (original `DO M=MF,MF+3`)
- detail: scans x_index_lo..x_index_lo+3 reading eos_table(6,...) -- one
  past the X dimension -- on every 1995-EOS call with X>0.6 (all solar
  calls). Benign today only because the aliased element (variable 2 at X=0)
  shares the same hole structure. The 01/06 vintages bound the loop at mf2.

## eos/opal (state) -- shared table_loaded_flag across the three vintages
- class: logical
- severity: medium
- confidence: high
- provenance: inherited (shared common/lreadco/)
- detail: one table_loaded_flag serves OPAL95/2001/2006; enabling two
  vintages in one run silently evaluates a vintage whose tables were never
  loaded.

## Low findings (agent-reported, not yet independently verified)
- t6rinterp.f90: reads an undefined `dix` (SAVE dropped in modernization;
  value provably dead, but the in-file comment falsely claims SAVE
  persists). provenance: modernization (comment) / dead either way.
- esac*.f90 (all three): the >1e20 overflow diagnostic tests the wrong
  array (inherited ESK-vs-ESK2 slip).
- rhoofp: secant can divide by zero when trial pressures coincide
  (self-limits to the -999 sentinel). inherited.
- eqbound: "edge by linear interpolation" comment describes an abandoned
  scheme (dead t_fraction). inherited.
- esac01/06: Cv scaled with 83.14511 on the no-rad path but 83.14510 in
  radsub01/06 (inconsistent gas constant at the 7th digit). inherited.

## Weak/uncertain observations
- (chased and cleared by the agent: X<1e-6 index block, iqu/ipu prechecks,
  boundary resets, rhoofp06 deviations, gmass electron term, sp literals --
  all faithful to originals or documented as deliberate.)
