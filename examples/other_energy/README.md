# An extra energy source

This example adds a heating term to the luminosity equation without
touching YREC's source. It is the template for tidal heating, magnetic
heating, an injected energy, or any other source that is not nuclear,
neutrino or gravothermal.

## How it works

`src/core/structure/energy_sources.f90` keeps a list of extra energy
sources. A source is a subroutine with this argument list:

```fortran
subroutine my_source(pt, eps, deps_dlnrho, deps_dlnt, ierr)
   use energy_sources, only: energy_point
   type(energy_point), intent(in) :: pt   ! m, dm (g), log10 T, P, rho, r (cgs), L (Lsun), xa, dt (s)
   double precision, intent(out) :: eps          ! erg/g/s; positive heats
   double precision, intent(out) :: deps_dlnrho  ! d eps/d ln rho at constant T
   double precision, intent(out) :: deps_dlnt    ! d eps/d ln T at constant rho
   integer, intent(out) :: ierr                  ! nonzero stops the solve
end subroutine
```

A program registers it with `add_energy_source('name', my_source, ierr)`
and then calls `run_yrec`, which reads the deck named on the command line
as the `yrec` executable does. At every shell of the interior model YREC
adds the sources' sum to dL/dm and their partials to the Newton
Jacobian. With no source registered nothing is evaluated, and every
result is bit-for-bit that of the plain `yrec`.

`yrec_extra_heat.f90` is such a program. Its source, `uniform_heat`,
heats at a constant rate `eps_extra` (erg/g/s) between two mass
coordinates, read from `extra_heat.nml`.

## Running it

From `src/`:

```
make other_energy
```

builds `examples/other_energy/yrec_extra_heat`. Then, from this
directory:

```
./yrec_extra_heat solar_heated.inlist
```

`solar_heated.inlist` is the testsuite's base solar model
(`testsuite/Test_solarGS98_base.inlist`) with its output renamed. The
history file `output/history.data` gets the columns of
`history_columns.list`, among them `lum_extra`, the luminosity of the
extra sources (Lsun). With the shipped `extra_heat.nml`, 0.5 erg/g/s
between 0.3 and 0.6 Msun adds about 0.08 Lsun. Set `eps_extra = 0` and
the run reproduces `../../src/yrec solar_heated.inlist` exactly.

## Writing your own

1. Copy `yrec_extra_heat.f90` and replace `uniform_heat`. Keep it a
   module procedure: an internal procedure cannot safely be registered.
2. Give honest partials. They do not change the converged model, but
   wrong ones slow or stop the Newton iteration.
3. Keep the source inside the interior model. Above the fitting point
   the envelope is integrated at constant L, so a source there has no
   effect.
4. Add a target for your program next to `other_energy` in
   `src/Makefile`.

The profile column `eps_extra` gives the sources' sum zone by zone; add
it to a profile column list to write it.
