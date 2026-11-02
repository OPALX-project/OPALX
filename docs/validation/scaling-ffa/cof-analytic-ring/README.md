# Analytic COF acceptance fixture

The 16 overlapping tanh sectors sum to a radial field B=B0*(r/R0)^(-0.5)
in the median plane. Momentum is chosen so r=R0=4.42 m is a circular orbit.
The field index -0.5 gives weak focusing in both planes, so the 4D fixed point
is nonsingular. This is independent of the school optics and probes actual
COF command acceptance, occurrence initialization and field evaluation.

Copy both input files into a fresh run directory, then run OPALX on `ring.in`.
For this machine:

```sh
OMP_NUM_THREADS=1 OMP_PROC_BIND=false mpiexec --host localhost:4 --map-by slot --bind-to none -n 1 /Users/adelmann/git/opalx/build-omp/src/opalx --info 1 ring.in
```

The command writes `orbit.json` and must converge from its displaced 4D seed.
All four section coordinates should be within 1e-9 of zero in their respective
units (metres and p/(mc)); the recorded check gave <3e-14 m. This loose check
bound includes the configured Newton and return-map rounding errors and does
not modify the command's convergence tolerances. No TRACK is executed, so the
FROMFILE particle data is unused by the empty reference container.
