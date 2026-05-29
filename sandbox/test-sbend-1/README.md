# test-sbend-1: single analytic SBEND placement consistency test

This directory is a minimal regression seed for one explicitly placed analytic
`SBEND` with a `0.3 m` upstream drift.  It starts with a one-sided entrance
pole-face case: the exported entrance port is centered at `(0, 0, 0.3)`, the
single particle is started exactly at `(0, 0, 0)`, and the exit monitor remains
at the hard SBEND exit.  Since the particle is on the reference point, the
`E1 = 45 deg` pole face should not move the particle away from the design exit.

## Analytic configuration

| quantity | value |
| --- | ---: |
| upstream drift | `0.3 m` |
| `L` | `1.0 m` |
| `ANGLE` | `45 deg` |
| `E1` | `45 deg` |
| `E2` | `0 deg` |
| `HGAP` | `0.02 m` |
| `FINT` | `0.5` |

The input uses explicit element pose attributes `X`, `Y`, `Z`, `THETA`, `PHI`,
and `PSI`.  The drift starts at the input start, the SBEND entrance edge is at
the drift exit, and the entry fringe field-begin point lies one fringe-support
length upstream of the SBEND edge.  The monitor is placed at the exported SBEND
exit edge.

The analytic fringe profile follows the OPAL `FM1DProfile1` default Enge model.
With no explicit `GAP`, OPALX uses `GAP = 2 * HGAP = 0.04 m`, so the default
profile half width is `5 * GAP = 0.2 m`.  For this `E1 = 45 deg`, `E2 = 0`
case the local SBEND field support is approximately
`[-0.2828427, 1.2261722] m` in the bend chart.  In lab coordinates the field
starts at `(0, 0, 0.0171573) m` and ends near
`(-0.5241048, 0, 1.3653009) m`.  The integrated Enge effective length remains
close to the old-OPAL chord-derived arc length, `1.0261729 m`.

## How to run

From this directory:

```sh
/Users/adelmann/git/opalx/build/src/opalx test-sbend-1.in
```

The useful diagnostic for this placement check is `test-sbend-1_Monitors.stat`.
The `.loss` file can be header-only for this one-particle monitor case, so it
is not the primary acceptance file.

## Practical acceptance check

Use the monitor row with the smallest local in-plane residual.  In the current
monitor statistics format, `mean_x`, `mean_y`, and `mean_s` are columns 15, 16,
and 17:

```sh
awk '$1=="test-sbend-1_exit_monitor" {
  r=sqrt($15*$15 + $17*$17)
  if (n++==0 || r<best) { best=r; mx=$15; my=$16; ms=$17; t=$3 }
}
END {
  print "best monitor residual =", best, \
        "mean_x =", mx, "mean_y =", my, "mean_s =", ms, "t =", t
}' test-sbend-1_Monitors.stat
```

The placement check is considered good when:

- OPALX exits with status 0.
- `test-sbend-1_Monitors.stat` is produced.
- The best monitor-plane residual `sqrt(mean_x^2 + mean_s^2)` is below
  `5.0e-5 m`.
- The output contains no `nan` and no OPALX error.
- `data/test-sbend-1_ElementPositions.txt` shows `ENTRY EDGE: B1` at the
  drift exit.

The expected exit monitor pose for this on-axis entrance-edge case is
`X = -0.3826834324 m`, `Z = 1.2238795325 m`, and `THETA = -45 deg`.

## Old OPAL comparison

`test-sbend-1_opal.in` is the OPAL 2022.1 comparison deck.  It uses the legacy
`ELEMEDGE` bend placement syntax and writes `data/test-sbend-1_opal_DesignPath.dat`
plus `test-sbend-1_opal.stat`.  The old deck keeps `DESIGNENERGY = 1000 MeV`
to match the `1 GeV` electron beam momentum; using `1 MeV` makes the old OPAL
field about three orders of magnitude too weak for this comparison.
The old OPAL deck starts at `S_OFFSET = 1.0 m`, tracks through the same `0.3 m`
drift, and places the SBEND at `ELEMEDGE = S_OFFSET + 0.3 m`.

OPALX and old OPAL both write into a directory named `data`, so run them in
separate working directories before comparing.  The helper script does that:

```sh
./run_opal_opalx_compare.sh
```

It defaults to `/tmp/opalx-sbend-compare`, runs OPALX in `opalx/`, runs OPAL
2022.1 in `opal-2022/`, and writes CSV plus Markdown summaries in
`comparison/`, including PNG plots for the design path, momentum, and sampled
`B_y` field.  The comparison intentionally uses only `*_DesignPath.dat` and
`.stat`; `.h5` and `.lbal` are ignored.

The report includes a dedicated sampled-field check for the reference path:
`design_path_bfield_summary.csv` compares the `Bfx/Bfy/Bfz` columns from
`*_DesignPath.dat`, and `stat_bfield_summary.csv` compares the `.stat`
`Bx_ref/By_ref/Bz_ref` columns.  These summaries include field peaks,
trapezoidal field integrals over the normalized path coordinate, and the
sampled field support interval.

Generated output is intentionally not kept in the repository unless this
placement check is later promoted to a reference regression with frozen
tolerances.
