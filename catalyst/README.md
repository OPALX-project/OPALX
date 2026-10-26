# OPALX Catalyst in-situ visualization

ParaView/Catalyst state scripts and helper modules for OPALX in-situ
visualization. The state scripts in this directory are loaded by the
simulation through IPPL's Catalyst adaptor (`catalyst_execute` entry point).

## Layout

- `catalyst-muon-test-1.py` — single-bunch MUON test pipeline (particles +
  geometry + transform filter)
- `catalyst-muon-fields.py` — self-field visualization state script
- `catalyst-muon-combined.py` — combined particles + fields state script
- `paraview/bunch_transform_filter.py` — programmable filter transforming
  geometry/particles between the bunch reference frame and the lab frame
- `paraview/*.pvsm` — ParaView state files for offline reproduction
- `docs/`, `*.md` — design notes and review documents

## Environment variables

| Variable | Default | Effect |
|---|---|---|
| `BUNCH_TRANSFORM_INVERSE` | `1` | `0` transforms lab -> bunch, `1` bunch -> lab |
| `BUNCH_TRANSFORM_FILTER` | repo path | Path to `bunch_transform_filter.py`, used by state-script loaders |
| `CATALYST_SNAPSHOT_SCRIPT` | repo path | Path to `snapshot_extractor.py`, used by state-script loaders |
| `BUNCH_TRANSFORM_DEBUG` | unset | `1` appends per-execution diagnostics to `/tmp/bunch_transform_filter_debug.txt` |
| `OPALX_CATALYST_DUMP_SNAPSHOT` | `0` | `N` creates a VTPC snapshot extractor writing on iteration `N` and every `N`th after |
| `OPALX_CATALYST_LIVE_FREQ` | `1` | Throttles Catalyst live updates to every `N`th tracking step |
| `OPALX_CATALYST_FIELDS` | unset | `ON` exports field channels on the bunch mesh (see `ParallelTracker.cpp`) |

Set `BUNCH_TRANSFORM_DEBUG=1` before a run to append filter diagnostics
(dataset traversal, transform origins, leaf point counts) to
`/tmp/bunch_transform_filter_debug.txt`.
