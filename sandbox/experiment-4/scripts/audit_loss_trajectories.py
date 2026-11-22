#!/usr/bin/env python3
"""Audit immutable H5 witness IDs across losses; never reconnect nearest neighbours.

Positions use H5's straight-line reference translation (RefPartR + local R).
The causal displacement test allows 1 pm for coordinate roundoff, not for
physical motion. Momentum reversals are diagnostics, not a convergence test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

C = 299792458.0
IP_TIME = 0.17 / C


def scalar(group, name):
    return float(np.asarray(group.attrs[name]).reshape(-1)[0])


def audit(path, ballistic_after=None):
    previous = {}
    ever_seen, lost = set(), set()
    times, populations, events, reversals = [], [], [], {}
    max_ratio = 0.0
    violations, duplicates, reappeared, late_births, nonfinite = 0, 0, 0, 0, 0
    worst = None
    ballistic_dp, ballistic_residual = 0.0, 0.0
    ballistic_loss_intervals = 0
    outside_intervals, outside_dp, outside_residual = 0, 0.0, 0.0
    outside_counts = []
    last_time = None
    with h5py.File(path, "r") as h5:
        groups = sorted(h5.values(), key=lambda g: int(g.name.rsplit("#", 1)[1]))
        for group in groups:
            if "id" not in group:
                continue
            t = scalar(group, "TIME")
            ids = np.asarray(group["id"], dtype=np.int64)
            r = np.column_stack([group[k][:] for k in ("x", "y", "z")])
            ref = np.asarray(group.attrs["RefPartR"]).reshape(3)
            r += ref
            p = np.column_stack([group[k][:] for k in ("px", "py", "pz")])
            angles = np.asarray(group.attrs.get("TaitBryantAngles", [0, 0, 0]))
            if np.any(angles != 0):
                raise ValueError("Nonzero reference rotation needs an explicit lab-frame transform")
            nonfinite += int(np.count_nonzero(~np.isfinite(r)))
            nonfinite += int(np.count_nonzero(~np.isfinite(p)))
            current = {int(i): (ri, pi) for i, ri, pi in zip(ids, r, p)}
            duplicates += len(ids) - len(current)
            new_ids = set(current) - ever_seen
            if ever_seen:
                late_births += len(new_ids)  # This deck injects all witnesses at one T0.
            reappeared += len(set(current) & lost)
            removed = set(previous) - set(current)
            loss_ratios = []
            loss_dp, loss_residual = 0.0, 0.0
            for i, (ri, pi) in current.items():
                record = reversals.setdefault(i, {"id": i, "px_reversals": 0, "x_crossings": 0,
                                                  "max_displacement_over_cdt": 0.0})
                if i not in previous:
                    continue
                r0, p0 = previous[i]
                dt = t - last_time
                if dt <= 0:
                    raise ValueError("Non-increasing H5 times")
                displacement = float(np.linalg.norm(ri - r0))
                ratio = displacement / (C * dt)
                # In this straight drift, two samples outside on the same z side
                # must retain momentum under the hard force cutoff. No external
                # elements act in this particular Experiment4 deck.
                same_outside_side = (ri[2] < .16 and r0[2] < .16) or (ri[2] > .18 and r0[2] > .18)
                if same_outside_side:
                    outside_intervals += 1
                    outside_dp = max(outside_dp, float(np.linalg.norm(pi - p0)))
                    velocity = C * p0 / np.sqrt(1 + np.dot(p0, p0))
                    outside_residual = max(outside_residual, float(np.linalg.norm(ri - r0 - velocity * dt)))
                if removed:
                    loss_ratios.append(ratio)
                    if ballistic_after is not None and last_time >= ballistic_after:
                        dp = float(np.linalg.norm(pi - p0))
                        v0 = C * p0 / np.sqrt(1 + np.dot(p0, p0))
                        residual = float(np.linalg.norm(ri - r0 - v0 * dt))
                        loss_dp = max(loss_dp, dp)
                        loss_residual = max(loss_residual, residual)
                record["max_displacement_over_cdt"] = max(record["max_displacement_over_cdt"], ratio)
                if displacement > C * dt + 1.0e-12:
                    violations += 1
                if ratio > max_ratio:
                    max_ratio = ratio
                    worst = {"id": i, "time_s": t, "step": int(scalar(group, "GlobalTrackStep")),
                             "displacement_m": displacement, "dt_s": dt, "ratio": ratio}
                if abs(t - IP_TIME) <= 20e-12:
                    # Exclude exact zeros and tiny roundoff-level normalized momenta.
                    record["px_reversals"] += int(pi[0] * p0[0] < 0 and min(abs(pi[0]), abs(p0[0])) > 1e-12)
                    record["x_crossings"] += int(ri[0] * r0[0] < 0)
            if removed:
                is_ballistic = ballistic_after is not None and last_time >= ballistic_after
                if is_ballistic:
                    ballistic_loss_intervals += 1
                    ballistic_dp = max(ballistic_dp, loss_dp)
                    ballistic_residual = max(ballistic_residual, loss_residual)
                events.append({"time_s": t, "step": int(scalar(group, "GlobalTrackStep")),
                               "removed": len(removed), "survivors": len(current),
                               "ballistic_interval": is_ballistic,
                               "max_ballistic_delta_p_mc": loss_dp if is_ballistic else None,
                               "max_ballistic_position_residual_m": loss_residual if is_ballistic else None,
                               "max_survivor_displacement_over_cdt": max(loss_ratios, default=0.0)})
            lost |= removed
            ever_seen |= set(current)
            times.append(t)
            populations.append(len(current))
            outside_counts.append(int(np.count_nonzero((r[:, 2] < .16) | (r[:, 2] > .18))))
            previous, last_time = current, t
    summary = {"file": str(path), "snapshots": len(times), "initial_count": populations[0],
               "final_count": populations[-1], "unique_ids_seen": len(ever_seen),
               "lost_ids": len(lost), "loss_intervals": len(events), "duplicate_ids": duplicates,
               "reappearing_ids": reappeared, "new_ids_after_initial_injection": late_births,
               "nonfinite_components": nonfinite, "causal_displacement_violations": violations,
               "max_displacement_over_cdt": max_ratio, "worst_interval": worst,
               "last_time_s": times[-1]}
    summary.update(ballistic_after_s=ballistic_after,
                   ballistic_loss_intervals=ballistic_loss_intervals,
                   max_ballistic_delta_p_on_loss_mc=ballistic_dp if ballistic_loss_intervals else None,
                   max_ballistic_position_residual_on_loss_m=ballistic_residual if ballistic_loss_intervals else None,
                   loss_validation_exercised=bool(events))
    summary.update(outside_same_side_intervals=outside_intervals,
                   max_delta_p_outside_window_mc=outside_dp if outside_intervals else None,
                   max_drift_residual_outside_window_m=outside_residual if outside_intervals else None,
                   final_outside_count=outside_counts[-1])
    particles = pd.DataFrame(reversals.values())
    summary["particles_with_px_reversals_near_overlap"] = int((particles.px_reversals > 0).sum())
    summary["max_px_reversals_near_overlap"] = int(particles.px_reversals.max())
    return summary, pd.DataFrame(events), particles, (times, populations)


def plot_tracks(path, ids):
    rows = []
    with h5py.File(path, "r") as h5:
        for group in sorted(h5.values(), key=lambda g: int(g.name.rsplit("#", 1)[1])):
            if "id" not in group:
                continue
            all_ids = group["id"][:]
            selection = np.isin(all_ids, ids)
            ref = np.asarray(group.attrs["RefPartR"]).reshape(3)
            for i, x, px in zip(all_ids[selection], group["x"][:][selection] + ref[0],
                                 group["px"][:][selection]):
                rows.append((int(i), scalar(group, "TIME"), x, px))
    return pd.DataFrame(rows, columns=["id", "time_s", "x_m", "px_mc"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--ranks", type=int, choices=(1, 2), required=True)
    parser.add_argument("--track12-reference", type=Path,
                        help="Optional accepted Track12 pointwise CSV; NOT a matched-physics control")
    args = parser.parse_args()
    run = args.run.resolve()
    output = run / "loss-audit"
    output.mkdir(exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    summaries = {}
    # Lifecycle messages have no time stamp; use the preceding tracker time
    # plus two coarse steps to conservatively exclude the transition interval.
    last_time, ballistic_after, failure = None, None, None
    with (run / "opalx.log").open() as log:
        for line in log:
            match = re.search(r"t=\s*([0-9.eE+-]+)\s*\[ps\]", line)
            if match:
                last_time = float(match[1]) * 1e-12
            if "CIC stencil lies outside" in line:
                failure = "Witness CIC stencil outside solved local field storage"
            if "BB-state=Completed" in line and last_time is not None:
                ballistic_after = last_time + 2e-12
                break
    for c, species in ((1, "electron"), (2, "positron")):
        path, = run.glob(f"*_c{c}.h5")
        summary, events, particles, counts = audit(path, ballistic_after)
        summaries[species] = summary
        events.to_csv(output / f"{species}_loss_intervals.csv", index=False)
        particles.to_csv(output / f"{species}_identity_summary.csv", index=False)
        axes[0, 0].plot((np.asarray(counts[0]) - IP_TIME) * 1e12, counts[1], label=species)
        if len(events):
            axes[0, 1].plot((events.time_s - IP_TIME) * 1e12,
                            events.max_survivor_displacement_over_cdt, ".", label=species, ms=3)
        if c == 2:
            selected = particles.sort_values(["px_reversals", "max_displacement_over_cdt"],
                                             ascending=False).head(8).id.tolist()
            tracks = plot_tracks(path, selected)
            tracks.to_csv(output / "selected_positron_tracks.csv", index=False)
            for i, track in tracks.groupby("id"):
                t = (track.time_s - IP_TIME) * 1e12
                near = abs(t) <= 20
                axes[1, 0].plot(t[near], track.x_m[near] * 1e6, lw=1, label=f"ID {i}")
                axes[1, 1].plot(t, track.x_m * 1e3, lw=1)
    axes[0, 0].set(xlabel="Time from nominal overlap [ps]", ylabel="Witness population")
    axes[0, 1].axhline(1, color="k", ls="--", lw=1, label="causal limit")
    if not any(s["loss_intervals"] for s in summaries.values()):
        axes[0, 1].text(.5, .5, "No loss interval reached\nLoss validation NOT exercised",
                        ha="center", va="center", transform=axes[0, 1].transAxes)
    axes[0, 1].set(xlabel="Time from nominal overlap [ps]", ylabel=r"Max survivor $|\Delta r|/(c\Delta t)$ at losses")
    axes[1, 0].set(xlabel="Time from nominal overlap [ps]", ylabel=r"Positron $x$ [$\mu$m]", title="Interaction region; fixed IDs")
    axes[1, 1].set(xlabel="Time from nominal overlap [ps]", ylabel="Positron x [mm]", title="Same particles: full saved interval")
    for ax in axes.flat:
        ax.grid(alpha=.2)
    for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(run.name + ": immutable-ID audit" + (" — ABORTED before losses" if failure else ""))
    fig.savefig(output / "loss_trajectory_audit.png", dpi=200)
    plt.close(fig)

    # Sparse primary output, where enabled, independently diagnoses source expansion.
    primary, = run.glob("*_c0.h5")
    primary_rows = []
    with h5py.File(primary, "r") as h5:
        for group in h5.values():
            if "x" in group:
                primary_rows.append({"time_s": scalar(group, "TIME"),
                                     "sigma_x_m": float(np.std(group["x"][:])),
                                     "sigma_y_m": float(np.std(group["y"][:]))})
    if primary_rows:
        primary_frame = pd.DataFrame(primary_rows).sort_values("time_s")
        primary_frame.to_csv(output / "primary_envelope.csv", index=False)
        nearest = primary_frame.iloc[(primary_frame.time_s - IP_TIME).abs().argmin()]
        summaries["primary_nearest_nominal_overlap"] = nearest.to_dict()
    if args.track12_reference:
        reference = pd.read_csv(args.track12_reference)
        rows = []
        for pair, track in reference[reference.species == "positron"].groupby("pair"):
            p = track.sort_values("time_s").px_opalx.to_numpy()
            rows.append({"pair": int(pair), "px_reversals": int(np.count_nonzero(p[1:] * p[:-1] < 0))})
        pd.DataFrame(rows).to_csv(output / "reference_track12_positron_reversals.csv", index=False)
        summaries["reference_track12"] = {
            "file": str(args.track12_reference.resolve()), "positron_reversals": rows,
            "caveat": "Different witnesses, timing, source evolution, and mesh; not a controlled equivalence test",
        }
    preparation = json.loads((run / "run_manifest.json").read_text())
    summaries["execution"] = {
        "failure": failure,
        "loss_validation_complete": failure is None and preparation["status"] == "completed"
        and all(summaries[s]["loss_validation_exercised"] for s in ("electron", "positron")),
    }
    summaries["provenance"] = {
        "git_revision": preparation["git_revision"],
        "executable_sha256": preparation["executable_sha256"],
        "mpi_ranks": args.ranks,
        "launch": "mpiexec --oversubscribe --bind-to none -n N OPALX DECK --info 4",
        "analysis_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "deck_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in run.glob("*.in")},
        "clock": "nominal IP/c; exact primary overlap may differ by beta and source offset",
        "scope": "No ID remapping. Causality is necessary, not sufficient, for trajectory accuracy.",
    }
    (output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
