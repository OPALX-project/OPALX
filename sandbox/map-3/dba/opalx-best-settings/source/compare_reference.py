#!/usr/bin/env python3
"""Compare saved OPALX maps with independent field and thin-edge references.

No OPALX runs. The CAS map is a nominal-design thin-edge approximation, not
an exact solution through the distributed Enge fields. Its difference includes
model/reference-orbit differences and must not be called integration error.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from variational_reference import ROOT, REPO, load_parameters, np, pd, particle_constants
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

CAS_SOURCE = "https://indico.cern.ch/event/1462429/contributions/6157280/attachments/2937404/5159798/Numerical_methods_tracking-2.pdf"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def analytic_thin_edge(p):
    """CAS slide 29, E1=E2=0, full gap g=2*HGAP; nominal hard-edge body.

    Coordinates (x,x',y,y',zeta,delta); entrance/exit fringe support lengths
    are replaced by drifts. Thus body positions and total design length stay
    unchanged. This model uses the nominal field-normalization rigidity.
    """
    u0, _ = particle_constants()
    gamma2 = 1 + u0*u0

    def drift(length):
        m = np.eye(6)
        m[0, 1] = m[2, 3] = length
        m[4, 5] = length / gamma2
        return m

    def plane(k, length):
        if k == 0:
            return np.array([[1., length], [0., 1.]])
        r = math.sqrt(abs(k))
        c, s = (math.cos(r*length), math.sin(r*length)) if k > 0 else (math.cosh(r*length), math.sinh(r*length))
        return np.array([[c, s/r], [-math.copysign(r, k)*s, c]])

    rho, angle = p.radius_m, math.radians(p.body_angle_deg)
    c, s = math.cos(angle), math.sin(angle)
    bend = np.eye(6)
    bend[:2, :2] = [[c, rho*s], [-s/rho, c]]
    bend[:2, 5] = [rho*(1-c), s]
    bend[2, 3] = p.body_length
    bend[4, :2] = [-s, -rho*(1-c)]
    bend[4, 5] = rho*(s-angle) + p.body_length/gamma2
    edge = np.eye(6)
    edge[3, 2] = math.tan(p.gap * p.fint / rho) / rho
    quad = np.eye(6)
    quad[:2, :2] = plane(p.charge_e * p.quadrupole_k1_seed_m2, p.quadrupole_length_m)
    quad[2:4, 2:4] = plane(-p.charge_e * p.quadrupole_k1_seed_m2, p.quadrupole_length_m)
    quad[4, 5] = p.quadrupole_length_m/gamma2
    result = np.eye(6)
    for m in [drift(p.half_width), edge, bend, edge, drift(p.half_width+p.field_free_drift_m),
              quad, drift(p.half_width+p.field_free_drift_m), edge, bend, edge, drift(p.half_width)]:
        result = m @ result
    return result


def plot_differences(differences, target):
    # Scaling positions by L*=1 m makes all entries dimensionless; numerical
    # values equal SI-unit matrix entries, but this is not a magnet length.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelsize": 11, "axes.titlesize": 12,
                         "savefig.dpi": 180, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 5.0), layout="constrained")
    limit = max(float(abs(d).max()) for d in differences.values())
    norm = SymLogNorm(linthresh=1e-10, vmin=-limit, vmax=limit, base=10)
    cmap = plt.get_cmap("RdBu_r")
    for ax, (name, diff) in zip(axes, differences.items()):
        im = ax.imshow(diff, cmap=cmap, norm=norm)
        ax.set(xticks=range(6), yticks=range(6), xticklabels=range(1, 7), yticklabels=range(1, 7),
               xlabel="Input coordinate j", ylabel="Output coordinate i", title=name)
        ax.tick_params(length=0)
        for i, j in np.ndindex(6, 6):
            value = diff[i, j]
            label = "0" if value == 0 else ("<1e-11" if abs(value) < 1e-11 else f"{value:+.0e}".replace("e-0", "e-"))
            rgba = cmap(norm(value))
            luminance = .2126*rgba[0] + .7152*rgba[1] + .0722*rgba[2]
            ax.text(j, i, label, ha="center", va="center", fontsize=8.2,
                    color="white" if luminance < .5 else "#111111")
    bar = fig.colorbar(im, ax=axes, location="bottom", shrink=.65, pad=.04, aspect=50)
    bar.set_label("Signed dimensionless difference (symmetric-log colour scale)")
    fig.suptitle("map-3: OPALX minus independent distributed-field reference", fontsize=15)
    fig.supxlabel("Coordinates: (x, x′, y, y′, ζ, δ); positions scaled by 1 m. Not the CAS thin-edge approximation.", fontsize=10)
    fig.savefig(target / "map-differences.png")
    fig.savefig(target / "map-differences.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=ROOT / "variational-reference-validated")
    parser.add_argument("--output", type=Path, default=ROOT / "map-comparison")
    args = parser.parse_args()
    provenance = json.loads((args.reference / "provenance.json").read_text())
    for name, sha in provenance["sources_sha256"].items():
        assert digest(REPO / name) == sha, name
    reference = np.loadtxt(args.reference / "reference-matrix.txt")
    np.testing.assert_array_equal(reference, np.loadtxt(args.reference / "dop853-finer-matrix.txt"))
    uncertainty = max(r["max_entry_difference_from_finer"] for r in provenance["refinements"]
                      if r["name"] in ("dop853-fine", "rk45-crosscheck"))
    assert uncertainty < 1e-9, uncertainty
    args.output.mkdir(parents=True, exist_ok=True)
    cas = analytic_thin_edge(load_parameters())
    np.savetxt(args.output / "cas-thin-edge-matrix.txt", cas, fmt="%.16e",
               header="CAS thin-edge approximation on nominal design orbit; NOT the distributed-field reference")
    rows, differences, sources = [], {}, {}
    for name in ("BORIS", "RK4", "DOP853"):
        case = ROOT / "opalx-best-settings" / name.lower()
        record = json.loads((case / "result.json").read_text())
        assert record["status"] == "OK" and record["ranks"] == 1
        for artifact, sha in record["artifacts_sha256"].items():
            assert digest(case / artifact) == sha, (name, artifact)
        matrix = np.loadtxt(case / "matrix.txt")
        sources[str((case / "matrix.txt").relative_to(REPO))] = digest(case / "matrix.txt")
        diff = matrix-reference
        differences[name] = diff
        i, j = np.unravel_index(abs(diff).argmax(), diff.shape)
        np.savetxt(args.output / f"{name.lower()}-minus-field-reference.txt", diff, fmt="%+.16e")
        np.savetxt(args.output / f"{name.lower()}-minus-cas-thin-edge.txt", matrix-cas, fmt="%+.16e")
        rows.append(dict(integrator=name, signed_determinant_error_from_saved_matrix=float(np.linalg.det(matrix)-1),
                         absolute_determinant_error_stdout=record["determinant_error"],
                         maximum_scaled_entry_error=float(abs(diff).max()), largest_entry=f"R{i+1}{j+1}",
                         delta_R16_m=diff[0, 5], delta_R26=diff[1, 5],
                         cas_model_difference_max=float(abs(matrix-cas).max()),
                         cas_delta_R16_m=matrix[0, 5]-cas[0, 5], cas_delta_R26=matrix[1, 5]-cas[1, 5]))
    table = pd.DataFrame(rows)
    table.to_csv(args.output / "differences.csv", index=False, float_format="%.16e")
    plot_differences(differences, args.output)
    source_paths = [Path(__file__), args.reference / "reference-matrix.txt", args.reference / "provenance.json"]
    sources.update({str(path.relative_to(REPO)): digest(path) for path in source_paths})
    (args.output / "provenance.json").write_text(json.dumps(dict(
        sources_sha256=sources, numerical_reference_refinement_difference=uncertainty,
        reference_is_closed_form=False, position_scale_m=1.0, cas_source=CAS_SOURCE,
        cas_slide=29, cas_is_distributed_field_solution=False), indent=2)+"\n")
    print(table.to_string(index=False))
    print("Reference refinement / cross-method difference:", uncertainty)
    print("CAS thin-edge map (different model):\n", cas)


if __name__ == "__main__":
    main()
