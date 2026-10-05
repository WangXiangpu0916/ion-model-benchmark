"""Analyze saved ion-water distances; never infer exchange rates from saved frames.

Usage: python analyze.py --input results/hydration --output results/analysis
       python analyze.py --input results/hydration --amide results/amide --output results/analysis
       python analyze.py --self-check

Distances must already use minimum-image geometry; water columns must retain their
identities. RDF range requires a shortest periodic box edge greater than 1.8 nm.
The fixed shell cutoffs are diagnostics, not exact experimental shell boundaries.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


EDGES = np.arange(451, dtype=float) * 0.002
RADII = (EDGES[:-1] + EDGES[1:]) / 2
SHELL = {"Mg": 0.30, "Ca": 0.35}
PEAK_WINDOW = {"Mg": (0.17, 0.26), "Ca": (0.20, 0.32)}
TARGET_DISTANCE = {"Mg": 0.209, "Ca": 0.246}
TARGET_CN = {"Mg": 6.0, "Ca": 8.0}
COLORS = {"CM": "#2878A6", "IOD": "#D36E2D"}
MODEL_ORDER = ("CM", "IOD")
ION_ORDER = ("Mg", "Ca")
NOTES = [
    "Observed trajectories, not independent frames, are replication units; SD is descriptive between-replica variation.",
    "Fitting targets were used to develop these models; agreement is an implementation/transfer check, not independent validation.",
    "Mg 0.209 +/- 0.004 nm is a supplied fitting-reference range, not a simulation confidence interval; Ca reference is 0.246 nm.",
    "Coordination cutoffs Mg 0.30 nm and Ca 0.35 nm, with +/-0.02 nm sensitivity, are diagnostic definitions.",
    "RDFs average frame histograms normalized by that frame's water number density; cumulative CN averages unweighted counts.",
    "Mean shell distance pools all in-shell distances; it is different from the RDF peak position.",
    "No automatic equilibration trimming or convergence pass is applied; input must contain the designated production segment.",
    "Shell-identity changes are observations between saved frames, not committed exchange events or kinetic rates.",
    "Minimum-image distance calculation and invariant water-column identities are responsibilities of the simulation writer.",
    "RDF use requires shortest periodic box edge >1.8 nm; volume alone cannot verify this condition for a general box.",
    "An amide initial-state difference indicates unresolved sampling; these traces do not estimate affinity or NMDAR function.",
]


def radial_distribution(distances, volumes, edges=EDGES):
    """Return frame-density-normalized g(r) and mean cumulative water counts."""
    nframes, nwater = distances.shape
    weighted = np.zeros(len(edges) - 1)
    counts = np.zeros_like(weighted)
    # Bound temporary histogram allocations for large distance arrays.
    for start in range(0, nframes, 128):
        part = distances[start:start + 128]
        indices = np.searchsorted(edges, part, side="right") - 1
        valid = (indices >= 0) & (indices < len(weighted))
        inverse_density = np.broadcast_to(volumes[start:start + len(part), None] / nwater, part.shape)
        weighted += np.bincount(indices[valid], weights=inverse_density[valid], minlength=len(weighted))
        counts += np.bincount(indices[valid], minlength=len(weighted))
    shell_volumes = 4 * np.pi / 3 * np.diff(edges ** 3)
    return weighted / nframes / shell_volumes, np.cumsum(counts) / nframes


def membership_changes(mask):
    changes = np.count_nonzero(mask[1:] != mask[:-1], axis=1)
    return {
        "shell_identity_changed_intervals": int(np.count_nonzero(changes)),
        "shell_identity_changed_interval_fraction": float(np.mean(changes > 0)),
        "total_membership_flips": int(changes.sum()),
        "distinct_waters_seen_in_shell": int(np.any(mask, axis=0).sum()),
        "identity_observation_intervals": len(changes),
    }


def conditional_distance(distances, cutoff):
    selected = distances[distances < cutoff]
    return float(selected.mean()) if selected.size else None


def read_case(metadata_path):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("status") != "complete":
        raise ValueError(f"{metadata_path}: case is not complete")
    required = ("campaign", "model", "ion", "replica", "initial_state", "temperature_K",
                "production_ps", "sample_interval_ps", "seed")
    missing = [key for key in required if key not in metadata]
    if missing:
        raise ValueError(f"{metadata_path}: missing metadata {missing}")
    if metadata["model"] not in MODEL_ORDER or metadata["ion"] not in ION_ORDER:
        raise ValueError(f"{metadata_path}: expected model CM/IOD and ion Mg/Ca")
    if metadata["campaign"] not in ("hydration", "amide", "pilot"):
        raise ValueError(f"{metadata_path}: invalid campaign")
    if type(metadata["replica"]) is not int or metadata["replica"] not in (1, 2, 3):
        raise ValueError(f"{metadata_path}: replica must be integer 1, 2, or 3")
    if metadata["initial_state"] not in ("hydrated", "contact", "separated"):
        raise ValueError(f"{metadata_path}: invalid initial_state")
    for key in ("temperature_K", "production_ps", "sample_interval_ps"):
        if not isinstance(metadata[key], (int, float)) or not np.isfinite(metadata[key]) or metadata[key] <= 0:
            raise ValueError(f"{metadata_path}: {key} must be finite and positive")
    samples_path = metadata_path.with_name("samples.npz")
    if hashlib.sha256(samples_path.read_bytes()).hexdigest() != metadata.get("samples_sha256"):
        raise ValueError(f"{samples_path}: checksum differs from the completed run manifest")
    with np.load(samples_path, allow_pickle=False) as archive:
        arrays = {key: np.asarray(archive[key], dtype=float) for key in (
            "distances_nm", "time_ps", "volume_nm3", "temperature_K", "potential_kj_mol",
            "min_chloride_nm", "amide_distance_nm")}
    distances = arrays["distances_nm"]
    if distances.ndim != 2 or distances.shape[0] < 4 or distances.shape[1] < 1:
        raise ValueError(f"{samples_path}: distances must have >=4 frames and >=1 water")
    if not np.all(np.isfinite(distances)) or np.any(distances < 0):
        raise ValueError(f"{samples_path}: distances must be finite and nonnegative")
    frames = len(distances)
    expected_frames = round(metadata["production_ps"] / metadata["sample_interval_ps"])
    if frames != expected_frames:
        raise ValueError(f"{samples_path}: expected {expected_frames} frames, found {frames}")
    if metadata.get("min_box_length_nm") is not None and metadata["min_box_length_nm"] <= 2*EDGES[-1]:
        raise ValueError(f"{samples_path}: box too short for requested RDF radius")
    for key, value in arrays.items():
        if key == "distances_nm":
            continue
        if value.shape != (frames,):
            raise ValueError(f"{samples_path}: {key} must have one entry per frame")
        if key not in ("min_chloride_nm", "amide_distance_nm") and not np.all(np.isfinite(value)):
            raise ValueError(f"{samples_path}: {key} contains nonfinite values")
    if np.any(arrays["volume_nm3"] <= 0) or np.any(arrays["temperature_K"] <= 0):
        raise ValueError(f"{samples_path}: volume and temperature must be positive")
    cadence = np.diff(arrays["time_ps"])
    if not np.isclose(arrays["time_ps"][0], metadata["sample_interval_ps"]) or not np.isclose(arrays["time_ps"][-1], metadata["production_ps"]):
        raise ValueError(f"{samples_path}: saved times do not cover the planned production segment")
    if np.any(cadence <= 0) or not np.allclose(cadence, metadata["sample_interval_ps"], rtol=1e-4, atol=1e-6):
        raise ValueError(f"{samples_path}: nonincreasing time or cadence differs from metadata")
    for key in ("min_chloride_nm", "amide_distance_nm"):
        values = arrays[key]
        if np.any(values[np.isfinite(values)] < 0) or np.any(np.isneginf(values)):
            raise ValueError(f"{samples_path}: invalid {key}")
    if metadata["campaign"] == "amide" and not np.all(np.isfinite(arrays["amide_distance_nm"])):
        raise ValueError(f"{samples_path}: amide campaign requires finite amide distances")
    if metadata["campaign"] != "amide" and not np.all(np.isnan(arrays["amide_distance_nm"])):
        raise ValueError(f"{samples_path}: hydration/pilot amide distances must all be NaN")
    case = analyze_arrays(metadata, arrays)
    case["summary"]["case"] = metadata_path.parent.name
    case["summary"]["source"] = str(metadata_path.resolve())
    return case


def analyze_arrays(metadata, arrays):
    distances = arrays["distances_nm"]
    frames, waters = distances.shape
    ion = metadata["ion"]
    cutoff = SHELL[ion]
    membership = distances < cutoff
    cn = membership.sum(axis=1)
    density = waters / arrays["volume_nm3"]
    rdf, cumulative_cn = radial_distribution(distances, arrays["volume_nm3"])
    peak_indices = np.flatnonzero((RADII >= PEAK_WINDOW[ion][0]) & (RADII <= PEAK_WINDOW[ion][1]))
    peak_index = peak_indices[np.argmax(rdf[peak_indices])]
    valid_peak = bool(np.max(rdf[peak_indices]) > 0)
    cn_values, cn_counts = np.unique(cn, return_counts=True)
    chloride = arrays["min_chloride_nm"]
    known_chloride = ~np.isnan(chloride)
    summary = {key: metadata[key] for key in ("campaign", "model", "ion", "replica", "initial_state",
               "production_ps", "sample_interval_ps", "seed")}
    summary.update({
        "target_temperature_K": metadata["temperature_K"],
        "frames": frames, "water_count": waters,
        "first_saved_time_ps": float(arrays["time_ps"][0]),
        "last_saved_time_ps": float(arrays["time_ps"][-1]),
        "shell_cutoff_nm": cutoff,
        "rdf_peak_nm": float(RADII[peak_index]) if valid_peak else None,
        "rdf_peak_at_search_boundary": bool(peak_index in (peak_indices[0], peak_indices[-1])) if valid_peak else None,
        "mean_shell_distance_nm": conditional_distance(distances, cutoff),
        "mean_CN": float(cn.mean()),
        "mean_temperature_K": float(arrays["temperature_K"].mean()),
        "mean_water_number_density_nm3": float(density.mean()),
        "mean_potential_kj_mol": float(arrays["potential_kj_mol"].mean()),
        "chloride_contact_fraction_lt_0p4nm": float(np.mean(chloride[known_chloride] < 0.4)) if np.any(known_chloride) else None,
        "chloride_distance_observed_frames": int(known_chloride.sum()),
        "CN_distribution": {str(int(n)): float(count / frames) for n, count in zip(cn_values, cn_counts)},
        "cutoff_sensitivity": [{"cutoff_nm": float(c), "mean_CN": float(np.mean(np.sum(distances < c, axis=1)))}
                               for c in (cutoff - 0.02, cutoff, cutoff + 0.02)],
        **membership_changes(membership),
    })
    series = {"CN": cn, "temperature_K": arrays["temperature_K"],
              "water_number_density_nm3": density, "potential_kj_mol": arrays["potential_kj_mol"]}
    half = frames // 2
    summary["halves"] = {}
    for name, values in series.items():
        first, second = float(values[:half].mean()), float(values[half:].mean())
        summary["halves"][name] = {"first": first, "second": second, "second_minus_first": second - first}
    summary["halves"]["mean_shell_distance_nm"] = {
        "first": conditional_distance(distances[:half], cutoff),
        "second": conditional_distance(distances[half:], cutoff)}
    summary["block_means"] = []
    for block in np.array_split(np.arange(frames), min(8, frames // 2)):
        summary["block_means"].append({"time_ps": float(arrays["time_ps"][block].mean()),
                                       "frames": len(block),
                                       **{name: float(values[block].mean()) for name, values in series.items()}})
    if metadata["campaign"] == "amide":
        amide = arrays["amide_distance_nm"]
        summary["mean_amide_distance_nm"] = float(amide.mean())
        summary["amide_fraction_within_diagnostic_water_cutoff"] = float(np.mean(amide < cutoff))
        summary["halves"]["amide_distance_nm"] = {
            "first": float(amide[:half].mean()), "second": float(amide[half:].mean())}
    return {"summary": summary, "rdf": rdf, "cumulative_CN": cumulative_cn,
            "time_ps": arrays["time_ps"], "amide_distance_nm": arrays["amide_distance_nm"],
            "CN": cn, "density": density, "temperature_K": arrays["temperature_K"]}


def load_cases(root, amide=False):
    paths = sorted(root.rglob("metadata.json"))
    if not paths:
        raise ValueError(f"No metadata.json files found in {root}")
    cases = [read_case(path) for path in paths]
    seen = set()
    for case in cases:
        s = case["summary"]
        if (s["campaign"] == "amide") != amide:
            raise ValueError(f"Unexpected campaign {s['campaign']} in {root}")
        key = tuple(s[k] for k in ("model", "ion", "replica", "initial_state"))
        if key in seen:
            raise ValueError(f"Duplicate case identity: {key}")
        seen.add(key)
    return cases


def groups(cases):
    for ion in ION_ORDER:
        for model in MODEL_ORDER:
            selected = [c for c in cases if c["summary"]["ion"] == ion and c["summary"]["model"] == model]
            yield ion, model, selected


def mean_sd(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return (float(values.mean()), float(values.std(ddof=1)) if len(values) > 1 else None) if len(values) else (None, None)


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_tables(output, cases, amide_cases):
    summaries = [c["summary"] for c in cases + amide_cases]
    rows = []
    for summary in summaries:
        rows.append({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value
                     for key, value in summary.items()})
    write_csv(output / "replica_summaries.csv", rows)
    aggregate = []
    for ion, model, selected in groups(cases):
        row = {"ion": ion, "model": model, "replicas": len(selected)}
        for metric in ("rdf_peak_nm", "mean_shell_distance_nm", "mean_CN",
                       "shell_identity_changed_interval_fraction", "chloride_contact_fraction_lt_0p4nm"):
            mean, sd = mean_sd([c["summary"][metric] for c in selected if c["summary"][metric] is not None])
            row[metric + "_mean"] = mean
            row[metric + "_SD"] = sd
        aggregate.append(row)
    payload = {"analysis_notes": NOTES, "replicas": summaries, "hydration_group_summaries": aggregate}
    (output / "replica_summaries.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_csv(output / "group_summaries.csv", aggregate)
    rdf_rows = []
    for case in cases + amide_cases:
        s = case["summary"]
        for radius, upper, g, cumulative in zip(RADII, EDGES[1:], case["rdf"], case["cumulative_CN"]):
            rdf_rows.append({**{key: s[key] for key in ("campaign", "case", "model", "ion", "replica", "initial_state")},
                             "radius_nm": float(radius), "bin_upper_nm": float(upper),
                             "g_M_waterO": float(g), "mean_CN_through_bin_upper": float(cumulative)})
    write_csv(output / "rdf_curves.csv", rdf_rows)
    (output / "analysis_notes.txt").write_text("\n".join(NOTES) + "\n", encoding="utf-8")


def make_figures(output, cases, amide_cases):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 180, "pdf.fonttype": 42})

    def finish(fig, name, note, top=0.94):
        fig.text(0.02, 0.01, note, fontsize=8, va="bottom")
        fig.tight_layout(rect=(0, 0.06, 1, top))
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"{name}.{suffix}")
        plt.close(fig)

    def scatter_replicas(ax, selected, metric, x, color, multiplier=1):
        values = [c["summary"][metric] for c in selected if c["summary"][metric] is not None]
        if not values:
            return
        values = np.asarray(values) * multiplier
        jitter = np.linspace(-0.10, 0.10, len(values)) if len(values) > 1 else [0]
        ax.scatter(x + np.asarray(jitter), values, color=color, s=25, alpha=0.85, zorder=3)
        mean, sd = mean_sd(values)
        ax.errorbar(x, mean, yerr=sd, color="black", fmt="_", ms=13, capsize=4, zorder=4)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.7))
    fig.suptitle("Hydration RDFs: individual replicas and equally weighted means")
    for ax, ion in zip(axes, ION_ORDER):
        for _, model, selected in (g for g in groups(cases) if g[0] == ion):
            for case in selected:
                ax.plot(RADII, case["rdf"], color=COLORS[model], alpha=0.25, lw=0.8)
            if selected:
                ax.plot(RADII, np.mean([c["rdf"] for c in selected], axis=0), color=COLORS[model],
                        lw=1.8, label=f"{model}, n={len(selected)}")
        ax.set(title=ion, xlabel="Metal-water O distance (nm)", ylabel="g(r)", xlim=(0.12, 0.65))
        ax.axvline(SHELL[ion], color="gray", ls=":", lw=1, label="Diagnostic shell cutoff")
        ax.legend(fontsize=8)
    finish(fig, "figure1_rdf", "Each frame is normalized by its own water density. Curves beyond 0.65 nm remain in rdf_curves.csv.")

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    fig.suptitle("Hydration metrics versus model-fitting references")
    for row, ion in enumerate(ION_ORDER):
        for col, metric in enumerate(("rdf_peak_nm", "mean_CN")):
            ax = axes[row, col]
            for x, model in enumerate(MODEL_ORDER):
                selected = [c for c in cases if c["summary"]["ion"] == ion and c["summary"]["model"] == model]
                scatter_replicas(ax, selected, metric, x, COLORS[model])
            target = TARGET_DISTANCE[ion] if col == 0 else TARGET_CN[ion]
            ax.axhline(target, color="#656565", ls="--", lw=1, label="Fitting reference")
            if col == 0 and ion == "Mg":
                ax.axhspan(target - 0.004, target + 0.004, color="gray", alpha=0.15, label="Reference range")
            ax.set(title=ion, xticks=[0, 1], xticklabels=MODEL_ORDER, xlim=(-0.4, 1.4),
                   ylabel="RDF peak (nm)" if col == 0 else "Mean water coordination")
            ax.legend(fontsize=8)
    finish(fig, "figure2_hydration_metrics", "Colored points: replicas. Black marks: mean +/- SD, not confidence intervals. Agreement is not independent validation.")

    fig, axes = plt.subplots(2, 3, figsize=(9, 6.1))
    fig.suptitle("Production diagnostics: consecutive block means", fontsize=11)
    metrics = (("CN", "Water coordination", "coordination"),
               ("water_number_density_nm3", "Water O density (nm$^{-3}$)", "water density"),
               ("temperature_K", "Temperature (K)", "temperature"))
    styles = {1: "-", 2: "--", 3: ":"}
    diagnostic_legend = {}
    for row, ion in enumerate(ION_ORDER):
        for col, (metric, label, title) in enumerate(metrics):
            ax = axes[row, col]
            for case in cases:
                s = case["summary"]
                if s["ion"] != ion:
                    continue
                blocks = s["block_means"]
                ax.plot([b["time_ps"] for b in blocks], [b[metric] for b in blocks],
                        color=COLORS[s["model"]], ls=styles[s["replica"]], lw=1.2,
                        label=f"{s['model']} r{s['replica']}")
            ax.set_title(f"{ion}: {title}", fontsize=11)
            ax.set_xlabel("Production time (ps)", fontsize=10)
            ax.set_ylabel(label, fontsize=10)
            ax.tick_params(labelsize=10)
            handles, labels = ax.get_legend_handles_labels()
            diagnostic_legend.update(zip(labels, handles))
    fig.legend(diagnostic_legend.values(), diagnostic_legend.keys(), loc="upper center",
               bbox_to_anchor=(0.5, 0.94), ncol=6, fontsize=9, frameon=False)
    finish(fig, "figure3_sampling_diagnostics",
           "Blocks are correlated diagnostics, not independent replicas. No automatic trimming.\n"
           "First/second-half changes are in replica_summaries.json.", top=0.86)

    fig, axes = plt.subplots(2, 2, figsize=(9, 6.2))
    fig.suptitle("Shell sensitivity and " + ("amide initial-state dependence" if amide_cases else "observed shell-identity changes"), fontsize=11)
    amide_legend = {}
    for row, ion in enumerate(ION_ORDER):
        ax = axes[row, 0]
        for _, model, selected in (g for g in groups(cases) if g[0] == ion):
            if not selected:
                continue
            curves = np.array([[v["mean_CN"] for v in c["summary"]["cutoff_sensitivity"]] for c in selected])
            cutoffs = [v["cutoff_nm"] for v in selected[0]["summary"]["cutoff_sensitivity"]]
            for curve in curves:
                ax.plot(cutoffs, curve, color=COLORS[model], alpha=0.25, lw=0.8)
            spread = curves.std(axis=0, ddof=1) if len(curves) > 1 else None
            ax.errorbar(cutoffs, curves.mean(axis=0), yerr=spread, color=COLORS[model], fmt="o-", capsize=3, label=model)
        ax.set(title=f"{ion}: shell-cutoff sensitivity", xlabel="Cutoff (nm)", ylabel="Mean water coordination")
        ax.set_ylim(min(ax.get_ylim()[0], TARGET_CN[ion] - .01), max(ax.get_ylim()[1], TARGET_CN[ion] + .01))
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.legend(fontsize=9)
        ax = axes[row, 1]
        if amide_cases:
            for case in amide_cases:
                s = case["summary"]
                if s["ion"] != ion:
                    continue
                indices = np.linspace(0, len(case["time_ps"]) - 1, min(2000, len(case["time_ps"])), dtype=int)
                group_label = f"{s['model']} {s['initial_state']}"
                line, = ax.plot(case["time_ps"][indices], case["amide_distance_nm"][indices], color=COLORS[s["model"]],
                                ls="-" if s["initial_state"] == "contact" else "--", alpha=0.65, lw=0.8)
                amide_legend.setdefault(group_label, line)
            ax.set(title=f"{ion}: carbonyl-O distance", xlabel="Production time (ps)", ylabel="Metal-amide O (nm)")
        else:
            for x, model in enumerate(MODEL_ORDER):
                selected = [c for c in cases if c["summary"]["ion"] == ion and c["summary"]["model"] == model]
                scatter_replicas(ax, selected, "shell_identity_changed_interval_fraction", x, COLORS[model], 100)
            ax.set(title=f"{ion}: adjacent saved-frame comparison", xticks=[0, 1], xticklabels=MODEL_ORDER,
                   xlim=(-0.4, 1.4), ylabel="Intervals with membership changes (%)")
    for ax in axes.flat:
        ax.title.set_fontsize(11)
        ax.xaxis.label.set_fontsize(10)
        ax.yaxis.label.set_fontsize(10)
        ax.tick_params(labelsize=10)
    if amide_cases:
        fig.legend(amide_legend.values(), amide_legend.keys(), loc="upper center",
                   bbox_to_anchor=(0.5, 0.94), ncol=4, fontsize=9, frameon=False)
    note = ("All replica curves are retained: complete design has 3 curves per model/start/ion.\n"
            "Initial-state dependence signals unresolved sampling; contact is diagnostic, not affinity."
            if amide_cases else "Identity changes depend on frame cadence and shell cutoff; they are not committed exchanges or kinetic rates. Error bars: SD.")
    finish(fig, "figure4_sensitivity", note, top=0.86 if amide_cases else 0.94)


def self_check():
    # Exact two-frame calculation: varying density must be normalized per frame.
    distances = np.array([[0.10, 0.30], [0.10, 0.15]])
    volumes = np.array([1.0, 2.0])
    edges = np.array([0.0, 0.2, 0.4])
    g, cumulative = radial_distribution(distances, volumes, edges)
    shell_volumes = 4 * np.pi / 3 * np.diff(edges ** 3)
    assert np.allclose(g * shell_volumes, [1.25, 0.25])
    assert np.allclose(cumulative, [1.5, 2.0])
    assert not np.allclose(g * shell_volumes, np.array([1.5, 0.5]) / np.mean(2 / volumes))
    changes = membership_changes(distances < 0.2)
    assert changes["total_membership_flips"] == 1 and changes["shell_identity_changed_interval_fraction"] == 1
    # Equal CN can hide a replaced water identity.
    assert membership_changes(np.array([[True, False], [False, True]]))["total_membership_flips"] == 2
    assert np.isclose(conditional_distance(distances, 0.2), 0.35 / 3)
    assert mean_sd([6.0]) == (6.0, None)
    assert mean_sd([5.0, 6.0, 7.0]) == (6.0, 1.0)
    print("Self-check passed: per-frame RDF normalization, CN integral, identity observations and replica SD.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--amide", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    if args.input is None or args.output is None:
        parser.error("--input and --output are required unless --self-check is used")
    try:
        cases = load_cases(args.input)
        amide_cases = load_cases(args.amide, amide=True) if args.amide else []
        for ion, model, selected in groups(cases):
            if len(selected) != 3:
                print(f"NOTE: {ion}/{model} has {len(selected)} replicas; expected 3. Summary remains descriptive.")
        args.output.mkdir(parents=True, exist_ok=True)
        write_tables(args.output, cases, amide_cases)
        make_figures(args.output, cases, amide_cases)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"Analysis failed: {exc}\n")
    print(f"Analyzed {len(cases)} hydration/pilot and {len(amide_cases)} amide cases into {args.output.resolve()}")
    print("Four PNG/PDF figures, replica and group summaries, RDF curves and interpretation notes written.")


if __name__ == "__main__":
    main()
