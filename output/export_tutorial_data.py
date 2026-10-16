"""Export checked, unrounded local measurements for the offline beginner tutorial.

Run: .venv/Scripts/python.exe output/export_tutorial_data.py
The exporter checks every sample hash and recomputes summaries before writing.
"""

import csv
import json
from pathlib import Path
import sys

import numpy as np
import openmm as mm
from openmm import unit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import analyze
from verify_project import digest, expected_cases, require, same_value, verify_analysis


def main():
    protocol_path = ROOT / "protocol.json"
    analysis_dir = ROOT / "results/analysis"
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    with (analysis_dir / "group_summaries.csv").open(encoding="utf-8", newline="") as stream:
        groups = [{k: v if k in ("ion", "model") else int(v) if k == "replicas" else float(v)
                   for k, v in row.items()} for row in csv.DictReader(stream)]
    payload = {"protocol": protocol, "groups": groups, "radii": analyze.RADII.tolist(),
               "time_ps": None, "hydration": [], "amide": []}
    summaries, metadata = {}, []
    sources = [protocol_path, ROOT / "analyze.py", ROOT / "verify_project.py", Path(__file__),
               analysis_dir / "replica_summaries.json", analysis_dir / "replica_summaries.csv",
               analysis_dir / "group_summaries.csv"]
    for campaign in ("hydration", "amide"):
        for name in sorted(expected_cases(campaign)):
            directory = ROOT / "results" / campaign / name
            meta_path, samples_path = directory / "metadata.json", directory / "samples.npz"
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            case = analyze.read_case(meta_path)
            summary = case["summary"]
            summaries[campaign, name] = summary
            metadata.append(meta)
            sources.extend((meta_path, samples_path))
            if payload["time_ps"] is None:
                payload["time_ps"] = case["time_ps"].tolist()
            require(np.array_equal(case["time_ps"], payload["time_ps"]), "Cases have different time axes")
            entry = {"id": name, "summary": summary, "cn": case["CN"].tolist(),
                     "run": {k: meta[k] for k in ("seed", "n_water", "platform", "min_box_length_nm", "samples_sha256")},
                     "mean_cation_molar": float(np.mean(case["density"] / summary["water_count"]) / 0.602214076),
                     "temperature_K": case["temperature_K"].tolist(),
                     "water_number_density_nm3": case["density"].tolist()}
            if campaign == "hydration":
                with np.load(samples_path, allow_pickle=False) as archive:
                    distances = archive["distances_nm"].astype(float)
                snapshots = []
                for frame in (0, 499, 999, 1499, 1999, 166, 167):
                    ids = np.flatnonzero(distances[frame] < 0.60)
                    snapshots.append({"time_ps": float(case["time_ps"][frame]),
                                      "waters": [[int(i), float(distances[frame, i])] for i in ids]})
                cutoff = summary["shell_cutoff_nm"]
                scan = [{"cutoff_nm": float(c), "mean_CN": float((distances < c).sum(axis=1).mean())}
                        for c in np.linspace(cutoff - .02, cutoff + .02, 9)]
                entry.update(rdf=case["rdf"].tolist(), cumulative_CN=case["cumulative_CN"].tolist(),
                             snapshots=snapshots, cutoff_scan=scan)
            else:
                entry["distance"] = case["amide_distance_nm"].tolist()
            payload[campaign].append(entry)
    audit = verify_analysis(analysis_dir, summaries)
    saved_groups = json.loads((analysis_dir / "replica_summaries.json").read_text(encoding="utf-8"))["hydration_group_summaries"]
    require(same_value(groups, saved_groups), "Group CSV differs from audited group JSON")
    production_ns = sum(m["production_ps"] for m in metadata) / 1000
    production_wall = sum(m["production_wall_seconds"] for m in metadata)
    payload["timing"] = {
        "cases": len(metadata), "production_ns": production_ns,
        "production_wall_seconds": production_wall,
        "case_wall_seconds": sum(m["wall_seconds"] for m in metadata),
        "production_ns_per_day": production_ns / production_wall * 86400,
        "min_ns_per_day": min(m["ns_per_day"] for m in metadata),
        "max_ns_per_day": max(m["ns_per_day"] for m in metadata),
        "devices": sorted({m["platform_properties"]["DeviceName"] for m in metadata}),
    }
    payload["provenance"] = {
        "audit": audit,
        "sources": [{"path": p.relative_to(ROOT).as_posix(), "sha256": digest(p)} for p in sources],
        "notes": ["Formal hydration and amide cases only; pilots excluded.",
                  "All series retain every saved production frame, at 1 ps cadence; no rounding or subsampling.",
                  "Snapshot water IDs are zero-based original distances_nm column IDs, not atom indices.",
                  "Snapshots show distances below 0.60 nm only; no coordinates or angles are exported; 167/168 ps identity example is included.",
                  "Hydration cutoff scans use all waters and all 2000 frames, not the five snapshots.",
                  "RDFs and summaries are recomputed with analyze.read_case and checked against existing JSON/CSV.",
                  "The same integer water ID across different cases does not imply the same physical water.",
                  "Timing reflects concurrent desktop/GPU use, not an isolated hardware benchmark."],
    }
    system_path = ROOT / "results/hydration/CM_Mg_r1_hydrated/system.xml"
    first_meta = json.loads(system_path.with_name("metadata.json").read_text())
    system = mm.XmlSerializer.deserialize(system_path.read_text())
    nonbonded = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    q, sigma, epsilon = nonbonded.getParticleParameters(first_meta["water_oxygen_indices"][0])
    payload["water_oxygen"] = {"charge_e": q.value_in_unit(unit.elementary_charge),
                               "sigma_nm": sigma.value_in_unit(unit.nanometer),
                               "epsilon_kj_mol": epsilon.value_in_unit(unit.kilojoule_per_mole)}
    payload["conversion"] = json.loads((ROOT / "parameters/freesolv/conversion_verification.json").read_text())
    for p in (system_path, ROOT / "parameters/freesolv/conversion_verification.json", ROOT / "parameters/source_manifest.json"):
        payload["provenance"]["sources"].append({"path": p.relative_to(ROOT).as_posix(), "sha256": digest(p)})
    destination = Path(__file__).with_name("tutorial-data.json")
    destination.write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(",", ":")), encoding="utf-8")
    print(f"Checked 36 cases; wrote {destination} ({destination.stat().st_size:,} bytes).")


if __name__ == "__main__":
    main()
