"""Read-only completion audit; missing work exits 1, and --self-check uses temporary files.

Run with the project Python. If pypdf is only in a bundled environment, provide
--pdf-python PATH. Use --data-only before building the report. These checks
establish artifact consistency, not scientific
correctness, equilibration, or the quality of the report's interpretation.
"""

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np
import openmm as mm
from openmm import app, unit

import analyze


ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify_hash(path, expected):
    require(digest(path) == expected, f"SHA-256 mismatch: {path}")


def expected_cases(campaign):
    states = ("hydrated",) if campaign == "hydration" else ("contact", "separated")
    return {f"{model}_{ion}_r{rep}_{state}": (model, ion, rep, state)
            for model, ion, rep, state in itertools.product(("CM", "IOD"), ("Mg", "Ca"), (1, 2, 3), states)}


def verify_sources():
    manifest = json.loads((ROOT / "parameters/source_manifest.json").read_text())
    names = [entry["file"] for entry in manifest["files"]]
    require(len(names) == len(set(names)) >= 7, "Missing or duplicate source manifest entries")
    for entry in manifest["files"]:
        path = ROOT / entry["file"]
        require(path.stat().st_size == entry["bytes"], f"Source size differs: {path}")
        verify_hash(path, entry["sha256"])
        require(entry.get("source_url") and entry.get("source_commit"), f"Source provenance missing: {path}")
    conversion_dir = ROOT / "parameters/freesolv"
    conversion = json.loads((conversion_dir / "conversion_verification.json").read_text())
    require(conversion["passed"] and len(conversion["geometries"]) >= 2, "Acetamide conversion was not verified")
    for name, checksum in conversion["sha256"].items():
        verify_hash(conversion_dir / name, checksum)
    for geometry in conversion["geometries"]:
        for errors in geometry["by_force"].values():
            require(errors["absolute_energy_error_kj_mol"] <= conversion["tolerances"]["absolute_energy_kj_mol"],
                    "Recorded acetamide energy discrepancy exceeds tolerance")
            require(errors["max_absolute_force_error_kj_mol_nm"] <= conversion["tolerances"]["max_absolute_force_kj_mol_nm"],
                    "Recorded acetamide force discrepancy exceeds tolerance")
    return f"{len(names)} upstream files and acetamide conversion hashes/errors checked"


def verify_case(path, campaign, identity, protocol):
    meta = json.loads(path.read_text())
    require(meta.get("status") == "complete", f"status={meta.get('status', 'missing')}")
    require(meta["campaign"] == campaign, "Campaign identity differs")
    require(tuple(meta[k] for k in ("model", "ion", "replica", "initial_state")) == identity, "Case identity differs")
    for key, value in protocol.items():
        if key not in ("status", "notes", "version"):
            require(meta["protocol"][key] == value, f"Numerical protocol differs: {key}")
    require(re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", meta["git_commit"]) is not None, "Invalid recorded Git commit")
    for name, hash_key in (("simulate.py", "code_sha256"), ("protocol.json", "protocol_sha256")):
        content = subprocess.check_output(["git", "show", f"{meta['git_commit']}:{name}"], cwd=ROOT, stderr=subprocess.PIPE)
        require(hashlib.sha256(content).hexdigest() == meta[hash_key], f"{name} does not match recorded Git commit/hash")
    expected_seed = protocol["seeds"][protocol["replicas"].index(meta["replica"])]
    require(meta["seed"] == expected_seed + (50000 if identity[3] == "contact" else 0), "Seed differs from protocol")
    directory = path.parent
    verify_hash(directory / "samples.npz", meta["samples_sha256"])
    duration, cadence = protocol["production_ps"], protocol["sample_interval_ps"]
    frame_count = round(duration / cadence)
    require(np.isclose(frame_count * cadence, duration), "Production duration is not a whole number of sample intervals")
    require(meta["production_ps"] == duration and meta["sample_interval_ps"] == cadence, "Top-level duration/cadence differs")
    require(meta["temperature_K"] == protocol["temperature_K"], "Top-level temperature differs from protocol")
    require(meta["completed_frames"] == frame_count, "Completed frame count differs from protocol")
    require(np.isfinite(meta["min_box_length_nm"]) and meta["min_box_length_nm"] > 2 * protocol["cutoff_nm"],
            "Missing/invalid minimum periodic box length")
    require(meta["box_geometry"] == "cubic with isotropic barostat", "Unexpected box geometry")
    with np.load(directory / "samples.npz", allow_pickle=False) as arrays:
        require(arrays["distances_nm"].shape == (frame_count, meta["n_water"]), "Sample shape differs from frame/water counts")
        require(np.allclose(arrays["time_ps"], np.arange(1, frame_count + 1) * cadence), "Sample time axis differs from protocol")
        require(np.isfinite(arrays["min_chloride_nm"]).all(), "Missing/nonfinite chloride distances")
        require((np.cbrt(arrays["volume_nm3"]) >= meta["min_box_length_nm"] - 1e-8).all(), "Volumes contradict minimum box length")
        last_box = np.cbrt(arrays["volume_nm3"][-1])
    # Reuse the analysis reader for finite-data, units/cadence and shell-metric validation.
    summary = analyze.read_case(path)["summary"]
    system = mm.XmlSerializer.deserialize((directory / "system.xml").read_text())
    topology = app.PDBFile(str(directory / "initial.pdb")).topology
    atoms = list(topology.atoms())
    require(system.getNumParticles() == len(atoms) == meta["n_atoms"], "System/PDB/metadata atom counts differ")
    water_indices = [a.index for a in atoms if a.residue.name == "HOH" and a.element.symbol == "O"]
    ion_indices = [a.index for a in atoms if a.element.symbol == identity[1]]
    chloride_indices = [a.index for a in atoms if a.element.symbol == "Cl"]
    require(water_indices == meta["water_oxygen_indices"] and len(water_indices) == meta["n_water"], "Water identities/count differ")
    require(len(ion_indices) == 1 and len(chloride_indices) == meta["n_chloride"] == 2, "Expected one cation and two chlorides")
    nb = next(force for force in system.getForces() if isinstance(force, mm.NonbondedForce))
    charges = [nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())]
    require(abs(sum(charges)) < 1e-6 and np.isclose(charges[ion_indices[0]], 2), "System is not neutral with a +2 cation")
    require(all(np.isclose(charges[i], -1) for i in chloride_indices), "Chloride charge differs from -1")
    require(nb.getNonbondedMethod() == mm.NonbondedForce.PME and np.isclose(nb.getCutoffDistance().value_in_unit(unit.nanometer), protocol["cutoff_nm"]),
            "PME/cutoff differs from protocol")
    require(np.isclose(nb.getEwaldErrorTolerance(), protocol["ewald_tolerance"]), "PME error tolerance differs")
    _, sigma, epsilon = nb.getParticleParameters(ion_indices[0])
    parameter = protocol["models"][identity[0]][identity[1]]
    require(np.isclose(sigma.value_in_unit(unit.nanometer), .2 * parameter["rmin_half_angstrom"] / 2**(1/6), rtol=1e-10), "Cation sigma differs")
    require(np.isclose(epsilon.value_in_unit(unit.kilojoule_per_mole), 4.184 * parameter["epsilon_kcal_mol"], rtol=1e-10), "Cation epsilon differs")
    require(meta["degrees_of_freedom"] == 3 * system.getNumParticles() - system.getNumConstraints() - 3, "Temperature degrees of freedom differ")
    final = app.PDBFile(str(directory / "final.pdb"))
    correction = meta.get("final_pdb_box_correction")
    if correction:
        verify_hash(directory / "final.pdb", correction["corrected_pdb_sha256"])
        require(correction["source_samples_sha256"] == meta["samples_sha256"], "PDB correction references different samples")
    final_vectors = np.asarray(final.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
    require(np.allclose(final_vectors, np.eye(3) * last_box, atol=6e-5, rtol=0), "Final PDB box differs from final NPT sample")
    require(len(list(final.topology.atoms())) == len(atoms), "Final PDB atom count differs")
    require(np.isfinite(final.positions.value_in_unit(unit.nanometer)).all(), "Final PDB contains nonfinite coordinates")
    return summary


def same_value(value, expected):
    """Compare recomputed fields, including JSON-encoded nested CSV cells."""
    if isinstance(value, str) and not isinstance(expected, str):
        if value == "" and expected is None:
            return True
        if type(expected) is bool:
            return value == str(expected)  # csv.writer writes Python boolean spelling.
        try:
            value = json.loads(value)
        except ValueError:
            return False
    if isinstance(expected, dict):
        return isinstance(value, dict) and value.keys() == expected.keys() and all(same_value(value[k], v) for k, v in expected.items())
    if isinstance(expected, list):
        return isinstance(value, list) and len(value) == len(expected) and all(same_value(v, e) for v, e in zip(value, expected))
    if type(expected) in (int, float):
        return type(value) in (int, float) and np.isfinite(value) and np.isclose(value, expected, rtol=1e-9, atol=1e-12)
    return type(value) is type(expected) and value == expected


def verify_analysis(directory, summaries):
    expected = {(campaign, case) for campaign in ("hydration", "amide") for case in expected_cases(campaign)}
    payload = json.loads((directory / "replica_summaries.json").read_text())
    rows = payload["replicas"]
    keys = [(row["campaign"], row["case"]) for row in rows]
    require(len(keys) == len(set(keys)) and set(keys) == expected, "Analysis JSON must cover exactly 12 hydration and 24 amide cases")
    require(set(summaries) == expected, "Cannot verify analysis values until all raw cases pass")
    with (directory / "replica_summaries.csv").open(newline="", encoding="utf-8") as stream:
        csv_rows = list(csv.DictReader(stream))
    csv_keys = [(row["campaign"], row["case"]) for row in csv_rows]
    require(len(csv_keys) == len(set(csv_keys)) and set(csv_keys) == expected, "Analysis CSV case coverage differs")
    for row in rows + csv_rows:
        key = row["campaign"], row["case"]
        for metric, recomputed in summaries[key].items():
            if metric != "source":  # Absolute source paths can change when a repository is moved.
                require(metric in row and same_value(row[metric], recomputed), f"Analysis {metric} differs from raw data: {key}")
    groups = payload["hydration_group_summaries"]
    require(len(groups) == 4 and {(r["model"], r["ion"]) for r in groups} == set(itertools.product(("CM", "IOD"), ("Mg", "Ca"))),
            "Hydration group coverage differs")
    for group in groups:
        require(group["replicas"] == 3, "Hydration group does not contain three replicas")
        selected = [s for (campaign, _), s in summaries.items() if campaign == "hydration" and s["model"] == group["model"] and s["ion"] == group["ion"]]
        for metric in ("mean_CN", "rdf_peak_nm", "mean_shell_distance_nm"):
            mean, sd = analyze.mean_sd([s[metric] for s in selected])
            require(np.isclose(group[metric + "_mean"], mean) and np.isclose(group[metric + "_SD"], sd), f"Group {metric} differs")
    return "36 exact case identities; JSON/CSV metrics recomputed; four three-replica hydration groups checked"


def verify_figures(directory):
    from PIL import Image
    for name in ("figure1_rdf", "figure2_hydration_metrics", "figure3_sampling_diagnostics", "figure4_sensitivity"):
        with Image.open(directory / f"{name}.png") as image:
            pixels = np.asarray(image.convert("RGB"))
        require(min(pixels.shape[:2]) >= 300 and np.isfinite(pixels).all() and np.ptp(pixels[..., :3]) > .1,
                f"Figure is too small, invalid or blank: {name}")
    return "Four PNG figures decoded; visual interpretation/layout still requires human review"


def verify_report(path, python):
    require(path.is_file(), f"Missing report: {path}")
    code = "from pypdf import PdfReader; import json,sys; r=PdfReader(sys.argv[1]); print(json.dumps({'pages':len(r.pages),'text':'\\n'.join(p.extract_text() or '' for p in r.pages)}))"
    result = subprocess.run([str(python), "-c", code, str(path)], capture_output=True, text=True, encoding="utf-8")
    require(result.returncode == 0, f"PDF inspection unavailable/failed: {result.stderr.strip()[-500:]}")
    report = json.loads(result.stdout)
    require(4 <= report["pages"] <= 6, f"Report has {report['pages']} pages; required 4-6")
    text = report["text"]
    require(len(text.split()) >= 600, "Report has insufficient extractable text for content review")
    lowered = text.lower()
    for terms in (("magnesium", "mg"), ("calcium", "ca"), ("acetamide", "amide"), ("sampling",), ("references", "bibliography"), ("codex", "ai assistance")):
        require(any(term in lowered for term in terms), f"Report content/provenance term missing: {terms}")
    excerpt = " ".join(text.split())[:240]
    return f"{report['pages']} pages, {len(text.split())} words; content excerpt: {excerpt}"


def verify_presentation(path):
    text = path.read_text(encoding="utf-8")
    words = len(text.split())
    require(words >= 250, "Presentation script is empty or too short to inspect")
    require(all(re.search(rf"\b{term}\b", text, re.I) for term in ("CM", "IOD")), "Presentation does not name both compared models")
    require(any(term in text.lower() for term in ("acetamide", "amide")), "Presentation omits the extension")
    return f"{words} words, approximately {words/140:.1f} minutes at 140 words/minute; pacing and scientific claims need review"


def self_check():
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "source"
        path.write_bytes(b"original")
        expected = digest(path)
        verify_hash(path, expected)
        path.write_bytes(b"tampered")
        try:
            verify_hash(path, expected)
        except ValueError:
            pass
        else:
            raise AssertionError("Tampered input was accepted")
    assert len(expected_cases("hydration")) == 12 and len(expected_cases("amide")) == 24
    nested = {"halves": {"CN": {"second_minus_first": .25}}, "cutoffs": [{"mean_CN": 6.0}]}
    assert same_value(json.dumps(nested), nested) and same_value("", None)
    assert same_value("False", False) and same_value("True", True)
    assert not same_value("True", False)
    changed = json.loads(json.dumps(nested))
    changed["halves"]["CN"]["second_minus_first"] = .5
    assert not same_value(changed, nested) and not same_value("NaN", 6.0)
    print("Self-check passed: hash/nested-metric tampering rejected and exact campaign case counts checked.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--analysis", type=Path, help="Defaults to RESULTS/analysis")
    parser.add_argument("--data-only", action="store_true", help="Audit raw data, summaries and figures before report generation")
    parser.add_argument("--report", type=Path, default=ROOT / "output/pdf/ion_model_benchmark_report.pdf")
    parser.add_argument("--presentation", type=Path, default=ROOT / "output/pdf/presentation_script.md")
    parser.add_argument("--pdf-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return 0
    failed = 0

    def check(label, function):
        nonlocal failed
        try:
            result = function()
        except (OSError, ValueError, KeyError, TypeError, StopIteration, subprocess.CalledProcessError) as error:
            failed += 1
            print(f"FAIL {label}: {error}")
            return None
        print(f"PASS {label}: {result if isinstance(result, str) else 'checked'}")
        return result

    check("parameter provenance", verify_sources)
    protocol = check("protocol readable", lambda: json.loads((ROOT / "protocol.json").read_text()))
    summaries = {}
    for campaign in ("hydration", "amide"):
        directory = args.results / campaign
        expected = expected_cases(campaign)
        paths = list(directory.rglob("metadata.json"))
        actual = {p.parent.name for p in paths}
        check(f"{campaign} inventory", lambda: require(actual == set(expected) and len(paths) == len(expected),
              f"missing={sorted(set(expected)-actual)}; unexpected={sorted(actual-set(expected))}") or f"{len(expected)} cases")
        if protocol is not None:
            for name in sorted(actual & set(expected)):
                summary = check(f"{campaign}/{name}", lambda: verify_case(directory / name / "metadata.json", campaign, expected[name], protocol))
                if summary is not None:
                    summaries[campaign, name] = summary
    analysis = args.analysis or args.results / "analysis"
    check("analysis matches raw data", lambda: verify_analysis(analysis, summaries))
    check("quantitative figures", lambda: verify_figures(analysis))
    if not args.data_only:
        check("English report", lambda: verify_report(args.report, args.pdf_python))
        check("presentation script", lambda: verify_presentation(args.presentation))
    status = subprocess.run(["git", "status", "--short"], cwd=ROOT, capture_output=True, text=True)
    print("INFO Git worktree: " + (status.stdout.strip() or ("clean" if status.returncode == 0 else status.stderr.strip())))
    print(f"Audit: {failed} incomplete/failed checks. This is an artifact-consistency audit, not a scientific certification.")
    return int(failed != 0)


if __name__ == "__main__":
    raise SystemExit(main())
