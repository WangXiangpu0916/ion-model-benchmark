"""Build the final six-page PDF and five-minute talk from a complete campaign.

Requires reportlab for PDF output; the separate --project-python environment must
contain the project's dependencies for verify_project.py's mandatory data audit.
Use a Python environment containing reportlab, e.g. the Codex bundled runtime.
No partial-report option: exactly 12 hydration and 24 amide cases are required.

python build_report.py --check
python build_report.py --output output/pdf
python build_report.py --self-check

The caller must render and visually review the resulting PDF before delivery.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import io
from itertools import product
import json
import math
from pathlib import Path
import re
import statistics
import subprocess
import sys
from xml.sax.saxutils import escape
from zoneinfo import ZoneInfo


ROOT = Path(__file__).resolve().parent
MODELS = ("CM", "IOD")
IONS = ("Mg", "Ca")
FIGURES = ("figure1_rdf.png", "figure2_hydration_metrics.png",
           "figure3_sampling_diagnostics.png", "figure4_sensitivity.png")
REFERENCES = [
    ("Li P, Roberts BP, Chakravorty DK, Merz KM Jr. Rational Design of Particle Mesh Ewald Compatible Lennard-Jones Parameters for +2 Metal Cations in Explicit Solvent. J Chem Theory Comput. 2013;9:2733-2748.", "10.1021/ct400146w"),
    ("Li Z, Song LF, Li P, Merz KM Jr. Systematic Parametrization of Divalent Metal Ions for the OPC3, OPC, TIP3P-FB, and TIP4P-FB Water Models. J Chem Theory Comput. 2020;16:4429-4442.", "10.1021/acs.jctc.0c00194"),
    ("Mobley DL, Guthrie JP. FreeSolv: a database of experimental and calculated hydration free energies, with input files. J Comput Aided Mol Des. 2014;28:711-720.", "10.1007/s10822-014-9747-x"),
    ("Duarte Ramos Matos G et al. Approaches for Calculating Solvation Free Energies and Enthalpies Demonstrated with an Update of the FreeSolv Database. J Chem Eng Data. 2017;62:1559-1569.", "10.1021/acs.jced.7b00104"),
    ("Grossfield A et al. Best Practices for Quantification of Uncertainty and Sampling Quality in Molecular Simulations. Living J Comput Mol Sci. 2018;1:5067.", "10.33011/livecoms.1.1.5067"),
]


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_keys():
    return {(campaign, model, ion, replica, state)
            for campaign in ("hydration", "amide")
            for model, ion, replica in product(MODELS, IONS, (1, 2, 3))
            for state in (("hydrated",) if campaign == "hydration" else ("contact", "separated"))}


def row_key(row):
    return tuple(row[k] for k in ("campaign", "model", "ion", "replica", "initial_state"))


def require_complete_keys(rows):
    keys = [row_key(row) for row in rows]
    expected = expected_keys()
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate case identities in final analysis")
    missing, extra = expected - set(keys), set(keys) - expected
    if missing or extra:
        raise ValueError(f"Final report requires all 12 hydration and 24 amide cases; missing={sorted(missing)}, unexpected={sorted(extra)}")


def finite(value, description):
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Invalid or unavailable {description}: {value!r}")
    return float(value)


def select(rows, campaign, ion, model, state=None):
    return [r for r in rows if r["campaign"] == campaign and r["ion"] == ion and r["model"] == model
            and (state is None or r["initial_state"] == state)]


def stats(rows, metric):
    values = [finite(r[metric], metric) for r in rows]
    if len(values) != 3:
        raise ValueError(f"Expected exactly three replicas for {metric}")
    return statistics.mean(values), statistics.stdev(values)


def formatted_stats(rows, metric, digits=3, scale=1):
    mean, sd = stats(rows, metric)
    return f"{mean * scale:.{digits}f} +/- {sd * scale:.{digits}f}"


def interval(values, digits=2):
    values = [finite(x, "range endpoint") for x in values]
    return f"{min(values):.{digits}f}-{max(values):.{digits}f}"


def audit_data(root, analysis, project_python):
    result = subprocess.run([str(project_python), str(root / "verify_project.py"), "--data-only",
                             "--results", str(root / "results"), "--analysis", str(analysis)],
                            cwd=root, text=True, encoding="utf-8", capture_output=True)
    if result.returncode != 0:
        raise ValueError("Required raw-data/summary audit failed:\n" + (result.stdout + result.stderr)[-4000:])
    return result.stdout.strip()


def load_inputs(root, analysis, project_python):
    summary_path = analysis / "replica_summaries.json"
    payload = read_json(summary_path)
    rows = payload["replicas"]
    require_complete_keys(rows)
    metadata = []
    for row in rows:
        case_name = row["case"]
        if not isinstance(case_name, str) or Path(case_name).name != case_name or case_name in (".", ".."):
            raise ValueError("Invalid case directory name")
        directory = root / "results" / row["campaign"] / case_name
        meta_path, samples_path = directory / "metadata.json", directory / "samples.npz"
        meta = read_json(meta_path)
        if meta.get("status") != "complete" or row_key(meta) != row_key(row):
            raise ValueError(f"Case is not complete or its identity changed: {directory}")
        if not samples_path.is_file() or sha256(samples_path) != meta.get("samples_sha256"):
            raise ValueError(f"Missing/changed saved simulation samples: {directory}")
        expected_frames = round(meta["production_ps"] / meta["sample_interval_ps"])
        if row["frames"] != expected_frames or meta.get("completed_frames") != expected_frames:
            raise ValueError(f"Incomplete production frames: {directory}")
        if row["water_count"] != meta["n_water"]:
            raise ValueError(f"Water counts differ between analysis and metadata: {directory}")
        for key in ("seed", "production_ps", "sample_interval_ps"):
            if row[key] != meta[key]:
                raise ValueError(f"Summary and metadata disagree on {key}: {directory}")
        for key in ("production_ps", "sample_interval_ps", "temperature_K"):
            if meta[key] != meta["protocol"][key]:
                raise ValueError(f"Top-level {key} contradicts embedded protocol: {directory}")
        if row["target_temperature_K"] != meta["temperature_K"]:
            raise ValueError(f"Summary temperature target differs from metadata: {directory}")
        if not math.isclose(row["last_saved_time_ps"], meta["production_ps"], rel_tol=1e-9, abs_tol=1e-6):
            raise ValueError(f"Incomplete saved production interval: {directory}")
        if finite(meta["min_box_length_nm"], "shortest box edge") <= 1.8:
            raise ValueError(f"RDF range extends beyond supported box geometry: {directory}")
        if summary_path.stat().st_mtime + 1 < meta_path.stat().st_mtime:
            raise ValueError(f"Analysis predates completed metadata; rerun analyze.py: {directory}")
        for metric in ("rdf_peak_nm", "mean_shell_distance_nm", "mean_CN", "mean_temperature_K",
                       "mean_water_number_density_nm3", "chloride_contact_fraction_lt_0p4nm"):
            finite(row[metric], metric)
        if row["campaign"] == "amide":
            fraction = finite(row["amide_fraction_within_diagnostic_water_cutoff"], "amide contact fraction")
            if not 0 <= fraction <= 1:
                raise ValueError("Invalid amide contact fraction")
            finite(row["mean_amide_distance_nm"], "amide distance")
        for key in ("ns_per_day", "production_wall_seconds", "wall_seconds"):
            if finite(meta[key], key) <= 0:
                raise ValueError(f"Nonpositive completed-case timing: {directory}")
        metadata.append(meta)
    for name in FIGURES:
        path = analysis / name
        if not path.is_file() or path.stat().st_mtime + 1 < summary_path.stat().st_mtime:
            raise ValueError(f"Missing/stale figure; regenerate analysis including --amide: {path}")
    source_manifest = read_json(root / "parameters" / "source_manifest.json")
    for source in source_manifest["files"]:
        path = root / source["file"]
        if sha256(path) != source["sha256"]:
            raise ValueError(f"Upstream parameter source checksum changed: {path}")
    verification = read_json(root / "parameters" / "freesolv" / "conversion_verification.json")
    if verification.get("passed") is not True:
        raise ValueError("Acetamide conversion verification did not pass")
    for name, digest in verification["sha256"].items():
        if sha256(root / "parameters" / "freesolv" / name) != digest:
            raise ValueError(f"Acetamide conversion input/output changed: {name}")
    # Require one calibration definition throughout; timing may differ by campaign.
    ion_parameters = metadata[0]["protocol"]["models"]
    if any(m["protocol"]["models"] != ion_parameters for m in metadata):
        raise ValueError("Ion parameter definitions differ between completed cases")
    sources = (root / "SOURCES.md").read_text(encoding="utf-8")
    if any(doi not in sources for _, doi in REFERENCES[:4]):
        raise ValueError("SOURCES.md lacks the parameter references used by this report")
    audit = audit_data(root, analysis, project_python)
    return {"rows": rows, "metadata": metadata, "verification": verification,
            "analysis": analysis, "root": root, "source_manifest": source_manifest, "data_audit": audit}


def facts(data):
    rows, metadata = data["rows"], data["metadata"]
    hydration = [r for r in rows if r["campaign"] == "hydration"]
    gaps, peak_changes = {}, {}
    for ion, model in product(IONS, MODELS):
        contact = stats(select(rows, "amide", ion, model, "contact"), "amide_fraction_within_diagnostic_water_cutoff")[0]
        separated = stats(select(rows, "amide", ion, model, "separated"), "amide_fraction_within_diagnostic_water_cutoff")[0]
        gaps[ion, model] = 100 * (contact - separated)
    for ion in IONS:
        peak_changes[ion] = stats(select(rows, "hydration", ion, "IOD"), "rdf_peak_nm")[0] - stats(select(rows, "hydration", ion, "CM"), "rdf_peak_nm")[0]
    throughput = sum(m["production_ps"] / 1000 for m in metadata) / sum(m["production_wall_seconds"] for m in metadata) * 86400
    cn_drift = [abs(r["halves"]["CN"]["second_minus_first"]) for r in hydration]
    temp_drift = [abs(r["halves"]["temperature_K"]["second_minus_first"]) for r in hydration]
    sensitivity = [max(x["mean_CN"] for x in r["cutoff_sensitivity"]) - min(x["mean_CN"] for x in r["cutoff_sensitivity"]) for r in hydration]
    errors = [v for geometry in data["verification"]["geometries"] for v in geometry["by_force"].values()]
    return {
        "gaps": gaps, "peak_changes": peak_changes, "throughput": throughput,
        "total_ns": sum(m["production_ps"] for m in metadata) / 1000,
        "case_wall_hours": sum(m["wall_seconds"] for m in metadata) / 3600,
        "max_cn_half_change": max(cn_drift), "max_temperature_half_change": max(temp_drift),
        "max_cutoff_CN_span": max(sensitivity),
        "no_identity_changes": sum(r["shell_identity_changed_intervals"] == 0 for r in hydration),
        "max_conversion_energy_error": max(e["absolute_energy_error_kj_mol"] for e in errors),
        "max_conversion_force_error": max(e["max_absolute_force_error_kj_mol_nm"] for e in errors),
        # One ion per box: <1/V> follows exactly from fixed water count / frame density.
        "cation_M": [1e24 / 6.02214076e23 * r["mean_water_number_density_nm3"] / r["water_count"] for r in rows],
    }


def protocol_values(data, key, campaign=None, scale=1, digits=2):
    metas = [m for m in data["metadata"] if campaign is None or m["campaign"] == campaign]
    values = sorted({float(m["protocol"][key]) * scale for m in metas})
    return f"{values[0]:.{digits}f}" if len(values) == 1 else interval(values, digits)


def build_pdf(data, author, report_date):
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.pdfgen.canvas import Canvas
        from reportlab.platypus import Paragraph, Table, TableStyle
        from reportlab.lib.utils import ImageReader
    except ImportError as exc:
        raise RuntimeError("PDF generation requires reportlab. Use the documented bundled Python or an environment with reportlab installed.") from exc

    rows, metadata = data["rows"], data["metadata"]
    f = facts(data)
    buffer = io.BytesIO()
    canvas = Canvas(buffer, pagesize=A4, pageCompression=1, invariant=1)
    canvas.setTitle("Ion-model calibration and coordination: an NMDAR-motivated benchmark")
    canvas.setAuthor(author)
    width, height = A4
    margin, content_width = 42, width - 84
    navy, teal = colors.HexColor("#17364B"), colors.HexColor("#147D92")
    styles = {
        "body": ParagraphStyle("body", fontName="Helvetica", fontSize=9, leading=12, textColor=navy),
        "small": ParagraphStyle("small", fontName="Helvetica", fontSize=7.7, leading=10, textColor=navy),
        "title": ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=20, leading=24, textColor=navy),
        "heading": ParagraphStyle("heading", fontName="Helvetica-Bold", fontSize=12, leading=15, textColor=teal),
        "cell": ParagraphStyle("cell", fontName="Helvetica", fontSize=7.8, leading=10, textColor=navy),
    }
    page, y = 0, 0

    def check_space(h):
        if y - h < 46:
            raise ValueError(f"Page {page} content overflows: reduce content/figure height before final generation")

    def para(text, style="body", gap=7):
        nonlocal y
        item = Paragraph(text, styles[style])
        _, h = item.wrap(content_width, height)
        check_space(h + gap)
        item.drawOn(canvas, margin, y - h)
        y -= h + gap

    def new_page(title):
        nonlocal page, y
        if page:
            canvas.showPage()
        page += 1
        canvas.setFillColor(navy)
        canvas.setFont("Helvetica", 7.5)
        canvas.drawString(margin, height - 24, "LOCAL MOLECULAR SIMULATION STUDY")
        canvas.drawRightString(width - margin, height - 24, report_date)
        canvas.setStrokeColor(teal)
        canvas.setLineWidth(1.2)
        canvas.line(margin, height - 31, width - margin, height - 31)
        canvas.setFillColor(colors.HexColor("#65717C"))
        canvas.setFont("Helvetica", 7)
        canvas.drawString(margin, 26, "CM/IOD + TIP3P | Hydration and acetamide | Descriptive three-replica comparisons")
        canvas.drawRightString(width - margin, 26, f"{page} / 6")
        y = height - 46
        para(title, "title", 10)

    def table(header, body, fractions):
        nonlocal y
        cells = [[Paragraph(escape(str(cell)), styles["cell"]) for cell in line] for line in [header] + body]
        item = Table(cells, colWidths=[content_width * q for q in fractions], hAlign="LEFT")
        item.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E5F1F4")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F7F9")]),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 6), ("RIGHTPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, teal),
        ]))
        _, h = item.wrap(content_width, height)
        check_space(h + 9)
        item.drawOn(canvas, margin, y - h)
        y -= h + 9

    def figure(name, max_height=None):
        nonlocal y
        image = ImageReader(str(data["analysis"] / name))
        iw, ih = image.getSize()
        scale = min(content_width / iw, max_height / ih if max_height else content_width / iw)
        draw_width, draw_height = iw * scale, ih * scale
        check_space(draw_height + 5)
        canvas.drawImage(image, margin + (content_width - draw_width) / 2, y - draw_height,
                         width=draw_width, height=draw_height, mask="auto")
        y -= draw_height + 5

    new_page("Ion-model calibration and coordination")
    para("An NMDAR-motivated Ca2+/Mg2+ benchmark in TIP3P water and acetamide", "heading", 5)
    para(f"{escape(author)} | Completed local computational study | {f['total_ns']:.1f} ns of production across 36 cases", "small", 6)
    para("<b>Question.</b> How much do two published ion calibrations change local hydration and finite-time amide contact, before attempting an NMDAR pore calculation?", gap=5)
    para(f"<b>Observed.</b> IOD-minus-CM mean RDF-peak shifts were {f['peak_changes']['Mg']:+.4f} nm for Mg and {f['peak_changes']['Ca']:+.4f} nm for Ca. The table records the independent replica variation and initial-state contrasts; it makes no significance or affinity claim.", gap=6)
    body = []
    for ion, model in product(IONS, MODELS):
        group = select(rows, "hydration", ion, model)
        body.append([f"{ion} / {model}", formatted_stats(group, "rdf_peak_nm", 4),
                     formatted_stats(group, "mean_CN", 2), f"{f['gaps'][ion, model]:+.1f}"])
    table(["Group", "RDF peak, nm (mean +/- SD)", "Water CN (mean +/- SD)", "Amide contact-start minus separated-start fraction (percentage points)"],
          body, [0.13, 0.27, 0.23, 0.37])
    figure("figure2_hydration_metrics.png", max_height=397)
    para("Figure 1. Hydration metrics and fitting references [1,2]. Three points per group are separate simulations; SD is not a confidence interval. Amide contact uses the diagnostic water-shell cutoff and is defined on page 5.", "small", 0)

    new_page("Methods and parameter provenance")
    para("The same full-charge (+2 e) Ca2+ and Mg2+ models were compared in standard TIP3P water. CM is Li-Merz's compromise 12-6 set; IOD is the distance-targeted set [1]. Lorentz-Berthelot mixing and PME were retained. No charge scaling or 12-6-4 term was introduced.")
    protocol_rows = [
        ["Design", "2 ions x 2 models x 3 replicas; plus contact/separated amide starts (24 cases)"],
        ["Temperature / pressure", f"{protocol_values(data, 'temperature_K')} K / {protocol_values(data, 'pressure_bar')} bar"],
        ["Initial periodic cube / counterions", f"{protocol_values(data, 'box_nm')} nm / one cation + two chloride ions; no extra salt"],
        ["Finite-concentration system", f"Cation {interval(f['cation_M'], 4)} M across case means from actual volumes; chloride is twice this concentration"],
        ["Integrator", f"Langevin-middle; step {protocol_values(data, 'timestep_ps', scale=1000, digits=1)} fs; friction {protocol_values(data, 'friction_per_ps')} per ps"],
        ["Constraints / pressure control", "Rigid water and H-bond constraints; isotropic Monte Carlo barostat every 25 steps"],
        ["PME / nonbonded cutoff", f"Error tolerance {protocol_values(data, 'ewald_tolerance', digits=5)}; real-space cutoff {protocol_values(data, 'cutoff_nm')} nm"],
        ["Equilibration / saved cadence", f"{protocol_values(data, 'equilibration_ps')} ps / {protocol_values(data, 'sample_interval_ps')} ps"],
        ["Production per case", f"Hydration {protocol_values(data, 'production_ps', 'hydration', .001)} ns; amide {protocol_values(data, 'production_ps', 'amide', .001)} ns"],
    ]
    table(["Item", "Recorded setup"], protocol_rows, [.31, .69])
    parameters = metadata[0]["protocol"]["models"]
    parameter_rows = []
    for ion, model in product(IONS, MODELS):
        p = parameters[model][ion]
        parameter_rows.append([f"{ion} / {model}", f"{p['rmin_half_angstrom']:.3f}", f"{p['epsilon_kcal_mol']:.8f}"])
    table(["Group", "Amber Rmin/2 (angstrom)", "epsilon (kcal/mol)"], parameter_rows, [.24, .38, .38])
    para(f"Replicas share a solvent-packing template but use separately seeded velocities, barostat moves and chloride placement. They test stochastic repeatability of this preparation; they do not independently sample every possible initial hydration state. The approximately {statistics.mean(f['cation_M']):.3f} M mean cation concentration is finite, not infinite dilution.", "small")
    para("<b>Acetamide.</b> FreeSolv mobley_8048190 supplies GAFF/AM1-BCC parameters [3,4]. Its neutral, nine-atom topology was converted using OpenMM's Amber reader, preserving exclusions and 1-4 scaling. Carbonyl oxygen O1 was initially placed in direct contact or at 0.65 nm; each start was then minimized and equilibrated. These are preparation labels, not enforced production states.")
    para(f"<b>Conversion check.</b> The source geometry and four perturbed geometries passed independent source-topology versus XML comparisons. Maximum per-force energy and atomic-force differences were {f['max_conversion_energy_error']:.3g} kJ/mol and {f['max_conversion_force_error']:.3g} kJ/mol/nm. This verifies format conversion, not metal-amide accuracy.")
    para("Exact source commits, original files and SHA-256 checksums are retained in SOURCES.md and parameters/source_manifest.json. All completed cases carry the protocol, code hash, initial/final structures and runtime metadata.", "small", 0)

    new_page("Hydration structure: the measured comparison")
    figure("figure1_rdf.png")
    para("Figure 2. Individual-replica metal-water oxygen RDFs and equally weighted group means. Each frame's histogram is normalized by its own water number density. The CSV retains the full 0-0.9 nm range; histogram bins are 0.002 nm.", "small")
    table(["Group", "Pooled in-shell mean distance (nm)", "Mean water CN", "Counterion contact fraction (%)"],
          [[f"{ion} / {model}", formatted_stats(select(rows, "hydration", ion, model), "mean_shell_distance_nm", 4),
            formatted_stats(select(rows, "hydration", ion, model), "mean_CN", 2),
            formatted_stats(select(rows, "hydration", ion, model), "chloride_contact_fraction_lt_0p4nm", 1, 100)]
           for ion, model in product(IONS, MODELS)], [.16, .30, .24, .30])
    para("The pooled in-shell mean distance and the RDF peak are different estimators. The mean conditions on all water oxygens inside the declared shell; the peak is the maximum RDF bin within the predeclared ion-specific window. A peak value is therefore discrete at this bin width, and a small replica SD does not eliminate binning or model uncertainty.")
    para("The reference peak distances are 0.209 +/- 0.004 nm for Mg and 0.246 nm for Ca; reference coordination values are 6 and 8 [2]. These values belong to the calibration lineage. Agreement, especially for IOD, shows target recovery rather than an independent prediction. No uncertainty is invented for the Ca target or the coordination values.")
    para("Water-shell cutoffs are 0.30 nm for Mg and 0.35 nm for Ca. Counterion contact is the fraction of saved frames with the nearest chloride below 0.40 nm. A counterion contact can alter the local environment and is reported rather than silently discarded. Acetamide was absent from these hydration comparisons.")

    new_page("Sampling diagnostics and definition sensitivity")
    figure("figure3_sampling_diagnostics.png")
    para("Figure 3. Consecutive block means in the designated production segment. Blocks display time dependence; they are not additional independent replicas. No automatic trimming was used to improve agreement.", "small")
    diagnostic_rows = []
    for ion, model in product(IONS, MODELS):
        group = select(rows, "hydration", ion, model)
        diagnostic_rows.append([f"{ion} / {model}",
            interval([r["halves"]["CN"]["second_minus_first"] for r in group], 3),
            interval([r["halves"]["temperature_K"]["second_minus_first"] for r in group], 2),
            interval([max(v["mean_CN"] for v in r["cutoff_sensitivity"]) - min(v["mean_CN"] for v in r["cutoff_sensitivity"]) for r in group], 3)])
    table(["Group", "Second-half minus first-half CN (replica range)", "Temperature half-change (K; range)", "CN span across cutoff +/-0.02 nm (range)"], diagnostic_rows, [.14, .30, .28, .28])
    para(f"Across hydration replicas, the largest absolute half-to-half CN change was {f['max_cn_half_change']:.3f}; the largest cutoff-dependent CN span was {f['max_cutoff_CN_span']:.3f}. These are diagnostic magnitudes, not a convergence test with a universal acceptance threshold.")
    para(f"In {f['no_identity_changes']} of 12 hydration trajectories, no first-shell membership change was observed between saved frames. A fixed coordination number can conceal water replacement, and unchanged water identities may reflect an unsampled slow process. Saved-frame observations do not establish committed exchange events or residence times.")
    para("Three independently seeded preparations provide descriptive repeatability checks. All replica values and the sample SD are reported, without independent-frame confidence intervals or p-value-based rankings. Shared preparation procedures and unvisited states can still limit independence and coverage [5].")

    new_page("Acetamide: finite-time initial-state contrasts")
    figure("figure4_sensitivity.png", max_height=340)
    para("Figure 4. Left: hydration shell-cutoff sensitivity. Right: metal-carbonyl oxygen distance in the 24 amide simulations. Solid/dashed preparation labels denote contact/separated starts; replicas are separate trajectories.", "small", 5)
    amide_rows = []
    for ion, model, state in product(IONS, MODELS, ("contact", "separated")):
        group = select(rows, "amide", ion, model, state)
        amide_rows.append([f"{ion} / {model}", state,
            formatted_stats(group, "mean_amide_distance_nm", 3),
            formatted_stats(group, "amide_fraction_within_diagnostic_water_cutoff", 1, 100)])
    table(["Group", "Prepared start", "Mean carbonyl-O distance (nm)", "Diagnostic contact fraction (%)"], amide_rows, [.17, .19, .32, .32])
    para("Diagnostic contact is carbonyl oxygen below 0.30 nm (Mg) or 0.35 nm (Ca), reusing the water-shell cutoff as an explicit operational definition. Each table entry is the mean +/- SD across three replicas. Neither the cutoff nor a contact-start trajectory defines a calibrated metal-amide bound state.", "small")
    para("Contact-start minus separated-start fraction contrasts (percentage points): " + "; ".join(f"{ion}/{model} {f['gaps'][ion, model]:+.1f}" for ion, model in product(IONS, MODELS)) + ". Initial-state contrasts reveal finite-time sensitivity; similar fractions alone would not establish equilibration. No equilibrium affinity, free-energy difference or kinetic rate is inferred.", "small", 0)

    new_page("Interpretation, reproducibility and references")
    para("<b>Contribution.</b> This is a completed, reproducible small-system comparison of calibration choices and local coordination observables. It supplies measured hydration differences and an amide sampling diagnostic relevant to planning an Asn-associated ion model. It does not simulate an NMDAR receptor or establish its pore mechanism.")
    para("<b>Next discriminating calculation.</b> Establish reversible sampling and an external metal-amide benchmark before using these combinations for quantitative binding claims. A later pore study would additionally need construct-specific geometry, an adequate ion/protein interaction model and a converged multidimensional sampling strategy. The current data do not select a universal winning force field.")
    devices = sorted({m.get("platform_properties", {}).get("DeviceName", m["platform"]) for m in metadata})
    versions = sorted({m["openmm_version"] for m in metadata})
    para(f"<b>Measured feasibility.</b> Production throughput across all 36 cases was {f['throughput']:.1f} ns/day (total simulated time divided by summed production wall time); case-specific values spanned {interval([m['ns_per_day'] for m in metadata], 1)} ns/day. Summed case wall time was {f['case_wall_hours']:.2f} hours. Pilot throughput is excluded. Device: {escape(', '.join(devices))}. OpenMM: {escape(', '.join(versions))}.", "small")
    para("<b>Reproduction.</b> Run simulate.py for hydration and amide with the committed protocol, then analyze.py with both campaigns, and build_report.py. The builder needs ReportLab; its separate project Python runs the mandatory verify_project.py data audit. Simulation dependencies are in requirements.txt. Raw NPZ/DCD/checkpoint files remain local and Git-ignored; tables, figures, manifests and scripts retain the audit record. Final generation requires 36 complete cases, raw-sample hashes and summaries recomputed from raw data.", "small")
    para("<b>Generative-AI contribution.</b> Codex assisted project scoping, source checking, code authoring, analysis and report drafting. Numerical findings are computed from the saved molecular-dynamics outputs; generated text does not replace parameter provenance or the independent checks. This assistance should remain disclosed when the work is shown to a supervisor.", "small")
    para("References", "heading", 4)
    for index, (citation, doi) in enumerate(REFERENCES, 1):
        para(f"[{index}] {escape(citation)} <link href='https://doi.org/{doi}' color='#147D92'>{doi}</link>", "small", 5)
    if page != 6:
        raise AssertionError("Report must have exactly six pages")
    canvas.save()
    return buffer.getvalue()


def presentation(data, author, report_date):
    rows = data["rows"]
    f = facts(data)
    mg_cm = stats(select(rows, "hydration", "Mg", "CM"), "rdf_peak_nm")[0]
    mg_iod = stats(select(rows, "hydration", "Mg", "IOD"), "rdf_peak_nm")[0]
    ca_cm = stats(select(rows, "hydration", "Ca", "CM"), "rdf_peak_nm")[0]
    ca_iod = stats(select(rows, "hydration", "Ca", "IOD"), "rdf_peak_nm")[0]
    largest = max(f["gaps"], key=lambda key: abs(f["gaps"][key]))
    text = f"""# Five-minute presentation script

{author} | {report_date}

Suggested visuals: report pages 1, 2, 3, 4 and 5. The figures show actual saved simulations. Speak the numerical contrasts as descriptive findings; do not add a significance or affinity claim.

## 0:00-0:50 - The question

My white paper asked how calcium could modify magnesium block in NMDA receptors. I chose a smaller question that I could finish on my own computer: how much do published ion calibration choices change hydration and contact with a neutral amide?

Acetamide provides an isolated amide group related to asparagine chemistry. It does not reproduce an asparagine cage. The purpose was to examine model behaviour and sampling requirements before making a pore-mechanism claim.

## 0:50-1:45 - What I ran and checked

I compared the Li-Merz compromise, or CM, and ion-oxygen-distance, or IOD, twelve-six parameter sets. Both use full plus-two charges and the same TIP3P water and interaction conventions. Each box contained one divalent ion and two chlorides, at approximately {statistics.mean(f['cation_M']):.3f} molar cation concentration from actual volumes, rather than infinite dilution.

There were twelve hydration simulations: two ions, two models and three replicas. The solvent-packing template is shared, with separately seeded velocities, barostat moves and chloride placement. I then ran twenty-four acetamide simulations, adding direct-contact and separated starting preparations. These produced {f['total_ns']:.1f} nanoseconds of recorded production. The acetamide conversion was checked against the original Amber topology using energies and forces in five geometries; a successful file conversion alone would not establish equivalent physics.

## 1:45-2:40 - Hydration results

For magnesium, the mean first RDF peak was {mg_cm:.4f} nanometres with CM and {mg_iod:.4f} with IOD. For calcium, the values were {ca_cm:.4f} and {ca_iod:.4f}. I show every replica, the group mean and standard deviation. Full distributions are available as machine-readable results.

These reference distances helped define the model family. Recovering them, particularly with IOD, checks implementation and target recovery; it does not independently demonstrate better binding predictions. The RDF peak and the mean distance inside a shell are also reported separately because they are different measurements.

## 2:40-3:55 - Sampling and the amide extension

I checked production halves, consecutive blocks, counterion contact and sensitivity to the shell boundary. The largest absolute hydration coordination change between halves was {f['max_cn_half_change']:.3f}. Moving the cutoff by plus or minus 0.02 nanometres gave a largest coordination span of {f['max_cutoff_CN_span']:.3f}. These diagnostics do not provide a universal convergence certificate.

For the amide test, I measured the carbonyl oxygen distance and the fraction of saved frames inside a stated diagnostic cutoff. The largest absolute contact-start versus separated-start fraction contrast occurred for {largest[0]} with {largest[1]}: contact-start minus separated-start was {f['gaps'][largest]:+.1f} percentage points. The other contrasts are in the table, and all twenty-four trajectories enter the analysis.

Persistent initial-state dependence is a sampling warning. Similar contact fractions would also be insufficient to prove equilibration. These finite-time fractions do not establish an affinity, free-energy difference or residence time.

## 3:55-5:00 - What the project delivers

The final deliverable is a local repository containing parameter provenance, runnable simulation and analysis code, the recorded protocol, quantitative figures and this short report. Actual production throughput averaged {f['throughput']:.1f} nanoseconds per day across the completed cases. That is a measured property of this workflow; the earlier short pilot is excluded from the estimate.

The scientific outcome is a bounded comparison of hydration structure and amide sampling sensitivity. It does not establish a preferred NMDAR ion model, simultaneous calcium-magnesium occupancy, or a knock-on mechanism. My next step would be to establish reversible sampling and an external metal-amide benchmark before attempting a quantitative binding calculation or a pore model.

Codex assisted the source work, code and drafting, and that contribution is disclosed in the report. The numerical claims can be traced back to the saved trajectories and inspected independently. I would like to discuss whether this benchmark is a useful starting point for the group's ion-selectivity work, and which next observable would most strongly discriminate between models.
"""
    words = len(re.findall(r"\b[\w'-]+\b", text))
    return text + f"\nScript length: approximately {words} words including headings; rehearse near {words / 5:.0f} words/minute and adjust pauses to five minutes.\n"


def self_check():
    keys = expected_keys()
    assert len(keys) == 36
    assert sum(k[0] == "hydration" for k in keys) == 12
    rows = [dict(zip(("campaign", "model", "ion", "replica", "initial_state"), key)) for key in keys]
    require_complete_keys(rows)
    for invalid in (rows[:-1], rows + rows[:1]):
        try:
            require_complete_keys(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("Incomplete or duplicate final campaign was accepted")
    assert stats([{"x": 5.}, {"x": 6.}, {"x": 7.}], "x") == (6., 1.)
    assert formatted_stats([{"x": 5.}, {"x": 6.}, {"x": 7.}], "x", 1) == "6.0 +/- 1.0"
    print("Report self-check passed: complete-case gate, duplicate rejection and replica statistics; no report generated.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--analysis", type=Path, default=ROOT / "results" / "analysis")
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "pdf")
    local_python = ROOT / ".venv" / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    parser.add_argument("--project-python", type=Path, default=local_python if local_python.is_file() else Path(sys.executable),
                        help="Python with the project dependencies; used for the mandatory raw/summary audit")
    parser.add_argument("--author", default="Xiangpu Wang")
    parser.add_argument("--date", help="Report date, YYYY-MM-DD; defaults to the current Australia/Sydney date")
    parser.add_argument("--check", action="store_true", help="Validate complete inputs without creating artifacts")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    try:
        data = load_inputs(ROOT, args.analysis.resolve(), args.project_python)
        if args.check:
            print("All 12 hydration and 24 amide cases, sample hashes, converted parameters and final figures verified. No artifacts generated.")
            return
        report_date = args.date or datetime.now(timezone.utc).astimezone(ZoneInfo("Australia/Sydney")).date().isoformat()
        datetime.strptime(report_date, "%Y-%m-%d")
        pdf = build_pdf(data, args.author, report_date)
        script = presentation(data, args.author, report_date)
        args.output.mkdir(parents=True, exist_ok=True)
        pdf_path = args.output / "ion_model_benchmark_report.pdf"
        script_path = args.output / "presentation_script.md"
        pdf_path.write_bytes(pdf)
        script_path.write_text(script, encoding="utf-8")
        try:
            commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        except (OSError, subprocess.SubprocessError):
            commit = "unavailable"
        manifest = {
            "report_date_Australia_Sydney": report_date, "author": args.author,
            "git_commit_at_build": commit, "case_counts": {"hydration": 12, "amide": 24},
            "summary_sha256": sha256(args.analysis / "replica_summaries.json"),
            "figure_sha256": {name: sha256(args.analysis / name) for name in FIGURES},
            "SOURCES_md_sha256": sha256(ROOT / "SOURCES.md"),
            "pdf_sha256": hashlib.sha256(pdf).hexdigest(),
            "presentation_sha256": sha256(script_path),
            "data_audit": data["data_audit"],
            "review_status": "Generated; rendering and visual QA required before delivery",
        }
        (args.output / "report_build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Final report refused: {exc}\n")
    print(f"Created six-page report: {pdf_path.resolve()}")
    print(f"Created five-minute presentation script: {script_path.resolve()}")
    print("Render all six pages and visually review before delivery.")


if __name__ == "__main__":
    main()
