# Divalent ion model benchmarks for NMDAR research

A small, reproducible molecular dynamics project for discussion with Ben Corry.
The project tests published Mg2+ and Ca2+ ion-water parameter combinations before
their use in more complex NMDA receptor simulations.

## Intended deliverables

- Runnable simulation and analysis code with documented parameter provenance.
- Two published ion-water model combinations, two ions and three seeded
  replicas for the hydration benchmark.
- A small amide-coordination extension after evaluating the hydration baseline.
- Quantitative figures, machine-readable summaries and measured local throughput.
- A four-to-six-page English report and a five-minute presentation script.

## Scientific scope

The study compares hydration structure, local coordination and sampling
sensitivity. It does not establish NMDAR binding affinity, ion permeation,
calcium-magnesium knock-on, or physiological unblock rates. CM and IOD use the
same TIP3P water and interaction conventions, so the comparison isolates these
two published ion-parameter choices within that fixed setup.

## Version control

This directory is an independent Git repository on branch `main`. Code, protocols,
parameter sources, run manifests, summary data, figures and report sources are
versioned. Virtual environments, caches and large trajectories are excluded from
Git and retained locally when needed. Each completed, checked stage receives a
commit. No remote repository has been created or configured.

## Status

All twelve hydration trajectories are complete and have passed the recorded
quality gate. Their analysis and four figures are saved in `results/analysis`;
the final combined analysis will replace this hydration-stage snapshot after
the twenty-four planned amide trajectories finish. Acetamide parameters have
passed source-versus-conversion energy and force checks. The final report and
presentation remain pending completion of that extension.

The local host has an NVIDIA GeForce RTX 3050 Laptop GPU with 4 GB VRAM and 16
logical CPU processors. A short pilot measured 392.6 ns/day; completed hydration
production rates span 112.8-508.4 ns/day under variable desktop/GPU load. Use
actual run metadata and logs for timing, not the short-pilot estimate.

## Run locally

From this directory in PowerShell, the existing environment can be used directly:

```powershell
.venv/Scripts/python.exe simulate.py --selfcheck
.venv/Scripts/python.exe analyze.py --self-check
.venv/Scripts/python.exe simulate.py --campaign hydration
.venv/Scripts/python.exe analyze.py --input results/hydration --output results/analysis
```

Only after evaluating the hydration diagnostics, run the amide extension and
combined analysis:

```powershell
.venv/Scripts/python.exe simulate.py --campaign amide
.venv/Scripts/python.exe analyze.py --input results/hydration --amide results/amide --output results/analysis
```

Do not start a duplicate campaign while the current one is live. The runner
checks existing completed cases using the recorded protocol and sample hash.
An incomplete directory requires inspection; it is never silently overwritten.

For a fresh environment, install Python 3.12 and use `uv venv --python 3.12 .venv`
followed by `uv pip install --python .venv/Scripts/python.exe -r requirements.txt`.
Core requirements support OpenCL/CPU. `requirements-local-lock.txt` records the
exact local environment including optional CUDA components. The Windows CUDA
launcher `run_gpu.py` exposes locally installed NVIDIA DLLs only to its process;
it accepts the same arguments as `simulate.py`, including `--platform CUDA`.

After both campaigns and their combined analysis finish, generate the report:

```powershell
python build_report.py --project-python .venv/Scripts/python.exe --output output/pdf
.venv/Scripts/python.exe verify_project.py --pdf-python (Get-Command python).Source
```

Here `python` denotes a document environment containing ReportLab 4.4.9 and
`pypdf`; the local Codex bundled Python provides both. The report builder invokes
the project Python to recompute and audit saved-data metrics before writing a
six-page PDF, an English presentation script and their checksum manifest.
The final PDF also requires visual review of every rendered page.

## What is measured

Read [PROTOCOL.md](PROTOCOL.md) before interpreting results and [SOURCES.md](SOURCES.md)
for original parameter files, fitting targets and the acetamide source. Water
frames are correlated; the three seeded trajectories are the replication units.
Agreement with an IOD fitting target is not independent predictive validation.
Amide starting-state differences measure finite-time sampling dependence, not
an equilibrium binding affinity.

All raw samples and trajectories stay in their local case directories. To share
a complete research archive, include those ignored data files separately from
the Git repository. The versioned manifests contain SHA-256 hashes for checking
the corresponding sample files. Git by itself is not a backup of ignored data.
