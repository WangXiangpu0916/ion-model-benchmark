# Divalent ion model benchmarks for NMDAR research

A small, reproducible molecular dynamics project for discussion with Ben Corry.
The project tests published Mg2+ and Ca2+ ion-water parameter combinations before
their use in more complex NMDA receptor simulations.

## Intended deliverables

- Runnable simulation and analysis code with documented parameter provenance.
- Two published ion-water model combinations, two ions and three independent
  replicas for the hydration benchmark.
- A small amide-coordination extension after evaluating the hydration baseline.
- Quantitative figures, machine-readable summaries and measured local throughput.
- A four-to-six-page English report and a five-minute presentation script.

## Scientific scope

The study compares hydration structure, local coordination and sampling
sensitivity. It does not establish NMDAR binding affinity, ion permeation,
calcium-magnesium knock-on, or physiological unblock rates. Model comparisons
include the prescribed water model; differences cannot automatically be assigned
to ion parameters alone.

## Version control

This directory is an independent Git repository on branch `main`. Code, protocols,
parameter sources, run manifests, summary data, figures and report sources are
versioned. Virtual environments, caches and large trajectories are excluded from
Git and retained locally when needed. Each completed, checked stage receives a
commit. No remote repository has been created or configured.

## Status

The source-verified protocol, simulation runner and quantitative analysis are
implemented. An OpenCL hydration pilot has completed; the full hydration
campaign is running. Acetamide parameters have passed source-versus-conversion
energy and force checks. Final biological conclusions and a final report are
pending completion and analysis of the production campaigns.

The local host has an NVIDIA GeForce RTX 3050 Laptop GPU with 4 GB VRAM and 16
logical CPU processors. A short pilot measured 392.6 ns/day; sustained production
subsequently measured about 80-90 ns/day. Use actual run metadata and logs for
timing, not the short-pilot estimate.

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
