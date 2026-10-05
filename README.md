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

Repository initialized. Parameter selection and local runtime preparation are in
progress. No molecular dynamics results have been generated yet.

The local host has an NVIDIA GeForce RTX 3050 Laptop GPU with 4 GB VRAM and 16
logical CPU processors. Actual simulation throughput remains to be measured.
