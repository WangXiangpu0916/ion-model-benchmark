# Work log

## 5 October 2026

- Initialized a standalone Git repository before implementation as requested.
  Initial commit: `e1a44b8`.
- Created a project-local Python 3.12 environment; installed OpenMM 8.6.1,
  NumPy 2.5.3, SciPy 1.18.1 and Matplotlib 3.11.2. No global Python environment
  was changed. CUDA dependencies are installed locally but OpenMM exposed
  Reference, CPU and OpenCL; OpenCL uses the NVIDIA RTX 3050 Laptop GPU.
- OpenMM installation test passed; median Reference/OpenCL force difference
  6.74321e-06 in the installation test's relative-error metric.
- Verified Li-Merz 2013 CM and IOD parameters from pinned source files and
  matched native OpenMM CM values. Retained original FreeSolv acetamide inputs.
- Completed a 10 ps equilibration + 50 ps production pilot: 2640 atoms,
  879 waters, Mg2+ and two Cl-. Production throughput 392.6 ns/day; mean
  temperature 299.8 K; diagnostic Mg coordination number 6.0. This is a pilot,
  not evidence that slow exchange has equilibrated.
- Committed source parameters, numerical protocol, runner and pilot as
  `285d52e`. Began the twelve-case hydration campaign on this code revision.
- Added analysis of RDFs, coordination distributions, replica variation, block
  diagnostics, cutoff sensitivity and amide distances; analytical self-checks
  passed. Only completed, checksummed runs enter analysis.
- Converted the pinned FreeSolv acetamide topology to OpenMM XML. Original and
  four perturbed geometries passed source-versus-conversion energy/force checks.
- Enabled the optional CUDA plugin using a process-local DLL search path. A
  concurrent CUDA pilot completed; its throughput is not an isolated benchmark.
  The hydration production process continues on OpenCL.
- Added a read-only completion audit. Code review identified a final-PDB box
  header bug and missing intermediate sample saves in the simulation runner.
  The live runner remains unchanged during hydration to preserve recorded source
  identity. Repairs will be documented and applied before the amide campaign.
- Prepared a report generator that requires the complete 12+24-case dataset.
  Formal report generation remains pending; temporary layout fixtures do not
  represent simulation results.
- Completed the first 2 ns hydration case, CM/Mg replica 1: 127.4 ns/day,
  mean temperature 298.30 K, mean diagnostic CN 6.0 and first RDF peak 0.203 nm.
  This single case is not the completed model comparison. Its final PDB CRYST1
  edges were corrected from the cube root of the last saved NPT volume; atom
  coordinates and raw samples were unchanged. The metadata records the original
  and corrected headers and file hashes. The completed-case audit then passed.
- Reviewed that completed trajectory's RDF normalization: the shell-volume-
  weighted mean over 0.8-0.9 nm was 1.00759. Six water identities persisted
  throughout the saved frames. Its 0.203 nm peak is below the quoted Mg reference
  interval; this is an observed single-case result, not a model ranking.
- Applied the same documented CRYST1-only correction to both completed pilots,
  checking the serialized cubic box and isotropic barostat. Verified corrected
  PDB dimensions and hashes. The original OpenCL pilot predates the first code
  commit and has less provenance metadata; it remains excluded from the formal
  production comparison. Partial pilot plots are retained only under `tmp/`.

## Assistance and ownership

OpenAI Codex assisted with literature/source retrieval, experimental design,
code generation, execution, verification, analysis and document preparation.
Subagents performed scoped parameter, analysis and code-review tasks in the
same project. These checks are not independent experimental replication or
supervisory approval. Literature values, simulation output and interpretation
must remain distinguishable in the final report.

The user owns the scientific decisions and should understand and review the
code, source references and results before presenting or submitting the work.
