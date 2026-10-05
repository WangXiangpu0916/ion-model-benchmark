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
- Completed and audited CM/Mg replica 2 (118.3 ns/day). Its mean temperature
  was 298.36 K, first RDF peak 0.203 nm and mean CN 6.0, matching the first
  replica's local geometry at the recorded resolution. Six water identities
  persisted through its saved frames. This repeatability does not establish
  equilibrium exchange or the still-pending CM/IOD comparison. Its final PDB
  header received the same provenance-recorded correction.

## 6 October 2026

- Completed and audited CM/Mg replica 3 (141.3 ns/day), including the documented
  final-PDB correction. All three CM/Mg replicas now have their planned 2000
  production frames. Their RDF peaks are all 0.203 nm and water CN is 6.0.
  The pooled in-shell distance, averaged across replicas, is
  0.204269 +/- 0.000031 nm (sample SD); mean temperature is
  298.40 +/- 0.13 K. Outer RDF averages over 0.8-0.9 nm span 1.00545-1.00759.
  No saved-frame shell membership changes were observed in these three runs;
  local repeatability does not resolve slow exchange or establish equilibrium.
- Runtime inspection found concurrent interactive GPU use. No other process or
  hardware setting was changed. Report timing is consequently described as
  observed workflow throughput rather than isolated hardware performance.
- Extended the existing audit to verify recorded corrected-PDB hashes against
  their referenced sample hashes. Confirmed Git preserves those PDB bytes.
- Completed and audited CM/Ca replica 1 (140.6 ns/day): RDF peak 0.247 nm,
  mean shell distance 0.252668 nm, mean CN 8.0515 and temperature 298.20 K.
  Diagnostic CN probabilities were 0.0055, 0.9375 and 0.0570 for CN 7, 8 and 9.
  The outer RDF average over 0.8-0.9 nm was 1.00010. There were 82 distinct
  in-shell water identities and 178 membership flips between saved frames;
  these are descriptive observations, not committed exchanges or a rate.
  No chloride distance below the diagnostic 0.4 nm cutoff was observed.
  Two further CM/Ca replicas remain pending.
- Completed and audited CM/Ca replica 2 (112.8 ns/day): RDF peak 0.249 nm,
  mean shell distance 0.252466 nm, CN 8.0335 and temperature 298.57 K.
  CN 8 occupied 96.15% of saved frames; 65 distinct shell waters and 130
  membership flips were observed. The half-to-half CN difference was -0.027;
  no diagnostic chloride contacts occurred. The two Ca peaks differ by one
  0.002 nm histogram bin; the third replica remains pending.
- Increased diagnostic-figure label sizes and replaced repeated amide legends
  with four model/start entries while retaining every replica curve. Numerical
  analysis is unchanged. Self-checks and a clearly marked synthetic six-page
  layout fixture passed; report pages 4 and 5 were visually inspected. Final
  real-data figures and report still require their own visual review.
- Completed and audited CM/Ca replica 3 (126.3 ns/day); the CM hydration
  baseline now contains all six planned trajectories. Across the three Ca
  replicas, mean CN is 8.042 +/- 0.00904, RDF peak 0.24833 +/- 0.00115 nm,
  and mean in-shell distance 0.252530 +/- 0.000119 nm (between-replica sample
  SDs). Each trajectory was dominated by CN 8, and all three showed saved-frame
  water-identity changes. The six IOD hydration trajectories are next; these
  CM-only results do not yet establish the parameter-set comparison.

- Completed and audited IOD/Mg replica 1 (125.7 ns/day), including the recorded
  final-PDB header correction. Its RDF peak is 0.207 nm, mean shell distance
  0.208576 nm, mean CN 6.0 and temperature 298.46 K. The peak is 0.004 nm
  above the completed CM/Mg group and inside the quoted Mg reference interval.
  Six water identities persisted throughout the saved frames. This is one IOD
  trajectory; the remaining replicas are needed for the planned comparison.

## Assistance and ownership

OpenAI Codex assisted with literature/source retrieval, experimental design,
code generation, execution, verification, analysis and document preparation.
Subagents performed scoped parameter, analysis and code-review tasks in the
same project. These checks are not independent experimental replication or
supervisory approval. Literature values, simulation output and interpretation
must remain distinguishable in the final report.

The user owns the scientific decisions and should understand and review the
code, source references and results before presenting or submitting the work.
