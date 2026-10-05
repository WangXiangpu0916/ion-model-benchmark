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

- Focused review of all seven then-complete trajectories found no new issue
  requiring a rerun. IOD/Mg replica 1 has half-trajectory RDF maxima of 0.205
  and 0.209 nm versus its pooled 0.207 nm maximum, with similar heights in
  neighboring bins. Peak-bin movement must be interpreted alongside shell
  distances and replicate variation; the outer RDF average was 1.0069.
- Completed and audited IOD/Mg replica 2 (119.3 ns/day), including the
  documented final-PDB correction. Its RDF peak is also 0.207 nm, mean shell
  distance 0.208636 nm, mean CN 6.0 and mean temperature 298.52 K. Six waters
  persisted in the shell at saved frames; the chloride-contact fraction was
  0.003. The final IOD/Mg replica and three IOD/Ca replicas remain pending.

- Completed and audited IOD/Mg replica 3 (138.7 ns/day), including the
  documented final-PDB correction. All six planned Mg trajectories are now
  complete. IOD/Mg RDF peaks are all 0.207 nm; its group mean shell distance
  is 0.208612 +/- 0.000032 nm and temperature 298.52 +/- 0.06 K (sample SDs).
  Compared with CM, the mean shell distance is larger by 0.004343 nm. Both
  models retained CN 6 and the same six water identities in each trajectory.
  These repeatable local structural differences do not establish water-exchange
  kinetics. Three IOD/Ca trajectories remain before the hydration quality gate.

- Completed and audited IOD/Ca replica 1 (115.2 ns/day), including the
  documented final-PDB correction. Its RDF peak is 0.243 nm, mean shell
  distance 0.247978 nm, mean CN 8.013 and temperature 298.46 K. CN 8 accounts
  for 97.5% of saved frames; 53 distinct shell-water identities and 98
  membership flips were recorded. No diagnostic chloride contacts occurred.
  The two remaining IOD/Ca replicas are needed before interpreting this group.

- Completed and audited IOD/Ca replica 2 (494.5 ns/day), including the
  documented final-PDB correction. Its RDF peak is 0.245 nm, mean shell
  distance 0.247859 nm, CN 7.9955 and temperature 298.23 K. CN 8 occupies
  97.15% of saved frames; 46 distinct shell-water identities and 82 membership
  flips were recorded. No diagnostic chloride contacts occurred. The faster
  throughput is an observed runtime, not a change to the numerical protocol.
- Revised report and presentation wording to make calibration-reference
  agreement conditional on actual results and describe between-replica
  variation without implying independently sampled initial hydration states.

- Completed and audited IOD/Ca replica 3 (508.4 ns/day); the hydration process
  exited successfully with all twelve planned trajectories. Its RDF peak is
  0.245 nm, mean shell distance 0.247961 nm and CN 7.999. The IOD/Ca group has
  a mean peak of 0.24433 +/- 0.00115 nm and CN 8.0025 +/- 0.00926. Final-PDB
  corrections are recorded for every affected old-runner case.
- Hydration quality gate passed: all saved values are finite; mean temperature
  spans 298.19-298.58 K, mean water density 32.9165-32.9542 per nm3, and outer
  RDF averages 0.9940-1.0090. Maximum absolute half-to-half changes are 0.571 K,
  0.0346 water oxygens per nm3, 22.93 kJ/mol potential energy and 0.027 water CN.
  No runaway behavior or source-conversion failure warrants holding the amide
  extension. A focused second review agreed. This gate does not establish
  equilibration of slow ligand or hydration-state changes.
- Generated and visually reviewed the four actual hydration-stage figures.
  The tiny Mg shell-cutoff effect is shown against an explicit coordination
  scale with a minimum reference-centered range, avoiding an easily missed
  additive-axis offset. Numerical results are unchanged.
- Applied the previously reviewed runner reliability patch only after the
  hydration process exited. Source/protocol/Git provenance is frozen once per
  invocation; initialized partial samples and metadata are saved atomically;
  failures are recorded; final PDB uses the final NPT box. System construction,
  periodic distances and LJ conversion are unchanged by AST comparison. The
  runnable self-check passed injected-save/failure checks and twelve analytic
  ion-water energy/force comparisons. Scientific settings remain frozen.

- Started the preplanned twenty-four-case amide campaign from clean commit
  39f0175. The first live partial snapshot contained 200 finite, initialized
  frames and its hash matched the metadata. The first complete CM/Mg/contact
  trajectory then passed the case audit, including final PDB dimensions and
  frozen invocation provenance, without a PDB-header correction.
  It ran at 441.3 ns/day; its diagnostic amide-contact fraction was 1.0, mean
  metal-carbonyl-O distance 0.200781 nm, water CN 5.0 and temperature 298.45 K.
  This describes a single contact-start trajectory; its separated counterpart
  and remaining replicas are still required for the planned comparison.

- The first CM/Mg/separated amide trajectory passed the case audit at
  452.3 ns/day. Its diagnostic contact fraction was 0.0, mean ion-amide-O
  distance 1.46408 nm, water CN 6.0 and temperature 298.43 K. Together with
  the first contact-start trajectory (fraction 1.0), this provides an initial
  example of finite-time preparation dependence. No affinity or equilibration
  claim is made; both remaining replica pairs are still planned.

- All six CM/Mg amide trajectories are complete and audited. In all three
  contact-start runs the diagnostic contact fraction is 1.0 and water CN is
  5.0; all three separated-start runs have contact fraction 0.0 and water CN
  6.0. Contact-start mean ion-amide-O distances span 0.200769-0.200929 nm;
  separated-start means span 1.37329-1.46408 nm. Temperatures span
  298.34-298.55 K, and diagnostic chloride-contact fractions are at most
  0.0015. This repeated starting-state contrast is retained as a finite-time
  sampling result. CM/Ca and both IOD groups remain in progress/pending.

- All six CM/Ca amide trajectories are complete and audited. Contact-start
  fractions are 0.9990, 0.8420 and 0.1485; all three separated-start fractions
  are 0.0. Contact-start mean water CN values are 7.0465, 7.2010 and 7.8965,
  versus 8.0355, 8.0435 and 8.0425 for separated starts. The large spread
  across contact-start trajectories is retained rather than reduced to a
  claim of a unique binding state. Mean temperatures span 298.17-298.76 K,
  and the largest diagnostic chloride-contact fraction is 0.0005.
  All twelve CM amide cases are now complete; twelve IOD cases remain.

- All six IOD/Mg amide trajectories are complete and audited. As with CM/Mg,
  contact-start fractions are 1.0 with water CN 5.0 in all three replicas;
  separated-start fractions are 0.0 with water CN 6.0. Contact-start mean
  ion-amide-O distances span 0.204965-0.205057 nm, compared with
  0.200769-0.200929 nm for CM/Mg. Mean temperatures span 298.29-298.53 K
  and diagnostic chloride contacts occupy at most 0.003 of saved frames.
  The final six IOD/Ca amide trajectories remain in progress/pending.

- All six IOD/Ca amide trajectories are complete and audited; the campaign
  exited successfully, completing all thirty-six formal cases and 72 ns of
  production. IOD/Ca contact-start fractions are 0.0050, 0.2210 and 0.3925;
  separated-start fractions are all 0.0. Contact-start mean water CN values
  are 7.9900, 7.7870 and 7.6015, versus 8.0040, 8.0000 and 8.0020 for
  separated starts. No diagnostic chloride contacts occurred in this group.
- Generated the combined analysis from all 36 cases. The report builder's
  mandatory data audit passed before generating the six-page English PDF and
  presentation. A focused scientific review found no new thermodynamic or
  counterion issue. Across amide trajectories, mean temperature spans
  298.174-298.758 K; the largest absolute half-to-half temperature change is
  0.614 K, and diagnostic chloride contacts occupy at most 0.3% of frames.
- The principal amide result is preparation dependence: all twelve separated
  starts have zero saved-frame contacts; all six Mg contact starts retain
  contact and five shell waters. Ca contact-start fractions average
  66.3 +/- 45.3% for CM and 20.6 +/- 19.4% for IOD (replica sample SD).
  Report and presentation now state these findings explicitly. Displayed
  shell/amide distance precision was increased to avoid rounding nonzero
  replica standard deviations to zero. Initial six-page visual review passed;
  the revised final PDF still requires its final render and audit.

## Assistance and ownership

OpenAI Codex assisted with literature/source retrieval, experimental design,
code generation, execution, verification, analysis and document preparation.
Subagents performed scoped parameter, analysis and code-review tasks in the
same project. These checks are not independent experimental replication or
supervisory approval. Literature values, simulation output and interpretation
must remain distinguishable in the final report.

The user owns the scientific decisions and should understand and review the
code, source references and results before presenting or submitting the work.
