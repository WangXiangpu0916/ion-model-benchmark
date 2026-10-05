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

## Assistance and ownership

OpenAI Codex assisted with literature/source retrieval, experimental design,
code generation, execution, verification, analysis and document preparation.
Subagents performed scoped parameter, analysis and code-review tasks in the
same project. These checks are not independent experimental replication or
supervisory approval. Literature values, simulation output and interpretation
must remain distinguishable in the final report.

The user owns the scientific decisions and should understand and review the
code, source references and results before presenting or submitting the work.
