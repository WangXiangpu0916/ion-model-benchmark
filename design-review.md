# Scientific design review

**Status: advisory, 5 October 2026.** This is an independent review of the proposed benchmark, not a frozen protocol or a report of simulation results. Model identities, implementation, reference targets and the measured computing budget must be fixed before production.

## Question and defensible endpoint

Compare two published, internally compatible ion/water parameter combinations for Ca2+ and Mg2+, with three separately prepared replicas per combination: 12 small aqueous systems. The primary endpoint is the reproducibility of local hydration structure, not ion exchange rates or NMDAR function. Differences between combinations using different waters belong to the entire combination; they cannot be attributed uniquely to the ion potential.

The strongest attainable claim is: **under the stated simulation conditions and budget, these combinations reproduce, differ from, or do not resolve selected published hydration observables; the report identifies which additional sampling or parameter validation is needed before an amide or pore calculation.** A stable hydration shell alone cannot establish binding affinity, selectivity, or knock-on.

## Implementation gate before spending the budget

- Record parameter source/table/version, water variant, charge and mass distribution, pair mixing rules, all extra interactions, restraints, boundary conditions, electrostatics, thermostat, barostat and time step. Independently evaluate a few ion-water pair energies/forces to check units and engine implementation. A short run without numerical failure is not this check.
- If using Duarte's dummy model, its published A/B coefficients use geometric mixing. Converting them to per-particle Lennard-Jones parameters and silently using Lorentz-Berthelot can change the model. Its original solvent sphere/SCAAS/LRF protocol also differs from PBC/PME; an implementation with changed conditions is a transfer benchmark, not an exact replication. Source definitions of metal-oxygen distance differ between its Mg and Ca reporting. [Duarte et al., 2014, sections II.2 and III.1; Table 4 footnote](https://pmc.ncbi.nlm.nih.gov/articles/PMC4180081/).
- If using a 12-6-4 model, audit the actual ion-water C4 interaction and its cutoff treatment, not merely the presence of a C4 parameter in a file. Preserve the published water-specific parameter pairing. [Li et al., 2020, Table 5](https://pmc.ncbi.nlm.nih.gov/articles/PMC8173364/).
- Fix the ion concentration, box size and counterion/background-charge convention. If counterions are present, track their contact with the central ion. Salt association must not silently become a difference in the hydration comparison.

## Minimum measurements

| Observable | Exact treatment | Interpretation |
|---|---|---|
| Metal-water oxygen RDF, g(r) | Minimum-image distances to water oxygens only; exclude dummy/virtual sites. Normalize shell counts by shell volume and water number density, consistently with the ensemble. Stay below half the shortest box length. Report bin width. | Plot each replica as well as the mean; compare first-peak position and the complete first-shell shape. Peak height alone is sensitive to binning. |
| Coordination number, N(t) | Count waters within a declared first-shell cutoff. Also integrate the RDF and check agreement with the direct count. Report mean and the full discrete distribution P(N). | An average can conceal a mixture of coordination states. A frozen count is not proof that their equilibrium populations are sampled. |
| First-shell distances | Record the per-frame shell-distance mean and RDF peak separately. Match the chosen published target's definition. | A distribution's mean and mode are different observables. Neither is automatically an experimental metal-oxygen distance. |
| Sampling diagnostics | Temperature, density/volume, potential energy, N(t), cumulative distance mean, and first-shell water identities. | Thermodynamic stability and local geometric stability answer different questions from shell renewal. |

Choose the cutoff from the source or a separately labelled pilot, then freeze it. For comparisons, report an ion-specific common cutoff plus a small documented cutoff sensitivity; if model-specific first minima are used, show those too. Do not select cutoffs separately to force equal coordination numbers. Verify that the RDF's distant region approaches its expected normalization where the finite box permits it.

Do not add diffusion, absolute hydration free energy or a pore PMF to the minimum deliverable. Those introduce separate size, electrostatic, sampling and reference-state issues. Raw mean interaction energy is not hydration free energy. The selected hydration metrics may have been used to fit the models: reproducing them is an implementation/transfer check, not independent predictive validation.

## Run and convergence decisions

1. **Pilot separately.** Measure throughput and output/storage costs; check whether time step and constraints are numerically appropriate. Use independently solvated/relaxed systems and independent velocities, not only three seeds on one trapped shell. Record any shared initial geometry.
2. **Freeze the production budget and target precision.** Use equal planned production lengths within each comparison. Set meaningful precision in distance and coordination units using the reference's resolution/uncertainty or the model difference the project aims to resolve. If those sources do not support a numeric tolerance, report estimated precision rather than inventing a pass threshold.
3. **Determine equilibration from relevant traces.** Inspect density/energy and local hydration jointly. Examine sensitivity to removing a longer initial interval. Retain the full trace so trimming can be audited; do not trim solely to obtain agreement. Correlation-based equilibration selection can mistake a final trapped state for an equilibrium ensemble. [Grossfield et al., uncertainty and sampling guidance](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067).
4. **Assess each observable, not the trajectory as a whole.** Compare early/late production blocks, cumulative estimates and replicas. Estimate correlation time and block uncertainty for fluctuating observables, increasing block length until the uncertainty stabilizes while enough blocks remain for an estimate. No plateau, few effective samples, continuing drift or dependence on the equilibration cut means unresolved precision. A constant coordination trace does not have a meaningfully estimated correlation time.
5. **Separate three verdicts.** Implementation verified; local structural estimates sufficiently resolved for the declared comparison; kinetic/state-population sampling unresolved. These can coexist. Neither overlapping error bars nor a change below an arbitrary 5% establishes convergence or equivalence. If the budget ends first, publish the unresolved comparison and its diagnostics rather than extending until a desired ranking appears.

## Uncertainty with three replicas

Treat each independently prepared trajectory as the replication unit. Show its estimate, then the equally weighted mean and between-replica standard deviation. Within-trajectory blocks diagnose autocorrelation and drift; they do not turn three replicas into hundreds of independent systems. This follows the distinction between correlated frames and independent measurements in [Grossfield et al.](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067).

If a confidence interval is needed for a scalar mean, an approximate normal-replica-means t interval is mean +/- 4.303 SD/sqrt(3), with two degrees of freedom. State its assumptions and small-sample instability. Display the three values regardless; avoid a p-value-based model ranking or a seemingly precise bootstrap interval obtained by resampling individual frames. For RDFs and P(N), show separate replica curves; any shaded replica range is a range, not a 95% confidence band. Retain both between-replica variation and within-replica time dependence in the report.

Reference uncertainty, simulation sampling uncertainty and model discrepancy are distinct. A published range without a stated statistical meaning must not be treated as a confidence interval. Report signed deviation and method compatibility; call an uncertain comparison unresolved.

## Water exchange: an optional descriptive diagnostic

Track water identities, not just changes in N. A water can exchange while coordination number stays constant. Define entry/exit and a recrossing/commitment rule before counting events; give the distance cutoffs, commitment duration and observation cadence. Test event detection on known synthetic entries/exits and on a more finely sampled pilot. An event shorter than the output interval may be missed. Continuous survival and intermittent occupancy correlations answer different questions and must not share the label "residence time" without a definition.

Report completed exchanges, trajectories containing them, exposure time and censored episodes. Initial-shell residence episodes are left-censored and surviving final episodes are right-censored; completed lifetimes alone are biased toward short events. The six waters around one Mg are not six independent simulations. With no exchanges, report "none observed in this duration" and leave a fitted exchange rate/lifetime unavailable. Do not report zero rate, infinite lifetime or use lack of exchange as proof of an incorrect model. Grotz and Schwierz explicitly distinguish experimentally matched microsecond exchange from deliberately accelerated exchange models. [Grotz and Schwierz, 2022](https://doi.org/10.1063/5.0087292).

## Conditional amide extension

Proceed only after the implementation and local-hydration checks are interpretable, and only if compatible published amide parameters and a meaningful external comparison are available. Use a declared neutral small-molecule analogue such as acetamide; it represents an isolated Asn-like chemical group, not an Asn cage. Ion-amide cross interactions require their own audit, especially extra polarization/C4 terms. Do not invent a missing cross parameter.

Keep this a contact/solvent-separated initial-condition test at a declared concentration. Compare the metal-carbonyl oxygen distance, water and amide coordination separately, and starting-state dependence. A bound start remaining bound while a separated start stays separated reveals hysteresis, not two independently established equilibrium populations. Without repeated interconversion or a properly converged enhanced-sampling calculation, no binding free energy, affinity ranking or residence-rate claim follows. Passing the hydration check does not validate amide interactions.

## Compact deliverable

Preserve input/provenance metadata, a parameter/condition comparison table, all replica scalar traces and aggregate hydration data. Produce four figures: replica RDFs with compatible references; P(N) and scalar replica estimates; time/block diagnostics; observed exchange/initial-state limitations. The final short report should distinguish actual results from unmet endpoints and state the next discriminating calculation. A reproducible bounded negative or unresolved result is scientifically useful.
