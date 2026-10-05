# Five-minute presentation script

Xiangpu Wang | 2026-10-06

Suggested visuals: report pages 1, 2, 3, 4 and 5. The figures show actual saved simulations. Speak the numerical contrasts as descriptive findings; do not add a significance or affinity claim.

## 0:00-0:50 - The question

My white paper asked how calcium could modify magnesium block in NMDA receptors. I chose a smaller question that I could finish on my own computer: how much do published ion calibration choices change hydration and contact with a neutral amide?

Acetamide provides an isolated amide group related to asparagine chemistry. It does not reproduce an asparagine cage. The purpose was to examine model behaviour and sampling requirements before making a pore-mechanism claim.

## 0:50-1:45 - What I ran and checked

I compared the Li-Merz compromise, or CM, and ion-oxygen-distance, or IOD, twelve-six parameter sets. Both use full plus-two charges and the same TIP3P water and interaction conventions. Each box contained one divalent ion and two chlorides, at approximately 0.063 molar cation concentration from actual volumes, rather than infinite dilution.

There were twelve hydration simulations: two ions, two models and three replicas. The solvent-packing template is shared, with separately seeded velocities, barostat moves and chloride placement. I then ran twenty-four acetamide simulations, adding direct-contact and separated starting preparations. These produced 72.0 nanoseconds of recorded production. The acetamide conversion was checked against the original Amber topology using energies and forces in five geometries; a successful file conversion alone would not establish equivalent physics.

## 1:45-2:40 - Hydration results

For magnesium, the mean first RDF peak was 0.2030 nanometres with CM and 0.2070 with IOD. For calcium, the values were 0.2483 and 0.2443. I show every replica, the group mean and standard deviation. Full distributions are available as machine-readable results.

These reference distances helped define the model family. Any agreement with these calibration references constitutes target recovery, not independent predictive validation. The RDF peak and the mean distance inside a shell are also reported separately because they are different measurements.

## 2:40-3:55 - Sampling and the amide extension

I checked production halves, consecutive blocks, counterion contact and sensitivity to the shell boundary. The largest absolute hydration coordination change between halves was 0.027. Moving the cutoff by plus or minus 0.02 nanometres gave a largest coordination span of 0.043. These diagnostics do not provide a universal convergence certificate.

The amide trajectories showed a clear preparation effect: 12 of twelve separated starts had no saved-frame contacts. For magnesium, contact-start minus separated-start fractions differed by 100.0 percentage points with CM and 100.0 with IOD. Calcium contact-start fractions averaged 66.3 percent with CM and 20.6 with IOD, with large replica standard deviations of 45.3 and 19.4 percentage points. That variability is part of the result.

Persistent initial-state dependence is a sampling warning. Similar contact fractions would also be insufficient to prove equilibration. These finite-time fractions do not establish an affinity, free-energy difference or residence time.

## 3:55-5:00 - What the project delivers

The final deliverable is a local repository containing parameter provenance, runnable simulation and analysis code, the recorded protocol, quantitative figures and this short report. Actual production throughput averaged 262.4 nanoseconds per day across the completed cases. That is a measured property of this workflow; the earlier short pilot is excluded from the estimate.

The scientific outcome is a bounded comparison of hydration structure and amide sampling sensitivity. It does not establish a preferred NMDAR ion model, simultaneous calcium-magnesium occupancy, or a knock-on mechanism. My next step would be to establish reversible sampling and an external metal-amide benchmark before attempting a quantitative binding calculation or a pore model.

Codex assisted the source work, code and drafting, and that contribution is disclosed in the report. The numerical claims can be traced back to the saved trajectories and inspected independently. I would like to discuss whether this benchmark is a useful starting point for the group's ion-selectivity work, and which next observable would most strongly discriminate between models.

Script length: approximately 669 words including headings; rehearse near 134 words/minute and adjust pauses to five minutes.
