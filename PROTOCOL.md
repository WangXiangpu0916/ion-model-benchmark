# Ion hydration and amide coordination protocol

Protocol frozen 5 October 2026 after an OpenCL throughput pilot, before hydration
production. This is a computational pilot study, not a preregistered clinical or
experimental study. The numerical settings and seeds are in `protocol.json`.

## Question and comparison

How do two published calibration choices for divalent ions change hydration
geometry and finite-time amide coordination under otherwise matched conditions?
The hydration phase compares Li-Merz 2013 CM and IOD 12-6 parameters in the same
TIP3P water model for Mg2+ and Ca2+. There are three seeded trajectories per ion
and model, twelve trajectories in total. RDF distance targets were used during
parameter fitting; matching them is an implementation and target-recovery check.

## Simulation budget

Each hydration system contains one divalent cation, two chloride counterions,
and water in an initially 3 nm cubic periodic box. Actual water count, fluctuating
volume, concentration and chloride contacts are reported. No extra salt is added.
OpenMM PME electrostatics, a 1 nm real-space cutoff, rigid TIP3P water, constrained
hydrogen bonds, a 2 fs Langevin-middle integrator, and isotropic 1 bar pressure
control are used at 298.15 K. After minimization and 100 ps equilibration, each
trajectory produces 2000 ps of data at 1 ps intervals. DCD coordinates are saved
every 10 ps. Independent seeds control velocities, thermostat, barostat and
counterion placement; a common water template is used for initial packing.

The 50 ps production pilot ran at 392.6 ns/day on the local NVIDIA RTX 3050 Laptop
GPU using OpenCL mixed precision. The twelve hydration production trajectories
therefore require about 1.5 hours at that measured rate, plus setup/equilibration;
thermal throttling and competing workloads can change the timing. A per-case
one-hour safety ceiling catches an unexpectedly slow or stalled run. Completed
cases are reused only when their protocol and sample checksum match.

**Runtime update during production:** the 1.5-hour estimate above was extrapolated
from the short pilot. Sustained throughput initially fell to about 80-90 ns/day
and later recovered; timing must be taken from completed-case metadata. The
planned trajectory lengths and numerical protocol remain unchanged.

## Observables and interpretation

- Ion-water oxygen RDF, normalized by each frame's water density. The analysis
  radius must stay below half the shortest box edge.
- First RDF peak and conditional mean first-shell distance, reported separately.
- Coordination at fixed diagnostic cutoffs: Mg 0.30 nm, Ca 0.35 nm, with
  +/-0.02 nm sensitivity. These choices do not make source CN values exact tests.
- Per-replica coordination distributions, block/half comparisons, temperature,
  density and ion-chloride contacts. Show all three replicates and their spread.
- Saved-frame water-identity changes are descriptive. No exchange rate, lifetime,
  equilibrium binding constant or receptor function will be estimated.

Frames are correlated and are not independent replicates. Parameter alternatives
are scientific scenarios, not statistical replicates. We will not lengthen runs
until significance or select a model solely because its training target matches.
If a metric drifts, contradicts a control, or depends on initialization, report
that finding and restrict the corresponding conclusion. Structural stability does
not certify equilibration of slow hydration or ligand-exchange processes.

## Amide extension and stopping rule

The extension starts after validating source parameter conversion, finite energies,
solvent structure and stable thermodynamic behavior in the hydration runs.
FreeSolv acetamide uses published GAFF/AM1-BCC parameters. Format conversion must
reproduce source Amber energies and forces before it is used in water.

Each ion/model is simulated from a direct-contact and a separated start, each
with three seeds. Settings and duration otherwise match the hydration phase.
The amide is unrestrained in production. Initial separation is not an imposed
binding condition. The primary result is the ion-carbonyl-oxygen distance and its
sensitivity to the starting state, alongside hydration counts and counterion
contacts. Persistent separation of outcomes diagnoses finite-time sampling
dependence; it does not measure an equilibrium affinity. A neutral amide in water
is a component model, not an NMDA receptor pore or a calibrated ion-amide model.

## Reproducibility and deliverables

Each case records the source-code hash, Git commit, complete numerical protocol,
seed, platform properties, periodic system XML, coordinates, thermodynamic log
and sample checksum. Large raw samples, checkpoints and DCD trajectories are
retained locally but ignored by Git. Derived summaries and figures are versioned.
The final report must distinguish literature values, measured simulations and
recommendations, and acknowledge AI assistance in code, analysis and writing.
