# Parameter provenance and reference values

Retrieved 5 October 2026 (Australia/Sydney). Exact retrieval URLs, repository commits, archive members, sizes, and SHA-256 checksums are recorded in `parameters/source_manifest.json`. Files under `parameters/amber/` and `parameters/freesolv/` are unchanged upstream inputs. Only the named acetamide members were retained from the FreeSolv archives (archive sizes 1,036,087 and 247,279 bytes).

## Ion models: same TIP3P water, different fitting objectives

**Primary parameter paper:** Li, P.; Roberts, B. P.; Chakravorty, D. K.; Merz, K. M., Jr. *Rational Design of Particle Mesh Ewald Compatible Lennard-Jones Parameters for +2 Metal Cations in Explicit Solvent.* Journal of Chemical Theory and Computation **9** (2013), 2733–2748. [DOI: 10.1021/ct400146w](https://doi.org/10.1021/ct400146w); [PubMed](https://pubmed.ncbi.nlm.nih.gov/23914143/); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3728907/).

The benchmark compares the published **compromise (CM)** and **ion–oxygen-distance (IOD)** 12–6 Lennard-Jones parameter sets. CM is the standard divalent set included in OpenMM `amber14/tip3p.xml`; IOD is fitted to the hydration-shell distance. Both use an ionic charge of **+2 e**, standard TIP3P water, Lorentz–Berthelot mixing, and PME electrostatics. No charge scaling or 12–6–4 term is introduced. The comparison tests sensitivity to the calibration objective; it does not establish that either set is best for every property.

The source values are in the Amber NONBON convention: atomic **Rmin/2 in angstroms**, epsilon in **kcal/mol**. Convert to OpenMM via `sigma_nm = 0.1 * 2 * Rmin_half_A / 2**(1/6)` and `epsilon_kJmol = 4.184 * epsilon_kcalmol`. Lorentz–Berthelot uses the arithmetic mean of atomic sigma and geometric mean of atomic epsilon.

| Set | Ion | Rmin/2 (Å) | epsilon (kcal/mol) | sigma (nm) | epsilon (kJ/mol) |
| --- | --- | ---: | ---: | ---: | ---: |
| CM | Mg2+ | 1.360 | 0.01020237 | 0.2423244513341723 | 0.04268671608 |
| CM | Ca2+ | 1.649 | 0.10592870 | 0.2938183972426839 | 0.4432056808 |
| IOD | Mg2+ | 1.395 | 0.01491700 | 0.2485607423611547 | 0.062412728 |
| IOD | Ca2+ | 1.608 | 0.08337961 | 0.2865130277539331 | 0.34886028824 |

Sources: [CM Amber parameter file](https://raw.githubusercontent.com/choderalab/ambermini/cfddd43ea945221350b5a4c31b2ba56d25958f6e/share/amber/dat/leap/parm/frcmod.ions234lm_126_tip3p), [IOD Amber parameter file](https://raw.githubusercontent.com/choderalab/ambermini/cfddd43ea945221350b5a4c31b2ba56d25958f6e/share/amber/dat/leap/parm/frcmod.ions234lm_iod_tip3p), and [OpenMM source for the native baseline](https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/data/amber14/tip3p.xml). Verify the installed OpenMM values in the runtime rather than depending on the moving upstream branch. The complete Amber files also contain ions outside this project's scope.

## Hydration reference values

The following are the target values compiled in Table 1 of Li, Z.; Song, L. F.; Li, P.; Merz, K. M., Jr. *Systematic Parametrization of Divalent Metal Ions for the OPC3, OPC, TIP3P-FB, and TIP4P-FB Water Models.* JCTC **16** (2020), 4429–4442. [DOI: 10.1021/acs.jctc.0c00194](https://doi.org/10.1021/acs.jctc.0c00194); [Table 1 and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8173364/). That paper attributes these targets to the earlier 2013 study and its underlying experimental references.

| Ion | First ion–water-oxygen RDF peak | First-shell coordination number |
| --- | --- | ---: |
| Mg2+ | 2.09 ± 0.04 Å | 6 |
| Ca2+ | 2.46 Å | 8 |

The Mg interval is the source's quoted uncertainty, not a measured uncertainty from this project or necessarily a standard error. The table provides no uncertainty for the Ca distance or either coordination number; do not invent one or treat those values as exact universal constants. Coordination depends on the operational shell definition and the experimental interpretation. Report the cutoff used, its sensitivity, and replica variation. IOD peak agreement is recovery of a fitting target, not an independent validation of the force field. Bulk hydration geometry also does not establish a binding free energy, water-exchange rate, or NMDAR selectivity.

## Acetamide extension

Use FreeSolv compound **mobley_8048190**, acetamide, SMILES `CC(=O)N`: nine atoms, neutral total charge. The supplied vacuum `prmtop` and `inpcrd` contain the complete preparameterized molecule; the `mol2` records atom names, coordinates, GAFF types, and AM1-BCC charges. The `frcmod` is a supplement, not a complete force field on its own. Preserve the original atom order, exclusions, torsions, and 1–4 scaling when converting formats. Amber coordinates are angstroms and Amber topology charge fields use Amber's internal scaling; use OpenMM's Amber reader rather than interpreting raw charge fields as electron charges.

FreeSolv's published calculations use **GAFF, AM1-BCC charges, and TIP3P water**. This matches the small-molecule/water conventions needed here without constructing charges from peptide caps. Source: [FreeSolv repository and methods](https://github.com/MobleyLab/FreeSolv/tree/6c7d19b4b565537365ffd22006aa2cd4643200c6), [database entry](https://raw.githubusercontent.com/MobleyLab/FreeSolv/6c7d19b4b565537365ffd22006aa2cd4643200c6/database.txt). FreeSolv lists acetamide hydration free energy as −9.71 kcal/mol experimentally and −8.82 ± 0.02 kcal/mol for GAFF. Its experimental uncertainty of 0.60 kcal/mol is a database-assigned default, not a reported experimental error. These are provenance/context values, not results of this project's simulations.

Cite: Mobley, D. L.; Guthrie, J. P. *FreeSolv: a database of experimental and calculated hydration free energies, with input files.* JCAMD **28** (2014), 711–720. [DOI: 10.1007/s10822-014-9747-x](https://doi.org/10.1007/s10822-014-9747-x). For the supplied calculation inputs, also cite Duarte Ramos Matos, G. et al. *Approaches for Calculating Solvation Free Energies and Enthalpies Demonstrated with an Update of the FreeSolv Database.* JCED **62** (2017), 1559–1569. [DOI: 10.1021/acs.jced.7b00104](https://doi.org/10.1021/acs.jced.7b00104). The upstream CC BY 4.0 license is retained in `parameters/freesolv/LICENSE`.

This extension is a transferability and sampling test of a neutral amide, not a separately calibrated metal–acetamide binding model. The source validates molecule hydration, not its Mg2+/Ca2+ affinity. Do not infer equilibrium preference merely because direct-contact and water-separated starts remain different over short trajectories.
