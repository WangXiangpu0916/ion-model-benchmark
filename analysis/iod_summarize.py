"""Check saved GROMACS analyses and write one English summary; no XTC reader."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import sys
import numpy as np

root=Path(__file__).resolve().parents[1]
run=Path(sys.argv[1]).resolve()
run.relative_to(root/'systems')
ion,rep=sys.argv[2],int(sys.argv[3]); cutoffs=list(map(float,sys.argv[4:]))
assert ion in ('mg','ca') and rep in (1,2,3)
assert len(cutoffs)==3 and cutoffs==sorted(set(cutoffs))
a=run/'analysis'; name=f'{ion}_iod_rep{rep:02d}'; label=ion.capitalize()

def load(filename,shape):
    x=np.loadtxt(a/filename,comments=('#','@'))
    assert x.shape==shape and np.isfinite(x).all(), filename
    assert np.all(np.diff(x[:,0])>0), filename
    return x

thermo=load('md_thermo.xvg',(20001,6))
rdf=load(f'{ion}_ow_rdf.xvg',(600,2)); cn=load(f'{ion}_ow_cn.xvg',(600,2))
counts=load(f'{ion}_cn_time.xvg',(2001,4)); cl=load(f'{ion}_cl_distances.xvg',(2001,3))
assert np.allclose(thermo[:,0],np.arange(20001)*.1)
assert np.allclose(counts[:,0],np.arange(2001)) and np.allclose(cl[:,0],counts[:,0])
assert np.allclose(counts[:,1:],np.rint(counts[:,1:]))
assert np.all(np.diff(counts[:,1:],axis=1)>=0)
assert np.all(np.diff(cn[:,1])>=-1e-6)
assert json.loads((a/'md_check.json').read_text())['end_ps']==2000
lines=(run/'md.gro').read_text().splitlines()[2:-1]
water_ids={int(l[15:20]) for l in lines if l[5:10].strip()=='WAT' and l[10:15].strip()=='O'}
assert len(water_ids)==859
shells=[]
for i,line in enumerate((a/f'{ion}_shell_ids.dat').read_text().splitlines()):
    parts=line.split();time=float(parts[0]);n=int(parts[1]);ids=frozenset(map(int,parts[2:]))
    assert time==i and n==len(ids)==counts[i,2] and ids<=water_ids
    shells.append(ids)
assert len(shells)==2001

def peak(x):
    first=x[(x[:,0]>=.15)&(x[:,0]<=cutoffs[1])]
    return float(first[np.argmax(first[:,1]),0])

seeds={}
for stage in ['nvt','npt','md']:
    seeds[stage]=int(re.search(r'^ld-seed\s*=\s*(\d+)',(run/f'{stage}.mdp').read_text(),re.M)[1])
with (run/'md.xtc').open('rb') as f: trajectory_hash=hashlib.file_digest(f,'sha256').hexdigest()
stats={'condition':f'{ion}_iod','replicate':rep,'run':str(run),'cutoffs_nm':cutoffs,
       'trajectory_sha256':trajectory_hash,'seeds':seeds,'rdf_peak_nm':peak(rdf),
       'temperature_K':float(thermo[:,2].mean()),'density_g_ml':float(thermo[:,5].mean()/1000),
       'first500_density_g_ml':float(thermo[thermo[:,0]<500,5].mean()/1000),
       'last500_density_g_ml':float(thermo[thermo[:,0]>=1500,5].mean()/1000),
       'coordination':[], 'unique_shell_waters':len(set.union(*(set(s) for s in shells))),
       'distinct_shell_sets':len(set(shells)),
       'persistent_water_ids':sorted(set.intersection(*(set(s) for s in shells))),
       'adjacent_set_changes':sum(x!=y for x,y in zip(shells,shells[1:])),
       'same_cn_set_changes':sum(x!=y and len(x)==len(y) for x,y in zip(shells,shells[1:])),
       'half_peaks_nm':[peak(load(f'{ion}_ow_rdf_half{i}.xvg',(600,2))) for i in (1,2)],
       'chloride_minima':[{'distance_nm':float(cl[:,j].min()),'time_ps':float(cl[cl[:,j].argmin(),0])} for j in (1,2)]}
for j,r in enumerate(cutoffs,1):
    avg=float(counts[:,j].mean()); cumulative=float(np.interp(r,cn[:,0],cn[:,1]))
    assert abs(avg-cumulative)<.02
    stats['coordination'].append({'cutoff_nm':r,'mean':avg,'rdf_cumulative':cumulative,
                                 'distribution':{str(int(k)):int(v) for k,v in sorted(Counter(counts[:,j]).items())}})
coord='\n'.join(f"   {x['cutoff_nm']:.2f} nm: mean {x['mean']:.6f}; CN counts {x['distribution']}." for x in stats['coordination'])
cltext='\n'.join(f"   Cl{i}: minimum {x['distance_nm']:.3f} nm at {x['time_ps']:.0f} ps." for i,x in enumerate(stats['chloride_minima'],1))
summary=f'''{label} IOD hydration replicate {rep:02d}: summary
Record date: 2026-11-02
Project root: {root}
Run directory: {run.relative_to(root)}
Status: completed 2 ns IOD trajectory with checked basic hydration analyses.

SYSTEM AND METHODS
- Composition: 1 {label}2+, 2 Cl-, 859 waters; 2580 atoms; total charge 0.
- Models: Li-Merz IOD 12-6 cation, standard Amber TIP3P, Joung-Cheatham Cl-.
- Sources: Li et al. (2013), DOI 10.1021/ct400146w, Table 7; AmberTools 26 source parameter files data/raw/frcmod.ions234lm_iod_tip3p and frcmod.ionsjc_tip3p.
- IOD raw-file SHA-256: 081a17684c4a1589d93af3c45b196ff6ca182912e443887ae799bbe3b633e16e.
- Build: systems/{ion}_iod/amber/build_matched.leap and system_matched.prmtop/inpcrd, using the corresponding valid unminimized CM PDB as the coordinate template; only cation LJ parameters change.
- tleap build, ParmEd conversion, composition, atom order, charges, water/Cl parameters, box, all 16 LJ pairs, effective TPR coefficients, SETTLE and exclusions were checked.
- Each IOD condition has its own minimization. Repeats within a condition share that minimized configuration and use new velocities/seeds; they are not independent initial configurations.
- Engine: GROMACS 2025.4, /opt/gromacs-2025.4-cuda/bin/gmx; 1 MPI rank, 4 OpenMP threads, GPU nonbonded/PME for dynamics.
- Preparation: 20 ps NVT -> 200 ps NPT -> checkpoint continuation to a total 1000 ps NPT.
- Production: 2 ns NPT, 1 fs timestep, 298 K, 1 bar; production time resets to 0 ps.
- V-rescale thermostat (tau-t 0.1 ps), isotropic C-rescale barostat (tau-p 5.0 ps).
- PME electrostatics, 1.0 nm real-space cutoff; LJ 1.0 nm cutoff, Potential-shift, DispCorr=EnerPres; xyz periodic boundaries.
- Rigid water uses SETTLE even though the MDP constraints setting is none.
- NVT gen-seed/ld-seed: {seeds['nvt']}; NPT ld-seed: {seeds['npt']}; MD ld-seed: {seeds['md']}. NPT and MD inherit velocities.
- Production: 2001 frames, 0-2000 ps, every 1 ps; 20001 energy records, every 0.1 ps. Finished normally; no detected dynamics numerical warnings.
- Trajectory SHA-256: {trajectory_hash}.
- Analyses use GROMACS and text-only Python checks, RDF bin width 0.002 nm, periodic boundaries. No MDAnalysis trajectory cache was generated.

OBSERVED RESULTS
1. Thermodynamics (arithmetic means of exported XVG records)
   Temperature: {stats['temperature_K']:.3f} K; density: {stats['density_g_ml']:.6f} g/mL.
   Density in 0-<500 ps and 1500-2000 ps: {stats['first500_density_g_ml']:.6f} and {stats['last500_density_g_ml']:.6f} g/mL.
   GROMACS accumulated energy averages may differ slightly from these exported-record means.
2. {label}-water oxygen RDF
   First peak: {stats['rdf_peak_nm']:.3f} nm. Provisional shell boundary: {cutoffs[1]:.2f} nm.
   RDF and cumulative-coordination distance grids are offset; each uses its own distance column.
3. Per-frame coordination and boundary sensitivity
{coord}
4. Shell identities at {cutoffs[1]:.2f} nm
   Counts match per-frame coordination in all 2001 frames, and all IDs identify WAT/O atoms.
   Distinct waters seen: {stats['unique_shell_waters']}; distinct ID sets: {stats['distinct_shell_sets']}.
   Adjacent saved-frame set changes: {stats['adjacent_set_changes']}; changes with unchanged CN: {stats['same_cn_set_changes']}.
   IDs present in every saved frame: {stats['persistent_water_ids']}.
5. Time-segment RDF
   0-999 ps (1000 frames) peak: {stats['half_peaks_nm'][0]:.3f} nm.
   1000-2000 ps (1001 frames) peak: {stats['half_peaks_nm'][1]:.3f} nm.
6. {label}-Cl distances
{cltext}
7. Data and figures
   analysis/md_thermo.xvg/png; {ion}_ow_rdf.xvg; {ion}_ow_cn.xvg; {ion}_ow_rdf_cn.png;
   {ion}_cn_selections.txt; {ion}_cn_time.xvg; {ion}_shell_ids.dat;
   {ion}_ow_rdf_half1.xvg; {ion}_ow_rdf_half2.xvg; {ion}_ow_rdf_halves.png;
   {ion}_cl_distances.xvg/png; {name}_statistics.json. Commands and checked figures are recorded in the process document.

LIMITATIONS AND REMAINING WORK
- This summary describes one seed. Report each planned replicate separately; shared initial configurations do not exclude initial-state dependence.
- Similar RDFs, stable temperature, or completed repeats do not prove full convergence of slow processes. Time halves are not independent replicates and adjacent frames are correlated.
- Water-ID set changes are not independent exchange-event counts or exchange rates. No observed change at 1 ps spacing cannot exclude brief departures and returns between frames.
- No dedicated ion-pair criterion has been applied; the water-shell cutoff is not a cation-Cl contact criterion.
- CM/IOD comparisons assess model dependence. Recovery of the fitted IOD target is not independent predictive validation or proof that one model is universally better.
- No reliable exchange rate or robust between-replicate uncertainty is inferred here.
- Full cross-engine periodic-system single-point energy/force comparison remains unfinished.
- Water-box results do not establish NMDAR binding, permeation, or blockade mechanisms; published structures require separate source-based evaluation.
'''
for p in [a/f'{name}_statistics.json',a/f'{name}_summary.txt']:assert not p.exists(),p
with (a/f'{name}_statistics.json').open('x') as f:json.dump(stats,f,indent=2)
with (a/f'{name}_summary.txt').open('x',encoding='utf-8') as f:f.write(summary)
print(json.dumps(stats,indent=2))
print('Summary:',a/f'{name}_summary.txt')