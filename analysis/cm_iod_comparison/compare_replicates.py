"""Compare completed CM, IOD, or both sets, retaining every seed and cutoff."""
from pathlib import Path
import hashlib
import json
import re
import sys
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

project=Path(__file__).resolve().parents[2]
mode=sys.argv[1] if len(sys.argv)>1 else 'cm_iod'
assert mode in ('cm','iod','cm_iod')
models=('cm','iod') if mode=='cm_iod' else (mode,)
folder=project/'analysis'/f'{mode}_comparison'
label='CM and IOD' if mode=='cm_iod' else mode.upper()
rows=[];curves={}
for ion in ('mg','ca'):
    for model in models:
        for rep in (1,2,3):
            run=project/'systems'/f'{ion}_{model}'
            run=run/'gromacs' if rep==1 else run/f'rep{rep:02d}'/'gromacs'
            a=run/'analysis'
            def load(name,shape):
                x=np.loadtxt(a/name,comments=('#','@'))
                assert x.shape==shape and np.isfinite(x).all(),(run,name)
                return x
            rdf=load(f'{ion}_ow_rdf.xvg',(600,2));cn=load(f'{ion}_cn_time.xvg',(2001,4))
            thermo=load('md_thermo.xvg',(20001,6));cl=load(f'{ion}_cl_distances.xvg',(2001,3))
            assert np.allclose(cn[:,0],np.arange(2001)) and np.allclose(thermo[:,0],np.arange(20001)*.1)
            assert 'Finished mdrun' in (run/'md.log').read_text()
            cutoffs=list(map(float,re.findall(r'within\s+([\d.]+)',(a/f'{ion}_cn_selections.txt').read_text())))
            assert len(cutoffs)==3
            shell=[]
            for i,line in enumerate((a/f'{ion}_shell_ids.dat').read_text().splitlines()):
                p=line.split();s=frozenset(map(int,p[2:]));assert float(p[0])==i and len(s)==int(p[1])==cn[i,2]
                shell.append(s)
            assert len(shell)==2001
            first=rdf[(rdf[:,0]>=.15)&(rdf[:,0]<=cutoffs[1])]
            with (run/'md.xtc').open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            row={'ion':ion,'model':model,'replicate':rep,'run':str(run),'trajectory_sha256':digest,
                 'rdf_peak_nm':float(first[first[:,1].argmax(),0]),'cutoffs_nm':cutoffs,
                 'mean_cn':cn[:,1:].mean(axis=0).tolist(),
                 'cn_distributions':[{str(int(k)):int(v) for k,v in sorted(Counter(cn[:,j]).items())} for j in (1,2,3)],
                 'central_cn_500ps_blocks':[float(cn[(cn[:,0]>=b)&((cn[:,0]<b+500) if b<1500 else (cn[:,0]<=2000)),2].mean()) for b in (0,500,1000,1500)],
                 'temperature_K':float(thermo[:,2].mean()),
                 'density_g_ml':float(thermo[:,5].mean()/1000),
                 'unique_shell_waters':len(set().union(*shell)),
                 'adjacent_set_changes':sum(x!=y for x,y in zip(shell,shell[1:])),
                 'chloride_minima_nm':cl[:,1:].min(axis=0).tolist()}
            rows.append(row);curves[ion,model,rep]=rdf
assert len(rows)==6*len(models) and len({r['trajectory_sha256'] for r in rows})==len(rows)
# All source data must be complete before creating any result.
outputs=['replicate_statistics.json',f'{mode}_rdf_replicates.png',f'{mode}_cn_replicates.png',f'{mode}_cn_distributions.png',f'{mode}_comparison_summary.txt']
folder.mkdir(exist_ok=True)
for name in outputs:assert not (folder/name).exists(),name
fig,axes=plt.subplots(2,2,figsize=(12,8))
for row,ion in enumerate(('mg','ca')):
    for col in (0,1):
        ax=axes[row,col]
        for model in models:
            color={'cm':'navy','iod':'darkorange'}[model]
            for rep,style in [(1,'-'),(2,'--'),(3,':')]:
                x=curves[ion,model,rep];ax.plot(x[:,0],x[:,1],color=color,linestyle=style,linewidth=1.3,label=f'{model.upper()} {rep:02d}')
        ax.set_xlim((.15,.6) if col==0 else ((.19,.225) if ion=='mg' else (.22,.28)))
        ax.set_xlabel(f'{ion.capitalize()}-O distance (nm)');ax.set_ylabel('RDF g(r)')
        ax.set_title(f'{ion.capitalize()}: '+('hydration shells' if col==0 else 'first peak'))
        ax.grid(alpha=.2);ax.legend(ncol=2,fontsize=8)
fig.suptitle(f'{label} hydration: three 2 ns trajectories per condition')
fig.tight_layout();fig.savefig(folder/f'{mode}_rdf_replicates.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4.5))
for ion,ax in zip(('mg','ca'),axes):
    subset=[r for r in rows if r['ion']==ion]
    for i,r in enumerate(subset):
        for j,(offset,color) in enumerate(zip((-.18,0,.18),('navy','seagreen','darkorange'))):
            series_label=f"{r['model'].upper()} {r['cutoffs_nm'][j]:.2f} nm" if r['replicate']==1 else None
            ax.scatter(i+offset,r['mean_cn'][j],color=color,marker='o' if r['model']=='cm' else 's',s=35,label=series_label)
    ax.set_xticks(range(len(subset)),[f"{r['model'].upper()} {r['replicate']:02d}" for r in subset],rotation=30)
    ax.set_ylabel('Mean water coordination number');ax.set_title(ion.capitalize())
    ax.grid(axis='y',alpha=.2);ax.legend(ncol=2,fontsize=8)
fig.suptitle('Coordination by run and shell boundary')
fig.tight_layout();fig.savefig(folder/f'{mode}_cn_replicates.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4.8))
for ion,ax in zip(('mg','ca'),axes):
    subset=[r for r in rows if r['ion']==ion]
    states=sorted({int(k) for r in subset for k in r['cn_distributions'][1]})
    bottom=np.zeros(len(subset))
    for state in states:
        values=np.array([r['cn_distributions'][1].get(str(state),0)/2001 for r in subset])
        ax.bar(range(len(subset)),values,bottom=bottom,label=f'CN {state}')
        bottom+=values
    assert np.allclose(bottom,1)
    ax.set_xticks(range(len(subset)),[f"{r['model'].upper()} {r['replicate']:02d}" for r in subset],rotation=30)
    ax.set_ylim(0,1.05);ax.set_ylabel('Fraction of saved frames')
    ax.set_title(f"{ion.capitalize()}: {subset[0]['cutoffs_nm'][1]:.2f} nm shell");ax.legend()
fig.suptitle(f'{label}: coordination-number distributions by run')
fig.tight_layout();fig.savefig(folder/f'{mode}_cn_distributions.png',dpi=180);plt.close(fig)
lines=[f'{label} hydration benchmark: per-replicate comparison','Record date: 2026-11-02',f'Project root: {project}',
       f'Analysis directory: {folder}',f'Status: all {len(rows)} {label} 2 ns trajectories and per-run analyses complete.','',
       'SYSTEM AND METHODS',
       f'1 divalent cation, 2 Cl-, 859 standard Amber TIP3P waters; Li-Merz 2013 {label} 12-6 cations and Joung-Cheatham Cl-.',
       'Each condition uses one minimized configuration and three distinct velocity/stochastic seeds; 20 ps NVT, 1000 ps NPT, then 2 ns NPT at 298 K and 1 bar.',
       'All replicas have 2001 trajectory frames and 20001 exported energy records. RDF bin width: 0.002 nm.',
       'Source: Li et al. 2013, DOI 10.1021/ct400146w, Tables 7-8. Parameter/conversion and sampled-frame crosschecks are stored separately.','',
       'OBSERVED RESULTS',
       'Condition Rep Peak_nm Cutoff_nm Mean_CN Temperature_K Density_g_ml Distinct_shell_waters Adjacent_set_changes']
for r in rows:lines.append(f"{r['ion']}_{r['model']} {r['replicate']:02d} {r['rdf_peak_nm']:.3f} {r['cutoffs_nm'][1]:.2f} {r['mean_cn'][1]:.6f} {r['temperature_K']:.3f} {r['density_g_ml']:.6f} {r['unique_shell_waters']} {r['adjacent_set_changes']}")
lines.extend(['','Per-condition statistics across three seeds (sample SD, not frame-based standard errors):'])
for ion in ('mg','ca'):
    for model in models:
        rs=[r for r in rows if r['ion']==ion and r['model']==model]
        peaks=[r['rdf_peak_nm'] for r in rs];means=[r['mean_cn'][1] for r in rs]
        lines.append(f'{ion}_{model}: peak {min(peaks):.3f}-{max(peaks):.3f} nm; central-cutoff mean CN {min(means):.6f}-{max(means):.6f}; CN mean +/- SD {np.mean(means):.6f} +/- {np.std(means,ddof=1):.6f}; peak mean +/- SD {np.mean(peaks):.6f} +/- {np.std(peaks,ddof=1):.6f} nm.')
lines.extend(['','Coordination distributions and boundary sensitivity (counts out of 2001):'])
for r in rows:
    for cutoff,mean,counts in zip(r['cutoffs_nm'],r['mean_cn'],r['cn_distributions']):
        lines.append(f"{r['ion']}_{r['model']} {r['replicate']:02d}: cutoff {cutoff:.2f} nm; mean {mean:.6f}; counts {counts}.")
    lines.append(f"  Central-cutoff CN in consecutive 500 ps blocks: {r['central_cn_500ps_blocks']}.")
lines.extend(['','The reported seed SD and range describe the three tested runs. RDF maxima are discretized at 0.002 nm; zero SD at this resolution is not zero uncertainty.',
              'The 500 ps block means describe time variation within each trajectory; the final block includes the 2000 ps endpoint. They are not additional independent replicas.'])
if mode=='cm_iod':
    lines.extend(['','Parameter sensitivity: IOD minus CM (difference of three-seed arithmetic means)'])
    for ion in ('mg','ca'):
        cm=[r for r in rows if r['ion']==ion and r['model']=='cm']
        iod=[r for r in rows if r['ion']==ion and r['model']=='iod']
        assert len(cm)==len(iod)==3
        assert all(r['cutoffs_nm']==cm[0]['cutoffs_nm'] for r in cm+iod)
        delta_peak=np.mean([r['rdf_peak_nm'] for r in iod])-np.mean([r['rdf_peak_nm'] for r in cm])
        delta_cn=np.mean([r['mean_cn'] for r in iod],axis=0)-np.mean([r['mean_cn'] for r in cm],axis=0)
        lines.append(f'{ion}: RDF peak shift {delta_peak:+.6f} nm.')
        for cutoff,delta in zip(cm[0]['cutoffs_nm'],delta_cn):
            lines.append(f'  Mean coordination change at {cutoff:.2f} nm: {delta:+.6f}.')
    lines.append('These differences compare the sampled local hydration structures under the two cation parameter sets with the same water/Cl models and protocol; seed-level SD and cutoff dependence are reported above.')
lines.extend(['','LIMITATIONS AND REMAINING WORK',
              'Seeds within each condition share the minimized configuration. Local structural agreement does not exclude initial-state dependence or demonstrate convergence of slow processes.',
              'Adjacent saved frames are correlated; time halves are not independent replicas. Three seeds do not provide precise uncertainty for rare events.',
              'Water-ID set changes are not independent exchange events or exchange rates. A 1 ps saving interval can miss brief departures and returns.',
              'CM/IOD differences describe model dependence; recovering a fitted IOD is not independent predictive validation.',
              'No dedicated ion-pair classification was applied. Water-shell cutoffs are not ion-Cl contact criteria.',
              'A full cross-engine periodic single-point energy/force comparison remains unfinished.',
              'Water-box results alone do not establish NMDAR binding, permeation or block; structure evidence is evaluated separately.'])
with (folder/'replicate_statistics.json').open('x') as f:json.dump(rows,f,indent=2)
with (folder/f'{mode}_comparison_summary.txt').open('x',encoding='utf-8') as f:f.write('\n'.join(lines)+'\n')
print('\n'.join(lines[12:]))
