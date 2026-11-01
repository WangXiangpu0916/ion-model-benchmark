"""Independent minimum-image counts/distances on 11 GROMACS-exported frames."""
from pathlib import Path
import re
import sys
import numpy as np

run=Path(sys.argv[1]).resolve();ion=sys.argv[2];cutoffs=np.array(list(map(float,sys.argv[3:])))
assert ion in ('mg','ca') and len(cutoffs)==3
folder=run/'analysis'
lines=(folder/'crosscheck_samples.gro').read_text().splitlines()
counts=np.loadtxt(folder/f'{ion}_cn_time.xvg',comments=('#','@'))
chlorides=np.loadtxt(folder/f'{ion}_cl_distances.xvg',comments=('#','@'))
ids={float(p[0]):set(map(int,p[2:])) for line in (folder/f'{ion}_shell_ids.dat').read_text().splitlines() if (p:=line.split())}
pos=0;times=[];errors=[]
while pos<len(lines):
    time=float(re.search(r't=\s*([\d.]+)',lines[pos])[1]);n=int(lines[pos+1]);assert n==2580
    atoms=lines[pos+2:pos+2+n]
    xyz=np.array([[float(v) for v in line[20:].split()[:3]] for line in atoms])
    box=np.array(list(map(float,lines[pos+2+n].split())))
    assert box.shape==(3,) and np.isfinite(xyz).all() and (box>0).all()
    assert atoms[0][5:10].strip()==ion.upper()
    oxy=np.array([i for i,line in enumerate(atoms) if line[5:10].strip()=='WAT' and line[10:15].strip()=='O'])
    cl=np.array([i for i,line in enumerate(atoms) if line[5:10].strip()=='Cl-'])
    assert len(oxy)==859 and len(cl)==2
    d=xyz-xyz[0];d-=box*np.rint(d/box);dist=np.linalg.norm(d,axis=1)
    row=counts[np.isclose(counts[:,0],time)][0]
    assert np.array_equal([(dist[oxy]<r).sum() for r in cutoffs],row[1:])
    assert set((oxy[dist[oxy]<cutoffs[1]]+1).tolist())==ids[time]
    reference=chlorides[np.isclose(chlorides[:,0],time)][0,1:]
    errors.extend(abs(dist[cl]-reference))
    times.append(time);pos+=n+3
assert times==list(np.arange(0,2001,200,dtype=float))
# GRO/XVG decimal output rounds distances; 0.001 nm tolerance covers that rounding.
assert max(errors)<.001
print('PASS: 11 frames (0, 200, ..., 2000 ps); all three CN cutoffs and shell IDs agree.')
print('PASS: independent minimum-image cation-Cl distances; max difference_nm:',max(errors))
print('Sources: GROMACS-exported GRO frames and original select/pairdist text results.')