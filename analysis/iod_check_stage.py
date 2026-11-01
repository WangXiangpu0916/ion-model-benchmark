"""Check one completed IOD dynamics stage using text data only; never opens XTC."""
from pathlib import Path
import json
import re
import sys
import numpy as np

run = Path(sys.argv[1]).resolve()
stage, label, end = sys.argv[2], sys.argv[3], float(sys.argv[4])
assert stage in ('nvt', 'npt', 'md')
xvg = run / 'analysis' / f'{label}_thermo.xvg'
data = np.loadtxt(xvg, comments=('#', '@'))
legends = re.findall(r'@ s\d+ legend "([^"]+)"', xvg.read_text())
assert legends == (['Potential', 'Temperature', 'Pressure'] if stage == 'nvt' else
                   ['Potential', 'Temperature', 'Pressure', 'Volume', 'Density'])
assert data.shape == (round(end*10)+1, len(legends)+1)
assert np.isfinite(data).all() and np.allclose(np.diff(data[:,0]), .1)
assert np.allclose(data[[0,-1],0], [0,end])
log = (run / f'{stage}.log').read_text()
assert 'Finished mdrun' in log
assert not re.search(r'LINCS WARNING|can not be settled|Fatal error', log, re.I)
assert re.search(rf'^\s*{round(end*1000)}\s+{end:.5f}\s*$', log, re.M)
lines = (run / f'{stage}.gro').read_text().splitlines()
assert int(lines[1]) == 2580
coords = np.array([[float(v) for v in line[20:].split()] for line in lines[2:-1]])
box = np.array([float(v) for v in lines[-1].split()])
assert coords.shape == (2580,6) and np.isfinite(coords).all()
assert box.shape == (3,) and np.isfinite(box).all() and (box>2).all()
means = dict(zip(legends, data[:,1:].mean(axis=0)))
late = dict(zip(legends, data[data[:,0]>=end/2,1:].mean(axis=0)))
blocks = []
if stage != 'nvt':
    for start in np.arange(0,end,end/5):
        rows=data[(data[:,0]>=start)&(data[:,0]<start+end/5)]
        blocks.append({'start_ps':float(start),'end_ps':float(start+end/5),
                       'temperature_K':float(rows[:,2].mean()),
                       'density_g_ml':float(rows[:,5].mean()/1000)})
report={'run':str(run),'stage':stage,'label':label,'end_ps':end,'records':len(data),
        'means':means,'last_half_means':late,'fifths':blocks,
        'box_nm':box.tolist(),'numerical_warning_patterns':0,
        'wall_seconds':float(re.findall(r'Time:\s+\S+\s+(\S+)',log)[-1])}
with (run/'analysis'/f'{label}_check.json').open('x') as f:json.dump(report,f,indent=2)
print(json.dumps(report,indent=2))