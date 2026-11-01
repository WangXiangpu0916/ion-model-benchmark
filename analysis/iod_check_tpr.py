"""Compare the effective TPR LJ matrix to the source Amber matrix."""
from pathlib import Path
import re
import subprocess
import sys
import numpy as np
import parmed as pmd

amber, tpr = map(Path, sys.argv[1:])
s = pmd.load_file(str(amber))
dump = subprocess.check_output(['/opt/gromacs-2025.4-cuda/bin/gmx','dump','-s',str(tpr)],text=True,stderr=subprocess.DEVNULL)
lj = np.array(re.findall(r'functype\[\d+\]=LJ_SR, c6=\s*(\S+), c12=\s*(\S+)',dump),dtype=float)
n = s.ptr('ntypes')
assert n == 4 and lj.shape == (n*n,2)
# In these ion/water systems, Amber and GROMACS preserve first-appearance type order.
types = list(dict.fromkeys(a.nb_idx for a in s.atoms))
assert types == list(range(1,n+1))
expected=[]
for i in types:
    for j in types:
        k=s.parm_data['NONBONDED_PARM_INDEX'][n*(i-1)+j-1]-1
        expected.append([s.parm_data['LENNARD_JONES_BCOEF'][k]*4.184/10**6,
                         s.parm_data['LENNARD_JONES_ACOEF'][k]*4.184/10**12])
expected=np.array(expected)
assert np.allclose(lj,expected,rtol=2e-6,atol=1e-15)
assert 'functype[16]=SETTLE' in dump
mask=expected!=0
print('PASS: all 16 effective TPR LJ pairs match Amber A/B coefficients after unit conversion.')
print('Maximum relative error:',np.max(np.abs(lj[mask]/expected[mask]-1)))
print('PASS: effective SETTLE term present.')
print('Not an energy/force comparison between engines.')