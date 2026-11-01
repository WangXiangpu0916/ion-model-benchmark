from pathlib import Path
from collections import Counter
import math
import numpy as np
import parmed as pmd

out = Path("../gromacs")
assert not (out / "system.top").exists()
assert not (out / "system.gro").exists()
s = pmd.load_file("system_matched.prmtop", xyz="system_matched.inpcrd")
old = pmd.load_file("../../mg_cm/amber/system.prmtop",
                    xyz="../../mg_cm/amber/system.inpcrd")
assert len(s.atoms) == 2580
assert Counter(r.name for r in s.residues) == {"MG": 1, "Cl-": 2, "WAT": 859}
assert abs(sum(a.charge for a in s.atoms)) < 1e-8
assert [(a.name, a.type, a.residue.name) for a in s.atoms] == [
    (a.name, a.type, a.residue.name) for a in old.atoms]
for a, b in zip(s.atoms, old.atoms):
    assert math.isclose(a.charge, b.charge, abs_tol=1e-10)
    assert math.isclose(a.mass, b.mass, abs_tol=1e-10)
    if a.idx:
        assert math.isclose(a.rmin, b.rmin, abs_tol=1e-10)
        assert math.isclose(a.epsilon, b.epsilon, abs_tol=1e-10)
assert math.isclose(s.atoms[0].rmin, 1.395, abs_tol=1e-9)
assert math.isclose(s.atoms[0].epsilon, 0.014917, abs_tol=1e-9)
delta = s.coordinates - old.coordinates
# PDB rounds each coordinate to 0.001 A; both Amber coordinate sets use the same box origin.
relative_error = np.max(np.abs(delta))
assert relative_error <= 0.0005001
assert np.allclose(s.box, old.box, rtol=0, atol=1e-7)
s.save(str(out / "system.top"), format="gromacs", overwrite=False)
s.save(str(out / "system.gro"), format="gro", precision=6, overwrite=False)
g = pmd.load_file(str(out / "system.top"), xyz=str(out / "system.gro"))
assert len(g.atoms) == 2580
assert g.defaults.comb_rule == 2
assert np.max(np.abs(g.coordinates - s.coordinates)) < 0.000006
for a, b in zip(s.atoms, g.atoms):
    assert (a.name, a.residue.name) == (b.name, b.residue.name)
    assert math.isclose(a.charge, b.charge, abs_tol=1e-8)
    assert math.isclose(a.rmin, b.rmin, abs_tol=1e-7)
    assert math.isclose(a.epsilon, b.epsilon, abs_tol=1e-8)
# Compare every Amber LJ type pair against the converted Lorentz-Berthelot values.
types = {a.nb_idx: a for a in s.atoms}
gtypes = {a.type: a for a in g.atoms}
n = s.ptr("ntypes")
errors = []
for i, ai in types.items():
    for j, aj in types.items():
        k = s.parm_data["NONBONDED_PARM_INDEX"][n*(i-1)+j-1] - 1
        gi, gj = gtypes[ai.type], gtypes[aj.type]
        epsilon = math.sqrt(gi.epsilon * gj.epsilon)
        rmin = gi.rmin + gj.rmin
        for key, expected in [("LENNARD_JONES_ACOEF", epsilon*rmin**12),
                              ("LENNARD_JONES_BCOEF", 2*epsilon*rmin**6)]:
            actual = s.parm_data[key][k]
            assert math.isclose(actual, expected, rel_tol=2e-6, abs_tol=1e-8)
            if actual:
                errors.append(abs(actual-expected)/abs(actual))
top = (out / "system.top").read_text()
assert "[ settles ]" in top and "[ exclusions ]" in top
print("PASS: 1 MG + 2 Cl- + 859 WAT; 2580 atoms; charge 0.")
print("PASS: atom order, charge, mass, water/Cl LJ, box match CM.")
print("Mg IOD Rmin/2_A=1.395 epsilon_kcal/mol=0.014917")
print("Mg IOD sigma_nm=", s.atoms[0].sigma/10,
      "epsilon_kJ/mol=", s.atoms[0].epsilon*4.184)
print("Maximum absolute coordinate rounding error_A=", relative_error)
print("PASS: all 16 LJ type pairs; max relative coefficient error=", max(errors))
print("PASS: GROMACS comb-rule=2; SETTLE and exclusions present.")
print("Cross-engine periodic energy/force comparison: NOT performed.")
