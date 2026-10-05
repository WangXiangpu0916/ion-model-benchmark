"""Convert the pinned FreeSolv acetamide inputs and verify the generated OpenMM XML.

Run with the project's .venv Python. No fitting or geometry minimization is done.
"""

from collections import defaultdict
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET

import numpy as np
import openmm as mm
from openmm import app, unit


HERE = Path(__file__).resolve().parent
ENERGY = unit.kilojoule_per_mole
LENGTH = unit.nanometer


def scalar(value):
    return value.value_in_unit_system(unit.md_unit_system) if unit.is_quantity(value) else value


def node(parent, tag, **attrs):
    return ET.SubElement(parent, tag, {key: str(scalar(value)) for key, value in attrs.items()})


def convert():
    source = app.AmberPrmtopFile(str(HERE / "mobley_8048190.prmtop"))
    coordinates = app.AmberInpcrdFile(str(HERE / "mobley_8048190.inpcrd"))
    system = source.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
    atoms = list(source.topology.atoms())
    assert len(atoms) == 9 and source.topology.getNumResidues() == 1
    atoms[0].residue.name = "ACM"
    atom_types = ["freesolv-acetamide-" + atom.name for atom in atoms]
    forces = {type(force).__name__: force for force in system.getForces()}
    assert set(forces) == {"HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "NonbondedForce", "CMMotionRemover"}
    nonbonded = forces["NonbondedForce"]
    root = ET.Element("ForceField")
    root.append(ET.Comment("Generated from pinned FreeSolv mobley_8048190.prmtop; see SOURCES.md and source_manifest.json."))
    types = node(root, "AtomTypes")
    residue = node(node(root, "Residues"), "Residue", name="ACM")
    for atom, atom_type in zip(atoms, atom_types):
        node(types, "Type", name=atom_type, element=atom.element.symbol, mass=system.getParticleMass(atom.index), **{"class": atom_type})
        node(residue, "Atom", name=atom.name, type=atom_type)
    edges = {frozenset((a.index, b.index)) for a, b in source.topology.bonds()}
    for a, b in source.topology.bonds():
        node(residue, "Bond", atomName1=a.name, atomName2=b.name)
    bonds = node(root, "HarmonicBondForce")
    source_bonds = forces["HarmonicBondForce"]
    for index in range(source_bonds.getNumBonds()):
        a, b, length, k = source_bonds.getBondParameters(index)
        node(bonds, "Bond", type1=atom_types[a], type2=atom_types[b], length=length, k=k)
    angles = node(root, "HarmonicAngleForce")
    source_angles = forces["HarmonicAngleForce"]
    for index in range(source_angles.getNumAngles()):
        a, b, c, angle, k = source_angles.getAngleParameters(index)
        node(angles, "Angle", type1=atom_types[a], type2=atom_types[b], type3=atom_types[c], angle=angle, k=k)
    torsions = node(root, "PeriodicTorsionForce", ordering="amber")
    grouped = defaultdict(list)
    source_torsions = forces["PeriodicTorsionForce"]
    for index in range(source_torsions.getNumTorsions()):
        a, b, c, d, periodicity, phase, k = source_torsions.getTorsionParameters(index)
        if all(frozenset(edge) in edges for edge in ((a, b), (b, c), (c, d))):
            tag, ordered = "Proper", (a, b, c, d)
        else:
            # Amber's improper central atom is third; XML lists it first.
            assert all(frozenset((c, other)) in edges for other in (a, b, d))
            tag, ordered = "Improper", (c, a, b, d)
        grouped[tag, ordered].append((periodicity, phase, k))
    for (tag, indices), terms in grouped.items():
        attrs = {f"type{i + 1}": atom_types[index] for i, index in enumerate(indices)}
        for i, (periodicity, phase, k) in enumerate(terms, 1):
            attrs.update({f"periodicity{i}": periodicity, f"phase{i}": phase, f"k{i}": k})
        node(torsions, tag, **attrs)
    nb = node(root, "NonbondedForce", coulomb14scale=1 / 1.2, lj14scale=0.5)
    for i, atom_type in enumerate(atom_types):
        q, sigma, epsilon = nonbonded.getParticleParameters(i)
        node(nb, "Atom", type=atom_type, charge=q, sigma=sigma, epsilon=epsilon)
    ET.indent(root)
    ET.ElementTree(root).write(HERE / "acetamide.xml", encoding="utf-8", xml_declaration=True)
    with (HERE / "acetamide.pdb").open("w", encoding="ascii", newline="\n") as handle:
        handle.write("REMARK   1 FREESOLV MOBLEY_8048190 ACETAMIDE; ORIGINAL COORDINATES\n")
        app.PDBFile.writeModel(source.topology, coordinates.positions, handle)
        app.PDBFile.writeFooter(source.topology, handle)
    return system, coordinates.positions


def verify(original, positions):
    pdb = app.PDBFile(str(HERE / "acetamide.pdb"))
    ff = app.ForceField("amber14/tip3p.xml", str(HERE / "acetamide.xml"))
    converted = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    assert original.getNumParticles() == converted.getNumParticles() == 9
    assert pdb.topology.getNumBonds() == 8
    groups = {"HarmonicBondForce": 0, "HarmonicAngleForce": 1, "PeriodicTorsionForce": 2, "NonbondedForce": 3, "CMMotionRemover": 4}
    counts = []
    exception_tables = []
    for system in (original, converted):
        count = {}
        for force in system.getForces():
            name = type(force).__name__
            force.setForceGroup(groups[name])
            if name == "HarmonicBondForce":
                count["bonds"] = force.getNumBonds()
            elif name == "HarmonicAngleForce":
                count["angles"] = force.getNumAngles()
            elif name == "PeriodicTorsionForce":
                count["nonzero_torsion_terms"] = sum(scalar(force.getTorsionParameters(i)[6]) != 0 for i in range(force.getNumTorsions()))
            elif name == "NonbondedForce":
                table = {}
                for i in range(force.getNumExceptions()):
                    a, b, q, sigma, epsilon = force.getExceptionParameters(i)
                    table[tuple(sorted((a, b)))] = np.array([scalar(q), scalar(sigma) if scalar(epsilon) else 0.0, scalar(epsilon)])
                exception_tables.append(table)
                count["nonbonded_exceptions"] = len(table)
        counts.append(count)
    assert counts[0] == counts[1] == {"bonds": 8, "angles": 12, "nonzero_torsion_terms": 14, "nonbonded_exceptions": 30}
    assert exception_tables[0].keys() == exception_tables[1].keys()
    exception_error = np.max([np.abs(exception_tables[0][key] - exception_tables[1][key]) for key in exception_tables[0]], axis=0)
    # Source LJ pair coefficients are stored to 8 significant digits; XML
    # recombines atomic LJ values. Verify that this roundoff stays below 1e-8.
    assert np.all(exception_error < 1e-8), exception_error
    contexts = []
    integrators = []
    for system in (original, converted):
        integrator = mm.VerletIntegrator(0.001)
        integrators.append(integrator)
        contexts.append(mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference")))
    rng = np.random.default_rng(20261005)
    base = np.array(positions.value_in_unit(LENGTH))
    assert np.max(np.abs(np.array(pdb.positions.value_in_unit(LENGTH)) - base)) < 1e-12
    samples = []
    for sample, perturbation_nm in enumerate((0.0, 0.002, 0.005, 0.01, 0.02)):
        xyz = base + rng.normal(0, perturbation_nm, base.shape)
        for context in contexts:
            context.setPositions(xyz)
        by_force = {}
        for name, group in groups.items():
            if name == "CMMotionRemover":
                continue
            states = [c.getState(getEnergy=True, getForces=True, groups=1 << group) for c in contexts]
            energy_error = abs(scalar(states[0].getPotentialEnergy() - states[1].getPotentialEnergy()))
            force_error = float(np.max(np.abs(states[0].getForces(asNumpy=True).value_in_unit(ENERGY / LENGTH) - states[1].getForces(asNumpy=True).value_in_unit(ENERGY / LENGTH))))
            assert energy_error < 1e-6, (sample, name, energy_error)
            assert force_error < 1e-5, (sample, name, force_error)
            by_force[name] = {"absolute_energy_error_kj_mol": energy_error, "max_absolute_force_error_kj_mol_nm": force_error}
        samples.append({"perturbation_std_nm": perturbation_nm, "by_force": by_force})
    result = {"passed": True, "openmm_version": mm.version.version, "platform": "Reference", "random_seed": 20261005, "atom_count": 9, "force_counts": counts[0], "improper_terms": 2, "source_zero_k_torsions_omitted_by_openmm_xml_loader": 3, "max_exception_errors_chargeprod_sigma_epsilon": exception_error.tolist(), "tolerances": {"absolute_energy_kj_mol": 1e-6, "max_absolute_force_kj_mol_nm": 1e-5}, "geometries": samples, "sha256": {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in ("mobley_8048190.prmtop", "mobley_8048190.inpcrd", "acetamide.xml", "acetamide.pdb")}}
    (HERE / "conversion_verification.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": True, "geometries_checked": len(samples), "force_counts": counts[0], "max_energy_error_kj_mol": max(v["absolute_energy_error_kj_mol"] for s in samples for v in s["by_force"].values()), "max_force_error_kj_mol_nm": max(v["max_absolute_force_error_kj_mol_nm"] for s in samples for v in s["by_force"].values())}, indent=2))


if __name__ == "__main__":
    verify(*convert())
