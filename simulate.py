"""Small explicit-solvent ion benchmarks. See protocol.json and SOURCES.md."""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import platform as host_platform
import random
import subprocess
import time

import numpy as np
import openmm as mm
from openmm import app, unit

ROOT = Path(__file__).resolve().parent


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def lj_parameters(parameters):
    return (0.2 * parameters["rmin_half_angstrom"] / 2**(1/6),
            4.184 * parameters["epsilon_kcal_mol"])


def distances(positions, index, others, box):
    """Minimum-image distances for the orthorhombic boxes used in this study."""
    delta = positions[others] - positions[index]
    delta -= box * np.round(delta / box)
    return np.linalg.norm(delta, axis=-1)


def build_system(protocol, model, ion, seed, campaign="hydration", initial_state="hydrated"):
    random.seed(seed)  # OpenMM Modeller uses Python's RNG when choosing counterion sites.
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue(ion.upper(), chain)
    topology.addAtom(ion.upper(), app.element.get_by_symbol(ion), residue)
    modeller = app.Modeller(topology, [mm.Vec3(0, 0, 0)] * unit.nanometer)
    ff = app.ForceField("amber14/tip3p.xml")
    if campaign == "amide":
        ligand_ff = ROOT / "parameters" / "freesolv" / "acetamide.xml"
        ligand_pdb = ROOT / "parameters" / "freesolv" / "acetamide.pdb"
        if not ligand_ff.exists() or not ligand_pdb.exists():
            raise FileNotFoundError("Verified acetamide XML and PDB are required for the extension")
        ff.loadFile(str(ligand_ff))
        ligand = app.PDBFile(str(ligand_pdb))
        atoms = list(ligand.topology.atoms())
        oxygen = next(a.index for a in atoms if a.element == app.element.oxygen)
        coords = ligand.positions.value_in_unit(unit.nanometer)
        coords = np.asarray(coords) - np.asarray(coords[oxygen])
        # Point the carbonyl oxygen toward the ion; place all other ligand atoms outward.
        bonded = next(b if a.index == oxygen else a for a, b in ligand.topology.bonds()
                      if oxygen in (a.index, b.index))
        direction = coords[bonded.index] / np.linalg.norm(coords[bonded.index])
        separation = ({"Mg": 0.21, "Ca": 0.245}[ion]
                      if initial_state == "contact" else 0.65)
        coords += direction * separation
        modeller.add(ligand.topology, coords * unit.nanometer)
    modeller.addSolvent(ff, model="tip3p", boxSize=mm.Vec3(*([protocol["box_nm"]]*3)),
                       neutralize=True, ionicStrength=0*unit.molar)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=protocol["cutoff_nm"]*unit.nanometer,
                             constraints=app.HBonds, rigidWater=True,
                             ewaldErrorTolerance=protocol["ewald_tolerance"])
    atoms = list(modeller.topology.atoms())
    ion_index = next(a.index for a in atoms if a.element.symbol == ion)
    water = [a.index for a in atoms if a.residue.name == "HOH" and a.element.symbol == "O"]
    chloride = [a.index for a in atoms if a.element.symbol == "Cl"]
    amide = [a.index for a in atoms if a.element.symbol == "O" and a.residue.name != "HOH"]
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    sigma, epsilon = lj_parameters(protocol["models"][model][ion])
    charge, old_sigma, old_epsilon = nb.getParticleParameters(ion_index)
    if abs(charge.value_in_unit(unit.elementary_charge)-2) > 1e-10:
        raise ValueError("Divalent ion charge must be +2")
    if model == "CM" and (abs(old_sigma.value_in_unit(unit.nanometer)-sigma) > 1e-10
                          or abs(old_epsilon.value_in_unit(unit.kilojoule_per_mole)-epsilon) > 1e-9):
        raise ValueError("Native OpenMM CM parameters no longer match the source-verified baseline")
    nb.setParticleParameters(ion_index, charge, sigma, epsilon)
    total_charge = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                       for i in range(system.getNumParticles()))
    if abs(total_charge) > 1e-6 or len(chloride) != 2:
        raise ValueError(f"Expected neutral +2/2Cl system, found net charge {total_charge}")
    barostat = mm.MonteCarloBarostat(protocol["pressure_bar"]*unit.bar,
                                   protocol["temperature_K"]*unit.kelvin, 25)
    barostat.setRandomNumberSeed(seed + 10000)
    system.addForce(barostat)
    return modeller, system, ion_index, water, chloride, amide


def run_case(protocol, out, model, ion, replica, seed, campaign, initial_state, selected_platform):
    out.mkdir(parents=True, exist_ok=False)
    meta = {"campaign": campaign, "model": model, "ion": ion, "replica": replica,
            "initial_state": initial_state, "seed": seed, "status": "preparing",
            "temperature_K": protocol["temperature_K"], "production_ps": protocol["production_ps"],
            "sample_interval_ps": protocol["sample_interval_ps"], "protocol": protocol,
            "openmm_version": mm.version.full_version, "python": host_platform.python_version(),
            "host_os": host_platform.platform(), "process_id": os.getpid(),
            "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "protocol_sha256": hashlib.sha256((ROOT/"protocol.json").read_bytes()).hexdigest(),
            "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "git_dirty_at_start": bool(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True).strip()),
            "started_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    write_json(out/"metadata.json", meta)
    start = time.perf_counter()
    modeller, system, ion_index, water, chloride, amide = build_system(
        protocol, model, ion, seed, campaign, initial_state)
    integrator = mm.LangevinMiddleIntegrator(protocol["temperature_K"]*unit.kelvin,
                                            protocol["friction_per_ps"]/unit.picosecond,
                                            protocol["timestep_ps"]*unit.picosecond)
    integrator.setRandomNumberSeed(seed)
    plat = mm.Platform.getPlatformByName(selected_platform)
    properties = {"Precision": "mixed"} if selected_platform in ("CUDA", "OpenCL") else {}
    if selected_platform == "CPU":
        properties["Threads"] = "4"
    simulation = app.Simulation(modeller.topology, system, integrator, plat, properties)
    meta.update(platform=selected_platform, platform_properties={
        name: plat.getPropertyValue(simulation.context, name) for name in plat.getPropertyNames()},
        n_atoms=system.getNumParticles(), n_water=len(water), n_chloride=len(chloride),
        water_oxygen_indices=water, box_geometry="cubic with isotropic barostat", min_box_length_nm=None,
        ion_lj_sigma_nm=lj_parameters(protocol["models"][model][ion])[0],
        ion_lj_epsilon_kj_mol=lj_parameters(protocol["models"][model][ion])[1])
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=2000)
    simulation.context.setVelocitiesToTemperature(protocol["temperature_K"]*unit.kelvin, seed+20000)
    with (out/"initial.pdb").open("w") as stream:
        app.PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True).getPositions(), stream)
    (out/"system.xml").write_text(mm.XmlSerializer.serialize(system), encoding="utf-8")
    meta["status"] = "equilibrating"
    write_json(out/"metadata.json", meta)
    eq_steps = round(protocol["equilibration_ps"]/protocol["timestep_ps"])
    simulation.step(eq_steps)
    sample_steps = round(protocol["sample_interval_ps"]/protocol["timestep_ps"])
    n_frames = round(protocol["production_ps"]/protocol["sample_interval_ps"])
    simulation.reporters.append(app.DCDReporter(str(out/"trajectory.dcd"), sample_steps*10))
    simulation.reporters.append(app.StateDataReporter(str(out/"thermodynamics.csv"), sample_steps*10,
                                step=True, time=True, potentialEnergy=True, temperature=True,
                                density=True, volume=True, speed=True))
    arrays = {"distances_nm": np.empty((n_frames, len(water)), dtype=np.float32)}
    for field in ("time_ps", "volume_nm3", "temperature_K", "potential_kj_mol", "min_chloride_nm", "amide_distance_nm"):
        arrays[field] = np.empty(n_frames, dtype=np.float64)
    dof = 3*system.getNumParticles()-system.getNumConstraints()-3
    gas_constant = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin)
    meta.update(status="running", equilibration_wall_seconds=time.perf_counter()-start,
                degrees_of_freedom=dof)
    write_json(out/"metadata.json", meta)
    production_start = time.perf_counter()
    for frame in range(n_frames):
        simulation.step(sample_steps)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        box = np.diag(vectors)
        if not np.allclose(vectors, np.diag(box)):
            raise ValueError("Analysis currently requires an orthorhombic box")
        if min(box) <= 2*protocol["cutoff_nm"]:
            raise ValueError("Box became too short for the configured cutoff")
        meta["min_box_length_nm"] = min(float(min(box)), meta["min_box_length_nm"] or float(min(box)))
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(xyz).all() or not np.isfinite(energy):
            raise ArithmeticError("Nonfinite simulation state")
        arrays["distances_nm"][frame] = distances(xyz, ion_index, water, box)
        arrays["time_ps"][frame] = (frame+1)*protocol["sample_interval_ps"]
        arrays["volume_nm3"][frame] = np.prod(box)
        arrays["temperature_K"][frame] = 2*state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(dof*gas_constant)
        arrays["potential_kj_mol"][frame] = energy
        arrays["min_chloride_nm"][frame] = min(distances(xyz, ion_index, chloride, box))
        arrays["amide_distance_nm"][frame] = distances(xyz, ion_index, amide, box)[0] if amide else np.nan
        if (frame+1) % max(1, n_frames//10) == 0 or frame == n_frames-1:
            meta.update(completed_frames=frame+1, last_sample_ps=arrays["time_ps"][frame],
                        wall_seconds=time.perf_counter()-start)
            write_json(out/"metadata.json", meta)
            simulation.saveCheckpoint(str(out/"checkpoint.chk"))
            print(f"{out.name}: {frame+1}/{n_frames} frames", flush=True)
        if time.perf_counter()-start > 3600:
            raise TimeoutError("One-hour per-case safety ceiling exceeded; inspect the saved checkpoint")
    wall_production = time.perf_counter()-production_start
    np.savez_compressed(out/"samples.npz", **arrays)
    final = simulation.context.getState(getPositions=True)
    with (out/"final.pdb").open("w") as stream:
        app.PDBFile.writeFile(modeller.topology, final.getPositions(), stream)
    meta.update(status="complete", production_wall_seconds=wall_production,
                ns_per_day=protocol["production_ps"]/1000/wall_production*86400,
                wall_seconds=time.perf_counter()-start,
                completed_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                samples_sha256=hashlib.sha256((out/"samples.npz").read_bytes()).hexdigest())
    write_json(out/"metadata.json", meta)
    print(f"Completed {out.name}: {meta['ns_per_day']:.1f} ns/day", flush=True)
    del simulation, integrator, system
    gc.collect()


def selfcheck():
    assert np.allclose(distances(np.array([[0.,0,0],[2.9,0,0],[0,0.2,0]]), 0, [1,2], np.array([3.,3,3])), [.1,.2])
    sigma, epsilon = lj_parameters({"rmin_half_angstrom": 1.36, "epsilon_kcal_mol": .01020237})
    assert abs(sigma-.2423244513341723) < 1e-12
    assert abs(epsilon-.04268671608) < 1e-12
    # At the LJ minimum, the energy equals -epsilon; verify the Amber convention.
    r = 2*.136
    assert abs(4*epsilon*((sigma/r)**12-(sigma/r)**6)+epsilon) < 1e-12
    print("Simulation self-check passed: periodic distance and source-unit/LJ conversion.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", choices=["pilot", "hydration", "amide"], default="pilot")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--platform", choices=["OpenCL", "CUDA", "CPU", "Reference"], default="OpenCL")
    parser.add_argument("--selfcheck", action="store_true")
    args = parser.parse_args()
    selfcheck()
    if args.selfcheck:
        return
    protocol = json.loads((ROOT/"protocol.json").read_text(encoding="utf-8"))
    out = args.output or ROOT/"results"/args.campaign
    if args.campaign == "pilot":
        protocol.update(equilibration_ps=10.0, production_ps=50.0)
        cases = [("CM", "Mg", 1, protocol["seeds"][0], "hydrated")]
    else:
        states = ["contact", "separated"] if args.campaign == "amide" else ["hydrated"]
        cases = [(model, ion, rep, seed+(50000 if state == "contact" else 0), state)
                 for model in protocol["models"] for ion in ("Mg", "Ca")
                 for rep, seed in zip(protocol["replicas"], protocol["seeds"]) for state in states]
    for model, ion, replica, seed, initial in cases:
        destination = out/f"{model}_{ion}_r{replica}_{initial}"
        if destination.exists():
            meta_file = destination/"metadata.json"
            meta = json.loads(meta_file.read_text()) if meta_file.exists() else {}
            samples = destination/"samples.npz"
            if (meta.get("status") == "complete" and meta.get("protocol") == protocol
                    and samples.exists() and hashlib.sha256(samples.read_bytes()).hexdigest() == meta.get("samples_sha256")):
                print(f"Verified completed case, skipping {destination.name}", flush=True)
                continue
            raise FileExistsError(f"Incomplete or changed case {destination}; inspect before choosing a new output")
        run_case(protocol, destination, model, ion, replica, seed, args.campaign, initial, args.platform)


if __name__ == "__main__":
    main()
