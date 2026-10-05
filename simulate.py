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
import tempfile
import time

import numpy as np
import openmm as mm
from openmm import app, unit

ROOT = Path(__file__).resolve().parent


def write_bytes(path, content):
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("wb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_json(path, value):
    write_bytes(path, (json.dumps(value, indent=2, allow_nan=False) + "\n").encode("utf-8"))


def freeze_inputs(source_path, protocol_path, repository):
    """Read source, protocol, and Git identity once for the whole invocation."""
    source = source_path.read_bytes()
    protocol = protocol_path.read_bytes()
    metadata = {
        "code_sha256": hashlib.sha256(source).hexdigest(),
        "protocol_sha256": hashlib.sha256(protocol).hexdigest(),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip(),
        "git_dirty_at_start": bool(subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"], cwd=repository, text=True).strip()),
        "invocation_started_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    return json.loads(protocol), metadata


def write_samples(path, arrays, completed_frames, metadata=None):
    """Atomically replace an NPZ containing only fully initialized frames."""
    if not arrays or not 0 <= completed_frames <= min(len(value) for value in arrays.values()):
        raise ValueError("Invalid completed frame count")
    snapshot = {name: value[:completed_frames] for name, value in arrays.items()}
    if metadata is not None:
        snapshot["metadata_json"] = np.array(json.dumps(metadata, allow_nan=False))
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("wb") as stream:
            np.savez_compressed(stream, **snapshot)
            stream.flush()
            os.fsync(stream.fileno())
        checksum = hashlib.sha256(temporary.read_bytes()).hexdigest()
        os.replace(temporary, path)
        return checksum
    finally:
        temporary.unlink(missing_ok=True)


def save_partial(out, arrays, completed_frames, meta, start):
    meta.update(completed_frames=completed_frames,
                last_sample_ps=float(arrays["time_ps"][completed_frames-1]) if completed_frames else None,
                wall_seconds=time.perf_counter()-start)
    # Embed metadata so the NPZ is interpretable even if interruption falls
    # between replacing it and replacing its external hash/metadata record.
    snapshot_meta = {key: value for key, value in meta.items() if not key.startswith("partial_samples_")}
    meta["partial_samples_sha256"] = write_samples(
        out/"samples.partial.npz", arrays, completed_frames, snapshot_meta)
    meta["partial_samples_frames"] = completed_frames
    write_json(out/"metadata.json", meta)


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


def run_case(protocol, out, model, ion, replica, seed, campaign, initial_state, selected_platform, provenance):
    out.mkdir(parents=True, exist_ok=False)
    meta = {"campaign": campaign, "model": model, "ion": ion, "replica": replica,
            "initial_state": initial_state, "seed": seed, "status": "preparing",
            "temperature_K": protocol["temperature_K"], "production_ps": protocol["production_ps"],
            "sample_interval_ps": protocol["sample_interval_ps"], "protocol": protocol,
            "openmm_version": mm.version.full_version, "python": host_platform.python_version(),
            "host_os": host_platform.platform(), "process_id": os.getpid(), **provenance,
            "started_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    write_json(out/"metadata.json", meta)
    start = time.perf_counter()
    arrays = None
    completed_frames = 0
    simulation = integrator = system = None
    try:
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
            completed_frames = frame+1
            if completed_frames % max(1, n_frames//10) == 0 or frame == n_frames-1:
                save_partial(out, arrays, completed_frames, meta, start)
                write_bytes(out/"checkpoint.chk", simulation.context.createCheckpoint())
                print(f"{out.name}: {frame+1}/{n_frames} frames", flush=True)
            if time.perf_counter()-start > 3600:
                raise TimeoutError("One-hour per-case safety ceiling exceeded; inspect the saved checkpoint")
        wall_production = time.perf_counter()-production_start
        samples_sha256 = write_samples(out/"samples.npz", arrays, completed_frames)
        final = simulation.context.getState(getPositions=True)
        modeller.topology.setPeriodicBoxVectors(final.getPeriodicBoxVectors())
        with (out/"final.pdb").open("w") as stream:
            app.PDBFile.writeFile(modeller.topology, final.getPositions(), stream)
        meta.update(status="complete", production_wall_seconds=wall_production,
                    ns_per_day=protocol["production_ps"]/1000/wall_production*86400,
                    wall_seconds=time.perf_counter()-start,
                    completed_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    samples_sha256=samples_sha256)
        write_json(out/"metadata.json", meta)
        print(f"Completed {out.name}: {meta['ns_per_day']:.1f} ns/day", flush=True)
    except BaseException as error:
        meta.update(status="failed", error_type=type(error).__name__, error_message=str(error),
                    completed_frames=completed_frames, wall_seconds=time.perf_counter()-start,
                    failed_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        try:
            if arrays is not None:
                save_partial(out, arrays, completed_frames, meta, start)
            else:
                write_json(out/"metadata.json", meta)
        except Exception as save_error:
            error.add_note(f"Could not update failure/partial-sample record: {save_error}")
            meta["partial_save_error"] = str(save_error)
            try:
                write_json(out/"metadata.json", meta)
            except Exception as metadata_error:
                error.add_note(f"Could not mark the case failed: {metadata_error}")
        raise
    finally:
        del simulation, integrator, system
        gc.collect()



def check_ion_water_pairs():
    """Check native TIP3P oxygen mixing and forces against the analytic pair law."""
    topology = app.Topology()
    residue = topology.addResidue("HOH", topology.addChain())
    oxygen = topology.addAtom("O", app.element.oxygen, residue)
    for name in ("H1", "H2"):
        topology.addBond(oxygen, topology.addAtom(name, app.element.hydrogen, residue))
    water_system = app.ForceField("amber14/tip3p.xml").createSystem(topology)
    water_force = next(f for f in water_system.getForces() if isinstance(f, mm.NonbondedForce))
    oxygen_parameters = water_force.getParticleParameters(0)
    oxygen_q = oxygen_parameters[0].value_in_unit(unit.elementary_charge)
    oxygen_sigma = oxygen_parameters[1].value_in_unit(unit.nanometer)
    oxygen_epsilon = oxygen_parameters[2].value_in_unit(unit.kilojoule_per_mole)
    # Published Rmin/2 (angstrom) and epsilon (kcal/mol), independent of protocol.json.
    cases = (("CM", "Mg", 1.360, .01020237), ("CM", "Ca", 1.649, .10592870),
             ("IOD", "Mg", 1.395, .01491700), ("IOD", "Ca", 1.608, .08337961))
    for model, ion, radius, epsilon in cases:
        system = mm.System()
        system.addParticle(24.305 if ion == "Mg" else 40.08)
        system.addParticle(15.99943)
        force = mm.NonbondedForce()
        force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        force.addParticle(2, *lj_parameters({"rmin_half_angstrom": radius, "epsilon_kcal_mol": epsilon}))
        force.addParticle(*oxygen_parameters)
        system.addForce(force)
        integrator = mm.VerletIntegrator(.001)
        context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
        mixed_rmin = .1*radius + oxygen_sigma*2**(1/6)/2
        mixed_epsilon = (4.184*epsilon*oxygen_epsilon)**.5
        coulomb = 138.935457644382 * 2 * oxygen_q  # kJ mol^-1 nm, OpenMM units.
        for r in (.21, .30, .45):
            context.setPositions([[0, 0, 0], [r, 0, 0]])
            state = context.getState(getEnergy=True, getForces=True)
            sixth = (mixed_rmin/r)**6
            expected_energy = coulomb/r + mixed_epsilon*(sixth**2-2*sixth)
            oxygen_fx = coulomb/r**2 + 12*mixed_epsilon*(sixth**2-sixth)/r
            assert np.isclose(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
                              expected_energy, rtol=1e-10, atol=1e-7), (model, ion, r, "energy")
            assert np.allclose(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer),
                               [[-oxygen_fx, 0, 0], [oxygen_fx, 0, 0]], rtol=1e-10, atol=1e-6), (model, ion, r, "force")
        del context, integrator


def selfcheck():
    assert np.allclose(distances(np.array([[0.,0,0],[2.9,0,0],[0,0.2,0]]), 0, [1,2], np.array([3.,3,3])), [.1,.2])
    sigma, epsilon = lj_parameters({"rmin_half_angstrom": 1.36, "epsilon_kcal_mol": .01020237})
    assert abs(sigma-.2423244513341723) < 1e-12
    assert abs(epsilon-.04268671608) < 1e-12
    # At the LJ minimum, the energy equals -epsilon; verify the Amber convention.
    r = 2*.136
    assert abs(4*epsilon*((sigma/r)**12-(sigma/r)**6)+epsilon) < 1e-12
    with tempfile.TemporaryDirectory() as directory:
        temporary = Path(directory)
        arrays = {"time_ps": np.arange(1., 6.), "distances_nm": np.arange(10.).reshape(5, 2)}
        arrays["distances_nm"][2:] = np.nan  # Represents the uninitialized tail.
        meta = {"status": "running"}
        save_partial(temporary, arrays, 2, meta, time.perf_counter())
        path = temporary/"samples.partial.npz"
        before = path.read_bytes()
        assert hashlib.sha256(before).hexdigest() == meta["partial_samples_sha256"]
        with np.load(path, allow_pickle=False) as saved:
            assert saved["distances_nm"].shape == (2, 2) and np.isfinite(saved["distances_nm"]).all()
            assert json.loads(str(saved["metadata_json"]))["completed_frames"] == 2
        original_save = np.savez_compressed
        def interrupted_save(stream, **values):
            stream.write(b"incomplete ZIP")
            raise OSError("injected write failure")
        np.savez_compressed = interrupted_save
        try:
            try:
                write_samples(path, arrays, 2)
            except OSError:
                pass
            else:
                raise AssertionError("Injected write failure was swallowed")
        finally:
            np.savez_compressed = original_save
        assert path.read_bytes() == before and not path.with_name(path.name+".tmp").exists()
        code_path, protocol_path = temporary/"source.py", temporary/"protocol.json"
        code_path.write_bytes(b"original source\n")
        protocol_path.write_bytes(b'{"production_ps": 5}\n')
        original_git = subprocess.check_output
        subprocess.check_output = lambda command, **kwargs: "0"*40 if "rev-parse" in command else ""
        try:
            protocol, provenance = freeze_inputs(code_path, protocol_path, temporary)
        finally:
            subprocess.check_output = original_git
        code_path.write_bytes(b"changed source\n")
        protocol_path.write_bytes(b'{"production_ps": 99}\n')
        assert protocol["production_ps"] == 5
        assert provenance["code_sha256"] == hashlib.sha256(b"original source\n").hexdigest()
        assert provenance["protocol_sha256"] == hashlib.sha256(b'{"production_ps": 5}\n').hexdigest()
        original_build = build_system
        def failed_build(*args, **kwargs):
            raise RuntimeError("injected preparation failure")
        globals()["build_system"] = failed_build
        try:
            try:
                run_case(dict(protocol, temperature_K=298.15, sample_interval_ps=1.), temporary/"failed-case",
                         "CM", "Mg", 1, 1, "pilot", "hydrated", "Reference", provenance)
            except RuntimeError:
                pass
            else:
                raise AssertionError("Injected case failure was swallowed")
        finally:
            globals()["build_system"] = original_build
        failed = json.loads((temporary/"failed-case/metadata.json").read_text())
        assert failed["status"] == "failed" and failed["completed_frames"] == 0
    check_ion_water_pairs()
    print("Simulation self-check passed: distances/LJ units, 12 ion-water energy/force checks, atomic initialized samples, and frozen input provenance.")


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
    protocol, provenance = freeze_inputs(Path(__file__), ROOT/"protocol.json", ROOT)
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
        run_case(protocol, destination, model, ion, replica, seed, args.campaign, initial, args.platform, provenance)


if __name__ == "__main__":
    main()
