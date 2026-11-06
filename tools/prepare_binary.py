#!/usr/bin/env python3
"""Prepare a homogeneous organic-solvent/water control, without PLGA or drug.

Run structures in the existing RDKit preparation environment, parameters in the
Amber environment, and assemble in the preparation environment. The plan command
uses only the Python standard library. All geometry is an initial condition, not
an assumption that the force field reproduces experimental miscibility.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

SCHEMA_VERSION = 1
AVOGADRO_NM3 = 602.214076  # molecules nm^-3 per (g cm^-3 / g mol^-1)
# Room-temperature pure-liquid estimates, used ONLY to choose molecule counts.
# NPT relaxes the non-additive mixture volume; no density is imposed during MD.
SOLVENTS = {
    "DMSO": {"smiles": "CS(=O)C", "resname": "DMS", "density": 1.095,
             "mass_g_mol": 78.134, "atoms_per_molecule": 10},
    "DMF": {"smiles": "CN(C)C=O", "resname": "DMF", "density": 0.944,
            "mass_g_mol": 73.095, "atoms_per_molecule": 12},
    "Acetone": {"smiles": "CC(=O)C", "resname": "ACE", "density": 0.785,
                "mass_g_mol": 58.080, "atoms_per_molecule": 10},
    "Acetonitrile": {"smiles": "CC#N", "resname": "ACN", "density": 0.786,
                     "mass_g_mol": 41.053, "atoms_per_molecule": 6},
}
WATER = {"smiles": "O", "resname": "WAT", "density": 0.997,
         "mass_g_mol": 18.01528, "atoms_per_molecule": 3}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def canonical_solvent(value):
    aliases = {name.lower(): name for name in SOLVENTS}
    aliases["acn"] = "Acetonitrile"
    try:
        return aliases[value.strip().lower()]
    except KeyError:
        raise argparse.ArgumentTypeError("Choose DMSO, DMF, Acetone, or Acetonitrile (ACN)") from None


def positive_number(value):
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("Expected a finite positive number")
    return result


def steps_for(duration_ps, dt_ps):
    if not math.isfinite(duration_ps) or not math.isfinite(duration_ps / dt_ps):
        raise ValueError("Stage duration is too large")
    count = round(duration_ps / dt_ps)
    if count < 1 or not math.isclose(count * dt_ps, duration_ps, rel_tol=1e-10, abs_tol=1e-9):
        raise ValueError(f"Stage duration {duration_ps} ps must be a positive multiple of {dt_ps} ps")
    return count


def validate_config(cfg):
    if cfg.get("schema_version") != SCHEMA_VERSION or cfg.get("solvent") not in SOLVENTS:
        raise ValueError("Unsupported binary configuration schema or solvent")
    box = cfg["box_nm"]
    if len(box) != 3 or any(not math.isfinite(v) or v < 3 for v in box) or len(set(box)) != 1:
        raise ValueError("Use a cubic box with an edge of at least 3 nm")
    for value in (cfg["temperature_K"], cfg["organic_volume"], cfg["water_volume"],
                  cfg["run"]["nvt_ps"], cfg["run"]["npt_ps"], cfg["run"]["production_ns"]):
        if not math.isfinite(value) or value <= 0:
            raise ValueError("Temperature, volumes and durations must be finite and positive")
    if type(cfg["seed"]) is not int or not 1 <= cfg["seed"] <= 2147483646:
        raise ValueError("Seed must be an integer from 1 to 2147483646")
    for stage, duration, dt in (("nvt", cfg["run"]["nvt_ps"], .001),
                                ("npt", cfg["run"]["npt_ps"], .002),
                                ("production", cfg["run"]["production_ns"] * 1000, .002)):
        try:
            steps_for(duration, dt)
        except ValueError as error:
            raise ValueError(f"{stage}: {error}") from error
    return cfg


def make_plan(cfg):
    validate_config(cfg)
    volume_nm3 = math.prod(cfg["box_nm"])
    if not math.isfinite(volume_nm3) or volume_nm3 > 1e6:
        raise ValueError("Box is too large for this at-most-99,999-atom packing workflow")
    volume_scale = max(cfg["organic_volume"], cfg["water_volume"])
    vo, vw = cfg["organic_volume"] / volume_scale, cfg["water_volume"] / volume_scale
    organic_fraction = vo / (vo + vw)
    organic = SOLVENTS[cfg["solvent"]]
    nw = round(AVOGADRO_NM3 * WATER["density"] * volume_nm3 * (1 - organic_fraction) / WATER["mass_g_mol"])
    no = round(AVOGADRO_NM3 * organic["density"] * volume_nm3 * organic_fraction / organic["mass_g_mol"])
    if min(nw, no) < 1:
        raise ValueError("This box and composition must contain at least one molecule of each component")
    count = nw * 3 + no * organic["atoms_per_molecule"]
    if count > 99999:
        raise ValueError("This packing workflow supports at most 99,999 atoms; reduce box size")
    actual_ratio = (no * organic["mass_g_mol"] / organic["density"]) / (nw * WATER["mass_g_mol"] / WATER["density"])
    return {"solvent": cfg["solvent"], "box_nm": cfg["box_nm"], "temperature_K": cfg["temperature_K"],
            "water_molecules": nw, "organic_molecules": no, "atom_count": count,
            "organic_mole_fraction": no / (nw + no),
            "requested_organic_volume_fraction": organic_fraction,
            "realized_organic_to_water_volume_ratio_from_counts": actual_ratio,
            "initial_density_g_cm3": (nw * WATER["mass_g_mol"] + no * organic["mass_g_mol"]) / (AVOGADRO_NM3 * volume_nm3),
            "count_basis": "pure-liquid additive volumes; integer-rounded counts; NPT volume is free",
            "force_field": "GAFF2 / AM1-BCC organic solvent + TIP3P water; provisional mixture model",
            "scope": "salt-free binary control; no PLGA, drug, poloxamer, injection, or evaporation"}


def read_prepared(out):
    cfg = validate_config(json.loads((out / "binary_config.json").read_text()))
    manifest = json.loads((out / "manifest.json").read_text())
    if manifest.get("configuration_sha256") != digest(cfg):
        raise ValueError("Prepared configuration changed; create a new run directory")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported manifest schema")
    return cfg, manifest


def structures(cfg, out):
    from rdkit import Chem
    from prepare_case import _small_molecule, _write_molecule

    plan = make_plan(cfg)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Refusing to replace nonempty preparation directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "binary_config.json", cfg)
    species, offset = [], 0
    for species_id, kind, data, count in (("WATER", "water", WATER, plan["water_molecules"]),
                                          ("ORG", "solvent", SOLVENTS[cfg["solvent"]], plan["organic_molecules"])):
        is_water = kind == "water"
        mol = _small_molecule(data["smiles"], cfg["seed"], water=is_water)
        if mol.GetNumAtoms() != data["atoms_per_molecule"]:
            raise RuntimeError(f"Unexpected atom count for {species_id}")
        folder = out / "molecules" / species_id
        _write_molecule(folder, mol, data["resname"], water=is_water)
        names = [line[12:16].strip() for line in (folder / "input.pdb").read_text().splitlines()
                 if line.startswith(("ATOM  ", "HETATM"))]
        stop = offset + count * mol.GetNumAtoms()
        species.append({"id": species_id, "kind": kind, "resname": data["resname"],
                        "count": count, "atoms_per_molecule": mol.GetNumAtoms(),
                        "atom_names": names, "atomic_numbers": [a.GetAtomicNum() for a in mol.GetAtoms()],
                        "formal_charge": Chem.GetFormalCharge(mol), "mass_g_mol": data["mass_g_mol"],
                        "initialization_density_g_cm3": data["density"],
                        "atom_start_0based": offset, "atom_stop_0based": stop,
                        "representative_atom_indices_0based": list(range(offset, stop, mol.GetNumAtoms()))})
        offset = stop
    manifest = {"schema_version": SCHEMA_VERSION, "configuration_sha256": digest(cfg),
                "solvent": cfg["solvent"], "temperature_K": cfg["temperature_K"],
                "box_nm": cfg["box_nm"], "atom_count": offset,
                "molecule_count": sum(s["count"] for s in species), "species": species,
                "initialization": "homogeneous random packing throughout the periodic cube",
                "plan": plan,
                "experimental_reference": {"doi": "10.1002/mame.202500269", "section": "4.3",
                                           "organic_ml": 1.5, "aqueous_ml": 10.0,
                                           "aqueous_poloxamer_percent": 0.5,
                                           "binary_control_omits_poloxamer": True},
                "source_sha256": {name: hashlib.sha256((Path(__file__).parent / name).read_bytes()).hexdigest()
                                  for name in ("prepare_binary.py", "prepare_case.py", "mdprep_forcefield.py")}}
    write_json(out / "manifest.json", manifest)
    print(json.dumps(plan, indent=2), flush=True)


def parameters(out):
    from mdprep_forcefield import _antechamber, _topology

    _, manifest = read_prepared(out)
    (out / "READY.json").unlink(missing_ok=True)
    for command in ("antechamber", "parmchk2", "tleap"):
        if shutil.which(command) is None:
            raise RuntimeError(f"{command} is missing; load the existing Amber module environment")
    for species in manifest["species"]:
        folder = out / "molecules" / species["id"]
        if species["formal_charge"] != 0:
            raise ValueError("Expected neutral water and solvent molecules")
        if species["kind"] == "water":
            _topology(folder, water=True)
        else:
            _antechamber(folder)
            _topology(folder)
    write_json(out / "parameterization.json", {
        "configuration_sha256": manifest["configuration_sha256"],
        "organic_method": "GAFF2 / AM1-BCC; existing small-molecule helpers",
        "water_method": "Amber TIP3P",
        "status": "parameters generated; mixture thermodynamics not validated",
        "species": [s["id"] for s in manifest["species"]]})


def pack(cfg, out, manifest):
    import numpy as np

    if shutil.which("packmol") is None:
        raise RuntimeError("Packmol is missing from the existing preparation environment")
    box = np.array(cfg["box_nm"]) * 10
    limits = " ".join(f"{v:.8f}" for v in box)
    lines = ["tolerance 2.0", "filetype pdb", "output packed.pdb", f"seed {cfg['seed']}",
             "pbc " + limits, "nloop 1000", "maxit 30"]
    for species in manifest["species"]:
        lines += [f"structure molecules/{species['id']}/input.pdb", f"  number {species['count']}",
                  "  inside box 0. 0. 0. " + limits, "end structure", ""]
    (out / "packmol.inp").write_text("\n".join(lines) + "\n")
    (out / "packed.pdb").unlink(missing_ok=True)
    print("Packing water and organic solvent throughout the same periodic cube...", flush=True)
    with (out / "packmol.inp").open() as inp, (out / "packmol.log").open("w") as log:
        result = subprocess.run(["packmol"], stdin=inp, stdout=log, stderr=subprocess.STDOUT, cwd=out)
    if result.returncode or "Success!" not in (out / "packmol.log").read_text() or not (out / "packed.pdb").is_file():
        raise RuntimeError(f"Packmol did not converge; inspect {out / 'packmol.log'}")
    records = [line for line in (out / "packed.pdb").read_text().splitlines()
               if line.startswith(("ATOM  ", "HETATM"))]
    if len(records) != manifest["atom_count"]:
        raise RuntimeError("Packmol atom count differs from the manifest")
    expected_names = [name for s in manifest["species"] for _ in range(s["count"]) for name in s["atom_names"]]
    if [line[12:16].strip() for line in records] != expected_names:
        raise RuntimeError("Packmol changed molecule or atom order")
    xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in records])
    if not np.isfinite(xyz).all():
        raise RuntimeError("Non-finite packed coordinates")
    return xyz


def write_mdps(cfg, out):
    folder = out / "mdp"
    folder.mkdir(exist_ok=True)
    common = """cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
vdw-modifier = Potential-shift
DispCorr = EnerPres
pbc = xyz
fourierspacing = 0.12
pme-order = 4
ewald-rtol = 1e-5
"""
    (folder / "em.mdp").write_text("integrator = steep\nnsteps = 50000\nemtol = 1000\nemstep = 0.001\nconstraints = none\n" + common)
    for stage in ("nvt", "npt", "production"):
        dt = .001 if stage == "nvt" else .002
        duration = cfg["run"][stage + "_ps"] if stage != "production" else cfg["run"]["production_ns"] * 1000
        tinit = 0 if stage == "nvt" else cfg["run"]["nvt_ps"]
        if stage == "production":
            tinit += cfg["run"]["npt_ps"]
        stride = round((1 if stage == "nvt" else 2 if stage == "npt" else 10) / dt)
        text = f"integrator = md\ndt = {dt}\nnsteps = {steps_for(duration, dt)}\ntinit = {tinit}\ninit-step = 0\n" + common
        text += f"""constraints = h-bonds
constraint-algorithm = lincs
lincs-order = 4
tcoupl = v-rescale
tc-grps = System
tau-t = 0.5
ref-t = {cfg['temperature_K']}
nstenergy = 500
nstlog = 500
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = {stride}
compressed-x-precision = 1000
comm-mode = Linear
nstcomm = 100
"""
        if stage == "nvt":
            text += f"continuation = no\ngen-vel = yes\ngen-temp = {cfg['temperature_K']}\ngen-seed = {cfg['seed']}\npcoupl = no\n"
        else:
            text += "continuation = yes\ngen-vel = no\npcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 5.0\nref-p = 1.0\ncompressibility = 4.5e-5\n"
        (folder / f"{stage}.mdp").write_text(text)


def write_index(out, manifest):
    groups = [("System", 0, manifest["atom_count"])]
    groups += [("Water" if s["kind"] == "water" else "Organic", s["atom_start_0based"], s["atom_stop_0based"])
               for s in manifest["species"]]
    with (out / "index.ndx").open("w") as handle:
        for name, start, stop in groups:
            handle.write(f"[ {name} ]\n")
            for first in range(start + 1, stop + 1, 15):
                handle.write(" ".join(str(i) for i in range(first, min(first + 15, stop + 1))) + "\n")


def assemble(out):
    import numpy as np
    import parmed as pmd
    from prepare_case import _check_periodic_contacts

    cfg, manifest = read_prepared(out)
    report = json.loads((out / "parameterization.json").read_text())
    if report["configuration_sha256"] != manifest["configuration_sha256"]:
        raise ValueError("Parameterization does not match the prepared configuration")
    (out / "READY.json").unlink(missing_ok=True)
    system, molecule_index, molecule_number = None, [], 0
    for species in manifest["species"]:
        folder = out / "molecules" / species["id"]
        unit = pmd.load_file(str(folder / "parameters.prmtop"), xyz=str(folder / "parameters.inpcrd"))
        if [a.name for a in unit.atoms] != species["atom_names"]:
            raise RuntimeError(f"Amber changed atom names/order for {species['id']}")
        if [a.atomic_number for a in unit.atoms] != species["atomic_numbers"]:
            raise RuntimeError(f"Amber changed atom identity for {species['id']}")
        if abs(sum(a.charge for a in unit.atoms)) > 1e-5:
            raise RuntimeError(f"Non-neutral topology for {species['id']}")
        block = unit * species["count"]
        system = block if system is None else system + block
        for _ in range(species["count"]):
            molecule_index.extend([molecule_number] * len(unit.atoms))
            molecule_number += 1
    xyz = pack(cfg, out, manifest)
    if len(system.atoms) != manifest["atom_count"]:
        raise RuntimeError("Topology and coordinate atom counts differ")
    box = np.array(cfg["box_nm"]) * 10
    _check_periodic_contacts(np.mod(xyz, box), box, molecule_index)
    system.coordinates = xyz
    system.box = [*box, 90., 90., 90.]
    system.save(str(out / "system.gro"), overwrite=True)
    system.save(str(out / "topol.top"), overwrite=True)
    write_mdps(cfg, out)
    write_index(out, manifest)
    ready = {"configuration_sha256": manifest["configuration_sha256"],
             "atom_count": len(system.atoms), "molecule_count": molecule_number,
             "net_charge_e": sum(a.charge for a in system.atoms),
             "status": "binary initial inputs assembled; minimization and MD still required",
             "input_sha256": {name: hashlib.sha256((out / name).read_bytes()).hexdigest()
                              for name in ("system.gro", "topol.top", "index.ndx", "manifest.json", "binary_config.json")},
             "notes": ["Both species initially occupy the same full cube without restraints.",
                       "Organic:water input is a pure-volume ratio, not a mole ratio.",
                       "Initial box volume is an additive-volume estimate; isotropic NPT relaxes it.",
                       "GAFF2/TIP3P binary miscibility requires validation against experiment."]}
    write_json(out / "READY.json", ready)
    print(f"Prepared {len(system.atoms)} atoms in {out}; stages em, nvt, npt, production", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    for stage in ("plan", "structures"):
        child = subparsers.add_parser(stage, allow_abbrev=False)
        child.add_argument("--solvent", type=canonical_solvent, required=True)
        child.add_argument("--output", type=Path, required=(stage == "structures"))
        child.add_argument("--box-nm", type=positive_number, default=6.0)
        child.add_argument("--organic-volume", type=positive_number, default=1.5,
                           help="Pure organic input volume; same units as --water-volume")
        child.add_argument("--water-volume", type=positive_number, default=10.0)
        child.add_argument("--temperature", type=positive_number, default=298.15)
        child.add_argument("--production-ns", type=positive_number, default=20.0)
        child.add_argument("--nvt-ps", type=positive_number, default=100.0)
        child.add_argument("--npt-ps", type=positive_number, default=1000.0)
        child.add_argument("--seed", type=int, default=1729)
    for stage in ("parameters", "assemble"):
        subparsers.add_parser(stage, allow_abbrev=False).add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.stage in ("plan", "structures"):
            cfg = {"schema_version": SCHEMA_VERSION, "solvent": args.solvent,
                   "box_nm": [args.box_nm] * 3, "organic_volume": args.organic_volume,
                   "water_volume": args.water_volume, "temperature_K": args.temperature,
                   "seed": args.seed, "run": {"nvt_ps": args.nvt_ps, "npt_ps": args.npt_ps,
                                               "production_ns": args.production_ns}}
            validate_config(cfg)
            if args.stage == "plan":
                print(json.dumps(make_plan(cfg), indent=2, allow_nan=False))
            else:
                structures(cfg, args.output.resolve())
        elif args.stage == "parameters":
            parameters(args.output.resolve())
        else:
            assemble(args.output.resolve())
    except (ValueError, RuntimeError, OSError, KeyError, ImportError) as error:
        print(f"Binary preparation failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
