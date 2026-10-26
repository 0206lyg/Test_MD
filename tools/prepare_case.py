#!/usr/bin/env python3
"""Generate auditable AA PLGA / drug / mixed-solvent PoC inputs.

Stages deliberately use separate environments: structures (RDKit), parameters
(AmberTools), assemble (ParmEd + Packmol). No calculations run on import.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

SCHEMA = 1
# Densities are room-temperature pure-liquid initialization estimates only.
# Molecular counts use ideal additive volumes; NPT adjusts the mixture density.
SOLVENTS = {
    "DMSO": {"smiles": "CS(=O)C", "resname": "DMS", "density": 1.095},
    "DMF": {"smiles": "CN(C)C=O", "resname": "DMF", "density": 0.944},
    "acetone": {"smiles": "CC(=O)C", "resname": "ACE", "density": 0.785},
    "ACN": {"smiles": "CC#N", "resname": "ACN", "density": 0.786},
}
DRUGS = {
    "lawsone": {"smiles": "O=C1C(O)=CC(=O)c2ccccc12", "resname": "LAW", "formula": "C10H6O3"},
    "menadione": {"smiles": "CC1=CC(=O)c2ccccc2C1=O", "resname": "MEN", "formula": "C11H8O2"},
    "plumbagin": {"smiles": "CC1=CC(=O)C2=C(C1=O)C=CC=C2O", "resname": "PLU", "formula": "C11H8O3"},
}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_config(path):
    cfg = json.loads(Path(path).read_text())
    required = {"schema_version", "box_nm", "temperature_K", "polymer", "mixture", "drug", "run"}
    if set(cfg) != required or cfg["schema_version"] != SCHEMA:
        raise ValueError(f"Expected configuration keys {sorted(required)} and schema_version={SCHEMA}")
    if len(cfg["box_nm"]) != 3 or min(cfg["box_nm"]) <= 2.4:
        raise ValueError("Use three box lengths greater than 2.4 nm")
    poly, mix, drug = cfg["polymer"], cfg["mixture"], cfg["drug"]
    if poly["repeat_count"] < 3 or poly["chains"] < 1:
        raise ValueError("At least one PLGA chain with three repeats is required")
    if poly["end_groups"] != "hydroxyl_carboxyl":
        raise ValueError("This force-field implementation supports neutral hydroxyl/carboxyl ends")
    if mix["solvent"] not in SOLVENTS or drug["name"] not in DRUGS:
        raise ValueError(f"Solvents: {list(SOLVENTS)}; drugs: {list(DRUGS)}")
    if not 0 < mix["water_mole_fraction"] < 1:
        raise ValueError("water_mole_fraction must be strictly between zero and one")
    if not 0 < mix["initial_solvent_fill_fraction"] <= 1 or mix["initial_polymer_density_g_cm3"] <= 0:
        raise ValueError("Invalid initial density/fill fraction")
    if drug["count"] < 1 or int(drug["count"]) != drug["count"]:
        raise ValueError("drug.count must be a positive integer")
    if drug["initial_region"] not in ("solvent", "polymer", "anywhere"):
        raise ValueError("drug.initial_region must be solvent, polymer, or anywhere")
    if cfg["temperature_K"] <= 0 or min(cfg["run"].values()) <= 0:
        raise ValueError("Temperature and stage lengths must be positive")
    return cfg


def _pdb(mol, path, resname, water=False):
    from rdkit import Chem
    from mdprep_forcefield import _atom_name
    named = Chem.Mol(mol)
    for i, atom in enumerate(named.GetAtoms()):
        name = (["O", "H1", "H2"][i] if water else _atom_name(i, atom.GetSymbol()))
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{name:>4}")
        info.SetResidueName(resname)
        info.SetResidueNumber(1)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    Chem.MolToPDBFile(named, str(path))


def _small_molecule(smiles, seed, water=False):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if water:
        conf = Chem.Conformer(3)
        angle = math.radians(104.52)
        for i, p in enumerate([(0., 0., 0.), (.9572, 0., 0.), (.9572*math.cos(angle), .9572*math.sin(angle), 0.)]):
            conf.SetAtomPosition(i, p)
        mol.AddConformer(conf)
    else:
        settings = AllChem.ETKDGv3()
        settings.randomSeed = seed
        if AllChem.EmbedMolecule(mol, settings) < 0:
            raise RuntimeError(f"Cannot embed {smiles}")
        if not AllChem.MMFFHasAllMoleculeParams(mol) or AllChem.MMFFOptimizeMolecule(mol, maxIters=3000) != 0:
            raise RuntimeError(f"Cannot obtain a converged small-molecule geometry for {smiles}")
        xyz = mol.GetConformer().GetPositions()
        xyz -= xyz.mean(axis=0)
        for i, p in enumerate(xyz):
            mol.GetConformer().SetAtomPosition(i, p)
    return mol


def _write_molecule(folder, mol, resname, water=False, mapping=None):
    from rdkit import Chem
    from mdprep_forcefield import export_mol2
    folder.mkdir(parents=True, exist_ok=True)
    export_mol2(mol, folder / "input.mol2", resname)
    with Chem.SDWriter(str(folder / "input.sdf")) as writer:
        writer.SetForceV3000(mol.GetNumAtoms() > 999)
        writer.write(mol)
    _pdb(mol, folder / "input.pdb", resname, water)
    if mapping is not None:
        write_json(folder / "atom_map.json", mapping)


def _reference_polymer(md_root, configuration_key):
    """Reuse a prepared case's portable PDB coordinates across RDKit versions."""
    from rdkit import Chem
    from mdprep_polymer import build_graph
    for marker in sorted((md_root / "cases").glob("*_input/manifest.json")):
        source = marker.parent
        if not (source / "READY.json").is_file():
            continue
        manifest = json.loads(marker.read_text())
        if manifest.get("polymer_configuration_sha256") != configuration_key:
            continue
        metadata = json.loads((source / "polymer.json").read_text())
        polymers = [s for s in manifest["species"] if s["kind"] == "polymer"]
        if len(polymers) != len(metadata["chains"]):
            raise RuntimeError(f"Invalid polymer reference metadata: {source}")
        combined = None
        for species, chain in zip(polymers, metadata["chains"]):
            mol, _ = build_graph(chain["sequence"], chain["stereo"], chain["end_groups"])
            xyz = []
            for line in (source / "molecules" / species["id"] / "input.pdb").read_text().splitlines():
                if line.startswith(("ATOM  ", "HETATM")):
                    xyz.append(tuple(float(line[k:k+8]) for k in (30, 38, 46)))
            if len(xyz) != mol.GetNumAtoms() or len(xyz) != chain["atom_count"]:
                raise RuntimeError(f"Invalid polymer reference coordinates: {source}")
            conf = Chem.Conformer(len(xyz))
            for i, point in enumerate(xyz):
                conf.SetAtomPosition(i, point)
            mol.AddConformer(conf)
            combined = mol if combined is None else Chem.CombineMols(combined, mol)
        print(f"Reusing the exact PLGA PDB coordinates from {source}", flush=True)
        metadata["coordinate_reference_case"] = str(source.relative_to(md_root))
        return combined, metadata
    return None


def structures(cfg, out):
    from rdkit import Chem, rdBase
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from mdprep_polymer import build_graph, build_polymer_system
    from mdprep_forcefield import prepare_fragments
    expected = digest(cfg)
    marker = out / "manifest.json"
    if marker.exists():
        old = json.loads(marker.read_text())
        if old["configuration_sha256"] != expected:
            raise ValueError("Output belongs to a different configuration. Use a new output directory.")
        print(f"Reusing matching molecular structures: {out}", flush=True)
        return
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "case.json", cfg)
    kwargs = dict(cfg["polymer"])
    kwargs["box_nm"] = cfg["box_nm"]
    polymer_configuration_key = digest({"polymer": cfg["polymer"], "box_nm": cfg["box_nm"]})
    builder_version = hashlib.sha256((Path(__file__).parent / "mdprep_polymer.py").read_bytes()).hexdigest()
    polymer_key = digest({"settings": kwargs, "rdkit": rdBase.rdkitVersion, "builder": builder_version})
    md_root = Path(os.environ.get("MD_ROOT", Path(__file__).resolve().parents[1]))
    cache = md_root / "cache" / "mdprep-polymer" / polymer_key
    if (cache / "polymer.rdmol").exists() and (cache / "metadata.json").exists():
        polymer = Chem.Mol((cache / "polymer.rdmol").read_bytes())
        metadata = json.loads((cache / "metadata.json").read_text())
        if polymer is None or polymer.GetNumAtoms() != metadata["atom_count"]:
            raise RuntimeError(f"Invalid polymer cache: {cache}")
        print("Reusing identical PLGA coordinates from the molecular cache.", flush=True)
    else:
        reference = _reference_polymer(md_root, polymer_configuration_key)
        if reference is not None:
            polymer, metadata = reference
        else:
            print("Growing PLGA chains (no full-chain QM or distance-geometry embedding)...", flush=True)
            polymer, metadata = build_polymer_system(**kwargs)
        cache.mkdir(parents=True, exist_ok=True)
        temporary = cache / f"polymer.{os.getpid()}.tmp"
        temporary.write_bytes(polymer.ToBinary(Chem.PropertyPickleOptions.CoordsAsDouble))
        temporary.replace(cache / "polymer.rdmol")
        write_json(cache / "metadata.json", metadata)
    write_json(out / "polymer.json", metadata)
    species = []
    polymer_xyz = polymer.GetConformer().GetPositions()
    for i, chain in enumerate(metadata["chains"]):
        # Rebuild each known chemical graph rather than asking RDKit to split
        # and reassign stereochemistry on a many-thousand-atom combined Mol.
        mol, _ = build_graph(chain["sequence"], chain["stereo"], chain["end_groups"])
        offset, count = chain["atom_offset"], chain["atom_count"]
        if mol.GetNumAtoms() != count or any(mol.GetAtomWithIdx(j).GetAtomicNum() != polymer.GetAtomWithIdx(offset+j).GetAtomicNum() for j in range(count)):
            raise RuntimeError("Polymer graph and coordinate atom order differ")
        conf = Chem.Conformer(count)
        for j, point in enumerate(polymer_xyz[offset:offset+count]):
            conf.SetAtomPosition(j, point)
        mol.AddConformer(conf)
        name = f"PLG{i:03d}"
        _write_molecule(out / "molecules" / name, mol, f"P{i:02d}", mapping=chain)
        species.append({"id": name, "kind": "polymer", "count": 1, "formal_charge": 0,
                        "atoms_per_molecule": mol.GetNumAtoms(), "mass_g_mol": Descriptors.MolWt(mol)})
    drug = DRUGS[cfg["drug"]["name"]]
    solvent = SOLVENTS[cfg["mixture"]["solvent"]]
    drug_mol = _small_molecule(drug["smiles"], 4281)
    if rdMolDescriptors.CalcMolFormula(drug_mol) != drug["formula"]:
        raise RuntimeError("Drug molecular formula does not match the declared identity")
    solvent_mol = _small_molecule(solvent["smiles"], 4282)
    water_mol = _small_molecule("O", 4283, water=True)
    # 1 g/cm3 = 602.214076 u/nm3. Water/solvent mole fraction excludes PLGA/drug.
    poly_volume = sum(s["mass_g_mol"] for s in species) / (602.214076 * cfg["mixture"]["initial_polymer_density_g_cm3"])
    total_volume = math.prod(cfg["box_nm"])
    available = (total_volume - poly_volume) * cfg["mixture"]["initial_solvent_fill_fraction"]
    if available <= 0:
        raise ValueError("Polymer mass exceeds the configured box volume at the initialization density")
    xw = cfg["mixture"]["water_mole_fraction"]
    vw = Descriptors.MolWt(water_mol) / (602.214076 * .997)
    vs = Descriptors.MolWt(solvent_mol) / (602.214076 * solvent["density"])
    count = round(available / (xw*vw + (1-xw)*vs))
    nw = round(xw*count)
    ns = count-nw
    if min(nw, ns) < 1:
        raise ValueError("The requested mixture must contain both water and organic solvent")
    for name, kind, mol, n, res in [("DRUG", "drug", drug_mol, cfg["drug"]["count"], drug["resname"]),
                                  ("ORG", "solvent", solvent_mol, ns, solvent["resname"]),
                                  ("WAT", "water", water_mol, nw, "WAT")]:
        _write_molecule(out / "molecules" / name, mol, res, water=(kind == "water"))
        species.append({"id": name, "kind": kind, "count": n, "formal_charge": 0,
                        "atoms_per_molecule": mol.GetNumAtoms(), "mass_g_mol": Descriptors.MolWt(mol),
                        "smiles": Chem.MolToSmiles(Chem.RemoveHs(mol))})
    print("Preparing short charge-transfer reference oligomers...", flush=True)
    reference_case = metadata.get("coordinate_reference_case")
    prepare_fragments(out, reference_case=(md_root / reference_case) if reference_case else None)
    manifest = {"schema_version": SCHEMA, "configuration_sha256": expected,
                "polymer_configuration_sha256": polymer_configuration_key,
                "polymer_cache_key": polymer_key,
                "species": species, "atom_count": sum(s["count"]*s["atoms_per_molecule"] for s in species),
                "water_mole_fraction_actual": nw/(nw+ns),
                "initial_polymer_volume_nm3": poly_volume,
                "status": "unequilibrated neutral-species PoC; fixed composition, no solvent removal",
                "force_field_status": "GAFF2/AM1-BCC with approximate polymer charge transfer; TIP3P water"}
    write_json(marker, manifest)
    print(f"Structures ready: {manifest['atom_count']} atoms; water {nw}; organic solvent {ns}", flush=True)


def _pack(cfg, out, manifest):
    import numpy as np
    if shutil.which("packmol") is None:
        raise RuntimeError("Packmol is missing from the preparation environment")
    box = np.array(cfg["box_nm"])*10
    slab = np.array(cfg["polymer"]["slab_nm"])*10
    low, high = (box-slab)/2, (box+slab)/2
    mix = cfg["mixture"]
    lines = [f"tolerance {mix['packing_tolerance_angstrom']}", "filetype pdb", "output packed.pdb",
             f"seed {mix['packing_seed']}", "pbc " + " ".join(map(str, box)), "nloop 200", "maxit 30"]
    for s in manifest["species"]:
        lines += [f"structure molecules/{s['id']}/input.pdb", f"  number {s['count']}"]
        if s["kind"] == "polymer":
            # Coordinates already positioned; Packmol rotates neither chains nor their monomers.
            lines += ["  fixed 0. 0. 0. 0. 0. 0."]
        else:
            lines += ["  inside box 0. 0. 0. " + " ".join(map(str, box))]
            if s["kind"] == "drug" and cfg["drug"]["initial_region"] != "anywhere":
                word = "outside" if cfg["drug"]["initial_region"] == "solvent" else "inside"
                lines += [f"  {word} box " + " ".join(map(str, np.r_[low, high]))]
        lines += ["end structure", ""]
    (out / "packmol.inp").write_text("\n".join(lines))
    packing_key = digest([lines] + [(out / "molecules" / s["id"] / "input.pdb").read_text() for s in manifest["species"]])
    cached = json.loads((out / "packing.json").read_text()) if (out / "packing.json").exists() else {}
    packed = out / "packed.pdb"
    valid_cache = packed.exists() and cached.get("input_sha256") == packing_key and cached.get("output_sha256") == hashlib.sha256(packed.read_bytes()).hexdigest()
    if not valid_cache:
        packed.unlink(missing_ok=True)
        print("Packing drug and mixed solvent around the fixed initial PLGA coordinates...", flush=True)
        with (out / "packmol.inp").open() as inp, (out / "packmol.log").open("w") as log:
            result = subprocess.run(["packmol"], stdin=inp, stdout=log, stderr=subprocess.STDOUT, cwd=out)
        log = (out / "packmol.log").read_text()
        if result.returncode or "Success!" not in log or not packed.is_file():
            raise RuntimeError(f"Packmol did not converge. See {out / 'packmol.log'}; no runnable inputs were marked ready.")
        write_json(out / "packing.json", {"input_sha256": packing_key, "output_sha256": hashlib.sha256(packed.read_bytes()).hexdigest()})
    coords = []
    for line in (out / "packed.pdb").read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    xyz = np.array(coords)
    if xyz.shape != (manifest["atom_count"], 3) or not np.isfinite(xyz).all():
        raise RuntimeError("Packed coordinates do not match the manifest atom count")
    # Fixed chains must be exactly preserved, up to PDB rounding.
    offset = 0
    for s in manifest["species"]:
        if s["kind"] != "polymer":
            break
        original = []
        for line in (out / "molecules" / s["id"] / "input.pdb").read_text().splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                original.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if not np.allclose(xyz[offset:offset+len(original)], original, atol=.002):
            raise RuntimeError("Packmol changed the supposedly fixed PLGA coordinates")
        offset += len(original)
    return xyz


def _mdps(cfg, out):
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
    for stage in ("nvt", "npt", "benchmark"):
        dt = .001 if stage == "nvt" else .002
        steps = round(cfg["run"][stage+"_ps"] / dt)
        if steps < 1:
            raise ValueError("Each MD stage must contain at least one timestep")
        text = f"integrator = md\ndt = {dt}\nnsteps = {steps}\n" + common
        text += f"""constraints = h-bonds
constraint-algorithm = lincs
lincs-order = 4
tcoupl = v-rescale
tc-grps = System
tau-t = 0.5
ref-t = {cfg['temperature_K']}
nstenergy = 1000
nstlog = 1000
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = 1000
compressed-x-precision = 1000
comm-mode = Linear
nstcomm = 100
"""
        if stage == "nvt":
            text += f"continuation = no\ngen-vel = yes\ngen-temp = {cfg['temperature_K']}\ngen-seed = {cfg['polymer']['seed']+30000}\npcoupl = no\n"
        else:
            # Keep slab lateral area, allow normal dimension/density to relax.
            text += "continuation = yes\ngen-vel = no\npcoupl = C-rescale\npcoupltype = semiisotropic\ntau-p = 5.0\nref-p = 1.0 1.0\ncompressibility = 0.0 4.5e-5\n"
        (folder / (stage+".mdp")).write_text(text)


def assemble(cfg, out):
    import numpy as np
    import parmed as pmd
    manifest = json.loads((out / "manifest.json").read_text())
    if manifest["configuration_sha256"] != digest(cfg):
        raise ValueError("Configuration differs from prepared structures")
    (out / "READY.json").unlink(missing_ok=True)
    system = None
    molecule_index = []
    mol_number = 0
    for species in manifest["species"]:
        folder = out / "molecules" / species["id"]
        unit = pmd.load_file(str(folder / "parameters.prmtop"), xyz=str(folder / "parameters.inpcrd"))
        if len(unit.atoms) != species["atoms_per_molecule"]:
            raise RuntimeError(f"Amber changed the number of atoms for {species['id']}")
        if abs(sum(a.charge for a in unit.atoms)) > 1e-5:
            raise RuntimeError(f"Non-neutral topology for {species['id']}")
        block = unit * species["count"]
        system = block if system is None else system + block
        for _ in range(species["count"]):
            molecule_index.extend([mol_number]*len(unit.atoms))
            mol_number += 1
    xyz = _pack(cfg, out, manifest)
    if len(system.atoms) != len(xyz):
        raise RuntimeError("Topology/coordinate atom counts differ")
    box_a = np.array(cfg["box_nm"])*10
    xyz_periodic = np.mod(xyz, box_a)
    # Detect severe intermolecular contacts, including opposite periodic faces.
    _check_periodic_contacts(xyz_periodic, box_a, molecule_index)
    system.coordinates = xyz
    system.box = [*box_a, 90., 90., 90.]
    system.save(str(out / "system.gro"), overwrite=True)
    system.save(str(out / "topol.top"), overwrite=True)
    _mdps(cfg, out)
    ready = {"atom_count": len(system.atoms), "molecule_count": mol_number,
             "net_charge_e": sum(a.charge for a in system.atoms),
             "configuration_sha256": manifest["configuration_sha256"],
             "status": "initial inputs assembled; requires energy minimization and relaxation",
             "notes": ["All species neutral; this is a declared initial chemical-state assumption.",
                       "Initial polymer conformation and mixture volume are not equilibrated states.",
                       "Short MD is not a converged partition/free-energy calculation."]}
    write_json(out / "READY.json", ready)
    print(f"Prepared {len(system.atoms)} atoms: {out / 'system.gro'} and {out / 'topol.top'}", flush=True)


def _check_periodic_contacts(xyz, box, molecule_index, cutoff=1.2):
    """Linear cell-list check using existing NumPy, without a SciPy dependency."""
    import itertools
    import numpy as np
    bins = np.floor(box/cutoff).astype(int)
    width = box/bins
    keys = np.floor(xyz/width).astype(int)
    cells = {}
    positions = xyz.tolist()
    for i, key in enumerate(keys):
        p = positions[i]
        for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3):
            cell = ((int(key[0])+dx)%int(bins[0]), (int(key[1])+dy)%int(bins[1]), (int(key[2])+dz)%int(bins[2]))
            for j in cells.get(cell, ()):
                if molecule_index[i] == molecule_index[j]:
                    continue
                distance2 = 0.
                for k in range(3):
                    delta = abs(p[k]-positions[j][k])
                    delta = min(delta, float(box[k])-delta)
                    distance2 += delta*delta
                if distance2 < cutoff*cutoff:
                    raise RuntimeError(f"Intermolecular contact {i+1}/{j+1} below {cutoff} A across periodic box; inspect packed.pdb")
        cells.setdefault(tuple(map(int, key)), []).append(i)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("structures", "parameters", "assemble"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg, out = read_config(args.config), args.output.resolve()
    if args.stage == "structures":
        structures(cfg, out)
    elif args.stage == "parameters":
        from mdprep_forcefield import parameterize
        manifest = json.loads((out / "manifest.json").read_text())
        if manifest["configuration_sha256"] != digest(cfg):
            raise ValueError("Configuration differs from prepared structures")
        parameterize(out)
    else:
        assemble(cfg, out)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")
