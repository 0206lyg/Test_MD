#!/usr/bin/env python3
"""Provisional GAFF2/AM1-BCC PoC parameters, with auditable polymer charge transfer.

Amber commands run in the Amber module environment; RDKit is imported only by the
coordinate/fragment preparation helpers. No QM calculation is run on a long chain.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor


METHOD = "gaff2-am1bcc-pentamer-neighbor-transfer-v1"
SMALL_MOLECULE_METHOD = "gaff2-am1bcc-neutral-roundoff-v1"


def _atom_name(index, element):
    alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    width = 4-len(element)
    if index >= 36 ** width:
        raise ValueError("Molecule exceeds the unique four-character atom-name limit")
    n = index
    digits = ""
    for _ in range(width):
        digits = alphabet[n % 36] + digits
        n //= 36
    return element + digits


def export_mol2(mol, path, resname="MOL", charges=None):
    """Write explicit connectivity and deterministic names in RDKit atom order."""
    from rdkit import Chem
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if mol.GetNumConformers() != 1:
        raise ValueError("A molecule must have exactly one coordinate conformer")
    if charges is None:
        charges = [0.0] * mol.GetNumAtoms()
    if len(charges) != mol.GetNumAtoms():
        raise ValueError("Charge and atom counts differ")
    xyz = mol.GetConformer()
    lines = ["@<TRIPOS>MOLECULE", resname,
             f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0", "SMALL", "USER_CHARGES", "", "@<TRIPOS>ATOM"]
    for i, atom in enumerate(mol.GetAtoms()):
        element = atom.GetSymbol()
        if element in ("C", "N"):
            suffix = "ar" if atom.GetIsAromatic() else {Chem.HybridizationType.SP: "1", Chem.HybridizationType.SP2: "2"}.get(atom.GetHybridization(), "3")
            atomtype = element + "." + suffix
            if element == "N" and any(b.GetOtherAtom(atom).GetSymbol() == "C" and any(b2.GetBondTypeAsDouble() == 2 and b2.GetOtherAtom(b.GetOtherAtom(atom)).GetSymbol() == "O" for b2 in b.GetOtherAtom(atom).GetBonds()) for b in atom.GetBonds()):
                atomtype = "N.am"
        elif element == "O":
            atomtype = "O.2" if any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds()) else "O.3"
        elif element == "S":
            oxo = sum(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(atom).GetSymbol() == "O" for b in atom.GetBonds())
            atomtype = "S.O2" if oxo == 2 else "S.O" if oxo == 1 else "S.3"
        else:
            atomtype = element
        p = xyz.GetAtomPosition(i)
        lines.append(f"{i+1:7d} {_atom_name(i, element):4s} {p.x:14.7f} {p.y:14.7f} {p.z:14.7f} {atomtype:6s} 1 {resname:3s} {charges[i]: .10f}")
    lines.append("@<TRIPOS>BOND")
    for i, bond in enumerate(mol.GetBonds()):
        order = "ar" if bond.GetIsAromatic() else str(int(bond.GetBondTypeAsDouble()))
        lines.append(f"{i+1:7d} {bond.GetBeginAtomIdx()+1:7d} {bond.GetEndAtomIdx()+1:7d} {order}")
    lines.extend(["@<TRIPOS>SUBSTRUCTURE", f"1 {resname} 1 GROUP 0 **** 0 ROOT", ""])
    path.write_text("\n".join(lines))


def _fragment_molecule(sequence, seed):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from mdprep_polymer import build_graph
    mol, mapping = build_graph(sequence, ["L" if x == "LA" else None for x in sequence])
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.maxIterations = 2000
    if AllChem.EmbedMolecule(mol, params) < 0:
        raise RuntimeError("Cannot embed the reduced PLGA fragment")
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise RuntimeError("MMFF geometry parameters are missing for a reduced PLGA fragment")
    converged = AllChem.MMFFOptimizeMolecule(mol, maxIters=3000)
    if converged != 0:
        raise RuntimeError("Reduced PLGA fragment MMFF geometry did not converge")
    mapping["stereochemistry_note"] = "all LA centres S; transferred charges shared by D/L mirror centres"
    return mol, mapping


def prepare_fragments(case_dir, reference_case=None):
    """Generate neutral pentamers covering all L/G nearest-neighbor triplets."""
    from rdkit import Chem
    base = Path(case_dir) / "charge_fragments"
    base.mkdir(parents=True, exist_ok=True)
    if reference_case is not None:
        reference = Path(reference_case) / "charge_fragments"
        specification = json.loads((reference / "specification.json").read_text())
        if specification.get("method") == METHOD:
            expected = [f"fragment{i:02d}" for i in range(10)]
            if specification.get("fragments") != expected:
                raise RuntimeError("Invalid reference charge-fragment specification")
            for name in expected:
                destination = base / name
                destination.mkdir(exist_ok=True)
                for filename in ("input.mol2", "input.sdf", "atom_map.json"):
                    shutil.copy2(reference / name / filename, destination / filename)
            shutil.copy2(reference / "specification.json", base / "specification.json")
            print("Reusing the same polymer charge-reference geometries.", flush=True)
            return
    sequences = ["G"+"".join(middle)+"G" for middle in itertools.product("GL", repeat=3)] + ["LLLLL", "LGLGL"]
    specification = {"method": METHOD, "interior_repeat": 2, "fragments": []}
    for index, letters in enumerate(sequences):
        name = f"fragment{index:02d}"
        folder = base / name
        folder.mkdir(exist_ok=True)
        seq = ["GA" if letter == "G" else "LA" for letter in letters]
        mol, mapping = _fragment_molecule(seq, 7129 + index)
        export_mol2(mol, folder / "input.mol2", "FRG")
        Chem.MolToMolFile(mol, str(folder / "input.sdf"))
        (folder / "atom_map.json").write_text(json.dumps(mapping, indent=2)+"\n")
        specification["fragments"].append(name)
    (base / "specification.json").write_text(json.dumps(specification, indent=2)+"\n")


def _mol2_atoms(path):
    atoms = []
    section = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line[9:]
        elif section == "ATOM" and line.strip():
            fields = line.split()
            if len(fields) < 9:
                raise ValueError(f"Incomplete MOL2 atom record in {path}")
            atoms.append({"name": fields[1], "charge": float(fields[8]), "type": fields[5]})
    if not atoms:
        raise ValueError(f"No MOL2 atoms in {path}")
    return atoms


def _run(argv, folder, log_name):
    with (folder / log_name).open("w") as log:
        log.write("COMMAND: " + " ".join(map(str, argv)) + "\n")
        log.flush()
        result = subprocess.run(list(map(str, argv)), cwd=folder, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{argv[0]} failed ({result.returncode}); see {folder / log_name}")


def _cache_signature(folder, polymer_charges=None):
    signature = hashlib.sha256((folder / "input.mol2").read_bytes())
    signature.update((METHOD if polymer_charges is not None else SMALL_MOLECULE_METHOD).encode())
    signature.update(json.dumps(polymer_charges).encode())
    amber_home = Path(os.environ.get("AMBERHOME", "/missing"))
    for relative in ("dat/leap/parm/gaff2.dat", "dat/antechamber/BCCPARM.DAT"):
        parameter_file = amber_home / relative
        if parameter_file.is_file():
            signature.update(parameter_file.read_bytes())
    executable = Path(shutil.which("antechamber")).resolve()
    signature.update(str(executable).encode())
    signature.update(executable.read_bytes())
    return signature.hexdigest()


def _cache_root():
    return Path(os.environ.get("MD_ROOT", Path(__file__).resolve().parents[1])) / "cache" / "mdprep-charges"


def _save_charge_cache(folder, signature):
    typed_file = folder / "typed.mol2"
    metadata = {"input_signature": signature, "typed_sha256": hashlib.sha256(typed_file.read_bytes()).hexdigest()}
    (folder / "antechamber_cache.json").write_text(json.dumps(metadata, indent=2)+"\n")
    destination = _cache_root() / signature
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        with tempfile.TemporaryDirectory(prefix="pending-", dir=destination.parent) as temporary:
            staging = Path(temporary) / "entry"
            staging.mkdir()
            for file in folder.iterdir():
                if file.is_file():
                    shutil.copy2(file, staging / file.name)
            try:
                staging.rename(destination)
            except OSError:
                if not destination.is_dir():
                    raise
                # Another job published the same completed calculation.


def _quantized_neutral(charges):
    units = [int(round(q*1_000_000)) for q in charges]
    remainder = -sum(units)
    if remainder:
        sign = 1 if remainder > 0 else -1
        ranked = sorted(range(len(units)), key=lambda i: sign*(charges[i]*1_000_000-units[i]), reverse=True)
        if abs(remainder) > len(ranked):
            raise RuntimeError("Unexpected charge quantization residual")
        for i in ranked[:abs(remainder)]:
            units[i] += sign
    return [u/1_000_000 for u in units]


def _correct_sqm_charge_roundoff(folder):
    """Correct only the neutrality loss from SQM's 3-decimal printed charges."""
    typed_file = folder / "typed.mol2"
    atoms = _mol2_atoms(typed_file)
    total = sum(a["charge"] for a in atoms)
    sqm = (folder / "sqm.out").read_text()
    reported = re.findall(r"Total Mulliken Charge\s*=\s*([-+0-9.eE]+)", sqm)
    if not reported or abs(float(reported[-1])) > 0.000501:
        raise RuntimeError(f"SQM did not confirm a neutral QM molecule in {folder}")
    if abs(total) > len(atoms)*0.000501:
        raise RuntimeError(f"Charge residual exceeds the SQM printed-charge rounding bound in {folder}")
    correction = -total/len(atoms)
    values = _quantized_neutral([a["charge"]+correction for a in atoms])
    shutil.copy2(typed_file, folder / "typed_before_roundoff_correction.mol2")
    output, section, index = [], None, 0
    for line in typed_file.read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line[9:]
        elif section == "ATOM" and line.strip():
            fields = line.split()
            fields[8] = f"{values[index]:.6f}"
            line = " ".join(fields)
            index += 1
        output.append(line)
    typed_file.write_text("\n".join(output)+"\n")
    (folder / "charge_roundoff_correction.json").write_text(json.dumps({
        "reason": "SQM prints individual Mulliken charges to 3 decimals; AM1-BCC preserves their rounded net charge",
        "raw_net_charge_e": total, "reported_qm_net_charge_e": float(reported[-1]),
        "uniform_per_atom_correction_e": correction, "quantized_charge_precision_e": 1e-6,
        "corrected_net_charge_e": sum(values),
    }, indent=2)+"\n")


def _antechamber(folder, polymer_charges=None):
    signature = _cache_signature(folder, polymer_charges)
    cache_file = folder / "antechamber_cache.json"
    typed_file = folder / "typed.mol2"
    shared_cache = _cache_root() / signature
    if not cache_file.is_file() and (shared_cache / "antechamber_cache.json").is_file():
        for file in shared_cache.iterdir():
            if file.is_file() and file.name not in {"input.mol2", "input.sdf", "input.pdb", "atom_map.json"}:
                shutil.copy2(file, folder / file.name)
    cached = False
    if cache_file.is_file() and typed_file.is_file():
        cache = json.loads(cache_file.read_text())
        cached = (cache.get("input_signature") == signature and
                  cache.get("typed_sha256") == hashlib.sha256(typed_file.read_bytes()).hexdigest())
    argv = ["antechamber", "-i", "input.mol2", "-fi", "mol2", "-o", "typed.mol2", "-fo", "mol2", "-at", "gaff2", "-nc", "0", "-s", "2", "-an", "n"]
    if polymer_charges is None:
        argv += ["-c", "bcc"]
    else:
        (folder / "charges.txt").write_text("\n".join(f"{q:.10f}" for q in polymer_charges)+"\n")
        argv += ["-c", "rc", "-cf", "charges.txt"]
    if not cached:
        print(f"Parameterizing {folder.name}: {'GAFF2 with transferred charges' if polymer_charges is not None else 'GAFF2 / AM1-BCC'}", flush=True)
        typed_file.unlink(missing_ok=True)
        if polymer_charges is None:
            (folder / "sqm.out").unlink(missing_ok=True)
        _run(argv, folder, "antechamber.log")
        if polymer_charges is None:
            _correct_sqm_charge_roundoff(folder)
    else:
        print(f"Reusing checked atom types/charges: {folder.name}", flush=True)
    original, typed = _mol2_atoms(folder / "input.mol2"), _mol2_atoms(folder / "typed.mol2")
    if [a["name"] for a in original] != [a["name"] for a in typed]:
        raise RuntimeError(f"Atom order/names changed in {folder}; refusing charge transfer")
    if abs(sum(a["charge"] for a in typed)) > 1e-4:
        raise RuntimeError(f"Unexpected non-neutral charge in {folder}")
    _save_charge_cache(folder, signature)
    return typed


def _charge_key(mapping, atom):
    sequence = mapping["sequence"]
    i = atom["repeat"]
    role = atom["role"]
    if atom.get("parent_role"):
        role = "H:" + atom["parent_role"]
    elif re.fullmatch(r"[A-Z]+:H[0-9]+", role):
        role = "H:" + role.split(":")[0]
    if i == 0:
        return ("first", sequence[i], role)
    if i == len(sequence)-1:
        return ("last", sequence[i], role)
    return ("middle", sequence[i-1], sequence[i], sequence[i+1], role)


def _charge_library(case_dir):
    samples = {}
    base = case_dir / "charge_fragments"
    specification = json.loads((base / "specification.json").read_text())
    if specification["method"] != METHOD:
        raise RuntimeError("Charge fragment specification is outdated; regenerate structures")
    folders = [base / name for name in specification["fragments"]]
    # Independent SQM jobs each use one thread in prepare_poc.sbatch.
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(_antechamber, folders))
    for folder, typed in zip(folders, results):
        mapping = json.loads((folder / "atom_map.json").read_text())
        if len(typed) != len(mapping["atoms"]):
            raise RuntimeError("Fragment atom metadata do not match the charged molecule")
        for atom in mapping["atoms"]:
            if atom["repeat"] not in {0, specification["interior_repeat"], len(mapping["sequence"])-1}:
                continue
            samples.setdefault(_charge_key(mapping, atom), []).append(typed[atom["index"]]["charge"])
    if not samples:
        raise RuntimeError("Missing charge fragments: run the structures preparation stage first")
    library = {key: sum(values)/len(values) for key, values in samples.items()}
    (case_dir / "charge_fragments" / "library.json").write_text(json.dumps({
        "method": METHOD,
        "status": "provisional; not experimentally validated for PLGA/drug partitioning",
        "fragment_selection": "central residue of ten pentamers for interior charges; first/last residues for neutral OH/COOH termini",
        "source_method": "Reduction-charge-transfer concept: doi:10.1021/acs.jcim.3c01691; this implementation uses Amber AM1-BCC, not AM1-BCC-ELF10",
        "templates": [{"key": list(k), "mean_charge_e": library[k], "sample_count": len(samples[k])} for k in sorted(library)],
    }, indent=2)+"\n")
    return library


def _polymer_parameters(case_dir):
    """Generate local GAFF2 additions on complete pentamers, checking duplicates."""
    base = case_dir / "charge_fragments"
    specification = json.loads((base / "specification.json").read_text())
    folders = [base / name for name in specification["fragments"]]

    def check(folder):
        _run(["parmchk2", "-i", "typed.mol2", "-f", "mol2", "-o", "missing.frcmod", "-s", "gaff2"], folder, "parmchk2.log")
        if re.search(r"ATTN|need revision", (folder / "missing.frcmod").read_text(), re.I):
            raise RuntimeError(f"Unresolved pentamer force-field parameter in {folder}")
        return folder / "missing.frcmod"

    with ThreadPoolExecutor(max_workers=2) as pool:
        files = list(pool.map(check, folders))
    sections = ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON")
    records = {section: {} for section in sections}
    origins = {}
    for file in files:
        section = None
        for line in file.read_text().splitlines():
            text = line.strip()
            if text in sections:
                section = text
                continue
            if not text or section is None:
                continue
            count = {"MASS": 1, "BOND": 2, "ANGLE": 3, "DIHE": 4, "IMPROPER": 4, "NONBON": 1}[section]
            pattern = r"^\s*" + r"\s*-\s*".join([r"([A-Za-z0-9*+]+)"]*count) + r"\s+(.+)$"
            match = re.match(pattern, line)
            if match is None:
                raise RuntimeError(f"Cannot parse {section} parameter in {file}: {line}")
            atomtypes = tuple(match.groups()[:count])
            values = []
            for token in match.groups()[count].split():
                try:
                    values.append(float(token))
                except ValueError:
                    break
            expected = {"MASS": (1, 2), "BOND": (2,), "ANGLE": (2,), "DIHE": (4,), "IMPROPER": (3,), "NONBON": (2,)}[section]
            if len(values) not in expected:
                raise RuntimeError(f"Unexpected numerical fields in {file}: {line}")
            # A Fourier torsion may contain several distinct periodicities.
            term = abs(values[-1]) if section in {"DIHE", "IMPROPER"} else None
            key = (atomtypes, term)
            old = records[section].get(key)
            if old is not None and old[0] != tuple(values):
                raise RuntimeError(f"Conflicting fragment parameters for {section} {atomtypes}: {origins[(section,key)]} and {file}")
            if old is None:
                records[section][key] = (tuple(values), line)
                origins[(section,key)] = str(file.relative_to(case_dir))
    merged = base / "polymer.frcmod"
    output = ["GAFF2 local additions from ten complete PLGA pentamers; duplicate values checked"]
    provenance = []
    for section in sections:
        output.extend([section])
        for key, (_, line) in records[section].items():
            output.append(line)
            provenance.append({"section": section, "atom_types": key[0], "source": origins[(section,key)], "parameter_line": line})
        output.append("")
    merged.write_text("\n".join(output)+"\n")
    (base / "polymer_parameter_sources.json").write_text(json.dumps(provenance, indent=2)+"\n")
    return merged


def _polymer_charges(mapping, library):
    if len(mapping["sequence"]) < 3:
        raise ValueError("PLGA charge transfer requires at least three residues")
    charges = [0.0] * len(mapping["atoms"])
    residues = {}
    for atom in mapping["atoms"]:
        key = _charge_key(mapping, atom)
        if key not in library:
            raise ValueError(f"No fragment charge available for {key}")
        charges[atom["index"]] = library[key]
        residues.setdefault(atom["repeat"], []).append(atom["index"])
    report = []
    for repeat, indices in sorted(residues.items()):
        before = sum(charges[i] for i in indices)
        adjustment = -before/len(indices)
        for i in indices:
            charges[i] += adjustment
        # Amber MOL2 writers commonly retain six decimals. Quantize each neutral
        # residue consistently so roundoff cannot accumulate along a long chain.
        units = {i: int(round(charges[i]*1_000_000)) for i in indices}
        remainder = -sum(units.values())
        if remainder:
            sign = 1 if remainder > 0 else -1
            ranked = sorted(indices, key=lambda i: sign*(charges[i]*1_000_000-units[i]), reverse=True)
            if abs(remainder) > len(ranked):
                raise RuntimeError("Unexpected charge quantization residual")
            for i in ranked[:abs(remainder)]:
                units[i] += sign
        for i in indices:
            charges[i] = units[i]/1_000_000
        report.append({"repeat": repeat, "uncorrected_charge_e": before,
                       "per_atom_correction_e": adjustment, "atom_count": len(indices),
                       "quantized_charge_precision_e": 1e-6})
    if abs(sum(charges)) > 1e-8:
        raise RuntimeError("Polymer neutrality correction failed")
    return charges, {"method": METHOD, "per_residue": report,
                     "maximum_absolute_per_atom_correction_e": max(abs(x["per_atom_correction_e"]) for x in report),
                     "net_charge_e": sum(charges)}


def _topology(folder, water=False, polymer_parameters=None):
    if water:
        commands = "source leaprc.water.tip3p\nmol = loadPdb input.pdb\nsaveAmberParm mol parameters.prmtop parameters.inpcrd\nquit\n"
    else:
        if polymer_parameters is None:
            _run(["parmchk2", "-i", "typed.mol2", "-f", "mol2", "-o", "missing.frcmod", "-s", "gaff2"], folder, "parmchk2.log")
        else:
            # LEaP below still checks every actual whole-chain parameter. Only
            # the redundant search for local GAFF2 additions is shared.
            shutil.copy2(polymer_parameters, folder / "missing.frcmod")
        if re.search(r"ATTN|need revision", (folder / "missing.frcmod").read_text(), re.I):
            raise RuntimeError(f"Unresolved force-field parameter in {folder / 'missing.frcmod'}")
        commands = "source leaprc.gaff2\nloadAmberParams missing.frcmod\nmol = loadMol2 typed.mol2\ncheck mol\nsaveAmberParm mol parameters.prmtop parameters.inpcrd\nquit\n"
    (folder / "tleap.in").write_text(commands)
    for name in ("parameters.prmtop", "parameters.inpcrd"):
        (folder / name).unlink(missing_ok=True)
    _run(["tleap", "-f", "tleap.in"], folder, "tleap.log")
    log = (folder / "tleap.log").read_text()
    if re.search(r"Errors\s*=\s*[1-9]|FATAL|Could not find", log, re.I):
        raise RuntimeError(f"LEaP parameter assignment failed; see {folder / 'tleap.log'}")
    for name in ("parameters.prmtop", "parameters.inpcrd"):
        if not (folder / name).is_file() or (folder / name).stat().st_size == 0:
            raise RuntimeError(f"LEaP did not produce {folder / name}")


def parameterize(case_dir):
    """Run in the Amber module environment. Writes one Amber topology per species."""
    case_dir = Path(case_dir).resolve()
    for command in ("antechamber", "parmchk2", "tleap"):
        if shutil.which(command) is None:
            raise RuntimeError(f"{command} is missing; load the configured Amber module")
    manifest = json.loads((case_dir / "manifest.json").read_text())
    library = _charge_library(case_dir)
    polymer_parameters = _polymer_parameters(case_dir)
    results = []
    for species in manifest["species"]:
        folder = case_dir / "molecules" / species["id"]
        if species.get("formal_charge", 0) != 0:
            raise ValueError("This initial PoC supports neutral chemical species only")
        kind = species["kind"].lower()
        if kind in ("water", "tip3p"):
            _topology(folder, water=True)
        elif kind in ("polymer", "plga"):
            mapping = json.loads((folder / "atom_map.json").read_text())
            charges, report = _polymer_charges(mapping, library)
            (folder / "charge_transfer.json").write_text(json.dumps(report, indent=2)+"\n")
            _antechamber(folder, charges)
            _topology(folder, polymer_parameters=polymer_parameters)
        else:
            _antechamber(folder)
            _topology(folder)
        result = {"id": species["id"], "prmtop": str((folder / "parameters.prmtop").relative_to(case_dir)),
                  "inpcrd": str((folder / "parameters.inpcrd").relative_to(case_dir))}
        results.append(result)
    report = {"method": METHOD, "status": "provisional PoC force field", "species": results}
    (case_dir / "parameterization.json").write_text(json.dumps(report, indent=2)+"\n")
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--fragments", action="store_true", help="Generate reduced polymers in the RDKit preparation environment")
    args = parser.parse_args()
    prepare_fragments(args.case_dir) if args.fragments else parameterize(args.case_dir)
