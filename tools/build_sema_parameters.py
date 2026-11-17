#!/usr/bin/env python3
"""Build an explicitly PROVISIONAL ff14SB semaglutide pilot parameter set.

Canonical residues use ff14SB. Noncanonical residues use ACE-X-NME AM1-BCC
charges, AMBER atom types, a CX alpha carbon, and parmchk2 protein-parameter
analogies. This is NOT the CHARMM36m/CGenFF model of Colston et al. (2025),
and a stable trajectory does not validate these custom parameters.

Requires Python >=3.9, RDKit, NumPy, ParmEd, AmberTools >=24.
"""
import argparse
import collections
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import parmed as pmd


def run(cmd, cwd):
    with (cwd / "commands.log").open("a") as out:
        out.write("\n$ " + " ".join(map(str, cmd)) + "\n")
        out.flush()
        child_env = os.environ.copy()
        # These small semiempirical matrices run faster without BLAS/OpenMP
        # oversubscription; production MD has a separate thread allocation.
        child_env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
        result = subprocess.run(list(map(str, cmd)), cwd=cwd, env=child_env, stdout=out, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("Command failed: %s; inspect %s" % (cmd, cwd / "commands.log"))


def read_reference(path):
    residues = collections.OrderedDict()
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            num = int(line[22:26])
            residues.setdefault(num, {"name": line[17:20].strip(), "atoms": {}})
            residues[num]["atoms"][line[12:16].strip()] = (line[76:78].strip(), np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    assert len(residues) == 31 and residues[2]["name"] == "AIB" and residues[20]["name"] == "SLY"
    assert sum(len(v["atoms"]) for v in residues.values()) == 291
    return residues


def capped_residue(residues, num, allow_repair=True):
    record = residues[num]
    name = record["name"]
    atoms = dict(record["atoms"])
    prev = residues[num - 1]["atoms"]
    nxt = residues[num + 1]["atoms"]
    # Neighbor backbone coordinates provide chemically faithful ACE/NME caps.
    atoms.update({"C0": prev["C"], "O0": prev["O"], "CM0": prev["CA"], "N0": nxt["N"], "CM1": nxt["CA"]})
    bonds = [("N", "CA", 1), ("CA", "C", 1), ("C", "O", 2), ("C0", "O0", 2), ("C0", "CM0", 1), ("C0", "N", 1), ("C", "N0", 1), ("N0", "CM1", 1)]
    if name == "AIB":
        bonds += [("CA", "CB1", 1), ("CA", "CB2", 1)]
    else:
        chain = "CA C4 C5 C6 C7 N2 C8 C9 O1 C10 C11 O2 C12 C13 N3 C14 C15 O3 C16 C17 O4 C18 C19 N4 C20 C21 C22 C23 N5 C24 C25 C26 C27 C28 C29 C30 C31 C32 C33 C34 C35 C36 C37 C38 C39 C40 C41".split()
        bonds += [(a, b, 1) for a, b in zip(chain, chain[1:])]
        bonds += [("C8", "O12", 2), ("C14", "O11", 2), ("C20", "O10", 2), ("C24", "O7", 2), ("C23", "C42", 1), ("C42", "O8", 1), ("C42", "O9", 2), ("C41", "O5", 1), ("C41", "O6", 2)]
    names = list(atoms)
    index = {a: i for i, a in enumerate(names)}
    rw = Chem.RWMol()
    for atom_name in names:
        atom = Chem.Atom(atoms[atom_name][0])
        atom.SetProp("atom_name", atom_name)
        atom.SetProp("part", "cap" if atom_name in {"C0", "O0", "CM0", "N0", "CM1"} else "residue")
        if name == "SLY" and atom_name in {"O5", "O8"}:
            atom.SetFormalCharge(-1)
        rw.AddAtom(atom)
    for a, b, order in bonds:
        rw.AddBond(index[a], index[b], Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
        distance = np.linalg.norm(atoms[a][1] - atoms[b][1])
        if not 0.95 < distance < 1.9:
            raise ValueError("Unphysical reference bond %s-%s: %.3f A" % (a, b, distance))
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conformer = Chem.Conformer(len(names))
    for i, atom_name in enumerate(names):
        conformer.SetAtomPosition(i, atoms[atom_name][1].tolist())
    mol.AddConformer(conformer)
    mol = Chem.AddHs(mol, addCoords=True)
    hcount = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            parent = atom.GetNeighbors()[0]
            hcount += 1
            hn = "H" if parent.GetProp("atom_name") == "N" else "HA" if parent.GetProp("atom_name") == "CA" else "H%03d" % hcount
            atom.SetProp("atom_name", hn)
            atom.SetProp("part", parent.GetProp("part"))
    # Only relax generated hydrogens; retain the published heavy-atom geometry.
    ff = AllChem.UFFGetMoleculeForceField(mol)
    for i in range(len(names)):
        ff.AddFixedPoint(i)
    ff.Minimize(maxIts=1000)
    Chem.AssignStereochemistryFrom3D(mol)
    if name == "SLY":
        for target in ("CA", "C23"):
            a = mol.GetAtomWithIdx(index[target])
            if not a.HasProp("_CIPCode") or a.GetProp("_CIPCode") != "S":
                if target == "C23" and allow_repair and a.HasProp("_CIPCode") and a.GetProp("_CIPCode") == "R":
                    # Official semaglutide contains L-gamma-Glu (S). The source
                    # supplementary geometry has this center inverted. Reflect
                    # its carboxyl branch through C22-C23-N5, preserving bond
                    # lengths/angles, then regenerate H and verify S explicitly.
                    center = record["atoms"]["C23"][1]
                    normal = np.cross(record["atoms"]["C22"][1]-center, record["atoms"]["N5"][1]-center)
                    normal /= np.linalg.norm(normal)
                    for branch_atom in ("C42", "O8", "O9"):
                        element, pos = record["atoms"][branch_atom]
                        record["atoms"][branch_atom] = (element, pos - 2*np.dot(pos-center, normal)*normal)
                    record["gamma_glu_chirality_repaired"] = "R-to-S carboxyl-branch reflection"
                    return capped_residue(residues, num, allow_repair=False)
                raise ValueError("Expected L/S stereochemistry at SLY %s" % target)
    return mol


def write_mol2(mol, name, path):
    xyz = mol.GetConformer()
    rows = ["@<TRIPOS>MOLECULE", name, "%d %d 1 0 0" % (mol.GetNumAtoms(), mol.GetNumBonds()), "SMALL", "USER_CHARGES", "", "@<TRIPOS>ATOM"]
    for a in mol.GetAtoms():
        pos = xyz.GetAtomPosition(a.GetIdx())
        if a.GetSymbol() == "H":
            typ = "H"
        elif a.GetSymbol() == "C":
            typ = "C.2" if any(b.GetBondTypeAsDouble() == 2 for b in a.GetBonds()) else "C.3"
        elif a.GetSymbol() == "N":
            typ = "N.am"
        elif a.GetSymbol() == "O":
            parent = a.GetNeighbors()[0]
            carboxyl = sum(n.GetSymbol() == "O" for n in parent.GetNeighbors()) == 2
            typ = "O.co2" if carboxyl else "O.2" if a.GetDegree() == 1 else "O.3"
        else:
            raise ValueError("Unexpected element")
        rows.append("%7d %-4s %11.6f %11.6f %11.6f %-5s 1 %-3s %10.6f" % (a.GetIdx() + 1, a.GetProp("atom_name"), pos.x, pos.y, pos.z, typ, name, float(a.GetFormalCharge())))
    rows.append("@<TRIPOS>BOND")
    for i, b in enumerate(mol.GetBonds(), 1):
        rows.append("%6d %6d %6d %d" % (i, b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1, int(b.GetBondTypeAsDouble())))
    rows += ["@<TRIPOS>SUBSTRUCTURE", "1 %s 1" % name]
    path.write_text("\n".join(rows) + "\n")


def preopt_mmff(mol):
    """Gas-phase MMFF94s preoptimization before the unchanged AM1-BCC stage."""
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise RuntimeError("MMFF94s preoptimization parameters are incomplete")
    if AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=10000):
        raise RuntimeError("MMFF94s preoptimization did not converge")
    Chem.AssignStereochemistryFrom3D(mol)
    for a in mol.GetAtoms():
        if a.HasProp("_CIPCode") and a.GetProp("atom_name") in ("CA", "C23") and a.GetProp("_CIPCode") != "S":
            raise RuntimeError("MMFF94s preoptimization changed stereochemistry")
    return mol


def build_fragment(residues, num, output, use_mmff=False, reuse=False, single_point=False):
    name = residues[num]["name"]
    charge = 0 if name == "AIB" else -2
    work = output / name
    work.mkdir(exist_ok=True)
    mol = capped_residue(residues, num)
    if use_mmff:
        mol = preopt_mmff(mol)
    prior_input = (work / "capped.mol2").read_bytes() if (work / "capped.mol2").exists() else None
    write_mol2(mol, name, work / "capped.mol2")
    convergence_marker = "Calculation Completed" if single_point else "geometry converged"
    if reuse:
        sqm_input = (work / "sqm.in").read_text() if (work / "sqm.in").exists() else ""
        stored_single_point = "maxcyc=0" in sqm_input.replace(" ", "")
        if stored_single_point != single_point:
            raise RuntimeError("Reused charge geometry method differs from requested method")
        sqm_path = work / "sqm.out"
        if prior_input != (work / "capped.mol2").read_bytes() or not (work / "charged.ac").exists() or not sqm_path.exists() or convergence_marker not in sqm_path.read_text() or "Calculation Completed" not in sqm_path.read_text():
            raise RuntimeError("Cannot reuse charges: identical input and completed/converged SQM output required for " + name)
    else:
        command = ["antechamber", "-i", "capped.mol2", "-fi", "mol2", "-o", "charged.ac", "-fo", "ac", "-at", "amber", "-c", "bcc", "-nc", str(charge), "-rn", name, "-an", "n", "-pf", "y", "-dr", "no"]
        if single_point:
            command += ["-ek", "qm_theory='AM1', maxcyc=0, scfconv=1.d-10, ndiis_attempts=700,"]
        run(command, work)
    sqm_text = (work / "sqm.out").read_text()
    if convergence_marker not in sqm_text or "Calculation Completed" not in sqm_text or "SCF CONVERGENCE FAILURE" in sqm_text.upper():
        raise RuntimeError("SQM did not fully converge for " + name)
    lines = []
    for line in (work / "charged.ac").read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            pieces = line.split()
            # AC format has the AMBER atom type in its final field.
            pos = line.rfind(pieces[-1])
            line = line[:pos] + "CX"
        lines.append(line)
    (work / "typed.ac").write_text("\n".join(lines) + "\n")
    cap_names = [a.GetProp("atom_name") for a in mol.GetAtoms() if a.GetProp("part") == "cap"]
    mc = ["HEAD_NAME N", "TAIL_NAME C", "MAIN_CHAIN CA"] + ["OMIT_NAME " + n for n in cap_names] + ["PRE_HEAD_TYPE C", "POST_TAIL_TYPE N", "CHARGE %.1f" % charge]
    (work / "residue.mc").write_text("\n".join(mc) + "\n")
    run(["prepgen", "-i", "typed.ac", "-o", name + ".prepi", "-m", "residue.mc", "-rn", name], work)
    # CX is the ff14SB backbone sp3 carbon type. For genuinely absent
    # AIB angles (CT-CX-CT), use the existing CT carbon-angle analogue;
    # retain every exact ff14SB backbone torsion when it is present.
    database = Path(os.environ["AMBERHOME"]) / "dat/antechamber/PARMCHK.DAT"
    atc = database.read_text()
    atc = atc.replace("PARM\tCX\t0\t\t0\t12.01\t0\t6\n", "PARM\tCX\t0\t\t0\t12.01\t0\t6\nEQUA\tCT\n")
    (work / "PARMCHK_CX_CT.DAT").write_text(atc)
    run(["parmchk2", "-i", name + ".prepi", "-f", "prepi", "-o", name + ".frcmod", "-s", "4", "-frc", "ff14SB", "-c", "PARMCHK_CX_CT.DAT"], work)
    frc = (work / (name + ".frcmod")).read_text()
    if "ATTN" in frc or "need revision" in frc.lower():
        raise RuntimeError("Unresolved force field terms in " + str(work / (name + ".frcmod")))
    original_charges = {}
    for line in (work / "typed.ac").read_text().splitlines():
        if line.startswith("ATOM"):
            fields = line.split()
            original_charges[line[13:17].strip()] = float(fields[-2])
    residue_charges = {}
    for line in (work / (name + ".prepi")).read_text().splitlines():
        fields = line.split()
        if len(fields) == 11 and fields[0].isdigit() and fields[2] != "DU":
            residue_charges[fields[1]] = float(fields[-1])
    equivalence_groups = [["CB1", "CB2"], ["H002", "H003", "H004", "H005", "H006", "H007"]] if name == "AIB" else [["O5", "O6"], ["O8", "O9"]]
    symmetry = [{"atoms": group, "final_charge_spread_e": max(residue_charges[n] for n in group)-min(residue_charges[n] for n in group)} for group in equivalence_groups]
    if any(check["final_charge_spread_e"] > 1e-6 for check in symmetry):
        raise RuntimeError("Expected chemically equivalent charge groups differ in " + name)
    sqm_audit = {"method": "AM1", "scf_tolerance": 1e-10, "completed": "Calculation Completed" in sqm_text, "geometry_minimization_requested": not single_point, "AM1_geometry_converged": None if single_point else "geometry converged" in sqm_text, "MMFF94s_preoptimization_converged": bool(use_mmff), "charge_variant": "AM1-BCC single point on MMFF94s geometry" if single_point else "AM1-BCC with AM1 geometry optimization"}
    return {"sqm_audit": sqm_audit, "symmetry_handling": "Antechamber AM1-BCC equivalent-atom charge handling retained; post-prepgen equivalence explicitly checked", "symmetry_checks": symmetry, "residue": name, "formal_charge": charge, "capped_formula": rdMolDescriptors.CalcMolFormula(mol), "cap_atoms_removed": cap_names, "supplemental_parameters": frc, "charge_sum_before_cap_correction": sum(original_charges[n] for n in residue_charges), "charge_sum_after_cap_correction": sum(residue_charges.values()), "atom_charge_changes": [{"atom": n, "capped_AM1BCC_e": original_charges[n], "final_polymer_residue_e": q, "correction_e": q-original_charges[n]} for n, q in residue_charges.items()], "CX_CT_analogy": "Only absent parameters may use existing CT analogue for CX; exact ff14SB entries retain priority", "provisional": True}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--charge-geometry", choices=("mmff-sp", "mmff-am1opt", "reference-am1opt"), default="mmff-sp", help="Default mmff-sp: converged MMFF94s geometry and single-point AM1-BCC; alternatives perform full AM1 geometry optimization")
    ap.add_argument("--preopt-mmff", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--reuse-completed-charges", action="store_true", help="Reuse only identical capped inputs with completed/converged SQM outputs")
    ap.add_argument("--single-point-am1", action="store_true", help="Documented Amber AM1-BCC variant: SCF on converged MMFF94s geometry, no gas-phase AM1 geometry minimization")
    args = ap.parse_args()
    args.preopt_mmff = args.charge_geometry in ("mmff-sp", "mmff-am1opt") or args.preopt_mmff
    args.single_point_am1 = args.charge_geometry == "mmff-sp" or args.single_point_am1
    if args.single_point_am1 and not args.preopt_mmff:
        ap.error("--single-point-am1 requires --preopt-mmff")
    for tool in ("antechamber", "prepgen", "parmchk2", "tleap"):
        if not shutil.which(tool):
            raise RuntimeError("Required AmberTools executable not found: " + tool)
    output = Path(args.out).resolve()
    output.mkdir(parents=True, exist_ok=True)
    reference = Path(args.reference).resolve()
    residues = read_reference(reference)
    audit = {"status": "PROVISIONAL_PARAMETERS_NOT_EXPERIMENTALLY_VALIDATED", "protein_force_field": "ff14SB", "custom_residues": [], "charge_method": "ACE-X-NME AM1-BCC, prepgen charge correction after cap removal", "limitations": ["Not the Nature Communications 2025 CHARMM36m/CGenFF parameter set", "Custom AM1-BCC charges are not multiconformer RESP charges", "AIB alpha-carbon and torsions use ff14SB AMBER-type analogies without a bespoke QM torsion fit", "Short MD stability is not force-field validation"]}
    audit["charge_geometry_option"] = args.charge_geometry
    audit["geometry_for_charge_calculation"] = "MMFF94s preoptimization then default SQM AM1 optimization (grms_tol=0.0005, scfconv=1e-10)" if args.preopt_mmff else "Reference heavy atoms/H relaxation then default SQM AM1 optimization (grms_tol=0.0005, scfconv=1e-10)"
    if args.single_point_am1:
        audit["charge_method"] = "Single-point AM1-BCC on converged MMFF94s ACE-X-NME geometries; fully converged AM1 SCF at 1e-10, maxcyc=0; prepgen cap-removal charge correction"
        audit["geometry_for_charge_calculation"] = "Converged MMFF94s geometry, then single-point AM1 SCF (maxcyc=0, scfconv=1e-10) and BCC; NOT AM1 geometry-optimized charges"
        audit["charge_method_source"] = "Amber25 Reference Manual sec16.2 tip7 (maxcyc=0 for trusted geometries; caveat about gas-phase anion conformations)"
    for num in (2, 20):
        audit["custom_residues"].append(build_fragment(residues, num, output, use_mmff=args.preopt_mmff, reuse=args.reuse_completed_charges, single_point=args.single_point_am1))
    # Convert CHARMM heavy-atom conventions to Amber; no heavy atoms are dropped.
    pdb_lines = []
    for line in reference.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = line[17:20].strip()
        atom = line[12:16].strip()
        if residue == "HIS":
            residue = "HIE"
        if residue == "ILE" and atom == "CD":
            atom = "CD1"
        if int(line[22:26]) == 31:
            atom = {"OT1": "O", "OT2": "OXT"}.get(atom, atom)
        pos = residues[int(line[22:26])]["atoms"][line[12:16].strip()][1]
        pdb_lines.append(line[:12] + ("%4s" % atom) + line[16:17] + ("%3s" % residue) + line[20:30] + ("%8.3f%8.3f%8.3f" % tuple(pos)) + line[54:])
    (output / "solute_heavy.pdb").write_text("\n".join(pdb_lines) + "\nTER\nEND\n")
    leap = ["source leaprc.protein.ff14SB"]
    for name in ("AIB", "SLY"):
        leap += ["loadamberprep %s/%s.prepi" % (name, name), "loadamberparams %s/%s.frcmod" % (name, name)]
    leap += ["mol = loadpdb solute_heavy.pdb", "check mol", "charge mol", "saveamberparm mol solute.prmtop solute.inpcrd", "savepdb mol solute.pdb", "quit"]
    (output / "build.leap").write_text("\n".join(leap) + "\n")
    run(["tleap", "-f", "build.leap"], output)
    if not (output / "solute.prmtop").exists():
        raise RuntimeError("LEaP did not write topology; inspect commands.log")
    structure = pmd.load_file(str(output / "solute.prmtop"), xyz=str(output / "solute.inpcrd"))
    formula = collections.Counter(a.element_name for a in structure.atoms)
    expected = {"C": 187, "H": 287, "N": 45, "O": 59}
    if dict(formula) != expected or abs(sum(a.charge for a in structure.atoms) + 4) > 1e-4:
        raise ValueError("Wrong composition/charge: %s / %.6f" % (formula, sum(a.charge for a in structure.atoms)))
    if len(structure.residues) != 31:
        raise ValueError("Expected 31 residues")
    audit["gamma_glu_geometry_correction"] = residues[20].get("gamma_glu_chirality_repaired", "none; S confirmed")
    audit["gamma_glu_identity_source"] = "https://precision.fda.gov/uniisearch/srs/unii/53axn4nnhx"
    audit["formula_at_selected_protonation"] = dict(formula)
    audit["net_charge_e"] = sum(a.charge for a in structure.atoms)
    audit["atom_count"] = len(structure.atoms)
    audit["residue_charges"] = [{"residue": r.idx+1, "name": r.name, "charge": sum(a.charge for a in r.atoms)} for r in structure.residues]
    structure.save(str(output / "solute.top"), overwrite=True)
    structure.save(str(output / "solute.gro"), overwrite=True)
    (output / "parameter_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"out": str(output), "charge": audit["net_charge_e"], "atoms": len(structure.atoms), "status": audit["status"]}, indent=2))


if __name__ == "__main__":
    main()
