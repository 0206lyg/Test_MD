#!/usr/bin/env python3
"""Reproducible PLGA graphs and unequilibrated, self-avoiding starting coils.

Coordinates are in Angstrom. This builder supplies geometry, not a force field
or an equilibrated polymer phase. Long chains use local internal coordinates;
whole-chain distance-geometry embedding is deliberately unnecessary.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np


def build_graph(sequence: Sequence[str], stereo: Sequence[str | None] | None = None,
                end_groups: str = "hydroxyl_carboxyl", explicit_hydrogens: bool = True):
    """Return (RDKit Mol, atom/repeat metadata), without a conformer.

    LA means a lactyl residue, GA a glycolyl residue. S=L and R=D at LA.
    The default chain is HO-CH(R)-C(=O)-...-O-CH(R)-C(=O)-OH.
    ``hydroxyl_methyl_ester`` changes only the terminal acid to its methyl ester.
    """
    from rdkit import Chem
    sequence = list(sequence)
    if not sequence or any(x not in {"LA", "GA"} for x in sequence):
        raise ValueError("sequence must contain one or more LA/GA residues")
    if end_groups not in {"hydroxyl_carboxyl", "hydroxyl_methyl_ester"}:
        raise ValueError("unsupported end_groups")
    stereo = list(stereo) if stereo is not None else ["L" if x == "LA" else None for x in sequence]
    if len(stereo) != len(sequence):
        raise ValueError("stereo and sequence lengths differ")
    rw = Chem.RWMol()
    atom_info, repeats = [], []

    def add(element, role, i):
        idx = rw.AddAtom(Chem.Atom(element))
        atom_info.append({"index": idx, "repeat": i, "monomer": sequence[i],
                          "role": role, "parent_role": None})
        return idx

    previous_c = None
    for i, kind in enumerate(sequence):
        if kind == "LA" and stereo[i] not in {"L", "D"}:
            raise ValueError("LA stereochemistry must be L or D")
        if kind == "GA" and stereo[i] is not None:
            raise ValueError("GA is achiral; stereo must be None")
        ids = {role: add(element, role, i) for role, element in
               [("O", "O"), ("CA", "C"), ("C", "C"), ("OC", "O")]}
        rw.AddBond(ids["O"], ids["CA"], Chem.BondType.SINGLE)
        rw.AddBond(ids["CA"], ids["C"], Chem.BondType.SINGLE)
        rw.AddBond(ids["C"], ids["OC"], Chem.BondType.DOUBLE)
        if previous_c is not None:
            rw.AddBond(previous_c, ids["O"], Chem.BondType.SINGLE)
        if kind == "LA":
            ids["CM"] = add("C", "CM", i)
            rw.AddBond(ids["CA"], ids["CM"], Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(ids["CA"]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        repeats.append({"repeat": i, "monomer": kind, "stereo": stereo[i], "atoms": ids})
        previous_c = ids["C"]
    terminal = add("O", "OT", len(sequence) - 1)
    rw.AddBond(previous_c, terminal, Chem.BondType.SINGLE)
    repeats[-1]["atoms"]["OT"] = terminal
    if end_groups == "hydroxyl_methyl_ester":
        cap = add("C", "CT", len(sequence) - 1)
        rw.AddBond(terminal, cap, Chem.BondType.SINGLE)
        repeats[-1]["atoms"]["CT"] = cap
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for rep in repeats:
        if rep["monomer"] == "LA":
            atom = mol.GetAtomWithIdx(rep["atoms"]["CA"])
            if atom.GetProp("_CIPCode") != ("S" if rep["stereo"] == "L" else "R"):
                atom.InvertChirality()
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if explicit_hydrogens:
        mol = Chem.AddHs(mol)
        counts = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            parent = atom.GetNeighbors()[0].GetIdx()
            info = atom_info[parent]
            counts[parent] = counts.get(parent, 0) + 1
            atom_info.append({"index": atom.GetIdx(), "repeat": info["repeat"],
                              "monomer": info["monomer"],
                              "role": f"{info['role']}:H{counts[parent]}",
                              "parent_role": info["role"], "parent_index": parent})
    return mol, {"sequence": sequence, "stereo": stereo, "end_groups": end_groups,
                 "atoms": atom_info, "repeats": repeats}


def _unit(v):
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        raise ValueError("degenerate internal coordinates")
    return v / norm


def _place(a, b, c, length, angle, dihedral):
    """Place atom d using distance(c,d), angle(b,c,d), dihedral(a,b,c,d)."""
    axis = _unit(c - b)
    normal = _unit(np.cross(b - a, axis))
    side = np.cross(normal, axis)
    theta, phi = math.radians(angle), math.radians(dihedral)
    return c + length * (-math.cos(theta) * axis + math.sin(theta) *
                         (math.cos(phi) * side + math.sin(phi) * normal))


class _Cells:
    def __init__(self, cutoff):
        self.cutoff = cutoff
        self.squared = cutoff * cutoff
        self.cells = {}
        self.offsets = list(itertools.product((-1, 0, 1), repeat=3))

    def key(self, xyz):
        return tuple(np.floor(xyz / self.cutoff).astype(int))

    def add(self, points, chain, repeat):
        records = []
        for xyz in points:
            key = self.key(xyz)
            record = (xyz, chain, repeat)
            self.cells.setdefault(key, []).append(record)
            records.append((key, record))
        return records

    def remove(self, records):
        for key, record in reversed(records):
            bucket = self.cells[key]
            for i in range(len(bucket) - 1, -1, -1):
                if bucket[i] is record:
                    bucket.pop(i)
                    break
            if not bucket:
                del self.cells[key]

    def clashes(self, points, chain, repeat):
        for xyz in points:
            base = self.key(xyz)
            for off in self.offsets:
                key = (base[0] + off[0], base[1] + off[1], base[2] + off[2])
                for other, other_chain, other_repeat in self.cells.get(key, ()):
                    # Same and adjacent repeats contain many bonded/1-4 pairs.
                    if other_chain == chain and abs(other_repeat - repeat) <= 1:
                        continue
                    delta = other - xyz
                    if np.dot(delta, delta) < self.squared:
                        return True
        return False


def _local_clash(points, previous, cutoff):
    """Check within/adjacent-residue pairs, excluding only 1-2 and 1-3 pairs."""
    # In O, CA, C, OC, (CM) order the only intrarepeat pairs separated
    # by three bonds are O--OC and (for lactyl residues) CM--OC.
    if np.linalg.norm(points[0] - points[3]) < cutoff:
        return True
    if len(points) == 5 and np.linalg.norm(points[4] - points[3]) < cutoff:
        return True
    if previous is not None:
        from_previous_c = np.array([2, 1, 0, 1, 2][:len(previous)])
        from_current_o = np.array([0, 1, 2, 3, 2][:len(points)])
        separation = from_previous_c[:, None] + 1 + from_current_o[None, :]
        distances = np.sum((previous[:, None, :] - points[None, :, :]) ** 2, axis=2)
        if np.any((separation > 2) & (distances < cutoff * cutoff)):
            return True
    return False


def _hydrogens(points, first, rng):
    """Tetrahedral starting H geometry; returns stable role-to-coordinate map."""
    oxygen, ca, carbon = points[:3]
    u, v = _unit(oxygen - ca), _unit(carbon - ca)
    parallel = -(u + v) / (3. * (1. + np.dot(u, v)))
    normal = _unit(np.cross(u, v))
    out = math.sqrt(1. - np.dot(parallel, parallel)) * normal
    result = {}
    if len(points) == 5:
        methyl = points[4]
        methyl_direction = _unit(methyl - ca)
        branch = 1 if np.dot(methyl_direction, normal) > 0 else -1
        result["CA:H1"] = ca + 1.09 * (parallel - branch * out)
        toward_ca = _unit(ca - methyl)
        axis1 = _unit(np.cross(toward_ca, normal))
        axis2 = np.cross(toward_ca, axis1)
        phase = rng.uniform(0., 2. * math.pi)
        for j in range(3):
            angle = phase + j * 2. * math.pi / 3.
            direction = -toward_ca / 3. + math.sqrt(8. / 9.) * (
                math.cos(angle) * axis1 + math.sin(angle) * axis2)
            result[f"CM:H{j+1}"] = methyl + 1.09 * direction
    else:
        result["CA:H1"] = ca + 1.09 * (parallel + out)
        result["CA:H2"] = ca + 1.09 * (parallel - out)
    if first:
        toward_ca = _unit(ca - oxygen)
        perpendicular = _unit(np.cross(toward_ca, carbon - ca))
        angle = rng.uniform(0., 2. * math.pi)
        radial = math.cos(angle) * perpendicular + math.sin(angle) * np.cross(toward_ca, perpendicular)
        result["O:H1"] = oxygen + .97 * (math.cos(math.radians(108.)) * toward_ca +
                                         math.sin(math.radians(108.)) * radial)
    return result


def _adjacent_h_clash(points, roles, previous):
    if previous is None:
        return False
    dc = {"O": 2, "CA": 1, "C": 0, "OC": 1, "CM": 2}
    do = {"O": 0, "CA": 1, "C": 2, "OC": 3, "CM": 2, "OT": 3}
    a = np.array([dc[r.split(":")[0]] + (":" in r) for r in previous["all_roles"]])
    b = np.array([do[r.split(":")[0]] + (":" in r) for r in roles])
    separation = a[:, None] + 1 + b[None, :]
    distances = np.sum((previous["all_points"][:, None, :] - points[None, :, :]) ** 2, axis=2)
    return bool(np.any((separation > 2) & (distances < 1.2 ** 2)))


def _same_repeat_clash(points, roles):
    names = ("O", "CA", "C", "OC", "CM", "OT")
    graph_distance = np.array([[0,1,2,3,2,3], [1,0,1,2,1,2],
        [2,1,0,1,2,1], [3,2,1,0,3,2], [2,1,2,3,0,3], [3,2,1,2,3,0]])
    parents = np.array([names.index(r.split(":")[0]) for r in roles])
    h = np.array([":" in r for r in roles], dtype=int)
    separation = graph_distance[parents[:, None], parents[None, :]] + h[:, None] + h[None, :]
    distances = np.sum((points[:, None, :] - points[None, :, :]) ** 2, axis=2)
    return bool(np.any((separation > 2) & (distances < 1.2 ** 2)))


def _grow_chain_iter(sequence, stereo, rng, lower, upper, cells, chain,
                max_proposals=500000, trials_per_growth=96, all_cells=None):
    if all_cells is None:
        all_cells = _Cells(1.2)
    records, proposal_count, backtracks = [], 0, 0
    next_restart = 25000
    while len(records) < len(sequence):
        if proposal_count >= next_restart:
            for rec in reversed(records):
                cells.remove(rec["cells"])
                all_cells.remove(rec["all_cells"])
            records.clear()
            next_restart += 25000
        i = len(records)
        start_candidates = None
        if i == 0:
            start_candidates = rng.uniform(lower + 3., upper - 3., size=(128, 3))
            occupied = np.asarray([entry[0] for bucket in cells.cells.values() for entry in bucket])
            if len(occupied):
                scores = []
                for point in start_candidates:
                    nearest = float(np.min(np.sum((occupied - point) ** 2, axis=1)))
                    boundary = float(np.min(np.minimum(point - lower, upper - point))) ** 2
                    scores.append(min(nearest, boundary))
                start_candidates = start_candidates[np.argsort(scores)[-16:]]
        accepted = False
        for _ in range(trials_per_growth):
            proposal_count += 1
            if proposal_count > max_proposals:
                for rec in reversed(records):
                    cells.remove(rec["cells"])
                    all_cells.remove(rec["all_cells"])
                raise RuntimeError(
                    f"Chain {chain + 1} could not be packed after {max_proposals} proposals. "
                    "No overlapping structure was written. Increase initial slab thickness "
                    "explicitly, then equilibrate/contract it with the chosen force field.")
            torsion1 = float(rng.uniform(-180., 180.))
            torsion2 = float(rng.uniform(-180., 180.))
            if i:
                prev_ca, prev_c, oxygen = records[-1]["tail"]
                ca = _place(prev_ca, prev_c, oxygen, 1.45, 116., 180.)
                carbon = _place(prev_c, oxygen, ca, 1.52, 110., torsion1)
            else:
                oxygen = start_candidates[int(rng.integers(len(start_candidates)))]
                direction = _unit(rng.normal(size=3))
                ca = oxygen + 1.45 * direction
                virtual = oxygen + _unit(np.cross(direction, rng.normal(size=3)))
                carbon = _place(virtual, oxygen, ca, 1.52, 110., torsion1)
            next_o = _place(oxygen, ca, carbon, 1.36, 113., torsion2)
            carbonyl = carbon - 1.22 * _unit(_unit(ca - carbon) + _unit(next_o - carbon))
            points = [oxygen, ca, carbon, carbonyl]
            if sequence[i] == "LA":
                u, v = _unit(oxygen - ca), _unit(carbon - ca)
                parallel = -(u + v) / (3. * (1. + np.dot(u, v)))
                out = math.sqrt(1. - np.dot(parallel, parallel)) * _unit(np.cross(u, v))
                # With this neighbor order, L is the positive branch (verified
                # below against RDKit's coordinate-derived CIP assignment).
                methyl = ca + 1.53 * (parallel + (1 if stereo[i] == "L" else -1) * out)
                points.append(methyl)
            candidate = np.asarray(points)
            hs = _hydrogens(candidate, i == 0, rng)
            roles = ["O", "CA", "C", "OC"] + (["CM"] if len(points) == 5 else [])
            all_points = list(points) + list(hs.values())
            all_roles = roles + list(hs)
            if i == len(sequence) - 1:
                terminal_h = _place(ca, carbon, next_o, .97, 108., 180.)
                hs["OT:H1"] = terminal_h
                all_points += [next_o, terminal_h]
                all_roles += ["OT", "OT:H1"]
            all_points = np.asarray(all_points)
            if _same_repeat_clash(all_points, all_roles):
                continue
            bounded = np.vstack((all_points, next_o))
            if np.any(bounded < lower) or np.any(bounded > upper):
                continue
            if _local_clash(candidate, records[-1]["points"] if i else None, cells.cutoff):
                continue
            if cells.clashes(candidate, chain, i) or cells.clashes([next_o], chain, i + 1):
                continue
            if all_cells.clashes(all_points, chain, i):
                continue
            if _adjacent_h_clash(all_points, all_roles, records[-1] if i else None):
                continue
            # Bias the *initial packing* toward persistent walks, reflecting
            # the heading near a boundary. This is not Boltzmann sampling;
            # subsequent force-field minimization/equilibration is required.
            if i >= 3 and _ < 64:
                direction = _unit(records[-1]["tail"][0] - records[max(0, i - 6)]["tail"][0])
                current = records[-1]["tail"][2]
                for axis in range(3):
                    if current[axis] < lower[axis] + 8. and direction[axis] < 0:
                        direction[axis] *= -1.
                    if current[axis] > upper[axis] - 8. and direction[axis] > 0:
                        direction[axis] *= -1.
                alignment = np.dot(_unit(next_o - current), direction)
                if rng.random() > math.exp(2. * (alignment - 1.)):
                    continue
            records.append({"points": candidate, "tail": (ca, carbon, next_o),
                            "cells": cells.add(np.vstack((candidate, next_o))
                                if i == len(sequence) - 1 else candidate, chain, i),
                            "all_cells": all_cells.add(all_points, chain, i),
                            "all_points": all_points, "all_roles": all_roles, "hydrogens": hs})
            accepted = True
            yield len(records)
            break
        if not accepted:
            backtracks += 1
            n_remove = min(len(records), int(rng.integers(3, 16)))
            if backtracks % 20 == 0:
                n_remove = min(len(records), 50)
            for _ in range(n_remove):
                record = records.pop()
                cells.remove(record["cells"])
                all_cells.remove(record["all_cells"])
    xyz = np.vstack([rec["points"] for rec in records] + [[records[-1]["tail"][2]]])
    hydrogen_coords = {(i, role): point for i, rec in enumerate(records)
                       for role, point in rec["hydrogens"].items()}
    return xyz, {"proposals": proposal_count, "backtracks": backtracks}, hydrogen_coords


def _grow_chain(*args, **kwargs):
    """Single-chain convenience wrapper around the incremental growth engine."""
    iterator = _grow_chain_iter(*args, **kwargs)
    while True:
        try:
            next(iterator)
        except StopIteration as done:
            return done.value


def _coordinate_chirality(mol, metadata, xyz):
    """Validate coordinate stereochemistry, correcting the initial branch convention."""
    from rdkit import Chem
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(xyz):
        conf.SetAtomPosition(i, tuple(map(float, point)))
    mol.AddConformer(conf, assignId=True)
    probe = Chem.Mol(mol)
    Chem.AssignAtomChiralTagsFromStructure(probe, replaceExistingTags=True)
    Chem.AssignStereochemistry(probe, cleanIt=True, force=True)
    corrections = 0
    for rep in metadata["repeats"]:
        if rep["monomer"] != "LA":
            continue
        ids = rep["atoms"]
        expected = "S" if rep["stereo"] == "L" else "R"
        atom = probe.GetAtomWithIdx(ids["CA"])
        if not atom.HasProp("_CIPCode"):
            raise RuntimeError("Could not infer alpha-carbon chirality from coordinates")
        if atom.GetProp("_CIPCode") != expected:
            ca = xyz[ids["CA"]]
            normal = _unit(np.cross(xyz[ids["O"]] - ca, xyz[ids["C"]] - ca))
            delta = xyz[ids["CM"]] - ca
            xyz[ids["CM"]] -= 2 * np.dot(delta, normal) * normal
            corrections += 1
    # A reflected methyl would invalidate collision checks. The growth formula
    # must be corrected in code if its convention differs from RDKit.
    if corrections:
        raise RuntimeError(f"Internal chirality convention error ({corrections} centers); no output written")
    return mol


def _soft_push_off(mol, lower, upper, heavy_min=1.4, all_min=1.2, max_cycles=30):
    """Geometric preparation, not physical energy minimization.

    Finite overlap penalties move initially overlapping independent chains.
    Harmonic constraints preserve the initial 1-2/1-3 distances and the signed
    volume of each tetrahedral center. A GAFF energy minimization is still
    mandatory before any MD. SciPy is imported only for this preparation step.
    """
    from scipy.spatial import cKDTree
    from scipy.optimize import minimize
    from itertools import combinations
    x0 = mol.GetConformer().GetPositions().copy()
    n = len(x0)
    heavy = np.array([a.GetAtomicNum() != 1 for a in mol.GetAtoms()])
    pairs = set()
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        neighbors = [a.GetIdx() for a in atom.GetNeighbors()]
        pairs.update(tuple(sorted((i, j))) for j in neighbors)
        pairs.update(tuple(sorted(p)) for p in combinations(neighbors, 2))
    geometry = np.array(sorted(pairs), dtype=int)
    a, b = geometry.T
    reference = np.linalg.norm(x0[a] - x0[b], axis=1)
    excluded = np.array(sorted(i * n + j for i, j in pairs))
    bounds = np.tile(np.c_[lower, upper], (n, 1))
    chiral = []
    for atom in mol.GetAtoms():
        if str(atom.GetChiralTag()) != "CHI_UNSPECIFIED":
            neighbors = [v.GetIdx() for v in atom.GetNeighbors() if v.GetAtomicNum() != 1]
            if len(neighbors) == 3:
                chiral.append([atom.GetIdx()] + neighbors)
    chiral = np.asarray(chiral, dtype=int).reshape(-1, 4)
    c, u, v, w = chiral.T

    def volumes(x):
        uu, vv, ww = x[u] - x[c], x[v] - x[c], x[w] - x[c]
        return np.einsum("ij,ij->i", uu, np.cross(vv, ww)), uu, vv, ww

    reference_volume = volumes(x0)[0]

    def contacts(x, radius):
        pp = cKDTree(x).query_pairs(radius, output_type="ndarray")
        return pp[~np.isin(pp[:, 0] * n + pp[:, 1], excluded, assume_unique=False)]

    def objective(flat, q, desired):
        x = flat.reshape(n, 3)
        gradient = np.zeros_like(x)
        delta = x[a] - x[b]
        r = np.maximum(np.linalg.norm(delta, axis=1), 1e-10)
        difference = r - reference
        energy = 100. * np.dot(difference, difference)
        dg = 200. * difference[:, None] * delta / r[:, None]
        np.add.at(gradient, a, dg)
        np.add.at(gradient, b, -dg)
        p1, p2 = q.T
        delta = x[p1] - x[p2]
        r = np.maximum(np.linalg.norm(delta, axis=1), 1e-10)
        overlap = np.minimum(r - desired, 0.)
        energy += 10. * np.dot(overlap, overlap)
        dg = 20. * overlap[:, None] * delta / r[:, None]
        np.add.at(gradient, p1, dg)
        np.add.at(gradient, p2, -dg)
        vv, uu, uv, uw = volumes(x)
        dv = vv - reference_volume
        energy += 10. * np.dot(dv, dv)
        gu = 20. * dv[:, None] * np.cross(uv, uw)
        gv = 20. * dv[:, None] * np.cross(uw, uu)
        gw = 20. * dv[:, None] * np.cross(uu, uv)
        np.add.at(gradient, u, gu)
        np.add.at(gradient, v, gv)
        np.add.at(gradient, w, gw)
        np.add.at(gradient, c, -gu - gv - gw)
        return energy, gradient.ravel()

    x = x0.copy()
    for cycle in range(max_cycles):
        q = contacts(x, max(heavy_min, all_min) + .8)
        desired = np.where(heavy[q[:, 0]] & heavy[q[:, 1]], heavy_min + .18, all_min + .18)
        result = minimize(objective, x.ravel(), args=(q, desired), jac=True,
            method="L-BFGS-B", bounds=bounds,
            options={"maxiter": 150, "maxcor": 5, "ftol": 1e-10, "gtol": .002, "maxls": 30})
        x = result.x.reshape(n, 3)
        q = contacts(x, max(heavy_min, all_min))
        threshold = np.where(heavy[q[:, 0]] & heavy[q[:, 1]], heavy_min, all_min)
        distances = np.linalg.norm(x[q[:, 0]] - x[q[:, 1]], axis=1)
        violations = int(np.sum(distances < threshold))
        max_error = float(np.max(np.abs(np.linalg.norm(x[a] - x[b], axis=1) - reference)))
        flipped = int(np.sum(volumes(x)[0] * reference_volume <= 0))
        if violations == 0 and max_error < .06 and flipped == 0:
            conf = mol.GetConformer()
            for i, point in enumerate(x):
                conf.SetAtomPosition(i, tuple(map(float, point)))
            return {"method": "finite geometric overlap penalty with covalent and chirality constraints",
                    "cycles": cycle + 1, "max_geometry_distance_error_angstrom": max_error,
                    "chirality_flips": 0, "nonbonded_distance_violations": 0}
    raise RuntimeError("Coordinate push-off did not reach its distance/geometry requirements; no usable system written")


def build_polymer_system(chains=30, repeat_count=90, la_fraction=0.5, d_fraction=0.5,
                         seed=1729, box_nm=(8., 8., 14.), slab_nm=(8., 8., 8.),
                         end_groups="hydroxyl_carboxyl", min_nonbonded_angstrom=1.4,
                         max_proposals_per_chain=50000, progress_callback=None):
    """Return an explicit-H combined Mol and metadata for bounded initial coils.

    LA/GA counts and D/L counts are balanced per chain (nearest integer), then
    independently shuffled. Spatial RNG streams are independent of sequence
    generation, so increasing packing effort cannot change polymer chemistry.
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    if chains < 1 or repeat_count < 1 or not 0 <= la_fraction <= 1 or not 0 <= d_fraction <= 1:
        raise ValueError("invalid chain count, repeat count or composition")
    if not np.isfinite(min_nonbonded_angstrom) or min_nonbonded_angstrom < 1.2:
        raise ValueError("minimum nonbonded heavy distance must be at least 1.2 Angstrom")
    if max_proposals_per_chain < 1:
        raise ValueError("max_proposals_per_chain must be positive")
    box = np.asarray(box_nm, dtype=float) * 10.
    slab = np.asarray(slab_nm, dtype=float) * 10.
    if box.shape != (3,) or slab.shape != (3,) or np.any(slab <= 10.) or np.any(slab > box):
        raise ValueError("slab must fit inside box and all slab dimensions must exceed 1 nm")
    if end_groups != "hydroxyl_carboxyl":
        raise ValueError("Coordinate builder currently supports hydroxyl_carboxyl; build_graph also supports methyl ester")
    lower = (box - slab) / 2. + 1.5
    upper = (box + slab) / 2. - 1.5
    streams = np.random.SeedSequence(seed).spawn(2 * chains)
    all_mols, chain_metadata = [], []
    offset = 0
    specs, iterators, grown = [], [], {}
    for chain in range(chains):
        rng_seq = np.random.default_rng(streams[2 * chain])
        rng_xyz = np.random.default_rng(streams[2 * chain + 1])
        n_la = round(repeat_count * la_fraction)
        sequence = ["LA"] * n_la + ["GA"] * (repeat_count - n_la)
        rng_seq.shuffle(sequence)
        labels = ["D"] * round(n_la * d_fraction) + ["L"] * (n_la - round(n_la * d_fraction))
        rng_seq.shuffle(labels)
        label_iter = iter(labels)
        stereo = [next(label_iter) if kind == "LA" else None for kind in sequence]
        heavy, info = build_graph(sequence, stereo, end_groups, explicit_hydrogens=False)
        specs.append((sequence, stereo, heavy, info))
        iterators.append(_grow_chain_iter(sequence, stereo, rng_xyz, lower, upper,
                          _Cells(min_nonbonded_angstrom), chain, max_proposals=max_proposals_per_chain,
                          all_cells=_Cells(1.2)))
    # Independent coils first; movable-chain overlap removal follows below.
    # Keeping spatial RNG streams independent also keeps chemistry reproducible.
    pending = list(range(chains))
    while pending:
        still_pending = []
        for chain in pending:
            try:
                repeat = next(iterators[chain])
                if progress_callback is not None:
                    progress_callback(chain, repeat)
                still_pending.append(chain)
            except StopIteration as done:
                grown[chain] = done.value
        pending = still_pending
    for chain, (sequence, stereo, heavy, info) in enumerate(specs):
        xyz, statistics, hydrogen_coords = grown[chain]
        heavy = _coordinate_chirality(heavy, info, xyz)
        mol = Chem.AddHs(heavy, addCoords=True)
        _, complete_info = build_graph(sequence, stereo, end_groups, explicit_hydrogens=True)
        info = complete_info
        conf = mol.GetConformer()
        for atom in info["atoms"]:
            if atom["parent_role"] is not None:
                conf.SetAtomPosition(atom["index"], tuple(map(float,
                    hydrogen_coords[(atom["repeat"], atom["role"]) ])))
        info.update({"chain": chain, "atom_offset": offset, "atom_count": mol.GetNumAtoms(),
                     "molar_mass_g_mol": Descriptors.MolWt(mol), **statistics})
        info["counts"] = {"LA": sequence.count("LA"), "GA": sequence.count("GA"),
                          "D": stereo.count("D"), "L": stereo.count("L")}
        # Individual-chain metadata indices remain local for antechamber/MOL2.
        for atom in mol.GetAtoms():
            entry = info["atoms"][atom.GetIdx()]
            residue = Chem.AtomPDBResidueInfo()
            role = entry["role"].replace(":", "")
            residue.SetName(f"{role[:4]:>4}")
            residue.SetResidueName("LAC" if entry["monomer"] == "LA" else "GLY")
            residue.SetResidueNumber(entry["repeat"] + 1)
            residue.SetChainId(chr(ord("A") + chain % 26))
            residue.SetIsHeteroAtom(True)
            atom.SetMonomerInfo(residue)
        all_mols.append(mol)
        chain_metadata.append(info)
        offset += mol.GetNumAtoms()
    combined = all_mols[0]
    for mol in all_mols[1:]:
        combined = Chem.CombineMols(combined, mol)
    packing = (_soft_push_off(combined, lower, upper, heavy_min=min_nonbonded_angstrom)
               if chains > 1 else {"method": "single self-avoiding coil; no interchain packing needed"})
    mass = sum(x["molar_mass_g_mol"] for x in chain_metadata)
    metadata = {"schema": "mdprep-polymer-1", "seed": seed, "box_nm": list(box_nm),
                "initial_slab_nm": list(slab_nm), "chains": chain_metadata,
                "atom_count": combined.GetNumAtoms(), "repeat_count_per_chain": repeat_count,
                "la_fraction_requested": la_fraction, "d_fraction_requested": d_fraction,
                "min_nonlocal_heavy_distance_angstrom": min_nonbonded_angstrom,
                "min_nonlocal_all_atom_distance_angstrom": 1.2,
                "coordinate_push_off": packing,
                "nominal_polymer_density_g_cm3": mass / (602.214076 * float(np.prod(slab_nm))),
                "status": "unequilibrated initial coordinates; minimize before MD",
                "end_groups": end_groups}
    return combined, metadata


def write_polymer_system(mol, metadata, output_directory):
    """Write complete PDB, and individual chain SDF/atom-map JSON files."""
    from rdkit import Chem
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    Chem.MolToPDBFile(mol, str(directory / "polymer.pdb"))
    (directory / "polymer.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (directory / "polymer.rdkit").write_bytes(mol.ToBinary(Chem.PropertyPickleOptions.CoordsAsDouble))
    xyz = mol.GetConformer().GetPositions()
    for i, info in enumerate(metadata["chains"], 1):
        fragment, _ = build_graph(info["sequence"], info["stereo"], metadata["end_groups"])
        conf = Chem.Conformer(fragment.GetNumAtoms())
        offset = info["atom_offset"]
        for j, point in enumerate(xyz[offset:offset + fragment.GetNumAtoms()]):
            conf.SetAtomPosition(j, tuple(map(float, point)))
        fragment.AddConformer(conf)
        stem = directory / f"chain_{i:03d}"
        with Chem.SDWriter(str(stem) + ".sdf") as writer:
            writer.SetForceV3000(True)
            writer.write(fragment)
        Path(str(stem) + ".json").write_text(json.dumps(info, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chains", type=int, default=30)
    parser.add_argument("--repeat-count", type=int, default=90)
    parser.add_argument("--la-fraction", type=float, default=.5)
    parser.add_argument("--d-fraction", type=float, default=.5)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--box-nm", type=float, nargs=3, default=(8., 8., 14.))
    parser.add_argument("--slab-nm", type=float, nargs=3, default=(8., 8., 8.))
    parser.add_argument("--minimum-heavy-distance", type=float, default=1.4)
    args = parser.parse_args()
    mol, metadata = build_polymer_system(chains=args.chains, repeat_count=args.repeat_count,
        la_fraction=args.la_fraction, d_fraction=args.d_fraction, seed=args.seed,
        box_nm=args.box_nm, slab_nm=args.slab_nm,
        min_nonbonded_angstrom=args.minimum_heavy_distance)
    write_polymer_system(mol, metadata, args.output)
    print(json.dumps({key: metadata[key] for key in
                      ("atom_count", "nominal_polymer_density_g_cm3", "status")}, indent=2))


if __name__ == "__main__":
    main()
