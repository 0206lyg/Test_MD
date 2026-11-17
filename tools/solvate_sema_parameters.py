#!/usr/bin/env python3
"""Rebuild the bundled OPC/NaCl box from the dry Amber parameter directory.

Developer/reproducibility utility, NOT used at job submission. Requires
AmberTools 24 and ParmEd. The shipped files avoid rebuilding parameters on HPC.
"""
import argparse
import collections
import hashlib
import json
import re
from pathlib import Path
import shutil
import subprocess

import numpy as np
import parmed as pmd


def volume_nm3(structure):
    vectors = structure.box_vectors.value_in_unit(pmd.unit.angstrom)
    return abs(float(np.linalg.det(np.asarray(vectors,dtype=float)))) / 1000.


def leap(commands, directory, filename):
    path = directory / filename
    path.write_text('\n'.join(commands + ['quit'])+'\n')
    with (directory / (filename+'.log')).open('w') as log:
        subprocess.run(['tleap','-f',filename],cwd=directory,stdout=log,
                       stderr=subprocess.STDOUT,check=True)


def ndx_group(name, indices):
    values = [str(i+1) for i in indices]
    return '[ %s ]\n' % name + '\n'.join(' '.join(values[i:i+15])
                for i in range(0,len(values),15))+'\n\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parameters', required=True, type=Path)
    parser.add_argument('--assets', required=True, type=Path)
    parser.add_argument('--provenance', required=True, type=Path)
    parser.add_argument('--padding-angstrom',type=float,default=10.)
    parser.add_argument('--salt-reference-volume-nm3',type=float,
                        help='Representative equilibrated volume from a pilot, to avoid initial packing-volume bias')
    args = parser.parse_args()
    work, assets = args.parameters.resolve(), args.assets.resolve()
    assets.mkdir(parents=True,exist_ok=True)
    start = ['source leaprc.protein.ff14SB', 'source leaprc.water.opc']
    for name in ('AIB','SLY'):
        start += ['loadamberprep %s/%s.prepi' % (name,name),
                  'loadamberparams %s/%s.frcmod' % (name,name)]
    start += ['mol = loadpdb solute_heavy.pdb',
              'solvateOct mol OPCBOX %.6f' % args.padding_angstrom]
    leap(start+['saveamberparm mol water_box.prmtop water_box.inpcrd'],work,'water_box.leap')
    water_box = pmd.load_file(str(work/'water_box.prmtop'),xyz=str(work/'water_box.inpcrd'))
    volume = volume_nm3(water_box)
    salt_volume = args.salt_reference_volume_nm3 or volume
    if salt_volume<=0 or not np.isfinite(salt_volume):
        raise ValueError('Salt reference volume must be positive and finite')
    pairs = int(round(.150 * .602214076 * salt_volume))
    leap(start+['addIonsRand mol Na+ 0',
                'addIonsRand mol Na+ %d Cl- %d' % (pairs,pairs),
                'check mol', 'charge mol',
                'saveamberparm mol system.prmtop system.inpcrd',
                'savepdb mol system.pdb'],work,'system.leap')
    full = pmd.load_file(str(work/'system.prmtop'),xyz=str(work/'system.inpcrd'))
    dry = pmd.load_file(str(work/'solute.prmtop'),xyz=str(work/'solute.inpcrd'))
    n = len(dry.atoms)
    atoms = full.atoms[:n]
    assert [a.name for a in atoms] == [a.name for a in dry.atoms]
    formula = dict(collections.Counter(a.element_name for a in atoms))
    assert formula == {'C':187,'H':287,'N':45,'O':59}, formula
    expected_names = 'HIE AIB GLU GLY THR PHE THR SER ASP VAL SER SER TYR LEU GLU GLY GLN ALA ALA SLY GLU PHE ILE ALA TRP LEU VAL ARG GLY ARG GLY'.split()
    assert [r.name for r in dry.residues] == expected_names
    q = sum(a.charge for a in atoms)
    assert abs(q+4)<1e-4 and abs(sum(a.charge for a in full.atoms))<1e-4
    assert all(b.type is not None and b.type.k>0 for b in dry.bonds)
    assert all(a.type is not None and a.type.k>0 for a in dry.angles)
    assert all(d.type is not None for d in full.dihedrals)
    assert all(a.atom_type is not None for a in full.atoms)
    full.save(str(assets/'topol.top'),overwrite=True)
    full.save(str(assets/'system.gro'),overwrite=True)
    # A whole-solute harmonic restraint during NVT/NPT only.
    top = (assets/'topol.top').read_text()
    sections = list(re.finditer(r'^\[ moleculetype \]',top,flags=re.M))
    if len(sections)<2:
        raise RuntimeError('Expected separate solute and solvent molecule types')
    insertion = '\n#ifdef SEMA_POSRES\n#include "posre_sema.itp"\n#endif\n\n'
    top = top[:sections[1].start()] + insertion + top[sections[1].start():]
    (assets/'topol.top').write_text(top)
    heavy = [a.idx for a in atoms if a.atomic_number>1]
    (assets/'posre_sema.itp').write_text('[ position_restraints ]\n; atom type kx ky kz (kJ/mol/nm^2)\n' + ''.join('%d 1 1000 1000 1000\n' % (i+1) for i in heavy))
    tail_names = {'C%d'%i for i in range(24,42)} | {'O5','O6','O7'}
    peptide_sly = {'N','CA','C','O','C4','C5','C6','C7','N2'}
    tail = [a.idx for a in atoms if a.residue.idx==19 and a.name in tail_names]
    peptide = [a.idx for a in atoms if a.atomic_number>1 and
               (a.residue.idx!=19 or a.name in peptide_sly)]
    backbone = [a.idx for a in atoms if a.name in ('N','CA','C')]
    groups = [('System',range(len(full.atoms))),('Semaglutide',range(n)),
              ('Backbone',backbone),('Peptide',peptide),('C18',tail),
              ('Water_and_ions',range(n,len(full.atoms)))]
    (assets/'index.ndx').write_text(''.join(ndx_group(*group) for group in groups))
    parameters = json.loads((work/'parameter_audit.json').read_text())
    solvent_counts = collections.Counter(r.name for r in full.residues[31:])
    na,cl = solvent_counts.get('Na+',0),solvent_counts.get('Cl-',0)
    waters = solvent_counts.get('WAT',0)
    assert na==pairs+4 and cl==pairs and waters>0,solvent_counts
    source = json.loads(args.provenance.read_text())
    source['force_field'] = ('ff14SB canonical residues; custom AMBER-type AIB and lipidated Lys with '+parameters['charge_method']+'; parmchk2 parameter analogies; OPC water; Li-Merz 12-6 OPC ions')
    source['parameter_status'] = 'provisional; not the published CHARMM36m/CGenFF force field'
    source['conversion'] = 'AmberTools 24.8 -> ParmEd 4.3.1 -> GROMACS explicit topology'
    stereocenters = []
    def stereo(residue,center,neighbors,expected):
        names = {a.name:a.idx for a in residue.atoms}
        i,near = names[center],[names[name] for name in neighbors]
        signed_volume = float(np.linalg.det((full.coordinates[near]-full.coordinates[i])/10.))
        actual = 'S' if signed_volume>0 else 'R'
        if actual!=expected:
            raise ValueError('Wrong initial stereochemistry: residue %d %s' % (residue.idx+1,center))
        stereocenters.append(dict(resnr=residue.idx+1,name=center,center=i,
             neighbors=near,expected=expected,initial_signed_volume_nm3=signed_volume))
    for residue in dry.residues:
        if residue.name not in ('GLY','AIB'):
            stereo(residue,'CA',['N','C','C4' if residue.name=='SLY' else 'CB'],'S')
        if residue.name=='THR':
            stereo(residue,'CB',['OG1','CA','CG2'],'R')
        if residue.name=='ILE':
            stereo(residue,'CB',['CA','CG1','CG2'],'S')
        if residue.name=='SLY':
            stereo(residue,'C23',['N5','C42','C22'],'S')
    manifest = dict(schema_version=1,atoms=[dict(index=a.idx,name=a.name,
        resnr=a.residue.idx+1,resname=a.residue.name,element=a.element_name,
        mass=a.mass,charge=a.charge) for a in atoms],
        bonds=[[b.atom1.idx,b.atom2.idx] for b in dry.bonds],
        peptide_backbone=[dict(resnr=r.idx+1,**{key:next((a.idx for a in r.atoms if a.name==key),None) for key in ('N','CA','C','O')}) for r in dry.residues],
        tail_atom_indices=tail,peptide_atom_indices=peptide,expected_formal_charge=-4,
        stereocenters=stereocenters,
        chemical_audit=dict(status='pass',residue_count=31,heavy_atoms=291,
            sequence=expected_names,formula_at_selected_state=formula,
            checked_stereocenters=len(stereocenters),
            notes=['Full terminal Gly retained; N terminus NH3+, C terminus COO-.',
                   'HIE neutral histidine is an assumed tautomer at pH 7.4.',
                   'Asp/Glu, gamma-Glu and C18 terminal carboxylates deprotonated; Arg protonated; acylated Lys neutral.',
                   'Gamma-Glu geometry correction: '+parameters.get('gamma_glu_geometry_correction','unknown')]),
        parameter_audit=dict(coverage_complete=True,status='provisional_physics_unvalidated',details=parameters),
        protonation=dict(target_pH=7.4,method='fixed-state assumption, not constant-pH MD',histidine='HIE'),
        provenance=source,system_atom_count=len(full.atoms),water_molecules=waters,
        physical_atoms=sum(a.atomic_number>0 for a in full.atoms),
        virtual_sites=sum(a.atomic_number==0 for a in full.atoms),
        ions=dict(Na=na,Cl=cl,neutralizing_Na=4,NaCl_pairs=pairs,
                  requested_added_salt_M=.150,initial_added_salt_M=pairs/(.602214076*volume),
                  salt_count_reference_volume_nm3=salt_volume,
                  concentration_at_reference_volume_M=pairs/(.602214076*salt_volume),
                  count_basis='representative equilibrated pilot volume' if args.salt_reference_volume_nm3 else 'initial packing volume',
                  note='Salt-pair concentration excludes neutralizing Na; NPT volume changes concentration. Initial packing concentration is lower; report the sampled production concentration.'),
        box_volume_nm3=volume,box_shape='truncated octahedron',
        box_lengths_nm=(full.box[:3]/10).tolist(),box_angles_deg=full.box[3:].tolist(),
        initial_padding_nm=args.padding_angstrom/10,position_restraints=True)
    (assets/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    reproducibility = assets/'parameter_source'
    reproducibility.mkdir(exist_ok=True)
    for name in ('solute_heavy.pdb','solute.pdb','solute.prmtop','solute.inpcrd','parameter_audit.json','build.leap','system.leap'):
        shutil.copy2(work/name,reproducibility/name)
    for name in ('AIB','SLY'):
        dest = reproducibility/name
        dest.mkdir(exist_ok=True)
        for item in (name+'.prepi',name+'.frcmod','capped.mol2','charged.ac','typed.ac','residue.mc','PARMCHK_CX_CT.DAT'):
            shutil.copy2(work/name/item,dest/item)
    checks = {str(p.relative_to(assets)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(assets.rglob('*')) if p.is_file() and p.name!='checksums.json'}
    (assets/'checksums.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ('system_atom_count','physical_atoms','water_molecules','ions','box_volume_nm3')},indent=2))


if __name__=='__main__':
    main()
