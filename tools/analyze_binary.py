#!/usr/bin/env python3
"""Molecule-based concentration diagnostics for a binary solvent trajectory.

Input: a concatenated, orthorhombic GRO trajectory, e.g. from
  printf '0\n' | gmx trjconv -s production.tpr -f production.xtc \
      -o analysis/frames.gro -dt 100 -pbc mol
Only NumPy is required. Atom counts/order must match prepare_binary.py's manifest.

For molecule centres r_i, species labels s_i=1 (organic) or 0 (water), and the
global organic mole fraction x, the reported concentration structure factor is
Scc(k)=|sum_i (s_i-x) exp(2*pi*j*k.r_i/L)|^2/[N*x*(1-x)].
Integer k is sampled in three lowest shells |k|^2=1,2,3. Each +/- pair is used
once. Independent uniformly distributed molecular positions have E[Scc]=1.
Real homogeneous liquids can differ from 1 because of molecular correlations;
this baseline is a reference, not a phase-separation decision threshold.
The Fourier magnitude is invariant under a common periodic translation, unlike
a time-averaged laboratory-frame density profile that can hide moving domains.
All concentration measures count molecules, not atoms or volume fractions.
The random-corrected mixing index is 1 for random species labels on the observed
centres and 0 for complete cell-wise separation at the chosen grid scale. It is
not clipped and is not a thermodynamic phase label. Outputs include 2D projected
maps and ParaView 3D cell fields / organic-molecule point time series in nm and ns.
"""
from __future__ import annotations

import argparse
import csv
import html
import itertools
import json
import math
from pathlib import Path
import re

import numpy as np


def read_layout(path):
    manifest = json.loads(Path(path).read_text())
    species = manifest.get("species", [])
    if len(species) != 2 or {s.get("kind") for s in species} != {"water", "solvent"}:
        raise ValueError("Expected exactly one water and one organic-solvent species")
    offset = 0
    layout = []
    labels = []
    for species_record in species:
        count = species_record["count"]
        size = species_record["atoms_per_molecule"]
        if any(not isinstance(v, int) or isinstance(v, bool) or v <= 0 for v in (count, size)):
            raise ValueError("Manifest molecule and atom counts must be positive integers")
        stop = offset + count * size
        if species_record.get("atom_start_0based", offset) != offset or species_record.get("atom_stop_0based", stop) != stop:
            raise ValueError("Manifest atom offsets disagree with species order")
        layout.append((offset, stop, count, size))
        labels.extend([int(species_record["kind"] == "solvent")] * count)
        offset = stop
    if manifest.get("atom_count") != offset:
        raise ValueError("Manifest total atom count disagrees with species counts")
    return manifest, layout, np.asarray(labels, dtype=float)


def gro_frames(path, expected_atoms):
    """Stream standard or extended-coordinate-precision GRO frames in nm."""
    previous_time = None
    with Path(path).open() as stream:
        frame_number = 0
        while True:
            title = stream.readline()
            if not title:
                break
            count_line = stream.readline()
            if not count_line:
                raise ValueError("Truncated GRO before atom count")
            try:
                count = int(count_line)
            except ValueError as error:
                raise ValueError(f"Invalid atom count at frame {frame_number}") from error
            if count != expected_atoms:
                raise ValueError(f"Frame {frame_number}: {count} atoms; expected {expected_atoms}")
            match = re.search(r"\bt\s*=\s*([-+\d.eE]+)", title)
            if match is None:
                raise ValueError("GRO frame title lacks t= time; extract a trajectory using gmx trjconv")
            time = float(match.group(1))
            if not math.isfinite(time) or (previous_time is not None and time <= previous_time):
                raise ValueError("Trajectory times must be finite and strictly increasing")
            previous_time = time
            xyz = np.empty((count, 3), dtype=float)
            width = None
            for i in range(count):
                line = stream.readline()
                if width is None:
                    decimals = [m.start() for m in re.finditer(r"\.", line[20:])]
                    if len(decimals) < 3:
                        raise ValueError("Malformed GRO coordinate fields")
                    width = decimals[1] - decimals[0]
                    if width < 5 or decimals[2] - decimals[1] != width:
                        raise ValueError("Inconsistent GRO coordinate precision")
                try:
                    xyz[i] = [float(line[20 + j * width:20 + (j + 1) * width]) for j in range(3)]
                except ValueError as error:
                    raise ValueError(f"Invalid coordinates at frame {frame_number}, atom {i + 1}") from error
            box_values = [float(v) for v in stream.readline().split()]
            if len(box_values) not in (3, 9):
                raise ValueError("GRO box must contain three or nine values")
            if len(box_values) == 9 and not np.allclose(box_values[3:], 0, atol=1e-8):
                raise ValueError("Only orthorhombic boxes are supported")
            box = np.asarray(box_values[:3])
            if not np.isfinite(xyz).all() or not np.isfinite(box).all() or np.any(box <= 0):
                raise ValueError("Nonfinite coordinates or invalid box")
            yield time, xyz, box
            frame_number += 1
    if frame_number == 0:
        raise ValueError("No GRO frames were found")


def molecule_centres(xyz, box, layout):
    blocks = []
    for start, stop, count, size in layout:
        atoms = xyz[start:stop].reshape(count, size, 3)
        anchor = atoms[:, :1, :]
        displacement = atoms - anchor
        displacement -= box * np.rint(displacement / box)
        blocks.append(np.mod(anchor[:, 0, :] + displacement.mean(axis=1), box))
    return np.concatenate(blocks, axis=0)


def concentration_statistics(centres, box, labels, grid=4):
    fraction = labels.mean()
    contrast = labels - fraction
    fractional = np.mod(centres / box, 1.0)
    shells = {1: [], 2: [], 3: []}
    for mode in itertools.product((-1, 0, 1), repeat=3):
        nonzero = [v for v in mode if v]
        if not nonzero or nonzero[0] < 0:
            continue
        amplitude = np.sum(contrast * np.exp(2j * np.pi * (fractional @ np.asarray(mode))))
        shells[sum(v*v for v in mode)].append(float(abs(amplitude)**2 / (len(labels)*fraction*(1-fraction))))
    total, organic = grid_counts(centres, box, labels, grid)
    occupied = total > 0
    # The expected weighted variance for random species assignment at fixed
    # overall counts is x(1-x)*(occupied_cells-1)/(N-1).
    weighted_variance = float(np.sum((organic[occupied]-fraction*total[occupied])**2 / total[occupied]) / len(labels))
    random_variance = float(fraction*(1-fraction)*(occupied.sum()-1)/(len(labels)-1))
    result = {f"Scc_shell{shell}_mean": float(np.mean(values)) for shell, values in shells.items()}
    result["Scc_low_modes_max"] = max(max(values) for values in shells.values())
    result.update({"local_x_weighted_variance": weighted_variance,
                   "random_label_variance_reference": random_variance,
                   "local_variance_over_random": weighted_variance/random_variance if random_variance > 0 else 0.0})
    denominator = fraction*(1-fraction)-random_variance
    segregation = (weighted_variance-random_variance)/denominator if denominator > 1e-15 else None
    result.update({"segregation_excess": segregation,
                   "mixing_index_random_corrected": 1-segregation if segregation is not None else None,
                   "occupied_cells": int(occupied.sum()),
                   "mean_organic_molecules_per_cell": float(labels.sum()/grid**3)})
    return result


def grid_counts(centres, box, labels, grid):
    """Return arrays indexed [x,y,z]; VTK serialization separately uses x-fast order."""
    cell = np.minimum((np.mod(centres/box, 1.0)*grid).astype(int), grid-1)
    flat = np.ravel_multi_index(cell.T, (grid, grid, grid))
    total = np.bincount(flat, minlength=grid**3).reshape((grid,)*3)
    organic = np.bincount(flat, weights=labels, minlength=grid**3).reshape((grid,)*3)
    return total, organic


def molecular_fraction(total, organic):
    result = np.full(total.shape, np.nan, dtype=float)
    np.divide(organic, total, out=result, where=total > 0)
    return result


def vtk_values(values, integer=False):
    flat = np.asarray(values).ravel(order="F")
    return " ".join(str(int(v)) if integer else f"{float(v):.9g}" for v in flat)


def write_vti(path, total, organic, box, fraction):
    grid = total.shape[0]
    cell_volume = float(np.prod(box)/grid**3)
    xlocal = molecular_fraction(total, organic)
    fields = {"organic_mole_fraction": xlocal, "relative_enrichment": xlocal/fraction,
              "organic_number_density_nm3": organic/cell_volume,
              "water_number_density_nm3": (total-organic)/cell_volume,
              "organic_molecule_count": organic, "total_molecule_count": total,
              "occupied_cell": (total > 0).astype(int)}
    extent = f"0 {grid} 0 {grid} 0 {grid}"
    spacing = " ".join(f"{v/grid:.9g}" for v in box)
    lines = ['<?xml version="1.0"?>', '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">',
             f'<ImageData WholeExtent="{extent}" Origin="0 0 0" Spacing="{spacing}">',
             f'<Piece Extent="{extent}"><CellData Scalars="relative_enrichment">']
    for name, values in fields.items():
        integer = name in ("organic_molecule_count", "total_molecule_count", "occupied_cell")
        lines.append(f'<DataArray type="{"Int32" if integer else "Float64"}" Name="{name}" format="ascii">{vtk_values(values, integer)}</DataArray>')
    lines += ['</CellData><PointData/></Piece></ImageData></VTKFile>']
    Path(path).write_text("\n".join(lines)+"\n")


def write_vtp(path, organic_centres):
    count = len(organic_centres)
    # A point needs a vertex cell to be visible as points/glyphs in ParaView.
    coordinates = " ".join(f"{v:.9g}" for v in organic_centres.ravel())
    indices = " ".join(map(str, range(count)))
    offsets = " ".join(map(str, range(1, count+1)))
    Path(path).write_text(f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
<PolyData><Piece NumberOfPoints="{count}" NumberOfVerts="{count}" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">{coordinates}</DataArray></Points>
<PointData Scalars="organic_molecule_index"><DataArray type="Int32" Name="organic_molecule_index" format="ascii">{indices}</DataArray></PointData>
<Verts><DataArray type="Int32" Name="connectivity" format="ascii">{indices}</DataArray><DataArray type="Int32" Name="offsets" format="ascii">{offsets}</DataArray></Verts>
</Piece></PolyData></VTKFile>
''')


def write_pvd(path, frames):
    lines = ['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>']
    lines += [f'<DataSet timestep="{time/1000:.12g}" group="" part="0" file="{name}"/>' for time, name in frames]
    lines.append('</Collection></VTKFile>')
    Path(path).write_text("\n".join(lines)+"\n")


def projection_fields(total, organic, fraction):
    """Column-integrated composition, not a central slice or mean of cell fractions."""
    return [("XY", (0, 1), molecular_fraction(total.sum(axis=2), organic.sum(axis=2))/fraction),
            ("XZ", (0, 2), molecular_fraction(total.sum(axis=1), organic.sum(axis=1))/fraction),
            ("YZ", (1, 2), molecular_fraction(total.sum(axis=0), organic.sum(axis=0))/fraction)]


def enrichment_colour(value):
    if not math.isfinite(value):
        return "#c7cbd1"
    # Fixed scale across solvents, frames and projections. Exact un-clipped values
    # are retained in VTK/CSV; orange saturates at >=3-fold enrichment.
    stops = ((43, 105, 155), (247, 247, 242), (194, 70, 31))
    value = min(3.0, max(0.0, value))
    a, b, weight = (stops[0], stops[1], value) if value <= 1 else (stops[1], stops[2], (value-1)/2)
    rgb = [round(u+(v-u)*weight) for u, v in zip(a, b)]
    return "#"+"".join(f"{v:02x}" for v in rgb)


def projection_svg(path, total, organic, box, fraction, title, note):
    fields = projection_fields(total, organic, fraction)
    elements = ['<svg xmlns="http://www.w3.org/2000/svg" width="1110" height="560" viewBox="0 0 1110 560">',
                '<rect width="1110" height="560" fill="white"/>', '<g font-family="sans-serif" fill="#182333">',
                f'<text x="35" y="30" font-size="20">{html.escape(title)}</text>',
                '<text x="35" y="54" font-size="13">Organic enrichment = local molecular mole fraction / global molecular mole fraction</text>']
    for index, (name, axes, values) in enumerate(fields):
        left, top, width = 65+index*360, 102, 275
        count = len(values)
        elements.append(f'<text x="{left}" y="84" font-size="16">{name} projection</text>')
        for i in range(count):
            for j in range(count):
                elements.append(f'<rect x="{left+i*width/count:.3f}" y="{top+(count-1-j)*width/count:.3f}" width="{width/count+.01:.3f}" height="{width/count+.01:.3f}" fill="{enrichment_colour(values[i,j])}"/>')
        elements.append(f'<rect x="{left}" y="{top}" width="{width}" height="{width}" fill="none" stroke="#182333"/>')
        for tick in (0, .5, 1):
            elements.append(f'<text x="{left+tick*width}" y="{top+width+20}" text-anchor="middle" font-size="12">{tick*box[axes[0]]:.3g}</text>')
            elements.append(f'<text x="{left-9}" y="{top+(1-tick)*width+4}" text-anchor="end" font-size="12">{tick*box[axes[1]]:.3g}</text>')
        elements.append(f'<text x="{left+width/2}" y="{top+width+40}" text-anchor="middle" font-size="13">{"xyz"[axes[0]]} (nm)</text>')
        elements.append(f'<text x="{left-43}" y="{top+width/2}" text-anchor="middle" transform="rotate(-90 {left-43} {top+width/2})" font-size="13">{"xyz"[axes[1]]} (nm)</text>')
    for i in range(300):
        elements.append(f'<rect x="{385+i}" y="445" width="1.1" height="16" fill="{enrichment_colour(3*i/299)}"/>')
    for value in (0, 1, 2, 3):
        elements.append(f'<text x="{385+100*value}" y="480" text-anchor="middle" font-size="12">{"3+" if value==3 else value}</text>')
    maximum = max(float(np.nanmax(v)) for _, _, v in fields if np.any(np.isfinite(v)))
    elements += [f'<text x="35" y="509" font-size="12">{html.escape(note)}</text>',
                 f'<text x="35" y="533" font-size="12">Grid: {len(total)} per axis; max projected enrichment = {maximum:.3g}; grey = empty column. Projections can hide depth variation.</text>', '</g></svg>']
    Path(path).write_text("\n".join(elements)+"\n")


def dispersion_svg(path, rows):
    points = [(r["time_ps"]/1000, r["mixing_index_random_corrected"]) for r in rows if r["mixing_index_random_corrected"] is not None]
    if not points:
        return
    x0, x1 = points[0][0], points[-1][0]
    if x1 <= x0:
        x0 = max(0, x0-.05)
        x1 = x0+.1
    low = min(0, min(y for _, y in points))-.05
    high = max(1, max(y for _, y in points))+.05
    def xx(x): return 75+780*(x-x0)/(x1-x0)
    def yy(y): return 370-280*(y-low)/(high-low)
    lines = ['<svg xmlns="http://www.w3.org/2000/svg" width="920" height="480" viewBox="0 0 920 480">',
             '<rect width="920" height="480" fill="white"/><g font-family="sans-serif" fill="#182333">',
             '<text x="35" y="32" font-size="20">Molecular mixing index, corrected for random-label fluctuations</text>']
    for i in range(5):
        x, y = x0+(x1-x0)*i/4, low+(high-low)*i/4
        lines += [f'<path d="M75,{yy(y)}H855" stroke="#ddd"/>',
                  f'<text x="65" y="{yy(y)+4}" text-anchor="end" font-size="12">{y:.3g}</text>',
                  f'<text x="{xx(x)}" y="393" text-anchor="middle" font-size="12">{x:.3g}</text>']
    for y in (0, 1):
        lines.append(f'<path d="M75,{yy(y)}H855" stroke="#777" stroke-dasharray="5 4"/>')
    coordinates = " ".join(f"{xx(x):.3f},{yy(y):.3f}" for x, y in points)
    lines += [f'<polyline points="{coordinates}" fill="none" stroke="#176f9d" stroke-width="2"/>',
              f'<circle cx="{xx(points[-1][0])}" cy="{yy(points[-1][1])}" r="3" fill="#176f9d"/>',
              '<path d="M75,90V370H855" fill="none" stroke="#182333"/>',
              '<text x="465" y="419" text-anchor="middle" font-size="13">Trajectory time (ns)</text>',
              '<text x="26" y="230" transform="rotate(-90 26 230)" text-anchor="middle" font-size="13">Mixing index M</text>',
              '<text x="35" y="448" font-size="12">M = 1: random-label reference; M = 0: complete species separation at the selected grid scale.</text>',
              '<text x="35" y="468" font-size="12">Values above 1 are retained. These are spatial diagnostics, not phase-equilibrium thresholds.</text></g></svg>']
    Path(path).write_text("\n".join(lines)+"\n")


def profiles(centres, box, labels, bins):
    for axis in range(3):
        fractional = np.mod(centres[:, axis]/box[axis], 1.0)
        organic, _ = np.histogram(fractional[labels == 1], bins=bins, range=(0, 1))
        water, _ = np.histogram(fractional[labels == 0], bins=bins, range=(0, 1))
        for b in range(bins):
            total = int(organic[b] + water[b])
            yield {"axis": "xyz"[axis], "bin": b, "fractional_position": (b+.5)/bins,
                   "position_nm": (b+.5)*box[axis]/bins,
                   "organic_molecules": int(organic[b]), "water_molecules": int(water[b]),
                   "organic_mole_fraction": float(organic[b]/total) if total else None,
                   "organic_number_density_nm3": float(organic[b]/(np.prod(box)/bins)),
                   "water_number_density_nm3": float(water[b]/(np.prod(box)/bins))}


def diagnostic_svg(path, rows, final_profiles, label):
    """A dependency-free overview; full time-resolved profiles remain in CSV."""
    elements = ['<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="510" viewBox="0 0 1080 510">',
                '<rect width="1080" height="510" fill="white"/>',
                '<g font-family="sans-serif" fill="#182333">',
                f'<text x="55" y="32" font-size="21">{html.escape(label)} / water: concentration diagnostics</text>']
    colours = ["#176f9d", "#be5429", "#357d48"]

    def panel(left, title, curves, xlo, xhi, ymax, xlabel, ylabel, baseline):
        top, width, height = 95, 430, 300
        if xhi <= xlo:
            xlo = max(0, xlo-.05)
            xhi = xlo+.1
        ymax = max(ymax, 1e-12)
        elements.append(f'<text x="{left}" y="65" font-size="15">{title}</text>')
        for i in range(5):
            xp = left+width*i/4
            yp = top+height-height*i/4
            elements.append(f'<path d="M{left},{yp}h{width}" stroke="#e2e6eb"/>')
            elements.append(f'<text x="{left-8}" y="{yp+4}" text-anchor="end" font-size="11">{ymax*i/4:.3g}</text>')
            elements.append(f'<text x="{xp}" y="{top+height+20}" text-anchor="middle" font-size="11">{xlo+(xhi-xlo)*i/4:.3g}</text>')
        if 0 <= baseline <= ymax:
            yp = top+height-height*baseline/ymax
            elements.append(f'<path d="M{left},{yp}h{width}" stroke="#777" stroke-dasharray="5 4"/>')
        for colour, (name, xs, ys) in zip(colours, curves):
            points = " ".join(f"{left+width*(x-xlo)/(xhi-xlo):.2f},{top+height-height*y/ymax:.2f}" for x, y in zip(xs, ys))
            elements.append(f'<polyline points="{points}" fill="none" stroke="{colour}" stroke-width="1.6"/>')
            if len(xs) == 1:
                elements.append(f'<circle cx="{left+width*(xs[0]-xlo)/(xhi-xlo)}" cy="{top+height-height*ys[0]/ymax}" r="3" fill="{colour}"/>')
            index = colours.index(colour)
            elements.append(f'<text x="{left+index*140}" y="{top+height+63}" fill="{colour}" font-size="12">{name}</text>')
        elements.append(f'<path d="M{left},{top}v{height}h{width}" fill="none" stroke="#182333"/>')
        elements.append(f'<text x="{left+width/2}" y="{top+height+41}" text-anchor="middle" font-size="12">{xlabel}</text>')
        elements.append(f'<text x="{left-42}" y="{top+height/2}" text-anchor="middle" transform="rotate(-90 {left-42} {top+height/2})" font-size="12">{ylabel}</text>')

    times = [r["time_ps"]/1000 for r in rows]
    curves = [(f"|k|² = {s}", times, [r[f"Scc_shell{s}_mean"] for r in rows]) for s in (1,2,3)]
    panel(62, "Long-wavelength concentration modes", curves, times[0], times[-1],
          max(1.1, max(max(c[2]) for c in curves)*1.05), "Trajectory time (ns)", "Normalized Scc", 1)
    curves = []
    for axis in "xyz":
        selected = [p for p in final_profiles if p["axis"] == axis and p["organic_mole_fraction"] is not None]
        curves.append((f"{axis} direction", [p["fractional_position"] for p in selected], [p["organic_mole_fraction"] for p in selected]))
    panel(605, f"Last-frame profiles (t = {times[-1]:.3g} ns)", curves, 0, 1,
          min(1.0, max(max(c[2]) for c in curves)*1.15), "Position / box length", "Organic mole fraction", rows[-1]["organic_mole_fraction"])
    elements.append('<text x="55" y="493" font-size="12">Dashed lines: ideal random Scc = 1 (left); global composition (right). Neither is a demixing threshold.</text>')
    elements.append('</g></svg>')
    Path(path).write_text("\n".join(elements)+"\n")


def analyze(trajectory, manifest_path, output, bins=12, grid=4):
    manifest, layout, labels = read_layout(manifest_path)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    vtk_folder = output/"vtk"
    vtk_folder.mkdir(exist_ok=True)
    rows = []
    final_profiles = []
    fraction = float(labels.mean())
    field_frames, point_frames = [], []
    previous = None
    accumulated_total = np.zeros((grid,)*3, dtype=float)
    accumulated_organic = np.zeros((grid,)*3, dtype=float)
    accumulated_box = np.zeros(3, dtype=float)
    weight_sum = 0.0
    with (output/"composition_profiles.csv").open("w", newline="") as pf:
        profile_writer = None
        for frame_index, (time, xyz, box) in enumerate(gro_frames(trajectory, manifest["atom_count"])):
            centres = molecule_centres(xyz, box, layout)
            total, organic = grid_counts(centres, box, labels, grid)
            if previous is not None:
                previous_time, previous_total, previous_organic, previous_box = previous
                dt = time-previous_time
                accumulated_total += .5*dt*(previous_total+total)
                accumulated_organic += .5*dt*(previous_organic+organic)
                accumulated_box += .5*dt*(previous_box+box)
                weight_sum += dt
            previous = (time, total, organic, box)
            field_name = f"vtk/concentration_{frame_index:06d}.vti"
            point_name = f"vtk/organic_{frame_index:06d}.vtp"
            write_vti(output/field_name, total, organic, box, fraction)
            write_vtp(output/point_name, centres[labels == 1])
            field_frames.append((time, field_name))
            point_frames.append((time, point_name))
            row = {"time_ps": time, "box_x_nm": float(box[0]), "box_y_nm": float(box[1]), "box_z_nm": float(box[2]),
                   "organic_mole_fraction": fraction, **concentration_statistics(centres, box, labels, grid)}
            rows.append(row)
            final_profiles = list(profiles(centres, box, labels, bins))
            for profile in final_profiles:
                record = {"time_ps": time, **profile}
                if profile_writer is None:
                    profile_writer = csv.DictWriter(pf, fieldnames=list(record))
                    profile_writer.writeheader()
                profile_writer.writerow(record)
    with (output/"concentration_timeseries.csv").open("w", newline="") as tf:
        writer = csv.DictWriter(tf, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    midpoint = .5*(rows[0]["time_ps"]+rows[-1]["time_ps"])
    late = [r for r in rows if r["time_ps"] >= midpoint]
    late_means = {}
    for key in rows[0]:
        if key not in ("time_ps", "organic_mole_fraction"):
            values = [r[key] for r in late if r[key] is not None]
            late_means[key] = float(np.mean(values)) if values else None
    if weight_sum == 0:
        accumulated_total = total.astype(float)
        accumulated_organic = organic.copy()
        accumulated_box = box.copy()
        weight_sum = 1.0
    mean_box = accumulated_box/weight_sum
    write_pvd(output/"concentration_3d.pvd", field_frames)
    write_pvd(output/"organic_molecules_3d.pvd", point_frames)
    projection_svg(output/"projection_last.svg", total, organic, box, fraction,
                   f"{manifest.get('solvent', 'Organic solvent')} / water: last-frame composition",
                   f"t = {rows[-1]['time_ps']/1000:.5g} ns; each panel integrates through the full third dimension.")
    projection_svg(output/"projection_time_average.svg", accumulated_total, accumulated_organic, mean_box, fraction,
                   f"{manifest.get('solvent', 'Organic solvent')} / water: time-aggregated composition",
                   f"Window {rows[0]['time_ps']/1000:.5g}-{rows[-1]['time_ps']/1000:.5g} ns; fractional-coordinate bins, mean box axes. Moving domains can wash out.")
    with (output/"projection_values.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("window", "plane", "bin_axis1", "bin_axis2", "relative_enrichment", "organic_mole_fraction"))
        for window, counts, counts_org in (("last_frame", total, organic), ("time_aggregated", accumulated_total, accumulated_organic)):
            for name, _, enrichment in projection_fields(counts, counts_org, fraction):
                for i, j in np.ndindex(enrichment.shape):
                    value = enrichment[i,j]
                    writer.writerow((window, name, i, j, value if math.isfinite(value) else "", value*fraction if math.isfinite(value) else ""))
    summary = {"solvent": manifest.get("solvent", "Organic solvent"), "frames": len(rows),
               "time_start_ps": rows[0]["time_ps"], "time_end_ps": rows[-1]["time_ps"],
               "molecules": len(labels), "organic_molecules": int(labels.sum()),
               "water_molecules": int(len(labels)-labels.sum()), "organic_mole_fraction": float(labels.mean()),
               "profile_bins_per_axis": bins, "local_grid_cells_per_axis": grid,
               "late_window_start_ps": midpoint, "late_window_frames": len(late),
               "late_window_means": late_means,
               "spatial_outputs": {"field_series": "concentration_3d.pvd", "organic_molecule_series": "organic_molecules_3d.pvd",
                                   "length_unit": "nm", "PVD_time_unit": "ns", "field_location": "CellData; VTK x-index fastest",
                                   "molecular_points": "One unweighted center of geometry per organic molecule, wrapped into the box",
                                   "projection_colour_scale": "relative enrichment x_local/x_global, fixed 0..3 with >=3 saturated",
                                   "empty_cells": "organic_mole_fraction and enrichment are NaN; occupied_cell=0; number densities=0",
                                   "projection_average": "Trapezoidal integration of organic and total column counts, then their ratio; bins follow fractional box coordinates; axes use time-average box lengths"},
               "diagnostic_definitions": {"Scc": "abs(sum((species_label-x)*exp(2*pi*i*k.r/L)))^2/(N*x*(1-x)); organic label=1, water=0; shell s uses integer |k|^2=s and one member of every +/- pair",
                                          "Scc_ideal_reference": 1.0,
                                          "centres": "Unweighted molecular centres of geometry after intramolecular minimum-image reconstruction",
                                          "local_variance": "Molecule-count-weighted variance of organic mole fraction over the periodic cubic grid",
                                          "random_label_variance_reference": "x*(1-x)*(occupied_cells-1)/(N-1) at fixed species counts",
                                          "segregation_excess": "S=(weighted_variance-random_label_variance)/(x*(1-x)-random_label_variance); E[S]=0 for random species labels on these centers, S=1 for complete cell-wise separation; no clipping",
                                          "mixing_index_random_corrected": "M=1-S; E[M]=1 for random labels, M=0 for complete cell-wise separation; may exceed 1; undefined (blank/null) if every occupied cell holds one molecule"},
               "interpretation": ["No automated phase classification is made. Ideal-random baselines are not thresholds for real liquid mixtures.",
                                  "Mixing/segregation indices depend on grid size and origin. Compare equal physical coarse-graining scales and inspect multiple grid choices.",
                                  "At the default composition a 4x4x4 grid has only about 3-5 organic molecules per cell; random count fluctuations are explicitly included in the reference.",
                                  "Inspect each-frame 3D fields: projections and time averages can conceal an interface oriented along the integrated direction or a translating domain.",
                                  "Inspect persistence/growth of concentration domains, all directional profiles, independent seeds, and box-size dependence.",
                                  "Molecular clustering and correlated concentration fluctuations can occur in a macroscopically miscible liquid.",
                                  "The input usually contains production after equilibration; times are taken from GRO titles.",
                                  "Late-window averages are descriptive and do not establish equilibration or independent statistical uncertainty."]}
    (output/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    diagnostic_svg(output/"concentration_diagnostics.svg", rows, final_profiles, summary["solvent"])
    dispersion_svg(output/"dispersion_timeseries.svg", rows)
    print(f"Analyzed {len(rows)} frames ({rows[0]['time_ps']:g} to {rows[-1]['time_ps']:g} ps): {output}")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, allow_abbrev=False)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bins", type=int, default=12)
    parser.add_argument("--grid", type=int, default=4)
    args = parser.parse_args()
    if not 2 <= args.bins <= 100 or not 2 <= args.grid <= 20:
        parser.error("--bins must be 2..100 and --grid must be 2..20")
    try:
        analyze(args.trajectory, args.manifest, args.output, args.bins, args.grid)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"Binary analysis failed: {error}\n")


if __name__ == "__main__":
    main()
