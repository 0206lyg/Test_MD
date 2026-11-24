"""Render the initial phase boundaries in nm, using the geometry used for packing."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path


def _element(parent, tag, **attrs):
    return ET.SubElement(parent, tag, {k.replace("_", "-"): str(v) for k, v in attrs.items()})


def _text(parent, x, y, value, size=15, **attrs):
    el = _element(parent, "text", x=x, y=y, font_size=size, **attrs)
    el.text = str(value)
    return el


def _pocket_labels(parent, positions, bounds):
    """Keep nearby projected IDs readable; a leader points to each actual center."""
    placed = []
    xmin, ymin, xmax, ymax = bounds
    for index, (x, y), radius in positions:
        candidates = [(x, y)] if radius > 12 else []
        for dist in (19, 30, 42):
            for dx, dy in ((1, 0), (0, -1), (-1, 0), (0, 1),
                           (.707, -.707), (-.707, -.707), (-.707, .707), (.707, .707)):
                candidates.append((x + dist * dx, y + dist * dy))
        for lx, ly in candidates:
            if not (xmin + 10 <= lx <= xmax - 10 and ymin + 10 <= ly <= ymax - 10):
                continue
            if all(math.hypot(lx - px, ly - py) > 20 for px, py in placed):
                break
        if math.hypot(lx - x, ly - y) > 1:
            _element(parent, "line", x1=x, y1=y, x2=lx, y2=ly, stroke="#236777")
            _element(parent, "circle", cx=x, cy=y, r=2.3, fill="#236777")
        _element(parent, "circle", cx=lx, cy=ly, r=9, fill="#ffffff", fill_opacity=.91)
        _text(parent, lx, ly + 4, index, 11, fill="#174e5d",
              text_anchor="middle", font_weight="bold")
        placed.append((lx, ly))


def write_preview(path, centers, pocket_radius, wet_radius, box_length, occupancies, seed):
    """Write a fixed-scale orthographic projection, equatorial section and pocket table."""
    pockets = [(i, tuple(map(float, c)), int(n))
               for i, (c, n) in enumerate(zip(centers, occupancies), 1)]
    svg = ET.Element("svg", xmlns="http://www.w3.org/2000/svg",
                     width="1180", height="620", viewBox="0 0 1180 620",
                     role="img", attrib={"aria-labelledby": "title desc"})
    _element(svg, "title", id="title").text = "Initial W₁ pockets inside the PLGA/DCM droplet"
    _element(svg, "desc", id="desc").text = (
        "Left: all pocket footprints projected onto x–y. Middle: the exact section at z = 0. "
        "Right: signed z coordinate, shortest oil film to the outer water and drug count for each pocket.")
    plot = _element(svg, "g", font_family="DejaVu Sans, Arial, sans-serif", fill="#1d3448")
    _element(plot, "rect", width=1180, height=620, fill="white")
    _text(plot, 30, 46, "Initial W₁ pockets inside the PLGA/DCM droplet", 29, font_weight="bold")
    _text(plot, 31, 79, "Wet diameter {:.2f} nm   ·   {} W₁ pockets   ·   {} semaglutide molecules   ·   Seed {}".format(
        2 * wet_radius, len(pockets), sum(n for _, _, n in pockets), seed), 17, fill="#425e70")
    colors = [("#d4eced", "#267d84"), ("#e7def1", "#82699b"), ("#eaf3fa", "#7b9bb1")]
    for x, (fill, stroke), label in zip((31, 243, 462), colors,
                                        ("W₁  Internal water", "O  PLGA in DCM", "W₂  Outer water")):
        _element(plot, "rect", x=x, y=102, width=21, height=15, fill=fill, stroke=stroke)
        _text(plot, x + 30, 115, label)
    _element(plot, "path", d="M30 137 H1150", stroke="#cdd7df")
    side, top = 313., 229.
    scale = side / box_length
    for lft, projection in ((30, True), (402, False)):
        cx, cy = lft + side / 2, top + side / 2
        title = "a   All W₁ pockets" if projection else "b   Equatorial section"
        caption = "Projection along z onto the x–y plane" if projection else "Only pockets intersecting z = 0"
        _text(plot, lft, 173, title, 21, font_weight="bold")
        _text(plot, lft, 198, caption, 14, fill="#526c7b")
        _element(plot, "rect", x=lft, y=top, width=side, height=side,
                 fill=colors[2][0], stroke=colors[2][1])
        _element(plot, "circle", cx=cx, cy=cy, r=wet_radius * scale,
                 fill=colors[1][0], stroke=colors[1][1], stroke_width=1.6)
        positions = []
        for index, (x, y, z), _ in sorted(pockets, key=lambda p: p[1][2]):
            if not projection and abs(z) >= pocket_radius:
                continue
            radius = pocket_radius if projection else math.sqrt(pocket_radius**2 - z**2)
            px, py = cx + scale * x, cy - scale * y
            _element(plot, "circle", cx=px, cy=py, r=scale * radius,
                     fill=colors[0][0], fill_opacity=.58 if projection else 1,
                     stroke=colors[0][1], stroke_width=1.55)
            positions.append((index, (px, py), scale * radius))
        _pocket_labels(plot, sorted(positions), (lft, top, lft + side, top + side))
    table = _element(plot, "g")
    _text(table, 773, 173, "c   Pocket geometry", 21, font_weight="bold")
    _text(table, 773, 198, "Oil film = shortest W₁-to-W₂ surface gap", 14, fill="#526c7b")
    cols = (790, 887, 1011, 1116)
    for col, label in zip(cols, ("ID", "z (nm)", "Oil film (nm)", "Drug count")):
        _text(table, col, 248, label, 14, font_weight="bold", text_anchor="middle")
    for index, (x, y, z), count in pockets:
        yrow = 254 + index * 27
        if index % 2:
            _element(table, "rect", x=763, y=yrow - 19, width=392, height=27, fill="#f1f5f7")
        gap = wet_radius - math.sqrt(x*x + y*y + z*z) - pocket_radius
        for col, value in zip(cols, (index, "{:.2f}".format(z), "{:.2f}".format(gap), count)):
            _text(table, col, yrow, value, 15, text_anchor="middle")
    _text(plot, 30, 571, "Both maps share one scale. Box width {:.2f} nm. W₁ radius {:.2f} nm.".format(
        box_length, pocket_radius), 15)
    _text(plot, 30, 598, "Projected circles can overlap when pockets are separated in depth.",
          14, fill="#526c7b")
    Path(path).write_bytes(ET.tostring(svg, encoding="utf-8", xml_declaration=True))
