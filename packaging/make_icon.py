#!/usr/bin/env python3
"""Draws MIL-C's mark and application icon.

The suite has one animal: the cow's head, face on, white on a gradient tile. MIL-X is the head
alone on blue; MIL-Q is the head inside a ring on green — a targeted assay. MIL-C is the head
inside a **gate** on violet: the polygon a cytometrist draws around the population that matters,
with its handles at the corners and a few events left outside. The violet is the 405 nm laser.

The head itself is not redrawn here. It comes, shape by shape, from the suite's shared generator
in the MIL-X repository (`milx/tools/make_icon.py`), found beside this one or through $MILX_ROOT,
so the three applications cannot drift apart. The files written are kept, so nothing else needs
that repository.

    python3 packaging/make_icon.py                 # the chosen mark: every size, .icns, .ico, .svg
    python3 packaging/make_icon.py --sheet FILE    # the three studies side by side, for choosing
    python3 packaging/make_icon.py --variant plot  # another study as the icon

Variants: `gate` (the mark), `stream` (cells in single file around the head), `plot` (the head
on a pair of axes).
"""
import argparse
import importlib.util
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "icons"

FLAT = "#4F3A9E"                 # the mark on paper, as #234B8C is MIL-X's and #10685C MIL-Q's
TOP = (122, 98, 216)             # #7a62d8, the tile, lit
BOTTOM = (55, 37, 122)           # #37257a
OUTSIDE = (201, 189, 242)        # #c9bdf2, the events outside the gate
WHITE = (255, 255, 255)
SS = 4

SHRINK = 0.70                    # the head inside its gate (0.74 inside MIL-Q's ring)
STROKE = 2.7                     # of the gate, in grid units (the ring's is 3)
HANDLE = 5.4                     # a vertex handle's side

# The gate, in grid units from the tile's centre: six vertices, drawn the way a hand draws one —
# no two sides alike, nothing parallel to the tile.
GATE = [(-23, -17), (6, -25), (27, -10), (25, 14), (3, 26), (-25, 17)]
# Events that fell outside it: (x, y, radius)
EVENTS = [(-33, -30, 1.9), (-27, -36, 1.4), (-36, -20, 1.3), (30, -29, 1.7), (36, -22, 1.3),
          (23, -35, 1.3), (35, 24, 1.8), (31, 32, 1.4), (38, 14, 1.2), (-33, 30, 1.6),
          (-26, 36, 1.3), (-37, 21, 1.2), (14, 36, 1.3), (-8, -36, 1.2)]
STREAM_N, STREAM_R, STREAM_DOT = 18, 29.5, 2.5
AXES = (-27, 27, 3.0)            # where the two axes meet (x = y mirrored), and their stroke


def shared():
    """The suite's generator: the head's shapes and the tile's geometry."""
    candidates = []
    if os.environ.get("MILX_ROOT"):
        candidates.append(Path(os.environ["MILX_ROOT"]) / "tools" / "make_icon.py")
    for sibling in ("OpenDIAL/milx", "mil-x/milx", "milx"):
        candidates.append(ROOT.parent / sibling / "tools" / "make_icon.py")
    for c in candidates:
        if c.is_file():
            spec = importlib.util.spec_from_file_location("suite_icon", c)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    sys.exit("the suite's icon generator was not found; set MILX_ROOT to the milx folder of the "
             "MIL-X repository (looked in: " + ", ".join(map(str, candidates)) + ")")


# ---- raster ---------------------------------------------------------------------------------

def _decorations(draw, variant, cx, cy, unit, colour=255, outside=255, small=False):
    """What stands around the head, on a mask or on a colour layer. `small` is the icon at 32
    pixels and under: the gate alone and a little heavier, because handles and stray events
    two pixels wide are noise, not drawing."""
    def P(x, y):
        return (cx + x * unit, cy + y * unit)

    def dot(x, y, r, fill):
        px, py = P(x, y)
        draw.ellipse((px - r * unit, py - r * unit, px + r * unit, py + r * unit), fill=fill)

    if variant == "gate":
        pts = [P(*v) for v in GATE]
        stroke = STROKE * (1.3 if small else 1.0)
        draw.line(pts + pts[:2], fill=colour, width=max(1, round(stroke * unit)), joint="curve")
        if small:
            return
        for x, y, r in EVENTS:
            dot(x, y, r, outside)
        h = HANDLE * unit / 2
        for px, py in pts:
            draw.rounded_rectangle((px - h, py - h, px + h, py + h), radius=h * 0.28, fill=colour)
    elif variant == "stream":
        for k in range(STREAM_N):
            a = 2 * math.pi * k / STREAM_N - math.pi / 2
            dot(STREAM_R * math.cos(a), STREAM_R * math.sin(a), STREAM_DOT, colour)
    elif variant == "plot":
        lo, hi, w = AXES
        for x, y, r in EVENTS:
            if x > lo + 4 and y < hi - 4:
                dot(x * 0.8, y * 0.8, r, outside)
        width = max(1, round(w * unit))
        draw.line([P(lo, -hi), P(lo, hi), P(hi, hi)], fill=colour, width=width, joint="curve")
        for end in (P(lo, -hi), P(hi, hi)):
            draw.ellipse((end[0] - width / 2, end[1] - width / 2,
                          end[0] + width / 2, end[1] + width / 2), fill=colour)


def render(size, variant="gate"):
    """The icon: the mark in white on the violet tile."""
    g = shared()
    s = size * SS
    img = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    box, side, radius, *_ = g.geometry(s)

    grad = Image.new("RGBA", (1, side))
    px = grad.load()
    for y in range(side):
        t = y / max(1, side - 1)
        px[0, y] = tuple(int(TOP[i] + (BOTTOM[i] - TOP[i]) * t) for i in range(3)) + (255,)
    grad = grad.resize((side, side), Image.BILINEAR)
    tile = Image.new("L", (s, s), 0)
    ImageDraw.Draw(tile).rounded_rectangle(box, radius=radius, fill=255)
    img.paste(grad, (box[0], box[1]), tile.crop(box))

    clear = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    hi = clear.copy()
    ImageDraw.Draw(hi).rounded_rectangle((box[0], box[1], box[2], box[1] + int(side * 0.42)),
                                         radius=radius, fill=(255, 255, 255, 26))
    hi = hi.filter(ImageFilter.GaussianBlur(side * 0.05))
    img.alpha_composite(Image.composite(hi, clear, tile))

    unit = side * 0.66 / (g.MARK_BOX[2] - g.MARK_BOX[0])
    cx, cy = box[0] + side / 2, box[1] + side / 2
    if variant == "plot":
        cx, cy = cx + 2.5 * unit, cy - 2.5 * unit     # the head sits in the plot, off the axes

    # the events first and fainter, then the gate and the head over them
    faint = Image.new("L", (s, s), 0)
    solid = Image.new("L", (s, s), 0)
    _decorations(ImageDraw.Draw(faint), variant, box[0] + side / 2, box[1] + side / 2, unit,
                 colour=0, outside=150, small=size <= 32)
    _decorations(ImageDraw.Draw(solid), variant, box[0] + side / 2, box[1] + side / 2, unit,
                 colour=255, outside=0, small=size <= 32)
    head_box = (box[0] + (cx - box[0] - side / 2), box[1] + (cy - box[1] - side / 2))
    head = g.cow_mask(s, (head_box[0], head_box[1], head_box[0] + side, head_box[1] + side),
                      side, shrink=SHRINK)
    solid = Image.composite(Image.new("L", (s, s), 255), solid, head)
    for mask, colour in ((faint, OUTSIDE), (solid, WHITE)):
        layer = Image.new("RGBA", (s, s), colour + (255,))
        layer.putalpha(Image.composite(mask, Image.new("L", (s, s), 0), tile))
        img.alpha_composite(layer)
    return img.resize((size, size), Image.LANCZOS)


# ---- vector ---------------------------------------------------------------------------------

def _head_svg(g, fill, cx=32.0, cy=32.0):
    """The head on the 64 grid, shrunk about (cx, cy): the same markup as the suite's drawings."""
    P, E, M = g.PATCH, g.EYES, g.MUZZLE
    return (
        '<defs><clipPath id="body"><path d="%s"/></clipPath>' % g.HEAD
        + '<mask id="knock" maskUnits="userSpaceOnUse" x="0" y="0" width="64" height="64">'
        + '<rect width="64" height="64" fill="#fff"/>'
        + '<ellipse cx="%s" cy="%s" rx="%s" ry="%s" transform="rotate(%s %s %s)" fill="#000" '
          'clip-path="url(#body)"/>' % (P[0], P[1], P[2], P[3], P[4], P[0], P[1])
        + "".join('<circle cx="%s" cy="%s" r="%s" fill="#000"/>' % e for e in E)
        + '<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="#000"/>' % M
        + "".join('<ellipse cx="%s" cy="%s" rx="%s" ry="%s" fill="#fff"/>' % n for n in g.NOSTRILS)
        + "</mask></defs>"
        + '<g transform="translate(%s,%s) scale(%s) translate(-32,-29)">' % (cx, cy, SHRINK)
        + '<g mask="url(#knock)" fill="%s">' % fill
        + "".join('<ellipse cx="%s" cy="%s" rx="%s" ry="%s" transform="rotate(%s %s %s)"/>'
                  % (e[0], e[1], e[2], e[3], e[4], e[0], e[1]) for e in (g.EAR_L, g.EAR_R))
        + '<path d="%s"/><path d="%s"/><path d="%s"/>' % (g.HORN_L, g.HORN_R, g.HEAD)
        + "</g></g>")


def _around_svg(variant, fill, outside, opacity):
    """The gate, the stream or the axes, on the 64 grid centred at (32, 32)."""
    out = []
    if variant == "gate":
        out += ['<circle cx="%.1f" cy="%.1f" r="%s" fill="%s" fill-opacity="%s"/>'
                % (32 + x, 32 + y, r, outside, opacity) for x, y, r in EVENTS]
        pts = " ".join("%.1f,%.1f" % (32 + x, 32 + y) for x, y in GATE)
        out.append('<polygon points="%s" fill="none" stroke="%s" stroke-width="%s" '
                   'stroke-linejoin="round"/>' % (pts, fill, STROKE))
        h = HANDLE / 2
        out += ['<rect x="%.2f" y="%.2f" width="%s" height="%s" rx="%.2f" fill="%s"/>'
                % (32 + x - h, 32 + y - h, HANDLE, HANDLE, h * 0.28, fill) for x, y in GATE]
    elif variant == "stream":
        for k in range(STREAM_N):
            a = 2 * math.pi * k / STREAM_N - math.pi / 2
            out.append('<circle cx="%.2f" cy="%.2f" r="%s" fill="%s"/>'
                       % (32 + STREAM_R * math.cos(a), 32 + STREAM_R * math.sin(a), STREAM_DOT, fill))
    elif variant == "plot":
        lo, hi, w = AXES
        out += ['<circle cx="%.1f" cy="%.1f" r="%s" fill="%s" fill-opacity="%s"/>'
                % (32 + x * 0.8, 32 + y * 0.8, r, outside, opacity)
                for x, y, r in EVENTS if x > lo + 4 and y < hi - 4]
        out.append('<polyline points="%s,%s %s,%s %s,%s" fill="none" stroke="%s" stroke-width="%s" '
                   'stroke-linecap="round" stroke-linejoin="round"/>'
                   % (32 + lo, 32 - hi, 32 + lo, 32 + hi, 32 + hi, 32 + hi, fill, w))
    return "".join(out)


def mark_svg(variant="gate", colour=FLAT):
    """The mark alone on the 64 grid, as `mil-x-front.svg` and `mil-q-front.svg` are. The events
    outside the gate reach past the grid, so the view is a little wider than 64."""
    g = shared()
    head = (34.5, 29.5) if variant == "plot" else (32, 32)
    view = "-8 -8 80 80" if variant in ("gate", "plot") else "0 0 64 64"
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="%s" width="320" height="320" '
            'role="img">\n' % view
            + _around_svg(variant, colour, colour, 0.45) + "\n"
            + _head_svg(g, colour, *head) + "\n</svg>\n")


def tile_svg(variant="gate", size=1024):
    g = shared()
    box, side, radius, *_ = g.geometry(size)
    unit = side * 0.66 / (g.MARK_BOX[2] - g.MARK_BOX[0])
    ox, oy = box[0] + side / 2 - 32 * unit, box[1] + side / 2 - 32 * unit
    head = (34.5, 29.5) if variant == "plot" else (32, 32)
    return "\n".join([
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d">'
        % (size, size, size, size),
        '<defs><linearGradient id="g" x1="0" y1="0" x2="0" y2="1"><stop offset="0" '
        'stop-color="#7a62d8"/><stop offset="1" stop-color="#37257a"/></linearGradient>',
        '<linearGradient id="hi" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#fff" '
        'stop-opacity="0.12"/><stop offset="1" stop-color="#fff" stop-opacity="0"/></linearGradient>',
        '<clipPath id="c"><rect x="%d" y="%d" width="%d" height="%d" rx="%d"/></clipPath></defs>'
        % (box[0], box[1], side, side, radius),
        '<rect x="%d" y="%d" width="%d" height="%d" rx="%d" fill="url(#g)"/>'
        % (box[0], box[1], side, side, radius),
        '<rect x="%d" y="%d" width="%d" height="%d" rx="%d" fill="url(#hi)"/>'
        % (box[0], box[1], side, side // 2, radius),
        '<g clip-path="url(#c)"><g transform="translate(%.2f %.2f) scale(%.4f)">' % (ox, oy, unit),
        _around_svg(variant, "#ffffff", "#c9bdf2", 0.6),
        _head_svg(g, "#ffffff", *head),
        "</g></g></svg>"])


# ---- files ------------------------------------------------------------------------------------

def export(variant="gate", name="MIL-C"):
    g = shared()
    OUT.mkdir(parents=True, exist_ok=True)
    master = render(1024, variant)
    (OUT / "png").mkdir(exist_ok=True)
    for size in g.PNG_SIZES:
        # the small ones are drawn at their own size, not shrunk from the large: a gate one
        # pixel wide is decided by the drawing, not by a resampling filter
        (master if size == 1024 else render(size, variant)).save(OUT / "png" / f"{name}-{size}.png")
    iconset = OUT / f"{name}.iconset"
    shutil.rmtree(iconset, ignore_errors=True)
    iconset.mkdir()
    for size, label in g.ICONSET:
        shutil.copy(OUT / "png" / f"{name}-{size}.png", iconset / f"icon_{label}.png")
    if shutil.which("iconutil"):
        subprocess.run(["iconutil", "--convert", "icns", str(iconset),
                        "--output", str(OUT / f"{name}.icns")], check=True)
    master.save(OUT / f"{name}.ico", sizes=[(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    (OUT / f"{name}.svg").write_text(tile_svg(variant), encoding="utf-8")
    (OUT / "mil-c-front.svg").write_text(mark_svg(variant), encoding="utf-8")
    render(256, variant).save(ROOT / "milc" / "icon.png")
    print(f"[icon] {variant} -> {OUT}: PNGs, {name}.iconset, {name}.icns, {name}.ico, {name}.svg, "
          "mil-c-front.svg; milc/icon.png")


def _font(size):
    for path in ("/System/Library/Fonts/Helvetica.ttc", "/System/Library/Fonts/SFNS.ttf"):
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def sheet(path):
    """The three studies, each large and at the sizes an icon is really seen at, beside the two
    marks they have to sit next to."""
    variants = [("gate", "A · gate"), ("stream", "B · stream"), ("plot", "C · plot")]
    W, H = 1500, 1010
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    big, small = _font(30), _font(20)
    for k, (v, label) in enumerate(variants):
        x = 60 + k * 480
        img.paste(render(400, v), (x, 50), render(400, v))
        d.text((x + 200, 470), label, fill="#56617a", font=big, anchor="mm")
        cursor = x + 40
        for size in (128, 64, 32, 16):
            icon = render(size, v)
            img.paste(icon, (cursor, 520 + (128 - size) // 2), icon)
            d.text((cursor + size / 2, 670), str(size), fill="#9aa1ac", font=small, anchor="mm")
            cursor += size + 26
    # the family: X, Q and the chosen C in a row, as in a Dock
    g = shared()
    d.text((60, 730), "the suite", fill="#56617a", font=big)
    row = [g.render(180, "milx"), g.render(180, "milq"), render(180, "gate")]
    for k, icon in enumerate(row):
        img.paste(icon, (60 + k * 190, 780), icon)
    for k, name in enumerate(("MIL-X", "MIL-Q", "MIL-C")):
        d.text((150 + k * 190, 975), name, fill="#9aa1ac", font=small, anchor="mm")
    img.save(path)
    print(f"[icon] sheet -> {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=["gate", "stream", "plot"], default="gate")
    ap.add_argument("--sheet", metavar="FILE", help="the studies side by side, and nothing else")
    a = ap.parse_args()
    if a.sheet:
        sheet(a.sheet)
    else:
        export(a.variant)


if __name__ == "__main__":
    main()
