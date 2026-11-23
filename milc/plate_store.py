"""
Plates kept with the project: what `plate.py` holds in memory, written so
that it is there again when the folder is opened.

A panel is stored as what it shows — its content, arrays and all — and
not as a way to make it again: the page of a paper must not change because
the analysis behind it was run once more. So a panel read back has no
plot behind it (`source`) and is not stale: it is what it is.

The words and the small numbers go in the project's JSON under `plates`;
the arrays go beside it in one
`.npz`, each under a name made from its own bytes, so an array that
several panels share (the same map coloured four ways, a plate and its
copy) is written once, and the file is written again only when the arrays
are no longer the same ones.

Nothing here knows the kinds of content. A content is a dataclass, and a
dataclass is stored as its fields: numbers, text, lists, tuples, dicts,
sets, enums, numpy arrays and other dataclasses, however nested. A kind
of content added elsewhere is stored without a word here. What cannot be
stored, or cannot be made again when it is read — its class is gone, its
arrays are missing — becomes a `Missing` panel that says so on the page
and keeps what was stored, so saving again loses nothing.

No Qt here.
"""

from __future__ import annotations

import dataclasses
import enum
import hashlib
import importlib
import math
import os
import sys
import textwrap
import zipfile
from dataclasses import dataclass, field

import numpy as np

from . import figures
from .figures import FigureStyle
from .plate import Panel, Plate

#: the arrays of the panels, beside the project's file
FILE = "MIL-C plates.npz"
VERSION = 1

#: an array of at most this many numbers is written in the JSON itself: a
#: contour is hundreds of short lines, and a zip of a thousand members is
#: slow to write and to read
INLINE = 64

_TAG = "$"
_PACKAGE = __name__.rsplit(".", 1)[0]


class Unstorable(TypeError):
    pass


@dataclass
class Missing:
    """In the place of a content that could not be stored or restored."""
    note: str
    #: what the project held for it, written back as it was
    stored: object = None
    arrays: dict = field(default_factory=dict)
    uses_colormap = False

    def category_names(self) -> list[tuple[str, str]]:
        return []


def _draw_missing(fig, c: Missing, style: FigureStyle) -> None:
    ax = fig.subplots()
    ax.set_axis_off()
    width = max(int(fig.bbox.width / fig.dpi * 72 / (style.font_size * 0.65)), 12)
    text = "This panel could not be restored.\n" + "\n".join(textwrap.wrap(c.note, width))
    ax.text(0.5, 0.5, text, ha="center", va="center", color="#8a8f98",
            fontsize=style.font_size, transform=ax.transAxes)


figures._DRAW[Missing] = _draw_missing


# --- any value → what JSON can hold ---------------------------------------------------

def _key(a: np.ndarray) -> str:
    h = hashlib.sha1(f"{a.dtype.str}{a.shape}".encode())
    h.update(np.ascontiguousarray(a).data)
    return h.hexdigest()[:20]


def _type_name(cls) -> str:
    return f"{cls.__module__}:{cls.__qualname__}"


def _array(a: np.ndarray, arrays: dict):
    if a.dtype.kind == "O":
        return {_TAG: "objects", "shape": list(a.shape),
                "v": [encode(x, arrays) for x in a.ravel().tolist()]}
    if a.dtype.kind not in "biufcU":
        raise Unstorable(f"an array of {a.dtype}")
    plain = a.dtype.kind in "biuU" or (a.dtype.kind == "f" and bool(np.isfinite(a).all()))
    if a.size <= INLINE and plain:
        return {_TAG: "array", "dtype": a.dtype.str, "shape": list(a.shape),
                "v": a.ravel().tolist()}
    key = _key(a)
    arrays[key] = a
    return {_TAG: "npz", "key": key}


def encode(value, arrays: dict):
    """`value` as JSON can hold it; the arrays too large for that are put
    in `arrays` under the names they are referred to by. `Unstorable`
    when something in it is none of the kinds this knows."""
    if value is None or isinstance(value, (bool, str)):
        return str(value) if isinstance(value, str) else value
    if isinstance(value, enum.Enum):
        return {_TAG: "enum", "type": _type_name(type(value)), "v": value.name}
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        value = float(value)
        return value if math.isfinite(value) else {_TAG: "float", "v": repr(value)}
    if isinstance(value, np.ndarray):
        return _array(value, arrays)
    if isinstance(value, np.generic):
        if value.dtype.kind not in "biuf":
            raise Unstorable(f"a {value.dtype}")
        return {_TAG: "scalar", "dtype": value.dtype.str, "v": encode(value.item(), arrays)}
    if isinstance(value, Missing):
        arrays.update(value.arrays)
        return value.stored
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {_TAG: "data", "type": _type_name(type(value)),
                "v": {f.name: encode(getattr(value, f.name), arrays)
                      for f in dataclasses.fields(value)}}
    if isinstance(value, list):
        return [encode(x, arrays) for x in value]
    if isinstance(value, (tuple, set, frozenset)):
        kind = "tuple" if isinstance(value, tuple) else type(value).__name__
        return {_TAG: kind, "v": [encode(x, arrays) for x in value]}
    if isinstance(value, dict):
        if all(isinstance(k, str) for k in value) and _TAG not in value:
            return {k: encode(v, arrays) for k, v in value.items()}
        return {_TAG: "map", "v": [[encode(k, arrays), encode(v, arrays)]
                                   for k, v in value.items()]}
    raise Unstorable(f"a {type(value).__name__}")


def _class(name: str, what):
    """The class stored as `module:Name`: one of this program's, or one
    from a module that has a kind of content among those that are drawn.
    Nothing else is imported, or made, on a file's word."""
    module, _, qualname = name.partition(":")
    ours = module == _PACKAGE or module.startswith(_PACKAGE + ".")
    if ours and module not in sys.modules:
        importlib.import_module(module)
    if not ours and module not in {c.__module__ for c in figures._DRAW}:
        raise LookupError(f"{name} is not a kind of content this program knows")
    obj = sys.modules[module]
    for part in qualname.split("."):
        obj = getattr(obj, part)
    if not (isinstance(obj, type) and what(obj)):
        raise LookupError(f"{name} is not what was stored")
    return obj


def decode(value, arrays):
    """What `encode` was given. `arrays` is anything indexed by name. An
    exception — of any kind — when it cannot be made again."""
    if isinstance(value, list):
        return [decode(x, arrays) for x in value]
    if not isinstance(value, dict):
        return value
    if _TAG not in value:
        return {k: decode(v, arrays) for k, v in value.items()}
    kind, v = value[_TAG], value.get("v")
    if kind == "npz":
        return arrays[value["key"]]
    if kind == "array":
        return np.array(v, dtype=np.dtype(value["dtype"])).reshape(value["shape"])
    if kind == "objects":
        out = np.empty(len(v), dtype=object)
        out[:] = [decode(x, arrays) for x in v]
        return out.reshape(value["shape"])
    if kind == "scalar":
        return np.dtype(value["dtype"]).type(decode(v, arrays))
    if kind == "float":
        return float(v)
    if kind in ("tuple", "set", "frozenset"):
        return {"tuple": tuple, "set": set, "frozenset": frozenset}[kind](
            decode(x, arrays) for x in v)
    if kind == "map":
        return {decode(k, arrays): decode(x, arrays) for k, x in v}
    if kind == "enum":
        return _class(value["type"], lambda c: issubclass(c, enum.Enum))[v]
    if kind == "data":
        cls = _class(value["type"], dataclasses.is_dataclass)
        given = {k: decode(x, arrays) for k, x in v.items()}
        known = {f.name: f for f in dataclasses.fields(cls)}
        obj = cls(**{k: x for k, x in given.items() if k in known and known[k].init})
        for k, x in given.items():
            if k in known and not known[k].init:
                object.__setattr__(obj, k, x)
        return obj
    raise ValueError(f"an unknown kind of value: {kind!r}")


# --- plates ↔ what the project holds ---------------------------------------------------

def _panel_record(panel: Panel, arrays: dict) -> dict:
    record = {"title": panel.title, "name": panel.name, "style": {},
              "colspan": panel.colspan, "rowspan": panel.rowspan,
              "kind": type(panel.content).__name__}
    own: dict = {}
    try:
        if isinstance(panel.content, Missing):
            stored = panel.content.stored or {}
            record.update({k: v for k, v in stored.items() if v is not None},
                          content=stored.get("content"))
            own = dict(panel.content.arrays)
        else:
            record["content"] = encode(panel.content, own)
        record["style"] = encode(panel.style, own)
    except Exception as exc:  # a panel that cannot be kept must not cost the others
        record.update(content=None, style={}, problem=f"{record['kind']} could not be "
                      f"stored: {exc}")
        own = {}
    arrays.update(own)
    return record


def dump(plates: list[Plate], current: int = 0) -> tuple[dict, dict]:
    """(what goes in the project's JSON, the arrays it refers to by name)."""
    arrays: dict = {}
    record = {"version": VERSION, "file": FILE, "current": int(current), "plates": [
        {"name": p.name, "columns": p.columns, "width_mm": p.width_mm,
         "height_mm": p.height_mm, "letters": p.letters, "style": p.style.to_dict(),
         "panels": [_panel_record(panel, arrays) for panel in p.panels]}
        for p in plates]}
    return record, arrays


def _panel(record: dict, arrays) -> Panel:
    stored = {k: record.get(k) for k in ("kind", "content", "problem")}
    kind = record.get("kind") or "panel"
    try:
        if record.get("content") is None:
            raise ValueError(record.get("problem") or f"{kind} was not stored")
        content = decode(record["content"], arrays)
        if type(content) not in figures._DRAW:
            raise LookupError(f"nothing here draws a {type(content).__name__}")
    except Exception as exc:
        if isinstance(exc, KeyError):
            note = f"The data of this {kind} are not in “{FILE}”."
        else:
            note = str(exc) if record.get("content") is None else \
                f"{kind} could not be read: {exc}"
        refs = _references(record.get("content"))
        content = Missing(note, stored, {k: arrays[k] for k in refs if k in arrays})
    try:
        style = decode(record.get("style"), arrays)
    except Exception:
        style = None
    return Panel(content, str(record.get("title") or ""),
                 dict(style) if isinstance(style, dict) else {},
                 max(int(record.get("colspan") or 1), 1),
                 max(int(record.get("rowspan") or 1), 1),
                 name=str(record.get("name") or ""))


def _references(value) -> set[str]:
    """The names of the arrays a stored value points to."""
    if isinstance(value, list):
        return set().union(*map(_references, value)) if value else set()
    if not isinstance(value, dict):
        return set()
    if value.get(_TAG) == "npz":
        return {value.get("key")}
    return set().union(*map(_references, value.values())) if value else set()


def restore(record, arrays) -> tuple[list[Plate], int]:
    """The plates of `dump`, and which was the current one. Whatever is
    wrong with a panel, the panel is there, as a `Missing`; whatever is
    wrong with a plate, the others are there."""
    plates = []
    if not isinstance(record, dict) or not isinstance(record.get("plates"), list):
        return plates, 0
    for k, r in enumerate(record["plates"]):
        try:
            height = r.get("height_mm")
            p = Plate(columns=max(int(r.get("columns") or 2), 1),
                      width_mm=float(r.get("width_mm") or 180.0),
                      height_mm=float(height) if height else None,
                      style=FigureStyle.from_dict(r.get("style")),
                      letters=bool(r.get("letters", True)),
                      name=str(r.get("name") or f"Plate {k + 1}"))
            for panel in r.get("panels") or []:
                try:
                    p.panels.append(_panel(panel, arrays))
                except Exception as exc:
                    p.panels.append(Panel(Missing(f"It could not be read: {exc}")))
            plates.append(p)
        except Exception:
            continue
    try:
        current = int(record.get("current") or 0)
    except (TypeError, ValueError):
        current = 0
    return plates, min(max(current, 0), max(len(plates) - 1, 0))


# --- the files ---------------------------------------------------------------------------

def _stored_keys(path: str) -> set[str] | None:
    try:
        with np.load(path, allow_pickle=False) as held:
            return set(held.files)
    except Exception:
        return None


def _write_arrays(path: str, arrays: dict) -> None:
    """An `.npz` as numpy reads it, whole or not at all. Packed lightly:
    the lightest packing is within a few percent of the heaviest on what
    a panel holds, in a fifth of the time, and this is written while the
    user waits."""
    partial = path + ".part"
    with zipfile.ZipFile(partial, "w", zipfile.ZIP_DEFLATED, compresslevel=1) as z:
        for key, a in arrays.items():
            with z.open(key + ".npy", "w", force_zip64=True) as member:
                np.lib.format.write_array(member, np.ascontiguousarray(a), allow_pickle=False)
    os.replace(partial, path)


def write(folder: str, plates: list[Plate], current: int = 0) -> dict | None:
    """Write the arrays beside the project (only when they are not the
    ones already there) and return what the project's JSON is to hold
    under `plates`; None, and no file, when there are no plates."""
    path = os.path.join(folder, FILE)
    if not plates:
        try:
            os.remove(path)
        except OSError:
            pass
        return None
    record, arrays = dump(plates, current)
    if set(arrays) != (_stored_keys(path) or set()) or (not arrays and os.path.exists(path)):
        if arrays:
            _write_arrays(path, arrays)
        else:
            os.remove(path)
    return record


class _Arrays:
    """The arrays of the side file, each read when first asked for."""

    def __init__(self, path: str):
        self._held, self._read = None, {}
        try:
            self._held = np.load(path, allow_pickle=False)
        except Exception:
            pass

    def __contains__(self, key) -> bool:
        return key in self._read or (self._held is not None and key in self._held.files)

    def __getitem__(self, key):
        if key not in self._read:
            if key not in self:
                raise KeyError(key)
            self._read[key] = self._held[key]
        return self._read[key]

    def close(self) -> None:
        if self._held is not None:
            self._held.close()


def read(folder: str, record) -> tuple[list[Plate], int]:
    """The plates a project holds (`record`: its `plates`), with which is
    the current one. Never an error: a project without plates has none."""
    if not record:
        return [], 0
    arrays = _Arrays(os.path.join(folder, FILE))
    try:
        return restore(record, arrays)
    except Exception:
        return [], 0
    finally:
        arrays.close()
