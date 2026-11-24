# PyInstaller build, shared by macOS, Windows and Linux.
#
# What is not obvious, and fails late — the bundle builds, and then the first
# UMAP does not run:
#
#   * numba compiles functions the first time they are called and keeps the
#     result on disk, keyed by the source file. A function marked cache=True
#     whose source is only inside the archive has nowhere to be keyed from.
#     So milc itself, umap and pynndescent are collected as .py files as well
#     as inside the archive; numba then keeps what it compiled in a cache
#     folder of the user's, and the second start is three times faster.
#   * matplotlib is used through three back ends chosen by name at run time
#     (Agg, SVG, PDF) and needs its own data folder for fonts.
#
# Build:  pyinstaller packaging/milc.spec --noconfirm

import re
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_submodules

MACOS = sys.platform == "darwin"
WINDOWS = sys.platform.startswith("win")

ROOT = Path(SPECPATH).parent
ICONS = Path(SPECPATH) / "icons"
# Linux has no icon inside an executable: the menu entry names one
ICON = (str(ICONS / "MIL-C.ico") if WINDOWS else
        str(ICONS / "MIL-C.icns") if MACOS else None)

# read rather than repeat: a version written twice disagrees with itself
VERSION = re.search(r'__version__ = "([^"]+)"',
                    (ROOT / "milc" / "__init__.py").read_text()).group(1)

# the whole tree goes, so the translations (help/pages/pt) travel with the
# English pages
datas = [(str(ROOT / "milc" / "icon.png"), "milc"),
         (str(ROOT / "milc" / "help" / "pages"), "milc/help/pages"),
         (str(ROOT / "LICENSE"), ".")]

hidden = [
    "matplotlib.backends.backend_agg",
    "matplotlib.backends.backend_svg",
    "matplotlib.backends.backend_pdf",
    "PyQt6.QtSvg",
    "sklearn.utils._typedefs",
    "sklearn.neighbors._partition_nodes",
    # openTSNE's compiled modules import one another from C, where
    # PyInstaller cannot see it: the bundle built and tSNE failed on
    # openTSNE._matrix_mul
    *collect_submodules("openTSNE"),
]

# nothing here reaches these; named so that one cannot creep in unnoticed
excluded = [
    "tkinter", "IPython", "jupyter", "notebook", "pytest", "pandas",
    "matplotlib.backends.backend_qtagg", "matplotlib.backends.backend_tkagg",
    "PyQt6.QtWebEngineCore", "PyQt6.QtWebEngineWidgets", "PyQt6.Qt3DCore",
    "PyQt6.QtBluetooth", "PyQt6.QtMultimedia", "PyQt6.QtQuick", "PyQt6.QtQml",
    "PyQt6.QtPdf", "PyQt6.QtSql", "PyQt6.QtDesigner",
]

analysis = Analysis(
    [str(ROOT / "run.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=hidden,
    excludes=excluded,
    module_collection_mode={
        "milc": "pyz+py",
        "umap": "pyz+py",
        "pynndescent": "pyz+py",
    },
    noarchive=False,
)

pyz = PYZ(analysis.pure)

executable = EXE(
    pyz,
    analysis.scripts,
    [],
    exclude_binaries=True,
    name="MIL-C",
    console=False,
    # an unsigned arm64 binary will not start at all; ad-hoc is enough
    codesign_identity=None,
    entitlements_file=None,
    icon=ICON,
)

collection = COLLECT(
    executable,
    analysis.binaries,
    analysis.datas,
    strip=False,
    upx=False,
    name="MIL-C",
)

if MACOS:
    app = BUNDLE(
        collection,
        name="MIL-C.app",
        icon=ICON,
        bundle_identifier="io.github.pedronobrez.milc",
        info_plist={
            "CFBundleName": "MIL-C",
            "CFBundleDisplayName": "MIL-C",
            "CFBundleShortVersionString": VERSION,
            "CFBundleVersion": VERSION,
            "NSHighResolutionCapable": True,
            # the window follows the system light or dark setting
            "NSRequiresAquaSystemAppearance": False,
            "NSHumanReadableCopyright": "GPL-3.0-or-later",
            "LSApplicationCategoryType": "public.app-category.medical",
        },
    )
