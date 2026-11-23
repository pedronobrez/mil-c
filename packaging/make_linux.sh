#!/usr/bin/env bash
# Wrap the Linux build of MIL-C — the folder PyInstaller leaves in
# dist/MIL-C — in the two forms a Linux user expects.
#
#     python -m PyInstaller packaging/milc.spec --noconfirm
#     bash packaging/make_linux.sh
#
#   MIL-C-<version>-linux-x86_64.deb      Debian, Ubuntu, Mint: installs to
#                                         /opt/MIL-C, with a menu entry and
#                                         `milc` on the PATH
#   MIL-C-<version>-linux-x86_64.tar.gz   any distribution, no administrator:
#                                         unpack, run ./install.sh (into
#                                         ~/.local) or MIL-C/MIL-C as it is
#
# The bundle holds its own Python, Qt and every compiled extension. What it
# takes from the system is the C library — so it is built on the oldest
# Ubuntu the runners offer, and runs on that or anything newer — and the few
# X libraries Qt loads by name, which the package lists as dependencies.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
dist="${1:-$root/dist}"
bundle="$dist/MIL-C"
[ -x "$bundle/MIL-C" ] || { echo "no $bundle/MIL-C: build the application first" >&2; exit 1; }

version="$(python3 -c "import re,pathlib;print(re.search(r'__version__ = \"([^\"]+)\"', pathlib.Path('$root/milc/__init__.py').read_text()).group(1))")"
machine="$(uname -m)"
case "$machine" in x86_64) arch=amd64 ;; aarch64) arch=arm64 ;; *) arch="$machine" ;; esac
name="MIL-C-$version-linux-$machine"

desktop() {   # desktop <Exec> <Icon>
    cat <<ENTRY
[Desktop Entry]
Type=Application
Name=MIL-C
GenericName=Flow cytometry analysis
Comment=Gate, balance, concatenate and cluster flow cytometry data
Exec=$1
Icon=$2
Terminal=false
Categories=Science;Biology;DataVisualization;
StartupWMClass=MIL-C
ENTRY
}

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

echo "==> $name.tar.gz"
stage="$work/tar/MIL-C-$version"
mkdir -p "$stage"
cp -a "$bundle" "$stage/MIL-C"
cp "$root/LICENSE" "$stage/"
cp "$here/icons/png/MIL-C-256.png" "$stage/mil-c.png"
cat > "$stage/install.sh" <<'INSTALL'
#!/usr/bin/env bash
# Install MIL-C for this user, with no administrator: the program into
# ~/.local/opt, `milc` into ~/.local/bin, an entry in the applications menu.
# ./install.sh --remove takes all of it away again.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data="${XDG_DATA_HOME:-$HOME/.local/share}"
opt="$HOME/.local/opt/MIL-C"
bin="$HOME/.local/bin/milc"
entry="$data/applications/mil-c.desktop"
icon="$data/icons/hicolor/256x256/apps/mil-c.png"
if [ "${1:-}" = "--remove" ]; then
    rm -rf "$opt" "$bin" "$entry" "$icon"
    echo "MIL-C was removed."
    exit 0
fi
mkdir -p "$(dirname "$opt")" "$(dirname "$bin")" "$(dirname "$entry")" "$(dirname "$icon")"
rm -rf "$opt"
cp -a "$here/MIL-C" "$opt"
ln -sf "$opt/MIL-C" "$bin"
cp "$here/mil-c.png" "$icon"
cat > "$entry" <<ENTRY
[Desktop Entry]
Type=Application
Name=MIL-C
GenericName=Flow cytometry analysis
Comment=Gate, balance, concatenate and cluster flow cytometry data
Exec=$opt/MIL-C
Icon=mil-c
Terminal=false
Categories=Science;Biology;DataVisualization;
StartupWMClass=MIL-C
ENTRY
command -v update-desktop-database >/dev/null 2>&1 && update-desktop-database "$data/applications" || true
echo "MIL-C is in $opt. Start it from the applications menu, or with: $bin"
case ":$PATH:" in *":$HOME/.local/bin:"*) ;; *) echo "($HOME/.local/bin is not on your PATH.)" ;; esac
INSTALL
chmod +x "$stage/install.sh"
tar -C "$work/tar" -czf "$dist/$name.tar.gz" "MIL-C-$version"

echo "==> $name.deb"
pkg="$work/deb"
mkdir -p "$pkg/DEBIAN" "$pkg/opt" "$pkg/usr/bin" "$pkg/usr/share/applications" \
         "$pkg/usr/share/doc/mil-c"
cp -a "$bundle" "$pkg/opt/MIL-C"
ln -s /opt/MIL-C/MIL-C "$pkg/usr/bin/milc"
desktop /opt/MIL-C/MIL-C mil-c > "$pkg/usr/share/applications/mil-c.desktop"
for size in 16 32 48 64 128 256 512; do
    icons="$pkg/usr/share/icons/hicolor/${size}x${size}/apps"
    mkdir -p "$icons"
    cp "$here/icons/png/MIL-C-$size.png" "$icons/mil-c.png"
done
cp "$root/LICENSE" "$pkg/usr/share/doc/mil-c/copyright"
cat > "$pkg/DEBIAN/control" <<CONTROL
Package: mil-c
Version: $version
Section: science
Priority: optional
Architecture: $arch
Installed-Size: $(du -sk "$pkg" | cut -f1)
Depends: libc6, libegl1, libgl1, libfontconfig1, libdbus-1-3, libglib2.0-0 | libglib2.0-0t64, libxkbcommon-x11-0, libxcb-cursor0, libxcb-icccm4, libxcb-image0, libxcb-keysyms1, libxcb-randr0, libxcb-render-util0, libxcb-shape0, libxcb-xinerama0
Maintainer: Pedro Nobre <pedronobrez@users.noreply.github.com>
Homepage: https://github.com/pedronobrez/mil-c
Description: Flow cytometry analysis: gate, balance, concatenate, cluster
 MIL-C reads FCS files and FlowJo workspaces, takes a balanced number of
 events from each sample, concatenates them and runs UMAP, tSNE, FlowSOM
 and PhenoGraph-style clustering, with statistics between groups and
 figures ready for a page.
CONTROL
find "$pkg" -type d -exec chmod 755 {} +
dpkg-deb --build --root-owner-group -Zxz "$pkg" "$dist/$name.deb" >/dev/null

ls -l "$dist/$name.tar.gz" "$dist/$name.deb"
