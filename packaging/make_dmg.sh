#!/usr/bin/env bash
# Build MIL-C.app and wrap it in a disk image.
#
#     bash packaging/make_dmg.sh [BUILD_DIR]
#
# BUILD_DIR (default: dist/ here) is where the application is built. Give one
# outside a synchronised folder — iCloud Drive, Dropbox — when the repository
# lives in one: a build is thousands of files, and a sync client reading them
# while they are written is slow at best. The image itself always lands in
# dist/ of the repository.
#
# The image is unsigned. That is not what makes macOS complain: Gatekeeper
# reacts to the com.apple.quarantine attribute, which is attached by whatever
# downloads a file, not by whatever builds it. An image built here and mounted
# here opens with no dialog at all. One that travels — mail, AirDrop, a
# browser — will ask, and the reader can right-click Open once, or run
#     xattr -dr com.apple.quarantine /Applications/MIL-C.app
#
# Apple Silicon does insist on *some* signature or the binary will not start;
# the ad-hoc one PyInstaller applies is enough and costs nothing.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
build="${1:-$root/dist}"
version="$(python3 -c "import re,pathlib;print(re.search(r'__version__ = \"([^\"]+)\"', pathlib.Path('$root/milc/__init__.py').read_text()).group(1))")"
mkdir -p "$root/dist" "$build"
dmg="$root/dist/MIL-C-$version-macos-$(uname -m).dmg"
app="$build/MIL-C.app"

echo "==> building the application in $build"
rm -rf "$app" "$build/MIL-C"
python3 -m PyInstaller "$here/milc.spec" --noconfirm --log-level WARN \
    --distpath "$build" --workpath "$build/.build"
[ -d "$app" ] || { echo "no $app was produced" >&2; exit 1; }

echo "==> checking it runs before shipping it"
# from an empty folder, so nothing is imported from the repository by accident
( cd "$(mktemp -d)" && "$app/Contents/MacOS/MIL-C" selftest )
# the self-test builds the window itself; this is the way a double click
# takes, start to finish, and the two have differed
"$app/Contents/MacOS/MIL-C" gui --quit-after 3
echo "    the window opens and closes"

echo "==> staging"
stage="$build/.dmg"
rm -rf "$stage" "$dmg"
mkdir -p "$stage"
cp -R "$app" "$stage/"
ln -s /Applications "$stage/Applications"

echo "==> creating $dmg"
hdiutil create -volname "MIL-C $version" -srcfolder "$stage" \
    -ov -format UDZO "$dmg" >/dev/null
rm -rf "$stage"

echo
echo "$dmg"
du -h "$dmg" | cut -f1 | xargs echo "size:"
du -sh "$app" | cut -f1 | xargs echo "application:"
