#!/usr/bin/env bash
# The macOS half of a release, from the Mac, for when the CI cannot build it
# (it does, on a tag or on request: this is the way round it).
#
#     bash packaging/release_mac.sh            # tests and disk image, in dist/
#     bash packaging/release_mac.sh v1.0.0     # the same, then attached to the
#                                              # release v1.0.0 (or build-N)
#
# 1. the whole suite, offscreen, as the CI ran it on macOS
# 2. make_dmg.sh: the application, its own self-test, the window opened and
#    closed, the disk image
# 3. with a TAG: the image attached to that release, and SHA256SUMS made
#    again over everything the release serves, the image included
#
# The build goes outside the repository (iCloud Drive syncs this folder, and
# a build is thousands of files); only the image lands in dist/.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
tag="${1:-}"
build="${MILC_BUILD_DIR:-${TMPDIR:-/tmp}/milc-macbuild}"
cd "$root"

if [ -n "$tag" ] && ! gh release view "$tag" >/dev/null 2>&1; then
    echo "there is no release $tag on GitHub: push the tag first, or run without one" >&2
    exit 1
fi

echo "==> the suite"
QT_QPA_PLATFORM=offscreen python3 -m pytest -q -p no:cacheprovider

echo "==> the disk image"
bash "$here/make_dmg.sh" "$build"
version="$(python3 -c "import milc; print(milc.__version__)")"
dmg="$root/dist/MIL-C-$version-macos-$(uname -m).dmg"
[ -f "$dmg" ] || { echo "no $dmg was made" >&2; exit 1; }

if [ -z "$tag" ]; then
    echo "$dmg"
    exit 0
fi

echo "==> attached to $tag"
gh release upload "$tag" "$dmg" --clobber

echo "==> SHA256SUMS over what $tag serves"
sums="$(mktemp -d)"
gh release download "$tag" --dir "$sums" \
    --pattern '*.dmg' --pattern '*.msi' --pattern '*.zip' \
    --pattern '*.deb' --pattern '*.tar.gz'
(cd "$sums" && shasum -a 256 MIL-C-* > SHA256SUMS && cat SHA256SUMS)
gh release upload "$tag" "$sums/SHA256SUMS" --clobber
rm -rf "$sums"
