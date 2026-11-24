---
title: Installation
---

MIL-C is published as installers on the releases page of its repository,
`https://github.com/pedronobrez/mil-c/releases`. Each release has five files, one or two per
system, and a `SHA256SUMS` file to check them against. Every installer carries its own Python,
Qt, numpy and compiled numerical code: nothing else has to be installed first, and none of them
needs R.

## The five installers

For version 1.0.0:

| File | System | What it does |
|---|---|---|
| `MIL-C-1.0.0-macos-arm64.dmg` | macOS on Apple Silicon | a disk image: drag **MIL-C** onto **Applications** |
| `MIL-C-1.0.0-windows-x64.msi` | Windows, 64-bit | installs for the current user, with a Start menu entry |
| `MIL-C-1.0.0-windows-x64.zip` | Windows, 64-bit | the same program as a folder: unpack anywhere and run `MIL-C.exe` |
| `MIL-C-1.0.0-linux-x86_64.deb` | Debian, Ubuntu, Mint | installs to `/opt/MIL-C`, with a menu entry and `milc` on the PATH |
| `MIL-C-1.0.0-linux-x86_64.tar.gz` | any Linux, x86-64 | unpack and run `install.sh`, or run the program where it lies |

Notes on each:

- **macOS.** Apple Silicon only: numba, with which UMAP and FlowSOM are compiled, no longer
  publishes builds for Intel Macs.
- **Windows `.msi`.** Installs into your own `AppData\Local\MIL-C` and asks for no administrator.
  A newer version replaces the older one rather than sitting beside it. Remove it from
  **Settings ▸ Apps** like any other program.
- **Windows `.zip`.** For a machine where an installer is not allowed to run. Nothing is written
  outside the folder you unpack it into, except the program's settings and caches.
- **Linux `.deb`.** `sudo apt install ./MIL-C-1.0.0-linux-x86_64.deb`. The package names the few
  system libraries Qt needs, so apt installs them too. `sudo apt remove mil-c` takes it away.
- **Linux `.tar.gz`.** For any distribution and no administrator:

```
tar -xzf MIL-C-1.0.0-linux-x86_64.tar.gz
MIL-C-1.0.0/install.sh            # into ~/.local, with a menu entry
MIL-C-1.0.0/install.sh --remove   # and out again
```

  `install.sh` puts the program in `~/.local/opt/MIL-C` and `milc` in `~/.local/bin`; it says so
  if that folder is not on your PATH. Without installing, `MIL-C-1.0.0/MIL-C/MIL-C` runs as it is.
- The Linux builds are made on Ubuntu 22.04 and run there or on anything as recent: a bundle
  takes the C library from the machine it runs on.

## The first start

None of the installers is signed, so the first start asks for a confirmation:

- **macOS**: right-click (or Control-click) **MIL-C** in Applications and choose **Open**, then
  **Open** again. If macOS only offers to move it to the bin, open **System Settings ▸ Privacy &
  Security**, find the line about MIL-C and choose **Open Anyway**. This is asked once.
- **Windows**: SmartScreen says it protected your PC. Choose **More info**, then **Run anyway**.

The first UMAP or FlowSOM after installing takes longer than the ones after it, because the
numerical code is compiled once and kept in your cache folder: measured in the macOS bundle, 16 s
for the first UMAP and 5 s for each after it. [[troubleshooting]] has more on first-launch
problems.

## Checking a download

`SHA256SUMS` lists the checksum of every installer, computed over the files the release serves.
To check yours:

```
shasum -a 256 -c SHA256SUMS --ignore-missing        # macOS
sha256sum -c SHA256SUMS --ignore-missing            # Linux
Get-FileHash MIL-C-1.0.0-windows-x64.msi            # Windows PowerShell: compare by eye
```

Every installer is larger than 50 MB; a much smaller file is a broken download.

## Checking the program: selftest

Every installed MIL-C can test itself. The self-test writes four invented FCS files to a
temporary folder, opens them as a session, concatenates them, runs UMAP, tSNE and FlowSOM, writes
a figure as SVG, PDF, PNG and TIFF, and then does the same through the window without showing it.
It prints each step with the time it took and ends with `MIL-C 1.0.0: self-test passed`.

```
/Applications/MIL-C.app/Contents/MacOS/MIL-C selftest      # macOS
milc selftest                                              # Linux, installed
MIL-C.exe selftest --log selftest.log                      # Windows, in the program's folder
```

On Windows the program has no console, so `--log FILE` writes the steps, and any failure, to a
file instead. `--no-window` skips the part that needs Qt. See [[command-line]].

## Running from source

MIL-C is a Python package (`pyproject.toml`), and runs from a copy of the repository on Python
3.10 or newer. It depends on numpy, scipy, scikit-learn, numba, umap-learn, openTSNE, matplotlib,
PyQt6 and pyqtgraph, which pip installs:

```
git clone https://github.com/pedronobrez/mil-c
cd mil-c
pip install -e .
python -m milc                     # the window
python -m milc gui path/to/folder  # the window, with a folder open
```

`pip install -e .` also puts a `milc` command on the PATH, the same as `python -m milc`. On
macOS, double-clicking `MIL-C.command` in the repository starts the window too.

Once it runs, [[opening-data]] is the next page.
