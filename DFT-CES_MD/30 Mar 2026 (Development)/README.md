# DFT-CES patch for LAMMPS (30 Mar 2026)

A patch and usage notes for applying **DFT-CES** functionality to a recent LAMMPS source.
The original DFT-CES patch targets LAMMPS `2 Aug 2023`, so it does not apply cleanly to the
latest source. This repository provides a re-based patch whose hunk positions and context
have been adjusted to match the newer LAMMPS version.

## Background

| Item | Detail |
|------|--------|
| LAMMPS source | Obtained on 2026-07-02 via `git clone https://github.com/lammps/lammps.git` |
| Cloned version | `30 Mar 2026 (Development)` (`src/version.h`) |
| Original patch | `stable_2Aug2023_update1_ext.patch` from [SeungJay/DFT-CES](https://github.com/SeungJay/DFT-CES/tree/main/DFT-CES_MD/stable_2Aug2023_update1) |
| Original patch base | LAMMPS `2 Aug 2023 (Update 1)` |
| New patch | `stable_2Aug2023_update1_ext_lammps30Mar2026.patch` |

The original patch (2023) and the cloned LAMMPS source (2026) differ, so the context and
line numbers no longer match. By comparing the 2023 reference tree
(`lammps-stable_2Aug2023_update1`) against the 2026 source, the patch was regenerated to fit
the 2026 version (done with the help of Claude Cowork).

## Files changed by the patch

**Modified existing files (5)**

- `src/domain.cpp`, `src/domain.h`
- `src/input.cpp`, `src/input.h`
- `src/run.cpp`

**Newly added files (4, DFT-CES core)**

- `src/grid.cpp`, `src/grid.h`
- `src/fix_gridforce.cpp`, `src/fix_gridforce.h`

## How to apply

> ⚠️ Use GNU `patch`. `git apply` may fail because it is strict about blank context lines.

```bash
# 1. Move to the LAMMPS root
cd /path/to/lammps

# 2. Dry-run to check whether it applies (does not modify any files)
patch -p1 --dry-run < ../stable_2Aug2023_update1_ext_lammps30Mar2026.patch

# 3. Apply for real
patch -p1 < ../stable_2Aug2023_update1_ext_lammps30Mar2026.patch
```

`-p1` strips the leading path component (`a/`, `b/`) and applies the changes to `src/...`.
If the dry-run reports all 9 files as `checking file ...` with no `FAILED` / `offset` / `fuzz`
warnings, everything is fine.

To revert:

```bash
patch -p1 -R < ../stable_2Aug2023_update1_ext_lammps30Mar2026.patch
```

## Build

The newly added `grid.cpp` and `fix_gridforce.cpp` are picked up automatically because both
make and cmake scan `src/*.cpp`; no manual registration is needed.

**traditional make**

```bash
cd src
make yes-<required-packages>
make mpi          # or: make serial
```

**cmake**

```bash
mkdir build && cd build
cmake ../cmake
cmake --build . -j
```

## Notes

- The new patch was verified against the cloned 2026 source with `patch -p1 --dry-run`
  (all files apply cleanly).
- If you later pull a newer LAMMPS version, the source may change again, so the patch may
  need to be re-based the same way (by comparing against the reference tree).
