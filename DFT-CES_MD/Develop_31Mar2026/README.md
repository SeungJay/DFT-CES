# DFT-CES patch for LAMMPS (develop, 30 Mar 2026) — MACE ML-IAP compatible

DFT-CES (grand canonical density functional theory in classical explicit solvent,
QM/MM electrostatic coupling) patch, ported to the LAMMPS **develop** branch
(`patch_30Mar2026`, commit `f0e9d2a9f8`) and adapted for use alongside the **MACE
ML-IAP** interface running under KOKKOS/CUDA.

This is a re-based and unit-adapted version of the original
`stable_2Aug2023_update1` DFT-CES patch.

## What this patch adds

New commands / classes:
- `grid` command (`grid.cpp`, `grid.h`) — reads a QM electrostatic potential cube file
  (Gaussian `.cube`, Rydberg units) onto a grid, and optionally writes the MD charge
  density back out.
- `fix gridforce` (`fix_gridforce.cpp`, `fix_gridforce.h`) — applies the force from the
  QM potential grid onto the (partial) charges of a group via trilinear interpolation,
  as a `POST_FORCE` fix.

Modified core files (hooks for the above): `domain.cpp/.h`, `input.cpp/.h`, `run.cpp`.

## Differences from the original `stable_2Aug2023_update1` patch

| Change | Reason |
|--------|--------|
| Re-based onto develop `patch_30Mar2026` line numbers | source tree moved since 2Aug2023 stable |
| **`units metal`** in `fix_gridforce`: dropped the `23.06092` (eV→kcal/mol) factor, now `weight*1.0` | original assumed `real` units (kcal/mol); MACE ML-IAP runs in `metal` (eV·Å). With the grid stored in eV and coordinates in Å, `q·V` is already eV and the force is already eV/Å — no conversion needed. |
| `abs(weight)` → **`fabs(weight)`** (8 sites, `fix_gridforce.cpp` charge spreading) | `abs` can bind to the integer `std::abs`, silently truncating a non-integer `weight`; `fabs` is unambiguously double. (Harmless when `weight` is integer-valued, e.g. `-1`, but safe for fractional weights.) |

## Applying the patch

From the LAMMPS source root (develop branch, ~`30 Mar 2026`):

```bash
patch -p1 < stable_2Aug2023_update1_ext_lammps30Mar2026.patch

# verify
ls src/grid.cpp src/grid.h src/fix_gridforce.cpp src/fix_gridforce.h
ls src/*.rej 2>/dev/null && echo "REJECTS!" || echo "clean"
grep -c "fabs(weight)" src/fix_gridforce.cpp   # 8
grep -c "weight\*1.0"  src/fix_gridforce.cpp   # 4
```

`grid.cpp`/`fix_gridforce.cpp` are plain C++ (CPU); they compile under the KOKKOS/CUDA
build via `nvcc_wrapper` as host code and are picked up automatically by the CMake source
GLOB (no extra `-D PKG_...` needed for DFT-CES itself).

## Usage (metal units)

`fix gridforce` needs per-atom charges, so use a charge-bearing atom style:

```
units           metal
atom_style      charge
read_data       data.lammps          # atoms as: id type q x y z
# ... set charges (from QM) ...

pair_style      mliap unified <model>.model-mliap_lammps.pt
pair_coeff      * * <elements>

group           mob type <...>
grid            ./V_ryd.cube yes 1 MDrho.cube     # QM potential (Ry) in; MD rho out
fix             gf mob gridforce -1 1 0           # weight sfactor cubeID
fix_modify      gf energy yes                     # add grid energy to total PE
fix             int all nvt temp 300 300 0.1      # REQUIRED: gridforce does not integrate
run             ...
```

`grid` command arguments: `grid <cube_in> <save yes|no> <ncubes> [out1 ...]`.
`fix gridforce` arguments: `fix ID group gridforce <weight> <sfactor> <cubeID>`
(6 args total; `weight = -1` follows the QM sign convention, `cubeID` in 0–4).

### Coupling with MACE ML-IAP

MACE (`pair_style mliap unified`) provides the interatomic PES; `fix gridforce` adds the
external QM-potential force on the charges. The two contributions are additive
(`F_total = F_MACE + F_grid` at a given configuration), verified to hold under the
single-GPU KOKKOS run (`-k on g 1 -sf kk`). `fix gridforce` runs on the host; LAMMPS
synchronizes the force array between device and host for the non-KOKKOS fix, so the grid
force is correctly included in the GPU integration.

## Notes

- `fix gridforce` is a **host (CPU)** fix. For a KOKKOS/CUDA run it is not ported to a
  device kernel; it relies on LAMMPS' automatic host↔device force synchronization. This is
  fine for typical system sizes (the grid interpolation is cheap relative to the MACE
  evaluation), but has not been performance-tuned.
- Charges must be physically meaningful (QM-derived partial charges) for the coupling to be
  correct — the grid force is `∝ q`.
- Applied to LAMMPS develop `f0e9d2a9f8` (2026.3.31); line numbers are matched to that
  snapshot. On a different develop commit the hunks may need re-basing.
