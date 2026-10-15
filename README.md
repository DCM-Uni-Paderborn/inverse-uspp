# Inverse USPP reference data

Original ultrasoft pseudopotentials supplemented with projector-indexed reference waves for inverse Gaussian representations in CP2K / Quickstep.

The collection contains 194 scalar UPFs for 95 elements from H through Am, with PBE and PBEsol data for every element. It includes all 86 distinct USPPs for 41 elements in the four SSSP v2.0 Efficiency and Precision collections. Another 108 PSLibrary potentials cover the remaining 54 elements. Those additional potentials are not part of the SSSP selection or its validation.

The UPFs occupy about 433 MB. Original potential data are stored once per distinct dataset. The catalog records membership in the SSSP collections.

## Reconstruction routes

- Route 3 covers 118 files. It repeats the source generator and accepts reference waves when the original dataset is reproduced within the recorded table-scaled tolerance of `1e-6`. The radial mesh and quadrature must agree exactly. This includes ten SSSP potentials and the 108 newly generated additional potentials.
- Route 4 covers 76 files. It defines an auxiliary reconstruction from complete exported waves when a compatible historical replay is not established. This includes 24 PSLibrary and 52 GBRV SSSP potentials. These waves are not claimed to be recovered historical original waves.

For both routes, the original local potential, projectors, nonlocal coefficients, augmentation and core density remain byte for byte unchanged. Only `has_wfc` is enabled and a `PP_FULL_WFC` block is added. A comment and the catalog identify each supplement.

An inverse calculation uses `E_inv[Psi] = E_UPF[T^-1 Psi]`. All density, normalization-metric, kinetic-energy and geometry corrections needed by the original functional must be retained. Invertibility alone does not justify the reduced `INVERSE_1C` identities or an all-electron interpretation of auxiliary waves.

The Vanderbilt export transforms reference waves into the final projector basis, then maps them onto the unchanged original UPF mesh. The catalog records the projector alignment and common outer continuation. For Re PBEsol, the generator exported the reference pairs and subsequently stopped at its original q-eigenvalue diagnostic. That file is explicitly an auxiliary reconstruction, not a successful historical potential replay.

## Use and scope

- `data/pbe/sssp/` and `data/pbesol/sssp/` contain the supplemented SSSP selections.
- `data/pbe/additional/` and `data/pbesol/additional/` contain the additional PSLibrary potentials.
- `catalog.json` records hashes, valence charge, source selection, reconstruction route and atomic conditioning.
- `sources.json` pins source archives, generator versions and the scope of the native CP2K checks.
- `inputs/gbrv/` contains the selected atomic generator inputs. QE generator inputs are embedded in the UPFs.

All 194 files passed the CP2K development reader and native atomic inverse initialization. The original fields were compared exactly, and atomic kinetic-energy and metric corrections were finite. Separate radial inverse round trips also passed. These are data and atomic-algebra checks. They do not establish multi-atom invertibility, SCF convergence, force/stress correctness or basis accuracy for every element. Use a CP2K version implementing the corresponding USPP representation and a suitable Gaussian basis.

No model-core kinetic-energy density was added. Complete valence reference pairs do not resolve the separate NLCC core-tau requirement of the chosen meta-GGA extension. PBE and PBEsol do not require that additional quantity.

## Preparation tools

`tools/prepare_upf.py` requires Python, NumPy and SciPy. It preserves and checks the original operator sections and refuses to overwrite an existing output.

```sh
python tools/prepare_upf.py original.upf generated.upf supplemented.upf
python tools/prepare_upf.py original.upf reference-waves.dat supplemented.upf --format vanderbilt
```

`generate_pslibrary.py` runs baseline and full-wave atomic calculations for additional elements. `generate_gbrv.py` runs the selected Vanderbilt inputs. `patch_vanderbilt.py` adds the reference export to a separate generator source copy. With recent GNU Fortran, link `date_compat.f90` when the legacy date routine is unavailable. Generator binaries, build trees and raw verification output are excluded.

## Sources

- [SSSP](https://www.materialscloud.org/discover/sssp)
- [PSLibrary](https://github.com/dalcorso/pslibrary)
- [GBRV](https://www.physics.rutgers.edu/gbrv/)
- [Vanderbilt generator](https://www.physics.rutgers.edu/~dhv/uspp/)

Original author notices and licenses remain in the datasets. PSLibrary and GBRV are distributed under the GNU GPL. The preparation tools are GPL-2.0-or-later.
