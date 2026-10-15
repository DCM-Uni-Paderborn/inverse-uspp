# Inverse USPP reference data

Original ultrasoft pseudopotentials supplemented with projector-indexed reference waves for inverse Gaussian representations in CP2K / Quickstep.

The first collection follows the SSSP v2.0 PBE and PBEsol Efficiency and Precision selections. Additional scalar PSLibrary and GBRV potentials extend element coverage. Files outside SSSP are identified separately in the catalog.

## Reconstruction routes

- Route 3 repeats the source generator and accepts reference waves when the original dataset is reproduced within the recorded compatibility bound.
- Route 4 defines an auxiliary reconstruction from complete exported waves when the source replay differs. The original local potential, projectors, nonlocal coefficients, augmentation and core density remain unchanged. These are not claimed to be recovered historical original waves.

An invertible reconstruction is required. All density, metric, kinetic-energy and geometry corrections needed by the original functional must be retained. Complete valence pairs do not supply a missing NLCC model-core kinetic density. Dataset checks establish data integrity and reconstruction algebra, not basis accuracy or complete CP2K feature coverage.

The machine-readable catalog records original and supplemented hashes, generator input, reconstruction route, compatibility and conditioning. Generator builds and intermediate output are not part of this repository.

## Sources

- [SSSP](https://www.materialscloud.org/discover/sssp)
- [PSLibrary](https://github.com/dalcorso/pslibrary)
- [GBRV](https://www.physics.rutgers.edu/gbrv/)
- [Vanderbilt generator](https://www.physics.rutgers.edu/~dhv/uspp/)

Original author notices and licenses remain in the datasets. PSLibrary and GBRV are distributed under the GNU GPL. The preparation tools are GPL-2.0-or-later.
