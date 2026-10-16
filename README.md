# Inverse USPP reference data

Original ultrasoft pseudopotentials supplemented with projector-indexed reference waves for inverse Gaussian representations in CP2K / Quickstep.

The collection contains 389 scalar UPFs for 95 elements from H through Am, with PBE and PBEsol data for every element and LDA data for 64 elements. It includes all 86 distinct USPPs for 41 elements in the four SSSP v2.0 Efficiency and Precision collections, 108 additional PSLibrary potentials, and the complete current GBRV USPP collection for PBE, PBEsol and LDA. Additional potentials are not part of the SSSP selection or its validation.

The UPFs occupy about 683 MB. Files shared by several SSSP selections are stored once, with membership recorded in the catalog. The full GBRV collection retains the upstream datasets and radial representations separately from the SSSP selection.

The additional PSLibrary files were generated to cover the 54 elements missing from the SSSP USPP selection. They are not the complete PSLibrary with all core partitions, cutoff variants and functionals.

## How the reference waves are obtained

- **Compatible generator replay**, recorded as route 3. We rerun the atomic pseudopotential generator with full reference-wave output enabled. The regenerated potential must reproduce the original tables within the recorded table-scaled tolerance of `1e-6`, with exactly the same radial mesh and quadrature. The resulting AE/pseudo reference pairs supplement that original potential.
- **Auxiliary reconstruction**, recorded as route 4. Complete exported reference pairs define an invertible reconstruction operator `T`, even when a compatible replay of the historical potential has not been established. We retain the original potential and evaluate its functional with the back-transformed soft orbitals. The reference waves define the change of orbital variables. They do not replace the original potential or claim to reproduce its historical AE reference waves.

The numbers are catalog identifiers. They do not indicate accuracy levels. The collection contains 118 compatible replays and 271 auxiliary reconstructions. The latter include all 195 upstream GBRV files.

For both approaches, the original local potential, projectors, nonlocal coefficients, augmentation and core density remain unchanged. In a UPF v2 file, only `has_wfc` is enabled and a `PP_FULL_WFC` block is added. Every other original data section remains byte for byte unchanged. A comment and the catalog identify each supplement.

An inverse calculation uses `E_inv[Psi] = E_UPF[T^-1 Psi]`. All density, normalization-metric, kinetic-energy and geometry corrections needed by the original functional must be retained. Invertibility alone does not justify the reduced `INVERSE_1C` identities or an all-electron interpretation of auxiliary waves.

The Vanderbilt export transforms reference waves into the final projector basis, then maps them onto the unchanged original UPF mesh. The catalog records the projector alignment and common outer continuation. For Re PBEsol, the generator exported the reference pairs and subsequently stopped at its original q-eigenvalue diagnostic. That file is explicitly an auxiliary reconstruction, not a successful historical potential replay.

The full GBRV collection contains 65 files per functional, including the Hf⁴⁺ alternative. The upstream aggregate LDA/PBEsol archives omit Bi, which is obtained from the individually published files. The Na v1.5 generator inputs are also obtained separately. `convert_gbrv.f90` uses QE's own reader, augmentation expansion and writer to convert UPF v1 to v2. Inner augmentation polynomials become angular-resolved radial tables on the original grid. Numeric fields and the expanded augmentation are checked against the v1 data to serialization roundoff. Reference supplementation then preserves every converted operator section exactly. For LDA reference generation, the atomic input selects Ceperley–Alder/Perdew–Zunger through `exfact=0`. The catalog records the Ta reference-input formatting repair and all reference-export diagnostics.

## Use and scope

- `data/<functional>/sssp/` contains the supplemented SSSP selections for PBE and PBEsol.
- `data/<functional>/additional/pslibrary/` contains the additional PSLibrary potentials for PBE and PBEsol.
- `data/<functional>/additional/gbrv/` contains the complete current GBRV UPFs for PBE, PBEsol and LDA.
- `catalog.json` records hashes, valence charge, source selection, reconstruction route and atomic conditioning.
- `sources.json` pins source archives, generator versions and the scope of the native CP2K checks.
- `inputs/gbrv/` contains the selected atomic generator inputs. QE generator inputs are embedded in the UPFs.

All 389 files passed the CP2K development reader and native atomic inverse initialization. The original v2 fields were compared exactly, and atomic kinetic-energy and metric corrections were finite. Separate radial inverse round trips also passed. These are data and atomic-algebra checks. They do not establish multi-atom invertibility, SCF convergence, force/stress correctness or basis accuracy for every element. Use a CP2K version implementing the corresponding USPP representation and a suitable Gaussian basis.

GBRV Ga q19, Pd q16 and Rh q15 require the additional closed-core configurations in the [CP2K development reader](https://github.com/DCM-Uni-Paderborn/cp2k/commit/c168ae38c55d03738a1c2b6607fb0bcdbc36b2ea).

No model-core kinetic-energy density was added. Complete valence reference pairs do not resolve the separate NLCC core-tau requirement of the chosen meta-GGA extension. PBE, PBEsol and LDA do not require that additional quantity.

## Preparation tools

`tools/prepare_upf.py` requires Python, NumPy and SciPy. It preserves and checks the original operator sections and refuses to overwrite an existing output.

```sh
python tools/prepare_upf.py original.upf generated.upf supplemented.upf
python tools/prepare_upf.py original.upf reference-waves.dat supplemented.upf --format vanderbilt
```

`generate_pslibrary.py` runs baseline and full-wave atomic calculations for additional elements. `generate_gbrv.py` runs the selected Vanderbilt inputs. `patch_vanderbilt.py` adds the reference export to a separate generator source copy. With recent GNU Fortran, link `date_compat.f90` when the legacy date routine is unavailable. Generator binaries, build trees and raw verification output are excluded.

`import_gbrv.py` imports the complete upstream PBE/PBEsol/LDA UPF collections into an existing preparation workspace. It checks pinned source hashes, invokes `convert_gbrv`, runs the atomic reference exports and prepares the auxiliary reconstructions. Build `convert_gbrv.f90` against the QE 6.3MaX modules and libraries, using the include and link settings of QE's `upf2upf2.x` target. `--retry-failed` retains successful exports. Native CP2K validation is a separate final check before publication.

## Other functionals and sources

Upstream [PSLibrary](https://github.com/dalcorso/pslibrary) also supplies generator recipes for PZ-LDA, PW91, BP, revPBE and WC, with scalar and fully relativistic variants. Their presence upstream does not mean those additional variants are included here.

Other USPP sources include the [original QE collection](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library) and the [legacy Vanderbilt library](https://www.physics.rutgers.edu/~dhv/uspp/). They require individual checks of pseudopotential type, licensing and available reconstruction data before inclusion.

## Sources

- [SSSP](https://www.materialscloud.org/discover/sssp)
- [PSLibrary](https://github.com/dalcorso/pslibrary)
- [GBRV](https://www.physics.rutgers.edu/gbrv/)
- [Vanderbilt generator](https://www.physics.rutgers.edu/~dhv/uspp/)

Original author notices and licenses remain in the datasets. PSLibrary and GBRV are distributed under the GNU GPL. The preparation tools are GPL-2.0-or-later.
