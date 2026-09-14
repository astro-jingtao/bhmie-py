# bhmiepy

**bhmiepy is a Python interface to the classic BHMIE Mie-scattering code,
based on the Fortran implementation maintained by Thomas P. Robitaille in
[hyperion-rt/bhmie](https://github.com/hyperion-rt/bhmie)** (BSD-2-Clause),
itself derived from the program in Appendix A of Bohren & Huffman,
*Absorption and Scattering of Light by Small Particles* (Wiley, 1983),
with extensive modifications by B. T. Draine. The core numerical
implementation is upstream's, not original work; bhmiepy contributes the
compiled Python binding, a NumPy-friendly vectorized API, and the
size-distribution-averaging pipeline as a library.

It computes scattering and absorption by homogeneous isotropic spheres
(efficiencies, asymmetry parameter, and amplitude scattering functions)
and, building on that, dust-population properties for one or more grain
materials with size distributions: cross-sections, extinction opacity,
albedo, and scattering-matrix elements S11/S12/S33/S34. The Fortran core
is compiled directly via f2py + Meson — no Python-level loops over
parameters.

**Status: early development** (not on PyPI yet; install from source).

## Installation

From a source checkout:

```bash
pip install .
```

Requirements: Python >= 3.10, numpy >= 1.26, and a C and a Fortran
compiler. The build uses meson-python + f2py; meson-python, ninja, and
numpy are pulled in automatically as build requirements.

On Linux and macOS, gfortran from your system package manager
(`sudo apt install gfortran` / `brew install gcc`) is sufficient.

On **Windows**, a known-good toolchain is conda-forge `flang` **plus
`flang-rt_win-64`** (the Flang runtime is a separate package) together
with the MSVC C compiler (e.g. from Visual Studio Build Tools):

```
conda create -n bhmiepy -c conda-forge python numpy meson-python ninja pytest pip flang flang-rt_win-64
```

Build from an environment that has run `vcvarsall.bat x64` and activated
the conda env, with the following additions (the conda flang packaging
needs them; see the project's CLAUDE.md for the rationale):
`AR=llvm-ar`, `RANLIB=llvm-ranlib`, the flang resource directory's
`lib\x86_64-pc-windows-msvc` appended to `LIB`, and
`FFLAGS=-static-libflangrt`.

### Development

```bash
git submodule update --init --recursive   # test fixtures + optional reference build
pip install -e . --no-build-isolation     # needs numpy, meson-python, ninja in the env
pytest                                    # fast test suite
pytest -m slow                            # golden-data + upstream-equivalence runs
```

The `upstream/` submodule provides refractive-index tables, recorded
reference outputs, and the sources for the optional `bhmie_ref`
executable (built automatically when the submodule is present;
`-Dupstream_ref=disabled` skips it). A plain `pip install .` without
the submodule works and runs the fast test suite.

## Usage — single spheres

```python
import bhmiepy

result = bhmiepy.bhmie(x=2.5, m=1.6 + 0.01j, nang=90)
result.qext, result.qsca, result.qback, result.g   # efficiencies
result.qabs, result.albedo                         # derived
result.s1, result.s2, result.angles                # amplitudes on 0..pi

# vectorized: x and/or m as arrays, computed in Fortran (no Python loop)
import numpy as np
result = bhmiepy.bhmie(np.logspace(-2, 2, 100), 1.6 + 0.01j, nang=90)
result.qext.shape                                  # (100,)

# physical quantities convenience
result = bhmiepy.compute(radius=0.1, wavelength=0.25, refractive_index=1.6 + 0.01j)
```

`x` is the size parameter `2*pi*a/lambda` and `m` the relative
refractive index (`Im(m) >= 0`); `nang` angles between 0 and 90 degrees
yield `2*nang - 1` amplitudes covering 0 to 180 degrees.

## Usage — dust populations

Size-distribution-averaged dust properties — the full capability of the
upstream CLI tool, as a library call (units follow upstream: micron,
g/cm^3, cm^2):

```python
import numpy as np
from bhmiepy import dust

silicate = dust.Material.from_file("upstream/examples/ri-data/silicate_ld93")
comp = dust.Component(
    material=silicate,
    distribution=dust.PowerLaw(0.005, 0.25, -3.5),  # n(a) ~ a^-3.5
    abundance_mass=0.627,
    density=3.3,
)
res = dust.compute_dust_properties(
    [comp],
    wavelengths=np.logspace(np.log10(0.01), np.log10(1000), 250),
    amin=0.005, amax=1.0, na=1000,
    n_angles=181, n_small_angles=10,
    gas_to_dust=141.84,
)
res.kappa_ext   # extinction opacity, cm^2/g (incl. gas)
res.albedo, res.g, res.s11, res.s12, res.s33, res.s34, res.angles
```

`PowerLaw`, `PeakedPowerLaw` (upstream 'ped'), and `TableDistribution`
(upstream 'table', also `TableDistribution.from_file`) are supported, and
`dust.read_parameter_file` parses upstream `.in` parameter files
directly — an existing upstream model can be recomputed as a library
call.

## Validation

- Cross-checked against the original F77 BHMIE (Bohren & Huffman 1983
  App. A / Draine), vendored as a single-precision oracle inside the
  extension, over fixed cases and random samples.
- Analytic limits: Rayleigh regime, non-absorbing spheres
  (`qsca = qext`), geometric limit, absorption non-negativity.
- The dust pipeline reproduces the upstream reference outputs shipped in
  `upstream/examples/` (mrn77 and kmh94_full): all quantities agree to
  better than the reference files' 5-significant-digit precision
  (~5e-5 relative). Those full-fidelity runs take ~30 s each and live
  under the `slow` pytest marker.
- CI runs the full test matrix (Ubuntu/macOS/Windows × Python 3.10–3.13
  × numpy 1.26/2.x).

### Benchmarks vs the pristine upstream code

Two independent A/B benchmarks (both separate from the pytest suite):

- `benchmarks/bench_upstream.py` — core routine: the package's optimized
  BHMIE copy vs a verbatim upstream copy compiled into the same extension.
  On a dust-shaped workload (2400 size x wavelength points) the results
  are bitwise identical and the package copy shows no measurable
  performance degradation; any speedup is compiler-dependent — parity
  with gfortran, ~2-3x with LLVM flang (whose codegen is more sensitive
  to the upstream copy's fixed 16 MB per-call workspace).
- `benchmarks/bench_dust_upstream.py` — dust pipeline (Phase 2): the
  library vs the upstream CLI itself (`bhmie_ref`), built verbatim from
  the `upstream/` submodule by meson (`-Dupstream_ref`, auto-enabled when
  the submodule is initialized). Same parameter files in, outputs
  compared at the upstream output precision, end-to-end wall-clock
  reported.

```
python benchmarks/bench_upstream.py
python benchmarks/bench_dust_upstream.py
```

## Repository layout

| Path | Contents |
| --- | --- |
| `src/bhmiepy/` | Python wrapper (MIT) |
| `src/bhmiepy/_fortran/bhmie.f90`, `types.f90` | Derived from upstream bhmie (BSD-2-Clause); see file headers |
| `src/bhmiepy/_fortran/bhmiepy_ext.f90` | New Fortran code written for bhmiepy (MIT) |
| `src/bhmiepy/_fortran/ref/` | Comparison-only sources: original F77 oracle + pristine upstream copy |
| `upstream/` | Git submodule of upstream hyperion-rt/bhmie pinned at `37c0729` (with its fortranlib submodule) — read-only reference & test fixtures |
| `tests/` | pytest suite (fast by default, `-m slow` for golden-data runs) |
| `benchmarks/` | A/B benchmarks against the pristine upstream code |
| `CHANGELOG.md` | Release notes |

## Citing

If you use bhmiepy in published research, please cite the underlying
method and implementation:

```bibtex
@book{bohren1983absorption,
  author    = {Bohren, Craig F. and Huffman, Donald R.},
  title     = {Absorption and Scattering of Light by Small Particles},
  publisher = {Wiley},
  address   = {New York},
  year      = {1983},
  isbn      = {978-0-471-29340-8},
  doi       = {10.1002/9783527618156}
}
```

and acknowledge the BHMIE implementation by Thomas P. Robitaille
(https://github.com/hyperion-rt/bhmie), which extends the original
program with modifications by B. T. Draine. An acknowledgment of
bhmiepy itself is welcome but not required.

## Credits and licenses

- Python wrapper code: MIT © 2026 Tao Jing — see [LICENSE](LICENSE).
- Vendored Fortran from hyperion-rt/bhmie: BSD-2-Clause © 2012-13
  Thomas P. Robitaille — see
  [LICENSES/BSD-2-Clause-bhmie.txt](LICENSES/BSD-2-Clause-bhmie.txt).
- Log-log interpolation/integration algorithms ported from
  astrofrog/fortranlib: BSD-2-Clause © 2009-13 Thomas P. Robitaille —
  see [LICENSES/BSD-2-Clause-fortranlib.txt](LICENSES/BSD-2-Clause-fortranlib.txt).
- Underlying BHMIE subroutine: Bohren & Huffman (1983), Appendix A,
  modified by B. T. Draine.

See [NOTICE](NOTICE) for the full attribution chain.
