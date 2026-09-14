# bhmiepy

Python interface to the BHMie Mie scattering code.

bhmiepy wraps the classic BHMIE Mie-scattering subroutine — computing
scattering and absorption by a homogeneous isotropic sphere — in a modern,
NumPy-friendly Python API, compiled directly from Fortran via f2py + Meson.

The wrapped Fortran implementation is the one maintained by Thomas P.
Robitaille in [hyperion-rt/bhmie](https://github.com/hyperion-rt/bhmie)
(BSD-2-Clause), itself derived from the program in Appendix A of
Bohren & Huffman, *Absorption and Scattering of Light by Small Particles*
(Wiley, 1983), with extensive modifications by B. T. Draine.

**Status: early development.**

## Usage (Phase 1: core routine)

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

Validated against the original F77 BHMIE (bundled as a single-precision
oracle in the extension), Rayleigh / geometric limits, and energy
conservation — see `tests/`.

## Usage (Phase 2: dust populations)

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
`dust.read_parameter_file` parses upstream `.in` parameter files directly.

Validated against the upstream reference outputs shipped in
`upstream/examples/` (mrn77 and kmh94_full): all quantities agree to better
than the reference files' 5-significant-digit precision (~5e-5 relative).
Those full-fidelity runs take ~30 s each and live under the `slow` pytest
marker (`pytest -m slow`).

## Installation

From a source checkout (not on PyPI):

```bash
pip install .
```

The build compiles Fortran, so a C and a Fortran compiler must be available
first (e.g. gfortran via your system package manager or conda-forge);
meson-python, ninja, and numpy are pulled in automatically as build
requirements. For development:

```bash
pip install -e . --no-build-isolation   # needs numpy, meson-python, ninja in the env
pytest                                  # fast test suite
pytest -m slow                          # golden-data + upstream-equivalence runs
```

On Windows, a known-good toolchain is conda-forge `flang` **plus
`flang-rt_win-64`** (the Flang runtime is a separate package) and the MSVC
C compiler — see CLAUDE.md for the recipe and toolchain notes.

Development from a git checkout also needs the upstream reference
submodule (refractive-index tables and reference outputs for the tests,
and sources for the optional `bhmie_ref` build):

```bash
git submodule update --init --recursive
```

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
  the `bhmie/` submodule by meson (`-Dupstream_ref`, auto-enabled when
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
| `upstream/` | Git submodule of upstream hyperion-rt/bhmie pinned at `37c0729` (with its fortranlib submodule) — read-only reference & test fixtures; `git submodule update --init --recursive` after clone |
| `tests/` | pytest suite |

## Credits and licenses

- Python wrapper code: MIT © 2026 Tao Jing — see [LICENSE](LICENSE).
- Vendored Fortran from hyperion-rt/bhmie: BSD-2-Clause © 2012-13
  Thomas P. Robitaille — see
  [LICENSES/BSD-2-Clause-bhmie.txt](LICENSES/BSD-2-Clause-bhmie.txt).
- Underlying BHMIE subroutine: Bohren & Huffman (1983), Appendix A,
  modified by B. T. Draine.

See [NOTICE](NOTICE) for the full attribution chain.

If you use bhmiepy in published research, please cite Bohren & Huffman
(1983) and acknowledge the BHMIE implementation by T. P. Robitaille
(https://github.com/hyperion-rt/bhmie).
