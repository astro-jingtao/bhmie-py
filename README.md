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

**bhmiepy is not on PyPI yet** — install from a source checkout. Besides
Python >= 3.10 you need a C and a Fortran compiler; numpy, meson-python,
and ninja are pulled in automatically as build requirements, so a plain
`pip install .` is all it takes once the compilers are set up.

### Linux

```bash
sudo apt install gfortran        # or your distribution's equivalent
pip install .
```

### macOS

```bash
brew install gcc                 # provides gfortran
pip install .
```

### Windows

Two toolchains are known to work; both are exercised in CI on every
push. Pick either one — set the environment up first, then run
`pip install .` from the same shell.

#### Option A — conda-forge flang + MSVC

1. Install the [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
   (or Visual Studio itself) with the *Desktop development with C++*
   workload. This provides the C compiler (`cl`) and the linker.

2. Create a conda environment with Fortran and the build tools (any
   prompt works for this step):

   ```
   conda create -n bhmiepy -c conda-forge python numpy meson-python ninja pip flang flang-rt_win-64
   ```

   `flang` is the LLVM Fortran compiler and `flang-rt_win-64` its
   runtime — the runtime is a separate package, and without it linking
   fails with `LNK1104: flang_rt.runtime.*.lib`.

3. Open an **x64 Native Tools Command Prompt for VS** from the Start
   menu (or run `vcvarsall.bat x64` in a `cmd` prompt), and activate
   the environment *in that prompt*:

   ```
   conda activate bhmiepy
   ```

4. Set the extra variables the conda flang packaging needs — one line
   each in the same prompt (typed interactively; inside a `.cmd`
   script write `%%i`. In PowerShell use `$env:AR = "llvm-ar"`-style
   assignments instead):

   ```
   set AR=llvm-ar
   set RANLIB=llvm-ranlib
   for /f %i in ('flang-new --print-resource-dir') do set "LIB=%LIB%;%i\lib\x86_64-pc-windows-msvc"
   set FFLAGS=-static-libflangrt
   ```

   Why each one: meson searches for `ar`/`ranlib`, but the conda
   environment only ships the `llvm-*` names; the MSVC linker does not
   know where the flang runtime library lives, so its directory is
   appended to `LIB`; and conda's flang ships no runtime DLL, so the
   Fortran code must be compiled against the static runtime
   (`FFLAGS`).

5. Build and install:

   ```
   pip install .
   ```

#### Option B — mingw-w64 gfortran

1. Install [MSYS2](https://www.msys2.org/), then the UCRT64 Fortran
   toolchain from an MSYS2 shell (or any mingw-w64 distribution such
   as [WinLibs](https://winlibs.com/)):

   ```
   pacman -S mingw-w64-ucrt-x86_64-gcc-fortran
   ```

2. In a `cmd` prompt with the toolchain's `bin` directory on `PATH`
   (for MSYS2 typically `set PATH=C:\msys64\ucrt64\bin;%PATH%`), set:

   ```
   set CC=gcc
   set FC=gfortran
   set FFLAGS=-static-libgfortran -static-libgcc
   set LDFLAGS=-static-libgfortran -static-libgcc -static-libquadmath -Wl,-Bstatic -lgfortran -lquadmath -lwinpthread -Wl,-Bdynamic
   ```

   The `LDFLAGS` line statically links the Fortran runtime into the
   extension and is **not optional**: since Python 3.8 the interpreter
   no longer resolves a module's DLL dependencies from `PATH`, so a
   dynamically linked `libgfortran-5.dll` makes `import bhmiepy` fail
   with an unhelpful "DLL load failed" error.

3. Build and install:

   ```
   pip install .
   ```

### Development

```bash
git clone --recurse-submodules https://github.com/astro-jingtao/bhmie-py
cd bhmie-py
pip install -e . --no-build-isolation     # needs numpy, meson-python, ninja in the env
pytest                                    # fast test suite
pytest -m slow                            # golden-data + upstream-equivalence runs
```

On Windows, build with either option's toolchain and environment
settings above (add `pytest` to the environment). The `upstream/`
submodule provides refractive-index tables, recorded reference outputs,
and the sources for the optional `bhmie_ref` executable (built
automatically when the submodule is present; `-Dupstream_ref=disabled`
skips it). A plain `pip install .` also works without the submodule —
the reference build is simply skipped — but two fast parameter-file
tests read fixtures from `upstream/examples/`, so run
`git submodule update --init --recursive` before `pytest` (the already
installed package does not need rebuilding for it).

The `upstream/` submodule provides refractive-index tables, recorded
reference outputs, and the sources for the optional `bhmie_ref`
executable (built automatically when the submodule is present;
`-Dupstream_ref=disabled` skips it). A plain `pip install .` also works
without the submodule — the reference build is simply skipped — but two
fast parameter-file tests read fixtures from `upstream/examples/`, so
run the one-line submodule command above before `pytest` (the already
installed package does not need rebuilding for it).

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
  × numpy 1.26/2.x), with Windows covered under both documented
  toolchains (conda flang + MSVC, and mingw-w64 gfortran).

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
program with modifications by B. T. Draine. Alongside these, an
acknowledgment of bhmiepy itself is also warmly welcomed.

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
