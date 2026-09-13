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

**Status: early development.** The roadmap lives in
[docs/DEV_PLAN.md](docs/DEV_PLAN.md).

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

Phase 2 will add the dust-population layer of the upstream code:
size-distribution averaging, refractive-index tables, scattering matrices,
and extinction opacities.

## Installation (development, Windows)

Requires a conda env with the Fortran/C toolchain (see CLAUDE.md for the
exact recipe: conda-forge `flang` + `flang-rt_win-64`, meson-python, ninja)
and Visual Studio Build Tools for the C compiler.

```powershell
cmd /c scripts\dev-install.cmd
conda run -n bhmiepy python -m pytest
```

## Repository layout

| Path | Contents |
| --- | --- |
| `src/bhmiepy/` | Python wrapper (MIT) |
| `src/bhmiepy/_fortran/bhmie.f90`, `types.f90` | Derived from upstream bhmie (BSD-2-Clause); see file headers |
| `src/bhmiepy/_fortran/bhmiepy_ext.f90` | New Fortran code written for bhmiepy (MIT) |
| `bhmie/` | Pristine snapshot of upstream hyperion-rt/bhmie @ `37c0729` — read-only reference & test fixtures |
| `tests/` | pytest suite |
| `docs/` | Development plan and design docs |

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
