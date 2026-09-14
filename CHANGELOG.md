# Changelog

All notable changes to bhmiepy are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [Semantic Versioning](https://semver.org/).

Version policy: the package develops under `0.x.y` until the first PyPI
release (deliberately not published yet — see the project roadmap);
`1.0.0` marks the first published release, after which breaking API
changes bump the major version. The version number lives in three
places — `pyproject.toml`, `meson.build`, and
`src/bhmiepy/__init__.py` — and a test enforces they stay in sync;
bump all three together. Dates below reflect when each version's
development finished; no version has been tagged or published yet, so
0.2.0 is the first publicly visible state of the repository.

## [Unreleased]

## [0.2.0] — 2026-09-15

Phase 2 (dust populations) plus publication-readiness work.

### Added

- Dust-population pipeline (`bhmiepy.dust`): `Material`
  (`from_file` / `interpolate`), size distributions `PowerLaw`,
  `PeakedPowerLaw` (upstream `ped`), and `TableDistribution`
  (`from_file`), `Component`, and `compute_dust_properties()` returning
  a `DustResult` (cross-sections, `kappa_ext`, `g`, scattering-matrix
  elements `s11/s12/s33/s34`, and derived `albedo` / `mu` /
  `polarization_90`). Units and numerical conventions follow the
  upstream CLI (micron, g/cm^3, cm^2).
- `dust.read_parameter_file()` parsing upstream `.in` parameter files
  (power/ped/table distributions, relative-path resolution).
- Golden-data validation against the reference outputs shipped in the
  `upstream/` submodule (mrn77, kmh94_full; agreement to the reference
  files' ~5 significant digits) under the `slow` pytest marker.
- `upstream/` read-only git submodule (hyperion-rt/bhmie, pinned) as
  test-fixture source, plus the optional verbatim `bhmie_ref` CLI
  build (`-Dupstream_ref`) for A/B benchmarking and the live
  upstream-equivalence slow test.
- Benchmarks (independent of pytest): `benchmarks/bench_upstream.py`
  (core routine: bitwise identical to pristine upstream, timings) and
  `benchmarks/bench_dust_upstream.py` (dust pipeline vs the upstream
  CLI end to end).
- CI: GitHub Actions test matrix over Ubuntu/macOS/Windows, Python
  3.10–3.13, and numpy 1.26/2.x.
- This changelog.

### Changed

- Performance: the BHMIE workspace is allocated at the series order
  `d(nmx)` instead of the fixed 16 MB `d(nmxx)` per call, and the
  dust-pipeline accumulation moved into the extension
  (`bhmie_dust_accum`) — results bitwise unchanged, fast test suite
  ~20x faster, dust pipeline ~21% faster end to end. Any speedup over
  the pristine upstream code is compiler-dependent (parity with
  gfortran, ~2–3x with LLVM flang).
- The nang-sized Fortran work arrays are heap-allocated, removing the
  practical upper bound on `nang` (stack overflow beyond ~20000 with
  stack-allocating compilers).

### Fixed

- Process-safety and silent-failure hardening from the pre-publication
  review: series-order and input guards on the dust path (oversized
  size parameters no longer reach the Fortran `stop`), descending and
  zero-integral table distributions rejected loudly instead of
  producing zeros/NaN, NaN table rows rejected at construction,
  non-overlapping size ranges raise instead of returning NaN `g`,
  weight skipping follows upstream's `not (weight > 0)` semantics.
- Empty (zero-element) `x`/`m` input raises a clear `ValueError`
  instead of a cryptic f2py error.
- `Material.interpolate()` rejects non-monotonic wavelength tables
  (previously silent garbage); single-row refractive-index and
  size-distribution files parse correctly (`ndmin=2`).
- The guard against oversized size parameters compares in float64
  (an int64 cast wrapped around beyond x ≈ 9.2e18 and bypassed it).

## [0.1.0] — 2026-09-13

Phase 1 (core routine): initial wrapper around the BHMIE subroutine.

### Added

- `bhmiepy.bhmie(x, m, nang)`: vectorized single-sphere Mie
  computation — the parameter loop runs in Fortran (f2py + meson
  build), NumPy broadcasting between `x` and `m`.
- `bhmiepy.compute(radius, wavelength, refractive_index)` convenience
  wrapper over physical quantities.
- `MieResult` dataclass (`qext/qsca/qback/g/s1/s2/angles` plus derived
  `qabs` and `albedo`), full input validation with readable errors
  (non-finite values, `Im(m) < 0`, series-order limit).
- Validation: cross-checks against the original single-precision F77
  BHMIE (vendored as an oracle in the extension), Rayleigh and
  geometric limits, energy conservation, and negative-input tests.
