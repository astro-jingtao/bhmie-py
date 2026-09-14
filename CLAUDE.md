# bhmiepy — project instructions

## What this is

bhmiepy is a Python wrapper (package name `bhmiepy`) around the BHMIE
Mie-scattering Fortran code from https://github.com/hyperion-rt/bhmie.
Development is phased: Phase 1 wraps the core single-sphere BHMIE routine;
Phase 2 wraps the upstream dust-population layer. The authoritative roadmap
is docs/DEV_PLAN.md — **local-only and untracked** (see the
authoritative-docs note below); keep it updated as phases complete.

## Authoritative docs — read before relying on any note

- `.my_note/` contains the project owner's personal notes. In particular,
  `.my_note/overall_guide.md` was used as a reference ONLY during the first
  two steps of the project: (1) repository initialization and (2) writing
  docs/DEV_PLAN.md. Do NOT consult it for any later work — it is not
  maintained and may silently go stale.
- Authoritative documents: docs/DEV_PLAN.md (roadmap — **local-only and
  untracked**, gitignored like `scripts/` because it is a personal working
  document), README.md (overview), this file (working rules).

## Code provenance and licensing rules (important)

The repository deliberately mixes code under two licenses:

- `src/bhmiepy/_fortran/bhmie.f90`, `src/bhmiepy/_fortran/types.f90` —
  copies of files from upstream hyperion-rt/bhmie, **BSD-2-Clause**,
  © 2012-13 Thomas P. Robitaille. When editing these files:
  - never delete or alter the license header block at the top of the file;
  - update the "Modifications for bhmiepy (relative to the upstream file)"
    line in the header whenever a change diverges from upstream.
- `src/bhmiepy/_fortran/ref/bhmie_f77.f` — verbatim copy of the original F77
  BHMIE (plus a header and one SAVE statement, both documented in the
  header). Same BSD-2-Clause terms; keep it unmodified. Lives in
  `_fortran/ref/` with the other comparison-only sources.
- `src/bhmiepy/_fortran/ref/bhmie_upstream.f90` — pristine copy of upstream
  `src/bhmie.f90` (module renamed `bhmie_routine_upstream` so it can
  coexist with the package copy in one build; a bridge subroutine was
  appended). BSD-2-Clause, same rules: keep the body verbatim, record any
  deviation in its header. It exists for benchmarking/auditing the
  optimized copy (`benchmarks/bench_upstream.py`, slow test
  `TestUpstreamEquivalence`).
- `src/bhmiepy/_loglog.py` — Python translation of algorithms from
  astrofrog/fortranlib, **BSD-2-Clause**, © 2009-13 Thomas P. Robitaille
  (see LICENSES/BSD-2-Clause-fortranlib.txt). It defines the upstream
  numerical conventions — do not "simplify" its formulas away from the
  fortranlib originals (golden tests depend on them).
- `src/bhmiepy/_fortran/bhmiepy_ext.f90` and everything under `src/bhmiepy/`
  except the upstream-/fortranlib-derived files listed above — new code
  written for bhmiepy, **MIT**.
- `upstream/` — git submodule (name `bhmie`, path `upstream`) of
  hyperion-rt/bhmie pinned at commit
  `37c072909a67d1027120762680a1c0b350875398`, with its own `fortranlib`
  submodule (`git submodule update --init --recursive` after clone). It is a
  **read-only** reference and test-fixture source (`upstream/examples/`
  holds refractive-index tables and reference outputs), and the source of
  the optional `bhmie_ref` executable that meson builds verbatim for A/B
  benchmarking (`-Dupstream_ref`; auto-enabled when the submodule is
  initialized, skipped otherwise so plain user builds never need it). Never
  edit it; the live copies are in `src/bhmiepy/_fortran/` (comparison-only
  copies in `src/bhmiepy/_fortran/ref/`).
- Attribution obligations live in `NOTICE`, `LICENSE` (MIT), and
  `LICENSES/BSD-2-Clause-bhmie.txt`. Keep them in sync when files move.

## Development environment (Windows)

- Dedicated conda env `bhmiepy` (conda-forge): Python 3.12, numpy,
  meson-python, ninja, pytest, flang **and `flang-rt_win-64`** (the flang
  runtime is a separate package; without it linking fails with
  `LNK1104: flang_rt.runtime.*.lib`):

  ```
  conda create -n bhmiepy -c conda-forge python=3.12 numpy meson-python ninja pytest pip flang
  conda install -n bhmiepy -c conda-forge flang-rt_win-64
  ```

- C compiler: Visual Studio 2022 BuildTools (`cl`), activated with
  `vcvarsall.bat x64`. Fortran: conda-forge flang (LLVM Flang 23).
- `scripts/` (dev-install.cmd, test.cmd) is **local-only and untracked**
  (listed in `.gitignore`): it hardcodes this machine's VS/conda paths and
  was purged from git history for that reason. Keep using it locally, but
  never commit it; README carries the general install instructions.
- Editable install: run `scripts\dev-install.cmd` from the repo root. It
  chains vcvarsall + conda activate and sets the extra environment pieces
  flang needs (`AR`/`RANLIB=llvm-*`, flang runtime dir on `LIB`,
  `FFLAGS=-static-libflangrt`); the script header documents each one. Do not
  inline this chain as a single `cmd /c '...'` line — `%LIB%` then expands
  before vcvarsall has run and breaks the MSVC library search.
- Run tests: `scripts\test.cmd` (fast suite) / `scripts\test.cmd -m slow`
  (golden-data + upstream-equivalence runs) / pass any pytest args through.
- Benchmark vs pristine upstream Fortran (both independent of the pytest
  suite): `python benchmarks\bench_upstream.py` (core routine: asserts
  identical results, reports timings) and
  `python benchmarks\bench_dust_upstream.py` (dust pipeline vs the
  `bhmie_ref` CLI built from the submodule; compares at upstream
  output-file precision, reports end-to-end timings).
- Rebuild after touching only Python files is unnecessary (editable); after
  touching Fortran or meson.build, rerun `scripts\dev-install.cmd` (a stale
  build makes the import-time rebuild fail and breaks `test.cmd`). If meson
  gets confused, delete the `build/` directory first.
- Build-system gotchas already handled (keep them handled):
  - the f2py `custom_target` must stay in the root `meson.build` (f2py writes
    outputs to the working directory);
  - every Fortran source containing subroutine implementations must be listed
    as an extension source — files only passed to f2py are parsed, not
    compiled;
  - Fortran sources must be pure ASCII (f2py's parser reads ASCII).

## Conventions

- Build system: meson + meson-python + f2py. The only file f2py parses is
  `src/bhmiepy/_fortran/bhmiepy_ext.f90` (the Python-facing wrapper); other
  Fortran files are compiled directly by meson.
- Tests: pytest, class-based style, lightweight and fast. Reference data for
  validation comes from `upstream/examples/` and from the vendored F77
  oracle `src/bhmiepy/_fortran/ref/bhmie_f77.f` (copied from
  `upstream/original/bhmie.f`).
- Prefer failing loudly over silent fallbacks. Do not hide build/runtime
  errors.
- Git-tracked content defaults to English. Where Chinese already exists,
  keep it and add an English rendering alongside — do not delete the
  Chinese.
- Temporary and throwaway artifacts stay in the system temp directory or
  project build dirs, never at drive roots or the home directory; clean them
  up when done.
- Do not stage or commit unless the user asks. When asked, include
  `.claude/knowledge/` changes from the knowledge-curator flow in the same
  round of commits.
