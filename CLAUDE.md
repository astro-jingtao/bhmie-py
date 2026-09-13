# bhmiepy — project instructions

## What this is

bhmiepy is a Python wrapper (package name `bhmiepy`) around the BHMIE
Mie-scattering Fortran code from https://github.com/hyperion-rt/bhmie.
Development is phased: Phase 1 wraps the core single-sphere BHMIE routine;
Phase 2 wraps the upstream dust-population layer. The authoritative roadmap
is [docs/DEV_PLAN.md](docs/DEV_PLAN.md) — keep it updated as phases complete.

## Authoritative docs — read before relying on any note

- `.my_note/` contains the project owner's personal notes. In particular,
  `.my_note/overall_guide.md` was used as a reference ONLY during the first
  two steps of the project: (1) repository initialization and (2) writing
  docs/DEV_PLAN.md. Do NOT consult it for any later work — it is not
  maintained and may silently go stale.
- Authoritative documents: docs/DEV_PLAN.md (roadmap), README.md (overview),
  this file (working rules).

## Code provenance and licensing rules (important)

The repository deliberately mixes code under two licenses:

- `src/bhmiepy/_fortran/bhmie.f90`, `src/bhmiepy/_fortran/types.f90` —
  copies of files from upstream hyperion-rt/bhmie, **BSD-2-Clause**,
  © 2012-13 Thomas P. Robitaille. When editing these files:
  - never delete or alter the license header block at the top of the file;
  - update the "Modifications for bhmiepy (relative to the upstream file)"
    line in the header whenever a change diverges from upstream.
- `src/bhmiepy/_fortran/bhmie_f77.f` — verbatim copy of the original F77
  BHMIE (plus a header and one SAVE statement, both documented in the
  header). Same BSD-2-Clause terms; keep it unmodified.
- `src/bhmiepy/_loglog.py` — Python translation of algorithms from
  astrofrog/fortranlib, **BSD-2-Clause**, © 2009-13 Thomas P. Robitaille
  (see LICENSES/BSD-2-Clause-fortranlib.txt). It defines the upstream
  numerical conventions — do not "simplify" its formulas away from the
  fortranlib originals (golden tests depend on them).
- `src/bhmiepy/_fortran/bhmiepy_ext.f90` and everything under `src/bhmiepy/`
  except the two files above — new code written for bhmiepy, **MIT**.
- `bhmie/` — pristine snapshot of upstream hyperion-rt/bhmie at commit
  `37c072909a67d1027120762680a1c0b350875398` (its `.git` was removed on
  purpose when it was folded into this repository). It is a **read-only**
  reference and test-fixture source (`bhmie/examples/` holds refractive-index
  tables and reference outputs). Never build from it and never edit it; the
  live copies are in `src/bhmiepy/_fortran/`.
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
- Editable install: run `scripts\dev-install.cmd` from the repo root. It
  chains vcvarsall + conda activate and sets the extra environment pieces
  flang needs (`AR`/`RANLIB=llvm-*`, flang runtime dir on `LIB`,
  `FFLAGS=-static-libflangrt`); the script header documents each one. Do not
  inline this chain as a single `cmd /c '...'` line — `%LIB%` then expands
  before vcvarsall has run and breaks the MSVC library search.
- Run tests: `conda run -n bhmiepy python -m pytest`
- Rebuild after touching only Python files is unnecessary (editable); after
  touching Fortran or meson.build, rerun `scripts\dev-install.cmd`. If meson
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
  validation comes from `bhmie/examples/` and from compiling
  `bhmie/original/bhmie.f` as an independent oracle.
- Prefer failing loudly over silent fallbacks. Do not hide build/runtime
  errors.
- Temporary and throwaway artifacts stay in the system temp directory or
  project build dirs, never at drive roots or the home directory; clean them
  up when done.
- Do not stage or commit unless the user asks. When asked, include
  `.claude/knowledge/` changes from the knowledge-curator flow in the same
  round of commits.
