---
name: f2py-meson-flang-windows
description: Build Python Fortran extension modules on Windows with f2py, meson-python and conda-forge flang (LLVM Flang). Use when creating or repairing a pyproject/meson setup that wraps Fortran for Python on Windows, when a conda env needs a Fortran compiler for f2py, when a build fails with flang/MSVC errors such as "LNK1104 flang_rt.runtime", "LNK1104 libcmt.lib", "Unknown linker(s): ar/gar", "LNK2001 unresolved external symbol <sub>_", "ascii codec can't decode byte" from crackfortran, or "runtime derived type info descriptor was not generated" aborting flang, or when calling or running legacy Fortran crashes with no traceback (exit code -1073741571 / 0xC00000FD stack overflow).
---

# f2py + meson + flang on Windows

## Role

You are a build engineer for Python packages that compile Fortran on
Windows via f2py + meson-python, using conda-forge flang (LLVM Flang) and
the MSVC toolchain. Your job is to:

1. Provision a working build environment (conda packages + VS BuildTools).
2. Write a correct build script and root `meson.build`.
3. Diagnose the characteristic link/parse failures and apply the known fix.

Working progression:

```text
env setup → build script → meson.build wiring → build + test → troubleshoot
```

---

# Hard Rules

## R1 — Build-script environment setup goes in a `.bat` file, never a one-line `cmd /c '...'` chain

A single `cmd /c 'call vcvarsall.bat x64 && set LIB=%LIB%;<dir> && pip
install -e .'` expands `%LIB%` when the line is parsed — before vcvarsall
has run — silently replacing it with an empty prefix, which later surfaces
as `LNK1104: libcmt.lib`. A `.bat` file expands per-line at execution time.
(If a one-liner is unavoidable, use delayed expansion `!VAR!` with
`cmd /v:on`.)

## R2 — Every Fortran file containing subroutine implementations must be listed as a `py.extension_module` source

Files passed only to the f2py `custom_target` are parsed for signatures but
never compiled — their symbols will be missing at link time
(`LNK2001: unresolved external symbol <name>_`), even though the symbols
flang emits look correct.

## R3 — The f2py `custom_target` is declared in the ROOT `meson.build`

f2py writes its generated files to the working directory (the build root),
while meson resolves declared outputs relative to the declaring
meson.build's build subdir. Declaring it in a nested meson.build fails with
`FileNotFoundError: <path>\<m>-f2pywrappers.f`.

## R4 — Fortran sources fed to f2py must be pure ASCII

f2py's crackfortran reads sources as ASCII; a UTF-8 em-dash or `©` in a
comment aborts with `'ascii' codec can't decode byte 0xe2`. Use `--` and
`(c)`.

## R5 — Keep the f2py-parsed file "F77-clean"; reach Fortran-90 modules through an external bridge

The file passed to the f2py `custom_target` must contain only bare
subroutines declared with plain `double precision` / `double complex`
types — no Fortran modules, no `use`, no interface blocks, no named kind
parameters (f2py's parser mishandles all of these for externally linked
code, and an unresolvable kind silently corrupts signatures). Reach
Fortran-90 module code by appending a small bare bridge subroutine after
`end module` in the module file; its external symbol is what the
f2py-parsed code calls. This also avoids the silent collision between a
module procedure name and a same-named F77 external — both would bind to
the `name_` symbol, calling the wrong routine.

## R6 — Call f2py entry points with keyword arguments after the leading arrays

f2py REORDERS the generated signature: an argument that owns another
argument's dimension (e.g. `nang` dimensioning `angles(nang)`) becomes
optional and is moved AFTER the array it dimensions — source
`(x, m, nang, angles)` generates `mod.f(x, m, angles, [n, nang])`.
Positional calls silently misbind (a scalar lands on the array slot and
fails with a confusing `shape(...) == n` error). Pass everything after
the leading arrays by keyword, and when in doubt print `ext.__doc__` to
see the generated signature.

## R7 — Guard EVERY entry into process-fatal Fortran; port predicates with their NaN semantics

Legacy Fortran routines report invalid input with `stop`/`error stop`,
which kills the whole Python process with no traceback. A guard on one
Python caller does not protect other callers of the same extension —
validate at every exported entry, or give the Fortran wrapper an `ierr`
out-argument that callers must check. When porting the routine's guard
formulas to Python: compare in float64 BEFORE any int cast (an int64
cast wraps at ~9.2e18 and silently bypasses range checks), and preserve
comparison polarity — upstream `if (w > 0)` skips NaN where a natural
rewrite `if w <= 0: pass-through` executes it; port as `if not (w > 0)`.

---

# Environment Setup (once per machine / env)

```text
conda create -n <env> -c conda-forge python=3.12 numpy meson-python ninja pytest pip flang
conda install -n <env> -c conda-forge flang-rt_win-64
```

Key facts (LLVM Flang 21–23, conda-forge win-64):

- `fortran_metapackage` does **not** exist on win-64; the compiler package
  is `flang`.
- The flang runtime is a **separate package** `flang-rt_win-64`. Without it,
  linking fails with `LNK1104: flang_rt.runtime.static.lib`.
- The env ships **no DLL** for the dynamic runtime, so Fortran must be
  compiled against the static runtime: `FFLAGS=-static-libflangrt`.
- C compiler: VS BuildTools `cl` via `vcvarsall.bat x64` (meson needs a C
  compiler for f2py's generated glue; flang uses the MSVC linker).
- flang on the windows-msvc target already emits gfortran-style `name_`
  symbols, matching f2py's generated C calls — `-funderscoring` is NOT
  needed (verify with `llvm-nm` on a compiled object before suspecting
  symbol mangling).

# Build Script

Copy `templates/dev-install.cmd`, adapt three paths (vcvarsall location,
conda root, env name). It sets everything the toolchain cannot discover:

| Variable | Value | Why |
| --- | --- | --- |
| `AR` / `RANLIB` | `llvm-ar` / `llvm-ranlib` | meson probes for `ar`/`gar`; the env only ships `llvm-*` |
| `LIB` | append `<resource-dir>\lib\x86_64-pc-windows-msvc` | link.exe must find `flang_rt.runtime.*.lib`; resource dir = `flang-new --print-resource-dir` |
| `FFLAGS` | `-static-libflangrt` | no dynamic-runtime DLL exists in the env |

# meson.build Wiring

Copy `templates/meson.build.txt` as the project root `meson.build`. Its
load-bearing details:

- `custom_target` runs `python -m numpy.f2py <wrapper.f90> -m <ext>` and
  declares the **exact** generated filenames. For bare (non-module)
  subroutines f2py generates `<ext>module.c` + `<ext>-f2pywrappers.f` (the
  wrapper file is written even when empty). `-f2pywrappers2.f90` appears
  only for Fortran-module subroutines. When unsure, probe:
  `python -m numpy.f2py <file>.f90 -m <ext>` in a scratch dir and list what
  appeared.
- `py.extension_module` sources = generated files + the wrapper `.f90`
  itself (see R2) + other implementation Fortran + `fortranobject.c` from
  `numpy.f2py.get_include()`, with numpy + f2py include dirs.
- Pure-Python package files are installed explicitly via
  `py.install_sources(..., subdir: '<pkg>')` — meson installs only what is
  declared.

# Troubleshooting

| Symptom | Cause | Fix |
| --- | --- | --- |
| `PackagesNotFoundError: fortran_metapackage` | wrong package name on win-64 | install `flang` |
| `LNK1104: flang_rt.runtime.static.lib` | runtime package missing | `conda install flang-rt_win-64` |
| `LNK1104: flang_rt.runtime.dynamic.lib` (final link) | resource dir not on `LIB` | append it (build script) |
| `.pyd` builds but fails to load | linked against dynamic runtime, no DLL ships | `FFLAGS=-static-libflangrt` |
| `Unknown linker(s): [['ar'], ['gar']]` | meson static-linker probe | `AR=llvm-ar RANLIB=llvm-ranlib` |
| `LNK1104: libcmt.lib` (sanity check) | `LIB` clobbered by parse-time `%LIB%` expansion | R1: use a `.bat` |
| `FileNotFoundError: <...>-f2pywrappers.f` | custom_target in nested meson.build | R3: move to root |
| `'ascii' codec can't decode byte ...` in crackfortran | non-ASCII chars in Fortran | R4: pure ASCII |
| `LNK2001: unresolved external symbol <sub>_` | implementation file only given to f2py, not compiled | R2: add to extension sources |
| Python process dies, no traceback, exit `-1073741571` (`0xC00000FD`) | legacy F77 local array larger than the 1 MB stack — flang stacks large locals, gfortran moves them to static storage | add `SAVE` to the vendored routine (restores F77 static semantics), record it in the file header |
| plain Fortran exe crashes `0xC00000FD` (works on Linux) | megabyte-scale **automatic arrays** (locals sized by dummy args, e.g. `s11(nwav, nang)`) go on the stack; Windows default stack reserve is 1 MB vs Linux 8 MB — distinct from the F77-local/`SAVE` row above | raise the reserve instead of editing pristine sources: `link_args : ['-Wl,/STACK:67108864']` guarded by `host_machine.system() == 'windows'` |
| flang aborts: `runtime derived type info descriptor was not generated and skipExternalRttiDefinition and ignoreMissingTypeDescriptors options are not set` (`LLVM ERROR: aborting`) | flang 23 codegen gap on derived types with allocatable components as `intent(out)` dummies | per-target `fortran_args : ['-mllvm', '--ignore-missing-type-desc']` — kebab spelling; the camelCase name the diagnostic prints is rejected by the option parser. Descriptors only feed runtime diagnostics; validate the binary's outputs afterwards |
| `import` triggers a ninja rebuild that fails (`WinError 2`, `LNK1104`, or `Regenerating build files` exits 1) | meson-python editable installs rebuild-on-import; the rebuild subprocess runs plain ninja/meson without vcvars/`LIB`/`FFLAGS` — fires after editing Fortran, after any meson.build change, and after `meson configure -D...` | rerun the full dev-install script after touching Fortran or meson.build, before invoking Python; run Python through an activated env |
| need a built non-installed target (auxiliary exe/lib) from an editable install | meson-python editable artifacts live in `build/<python-tag>/` (e.g. `build/cp312/`), not `build/` directly | locate by scanning `build/` plus one level of subdirs (suffix-agnostic: no suffix on POSIX, `.exe` on Windows) |
| `Successfully installed` but `import` fails on missing symbol | stale `build/` dir from earlier attempt | delete `build/`, rebuild |

---

# Out of Scope

- Linux/macOS builds (gfortran there; none of these workarounds apply).
- setuptools / numpy.distutils legacy builds (deprecated path).
- Producing distributable wheels / cibuildwheel / CI matrices — the recipe
  here is for local dev installs (`pip install -e . --no-build-isolation`).
- Intel oneAPI (ifx) or MSYS2 gfortran toolchains — fallbacks if flang
  fails, not covered here.

---

# Supporting Files

- `templates/dev-install.cmd` — the one-shot build script; copy, adapt the
  three paths, run from the repo root. The header comment documents each
  environment tweak.
- `templates/meson.build.txt` — root meson.build pattern for an f2py
  extension inside a `src/`-layout package; copy and rename the extension
  and sources.
