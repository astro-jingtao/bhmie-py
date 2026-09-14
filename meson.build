project('bhmiepy', 'c', 'fortran',
  version : '0.2.0',
  meson_version : '>=1.1.0',
  default_options : ['buildtype=release', 'c_std=c99'])

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

py.install_sources(
  'src/bhmiepy/__init__.py',
  'src/bhmiepy/_loglog.py',
  'src/bhmiepy/_upstream_cli.py',
  'src/bhmiepy/mie.py',
  'src/bhmiepy/dust.py',
  subdir : 'bhmiepy',
)

# Fortran sources compiled and linked into the extension. types.f90 and
# bhmie.f90 are derived from upstream hyperion-rt/bhmie (BSD-2-Clause);
# bhmiepy_ext.f90 is new code written for bhmiepy (MIT). bhmiepy_ext.f90 is
# both parsed by f2py (via custom_target above) and compiled directly -- the
# f2py step only generates the C API layer, the subroutine implementations
# come from compiling the raw source.

fortran_sources = files(
  'src/bhmiepy/_fortran/bhmiepy_ext.f90',
  'src/bhmiepy/_fortran/types.f90',
  'src/bhmiepy/_fortran/bhmie.f90',
  'src/bhmiepy/_fortran/ref/bhmie_upstream.f90',
  'src/bhmiepy/_fortran/ref/bhmie_f77.f',
)

# For bare subroutines f2py generates the F77-style wrapper
# (-f2pywrappers.f, possibly empty); the -f2pywrappers2.f90 variant only
# appears for Fortran module subroutines. Verified with:
#   python -m numpy.f2py bhmiepy_ext.f90 -m _bhmiepy_ext
# The custom_target must be declared in the root meson.build: f2py writes
# its output files to the working directory, which for meson is the build
# root, and meson looks for declared outputs relative to the declaring
# meson.build's build subdir.

f2py_ext_sources = custom_target('bhmiepy-ext-f2py',
  input : ['src/bhmiepy/_fortran/bhmiepy_ext.f90'],
  output : ['_bhmiepy_extmodule.c', '_bhmiepy_ext-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_bhmiepy_ext'],
)

py.extension_module('_bhmiepy_ext',
  f2py_ext_sources,
  fortran_sources,
  incdir_f2py / 'fortranobject.c',
  include_directories : inc_np,
  dependencies : py_dep,
  install : true,
  subdir : 'bhmiepy',
)

# Optional pristine upstream CLI ("bhmie_ref"), compiled verbatim from the
# read-only upstream/ git submodule (hyperion-rt/bhmie) for A/B
# benchmarking (benchmarks/bench_dust_upstream.py) and the slow
# upstream-equivalence test. Source order mirrors upstream Makefile.in.
# The compiler is whatever meson picks for this platform (flang in the
# dev env, gfortran elsewhere) -- no per-platform commands here.
# Auto-skips when the submodule (or its nested fortranlib) is not
# initialized, so plain user builds from a submodule-less checkout never
# require it.

fs = import('fs')
upstream_ref = get_option('upstream_ref')
upstream_ref_ready = fs.is_file('upstream/src/fortranlib/src/lib_array.f90')
if upstream_ref.enabled() and not upstream_ref_ready
  error(
    'upstream_ref=enabled but the bhmie/ submodule is not initialized.\n' +
    'Run: git submodule update --init --recursive',
  )
endif
if upstream_ref.disabled() or not upstream_ref_ready
  summary({'bhmie_ref (upstream CLI)' : 'not built'}, section : 'upstream reference')
else
  upstream_ref_link_args = []
  if host_machine.system() == 'windows'
    # upstream's compute_dust_properties keeps megabytes of automatic
    # arrays (s11..s34, sized by wavelengths x angles) on the stack; the
    # Windows default 1 MB stack overflows (0xC00000FD) where Linux's
    # 8 MB does not. Raise the stack reserve instead of touching the
    # pristine sources.
    upstream_ref_link_args = ['-Wl,/STACK:67108864']
  endif
  executable(
    'bhmie_ref',
    'upstream/src/types.f90',
    'upstream/src/fortranlib/src/lib_array.f90',
    'upstream/src/material.f90',
    'upstream/src/distributions.f90',
    'upstream/src/bhmie.f90',
    'upstream/src/bhmie_caller.f90',
    'upstream/src/main.f90',
    # flang 23 aborts on upstream material.f90/distributions.f90 with
    # "runtime derived type info descriptor was not generated" (derived
    # types with allocatable components as intent(out) dummies). The flag
    # the diagnostic itself names skips those descriptors; they only feed
    # runtime diagnostics, and the A/B outputs (compared against the
    # recorded golden files) verify the binary stays correct. No source
    # changes to the pristine upstream code.
    fortran_args : ['-mllvm', '--ignore-missing-type-desc'],
    link_args : upstream_ref_link_args,
    install : false,
  )
  summary({'bhmie_ref (upstream CLI)' : 'built (from bhmie/ submodule)'}, section : 'upstream reference')
endif
