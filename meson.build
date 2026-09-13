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
  'src/bhmiepy/_fortran/bhmie_upstream.f90',
  'src/bhmiepy/_fortran/bhmie_f77.f',
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
