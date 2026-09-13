@echo off
rem ===========================================================================
rem One-shot development build for a f2py + meson-python Fortran extension
rem on Windows (conda-forge flang + VS BuildTools).
rem
rem Adapt these three paths, then run from the repo root:
rem   1. vcvarsall.bat location (VS 2022 BuildTools or Community)
rem   2. conda root (...\Scripts\activate.bat)
rem   3. conda env name
rem
rem Why the environment tweaks below are needed (flang 21-23 / LLVM,
rem conda-forge; see the skill's SKILL.md for the full table):
rem   AR/RANLIB  - meson looks for "ar"/"gar" as the Fortran static linker,
rem                the conda env only ships llvm-ar/llvm-ranlib.
rem   LIB        - flang embeds a reference to its runtime library
rem                (flang_rt.runtime.*.lib) in the objects it produces, but
rem                link.exe does not know where that library lives.
rem   FFLAGS     - the conda flang-rt package ships no DLL for the dynamic
rem                runtime, so Fortran code must be compiled against the
rem                static runtime (-static-libflangrt) or the built .pyd
rem                would fail to load.
rem
rem Do NOT inline this as a single `cmd /c 'a && set LIB=%LIB%;x && b'`
rem line: %LIB% then expands at parse time, BEFORE vcvarsall has run, and
rem silently destroys the MSVC library search path.
rem ===========================================================================
setlocal

call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
if errorlevel 1 exit /b 1

call C:\Users\<user>\anaconda3\Scripts\activate.bat <env>
if errorlevel 1 exit /b 1

set AR=llvm-ar
set RANLIB=llvm-ranlib

for /f "usebackq delims=" %%i in (`flang-new --print-resource-dir`) do set "FLANG_RESOURCE_DIR=%%i"
set "LIB=%LIB%;%FLANG_RESOURCE_DIR%\lib\x86_64-pc-windows-msvc"

set "FFLAGS=-static-libflangrt"

pip install -e . --no-build-isolation %*
set "RC=%ERRORLEVEL%"
endlocal & exit /b %RC%
