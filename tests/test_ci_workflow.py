"""Guards the CI workflow file against package-spec composition bugs.

The first real CI run (2026-09-15) failed on every Windows/numpy-2 job
during environment creation with::

    critical libmamba Error parsing version ">=2". Version contains
    invalid characters in >=2.

The workflow rendered ``numpy=${{ ... && '=1.26.*' || '>=2,<3' }}`` into
``numpy=>=2,<3``: the literal ``numpy=`` prefix was composed with a branch
that did not carry the package name. The numpy-1.26 branch happened to
produce the *valid* spec ``numpy==1.26.*``, which is why those jobs got
past environment creation and hid the bug. These tests keep every numpy
spec branch a complete, self-contained micromamba matchspec.
"""

import glob
import importlib.util
import platform
import re
from pathlib import Path

import pytest

CI_FILE = Path(__file__).parents[1] / ".github" / "workflows" / "ci.yml"
PE_IMPORTS = Path(__file__).parents[1] / ".github" / "scripts" / "pe_imports.py"


class TestCreateArgsSpecs:
    """The numpy ternary in ci.yml must render to valid micromamba specs."""

    def _numpy_branches(self):
        text = CI_FILE.read_text(encoding="utf-8")
        m = re.search(
            r"matrix\.numpy == '1\.26' && '([^']+)' \|\| '([^']+)'", text
        )
        assert m, (
            "numpy version ternary not found in .github/workflows/ci.yml "
            "(workflow restructured? update this test)"
        )
        return m.group(1), m.group(2)

    def test_branches_are_complete_specs(self):
        # Each branch must carry the package name itself, so that no
        # outer 'numpy=' prefix can compose with it into a doubled
        # operator like 'numpy=>=2,<3'.
        for spec in self._numpy_branches():
            assert re.fullmatch(
                r"numpy(==|=|>=|<=|<|>|!=)[0-9][0-9.*,\s<>!=]*", spec
            ), f"not a complete micromamba spec: {spec!r}"

    def test_no_doubled_operators(self):
        # The original bug: 'numpy=' + '>=2,<3' -> 'numpy=>=2,<3'.
        for spec in self._numpy_branches():
            assert "=>" not in spec and "=<" not in spec, (
                f"doubled operator in spec: {spec!r}"
            )


class TestWindowsGfortranCorner:
    """The experimental Windows gfortran job must stay a SEPARATE job.

    A matrix ``include:`` entry whose keys all match an existing matrix
    combination MERGES into that job instead of adding a new one. The
    single-value ``fortran: ['flang']`` dimension is what makes the
    gfortran include entry create a new job; deleting the dimension
    would silently turn the windows/py3.12/numpy-2 flang job INTO the
    gfortran corner, losing flang coverage for that cell.
    """

    def test_fortran_dimension_present(self):
        text = CI_FILE.read_text(encoding="utf-8")
        assert re.search(r"^\s*fortran:\s*\[\s*'flang'\s*\]", text, re.M), (
            "matrix.fortran dimension missing: the gfortran include "
            "would merge into an existing job instead of adding one"
        )

    def test_gfortran_corner_is_required(self):
        text = CI_FILE.read_text(encoding="utf-8")
        assert re.search(r"fortran:\s*'gfortran'", text), (
            "gfortran include entry not found in ci.yml"
        )
        assert "continue-on-error" not in text, (
            "the gfortran corner is regular blocking coverage now; "
            "re-adding continue-on-error must be a conscious decision, "
            "not a leftover"
        )


class TestPeImports:
    """The PE import-table dumper must stay working.

    It is the Windows mingw-gfortran corner's only DLL diagnostic (Git
    Bash ships neither ldd nor objdump). Validated here against the
    built extension itself: the fast suite already requires the editable
    build to import, so the .pyd is present whenever these tests run.
    """

    def _load(self):
        spec = importlib.util.spec_from_file_location("pe_imports", PE_IMPORTS)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    @pytest.mark.skipif(
        platform.system() != "Windows", reason="PE parser targets Windows .pyd"
    )
    def test_parses_built_extension(self):
        pyds = glob.glob(
            str(Path(__file__).parents[1] / "build" / "*" / "_bhmiepy_ext*.pyd")
        )
        assert pyds, "built _bhmiepy_ext .pyd not found under build/"
        imports = self._load().pe_imports(pyds[0])
        assert imports, "no DLL imports parsed from the built extension"
        assert any(dll.upper().startswith("PYTHON3") for dll in imports), imports
        assert any(dll.upper() == "KERNEL32.DLL" for dll in imports), imports
