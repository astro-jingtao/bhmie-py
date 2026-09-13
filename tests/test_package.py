"""Smoke tests for the bhmiepy package skeleton."""

import bhmiepy
from bhmiepy import _bhmiepy_ext


class TestPackageBasics:
    def test_version_is_set(self):
        assert isinstance(bhmiepy.__version__, str)
        assert bhmiepy.__version__ != ""

    def test_backend_reports_double_precision(self):
        prec_digits, range_digits = _bhmiepy_ext.backend_info()
        assert prec_digits >= 15
        assert range_digits >= 307
