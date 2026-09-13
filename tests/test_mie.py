"""Tests for the bhmiepy public API (Phase 1: core BHMIE routine)."""

import numpy as np
import pytest

from bhmiepy import MieResult, _bhmiepy_ext, bhmie, compute


class TestAgainstF77Oracle:
    """Cross-check the double-precision port against the original F77 routine.

    The F77 interface is single precision (inputs truncated on entry,
    amplitudes rounded on exit), so agreement is limited to single-
    precision level; the rtol values below are individually relaxed from
    numpy's assert_allclose default accordingly.
    """

    NANG = 45

    CASES = [
        (0.1, 1.5 + 0.0j),
        (1.0, 1.5 + 0.01j),
        (5.0, 1.33 + 0.0j),
        (10.0, 1.7 + 0.5j),
        (50.0, 1.5 + 0.001j),
        (100.0, 2.0 + 1.0j),
    ]

    @staticmethod
    def _oracle(x, m, nang):
        return _bhmiepy_ext.bhmie_f77_ref(x, m, nang)

    def test_efficiencies_fixed_cases(self):
        for x, m in self.CASES:
            res = bhmie(x, m, nang=self.NANG)
            qext, qsca, qback, g, _, _ = self._oracle(x, m, self.NANG)
            np.testing.assert_allclose(res.qext, qext, rtol=1e-5, atol=1e-9,
                                       err_msg=f"qext at x={x}, m={m}")
            np.testing.assert_allclose(res.qsca, qsca, rtol=1e-5, atol=1e-9,
                                       err_msg=f"qsca at x={x}, m={m}")
            np.testing.assert_allclose(res.qback, qback, rtol=1e-5, atol=1e-9,
                                       err_msg=f"qback at x={x}, m={m}")
            np.testing.assert_allclose(res.g, g, rtol=1e-5, atol=1e-9,
                                       err_msg=f"g at x={x}, m={m}")

    def test_amplitudes_fixed_cases(self):
        for x, m in self.CASES:
            res = bhmie(x, m, nang=self.NANG)
            _, _, _, _, s1o, s2o = self._oracle(x, m, self.NANG)
            for got, want, name in [
                (res.s1, s1o, "s1"),
                (res.s2, s2o, "s2"),
            ]:
                scale = np.max(np.abs(want))
                assert np.max(np.abs(got - want)) <= 1e-4 * scale, (
                    f"{name} at x={x}, m={m}: max diff "
                    f"{np.max(np.abs(got - want)):.3e} vs scale {scale:.3e}"
                )

    def test_random_sample_batched(self):
        # One batched call against per-point oracle calls: also exercises
        # the Fortran-side loop over parameter arrays.
        rng = np.random.default_rng(42)
        x = np.exp(rng.uniform(np.log(0.01), np.log(200.0), size=12))
        m = rng.uniform(1.1, 2.5, size=12) + 1j * rng.uniform(0.0, 1.0, size=12)

        res = bhmie(x, m, nang=self.NANG)
        assert res.qext.shape == (12,)
        for i in range(12):
            qext, qsca, qback, g, s1o, _ = self._oracle(x[i], m[i], self.NANG)
            np.testing.assert_allclose(res.qext[i], qext, rtol=1e-5, atol=1e-9)
            np.testing.assert_allclose(res.qsca[i], qsca, rtol=1e-5, atol=1e-9)
            np.testing.assert_allclose(res.qback[i], qback, rtol=1e-5, atol=1e-9)
            np.testing.assert_allclose(res.g[i], g, rtol=1e-5, atol=1e-9)
            scale = np.max(np.abs(s1o))
            assert np.max(np.abs(res.s1[:, i] - s1o)) <= 1e-4 * scale


class TestPhysicalLimits:
    def test_nonabsorbing_qsca_equals_qext(self):
        # For real m, qsca and qext agree analytically. Numerically qext
        # comes from the real part of the forward amplitude, which for
        # small x involves heavy cancellation (~8 digits lost at x = 0.1),
        # so rtol is individually relaxed to 1e-6.
        for x in [0.1, 1.0, 10.0, 50.0]:
            res = bhmie(x, 1.5 + 0.0j, nang=10)
            np.testing.assert_allclose(res.qsca, res.qext, rtol=1e-6,
                                       err_msg=f"x={x}")

    def test_rayleigh_small_x(self):
        x, m = 1.0e-3, 1.5 + 0.0j
        lorentz = (m**2 - 1.0) / (m**2 + 2.0)
        expected = (8.0 / 3.0) * x**4 * abs(lorentz) ** 2
        res = bhmie(x, m, nang=10)
        np.testing.assert_allclose(res.qsca, expected, rtol=1e-3)

    def test_rayleigh_g_tends_to_zero(self):
        res = bhmie(1.0e-3, 1.5 + 0.01j, nang=10)
        assert abs(res.g) < 1e-5

    def test_geometric_limit_large_x(self):
        res = bhmie(300.0, 1.5 + 0.01j, nang=10)
        assert 1.9 < float(res.qext) < 2.4

    def test_qabs_nonnegative_for_absorbing(self):
        for x, m in TestAgainstF77Oracle.CASES:
            if m.imag > 0:
                res = bhmie(x, m, nang=10)
                assert float(res.qabs) >= -1e-10 * float(res.qext), (
                    f"qabs < 0 at x={x}, m={m}"
                )


class TestValidation:
    def test_nang_below_two(self):
        with pytest.raises(ValueError, match="nang must be >= 2"):
            bhmie(1.0, 1.5j, nang=1)

    def test_nang_not_integer(self):
        with pytest.raises(TypeError):
            bhmie(1.0, 1.5j, nang=2.5)

    @pytest.mark.parametrize("x", [0.0, -1.0, np.nan, np.inf])
    def test_bad_x(self, x):
        with pytest.raises(ValueError):
            bhmie(x, 1.5 + 0.0j, nang=10)

    def test_negative_imaginary_m(self):
        with pytest.raises(ValueError, match="Im\\(m\\) must be >= 0"):
            bhmie(1.0, 1.5 - 0.1j, nang=10)

    def test_nonfinite_m(self):
        with pytest.raises(ValueError):
            bhmie(1.0, complex(np.nan, 0.5), nang=10)

    def test_x_beyond_series_limit(self):
        # nmx = nint(x + 4x**(1/3) + 2) + 15 must stay below nmxx = 1e6
        with pytest.raises(ValueError, match="series order"):
            bhmie(1.2e6, 1.5 + 0.0j, nang=10)

    def test_compute_bad_radius(self):
        with pytest.raises(ValueError, match="radius"):
            compute(-0.1, 0.25, 1.5 + 0.0j)

    def test_compute_bad_wavelength(self):
        with pytest.raises(ValueError, match="wavelength"):
            compute(0.1, 0.0, 1.5 + 0.0j)


class TestShapesAndBroadcasting:
    def test_scalar_inputs_give_scalar_result(self):
        res = bhmie(2.5, 1.6 + 0.01j, nang=90)
        assert isinstance(res, MieResult)
        assert res.qext.shape == ()
        assert res.s1.shape == (179,)
        assert res.s2.shape == (179,)

    def test_vector_x_scalar_m(self):
        res = bhmie(np.logspace(-1, 2, 5), 1.6 + 0.01j, nang=5)
        assert res.qext.shape == (5,)
        assert res.s1.shape == (9, 5)

    def test_broadcasting(self):
        x = np.logspace(-1, 1, 3)[:, None]
        m = (1.3 + 1j * np.linspace(0.0, 0.5, 4))[None, :]
        res = bhmie(x, m, nang=5)
        assert res.qext.shape == (3, 4)
        assert res.s1.shape == (9, 3, 4)

    def test_angles_grid(self):
        nang = 90
        res = bhmie(1.0, 1.5 + 0.0j, nang=nang)
        assert res.angles.shape == (2 * nang - 1,)
        assert res.angles[0] == 0.0
        np.testing.assert_allclose(res.angles[-1], np.pi)
        # uniform grid 0..pi with 2*nang-1 points: step pi / (2*(nang-1))
        np.testing.assert_allclose(np.diff(res.angles), np.pi / (2 * (nang - 1)))

    def test_derived_properties_shape(self):
        res = bhmie(np.array([1.0, 10.0]), 1.5 + 0.01j, nang=5)
        assert res.qabs.shape == (2,)
        assert res.albedo.shape == (2,)
        np.testing.assert_allclose(res.qabs + res.qsca, res.qext)


class TestCompute:
    def test_compute_matches_bhmie(self):
        radius, wavelength, m = 0.1, 0.25, 1.6 + 0.01j
        got = compute(radius, wavelength, m, nang=30)
        want = bhmie(2.0 * np.pi * radius / wavelength, m, nang=30)
        for field in ["qext", "qsca", "qback", "g"]:
            np.testing.assert_allclose(getattr(got, field), getattr(want, field))

    def test_compute_vectorized(self):
        got = compute(np.array([0.05, 0.1, 0.2]), 0.25, 1.6 + 0.01j, nang=5)
        assert got.qext.shape == (3,)
