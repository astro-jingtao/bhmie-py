"""Tests for the Phase 2 dust pipeline (distributions, materials, ported
log-log helpers, pipeline plumbing, and golden-data validation)."""

from pathlib import Path

import numpy as np
import pytest

import bhmiepy
from bhmiepy import dust
from bhmiepy._loglog import integral_loglog, integral_loglog_subset, interp1d_loglog, logspace

EXAMPLES = Path(__file__).resolve().parents[1] / "bhmie" / "examples"


class TestLoglogPort:
    """The ported fortranlib algorithms (exactness properties)."""

    def test_power_law_integral_exact(self):
        # the log-log interpolant through a power law IS the power law,
        # so the integral is analytic
        x = np.logspace(0, 2, 7)
        for p in (-3.5, -1.0, 0.5, 2.0):
            y = x**p
            got = integral_loglog(x, y)
            if abs(p + 1.0) < 1e-12:
                want = np.log(100.0)
            else:
                want = (100.0 ** (p + 1.0) - 1.0) / (p + 1.0)
            np.testing.assert_allclose(got, want, rtol=1e-13)

    def test_subset_matches_analytic_and_total(self):
        x = np.logspace(0, 2, 11)
        y = x**-2.0
        np.testing.assert_allclose(
            integral_loglog_subset(x, y, 2.0, 50.0), 0.5 - 50.0**-1.0, rtol=1e-13
        )
        np.testing.assert_allclose(
            integral_loglog_subset(x, y, x[0], x[-1]), integral_loglog(x, y), rtol=1e-13
        )
        # inverted limits yield the negated integral (upstream does not
        # guard against x1 > x2; the pipeline never passes inverted limits)
        np.testing.assert_allclose(
            integral_loglog_subset(x, y, 5.0, 2.0), -(0.5 - 0.2), rtol=1e-13
        )

    def test_interp_exact_for_power_law_and_zero_segments(self):
        x = np.array([1.0, 10.0, 100.0])
        y = x**2
        np.testing.assert_allclose(interp1d_loglog(x, y, 30.0), 900.0, rtol=1e-14)
        np.testing.assert_allclose(interp1d_loglog(x, y, 10.0), 100.0, rtol=0)
        y0 = np.array([1.0, 0.0, 4.0])
        assert interp1d_loglog(x, y0, 50.0) == 0.0

    def test_interp_out_of_bounds_raises(self):
        x = np.array([1.0, 10.0])
        with pytest.raises(ValueError, match="out of bounds"):
            interp1d_loglog(x, x, 11.0)

    def test_logspace_matches_numpy(self):
        np.testing.assert_allclose(
            logspace(0.01, 1000.0, 250), np.logspace(np.log10(0.01), np.log10(1000.0), 250)
        )


class TestDistributions:
    def test_power_law_weights(self):
        d = dust._make_distribution(dust.PowerLaw(0.005, 1.0, -3.5))
        assert d.weight_number(0.005, 1.0) == pytest.approx(1.0)
        assert d.weight_number(0.0, 2.0) == pytest.approx(1.0)  # full coverage
        assert d.weight_number(2.0, 3.0) == 0.0
        np.testing.assert_allclose(
            d.weight_number(0.01, 0.1),
            (0.1**-2.5 - 0.01**-2.5) / (1.0**-2.5 - 0.005**-2.5),
            rtol=1e-14,
        )
        # partial overlap is clamped to [amin, amax]
        np.testing.assert_allclose(
            d.weight_number(0.5, 2.0), d.weight_number(0.5, 1.0), rtol=1e-14
        )

    def test_power_law_average_volume_against_numeric(self):
        d = dust._make_distribution(dust.PowerLaw(0.005, 1.0, -3.5))
        a = np.logspace(np.log10(0.005), 0.0, 2_000_001)
        n = a**-3.5
        want = (
            np.trapezoid(4.0 / 3.0 * 3.1415926 * a**3 * n, a) / np.trapezoid(n, a)
        )
        np.testing.assert_allclose(d.average_volume(), want, rtol=1e-4)

    def test_peaked_grid_matches_upstream_construction(self):
        d = dust._make_distribution(dust.PeakedPowerLaw(0.005, 0.15, -3.5))
        i = np.arange(1, 1001)
        np.testing.assert_allclose(d.a, 10.0 ** (-5.0 + i / 1000.0 * 10.0), rtol=0)
        assert d.amax == pytest.approx(1e5)
        assert 0.0 < d.weight_number(0.005, 1e5) <= 1.0 + 1e-12

    def test_table_distribution_from_file(self, tmp_path):
        f = tmp_path / "dist.tbl"
        a = np.logspace(-2, 0, 50)
        f.write_text("\n".join(f"{ai} {ni}" for ai, ni in zip(a, a**-3.0)))
        d = dust._make_distribution(dust.TableDistribution.from_file(f))
        assert d.weight_number(a[0], a[-1]) == pytest.approx(1.0, abs=1e-12)
        assert d.average_volume() > 0.0


class TestMaterial:
    def test_interpolate_loglog(self):
        m = dust.Material(
            wavelengths=np.array([0.5, 1.0, 2.0, 4.0]),
            refractive_indices=(np.array([0.5, 1.0, 2.0, 4.0]) ** 1.5)
            + 0.0j,
        )
        mi = m.interpolate(np.array([1.0, 2.8284271247461903]))
        # log-log interpolation is exact for power laws
        np.testing.assert_allclose(
            mi.refractive_indices.real, [1.0, 2.8284271247461903**1.5], rtol=1e-12
        )
        np.testing.assert_allclose(mi.refractive_indices.imag, 0.0)

    def test_interpolate_range_error(self):
        m = dust.Material(
            wavelengths=np.array([1.0, 2.0]), refractive_indices=np.array([1.5, 1.6])
        )
        with pytest.raises(ValueError, match="only covers"):
            m.interpolate(np.array([0.5, 2.5]))

    def test_from_file_roundtrip(self, tmp_path):
        f = tmp_path / "ri.dat"
        f.write_text("1.0 1.5 0.1\n2.0 1.6 0.2\n")
        m = dust.Material.from_file(f)
        np.testing.assert_allclose(m.wavelengths, [1.0, 2.0])
        np.testing.assert_allclose(m.refractive_indices, [1.5 + 0.1j, 1.6 + 0.2j])


class TestDustPipelinePlumbing:
    """Small end-to-end run against a naive per-point reimplementation."""

    WAV = np.array([0.5, 1.0, 2.0, 10.0])
    NANG, NSMALL, NA = 5, 2, 10
    AMIN, AMAX = 0.01, 0.3

    def _material(self):
        return dust.Material(
            wavelengths=np.array([0.3, 1.0, 5.0, 50.0]),
            refractive_indices=np.array(
                [1.7 + 0.5j, 1.6 + 0.1j, 1.5 + 0.01j, 1.4 + 0.001j]
            ),
        )

    def test_matches_naive_accumulation(self):
        comp = dust.Component(self._material(), dust.PowerLaw(0.01, 0.3, -3.5), 1.0, 3.3)
        res = dust.compute_dust_properties(
            [comp], self.WAV, self.AMIN, self.AMAX, self.NA,
            self.NANG, self.NSMALL, gas_to_dust=100.0,
        )

        dist = dust._make_distribution(comp.distribution)
        mat = comp.material.interpolate(self.WAV)
        nwav = len(self.WAV)
        nang = self.NANG + self.NSMALL
        cext = np.zeros(nwav)
        csca = np.zeros(nwav)
        gsca = np.zeros(nwav)
        logamin = np.log10(self.AMIN)
        logastep = (np.log10(self.AMAX) - logamin) / self.NA
        for ia in range(1, self.NA + 1):
            a1 = 10.0 ** (logamin + logastep * (ia - 1))
            a = 10.0 ** (logamin + logastep * (ia - 0.5))
            a2 = 10.0 ** (logamin + logastep * ia)
            w = dist.weight_number(a1, a2)
            if w <= 0.0:
                continue
            cs = np.pi * a * a * 1.0e-8
            for iw in range(nwav):
                r = bhmiepy.bhmie(
                    2.0 * np.pi * a / self.WAV[iw],
                    mat.refractive_indices[iw],
                    nang=nang,
                )
                cext[iw] += r.qext * cs * w
                csca[iw] += r.qsca * cs * w
                gsca[iw] += r.g * r.qsca * cs * w

        np.testing.assert_allclose(res.cext, cext, rtol=1e-12)
        np.testing.assert_allclose(res.csca, csca, rtol=1e-12)
        np.testing.assert_allclose(res.g, gsca / csca, rtol=1e-10)

        nang2 = 2 * nang - 1
        assert res.s11.shape == (nwav, nang2)
        assert res.angles.shape == (nang2,)
        assert res.angles[0] == 0.0
        np.testing.assert_allclose(res.angles[-1], np.pi)
        # mirrored grid: theta_j = pi - theta_{N-1-j}
        np.testing.assert_allclose(res.angles, np.pi - res.angles[::-1], atol=1e-15)
        assert res.albedo.shape == (nwav,)
        assert res.polarization_90.shape == (nwav,)

    def test_vec_ang_uniform_grid_matches_vec(self):
        x = np.array([0.5, 2.0, 7.0])
        m = np.array([1.5 + 0.01j, 1.6 + 0.1j, 1.3 + 0.0j])
        uni = np.linspace(0.0, np.pi / 2, 9)
        uni[0] = 0.0
        uni[-1] = np.pi / 2
        out1 = bhmiepy._bhmiepy_ext.bhmie_vec(x, m, 9)
        out2 = bhmiepy._bhmiepy_ext.bhmie_vec_ang(x, m, angles=uni)
        for a, b in zip(out1, out2):
            np.testing.assert_allclose(a, b, rtol=1e-13)


class TestDustValidation:
    def test_bad_n_angles(self):
        comp = dust.Component(
            dust.Material(np.array([1.0, 2.0]), np.array([1.5, 1.6])),
            dust.PowerLaw(0.1, 0.3, -3.5), 1.0, 3.0,
        )
        with pytest.raises(ValueError, match="n_angles"):
            dust.compute_dust_properties([comp], [1.0], 0.1, 0.3, 5, 1, 0)

    def test_bad_size_range(self):
        comp = dust.Component(
            dust.Material(np.array([1.0, 2.0]), np.array([1.5, 1.6])),
            dust.PowerLaw(0.1, 0.3, -3.5), 1.0, 3.0,
        )
        with pytest.raises(ValueError, match="amin < amax"):
            dust.compute_dust_properties([comp], [1.0], 0.3, 0.1, 5, 5)

    def test_negative_gas_to_dust(self):
        comp = dust.Component(
            dust.Material(np.array([1.0, 2.0]), np.array([1.5, 1.6])),
            dust.PowerLaw(0.1, 0.3, -3.5), 1.0, 3.0,
        )
        with pytest.raises(ValueError, match="gas_to_dust"):
            dust.compute_dust_properties([comp], [1.0], 0.1, 0.3, 5, 5, 0, -1.0)

    def test_unknown_distribution_type(self):
        with pytest.raises(TypeError, match="unknown size distribution"):
            dust._make_distribution("not-a-distribution")


class TestParameterFile:
    def test_parse_mrn77(self):
        inp = dust.read_parameter_file(EXAMPLES / "mrn77" / "mrn77.in")
        assert inp.prefix == "mrn77"
        assert inp.output_format == 1
        assert (inp.amin, inp.amax, inp.na) == (0.005, 2.0, 1000)
        assert (inp.n_angles, inp.n_small_angles) == (181, 10)
        assert inp.gas_to_dust == pytest.approx(141.84)
        assert inp.wavelengths.size == 250
        np.testing.assert_allclose(inp.wavelengths[0], 0.01, rtol=1e-12)
        np.testing.assert_allclose(inp.wavelengths[-1], 1000.0, rtol=1e-12)
        assert len(inp.components) == 3
        assert all(isinstance(c.distribution, dust.PowerLaw) for c in inp.components)
        assert [c.distribution.amax for c in inp.components] == [1.0, 0.25, 0.25]
        assert inp.components[0].material.wavelengths.size > 0

    def test_parse_kmh94_full_table(self):
        inp = dust.read_parameter_file(EXAMPLES / "kmh94_full" / "kmh94_3.1_full.in")
        assert inp.output_format == 2
        assert all(isinstance(c.distribution, dust.TableDistribution) for c in inp.components)
        assert inp.components[0].distribution.a.size > 0


# ---------------------------------------------------------------------------
# Golden-data validation: full runs against the upstream reference outputs.
# These reproduce the upstream CLI computation (~700k Mie calls each, i.e.
# minutes) -- deselected by default, run explicitly with `pytest -m slow`.


def _compare(name, got, want, rtol=1e-3):
    got = np.asarray(got, dtype=float)
    want = np.asarray(want, dtype=float)
    np.testing.assert_allclose(got, want, rtol=rtol, atol=1e-4 * np.max(np.abs(want)),
                               err_msg=name)


@pytest.mark.slow
class TestGoldenUpstreamOutputs:
    def test_mrn77_summary(self):
        inp = dust.read_parameter_file(EXAMPLES / "mrn77" / "mrn77.in")
        res = dust.compute_dust_properties(
            inp.components, inp.wavelengths, inp.amin, inp.amax, inp.na,
            inp.n_angles, inp.n_small_angles, inp.gas_to_dust,
        )
        ref = np.loadtxt(EXAMPLES / "mrn77" / "mrn77.summary")
        _compare("wav", res.wavelengths, ref[:, 0])
        _compare("cext", res.cext, ref[:, 1])
        _compare("csca", res.csca, ref[:, 2])
        _compare("kappa", res.kappa_ext, ref[:, 3])
        _compare("g", res.g, ref[:, 4])
        _compare("pol90", res.polarization_90, ref[:, 5])

    def test_kmh94_full_files(self):
        inp = dust.read_parameter_file(EXAMPLES / "kmh94_full" / "kmh94_3.1_full.in")
        res = dust.compute_dust_properties(
            inp.components, inp.wavelengths, inp.amin, inp.amax, inp.na,
            inp.n_angles, inp.n_small_angles, inp.gas_to_dust,
        )
        d = str(EXAMPLES / "kmh94_full" / "kmh94_3.1_full" / "kmh94_3.1_full")
        _compare("wav", res.wavelengths, np.loadtxt(d + ".wav"))
        _compare("mu", res.mu, np.loadtxt(d + ".mu"))
        _compare("alb", res.albedo, np.loadtxt(d + ".alb"))
        _compare("chi", res.kappa_ext, np.loadtxt(d + ".chi"))
        _compare("g", res.g, np.loadtxt(d + ".g"))
        for f, attr in (("f11", "s11"), ("f12", "s12"), ("f33", "s33"), ("f34", "s34")):
            _compare(f, getattr(res, attr), np.loadtxt(d + f".{f}"))
