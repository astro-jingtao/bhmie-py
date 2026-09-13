"""Dust-population scattering properties (Phase 2).

Python port of the pipeline in hyperion-rt/bhmie (bhmie_caller.f90,
distributions.f90, material.f90): size-distribution averaging of Mie
results over one or more grain materials, producing cross-sections,
extinction opacity, albedo, asymmetry parameter, and scattering-matrix
elements S11/S12/S33/S34.

Numerical conventions follow upstream exactly (log-log interpolation and
integration from fortranlib, ported in ``_loglog.py``; units: micron for
sizes/wavelengths, g/cm^3 for grain density, cm^2 for cross-sections,
cm^2/g for opacity), so results match the upstream CLI tool; the golden
tests validate that against the reference outputs shipped in ``bhmie/examples``.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import _bhmiepy_ext
from ._loglog import integral_loglog, integral_loglog_subset, interp1d_loglog, logspace

# upstream pi used in average_volume (distributions.f90: 3.1415926) --
# kept truncated for bit-level fidelity of the port
_PI_UPSTREAM = 3.1415926


# ---------------------------------------------------------------------------
# Materials


@dataclass
class Material:
    """Tabulated refractive index m(lambda): wavelengths in micron (ascending
    or descending, as in the upstream ri-data files), complex indices."""

    wavelengths: np.ndarray
    refractive_indices: np.ndarray

    @classmethod
    def from_file(cls, path) -> "Material":
        path = Path(path)
        data = np.loadtxt(path)
        if data.ndim != 2 or data.shape[1] < 3:
            raise ValueError(
                f"refractive index file {path} must have three columns "
                "(wavelength, Re(m), Im(m))"
            )
        return cls(
            wavelengths=np.ascontiguousarray(data[:, 0], dtype=np.float64),
            refractive_indices=np.ascontiguousarray(
                data[:, 1] + 1j * data[:, 2], dtype=np.complex128
            ),
        )

    def interpolate(self, wav: np.ndarray) -> "Material":
        """New Material interpolated onto ``wav`` (log-log in Re and Im
        separately), mirroring upstream interpolate_material."""
        wav = np.asarray(wav, dtype=np.float64)
        wmin, wmax = self.wavelengths.min(), self.wavelengths.max()
        if wav.max() > wmax or wav.min() < wmin:
            raise ValueError(
                f"refractive index table only covers "
                f"[{wmin:.4f}, {wmax:.4f}] microns; requested "
                f"[{wav.min():.4f}, {wav.max():.4f}]"
            )
        ref_real = interp1d_loglog(self.wavelengths, self.refractive_indices.real, wav)
        ref_imag = interp1d_loglog(self.wavelengths, self.refractive_indices.imag, wav)
        return Material(
            wavelengths=np.array(wav, dtype=np.float64, copy=True),
            refractive_indices=ref_real + 1j * ref_imag,
        )


# ---------------------------------------------------------------------------
# Size distributions


@dataclass
class PowerLaw:
    """n(a) proportional to a**apower on [amin, amax] (upstream 'power')."""

    amin: float
    amax: float
    apower: float


@dataclass
class PeakedPowerLaw:
    """n(a) proportional to a**apower * exp(-a/aturn) for a >= amin
    (upstream 'ped'), tabulated internally exactly like upstream: 1000
    points at a_i = 10**(-5 + i/1000 * 10), i = 1..1000."""

    amin: float
    aturn: float
    apower: float


@dataclass
class TableDistribution:
    """Number distribution given as a two-column table (a in micron,
    relative n(a)), normalized internally like upstream."""

    a: np.ndarray
    n: np.ndarray

    @classmethod
    def from_file(cls, path) -> "TableDistribution":
        path = Path(path)
        data = np.loadtxt(path)
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(
                f"size distribution file {path} must have two columns "
                "(size in microns, relative number)"
            )
        return cls(a=data[:, 0], n=data[:, 1])


class _NumericalDistribution:
    """Internal: normalized tabulated distribution (upstream type 2)."""

    def __init__(self, a: np.ndarray, n_raw: np.ndarray, amin: float, amax: float):
        self.a = np.asarray(a, dtype=np.float64)
        n_raw = np.asarray(n_raw, dtype=np.float64)
        self.m_raw = n_raw * self.a**3
        self.n = n_raw / integral_loglog(self.a, n_raw)
        self.m = self.m_raw / integral_loglog(self.a, self.m_raw)
        self.amin = amin
        self.amax = amax

    def weight_number(self, a1: float, a2: float) -> float:
        if a2 < self.amin or a1 > self.amax:
            return 0.0
        if a1 < self.amin and a2 > self.amax:
            return 1.0
        return integral_loglog_subset(self.a, self.n, max(a1, self.amin), min(a2, self.amax))

    def average_volume(self) -> float:
        return integral_loglog(
            self.a, self.n * 4.0 / 3.0 * _PI_UPSTREAM * self.a**3
        ) / integral_loglog(self.a, self.n)


def _make_distribution(spec):
    """Normalize a public distribution spec into an object with
    weight_number / average_volume."""
    if isinstance(spec, PowerLaw):
        return _PowerLawImpl(spec.amin, spec.amax, spec.apower)
    if isinstance(spec, PeakedPowerLaw):
        i = np.arange(1, 1001)
        a = 10.0 ** (-5.0 + i / 1000.0 * 10.0)
        n = a ** spec.apower * np.exp(-a / spec.aturn)
        return _NumericalDistribution(a, n, spec.amin, a[-1])
    if isinstance(spec, TableDistribution):
        a = np.asarray(spec.a, dtype=np.float64)
        n = np.asarray(spec.n, dtype=np.float64)
        if np.any(np.diff(a) <= 0) and np.any(np.diff(a) >= 0):
            raise ValueError("table distribution must be sorted (ascending or descending)")
        if np.any(a <= 0) or np.any(n < 0):
            raise ValueError("table distribution requires a > 0 and n >= 0")
        return _NumericalDistribution(a, n, a[0], a[-1])
    raise TypeError(f"unknown size distribution spec: {type(spec).__name__}")


class _PowerLawImpl:
    def __init__(self, amin, amax, apower):
        self.amin = amin
        self.amax = amax
        self.apower = apower
        p1 = apower + 1.0
        self._norm = amax**p1 - amin**p1
        if self._norm == 0.0:
            raise ValueError("degenerate power-law normalization (apower = -1)")

    def weight_number(self, a1: float, a2: float) -> float:
        if a2 < self.amin or a1 > self.amax:
            return 0.0
        if a1 < self.amin and a2 > self.amax:
            return 1.0
        p1 = self.apower + 1.0
        return (min(a2, self.amax) ** p1 - max(a1, self.amin) ** p1) / self._norm

    def average_volume(self) -> float:
        p1 = self.apower + 1.0
        p4 = self.apower + 4.0
        return (
            4.0
            / 3.0
            * _PI_UPSTREAM
            * (self.amax**p4 - self.amin**p4)
            / self._norm
            * p1
            / p4
        )


# ---------------------------------------------------------------------------
# Dust model and computation


@dataclass
class Component:
    """One chemical component: refractive index table, size distribution,
    mass abundance (relative weights are renormalized internally), and
    grain material density in g/cm^3."""

    material: Material
    distribution: PowerLaw | PeakedPowerLaw | TableDistribution
    abundance_mass: float
    density: float


@dataclass
class DustResult:
    """Size-distribution-averaged dust properties per wavelength.

    cext/csca/cback: cross-sections in cm^2 (per normalized particle mix);
    kappa_ext: extinction opacity in cm^2/g including gas via gas_to_dust;
    g: <cos(theta)>; s11/s12/s33/s34: scattering matrix elements with shape
    (nwav, 2*(n_angles+n_small_angles)-1) on the mirrored 0..pi grid."""

    wavelengths: np.ndarray
    angles: np.ndarray
    cext: np.ndarray
    csca: np.ndarray
    cback: np.ndarray
    kappa_ext: np.ndarray
    g: np.ndarray
    s11: np.ndarray
    s12: np.ndarray
    s33: np.ndarray
    s34: np.ndarray

    @property
    def albedo(self) -> np.ndarray:
        return self.csca / self.cext

    @property
    def mu(self) -> np.ndarray:
        """cos(angles), matching the upstream '.mu' output file."""
        return np.cos(self.angles)

    @property
    def polarization_90(self) -> np.ndarray:
        """-S12/S11 at 90 degrees (upstream summary's last column)."""
        i90 = (len(self.angles) - 1) // 2  # index of pi/2 in the mirrored grid
        return -self.s12[:, i90] / self.s11[:, i90]


def compute_dust_properties(
    components,
    wavelengths,
    amin: float,
    amax: float,
    na: int,
    n_angles: int,
    n_small_angles: int = 0,
    gas_to_dust: float = 100.0,
) -> DustResult:
    """Compute size-distribution-averaged dust properties (port of upstream
    compute_dust_properties).

    Parameters
    ----------
    components : sequence of Component
    wavelengths : array_like
        Wavelengths in micron (used for every component).
    amin, amax, na : float, float, int
        Overall size range in micron and number of log-spaced size bins.
    n_angles : int
        Number of scattering angles between 0 and 90 degrees (>= 2).
    n_small_angles : int, optional
        Extra fine angles near 0 degrees (also mirrored near 180).
    gas_to_dust : float, optional
        Gas-to-dust mass ratio for the opacity (0 for dust-only).
    """
    n_angles = operator.index(n_angles)
    n_small_angles = operator.index(n_small_angles)
    na = operator.index(na)
    if n_angles < 2:
        raise ValueError("n_angles must be >= 2")
    if n_small_angles < 0:
        raise ValueError("n_small_angles must be >= 0")
    if na < 1:
        raise ValueError("na must be >= 1")
    if not 0.0 < amin < amax:
        raise ValueError("size range must satisfy 0 < amin < amax")
    if gas_to_dust < 0.0:
        raise ValueError("gas_to_dust must be >= 0")
    if len(components) == 0:
        raise ValueError("at least one component is required")

    wavelengths = np.asarray(wavelengths, dtype=np.float64)
    if np.any(wavelengths <= 0.0) or not np.all(np.isfinite(wavelengths)):
        raise ValueError("wavelengths must be positive and finite")
    nwav = wavelengths.size

    abundance_mass = np.array([c.abundance_mass for c in components], dtype=np.float64)
    if np.any(abundance_mass <= 0.0):
        raise ValueError("abundance_mass must be > 0")
    densities = np.array([c.density for c in components], dtype=np.float64)
    if np.any(densities <= 0.0):
        raise ValueError("density must be > 0")
    # upstream main.f90 renormalizes the mass abundances to sum to one
    abundance_mass = abundance_mass / abundance_mass.sum()

    nang = n_angles + n_small_angles
    nang2 = 2 * nang - 1

    # --- materials interpolated to the requested wavelength grid
    materials = [c.material.interpolate(wavelengths) for c in components]

    # --- distributions normalized like upstream
    dists = [_make_distribution(c.distribution) for c in components]

    # --- number abundances (upstream bhmie_caller.f90)
    average_particle_mass = np.array(
        [d.average_volume() / 1.0e12 * dens for d, dens in zip(dists, densities)]
    )
    if np.any(average_particle_mass <= 0.0):
        raise ValueError("average particle mass must be > 0")
    abundance_number = abundance_mass / average_particle_mass
    abundance_number = abundance_number / abundance_number.sum()

    # --- angle grid (exact upstream construction; the Fortran routine
    #     requires angles[0] == 0 and angles[-1] == pi/2 exactly)
    angles = np.empty(nang)
    angles[0] = 0.0
    k = np.arange(1, n_angles)
    angles[n_small_angles + 1 : nang] = k / (n_angles - 1) * (np.pi / 2.0)
    for ia in range(1, n_small_angles + 1):
        cur = n_small_angles + 1 - ia
        angles[cur] = angles[cur + 1] / np.sqrt(10.0)
    angles[-1] = np.pi / 2.0

    # --- log-spaced size bins (edges a1/a2, centers a)
    logamin = np.log10(amin)
    logastep = (np.log10(amax) - logamin) / na
    idx = np.arange(1, na + 1, dtype=np.float64)
    a_lo = 10.0 ** (logamin + logastep * (idx - 1.0))
    a_mid = 10.0 ** (logamin + logastep * (idx - 0.5))
    a_hi = 10.0 ** (logamin + logastep * idx)

    cext = np.zeros(nwav)
    csca = np.zeros(nwav)
    cback = np.zeros(nwav)
    gsca = np.zeros(nwav)
    s11 = np.zeros((nwav, nang2))
    s12 = np.zeros((nwav, nang2))
    s33 = np.zeros((nwav, nang2))
    s34 = np.zeros((nwav, nang2))

    for ic, (dist, mat) in enumerate(zip(dists, materials)):
        for ia in range(na):
            weight = dist.weight_number(a_lo[ia], a_hi[ia]) * abundance_number[ic]
            if weight <= 0.0:
                continue
            a = a_mid[ia]
            # microns -> cm^2: pi*a^2 * 1e-8
            cross_section = np.pi * a * a * 1.0e-8
            x = 2.0 * np.pi * a / wavelengths
            qext, qsca, qback, g, s1, s2 = _bhmiepy_ext.bhmie_vec_ang(
                x, mat.refractive_indices, angles=angles
            )
            cext += qext * cross_section * weight
            csca += qsca * cross_section * weight
            cback += qback * cross_section * weight
            gsca += g * qsca * cross_section * weight
            a1sq = np.abs(s1) ** 2
            a2sq = np.abs(s2) ** 2
            s11 += weight * (0.5 * (a1sq + a2sq)).T
            s12 += weight * (0.5 * (-a1sq + a2sq)).T
            s2cs1 = s2 * np.conj(s1)
            s33 += weight * s2cs1.real.T
            s34 += weight * s2cs1.imag.T

    kappa_ext = cext * np.sum(abundance_mass / average_particle_mass)
    kappa_ext = kappa_ext / (1.0 + gas_to_dust)

    with np.errstate(invalid="ignore", divide="ignore"):
        gsca = gsca / csca

    angles_full = np.concatenate([angles, np.pi - angles[-2::-1]])

    return DustResult(
        wavelengths=wavelengths,
        angles=angles_full,
        cext=cext,
        csca=csca,
        cback=cback,
        kappa_ext=kappa_ext,
        g=gsca,
        s11=s11,
        s12=s12,
        s33=s33,
        s34=s34,
    )


# ---------------------------------------------------------------------------
# Upstream parameter-file reader (compatibility + golden tests)


@dataclass
class DustInputFile:
    """Parsed upstream bhmie parameter file (.in)."""

    prefix: str
    output_format: int
    amin: float
    amax: float
    na: int
    n_angles: int
    n_small_angles: int
    gas_to_dust: float
    wavelengths: np.ndarray
    components: list = field(default_factory=list)


def _first_token(line: str) -> str:
    tok = line.split()[0]
    return tok.strip("'\"")


def _floats(line: str, n: int):
    vals = [float(t) for t in line.split()[:n]]
    if len(vals) != n:
        raise ValueError(f"expected {n} numbers on line: {line!r}")
    return vals


def read_parameter_file(path) -> DustInputFile:
    """Read an upstream bhmie parameter file (see bhmie/README.md for the
    format). Relative refractive-index and table paths are resolved against
    the parameter file's directory, as the upstream CLI expects to be run
    from there."""
    path = Path(path)
    lines = path.read_text().splitlines()
    base = path.parent

    def line(i):
        return lines[i]

    parsed = DustInputFile(
        prefix=_first_token(line(0)),
        output_format=int(_first_token(line(1))),
        amin=_floats(line(2), 1)[0],
        amax=_floats(line(3), 1)[0],
        na=int(_first_token(line(4))),
        n_angles=int(_first_token(line(5))),
        n_small_angles=int(_first_token(line(6))),
        gas_to_dust=_floats(line(8), 1)[0],
        wavelengths=logspace(*_floats(line(9), 3)[0:2], int(line(9).split()[2])),
        components=[],
    )
    n_components = int(_first_token(line(7)))

    i = 10
    for _ in range(n_components):
        if not lines[i].strip().startswith("---"):
            raise ValueError(f"expected component separator at line {i + 1}")
        i += 1
        abundance = _floats(line(i), 1)[0]
        i += 1
        density = _floats(line(i), 1)[0]
        i += 1
        material = Material.from_file(base / _first_token(line(i)))
        i += 1
        dist_type = _first_token(line(i)).lower()
        i += 1
        if dist_type == "power":
            a0, a1_, a2_ = _floats(line(i), 3)
            dist = PowerLaw(a0, a1_, a2_)
        elif dist_type == "ped":
            a0, aturn, apower = _floats(line(i), 3)
            dist = PeakedPowerLaw(a0, aturn, apower)
        elif dist_type == "table":
            dist = TableDistribution.from_file(base / _first_token(line(i)))
        else:
            raise ValueError(f"unknown distribution type {dist_type!r} at line {i + 1}")
        i += 1
        parsed.components.append(
            Component(
                material=material,
                distribution=dist,
                abundance_mass=abundance,
                density=density,
            )
        )
    return parsed
