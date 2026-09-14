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
tests validate that against the reference outputs shipped in ``upstream/examples``.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import _bhmiepy_ext
from ._loglog import integral_loglog, integral_loglog_subset, interp1d_loglog, logspace
from .mie import _check_series_order, _validate_inputs

# distributions.f90 writes its average_volume with the literal 3.1415926;
# in Fortran that literal is single precision (~3.14159250), so the port is
# NOT bit-faithful either way -- we use the double value of the written
# digits. The ~1e-7 relative difference sits far below the golden-test
# tolerance and upstream's own output precision.
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
        # ndmin=2: a single-row file must stay a (1, 3) table, not collapse
        # to a 1-D array that the column check below would misreport
        data = np.loadtxt(path, ndmin=2)
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
        diff = np.diff(self.wavelengths)
        if not (np.all(diff >= 0.0) or np.all(diff <= 0.0)):
            # also false for NaN entries, so non-finite tables are caught
            # here rather than bisecting to arbitrary intervals downstream
            raise ValueError(
                "material wavelengths must be sorted ascending or descending"
            )
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
        # ndmin=2: same single-row concern as Material.from_file
        data = np.loadtxt(path, ndmin=2)
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(
                f"size distribution file {path} must have two columns "
                "(size in microns, relative number)"
            )
        return cls(a=data[:, 0], n=data[:, 1])


class _NumericalDistribution:
    """Internal: normalized tabulated distribution (upstream type 2).
    Expects an ascending size grid (callers flip descending tables first:
    the log-log subset integral, like fortranlib's, assumes ascending)."""

    def __init__(self, a: np.ndarray, n_raw: np.ndarray, amin: float, amax: float):
        a = np.asarray(a, dtype=np.float64)
        n_raw = np.asarray(n_raw, dtype=np.float64)
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(n_raw)):
            raise ValueError("table distribution contains non-finite values")
        if not np.all(np.diff(a) >= 0.0):
            raise ValueError("internal: size table must be ascending")
        norm = integral_loglog(a, n_raw)
        if not (norm > 0.0):  # also false for NaN
            raise ValueError(
                "table distribution has zero (or non-finite) total integral; "
                "cannot normalize"
            )
        self.a = a
        self.n = n_raw / norm
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
        if not spec.aturn > 0.0:
            raise ValueError("PeakedPowerLaw requires aturn > 0")
        i = np.arange(1, 1001)
        a = 10.0 ** (-5.0 + i / 1000.0 * 10.0)
        n = a ** spec.apower * np.exp(-a / spec.aturn)
        return _NumericalDistribution(a, n, spec.amin, a[-1])
    if isinstance(spec, TableDistribution):
        a = np.asarray(spec.a, dtype=np.float64)
        n = np.asarray(spec.n, dtype=np.float64)
        if np.any(a <= 0) or np.any(n < 0):
            raise ValueError("table distribution requires a > 0 and n >= 0")
        diff = np.diff(a)
        if np.all(diff >= 0.0):
            pass  # ascending (duplicates allowed; zero-width bins contribute 0)
        elif np.all(diff <= 0.0):
            a, n = a[::-1], n[::-1]  # flip: the subset integral assumes ascending
        else:
            raise ValueError(
                "table distribution sizes must be sorted (ascending or descending)"
            )
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

    wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))
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

    # Guard every (bin, wavelength) pair BEFORE any Fortran call: an
    # oversized size parameter would hit the Fortran "stop" and kill the
    # whole Python process (same guard mie.bhmie applies, applied to the
    # full batch here so the loop below is safe).
    for mat in materials:
        x_all = 2.0 * np.pi * a_mid[:, None] / wavelengths[None, :]
        _validate_inputs(x_all, np.broadcast_to(mat.refractive_indices, x_all.shape))
        _check_series_order(x_all, np.broadcast_to(mat.refractive_indices, x_all.shape))

    cext = np.zeros(nwav)
    csca = np.zeros(nwav)
    cback = np.zeros(nwav)
    gsca = np.zeros(nwav)
    # Fortran-ordered: bhmie_dust_accum accumulates into them in place
    # (intent(inout)); converted back to C order in the result below.
    s11 = np.zeros((nwav, nang2), order="F")
    s12 = np.zeros((nwav, nang2), order="F")
    s33 = np.zeros((nwav, nang2), order="F")
    s34 = np.zeros((nwav, nang2), order="F")

    # Main loop, one extension call per component: the core computation and
    # the weighted accumulation both happen in Fortran (bhmie_dust_accum),
    # so the huge per-point s1/s2 amplitude arrays never cross into Python
    # (materializing and reducing them in NumPy dominated the pipeline's
    # non-core runtime; the summation order now matches upstream's own
    # sequential per-bin loop).
    for ic, (dist, mat) in enumerate(zip(dists, materials)):
        weights = (
            np.array(
                [dist.weight_number(lo, hi) for lo, hi in zip(a_lo, a_hi)]
            )
            * abundance_number[ic]
        )
        # upstream semantics: 'if (weight_number > 0)' -- false for NaN, so
        # a non-finite weight skips the bin instead of poisoning the
        # accumulators (construction-time checks reject the known causes)
        keep = np.nonzero(weights > 0.0)[0]
        if keep.size == 0:
            # no contributing bin: nothing to accumulate (upstream's loop
            # body simply never runs; the zero-scattering check below is
            # the loud failure). Also avoids a zero-size f2py call.
            continue
        a = a_mid[keep]
        w = weights[keep]
        # microns -> cm^2: pi*a^2 * 1e-8
        cross = np.pi * a * a * 1.0e-8
        x = 2.0 * np.pi * a[:, None] / wavelengths[None, :]
        m = np.broadcast_to(mat.refractive_indices, (keep.size, nwav))
        _bhmiepy_ext.bhmie_dust_accum(
            x=x, m=m, weights=w, cross=cross, angles=angles,
            cext=cext, csca=csca, cback=cback, gsca=gsca,
            s11=s11, s12=s12, s33=s33, s34=s34,
        )

    kappa_ext = cext * np.sum(abundance_mass / average_particle_mass)
    kappa_ext = kappa_ext / (1.0 + gas_to_dust)

    if np.any(csca <= 0.0):
        raise ValueError(
            "total scattering cross-section is zero at "
            f"{int(np.sum(csca <= 0.0))} of {nwav} wavelengths -- the size "
            "range [amin, amax] does not overlap any component distribution"
        )
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
        s11=np.ascontiguousarray(s11),
        s12=np.ascontiguousarray(s12),
        s33=np.ascontiguousarray(s33),
        s34=np.ascontiguousarray(s34),
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
    """Read an upstream bhmie parameter file (see upstream/README.md for the
    format). Blank lines are skipped and component separators are consumed
    like the upstream CLI's list-directed reads (any single record). Relative
    refractive-index and table paths are resolved against the parameter
    file's directory, as the upstream CLI expects to be run from there."""
    path = Path(path)
    records = [ln for ln in path.read_text().splitlines() if ln.strip()]
    base = path.parent

    def record(i):
        if i >= len(records):
            raise ValueError(f"parameter file ended early (wanted record {i + 1})")
        return records[i]

    wmin, wmax, nw = _floats(record(9), 3)
    parsed = DustInputFile(
        prefix=_first_token(record(0)),
        output_format=int(_first_token(record(1))),
        amin=_floats(record(2), 1)[0],
        amax=_floats(record(3), 1)[0],
        na=int(_first_token(record(4))),
        n_angles=int(_first_token(record(5))),
        n_small_angles=int(_first_token(record(6))),
        gas_to_dust=_floats(record(8), 1)[0],
        wavelengths=logspace(wmin, wmax, int(nw)),
        components=[],
    )
    n_components = int(_first_token(record(7)))

    i = 10
    for _ in range(n_components):
        i += 1  # separator record (consumed without interpretation, upstream read(32,*))
        abundance = _floats(record(i), 1)[0]
        i += 1
        density = _floats(record(i), 1)[0]
        i += 1
        material = Material.from_file(base / _first_token(record(i)))
        i += 1
        dist_type = _first_token(record(i)).lower()
        i += 1
        if dist_type == "power":
            a0, a1_, a2_ = _floats(record(i), 3)
            dist = PowerLaw(a0, a1_, a2_)
        elif dist_type == "ped":
            a0, aturn, apower = _floats(record(i), 3)
            dist = PeakedPowerLaw(a0, aturn, apower)
        elif dist_type == "table":
            dist = TableDistribution.from_file(base / _first_token(record(i)))
        else:
            raise ValueError(f"unknown distribution type {dist_type!r} in record {i + 1}")
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
