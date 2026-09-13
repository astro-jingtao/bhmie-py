"""Single-sphere Mie scattering via the BHMIE Fortran routine.

This module exposes the BHMIE subroutine (Bohren & Huffman 1983, Appendix A,
as modified by B. T. Draine and maintained by T. P. Robitaille) through a
NumPy-friendly vectorized interface.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass

import numpy as np

from . import _bhmiepy_ext

# Mirror of nmxx in _fortran/bhmie.f90: the upstream routine allocates its
# logarithmic-derivative array with this size and responds to an exceeded
# series order with a Fortran "stop", which would kill the whole Python
# process -- so the limit is enforced here, before the call.
_NMXX = 1_000_000


@dataclass
class MieResult:
    """Result of a single-sphere Mie computation.

    Attributes
    ----------
    x, m : array or scalar
        Size parameter(s) ``x = 2*pi*a/lambda`` and relative refractive
        index(es) ``m`` as given, broadcast to a common shape.
    angles : array
        Scattering angles in radians, ``linspace(0, pi, 2*nang - 1)``.
        ``s1``/``s2`` are evaluated at these angles.
    qext, qsca, qback, g : array or scalar
        Efficiencies for extinction, scattering, backscattering, and the
        asymmetry parameter ``<cos(theta)>``.
    s1, s2 : array
        Amplitude scattering functions, shape ``(2*nang - 1,) + x.shape``.
        Following the upstream convention, ``s1 = -i*f_22`` (incident and
        scattered E perpendicular to the scattering plane) and
        ``s2 = -i*f_11`` (parallel); see Bohren & Huffman (1983) eq. 3.12.
    """

    x: np.ndarray
    m: np.ndarray
    angles: np.ndarray
    qext: np.ndarray
    qsca: np.ndarray
    qback: np.ndarray
    g: np.ndarray
    s1: np.ndarray
    s2: np.ndarray

    @property
    def qabs(self) -> np.ndarray:
        """Absorption efficiency ``qext - qsca``."""
        return self.qext - self.qsca

    @property
    def albedo(self) -> np.ndarray:
        """Single-scattering albedo ``qsca / qext``."""
        return self.qsca / self.qext


def _validate_inputs(x_arr: np.ndarray, m_arr: np.ndarray, nang: int) -> None:
    if not np.all(np.isfinite(x_arr)):
        raise ValueError("x contains non-finite values")
    if not np.all(np.isfinite(m_arr)):
        raise ValueError("m contains non-finite values")
    if np.any(x_arr <= 0.0):
        raise ValueError("x must be positive (got x <= 0)")
    if np.any(m_arr.imag < 0.0):
        raise ValueError("Im(m) must be >= 0 (absorbing medium convention)")


def _check_series_order(x_b: np.ndarray, m_b: np.ndarray) -> None:
    """Reject inputs whose Mie series order would exceed the Fortran limit.

    Mirrors the upstream computation ``nmx = nint(max(x + 4*x**(1/3) + 2,
    |m|*x)) + 15``; the routine would otherwise abort the process.
    """
    xstop = x_b + 4.0 * x_b ** (1.0 / 3.0) + 2.0
    nmx = np.round(np.maximum(xstop, np.abs(m_b) * x_b)).astype(np.int64) + 15
    if np.any(nmx > _NMXX):
        i = int(np.argmax(nmx))
        raise ValueError(
            f"size parameter too large: series order nmx = {nmx.flat[i]} "
            f"exceeds the Fortran limit nmxx = {_NMXX} "
            f"(at x = {x_b.flat[i]}, m = {m_b.flat[i]})"
        )


def bhmie(x, m, nang: int = 90) -> MieResult:
    """Mie scattering efficiencies and amplitudes for homogeneous spheres.

    Parameters
    ----------
    x : float or array_like
        Size parameter ``2*pi*a/lambda`` (> 0). Broadcast against ``m``.
    m : complex or array_like of complex
        Relative refractive index (sphere / medium). ``Im(m) >= 0`` is
        required (absorption convention).
    nang : int, optional
        Number of angles between 0 and 90 degrees; amplitudes are returned
        at ``2*nang - 1`` angles covering 0 to 180 degrees inclusive
        (default 90, i.e. one-degree sampling).

    Returns
    -------
    MieResult
        See :class:`MieResult`. Scalar ``x``/``m`` give scalar efficiencies
        and ``s1``/``s2`` of shape ``(2*nang - 1,)``; array inputs give
        results of the broadcast shape (amplitudes: ``(2*nang - 1,) +
        shape``).

    Notes
    -----
    The heavy loop runs in Fortran. Each call allocates a workspace of
    ``nmxx`` complex doubles (~16 MB) inside the routine, so very large
    parameter batches are memory-bandwidth bound but not memory limited.
    """
    nang = operator.index(nang)
    if nang < 2:
        raise ValueError(f"nang must be >= 2 (got {nang})")

    x_arr = np.asarray(x, dtype=np.float64)
    m_arr = np.asarray(m, dtype=np.complex128)
    _validate_inputs(x_arr, m_arr, nang)

    x_b, m_b = np.broadcast_arrays(x_arr, m_arr)
    _check_series_order(x_b, m_b)

    shape = x_b.shape
    n = int(x_b.size)
    flat_x = np.ascontiguousarray(x_b).reshape(-1)
    flat_m = np.ascontiguousarray(m_b).reshape(-1)

    qext, qsca, qback, g, s1, s2 = _bhmiepy_ext.bhmie_vec(flat_x, flat_m, nang)

    nang2 = 2 * nang - 1
    angles = np.linspace(0.0, math.pi, nang2)

    return MieResult(
        x=x_b,
        m=m_b,
        angles=angles,
        qext=qext.reshape(shape),
        qsca=qsca.reshape(shape),
        qback=qback.reshape(shape),
        g=g.reshape(shape),
        s1=s1.reshape((nang2,) + shape),
        s2=s2.reshape((nang2,) + shape),
    )


def compute(radius, wavelength, refractive_index, nang: int = 90) -> MieResult:
    """Mie scattering from physical quantities.

    ``radius`` and ``wavelength`` may carry any (identical) unit; only their
    ratio enters through ``x = 2*pi*radius/wavelength``.

    Parameters
    ----------
    radius, wavelength : float or array_like
        Sphere radius and wavelength (both > 0, same unit, broadcast
        against each other and against ``refractive_index``).
    refractive_index : complex or array_like of complex
        Relative refractive index with ``Im >= 0``.

    Returns
    -------
    MieResult
        Same as :func:`bhmie`.
    """
    radius = np.asarray(radius, dtype=np.float64)
    wavelength = np.asarray(wavelength, dtype=np.float64)
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("radius must be finite and positive")
    if not np.all(np.isfinite(wavelength)) or np.any(wavelength <= 0.0):
        raise ValueError("wavelength must be finite and positive")

    x = 2.0 * math.pi * radius / wavelength
    return bhmie(x, refractive_index, nang=nang)
