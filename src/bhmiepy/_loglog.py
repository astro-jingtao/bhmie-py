"""Faithful NumPy ports of the fortranlib routines the upstream bhmie
pipeline depends on (astrofrog/fortranlib, lib_array module: integral_loglog,
interp1d_loglog, locate, logspace).

The upstream dust pipeline (material interpolation, size-distribution
weights, average volumes) is defined in terms of these exact algorithms;
the golden-data tests against upstream outputs only hold because this
module reproduces them operation by operation, including their special
cases (zero values, power-law slope -1, domain clamping). Do not "simplify"
the formulas away from the originals.

Upstream: https://github.com/astrofrog/fortranlib, Copyright (c) 2009-13
Thomas P. Robitaille, BSD 2-Clause License (see
LICENSES/BSD-2-Clause-fortranlib.txt in the repository). This file is a
Python translation of parts of its lib_array module and is therefore a
derivative work distributed under the same license terms.
"""

from __future__ import annotations

import numpy as np


def locate(x: np.ndarray, v: float) -> int:
    """Bisection index: x[j] <= v < x[j+1] for ascending x (mirrored for
    descending), with the upstream special cases v == x[0] -> 0 and
    v == x[-1] -> len(x) - 2. v must lie inside [min(x), max(x)]."""
    n = x.size
    if v == x[0]:
        return 0
    if v == x[-1]:
        return n - 2
    if v < min(x[0], x[-1]) or v > max(x[0], x[-1]):
        raise ValueError(f"locate: {v} outside tabulated range [{x[0]}, {x[-1]}]")
    ascnd = x[-1] >= x[0]
    jl, ju = -1, n
    while ju - jl > 1:
        jm = (ju + jl) // 2
        if (v >= x[jm]) == ascnd:
            jl = jm
        else:
            ju = jm
    return jl


def _interp1d_single_loglog(x1, y1, x2, y2, v):
    """10**(log10(y1) + frac*(log10(y2)-log10(y1))); 0 if either endpoint
    is 0 (upstream convention)."""
    if y1 == 0.0 or y2 == 0.0:
        return 0.0
    frac = (np.log10(v) - np.log10(x1)) / (np.log10(x2) - np.log10(x1))
    return 10.0 ** (np.log10(y1) + frac * (np.log10(y2) - np.log10(y1)))


def interp1d_loglog(x: np.ndarray, y: np.ndarray, v) -> np.ndarray:
    """Log-log interpolation of y(x) at v (scalar or array).

    Raises ValueError outside the tabulated range -- the caller is expected
    to range-check with domain knowledge first, mirroring upstream
    interpolate_material."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    vs = np.atleast_1d(np.asarray(v, dtype=np.float64))
    lo, hi = min(x[0], x[-1]), max(x[0], x[-1])
    if np.any(vs < lo) or np.any(vs > hi):
        raise ValueError(
            f"interpolation out of bounds: values in [{vs.min()}, {vs.max()}] "
            f"outside tabulated [{lo}, {hi}]"
        )
    if x.size == 1:
        # single-row table: the range check above pins v to the only point
        return np.full_like(vs, y[0]).reshape(np.shape(v)) if np.ndim(v) else y[0]
    out = np.empty_like(vs)
    for k, vk in enumerate(vs):
        i = locate(x, vk)
        out[k] = _interp1d_single_loglog(x[i], y[i], x[i + 1], y[i + 1], vk)
    return out.reshape(np.shape(v)) if np.ndim(v) else out[0]


def _trapezium_loglog(x1, y1, x2, y2):
    """Vectorized port of trapezium_loglog_dp: exact integral of the
    log-log interpolant (power law through the two points) over [x1, x2]."""
    x1 = np.asarray(x1, dtype=np.float64)
    y1 = np.asarray(y1, dtype=np.float64)
    x2 = np.asarray(x2, dtype=np.float64)
    y2 = np.asarray(y2, dtype=np.float64)
    out = np.zeros(x1.shape, dtype=np.float64)
    nz = (x1 != x2) & (y1 != 0.0) & (y2 != 0.0)
    if np.any(nz):
        b = np.log10(y1[nz] / y2[nz]) / np.log10(x1[nz] / x2[nz])
        near = np.abs(b + 1.0) < 1e-10
        x1n, y1n, x2n = x1[nz], y1[nz], x2[nz]
        res = np.empty(b.shape, dtype=np.float64)
        res[near] = x1n[near] * y1n[near] * np.log(x2n[near] / x1n[near])
        rn = ~near
        res[rn] = (
            y1n[rn]
            * (x2n[rn] * (x2n[rn] / x1n[rn]) ** b[rn] - x1n[rn])
            / (b[rn] + 1.0)
        )
        out[nz] = res
    return out


def integral_loglog(x: np.ndarray, y: np.ndarray) -> float:
    """Total integral of y(x) using the log-log interpolant on every
    interval (port of integral_loglog_dp)."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    return float(np.sum(_trapezium_loglog(x[:-1], y[:-1], x[1:], y[1:])))


def integral_loglog_subset(x: np.ndarray, y: np.ndarray, a1: float, a2: float) -> float:
    """Integral of y(x) over [a1, a2], with both limits clamped to the
    tabulated domain and log-log interpolation at the limits
    (port of integral_loglog_dp(x, y, x1, x2))."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = x.size
    if a1 > x[n - 1] or a2 < x[0]:
        return 0.0

    if a1 > x[0]:
        i1 = locate(x, a1)
        f1 = _interp1d_single_loglog(x[i1], y[i1], x[i1 + 1], y[i1 + 1], a1)
        xx1 = a1
    else:
        # clamp so no spurious area is added outside the domain
        i1 = -1
        f1 = y[0]
        xx1 = x[0]

    if a2 < x[n - 1]:
        i2 = locate(x, a2)
        f2 = _interp1d_single_loglog(x[i2], y[i2], x[i2 + 1], y[i2 + 1], a2)
        xx2 = a2
    else:
        i2 = n - 1
        f2 = y[n - 1]
        xx2 = x[n - 1]

    if i2 > i1:
        total = 0.0
        if i2 > i1 + 1:
            total = integral_loglog(x[i1 + 1 : i2 + 1], y[i1 + 1 : i2 + 1])
        total += float(_trapezium_loglog(xx1, f1, x[i1 + 1], y[i1 + 1]))
        total += float(_trapezium_loglog(x[i2], y[i2], xx2, f2))
        return total
    return float(_trapezium_loglog(xx1, f1, xx2, f2))


def logspace(wmin: float, wmax: float, n: int) -> np.ndarray:
    """10**linspace(log10(wmin), log10(wmax), n), matching upstream
    logspace(wav_min, wav_max, wavelengths) in main.f90, including its
    guard: n == 1 with wmin != wmax is rejected loudly instead of silently
    collapsing the grid to [wmin]."""
    if n == 1 and wmin != wmax:
        raise ValueError(
            f"cannot build a 1-point log grid with wmin != wmax "
            f"({wmin} != {wmax}) -- likely a misconfigured wavelength range"
        )
    return 10.0 ** np.linspace(np.log10(wmin), np.log10(wmax), n)
