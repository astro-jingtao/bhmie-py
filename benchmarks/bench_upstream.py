"""Benchmark: bhmiepy's optimized BHMIE copy vs the pristine upstream code.

Both variants are compiled into the extension:
  - bhmie_vec          -- the package copy (allocate(d(nmx)) workspace,
                          heap-allocated angle arrays; see
                          src/bhmiepy/_fortran/bhmie.f90 header)
  - bhmie_vec_upstream -- pristine upstream body incl. the original
                          allocate(d(nmxx)) 16 MB workspace per call
                          (src/bhmiepy/_fortran/ref/bhmie_upstream.f90)

The script asserts the two produce IDENTICAL results (same algorithm,
only allocation strategy differs) and reports wall-clock times.

Run from an activated bhmiepy environment:
    python benchmarks/bench_upstream.py [--repeats N]
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from bhmiepy import _bhmiepy_ext as ext

NANG = 91


def make_workload():
    # a dust-pipeline-shaped batch: sizes x wavelengths grid
    sizes = np.logspace(np.log10(0.01), np.log10(2.0), 60)
    wavelengths = np.logspace(np.log10(0.1), np.log10(10.0), 40)
    x = (2.0 * np.pi * sizes[:, None] / wavelengths[None, :]).reshape(-1)
    m = np.full_like(x, 1.6 + 0.01j, dtype=np.complex128)
    return x, m


def run_once(func, x, m):
    return func(x, m, NANG)


def best_time(func, x, m, repeats):
    out = run_once(func, x, m)  # warmup + correctness capture
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        run_once(func, x, m)
        best = min(best, time.perf_counter() - t0)
    return best, out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()

    x, m = make_workload()
    print(f"workload: {x.size} (size, wavelength) points, nang={NANG}")

    t_ours, out_ours = best_time(ext.bhmie_vec, x, m, args.repeats)
    t_upstream, out_up = best_time(ext.bhmie_vec_upstream, x, m, args.repeats)

    names = ["qext", "qsca", "qback", "g", "s1", "s2"]
    identical = all(
        np.array_equal(a, b) for a, b in zip(out_ours, out_up)
    )
    for name, a, b in zip(names, out_ours, out_up):
        max_diff = np.max(np.abs(a - b)) if a.size else 0.0
        print(f"  {name:5s}: identical={np.array_equal(a, b)} "
              f"(max |diff| = {max_diff:.3e})")

    print()
    print(f"  bhmiepy (optimized)   : {t_ours:8.3f} s")
    print(f"  upstream (pristine)   : {t_upstream:8.3f} s")
    print(f"  speedup               : {t_upstream / t_ours:8.2f}x")

    if not identical:
        raise SystemExit("FAIL: results differ between optimized and upstream copies")
    print("\nOK: results are identical; performance "
          f"{'matches or exceeds' if t_ours <= t_upstream else 'TRAILS'} "
          "the pristine upstream")


if __name__ == "__main__":
    main()
