"""Benchmark: bhmiepy's dust pipeline vs the pristine upstream CLI.

Unlike the golden tests (which compare against the reference outputs
RECORDED in upstream/examples/), this compiles the upstream code from the
upstream/ git submodule (meson target ``bhmie_ref``, zero source changes),
runs it live on the same parameter files, and compares the outputs at
the upstream file-format precision (~5 significant digits).

Timing: the upstream number is the CLI's end-to-end wall clock (process
start + parameter/ri-table reads + computation + output writes -- that is
upstream's real usage form); the bhmiepy number is a single in-process
compute_dust_properties call. The core-routine-only A/B lives in
bench_upstream.py.

Run from the repo root of an activated environment:
    python benchmarks/bench_dust_upstream.py [--workloads mrn77,custom,kmh94]
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import time
from pathlib import Path

from bhmiepy import dust
from bhmiepy._upstream_cli import (
    compare_to_result,
    find_ref_binary,
    load_outputs,
    prepare_run_dir,
    run_ref,
)

EXAMPLES = Path(__file__).resolve().parents[1] / "upstream" / "examples"

# Dust-shaped single-component workload (same grid shape as mrn77):
# output format 2 so the full file set (angles, albedo, scattering
# matrix) is exercised.
CUSTOM_IN = """\
'bench_dust' = prefix for results
2 = output format (see README)
0.005 = minimum grain size
1.0 = maximum grain size
1000 = number of different grain sizes
181 = number of scattering angles
10 = number of extra fine scattering angles
1 = number of components
141.84 = gas-to-dust ratio
0.01 1000 250 = wavelength parameters
---
1.000 = mass fraction in this component
3.30 = grain density (g/cm^3)
'../ri-data/silicate_ld93' = refractive index file
power = type of size distribution
0.005 0.25 -3.5 = size distribution parameters
---
"""


def stage_custom(workdir: Path):
    """Write the custom workload with the layout the .in file expects."""
    model_dir = workdir / "custom"
    model_dir.mkdir()
    (model_dir / "bench_dust.in").write_text(CUSTOM_IN)
    ri = workdir / "ri-data"
    ri.mkdir()
    shutil.copy2(EXAMPLES / "ri-data" / "silicate_ld93", ri / "silicate_ld93")
    return model_dir, "bench_dust", 2, "bench_dust.in"


def run_workload(binary, name, staged, workdir):
    if name == "custom":
        model_dir, prefix, fmt, param_name = stage_custom(workdir)
    else:
        param_file = EXAMPLES / name / f"{name}.in"
        model_dir, prefix, fmt = prepare_run_dir(param_file, workdir)
        param_name = param_file.name

    t_exe = run_ref(binary, model_dir, param_name)
    outputs = load_outputs(model_dir, prefix, fmt)

    inp = dust.read_parameter_file(model_dir / param_name)
    t0 = time.perf_counter()
    res = dust.compute_dust_properties(
        inp.components, inp.wavelengths, inp.amin, inp.amax, inp.na,
        inp.n_angles, inp.n_small_angles, inp.gas_to_dust,
    )
    t_lib = time.perf_counter() - t0

    scaled = compare_to_result(outputs, res)
    return t_exe, t_lib, scaled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workloads", default="mrn77,custom",
        help="comma-separated subset of mrn77, custom, kmh94 "
             "(kmh94 adds ~minutes of upstream runtime)",
    )
    args = parser.parse_args()

    binary = find_ref_binary()
    print(f"upstream CLI : {binary}")
    print(f"workloads    : {args.workloads}\n")

    failures = []
    for name in [w.strip() for w in args.workloads.split(",") if w.strip()]:
        with tempfile.TemporaryDirectory(prefix="bhmiepy_ab_") as tmp:
            t_exe, t_lib, scaled = run_workload(
                binary, name, None, Path(tmp)
            )
        print(f"[{name}]")
        print(f"  upstream CLI (end-to-end) : {t_exe:8.2f} s")
        print(f"  bhmiepy (library call)    : {t_lib:8.2f} s")
        print(f"  end-to-end speedup        : {t_exe / t_lib:8.2f}x")
        worst = 0.0
        for qty, ratio in sorted(scaled.items()):
            flag = "ok" if ratio <= 1.0 else "FAIL"
            print(f"    {qty:16s} {ratio:8.3f}  {flag}")
            worst = max(worst, ratio)
        if worst > 1.0:
            failures.append(name)
        print()

    if failures:
        raise SystemExit(
            f"FAIL: upstream vs bhmiepy outside tolerance for: {failures}"
        )
    print("OK: all quantities within the golden tolerance "
          "(upstream output precision ~5 significant digits)")


if __name__ == "__main__":
    main()
