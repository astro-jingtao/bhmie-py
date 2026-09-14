"""Shared tooling for A/B runs against the pristine upstream bhmie CLI.

The upstream reference executable (``bhmie_ref``) is built verbatim from
the read-only upstream/ git submodule (hyperion-rt/bhmie; see meson.build,
option ``upstream_ref``) and speaks the upstream CLI contract: it takes a
parameter-file path as its only argument, resolves the refractive-index
and size-table paths recorded in that file relative to its working
directory, and writes results next to the prefix recorded in the file.

Used by the slow upstream-equivalence test (tests/test_dust_ref.py) and
by benchmarks/bench_dust_upstream.py. Not part of the public API.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

EXE_NAME = "bhmie_ref"
RUN_TIMEOUT = 3600  # seconds; golden-shaped workloads take minutes upstream

# upstream output format 2: file extension -> result attribute name
FORMAT2_FIELDS = {
    "wav": "wavelengths",
    "mu": "mu",
    "alb": "albedo",
    "chi": "kappa_ext",
    "g": "g",
    "f11": "s11",
    "f12": "s12",
    "f33": "s33",
    "f34": "s34",
}

# tolerance shape mirrored from the golden tests (tests/test_dust.py
# _compare): output files carry ~5 significant digits
RTOL = 1e-3
ATOL_SCALE = 1e-4


def find_ref_binary(repo_root: Path | None = None) -> Path:
    """Locate the bhmie_ref executable.

    No platform binding: the binary carries no suffix on POSIX and
    ``.exe`` on Windows; both are probed where relevant. Search order:
    the ``BHMIEPY_REF_BIN`` override, the meson-python editable build
    directory (``<repo>/build``), then ``$PATH``.
    """
    override = os.environ.get("BHMIEPY_REF_BIN")
    if override:
        path = Path(override)
        if not path.is_file():
            raise FileNotFoundError(
                f"BHMIEPY_REF_BIN is set to {override!r}, which does not exist"
            )
        return path

    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[2]
    names = [EXE_NAME] + ([EXE_NAME + ".exe"] if os.name == "nt" else [])
    build_dirs = [repo_root / "build"]
    if build_dirs[0].is_dir():
        # meson-python editable builds live in build/<python-tag>/ (e.g.
        # cp312); accept both layouts
        build_dirs += [d for d in build_dirs[0].iterdir() if d.is_dir()]
    for directory in build_dirs:
        for name in names:
            candidate = directory / name
            if candidate.is_file():
                return candidate

    found = shutil.which(EXE_NAME)
    if found:
        return Path(found)

    raise FileNotFoundError(
        f"{EXE_NAME} (pristine upstream CLI) not found. It is built from the "
        "upstream/ submodule: run 'git submodule update --init --recursive', "
        "then reinstall the package from source "
        "(editable: 'pip install -e . --no-build-isolation'). "
        "Alternatively point BHMIEPY_REF_BIN at the binary, or disable the "
        "target with -Dupstream_ref=disabled."
    )


def _records(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln.strip()]


def _first_token(line: str) -> str:
    return line.split()[0].strip("'\"")


def prepare_run_dir(param_file, dest_root) -> tuple[Path, str, int]:
    """Stage an upstream-CLI run inside ``dest_root``.

    Copies the parameter file and every file it references (refractive
    index tables, size-distribution tables) preserving the relative
    layout the parameter file expects (``../ri-data/...`` resolves
    against the model directory, exactly like the upstream CLI run from
    there). Record indexing mirrors bhmiepy.dust.read_parameter_file.

    Returns ``(model_dir, prefix, output_format)``.
    """
    param_file = Path(param_file)
    dest_root = Path(dest_root)
    records = _records(param_file.read_text())
    prefix = _first_token(records[0])
    output_format = int(_first_token(records[1]))
    n_components = int(_first_token(records[7]))

    refs: list[str] = []
    i = 10
    for _ in range(n_components):
        i += 1                                   # separator record
        i += 2                                   # abundance, density
        refs.append(_first_token(records[i]))    # refractive-index file
        i += 1
        dist_type = _first_token(records[i]).lower()
        i += 1
        if dist_type == "table":
            refs.append(_first_token(records[i]))  # size-table file
        elif dist_type not in ("power", "ped"):
            raise ValueError(f"unknown distribution type {dist_type!r}")
        i += 1                                   # parameters record

    model_dir = dest_root / param_file.parent.name
    model_dir.mkdir(parents=True)
    shutil.copy2(param_file, model_dir / param_file.name)
    for ref in refs:
        src = (param_file.parent / ref).resolve()
        dst = model_dir / ref
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    if "/" in prefix.replace("\\", "/"):
        (model_dir / prefix).parent.mkdir(parents=True, exist_ok=True)
    return model_dir, prefix, output_format


def run_ref(binary, model_dir, param_name, timeout: int = RUN_TIMEOUT) -> float:
    """Run the upstream CLI inside ``model_dir``; return wall-clock seconds.

    Raises on nonzero exit or timeout -- no silent failures.
    """
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(binary), param_name],
        cwd=model_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(
            f"{EXE_NAME} failed (exit {proc.returncode}): "
            f"{proc.stderr[-2000:] or proc.stdout[-2000:]}"
        )
    return elapsed


def load_outputs(model_dir, prefix, output_format: int) -> dict:
    """Load upstream CLI outputs as a dict keyed by bhmiepy result
    attribute names (wavelengths, cext, csca, kappa_ext, g,
    polarization_90 / mu, albedo, s11..s34)."""
    base = str(Path(model_dir) / prefix)
    out: dict = {}
    if output_format == 1:
        data = np.loadtxt(base + ".summary", ndmin=2)
        for name, col in (
            ("wavelengths", 0), ("cext", 1), ("csca", 2),
            ("kappa_ext", 3), ("g", 4), ("polarization_90", 5),
        ):
            out[name] = data[:, col]
    elif output_format == 2:
        for ext, name in FORMAT2_FIELDS.items():
            out[name] = np.loadtxt(f"{base}.{ext}", ndmin=1)
    else:
        raise ValueError(
            f"unsupported output format {output_format}; the A/B uses formats 1 and 2"
        )
    return out


def compare_to_result(outputs: dict, res, rtol: float = RTOL,
                      atol_scale: float = ATOL_SCALE) -> dict:
    """Compare upstream outputs against a compute_dust_properties result.

    Returns per-quantity deviations scaled by the golden-test tolerance
    ``atol + rtol*|want|`` (atol = atol_scale * max|want|): a value <= 1
    is within tolerance, exactly like np.testing.assert_allclose with the
    golden tests' parameters. Shape mismatches map to infinity.
    """
    report: dict = {}
    for name, want in outputs.items():
        got = np.asarray(getattr(res, name), dtype=float)
        want = np.asarray(want, dtype=float)
        if got.shape != want.shape:
            report[name] = float("inf")
            continue
        if want.size == 0:
            report[name] = 0.0
            continue
        atol = atol_scale * float(np.max(np.abs(want)))
        scaled = np.abs(got - want) / (atol + rtol * np.abs(want))
        report[name] = float(np.max(scaled))
    return report
