"""Slow A/B test: the pristine upstream CLI (bhmie_ref) vs bhmiepy's dust
pipeline, run live on the same parameter file.

Independent of the golden tests (tests/test_dust.py compares against the
RECORDED reference outputs shipped in upstream/examples/): here the upstream
code is compiled from the upstream/ submodule and executed on the fly, which
additionally validates the reference build itself. Both sides are compared
at the upstream output-file precision (~5 significant digits)."""

from pathlib import Path

import pytest

from bhmiepy import dust
from bhmiepy._upstream_cli import (
    compare_to_result,
    load_outputs,
    prepare_run_dir,
    run_ref,
)

EXAMPLES = Path(__file__).resolve().parents[1] / "upstream" / "examples"


@pytest.mark.slow
class TestDustUpstreamEquivalence:
    def test_mrn77_live_cli(self, tmp_path, bhmie_ref):
        model_dir, prefix, fmt = prepare_run_dir(
            EXAMPLES / "mrn77" / "mrn77.in", tmp_path
        )
        run_ref(bhmie_ref, model_dir, "mrn77.in")
        outputs = load_outputs(model_dir, prefix, fmt)

        inp = dust.read_parameter_file(model_dir / "mrn77.in")
        res = dust.compute_dust_properties(
            inp.components, inp.wavelengths, inp.amin, inp.amax, inp.na,
            inp.n_angles, inp.n_small_angles, inp.gas_to_dust,
        )

        scaled = compare_to_result(outputs, res)
        for name, ratio in sorted(scaled.items()):
            assert ratio <= 1.0, (
                f"{name}: exceeds the golden tolerance by {ratio:.2f}x"
            )
