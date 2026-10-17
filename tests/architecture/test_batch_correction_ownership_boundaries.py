from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest
from import_graph_support import build_import_graph

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
CAPABILITY_PREFIX = "phospy.preprocessing.batch_correction"


def test_batch_correction_capability_has_no_dataset_or_core_dependency() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    forbidden_prefixes = (
        "phospy.science.datasets",
        "phospy.workflows",
    )
    offenders = [
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith(CAPABILITY_PREFIX)
        if record.target.startswith(forbidden_prefixes)
    ]

    assert offenders == []


def test_generic_provenance_has_no_batch_correction_capability_dependency() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    offenders = [
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.provenance")
        if record.target.startswith(CAPABILITY_PREFIX)
    ]

    assert offenders == []


def test_retained_core_stage_is_only_an_adapter_and_contract() -> None:
    path = PACKAGE_ROOT / "preprocessing" / "stages" / "batch_correction.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    public_definitions = {
        node.name
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
        and not node.name.startswith("_")
    }

    assert public_definitions == {"BatchCorrectionStage"}
    assert not any(
        isinstance(node, ast.ClassDef)
        and node.name.endswith(("Engine", "Executor", "Provenance"))
        for node in tree.body
    )


def test_core_output_module_does_not_forward_the_passive_output_type() -> None:
    import phospy.preprocessing.execution.corrected_output as core_output
    from phospy.preprocessing.batch_correction.contracts.output import (
        CorrectedPreprocessingOutput,
    )

    assert not hasattr(core_output, "CorrectedPreprocessingOutput")
    assert CorrectedPreprocessingOutput.__module__.endswith(".contracts.output")


@pytest.mark.parametrize(
    "module_order,identity_assertions",
    (
        (
            (
                "phospy.provenance.serialization.batch_correction",
                "phospy.preprocessing.batch_correction.provenance.assembly",
            ),
            "provenance",
        ),
        (
            (
                "phospy.preprocessing.batch_correction.provenance.assembly",
                "phospy.provenance.serialization.batch_correction",
            ),
            "provenance",
        ),
        (
            (
                "phospy.preprocessing.batch_correction.adapters.preprocessing",
                "phospy.preprocessing.stages.batch_correction",
            ),
            "adapter",
        ),
        (
            (
                "phospy.preprocessing.stages.batch_correction",
                "phospy.preprocessing.batch_correction.adapters.preprocessing",
            ),
            "adapter",
        ),
    ),
)
def test_batch_correction_boundaries_import_in_both_orders(
    module_order: tuple[str, str],
    identity_assertions: str,
    tmp_path: Path,
) -> None:
    code = f"""
import importlib
for module_name in {module_order!r}:
    importlib.import_module(module_name)

if {identity_assertions!r} == 'provenance':
    from phospy.provenance.models import BatchCorrectionProvenance as aggregate
    from phospy.provenance.models.batch_correction import BatchCorrectionProvenance as owner
    import phospy.preprocessing.provenance.models as preprocessing_provenance
    assert aggregate is owner
    assert preprocessing_provenance._BatchCorrectionProvenance is owner
    assert preprocessing_provenance.PreprocessingStageProvenance.__annotations__[
        'batch_correction_provenance'
    ] == '_BatchCorrectionProvenance | None'
else:
    from phospy.preprocessing.batch_correction.adapters.preprocessing import (
        SpsRuvStyleBatchCorrectionWorkflowRunner,
    )
    from phospy.preprocessing.stages.batch_correction import (
        BatchCorrectionStage,
    )
    runner = SpsRuvStyleBatchCorrectionWorkflowRunner()
    stage = BatchCorrectionStage(sps_ruv_runner=runner)
    assert callable(runner.run)
    assert stage._executor._sps_ruv_runner is runner
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(SRC_ROOT)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
