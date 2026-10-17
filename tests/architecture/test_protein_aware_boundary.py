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


def test_protein_aware_contracts_have_one_owner_and_expected_modules() -> None:
    from phospy.datasets.protein_aware import (
        ProteinAwarePreparationReport,
        ProteinAwarePreparationResult,
    )
    from phospy.preprocessing.results.protein_aware import (
        ProteinAwarePreparationStageOutput,
    )

    assert ProteinAwarePreparationStageOutput.__module__ == (
        "phospy.preprocessing.results.protein_aware"
    )
    assert ProteinAwarePreparationReport.__module__ == ("phospy.datasets.protein_aware")
    assert ProteinAwarePreparationResult.__module__ == ("phospy.datasets.protein_aware")


def test_stage_output_is_absent_from_public_and_obsolete_aggregates() -> None:
    import phospy
    import phospy.advanced as advanced
    import phospy.api as ordinary

    for module in (phospy, ordinary, advanced):
        assert not hasattr(module, "ProteinAwarePreparationStageOutput")


def test_preprocessing_has_no_dataset_or_run_lineage_import() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.preprocessing")
        if record.target.startswith(
            ("phospy.datasets", "phospy.science.datasets", "phospy.run_lineage")
        )
    )

    assert offenders == []


def test_dataset_leaf_has_no_preprocessing_target() -> None:
    path = PACKAGE_ROOT / "datasets" / "protein_aware.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    targets = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        and node.module is not None
        and node.module.startswith("phospy.preprocessing")
    }
    targets.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
        if alias.name.startswith("phospy.preprocessing")
    )

    assert targets == set()


def test_differential_uses_dataset_boundary_without_stage_implementation() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith(
            ("phospy.differential", "phospy.differential")
        )
        if record.target
        in {
            "phospy.preprocessing.results.protein_aware",
            "phospy.preprocessing.stages.protein_aware",
        }
    )

    assert offenders == []


@pytest.mark.parametrize(
    "module_name",
    (
        "phospy.science.datasets.preprocessing.protein_aware_models",
        "phospy.contracts.results.preprocessing",
    ),
)
def test_removed_protein_aware_routes_fail(module_name: str, tmp_path: Path) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(SRC_ROOT)
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            f"import importlib; importlib.import_module({module_name!r})",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "ModuleNotFoundError" in completed.stderr


_RUN_AND_COMPOSE = r"""
import pandas as pd

from phospy.preprocessing.contracts.configuration.total_protein import (
    DatasetProteinAwarePreparationConfig,
)
from phospy.preprocessing.protein_aware.mapping import ProteinMappingConfig
from phospy.preprocessing.stages.protein_aware import ProteinAwarePreparationStage
from phospy.quantitative.models import IntensityScaleState, MatrixIntensityScaleState

phospho = pd.DataFrame(
    {"sample_a": [1.0], "sample_b": [2.0]},
    index=pd.Index(["MAPK14;Y182;"], name="site_id"),
)
site_metadata = pd.DataFrame(
    {"protein_id": ["P53778"], "gene_symbol": ["MAPK14"], "site": ["Y182"]},
    index=phospho.index.copy(),
)
total = pd.DataFrame(
    {"sample_a": [10.0], "sample_b": [11.0]},
    index=pd.Index(["P53778"], name="protein_id"),
)
state = IntensityScaleState(
    phospho=MatrixIntensityScaleState.log2(established_by="protein-aware.phospho"),
    total=MatrixIntensityScaleState.log2(established_by="protein-aware.total"),
)
stage_output = ProteinAwarePreparationStage().run(
    phospho=phospho,
    site_metadata=site_metadata,
    total=total,
    transformation_state=state,
    config=DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs"),
    mapping_config=ProteinMappingConfig(protein_identifier_columns=("protein_id",)),
)
assert stage_output is not None
"""


@pytest.mark.parametrize("dataset_first", (False, True))
def test_protein_aware_boundary_imports_and_composes_in_both_orders(
    dataset_first: bool,
    tmp_path: Path,
) -> None:
    imports = ""
    stage_module_assertion = ""
    if dataset_first:
        imports = """
from phospy.run_lineage.models import RunProvenance
from phospy.datasets.construction.protein_aware import compose_protein_aware_preparation_result
from phospy.datasets.internal_view import DatasetInternalView, ProteinAwarePreparationInternalView
from phospy.datasets.protein_aware import ProteinAwarePreparationReport, ProteinAwarePreparationResult
"""
    else:
        stage_module_assertion = """
import sys
assert not any(name == 'phospy.run_lineage' or name.startswith('phospy.run_lineage.') for name in sys.modules)
assert not any(name == 'phospy.science.datasets' or name.startswith('phospy.science.datasets.') for name in sys.modules)
"""
    late_imports = """
from phospy.run_lineage.models import RunProvenance
from phospy.datasets.construction.protein_aware import compose_protein_aware_preparation_result
from phospy.datasets.internal_view import DatasetInternalView, ProteinAwarePreparationInternalView
from phospy.datasets.protein_aware import ProteinAwarePreparationReport, ProteinAwarePreparationResult
"""
    code = (
        imports
        + _RUN_AND_COMPOSE
        + stage_module_assertion
        + late_imports
        + """
result = compose_protein_aware_preparation_result(stage_output, provenance=None)
assert result.__class__ is ProteinAwarePreparationResult
assert result.report.__class__ is ProteinAwarePreparationReport
assert result.report.provenance is None
assert result.to_payload()["report"]["schema_version"] == 1
protein_view = ProteinAwarePreparationInternalView(result)
assert protein_view.matched_pairs.shape[0] == 1
assert RunProvenance.__module__ == "phospy.run_lineage.models"
assert DatasetInternalView.__module__ == "phospy.datasets.internal_view"
"""
    )
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
