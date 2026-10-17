from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest
from import_graph_support import build_import_graph
from logical_ownership_policy import PREPARATION_RECORD_ALLOWED_SYMBOLS

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"

_OLD_DATASET_MODULES = (
    "phospy.contracts.dataset_build",
    "phospy.preprocessing.validation.sample_groups",
    "phospy.provenance.derived_quantitative",
    "phospy.provenance.models.trusted_assertions",
    "phospy.science.datasets._processing_state.json_contracts",
    "phospy.science.datasets._processing_state.missing_data",
    "phospy.science.datasets._processing_state.models",
    "phospy.science.datasets._processing_state.total_protein",
    "phospy.science.datasets.construction.analysis_ready",
    "phospy.science.datasets.construction.fingerprints",
    "phospy.science.datasets.construction.protein_aware_binding",
    "phospy.science.datasets.construction.service",
    "phospy.science.datasets.construction.trusted_assertions",
    "phospy.science.datasets.construction.validation",
    "phospy.science.datasets.derived_quantitative",
    "phospy.science.datasets.direct_construction",
    "phospy.science.datasets.imputation_metadata",
    "phospy.science.datasets.internal_frame_store",
    "phospy.science.datasets.internal_view",
    "phospy.science.datasets.models",
    "phospy.science.datasets.organism_coherence",
    "phospy.science.datasets.processing_state",
    "phospy.science.datasets.protein_aware",
    "phospy.science.tables.datasets",
    "phospy.science.datasets.builders.contracts",
    "phospy.science.datasets.builders.executor",
    "phospy.science.datasets.builders.interpreter",
    "phospy.science.datasets.builders.interpreter_collaborators",
    "phospy.science.datasets.builders.normalization_reporter",
    "phospy.science.datasets.builders.normalizer",
    "phospy.science.datasets.builders.preprocessing",
    "phospy.science.datasets.builders.protein_aware_binding",
    "phospy.science.datasets.builders.provenance_assembler",
    "phospy.science.datasets.builders.public",
    "phospy.science.datasets.builders.reader",
    "phospy.science.datasets.builders.report_assembler",
    "phospy.science.datasets.builders.sample_metadata_normalizer",
    "phospy.science.datasets.builders.sequence_derivation",
    "phospy.science.datasets.builders.site_identity_derivation",
    "phospy.science.datasets.builders.site_metadata_normalizer",
    "phospy.science.datasets.builders.site_sequence_boundary",
    "phospy.science.datasets.builders.state_builder",
    "phospy.science.datasets.builders.total_matrix_normalizer",
    "phospy.science.datasets.builders.transformation_resolver",
    "phospy.science.datasets.builders.transformation_state",
    "phospy.validation.datasets.builder_request",
    "phospy.validation.datasets.inputs",
    "phospy.workflows.result_caveat_helpers",
)


def _graph():
    return build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)


def _top_level_declarations(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    declarations: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            declarations.add(node.name)
        elif isinstance(node, ast.Assign):
            declarations.update(
                target.id for target in node.targets if isinstance(target, ast.Name)
            )
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            declarations.add(node.target.id)
    declarations.discard("__all__")
    return declarations


def _is_dataset_model_or_durable_module(module_name: str) -> bool:
    if module_name == "phospy.datasets":
        return True
    if not module_name.startswith("phospy.datasets."):
        return False
    lifecycle_roots = (
        "phospy.datasets.builder",
        "phospy.datasets.contracts",
        "phospy.datasets.construction",
        "phospy.datasets.validation",
    )
    return not any(
        module_name == root or module_name.startswith(f"{root}.")
        for root in lifecycle_roots
    )


def test_dataset_model_and_durable_roles_do_not_import_preprocessing_or_builders() -> (
    None
):
    graph = _graph()
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if _is_dataset_model_or_durable_module(record.source_module)
        if record.target.startswith("phospy.preprocessing")
        or record.target.startswith("phospy.science.datasets.builders")
        or record.target in {"phospy", "phospy.api", "phospy.advanced"}
    )

    assert offenders == []


def test_preprocessing_has_no_dataset_or_run_lineage_dependency() -> None:
    graph = _graph()
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module == "phospy.preprocessing"
        or record.source_module.startswith("phospy.preprocessing.")
        if record.target == "phospy.datasets"
        or record.target.startswith("phospy.datasets.")
        or record.target == "phospy.run_lineage"
        or record.target.startswith("phospy.run_lineage.")
    )

    assert offenders == []


def test_preparation_records_are_passive_closed_leaves_with_exact_declarations() -> (
    None
):
    graph = _graph()
    actual_targets = {
        record.target
        for record in graph.records
        if record.source_module == "phospy.preparation_records"
        or record.source_module.startswith("phospy.preparation_records.")
    }
    assert actual_targets == {
        "phospy.errors.input",
        "phospy.errors.validation",
        "phospy.frames.comparison",
        "phospy.frames.ownership",
        "phospy.frames.validation",
        "phospy.policies",
        "phospy.provenance.hashing",
    }

    for module_name, expected in PREPARATION_RECORD_ALLOWED_SYMBOLS.items():
        path = SRC_ROOT.joinpath(*module_name.split(".")).with_suffix(".py")
        assert _top_level_declarations(path) == set(expected)

    for relative in (
        Path("preparation_records/__init__.py"),
        Path("preparation_records/validation/__init__.py"),
    ):
        tree = ast.parse((PACKAGE_ROOT / relative).read_text(encoding="utf-8"))
        assert not any(
            isinstance(node, (ast.Import, ast.ImportFrom)) for node in tree.body
        )


def test_dataset_run_lineage_composition_is_published_and_one_way() -> None:
    graph = _graph()
    dataset_to_lineage = {
        record.target
        for record in graph.records
        if record.source_module.startswith("phospy.datasets.")
        if record.target == "phospy.run_lineage"
        or record.target.startswith("phospy.run_lineage.")
    }
    lineage_to_dataset = {
        record.target
        for record in graph.records
        if record.source_module == "phospy.run_lineage"
        or record.source_module.startswith("phospy.run_lineage.")
        if record.target == "phospy.datasets"
        or record.target.startswith("phospy.datasets.")
    }

    assert dataset_to_lineage == {
        "phospy.run_lineage.models",
        "phospy.run_lineage.serialization",
    }
    assert lineage_to_dataset == set()
    assert not any(
        record.source_module in {"phospy.datasets", "phospy.science.datasets"}
        and record.target.startswith("phospy.run_lineage")
        for record in graph.records
    )


def test_builder_composition_uses_only_published_dataset_and_preprocessing_paths() -> (
    None
):
    graph = _graph()
    composition_records = [
        record
        for record in graph.records
        if record.source_module == "phospy.datasets.builder"
        or record.source_module.startswith("phospy.datasets.construction.")
    ]
    assert any(
        record.target == "phospy.datasets.models" for record in composition_records
    )
    assert any(
        record.target.startswith("phospy.preprocessing.")
        for record in composition_records
    )
    assert not any(
        record.target.startswith("phospy.datasets._internal")
        for record in composition_records
    )


def test_reader_edge_is_exact_and_symbol_scoped() -> None:
    graph = _graph()
    dataset_to_io = [
        record
        for record in graph.records
        if record.source_module == "phospy.datasets"
        or record.source_module.startswith("phospy.datasets.")
        if record.target == "phospy.io" or record.target.startswith("phospy.io.")
    ]
    assert {(record.source_module, record.target) for record in dataset_to_io} == {
        ("phospy.datasets.builder", "phospy.io.readers.dataset_inputs")
    }

    builder_path = PACKAGE_ROOT / "datasets/builder.py"
    tree = ast.parse(
        builder_path.read_text(encoding="utf-8"), filename=str(builder_path)
    )
    imports = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        if node.module == "phospy.io.readers.dataset_inputs"
    ]
    assert len(imports) == 1
    assert [(alias.name, alias.asname) for alias in imports[0].names] == [
        ("DatasetPathTableReader", None),
    ]


def test_reader_initialization_closure_is_dataset_free() -> None:
    graph = _graph()
    eager_edges: dict[str, set[str]] = {}
    for module_name, path in graph.module_paths.items():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        targets: set[str] = set()
        for node in tree.body:
            imported_names: list[str] = []
            if isinstance(node, ast.Import):
                imported_names.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                imported_names.append(node.module)
                imported_names.extend(
                    f"{node.module}.{alias.name}" for alias in node.names
                )
            for imported_name in imported_names:
                parts = imported_name.split(".")
                for index in range(len(parts), 0, -1):
                    candidate = ".".join(parts[:index])
                    if candidate in graph.modules:
                        targets.add(candidate)
                        break
        eager_edges[module_name] = targets

    roots = {
        "phospy.io",
        "phospy.io.readers",
        "phospy.io.readers.dataset_inputs",
    }
    closure = set(roots)
    pending = list(roots)
    while pending:
        source = pending.pop()
        for target in eager_edges.get(source, set()):
            if target in closure:
                continue
            closure.add(target)
            pending.append(target)

    assert not any(
        module == "phospy.datasets" or module.startswith("phospy.datasets.")
        for module in closure
    )


@pytest.mark.parametrize(
    "initializer",
    (
        "phospy.datasets",
        "phospy.datasets.construction",
        "phospy.datasets.validation",
    ),
)
def test_dataset_initializers_do_not_eagerly_load_builder(initializer: str) -> None:
    code = f"""
import importlib
import sys

importlib.import_module({initializer!r})
assert "phospy.datasets.builder" not in sys.modules
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(SRC_ROOT)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize("reader_first", (True, False))
def test_builder_reader_and_facade_import_orders_are_stable(
    reader_first: bool,
    tmp_path: Path,
) -> None:
    first, second = (
        ("phospy.io.readers.dataset_inputs", "phospy.datasets.builder")
        if reader_first
        else ("phospy.datasets.builder", "phospy.io.readers.dataset_inputs")
    )
    code = f"""
import importlib
import sys

first = importlib.import_module({first!r})
if {reader_first!r}:
    assert not any(
        name == "phospy.datasets" or name.startswith("phospy.datasets.")
        for name in sys.modules
    )
second = importlib.import_module({second!r})
owner = importlib.import_module("phospy.datasets.builder")
api = importlib.import_module("phospy.api")
root = importlib.import_module("phospy")
assert importlib.import_module({first!r}) is first
assert importlib.import_module({second!r}) is second
assert root.AnalysisReadyDatasetBuilder is api.AnalysisReadyDatasetBuilder
assert api.AnalysisReadyDatasetBuilder is owner.AnalysisReadyDatasetBuilder
assert owner.AnalysisReadyDatasetBuilder.__module__ == "phospy.datasets.builder"
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


def test_io_reconstructs_through_trusted_factory_without_builder_or_internal_view() -> (
    None
):
    graph = _graph()
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module == "phospy.io"
        or record.source_module.startswith("phospy.io.")
        if record.target.startswith("phospy.science.datasets.builders")
        or record.target.startswith("phospy.datasets._construction")
        or record.target.startswith("phospy.datasets._internal")
        or record.target == "phospy.datasets.internal_view"
    )
    assert offenders == []

    kinase = (PACKAGE_ROOT / "io/bundles/_kinase/reconstruction.py").read_text(
        encoding="utf-8"
    )
    signalome = (PACKAGE_ROOT / "io/bundles/_signalome/reconstruction.py").read_text(
        encoding="utf-8"
    )
    assert "AnalysisReadyPhosphoDataset.from_trusted_tables(" in kinase
    assert "AnalysisReadyPhosphoDataset.from_trusted_tables(" in signalome
    assert "dataset.site_metadata" in signalome
    assert "DatasetInternalView" not in signalome


def test_removed_dataset_routes_fail_in_fresh_process(tmp_path: Path) -> None:
    code = f"""
import importlib

for module_name in {_OLD_DATASET_MODULES!r}:
    try:
        importlib.import_module(module_name)
    except ModuleNotFoundError:
        continue
    raise AssertionError(f"obsolete Dataset route resolved: {{module_name}}")

for symbol_name in ("DatasetBuildRequest", "DatasetInput"):
    try:
        exec(f"from phospy.contracts.requests import {{symbol_name}}")
    except ImportError:
        continue
    raise AssertionError(f"obsolete mixed-request route resolved: {{symbol_name}}")
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


@pytest.mark.parametrize("dataset_first", (True, False))
def test_dataset_and_preprocessing_fresh_import_orders_preserve_identity(
    dataset_first: bool,
    tmp_path: Path,
) -> None:
    first, second = (
        ("phospy.datasets.models", "phospy.preprocessing")
        if dataset_first
        else ("phospy.preprocessing", "phospy.datasets.models")
    )
    code = f"""
import importlib

first = importlib.import_module({first!r})
second = importlib.import_module({second!r})
models = importlib.import_module("phospy.datasets.models")
preprocessing = importlib.import_module("phospy.preprocessing")
assert importlib.import_module({first!r}) is first
assert importlib.import_module({second!r}) is second
from phospy.api import AnalysisReadyPhosphoDataset as public_model
from phospy.datasets.models import AnalysisReadyPhosphoDataset as owner_model
from phospy.preparation_records.batch_correction import BatchCorrectionReport as owner_report
from phospy.preprocessing.batch_correction.contracts import BatchCorrectionReport as consumer_report
assert public_model is owner_model
assert consumer_report is owner_report
assert owner_model.__module__ == "phospy.datasets.models"
assert owner_report.__module__ == "phospy.preparation_records.batch_correction"
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
