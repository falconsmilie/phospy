from __future__ import annotations

import ast
import csv
import json
import os
import subprocess
import sys
from collections.abc import Iterable
from pathlib import Path

import pytest
from import_graph_support import (
    ImportGraph,
    ImportRecord,
    build_eager_import_graph,
    build_import_graph,
)
from import_graph_support import (
    imported_modules as _imported_modules,
)
from import_graph_support import (
    package_name as _package_name,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"


def _kinase_owned_modules() -> frozenset[str]:
    owner_map = PROJECT_ROOT / "architecture/module-owners.csv"
    with owner_map.open(encoding="utf-8", newline="") as handle:
        return frozenset(
            row["module"] for row in csv.DictReader(handle) if row["owner"] == "kinase"
        )


_KINASE_OWNED_MODULES = _kinase_owned_modules()

ALLOWED_PACKAGE_EDGES = frozenset(
    {
        ("phospy.activities", "phospy.errors"),
        ("phospy.activities", "phospy.frames"),
        ("phospy.activities", "phospy.quantitative"),
        ("phospy.advanced", "phospy._api_inventory"),
        ("phospy.advanced", "phospy.differential"),
        ("phospy.advanced", "phospy.evidence"),
        ("phospy.advanced", "phospy.enrichment"),
        ("phospy.advanced", "phospy.io"),
        ("phospy.advanced", "phospy.kinase"),
        ("phospy.advanced", "phospy.preprocessing"),
        ("phospy.advanced", "phospy.references"),
        ("phospy.advanced", "phospy.signalome"),
        ("phospy.advanced", "phospy.statistics"),
        ("phospy.api", "phospy._api_inventory"),
        ("phospy.api", "phospy.io"),
        ("phospy.api", "phospy.datasets"),
        ("phospy.api", "phospy.differential"),
        ("phospy.api", "phospy.errors"),
        ("phospy.api", "phospy.enrichment"),
        ("phospy.api", "phospy.kinase"),
        ("phospy.api", "phospy.preprocessing"),
        ("phospy.api", "phospy.provenance"),
        ("phospy.api", "phospy.quantitative"),
        ("phospy.api", "phospy.references"),
        ("phospy.api", "phospy.taxonomy"),
        ("phospy.api", "phospy.signalome"),
        ("phospy.frames", "phospy.errors"),
        ("phospy.evidence", "phospy.errors"),
        ("phospy.evidence", "phospy.frames"),
        ("phospy.evidence", "phospy.phosphosites"),
        ("phospy.evidence", "phospy.policies"),
        ("phospy.evidence", "phospy.quantitative"),
        ("phospy.enrichment", "phospy.errors"),
        ("phospy.enrichment", "phospy.frames"),
        ("phospy.enrichment", "phospy.provenance"),
        ("phospy.enrichment", "phospy.run_lineage"),
        ("phospy.enrichment", "phospy.quantitative"),
        ("phospy.enrichment", "phospy.statistics"),
        ("phospy.enrichment", "phospy.validation"),
        ("phospy.datasets", "phospy.errors"),
        ("phospy.datasets", "phospy.evidence"),
        ("phospy.datasets", "phospy.frames"),
        ("phospy.datasets", "phospy.phosphosites"),
        ("phospy.datasets", "phospy.preparation_records"),
        ("phospy.datasets", "phospy.provenance"),
        ("phospy.datasets", "phospy.preprocessing"),
        ("phospy.datasets", "phospy.quantitative"),
        ("phospy.datasets", "phospy.references"),
        ("phospy.datasets", "phospy.run_lineage"),
        ("phospy.datasets", "phospy.taxonomy"),
        ("phospy.differential", "phospy.datasets"),
        ("phospy.differential", "phospy.errors"),
        ("phospy.differential", "phospy.evidence"),
        ("phospy.differential", "phospy.frames"),
        ("phospy.differential", "phospy.phosphosites"),
        ("phospy.differential", "phospy.policies"),
        ("phospy.differential", "phospy.preparation_records"),
        ("phospy.differential", "phospy.provenance"),
        ("phospy.differential", "phospy.quantitative"),
        ("phospy.differential", "phospy.references"),
        ("phospy.differential", "phospy.run_lineage"),
        ("phospy.differential", "phospy.statistics"),
        ("phospy.differential", "phospy.validation"),
        ("phospy.io", "phospy.datasets"),
        ("phospy.io", "phospy.enrichment"),
        ("phospy.io", "phospy.errors"),
        ("phospy.io", "phospy.evidence"),
        ("phospy.io", "phospy.frames"),
        ("phospy.io", "phospy.kinase"),
        ("phospy.io", "phospy.phosphosites"),
        ("phospy.io", "phospy.preprocessing"),
        ("phospy.io", "phospy.preparation_records"),
        ("phospy.io", "phospy.provenance"),
        ("phospy.io", "phospy.quantitative"),
        ("phospy.io", "phospy.references"),
        ("phospy.io", "phospy.run_lineage"),
        ("phospy.io", "phospy.taxonomy"),
        ("phospy.io", "phospy.validation"),
        ("phospy.io", "phospy.activities"),
        ("phospy.io", "phospy.signalome"),
        ("phospy.policies", "phospy.errors"),
        ("phospy.preparation_records", "phospy.errors"),
        ("phospy.preparation_records", "phospy.frames"),
        ("phospy.preparation_records", "phospy.policies"),
        ("phospy.preparation_records", "phospy.provenance"),
        ("phospy.kinase", "phospy.errors"),
        ("phospy.kinase", "phospy.activities"),
        ("phospy.kinase", "phospy.datasets"),
        ("phospy.kinase", "phospy.evidence"),
        ("phospy.kinase", "phospy.frames"),
        ("phospy.kinase", "phospy.phosphosites"),
        ("phospy.kinase", "phospy.policies"),
        ("phospy.kinase", "phospy.provenance"),
        ("phospy.kinase", "phospy.quantitative"),
        ("phospy.kinase", "phospy.references"),
        ("phospy.kinase", "phospy.run_lineage"),
        ("phospy.kinase", "phospy.taxonomy"),
        ("phospy.kinase", "phospy.statistics"),
        ("phospy.provenance", "phospy.errors"),
        ("phospy.preprocessing", "phospy.errors"),
        ("phospy.preprocessing", "phospy.evidence"),
        ("phospy.preprocessing", "phospy.frames"),
        ("phospy.preprocessing", "phospy.phosphosites"),
        ("phospy.preprocessing", "phospy.policies"),
        ("phospy.preprocessing", "phospy.provenance"),
        ("phospy.preprocessing", "phospy.preparation_records"),
        ("phospy.preprocessing", "phospy.quantitative"),
        ("phospy.preprocessing", "phospy.references"),
        ("phospy.preprocessing", "phospy.taxonomy"),
        ("phospy.phosphosites", "phospy.errors"),
        ("phospy.phosphosites", "phospy.frames"),
        ("phospy.phosphosites", "phospy.taxonomy"),
        ("phospy.quantitative", "phospy.errors"),
        ("phospy.quantitative", "phospy.provenance"),
        ("phospy.references", "phospy.errors"),
        ("phospy.references", "phospy.frames"),
        ("phospy.references", "phospy.phosphosites"),
        ("phospy.references", "phospy.policies"),
        ("phospy.references", "phospy.provenance"),
        ("phospy.references", "phospy.taxonomy"),
        ("phospy.run_lineage", "phospy.errors"),
        ("phospy.run_lineage", "phospy.evidence"),
        ("phospy.run_lineage", "phospy.preprocessing"),
        ("phospy.run_lineage", "phospy.provenance"),
        ("phospy.run_lineage", "phospy.references"),
        ("phospy.signalome", "phospy.datasets"),
        ("phospy.signalome", "phospy.errors"),
        ("phospy.signalome", "phospy.evidence"),
        ("phospy.signalome", "phospy.frames"),
        ("phospy.signalome", "phospy.kinase"),
        ("phospy.signalome", "phospy.phosphosites"),
        ("phospy.signalome", "phospy.policies"),
        ("phospy.signalome", "phospy.provenance"),
        ("phospy.signalome", "phospy.quantitative"),
        ("phospy.signalome", "phospy.references"),
        ("phospy.signalome", "phospy.run_lineage"),
        ("phospy.signalome", "phospy.validation"),
        ("phospy.validation", "phospy.errors"),
        ("phospy.validation", "phospy.frames"),
    }
)


def test_sps_discovery_request_is_passive_at_the_public_boundary() -> None:
    source = (
        PACKAGE_ROOT / "preprocessing/control_sites/discovery_contracts.py"
    ).read_text(encoding="utf-8")
    tree = ast.parse(source)
    request_class = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "SpsDiscoveryRequest"
    )

    assert not any(
        isinstance(node, ast.FunctionDef) and node.name == "__post_init__"
        for node in request_class.body
    )


def test_advanced_workflow_dependency_is_limited_to_sps_public_boundary() -> None:
    graph = _build_import_graph()
    advanced_workflow_targets = {
        record.target
        for record in graph.records
        if record.source_module.startswith("phospy.advanced")
        if record.target.startswith("phospy.workflows")
    }

    assert advanced_workflow_targets == set()


FORBIDDEN_PACKAGE_EDGES = frozenset(
    {
        ("phospy.errors", "phospy.api"),
        ("phospy.errors", "phospy.contracts"),
        ("phospy.errors", "phospy.io"),
        ("phospy.errors", "phospy.science"),
        ("phospy.errors", "phospy.tables"),
        ("phospy.errors", "phospy.validation"),
        ("phospy.errors", "phospy.workflows"),
        ("phospy.io", "phospy.api"),
        ("phospy.phosphosites", "phospy.contracts"),
        ("phospy.phosphosites", "phospy.io"),
        ("phospy.phosphosites", "phospy.science"),
        ("phospy.phosphosites", "phospy.workflows"),
        ("phospy.references", "phospy.evidence"),
        ("phospy.references", "phospy.io"),
        ("phospy.references", "phospy.science"),
        ("phospy.evidence", "phospy.references"),
        ("phospy.science", "phospy.contracts"),
        ("phospy.science", "phospy.io"),
        ("phospy.science", "phospy.tables"),
        ("phospy.science", "phospy.validation"),
        ("phospy.science", "phospy.workflows"),
        ("phospy.validation", "phospy.workflows"),
    }
)


def test_package_dependency_graph_is_acyclic() -> None:
    graph = _build_import_graph()
    package_edges = _architectural_package_edges(graph)
    cycles = [
        tuple(sorted(component))
        for component in _strongly_connected_components(
            _package_adjacency(package_edges)
        )
        if len(component) > 1
    ]

    assert cycles == [], _format_cycles(cycles)


def test_complete_eager_module_graph_is_acyclic_without_classification() -> None:
    graph = build_eager_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    cycles = [
        tuple(sorted(component))
        for component in _strongly_connected_components(graph.module_edges)
        if len(component) > 1
    ]

    assert cycles == [], _format_cycles(cycles)


def test_governed_condensation_classifications_match_source() -> None:
    policy = json.loads(
        (
            PROJECT_ROOT / "architecture/permanent-condensation-classifications.json"
        ).read_text(encoding="utf-8")
    )
    observed = tuple(
        (
            item["classification_id"],
            item["source_path"],
            item["source_module"],
            item["target_module"],
            item["symbol"],
        )
        for item in policy["classifications"]
    )

    assert observed == _GOVERNED_PACKAGE_CONDENSATION_CLASSIFICATIONS
    for classification in observed:
        _classification_id, source_path, source_module, target_module, symbol = (
            classification
        )
        source = (PROJECT_ROOT / source_path).read_text(encoding="utf-8")
        assert (
            _governed_import_form_problems(
                source,
                source_module=source_module,
                target_module=target_module,
                symbol=symbol,
            )
            == []
        )


@pytest.mark.parametrize(
    "source",
    (
        "from phospy.validation.common.input_sources import "
        "validate_dataframe_or_path as validate_source\n",
        "from phospy.validation.common.input_sources import *\n",
        "from importlib import import_module\n"
        "validate_source = import_module("
        "'phospy.validation.common.input_sources').validate_dataframe_or_path\n",
        "register('phospy.validation.common.input_sources.validate_dataframe_or_path')\n",
        "from phospy.validation.common.input_sources import "
        "validate_dataframe_or_path\n"
        "from phospy.validation.common.input_sources import "
        "validate_dataframe_or_path\n",
    ),
    ids=("alias", "wildcard", "dynamic", "registration", "duplicate"),
)
def test_governed_condensation_classification_rejects_non_exact_forms(
    source: str,
) -> None:
    assert _governed_import_form_problems(
        source,
        source_module="phospy.datasets.validation.inputs",
        target_module="phospy.validation.common.input_sources",
        symbol="validate_dataframe_or_path",
    )


def test_package_dependency_graph_matches_allowed_edge_table() -> None:
    graph = _build_import_graph()
    package_edges = _architectural_package_edges(graph)

    unexpected = sorted(package_edges - ALLOWED_PACKAGE_EDGES)
    stale_allowed = sorted(ALLOWED_PACKAGE_EDGES - package_edges)

    assert unexpected == [], _format_package_edges(unexpected, graph.records)
    assert stale_allowed == [], "\n".join(f"{s} -> {t}" for s, t in stale_allowed)


def test_forbidden_package_edges_are_absent() -> None:
    graph = _build_import_graph()
    offenders = sorted(_architectural_package_edges(graph) & FORBIDDEN_PACKAGE_EDGES)

    assert offenders == [], _format_package_edges(offenders, graph.records)


def test_import_extractor_includes_all_static_import_forms() -> None:
    source = """
from typing import TYPE_CHECKING
import phospy.datasets.models
from phospy.datasets.validation import inputs
from .datasets import validation
import importlib
importlib.import_module("phospy.io.bundles")
__import__("phospy.kinase.workflow")
if TYPE_CHECKING:
    from phospy.references.models import ReferenceBundle
def _local_imports():
    import phospy.preprocessing.contracts.configuration.dataset
    from .signalome import models
    importlib.import_module("phospy.io.bundles.signalome")
"""
    path = PACKAGE_ROOT / "__init__.py"
    imports = tuple(
        imported
        for imported, _line in _imported_modules("phospy", path, ast.parse(source))
    )

    assert "phospy.datasets.models" in imports
    assert "phospy.datasets.validation" in imports
    assert "phospy.datasets.validation.inputs" in imports
    assert "phospy.datasets" in imports
    assert "phospy.references.models" in imports
    assert "phospy.references.models.ReferenceBundle" in imports
    assert "phospy.preprocessing.contracts.configuration.dataset" in imports
    assert "phospy.signalome" in imports
    assert "phospy.signalome.models" in imports
    assert "phospy.io.bundles" in imports
    assert "phospy.io.bundles.signalome" in imports
    assert "phospy.kinase.workflow" in imports


def test_import_extractor_records_ast_line_numbers() -> None:
    source = """
import phospy.datasets.models
from phospy.datasets.validation import inputs
import importlib
importlib.import_module("phospy.io.bundles")
"""
    path = PACKAGE_ROOT / "__init__.py"
    imports = tuple(_imported_modules("phospy", path, ast.parse(source)))

    assert imports == (
        ("phospy.datasets.models", 2),
        ("phospy.datasets.validation", 3),
        ("phospy.datasets.validation.inputs", 3),
        ("importlib", 4),
        ("phospy.io.bundles", 5),
    )


def test_public_import_routes_load_in_clean_interpreters() -> None:
    code = "\n".join(
        [
            "from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow",
            "from phospy.api import ReferenceBundleBuilder",
            "from phospy.datasets.contracts import DatasetBuildRequest",
            "from phospy.errors import PhosPyError",
            "from phospy.references.models import ReferenceBundle",
            "from phospy.datasets.builder import "
            "AnalysisReadyDatasetBuilder as ScienceDatasetBuilder",
            "from phospy.references.validation.bundle import ReferenceBundleValidator",
            "assert AnalysisReadyDatasetBuilder",
            "assert ReferenceBundleBuilder",
            "assert DatasetBuildRequest",
            "assert PhosPyError",
            "assert ReferenceBundle",
            "assert ScienceDatasetBuilder",
            "assert ReferenceBundleValidator",
            "assert KinaseWorkflow",
        ]
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [
            str(SRC_ROOT),
            environment.get("PYTHONPATH", ""),
        ]
    )

    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "shared_imports",
    [
        ("phospy.quantitative", "phospy.statistics"),
        ("phospy.statistics", "phospy.quantitative"),
    ],
)
def test_shared_quantitative_and_statistics_imports_add_no_capability_or_dataset(
    shared_imports: tuple[str, str],
) -> None:
    code = "\n".join(
        [
            "import importlib",
            "import sys",
            "import phospy",
            "isolated_prefixes = ('phospy.quantitative', 'phospy.statistics', "
            "'phospy.science.activities', 'phospy.science.datasets')",
            "for module_name in tuple(sys.modules):",
            "    if module_name.startswith(isolated_prefixes):",
            "        sys.modules.pop(module_name)",
            *(f"importlib.import_module({name!r})" for name in shared_imports),
            "forbidden = sorted(name for name in sys.modules if "
            "name.startswith(('phospy.science.activities', "
            "'phospy.science.datasets')))",
            "assert forbidden == [], forbidden",
        ]
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(SRC_ROOT), environment.get("PYTHONPATH", "")]
    )

    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "owner_imports",
    [
        ("references", "evidence"),
        ("evidence", "references"),
    ],
)
def test_references_and_evidence_owner_leaves_load_in_both_fresh_process_orders(
    owner_imports: tuple[str, str],
) -> None:
    owner_modules = {
        "references": (
            "phospy.references.builder",
            "phospy.references.models",
            "phospy.references.tables",
            "phospy.references.sequences.repository",
        ),
        "evidence": (
            "phospy.evidence.localisation",
            "phospy.evidence.sequence_resolution",
            "phospy.evidence.quantification_depth",
        ),
    }
    imports = tuple(
        module_name for owner in owner_imports for module_name in owner_modules[owner]
    )
    code = "\n".join(
        [
            "import importlib",
            "import inspect",
            "import pathlib",
            "import tempfile",
            "import pandas as pd",
            *(f"importlib.import_module({name!r})" for name in imports),
            "from phospy.references.builder import ReferenceBundleBuilder",
            "from phospy.references.tables import KinaseSubstrateReference",
            "from phospy.references.sequences.repository import "
            "FastaProteinSequenceRepository",
            "from phospy.evidence.localisation import "
            "normalise_localisation_confidence_series",
            "from phospy.evidence.quantification_depth import "
            "validate_quantification_depth_series",
            "from phospy.evidence.sequence_resolution import "
            "PhosphositeSequenceResolutionRequest, PhosphositeSequenceResolver",
            "assert tuple(inspect.signature(ReferenceBundleBuilder).parameters) == ()",
            "assert isinstance(ReferenceBundleBuilder(), ReferenceBundleBuilder)",
            "table = KinaseSubstrateReference(frame=pd.DataFrame({"
            "'kinase': ['akt1'], 'substrate_site': ['Mapk1;s123;']}))",
            "assert table.frame.loc[0, 'kinase'] == 'AKT1'",
            "assert table.frame.loc[0, 'substrate_site'] == 'MAPK1;S123;'",
            "with tempfile.TemporaryDirectory() as tmp_dir:",
            "    fasta = pathlib.Path(tmp_dir) / 'proteins.fasta'",
            "    fasta.write_text('>P1\\nAASAA\\n', encoding='utf-8')",
            "    repository = FastaProteinSequenceRepository.from_path(fasta)",
            "    resolution = PhosphositeSequenceResolver().run("
            "PhosphositeSequenceResolutionRequest('P1', 'S3', 2), repository)",
            "    assert resolution.status == 'resolved'",
            "    assert resolution.site_sequence == 'AASAA'",
            "normalised, report, warnings = normalise_localisation_confidence_series("
            "pd.Series([75.0], index=['site']), source_column='probability', "
            "scale='percent')",
            "assert normalised.tolist() == [0.75]",
            "assert report.to_payload()['invalid_count'] == 0",
            "assert warnings == ()",
            "depth = validate_quantification_depth_series("
            "pd.Series([2.0], index=['site']), field_name='depth')",
            "assert depth.name == 'quantification_depth'",
            "assert depth.tolist() == [2.0]",
        ]
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(SRC_ROOT), environment.get("PYTHONPATH", "")]
    )

    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


def _build_import_graph() -> ImportGraph:
    return build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)


_GOVERNED_PACKAGE_CONDENSATION_CLASSIFICATIONS = (
    (
        "dataset-builder-path-reader",
        "src/phospy/datasets/builder.py",
        "phospy.datasets.builder",
        "phospy.io.readers.dataset_inputs",
        "DatasetPathTableReader",
    ),
    (
        "dataset-input-source-validation",
        "src/phospy/datasets/construction/input_resolution.py",
        "phospy.datasets.construction.input_resolution",
        "phospy.validation.common.input_sources",
        "validate_dataframe_or_path",
    ),
    (
        "dataset-request-source-validation",
        "src/phospy/datasets/validation/inputs.py",
        "phospy.datasets.validation.inputs",
        "phospy.validation.common.input_sources",
        "validate_dataframe_or_path",
    ),
)

_ROLE_SCOPED_IMPORTS_REMOVED_BEFORE_PACKAGE_CONDENSATION = frozenset(
    (source_module, target_module)
    for (
        _classification_id,
        _source_path,
        source_module,
        target_module,
        _symbol,
    ) in _GOVERNED_PACKAGE_CONDENSATION_CLASSIFICATIONS
)


def _governed_import_form_problems(
    source: str,
    *,
    source_module: str,
    target_module: str,
    symbol: str,
) -> list[str]:
    tree = ast.parse(source)
    exact_import_count = 0
    problems: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == target_module:
            names = [(alias.name, alias.asname) for alias in node.names]
            if node.level == 0 and names == [(symbol, None)]:
                exact_import_count += 1
            else:
                problems.append(f"non-exact import form at line {node.lineno}: {names}")
        elif isinstance(node, ast.Import):
            if any(
                alias.name == target_module
                or alias.name.startswith(f"{target_module}.")
                for alias in node.names
            ):
                problems.append(f"module import at line {node.lineno}")
        elif isinstance(node, ast.Call):
            if _is_dynamic_import_call(node, target_module=target_module):
                problems.append(f"dynamic import at line {node.lineno}")
            if _is_registration_call(
                node,
                target_module=target_module,
                symbol=symbol,
            ):
                problems.append(f"registration import at line {node.lineno}")
    if exact_import_count != 1:
        problems.append(
            f"{source_module} requires exactly one direct import of "
            f"{target_module}.{symbol}; observed={exact_import_count}"
        )
    return problems


def _is_dynamic_import_call(node: ast.Call, *, target_module: str) -> bool:
    function_name: str | None = None
    if isinstance(node.func, ast.Name):
        function_name = node.func.id
    elif isinstance(node.func, ast.Attribute):
        function_name = node.func.attr
    if function_name not in {"import_module", "__import__"} or not node.args:
        return False
    target = node.args[0]
    return isinstance(target, ast.Constant) and target.value == target_module


def _is_registration_call(
    node: ast.Call,
    *,
    target_module: str,
    symbol: str,
) -> bool:
    if isinstance(node.func, ast.Name):
        function_name = node.func.id
    elif isinstance(node.func, ast.Attribute):
        function_name = node.func.attr
    else:
        return False
    if "register" not in function_name.lower():
        return False
    governed_values = {
        target_module,
        symbol,
        f"{target_module}.{symbol}",
    }
    arguments = (*node.args, *(keyword.value for keyword in node.keywords))
    return any(
        isinstance(argument, ast.Constant) and argument.value in governed_values
        for argument in arguments
    )


def _architectural_package_edges(
    graph: ImportGraph,
) -> frozenset[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    targets_by_source_line: dict[tuple[str, int], set[str]] = {}
    for record in graph.records:
        targets_by_source_line.setdefault(
            (record.source_module, record.line), set()
        ).add(record.target)
    for record in graph.records:
        if (record.source_module, record.target) in (
            _ROLE_SCOPED_IMPORTS_REMOVED_BEFORE_PACKAGE_CONDENSATION
        ):
            continue
        same_line_targets = targets_by_source_line[(record.source_module, record.line)]
        if any(target.startswith(f"{record.target}.") for target in same_line_targets):
            continue
        if (
            record.source_module in _KINASE_OWNED_MODULES
            and record.target in _KINASE_OWNED_MODULES
        ):
            continue
        source = _architectural_package_name(record.source_module)
        target = _architectural_package_name(record.target)
        if source is not None and target is not None and source != target:
            edges.add((source, target))
    return frozenset(edges)


def _architectural_package_name(module_name: str) -> str | None:
    return _package_name(module_name)


def _package_adjacency(
    package_edges: Iterable[tuple[str, str]],
) -> dict[str, frozenset[str]]:
    adjacency: dict[str, set[str]] = {}
    for source, target in package_edges:
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set())
    return {source: frozenset(targets) for source, targets in adjacency.items()}


def _strongly_connected_components(
    graph: dict[str, frozenset[str]],
) -> list[tuple[str, ...]]:
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    components: list[tuple[str, ...]] = []

    def strong_connect(node: str) -> None:
        nonlocal index
        indices[node] = index
        lowlinks[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)

        for target in graph.get(node, ()):
            if target not in indices:
                strong_connect(target)
                lowlinks[node] = min(lowlinks[node], lowlinks[target])
            elif target in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[target])

        if lowlinks[node] != indices[node]:
            return
        component: list[str] = []
        while True:
            target = stack.pop()
            on_stack.remove(target)
            component.append(target)
            if target == node:
                break
        components.append(tuple(component))

    for node in graph:
        if node not in indices:
            strong_connect(node)
    return components


def _format_cycles(cycles: list[tuple[str, ...]]) -> str:
    if not cycles:
        return ""
    return "\n\n".join("\n".join(cycle) for cycle in cycles)


def _format_package_edges(
    edges: Iterable[tuple[str, str]],
    records: tuple[ImportRecord, ...],
) -> str:
    lines: list[str] = []
    for source, target in edges:
        lines.append(f"{source} -> {target}")
        examples = [
            record
            for record in records
            if _package_name(record.source_module) == source
            and _package_name(record.target) == target
        ][:8]
        for record in examples:
            lines.append(
                f"  {record.source_path.relative_to(PROJECT_ROOT)}:{record.line}"
            )
    return "\n".join(lines)
