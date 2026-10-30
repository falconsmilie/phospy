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

from tools.architecture.import_graph import (
    ImportGraph,
    ImportRecord,
    build_eager_import_graph,
    build_import_graph,
    strongly_connected_components,
)
from tools.architecture.import_graph import (
    imported_modules as _imported_modules,
)
from tools.architecture.logical_ownership import (
    CAPABILITY_OWNERS,
    SHARED_ALLOWED_TARGETS,
    SHARED_OWNERS,
    architectural_package_unit,
    owner_namespace_prefixes,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"


def _module_owners() -> dict[str, str]:
    owner_map = PROJECT_ROOT / "architecture/module-owners.csv"
    with owner_map.open(encoding="utf-8", newline="") as handle:
        return {row["module"]: row["owner"] for row in csv.DictReader(handle)}


_MODULE_OWNERS = _module_owners()
_KINASE_OWNED_MODULES = frozenset(
    module for module, owner in _MODULE_OWNERS.items() if owner == "kinase"
)


def _package_roots_for_owners(owners: Iterable[str]) -> frozenset[str]:
    return owner_namespace_prefixes(
        owners,
        owners_by_module=_MODULE_OWNERS,
    )


def _shared_owner_forbidden_roots(source_owner: str) -> frozenset[str]:
    allowed_owners = SHARED_ALLOWED_TARGETS[source_owner]
    forbidden_shared_owners = SHARED_OWNERS - allowed_owners - {source_owner}
    forbidden_owners = CAPABILITY_OWNERS | forbidden_shared_owners
    forbidden_roots = _package_roots_for_owners(forbidden_owners)
    assert forbidden_roots, f"empty forbidden-owner classification for {source_owner}"
    return forbidden_roots


def _run_shared_owner_import_probe(
    *,
    source_module: str,
    forbidden_roots: frozenset[str],
    injected_module: str | None = None,
) -> tuple[str, ...]:
    assert forbidden_roots, f"empty forbidden roots for {source_module}"
    code = """
import importlib
import json
import sys

source_module = sys.argv[1]
forbidden_roots = tuple(json.loads(sys.argv[2]))
injected_module = sys.argv[3] or None
assert forbidden_roots, f"empty forbidden roots for {source_module}"

importlib.import_module(source_module)
if injected_module is not None:
    importlib.import_module(injected_module)

forbidden = sorted(
    module_name
    for module_name in sys.modules
    if any(
        module_name == root or module_name.startswith(f"{root}.")
        for root in forbidden_roots
    )
)
print(json.dumps(forbidden))
"""
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(SRC_ROOT), environment.get("PYTHONPATH", "")]
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            source_module,
            json.dumps(sorted(forbidden_roots)),
            injected_module or "",
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    return tuple(json.loads(completed.stdout))


ALLOWED_PACKAGE_EDGES = frozenset(
    {
        ("phospy.core.activities", "phospy.core.errors"),
        ("phospy.core.activities", "phospy.core.frames"),
        ("phospy.core.activities", "phospy.core.quantitative"),
        ("phospy.api.advanced", "phospy.api"),
        ("phospy.api.advanced", "phospy.differential"),
        ("phospy.api.advanced", "phospy.core.evidence"),
        ("phospy.api.advanced", "phospy.enrichment"),
        ("phospy.api.advanced", "phospy.io"),
        ("phospy.api.advanced", "phospy.kinase"),
        ("phospy.api.advanced", "phospy.preprocessing"),
        ("phospy.api.advanced", "phospy.core.references"),
        ("phospy.api.advanced", "phospy.signalome"),
        ("phospy.api.advanced", "phospy.core.statistics"),
        ("phospy.api", "phospy.io"),
        ("phospy.api", "phospy.datasets"),
        ("phospy.api", "phospy.differential"),
        ("phospy.api", "phospy.core.errors"),
        ("phospy.api", "phospy.enrichment"),
        ("phospy.api", "phospy.kinase"),
        ("phospy.api", "phospy.preprocessing"),
        ("phospy.api", "phospy.core.provenance"),
        ("phospy.api", "phospy.core.quantitative"),
        ("phospy.api", "phospy.core.references"),
        ("phospy.api", "phospy.core.taxonomy"),
        ("phospy.api", "phospy.signalome"),
        ("phospy.core.frames", "phospy.core.errors"),
        ("phospy.core.evidence", "phospy.core.errors"),
        ("phospy.core.evidence", "phospy.core.frames"),
        ("phospy.core.evidence", "phospy.core.phosphosites"),
        ("phospy.core.evidence", "phospy.core.policies"),
        ("phospy.core.evidence", "phospy.core.quantitative"),
        ("phospy.enrichment", "phospy.core.errors"),
        ("phospy.enrichment", "phospy.core.frames"),
        ("phospy.enrichment", "phospy.core.provenance"),
        ("phospy.enrichment", "phospy.core.run_lineage"),
        ("phospy.enrichment", "phospy.core.quantitative"),
        ("phospy.enrichment", "phospy.core.statistics"),
        ("phospy.datasets", "phospy.core.errors"),
        ("phospy.datasets", "phospy.core.evidence"),
        ("phospy.datasets", "phospy.core.frames"),
        ("phospy.datasets", "phospy.core.phosphosites"),
        ("phospy.datasets", "phospy.core.preparation_records"),
        ("phospy.datasets", "phospy.core.provenance"),
        ("phospy.datasets", "phospy.preprocessing"),
        ("phospy.datasets", "phospy.core.quantitative"),
        ("phospy.datasets", "phospy.core.references"),
        ("phospy.datasets", "phospy.core.run_lineage"),
        ("phospy.datasets", "phospy.core.taxonomy"),
        ("phospy.differential", "phospy.datasets"),
        ("phospy.differential", "phospy.core.errors"),
        ("phospy.differential", "phospy.core.evidence"),
        ("phospy.differential", "phospy.core.frames"),
        ("phospy.differential", "phospy.core.phosphosites"),
        ("phospy.differential", "phospy.core.policies"),
        ("phospy.differential", "phospy.core.preparation_records"),
        ("phospy.differential", "phospy.core.provenance"),
        ("phospy.differential", "phospy.core.quantitative"),
        ("phospy.differential", "phospy.core.references"),
        ("phospy.differential", "phospy.core.run_lineage"),
        ("phospy.differential", "phospy.core.statistics"),
        ("phospy.io", "phospy.datasets"),
        ("phospy.io", "phospy.enrichment"),
        ("phospy.io", "phospy.core.errors"),
        ("phospy.io", "phospy.core.evidence"),
        ("phospy.io", "phospy.core.frames"),
        ("phospy.io", "phospy.kinase"),
        ("phospy.io", "phospy.core.phosphosites"),
        ("phospy.io", "phospy.preprocessing"),
        ("phospy.io", "phospy.core.preparation_records"),
        ("phospy.io", "phospy.core.provenance"),
        ("phospy.io", "phospy.core.quantitative"),
        ("phospy.io", "phospy.core.references"),
        ("phospy.io", "phospy.core.run_lineage"),
        ("phospy.io", "phospy.core.taxonomy"),
        ("phospy.io", "phospy.core.validation"),
        ("phospy.io", "phospy.core.activities"),
        ("phospy.io", "phospy.signalome"),
        ("phospy.core.policies", "phospy.core.errors"),
        ("phospy.core.preparation_records", "phospy.core.errors"),
        ("phospy.core.preparation_records", "phospy.core.frames"),
        ("phospy.core.preparation_records", "phospy.core.policies"),
        ("phospy.core.preparation_records", "phospy.core.provenance"),
        ("phospy.kinase", "phospy.core.errors"),
        ("phospy.kinase", "phospy.core.activities"),
        ("phospy.kinase", "phospy.datasets"),
        ("phospy.kinase", "phospy.core.evidence"),
        ("phospy.kinase", "phospy.core.frames"),
        ("phospy.kinase", "phospy.core.phosphosites"),
        ("phospy.kinase", "phospy.core.policies"),
        ("phospy.kinase", "phospy.core.provenance"),
        ("phospy.kinase", "phospy.core.quantitative"),
        ("phospy.kinase", "phospy.core.references"),
        ("phospy.kinase", "phospy.core.run_lineage"),
        ("phospy.kinase", "phospy.core.taxonomy"),
        ("phospy.kinase", "phospy.core.statistics"),
        ("phospy.core.provenance", "phospy.core.errors"),
        ("phospy.preprocessing", "phospy.core.errors"),
        ("phospy.preprocessing", "phospy.core.evidence"),
        ("phospy.preprocessing", "phospy.core.frames"),
        ("phospy.preprocessing", "phospy.core.phosphosites"),
        ("phospy.preprocessing", "phospy.core.policies"),
        ("phospy.preprocessing", "phospy.core.provenance"),
        ("phospy.preprocessing", "phospy.core.preparation_records"),
        ("phospy.preprocessing", "phospy.core.quantitative"),
        ("phospy.preprocessing", "phospy.core.references"),
        ("phospy.preprocessing", "phospy.core.taxonomy"),
        ("phospy.core.phosphosites", "phospy.core.errors"),
        ("phospy.core.phosphosites", "phospy.core.frames"),
        ("phospy.core.phosphosites", "phospy.core.taxonomy"),
        ("phospy.core.quantitative", "phospy.core.errors"),
        ("phospy.core.quantitative", "phospy.core.provenance"),
        ("phospy.core.references", "phospy.core.errors"),
        ("phospy.core.references", "phospy.core.frames"),
        ("phospy.core.references", "phospy.core.phosphosites"),
        ("phospy.core.references", "phospy.core.policies"),
        ("phospy.core.references", "phospy.core.provenance"),
        ("phospy.core.references", "phospy.core.taxonomy"),
        ("phospy.core.run_lineage", "phospy.core.errors"),
        ("phospy.core.run_lineage", "phospy.core.evidence"),
        ("phospy.core.run_lineage", "phospy.preprocessing"),
        ("phospy.core.run_lineage", "phospy.core.provenance"),
        ("phospy.core.run_lineage", "phospy.core.references"),
        ("phospy.signalome", "phospy.datasets"),
        ("phospy.signalome", "phospy.core.errors"),
        ("phospy.signalome", "phospy.core.evidence"),
        ("phospy.signalome", "phospy.core.frames"),
        ("phospy.signalome", "phospy.kinase"),
        ("phospy.signalome", "phospy.core.phosphosites"),
        ("phospy.signalome", "phospy.core.policies"),
        ("phospy.signalome", "phospy.core.provenance"),
        ("phospy.signalome", "phospy.core.quantitative"),
        ("phospy.signalome", "phospy.core.references"),
        ("phospy.signalome", "phospy.core.run_lineage"),
        ("phospy.signalome", "phospy.core.validation"),
        ("phospy.core.validation", "phospy.core.errors"),
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
        if record.source_module.startswith("phospy.api.advanced")
        if record.target.startswith("phospy.workflows")
    }

    assert advanced_workflow_targets == set()


FORBIDDEN_PACKAGE_EDGES = frozenset(
    {
        ("phospy.core.errors", "phospy.api"),
        ("phospy.core.errors", "phospy.contracts"),
        ("phospy.core.errors", "phospy.io"),
        ("phospy.core.errors", "phospy.science"),
        ("phospy.core.errors", "phospy.tables"),
        ("phospy.core.errors", "phospy.core.validation"),
        ("phospy.core.errors", "phospy.workflows"),
        ("phospy.io", "phospy.api"),
        ("phospy.core.phosphosites", "phospy.contracts"),
        ("phospy.core.phosphosites", "phospy.io"),
        ("phospy.core.phosphosites", "phospy.science"),
        ("phospy.core.phosphosites", "phospy.workflows"),
        ("phospy.core.references", "phospy.core.evidence"),
        ("phospy.core.references", "phospy.io"),
        ("phospy.core.references", "phospy.science"),
        ("phospy.core.evidence", "phospy.core.references"),
        ("phospy.science", "phospy.contracts"),
        ("phospy.science", "phospy.io"),
        ("phospy.science", "phospy.tables"),
        ("phospy.science", "phospy.core.validation"),
        ("phospy.science", "phospy.workflows"),
        ("phospy.core.validation", "phospy.workflows"),
    }
)


def test_package_dependency_graph_is_acyclic() -> None:
    graph = _build_import_graph()
    package_edges = _architectural_package_edges(graph)
    cycles = [
        tuple(sorted(component))
        for component in strongly_connected_components(
            _package_adjacency(package_edges)
        )
        if len(component) > 1
    ]

    assert cycles == [], _format_cycles(cycles)


def test_complete_eager_module_graph_is_acyclic_without_classification() -> None:
    graph = build_eager_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    cycles = [
        tuple(sorted(component))
        for component in strongly_connected_components(graph.module_edges)
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
            item["rationale"],
        )
        for item in policy["classifications"]
    )

    assert observed == _GOVERNED_PACKAGE_CONDENSATION_CLASSIFICATIONS
    for classification in observed:
        (
            _classification_id,
            source_path,
            source_module,
            target_module,
            symbol,
            rationale,
        ) = classification
        assert rationale.strip()
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
        "from phospy.core.validation.common.input_sources import "
        "validate_dataframe_or_path as validate_source\n",
        "from phospy.core.validation.common.input_sources import *\n",
        "from importlib import import_module\n"
        "validate_source = import_module("
        "'phospy.core.validation.common.input_sources').validate_dataframe_or_path\n",
        "register('phospy.core.validation.common.input_sources.validate_dataframe_or_path')\n",
        "from phospy.core.validation.common.input_sources import "
        "validate_dataframe_or_path\n"
        "from phospy.core.validation.common.input_sources import "
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
        target_module="phospy.core.validation.common.input_sources",
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
    from phospy.core.references.models import ReferenceBundle
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
    assert "phospy.core.references.models" in imports
    assert "phospy.core.references.models.ReferenceBundle" in imports
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


def _assert_imports_load_in_clean_interpreter(code_lines: Iterable[str]) -> None:
    code = "\n".join(code_lines)
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


def test_supported_aggregate_facades_load_in_clean_interpreter() -> None:
    _assert_imports_load_in_clean_interpreter(
        [
            "from phospy.api import AnalysisReadyDatasetBuilder, "
            "DatasetBuildRequest, KinaseWorkflow, PhosPyError, "
            "ReferenceBundle, ReferenceBundleBuilder",
            "from phospy.api.advanced import SpsDiscoveryRequest",
            "assert AnalysisReadyDatasetBuilder",
            "assert KinaseWorkflow",
            "assert DatasetBuildRequest",
            "assert PhosPyError",
            "assert ReferenceBundle",
            "assert ReferenceBundleBuilder",
            "assert SpsDiscoveryRequest",
        ]
    )


def test_internal_owner_modules_load_after_facades_in_clean_interpreter() -> None:
    _assert_imports_load_in_clean_interpreter(
        [
            "from phospy.api import AnalysisReadyDatasetBuilder, KinaseWorkflow, "
            "ReferenceBundleBuilder",
            "from phospy.datasets.contracts import DatasetBuildRequest",
            "from phospy.core.errors import PhosPyError",
            "from phospy.core.references.models import ReferenceBundle",
            "from phospy.datasets.builder import "
            "AnalysisReadyDatasetBuilder as ScienceDatasetBuilder",
            "from phospy.core.references.validation.bundle import ReferenceBundleValidator",
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


@pytest.mark.parametrize("source_owner", ["quantitative", "statistics"])
def test_shared_owner_forbidden_roots_follow_current_ownership_policy(
    source_owner: str,
) -> None:
    forbidden_roots = _shared_owner_forbidden_roots(source_owner)

    assert forbidden_roots
    assert "phospy.datasets" in forbidden_roots
    assert f"phospy.core.{source_owner}" not in forbidden_roots
    assert _package_roots_for_owners(SHARED_ALLOWED_TARGETS[source_owner]).isdisjoint(
        forbidden_roots
    )


@pytest.mark.parametrize("source_owner", ["quantitative", "statistics"])
def test_shared_owner_imports_load_no_forbidden_current_owner(
    source_owner: str,
) -> None:
    forbidden = _run_shared_owner_import_probe(
        source_module=f"phospy.core.{source_owner}",
        forbidden_roots=_shared_owner_forbidden_roots(source_owner),
    )

    assert forbidden == ()


@pytest.mark.parametrize(
    "injected_module",
    ["phospy.datasets", "phospy.datasets.models"],
)
def test_shared_owner_import_probe_detects_forbidden_root_and_descendant(
    injected_module: str,
) -> None:
    forbidden = _run_shared_owner_import_probe(
        source_module="phospy.core.quantitative",
        forbidden_roots=_shared_owner_forbidden_roots("quantitative"),
        injected_module=injected_module,
    )

    assert injected_module in forbidden


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
            "phospy.core.references.builder",
            "phospy.core.references.models",
            "phospy.core.references.tables",
            "phospy.core.references.sequences.repository",
        ),
        "evidence": (
            "phospy.core.evidence.localisation",
            "phospy.core.evidence.sequence_resolution",
            "phospy.core.evidence.quantification_depth",
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
            "from phospy.core.references.builder import ReferenceBundleBuilder",
            "from phospy.core.references.tables import KinaseSubstrateReference",
            "from phospy.core.references.sequences.repository import "
            "FastaProteinSequenceRepository",
            "from phospy.core.evidence.localisation import "
            "normalise_localisation_confidence_series",
            "from phospy.core.evidence.quantification_depth import "
            "validate_quantification_depth_series",
            "from phospy.core.evidence.sequence_resolution import "
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
        "ADR-0053 D05B-1 authorises this exact Dataset builder acquisition of "
        "the IO-owned path reader while IO retains filesystem ownership.",
    ),
    (
        "dataset-input-source-validation",
        "src/phospy/datasets/construction/input_resolution.py",
        "phospy.datasets.construction.input_resolution",
        "phospy.core.validation.common.input_sources",
        "validate_dataframe_or_path",
        "ADR-0053 D05B-2 authorises this exact Dataset input-resolution use of "
        "the shared DataFrame-or-path validation helper.",
    ),
    (
        "dataset-request-source-validation",
        "src/phospy/datasets/validation/inputs.py",
        "phospy.datasets.validation.inputs",
        "phospy.core.validation.common.input_sources",
        "validate_dataframe_or_path",
        "ADR-0053 D05B-2 authorises this exact Dataset request-validation use "
        "of the shared DataFrame-or-path validation helper.",
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
        _rationale,
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
    return architectural_package_unit(
        module_name,
        owners_by_module=_MODULE_OWNERS,
    )


def _package_adjacency(
    package_edges: Iterable[tuple[str, str]],
) -> dict[str, frozenset[str]]:
    adjacency: dict[str, set[str]] = {}
    for source, target in package_edges:
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set())
    return {source: frozenset(targets) for source, targets in adjacency.items()}


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
            if _architectural_package_name(record.source_module) == source
            and _architectural_package_name(record.target) == target
        ][:8]
        for record in examples:
            lines.append(
                f"  {record.source_path.relative_to(PROJECT_ROOT)}:{record.line}"
            )
    return "\n".join(lines)
