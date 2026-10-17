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
ENRICHMENT_ROOT = PACKAGE_ROOT / "enrichment"


def test_enrichment_owner_has_no_forbidden_outward_dependency() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    forbidden_prefixes = (
        "phospy.io",
        "phospy.science.datasets",
        "phospy.activities",
        "phospy.kinase.activity",
        "phospy.differential",
        "phospy.science.prediction",
        "phospy.signalome",
        "phospy.differential",
        "phospy.kinase",
        "phospy.signalome",
    )
    offenders = [
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.enrichment")
        if record.target.startswith(forbidden_prefixes)
    ]
    reverse_provenance = [
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.provenance")
        if record.target.startswith("phospy.enrichment")
    ]

    assert offenders == []
    assert reverse_provenance == []


def test_enrichment_request_does_not_import_dataset_input_alias() -> None:
    request_path = ENRICHMENT_ROOT / "contracts" / "requests.py"
    source = request_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(request_path))
    imported_targets = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }

    assert "DatasetInput" not in source
    assert all("datasets" not in target for target in imported_targets)
    assert "from pathlib import Path" in source
    assert "from os import PathLike" in source
    assert "import pandas as pd" in source


def test_io_enrichment_reader_imports_only_the_two_passive_owner_leaves() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    targets = {
        record.target
        for record in graph.records
        if record.source_module == "phospy.io.readers.enrichment_sets"
        if record.target.startswith("phospy.enrichment")
    }

    assert targets == {
        "phospy.enrichment.contracts.collections",
        "phospy.enrichment.contracts.vocabulary",
    }


@pytest.mark.parametrize(
    "module_order",
    (
        (
            "phospy.io.readers.enrichment_sets",
            "phospy.enrichment.contracts.collections",
            "phospy.enrichment.contracts.requests",
        ),
        (
            "phospy.enrichment.contracts.requests",
            "phospy.enrichment.contracts.collections",
            "phospy.io.readers.enrichment_sets",
        ),
        (
            "phospy.api",
            "phospy.advanced",
            "phospy.enrichment.workflow",
            "phospy.enrichment.contracts.results",
        ),
        (
            "phospy.enrichment.contracts.results",
            "phospy.enrichment.workflow",
            "phospy.advanced",
            "phospy.api",
        ),
        (
            "phospy.provenance",
            "phospy.enrichment.provenance",
            "phospy.enrichment.contracts.results",
        ),
        (
            "phospy.enrichment.contracts.results",
            "phospy.enrichment.provenance",
            "phospy.provenance",
        ),
    ),
)
def test_enrichment_import_orders_construct_and_execute_real_objects(
    module_order: tuple[str, ...],
    tmp_path: Path,
) -> None:
    code = f"""
import importlib
for module_name in {module_order!r}:
    importlib.import_module(module_name)

import phospy.api as api
from phospy.enrichment.contracts.collections import GeneSetCollection, PtmSetCollection
from phospy.enrichment.contracts.config import EnrichmentConfig
from phospy.enrichment.contracts.requests import EnrichmentWorkflowRequest
from phospy.enrichment.contracts.results import EnrichmentWorkflowResult
from phospy.enrichment.contracts.vocabulary import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
)
from phospy.enrichment.workflow import EnrichmentWorkflow

assert api.EnrichmentConfig is EnrichmentConfig
assert api.EnrichmentWorkflowRequest is EnrichmentWorkflowRequest
assert api.EnrichmentWorkflowResult is EnrichmentWorkflowResult
assert api.EnrichmentWorkflow is EnrichmentWorkflow
assert tuple(EnrichmentWorkflow.__init__.__annotations__) == ('return',)
workflow = EnrichmentWorkflow()
gene = workflow.run(EnrichmentWorkflowRequest(
    identifier_column='gene_symbol',
    identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    set_collection=GeneSetCollection(
        sets={{'B': ('B',), 'A': ('A', 'B')}},
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    source_name='enrichment_import_order',
    ),
    selected_identifiers=('A',),
    background_universe=('A', 'B'),
))
ptm = workflow.run(EnrichmentWorkflowRequest(
    identifier_column='site_key',
    identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    set_collection=PtmSetCollection(
        sets={{'M': ('rat|P1|S1',), 'N': ('rat|P1|T2',)}},
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    source_name='enrichment_import_order',
    ),
    selected_identifiers=('rat|P1|S1',),
    background_universe=('rat|P1|S1', 'rat|P1|T2'),
))
for result in (gene, ptm):
    assert result.provenance is not None
    assert result.provenance.workflow_parameters['correction_owner'] == 'ora_engine'
    assert result.provenance.workflow_parameters['number_of_tests'] == 2
assert [record.term_id for record in gene.records] == ['A', 'B']
assert [record.term_id for record in ptm.records] == ['M', 'N']
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
