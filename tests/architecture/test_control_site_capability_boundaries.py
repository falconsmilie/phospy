from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest
from import_graph_support import build_import_graph
from logical_ownership_policy import PUBLIC_FACADE_ALLOWED_TARGETS

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
CONTROL_SITE_ROOT = PACKAGE_ROOT / "preprocessing" / "control_sites"


def test_control_site_discovery_module_graph_is_acyclic() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    modules = {
        module
        for module in graph.modules
        if module.startswith("phospy.preprocessing.control_sites")
    }
    dependencies = {
        module: set(graph.module_edges[module]).intersection(modules)
        for module in modules
    }

    while dependencies:
        leaves = {module for module, targets in dependencies.items() if not targets}
        assert leaves, f"cyclic control-site modules: {sorted(dependencies)!r}"
        dependencies = {
            module: targets.difference(leaves)
            for module, targets in dependencies.items()
            if module not in leaves
        }


def test_passive_discovery_boundaries_do_not_import_numerical_execution() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    passive_modules = {
        "phospy.preprocessing.control_sites.discovery_contracts",
        "phospy.preprocessing.control_sites.discovery_validation",
    }
    forbidden_targets = {
        "phospy.preprocessing.control_sites.discovery_execution",
        "phospy.preprocessing.control_sites.discovery_ranking",
        "phospy.preprocessing.control_sites.discovery_results",
    }

    assert {
        record.target
        for record in graph.records
        if record.source_module in passive_modules
        if record.target in forbidden_targets
    } == set()


def test_public_facade_rejects_internal_control_site_planning_and_validation() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    internal_modules = {
        "phospy.preprocessing.control_sites.discovery_planning",
        "phospy.preprocessing.control_sites.discovery_validation",
        "phospy.preprocessing.control_sites.validation",
    }

    assert internal_modules.isdisjoint(PUBLIC_FACADE_ALLOWED_TARGETS["preprocessing"])
    assert {
        record.target
        for record in graph.records
        if record.source_module == "phospy.advanced"
        if record.target in internal_modules
    } == set()


def test_control_site_owner_has_no_dataset_io_or_analytical_dependency() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    forbidden_prefixes = (
        "phospy.io",
        "phospy.science.datasets",
        "phospy.differential",
        "phospy.activities",
        "phospy.kinase.activity",
        "phospy.science.prediction",
        "phospy.signalome",
    )
    offenders = [
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.preprocessing.control_sites")
        if record.target.startswith(forbidden_prefixes)
    ]
    reverse_provenance = [
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.provenance")
        if record.target.startswith("phospy.preprocessing.control_sites")
    ]

    assert offenders == []
    assert reverse_provenance == []


def test_quantitative_diagnostic_leaf_has_no_capability_dependency() -> None:
    path = PACKAGE_ROOT / "quantitative" / "declared_scale_diagnostics.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    targets = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }

    assert all("preprocessing" not in target for target in targets)
    assert all("datasets" not in target for target in targets)


def test_control_site_owner_does_not_restore_dataset_resolver_edge() -> None:
    sources = "\n".join(
        path.read_text(encoding="utf-8") for path in CONTROL_SITE_ROOT.glob("*.py")
    )

    assert "DatasetIntensityScaleResolver" not in sources
    assert "phospy.science.datasets.builders.transformation_resolver" not in (
        "\n".join(
            line
            for line in sources.splitlines()
            if not line.strip().startswith(("#", '"'))
        )
    )


def test_correction_consumes_only_published_control_site_contract_boundaries() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    correction_prefixes = (
        "phospy.preprocessing.batch_correction.execution",
        "phospy.preprocessing.batch_correction",
        "phospy.validation.workflows.batch_correction",
    )
    targets = {
        record.target
        for record in graph.records
        if record.source_module.startswith(correction_prefixes)
        if record.target.startswith("phospy.preprocessing.control_sites")
    }

    assert targets <= {
        "phospy.preprocessing.control_sites.contracts",
        "phospy.preprocessing.control_sites.validation",
    }


def test_advanced_exports_are_direct_control_site_owner_identities() -> None:
    import phospy.advanced as advanced
    from phospy.preprocessing.control_sites.contracts import ControlSiteSet
    from phospy.preprocessing.control_sites.discovery_contracts import (
        SpsDiscoveryRequest,
    )
    from phospy.preprocessing.control_sites.discovery_results import (
        SpsDiscoveryResult,
    )
    from phospy.preprocessing.control_sites.workflow import SpsDiscoveryWorkflow

    assert advanced.ControlSiteSet is ControlSiteSet
    assert advanced.SpsDiscoveryRequest is SpsDiscoveryRequest
    assert advanced.SpsDiscoveryResult is SpsDiscoveryResult
    assert advanced.SpsDiscoveryWorkflow is SpsDiscoveryWorkflow


@pytest.mark.parametrize(
    "module_order",
    (
        (
            "phospy.preprocessing.control_sites.workflow",
            "phospy.preprocessing.control_sites.discovery_results",
            "phospy.preprocessing.batch_correction.contracts.workflow",
            "phospy.preprocessing.batch_correction.planning",
        ),
        (
            "phospy.preprocessing.batch_correction.planning",
            "phospy.preprocessing.batch_correction.contracts.workflow",
            "phospy.preprocessing.control_sites.discovery_results",
            "phospy.preprocessing.control_sites.workflow",
        ),
    ),
)
def test_discovery_restore_and_correction_work_in_both_import_orders(
    module_order: tuple[str, ...],
    tmp_path: Path,
) -> None:
    code = f"""
import importlib
for module_name in {module_order!r}:
    importlib.import_module(module_name)

import pandas as pd
from phospy.advanced import (
    ControlSiteSet,
    CorrectionMissingnessPolicy,
    ObservationMask,
    OriginallyMissingCellTracking,
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryResult,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
)
from phospy.preprocessing.batch_correction.contracts import (
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
)
from phospy.preprocessing.batch_correction import (
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
)

sites = pd.Index((
    'phospy:v1|organism=rat|protein_namespace=protein_id|protein_identifier=P1|residue=S|position=1',
    'phospy:v1|organism=rat|protein_namespace=protein_id|protein_identifier=P2|residue=S|position=2',
    'phospy:v1|organism=rat|protein_namespace=protein_id|protein_identifier=P3|residue=S|position=3',
), name='site_key')
conditions = {{'a1': 'a', 'a2': 'a', 'b1': 'b', 'b2': 'b'}}
references = tuple(
    SpsReferenceDataset.from_condition_relative_log2(
        dataset_id=dataset_id,
        intensities=pd.DataFrame(values, index=sites),
        condition_by_sample=conditions,
        log2_scale_established_by=f'{{dataset_id}} scale',
        baseline_centering_established_by=f'{{dataset_id}} baseline',
        organism='rat',
        baseline_context='a',
        reference_context='fresh interpreter fixture',
        source_name=dataset_id,
        source_version='1',
        source_uri=f'https://example.test/{{dataset_id}}',
    )
    for dataset_id, values in (
        ('r1', {{'a1': (0.0, 0.0, 0.0), 'a2': (0.0, 0.0, 0.0),
                'b1': (0.1, 0.2, 2.0), 'b2': (0.1, 0.2, 2.0)}}),
        ('r2', {{'a1': (0.0, 0.0, 0.0), 'a2': (0.0, 0.0, 0.0),
                'b1': (0.2, 0.3, 3.0), 'b2': (0.2, 0.3, 3.0)}}),
    )
)
discovery = SpsDiscoveryWorkflow().run(SpsDiscoveryRequest(
    reference_datasets=references,
    config=SpsDiscoveryConfig(top_n=2, minimum_shared_sites=2),
))
restored_result = SpsDiscoveryResult.from_payload(discovery.to_payload())
controls = ControlSiteSet.from_payload(restored_result.control_site_set.to_payload())
assert controls.to_payload() == discovery.control_site_set.to_payload()

samples = ('s1', 's2', 's3', 's4')
target = pd.DataFrame(
    {{'s1': (1.0, 2.0, 5.0), 's2': (2.0, 4.0, 5.5),
     's3': (3.0, 1.0, 6.0), 's4': (4.0, 3.0, 7.0)}},
    index=sites,
)
config = InternalBatchCorrectionRequest(
    method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
    batch_column='batch',
    condition_columns=('condition',),
    replicate_column=None,
    control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
    control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
    missing_value_policy=InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING,
    imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
    n_unwanted_factors=1,
    stage_order=InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM,
    diagnostics_enabled=True,
)
mask = ObservationMask(feature_ids=tuple(sites), sample_ids=samples,
                       originally_missing_cells=())
result = BatchCorrectionWorkflow().run(BatchCorrectionWorkflowRequest(
    phospho=target,
    config=config,
    sample_metadata=pd.DataFrame(
        {{'batch': ('b1', 'b1', 'b2', 'b2'),
         'condition': ('a', 'b', 'a', 'b')}},
        index=samples,
    ),
    control_site_set=controls,
    missingness_policy=CorrectionMissingnessPolicy(
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,
        observation_mask=mask,
    ),
))
assert result.corrected_matrix.shape == target.shape
assert result.provenance.selected_site_key_rows == discovery.selected_site_keys
assert result.provenance.control_site_source['sps_discovery_identity'] == discovery.discovery_identity
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
