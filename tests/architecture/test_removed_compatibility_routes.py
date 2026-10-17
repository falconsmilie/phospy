from __future__ import annotations

import importlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

import phospy.kinase.scoring.rank_fusion as scoring_rank_fusion
import phospy.preprocessing.planning.registry as stage_registry
import phospy.signalome.clustering.exact_python as exact_python

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RETIRED_SIGNALOME_MODULES = (
    "phospy.contracts.configs.signalome",
    "phospy.contracts.results.signalome",
    "phospy.science.signalomes",
    "phospy.science.signalomes._result_validation",
    "phospy.science.signalomes.assignments",
    "phospy.science.signalomes.clustering",
    "phospy.science.signalomes.clustering.backend_dispatch",
    "phospy.science.signalomes.clustering.backends",
    "phospy.science.signalomes.clustering.backends.exact_python",
    "phospy.science.signalomes.clustering.backends.scipy_hierarchical",
    "phospy.science.signalomes.clustering.candidate_scoring",
    "phospy.science.signalomes.clustering.candidate_selection",
    "phospy.science.signalomes.clustering.contracts",
    "phospy.science.signalomes.clustering.diagnostic_schemas",
    "phospy.science.signalomes.clustering.diagnostics",
    "phospy.science.signalomes.clustering.exact_python",
    "phospy.science.signalomes.clustering.models",
    "phospy.science.signalomes.clustering.module_selection",
    "phospy.science.signalomes.clustering.orchestration",
    "phospy.science.signalomes.clustering.policies",
    "phospy.science.signalomes.clustering.protein_modules",
    "phospy.science.signalomes.clustering.scale_guards",
    "phospy.science.signalomes.clustering.scientific_policies",
    "phospy.science.signalomes.clustering.scipy_hierarchical",
    "phospy.science.signalomes.clustering.scoring",
    "phospy.science.signalomes.clustering.selection",
    "phospy.science.signalomes.clustering.stability",
    "phospy.science.signalomes.clustering.tree_building",
    "phospy.science.signalomes.clustering.tree_engine_adapter",
    "phospy.science.signalomes.clustering.validation",
    "phospy.science.signalomes.constants",
    "phospy.science.signalomes.context",
    "phospy.science.signalomes.expanded",
    "phospy.science.signalomes.models",
    "phospy.science.signalomes.modules",
    "phospy.science.signalomes.network",
    "phospy.science.signalomes.science",
    "phospy.science.tables.signalome",
    "phospy.science.tables.signalome.assignments",
    "phospy.science.tables.signalome.common",
    "phospy.science.tables.signalome.context",
    "phospy.science.tables.signalome.modules",
    "phospy.science.tables.signalome.network",
    "phospy.workflows.signalome",
    "phospy.workflows.signalome.alignment_diagnostics",
    "phospy.workflows.signalome.boundary_errors",
    "phospy.workflows.signalome.caveats",
    "phospy.workflows.signalome.clustering_runner",
    "phospy.workflows.signalome.component_helpers",
    "phospy.workflows.signalome.component_models",
    "phospy.workflows.signalome.constants",
    "phospy.workflows.signalome.context_tables",
    "phospy.workflows.signalome.contracts",
    "phospy.workflows.signalome.executor",
    "phospy.workflows.signalome.interpreter",
    "phospy.workflows.signalome.matrix_alignment",
    "phospy.workflows.signalome.module_tables",
    "phospy.workflows.signalome.network_builder",
    "phospy.workflows.signalome.protein_resolution",
    "phospy.workflows.signalome.provenance",
    "phospy.workflows.signalome.public",
    "phospy.workflows.signalome.result_assembly",
    "phospy.workflows.signalome.row_attrition",
    "phospy.workflows.signalome.scientific_policies",
    "phospy.workflows.signalome.score_matrix_selection",
    "phospy.workflows.signalome.score_preconditioning",
    "phospy.workflows.signalome.validator",
)

RETIRED_HORIZONTAL_MODULES = (
    "phospy.contracts",
    "phospy.contracts.configs",
    "phospy.contracts.requests",
    "phospy.contracts.results",
    "phospy.contracts.results.base",
    "phospy.science",
    "phospy.science.configs",
    "phospy.science.datasets",
    "phospy.science.references",
    "phospy.science.references.validation",
    "phospy.science.tables",
    "phospy.workflows",
    "phospy.validation.configs",
    "phospy.validation.configs.preprocessing",
    "phospy.validation.datasets",
    "phospy.validation.datasets.analysis_ready",
    "phospy.validation.datasets.fragpipe",
    "phospy.validation.datasets.group_coverage_filter",
    "phospy.validation.datasets.importers",
    "phospy.validation.datasets.maxquant",
    "phospy.validation.references",
    "phospy.validation.workflows",
)

COMPATIBILITY_ONLY_MODULES = (
    (
        "phospy._deprecations",
        "phospy.api._compat",
        "phospy.api.datasets",
        "phospy.api.enums",
        "phospy.api.requests",
        "phospy.api.results",
        "phospy.api.workflows",
        "phospy.advanced.configs",
        "phospy.advanced.results",
        "phospy.contracts.configs._validation",
        "phospy.contracts.configs.common",
        "phospy.contracts.configs.enrichment",
        "phospy.contracts.configs.preprocessing.control_sites",
        "phospy.contracts.enrichment_identifier_sets",
        "phospy.contracts.result_caveats",
        "phospy.contracts.results.caveats",
        "phospy.contracts.results.enrichment",
        "phospy.provenance.organisms",
        "phospy.science.configs._validation",
        "phospy.science.configs.common",
        "phospy.science.configs.enrichment",
        "phospy.science.configs.preprocessing.control_sites",
        "phospy.science.batch_correction.sps_discovery",
        "phospy.science.datasets.preprocessing.control_sites",
        "phospy.science.enrichment.collections",
        "phospy.science.enrichment.models",
        "phospy.science.enrichment.ora",
        "phospy.science.result_caveats",
        "phospy.science.sites",
        "phospy.science.sites.identifiers",
        "phospy.science.sites.identity",
        "phospy.science.sites.identity_columns",
        "phospy.science.sites.identity_contracts",
        "phospy.science.sites.identity_rules",
        "phospy.science.sites.identity_rules.contracts",
        "phospy.science.sites.identity_rules.dataset_identity",
        "phospy.science.sites.identity_rules.messages",
        "phospy.science.sites.identity_rules.parsing",
        "phospy.science.sites.identity_rules.reference_context",
        "phospy.science.sites.identity_rules.result_identity",
        "phospy.science.sites.metadata_validation",
        "phospy.science.sites.organisms",
        "phospy.science.sites.sequence_context",
        "phospy.science.sites.site_keys",
        "phospy.science.sites.validation",
        "phospy.signalome.clustering.protocol",
        "phospy.science.statistics.multiple_testing",
        "phospy.differential.multiple_testing",
        "phospy.science.quantitative_method_contracts",
        "phospy.science.transformations._authority",
        "phospy.science.transformations._payloads",
        "phospy.science.transformations.contracts",
        "phospy.science.transformations.models",
        "phospy.science.transformations.policy",
        "phospy.science.transformations.provenance",
        "phospy.science.transformations.quantitative_contracts",
        "phospy.science.transformations.scale_state",
        "phospy.science.transformations.scale_values",
        "phospy.science.transformations.state_coherence",
        "phospy.science.transformations.transformers.identity",
        "phospy.science.transformations.transformers.log2",
        "phospy.science.tables.base",
        "phospy.tables.activity",
        "phospy.tables.base",
        "phospy.tables.datasets",
        "phospy.tables.differential",
        "phospy.tables.kinase",
        "phospy.tables.references",
        "phospy.validation.workflows.configs",
        "phospy.validation.identity_contracts",
        "phospy.validation.datasets.display_site_identity",
        "phospy.validation.datasets.protein_scoped_site_identity",
        "phospy.validation.transformations.state",
        "phospy.validation.workflows.identity",
        "phospy.validation.workflows.enrichment",
        "phospy.validation.workflows.batch_correction.control_site_workflow",
        "phospy.validation.workflows.batch_correction.control_sites",
        "phospy.validation.workflows.batch_correction.sps_discovery",
        "phospy.validation.workflows.method_quantitative",
        "phospy.validation.workflows.quantitative",
        "phospy.workflows._pandas_typing",
        "phospy.workflows._row_attrition",
        "phospy.workflows.intensity_scale_evidence",
        "phospy.workflows.enrichment.caveats",
        "phospy.workflows.enrichment.executor",
        "phospy.workflows.enrichment.interpreter",
        "phospy.workflows.enrichment.models",
        "phospy.workflows.enrichment.ora_execution",
        "phospy.workflows.enrichment.public",
        "phospy.workflows.enrichment.provenance",
        "phospy.workflows.enrichment.result_assembly",
        "phospy.workflows.enrichment.set_filtering",
        "phospy.workflows.enrichment.validator",
        "phospy.workflows.enrichment.workflow",
        "phospy.workflows.batch_correction.sps_discovery",
        "phospy.contracts.configs.preprocessing.batch_correction",
        "phospy.contracts.configs.preprocessing.correction_missingness",
        "phospy.contracts.configs.preprocessing.internal_batch_correction",
        "phospy.science.batch_correction.executor",
        "phospy.science.batch_correction.ruv_iii",
        "phospy.science.batch_correction.ruv_iii_executor",
        "phospy.science.configs.preprocessing.batch_correction",
        "phospy.science.configs.preprocessing.correction_missingness",
        "phospy.science.configs.preprocessing.internal_batch_correction",
        "phospy.science.datasets.preprocessing.batch_correction",
        "phospy.science.datasets.preprocessing.batch_correction_engine",
        "phospy.science.datasets.preprocessing.batch_correction_metadata",
        "phospy.science.datasets.preprocessing.batch_correction_models",
        "phospy.science.datasets.preprocessing.batch_correction_provenance",
        "phospy.science.datasets.preprocessing.batch_correction_provenance_payloads",
        "phospy.science.datasets.preprocessing.batch_correction_provenance_validation",
        "phospy.validation.datasets.batch_correction",
        "phospy.validation.datasets.batch_correction_controls",
        "phospy.validation.datasets.batch_correction_design",
        "phospy.validation.datasets.batch_correction_provenance",
        "phospy.validation.workflows.batch_correction.design",
        "phospy.validation.workflows.batch_correction.missingness",
        "phospy.validation.workflows.batch_correction.protocols",
        "phospy.validation.workflows.batch_correction.request",
        "phospy.validation.workflows.batch_correction.stage_order",
        "phospy.workflows.batch_correction.contracts",
        "phospy.workflows.batch_correction.interpreter",
        "phospy.workflows.batch_correction.preprocessing_adapter",
        "phospy.workflows.batch_correction.provenance",
        "phospy.workflows.batch_correction.workflow",
        "phospy.io.bundles._shared.intensity_scale_state",
        "phospy.contracts.kinase_reference_projection",
        "phospy.contracts.configs.prediction",
        "phospy.science.configs.prediction",
        "phospy.science.prediction",
        "phospy.science.prediction.candidates",
        "phospy.science.prediction.execution",
        "phospy.science.prediction.internal_view",
        "phospy.science.prediction.models",
        "phospy.science.prediction.policies",
        "phospy.science.prediction.sampling_core",
        "phospy.science.prediction.sampling_runtime",
        "phospy.science.prediction.scientific_policies",
        "phospy.science.prediction.sequence_validation",
        "phospy.science.prediction.svm",
        "phospy.science.prediction.motif_scoring",
        "phospy.science.prediction.motif_scoring.frequency_matrices",
        "phospy.science.prediction.motif_scoring.kinase_library",
        "phospy.science.prediction.motif_scoring.library_validation",
        "phospy.science.prediction.motif_scoring.models",
        "phospy.science.prediction.motif_scoring.scaling",
        "phospy.science.prediction.motif_scoring.scoring",
        "phospy.science.prediction.scoring",
        "phospy.science.scoring",
        "phospy.science.scoring.policy_models",
        "phospy.science.references.kinase_library",
        "phospy.science.references.kinase_library_loading",
        "phospy.science.references.kinase_library_models",
        "phospy.science.references.validation.kinase_library",
        "phospy.science.tables.kinase",
        "phospy.validation.references.kinase_library",
        "phospy.workflows.kinase.contributions",
        "phospy.workflows.kinase.kinase_library_scoring",
        "phospy.workflows.kinase.reference_projection",
        "phospy.workflows.kinase.scoring_mode_contracts",
        "phospy.workflows.kinase.scoring_runner",
        "phospy.workflows.kinase.scoring_transforms",
        "phospy.workflows.kinase.component_models",
        "phospy.workflows.kinase.prediction_runner",
        "phospy.workflows.kinase.science",
        "phospy.science.activities",
        "phospy.science.activities.diagnostics",
        "phospy.science.activities.inputs",
        "phospy.science.activities.membership",
        "phospy.science.activities.method_contracts",
        "phospy.science.activities.method_models",
        "phospy.science.activities.method_quantitative_validation",
        "phospy.science.activities.methods",
        "phospy.science.activities.methods.ksea_zscore",
        "phospy.science.activities.methods.ssgsea_substrate_enrichment",
        "phospy.science.activities.methods.weighted_substrate_activity",
        "phospy.science.activities.models",
        "phospy.science.activities.quantitative_method_contracts",
        "phospy.science.activities.result_validation",
        "phospy.science.activities.results",
        "phospy.science.activities.scientific_policies",
        "phospy.science.activities.scoring",
        "phospy.science.activities.semantics",
        "phospy.science.activities.statistics",
        "phospy.science.activities.threshold_membership",
        "phospy.science.tables.activity",
        "phospy.validation.workflows.activity",
        "phospy.workflows.kinase.activity_runner",
        "phospy.workflows.kinase.membership",
        "phospy.workflows.kinase",
        "phospy.contracts.results.kinase",
        "phospy.contracts.results.kinase_boundary",
        "phospy.contracts.configs.kinase",
        "phospy.science.configs.kinase",
    )
    + RETIRED_SIGNALOME_MODULES
    + RETIRED_HORIZONTAL_MODULES
)


@pytest.mark.parametrize("module_name", COMPATIBILITY_ONLY_MODULES)
def test_compatibility_only_modules_are_physically_absent(module_name: str) -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module(module_name)


def test_deleted_compatibility_package_directories_have_no_source_modules() -> None:
    for relative in (
        "api/configs",
        "science/statistics",
        "science/signalomes",
        "science/tables/signalome",
        "science/activities",
        "science/enrichment",
        "science/prediction",
        "science/scoring",
        "science/transformations",
        "tables",
        "tables/signalome",
        "validation/transformations",
        "workflows/enrichment",
        "workflows/kinase",
        "workflows/signalome",
    ):
        assert list((PROJECT_ROOT / "src" / "phospy" / relative).rglob("*.py")) == []


def test_every_retired_signalome_route_fails_in_a_fresh_process() -> None:
    script = """
import importlib
import json
import sys

survivors = []
for module_name in json.loads(sys.argv[1]):
    try:
        importlib.import_module(module_name)
    except ModuleNotFoundError:
        continue
    survivors.append(module_name)
if survivors:
    raise SystemExit(f"retired Signalome routes still importable: {survivors}")
"""
    completed = subprocess.run(
        [sys.executable, "-c", script, json.dumps(RETIRED_SIGNALOME_MODULES)],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


@pytest.mark.parametrize(
    "module_order",
    [
        (
            "phospy.signalome.requests",
            "phospy.signalome.config",
            "phospy.signalome.results",
            "phospy.signalome.errors",
            "phospy.api",
            "phospy.advanced",
            "phospy",
        ),
        (
            "phospy",
            "phospy.advanced",
            "phospy.api",
            "phospy.signalome.errors",
            "phospy.signalome.results",
            "phospy.signalome.config",
            "phospy.signalome.requests",
        ),
    ],
)
def test_signalome_facade_identities_are_stable_in_opposing_fresh_import_orders(
    module_order: tuple[str, ...],
) -> None:
    script = """
import importlib
import json
import sys

for module_name in json.loads(sys.argv[1]):
    importlib.import_module(module_name)

import phospy
import phospy.advanced as advanced
import phospy.api as api
from phospy.signalome.config import SignalomeConfig
from phospy.signalome.errors import SignalomeScaleError
from phospy.signalome.public import SignalomeWorkflow
from phospy.signalome.requests import SignalomeWorkflowRequest
from phospy.signalome.results import SignalomeWorkflowResult

assert phospy.SignalomeWorkflow is api.SignalomeWorkflow is SignalomeWorkflow
assert api.SignalomeWorkflowRequest is SignalomeWorkflowRequest
assert advanced.SignalomeConfig is SignalomeConfig
assert api.SignalomeWorkflowResult is SignalomeWorkflowResult
assert api.SignalomeScaleError is SignalomeScaleError
"""
    completed = subprocess.run(
        [sys.executable, "-c", script, json.dumps(module_order)],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_former_semi_public_exports_are_internal_only() -> None:
    assert not hasattr(stage_registry, "PreprocessingStageMetadata")
    assert not hasattr(
        scoring_rank_fusion, "fuse_profile_and_motif_scores_by_rank_weight"
    )
    assert "__all__" not in vars(exact_python)
    for name in ("_WardClusterTree", "_build_cluster_tree", "cluster_sites"):
        assert not hasattr(exact_python, name)


def test_removed_registry_and_dispatch_tokens_do_not_survive_in_runtime_source() -> (
    None
):
    package_root = PROJECT_ROOT / "src" / "phospy"
    combined = "\n".join(
        path.read_text(encoding="utf-8") for path in sorted(package_root.rglob("*.py"))
    )
    for token in (
        "PhosPyDeprecationWarning",
        "warn_deprecated",
        "deprecated_advanced_export",
        "deprecated_config_export",
        "deprecated_request_export",
    ):
        assert token not in combined
