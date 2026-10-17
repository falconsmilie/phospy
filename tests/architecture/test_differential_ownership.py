from __future__ import annotations

import ast
import importlib
from pathlib import Path

from import_graph_support import build_import_graph

import phospy
import phospy.advanced
import phospy.api

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
DIFFERENTIAL_ROOT = PACKAGE_ROOT / "differential"

OLD_DIFFERENTIAL_MODULES = (
    "phospy.contracts.configs.differential",
    "phospy.contracts.results.differential",
    "phospy.science.configs.differential",
    "phospy.science.design",
    "phospy.science.design.contrast_helpers",
    "phospy.science.design.matrix_builder",
    "phospy.science.design.models",
    "phospy.science.differential",
    "phospy.science.differential.aggregation",
    "phospy.science.differential.aggregation.executor",
    "phospy.science.differential.aggregation.experimental",
    "phospy.science.differential.aggregation.models",
    "phospy.science.differential.aggregation.scientific_policies",
    "phospy.science.differential.compound_symmetry_gls",
    "phospy.science.differential.duplicate_correlation",
    "phospy.science.differential.empirical_bayes",
    "phospy.science.differential.executor",
    "phospy.science.differential.internal_view",
    "phospy.science.differential.linear_model",
    "phospy.science.differential.models",
    "phospy.science.differential.models.design",
    "phospy.science.differential.models.diagnostics",
    "phospy.science.differential.models.duplicate_correlation",
    "phospy.science.differential.models.empirical_bayes_config",
    "phospy.science.differential.models.fit",
    "phospy.science.differential.models.protein_aware",
    "phospy.science.differential.models.provenance",
    "phospy.science.differential.models.results",
    "phospy.science.differential.models.tables",
    "phospy.science.differential.policy_models",
    "phospy.science.differential.protein_covariate_adjusted",
    "phospy.science.differential.quantification_depth",
    "phospy.science.differential.validator",
    "phospy.science.tables._differential_validation",
    "phospy.science.tables.differential",
    "phospy.validation.workflows.differential",
    "phospy.validation.workflows.differential_design_rules",
    "phospy.validation.workflows.differential_result_tables",
    "phospy.workflows.differential",
    "phospy.workflows.differential.caveats",
    "phospy.workflows.differential.derived_quantitative_lineage",
    "phospy.workflows.differential.design_assembly",
    "phospy.workflows.differential.eligibility",
    "phospy.workflows.differential.executor",
    "phospy.workflows.differential.fitting",
    "phospy.workflows.differential.imputation_inference",
    "phospy.workflows.differential.interpreter",
    "phospy.workflows.differential.models",
    "phospy.workflows.differential.protein_aware_eligibility",
    "phospy.workflows.differential.protein_aware_inputs",
    "phospy.workflows.differential.provenance",
    "phospy.workflows.differential.public",
    "phospy.workflows.differential.reliability",
    "phospy.workflows.differential.replicates",
    "phospy.workflows.differential.result_assembly",
    "phospy.workflows.differential.validator",
)


def test_old_differential_routes_are_absent() -> None:
    loaded = []
    for module_name in OLD_DIFFERENTIAL_MODULES:
        try:
            importlib.import_module(module_name)
        except ModuleNotFoundError:
            continue
        loaded.append(module_name)

    assert loaded == []


def test_public_inventories_and_owner_identities_are_exact() -> None:
    assert tuple(
        map(len, (phospy.__all__, phospy.api.__all__, phospy.advanced.__all__))
    ) == (
        3,
        66,
        124,
    )
    owner_routes = {
        "BatchCovariate": "phospy.differential.design.models",
        "CategoricalCovariate": "phospy.differential.design.models",
        "ContinuousCovariate": "phospy.differential.design.models",
        "FixedEffectCovariate": "phospy.differential.design.models",
        "ExperimentalDesign": "phospy.differential.design.models",
        "SampleDesignRecord": "phospy.differential.design.models",
        "Contrast": "phospy.differential.design.models",
        "DifferentialAnalysisRequest": "phospy.differential.contracts",
        "all_pairwise_contrasts": "phospy.differential.design.contrast_helpers",
        "contrasts_vs_control": "phospy.differential.design.contrast_helpers",
        "DifferentialAnalysisWorkflow": "phospy.differential.public",
        "DifferentialAnalysisResult": "phospy.differential.results",
    }
    for name, module_name in owner_routes.items():
        assert getattr(phospy.api, name) is getattr(
            importlib.import_module(module_name), name
        )
    advanced_owner_routes = {
        "DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1": "phospy.differential.policies",
        "DifferentialAnalysisConfig": "phospy.differential.config",
        "DifferentialImputedValuePolicy": "phospy.differential.policies",
        "DifferentialModelDiagnostics": "phospy.differential.diagnostics",
        "DifferentialProteinAwareModelConfig": "phospy.differential.config",
        "DifferentialReliabilityProfile": "phospy.differential.policies",
        "EmpiricalBayesConfig": "phospy.differential.moderation.config",
        "MultipleTestingConfig": "phospy.differential.config",
        "MultipleTestingMethod": "phospy.differential.policies",
        "PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION": "phospy.differential.policies",
        "PAIRED_DESIGN_POLICY_FIXED_BLOCK": "phospy.differential.policies",
        "PAIRED_DESIGN_POLICY_REJECT": "phospy.differential.policies",
        "PairedDesignPolicy": "phospy.differential.policies",
        "QuantificationDepthTrendDiagnostics": "phospy.differential.diagnostics",
        "TechnicalReplicatePolicy": "phospy.differential.policies",
        "filter_differential_results": "phospy.differential.tables",
        "rank_differential_results": "phospy.differential.tables",
    }
    for name, module_name in advanced_owner_routes.items():
        assert getattr(phospy.advanced, name) is getattr(
            importlib.import_module(module_name), name
        )


def test_owner_initializers_are_narrow() -> None:
    offenders = []
    for path in sorted(DIFFERENTIAL_ROOT.rglob("__init__.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if any(
            isinstance(node, ast.Import | ast.ImportFrom) for node in ast.walk(tree)
        ):
            offenders.append(path.relative_to(PROJECT_ROOT).as_posix())

    assert offenders == []


def test_dependency_boundaries_are_exact() -> None:
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    records = tuple(graph.records)
    upstream_prefixes = ("phospy.datasets", "phospy.preprocessing", "phospy.statistics")
    reverse = sorted(
        (record.source_module, record.target)
        for record in records
        if record.source_module.startswith(upstream_prefixes)
        if record.target.startswith("phospy.differential")
    )
    differential_targets = {
        record.target
        for record in records
        if record.source_module.startswith("phospy.differential")
        if record.target.startswith(
            (
                "phospy.datasets",
                "phospy.preprocessing",
                "phospy.preparation_records",
                "phospy.io",
                "phospy.statistics",
            )
        )
    }

    assert reverse == []
    assert differential_targets == {
        "phospy.datasets.derived_quantitative",
        "phospy.datasets.internal_view",
        "phospy.datasets.models",
        "phospy.datasets.processing_state",
        "phospy.datasets.trusted_construction_caveats",
        "phospy.preparation_records.batch_correction",
        "phospy.preparation_records.policies",
        "phospy.preparation_records.protein_aware",
        "phospy.statistics.multiple_testing",
    }


def test_result_table_validator_definitions_are_single_owner() -> None:
    names = {
        "require_boolean",
        "require_column_name",
        "require_differential_result_columns",
        "require_na_position",
        "require_non_negative_threshold",
        "require_numeric_result_column",
        "require_probability_threshold",
    }
    definitions = {name: [] for name in names}
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in names:
                definitions[node.name].append(path.relative_to(PROJECT_ROOT).as_posix())

    assert definitions == {
        name: ["src/phospy/differential/validation/result_tables.py"] for name in names
    }


def test_differential_has_no_legacy_route_registration_or_forwarder() -> None:
    offenders = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        source = path.read_text(encoding="utf-8")
        for old_module in OLD_DIFFERENTIAL_MODULES:
            if old_module in source:
                offenders.append(
                    f"{path.relative_to(PROJECT_ROOT).as_posix()} -> {old_module}"
                )

    assert offenders == []
