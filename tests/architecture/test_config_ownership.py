from __future__ import annotations

import ast
import importlib
from pathlib import Path
from typing import get_type_hints

from phospy.advanced import DifferentialAnalysisConfig
from phospy.workflows.differential.interpreter import _resolve_execution_config
from phospy.workflows.differential.models import ResolvedDifferentialExecutionConfig
from phospy.workflows.kinase.contracts import ResolvedKinaseExecutionConfig
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.signalome.contracts import ResolvedSignalomeExecutionConfig
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PHOSPY_ROOT = SRC_ROOT / "phospy"
CONTRACT_CONFIG_ROOT = PHOSPY_ROOT / "contracts" / "configs"
SCIENCE_CONFIG_ROOT = PHOSPY_ROOT / "science" / "configs"
SCIENCE_ROOT = PHOSPY_ROOT / "science"

PRIVATE_DUPLICATE_ALLOWLIST_WITH_REASONS: dict[str, str] = {}

FULL_MODULE_REEXPORTS = (
    ("phospy.contracts.configs._validation", "phospy.policies.config_values"),
    ("phospy.science.configs._validation", "phospy.policies.config_values"),
    ("phospy.contracts.configs.dataset", "phospy.science.configs.dataset"),
    ("phospy.contracts.configs.localisation", "phospy.science.configs.localisation"),
    (
        "phospy.contracts.configs.reference_context",
        "phospy.science.configs.reference_context",
    ),
    (
        "phospy.contracts.configs.preprocessing._validation",
        "phospy.science.configs.preprocessing.validation",
    ),
    (
        "phospy.science.configs.preprocessing._validation",
        "phospy.science.configs.preprocessing.validation",
    ),
    (
        "phospy.contracts.configs.preprocessing.batch_correction",
        "phospy.science.configs.preprocessing.batch_correction",
    ),
    (
        "phospy.contracts.configs.preprocessing.comparisons",
        "phospy.science.configs.preprocessing.comparisons",
    ),
    (
        "phospy.contracts.configs.preprocessing.control_sites",
        "phospy.science.configs.preprocessing.control_sites",
    ),
    (
        "phospy.contracts.configs.preprocessing.correction_missingness",
        "phospy.science.configs.preprocessing.correction_missingness",
    ),
    (
        "phospy.contracts.configs.preprocessing.coverage_filter",
        "phospy.science.configs.preprocessing.coverage_filter",
    ),
    (
        "phospy.contracts.configs.preprocessing.intensity_transform",
        "phospy.science.configs.preprocessing.intensity_transform",
    ),
    (
        "phospy.contracts.configs.preprocessing.internal_batch_correction",
        "phospy.science.configs.preprocessing.internal_batch_correction",
    ),
    (
        "phospy.contracts.configs.preprocessing.localisation",
        "phospy.science.configs.preprocessing.localisation",
    ),
    (
        "phospy.contracts.configs.preprocessing.missing_data",
        "phospy.science.configs.preprocessing.missing_data",
    ),
    (
        "phospy.contracts.configs.preprocessing.normalisation",
        "phospy.science.configs.preprocessing.normalisation",
    ),
    (
        "phospy.contracts.configs.preprocessing.site_matrix",
        "phospy.science.configs.preprocessing.site_matrix",
    ),
    (
        "phospy.contracts.configs.preprocessing.site_sequence",
        "phospy.science.configs.preprocessing.site_sequence",
    ),
    (
        "phospy.contracts.configs.preprocessing.total_protein",
        "phospy.science.configs.preprocessing.total_protein",
    ),
)

COMMON_NAME_REEXPORTS = (
    ("phospy.contracts.configs.differential", "phospy.science.configs.differential"),
    ("phospy.contracts.configs.enrichment", "phospy.science.configs.enrichment"),
    ("phospy.contracts.configs.kinase", "phospy.science.configs.kinase"),
    ("phospy.contracts.configs.prediction", "phospy.science.configs.prediction"),
    ("phospy.contracts.configs.signalome", "phospy.science.configs.signalome"),
)

PUBLIC_WORKFLOW_CONFIG_DTO_NAMES = frozenset(
    {
        "DifferentialAnalysisConfig",
        "EnrichmentConfig",
        "KinaseActivityConfig",
        "KinaseAttritionPolicy",
        "KinasePredictionConfig",
        "KinaseScoringConfig",
        "MultipleTestingConfig",
        "SignalomeClusteringConfig",
        "SignalomeConfig",
        "SignalomeOutputConfig",
        "SignalomePerformanceConfig",
        "SignalomeScientificConfig",
        "SignalomeValidationConfig",
    }
)


def test_config_trees_do_not_define_same_exported_class_enum_or_function() -> None:
    definitions: dict[str, list[tuple[str, Path]]] = {}
    for owner, root in (
        ("contracts", CONTRACT_CONFIG_ROOT),
        ("science", SCIENCE_CONFIG_ROOT),
    ):
        for path in sorted(root.rglob("*.py")):
            for symbol_name in _exported_definition_names(path):
                definitions.setdefault(symbol_name, []).append((owner, path))

    offenders = {
        symbol_name: locations
        for symbol_name, locations in definitions.items()
        if {owner for owner, _ in locations} == {"contracts", "science"}
        if symbol_name not in PRIVATE_DUPLICATE_ALLOWLIST_WITH_REASONS
    }

    assert offenders == {}, _format_duplicate_definitions(offenders)
    assert all(PRIVATE_DUPLICATE_ALLOWLIST_WITH_REASONS.values())


def test_deliberate_science_owned_reexports_preserve_object_identity() -> None:
    for reexport_module_name, owner_module_name in FULL_MODULE_REEXPORTS:
        reexport_module = importlib.import_module(reexport_module_name)
        owner_module = importlib.import_module(owner_module_name)
        for symbol_name in reexport_module.__all__:
            assert getattr(reexport_module, symbol_name) is getattr(
                owner_module,
                symbol_name,
            ), symbol_name

    for reexport_module_name, owner_module_name in COMMON_NAME_REEXPORTS:
        reexport_module = importlib.import_module(reexport_module_name)
        owner_module = importlib.import_module(owner_module_name)
        common_names = sorted(set(reexport_module.__all__) & set(owner_module.__all__))
        assert common_names
        for symbol_name in common_names:
            assert getattr(reexport_module, symbol_name) is getattr(
                owner_module,
                symbol_name,
            ), symbol_name


def test_science_package_does_not_import_contracts() -> None:
    offenders: list[str] = []
    for path in sorted(SCIENCE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported_module in _imported_modules(tree):
            if imported_module == "phospy.contracts" or imported_module.startswith(
                "phospy.contracts."
            ):
                offenders.append(f"{_module_name(path)} -> {imported_module}")

    assert offenders == []


def test_public_api_config_imports_remain_available() -> None:
    public_configs = importlib.import_module("phospy.api.configs")

    missing = [
        symbol_name
        for symbol_name in public_configs.__all__
        if not hasattr(public_configs, symbol_name)
    ]

    assert missing == []


def test_workflow_interpreters_return_distinct_resolved_config_models() -> None:
    kinase_hints = get_type_hints(KinaseWorkflowInterpreter._resolve_execution_config)
    signalome_hints = get_type_hints(
        SignalomeWorkflowInterpreter._resolve_execution_config
    )

    assert kinase_hints["return"] is ResolvedKinaseExecutionConfig
    assert signalome_hints["return"] is ResolvedSignalomeExecutionConfig

    differential_config = DifferentialAnalysisConfig()
    resolved = _resolve_execution_config(differential_config)

    assert isinstance(resolved, ResolvedDifferentialExecutionConfig)
    assert resolved is not differential_config
    assert resolved.empirical_bayes is differential_config.empirical_bayes
    assert (
        resolved.multiple_testing_method == differential_config.multiple_testing.method
    )


def test_numerical_science_functions_do_not_accept_public_workflow_config_dtos() -> (
    None
):
    offenders: list[str] = []
    for path in sorted(SCIENCE_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
                continue
            for argument in (
                *node.args.posonlyargs,
                *node.args.args,
                *node.args.kwonlyargs,
            ):
                if argument.annotation is None:
                    continue
                annotation = ast.unparse(argument.annotation)
                matched_names = _public_dto_names_in_annotation(annotation)
                for matched_name in sorted(matched_names):
                    offenders.append(
                        f"{_module_name(path)}.{node.name}({argument.arg}: "
                        f"{matched_name})"
                    )

    assert offenders == []


def _exported_definition_names(path: Path) -> tuple[str, ...]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    names: list[str] = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef | ast.FunctionDef):
            if node.name.startswith("_"):
                continue
            names.append(node.name)
    return tuple(names)


def _imported_modules(tree: ast.AST) -> tuple[str, ...]:
    modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            modules.append(node.module)
    return tuple(modules)


def _public_dto_names_in_annotation(annotation: str) -> set[str]:
    return {
        name
        for name in PUBLIC_WORKFLOW_CONFIG_DTO_NAMES
        if name in annotation.replace("'", "").replace('"', "")
    }


def _module_name(path: Path) -> str:
    return "phospy." + ".".join(path.relative_to(PHOSPY_ROOT).with_suffix("").parts)


def _format_duplicate_definitions(
    offenders: dict[str, list[tuple[str, Path]]],
) -> str:
    lines: list[str] = []
    for symbol_name, locations in sorted(offenders.items()):
        lines.append(symbol_name)
        for owner, path in locations:
            lines.append(f"  {owner}: {path.relative_to(PROJECT_ROOT)}")
    return "\n".join(lines)
