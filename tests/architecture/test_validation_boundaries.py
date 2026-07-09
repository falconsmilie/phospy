from __future__ import annotations

import ast
import inspect
from pathlib import Path

import phospy
import phospy.api as public_api
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.resolved_validator import (
    ResolvedKinaseEligibilityValidator,
)
from phospy.workflows.kinase.scoring_mode_contracts import (
    KinaseScoringModeInputContract,
)
from phospy.workflows.kinase.sequence_contracts import (
    kinase_sequence_context_contract,
)
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
API_ROOT = SRC_ROOT / "phospy" / "api"
DOCS_ROOT = PROJECT_ROOT / "docs"
REQUEST_CONTRACTS_PATH = SRC_ROOT / "phospy" / "contracts" / "requests.py"
DATASET_VALIDATION_PREFIX = "phospy.validation.datasets"

DATASET_VALIDATOR_EXPORT_NAMES = frozenset(
    {
        "AnalysisReadyDatasetModelBoundaryValidator",
        "BatchCorrectionAdequacyValidator",
        "DISPLAY_SITE_CONTEXT_COLUMNS",
        "DatasetInputSourceValidator",
        "DatasetPreprocessingConfigValidator",
        "PhosphositeImportRequestValidator",
        "enforce_display_id_column",
        "enforce_site_key_column",
        "enforce_site_key_index",
        "enforce_site_key_matches_metadata",
        "enforce_unique_display_site_identity_rows",
        "enforce_unique_site_key_identity",
    }
)

EXPECTED_REQUEST_DTO_NAMES = frozenset(
    {
        "DatasetBuildRequest",
        "DifferentialAnalysisRequest",
        "EnrichmentWorkflowRequest",
        "KinaseWorkflowRequest",
        "PhosphositeImportRequest",
        "SignalomeWorkflowRequest",
    }
)

ALLOWED_REQUEST_POST_INIT_CALLS = frozenset(
    {
        "bool",
        "dict",
        "float",
        "frozenset",
        "int",
        "isinstance",
        "list",
        "object.__setattr__",
        "set",
        "str",
        "tuple",
    }
)

NONTRIVIAL_REQUEST_POST_INIT_NODES = (
    ast.Assert,
    ast.AsyncFor,
    ast.AsyncWith,
    ast.For,
    ast.Import,
    ast.ImportFrom,
    ast.Match,
    ast.Raise,
    ast.Try,
    ast.While,
    ast.With,
)


def _imported_modules(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.extend(alias.name for alias in node.names)
        if isinstance(node, ast.ImportFrom) and node.module is not None:
            modules.append(node.module)
    return modules


def _dotted_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        value_name = _dotted_name(node.value)
        if value_name is None:
            return None
        return f"{value_name}.{node.attr}"
    return None


def _has_dataclass_decorator(node: ast.ClassDef) -> bool:
    for decorator in node.decorator_list:
        decorated = decorator.func if isinstance(decorator, ast.Call) else decorator
        name = _dotted_name(decorated)
        if name == "dataclass" or (name is not None and name.endswith(".dataclass")):
            return True
    return False


def _request_dto_classes(tree: ast.Module) -> list[ast.ClassDef]:
    return [
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef)
        and node.name.endswith("Request")
        and _has_dataclass_decorator(node)
    ]


def _method_named(node: ast.ClassDef, method_name: str) -> ast.FunctionDef | None:
    for item in node.body:
        if isinstance(item, ast.FunctionDef) and item.name == method_name:
            return item
    return None


def _self_attribute_path(node: ast.Attribute) -> tuple[str, ...]:
    parts: list[str] = []
    current: ast.AST = node
    while isinstance(current, ast.Attribute):
        parts.append(current.attr)
        current = current.value
    if isinstance(current, ast.Name) and current.id == "self":
        return tuple(reversed(parts))
    return ()


def _is_allowed_request_post_init_call(call_name: str | None) -> bool:
    if call_name in ALLOWED_REQUEST_POST_INIT_CALLS:
        return True
    if call_name is None or "." in call_name:
        return False

    return call_name[:1].isupper() and not call_name.endswith("Validator")


def _nontrivial_request_post_init_reasons(
    class_name: str, method: ast.FunctionDef
) -> list[str]:
    reasons: set[str] = set()

    for node in ast.walk(method):
        if isinstance(node, NONTRIVIAL_REQUEST_POST_INIT_NODES):
            reasons.add(f"{class_name}.__post_init__ uses {type(node).__name__}")
        if isinstance(node, ast.Call):
            call_name = _dotted_name(node.func)
            if not _is_allowed_request_post_init_call(call_name):
                reasons.add(
                    f"{class_name}.__post_init__ calls {call_name or '<dynamic call>'}"
                )
        if isinstance(node, ast.Compare) and any(
            not isinstance(operator, (ast.Is, ast.IsNot)) for operator in node.ops
        ):
            reasons.add(f"{class_name}.__post_init__ compares domain values")
        if isinstance(node, ast.Attribute):
            path = _self_attribute_path(node)
            if len(path) > 1:
                reasons.add(
                    f"{class_name}.__post_init__ reads nested scientific state "
                    f"self.{'.'.join(path)}"
                )

    return sorted(reasons)


def test_api_package_does_not_import_private_dataset_validation_modules() -> None:
    offenders = sorted(
        str(path.relative_to(PROJECT_ROOT))
        for path in API_ROOT.rglob("*.py")
        if any(
            module == DATASET_VALIDATION_PREFIX
            or module.startswith(f"{DATASET_VALIDATION_PREFIX}.")
            for module in _imported_modules(path)
        )
    )

    assert offenders == []


def test_api_all_has_no_dataset_validator_exports() -> None:
    explicit_name_leaks = sorted(
        DATASET_VALIDATOR_EXPORT_NAMES.intersection(public_api.__all__)
    )
    validator_name_leaks = sorted(
        name for name in public_api.__all__ if name.endswith("Validator")
    )
    dataset_module_leaks = sorted(
        name
        for name in public_api.__all__
        if getattr(getattr(public_api, name), "__module__", "").startswith(
            DATASET_VALIDATION_PREFIX
        )
    )

    assert explicit_name_leaks == []
    assert validator_name_leaks == []
    assert dataset_module_leaks == []


def test_workflow_validators_compose_shared_and_domain_validation() -> None:
    differential_source = inspect.getsource(DifferentialAnalysisValidator.run)
    kinase_source = inspect.getsource(KinaseWorkflowValidator.run)
    signalome_source = inspect.getsource(
        SignalomeWorkflowValidator._require_site_identity_and_protein_grouping_metadata
    )
    signalome_grouping_source = inspect.getsource(
        SignalomeWorkflowValidator._require_signalome_protein_grouping_metadata
    )

    for source in (differential_source, kinase_source, signalome_source):
        assert "enforce_workflow_site_identity_contract(" in source

    assert "validate_reference_context_compatibility(" in kinase_source
    assert "self._dataset_eligibility_validator.run(" in differential_source
    assert "self._design_validator.run(" in differential_source
    assert "self._technical_replicate_planner.run(" in differential_source
    assert "enforce_localisation_requirement(" in kinase_source
    assert "enforce_localisation_requirement(" in signalome_source
    assert "enforce_required_non_empty_string_column(" in signalome_grouping_source


def test_kinase_scoring_mode_input_requirements_are_contract_owned() -> None:
    validator_source = inspect.getsource(KinaseWorkflowValidator.run)
    resolved_validator_source = inspect.getsource(
        ResolvedKinaseEligibilityValidator.run
    )
    sequence_contract_source = inspect.getsource(kinase_sequence_context_contract)

    assert "kinase_scoring_mode_input_contract(" in validator_source
    assert "kinase_scoring_mode_input_contract(" in resolved_validator_source
    assert "kinase_scoring_mode_input_contract(" in sequence_contract_source
    assert "KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY" not in validator_source


def test_kinase_scoring_mode_input_contract_remains_internal() -> None:
    assert not hasattr(public_api, "KinaseScoringModeInputContract")
    assert not hasattr(phospy, "KinaseScoringModeInputContract")
    assert KinaseScoringModeInputContract.__module__.startswith(
        "phospy.workflows.kinase."
    )


def test_workflow_request_dtos_do_not_own_domain_validation() -> None:
    tree = ast.parse(
        REQUEST_CONTRACTS_PATH.read_text(encoding="utf-8"),
        filename=str(REQUEST_CONTRACTS_PATH),
    )
    request_classes = _request_dto_classes(tree)
    discovered_request_names = {node.name for node in request_classes}
    offenders = [
        reason
        for class_node in request_classes
        if (post_init := _method_named(class_node, "__post_init__")) is not None
        for reason in _nontrivial_request_post_init_reasons(class_node.name, post_init)
    ]

    assert EXPECTED_REQUEST_DTO_NAMES.issubset(discovered_request_names)
    assert offenders == []


def test_dataset_request_dtos_do_not_call_private_dataset_validators() -> None:
    request_source = REQUEST_CONTRACTS_PATH.read_text(encoding="utf-8")

    assert DATASET_VALIDATION_PREFIX not in request_source
    for symbol_name in DATASET_VALIDATOR_EXPORT_NAMES:
        assert symbol_name not in request_source


def test_docs_do_not_present_dataset_validators_as_user_import_path() -> None:
    docs_with_private_imports = sorted(
        str(path.relative_to(PROJECT_ROOT))
        for path in DOCS_ROOT.rglob("*.md")
        if "from phospy.validation.datasets" in path.read_text(encoding="utf-8")
        or "import phospy.validation.datasets" in path.read_text(encoding="utf-8")
    )
    ownership_doc = (DOCS_ROOT / "validation-ownership.md").read_text(encoding="utf-8")

    assert docs_with_private_imports == []
    assert "Users validate dataset inputs by constructing datasets" in ownership_doc
    assert "not supported user entrypoints" in ownership_doc
