from __future__ import annotations

import ast
import inspect
import textwrap
from pathlib import Path

import pandas as pd
import pytest

import phospy
import phospy.advanced as advanced_api
import phospy.api as public_api
import phospy.phosphosites.contracts as identity_contracts
from phospy.activities.validation import (
    MethodQuantitativeInputValidator,
)
from phospy.datasets.construction.executor import DatasetBuildExecutor
from phospy.datasets.schemas import (
    SiteMetadataTable,
)
from phospy.differential.validation.dataset import (
    DifferentialDatasetEligibilityValidator,
)
from phospy.differential.validation.workflow import DifferentialAnalysisValidator
from phospy.enrichment.validation import EnrichmentWorkflowValidator
from phospy.errors.validation import DatasetValidationError
from phospy.kinase.resolved_validator import (
    ResolvedKinaseEligibilityValidator,
)
from phospy.kinase.scoring.contracts import (
    KinaseScoringModeInputContract,
)
from phospy.kinase.sequence_contracts import (
    kinase_sequence_context_contract,
)
from phospy.kinase.validator import (
    KinaseWorkflowValidator,
    kinase_profile_scoring_workflow_input_contract,
)
from phospy.phosphosites.contracts import (
    ANALYSIS_READY_DATASET_BASE_IDENTITY_CONTRACT,
    ANALYSIS_READY_IDENTITY_COLUMNS,
    SITE_SEQUENCE_COLUMN,
    SequenceContextRequirement,
)
from phospy.preprocessing.validation.quantitative_evidence import (
    AdditivePreprocessingScaleGuard,
)
from phospy.quantitative.validation import (
    WorkflowQuantitativeInputContract,
    WorkflowQuantitativeInputValidator,
)
from phospy.references.context_compatibility import (
    validate_reference_context_compatibility,
)
from phospy.signalome.validator import (
    SignalomeWorkflowValidator,
    signalome_workflow_input_contract,
)
from phospy.taxonomy import Organism
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
API_ROOT = SRC_ROOT / "phospy" / "api"
CONTRACTS_ROOT = SRC_ROOT / "phospy" / "contracts"
ENRICHMENT_ROOT = SRC_ROOT / "phospy" / "enrichment"
DOCS_ROOT = PROJECT_ROOT / "docs"
WORKFLOW_CAPABILITY_ROOTS = tuple(
    SRC_ROOT / "phospy" / name
    for name in ("differential", "enrichment", "kinase", "signalome")
)
DATASET_VALIDATION_ROOT = SRC_ROOT / "phospy" / "datasets" / "validation"
IO_IMPORTER_CONTRACTS_PATH = SRC_ROOT / "phospy" / "io" / "importers" / "contracts.py"
ENRICHMENT_REQUEST_CONTRACTS_PATH = ENRICHMENT_ROOT / "contracts" / "requests.py"
DATASET_BUILD_CONTRACTS_PATH = SRC_ROOT / "phospy" / "datasets" / "contracts.py"
DIFFERENTIAL_REQUEST_CONTRACTS_PATH = (
    SRC_ROOT / "phospy" / "differential" / "contracts.py"
)
KINASE_REQUEST_CONTRACTS_PATH = SRC_ROOT / "phospy" / "kinase" / "requests.py"
SIGNALOME_REQUEST_CONTRACTS_PATH = SRC_ROOT / "phospy" / "signalome" / "requests.py"
REQUEST_DTO_CONTRACT_PATHS = (
    IO_IMPORTER_CONTRACTS_PATH,
    ENRICHMENT_REQUEST_CONTRACTS_PATH,
    DATASET_BUILD_CONTRACTS_PATH,
    DIFFERENTIAL_REQUEST_CONTRACTS_PATH,
    KINASE_REQUEST_CONTRACTS_PATH,
    SIGNALOME_REQUEST_CONTRACTS_PATH,
)
DATASET_VALIDATION_PREFIX = "phospy.datasets.validation"
VALIDATION_PREFIX = "phospy.validation"
PRIVATE_IMPLEMENTATION_CLASS_SUFFIXES = ("Validator", "Interpreter", "Executor")
IDENTITY_CONTRACT_CLASS_DEFINITION_OWNERS = {
    "PhosphositeIdentityContract": (
        SRC_ROOT / "phospy" / "phosphosites" / "contracts.py"
    ),
    "SequenceContextContract": (
        SRC_ROOT / "phospy" / "phosphosites" / "sequence_context.py"
    ),
}

PUBLIC_BOUNDARY_MODULES = (
    ("phospy", phospy),
    ("phospy.api", public_api),
)

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


def _defined_public_validation_names(root: Path) -> set[str]:
    names: set[str] = set()
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and node.name.endswith("Validator"):
                names.add(node.name)
            if isinstance(node, ast.FunctionDef) and node.name.startswith(
                ("enforce_", "validate_")
            ):
                names.add(node.name)
    return names


def _defined_private_implementation_class_names(root: Path) -> set[str]:
    names: set[str] = set()
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name.endswith(
                PRIVATE_IMPLEMENTATION_CLASS_SUFFIXES
            ):
                names.add(node.name)
    return names


def _called_function_names(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return {
        call_name
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        if (call_name := _dotted_name(node.func)) is not None
    }


def _dotted_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        value_name = _dotted_name(node.value)
        if value_name is None:
            return None
        return f"{value_name}.{node.attr}"
    return None


def _node_mentions_any_name(node: ast.AST, names: set[str]) -> bool:
    dotted_name = _dotted_name(node)
    if dotted_name in names or (
        dotted_name is not None and dotted_name.rsplit(".", maxsplit=1)[-1] in names
    ):
        return True
    return any(
        _node_mentions_any_name(child, names) for child in ast.iter_child_nodes(node)
    )


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


def _assert_public_boundary_import_fails(module_name: str, symbol_name: str) -> None:
    with pytest.raises(ImportError):
        exec(f"from {module_name} import {symbol_name}", {})


def _public_boundary_leaks(symbol_names: set[str]) -> list[str]:
    leaks: list[str] = []
    for module_name, module in PUBLIC_BOUNDARY_MODULES:
        exported_names = set(getattr(module, "__all__", ()))
        for symbol_name in sorted(symbol_names):
            if symbol_name in exported_names or hasattr(module, symbol_name):
                leaks.append(f"{module_name}.{symbol_name}")
    return leaks


def _validation_owned_symbol_module_name(symbol: object) -> str:
    module = inspect.getmodule(symbol)
    assert module is not None
    return module.__name__


def _validation_owned_symbol_path(symbol: object) -> Path:
    return Path(inspect.getfile(symbol)).resolve()


def _class_definition_paths(class_name: str) -> set[Path]:
    paths: set[Path] = set()
    for path in PACKAGE_ROOT.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if any(
            isinstance(node, ast.ClassDef) and node.name == class_name
            for node in tree.body
        ):
            paths.add(path.resolve())
    return paths


def _minimal_analysis_ready_payload() -> dict[str, object]:
    site_index = protein_site_key_index(
        protein_identifiers=["MAPK14"],
        sites=["Y182"],
    )
    return {
        "phospho": pd.DataFrame(
            {"sample_a": [1.0], "sample_b": [1.5]},
            index=site_index.copy(),
        ),
        "site_metadata": pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": ["MAPK14;Y182;"],
                **site_key_context_columns(site_index),
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"],
            },
            index=site_index.copy(),
        ),
        "sample_metadata": None,
        "total": None,
        "comparisons": None,
        "organism": Organism.RAT,
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=False,
        ),
        "processing_state": supported_linear_processing_state(
            has_total_matrix=False,
        ),
    }


def test_only_builder_adapter_imports_private_dataset_validation_modules() -> None:
    offenders = sorted(
        str(path.relative_to(PROJECT_ROOT))
        for path in API_ROOT.rglob("*.py")
        if any(
            module == DATASET_VALIDATION_PREFIX
            or module.startswith(f"{DATASET_VALIDATION_PREFIX}.")
            for module in _imported_modules(path)
        )
        if path.name != "builders.py"
    )

    assert offenders == []


def test_public_boundaries_do_not_export_private_dataset_validation_names() -> None:
    private_dataset_validation_names = _defined_public_validation_names(
        DATASET_VALIDATION_ROOT
    )

    assert private_dataset_validation_names
    assert _public_boundary_leaks(private_dataset_validation_names) == []


def test_public_boundary_exports_do_not_resolve_to_validation_modules() -> None:
    offenders = sorted(
        f"{module_name}.{name}"
        for module_name, module in PUBLIC_BOUNDARY_MODULES
        for name in getattr(module, "__all__", ())
        if getattr(getattr(module, name, None), "__module__", "").startswith(
            VALIDATION_PREFIX
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


def test_public_boundaries_do_not_export_validators_interpreters_or_executors() -> None:
    offenders = sorted(
        f"{module_name}.{name}"
        for module_name, module in PUBLIC_BOUNDARY_MODULES
        for name in getattr(module, "__all__", ())
        if name.endswith(PRIVATE_IMPLEMENTATION_CLASS_SUFFIXES)
    )

    assert offenders == []


def test_workflow_implementation_classes_are_not_public_boundary_imports() -> None:
    private_implementation_names = set().union(
        *(
            _defined_private_implementation_class_names(root)
            for root in WORKFLOW_CAPABILITY_ROOTS
        )
    )

    assert private_implementation_names
    assert _public_boundary_leaks(private_implementation_names) == []


def test_selected_private_validation_names_cannot_be_imported_from_public_modules() -> (
    None
):
    selected_private_names = DATASET_VALIDATOR_EXPORT_NAMES | {
        "DifferentialAnalysisExecutor",
        "DifferentialAnalysisInterpreter",
        "DifferentialAnalysisValidator",
        "EnrichmentWorkflowExecutor",
        "EnrichmentWorkflowInterpreter",
        "EnrichmentWorkflowValidator",
        "KinaseWorkflowExecutor",
        "KinaseWorkflowInterpreter",
        "KinaseWorkflowValidator",
        "SignalomeWorkflowExecutor",
        "SignalomeWorkflowInterpreter",
        "SignalomeWorkflowValidator",
        "WorkflowQuantitativeInputContract",
        "WorkflowQuantitativeInputValidator",
        "kinase_profile_scoring_workflow_input_contract",
        "signalome_workflow_input_contract",
        "validate_reference_context_compatibility",
    }

    for module_name, _module in PUBLIC_BOUNDARY_MODULES:
        for symbol_name in selected_private_names:
            _assert_public_boundary_import_fails(module_name, symbol_name)


def test_workflow_validators_compose_shared_and_domain_validation() -> None:
    differential_source = inspect.getsource(DifferentialAnalysisValidator.run)
    differential_eligibility_source = inspect.getsource(
        DifferentialDatasetEligibilityValidator.run
    )
    enrichment_source = inspect.getsource(EnrichmentWorkflowValidator.run)
    kinase_source = inspect.getsource(KinaseWorkflowValidator.run)
    signalome_run_source = inspect.getsource(SignalomeWorkflowValidator.run)
    signalome_source = inspect.getsource(
        SignalomeWorkflowValidator._require_site_identity_and_protein_grouping_metadata
    )
    signalome_grouping_source = inspect.getsource(
        SignalomeWorkflowValidator._require_signalome_protein_grouping_metadata
    )

    for source in (differential_source, kinase_source, signalome_source):
        assert "enforce_workflow_site_identity_contract(" in source

    assert "validate_reference_context_compatibility(" in kinase_source
    assert "dataset_eligibility_validator.run(" in differential_source
    assert "design_validator.run(" in differential_source
    assert "technical_replicate_planner.run(" in differential_source
    assert "self._contract_validator.run(" in enrichment_source
    assert "self._quantitative_input_validator.run(" in differential_eligibility_source
    assert (
        "DIFFERENTIAL_LOG_ABUNDANCE_INPUT_CONTRACT" in differential_eligibility_source
    )
    assert "self._method_quantitative_input_validator.run(" in kinase_source
    assert "kinase_scoring_method_quantitative_input_contract(" in kinase_source
    assert "kinase_activity_method_quantitative_input_contract(" in kinase_source
    assert "kinase_profile_scoring_workflow_input_contract(" not in kinase_source
    assert "self._quantitative_input_validator.run(" in signalome_run_source
    assert "signalome_workflow_input_contract(" in signalome_run_source
    assert "enforce_localisation_requirement(" in kinase_source
    assert "enforce_localisation_requirement(" in signalome_source
    assert "enforce_required_non_empty_string_column(" in signalome_grouping_source


def test_differential_validator_collaborators_are_protocol_dispatched() -> None:
    source = textwrap.dedent(inspect.getsource(DifferentialAnalysisValidator.run))
    tree = ast.parse(source)
    concrete_collaborators = {
        "DifferentialDatasetEligibilityValidator",
        "ExperimentalDesignContractValidator",
        "TechnicalReplicateAggregationPlanner",
    }
    concrete_branches: list[str] = []

    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and _dotted_name(node.func) == "isinstance"
            and len(node.args) >= 2
            and _node_mentions_any_name(node.args[1], concrete_collaborators)
        ):
            concrete_branches.append(f"line {node.lineno}")

    assert concrete_branches == []
    assert "cast(Any" not in source


def test_additive_preprocessing_scale_policy_is_private_builder_owned() -> None:
    executor_source = inspect.getsource(DatasetBuildExecutor.run)
    executor_guard_source = inspect.getsource(
        DatasetBuildExecutor._validate_additive_preprocessing_scale_policy
    )
    workflow_sources = (
        inspect.getsource(DifferentialAnalysisValidator),
        inspect.getsource(KinaseWorkflowValidator),
        inspect.getsource(SignalomeWorkflowValidator),
    )

    assert (
        AdditivePreprocessingScaleGuard.__module__
        == "phospy.preprocessing.validation.quantitative_evidence"
    )
    assert not hasattr(public_api, "AdditivePreprocessingScaleGuard")
    assert "_validate_additive_preprocessing_scale_policy(" in executor_source
    assert "self._additive_preprocessing_scale_guard.run(" in executor_guard_source
    assert all(
        "AdditivePreprocessingScaleGuard" not in source for source in workflow_sources
    )
    assert all(
        "additive_preprocessing_scale" not in source for source in workflow_sources
    )


def test_signalome_result_identity_validation_is_science_owned_and_workflow_composed() -> (
    None
):
    from phospy.signalome import result_validation
    from phospy.signalome.result_assembly import SignalomeResultAssembler

    assembler_source = inspect.getsource(SignalomeResultAssembler.run)
    validation_source = inspect.getsource(
        result_validation.validate_signalome_result_site_level_identity
    )

    assert (
        result_validation.validate_signalome_result_site_level_identity.__module__
        == "phospy.signalome.result_validation"
    )
    assert "validate_signalome_result_site_level_identity(" in assembler_source
    assert "phospy.workflows" not in validation_source


def test_quantitative_validators_follow_domain_ownership() -> None:
    assert MethodQuantitativeInputValidator.__module__ == (
        "phospy.activities.validation"
    )
    assert WorkflowQuantitativeInputContract.__module__ == (
        "phospy.quantitative.validation"
    )
    assert WorkflowQuantitativeInputValidator.__module__ == (
        "phospy.quantitative.validation"
    )
    assert kinase_profile_scoring_workflow_input_contract.__module__ == (
        "phospy.kinase.validator"
    )
    assert signalome_workflow_input_contract.__module__ == (
        "phospy.signalome.validator"
    )


def test_phosphosite_identity_contracts_have_direct_owner_leaves() -> None:
    from phospy.phosphosites import contracts, dataset_identity, result_identity

    assert (
        validate_reference_context_compatibility.__module__
        == "phospy.references.context_compatibility"
    )
    assert (
        identity_contracts.PhosphositeIdentityContract
        is contracts.PhosphositeIdentityContract
    )
    assert (
        dataset_identity.enforce_site_key_column.__module__
        == "phospy.phosphosites.dataset_identity"
    )
    assert (
        result_identity.enforce_result_identity_metadata_coherence.__module__
        == "phospy.phosphosites.result_identity"
    )
    for name in (
        "REFERENCE_CONTEXT_IDENTITY_FIELDS",
        "REFERENCE_CONTEXT_UNKNOWN_CAVEAT_CODE",
        "ReferenceContextCompatibilityWarning",
        "validate_reference_context_compatibility",
    ):
        assert not hasattr(identity_contracts, name)


def test_former_phosphosite_identity_forwarding_routes_are_absent() -> None:
    for relative_path in (
        "science/sites/identity_columns.py",
        "science/sites/identity_contracts.py",
        "validation/identity_contracts.py",
        "validation/datasets/protein_scoped_site_identity.py",
    ):
        assert not (SRC_ROOT / "phospy" / relative_path).exists()


def test_identity_validators_are_not_exposed_from_advanced_api() -> None:
    validator_names = {
        "enforce_analysis_ready_site_key_index",
        "enforce_display_id_column",
        "enforce_phosphosite_identity_contract",
        "enforce_result_identity_metadata_coherence",
        "enforce_site_key_column",
        "enforce_site_key_matches_metadata",
        "validate_reference_context_compatibility",
    }
    exported_names = set(getattr(advanced_api, "__all__", ()))

    assert validator_names.isdisjoint(exported_names)
    assert not any(hasattr(advanced_api, name) for name in validator_names)


def test_identity_contract_classes_have_single_concrete_owner() -> None:
    observed = {
        class_name: _class_definition_paths(class_name)
        for class_name in IDENTITY_CONTRACT_CLASS_DEFINITION_OWNERS
    }

    expected = {
        class_name: {owner.resolve()}
        for class_name, owner in IDENTITY_CONTRACT_CLASS_DEFINITION_OWNERS.items()
    }
    assert observed == expected


def test_reference_context_policy_is_not_exported_from_private_validation_api() -> None:
    assert "ReferenceContextCompatibilityPolicy" not in identity_contracts.__all__
    assert not hasattr(identity_contracts, "ReferenceContextCompatibilityPolicy")


def test_request_dto_modules_do_not_import_validation_domains() -> None:
    offenders = sorted(
        f"{path.relative_to(PROJECT_ROOT)} imports {module}"
        for path in REQUEST_DTO_CONTRACT_PATHS
        for module in _imported_modules(path)
        if module == VALIDATION_PREFIX or module.startswith(f"{VALIDATION_PREFIX}.")
    )

    assert offenders == []


def test_contract_modules_do_not_import_validation_domains() -> None:
    offenders = sorted(
        f"{path.relative_to(PROJECT_ROOT)} imports {module}"
        for path in CONTRACTS_ROOT.rglob("*.py")
        for module in _imported_modules(path)
        if module == VALIDATION_PREFIX or module.startswith(f"{VALIDATION_PREFIX}.")
    )

    assert offenders == []


def test_public_constructors_do_not_raise_workflow_boundary_exceptions() -> None:
    offenders = sorted(
        f"{name} ({obj.__module__})"
        for name in public_api.__all__
        if inspect.isclass(obj := getattr(public_api, name, None))
        and "__post_init__" in getattr(obj, "__dict__", {})
        if "WorkflowValidationError" in inspect.getsource(obj.__post_init__)
        or "WorkflowBoundaryError" in inspect.getsource(obj.__post_init__)
    )

    assert offenders == []


def test_request_dto_modules_do_not_call_dataset_or_workflow_validators() -> None:
    offenders = sorted(
        f"{path.relative_to(PROJECT_ROOT)} calls {call_name}"
        for path in REQUEST_DTO_CONTRACT_PATHS
        for call_name in _called_function_names(path)
        if call_name.endswith("Validator")
        or call_name.startswith(("enforce_", "validate_"))
        or ".enforce_" in call_name
        or ".validate_" in call_name
    )

    assert offenders == []


def test_analysis_ready_dataset_site_sequence_requirement_is_boundary_owned() -> None:
    site_metadata_source = inspect.getsource(SiteMetadataTable._validate_frame)

    assert SITE_SEQUENCE_COLUMN in ANALYSIS_READY_IDENTITY_COLUMNS
    assert (
        ANALYSIS_READY_DATASET_BASE_IDENTITY_CONTRACT.sequence_context
        is SequenceContextRequirement.PRESENT
    )
    assert 'column_name="site_sequence"' in site_metadata_source
    assert "validate_site_sequence_column(" in site_metadata_source


def test_analysis_ready_dataset_still_requires_site_sequence() -> None:
    payload = _minimal_analysis_ready_payload()
    site_metadata = payload["site_metadata"]
    assert isinstance(site_metadata, pd.DataFrame)
    payload["site_metadata"] = site_metadata.drop(columns=["site_sequence"])

    with pytest.raises(DatasetValidationError, match="site_sequence"):
        trusted_analysis_ready_dataset_from_tables(**payload)


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
    assert (
        KinaseScoringModeInputContract.__module__ == "phospy.kinase.scoring.contracts"
    )


def test_workflow_request_dtos_do_not_own_domain_validation() -> None:
    request_classes = [
        class_node
        for path in REQUEST_DTO_CONTRACT_PATHS
        for class_node in _request_dto_classes(
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        )
    ]
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
    request_source = "\n".join(
        path.read_text(encoding="utf-8") for path in REQUEST_DTO_CONTRACT_PATHS
    )

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
