from __future__ import annotations

import dataclasses
import inspect
from collections.abc import Sequence

import phospy
import phospy.api as public_api
import phospy.api.results as api_results
import phospy.contracts.results as contract_results
import phospy.science.datasets.preprocessing.batch_correction as batch_correction
import phospy.science.transformations as transformations

_FORBIDDEN_PUBLIC_PARAMETER_EXACT = frozenset(
    {
        "_assume_owned",
        "assume_owned",
        "_owned",
        "owned",
        "copy",
        "copy_input",
        "copy_inputs",
        "copy_data",
        "no_copy",
        "skip_copy",
        "ownership_token",
        "transfer_token",
        "skip_validation",
        "_skip_validation",
        "disable_validation",
        "bypass_validation",
        "skip_fingerprint",
        "disable_fingerprint",
        "bypass_fingerprint",
        "trust_fingerprint",
        "suppress_warnings",
        "ignore_warnings",
        "disable_warnings",
        "silence_warnings",
        "validator",
        "request_validator",
        "source_validator",
        "config_validator",
        "interpreter",
        "executor",
        "path_reader",
        "source_reader",
        "batch_correction_runner",
        "internal_view",
    }
)
_FORBIDDEN_PUBLIC_PARAMETER_FRAGMENTS = (
    "assume_own",
    "owned_input",
    "ownership",
    "copy_bypass",
    "validation_bypass",
    "fingerprint_bypass",
    "suppress_warning",
    "ignore_warning",
    "disable_warning",
    "silence_warning",
    "internal_view",
)
_FORBIDDEN_PUBLIC_EXPORT_FRAGMENTS = (
    "Validator",
    "Interpreter",
    "Executor",
    "InternalView",
    "Ownership",
    "OwnershipToken",
    "OwnedFactory",
    "AssumeOwned",
)

_JSON_FIELD_RUNTIME_COVERAGE = {
    "phospy.contracts.results.base.ImporterMissingIntensitySummary."
    "missing_values_by_sample_id": (
        "recursive JSON freezing covered by exported JSON state tests"
    ),
    "phospy.contracts.results.base.ImporterMissingIntensitySummary."
    "missing_values_by_source_column": (
        "recursive JSON freezing covered by exported JSON state tests"
    ),
    "phospy.contracts.results.base.ImporterQualityReport.format_specific": (
        "recursive JSON freezing covered by exported JSON state tests"
    ),
    "phospy.contracts.results.base.PhosphositeImportResult.diagnostics": (
        "runtime public-boundary adversarial registry covers importer diagnostics"
    ),
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult.background_summary": (
        "runtime public-boundary adversarial registry covers enrichment summaries"
    ),
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult.diagnostics": (
        "runtime public-boundary adversarial registry covers enrichment diagnostics"
    ),
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult.method_metadata": (
        "runtime public-boundary adversarial registry covers enrichment metadata"
    ),
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
    "set_collection_summary": (
        "runtime public-boundary adversarial registry covers enrichment summaries"
    ),
    "phospy.contracts.results.kinase.KinaseWorkflowAttritionProvenance.metrics": (
        "runtime public-boundary adversarial registry covers kinase attrition evidence"
    ),
    "phospy.contracts.results.kinase.KinaseWorkflowAttritionProvenance.policy": (
        "runtime public-boundary adversarial registry covers kinase attrition policy"
    ),
    "phospy.contracts.results.kinase.KinaseWorkflowAttritionProvenance."
    "policy_violations": (
        "runtime public-boundary adversarial registry covers nested violations"
    ),
    "phospy.science.datasets.preprocessing.batch_correction."
    "BatchCorrectionResult.diagnostics": (
        "runtime public-boundary adversarial registry covers processing diagnostics"
    ),
    "phospy.science.datasets.preprocessing.protein_aware_preparation."
    "ProteinAwarePreparationReport.policy_parameters": (
        "runtime public-boundary adversarial registry covers processing policy JSON"
    ),
    "phospy.science.prediction.models.KinaseScoringResult.score_scale_metadata": (
        "runtime public-boundary adversarial registry covers intensity score metadata"
    ),
    "phospy.science.result_caveats.ResultCaveat.details": (
        "runtime public-boundary adversarial registry covers result caveat details"
    ),
    "phospy.science.transformations.models."
    "IntensityScaleEstablishmentProvenance.parameters": (
        "runtime public-boundary adversarial registry covers intensity-scale evidence"
    ),
}

_JSON_FIELD_EXEMPTIONS = {
    "phospy.contracts.requests.PhosphositeImportRequest.sample_intensity_columns": (
        "request input mapping, not retained result JSON state; importer validation "
        "normalizes it before constructing public results"
    ),
    "phospy.contracts.results.base.PhosphositeImportResult._sample_column_mapping": (
        "private scalar str->str map; public sample_column_mapping returns a copy"
    ),
    "phospy.science.configs.preprocessing.correction_missingness."
    "TemporaryImputationPolicy.method_parameters": (
        "normalized immutable tuple of scalar JSON pairs with method-specific "
        "validation"
    ),
    "phospy.science.datasets.preprocessing.control_sites.ControlSiteAnnotation."
    "metadata_missing_reason": (
        "normalized scalar str->str reason map; nested JSON state is not accepted"
    ),
    "phospy.science.datasets.preprocessing.control_sites.ControlSiteSourceMetadata."
    "metadata_missing_reason": (
        "normalized scalar str->str reason map; nested JSON state is not accepted"
    ),
    "phospy.science.design.models.SampleDesignRecord.covariates": (
        "normalized scalar covariate map; nested JSON state is not accepted"
    ),
    "phospy.science.differential.models.results.DifferentialAnalysisResult."
    "_contrast_tables": (
        "private DataFrame mapping; public table exports are covered by ownership "
        "probes"
    ),
    "phospy.science.differential.models.results.DifferentialAnalysisResult."
    "workflow_provenance": (
        "optional workflow-owned provenance payload; covered by provenance golden "
        "and differential workflow tests"
    ),
}


def _public_boundary_exports() -> dict[str, object]:
    exports: dict[str, object] = {}
    for module_name, module in (("phospy", phospy), ("phospy.api", public_api)):
        names = getattr(module, "__all__", ())
        assert isinstance(names, Sequence)
        for name in names:
            assert isinstance(name, str)
            assert not name.startswith("_"), f"{module_name}.{name}"
            assert hasattr(module, name), f"{module_name}.{name}"
            exports[f"{module_name}.{name}"] = getattr(module, name)
    return exports


def _is_forbidden_public_parameter(name: str) -> bool:
    normalized = name.lower()
    if name.startswith("_") or normalized in _FORBIDDEN_PUBLIC_PARAMETER_EXACT:
        return True
    if any(
        fragment in normalized for fragment in _FORBIDDEN_PUBLIC_PARAMETER_FRAGMENTS
    ):
        return True
    if any(
        fragment in normalized for fragment in ("validator", "interpreter", "executor")
    ):
        return True
    disabling_tokens = ("skip", "disable", "bypass", "ignore", "suppress", "silence")
    if "warning" in normalized and any(
        token in normalized for token in disabling_tokens
    ):
        return True
    if "validation" in normalized and any(
        token in normalized for token in disabling_tokens
    ):
        return True
    if "fingerprint" in normalized and any(
        token in normalized for token in (*disabling_tokens, "trust")
    ):
        return True
    if "check" in normalized and any(token in normalized for token in disabling_tokens):
        return True
    return False


def _forbidden_public_parameters(owner: object) -> list[str]:
    try:
        parameters = inspect.signature(owner).parameters
    except (TypeError, ValueError):
        return []
    return [
        name
        for name in parameters
        if name not in {"self", "cls"} and _is_forbidden_public_parameter(name)
    ]


def _json_like_public_model_field_ids() -> set[str]:
    classes: set[type[object]] = set()
    for module in (
        phospy,
        public_api,
        api_results,
        contract_results,
        batch_correction,
        transformations,
    ):
        for name in getattr(module, "__all__", ()):
            value = getattr(module, name, None)
            if inspect.isclass(value):
                classes.add(value)

    field_ids: set[str] = set()
    for owner in classes:
        if not dataclasses.is_dataclass(owner):
            continue
        for field in dataclasses.fields(owner):
            annotation = str(field.type)
            if not (
                "Mapping[" in annotation
                or "dict[" in annotation
                or "Dict[" in annotation
                or "Json" in annotation
            ):
                continue
            field_ids.add(f"{owner.__module__}.{owner.__qualname__}.{field.name}")
    return field_ids


def test_root_and_api_exported_signatures_have_no_private_boundary_controls() -> None:
    offenders: dict[str, list[str]] = {}
    for exported_name, exported in sorted(_public_boundary_exports().items()):
        if not (inspect.isclass(exported) or inspect.isfunction(exported)):
            continue
        forbidden = _forbidden_public_parameters(exported)
        if forbidden:
            offenders[exported_name] = forbidden

    assert offenders == {}


def test_root_and_api_exports_do_not_expose_private_boundary_roles() -> None:
    offenders: dict[str, str] = {}
    for exported_name in _public_boundary_exports():
        symbol_name = exported_name.rsplit(".", maxsplit=1)[1]
        for fragment in _FORBIDDEN_PUBLIC_EXPORT_FRAGMENTS:
            if fragment in symbol_name:
                offenders[exported_name] = fragment
                break

    assert offenders == {}


def test_exported_json_like_public_model_fields_are_inventoried() -> None:
    expected = set(_JSON_FIELD_RUNTIME_COVERAGE) | set(_JSON_FIELD_EXEMPTIONS)

    assert _json_like_public_model_field_ids() == expected
    assert all(reason.strip() for reason in _JSON_FIELD_RUNTIME_COVERAGE.values())
    assert all(reason.strip() for reason in _JSON_FIELD_EXEMPTIONS.values())
    assert len(_JSON_FIELD_EXEMPTIONS) < len(_JSON_FIELD_RUNTIME_COVERAGE)
