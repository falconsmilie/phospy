"""Package-private deprecation governance for retained PhosPy compatibility."""

from __future__ import annotations

import importlib
import warnings
from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

from phospy._api_inventory import (
    ADVANCED_CONFIG_API,
    ADVANCED_PUBLIC_API,
    ADVANCED_RESULT_API,
    API_COMPATIBILITY_INTRODUCED_VERSION,
    API_COMPATIBILITY_PLANNED_REMOVAL_VERSION,
    COMPATIBILITY_CONFIG_MODULES,
    CONFIG_COMPATIBILITY_ADVANCED_ROUTE_OVERRIDES,
    REQUEST_COMPATIBILITY_ADVANCED_API,
    STABLE_CONFIG_API,
    STABLE_REQUEST_API,
)


class PhosPyDeprecationWarning(DeprecationWarning):
    """Warning category for retained user-facing PhosPy deprecations."""


DeprecationKind = Literal[
    "argument-alias",
    "class-alias",
    "classmethod-alias",
    "constructor-behaviour",
    "function-alias",
    "import-route",
    "method-alias",
    "property-alias",
    "value-alias",
]
DeprecationStability = Literal["stable", "advanced", "unsupported", "internal"]


@dataclass(frozen=True, slots=True)
class RetainedDeprecation:
    """Metadata for one retained deprecation during its compatibility window."""

    identifier: str
    kind: DeprecationKind
    owner_module: str
    deprecated: str
    replacement: str
    introduced_version: str
    planned_removal_version: str
    stability: DeprecationStability
    replacement_module: str
    replacement_name: str
    summary: str = ""
    deprecated_module: str | None = None
    deprecated_name: str | None = None


_DEFAULT_INTRODUCED_VERSION = "1.6.0"
_DEFAULT_PLANNED_REMOVAL_VERSION = "2.0.0"


def api_compatibility_deprecation_id(*, old_module: str, name: str) -> str:
    """Return the shared deprecation identifier for an API compatibility route."""

    return f"api-import:{old_module}.{name}"


def deprecation_record(identifier: str) -> RetainedDeprecation | None:
    """Return registered deprecation metadata by identifier."""

    return _deprecation_map().get(identifier)


def retained_deprecations() -> tuple[RetainedDeprecation, ...]:
    """Return all retained package-level deprecations."""

    return tuple(_deprecation_map().values())


def compatibility_deprecation_record(
    *,
    old_module: str,
    name: str,
) -> RetainedDeprecation | None:
    """Return shared metadata for one deprecated API compatibility export."""

    return deprecation_record(
        api_compatibility_deprecation_id(old_module=old_module, name=name)
    )


def compatibility_deprecation_records() -> tuple[RetainedDeprecation, ...]:
    """Return all retained deprecated API compatibility export records."""

    return tuple(
        record for record in retained_deprecations() if record.kind == "import-route"
    )


def warn_deprecated(identifier: str, *, stacklevel: int) -> None:
    """Emit a reviewed package-specific deprecation warning."""

    if stacklevel < 1:
        raise ValueError("deprecation warning stacklevel must be at least 1")
    record = deprecation_record(identifier)
    if record is None:
        raise KeyError(f"unknown PhosPy deprecation identifier: {identifier}")
    warnings.warn(
        _format_deprecation_warning(record),
        PhosPyDeprecationWarning,
        stacklevel=stacklevel,
    )


def _format_deprecation_warning(record: RetainedDeprecation) -> str:
    summary = record.summary or f"{record.deprecated} is deprecated"
    summary = summary.rstrip()
    separator = " " if summary.endswith(".") else "; "
    return (
        f"{summary}{separator}use {record.replacement}. "
        f"This deprecation was introduced in PhosPy {record.introduced_version} "
        f"and is planned for removal in PhosPy {record.planned_removal_version}. "
        f"Owner: {record.owner_module}; stability: {record.stability}."
    )


@lru_cache(maxsize=1)
def _deprecation_map() -> dict[str, RetainedDeprecation]:
    entries: dict[str, RetainedDeprecation] = {}
    for record in (*_static_deprecations(), *_api_compatibility_deprecations()):
        if record.identifier in entries:
            raise RuntimeError(
                f"duplicate PhosPy deprecation identifier: {record.identifier}"
            )
        entries[record.identifier] = record
    return entries


def _static_deprecations() -> tuple[RetainedDeprecation, ...]:
    records: list[RetainedDeprecation] = [
        _record(
            identifier="science.differential.DifferentialAnalysis",
            kind="class-alias",
            owner_module="phospy.science.differential.public",
            deprecated="DifferentialAnalysis",
            replacement="`from phospy import DifferentialAnalysisWorkflow`",
            stability="stable",
            replacement_module="phospy",
            replacement_name="DifferentialAnalysisWorkflow",
            summary=(
                "DifferentialAnalysis is deprecated; use "
                "DifferentialAnalysisWorkflow from top-level phospy"
            ),
        ),
        _record(
            identifier="contracts.kinase.KinaseScoringConfig.default",
            kind="classmethod-alias",
            owner_module="phospy.contracts.configs.kinase",
            deprecated="KinaseScoringConfig.default()",
            replacement="KinaseScoringConfig.exploratory()",
            stability="advanced",
            replacement_module="phospy.advanced",
            replacement_name="KinaseScoringConfig",
            summary=(
                "KinaseScoringConfig.default() is deprecated because the name "
                "is ambiguous"
            ),
        ),
        _record(
            identifier="science.kinase.scoring_mode.kinase_library_motif",
            kind="value-alias",
            owner_module="phospy.science.configs.kinase",
            deprecated="'kinase_library_motif'",
            replacement="'kinase_library_contextual_motif'",
            stability="advanced",
            replacement_module="phospy.advanced.configs",
            replacement_name="KinaseScoringMode",
            summary=(
                "'kinase_library_motif' is deprecated because it requires "
                "contextual profile/reference-substrate eligibility. Use "
                "'kinase_library_contextual_motif' for the current behavior."
            ),
        ),
        _record(
            identifier="prediction.motif_library.bare_sequence",
            kind="value-alias",
            owner_module=("phospy.science.prediction.motif_scoring.library_validation"),
            deprecated="bare motif sequence strings in motif_sequences",
            replacement=(
                "ExplicitMotifSequence values or mapping entries with "
                "reference_id, site_id, kinase, and sequence fields"
            ),
            stability="unsupported",
            replacement_module=("phospy.science.prediction.motif_scoring.models"),
            replacement_name="ExplicitMotifSequence",
            summary=(
                "Bare motif sequence strings in motif_sequences are deprecated "
                "because they omit stable reference and site identity metadata "
                "needed for reproducible motif-library validation"
            ),
        ),
        _record(
            identifier="activities.inputs.missing_activity_input",
            kind="constructor-behaviour",
            owner_module="phospy.science.activities.inputs",
            deprecated=(
                "KinaseActivityInputs construction without typed activity_input "
                "semantics"
            ),
            replacement="ActivityInputMatrix.sample_level_abundance(...)",
            stability="unsupported",
            replacement_module="phospy.science.activities.semantics",
            replacement_name="ActivityInputMatrix",
            summary=(
                "KinaseActivityInputs constructed without typed activity_input "
                "semantics is deprecated"
            ),
        ),
        _record(
            identifier="activities.ssgsea.effect_matrix_dataframe",
            kind="argument-alias",
            owner_module=(
                "phospy.science.activities.methods.ssgsea_substrate_enrichment"
            ),
            deprecated="raw DataFrame effect_matrix",
            replacement=(
                "ActivityInputMatrix.contrast_log_fold_change(...) or "
                "ActivityInputMatrix.standardised_effect(...)"
            ),
            stability="unsupported",
            replacement_module="phospy.science.activities.semantics",
            replacement_name="ActivityInputMatrix",
            summary="Passing a raw DataFrame as effect_matrix is deprecated",
        ),
        _record(
            identifier="activities.result.missing_semantics",
            kind="constructor-behaviour",
            owner_module="phospy.science.activities.result_validation",
            deprecated=(
                "KinaseActivityResult construction without explicit activity "
                "input semantics"
            ),
            replacement=(
                "KinaseActivityResult(..., input_semantics=..., profile_metadata=...)"
            ),
            stability="stable",
            replacement_module="phospy.api",
            replacement_name="KinaseActivityResult",
            summary=(
                "KinaseActivityResult constructed without explicit activity "
                "input semantics is deprecated"
            ),
        ),
        _record(
            identifier="activities.result.activity_scores",
            kind="property-alias",
            owner_module="phospy.science.activities.results",
            deprecated="KinaseActivityResult.activity_scores",
            replacement="KinaseActivityResult.activity_matrix",
            stability="stable",
            replacement_module="phospy.api",
            replacement_name="KinaseActivityResult",
        ),
        _record(
            identifier="activities.result.weighted_activity",
            kind="property-alias",
            owner_module="phospy.science.activities.results",
            deprecated="KinaseActivityResult.weighted_activity",
            replacement="KinaseActivityResult.activity_matrix",
            stability="stable",
            replacement_module="phospy.api",
            replacement_name="KinaseActivityResult",
        ),
        _record(
            identifier="activities.result.legacy_condition_statistics_table",
            kind="method-alias",
            owner_module="phospy.science.activities.results",
            deprecated=(
                "KinaseActivityResult.legacy_condition_statistics_table_dataframe()"
            ),
            replacement=(
                "statistics_table_dataframe() or statistics_table and the "
                "canonical profile_id column"
            ),
            stability="stable",
            replacement_module="phospy.api",
            replacement_name="KinaseActivityResult",
            summary=(
                "KinaseActivityResult."
                "legacy_condition_statistics_table_dataframe() is deprecated "
                "and does not establish a biological condition contract"
            ),
        ),
        _record(
            identifier="preprocessing.pipeline.stage_metadata_registry",
            kind="argument-alias",
            owner_module="phospy.science.datasets.preprocessing.pipeline",
            deprecated="stage_metadata_registry",
            replacement="stage_contract_registry",
            stability="unsupported",
            replacement_module="phospy.science.datasets.preprocessing.pipeline",
            replacement_name="PreprocessingPipeline",
            planned_removal_version="1.8.0",
            summary=(
                "PreprocessingPipeline(stage_metadata_registry=...) is "
                "deprecated because stage_metadata_registry is a legacy alias "
                "for stage_contract_registry"
            ),
        ),
        _record(
            identifier="workflows.differential.TechnicalReplicateResolver",
            kind="class-alias",
            owner_module="phospy.workflows.differential.replicates",
            deprecated="TechnicalReplicateResolver",
            replacement=(
                "TechnicalReplicateAggregationPlanner and TechnicalReplicateAggregator"
            ),
            stability="unsupported",
            replacement_module="phospy.workflows.differential.replicates",
            replacement_name="TechnicalReplicateAggregationPlanner",
            summary=(
                "TechnicalReplicateResolver is deprecated; use "
                "TechnicalReplicateAggregationPlanner and "
                "TechnicalReplicateAggregator"
            ),
        ),
    ]
    records.extend(_enrichment_alias_deprecations())
    records.extend(_activity_summary_alias_deprecations())
    return tuple(records)


def _enrichment_alias_deprecations() -> tuple[RetainedDeprecation, ...]:
    return tuple(
        _record(
            identifier=f"io.enrichment.{deprecated_name}",
            kind="function-alias",
            owner_module="phospy.io.readers.enrichment_sets",
            deprecated=deprecated_name,
            replacement=replacement_name,
            stability="unsupported",
            replacement_module="phospy.io.readers.enrichment_sets",
            replacement_name=replacement_name,
        )
        for deprecated_name, replacement_name in (
            ("load_enrichment_sets_gmt", "read_enrichment_sets_gmt"),
            ("load_enrichment_sets_table", "read_enrichment_sets_table"),
            ("load_enrichment_sets_csv", "read_enrichment_sets_csv"),
            ("load_enrichment_sets_tsv", "read_enrichment_sets_tsv"),
        )
    )


def _activity_summary_alias_deprecations() -> tuple[RetainedDeprecation, ...]:
    return tuple(
        _record(
            identifier=f"activities.method_summary.{deprecated_name}",
            kind="property-alias",
            owner_module="phospy.science.activities.method_models",
            deprecated=f"ActivityMethodSummary.{deprecated_name}",
            replacement=replacement_name,
            stability="unsupported",
            replacement_module="phospy.science.activities.method_models",
            replacement_name="ActivityMethodSummary",
            summary=(
                f"ActivityMethodSummary.{deprecated_name} is deprecated; use "
                f"{replacement_name}. The deprecated condition-named alias is "
                "only a compatibility counter name and does not define "
                "biological condition semantics."
            ),
        )
        for deprecated_name, replacement_name in (
            ("kinase_condition_pairs_evaluated", "kinase_profile_pairs_evaluated"),
            ("kinase_condition_pairs_computed", "kinase_profile_pairs_computed"),
            (
                "kinase_condition_pairs_insufficient_substrates",
                "kinase_profile_pairs_insufficient_substrates",
            ),
            (
                "kinase_condition_pairs_invalid_background_variance",
                "kinase_profile_pairs_invalid_background_variance",
            ),
            (
                "kinase_condition_pairs_no_finite_background_values",
                "kinase_profile_pairs_no_finite_background_values",
            ),
            (
                "kinase_condition_pairs_no_finite_substrate_values",
                "kinase_profile_pairs_no_finite_substrate_values",
            ),
        )
    )


def _api_compatibility_deprecations() -> tuple[RetainedDeprecation, ...]:
    entries: dict[tuple[str, str], RetainedDeprecation] = {}

    for name in ADVANCED_PUBLIC_API:
        _register_api_compatibility(
            entries,
            old_module="phospy.api",
            name=name,
            owner_module="phospy.advanced",
            replacement_module="phospy.advanced",
            stability="advanced",
        )

    for name in REQUEST_COMPATIBILITY_ADVANCED_API:
        _register_api_compatibility(
            entries,
            old_module="phospy.api.requests",
            name=name,
            owner_module="phospy.advanced",
            replacement_module="phospy.advanced",
            stability="advanced",
        )

    advanced_config_routes = _advanced_config_routes()
    for name in ADVANCED_CONFIG_API:
        for old_module in advanced_config_routes.get(name, ()):
            _register_api_compatibility(
                entries,
                old_module=old_module,
                name=name,
                owner_module="phospy.advanced.configs",
                replacement_module="phospy.advanced.configs",
                stability="advanced",
            )

    for name in ADVANCED_RESULT_API:
        _register_api_compatibility(
            entries,
            old_module="phospy.api.results",
            name=name,
            owner_module="phospy.advanced.results",
            replacement_module="phospy.advanced.results",
            stability="advanced",
        )

    _register_contract_request_compatibility(entries)
    _register_contract_config_compatibility(entries)
    return tuple(entries.values())


def _register_contract_request_compatibility(
    entries: dict[tuple[str, str], RetainedDeprecation],
) -> None:
    owner_module = "phospy.contracts.requests"
    owner = importlib.import_module(owner_module)
    for name in getattr(owner, "__all__", ()):
        if name in STABLE_REQUEST_API:
            continue
        _register_api_compatibility(
            entries,
            old_module="phospy.api.requests",
            name=name,
            owner_module=owner_module,
            replacement_module=owner_module,
            stability="unsupported",
        )


def _register_contract_config_compatibility(
    entries: dict[tuple[str, str], RetainedDeprecation],
) -> None:
    for old_module in COMPATIBILITY_CONFIG_MODULES:
        owner_module = _config_owner_module(old_module)
        owner = importlib.import_module(owner_module)
        for name in getattr(owner, "__all__", ()):
            if name in STABLE_CONFIG_API or name in ADVANCED_CONFIG_API:
                continue
            _register_api_compatibility(
                entries,
                old_module=old_module,
                name=name,
                owner_module=owner_module,
                replacement_module=owner_module,
                stability="unsupported",
            )


def _advanced_config_routes() -> dict[str, tuple[str, ...]]:
    routes: dict[str, list[str]] = {name: [] for name in ADVANCED_CONFIG_API}
    for old_module in COMPATIBILITY_CONFIG_MODULES:
        owner_module = _config_owner_module(old_module)
        owner = importlib.import_module(owner_module)
        owner_exports = set(getattr(owner, "__all__", ()))
        for name in ADVANCED_CONFIG_API:
            if old_module == "phospy.api.configs" or name in owner_exports:
                routes[name].append(old_module)
    for name, route_overrides in CONFIG_COMPATIBILITY_ADVANCED_ROUTE_OVERRIDES.items():
        routes.setdefault(name, [])
        for old_module in route_overrides:
            if old_module not in routes[name]:
                routes[name].append(old_module)
    return {name: tuple(route_names) for name, route_names in routes.items()}


def _config_owner_module(old_module: str) -> str:
    return old_module.replace(
        "phospy.api.configs",
        "phospy.contracts.configs",
        1,
    )


def _register_api_compatibility(
    entries: dict[tuple[str, str], RetainedDeprecation],
    *,
    old_module: str,
    name: str,
    owner_module: str,
    replacement_module: str,
    stability: DeprecationStability,
) -> None:
    entries.setdefault(
        (old_module, name),
        _record(
            identifier=api_compatibility_deprecation_id(
                old_module=old_module,
                name=name,
            ),
            kind="import-route",
            owner_module=owner_module,
            deprecated=f"{old_module}.{name}",
            replacement=f"`from {replacement_module} import {name}`",
            introduced_version=API_COMPATIBILITY_INTRODUCED_VERSION,
            planned_removal_version=API_COMPATIBILITY_PLANNED_REMOVAL_VERSION,
            stability=stability,
            replacement_module=replacement_module,
            replacement_name=name,
            summary=(
                f"{old_module}.{name} is deprecated as {_api_stability_note(stability)}"
            ),
            deprecated_module=old_module,
            deprecated_name=name,
        ),
    )


def _api_stability_note(stability: DeprecationStability) -> str:
    if stability == "advanced":
        return "an advanced API route"
    return "an unsupported compatibility route"


def _record(
    *,
    identifier: str,
    kind: DeprecationKind,
    owner_module: str,
    deprecated: str,
    replacement: str,
    stability: DeprecationStability,
    replacement_module: str,
    replacement_name: str,
    introduced_version: str = _DEFAULT_INTRODUCED_VERSION,
    planned_removal_version: str = _DEFAULT_PLANNED_REMOVAL_VERSION,
    summary: str = "",
    deprecated_module: str | None = None,
    deprecated_name: str | None = None,
) -> RetainedDeprecation:
    return RetainedDeprecation(
        identifier=identifier,
        kind=kind,
        owner_module=owner_module,
        deprecated=deprecated,
        replacement=replacement,
        introduced_version=introduced_version,
        planned_removal_version=planned_removal_version,
        stability=stability,
        replacement_module=replacement_module,
        replacement_name=replacement_name,
        summary=summary,
        deprecated_module=deprecated_module,
        deprecated_name=deprecated_name,
    )


__all__ = [
    "PhosPyDeprecationWarning",
    "RetainedDeprecation",
    "api_compatibility_deprecation_id",
    "compatibility_deprecation_record",
    "compatibility_deprecation_records",
    "deprecation_record",
    "retained_deprecations",
    "warn_deprecated",
]
