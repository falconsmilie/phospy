"""Run-level provenance payload serialization."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from phospy.provenance.models import RunProvenance
from phospy.provenance.reference_context import ReferenceContext
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.provenance.serialization._payload import (
    optional_int,
    optional_mapping,
    optional_str,
    require_mapping,
    require_sequence,
    to_json_safe,
    to_json_value,
)
from phospy.provenance.serialization.environment import (
    environment_from_payload,
    environment_to_payload,
)
from phospy.provenance.serialization.references import (
    reference_from_payload,
    reference_to_payload,
)
from phospy.provenance.serialization.stages import (
    stage_from_payload,
    stage_to_payload,
)
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)

_LEGACY_PEPTIDE_RESOLUTION_AGGREGATION_POLICY = "mapping_weighted_mean"

_LEGACY_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_POLICY = "sum_to_one_per_peptide_row"

_LEGACY_DUPLICATE_PEPTIDE_POLICY = "retain_all_peptide_rows_as_independent_observations"

_LEGACY_MIXED_AMBIGUITY_POLICY = (
    "mixed_ambiguous_and_unambiguous_rows_share_same_weighted_mean_aggregation"
)

_CURRENT_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY = (
    "explicit_mapping_weight_when_supplied_else_equal_fraction_per_resolved_site"
)

_CURRENT_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY = (
    "sum_to_one_per_peptide_evidence_row"
)

_CURRENT_PEPTIDE_SIGNAL_ALLOCATION_POLICY = (
    "multiply_peptide_signal_by_mapping_fraction"
)

_CURRENT_PEPTIDE_SITE_SUMMARISATION_POLICY = "arithmetic_mean_of_allocated_signals"

_CURRENT_PEPTIDE_DUPLICATE_EVIDENCE_POLICY = (
    "retain_duplicate_peptide_evidence_rows_as_separate_observations"
)

_CURRENT_PEPTIDE_MIXED_AMBIGUITY_POLICY = (
    "combine_ambiguous_and_unambiguous_allocated_signals_in_site_mean"
)

_CURRENT_PEPTIDE_LOCALISATION_AGGREGATION_POLICY = (
    "arithmetic_mean_of_finite_reported_localisation_values"
)

_CURRENT_PEPTIDE_LEGACY_AGGREGATION_ALIAS = (
    "legacy_alias_for_arithmetic_mean_of_allocated_signals"
)


def to_payload(provenance: RunProvenance) -> dict[str, object]:
    """Serialize run provenance to a JSON-safe payload."""

    return {
        "environment": environment_to_payload(provenance.environment),
        "input_tables": [
            table_fingerprint_to_payload(item) for item in provenance.input_tables
        ],
        "preprocessing_stages": [
            stage_to_payload(item) for item in provenance.preprocessing_stages
        ],
        "reference": (
            None
            if provenance.reference is None
            else reference_to_payload(provenance.reference)
        ),
        "reference_context": (
            None
            if provenance.reference_context is None
            else provenance.reference_context.to_payload()
        ),
        "workflow_name": provenance.workflow_name,
        "workflow_parameters": to_json_safe(provenance.workflow_parameters),
        "random_state": provenance.random_state,
        "random_seed_policy": provenance.random_seed_policy,
        "output_tables": [
            table_fingerprint_to_payload(item) for item in provenance.output_tables
        ],
        "scientific_policies": [
            item.to_payload() for item in provenance.scientific_policies
        ],
    }


def from_payload(payload: Mapping[str, object]) -> RunProvenance:
    """Deserialize run provenance from a decoded payload."""

    payload = require_mapping(payload, field_name="provenance")
    environment_payload = require_mapping(
        payload.get("environment"),
        field_name="provenance.environment",
    )
    input_tables_payload = require_sequence(
        payload.get("input_tables"),
        field_name="provenance.input_tables",
    )
    stages_payload = require_sequence(
        payload.get("preprocessing_stages"),
        field_name="provenance.preprocessing_stages",
    )
    output_tables_payload = require_sequence(
        payload.get("output_tables"),
        field_name="provenance.output_tables",
    )
    scientific_policies_payload = require_sequence(
        payload.get("scientific_policies", []),
        field_name="provenance.scientific_policies",
    )
    reference_raw = payload.get("reference")
    if reference_raw is None:
        reference = None
    else:
        reference = reference_from_payload(
            require_mapping(reference_raw, field_name="provenance.reference")
        )
    reference_context_payload = optional_mapping(
        payload.get("reference_context"),
        field_name="provenance.reference_context",
    )
    workflow_parameters = require_mapping(
        payload.get("workflow_parameters"),
        field_name="provenance.workflow_parameters",
    )
    workflow_parameters = _normalise_workflow_parameters(workflow_parameters)
    return RunProvenance(
        environment=environment_from_payload(environment_payload),
        input_tables=tuple(
            table_fingerprint_from_payload(
                require_mapping(
                    item,
                    field_name=f"provenance.input_tables[{position}]",
                )
            )
            for position, item in enumerate(input_tables_payload)
        ),
        preprocessing_stages=tuple(
            stage_from_payload(
                require_mapping(
                    item,
                    field_name=f"provenance.preprocessing_stages[{position}]",
                )
            )
            for position, item in enumerate(stages_payload)
        ),
        reference=reference,
        workflow_name=optional_str(
            payload.get("workflow_name"),
            field_name="provenance.workflow_name",
        ),
        workflow_parameters={
            key: to_json_value(value) for key, value in workflow_parameters.items()
        },
        random_state=optional_int(
            payload.get("random_state"),
            field_name="provenance.random_state",
        ),
        random_seed_policy=optional_str(
            payload.get("random_seed_policy"),
            field_name="provenance.random_seed_policy",
        ),
        output_tables=tuple(
            table_fingerprint_from_payload(
                require_mapping(
                    item,
                    field_name=f"provenance.output_tables[{position}]",
                )
            )
            for position, item in enumerate(output_tables_payload)
        ),
        scientific_policies=tuple(
            ScientificPolicyRecord.from_payload(
                require_mapping(
                    item,
                    field_name=f"provenance.scientific_policies[{position}]",
                )
            )
            for position, item in enumerate(scientific_policies_payload)
        ),
        reference_context=None
        if reference_context_payload is None
        else ReferenceContext.from_payload(reference_context_payload),
    )


def _normalise_workflow_parameters(
    workflow_parameters: Mapping[str, object],
) -> dict[str, object]:
    normalized = dict(workflow_parameters)
    peptide_resolution = normalized.get("peptide_evidence_resolution")
    if isinstance(peptide_resolution, Mapping):
        peptide_resolution_payload = cast(Mapping[str, object], peptide_resolution)
        normalized["peptide_evidence_resolution"] = (
            _normalise_peptide_evidence_resolution_workflow_parameter(
                peptide_resolution_payload
            )
        )
    return normalized


def _normalise_peptide_evidence_resolution_workflow_parameter(
    payload: Mapping[str, object],
) -> dict[str, object]:
    normalized = dict(payload)
    normalized["mapping_weight_source_policy"] = str(
        normalized.get("mapping_weight_source_policy")
        or _CURRENT_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY
    )

    legacy_normalization = normalized.get("mapping_weight_normalisation")
    if "mapping_weight_normalization_policy" not in normalized:
        normalized["mapping_weight_normalization_policy"] = (
            _CURRENT_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY
            if legacy_normalization
            in (
                None,
                _LEGACY_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_POLICY,
                _CURRENT_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY,
            )
            else str(legacy_normalization)
        )
    normalized["mapping_weight_normalisation"] = str(
        normalized.get("mapping_weight_normalisation")
        or normalized["mapping_weight_normalization_policy"]
    )
    if (
        normalized["mapping_weight_normalisation"]
        == _LEGACY_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_POLICY
    ):
        normalized["mapping_weight_normalisation"] = str(
            normalized["mapping_weight_normalization_policy"]
        )

    normalized["signal_allocation_policy"] = str(
        normalized.get("signal_allocation_policy")
        or _CURRENT_PEPTIDE_SIGNAL_ALLOCATION_POLICY
    )
    normalized["site_summarisation_policy"] = str(
        normalized.get("site_summarisation_policy")
        or _CURRENT_PEPTIDE_SITE_SUMMARISATION_POLICY
    )

    legacy_duplicate_policy = normalized.get("duplicate_peptide_policy")
    if "duplicate_evidence_policy" not in normalized:
        normalized["duplicate_evidence_policy"] = (
            _CURRENT_PEPTIDE_DUPLICATE_EVIDENCE_POLICY
            if legacy_duplicate_policy
            in (
                None,
                _LEGACY_DUPLICATE_PEPTIDE_POLICY,
                _CURRENT_PEPTIDE_DUPLICATE_EVIDENCE_POLICY,
            )
            else str(legacy_duplicate_policy)
        )
    normalized["duplicate_peptide_policy"] = str(
        normalized.get("duplicate_peptide_policy")
        or normalized["duplicate_evidence_policy"]
    )
    if normalized["duplicate_peptide_policy"] == _LEGACY_DUPLICATE_PEPTIDE_POLICY:
        normalized["duplicate_peptide_policy"] = str(
            normalized["duplicate_evidence_policy"]
        )

    mixed_ambiguity_policy = normalized.get("mixed_ambiguity_policy")
    normalized["mixed_ambiguity_policy"] = (
        _CURRENT_PEPTIDE_MIXED_AMBIGUITY_POLICY
        if mixed_ambiguity_policy
        in (
            None,
            _LEGACY_MIXED_AMBIGUITY_POLICY,
            _CURRENT_PEPTIDE_MIXED_AMBIGUITY_POLICY,
        )
        else str(mixed_ambiguity_policy)
    )
    normalized["localisation_aggregation_policy"] = str(
        normalized.get("localisation_aggregation_policy")
        or _CURRENT_PEPTIDE_LOCALISATION_AGGREGATION_POLICY
    )

    aggregation_policy = normalized.get("aggregation_policy")
    normalized["aggregation_policy"] = (
        _CURRENT_PEPTIDE_LEGACY_AGGREGATION_ALIAS
        if aggregation_policy
        in (
            None,
            _LEGACY_PEPTIDE_RESOLUTION_AGGREGATION_POLICY,
            _CURRENT_PEPTIDE_LEGACY_AGGREGATION_ALIAS,
        )
        else str(aggregation_policy)
    )
    normalized["aggregation_formula"] = str(
        normalized.get("aggregation_formula")
        or (
            "a[p,s,j] = mapping_fraction[p,s] * peptide_signal[p,j]; "
            "site_signal[s,j] = arithmetic_mean(a[p,s,j] for retained peptide "
            "rows p mapped to s)"
        )
    )
    return normalized
