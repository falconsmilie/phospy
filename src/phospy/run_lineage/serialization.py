"""Run-level provenance payload serialization."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.evidence.serialization.historical_peptide import (
    normalise_historical_workflow_parameters,
)
from phospy.preprocessing.provenance.serialization import (
    stage_from_payload,
    stage_to_payload,
)
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
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.references.provenance.context import ReferenceContext
from phospy.references.provenance.serialization import (
    reference_from_payload,
    reference_to_payload,
)
from phospy.run_lineage.models import RunProvenance


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
    workflow_parameters = normalise_historical_workflow_parameters(workflow_parameters)
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
