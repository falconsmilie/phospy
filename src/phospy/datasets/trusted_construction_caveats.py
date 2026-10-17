"""Dataset trusted-construction result caveat assembly."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

from phospy.datasets.direct_construction import (
    DIRECT_CONSTRUCTION_SOURCE,
    DIRECT_CONSTRUCTION_WORKFLOW_NAME,
)
from phospy.datasets.trusted_assertions import (
    TrustedDatasetConstructionAssertions,
)
from phospy.provenance.result_caveats import ResultCaveat, ResultCaveatSeverity

if TYPE_CHECKING:
    from phospy.datasets.models import AnalysisReadyPhosphoDataset

_TRUSTED_CONSTRUCTION_ASSERTION_FIELDS = (
    "identity_user_asserted",
    "quantitative_meaning_user_asserted",
    "localisation_user_asserted",
    "sequence_user_asserted",
    "reference_context_user_asserted",
)


def build_direct_trusted_dataset_construction_caveat(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    code: str,
    workflow_scope: str,
    workflow_label: str,
) -> ResultCaveat | None:
    """Build a workflow caveat for direct trusted dataset construction."""

    details = direct_trusted_dataset_construction_details(
        dataset=dataset,
        workflow_scope=workflow_scope,
    )
    if details is None:
        return None
    severity = direct_trusted_dataset_caveat_severity(details)
    message = (
        "Input dataset was directly constructed as trusted analysis-ready state; "
        f"{workflow_label} did not rerun dataset-building validation."
    )
    if severity == "warning":
        message = (
            message
            + " Trusted construction assertion metadata is absent or incomplete."
        )
    return ResultCaveat(
        code=code,
        severity=severity,
        message=message,
        details=details,
    )


def direct_trusted_dataset_construction_details(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    workflow_scope: str,
) -> dict[str, object] | None:
    """Return direct trusted-construction details for result caveats."""

    provenance = dataset.provenance
    construction = (
        {}
        if provenance is None
        else _construction_payload(provenance.workflow_parameters)
    )
    assertions = dataset.trusted_construction_assertions
    source = construction.get("source")
    if (
        assertions is None
        and provenance is not None
        and provenance.workflow_name != DIRECT_CONSTRUCTION_WORKFLOW_NAME
        and source != DIRECT_CONSTRUCTION_SOURCE
    ):
        return None

    assertion_payload = _trusted_assertion_payload(
        assertions=assertions,
        construction=construction,
    )
    missing_assertions = _missing_trusted_assertions(assertion_payload)
    metadata_provided = assertion_payload.get("assertion_metadata_provided")
    metadata_provided_bool = (
        bool(metadata_provided)
        if isinstance(
            metadata_provided,
            bool,
        )
        else False
    )
    details: dict[str, object] = {
        "workflow_scope": workflow_scope,
        "dataset_workflow_name": None
        if provenance is None
        else provenance.workflow_name,
        "construction_source": (
            DIRECT_CONSTRUCTION_SOURCE if source is None else str(source)
        ),
        "trusted_assertion_metadata_provided": metadata_provided_bool,
        "trusted_construction_assertions": assertion_payload,
        "missing_trusted_assertions": list(missing_assertions),
    }
    builder_used = construction.get("builder_used")
    if isinstance(builder_used, bool):
        details["builder_used"] = builder_used
    warning = construction.get("warning")
    if isinstance(warning, str) and warning.strip():
        details["construction_warning"] = warning.strip()
    assertion_warning = construction.get("assertion_warning")
    if isinstance(assertion_warning, str) and assertion_warning.strip():
        details["assertion_warning"] = assertion_warning.strip()
    return details


def direct_trusted_dataset_caveat_severity(
    details: Mapping[str, object],
) -> ResultCaveatSeverity:
    """Return warning when trusted construction assertion metadata is incomplete."""

    missing = details.get("missing_trusted_assertions")
    metadata_provided = details.get("trusted_assertion_metadata_provided")
    if metadata_provided is not True:
        return "warning"
    if (
        isinstance(missing, Sequence)
        and not isinstance(missing, (str, bytes, bytearray))
        and missing
    ):
        return "warning"
    return "info"


def _construction_payload(parameters: Mapping[str, object]) -> Mapping[str, object]:
    construction = parameters.get("construction")
    if isinstance(construction, Mapping):
        return construction
    return {}


def _trusted_assertion_payload(
    *,
    assertions: TrustedDatasetConstructionAssertions | None,
    construction: Mapping[str, object],
) -> dict[str, object]:
    if assertions is not None:
        return assertions.to_payload()
    raw_payload = construction.get("trusted_construction_assertions")
    if isinstance(raw_payload, Mapping):
        return {str(key): value for key, value in raw_payload.items()}
    return TrustedDatasetConstructionAssertions.missing().to_payload()


def _missing_trusted_assertions(
    assertion_payload: Mapping[str, object],
) -> tuple[str, ...]:
    raw_missing = assertion_payload.get("missing_assertions")
    if isinstance(raw_missing, Sequence) and not isinstance(
        raw_missing,
        (str, bytes, bytearray),
    ):
        return tuple(str(item) for item in raw_missing)
    return tuple(
        field_name
        for field_name in _TRUSTED_CONSTRUCTION_ASSERTION_FIELDS
        if assertion_payload.get(field_name) is not True
    )


__all__ = [
    "build_direct_trusted_dataset_construction_caveat",
    "direct_trusted_dataset_caveat_severity",
    "direct_trusted_dataset_construction_details",
]
