"""Trusted-table assertion checks for dataset reconstruction."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, cast

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.provenance.models import RunProvenance, TrustedDatasetConstructionAssertions

if TYPE_CHECKING:
    from phospy.science.datasets.construction.analysis_ready import (
        AnalysisReadyPhosphoDataset,
    )


def _resolve_trusted_construction_assertions(
    *,
    trusted_construction_assertions: TrustedDatasetConstructionAssertions | None,
    provenance: RunProvenance | None,
    assume_owned: bool,
) -> TrustedDatasetConstructionAssertions | None:
    del provenance, assume_owned
    return trusted_construction_assertions


def _require_complete_from_trusted_assertions(
    *,
    dataset: AnalysisReadyPhosphoDataset,
) -> None:
    assertions = dataset.trusted_construction_assertions
    required_message = (
        "AnalysisReadyPhosphoDataset.from_trusted_tables requires "
        "trusted_construction_assertions with typed evidence or an explicit "
        "waiver for identity, intensity scale, quantitative meaning, aligned "
        "structure, localisation, sequence, and reference context"
    )
    if assertions is None:
        raise DatasetValidationError(required_message)
    if not assertions.assertion_metadata_provided:
        raise DatasetValidationError(required_message)
    if not assertions.all_required_assertions_present:
        raise DatasetValidationError(
            required_message + "; missing: " + ", ".join(assertions.missing_assertions)
        )
    _require_assertions_linked_to_provenance(dataset=dataset)


def _require_assertions_linked_to_provenance(
    *,
    dataset: AnalysisReadyPhosphoDataset,
) -> None:
    assertions = dataset.trusted_construction_assertions
    provenance = dataset.provenance
    if assertions is None or provenance is None:
        raise DatasetValidationError(
            "dataset.trusted_construction_assertions must be linked to "
            "dataset.provenance"
        )
    construction = provenance.workflow_parameters.get("construction")
    if not isinstance(construction, Mapping):
        raise DatasetValidationError(
            "dataset.provenance.workflow_parameters['construction'] must record "
            "trusted construction assertion provenance"
        )
    construction_payload = cast(Mapping[str, object], construction)
    raw_assertion_payload = construction_payload.get("trusted_construction_assertions")
    if not isinstance(raw_assertion_payload, Mapping):
        raise DatasetValidationError(
            "dataset.provenance.workflow_parameters['construction'] must record "
            "trusted_construction_assertions as a self-contained assertion payload"
        )
    try:
        provenance_assertions = TrustedDatasetConstructionAssertions.from_payload(
            cast(Mapping[str, object], raw_assertion_payload),
            field_name=(
                "dataset.provenance.workflow_parameters['construction']"
                "['trusted_construction_assertions']"
            ),
        )
    except PhosPyInputError as exc:
        raise DatasetValidationError(
            "dataset.provenance.workflow_parameters['construction']"
            "['trusted_construction_assertions'] must be a current-schema trusted "
            f"construction assertion payload; {exc}"
        ) from exc
    if provenance_assertions.to_payload() != assertions.to_payload():
        raise DatasetValidationError(
            "dataset.trusted_construction_assertions must be provenance-linked; "
            "trusted_construction_assertions payload does not match"
        )
    observed = construction_payload.get("trusted_construction_assertion_fingerprint")
    if observed != provenance_assertions.assertion_fingerprint:
        raise DatasetValidationError(
            "dataset.trusted_construction_assertions must be provenance-linked; "
            "trusted_construction_assertion_fingerprint does not match "
            "trusted_construction_assertions payload"
        )
    if observed != assertions.assertion_fingerprint:
        raise DatasetValidationError(
            "dataset.trusted_construction_assertions must be provenance-linked; "
            "trusted_construction_assertion_fingerprint does not match supplied "
            "trusted_construction_assertions"
        )
