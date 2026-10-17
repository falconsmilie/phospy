"""Passive caller inputs and validation diagnostics for SPS discovery."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import ReferenceValidationError
from phospy.preprocessing.control_sites._discovery_payload import (
    require_int_at_least,
    require_non_empty_string,
    require_non_negative_int,
    require_payload_int,
    require_payload_mapping,
    require_payload_string,
    require_positive_int,
)
from phospy.provenance.hashing import (
    DEFAULT_STABLE_JSON_HASH_ALGORITHM,
    fingerprint_table_normalized_axes,
)
from phospy.provenance.immutability import freeze_json_mapping, thaw_json_mapping
from phospy.provenance.models import (
    TableFingerprint,
)
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.quantitative._authority import (
    sps_reference_establishment_authority,
    sps_reference_quantitative_meaning_transition_authority,
)
from phospy.quantitative.declared_scale_diagnostics import (
    _run_declared_scale_sanity_diagnostics,  # pyright: ignore[reportPrivateUsage] - The authorised shared helpers retain their historical private names.
)
from phospy.quantitative.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentProvenance,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningTransitionProvenance,
    establish_intensity_scale_state,
)
from phospy.taxonomy.organisms import Organism, normalize_optional_organism

SPS_DISCOVERY_SELECTION_METHOD_CONSENSUS_STABILITY = "consensus_stability"
SPS_DISCOVERY_ALGORITHM_ID = "phospy_sps_consensus_stability"
SPS_DISCOVERY_ALGORITHM_VERSION = "2.0.0"
SPS_DISCOVERY_LEGACY_ALGORITHM_VERSION = "1.0.0"
SPS_DISCOVERY_CROSS_COUNT_ORDERING_POLICY = (
    "contributing_dataset_count_descending_then_consensus_stability_descending"
)
SPS_DISCOVERY_LEGACY_CROSS_COUNT_ORDERING_POLICY = (
    "consensus_stability_descending_ignoring_contribution_count"
)
SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING = "site_key_ascending"
SPS_DISCOVERY_CONTROL_SOURCE_TYPE = "sps_discovery"
SPS_DISCOVERY_IDENTITY_SCHEMA_V1 = "phospy-sps-discovery-result-v1"
SPS_DISCOVERY_IDENTITY_SCHEMA = "phospy-sps-discovery-result-v2"
SPS_DISCOVERY_IDENTITY_PREFIX = f"{DEFAULT_STABLE_JSON_HASH_ALGORITHM}:"
SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER = "sps_reference_matrix_fingerprint"
SPS_REFERENCE_CENTERING_ASSERTION_OPERATION = (
    "phospy.sps_reference.external_reference_centering_assertion"
)
SPS_REFERENCE_DECLARED_SCALE_PRODUCER_ID = (
    "phospy.science.transformations.transformers.identity.IdentityTransformer"
)
SPS_REFERENCE_DATASET_OVERRIDE_SCOPE_WARNING = (
    "The suspicious declared input intensity scale was accepted for dataset "
    "construction only. Downstream workflows may independently validate this "
    "declaration and may require their own explicit override. Differential "
    "analysis requires "
    "DifferentialAnalysisConfig.allow_suspicious_declared_input_scale=True."
)

SpsDiscoverySelectionMethod = Literal["consensus_stability"]
SpsDiscoveryTieHandling = Literal["site_key_ascending"]


@dataclass(frozen=True, slots=True)
class SpsDiscoveryConfig:
    """Immutable scientific configuration for later SPS selection.

    ``minimum_datasets_per_site`` is evaluated independently for every site;
    reference datasets do not need to share samples, group labels, or a full
    common site intersection.
    """

    top_n: int = 100
    minimum_reference_datasets: int = 2
    minimum_datasets_per_site: int = 2
    minimum_shared_sites: int = 1
    tie_handling: SpsDiscoveryTieHandling = (
        SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING
    )
    selection_method: SpsDiscoverySelectionMethod = (
        SPS_DISCOVERY_SELECTION_METHOD_CONSENSUS_STABILITY
    )

    def __post_init__(self) -> None:
        require_positive_int(self.top_n, field_name="sps_discovery.top_n")
        require_int_at_least(
            self.minimum_reference_datasets,
            minimum=2,
            field_name="sps_discovery.minimum_reference_datasets",
        )
        require_int_at_least(
            self.minimum_datasets_per_site,
            minimum=2,
            field_name="sps_discovery.minimum_datasets_per_site",
        )
        require_positive_int(
            self.minimum_shared_sites,
            field_name="sps_discovery.minimum_shared_sites",
        )
        if self.tie_handling != SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING:
            raise PhosPyInputError(
                "sps_discovery.tie_handling must be 'site_key_ascending'"
            )
        if self.selection_method != SPS_DISCOVERY_SELECTION_METHOD_CONSENSUS_STABILITY:
            raise PhosPyInputError(
                "sps_discovery.selection_method must be 'consensus_stability'"
            )

    def to_payload(self) -> dict[str, object]:
        """Return a deterministic JSON-compatible configuration payload."""

        return {
            "top_n": self.top_n,
            "minimum_reference_datasets": self.minimum_reference_datasets,
            "minimum_datasets_per_site": self.minimum_datasets_per_site,
            "minimum_shared_sites": self.minimum_shared_sites,
            "tie_handling": self.tie_handling,
            "selection_method": self.selection_method,
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsDiscoveryConfig:
        """Restore a validated SPS configuration from its public payload."""

        resolved = require_payload_mapping(payload, field_name="sps_discovery_config")
        return cls(
            top_n=require_payload_int(
                resolved.get("top_n"), field_name="sps_discovery_config.top_n"
            ),
            minimum_reference_datasets=require_payload_int(
                resolved.get("minimum_reference_datasets"),
                field_name="sps_discovery_config.minimum_reference_datasets",
            ),
            minimum_datasets_per_site=require_payload_int(
                resolved.get("minimum_datasets_per_site"),
                field_name="sps_discovery_config.minimum_datasets_per_site",
            ),
            minimum_shared_sites=require_payload_int(
                resolved.get("minimum_shared_sites"),
                field_name="sps_discovery_config.minimum_shared_sites",
            ),
            tie_handling=cast(
                SpsDiscoveryTieHandling,
                require_payload_string(
                    resolved.get("tie_handling"),
                    field_name="sps_discovery_config.tie_handling",
                ),
            ),
            selection_method=cast(
                SpsDiscoverySelectionMethod,
                require_payload_string(
                    resolved.get("selection_method"),
                    field_name="sps_discovery_config.selection_method",
                ),
            ),
        )


@dataclass(frozen=True, slots=True, eq=False, init=False)
class SpsReferenceDataset:
    """One phosphoproteomics reference matrix and its sample groups.

    The matrix is defensively copied. Rows must ultimately be governed by its
    encoded ``site_key`` index; the private request validator owns that check.
    ``condition_by_sample`` is frozen recursively so caller mutation cannot
    alter a validated request. ``intensity_scale_state`` carries the governed
    scale, quantitative meaning, and establishment evidence; SPS validation
    requires it to describe pre-established condition-relative log2 values.
    Successful discovery additionally requires a typed organism, explicit
    biological baseline and reference context, and complete source identity.
    """

    dataset_id: str
    _owned_intensities: object = field(repr=False)
    condition_by_sample: Mapping[str, str]
    intensity_scale_state: IntensityScaleState
    source_name: str | None = None
    source_version: str | None = None
    source_uri: str | None = None
    organism: Organism | None = None
    baseline_context: str | None = None
    reference_context: str | None = None

    def __init__(
        self,
        dataset_id: str,
        intensities: pd.DataFrame,
        condition_by_sample: Mapping[str, str],
        intensity_scale_state: IntensityScaleState,
        source_name: str | None = None,
        source_version: str | None = None,
        source_uri: str | None = None,
        *,
        organism: Organism | str | None = None,
        baseline_context: str | None = None,
        reference_context: str | None = None,
    ) -> None:
        object.__setattr__(self, "dataset_id", dataset_id)
        object.__setattr__(
            self,
            "_owned_intensities",
            (
                intensities.copy(deep=True)
                if isinstance(cast(object, intensities), pd.DataFrame)
                else intensities
            ),
        )
        object.__setattr__(self, "condition_by_sample", condition_by_sample)
        object.__setattr__(self, "intensity_scale_state", intensity_scale_state)
        object.__setattr__(self, "organism", organism)
        object.__setattr__(self, "baseline_context", baseline_context)
        object.__setattr__(self, "reference_context", reference_context)
        object.__setattr__(self, "source_name", source_name)
        object.__setattr__(self, "source_version", source_version)
        object.__setattr__(self, "source_uri", source_uri)
        self.__post_init__()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dataset_id",
            require_non_empty_string(
                self.dataset_id, field_name="sps_reference_dataset.dataset_id"
            ),
        )
        if isinstance(cast(object, self.condition_by_sample), Mapping):
            try:
                frozen_conditions = freeze_json_mapping(
                    self.condition_by_sample,
                    field_name=(
                        f"sps_reference_dataset[{self.dataset_id!r}]."
                        "condition_by_sample"
                    ),
                )
            except PhosPyInputError as exc:
                validation = SpsDiscoveryValidationResult(
                    issues=(
                        SpsDiscoveryValidationIssue(
                            code="invalid_condition_metadata",
                            message=str(exc),
                            dataset_id=self.dataset_id,
                            field_name="condition_by_sample",
                        ),
                    ),
                    reference_dataset_count=1,
                    potential_overlap_site_count=0,
                )
                raise SpsDiscoveryValidationError(validation) from exc
            object.__setattr__(self, "condition_by_sample", frozen_conditions)
        object.__setattr__(
            self,
            "organism",
            normalize_optional_organism(
                self.organism,
                field_name=(f"sps_reference_dataset[{self.dataset_id!r}].organism"),
                error_type=PhosPyInputError,
            ),
        )
        for field_name in (
            "baseline_context",
            "reference_context",
            "source_name",
            "source_version",
            "source_uri",
        ):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(
                    self,
                    field_name,
                    require_non_empty_string(
                        value,
                        field_name=(
                            f"sps_reference_dataset[{self.dataset_id!r}].{field_name}"
                        ),
                    ),
                )

    @classmethod
    def from_condition_relative_log2(
        cls,
        dataset_id: str,
        intensities: pd.DataFrame,
        condition_by_sample: Mapping[str, str],
        *,
        log2_scale_established_by: str,
        baseline_centering_established_by: str,
        organism: Organism | str | None = None,
        baseline_context: str | None = None,
        reference_context: str | None = None,
        source_name: str | None = None,
        source_version: str | None = None,
        source_uri: str | None = None,
    ) -> SpsReferenceDataset:
        """Assert pre-established SPS input semantics without changing values.

        This supported advanced boundary records a caller assertion that the
        supplied matrix is already on the log2 scale and was already centred on
        the appropriate reference/control baseline. It does not choose a
        baseline, centre values, transform values, or infer semantics from the
        numeric matrix. Both evidence-source identifiers are required, and the
        resulting typed evidence is bound to the exact submitted matrix.
        Organism, biological contexts, and source identity remain caller
        assertions and are validated by ``SpsDiscoveryWorkflow``.
        """

        resolved_dataset_id = require_non_empty_string(
            dataset_id,
            field_name="sps_reference_dataset.dataset_id",
        )
        scale_source = require_non_empty_string(
            log2_scale_established_by,
            field_name=("sps_reference_dataset.log2_scale_established_by"),
        )
        baseline_source = require_non_empty_string(
            baseline_centering_established_by,
            field_name=("sps_reference_dataset.baseline_centering_established_by"),
        )
        if not isinstance(cast(object, intensities), pd.DataFrame):
            raise PhosPyInputError(
                "sps_reference_dataset.intensities must be a pandas DataFrame"
            )
        owned_intensities = intensities.copy(deep=True)
        matrix_fingerprint = sps_reference_intensity_fingerprint(
            dataset_id=resolved_dataset_id,
            intensities=owned_intensities,
        )
        fingerprint_payload = table_fingerprint_to_payload(matrix_fingerprint)
        declared_state = IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(),
            quantity=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )
        diagnostic_warnings = _run_declared_scale_sanity_diagnostics(
            declared_scale_kind=IntensityScaleKind.LOG2,
            phospho=owned_intensities,
            total=None,
        )
        if diagnostic_warnings:
            diagnostic_warnings = (
                *diagnostic_warnings,
                SPS_REFERENCE_DATASET_OVERRIDE_SCOPE_WARNING,
            )
        established_state = establish_intensity_scale_state(
            declared_state,
            established_via=SPS_REFERENCE_DECLARED_SCALE_PRODUCER_ID,
            establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
            evidence_level=IntensityScaleEvidenceLevel.DECLARED_BY_USER,
            transformer_name=None,
            input_declaration_source=scale_source,
            parameters={
                SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER: fingerprint_payload,
            },
            trace_id=None,
            diagnostic_warnings=diagnostic_warnings,
            _authority=sps_reference_establishment_authority(),
        )
        meaning_provenance = QuantitativeMeaningTransitionProvenance(
            source_quantity=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            target_quantity=QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
            operation_id=SPS_REFERENCE_CENTERING_ASSERTION_OPERATION,
            producer_id=baseline_source,
            evidence_mode=QuantitativeMeaningEvidenceMode.DECLARED_BY_CALLER,
            parameters={
                "baseline_centering_established_by": baseline_source,
                "asserted_semantics": (
                    "condition_relative_log2_reference_control_centered"
                ),
            },
            output_table_fingerprint=matrix_fingerprint,
        )
        quantitative_state = established_state.transition_quantitative_meaning(
            target_quantity=QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
            provenance=meaning_provenance,
            authority=(sps_reference_quantitative_meaning_transition_authority()),
        )
        return cls(
            dataset_id=resolved_dataset_id,
            intensities=owned_intensities,
            condition_by_sample=condition_by_sample,
            intensity_scale_state=quantitative_state,
            organism=organism,
            baseline_context=baseline_context,
            reference_context=reference_context,
            source_name=source_name,
            source_version=source_version,
            source_uri=source_uri,
        )

    @property
    def intensities(self) -> pd.DataFrame:
        """Return a mutation-isolated snapshot of the reference matrix."""

        if not isinstance(self._owned_intensities, pd.DataFrame):
            raise PhosPyInputError(
                f"SPS reference dataset {self.dataset_id!r} intensities must be a "
                "pandas DataFrame"
            )
        return self._owned_intensities.copy(deep=True)

    def _intensities_snapshot(self) -> object:
        """Return an owner-detached matrix for trusted package collaborators."""

        if isinstance(self._owned_intensities, pd.DataFrame):
            return self._owned_intensities.copy(deep=True)
        return self._owned_intensities

    @property
    def sample_conditions(self) -> dict[str, object]:
        """Return a fresh JSON-shaped sample-to-condition mapping."""

        return thaw_json_mapping(
            self.condition_by_sample,
            field_name=(
                f"sps_reference_dataset[{self.dataset_id!r}].condition_by_sample"
            ),
        )


@dataclass(frozen=True, slots=True)
class SpsDiscoveryValidationIssue:
    """One structured SPS reference-input validation failure."""

    code: str
    message: str
    dataset_id: str | None = None
    field_name: str | None = None


@dataclass(frozen=True, slots=True)
class SpsDiscoveryValidationResult:
    """Aggregated validation outcome for a discovery request."""

    issues: tuple[SpsDiscoveryValidationIssue, ...]
    reference_dataset_count: int
    potential_overlap_site_count: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "issues", tuple(self.issues))
        if not all(
            isinstance(cast(object, issue), SpsDiscoveryValidationIssue)
            for issue in self.issues
        ):
            raise PhosPyInputError(
                "sps_discovery_validation.issues must contain only "
                "SpsDiscoveryValidationIssue values"
            )
        require_non_negative_int(
            self.reference_dataset_count,
            field_name="sps_discovery_validation.reference_dataset_count",
        )
        require_non_negative_int(
            self.potential_overlap_site_count,
            field_name="sps_discovery_validation.potential_overlap_site_count",
        )

    @property
    def valid(self) -> bool:
        return not self.issues

    def to_payload(self) -> dict[str, object]:
        return {
            "valid": self.valid,
            "reference_dataset_count": self.reference_dataset_count,
            "potential_overlap_site_count": self.potential_overlap_site_count,
            "issues": [
                {
                    "code": issue.code,
                    "message": issue.message,
                    "dataset_id": issue.dataset_id,
                    "field_name": issue.field_name,
                }
                for issue in self.issues
            ],
        }


class SpsDiscoveryValidationError(ReferenceValidationError):
    """Invalid SPS discovery references with a structured validation result."""

    validation_result: SpsDiscoveryValidationResult

    def __init__(self, validation_result: SpsDiscoveryValidationResult) -> None:
        self.validation_result = validation_result
        codes = ", ".join(issue.code for issue in validation_result.issues)
        super().__init__(f"SPS discovery request validation failed: {codes}")


def require_sps_quantitative_state(
    state: object,
    *,
    field_name: str,
    intensity_fingerprint: TableFingerprint,
) -> IntensityScaleState:
    """Require the governed quantitative domain consumed by SPS discovery."""

    if not isinstance(state, IntensityScaleState):
        raise PhosPyInputError(f"{field_name} must be an IntensityScaleState")
    establishment = state.establishment_provenance
    if not state.is_established or establishment is None:
        raise PhosPyInputError(
            f"{field_name} must carry established intensity-scale provenance"
        )
    if establishment.evidence_level is IntensityScaleEvidenceLevel.UNKNOWN:
        raise PhosPyInputError(
            f"{field_name} must carry known intensity-scale establishment evidence"
        )
    if state.kind is not IntensityScaleKind.LOG2:
        raise PhosPyInputError(f"{field_name} must establish the log2 scale")
    if state.quantity is not QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE:
        raise PhosPyInputError(
            f"{field_name} must establish condition-relative contrast_log2_fold_change "
            "meaning"
        )
    meaning_provenance = state.quantitative_meaning_provenance
    if meaning_provenance is None:
        raise PhosPyInputError(
            f"{field_name} must carry quantitative-meaning establishment provenance"
        )
    if not is_accepted_sps_meaning_evidence(meaning_provenance):
        raise PhosPyInputError(
            f"{field_name} quantitative-meaning evidence must establish that "
            "reference/control centring occurred"
        )
    if meaning_provenance.target_quantity is not state.quantity:
        raise PhosPyInputError(
            f"{field_name} quantitative-meaning provenance must match its meaning"
        )
    if not sps_scale_evidence_matches_matrix(
        establishment,
        intensity_fingerprint=intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name} scale establishment evidence must be bound to the "
            "submitted SPS reference intensity matrix"
        )
    if not sps_meaning_evidence_matches_matrix(
        meaning_provenance,
        intensity_fingerprint=intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name} reference/control-centering evidence must be bound to "
            "the submitted SPS reference intensity matrix"
        )
    return state


def sps_reference_intensity_fingerprint(
    *,
    dataset_id: str,
    intensities: pd.DataFrame,
) -> TableFingerprint:
    return fingerprint_table_normalized_axes(
        intensities,
        name=f"sps_reference.{dataset_id}.intensities",
    )


def is_accepted_sps_meaning_evidence(
    provenance: QuantitativeMeaningTransitionProvenance,
) -> bool:
    evidence_mode = cast(
        QuantitativeMeaningEvidenceMode,
        provenance.evidence_mode,
    )
    if evidence_mode in {
        QuantitativeMeaningEvidenceMode.DERIVED_BY_PHOSPY_OPERATION,
        QuantitativeMeaningEvidenceMode.RESTORED_FROM_TRUSTED_SERIALIZED_PROVENANCE,
    }:
        return True
    return (
        evidence_mode is QuantitativeMeaningEvidenceMode.DECLARED_BY_CALLER
        and provenance.operation_id == SPS_REFERENCE_CENTERING_ASSERTION_OPERATION
    )


def sps_scale_evidence_matches_matrix(
    provenance: IntensityScaleEstablishmentProvenance,
    *,
    intensity_fingerprint: TableFingerprint,
) -> bool:
    payload = provenance.parameters.get(SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER)
    return sps_fingerprint_payload_matches(
        payload,
        intensity_fingerprint=intensity_fingerprint,
    )


def sps_meaning_evidence_matches_matrix(
    provenance: QuantitativeMeaningTransitionProvenance,
    *,
    intensity_fingerprint: TableFingerprint,
) -> bool:
    return sps_fingerprint_payload_matches(
        provenance.output_table_fingerprint,
        intensity_fingerprint=intensity_fingerprint,
    )


def sps_fingerprint_payload_matches(
    payload: object,
    *,
    intensity_fingerprint: TableFingerprint,
) -> bool:
    if not isinstance(payload, Mapping):
        return False
    try:
        resolved = table_fingerprint_from_payload(cast(Mapping[str, object], payload))
    except (PhosPyInputError, TypeError, ValueError):
        return False
    return resolved == intensity_fingerprint


@dataclass(frozen=True, slots=True, eq=False)
class SpsDiscoveryRequest:
    """Passive caller-intent payload for SPS discovery execution."""

    reference_datasets: Sequence[SpsReferenceDataset]
    config: SpsDiscoveryConfig


__all__ = [
    "SPS_DISCOVERY_ALGORITHM_ID",
    "SPS_DISCOVERY_ALGORITHM_VERSION",
    "SPS_DISCOVERY_CROSS_COUNT_ORDERING_POLICY",
    "SPS_DISCOVERY_CONTROL_SOURCE_TYPE",
    "SPS_DISCOVERY_SELECTION_METHOD_CONSENSUS_STABILITY",
    "SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING",
    "SpsDiscoveryConfig",
    "SpsDiscoveryRequest",
    "SpsDiscoverySelectionMethod",
    "SpsDiscoveryTieHandling",
    "SpsDiscoveryValidationError",
    "SpsDiscoveryValidationIssue",
    "SpsDiscoveryValidationResult",
    "SpsReferenceDataset",
]
