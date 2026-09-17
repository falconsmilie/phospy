"""Typed contracts and numerical science for reference-data SPS discovery.

SPS discovery is deliberately independent of native SPS/RUV-style correction:
it derives negative controls from explicit reference datasets and never treats
the target correction matrix as implicit reference evidence.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, cast

import pandas as pd
from scipy.stats import chi2

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.hashing import fingerprint_table_normalized_axes
from phospy.provenance.immutability import freeze_json_mapping, thaw_json_mapping
from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.science.sites.validation import require_site_key_index
from phospy.science.transformations._authority import (
    sps_reference_quantitative_meaning_transition_authority,
)
from phospy.science.transformations.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentProvenance,
    IntensityScaleEstablishmentSource,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningTransitionProvenance,
)
from phospy.science.transformations.transformers import IdentityTransformer

SPS_DISCOVERY_SELECTION_METHOD_CONSENSUS_STABILITY = "consensus_stability"
SPS_DISCOVERY_ALGORITHM_ID = "phospy_sps_consensus_stability"
SPS_DISCOVERY_ALGORITHM_VERSION = "1.0.0"
SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING = "site_key_ascending"
SPS_DISCOVERY_CONTROL_SOURCE_TYPE = "sps_discovery"
SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER = "sps_reference_matrix_fingerprint"
SPS_REFERENCE_CENTERING_ASSERTION_OPERATION = (
    "phospy.sps_reference.external_reference_centering_assertion"
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
        _require_positive_int(self.top_n, field_name="sps_discovery.top_n")
        _require_int_at_least(
            self.minimum_reference_datasets,
            minimum=2,
            field_name="sps_discovery.minimum_reference_datasets",
        )
        _require_int_at_least(
            self.minimum_datasets_per_site,
            minimum=2,
            field_name="sps_discovery.minimum_datasets_per_site",
        )
        _require_positive_int(
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

        resolved = _require_payload_mapping(payload, field_name="sps_discovery_config")
        return cls(
            top_n=_require_payload_int(
                resolved.get("top_n"), field_name="sps_discovery_config.top_n"
            ),
            minimum_reference_datasets=_require_payload_int(
                resolved.get("minimum_reference_datasets"),
                field_name="sps_discovery_config.minimum_reference_datasets",
            ),
            minimum_datasets_per_site=_require_payload_int(
                resolved.get("minimum_datasets_per_site"),
                field_name="sps_discovery_config.minimum_datasets_per_site",
            ),
            minimum_shared_sites=_require_payload_int(
                resolved.get("minimum_shared_sites"),
                field_name="sps_discovery_config.minimum_shared_sites",
            ),
            tie_handling=cast(
                SpsDiscoveryTieHandling,
                _require_payload_string(
                    resolved.get("tie_handling"),
                    field_name="sps_discovery_config.tie_handling",
                ),
            ),
            selection_method=cast(
                SpsDiscoverySelectionMethod,
                _require_payload_string(
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
    """

    dataset_id: str
    _owned_intensities: object = field(repr=False)
    condition_by_sample: Mapping[str, str]
    intensity_scale_state: IntensityScaleState
    source_name: str | None = None
    source_version: str | None = None
    source_uri: str | None = None

    def __init__(
        self,
        dataset_id: str,
        intensities: pd.DataFrame,
        condition_by_sample: Mapping[str, str],
        intensity_scale_state: IntensityScaleState,
        source_name: str | None = None,
        source_version: str | None = None,
        source_uri: str | None = None,
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
        object.__setattr__(self, "source_name", source_name)
        object.__setattr__(self, "source_version", source_version)
        object.__setattr__(self, "source_uri", source_uri)
        self.__post_init__()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dataset_id",
            _require_non_empty_string(
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
        for field_name in ("source_name", "source_version", "source_uri"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(
                    self,
                    field_name,
                    _require_non_empty_string(
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
        """

        resolved_dataset_id = _require_non_empty_string(
            dataset_id,
            field_name="sps_reference_dataset.dataset_id",
        )
        scale_source = _require_non_empty_string(
            log2_scale_established_by,
            field_name=("sps_reference_dataset.log2_scale_established_by"),
        )
        baseline_source = _require_non_empty_string(
            baseline_centering_established_by,
            field_name=("sps_reference_dataset.baseline_centering_established_by"),
        )
        if not isinstance(cast(object, intensities), pd.DataFrame):
            raise PhosPyInputError(
                "sps_reference_dataset.intensities must be a pandas DataFrame"
            )
        owned_intensities = intensities.copy(deep=True)
        matrix_fingerprint = _sps_reference_intensity_fingerprint(
            dataset_id=resolved_dataset_id,
            intensities=owned_intensities,
        )
        fingerprint_payload = table_fingerprint_to_payload(matrix_fingerprint)
        declared_state = IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(),
            quantity=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )
        established_state = (
            DatasetIntensityScaleResolver(transformer=IdentityTransformer())
            .run(
                phospho=owned_intensities,
                total=None,
                expected_scale_kind=IntensityScaleKind.LOG2,
                declared_input_scale_state=declared_state,
                declared_input_establishment_mode=(
                    IntensityScaleEstablishmentMode.DECLARED
                ),
                input_declaration_source=scale_source,
                scale_establishment_parameters={
                    SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER: fingerprint_payload,
                },
            )
            .intensity_scale_state
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

    def _to_provenance(
        self,
        *,
        sites_passing_required_data: int | None = None,
        sites_entering_consensus: int | None = None,
    ) -> SpsReferenceDatasetProvenance:
        """Build source provenance for the validated workflow boundary."""

        return _reference_provenance(
            self,
            sites_passing_required_data=sites_passing_required_data,
            sites_entering_consensus=sites_entering_consensus,
        )


@dataclass(frozen=True, slots=True)
class SpsSampleConditionAssignment:
    """Condition/group assignment for one reference sample."""

    sample_id: str
    condition: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "sample_id",
            _require_non_empty_string(
                self.sample_id, field_name="sps_sample_condition.sample_id"
            ),
        )
        object.__setattr__(
            self,
            "condition",
            _require_non_empty_string(
                self.condition, field_name="sps_sample_condition.condition"
            ),
        )

    def to_payload(self) -> dict[str, object]:
        return {"sample_id": self.sample_id, "condition": self.condition}

    @classmethod
    def from_payload(
        cls, payload: Mapping[str, object]
    ) -> SpsSampleConditionAssignment:
        """Restore a validated sample-condition assignment."""

        resolved = _require_payload_mapping(payload, field_name="sps_sample_condition")
        return cls(
            sample_id=_require_payload_string(
                resolved.get("sample_id"),
                field_name="sps_sample_condition.sample_id",
            ),
            condition=_require_payload_string(
                resolved.get("condition"),
                field_name="sps_sample_condition.condition",
            ),
        )


@dataclass(frozen=True, slots=True)
class SpsReferenceDatasetProvenance:
    """Immutable identity and fingerprint for one contributing reference."""

    dataset_id: str
    site_count: int
    sample_count: int
    sample_conditions: tuple[SpsSampleConditionAssignment, ...]
    intensity_fingerprint: TableFingerprint
    intensity_scale_kind: IntensityScaleKind
    quantitative_meaning: QuantitativeMeaning
    intensity_scale_establishment: IntensityScaleEstablishmentProvenance
    quantitative_meaning_establishment: QuantitativeMeaningTransitionProvenance
    source_name: str | None = None
    source_version: str | None = None
    source_uri: str | None = None
    sites_passing_required_data: int | None = None
    sites_entering_consensus: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dataset_id",
            _require_non_empty_string(
                self.dataset_id, field_name="sps_source_provenance.dataset_id"
            ),
        )
        _require_non_negative_int(
            self.site_count, field_name="sps_source_provenance.site_count"
        )
        _require_non_negative_int(
            self.sample_count, field_name="sps_source_provenance.sample_count"
        )
        object.__setattr__(self, "sample_conditions", tuple(self.sample_conditions))
        if not all(
            isinstance(cast(object, item), SpsSampleConditionAssignment)
            for item in self.sample_conditions
        ):
            raise PhosPyInputError(
                "sps_source_provenance.sample_conditions must contain only "
                "SpsSampleConditionAssignment values"
            )
        if len(self.sample_conditions) != self.sample_count:
            raise PhosPyInputError(
                "sps_source_provenance.sample_conditions length must equal sample_count"
            )
        sample_ids = tuple(item.sample_id for item in self.sample_conditions)
        if len(set(sample_ids)) != len(sample_ids):
            raise PhosPyInputError(
                "sps_source_provenance.sample_conditions must not contain "
                "duplicate sample_id values"
            )
        if not isinstance(cast(object, self.intensity_fingerprint), TableFingerprint):
            raise PhosPyInputError(
                "sps_source_provenance.intensity_fingerprint must be a TableFingerprint"
            )
        if (
            self.intensity_fingerprint.rows != self.site_count
            or self.intensity_fingerprint.columns != self.sample_count
        ):
            raise PhosPyInputError(
                "sps_source_provenance intensity fingerprint shape must match "
                "site_count and sample_count"
            )
        if self.intensity_fingerprint.index_name != "site_key":
            raise PhosPyInputError(
                "sps_source_provenance.intensity_fingerprint must describe a "
                "governed 'site_key' index"
            )
        _require_sps_provenance_quantitative_contract(self)
        for field_name in ("source_name", "source_version", "source_uri"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(
                    self,
                    field_name,
                    _require_non_empty_string(
                        value,
                        field_name=f"sps_source_provenance.{field_name}",
                    ),
                )
        for field_name in (
            "sites_passing_required_data",
            "sites_entering_consensus",
        ):
            value = getattr(self, field_name)
            if value is not None:
                _require_non_negative_int(
                    value,
                    field_name=f"sps_source_provenance.{field_name}",
                )
                if value > self.site_count:
                    raise PhosPyInputError(
                        f"sps_source_provenance.{field_name} must not exceed site_count"
                    )
        if (
            self.sites_passing_required_data is not None
            and self.sites_entering_consensus is not None
            and self.sites_entering_consensus > self.sites_passing_required_data
        ):
            raise PhosPyInputError(
                "sps_source_provenance.sites_entering_consensus must not exceed "
                "sites_passing_required_data"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "dataset_id": self.dataset_id,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "source_uri": self.source_uri,
            "site_count": self.site_count,
            "sample_count": self.sample_count,
            "sites_passing_required_data": self.sites_passing_required_data,
            "sites_entering_consensus": self.sites_entering_consensus,
            "sample_conditions": [item.to_payload() for item in self.sample_conditions],
            "intensity_fingerprint": table_fingerprint_to_payload(
                self.intensity_fingerprint
            ),
            "quantitative_state": {
                "scale": self.intensity_scale_kind.value,
                "quantitative_meaning": self.quantitative_meaning.value,
                "scale_establishment": self.intensity_scale_establishment.to_payload(),
                "quantitative_meaning_establishment": (
                    self.quantitative_meaning_establishment.to_payload()
                ),
            },
        }

    @classmethod
    def from_payload(
        cls, payload: Mapping[str, object]
    ) -> SpsReferenceDatasetProvenance:
        """Restore matrix-free SPS reference provenance from serialized evidence."""

        resolved = _require_payload_mapping(payload, field_name="sps_source_provenance")
        quantitative_state = _require_payload_mapping(
            resolved.get("quantitative_state"),
            field_name="sps_source_provenance.quantitative_state",
        )
        scale_establishment_payload = _require_payload_mapping(
            quantitative_state.get("scale_establishment"),
            field_name=("sps_source_provenance.quantitative_state.scale_establishment"),
        )
        meaning_establishment_payload = _require_payload_mapping(
            quantitative_state.get("quantitative_meaning_establishment"),
            field_name=(
                "sps_source_provenance.quantitative_state."
                "quantitative_meaning_establishment"
            ),
        )
        fingerprint_payload = _require_payload_mapping(
            resolved.get("intensity_fingerprint"),
            field_name="sps_source_provenance.intensity_fingerprint",
        )
        sample_condition_payloads = _require_payload_sequence(
            resolved.get("sample_conditions"),
            field_name="sps_source_provenance.sample_conditions",
        )
        try:
            intensity_scale_kind = IntensityScaleKind(
                _require_payload_string(
                    quantitative_state.get("scale"),
                    field_name="sps_source_provenance.quantitative_state.scale",
                )
            )
            quantitative_meaning = QuantitativeMeaning(
                _require_payload_string(
                    quantitative_state.get("quantitative_meaning"),
                    field_name=(
                        "sps_source_provenance.quantitative_state.quantitative_meaning"
                    ),
                )
            )
            scale_establishment = _scale_establishment_from_payload(
                scale_establishment_payload
            )
            meaning_establishment = (
                QuantitativeMeaningTransitionProvenance.from_payload(
                    meaning_establishment_payload
                )
            )
            intensity_fingerprint = table_fingerprint_from_payload(fingerprint_payload)
        except (TypeError, ValueError) as exc:
            raise PhosPyInputError(
                "sps_source_provenance contains invalid quantitative or "
                "fingerprint evidence"
            ) from exc
        return cls(
            dataset_id=_require_payload_string(
                resolved.get("dataset_id"),
                field_name="sps_source_provenance.dataset_id",
            ),
            source_name=_optional_payload_string(
                resolved.get("source_name"),
                field_name="sps_source_provenance.source_name",
            ),
            source_version=_optional_payload_string(
                resolved.get("source_version"),
                field_name="sps_source_provenance.source_version",
            ),
            source_uri=_optional_payload_string(
                resolved.get("source_uri"),
                field_name="sps_source_provenance.source_uri",
            ),
            site_count=_require_payload_int(
                resolved.get("site_count"),
                field_name="sps_source_provenance.site_count",
            ),
            sample_count=_require_payload_int(
                resolved.get("sample_count"),
                field_name="sps_source_provenance.sample_count",
            ),
            sites_passing_required_data=_optional_payload_int(
                resolved.get("sites_passing_required_data"),
                field_name="sps_source_provenance.sites_passing_required_data",
            ),
            sites_entering_consensus=_optional_payload_int(
                resolved.get("sites_entering_consensus"),
                field_name="sps_source_provenance.sites_entering_consensus",
            ),
            sample_conditions=tuple(
                SpsSampleConditionAssignment.from_payload(
                    _require_payload_mapping(
                        item,
                        field_name=(
                            f"sps_source_provenance.sample_conditions[{position}]"
                        ),
                    )
                )
                for position, item in enumerate(sample_condition_payloads)
            ),
            intensity_fingerprint=intensity_fingerprint,
            intensity_scale_kind=intensity_scale_kind,
            quantitative_meaning=quantitative_meaning,
            intensity_scale_establishment=scale_establishment,
            quantitative_meaning_establishment=meaning_establishment,
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
        _require_non_negative_int(
            self.reference_dataset_count,
            field_name="sps_discovery_validation.reference_dataset_count",
        )
        _require_non_negative_int(
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


@dataclass(frozen=True, slots=True)
class SpsDatasetSiteStatistic:
    """Reference-specific rank/stability evidence for one candidate site."""

    dataset_id: str
    rank: int
    stability_score: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dataset_id",
            _require_non_empty_string(
                self.dataset_id, field_name="sps_site_statistic.dataset_id"
            ),
        )
        _require_positive_int(self.rank, field_name="sps_site_statistic.rank")
        stability_score = _require_finite_float(
            self.stability_score, field_name="sps_site_statistic.stability_score"
        )
        if stability_score < 0.0:
            raise PhosPyInputError(
                "sps_site_statistic.stability_score must be non-negative"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "dataset_id": self.dataset_id,
            "rank": self.rank,
            "stability_score": self.stability_score,
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsDatasetSiteStatistic:
        """Restore one validated reference-specific site statistic."""

        resolved = _require_payload_mapping(payload, field_name="sps_site_statistic")
        return cls(
            dataset_id=_require_payload_string(
                resolved.get("dataset_id"),
                field_name="sps_site_statistic.dataset_id",
            ),
            rank=_require_payload_int(
                resolved.get("rank"), field_name="sps_site_statistic.rank"
            ),
            stability_score=_require_payload_float(
                resolved.get("stability_score"),
                field_name="sps_site_statistic.stability_score",
            ),
        )


@dataclass(frozen=True, slots=True)
class SpsSiteStabilityRecord:
    """Consensus ordering and contributing evidence for one SPS candidate."""

    site_key: str
    consensus_rank: int
    consensus_stability_score: float
    contributing_dataset_count: int
    dataset_statistics: tuple[SpsDatasetSiteStatistic, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "site_key",
            _require_site_key(self.site_key, field_name="sps_site_record.site_key"),
        )
        _require_positive_int(
            self.consensus_rank, field_name="sps_site_record.consensus_rank"
        )
        consensus_stability_score = _require_finite_float(
            self.consensus_stability_score,
            field_name="sps_site_record.consensus_stability_score",
        )
        if not 0.0 <= consensus_stability_score <= 1.0:
            raise PhosPyInputError(
                "sps_site_record.consensus_stability_score must be between 0 and 1"
            )
        _require_positive_int(
            self.contributing_dataset_count,
            field_name="sps_site_record.contributing_dataset_count",
        )
        object.__setattr__(self, "dataset_statistics", tuple(self.dataset_statistics))
        if not all(
            isinstance(cast(object, item), SpsDatasetSiteStatistic)
            for item in self.dataset_statistics
        ):
            raise PhosPyInputError(
                "sps_site_record.dataset_statistics must contain only "
                "SpsDatasetSiteStatistic values"
            )
        if self.contributing_dataset_count != len(self.dataset_statistics):
            raise PhosPyInputError(
                "sps_site_record.contributing_dataset_count must equal the "
                "number of dataset_statistics"
            )
        dataset_ids = tuple(item.dataset_id for item in self.dataset_statistics)
        if len(set(dataset_ids)) != len(dataset_ids):
            raise PhosPyInputError(
                "sps_site_record.dataset_statistics must not contain duplicate "
                "dataset_id values"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "site_key": self.site_key,
            "consensus_rank": self.consensus_rank,
            "consensus_stability_score": self.consensus_stability_score,
            "contributing_dataset_count": self.contributing_dataset_count,
            "dataset_statistics": [
                item.to_payload() for item in self.dataset_statistics
            ],
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsSiteStabilityRecord:
        """Restore one validated consensus-ranking record."""

        resolved = _require_payload_mapping(payload, field_name="sps_site_record")
        statistics = _require_payload_sequence(
            resolved.get("dataset_statistics"),
            field_name="sps_site_record.dataset_statistics",
        )
        return cls(
            site_key=_require_payload_string(
                resolved.get("site_key"), field_name="sps_site_record.site_key"
            ),
            consensus_rank=_require_payload_int(
                resolved.get("consensus_rank"),
                field_name="sps_site_record.consensus_rank",
            ),
            consensus_stability_score=_require_payload_float(
                resolved.get("consensus_stability_score"),
                field_name="sps_site_record.consensus_stability_score",
            ),
            contributing_dataset_count=_require_payload_int(
                resolved.get("contributing_dataset_count"),
                field_name="sps_site_record.contributing_dataset_count",
            ),
            dataset_statistics=tuple(
                SpsDatasetSiteStatistic.from_payload(
                    _require_payload_mapping(
                        item,
                        field_name=f"sps_site_record.dataset_statistics[{position}]",
                    )
                )
                for position, item in enumerate(statistics)
            ),
        )


@dataclass(frozen=True, slots=True)
class SpsSelectionBoundaryCounts:
    """Site counts at important SPS selection boundaries."""

    total_unique_sites: int
    sites_meeting_dataset_overlap: int
    sites_with_valid_stability: int
    sites_ranked: int
    sites_selected: int

    def __post_init__(self) -> None:
        values = (
            self.total_unique_sites,
            self.sites_meeting_dataset_overlap,
            self.sites_with_valid_stability,
            self.sites_ranked,
            self.sites_selected,
        )
        for field_name, value in zip(
            (
                "total_unique_sites",
                "sites_meeting_dataset_overlap",
                "sites_with_valid_stability",
                "sites_ranked",
                "sites_selected",
            ),
            values,
            strict=True,
        ):
            _require_non_negative_int(
                value, field_name=f"sps_selection_boundaries.{field_name}"
            )
        if any(left < right for left, right in zip(values, values[1:], strict=False)):
            raise PhosPyInputError(
                "sps_selection_boundaries counts must be monotonically non-increasing"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "total_unique_sites": self.total_unique_sites,
            "sites_meeting_dataset_overlap": self.sites_meeting_dataset_overlap,
            "sites_with_valid_stability": self.sites_with_valid_stability,
            "sites_ranked": self.sites_ranked,
            "sites_selected": self.sites_selected,
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsSelectionBoundaryCounts:
        """Restore validated SPS boundary counts."""

        resolved = _require_payload_mapping(
            payload, field_name="sps_selection_boundaries"
        )
        return cls(
            total_unique_sites=_require_payload_int(
                resolved.get("total_unique_sites"),
                field_name="sps_selection_boundaries.total_unique_sites",
            ),
            sites_meeting_dataset_overlap=_require_payload_int(
                resolved.get("sites_meeting_dataset_overlap"),
                field_name=("sps_selection_boundaries.sites_meeting_dataset_overlap"),
            ),
            sites_with_valid_stability=_require_payload_int(
                resolved.get("sites_with_valid_stability"),
                field_name="sps_selection_boundaries.sites_with_valid_stability",
            ),
            sites_ranked=_require_payload_int(
                resolved.get("sites_ranked"),
                field_name="sps_selection_boundaries.sites_ranked",
            ),
            sites_selected=_require_payload_int(
                resolved.get("sites_selected"),
                field_name="sps_selection_boundaries.sites_selected",
            ),
        )


@dataclass(frozen=True, slots=True)
class SpsDiscoveryProvenance:
    """Immutable, comparable provenance for one SPS discovery result."""

    source_datasets: tuple[SpsReferenceDatasetProvenance, ...]
    config: SpsDiscoveryConfig
    requested_control_count: int
    actual_control_count: int
    selection_boundaries: SpsSelectionBoundaryCounts
    algorithm_id: str = SPS_DISCOVERY_ALGORITHM_ID
    algorithm_version: str = SPS_DISCOVERY_ALGORITHM_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_datasets", tuple(self.source_datasets))
        if not self.source_datasets:
            raise PhosPyInputError(
                "sps_discovery_provenance.source_datasets must not be empty"
            )
        if not all(
            isinstance(cast(object, item), SpsReferenceDatasetProvenance)
            for item in self.source_datasets
        ):
            raise PhosPyInputError(
                "sps_discovery_provenance.source_datasets must contain only "
                "SpsReferenceDatasetProvenance values"
            )
        dataset_ids = tuple(item.dataset_id for item in self.source_datasets)
        if len(set(dataset_ids)) != len(dataset_ids):
            raise PhosPyInputError(
                "sps_discovery_provenance.source_datasets must not contain "
                "duplicate dataset_id values"
            )
        if not isinstance(cast(object, self.config), SpsDiscoveryConfig):
            raise PhosPyInputError(
                "sps_discovery_provenance.config must be an SpsDiscoveryConfig"
            )
        if len(self.source_datasets) < self.config.minimum_reference_datasets:
            raise PhosPyInputError(
                "sps_discovery_provenance.source_datasets does not meet "
                "config.minimum_reference_datasets"
            )
        if not isinstance(
            cast(object, self.selection_boundaries), SpsSelectionBoundaryCounts
        ):
            raise PhosPyInputError(
                "sps_discovery_provenance.selection_boundaries must be an "
                "SpsSelectionBoundaryCounts"
            )
        _require_positive_int(
            self.requested_control_count,
            field_name="sps_discovery_provenance.requested_control_count",
        )
        _require_non_negative_int(
            self.actual_control_count,
            field_name="sps_discovery_provenance.actual_control_count",
        )
        if self.requested_control_count != self.config.top_n:
            raise PhosPyInputError(
                "sps_discovery_provenance.requested_control_count must equal "
                "config.top_n"
            )
        if self.actual_control_count > self.requested_control_count:
            raise PhosPyInputError(
                "sps_discovery_provenance.actual_control_count must not exceed "
                "requested_control_count"
            )
        if self.selection_boundaries.sites_selected != self.actual_control_count:
            raise PhosPyInputError(
                "sps_discovery_provenance selection boundary sites_selected must "
                "equal actual_control_count"
            )
        object.__setattr__(
            self,
            "algorithm_id",
            _require_non_empty_string(
                self.algorithm_id,
                field_name="sps_discovery_provenance.algorithm_id",
            ),
        )
        object.__setattr__(
            self,
            "algorithm_version",
            _require_non_empty_string(
                self.algorithm_version,
                field_name="sps_discovery_provenance.algorithm_version",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "source_datasets": [item.to_payload() for item in self.source_datasets],
            "algorithm_id": self.algorithm_id,
            "algorithm_version": self.algorithm_version,
            "selection_method": self.config.selection_method,
            "requested_control_count": self.requested_control_count,
            "actual_control_count": self.actual_control_count,
            "algorithm_parameters": self.config.to_payload(),
            "selection_boundaries": self.selection_boundaries.to_payload(),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsDiscoveryProvenance:
        """Restore validated SPS discovery provenance without source matrices."""

        resolved = _require_payload_mapping(
            payload, field_name="sps_discovery_provenance"
        )
        source_payloads = _require_payload_sequence(
            resolved.get("source_datasets"),
            field_name="sps_discovery_provenance.source_datasets",
        )
        config = SpsDiscoveryConfig.from_payload(
            _require_payload_mapping(
                resolved.get("algorithm_parameters"),
                field_name="sps_discovery_provenance.algorithm_parameters",
            )
        )
        serialized_selection_method = _require_payload_string(
            resolved.get("selection_method"),
            field_name="sps_discovery_provenance.selection_method",
        )
        if serialized_selection_method != config.selection_method:
            raise PhosPyInputError(
                "sps_discovery_provenance.selection_method must match "
                "algorithm_parameters.selection_method"
            )
        return cls(
            source_datasets=tuple(
                SpsReferenceDatasetProvenance.from_payload(
                    _require_payload_mapping(
                        item,
                        field_name=(
                            f"sps_discovery_provenance.source_datasets[{position}]"
                        ),
                    )
                )
                for position, item in enumerate(source_payloads)
            ),
            config=config,
            requested_control_count=_require_payload_int(
                resolved.get("requested_control_count"),
                field_name="sps_discovery_provenance.requested_control_count",
            ),
            actual_control_count=_require_payload_int(
                resolved.get("actual_control_count"),
                field_name="sps_discovery_provenance.actual_control_count",
            ),
            selection_boundaries=SpsSelectionBoundaryCounts.from_payload(
                _require_payload_mapping(
                    resolved.get("selection_boundaries"),
                    field_name="sps_discovery_provenance.selection_boundaries",
                )
            ),
            algorithm_id=_require_payload_string(
                resolved.get("algorithm_id"),
                field_name="sps_discovery_provenance.algorithm_id",
            ),
            algorithm_version=_require_payload_string(
                resolved.get("algorithm_version"),
                field_name="sps_discovery_provenance.algorithm_version",
            ),
        )


@dataclass(frozen=True, slots=True)
class SpsDiscoveryResult:
    """Ranked SPS discovery output convertible to ``ControlSiteSet``."""

    selected_site_keys: tuple[str, ...]
    site_ranking: tuple[SpsSiteStabilityRecord, ...]
    provenance: SpsDiscoveryProvenance

    def __post_init__(self) -> None:
        object.__setattr__(self, "selected_site_keys", tuple(self.selected_site_keys))
        object.__setattr__(self, "site_ranking", tuple(self.site_ranking))
        if not isinstance(cast(object, self.provenance), SpsDiscoveryProvenance):
            raise PhosPyInputError(
                "sps_discovery_result.provenance must be an SpsDiscoveryProvenance"
            )
        if not all(
            isinstance(cast(object, record), SpsSiteStabilityRecord)
            for record in self.site_ranking
        ):
            raise PhosPyInputError(
                "sps_discovery_result.site_ranking must contain only "
                "SpsSiteStabilityRecord values"
            )
        for position, site_key in enumerate(self.selected_site_keys):
            _require_site_key(
                site_key,
                field_name=f"sps_discovery_result.selected_site_keys[{position}]",
            )
        if len(set(self.selected_site_keys)) != len(self.selected_site_keys):
            raise PhosPyInputError(
                "sps_discovery_result.selected_site_keys must be unique"
            )
        expected_ranks = tuple(range(1, len(self.site_ranking) + 1))
        observed_ranks = tuple(item.consensus_rank for item in self.site_ranking)
        if observed_ranks != expected_ranks:
            raise PhosPyInputError(
                "sps_discovery_result.site_ranking consensus ranks must be "
                "contiguous and ordered from 1"
            )
        for previous, current in zip(
            self.site_ranking,
            self.site_ranking[1:],
            strict=False,
        ):
            if previous.consensus_stability_score < current.consensus_stability_score:
                raise PhosPyInputError(
                    "sps_discovery_result.site_ranking consensus stability scores "
                    "must be non-increasing"
                )
        if self.provenance.config.tie_handling == (
            SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING
        ):
            for previous, current in zip(
                self.site_ranking,
                self.site_ranking[1:],
                strict=False,
            ):
                if (
                    previous.consensus_stability_score
                    == current.consensus_stability_score
                    and previous.site_key > current.site_key
                ):
                    raise PhosPyInputError(
                        "sps_discovery_result.site_ranking equal consensus scores "
                        "must use ascending site_key tie handling"
                    )
        ranked_site_keys = tuple(item.site_key for item in self.site_ranking)
        if len(set(ranked_site_keys)) != len(ranked_site_keys):
            raise PhosPyInputError(
                "sps_discovery_result.site_ranking site_key values must be unique"
            )
        if self.selected_site_keys != ranked_site_keys[: len(self.selected_site_keys)]:
            raise PhosPyInputError(
                "sps_discovery_result.selected_site_keys must be the leading "
                "consensus-ranked site_key values"
            )
        if self.provenance.actual_control_count != len(self.selected_site_keys):
            raise PhosPyInputError(
                "sps_discovery_result selected site count must equal provenance "
                "actual_control_count"
            )
        if self.provenance.selection_boundaries.sites_ranked != len(self.site_ranking):
            raise PhosPyInputError(
                "sps_discovery_result ranking length must equal provenance "
                "selection boundary sites_ranked"
            )
        expected_selected_count = min(
            self.provenance.config.top_n,
            len(self.site_ranking),
        )
        if len(self.selected_site_keys) != expected_selected_count:
            raise PhosPyInputError(
                "sps_discovery_result selected site count must equal exactly "
                "min(config.top_n, sites_ranked)"
            )
        if len(self.site_ranking) < self.provenance.config.minimum_shared_sites:
            raise PhosPyInputError(
                "sps_discovery_result ranking length must meet "
                "config.minimum_shared_sites"
            )
        source_ids = {item.dataset_id for item in self.provenance.source_datasets}
        statistic_counts_by_source: Counter[str] = Counter()
        for record in self.site_ranking:
            if (
                record.contributing_dataset_count
                < self.provenance.config.minimum_datasets_per_site
            ):
                raise PhosPyInputError(
                    "sps_discovery_result.site_ranking contains a site below "
                    "config.minimum_datasets_per_site"
                )
            statistic_ids = {item.dataset_id for item in record.dataset_statistics}
            unknown = sorted(statistic_ids.difference(source_ids))
            if unknown:
                raise PhosPyInputError(
                    "sps_discovery_result.site_ranking references unknown source "
                    f"dataset ids: {', '.join(unknown)}"
                )
            statistic_counts_by_source.update(statistic_ids)
        for source in self.provenance.source_datasets:
            statistic_count = statistic_counts_by_source[source.dataset_id]
            if (
                source.sites_entering_consensus is not None
                and source.sites_entering_consensus != statistic_count
            ):
                raise PhosPyInputError(
                    "sps_discovery_result source sites_entering_consensus must equal "
                    "the number of ranking records containing that dataset"
                )
        for record in self.site_ranking:
            for statistic in record.dataset_statistics:
                if statistic.rank > statistic_counts_by_source[statistic.dataset_id]:
                    raise PhosPyInputError(
                        "sps_discovery_result reference-specific rank must not exceed "
                        "that source dataset's consensus entry count"
                    )

    def to_control_site_set(self) -> ControlSiteSet:
        """Return selected SPS controls through the existing downstream contract."""

        return ControlSiteSet.from_site_keys(
            self.selected_site_keys,
            source_metadata=ControlSiteSourceMetadata(
                source_type=SPS_DISCOVERY_CONTROL_SOURCE_TYPE,
                identifier_namespace="site_key",
                source_name=self.provenance.algorithm_id,
                source_version=self.provenance.algorithm_version,
                selection_method=self.provenance.config.selection_method,
                metadata_missing_reason={
                    "organism": (
                        "SPS reference records do not carry a typed organism field; "
                        "organism coherence remains caller-audited reference context"
                    ),
                    "license": "retained in caller-owned reference source records",
                    "redistribution": (
                        "no external SPS reference data is bundled or redistributed"
                    ),
                },
            ),
            label="sps_discovered_control",
        )

    @property
    def control_site_set(self) -> ControlSiteSet:
        """Expose the existing SPS/RUV negative-control representation."""

        return self.to_control_site_set()

    def to_payload(self) -> dict[str, object]:
        return {
            "selected_site_keys": list(self.selected_site_keys),
            "site_ranking": [item.to_payload() for item in self.site_ranking],
            "provenance": self.provenance.to_payload(),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsDiscoveryResult:
        """Restore a validated discovery result from its matrix-free payload."""

        resolved = _require_payload_mapping(payload, field_name="sps_discovery_result")
        selected = _require_payload_sequence(
            resolved.get("selected_site_keys"),
            field_name="sps_discovery_result.selected_site_keys",
        )
        ranking = _require_payload_sequence(
            resolved.get("site_ranking"),
            field_name="sps_discovery_result.site_ranking",
        )
        return cls(
            selected_site_keys=tuple(
                _require_payload_string(
                    item,
                    field_name=f"sps_discovery_result.selected_site_keys[{position}]",
                )
                for position, item in enumerate(selected)
            ),
            site_ranking=tuple(
                SpsSiteStabilityRecord.from_payload(
                    _require_payload_mapping(
                        item,
                        field_name=f"sps_discovery_result.site_ranking[{position}]",
                    )
                )
                for position, item in enumerate(ranking)
            ),
            provenance=SpsDiscoveryProvenance.from_payload(
                _require_payload_mapping(
                    resolved.get("provenance"),
                    field_name="sps_discovery_result.provenance",
                )
            ),
        )


@dataclass(frozen=True, slots=True)
class _ReferenceStabilityEvidence:
    dataset_id: str
    scores_by_site: Mapping[str, float]


class SpsDiscoveryExecutor:
    """Derive a deterministic consensus SPS ranking from validated references.

    Replicates are averaged within condition.  The reference-specific change
    magnitude is the largest absolute condition mean for a site, matching the
    central ``getSPS`` stability statistic for condition-relative inputs.
    Smaller magnitudes are more stable.  Reference ranks are converted to the
    half-step empirical quantiles used by ``getSPS`` and combined with its
    Fisher-style chi-square survival score; larger consensus scores are more
    stable.

    Missing replicate cells are ignored only when their condition retains at
    least one finite observation.  A site with an entirely missing condition
    is not rankable in that reference.  It remains eligible only when it has
    valid evidence in ``minimum_datasets_per_site`` references.
    """

    def run(
        self,
        *,
        reference_datasets: Sequence[SpsReferenceDataset],
        config: SpsDiscoveryConfig,
    ) -> SpsDiscoveryResult:
        references = tuple(sorted(reference_datasets, key=lambda item: item.dataset_id))
        evidence = tuple(_reference_stability_evidence(item) for item in references)

        site_occurrences: Counter[str] = Counter()
        valid_occurrences: Counter[str] = Counter()
        for reference in references:
            intensities = reference.intensities
            site_occurrences.update(str(value) for value in intensities.index)
        for item in evidence:
            valid_occurrences.update(item.scores_by_site.keys())

        overlapping_sites = {
            site_key
            for site_key, count in site_occurrences.items()
            if count >= config.minimum_datasets_per_site
        }
        rankable_sites = tuple(
            sorted(
                site_key
                for site_key in overlapping_sites
                if valid_occurrences[site_key] >= config.minimum_datasets_per_site
            )
        )
        if len(rankable_sites) < config.minimum_shared_sites:
            raise SpsDiscoveryValidationError(
                SpsDiscoveryValidationResult(
                    issues=(
                        SpsDiscoveryValidationIssue(
                            code="insufficient_rankable_overlap",
                            message=(
                                "SPS discovery found "
                                f"{len(rankable_sites)} sites with valid condition-level "
                                "stability evidence in at least "
                                f"{config.minimum_datasets_per_site} reference datasets; "
                                f"at least {config.minimum_shared_sites} are required"
                            ),
                            field_name="reference_datasets",
                        ),
                    ),
                    reference_dataset_count=len(references),
                    potential_overlap_site_count=len(overlapping_sites),
                )
            )

        statistics_by_site: dict[str, list[SpsDatasetSiteStatistic]] = {
            site_key: [] for site_key in rankable_sites
        }
        quantiles_by_site: dict[str, list[float]] = {
            site_key: [] for site_key in rankable_sites
        }
        consensus_entries_by_reference: dict[str, int] = {}
        for item in evidence:
            scored_candidates = tuple(
                (site_key, item.scores_by_site[site_key])
                for site_key in rankable_sites
                if site_key in item.scores_by_site
            )
            ranks = _stable_ranks(scored_candidates)
            candidate_count = len(scored_candidates)
            consensus_entries_by_reference[item.dataset_id] = candidate_count
            for site_key, score in scored_candidates:
                competition_rank, average_rank = ranks[site_key]
                statistics_by_site[site_key].append(
                    SpsDatasetSiteStatistic(
                        dataset_id=item.dataset_id,
                        rank=competition_rank,
                        stability_score=score,
                    )
                )
                quantiles_by_site[site_key].append(
                    (float(candidate_count) - average_rank + 0.5)
                    / float(candidate_count)
                )

        consensus_scores = {
            site_key: _fisher_style_consensus_score(quantiles_by_site[site_key])
            for site_key in rankable_sites
        }
        ordered_sites = tuple(
            sorted(
                rankable_sites,
                key=lambda site_key: (-consensus_scores[site_key], site_key),
            )
        )
        ranking = tuple(
            SpsSiteStabilityRecord(
                site_key=site_key,
                consensus_rank=position,
                consensus_stability_score=consensus_scores[site_key],
                contributing_dataset_count=len(statistics_by_site[site_key]),
                dataset_statistics=tuple(statistics_by_site[site_key]),
            )
            for position, site_key in enumerate(ordered_sites, start=1)
        )
        selected = ordered_sites[: config.top_n]
        boundaries = SpsSelectionBoundaryCounts(
            total_unique_sites=len(site_occurrences),
            sites_meeting_dataset_overlap=len(overlapping_sites),
            sites_with_valid_stability=len(rankable_sites),
            sites_ranked=len(ranking),
            sites_selected=len(selected),
        )
        evidence_by_id = {item.dataset_id: item for item in evidence}
        provenance = SpsDiscoveryProvenance(
            source_datasets=tuple(
                _reference_provenance(
                    reference,
                    sites_passing_required_data=len(
                        evidence_by_id[reference.dataset_id].scores_by_site
                    ),
                    sites_entering_consensus=consensus_entries_by_reference[
                        reference.dataset_id
                    ],
                )
                for reference in references
            ),
            config=config,
            requested_control_count=config.top_n,
            actual_control_count=len(selected),
            selection_boundaries=boundaries,
        )
        return SpsDiscoveryResult(
            selected_site_keys=selected,
            site_ranking=ranking,
            provenance=provenance,
        )


def _reference_stability_evidence(
    reference: SpsReferenceDataset,
) -> _ReferenceStabilityEvidence:
    intensities = reference.intensities
    conditions = tuple(
        sorted({str(value) for value in reference.condition_by_sample.values()})
    )
    samples_by_condition = {
        condition: tuple(
            sorted(
                str(column)
                for column in intensities.columns
                if str(reference.condition_by_sample[str(column)]) == condition
            )
        )
        for condition in conditions
    }
    scores: dict[str, float] = {}
    for site_key in sorted(str(value) for value in intensities.index):
        condition_means: list[float] = []
        for condition in conditions:
            values = tuple(
                float(cast(float | int, intensities.at[site_key, sample_id]))
                for sample_id in samples_by_condition[condition]
                if not pd.isna(intensities.at[site_key, sample_id])
            )
            if not values or not all(math.isfinite(value) for value in values):
                condition_means = []
                break
            condition_means.append(math.fsum(values) / float(len(values)))
        if condition_means:
            scores[site_key] = max(abs(value) for value in condition_means)
    return _ReferenceStabilityEvidence(
        dataset_id=reference.dataset_id,
        scores_by_site=scores,
    )


def _reference_provenance(
    reference: SpsReferenceDataset,
    *,
    sites_passing_required_data: int | None = None,
    sites_entering_consensus: int | None = None,
) -> SpsReferenceDatasetProvenance:
    intensities = reference.intensities
    intensity_fingerprint = _sps_reference_intensity_fingerprint(
        dataset_id=reference.dataset_id,
        intensities=intensities,
    )
    quantitative_state = _require_sps_quantitative_state(
        reference.intensity_scale_state,
        field_name=(
            f"sps_reference_dataset[{reference.dataset_id!r}].intensity_scale_state"
        ),
        intensity_fingerprint=intensity_fingerprint,
    )
    scale_establishment = quantitative_state.establishment_provenance
    meaning_establishment = quantitative_state.quantitative_meaning_provenance
    quantitative_meaning = quantitative_state.quantity
    if (
        scale_establishment is None
        or meaning_establishment is None
        or quantitative_meaning is None
    ):
        raise AssertionError("validated SPS quantitative state lost required evidence")
    assignments = tuple(
        SpsSampleConditionAssignment(
            sample_id=sample_id,
            condition=str(condition),
        )
        for sample_id, condition in sorted(reference.condition_by_sample.items())
    )
    return SpsReferenceDatasetProvenance(
        dataset_id=reference.dataset_id,
        source_name=reference.source_name,
        source_version=reference.source_version,
        source_uri=reference.source_uri,
        site_count=int(intensities.shape[0]),
        sample_count=int(intensities.shape[1]),
        sample_conditions=assignments,
        intensity_fingerprint=intensity_fingerprint,
        intensity_scale_kind=quantitative_state.kind,
        quantitative_meaning=quantitative_meaning,
        intensity_scale_establishment=scale_establishment,
        quantitative_meaning_establishment=meaning_establishment,
        sites_passing_required_data=sites_passing_required_data,
        sites_entering_consensus=sites_entering_consensus,
    )


def _stable_ranks(
    scored_sites: Sequence[tuple[str, float]],
) -> dict[str, tuple[int, float]]:
    """Return deterministic competition ranks and PhosR-compatible midranks."""

    ordered = sorted(scored_sites, key=lambda item: (item[1], item[0]))
    ranks: dict[str, tuple[int, float]] = {}
    start = 0
    while start < len(ordered):
        stop = start + 1
        while stop < len(ordered) and ordered[stop][1] == ordered[start][1]:
            stop += 1
        competition_rank = start + 1
        average_rank = (float(start + 1) + float(stop)) / 2.0
        for site_key, _ in ordered[start:stop]:
            ranks[site_key] = (competition_rank, average_rank)
        start = stop
    return ranks


def _fisher_style_consensus_score(quantiles: Sequence[float]) -> float:
    contribution_count = len(quantiles)
    statistic = -2.0 * math.fsum(math.log(value) for value in quantiles)
    degrees_of_freedom = 2 * (contribution_count - 1)
    return float(chi2.sf(statistic, degrees_of_freedom))


def _require_sps_quantitative_state(
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
    if not _is_accepted_sps_meaning_evidence(meaning_provenance):
        raise PhosPyInputError(
            f"{field_name} quantitative-meaning evidence must establish that "
            "reference/control centring occurred"
        )
    if meaning_provenance.target_quantity is not state.quantity:
        raise PhosPyInputError(
            f"{field_name} quantitative-meaning provenance must match its meaning"
        )
    if not _sps_scale_evidence_matches_matrix(
        establishment,
        intensity_fingerprint=intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name} scale establishment evidence must be bound to the "
            "submitted SPS reference intensity matrix"
        )
    if not _sps_meaning_evidence_matches_matrix(
        meaning_provenance,
        intensity_fingerprint=intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name} reference/control-centering evidence must be bound to "
            "the submitted SPS reference intensity matrix"
        )
    return state


def _require_sps_provenance_quantitative_contract(
    provenance: SpsReferenceDatasetProvenance,
) -> None:
    """Validate explicit, comparable SPS quantitative provenance fields."""

    field_name = "sps_source_provenance"
    if not isinstance(
        cast(object, provenance.intensity_scale_kind), IntensityScaleKind
    ):
        raise PhosPyInputError(
            f"{field_name}.intensity_scale_kind must be an IntensityScaleKind"
        )
    if provenance.intensity_scale_kind is not IntensityScaleKind.LOG2:
        raise PhosPyInputError(
            f"{field_name}.intensity_scale_kind must establish the log2 scale"
        )
    if not isinstance(
        cast(object, provenance.quantitative_meaning), QuantitativeMeaning
    ):
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning must be a QuantitativeMeaning"
        )
    if provenance.quantitative_meaning is not (
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE
    ):
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning must establish condition-relative "
            "contrast_log2_fold_change meaning"
        )
    scale_establishment = provenance.intensity_scale_establishment
    if not isinstance(
        cast(object, scale_establishment),
        IntensityScaleEstablishmentProvenance,
    ):
        raise PhosPyInputError(
            f"{field_name}.intensity_scale_establishment must be an "
            "IntensityScaleEstablishmentProvenance"
        )
    if (
        scale_establishment.scale != IntensityScaleKind.LOG2.value
        or scale_establishment.evidence_level is IntensityScaleEvidenceLevel.UNKNOWN
    ):
        raise PhosPyInputError(
            f"{field_name}.intensity_scale_establishment must carry known log2 "
            "establishment evidence"
        )
    meaning_establishment = provenance.quantitative_meaning_establishment
    if not isinstance(
        cast(object, meaning_establishment),
        QuantitativeMeaningTransitionProvenance,
    ):
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning_establishment must be a "
            "QuantitativeMeaningTransitionProvenance"
        )
    if not _is_accepted_sps_meaning_evidence(meaning_establishment):
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning_establishment must establish "
            "that reference/control centring occurred"
        )
    if meaning_establishment.target_quantity is not provenance.quantitative_meaning:
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning_establishment must match "
            "quantitative_meaning"
        )
    if not _sps_scale_evidence_matches_matrix(
        scale_establishment,
        intensity_fingerprint=provenance.intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name}.intensity_scale_establishment must be bound to "
            "intensity_fingerprint"
        )
    if not _sps_meaning_evidence_matches_matrix(
        meaning_establishment,
        intensity_fingerprint=provenance.intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning_establishment must be bound to "
            "intensity_fingerprint"
        )


def _sps_reference_intensity_fingerprint(
    *,
    dataset_id: str,
    intensities: pd.DataFrame,
) -> TableFingerprint:
    return fingerprint_table_normalized_axes(
        intensities,
        name=f"sps_reference.{dataset_id}.intensities",
    )


def _is_accepted_sps_meaning_evidence(
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


def _sps_scale_evidence_matches_matrix(
    provenance: IntensityScaleEstablishmentProvenance,
    *,
    intensity_fingerprint: TableFingerprint,
) -> bool:
    payload = provenance.parameters.get(SPS_REFERENCE_MATRIX_FINGERPRINT_PARAMETER)
    return _sps_fingerprint_payload_matches(
        payload,
        intensity_fingerprint=intensity_fingerprint,
    )


def _sps_meaning_evidence_matches_matrix(
    provenance: QuantitativeMeaningTransitionProvenance,
    *,
    intensity_fingerprint: TableFingerprint,
) -> bool:
    return _sps_fingerprint_payload_matches(
        provenance.output_table_fingerprint,
        intensity_fingerprint=intensity_fingerprint,
    )


def _sps_fingerprint_payload_matches(
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


def _require_site_key(value: object, *, field_name: str) -> str:
    text = _require_non_empty_string(value, field_name=field_name)
    if value != text:
        raise PhosPyInputError(
            f"{field_name} must be an exact canonical site_key without "
            "surrounding whitespace"
        )
    require_site_key_index(
        pd.Index([text], name="site_key"),
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    return text


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def _require_positive_int(value: object, *, field_name: str) -> int:
    return _require_int_at_least(value, minimum=1, field_name=field_name)


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    return _require_int_at_least(value, minimum=0, field_name=field_name)


def _require_int_at_least(value: object, *, minimum: int, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise PhosPyInputError(f"{field_name} must be an integer >= {minimum}")
    return value


def _require_finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    resolved = float(value)
    if not math.isfinite(resolved):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    return resolved


def _scale_establishment_from_payload(
    payload: Mapping[str, object],
) -> IntensityScaleEstablishmentProvenance:
    """Restore repository-authoritative typed scale-establishment evidence."""

    resolved = _require_payload_mapping(
        payload,
        field_name="sps_source_provenance.quantitative_state.scale_establishment",
    )
    prefix = "sps_source_provenance.quantitative_state.scale_establishment"
    try:
        mode = IntensityScaleEstablishmentMode(
            _require_payload_string(
                resolved.get("establishment_mode"),
                field_name=f"{prefix}.establishment_mode",
            )
        )
        source = IntensityScaleEstablishmentSource(
            _require_payload_string(
                resolved.get("establishment_source"),
                field_name=f"{prefix}.establishment_source",
            )
        )
        evidence_level = IntensityScaleEvidenceLevel(
            _require_payload_string(
                resolved.get("evidence_level"),
                field_name=f"{prefix}.evidence_level",
            )
        )
    except ValueError as exc:
        raise PhosPyInputError(
            f"{prefix} contains an unsupported typed policy value"
        ) from exc
    warnings = _require_payload_sequence(
        resolved.get("diagnostic_warnings"),
        field_name=f"{prefix}.diagnostic_warnings",
    )
    return IntensityScaleEstablishmentProvenance(
        scale=_require_payload_string(
            resolved.get("scale"), field_name=f"{prefix}.scale"
        ),
        mode=mode,
        source=source,
        evidence_level=evidence_level,
        transformer_name=_optional_payload_string(
            resolved.get("transformer_name"),
            field_name=f"{prefix}.transformer_name",
        ),
        input_declaration_source=_optional_payload_string(
            resolved.get("input_declaration_source"),
            field_name=f"{prefix}.input_declaration_source",
        ),
        parameters=_require_payload_mapping(
            resolved.get("parameters"), field_name=f"{prefix}.parameters"
        ),
        trace_id=_optional_payload_string(
            resolved.get("trace_id"), field_name=f"{prefix}.trace_id"
        ),
        diagnostic_warnings=tuple(
            _require_payload_string(
                item, field_name=f"{prefix}.diagnostic_warnings[{position}]"
            )
            for position, item in enumerate(warnings)
        ),
    )


def _require_payload_mapping(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be an object")
    resolved: dict[str, object] = {}
    for key, item in cast(Mapping[object, object], value).items():
        if not isinstance(key, str):
            raise PhosPyInputError(f"{field_name} keys must be strings")
        resolved[key] = item
    return resolved


def _require_payload_sequence(value: object, *, field_name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must be an array")
    return tuple(cast(Sequence[object], value))


def _require_payload_string(value: object, *, field_name: str) -> str:
    return _require_non_empty_string(value, field_name=field_name)


def _optional_payload_string(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_payload_string(value, field_name=field_name)


def _require_payload_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an integer")
    return value


def _optional_payload_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return _require_payload_int(value, field_name=field_name)


def _require_payload_float(value: object, *, field_name: str) -> float:
    return _require_finite_float(value, field_name=field_name)


__all__ = [
    "SPS_DISCOVERY_ALGORITHM_ID",
    "SPS_DISCOVERY_ALGORITHM_VERSION",
    "SPS_DISCOVERY_CONTROL_SOURCE_TYPE",
    "SPS_DISCOVERY_SELECTION_METHOD_CONSENSUS_STABILITY",
    "SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING",
    "SpsDatasetSiteStatistic",
    "SpsDiscoveryConfig",
    "SpsDiscoveryProvenance",
    "SpsDiscoveryResult",
    "SpsDiscoverySelectionMethod",
    "SpsDiscoveryTieHandling",
    "SpsDiscoveryValidationError",
    "SpsDiscoveryValidationIssue",
    "SpsDiscoveryValidationResult",
    "SpsReferenceDataset",
    "SpsReferenceDatasetProvenance",
    "SpsSampleConditionAssignment",
    "SpsSelectionBoundaryCounts",
    "SpsSiteStabilityRecord",
]
