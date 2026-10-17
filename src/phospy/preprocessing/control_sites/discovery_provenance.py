"""Matrix-free provenance contracts and codecs for SPS discovery."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.preprocessing.control_sites._discovery_payload import (
    optional_payload_int,
    optional_payload_string,
    require_non_empty_string,
    require_non_negative_int,
    require_payload_int,
    require_payload_mapping,
    require_payload_sequence,
    require_payload_string,
    require_positive_int,
)
from phospy.preprocessing.control_sites.discovery_contracts import (
    SPS_DISCOVERY_ALGORITHM_ID,
    SPS_DISCOVERY_ALGORITHM_VERSION,
    SPS_DISCOVERY_CROSS_COUNT_ORDERING_POLICY,
    SPS_DISCOVERY_LEGACY_ALGORITHM_VERSION,
    SPS_DISCOVERY_LEGACY_CROSS_COUNT_ORDERING_POLICY,
    SpsDiscoveryConfig,
    SpsReferenceDataset,
    is_accepted_sps_meaning_evidence,
    require_sps_quantitative_state,
    sps_meaning_evidence_matches_matrix,
    sps_reference_intensity_fingerprint,
    sps_scale_evidence_matches_matrix,
)
from phospy.provenance.models import (
    TableFingerprint,
)
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.quantitative.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentProvenance,
    IntensityScaleEstablishmentSource,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    QuantitativeMeaning,
    QuantitativeMeaningTransitionProvenance,
)
from phospy.taxonomy.organisms import Organism, normalize_optional_organism


@dataclass(frozen=True, slots=True)
class SpsSampleConditionAssignment:
    """Condition/group assignment for one reference sample."""

    sample_id: str
    condition: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "sample_id",
            require_non_empty_string(
                self.sample_id, field_name="sps_sample_condition.sample_id"
            ),
        )
        object.__setattr__(
            self,
            "condition",
            require_non_empty_string(
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

        resolved = require_payload_mapping(payload, field_name="sps_sample_condition")
        return cls(
            sample_id=require_payload_string(
                resolved.get("sample_id"),
                field_name="sps_sample_condition.sample_id",
            ),
            condition=require_payload_string(
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
    source_name: str
    source_version: str
    source_uri: str
    sites_passing_required_data: int | None = None
    sites_entering_consensus: int | None = None
    organism: Organism = field(kw_only=True)
    baseline_context: str = field(kw_only=True)
    reference_context: str = field(kw_only=True)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dataset_id",
            require_non_empty_string(
                self.dataset_id, field_name="sps_source_provenance.dataset_id"
            ),
        )
        require_non_negative_int(
            self.site_count, field_name="sps_source_provenance.site_count"
        )
        require_non_negative_int(
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
        normalized_organism = normalize_optional_organism(
            self.organism,
            field_name="sps_source_provenance.organism",
            error_type=PhosPyInputError,
        )
        if normalized_organism is None:
            raise PhosPyInputError(
                "sps_source_provenance.organism must be a supported organism"
            )
        object.__setattr__(self, "organism", normalized_organism)
        for field_name in (
            "baseline_context",
            "reference_context",
            "source_name",
            "source_version",
            "source_uri",
        ):
            object.__setattr__(
                self,
                field_name,
                require_non_empty_string(
                    getattr(self, field_name),
                    field_name=f"sps_source_provenance.{field_name}",
                ),
            )
        for field_name in (
            "sites_passing_required_data",
            "sites_entering_consensus",
        ):
            value = getattr(self, field_name)
            if value is not None:
                require_non_negative_int(
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
            "organism": self.organism.value,
            "baseline_context": self.baseline_context,
            "reference_context": self.reference_context,
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

        resolved = require_payload_mapping(payload, field_name="sps_source_provenance")
        quantitative_state = require_payload_mapping(
            resolved.get("quantitative_state"),
            field_name="sps_source_provenance.quantitative_state",
        )
        scale_establishment_payload = require_payload_mapping(
            quantitative_state.get("scale_establishment"),
            field_name=("sps_source_provenance.quantitative_state.scale_establishment"),
        )
        meaning_establishment_payload = require_payload_mapping(
            quantitative_state.get("quantitative_meaning_establishment"),
            field_name=(
                "sps_source_provenance.quantitative_state."
                "quantitative_meaning_establishment"
            ),
        )
        fingerprint_payload = require_payload_mapping(
            resolved.get("intensity_fingerprint"),
            field_name="sps_source_provenance.intensity_fingerprint",
        )
        sample_condition_payloads = require_payload_sequence(
            resolved.get("sample_conditions"),
            field_name="sps_source_provenance.sample_conditions",
        )
        try:
            intensity_scale_kind = IntensityScaleKind(
                require_payload_string(
                    quantitative_state.get("scale"),
                    field_name="sps_source_provenance.quantitative_state.scale",
                )
            )
            quantitative_meaning = QuantitativeMeaning(
                require_payload_string(
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
            organism = Organism(
                require_payload_string(
                    resolved.get("organism"),
                    field_name="sps_source_provenance.organism",
                )
            )
        except (TypeError, ValueError) as exc:
            raise PhosPyInputError(
                "sps_source_provenance contains invalid quantitative or "
                "fingerprint evidence"
            ) from exc
        return cls(
            dataset_id=require_payload_string(
                resolved.get("dataset_id"),
                field_name="sps_source_provenance.dataset_id",
            ),
            organism=organism,
            baseline_context=require_payload_string(
                resolved.get("baseline_context"),
                field_name="sps_source_provenance.baseline_context",
            ),
            reference_context=require_payload_string(
                resolved.get("reference_context"),
                field_name="sps_source_provenance.reference_context",
            ),
            source_name=require_payload_string(
                resolved.get("source_name"),
                field_name="sps_source_provenance.source_name",
            ),
            source_version=require_payload_string(
                resolved.get("source_version"),
                field_name="sps_source_provenance.source_version",
            ),
            source_uri=require_payload_string(
                resolved.get("source_uri"),
                field_name="sps_source_provenance.source_uri",
            ),
            site_count=require_payload_int(
                resolved.get("site_count"),
                field_name="sps_source_provenance.site_count",
            ),
            sample_count=require_payload_int(
                resolved.get("sample_count"),
                field_name="sps_source_provenance.sample_count",
            ),
            sites_passing_required_data=optional_payload_int(
                resolved.get("sites_passing_required_data"),
                field_name="sps_source_provenance.sites_passing_required_data",
            ),
            sites_entering_consensus=optional_payload_int(
                resolved.get("sites_entering_consensus"),
                field_name="sps_source_provenance.sites_entering_consensus",
            ),
            sample_conditions=tuple(
                SpsSampleConditionAssignment.from_payload(
                    require_payload_mapping(
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
            require_non_negative_int(
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

        resolved = require_payload_mapping(
            payload, field_name="sps_selection_boundaries"
        )
        return cls(
            total_unique_sites=require_payload_int(
                resolved.get("total_unique_sites"),
                field_name="sps_selection_boundaries.total_unique_sites",
            ),
            sites_meeting_dataset_overlap=require_payload_int(
                resolved.get("sites_meeting_dataset_overlap"),
                field_name=("sps_selection_boundaries.sites_meeting_dataset_overlap"),
            ),
            sites_with_valid_stability=require_payload_int(
                resolved.get("sites_with_valid_stability"),
                field_name="sps_selection_boundaries.sites_with_valid_stability",
            ),
            sites_ranked=require_payload_int(
                resolved.get("sites_ranked"),
                field_name="sps_selection_boundaries.sites_ranked",
            ),
            sites_selected=require_payload_int(
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
    organism: Organism = field(init=False)
    cross_count_ordering_policy: str = field(init=False)

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
        organisms = {item.organism for item in self.source_datasets}
        if len(organisms) != 1:
            raise PhosPyInputError(
                "sps_discovery_provenance.source_datasets must declare one "
                "coherent organism"
            )
        object.__setattr__(self, "organism", next(iter(organisms)))
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
        require_positive_int(
            self.requested_control_count,
            field_name="sps_discovery_provenance.requested_control_count",
        )
        require_non_negative_int(
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
            require_non_empty_string(
                self.algorithm_id,
                field_name="sps_discovery_provenance.algorithm_id",
            ),
        )
        object.__setattr__(
            self,
            "algorithm_version",
            require_non_empty_string(
                self.algorithm_version,
                field_name="sps_discovery_provenance.algorithm_version",
            ),
        )
        policies_by_version = {
            SPS_DISCOVERY_LEGACY_ALGORITHM_VERSION: (
                SPS_DISCOVERY_LEGACY_CROSS_COUNT_ORDERING_POLICY
            ),
            SPS_DISCOVERY_ALGORITHM_VERSION: SPS_DISCOVERY_CROSS_COUNT_ORDERING_POLICY,
        }
        try:
            cross_count_policy = policies_by_version[self.algorithm_version]
        except KeyError as exc:
            raise PhosPyInputError(
                "sps_discovery_provenance.algorithm_version must be one of "
                f"{sorted(policies_by_version)!r}"
            ) from exc
        object.__setattr__(self, "cross_count_ordering_policy", cross_count_policy)

    def to_payload(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "organism": self.organism.value,
            "source_datasets": [item.to_payload() for item in self.source_datasets],
            "algorithm_id": self.algorithm_id,
            "algorithm_version": self.algorithm_version,
            "selection_method": self.config.selection_method,
            "requested_control_count": self.requested_control_count,
            "actual_control_count": self.actual_control_count,
            "algorithm_parameters": self.config.to_payload(),
            "selection_boundaries": self.selection_boundaries.to_payload(),
        }
        if self.algorithm_version == SPS_DISCOVERY_ALGORITHM_VERSION:
            payload["cross_count_ordering_policy"] = self.cross_count_ordering_policy
        return payload

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsDiscoveryProvenance:
        """Restore validated SPS discovery provenance without source matrices."""

        resolved = require_payload_mapping(
            payload, field_name="sps_discovery_provenance"
        )
        source_payloads = require_payload_sequence(
            resolved.get("source_datasets"),
            field_name="sps_discovery_provenance.source_datasets",
        )
        config = SpsDiscoveryConfig.from_payload(
            require_payload_mapping(
                resolved.get("algorithm_parameters"),
                field_name="sps_discovery_provenance.algorithm_parameters",
            )
        )
        serialized_selection_method = require_payload_string(
            resolved.get("selection_method"),
            field_name="sps_discovery_provenance.selection_method",
        )
        if serialized_selection_method != config.selection_method:
            raise PhosPyInputError(
                "sps_discovery_provenance.selection_method must match "
                "algorithm_parameters.selection_method"
            )
        result = cls(
            source_datasets=tuple(
                SpsReferenceDatasetProvenance.from_payload(
                    require_payload_mapping(
                        item,
                        field_name=(
                            f"sps_discovery_provenance.source_datasets[{position}]"
                        ),
                    )
                )
                for position, item in enumerate(source_payloads)
            ),
            config=config,
            requested_control_count=require_payload_int(
                resolved.get("requested_control_count"),
                field_name="sps_discovery_provenance.requested_control_count",
            ),
            actual_control_count=require_payload_int(
                resolved.get("actual_control_count"),
                field_name="sps_discovery_provenance.actual_control_count",
            ),
            selection_boundaries=SpsSelectionBoundaryCounts.from_payload(
                require_payload_mapping(
                    resolved.get("selection_boundaries"),
                    field_name="sps_discovery_provenance.selection_boundaries",
                )
            ),
            algorithm_id=require_payload_string(
                resolved.get("algorithm_id"),
                field_name="sps_discovery_provenance.algorithm_id",
            ),
            algorithm_version=require_payload_string(
                resolved.get("algorithm_version"),
                field_name="sps_discovery_provenance.algorithm_version",
            ),
        )
        serialized_cross_count_policy = resolved.get("cross_count_ordering_policy")
        if result.algorithm_version == SPS_DISCOVERY_ALGORITHM_VERSION:
            policy = require_payload_string(
                serialized_cross_count_policy,
                field_name=("sps_discovery_provenance.cross_count_ordering_policy"),
            )
            if policy != result.cross_count_ordering_policy:
                raise PhosPyInputError(
                    "sps_discovery_provenance.cross_count_ordering_policy must "
                    "match algorithm_version"
                )
        elif serialized_cross_count_policy is not None:
            raise PhosPyInputError(
                "legacy SPS discovery provenance must not contain "
                "cross_count_ordering_policy"
            )
        serialized_organism = resolved.get("organism")
        if serialized_organism is not None:
            normalized_organism = normalize_optional_organism(
                serialized_organism,
                field_name="sps_discovery_provenance.organism",
                error_type=PhosPyInputError,
            )
            if normalized_organism is not result.organism:
                raise PhosPyInputError(
                    "sps_discovery_provenance.organism must match all source "
                    "dataset organisms"
                )
        return result


def reference_provenance(
    reference: SpsReferenceDataset,
    *,
    sites_passing_required_data: int | None = None,
    sites_entering_consensus: int | None = None,
) -> SpsReferenceDatasetProvenance:
    intensities = reference.intensities
    intensity_fingerprint = sps_reference_intensity_fingerprint(
        dataset_id=reference.dataset_id,
        intensities=intensities,
    )
    quantitative_state = require_sps_quantitative_state(
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
        organism=cast(Organism, reference.organism),
        baseline_context=cast(str, reference.baseline_context),
        reference_context=cast(str, reference.reference_context),
        source_name=cast(str, reference.source_name),
        source_version=cast(str, reference.source_version),
        source_uri=cast(str, reference.source_uri),
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
    if not is_accepted_sps_meaning_evidence(meaning_establishment):
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning_establishment must establish "
            "that reference/control centring occurred"
        )
    if meaning_establishment.target_quantity is not provenance.quantitative_meaning:
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning_establishment must match "
            "quantitative_meaning"
        )
    if not sps_scale_evidence_matches_matrix(
        scale_establishment,
        intensity_fingerprint=provenance.intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name}.intensity_scale_establishment must be bound to "
            "intensity_fingerprint"
        )
    if not sps_meaning_evidence_matches_matrix(
        meaning_establishment,
        intensity_fingerprint=provenance.intensity_fingerprint,
    ):
        raise PhosPyInputError(
            f"{field_name}.quantitative_meaning_establishment must be bound to "
            "intensity_fingerprint"
        )


def _scale_establishment_from_payload(
    payload: Mapping[str, object],
) -> IntensityScaleEstablishmentProvenance:
    """Restore repository-authoritative typed scale-establishment evidence."""

    resolved = require_payload_mapping(
        payload,
        field_name="sps_source_provenance.quantitative_state.scale_establishment",
    )
    prefix = "sps_source_provenance.quantitative_state.scale_establishment"
    try:
        mode = IntensityScaleEstablishmentMode(
            require_payload_string(
                resolved.get("establishment_mode"),
                field_name=f"{prefix}.establishment_mode",
            )
        )
        source = IntensityScaleEstablishmentSource(
            require_payload_string(
                resolved.get("establishment_source"),
                field_name=f"{prefix}.establishment_source",
            )
        )
        evidence_level = IntensityScaleEvidenceLevel(
            require_payload_string(
                resolved.get("evidence_level"),
                field_name=f"{prefix}.evidence_level",
            )
        )
    except ValueError as exc:
        raise PhosPyInputError(
            f"{prefix} contains an unsupported typed policy value"
        ) from exc
    warnings = require_payload_sequence(
        resolved.get("diagnostic_warnings"),
        field_name=f"{prefix}.diagnostic_warnings",
    )
    return IntensityScaleEstablishmentProvenance(
        scale=require_payload_string(
            resolved.get("scale"), field_name=f"{prefix}.scale"
        ),
        mode=mode,
        source=source,
        evidence_level=evidence_level,
        transformer_name=optional_payload_string(
            resolved.get("transformer_name"),
            field_name=f"{prefix}.transformer_name",
        ),
        input_declaration_source=optional_payload_string(
            resolved.get("input_declaration_source"),
            field_name=f"{prefix}.input_declaration_source",
        ),
        parameters=require_payload_mapping(
            resolved.get("parameters"), field_name=f"{prefix}.parameters"
        ),
        trace_id=optional_payload_string(
            resolved.get("trace_id"), field_name=f"{prefix}.trace_id"
        ),
        diagnostic_warnings=tuple(
            require_payload_string(
                item, field_name=f"{prefix}.diagnostic_warnings[{position}]"
            )
            for position, item in enumerate(warnings)
        ),
    )


__all__ = [
    "SpsDiscoveryProvenance",
    "SpsReferenceDatasetProvenance",
    "SpsSampleConditionAssignment",
    "SpsSelectionBoundaryCounts",
]
