"""Typed contracts for reference-data SPS discovery.

This module deliberately defines contracts only.  SPS ranking mathematics is
owned by a later scientific implementation; constructing these values does not
select sites or change native SPS/RUV-style correction behaviour.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal, cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.immutability import freeze_json_mapping, thaw_json_mapping
from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization.tables import table_fingerprint_to_payload
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.science.sites.validation import require_site_key_index

SPS_DISCOVERY_SELECTION_METHOD_CONSENSUS_STABILITY = "consensus_stability"
SPS_DISCOVERY_ALGORITHM_ID = "phospy_sps_consensus_stability"
SPS_DISCOVERY_ALGORITHM_VERSION = "contract-v1"
SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING = "site_key_ascending"
SPS_DISCOVERY_CONTROL_SOURCE_TYPE = "sps_discovery"

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


@dataclass(frozen=True, slots=True, eq=False, init=False)
class SpsReferenceDataset:
    """One phosphoproteomics reference matrix and its sample groups.

    The matrix is defensively copied. Rows must ultimately be governed by its
    encoded ``site_key`` index; the private request validator owns that check.
    ``condition_by_sample`` is frozen recursively so caller mutation cannot
    alter a validated request.
    """

    dataset_id: str
    _owned_intensities: object = field(repr=False)
    condition_by_sample: Mapping[str, str]
    source_name: str | None = None
    source_version: str | None = None
    source_uri: str | None = None

    def __init__(
        self,
        dataset_id: str,
        intensities: pd.DataFrame,
        condition_by_sample: Mapping[str, str],
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

    def _to_provenance(self) -> SpsReferenceDatasetProvenance:
        """Build source provenance for the validated workflow boundary."""

        intensities = self._intensities_snapshot()
        if not isinstance(intensities, pd.DataFrame):
            raise PhosPyInputError(
                f"SPS reference dataset {self.dataset_id!r} intensities must be a "
                "pandas DataFrame"
            )
        assignments = tuple(
            SpsSampleConditionAssignment(
                sample_id=sample_id,
                condition=str(condition),
            )
            for sample_id, condition in sorted(self.condition_by_sample.items())
        )
        return SpsReferenceDatasetProvenance(
            dataset_id=self.dataset_id,
            source_name=self.source_name,
            source_version=self.source_version,
            source_uri=self.source_uri,
            site_count=int(intensities.shape[0]),
            sample_count=int(intensities.shape[1]),
            sample_conditions=assignments,
            intensity_fingerprint=fingerprint_table(
                intensities,
                name=f"sps_reference.{self.dataset_id}.intensities",
            ),
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


@dataclass(frozen=True, slots=True)
class SpsReferenceDatasetProvenance:
    """Immutable identity and fingerprint for one contributing reference."""

    dataset_id: str
    site_count: int
    sample_count: int
    sample_conditions: tuple[SpsSampleConditionAssignment, ...]
    intensity_fingerprint: TableFingerprint
    source_name: str | None = None
    source_version: str | None = None
    source_uri: str | None = None

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

    def to_payload(self) -> dict[str, object]:
        return {
            "dataset_id": self.dataset_id,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "source_uri": self.source_uri,
            "site_count": self.site_count,
            "sample_count": self.sample_count,
            "sample_conditions": [item.to_payload() for item in self.sample_conditions],
            "intensity_fingerprint": table_fingerprint_to_payload(
                self.intensity_fingerprint
            ),
        }


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
        _require_finite_float(
            self.stability_score, field_name="sps_site_statistic.stability_score"
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "dataset_id": self.dataset_id,
            "rank": self.rank,
            "stability_score": self.stability_score,
        }


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
        _require_finite_float(
            self.consensus_stability_score,
            field_name="sps_site_record.consensus_stability_score",
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
        source_ids = {item.dataset_id for item in self.provenance.source_datasets}
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
                        "SPS references may be independently sourced; organism "
                        "coherence is not inferred by the contract-only ticket"
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
