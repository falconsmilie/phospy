"""Current and historical result reconstruction for SPS discovery."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.phosphosites.site_keys import decode_site_key
from phospy.preprocessing.control_sites._discovery_payload import (
    require_finite_float,
    require_non_empty_string,
    require_payload_float,
    require_payload_int,
    require_payload_mapping,
    require_payload_sequence,
    require_payload_string,
    require_positive_int,
    require_site_key,
)
from phospy.preprocessing.control_sites.contracts import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.preprocessing.control_sites.discovery_contracts import (
    SPS_DISCOVERY_CONTROL_SOURCE_TYPE,
    SPS_DISCOVERY_CROSS_COUNT_ORDERING_POLICY,
    SPS_DISCOVERY_IDENTITY_PREFIX,
    SPS_DISCOVERY_IDENTITY_SCHEMA,
    SPS_DISCOVERY_IDENTITY_SCHEMA_V1,
    SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING,
)
from phospy.preprocessing.control_sites.discovery_provenance import (
    SpsDiscoveryProvenance,
)
from phospy.preprocessing.control_sites.discovery_ranking import (
    fisher_style_consensus_score,
    stable_ranks,
)
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.models import (
    JsonValue,
)


@dataclass(frozen=True, slots=True)
class SpsDatasetSiteStatistic:
    """Reference-specific ranking evidence for one candidate site.

    ``rank`` is the deterministic competition rank used for reporting, while
    ``rank_quantile`` is the half-step quantile derived from the estimator
    midrank and supplied directly to the Fisher-style consensus calculation.
    They differ for tied stability scores.
    """

    dataset_id: str
    rank: int
    stability_score: float
    rank_quantile: float = field(kw_only=True)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dataset_id",
            require_non_empty_string(
                self.dataset_id, field_name="sps_site_statistic.dataset_id"
            ),
        )
        require_positive_int(self.rank, field_name="sps_site_statistic.rank")
        stability_score = require_finite_float(
            self.stability_score, field_name="sps_site_statistic.stability_score"
        )
        if stability_score < 0.0:
            raise PhosPyInputError(
                "sps_site_statistic.stability_score must be non-negative"
            )
        rank_quantile = require_finite_float(
            self.rank_quantile, field_name="sps_site_statistic.rank_quantile"
        )
        if not 0.0 < rank_quantile < 1.0:
            raise PhosPyInputError(
                "sps_site_statistic.rank_quantile must be strictly between 0 and 1"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "dataset_id": self.dataset_id,
            "rank": self.rank,
            "stability_score": self.stability_score,
            "rank_quantile": self.rank_quantile,
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsDatasetSiteStatistic:
        """Restore one validated reference-specific site statistic."""

        resolved = require_payload_mapping(payload, field_name="sps_site_statistic")
        return cls(
            dataset_id=require_payload_string(
                resolved.get("dataset_id"),
                field_name="sps_site_statistic.dataset_id",
            ),
            rank=require_payload_int(
                resolved.get("rank"), field_name="sps_site_statistic.rank"
            ),
            stability_score=require_payload_float(
                resolved.get("stability_score"),
                field_name="sps_site_statistic.stability_score",
            ),
            rank_quantile=require_payload_float(
                resolved.get("rank_quantile"),
                field_name="sps_site_statistic.rank_quantile",
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
            require_site_key(self.site_key, field_name="sps_site_record.site_key"),
        )
        require_positive_int(
            self.consensus_rank, field_name="sps_site_record.consensus_rank"
        )
        consensus_stability_score = require_finite_float(
            self.consensus_stability_score,
            field_name="sps_site_record.consensus_stability_score",
        )
        if not 0.0 <= consensus_stability_score <= 1.0:
            raise PhosPyInputError(
                "sps_site_record.consensus_stability_score must be between 0 and 1"
            )
        require_positive_int(
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

        resolved = require_payload_mapping(payload, field_name="sps_site_record")
        statistics = require_payload_sequence(
            resolved.get("dataset_statistics"),
            field_name="sps_site_record.dataset_statistics",
        )
        return cls(
            site_key=require_payload_string(
                resolved.get("site_key"), field_name="sps_site_record.site_key"
            ),
            consensus_rank=require_payload_int(
                resolved.get("consensus_rank"),
                field_name="sps_site_record.consensus_rank",
            ),
            consensus_stability_score=require_payload_float(
                resolved.get("consensus_stability_score"),
                field_name="sps_site_record.consensus_stability_score",
            ),
            contributing_dataset_count=require_payload_int(
                resolved.get("contributing_dataset_count"),
                field_name="sps_site_record.contributing_dataset_count",
            ),
            dataset_statistics=tuple(
                SpsDatasetSiteStatistic.from_payload(
                    require_payload_mapping(
                        item,
                        field_name=f"sps_site_record.dataset_statistics[{position}]",
                    )
                )
                for position, item in enumerate(statistics)
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
            require_site_key(
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
        for record in self.site_ranking:
            if (
                record.contributing_dataset_count
                < self.provenance.config.minimum_datasets_per_site
            ):
                raise PhosPyInputError(
                    "sps_discovery_result.site_ranking contains a site below "
                    "config.minimum_datasets_per_site"
                )
        count_first_ordering = self.provenance.cross_count_ordering_policy == (
            SPS_DISCOVERY_CROSS_COUNT_ORDERING_POLICY
        )
        for previous, current in zip(
            self.site_ranking,
            self.site_ranking[1:],
            strict=False,
        ):
            if count_first_ordering and (
                previous.contributing_dataset_count < current.contributing_dataset_count
            ):
                raise PhosPyInputError(
                    "sps_discovery_result.site_ranking contributing dataset counts "
                    "must be non-increasing"
                )
            comparable_consensus = (
                not count_first_ordering
                or previous.contributing_dataset_count
                == current.contributing_dataset_count
            )
            if (
                comparable_consensus
                and previous.consensus_stability_score
                < current.consensus_stability_score
            ):
                raise PhosPyInputError(
                    "sps_discovery_result.site_ranking consensus stability scores "
                    "must be non-increasing within an evidence-count stratum"
                )
        if self.provenance.config.tie_handling == (
            SPS_DISCOVERY_TIE_HANDLING_SITE_KEY_ASCENDING
        ):
            for previous, current in zip(
                self.site_ranking,
                self.site_ranking[1:],
                strict=False,
            ):
                same_ordering_stratum = (
                    not count_first_ordering
                    or previous.contributing_dataset_count
                    == current.contributing_dataset_count
                )
                if (
                    same_ordering_stratum
                    and previous.consensus_stability_score
                    == current.consensus_stability_score
                    and previous.site_key > current.site_key
                ):
                    raise PhosPyInputError(
                        "sps_discovery_result.site_ranking equal consensus scores "
                        "must use ascending site_key tie handling"
                    )
        ranked_site_keys = tuple(item.site_key for item in self.site_ranking)
        for position, site_key in enumerate(ranked_site_keys):
            ranked_organism = decode_site_key(
                site_key,
                field_name=(f"sps_discovery_result.site_ranking[{position}].site_key"),
                error_type=PhosPyInputError,
            ).organism
            if ranked_organism != self.provenance.organism:
                raise PhosPyInputError(
                    "sps_discovery_result.site_ranking"
                    f"[{position}].site_key organism {ranked_organism.value!r} "
                    "must match sps_discovery_result.provenance.organism "
                    f"{self.provenance.organism.value!r}"
                )
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
        statistics_by_source: dict[str, list[tuple[str, SpsDatasetSiteStatistic]]] = {
            source_id: [] for source_id in source_ids
        }
        for record in self.site_ranking:
            for statistic in record.dataset_statistics:
                statistics_by_source[statistic.dataset_id].append(
                    (record.site_key, statistic)
                )
        for dataset_id, source_statistics in statistics_by_source.items():
            expected_ranks = stable_ranks(
                tuple(
                    (site_key, statistic.stability_score)
                    for site_key, statistic in source_statistics
                )
            )
            entry_count = len(source_statistics)
            for site_key, statistic in source_statistics:
                expected_rank, expected_midrank = expected_ranks[site_key]
                if statistic.rank != expected_rank:
                    raise PhosPyInputError(
                        "sps_discovery_result reference-specific rank must match "
                        f"the stability-score competition rank for {dataset_id!r}"
                    )
                expected_quantile = (
                    float(entry_count) - expected_midrank + 0.5
                ) / float(entry_count)
                if not math.isclose(
                    statistic.rank_quantile,
                    expected_quantile,
                    rel_tol=1e-12,
                    abs_tol=1e-15,
                ):
                    raise PhosPyInputError(
                        "sps_discovery_result reference-specific rank_quantile must "
                        "match the tie-derived midrank quantile"
                    )
        for record in self.site_ranking:
            expected_consensus = fisher_style_consensus_score(
                tuple(
                    statistic.rank_quantile for statistic in record.dataset_statistics
                )
            )
            if not math.isclose(
                record.consensus_stability_score,
                expected_consensus,
                rel_tol=1e-12,
                abs_tol=1e-15,
            ):
                raise PhosPyInputError(
                    "sps_discovery_result consensus_stability_score must match "
                    "the serialized reference rank_quantiles"
                )

    def to_control_site_set(self) -> ControlSiteSet:
        """Return selected SPS controls through the existing downstream contract."""

        return ControlSiteSet.from_site_keys(
            self.selected_site_keys,
            source_metadata=ControlSiteSourceMetadata(
                source_type=SPS_DISCOVERY_CONTROL_SOURCE_TYPE,
                organism=self.provenance.organism.value,
                identifier_namespace="site_key",
                source_name=self.provenance.algorithm_id,
                source_version=self.provenance.algorithm_version,
                selection_method=self.provenance.config.selection_method,
                sps_discovery_identity=self.discovery_identity,
                metadata_missing_reason={
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

    @property
    def discovery_identity(self) -> str:
        """Return the canonical immutable identity of this discovery evidence."""

        digest = hash_json_payload(cast(JsonValue, self._identity_payload()))
        return f"{SPS_DISCOVERY_IDENTITY_PREFIX}{digest}"

    def _identity_payload(self) -> dict[str, object]:
        return {
            "identity_schema": SPS_DISCOVERY_IDENTITY_SCHEMA,
            "selected_site_keys": list(self.selected_site_keys),
            "site_ranking": [item.to_payload() for item in self.site_ranking],
            "provenance": self.provenance.to_payload(),
        }

    def _legacy_v1_identity_payload(self) -> dict[str, object]:
        return {
            "identity_schema": SPS_DISCOVERY_IDENTITY_SCHEMA_V1,
            "selected_site_keys": list(self.selected_site_keys),
            "site_ranking": [
                {
                    **{
                        key: value
                        for key, value in record.to_payload().items()
                        if key != "dataset_statistics"
                    },
                    "dataset_statistics": [
                        {
                            key: value
                            for key, value in statistic.to_payload().items()
                            if key != "rank_quantile"
                        }
                        for statistic in record.dataset_statistics
                    ],
                }
                for record in self.site_ranking
            ],
            "provenance": self.provenance.to_payload(),
        }

    def _legacy_v1_discovery_identity(self) -> str:
        digest = hash_json_payload(cast(JsonValue, self._legacy_v1_identity_payload()))
        return f"{SPS_DISCOVERY_IDENTITY_PREFIX}{digest}"

    def to_payload(self) -> dict[str, object]:
        return {
            **self._identity_payload(),
            "discovery_identity": self.discovery_identity,
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> SpsDiscoveryResult:
        """Restore a validated discovery result from its matrix-free payload."""

        resolved = require_payload_mapping(payload, field_name="sps_discovery_result")
        selected = require_payload_sequence(
            resolved.get("selected_site_keys"),
            field_name="sps_discovery_result.selected_site_keys",
        )
        ranking = require_payload_sequence(
            resolved.get("site_ranking"),
            field_name="sps_discovery_result.site_ranking",
        )
        serialized_schema = resolved.get("identity_schema")
        schema = (
            SPS_DISCOVERY_IDENTITY_SCHEMA_V1
            if serialized_schema is None
            else require_payload_string(
                serialized_schema,
                field_name="sps_discovery_result.identity_schema",
            )
        )
        if schema not in {
            SPS_DISCOVERY_IDENTITY_SCHEMA_V1,
            SPS_DISCOVERY_IDENTITY_SCHEMA,
        }:
            raise PhosPyInputError(
                "sps_discovery_result.identity_schema is not supported"
            )
        ranking_payloads = tuple(
            require_payload_mapping(
                item,
                field_name=f"sps_discovery_result.site_ranking[{position}]",
            )
            for position, item in enumerate(ranking)
        )
        if schema == SPS_DISCOVERY_IDENTITY_SCHEMA_V1:
            ranking_payloads = _migrate_sps_discovery_v1_ranking(ranking_payloads)
        result = cls(
            selected_site_keys=tuple(
                require_payload_string(
                    item,
                    field_name=f"sps_discovery_result.selected_site_keys[{position}]",
                )
                for position, item in enumerate(selected)
            ),
            site_ranking=tuple(
                SpsSiteStabilityRecord.from_payload(item) for item in ranking_payloads
            ),
            provenance=SpsDiscoveryProvenance.from_payload(
                require_payload_mapping(
                    resolved.get("provenance"),
                    field_name="sps_discovery_result.provenance",
                )
            ),
        )
        serialized_identity = resolved.get("discovery_identity")
        if serialized_identity is not None:
            restored_identity = require_payload_string(
                serialized_identity,
                field_name="sps_discovery_result.discovery_identity",
            )
            expected_identity = (
                result._legacy_v1_discovery_identity()
                if schema == SPS_DISCOVERY_IDENTITY_SCHEMA_V1
                else result.discovery_identity
            )
            if restored_identity != expected_identity:
                raise PhosPyInputError(
                    "sps_discovery_result.discovery_identity does not match the "
                    "canonical discovery evidence"
                )
        return result


def _migrate_sps_discovery_v1_ranking(
    ranking_payloads: Sequence[Mapping[str, object]],
) -> tuple[Mapping[str, object], ...]:
    """Upgrade baseline v1 ranking evidence with reconstructible quantiles."""

    scored_sites_by_source: dict[str, list[tuple[str, float]]] = {}
    for record_position, record in enumerate(ranking_payloads):
        site_key = require_payload_string(
            record.get("site_key"),
            field_name=f"sps_discovery_result.site_ranking[{record_position}].site_key",
        )
        statistics = require_payload_sequence(
            record.get("dataset_statistics"),
            field_name=(
                "sps_discovery_result.site_ranking"
                f"[{record_position}].dataset_statistics"
            ),
        )
        for statistic_position, statistic_payload in enumerate(statistics):
            field_name = (
                "sps_discovery_result.site_ranking"
                f"[{record_position}].dataset_statistics[{statistic_position}]"
            )
            statistic = require_payload_mapping(
                statistic_payload,
                field_name=field_name,
            )
            if "rank_quantile" in statistic:
                raise PhosPyInputError(
                    "phospy-sps-discovery-result-v1 dataset statistics must not "
                    "contain rank_quantile; use the v2 identity schema"
                )
            dataset_id = require_payload_string(
                statistic.get("dataset_id"),
                field_name=f"{field_name}.dataset_id",
            )
            stability_score = require_payload_float(
                statistic.get("stability_score"),
                field_name=f"{field_name}.stability_score",
            )
            scored_sites_by_source.setdefault(dataset_id, []).append(
                (site_key, stability_score)
            )

    ranks_by_source = {
        dataset_id: stable_ranks(scored_sites)
        for dataset_id, scored_sites in scored_sites_by_source.items()
    }
    migrated: list[Mapping[str, object]] = []
    for record_position, record in enumerate(ranking_payloads):
        site_key = require_payload_string(
            record.get("site_key"),
            field_name=f"sps_discovery_result.site_ranking[{record_position}].site_key",
        )
        statistics = require_payload_sequence(
            record.get("dataset_statistics"),
            field_name=(
                "sps_discovery_result.site_ranking"
                f"[{record_position}].dataset_statistics"
            ),
        )
        migrated_statistics: list[dict[str, object]] = []
        for statistic_position, statistic_payload in enumerate(statistics):
            field_name = (
                "sps_discovery_result.site_ranking"
                f"[{record_position}].dataset_statistics[{statistic_position}]"
            )
            statistic = require_payload_mapping(
                statistic_payload,
                field_name=field_name,
            )
            dataset_id = require_payload_string(
                statistic.get("dataset_id"),
                field_name=f"{field_name}.dataset_id",
            )
            entry_count = len(scored_sites_by_source[dataset_id])
            _, midrank = ranks_by_source[dataset_id][site_key]
            migrated_statistics.append(
                {
                    **statistic,
                    "rank_quantile": (float(entry_count) - midrank + 0.5)
                    / float(entry_count),
                }
            )
        migrated.append(
            {
                **record,
                "dataset_statistics": migrated_statistics,
            }
        )
    return tuple(migrated)


__all__ = [
    "SpsDatasetSiteStatistic",
    "SpsDiscoveryResult",
    "SpsSiteStabilityRecord",
]
