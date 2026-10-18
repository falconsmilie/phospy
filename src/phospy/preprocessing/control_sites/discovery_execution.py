"""Numerical execution for validated SPS discovery plans."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.preprocessing.control_sites.discovery_contracts import (
    SpsDiscoveryValidationError,
    SpsDiscoveryValidationIssue,
    SpsDiscoveryValidationResult,
    SpsReferenceDataset,
)
from phospy.preprocessing.control_sites.discovery_planning import SpsDiscoveryPlan
from phospy.preprocessing.control_sites.discovery_provenance import (
    SpsDiscoveryProvenance,
    SpsSelectionBoundaryCounts,
    reference_provenance,
)
from phospy.preprocessing.control_sites.discovery_ranking import (
    fisher_style_consensus_score,
    stable_ranks,
)
from phospy.preprocessing.control_sites.discovery_results import (
    SpsDatasetSiteStatistic,
    SpsDiscoveryResult,
    SpsSiteStabilityRecord,
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
    Fisher-style chi-square survival score. Within one evidence-count stratum,
    larger consensus scores are more stable. Across strata, sites supported by
    more references rank first; this conservative evidence-count policy keeps
    absence of reference evidence from acting as favorable evidence.

    Missing replicate cells are ignored only when their condition retains at
    least one finite observation.  A site with an entirely missing condition
    is not rankable in that reference.  It remains eligible only when it has
    valid evidence in ``minimum_datasets_per_site`` references.
    """

    def run(self, plan: SpsDiscoveryPlan) -> SpsDiscoveryResult:
        if not isinstance(  # pyright: ignore[reportUnnecessaryIsInstance] - The runtime guard rejects invalid callers despite the precise internal annotation.
            plan, SpsDiscoveryPlan
        ):
            raise TypeError("plan must be an SpsDiscoveryPlan")
        reference_datasets = plan.reference_datasets
        config = plan.config
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
            ranks = stable_ranks(scored_candidates)
            candidate_count = len(scored_candidates)
            consensus_entries_by_reference[item.dataset_id] = candidate_count
            for site_key, score in scored_candidates:
                competition_rank, average_rank = ranks[site_key]
                rank_quantile = (float(candidate_count) - average_rank + 0.5) / float(
                    candidate_count
                )
                statistics_by_site[site_key].append(
                    SpsDatasetSiteStatistic(
                        dataset_id=item.dataset_id,
                        rank=competition_rank,
                        stability_score=score,
                        rank_quantile=rank_quantile,
                    )
                )
                quantiles_by_site[site_key].append(rank_quantile)

        consensus_scores = {
            site_key: fisher_style_consensus_score(quantiles_by_site[site_key])
            for site_key in rankable_sites
        }
        ordered_sites = tuple(
            sorted(
                rankable_sites,
                key=lambda site_key: (
                    -len(statistics_by_site[site_key]),
                    -consensus_scores[site_key],
                    site_key,
                ),
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
                reference_provenance(
                    reference,
                    intensity_fingerprint=plan.intensity_fingerprints_by_dataset.get(
                        reference.dataset_id
                    ),
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
    sample_values_by_condition = {
        condition: tuple(
            intensities[sample_id].to_numpy(copy=False)
            for sample_id in samples_by_condition[condition]
        )
        for condition in conditions
    }
    site_rows = tuple(
        sorted(
            (str(site_key), position)
            for position, site_key in enumerate(intensities.index)
        )
    )
    scores: dict[str, float] = {}
    for site_key, row_position in site_rows:
        condition_means: list[float] = []
        for condition in conditions:
            values: list[float] = []
            for sample_values in sample_values_by_condition[condition]:
                value = sample_values[row_position]
                if not pd.isna(value):
                    values.append(float(cast(float | int, value)))
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


__all__: list[str] = []
