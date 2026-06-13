"""Protein-aware sample alignment and site eligibility diagnostics."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import pandas as pd

from phospy.contracts.configs.preprocessing.total_protein import (
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    DatasetProteinAwarePreparationMappingPolicy,
)
from phospy.policies import PolicyEnum
from phospy.science.datasets.preprocessing.protein_mapping import (
    ProteinMappingRecord,
    ProteinMappingResult,
    ProteinMappingStatus,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from phospy.validation.configs.preprocessing import (
    validate_protein_aware_sample_alignment_config,
)

PROTEIN_AWARE_REASON_MATCHED_PROTEIN_AVAILABLE = "matched_protein_available"
PROTEIN_AWARE_REASON_MISSING_SITE_PROTEIN_IDENTIFIER = "missing_site_protein_identifier"
PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW = "missing_total_protein_row"
PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING = "ambiguous_protein_mapping"
PROTEIN_AWARE_REASON_SAMPLE_MISMATCH = "sample_mismatch"
PROTEIN_AWARE_REASON_INCOMPATIBLE_TRANSFORMATION_STATE = (
    "incompatible_transformation_state"
)


class ProteinAwarePreparationEligibility(PolicyEnum):
    """Per-site protein-aware preparation eligibility status."""

    ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION = "eligible_for_protein_aware_preparation"
    FALLBACK_TO_PHOSPHO_ONLY = "fallback_to_phospho_only"
    EXCLUDED_FROM_PREPARATION = "excluded_from_preparation"


@dataclass(frozen=True, slots=True)
class ProteinAwareAlignmentConfig:
    """Configuration for protein-aware alignment diagnostics.

    `allow_reordered_samples=True` means reordered columns may be diagnosed as
    compatible. This collaborator still does not reorder matrices.
    """

    protein_mapping_policy: DatasetProteinAwarePreparationMappingPolicy = (
        DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS
    )
    allow_reordered_samples: bool = False

    def __post_init__(self) -> None:
        validate_protein_aware_sample_alignment_config(
            protein_mapping_policy=self.protein_mapping_policy,
            allow_reordered_samples=self.allow_reordered_samples,
            supported_mapping_policies=DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
        )
        object.__setattr__(
            self,
            "protein_mapping_policy",
            str(self.protein_mapping_policy).strip(),
        )


@dataclass(frozen=True, slots=True)
class ProteinAwareSampleAlignmentDiagnostics:
    """Diagnostics describing phospho/total sample-column compatibility."""

    phospho_sample_columns: tuple[str, ...]
    total_protein_sample_columns: tuple[str, ...]
    exact_sample_order_match: bool
    sample_order_compatible: bool
    reordered_sample_columns: bool
    allow_reordered_samples: bool
    missing_total_protein_samples: tuple[str, ...]
    extra_total_protein_samples: tuple[str, ...]

    def to_payload(self) -> dict[str, object]:
        return {
            "phospho_sample_columns": list(self.phospho_sample_columns),
            "total_protein_sample_columns": list(self.total_protein_sample_columns),
            "exact_sample_order_match": self.exact_sample_order_match,
            "sample_order_compatible": self.sample_order_compatible,
            "reordered_sample_columns": self.reordered_sample_columns,
            "allow_reordered_samples": self.allow_reordered_samples,
            "missing_total_protein_samples": list(self.missing_total_protein_samples),
            "extra_total_protein_samples": list(self.extra_total_protein_samples),
        }


@dataclass(frozen=True, slots=True)
class ProteinAwareTransformationStateDiagnostics:
    """Diagnostics describing phospho/total transformation-state compatibility."""

    compatible: bool
    phospho_transformation_state: dict[str, object] | None
    total_protein_transformation_state: dict[str, object] | None

    def to_payload(self) -> dict[str, object]:
        return {
            "transformation_state_compatible": self.compatible,
            "phospho_transformation_state": self.phospho_transformation_state,
            "total_protein_transformation_state": (
                self.total_protein_transformation_state
            ),
        }


@dataclass(frozen=True, slots=True)
class ProteinAwareSiteEligibilityDiagnostic:
    """Per-site protein-aware preparation eligibility diagnostic."""

    site_key: str
    eligibility: ProteinAwarePreparationEligibility
    reasons: tuple[str, ...]
    mapping_status: ProteinMappingStatus
    protein_identifier: str | None
    total_protein_row_key: str | None

    def to_payload(self) -> dict[str, object]:
        return {
            "site_key": self.site_key,
            "eligibility": self.eligibility.value,
            "reasons": list(self.reasons),
            "mapping_status": self.mapping_status.value,
            "protein_identifier": self.protein_identifier,
            "total_protein_row_key": self.total_protein_row_key,
        }


@dataclass(frozen=True, slots=True)
class ProteinAwareAlignmentEligibilityDiagnostics:
    """Complete protein-aware sample alignment and eligibility diagnostics."""

    sample_alignment: ProteinAwareSampleAlignmentDiagnostics
    transformation_state: ProteinAwareTransformationStateDiagnostics
    site_eligibility: tuple[ProteinAwareSiteEligibilityDiagnostic, ...]
    global_reasons: tuple[str, ...]

    @property
    def eligible_for_protein_aware_preparation(self) -> tuple[str, ...]:
        return _sites_with_status(
            self.site_eligibility,
            ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION,
        )

    @property
    def fallback_to_phospho_only(self) -> tuple[str, ...]:
        return _sites_with_status(
            self.site_eligibility,
            ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY,
        )

    @property
    def excluded_from_preparation(self) -> tuple[str, ...]:
        return _sites_with_status(
            self.site_eligibility,
            ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION,
        )

    @property
    def eligibility_by_site(self) -> dict[str, ProteinAwarePreparationEligibility]:
        return {
            diagnostic.site_key: diagnostic.eligibility
            for diagnostic in self.site_eligibility
        }

    @property
    def reasons_by_site(self) -> dict[str, tuple[str, ...]]:
        return {
            diagnostic.site_key: diagnostic.reasons
            for diagnostic in self.site_eligibility
        }

    @property
    def mapping_status_by_site(self) -> dict[str, ProteinMappingStatus]:
        return {
            diagnostic.site_key: diagnostic.mapping_status
            for diagnostic in self.site_eligibility
        }

    def to_payload(self) -> dict[str, object]:
        return {
            **self.sample_alignment.to_payload(),
            **self.transformation_state.to_payload(),
            "global_reasons": list(self.global_reasons),
            "eligible_for_protein_aware_preparation": list(
                self.eligible_for_protein_aware_preparation
            ),
            "fallback_to_phospho_only": list(self.fallback_to_phospho_only),
            "excluded_from_preparation": list(self.excluded_from_preparation),
            "site_eligibility": [
                diagnostic.to_payload() for diagnostic in self.site_eligibility
            ],
            "reasons_by_site": {
                site_key: list(reasons)
                for site_key, reasons in self.reasons_by_site.items()
            },
            "mapping_status_by_site": {
                site_key: status.value
                for site_key, status in self.mapping_status_by_site.items()
            },
        }


class ProteinAwareAlignmentEligibilityResolver:
    """Diagnose aligned protein-aware preparation eligibility.

    This collaborator reports readiness only. It does not run a model, subtract
    total protein, or decide downstream differential-analysis behavior.
    """

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None,
        mapping_result: ProteinMappingResult,
        intensity_scale_state: IntensityScaleState | None,
        config: ProteinAwareAlignmentConfig | None = None,
    ) -> ProteinAwareAlignmentEligibilityDiagnostics:
        resolved_config = config or ProteinAwareAlignmentConfig()
        sample_alignment = _diagnose_sample_alignment(
            phospho=phospho,
            total=total,
            allow_reordered_samples=resolved_config.allow_reordered_samples,
        )
        transformation_state = _diagnose_transformation_state(
            intensity_scale_state=intensity_scale_state,
            has_total_matrix=total is not None,
        )
        global_reasons = _resolve_global_reasons(
            sample_alignment=sample_alignment,
            transformation_state=transformation_state,
        )
        site_eligibility = tuple(
            _diagnose_site_eligibility(
                record=record,
                protein_mapping_policy=str(resolved_config.protein_mapping_policy),
                global_reasons=global_reasons,
            )
            for record in mapping_result.records
        )
        return ProteinAwareAlignmentEligibilityDiagnostics(
            sample_alignment=sample_alignment,
            transformation_state=transformation_state,
            site_eligibility=site_eligibility,
            global_reasons=global_reasons,
        )


def _diagnose_sample_alignment(
    *,
    phospho: pd.DataFrame,
    total: pd.DataFrame | None,
    allow_reordered_samples: bool,
) -> ProteinAwareSampleAlignmentDiagnostics:
    phospho_samples = _column_labels(phospho.columns)
    total_samples = () if total is None else _column_labels(total.columns)
    exact_match = phospho_samples == total_samples
    same_samples = Counter(phospho_samples) == Counter(total_samples)
    reordered = same_samples and not exact_match
    sample_order_compatible = exact_match or (
        bool(allow_reordered_samples) and reordered
    )

    total_sample_counts = Counter(total_samples)
    phospho_sample_counts = Counter(phospho_samples)
    missing_total_samples = tuple(
        sample
        for sample in phospho_samples
        if total_sample_counts.get(sample, 0) < phospho_sample_counts[sample]
    )
    extra_total_samples = tuple(
        sample
        for sample in total_samples
        if phospho_sample_counts.get(sample, 0) < total_sample_counts[sample]
    )
    if exact_match:
        missing_total_samples = ()
        extra_total_samples = ()

    return ProteinAwareSampleAlignmentDiagnostics(
        phospho_sample_columns=phospho_samples,
        total_protein_sample_columns=total_samples,
        exact_sample_order_match=exact_match,
        sample_order_compatible=sample_order_compatible,
        reordered_sample_columns=reordered,
        allow_reordered_samples=bool(allow_reordered_samples),
        missing_total_protein_samples=missing_total_samples,
        extra_total_protein_samples=extra_total_samples,
    )


def _diagnose_transformation_state(
    *,
    intensity_scale_state: IntensityScaleState | None,
    has_total_matrix: bool,
) -> ProteinAwareTransformationStateDiagnostics:
    if intensity_scale_state is None:
        return ProteinAwareTransformationStateDiagnostics(
            compatible=False,
            phospho_transformation_state=None,
            total_protein_transformation_state=None,
        )
    phospho_state = intensity_scale_state.phospho
    total_state = intensity_scale_state.total
    if not has_total_matrix:
        total_state = None
    compatible = (
        total_state is not None
        and phospho_state.kind is total_state.kind
        and phospho_state.transformed is total_state.transformed
    )
    return ProteinAwareTransformationStateDiagnostics(
        compatible=compatible,
        phospho_transformation_state=_matrix_state_payload(phospho_state),
        total_protein_transformation_state=(
            None if total_state is None else _matrix_state_payload(total_state)
        ),
    )


def _resolve_global_reasons(
    *,
    sample_alignment: ProteinAwareSampleAlignmentDiagnostics,
    transformation_state: ProteinAwareTransformationStateDiagnostics,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if not sample_alignment.sample_order_compatible:
        reasons.append(PROTEIN_AWARE_REASON_SAMPLE_MISMATCH)
    if not transformation_state.compatible:
        reasons.append(PROTEIN_AWARE_REASON_INCOMPATIBLE_TRANSFORMATION_STATE)
    return tuple(reasons)


def _diagnose_site_eligibility(
    *,
    record: ProteinMappingRecord,
    protein_mapping_policy: str,
    global_reasons: tuple[str, ...],
) -> ProteinAwareSiteEligibilityDiagnostic:
    mapping_reason = _reason_for_mapping_status(record.status)
    reasons = _dedupe_reasons((mapping_reason, *global_reasons))
    if global_reasons:
        eligibility = ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
    elif record.status is ProteinMappingStatus.MATCHED:
        eligibility = (
            ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION
        )
    elif record.status in {
        ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER,
        ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
    } and (
        protein_mapping_policy
        == DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT
    ):
        eligibility = ProteinAwarePreparationEligibility.FALLBACK_TO_PHOSPHO_ONLY
    else:
        eligibility = ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION
    return ProteinAwareSiteEligibilityDiagnostic(
        site_key=record.site_key,
        eligibility=eligibility,
        reasons=reasons,
        mapping_status=record.status,
        protein_identifier=record.protein_identifier,
        total_protein_row_key=record.total_protein_row_key,
    )


def _reason_for_mapping_status(status: ProteinMappingStatus) -> str:
    if status is ProteinMappingStatus.MATCHED:
        return PROTEIN_AWARE_REASON_MATCHED_PROTEIN_AVAILABLE
    if status is ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER:
        return PROTEIN_AWARE_REASON_MISSING_SITE_PROTEIN_IDENTIFIER
    if status is ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW:
        return PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW
    if status in {
        ProteinMappingStatus.AMBIGUOUS_SITE_PROTEIN_MAPPING,
        ProteinMappingStatus.AMBIGUOUS_TOTAL_PROTEIN_MAPPING,
    }:
        return PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING
    return str(status.value)


def _matrix_state_payload(state: MatrixIntensityScaleState) -> dict[str, object]:
    return {
        "kind": state.kind.value,
        "transformed": state.transformed,
        "established_by": state.established_by,
    }


def _column_labels(columns: pd.Index) -> tuple[str, ...]:
    return tuple(str(column) for column in columns.tolist())


def _sites_with_status(
    diagnostics: tuple[ProteinAwareSiteEligibilityDiagnostic, ...],
    status: ProteinAwarePreparationEligibility,
) -> tuple[str, ...]:
    return tuple(
        diagnostic.site_key
        for diagnostic in diagnostics
        if diagnostic.eligibility is status
    )


def _dedupe_reasons(reasons: tuple[str, ...]) -> tuple[str, ...]:
    seen: set[str] = set()
    deduped: list[str] = []
    for reason in reasons:
        if reason in seen:
            continue
        seen.add(reason)
        deduped.append(reason)
    return tuple(deduped)


__all__ = [
    "PROTEIN_AWARE_REASON_AMBIGUOUS_PROTEIN_MAPPING",
    "PROTEIN_AWARE_REASON_INCOMPATIBLE_TRANSFORMATION_STATE",
    "PROTEIN_AWARE_REASON_MATCHED_PROTEIN_AVAILABLE",
    "PROTEIN_AWARE_REASON_MISSING_SITE_PROTEIN_IDENTIFIER",
    "PROTEIN_AWARE_REASON_MISSING_TOTAL_PROTEIN_ROW",
    "PROTEIN_AWARE_REASON_SAMPLE_MISMATCH",
    "ProteinAwareAlignmentConfig",
    "ProteinAwareAlignmentEligibilityDiagnostics",
    "ProteinAwareAlignmentEligibilityResolver",
    "ProteinAwarePreparationEligibility",
    "ProteinAwareSampleAlignmentDiagnostics",
    "ProteinAwareSiteEligibilityDiagnostic",
    "ProteinAwareTransformationStateDiagnostics",
]
