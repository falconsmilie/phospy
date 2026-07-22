"""Dataset-builder evidence resolution at the site-level boundary."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import Any, cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.models import PeptideEvidenceTable
from phospy.science.evidence.multi_site import (
    MULTI_SITE_POLICY_ERROR,
    MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    MULTI_SITE_POLICY_KEEP_JOINT,
    MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
    MultiSiteHandlingConfig,
    parse_phospho_site_tokens,
)
from phospy.science.sites.identifiers import parse_canonical_site_identifier

DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED = "site_level_resolved"
DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE = "peptide_evidence"
SUPPORTED_DATASET_SITE_RESOLUTION_MODES: tuple[str, ...] = (
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
)

DATASET_MULTI_SITE_POLICY_REJECT = "reject"
DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING = (
    "exclude_from_sequence_scoring"
)
DATASET_MULTI_SITE_POLICY_KEEP_JOINT = "keep_joint"
DATASET_MULTI_SITE_POLICY_SPLIT = "split"
SUPPORTED_DATASET_MULTI_SITE_POLICIES: tuple[str, ...] = (
    DATASET_MULTI_SITE_POLICY_REJECT,
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_MULTI_SITE_POLICY_SPLIT,
)

_POLICY_TO_MULTI_SITE_HANDLING_POLICY: dict[str, str] = {
    DATASET_MULTI_SITE_POLICY_REJECT: MULTI_SITE_POLICY_ERROR,
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING: (
        MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
    ),
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT: MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_MULTI_SITE_POLICY_SPLIT: MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
}

DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN = (
    "mapping_weighted_mean"
)
DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS = (
    "legacy_alias_for_arithmetic_mean_of_allocated_signals"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL = (
    "explicit_mapping_weight_when_supplied_else_equal_fraction_per_resolved_site"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT = "explicit_mapping_weight"
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL = (
    "derived_equal_weight_per_mapped_site"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW = (
    "sum_to_one_per_peptide_evidence_row"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE = (
    DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
)
DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION = (
    "multiply_peptide_signal_by_mapping_fraction"
)
DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS = (
    "arithmetic_mean_of_allocated_signals"
)
DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS = (
    "retain_duplicate_peptide_evidence_rows_as_separate_observations"
)
DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS = (
    DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
)
DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS = (
    "combine_ambiguous_and_unambiguous_allocated_signals_in_site_mean"
)
DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN = (
    DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
)
DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES = (
    "arithmetic_mean_of_finite_reported_localisation_values"
)
DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR = "validate_without_repair"
_LEGACY_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE = "sum_to_one_per_peptide_row"
_LEGACY_DUPLICATE_PEPTIDE_POLICY_RETAIN_ALL_ROWS = (
    "retain_all_peptide_rows_as_independent_observations"
)
_LEGACY_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN = (
    "mixed_ambiguous_and_unambiguous_rows_share_same_weighted_mean_aggregation"
)
_MAPPING_FRACTION_COLUMN = "mapping_fraction"
_MAPPING_WEIGHT_SUM_TOLERANCE = 1e-6
_SITE_SEQUENCE_SOURCE_PROVIDED = "provided_site_sequence"
_SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT = "peptide_sequence_site_string"
_SITE_SEQUENCE_SOURCE_MISSING = "missing"


class _MappingWeightSourcePolicy(str, Enum):
    EXPLICIT_OR_DERIVED_EQUAL = (
        DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
    )


class _MappingWeightNormalizationPolicy(str, Enum):
    SUM_TO_ONE_PER_PEPTIDE_ROW = (
        DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
    )


class _SignalAllocationPolicy(str, Enum):
    MULTIPLY_BY_MAPPING_FRACTION = (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )


class _SiteSummarisationPolicy(str, Enum):
    ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS = (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )


class _DuplicateEvidencePolicy(str, Enum):
    RETAIN_DUPLICATE_ROWS = (
        DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
    )


class _MixedAmbiguityPolicy(str, Enum):
    COMBINE_ALLOCATED_SIGNALS = (
        DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
    )


class _LocalisationAggregationPolicy(str, Enum):
    ARITHMETIC_MEAN_OF_FINITE_VALUES = (
        DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES
    )


@dataclass(frozen=True, slots=True)
class _PeptideEvidenceResolutionPolicies:
    mapping_weight_source_policy: _MappingWeightSourcePolicy
    mapping_weight_normalization_policy: _MappingWeightNormalizationPolicy
    signal_allocation_policy: _SignalAllocationPolicy
    site_summarisation_policy: _SiteSummarisationPolicy
    duplicate_evidence_policy: _DuplicateEvidencePolicy
    mixed_ambiguity_policy: _MixedAmbiguityPolicy
    localisation_aggregation_policy: _LocalisationAggregationPolicy
    legacy_aggregation_policy: str
    aggregation_formula: str

    def to_payload(self) -> dict[str, str]:
        return {
            "mapping_weight_source_policy": (self.mapping_weight_source_policy.value),
            "mapping_weight_normalization_policy": (
                self.mapping_weight_normalization_policy.value
            ),
            "signal_allocation_policy": self.signal_allocation_policy.value,
            "site_summarisation_policy": self.site_summarisation_policy.value,
            "duplicate_evidence_policy": self.duplicate_evidence_policy.value,
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy.value,
            "localisation_aggregation_policy": (
                self.localisation_aggregation_policy.value
            ),
            "aggregation_policy": self.legacy_aggregation_policy,
            "aggregation_formula": self.aggregation_formula,
        }


_CURRENT_RESOLUTION_POLICIES = _PeptideEvidenceResolutionPolicies(
    mapping_weight_source_policy=(_MappingWeightSourcePolicy.EXPLICIT_OR_DERIVED_EQUAL),
    mapping_weight_normalization_policy=(
        _MappingWeightNormalizationPolicy.SUM_TO_ONE_PER_PEPTIDE_ROW
    ),
    signal_allocation_policy=_SignalAllocationPolicy.MULTIPLY_BY_MAPPING_FRACTION,
    site_summarisation_policy=(
        _SiteSummarisationPolicy.ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    ),
    duplicate_evidence_policy=_DuplicateEvidencePolicy.RETAIN_DUPLICATE_ROWS,
    mixed_ambiguity_policy=_MixedAmbiguityPolicy.COMBINE_ALLOCATED_SIGNALS,
    localisation_aggregation_policy=(
        _LocalisationAggregationPolicy.ARITHMETIC_MEAN_OF_FINITE_VALUES
    ),
    legacy_aggregation_policy=(DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS),
    aggregation_formula=(
        "a[p,s,j] = mapping_fraction[p,s] * peptide_signal[p,j]; "
        "site_signal[s,j] = arithmetic_mean(a[p,s,j] for retained peptide rows "
        "p mapped to s)"
    ),
)


@dataclass(frozen=True, slots=True)
class PeptideEvidenceResolutionSummary:
    """Structured summary for peptide-to-site resolution provenance."""

    input_mode: str
    multi_site_policy: str | None
    peptide_observations_received: int
    unique_site_ids_produced: int
    ambiguous_observations: int
    excluded_observations: int
    split_observations: int
    mapping_weight_source_policy: str
    mapping_weight_normalization_policy: str
    signal_allocation_policy: str
    site_summarisation_policy: str
    duplicate_evidence_policy: str
    mixed_ambiguity_policy: str
    localisation_aggregation_policy: str
    aggregation_policy: str
    aggregation_formula: str
    mapping_weight_source: str
    mapping_weight_normalisation: str
    duplicate_peptide_policy: str
    duplicate_peptide_rows: int
    site_sequence_column_present: bool
    provided_site_sequence_count: int
    accepted_site_sequence_count: int
    rejected_site_sequence_count: int
    provided_site_sequence_used_count: int
    peptide_context_derived_site_sequence_count: int
    missing_site_sequence_count: int
    site_sequence_policy: str

    def to_payload(self) -> dict[str, object]:
        return {
            "input_mode": self.input_mode,
            "multi_site_policy": self.multi_site_policy,
            "peptide_observations_received": int(self.peptide_observations_received),
            "unique_site_ids_produced": int(self.unique_site_ids_produced),
            "ambiguous_observations": int(self.ambiguous_observations),
            "excluded_observations": int(self.excluded_observations),
            "split_observations": int(self.split_observations),
            "mapping_weight_source_policy": self.mapping_weight_source_policy,
            "mapping_weight_normalization_policy": (
                self.mapping_weight_normalization_policy
            ),
            "signal_allocation_policy": self.signal_allocation_policy,
            "site_summarisation_policy": self.site_summarisation_policy,
            "duplicate_evidence_policy": self.duplicate_evidence_policy,
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy,
            "localisation_aggregation_policy": self.localisation_aggregation_policy,
            "aggregation_policy": self.aggregation_policy,
            "aggregation_formula": self.aggregation_formula,
            "mapping_weight_source": self.mapping_weight_source,
            "mapping_weight_normalisation": self.mapping_weight_normalisation,
            "duplicate_peptide_policy": self.duplicate_peptide_policy,
            "duplicate_peptide_rows": int(self.duplicate_peptide_rows),
            "site_sequence_column_present": bool(self.site_sequence_column_present),
            "provided_site_sequence_count": int(self.provided_site_sequence_count),
            "accepted_site_sequence_count": int(self.accepted_site_sequence_count),
            "rejected_site_sequence_count": int(self.rejected_site_sequence_count),
            "provided_site_sequence_used_count": int(
                self.provided_site_sequence_used_count
            ),
            "peptide_context_derived_site_sequence_count": int(
                self.peptide_context_derived_site_sequence_count
            ),
            "missing_site_sequence_count": int(self.missing_site_sequence_count),
            "site_sequence_policy": self.site_sequence_policy,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> PeptideEvidenceResolutionSummary:
        """Deserialize current or legacy peptide-resolution summary payloads."""

        normalized = normalise_peptide_evidence_resolution_summary_payload(payload)
        return cls(
            input_mode=str(normalized["input_mode"]),
            multi_site_policy=(
                None
                if normalized.get("multi_site_policy") is None
                else str(normalized["multi_site_policy"])
            ),
            peptide_observations_received=_payload_int(
                normalized["peptide_observations_received"]
            ),
            unique_site_ids_produced=_payload_int(
                normalized["unique_site_ids_produced"]
            ),
            ambiguous_observations=_payload_int(normalized["ambiguous_observations"]),
            excluded_observations=_payload_int(normalized["excluded_observations"]),
            split_observations=_payload_int(normalized["split_observations"]),
            mapping_weight_source_policy=str(
                normalized["mapping_weight_source_policy"]
            ),
            mapping_weight_normalization_policy=str(
                normalized["mapping_weight_normalization_policy"]
            ),
            signal_allocation_policy=str(normalized["signal_allocation_policy"]),
            site_summarisation_policy=str(normalized["site_summarisation_policy"]),
            duplicate_evidence_policy=str(normalized["duplicate_evidence_policy"]),
            mixed_ambiguity_policy=str(normalized["mixed_ambiguity_policy"]),
            localisation_aggregation_policy=str(
                normalized["localisation_aggregation_policy"]
            ),
            aggregation_policy=str(normalized["aggregation_policy"]),
            aggregation_formula=str(normalized["aggregation_formula"]),
            mapping_weight_source=str(normalized["mapping_weight_source"]),
            mapping_weight_normalisation=str(
                normalized["mapping_weight_normalisation"]
            ),
            duplicate_peptide_policy=str(normalized["duplicate_peptide_policy"]),
            duplicate_peptide_rows=_payload_int(normalized["duplicate_peptide_rows"]),
            site_sequence_column_present=bool(
                normalized["site_sequence_column_present"]
            ),
            provided_site_sequence_count=_payload_int(
                normalized["provided_site_sequence_count"]
            ),
            accepted_site_sequence_count=_payload_int(
                normalized["accepted_site_sequence_count"]
            ),
            rejected_site_sequence_count=_payload_int(
                normalized["rejected_site_sequence_count"]
            ),
            provided_site_sequence_used_count=_payload_int(
                normalized["provided_site_sequence_used_count"]
            ),
            peptide_context_derived_site_sequence_count=_payload_int(
                normalized["peptide_context_derived_site_sequence_count"]
            ),
            missing_site_sequence_count=_payload_int(
                normalized["missing_site_sequence_count"]
            ),
            site_sequence_policy=str(normalized["site_sequence_policy"]),
        )


def _payload_int(value: object) -> int:
    return int(cast(Any, value))


def normalise_peptide_evidence_resolution_summary_payload(
    payload: Mapping[str, object],
) -> dict[str, object]:
    """Return a current-policy peptide-resolution payload from old or new data."""

    normalized = dict(payload)
    current = _CURRENT_RESOLUTION_POLICIES.to_payload()

    normalized["mapping_weight_source_policy"] = str(
        normalized.get("mapping_weight_source_policy")
        or current["mapping_weight_source_policy"]
    )

    legacy_normalization = normalized.get("mapping_weight_normalisation")
    if "mapping_weight_normalization_policy" not in normalized:
        normalized["mapping_weight_normalization_policy"] = (
            current["mapping_weight_normalization_policy"]
            if legacy_normalization
            in (
                None,
                _LEGACY_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE,
                current["mapping_weight_normalization_policy"],
            )
            else str(legacy_normalization)
        )
    normalized["mapping_weight_normalisation"] = str(
        normalized.get("mapping_weight_normalisation")
        or normalized["mapping_weight_normalization_policy"]
    )
    if (
        normalized["mapping_weight_normalisation"]
        == _LEGACY_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE
    ):
        normalized["mapping_weight_normalisation"] = str(
            normalized["mapping_weight_normalization_policy"]
        )

    normalized["signal_allocation_policy"] = str(
        normalized.get("signal_allocation_policy")
        or current["signal_allocation_policy"]
    )
    normalized["site_summarisation_policy"] = str(
        normalized.get("site_summarisation_policy")
        or current["site_summarisation_policy"]
    )

    legacy_duplicate_policy = normalized.get("duplicate_peptide_policy")
    if "duplicate_evidence_policy" not in normalized:
        normalized["duplicate_evidence_policy"] = (
            current["duplicate_evidence_policy"]
            if legacy_duplicate_policy
            in (
                None,
                _LEGACY_DUPLICATE_PEPTIDE_POLICY_RETAIN_ALL_ROWS,
                current["duplicate_evidence_policy"],
            )
            else str(legacy_duplicate_policy)
        )
    normalized["duplicate_peptide_policy"] = str(
        normalized.get("duplicate_peptide_policy")
        or normalized["duplicate_evidence_policy"]
    )
    if (
        normalized["duplicate_peptide_policy"]
        == _LEGACY_DUPLICATE_PEPTIDE_POLICY_RETAIN_ALL_ROWS
    ):
        normalized["duplicate_peptide_policy"] = str(
            normalized["duplicate_evidence_policy"]
        )

    mixed_ambiguity_policy = normalized.get("mixed_ambiguity_policy")
    normalized["mixed_ambiguity_policy"] = (
        current["mixed_ambiguity_policy"]
        if mixed_ambiguity_policy
        in (
            None,
            _LEGACY_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN,
            current["mixed_ambiguity_policy"],
        )
        else str(mixed_ambiguity_policy)
    )
    normalized["localisation_aggregation_policy"] = str(
        normalized.get("localisation_aggregation_policy")
        or current["localisation_aggregation_policy"]
    )

    aggregation_policy = normalized.get("aggregation_policy")
    normalized["aggregation_policy"] = (
        current["aggregation_policy"]
        if aggregation_policy
        in (
            None,
            DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN,
            current["aggregation_policy"],
        )
        else str(aggregation_policy)
    )
    normalized["aggregation_formula"] = str(
        normalized.get("aggregation_formula") or current["aggregation_formula"]
    )
    return normalized


@dataclass(frozen=True, slots=True)
class PeptideEvidenceResolutionResult:
    """Site-level matrices produced from peptide-level evidence."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    summary: PeptideEvidenceResolutionSummary


class PeptideEvidenceDatasetResolver:
    """Resolve peptide-level evidence into site-level dataset-builder tables."""

    def run(
        self,
        *,
        evidence: PeptideEvidenceTable,
        multi_site_policy: str,
    ) -> PeptideEvidenceResolutionResult:
        if not isinstance(evidence, PeptideEvidenceTable):
            raise PhosPyInputError(
                "dataset peptide evidence resolution requires a PeptideEvidenceTable"
            )
        _validate_dataset_multi_site_policy(
            multi_site_policy,
            field_name="dataset build request multi_site_policy",
        )
        evidence_frame = evidence.to_dataframe()
        mapping = evidence.site_mapping.to_dataframe()
        site_sequence_column_present = "site_sequence" in evidence_frame.columns
        provided_site_sequence_count = (
            _count_non_empty_strings(evidence_frame.loc[:, "site_sequence"])
            if site_sequence_column_present
            else 0
        )
        peptide_observations_received = int(evidence_frame.shape[0])
        ambiguous_observations = int(
            evidence_frame.loc[:, "multi_site"].astype(bool).sum()
        )
        excluded_observations = (
            ambiguous_observations
            if multi_site_policy
            == DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
            else 0
        )
        split_observations = (
            ambiguous_observations
            if multi_site_policy == DATASET_MULTI_SITE_POLICY_SPLIT
            else 0
        )
        duplicate_peptide_rows = int(
            evidence_frame.loc[:, "peptide_sequence"]
            .astype(str)
            .duplicated(keep=False)
            .sum()
        )

        if mapping.empty:
            raise PhosPyInputError(
                "dataset build request peptide_evidence resolved to zero mapped "
                "site rows after applying multi_site_policy"
            )
        mapped_rows = _join_peptide_rows_to_site_mapping(
            evidence_frame=evidence_frame,
            mapping=mapping,
            sample_columns=evidence.sample_intensity_columns,
        )
        if mapped_rows.empty:
            raise PhosPyInputError(
                "dataset build request peptide_evidence resolved to zero mapped "
                "site rows after joining peptide evidence and site mapping"
            )
        resolved_mapping = _resolve_and_validate_mapping_fractions(
            mapped_rows=mapped_rows
        )
        allocated_rows = _allocate_peptide_signals_to_resolved_sites(
            mapped_rows=resolved_mapping.rows,
            sample_columns=evidence.sample_intensity_columns,
        )
        phospho = _summarise_allocated_site_signals(
            allocated_rows=allocated_rows,
            sample_columns=evidence.sample_intensity_columns,
        )
        site_metadata, site_sequence_summary = (
            _aggregate_site_metadata_and_localisation(
                mapped_rows=allocated_rows,
                site_ids=phospho.index,
                multi_site_policy=multi_site_policy,
            )
        )
        accepted_site_sequence_count = _count_non_empty_strings(
            site_metadata.loc[:, "site_sequence"]
        )
        summary = _build_resolution_summary(
            multi_site_policy=multi_site_policy,
            peptide_observations_received=peptide_observations_received,
            unique_site_ids_produced=int(phospho.shape[0]),
            ambiguous_observations=ambiguous_observations,
            excluded_observations=excluded_observations,
            split_observations=split_observations,
            mapping_weight_source=resolved_mapping.mapping_weight_source,
            duplicate_peptide_rows=duplicate_peptide_rows,
            site_sequence_column_present=site_sequence_column_present,
            provided_site_sequence_count=provided_site_sequence_count,
            accepted_site_sequence_count=accepted_site_sequence_count,
            site_sequence_summary=site_sequence_summary,
        )
        return PeptideEvidenceResolutionResult(
            phospho=phospho,
            site_metadata=site_metadata,
            summary=summary,
        )


def build_multi_site_handling_config_for_dataset_policy(
    *,
    multi_site_policy: str,
) -> MultiSiteHandlingConfig:
    """Translate dataset-builder multi-site policy to evidence config."""

    _validate_dataset_multi_site_policy(
        multi_site_policy,
        field_name="dataset build request multi_site_policy",
    )
    resolved_policy = _POLICY_TO_MULTI_SITE_HANDLING_POLICY[multi_site_policy]
    return MultiSiteHandlingConfig(
        statistical_modeling_policy=resolved_policy,
        kinase_sequence_scoring_policy=resolved_policy,
    )


@dataclass(frozen=True, slots=True)
class _ResolvedMappingFractions:
    rows: pd.DataFrame
    mapping_weight_source: str


def _join_peptide_rows_to_site_mapping(
    *,
    evidence_frame: pd.DataFrame,
    mapping: pd.DataFrame,
    sample_columns: tuple[str, ...],
) -> pd.DataFrame:
    peptide_fields = [
        "peptide_row_id",
        "protein_accession",
        "site_string",
        "peptide_sequence",
        "multi_site",
    ]
    if "site_sequence" in evidence_frame.columns:
        peptide_fields.append("site_sequence")
    if "localisation_confidence" in evidence_frame.columns:
        peptide_fields.append("localisation_confidence")
    peptide_rows = evidence_frame.loc[:, peptide_fields + list(sample_columns)].copy(
        deep=True
    )
    return mapping.merge(peptide_rows, how="inner", on="peptide_row_id")


def _resolve_and_validate_mapping_fractions(
    *,
    mapped_rows: pd.DataFrame,
) -> _ResolvedMappingFractions:
    resolved = mapped_rows.copy(deep=True)
    mapping_weight_source = DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT
    if "mapping_weight" not in resolved.columns:
        counts = resolved.groupby("peptide_row_id", sort=False).size().astype(float)
        resolved.loc[:, "mapping_weight"] = resolved.loc[:, "peptide_row_id"].map(
            lambda peptide_row_id: float(1.0 / counts.loc[peptide_row_id])
        )
        mapping_weight_source = DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL
    mapping_fractions = pd.to_numeric(
        resolved.loc[:, "mapping_weight"],
        errors="coerce",
    )
    if mapping_fractions.isna().any() or (mapping_fractions <= 0.0).any():
        raise PhosPyInputError(
            "dataset build request peptide_evidence site mapping contains "
            "non-positive or non-numeric mapping_weight values"
        )
    per_peptide_weight_sum = mapping_fractions.groupby(
        resolved.loc[:, "peptide_row_id"]
    ).sum()
    invalid_weight_rows = per_peptide_weight_sum.loc[
        (per_peptide_weight_sum - 1.0).abs() > _MAPPING_WEIGHT_SUM_TOLERANCE
    ]
    if not invalid_weight_rows.empty:
        preview = ", ".join(
            f"{str(peptide_row_id)!r}={float(total_weight):.6f}"
            for peptide_row_id, total_weight in invalid_weight_rows.iloc[:5].items()
        )
        suffix = "" if int(invalid_weight_rows.shape[0]) <= 5 else " ..."
        raise PhosPyInputError(
            "dataset build request peptide_evidence mapping_weight values must sum "
            "to 1.0 per peptide_row_id; invalid totals: "
            f"{preview}{suffix}"
        )
    resolved.loc[:, _MAPPING_FRACTION_COLUMN] = mapping_fractions.to_numpy(dtype=float)
    return _ResolvedMappingFractions(
        rows=resolved,
        mapping_weight_source=mapping_weight_source,
    )


def _allocate_peptide_signals_to_resolved_sites(
    *,
    mapped_rows: pd.DataFrame,
    sample_columns: tuple[str, ...],
) -> pd.DataFrame:
    allocated = mapped_rows.copy(deep=True)
    mapping_fractions = allocated.loc[:, _MAPPING_FRACTION_COLUMN].to_numpy(dtype=float)
    for sample_column in sample_columns:
        allocated.loc[:, sample_column] = (
            pd.to_numeric(allocated.loc[:, sample_column], errors="coerce")
            * mapping_fractions
        )
    return allocated


def _summarise_allocated_site_signals(
    *,
    allocated_rows: pd.DataFrame,
    sample_columns: tuple[str, ...],
) -> pd.DataFrame:
    matrix = (
        allocated_rows.groupby("site_id", sort=True)[list(sample_columns)]
        .mean(numeric_only=True)
        .astype(float)
    )
    matrix.index = pd.Index(matrix.index.astype(str), name="site_id")
    return matrix


def _build_resolution_summary(
    *,
    multi_site_policy: str,
    peptide_observations_received: int,
    unique_site_ids_produced: int,
    ambiguous_observations: int,
    excluded_observations: int,
    split_observations: int,
    mapping_weight_source: str,
    duplicate_peptide_rows: int,
    site_sequence_column_present: bool,
    provided_site_sequence_count: int,
    accepted_site_sequence_count: int,
    site_sequence_summary: Mapping[str, int],
) -> PeptideEvidenceResolutionSummary:
    policy_payload = _CURRENT_RESOLUTION_POLICIES.to_payload()
    return PeptideEvidenceResolutionSummary(
        input_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        multi_site_policy=multi_site_policy,
        peptide_observations_received=peptide_observations_received,
        unique_site_ids_produced=unique_site_ids_produced,
        ambiguous_observations=ambiguous_observations,
        excluded_observations=excluded_observations,
        split_observations=split_observations,
        mapping_weight_source_policy=policy_payload["mapping_weight_source_policy"],
        mapping_weight_normalization_policy=(
            policy_payload["mapping_weight_normalization_policy"]
        ),
        signal_allocation_policy=policy_payload["signal_allocation_policy"],
        site_summarisation_policy=policy_payload["site_summarisation_policy"],
        duplicate_evidence_policy=policy_payload["duplicate_evidence_policy"],
        mixed_ambiguity_policy=policy_payload["mixed_ambiguity_policy"],
        localisation_aggregation_policy=(
            policy_payload["localisation_aggregation_policy"]
        ),
        aggregation_policy=policy_payload["aggregation_policy"],
        aggregation_formula=policy_payload["aggregation_formula"],
        mapping_weight_source=mapping_weight_source,
        mapping_weight_normalisation=(
            policy_payload["mapping_weight_normalization_policy"]
        ),
        duplicate_peptide_policy=policy_payload["duplicate_evidence_policy"],
        duplicate_peptide_rows=duplicate_peptide_rows,
        site_sequence_column_present=site_sequence_column_present,
        provided_site_sequence_count=provided_site_sequence_count,
        accepted_site_sequence_count=accepted_site_sequence_count,
        rejected_site_sequence_count=int(
            site_sequence_summary["rejected_provided_context_count"]
        ),
        provided_site_sequence_used_count=int(
            site_sequence_summary["provided_site_sequence_used_count"]
        ),
        peptide_context_derived_site_sequence_count=int(
            site_sequence_summary["peptide_context_derived_site_sequence_count"]
        ),
        missing_site_sequence_count=int(
            site_sequence_summary["missing_site_sequence_count"]
        ),
        site_sequence_policy=(
            DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
        ),
    )


@dataclass(frozen=True, slots=True)
class _SiteSequenceResolution:
    site_sequence: str | None
    source: str
    rejected_provided_context_count: int = 0


@dataclass(frozen=True, slots=True)
class _InvalidProvidedSiteSequence:
    value: str
    reason: str


@dataclass(frozen=True, slots=True)
class _PeptideContextSequenceDerivation:
    site_sequence: str | None
    distinct_sequences: tuple[str, ...]


def _aggregate_site_metadata_and_localisation(
    *,
    mapped_rows: pd.DataFrame,
    site_ids: pd.Index,
    multi_site_policy: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    grouped = mapped_rows.groupby("site_id", sort=True)
    include_localisation_confidence = "localisation_confidence" in mapped_rows.columns
    site_rows: list[dict[str, object]] = []
    site_sequence_source_counts = {
        _SITE_SEQUENCE_SOURCE_PROVIDED: 0,
        _SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT: 0,
        _SITE_SEQUENCE_SOURCE_MISSING: 0,
    }
    rejected_provided_context_count = 0
    for site_id in site_ids.astype(str).tolist():
        group = grouped.get_group(site_id)
        gene_symbol, site = parse_canonical_site_identifier(
            site_id,
            field_name="dataset peptide evidence site_id",
            error_type=PhosPyInputError,
        )
        protein_accession = _single_non_empty_string_or_error(
            group.loc[:, "protein_accession"],
            field_name="protein_accession",
            site_id=site_id,
        )
        site_sequence_resolution = _resolve_site_sequence_for_resolved_site(
            group=group,
            site_id=site_id,
            resolved_site_token=site,
            multi_site_policy=multi_site_policy,
        )
        site_sequence_source_counts[site_sequence_resolution.source] += 1
        rejected_provided_context_count += (
            site_sequence_resolution.rejected_provided_context_count
        )
        site_rows.append(
            {
                "site_id": site_id,
                "gene_symbol": gene_symbol,
                "site": site,
                "site_sequence": site_sequence_resolution.site_sequence,
                "protein_accession": protein_accession,
                "protein_namespace": (
                    "protein_accession" if protein_accession is not None else None
                ),
                "protein_identifier": protein_accession,
            }
        )
        if include_localisation_confidence:
            site_rows[-1]["localisation_confidence"] = (
                _aggregate_localisation_confidence(
                    group.loc[:, "localisation_confidence"]
                )
            )
    site_metadata = pd.DataFrame(site_rows).set_index("site_id", drop=True)
    site_metadata.index = pd.Index(site_metadata.index.astype(str), name="site_id")
    summary = {
        "rejected_provided_context_count": int(rejected_provided_context_count),
        "provided_site_sequence_used_count": int(
            site_sequence_source_counts[_SITE_SEQUENCE_SOURCE_PROVIDED]
        ),
        "peptide_context_derived_site_sequence_count": int(
            site_sequence_source_counts[_SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT]
        ),
        "missing_site_sequence_count": int(
            site_sequence_source_counts[_SITE_SEQUENCE_SOURCE_MISSING]
        ),
    }
    return site_metadata, summary


def _resolve_site_sequence_for_resolved_site(
    *,
    group: pd.DataFrame,
    site_id: str,
    resolved_site_token: str,
    multi_site_policy: str,
) -> _SiteSequenceResolution:
    supplied_values = (
        _non_empty_strings(group.loc[:, "site_sequence"])
        if "site_sequence" in group.columns
        else []
    )
    valid_sequences: set[str] = set()
    invalid_sequences: list[_InvalidProvidedSiteSequence] = []
    for supplied_value in supplied_values:
        try:
            normalized = _normalize_site_sequence_for_resolved_site(
                site_id=site_id,
                site_sequence=supplied_value,
                resolved_site_token=resolved_site_token,
            )
        except PhosPyInputError as exc:
            invalid_sequences.append(
                _InvalidProvidedSiteSequence(
                    value=supplied_value,
                    reason=" ".join(str(exc).split()),
                )
            )
            continue
        if normalized is not None:
            valid_sequences.add(normalized)

    distinct_valid_sequences = tuple(sorted(valid_sequences))
    if len(distinct_valid_sequences) > 1:
        _raise_conflicting_supplied_site_sequences(
            site_id=site_id,
            distinct_sequences=distinct_valid_sequences,
        )
    if len(distinct_valid_sequences) == 1 and invalid_sequences:
        _raise_mixed_supplied_site_sequence_evidence(
            site_id=site_id,
            valid_sequence=distinct_valid_sequences[0],
            invalid_sequences=tuple(invalid_sequences),
        )
    if len(distinct_valid_sequences) == 1:
        return _SiteSequenceResolution(
            site_sequence=distinct_valid_sequences[0],
            source=_SITE_SEQUENCE_SOURCE_PROVIDED,
        )

    split_multisite_context = _is_split_multisite_context(
        group=group,
        multi_site_policy=multi_site_policy,
    )
    if supplied_values:
        if split_multisite_context:
            derived = _derive_site_sequence_from_peptide_context(
                group=group,
                site_id=site_id,
                resolved_site_token=resolved_site_token,
            )
            if derived.site_sequence is not None:
                return _SiteSequenceResolution(
                    site_sequence=derived.site_sequence,
                    source=_SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT,
                    rejected_provided_context_count=len(invalid_sequences),
                )
            _raise_invalid_supplied_site_sequences(
                site_id=site_id,
                invalid_sequences=tuple(invalid_sequences),
                derived_sequences=derived.distinct_sequences,
                derivation_allowed=True,
            )
        _raise_invalid_supplied_site_sequences(
            site_id=site_id,
            invalid_sequences=tuple(invalid_sequences),
            derived_sequences=(),
            derivation_allowed=False,
        )

    if split_multisite_context:
        derived = _derive_site_sequence_from_peptide_context(
            group=group,
            site_id=site_id,
            resolved_site_token=resolved_site_token,
        )
        if derived.site_sequence is not None:
            return _SiteSequenceResolution(
                site_sequence=derived.site_sequence,
                source=_SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT,
            )

    return _SiteSequenceResolution(
        site_sequence=None,
        source=_SITE_SEQUENCE_SOURCE_MISSING,
    )


def _is_split_multisite_context(
    *,
    group: pd.DataFrame,
    multi_site_policy: str,
) -> bool:
    if multi_site_policy != DATASET_MULTI_SITE_POLICY_SPLIT:
        return False
    if "multi_site" not in group.columns:
        return False
    return bool(group.loc[:, "multi_site"].astype(bool).any())


def _derive_site_sequence_from_peptide_context(
    *,
    group: pd.DataFrame,
    site_id: str,
    resolved_site_token: str,
) -> _PeptideContextSequenceDerivation:
    derived_sequences: set[str] = set()
    for _, row in group.iterrows():
        derived = _derive_site_sequence_from_peptide_row(
            row=row,
            site_id=site_id,
            resolved_site_token=resolved_site_token,
        )
        if derived is not None:
            derived_sequences.add(derived)
    distinct = tuple(sorted(derived_sequences))
    if len(distinct) != 1:
        return _PeptideContextSequenceDerivation(
            site_sequence=None,
            distinct_sequences=distinct,
        )
    return _PeptideContextSequenceDerivation(
        site_sequence=distinct[0],
        distinct_sequences=distinct,
    )


def _derive_site_sequence_from_peptide_row(
    *,
    row: pd.Series,
    site_id: str,
    resolved_site_token: str,
) -> str | None:
    peptide_sequence = _optional_row_string(row, "peptide_sequence")
    site_string = _optional_row_string(row, "site_string")
    if peptide_sequence is None or site_string is None:
        return None
    sequence = peptide_sequence.strip().upper()
    if not sequence or not sequence.isalpha():
        return None
    try:
        resolved_tokens = parse_phospho_site_tokens(
            resolved_site_token,
            field_name="dataset peptide evidence resolved site token",
        )
        declared_tokens = parse_phospho_site_tokens(
            site_string,
            field_name="dataset peptide evidence site_string",
        )
    except PhosPyInputError:
        return None
    if len(resolved_tokens) != 1:
        return None
    resolved_token = resolved_tokens[0]
    possible_starts: set[int] | None = None
    for token in declared_tokens:
        token_starts = {
            int(token.position) - peptide_position + 1
            for peptide_position, residue in enumerate(sequence, start=1)
            if residue == token.residue
        }
        if not token_starts:
            return None
        possible_starts = (
            token_starts
            if possible_starts is None
            else possible_starts.intersection(token_starts)
        )
    if possible_starts is None or len(possible_starts) != 1:
        return None
    protein_start = next(iter(possible_starts))
    peptide_positions = [
        peptide_position
        for peptide_position, residue in enumerate(sequence, start=1)
        if residue == resolved_token.residue
        and protein_start + peptide_position - 1 == int(resolved_token.position)
    ]
    if len(peptide_positions) != 1:
        return None
    peptide_position = peptide_positions[0]
    flank = min(peptide_position - 1, len(sequence) - peptide_position)
    start = peptide_position - flank - 1
    end = peptide_position + flank
    derived = sequence[start:end]
    return _normalize_site_sequence_for_resolved_site(
        site_id=site_id,
        site_sequence=derived,
        resolved_site_token=resolved_site_token,
    )


def _optional_row_string(row: pd.Series, column_name: str) -> str | None:
    if column_name not in row.index:
        return None
    value = row.loc[column_name]
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _aggregate_localisation_confidence(values: pd.Series) -> float | None:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.loc[numeric.notna()]
    if finite.empty:
        return None
    return float(finite.mean())


def _normalize_site_sequence_for_resolved_site(
    *,
    site_id: str,
    site_sequence: str | None,
    resolved_site_token: str,
) -> str | None:
    if site_sequence is None:
        return None
    sequence = site_sequence.strip().upper()
    expected_residue = resolved_site_token.strip().upper()[:1]
    if expected_residue not in {"S", "T", "Y"}:
        return sequence
    if len(sequence) < 3:
        return sequence
    if not sequence.isalpha() or (len(sequence) % 2 == 0):
        return sequence
    centre = len(sequence) // 2
    observed_residue = sequence[centre]
    if observed_residue == expected_residue:
        return sequence
    raise PhosPyInputError(
        "dataset peptide evidence site_sequence centre residue mismatch for "
        f"site_id={site_id!r}: expected={expected_residue!r} from resolved site "
        f"token {resolved_site_token!r}, observed={observed_residue!r}. Do not "
        "provide peptide-evidence site_sequence values that disagree with "
        "resolved site identity; remove the sequence to enable reference "
        "derivation or correct the upstream evidence."
    )


def _single_non_empty_string_or_error(
    values: pd.Series,
    *,
    field_name: str,
    site_id: str,
) -> str | None:
    distinct = tuple(dict.fromkeys(_non_empty_strings(values)))
    if len(distinct) <= 1:
        return distinct[0] if distinct else None
    preview = ", ".join(repr(value) for value in distinct[:5])
    suffix = "" if len(distinct) <= 5 else " ..."
    raise PhosPyInputError(
        f"{field_name} must contain at most one distinct non-empty value per "
        f"resolved site_id; site_id={site_id!r}, conflicting_values=["
        f"{preview}{suffix}]. Suggested fix: disambiguate peptide-site mapping "
        "or split rows before building."
    )


def _non_empty_strings(values: pd.Series) -> list[str]:
    tokens: list[str] = []
    for value in values.tolist():
        try:
            if bool(pd.isna(value)):
                continue
        except (TypeError, ValueError):
            pass
        if value is None:
            continue
        text = str(value).strip()
        if text:
            tokens.append(text)
    return tokens


def _count_non_empty_strings(values: pd.Series) -> int:
    count = 0
    for value in values.tolist():
        try:
            if bool(pd.isna(value)):
                continue
        except (TypeError, ValueError):
            pass
        if value is None:
            continue
        if str(value).strip():
            count += 1
    return count


def _raise_conflicting_supplied_site_sequences(
    *,
    site_id: str,
    distinct_sequences: tuple[str, ...],
) -> None:
    preview = _preview_quoted_values(distinct_sequences)
    raise PhosPyInputError(
        "dataset peptide evidence site_sequence values conflict for resolved "
        f"site_id={site_id!r}: distinct_normalized_value_count="
        f"{len(distinct_sequences)}, values=[{preview}]. PhosPy rejects "
        "conflicting supplied site-sequence contexts instead of selecting by "
        "row order, frequency, or lexical order. Correct the source evidence "
        "or choose an explicit upstream reference-resolution policy before "
        "dataset building."
    )


def _raise_mixed_supplied_site_sequence_evidence(
    *,
    site_id: str,
    valid_sequence: str,
    invalid_sequences: tuple[_InvalidProvidedSiteSequence, ...],
) -> None:
    invalid_preview = _preview_invalid_site_sequence_values(invalid_sequences)
    raise PhosPyInputError(
        "dataset peptide evidence site_sequence values are inconsistent for "
        f"resolved site_id={site_id!r}: valid_normalized_value="
        f"{valid_sequence!r}, invalid_supplied_values=[{invalid_preview}]. "
        "Mixed valid and invalid supplied evidence must not be silently "
        "reduced to the valid value. Correct the source evidence or choose an "
        "explicit upstream reference-resolution policy before dataset building."
    )


def _raise_invalid_supplied_site_sequences(
    *,
    site_id: str,
    invalid_sequences: tuple[_InvalidProvidedSiteSequence, ...],
    derived_sequences: tuple[str, ...],
    derivation_allowed: bool,
) -> None:
    invalid_preview = _preview_invalid_site_sequence_values(invalid_sequences)
    if derivation_allowed:
        derived_preview = _preview_quoted_values(derived_sequences)
        derivation_detail = (
            "peptide-context derivation did not establish exactly one fallback "
            f"sequence; derived_candidate_count={len(derived_sequences)}, "
            f"derived_candidates=[{derived_preview}]"
        )
    else:
        derivation_detail = (
            "peptide-context derivation is only available for split multi-site "
            "context under multi_site_policy='split'"
        )
    raise PhosPyInputError(
        "dataset peptide evidence site_sequence values are invalid for resolved "
        f"site_id={site_id!r}: invalid_supplied_values=[{invalid_preview}]. "
        f"{derivation_detail}. Correct the source evidence or choose an explicit "
        "upstream reference-resolution policy before dataset building."
    )


def _preview_quoted_values(values: tuple[str, ...]) -> str:
    preview = ", ".join(repr(value) for value in values[:5])
    suffix = "" if len(values) <= 5 else " ..."
    return f"{preview}{suffix}"


def _preview_invalid_site_sequence_values(
    invalid_sequences: tuple[_InvalidProvidedSiteSequence, ...],
) -> str:
    distinct = tuple(
        sorted(
            {
                (invalid_sequence.value, invalid_sequence.reason)
                for invalid_sequence in invalid_sequences
            }
        )
    )
    preview = ", ".join(
        f"value={value!r}, reason={reason!r}" for value, reason in distinct[:5]
    )
    suffix = "" if len(distinct) <= 5 else " ..."
    return f"{preview}{suffix}"


def _validate_dataset_multi_site_policy(policy: object, *, field_name: str) -> None:
    if (
        not isinstance(policy, str)
        or policy not in SUPPORTED_DATASET_MULTI_SITE_POLICIES
    ):
        supported = ", ".join(
            repr(value) for value in SUPPORTED_DATASET_MULTI_SITE_POLICIES
        )
        raise PhosPyInputError(f"{field_name} must be one of: {supported}")


__all__ = [
    "DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING",
    "DATASET_MULTI_SITE_POLICY_KEEP_JOINT",
    "DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS",
    "DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS",
    "DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN",
    "DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION",
    "DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR",
    "DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS",
    "DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS",
    "DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN",
    "DATASET_MULTI_SITE_POLICY_REJECT",
    "DATASET_MULTI_SITE_POLICY_SPLIT",
    "DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE",
    "DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED",
    "SUPPORTED_DATASET_MULTI_SITE_POLICIES",
    "SUPPORTED_DATASET_SITE_RESOLUTION_MODES",
    "PeptideEvidenceDatasetResolver",
    "PeptideEvidenceResolutionResult",
    "PeptideEvidenceResolutionSummary",
    "build_multi_site_handling_config_for_dataset_policy",
    "normalise_peptide_evidence_resolution_summary_payload",
]
