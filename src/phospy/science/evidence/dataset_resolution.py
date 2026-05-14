"""Dataset-builder evidence resolution at the site-level boundary."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.models import PeptideEvidenceTable
from phospy.science.evidence.multi_site import (
    MULTI_SITE_POLICY_ERROR,
    MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    MULTI_SITE_POLICY_KEEP_JOINT,
    MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
    MultiSiteHandlingConfig,
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
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT = "explicit_mapping_weight"
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL = (
    "derived_equal_weight_per_mapped_site"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE = (
    "sum_to_one_per_peptide_row"
)
DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS = (
    "retain_all_peptide_rows_as_independent_observations"
)
DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN = (
    "mixed_ambiguous_and_unambiguous_rows_share_same_weighted_mean_aggregation"
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
    aggregation_policy: str
    aggregation_formula: str
    mapping_weight_source: str
    mapping_weight_normalisation: str
    duplicate_peptide_policy: str
    duplicate_peptide_rows: int
    mixed_ambiguity_policy: str

    def to_payload(self) -> dict[str, object]:
        return {
            "input_mode": self.input_mode,
            "multi_site_policy": self.multi_site_policy,
            "peptide_observations_received": int(self.peptide_observations_received),
            "unique_site_ids_produced": int(self.unique_site_ids_produced),
            "ambiguous_observations": int(self.ambiguous_observations),
            "excluded_observations": int(self.excluded_observations),
            "split_observations": int(self.split_observations),
            "aggregation_policy": self.aggregation_policy,
            "aggregation_formula": self.aggregation_formula,
            "mapping_weight_source": self.mapping_weight_source,
            "mapping_weight_normalisation": self.mapping_weight_normalisation,
            "duplicate_peptide_policy": self.duplicate_peptide_policy,
            "duplicate_peptide_rows": int(self.duplicate_peptide_rows),
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy,
        }


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
        mapped_rows, mapping_weight_source = _build_mapped_rows(
            evidence_frame=evidence_frame,
            mapping=mapping,
            sample_columns=evidence.sample_intensity_columns,
        )
        if mapped_rows.empty:
            raise PhosPyInputError(
                "dataset build request peptide_evidence resolved to zero mapped "
                "site rows after joining peptide evidence and site mapping"
            )
        phospho = _aggregate_site_matrix(
            mapped_rows=mapped_rows,
            sample_columns=evidence.sample_intensity_columns,
        )
        site_metadata = _build_site_metadata(
            mapped_rows=mapped_rows, site_ids=phospho.index
        )
        summary = PeptideEvidenceResolutionSummary(
            input_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            multi_site_policy=multi_site_policy,
            peptide_observations_received=peptide_observations_received,
            unique_site_ids_produced=int(phospho.shape[0]),
            ambiguous_observations=ambiguous_observations,
            excluded_observations=excluded_observations,
            split_observations=split_observations,
            aggregation_policy=(
                DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN
            ),
            aggregation_formula=(
                "site_intensity = mean(per_peptide_intensity * mapping_weight)"
            ),
            mapping_weight_source=mapping_weight_source,
            mapping_weight_normalisation=(
                DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE
            ),
            duplicate_peptide_policy=DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS,
            duplicate_peptide_rows=duplicate_peptide_rows,
            mixed_ambiguity_policy=(
                DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN
            ),
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


def _build_mapped_rows(
    *,
    evidence_frame: pd.DataFrame,
    mapping: pd.DataFrame,
    sample_columns: tuple[str, ...],
) -> tuple[pd.DataFrame, str]:
    peptide_fields = ["peptide_row_id", "protein_accession"]
    if "site_sequence" in evidence_frame.columns:
        peptide_fields.append("site_sequence")
    if "localisation_confidence" in evidence_frame.columns:
        peptide_fields.append("localisation_confidence")
    peptide_rows = evidence_frame.loc[:, peptide_fields + list(sample_columns)].copy(
        deep=True
    )
    merged = mapping.merge(peptide_rows, how="inner", on="peptide_row_id")
    if merged.empty:
        return merged, DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL
    mapping_weight_source = DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT
    if "mapping_weight" not in merged.columns:
        counts = merged.groupby("peptide_row_id", sort=False).size().astype(float)
        merged.loc[:, "mapping_weight"] = merged.loc[:, "peptide_row_id"].map(
            lambda peptide_row_id: float(1.0 / counts.loc[str(peptide_row_id)])
        )
        mapping_weight_source = DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL
    weights = pd.to_numeric(merged.loc[:, "mapping_weight"], errors="coerce")
    if weights.isna().any() or (weights <= 0.0).any():
        raise PhosPyInputError(
            "dataset build request peptide_evidence site mapping contains "
            "non-positive or non-numeric mapping_weight values"
        )
    per_peptide_weight_sum = weights.groupby(merged.loc[:, "peptide_row_id"]).sum()
    invalid_weight_rows = per_peptide_weight_sum.loc[
        (per_peptide_weight_sum - 1.0).abs() > 1e-6
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
    for sample_column in sample_columns:
        merged.loc[:, sample_column] = pd.to_numeric(
            merged.loc[:, sample_column], errors="coerce"
        ) * weights.to_numpy(dtype=float)
    return merged, mapping_weight_source


def _aggregate_site_matrix(
    *,
    mapped_rows: pd.DataFrame,
    sample_columns: tuple[str, ...],
) -> pd.DataFrame:
    matrix = (
        mapped_rows.groupby("site_id", sort=True)[list(sample_columns)]
        .mean(numeric_only=True)
        .astype(float)
    )
    matrix.index = pd.Index(matrix.index.astype(str), name="site_id")
    return matrix


def _build_site_metadata(
    *, mapped_rows: pd.DataFrame, site_ids: pd.Index
) -> pd.DataFrame:
    grouped = mapped_rows.groupby("site_id", sort=True)
    include_localisation_confidence = "localisation_confidence" in mapped_rows.columns
    site_rows: list[dict[str, object]] = []
    for site_id in site_ids.astype(str).tolist():
        group = grouped.get_group(site_id)
        gene_symbol, site = parse_canonical_site_identifier(
            site_id,
            field_name="dataset peptide evidence site_id",
            error_type=PhosPyInputError,
        )
        protein_id = _first_non_empty_string(group.loc[:, "protein_accession"])
        site_sequence = (
            _first_non_empty_string(group.loc[:, "site_sequence"])
            if "site_sequence" in group.columns
            else None
        )
        site_rows.append(
            {
                "site_id": site_id,
                "gene_symbol": gene_symbol,
                "site": site,
                "site_sequence": _site_sequence_for_resolved_site(
                    site_sequence=site_sequence,
                    resolved_site_token=site,
                ),
                "protein_id": protein_id,
            }
        )
        if include_localisation_confidence:
            site_rows[-1]["localisation_confidence"] = _mean_localisation_confidence(
                group.loc[:, "localisation_confidence"]
            )
    site_metadata = pd.DataFrame(site_rows).set_index("site_id", drop=True)
    site_metadata.index = pd.Index(site_metadata.index.astype(str), name="site_id")
    return site_metadata


def _mean_localisation_confidence(values: pd.Series) -> float | None:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.loc[numeric.notna()]
    if finite.empty:
        return None
    return float(finite.mean())


def _site_sequence_for_resolved_site(
    *,
    site_sequence: str | None,
    resolved_site_token: str,
) -> str | None:
    if site_sequence is None:
        return None
    sequence = site_sequence.strip().upper()
    if len(sequence) < 3:
        return sequence
    if not sequence.isalpha() or (len(sequence) % 2 == 0):
        return sequence
    expected_residue = resolved_site_token.strip().upper()[:1]
    if expected_residue not in {"S", "T", "Y"}:
        return sequence
    center = len(sequence) // 2
    observed_residue = sequence[center]
    if observed_residue == expected_residue:
        return sequence
    return f"{sequence[:center]}{expected_residue}{sequence[center + 1 :]}"


def _first_non_empty_string(values: pd.Series) -> str | None:
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
            return text
    return None


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
    "DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN",
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
]
