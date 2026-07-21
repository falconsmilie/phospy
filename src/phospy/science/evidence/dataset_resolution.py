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
DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR = "validate_without_repair"
_SITE_SEQUENCE_SOURCE_PROVIDED = "provided_site_sequence"
_SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT = "peptide_sequence_site_string"
_SITE_SEQUENCE_SOURCE_MISSING = "missing"


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
            "aggregation_policy": self.aggregation_policy,
            "aggregation_formula": self.aggregation_formula,
            "mapping_weight_source": self.mapping_weight_source,
            "mapping_weight_normalisation": self.mapping_weight_normalisation,
            "duplicate_peptide_policy": self.duplicate_peptide_policy,
            "duplicate_peptide_rows": int(self.duplicate_peptide_rows),
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy,
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
        site_metadata, site_sequence_summary = _build_site_metadata(
            mapped_rows=mapped_rows,
            site_ids=phospho.index,
            multi_site_policy=multi_site_policy,
        )
        accepted_site_sequence_count = _count_non_empty_strings(
            site_metadata.loc[:, "site_sequence"]
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


def _build_site_metadata(
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
            site_rows[-1]["localisation_confidence"] = _mean_localisation_confidence(
                group.loc[:, "localisation_confidence"]
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


def _mean_localisation_confidence(values: pd.Series) -> float | None:
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
    "DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN",
    "DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR",
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
