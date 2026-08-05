"""Site-sequence validation, derivation, and sequence-resolution diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.dataset_resolution.contracts import (
    DATASET_MULTI_SITE_POLICY_SPLIT,
    SITE_SEQUENCE_SOURCE_MISSING,
    SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT,
    SITE_SEQUENCE_SOURCE_PROVIDED,
)
from phospy.science.evidence.multi_site import parse_phospho_site_tokens


@dataclass(frozen=True, slots=True)
class SiteSequenceResolution:
    """Resolved site-sequence value and provenance for one resolved site."""

    site_sequence: str | None
    source: str
    rejected_provided_context_count: int = 0


@dataclass(frozen=True, slots=True)
class InvalidProvidedSiteSequence:
    """Invalid supplied site-sequence context plus deterministic error detail."""

    value: str
    reason: str


@dataclass(frozen=True, slots=True)
class PeptideContextSequenceDerivation:
    """Candidate site-sequence derivation from peptide sequence and site string."""

    site_sequence: str | None
    distinct_sequences: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionDiagnostics:
    """Sequence-resolution counts passed to summary assembly."""

    rejected_provided_context_count: int
    provided_site_sequence_used_count: int
    peptide_context_derived_site_sequence_count: int
    missing_site_sequence_count: int


def resolve_site_sequence_for_resolved_site(
    *,
    group: pd.DataFrame,
    site_id: str,
    resolved_site_token: str,
    multi_site_policy: str,
) -> SiteSequenceResolution:
    """Resolve the sequence context for one resolved site group."""

    supplied_values = (
        non_empty_strings(group.loc[:, "site_sequence"])
        if "site_sequence" in group.columns
        else []
    )
    valid_sequences: set[str] = set()
    invalid_sequences: list[InvalidProvidedSiteSequence] = []
    for supplied_value in supplied_values:
        try:
            normalized = normalize_site_sequence_for_resolved_site(
                site_id=site_id,
                site_sequence=supplied_value,
                resolved_site_token=resolved_site_token,
            )
        except PhosPyInputError as exc:
            invalid_sequences.append(
                InvalidProvidedSiteSequence(
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
        return SiteSequenceResolution(
            site_sequence=distinct_valid_sequences[0],
            source=SITE_SEQUENCE_SOURCE_PROVIDED,
        )

    split_multisite_context = is_split_multisite_context(
        group=group,
        multi_site_policy=multi_site_policy,
    )
    if supplied_values:
        if split_multisite_context:
            derived = derive_site_sequence_from_peptide_context(
                group=group,
                site_id=site_id,
                resolved_site_token=resolved_site_token,
            )
            if derived.site_sequence is not None:
                return SiteSequenceResolution(
                    site_sequence=derived.site_sequence,
                    source=SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT,
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
        derived = derive_site_sequence_from_peptide_context(
            group=group,
            site_id=site_id,
            resolved_site_token=resolved_site_token,
        )
        if derived.site_sequence is not None:
            return SiteSequenceResolution(
                site_sequence=derived.site_sequence,
                source=SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT,
            )

    return SiteSequenceResolution(
        site_sequence=None,
        source=SITE_SEQUENCE_SOURCE_MISSING,
    )


def is_split_multisite_context(
    *,
    group: pd.DataFrame,
    multi_site_policy: str,
) -> bool:
    if multi_site_policy != DATASET_MULTI_SITE_POLICY_SPLIT:
        return False
    if "multi_site" not in group.columns:
        return False
    return bool(group.loc[:, "multi_site"].astype(bool).any())


def derive_site_sequence_from_peptide_context(
    *,
    group: pd.DataFrame,
    site_id: str,
    resolved_site_token: str,
) -> PeptideContextSequenceDerivation:
    """Derive a resolved site sequence from peptide context, if unambiguous."""

    derived_sequences: set[str] = set()
    for _, row in group.iterrows():
        derived = derive_site_sequence_from_peptide_row(
            row=row,
            site_id=site_id,
            resolved_site_token=resolved_site_token,
        )
        if derived is not None:
            derived_sequences.add(derived)
    distinct = tuple(sorted(derived_sequences))
    if len(distinct) != 1:
        return PeptideContextSequenceDerivation(
            site_sequence=None,
            distinct_sequences=distinct,
        )
    return PeptideContextSequenceDerivation(
        site_sequence=distinct[0],
        distinct_sequences=distinct,
    )


def derive_site_sequence_from_peptide_row(
    *,
    row: pd.Series,
    site_id: str,
    resolved_site_token: str,
) -> str | None:
    peptide_sequence = optional_row_string(row, "peptide_sequence")
    site_string = optional_row_string(row, "site_string")
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
    return normalize_site_sequence_for_resolved_site(
        site_id=site_id,
        site_sequence=derived,
        resolved_site_token=resolved_site_token,
    )


def optional_row_string(row: pd.Series, column_name: str) -> str | None:
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


def normalize_site_sequence_for_resolved_site(
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


def non_empty_strings(values: pd.Series) -> list[str]:
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


def count_non_empty_strings(values: pd.Series) -> int:
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
    invalid_sequences: tuple[InvalidProvidedSiteSequence, ...],
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
    invalid_sequences: tuple[InvalidProvidedSiteSequence, ...],
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
    invalid_sequences: tuple[InvalidProvidedSiteSequence, ...],
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


__all__ = [
    "InvalidProvidedSiteSequence",
    "PeptideContextSequenceDerivation",
    "SiteSequenceResolution",
    "SiteSequenceResolutionDiagnostics",
    "count_non_empty_strings",
    "derive_site_sequence_from_peptide_context",
    "derive_site_sequence_from_peptide_row",
    "is_split_multisite_context",
    "non_empty_strings",
    "normalize_site_sequence_for_resolved_site",
    "optional_row_string",
    "resolve_site_sequence_for_resolved_site",
]
