from __future__ import annotations

import warnings
from collections import Counter
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.science.prediction.motif_scoring.frequency_matrices import (
    _build_frequency_matrix_from_windows,
)
from phospy.science.prediction.motif_scoring.models import (
    _MOTIF_LIBRARY_VALIDATION_ATTR,
    AMINO_ACIDS,
    DEFAULT_MOTIF_FLANK_SIZE,
    ExplicitMotifSequence,
    MotifLibraryValidationResult,
    MotifLibraryValidationRow,
    _extract_sequence_window,
    _is_missing_scalar,
    _LibraryCandidate,
    _normalize_sequence_value,
)
from phospy.science.prediction.sequence_validation import (
    SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID,
    SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE,
    SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH,
    SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER,
    SEQUENCE_VALIDATION_STATUS_VALID,
    MotifSequenceValidator,
    SequenceValidationInput,
)

_BARE_MOTIF_SEQUENCE_DEPRECATION_MESSAGE = (
    "Bare motif sequence strings in motif_sequences are deprecated and will be "
    "rejected in a future release because they omit stable reference and site "
    "identity metadata needed for reproducible motif-library validation. Pass "
    "ExplicitMotifSequence values or mapping entries with reference_id, site_id, "
    "kinase, and sequence fields."
)


def build_motif_library(
    *,
    kinase_substrate_map: pd.DataFrame,
    site_sequences: pd.Series,
    site_identities: Mapping[str, str] | pd.Series | None = None,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """Build per-kinase motif frequency matrices from reference sequences."""

    if flank_size < 0:
        raise ValueError("flank_size must be >= 0")

    sequence_lookup = {
        str(site_id): sequence for site_id, sequence in site_sequences.items()
    }
    identity_lookup = (
        {}
        if site_identities is None
        else {
            str(site_id): str(identity) for site_id, identity in site_identities.items()
        }
    )
    candidates: list[_LibraryCandidate] = []
    for kinase, grouped in kinase_substrate_map.groupby("kinase", sort=False):
        sites = list(dict.fromkeys(grouped.loc[:, "substrate_site"].astype(str)))
        for reference_id in sites:
            candidates.append(
                _LibraryCandidate(
                    reference_id=reference_id,
                    site_id=identity_lookup.get(reference_id, reference_id),
                    kinase=str(kinase),
                    sequence_input=_normalize_sequence_value(
                        sequence_lookup.get(reference_id, np.nan)
                    ),
                )
            )

    frequency_matrices, size_series, validation = _build_motif_library_from_candidates(
        candidates,
        flank_size=flank_size,
    )
    size_series.attrs[_MOTIF_LIBRARY_VALIDATION_ATTR] = validation
    return frequency_matrices, size_series


def build_motif_library_from_sequences(
    *,
    motif_sequences: Mapping[
        str, Sequence[str | Mapping[str, object] | ExplicitMotifSequence]
    ],
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """Build motif frequency matrices from explicit per-kinase sequences.

    Prefer structured metadata carrying `reference_id`, optional `site_id`,
    optional `kinase`, and `sequence`. Bare sequence strings remain accepted
    during the deprecation window but emit `DeprecationWarning`.
    """

    if flank_size < 0:
        raise ValueError("flank_size must be >= 0")

    candidates: list[_LibraryCandidate] = []
    saw_bare_sequence = False
    for kinase, sequences in motif_sequences.items():
        kinase_name = str(kinase)
        for index, entry in enumerate(sequences):
            saw_bare_sequence = saw_bare_sequence or isinstance(entry, str)
            explicit_entry = _normalize_explicit_motif_sequence(
                entry=entry,
                default_kinase=kinase_name,
                default_reference_id=f"{kinase_name}#{index}",
            )
            candidates.append(
                _LibraryCandidate(
                    reference_id=explicit_entry.reference_id,
                    site_id=explicit_entry.site_id,
                    kinase=explicit_entry.kinase,
                    sequence_input=_normalize_sequence_value(explicit_entry.sequence),
                )
            )
    if saw_bare_sequence:
        warnings.warn(
            _BARE_MOTIF_SEQUENCE_DEPRECATION_MESSAGE,
            DeprecationWarning,
            stacklevel=2,
        )
    frequency_matrices, size_series, validation = _build_motif_library_from_candidates(
        candidates,
        flank_size=flank_size,
    )
    size_series.attrs[_MOTIF_LIBRARY_VALIDATION_ATTR] = validation
    return frequency_matrices, size_series


def get_motif_library_validation(
    motif_sizes: pd.Series,
) -> MotifLibraryValidationResult | None:
    """Extract attached motif-library validation diagnostics from motif sizes."""

    validation = motif_sizes.attrs.get(_MOTIF_LIBRARY_VALIDATION_ATTR)
    if isinstance(validation, MotifLibraryValidationResult):
        return validation
    return None


def _build_motif_library_from_candidates(
    candidates: Sequence[_LibraryCandidate],
    *,
    flank_size: int,
) -> tuple[dict[str, pd.DataFrame], pd.Series, MotifLibraryValidationResult]:
    expected_window_size = (2 * flank_size) + 1
    validator = MotifSequenceValidator(
        expected_window_size=expected_window_size,
        require_phospho_centre_residue=True,
        enforce_site_identity_format=True,
        window_length_policy="centred_superset",
        residue_validation_scope="centre_window",
    )
    validation_inputs = [
        SequenceValidationInput(
            site_id=candidate.reference_id,
            site_sequence=candidate.sequence_input,
            site_identity=candidate.site_id,
        )
        for candidate in candidates
    ]
    validation_result = validator.run(rows=validation_inputs)
    status_counts = Counter(row.status for row in validation_result.rows)

    accepted_windows_by_kinase: dict[str, list[str]] = {}
    validation_rows: list[MotifLibraryValidationRow] = []
    excluded_reference_ids: list[str] = []
    for candidate, row in zip(candidates, validation_result.rows, strict=True):
        validation_rows.append(
            MotifLibraryValidationRow(
                reference_id=candidate.reference_id,
                site_id=candidate.site_id,
                kinase=candidate.kinase,
                sequence=row.sequence,
                status=row.status,
                reason=row.reason,
                expected_centre_residue=row.expected_centre_residue,
                observed_centre_residue=row.observed_centre_residue,
                sequence_length=row.sequence_length,
            )
        )
        if row.status != SEQUENCE_VALIDATION_STATUS_VALID:
            excluded_reference_ids.append(candidate.reference_id)
            continue
        if row.sequence is None:
            continue
        accepted_window = _extract_sequence_window(row.sequence, flank_size)
        if _is_missing_scalar(accepted_window):
            excluded_reference_ids.append(candidate.reference_id)
            continue
        accepted_windows_by_kinase.setdefault(candidate.kinase, []).append(
            str(accepted_window)
        )

    frequency_matrices: dict[str, pd.DataFrame] = {}
    motif_sizes: dict[str, float] = {}
    for kinase, windows in accepted_windows_by_kinase.items():
        if not windows:
            continue
        windows_series = pd.Series(windows, dtype=object)
        frequency_matrices[kinase] = _build_frequency_matrix_from_windows(
            windows_series,
            context=f"kinase={kinase}",
        )
        motif_sizes[kinase] = float(len(windows_series))

    size_series = pd.Series(motif_sizes, dtype=float, name="motif_size")
    size_series.index.name = "kinase"
    validation = MotifLibraryValidationResult(
        total_reference_sequences=validation_result.total_sequences,
        accepted_reference_sequences=validation_result.valid_sequences,
        excluded_reference_sequences=validation_result.invalid_sequences,
        missing_sequences=status_counts[SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE],
        short_sequences=status_counts[SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE],
        off_centre_sequences=status_counts[
            SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE
        ],
        site_residue_mismatches=status_counts[
            SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH
        ],
        invalid_site_ids=status_counts[SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID],
        non_phospho_centre_residues=status_counts[
            SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE
        ],
        unsupported_residue_characters=status_counts[
            SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER
        ],
        sequences_excluded_from_motif_profile_construction=(
            validation_result.invalid_sequences
        ),
        expected_window_size=expected_window_size,
        supported_amino_acids=tuple(AMINO_ACIDS),
        accepted_window_length_policy=(
            f"input sequences must be centred and odd-length with minimum length "
            f"{expected_window_size}; scoring windows are centre-extracted to exactly "
            f"{expected_window_size} residues (2 * flank_size + 1)"
        ),
        unsupported_residue_policy=(
            "exclude any window containing unsupported amino acids; "
            f"supported residues: {', '.join(AMINO_ACIDS)}"
        ),
        excluded_reference_ids=tuple(excluded_reference_ids),
        rows=tuple(validation_rows),
    )
    return frequency_matrices, size_series, validation


def _normalize_explicit_motif_sequence(
    *,
    entry: str | Mapping[str, object] | ExplicitMotifSequence,
    default_kinase: str,
    default_reference_id: str,
) -> ExplicitMotifSequence:
    if isinstance(entry, ExplicitMotifSequence):
        reference_id = _coerce_identifier(
            entry.reference_id, fallback=default_reference_id
        )
        kinase = _coerce_identifier(entry.kinase, fallback=default_kinase)
        if kinase != default_kinase:
            raise ValueError(
                "explicit motif sequence kinase metadata must match the parent "
                f"mapping key '{default_kinase}'"
            )
        return ExplicitMotifSequence(
            reference_id=reference_id,
            site_id=_coerce_optional_identifier(entry.site_id),
            kinase=kinase,
            sequence=entry.sequence,
        )
    if isinstance(entry, Mapping):
        reference_id = _coerce_identifier(
            entry.get("reference_id"),
            fallback=default_reference_id,
        )
        entry_kinase = _coerce_identifier(entry.get("kinase"), fallback=default_kinase)
        if entry_kinase != default_kinase:
            raise ValueError(
                "explicit motif sequence kinase metadata must match the parent "
                f"mapping key '{default_kinase}'"
            )
        return ExplicitMotifSequence(
            reference_id=reference_id,
            site_id=_coerce_optional_identifier(entry.get("site_id")),
            kinase=entry_kinase,
            sequence=entry.get("sequence"),
        )
    return ExplicitMotifSequence(
        reference_id=default_reference_id,
        site_id=None,
        kinase=default_kinase,
        sequence=entry,
    )


def _coerce_identifier(value: object, *, fallback: str) -> str:
    identifier = _coerce_optional_identifier(value)
    if identifier is None:
        return str(fallback)
    return identifier


def _coerce_optional_identifier(value: object) -> str | None:
    if value is None:
        return None
    if bool(pd.Series((value,), dtype="object").isna().iat[0]):
        return None
    text = str(value).strip()
    if text == "":
        return None
    return text
