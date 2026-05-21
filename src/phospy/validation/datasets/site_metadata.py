"""Phosphosite metadata validation and workflow policy enforcement."""

from __future__ import annotations

import math
import numbers
from dataclasses import dataclass
from typing import TypeVar

import pandas as pd

from phospy.contracts.configs.localisation import LocalisationRequirement
from phospy.science.sites.identifiers import ParsedSiteToken, try_parse_site_token
from phospy.science.sites.identity import (
    build_phosphosite_identity,
    validate_identity_optional_columns,
)

ErrorType = TypeVar("ErrorType", bound=Exception)
_PHOSPHORYLATABLE_RESIDUES = frozenset({"S", "T", "Y"})
_EXAMPLE_LIMIT = 5
_SITE_POSITION_CANDIDATE_COLUMNS = ("site_position", "position")
_SUPPORTED_GAP_SEQUENCE_CHARACTERS = frozenset({"_", "-"})


@dataclass(frozen=True, slots=True)
class LocalisationProbabilityAssessment:
    """Parsed localisation-probability assessment for one metadata column."""

    normalized: pd.Series
    missing_mask: pd.Series
    invalid_mask: pd.Series
    invalid_examples: tuple[str, ...]

    @property
    def missing_count(self) -> int:
        return int(self.missing_mask.sum())

    @property
    def invalid_count(self) -> int:
        return int(self.invalid_mask.sum())


def validate_site_identity_metadata(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    site_column: str = "site",
    site_sequence_column: str = "site_sequence",
    residue_column: str = "residue",
    allow_opaque_site_values: bool = False,
) -> None:
    """Validate row-level phosphosite identity coherence metadata."""

    validate_identity_optional_columns(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
    )

    malformed_site_values: list[str] = []
    inconsistent_residue_rows: list[str] = []
    inconsistent_position_rows: list[str] = []
    non_phospho_centre_rows: list[str] = []
    centre_mismatch_rows: list[str] = []
    invalid_explicit_residue_rows: list[str] = []
    invalid_explicit_position_rows: list[str] = []

    site_positions = _resolve_site_position_series(site_metadata)

    for site_id in site_metadata.index.tolist():
        site_value = site_metadata.at[site_id, site_column]
        parsed_site = try_parse_site_token(site_value)
        if parsed_site is None and not allow_opaque_site_values:
            malformed_site_values.append(f"{site_id!r}:{site_value!r}")

        explicit_residue = _resolve_optional_residue(
            site_metadata.at[site_id, residue_column]
            if residue_column in site_metadata.columns
            else None
        )
        if residue_column in site_metadata.columns and explicit_residue is None:
            raw_residue = site_metadata.at[site_id, residue_column]
            if not _is_missing(raw_residue):
                invalid_explicit_residue_rows.append(f"{site_id!r}:{raw_residue!r}")

        explicit_position = _resolve_optional_position(site_positions.at[site_id])
        if site_positions.name is not None and explicit_position is None:
            raw_position = site_positions.at[site_id]
            if not _is_missing(raw_position):
                invalid_explicit_position_rows.append(f"{site_id!r}:{raw_position!r}")

        if parsed_site is not None and explicit_residue is not None:
            if explicit_residue != parsed_site.residue:
                inconsistent_residue_rows.append(
                    f"{site_id!r}: site={parsed_site.residue!r}, "
                    f"residue_column={explicit_residue!r}"
                )
        if parsed_site is not None and explicit_position is not None:
            if explicit_position != parsed_site.position:
                inconsistent_position_rows.append(
                    f"{site_id!r}: site={parsed_site.position}, "
                    f"site_position_column={explicit_position}"
                )

        if site_sequence_column not in site_metadata.columns:
            continue
        parsed_sequence = _resolve_optional_sequence(
            site_metadata.at[site_id, site_sequence_column]
        )
        if parsed_sequence is None:
            continue
        if not _sequence_supports_central_residue_check(parsed_sequence):
            continue
        centre_residue = _resolve_central_residue(parsed_sequence)
        if centre_residue is None:
            continue
        if centre_residue not in _PHOSPHORYLATABLE_RESIDUES:
            non_phospho_centre_rows.append(
                f"{site_id!r}: centre={centre_residue!r}, sequence={parsed_sequence!r}"
            )
            continue
        expected_residue = _resolve_expected_residue(parsed_site, explicit_residue)
        if expected_residue is not None and centre_residue != expected_residue:
            centre_mismatch_rows.append(
                f"{site_id!r}: expected={expected_residue!r}, observed={centre_residue!r}"
            )

    details: list[str] = []
    if malformed_site_values:
        details.append(
            "site values must use strict 'S/T/Y<position>' tokens (example: "
            "'S123') unless opaque-site mode is explicitly enabled; "
            + _summarise_examples(malformed_site_values)
        )
    if invalid_explicit_residue_rows:
        details.append(
            f"{field_name}.residue must be one residue letter when provided; "
            + _summarise_examples(invalid_explicit_residue_rows)
        )
    if invalid_explicit_position_rows:
        details.append(
            f"{field_name}.{site_positions.name} must be an integer >= 1 when provided; "
            + _summarise_examples(invalid_explicit_position_rows)
        )
    if inconsistent_residue_rows:
        details.append(
            "residue column must match parsed site residue when both are present; "
            + _summarise_examples(inconsistent_residue_rows)
        )
    if inconsistent_position_rows:
        details.append(
            "site position column must match parsed site position when both are "
            "present; " + _summarise_examples(inconsistent_position_rows)
        )
    if non_phospho_centre_rows:
        details.append(
            f"{field_name}.site_sequence must contain a centred phosphorylatable "
            "residue (S/T/Y); " + _summarise_examples(non_phospho_centre_rows)
        )
    if centre_mismatch_rows:
        details.append(
            "site_sequence central residue must agree with site/residue metadata; "
            + _summarise_examples(centre_mismatch_rows)
        )

    if details:
        raise error_type(
            f"{field_name} phosphosite identity metadata validation failed; "
            + "; ".join(details)
        )


def validate_localisation_probability_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "localisation_probability",
) -> None:
    """Validate optional localisation probability values when the column exists."""

    assessment = assess_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=column_name,
    )
    if assessment is None or assessment.invalid_count == 0:
        return
    raise error_type(
        f"{field_name}.{column_name} must contain values in [0.0, 1.0] or missing; "
        f"invalid_row_count={assessment.invalid_count}; "
        f"examples={_summarise_examples(list(assessment.invalid_examples), limit=3)}"
    )


def enforce_localisation_requirement(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    workflow_name: str,
    requirement: LocalisationRequirement,
    error_type: type[ErrorType],
    column_name: str = "localisation_confidence",
) -> None:
    """Enforce workflow-level localisation policy using row-context diagnostics."""

    if not requirement.requires_probability_column:
        return
    resolved_column_name = _resolve_localisation_column_name(
        site_metadata=site_metadata,
        requested_column_name=column_name,
    )
    if resolved_column_name is None:
        site_examples = _site_id_examples(site_metadata.index)
        raise error_type(
            f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
            f"missing required column={field_name}.{column_name}; "
            f"affected_rows={int(site_metadata.shape[0])}; "
            f"example_site_ids={site_examples}"
        )

    assessment = assess_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=resolved_column_name,
    )
    if assessment is None:  # pragma: no cover - defensive guard
        return
    if assessment.invalid_count > 0:
        invalid_sites = _site_id_examples(
            site_metadata.index[assessment.invalid_mask.to_numpy()]
        )
        raise error_type(
            f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
            f"invalid values in {field_name}.{resolved_column_name}; "
            f"affected_rows={assessment.invalid_count}; "
            f"example_site_ids={invalid_sites}; "
            f"example_values={_summarise_examples(list(assessment.invalid_examples), limit=3)}"
        )
    if requirement.require_present and assessment.missing_count > 0:
        missing_sites = _site_id_examples(
            site_metadata.index[assessment.missing_mask.to_numpy()]
        )
        raise error_type(
            f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
            f"missing values in {field_name}.{resolved_column_name}; "
            f"affected_rows={assessment.missing_count}; "
            f"example_site_ids={missing_sites}"
        )
    if requirement.minimum_probability is None:
        return
    below_threshold = assessment.normalized.notna() & (
        assessment.normalized.astype("float64") < requirement.minimum_probability
    )
    below_threshold_count = int(below_threshold.sum())
    if below_threshold_count <= 0:
        return
    below_threshold_sites = _site_id_examples(
        site_metadata.index[below_threshold.to_numpy()]
    )
    threshold = float(requirement.minimum_probability)
    raise error_type(
        f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
        f"{field_name}.{resolved_column_name} must be >= {threshold:.3f}; "
        f"affected_rows={below_threshold_count}; "
        f"example_site_ids={below_threshold_sites}"
    )


def enforce_required_non_empty_string_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    workflow_name: str,
    column_name: str,
    error_type: type[ErrorType],
) -> None:
    """Require one site-metadata column to be non-missing, non-empty strings."""

    if column_name not in site_metadata.columns:
        site_examples = _site_id_examples(site_metadata.index)
        raise error_type(
            f"{field_name} is missing required columns: {column_name}; "
            f"{workflow_name} requires {field_name}.{column_name}; "
            f"missing required column={field_name}.{column_name}; "
            f"affected_rows={int(site_metadata.shape[0])}; "
            f"example_site_ids={site_examples}"
        )
    column = site_metadata.loc[:, column_name]
    invalid_mask = pd.Series(False, index=column.index.copy(), dtype="boolean")
    for site_id, raw_value in column.items():
        if _is_missing(raw_value):
            invalid_mask.at[site_id] = True
            continue
        if not isinstance(raw_value, str):
            invalid_mask.at[site_id] = True
            continue
        if raw_value.strip() == "":
            invalid_mask.at[site_id] = True
    invalid_count = int(invalid_mask.sum())
    if invalid_count == 0:
        return
    invalid_sites = _site_id_examples(site_metadata.index[invalid_mask.to_numpy()])
    raise error_type(
        f"{workflow_name} requires {field_name}.{column_name} to contain non-empty "
        f"string values; affected_rows={invalid_count}; "
        f"example_site_ids={invalid_sites}"
    )


def validate_site_sequence_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "site_sequence",
) -> None:
    """Validate site-sequence strings as plausible amino-acid contexts."""

    if column_name not in site_metadata.columns:
        return
    invalid_rows: list[str] = []
    values = site_metadata.loc[:, column_name]
    for site_id, raw_value in values.items():
        if not isinstance(raw_value, str):
            invalid_rows.append(f"{site_id!r}:{raw_value!r}")
            continue
        sequence = raw_value.strip().upper()
        if len(sequence) < 3:
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:sequence_too_short")
            continue
        if not any(character.isalpha() for character in sequence):
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:no_residue_letters")
            continue
        if any(
            (not character.isalpha()) and character != "_" for character in sequence
        ):
            invalid_rows.append(
                f"{site_id!r}:{raw_value!r}:unsupported_sequence_characters"
            )
    if not invalid_rows:
        return
    raise error_type(
        f"{field_name}.{column_name} must be plausible amino-acid context strings; "
        + _summarise_examples(invalid_rows)
    )


def enforce_centred_site_sequence_context(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    workflow_name: str,
    error_type: type[ErrorType],
    site_column: str = "site",
    site_sequence_column: str = "site_sequence",
    residue_column: str = "residue",
    allow_gapped_sequence_context: bool = False,
    allow_unknown_site_residue: bool = False,
) -> None:
    """Require centred residue context for sequence-aware workflow execution."""

    if site_sequence_column not in site_metadata.columns:
        site_examples = _site_id_examples(site_metadata.index)
        raise error_type(
            f"{workflow_name} requires centred sequence context in "
            f"{field_name}.{site_sequence_column}; "
            f"missing required column={field_name}.{site_sequence_column}; "
            f"affected_rows={int(site_metadata.shape[0])}; "
            f"example_site_ids={site_examples}"
        )

    missing_sequence_rows: list[str] = []
    even_length_rows: list[str] = []
    unsupported_character_rows: list[str] = []
    unknown_expected_residue_rows: list[str] = []
    non_phospho_centre_rows: list[str] = []
    centre_mismatch_rows: list[str] = []

    for site_id in site_metadata.index.tolist():
        raw_sequence = site_metadata.at[site_id, site_sequence_column]
        sequence = _resolve_optional_sequence(raw_sequence)
        if sequence is None:
            missing_sequence_rows.append(f"{site_id!r}:{raw_sequence!r}")
            continue

        sequence_length = len(sequence)
        if sequence_length % 2 == 0:
            even_length_rows.append(
                f"{site_id!r}:{sequence!r}:length={sequence_length}"
            )
            continue

        unsupported_characters = sorted(
            {character for character in sequence if not character.isalpha()}
        )
        if (
            allow_gapped_sequence_context
            and unsupported_characters
            and all(
                character in _SUPPORTED_GAP_SEQUENCE_CHARACTERS
                for character in unsupported_characters
            )
        ):
            unsupported_characters = []
        if unsupported_characters:
            joined = "".join(unsupported_characters)
            unsupported_character_rows.append(
                f"{site_id!r}:{sequence!r}:unsupported_characters={joined!r}"
            )
            continue

        centre_residue = sequence[sequence_length // 2]
        if centre_residue not in _PHOSPHORYLATABLE_RESIDUES:
            non_phospho_centre_rows.append(
                f"{site_id!r}:{sequence!r}:centre={centre_residue!r}"
            )
            continue

        parsed_site = None
        if site_column in site_metadata.columns:
            parsed_site = try_parse_site_token(site_metadata.at[site_id, site_column])
        explicit_residue = None
        if residue_column in site_metadata.columns:
            explicit_residue = _resolve_optional_residue(
                site_metadata.at[site_id, residue_column]
            )
        expected_residue = _resolve_expected_residue(parsed_site, explicit_residue)
        if expected_residue is None:
            if allow_unknown_site_residue:
                continue
            observed_site = (
                None
                if site_column not in site_metadata.columns
                else site_metadata.at[site_id, site_column]
            )
            observed_residue = (
                None
                if residue_column not in site_metadata.columns
                else site_metadata.at[site_id, residue_column]
            )
            unknown_expected_residue_rows.append(
                f"{site_id!r}:site={observed_site!r}:residue={observed_residue!r}"
            )
            continue
        if centre_residue != expected_residue:
            centre_mismatch_rows.append(
                f"{site_id!r}:expected={expected_residue!r}:"
                f"observed={centre_residue!r}:sequence={sequence!r}"
            )

    details: list[str] = []
    if missing_sequence_rows:
        details.append(
            f"missing or blank {field_name}.{site_sequence_column} values; "
            + _summarise_examples(missing_sequence_rows)
        )
    if even_length_rows:
        details.append(
            f"{field_name}.{site_sequence_column} must be odd length for centred "
            "context; " + _summarise_examples(even_length_rows)
        )
    if unsupported_character_rows:
        details.append(
            f"{field_name}.{site_sequence_column} contains unsupported non-letter "
            "characters under the configured centred-context policy; "
            + _summarise_examples(unsupported_character_rows)
        )
    if unknown_expected_residue_rows:
        details.append(
            "cannot resolve expected phosphosite residue from site/residue metadata "
            "under strict centred-context policy; "
            + _summarise_examples(unknown_expected_residue_rows)
        )
    if non_phospho_centre_rows:
        details.append(
            f"{field_name}.{site_sequence_column} centre residue must be one of "
            "S/T/Y; " + _summarise_examples(non_phospho_centre_rows)
        )
    if centre_mismatch_rows:
        details.append(
            f"{field_name}.{site_sequence_column} centre residue must match the "
            f"site token residue from {field_name}.{site_column}; "
            + _summarise_examples(centre_mismatch_rows)
        )

    if not details:
        return
    raise error_type(
        f"{workflow_name} requires centred sequence context in "
        f"{field_name}.{site_sequence_column}; " + "; ".join(details)
    )


def enforce_site_identity_rows(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    allow_opaque_site_values: bool = False,
) -> None:
    """Enforce row-level phosphosite identity parsing for workflow boundaries."""

    for row_position in range(len(site_metadata.index)):
        site_id = site_metadata.index[row_position]
        row = site_metadata.iloc[row_position]
        _ = build_phosphosite_identity(
            display_id=site_id,
            gene_symbol=row["gene_symbol"],
            site=row["site"],
            allow_opaque_site_values=allow_opaque_site_values,
            protein_id=(
                row["protein_id"] if "protein_id" in site_metadata.columns else None
            ),
            protein_accession=(
                row["protein_accession"]
                if "protein_accession" in site_metadata.columns
                else None
            ),
            field_name=f"{field_name}[{site_id!r}]",
            error_type=error_type,
        )


def assess_localisation_probability_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "localisation_probability",
) -> LocalisationProbabilityAssessment | None:
    """Parse optional localisation probability values with diagnostics."""

    if column_name not in site_metadata.columns:
        return None
    values = site_metadata.loc[:, column_name]
    normalized = pd.Series(pd.NA, index=values.index.copy(), dtype="Float64")
    missing_mask = pd.Series(False, index=values.index.copy(), dtype="boolean")
    invalid_mask = pd.Series(False, index=values.index.copy(), dtype="boolean")
    invalid_examples: list[str] = []

    for site_id, raw_value in values.items():
        parsed = _parse_localisation_probability(raw_value)
        if parsed is None:
            missing_mask.at[site_id] = True
            continue
        if isinstance(parsed, float):
            normalized.at[site_id] = parsed
            continue
        invalid_mask.at[site_id] = True
        if len(invalid_examples) < _EXAMPLE_LIMIT:
            invalid_examples.append(f"{site_id!r}:{raw_value!r}:{parsed}")

    if bool((missing_mask & invalid_mask).any()):
        raise error_type(
            f"{field_name}.{column_name} localisation parsing produced inconsistent "
            "missing/invalid masks"
        )
    return LocalisationProbabilityAssessment(
        normalized=normalized,
        missing_mask=missing_mask,
        invalid_mask=invalid_mask,
        invalid_examples=tuple(invalid_examples),
    )


def validate_localisation_confidence_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "localisation_confidence",
) -> None:
    """Validate optional localisation-confidence values when the column exists."""

    validate_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=column_name,
    )


def assess_localisation_confidence_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "localisation_confidence",
) -> LocalisationProbabilityAssessment | None:
    """Parse optional localisation-confidence values with diagnostics."""

    return assess_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=column_name,
    )


def _resolve_site_position_series(site_metadata: pd.DataFrame) -> pd.Series:
    for column_name in _SITE_POSITION_CANDIDATE_COLUMNS:
        if column_name in site_metadata.columns:
            series = site_metadata.loc[:, column_name].copy()
            series.name = column_name
            return series
    series = pd.Series(pd.NA, index=site_metadata.index.copy(), dtype="object")
    series.name = None
    return series


def _resolve_optional_residue(value: object) -> str | None:
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        return None
    token = value.strip().upper()
    if len(token) != 1:
        return None
    if not token.isalpha():
        return None
    return token


def _resolve_optional_position(value: object) -> int | None:
    if _is_missing(value):
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, numbers.Integral):
        integer_value = int(value)
        return integer_value if integer_value >= 1 else None
    if isinstance(value, numbers.Real):
        numeric_value = float(value)
        if not numeric_value.is_integer():
            return None
        integer_value = int(numeric_value)
        return integer_value if integer_value >= 1 else None
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "":
            return None
        try:
            parsed = int(stripped)
        except ValueError:
            return None
        return parsed if parsed >= 1 else None
    return None


def _resolve_optional_sequence(value: object) -> str | None:
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        return None
    stripped = value.strip().upper()
    if stripped == "":
        return None
    return stripped


def _resolve_central_residue(site_sequence: str) -> str | None:
    sequence_length = len(site_sequence)
    if sequence_length == 0 or sequence_length % 2 == 0:
        return None
    return site_sequence[sequence_length // 2]


def _sequence_supports_central_residue_check(site_sequence: str) -> bool:
    if not site_sequence.isalpha():
        return False
    return len(site_sequence) % 2 == 1


def _resolve_expected_residue(
    parsed_site: ParsedSiteToken | None,
    explicit_residue: str | None,
) -> str | None:
    if explicit_residue is not None:
        return explicit_residue
    if parsed_site is not None:
        return parsed_site.residue
    return None


def _parse_localisation_probability(value: object) -> float | str | None:
    if _is_missing(value):
        return None
    if isinstance(value, bool):
        return "bool_not_allowed"
    if isinstance(value, str):
        token = value.strip()
        if token == "":
            return None
        try:
            numeric_value = float(token)
        except ValueError:
            return "not_numeric"
    elif isinstance(value, (int, float)):
        numeric_value = float(value)
    else:
        return "unsupported_type"
    if not math.isfinite(numeric_value):
        return "not_finite"
    if numeric_value < 0.0 or numeric_value > 1.0:
        return "out_of_range"
    return float(numeric_value)


def _resolve_localisation_column_name(
    *,
    site_metadata: pd.DataFrame,
    requested_column_name: str,
) -> str | None:
    if requested_column_name in site_metadata.columns:
        return requested_column_name
    if (
        requested_column_name == "localisation_confidence"
        and "localisation_probability" in site_metadata.columns
    ):
        return "localisation_probability"
    return None


def _site_id_examples(index: pd.Index, *, limit: int = _EXAMPLE_LIMIT) -> str:
    labels = [str(value) for value in index.tolist()]
    if not labels:
        return "(none)"
    preview = ", ".join(repr(label) for label in labels[:limit])
    suffix = "" if len(labels) <= limit else f", +{len(labels) - limit} more"
    return f"[{preview}{suffix}]"


def _summarise_examples(values: list[str], *, limit: int = _EXAMPLE_LIMIT) -> str:
    if not values:
        return "(none)"
    preview = ", ".join(values[:limit])
    suffix = "" if len(values) <= limit else f", +{len(values) - limit} more"
    return f"[{preview}{suffix}]"


def _is_missing(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


__all__ = [
    "LocalisationProbabilityAssessment",
    "assess_localisation_confidence_column",
    "assess_localisation_probability_column",
    "enforce_centred_site_sequence_context",
    "enforce_site_identity_rows",
    "enforce_required_non_empty_string_column",
    "enforce_localisation_requirement",
    "validate_localisation_confidence_column",
    "validate_site_sequence_column",
    "validate_localisation_probability_column",
    "validate_site_identity_metadata",
]
