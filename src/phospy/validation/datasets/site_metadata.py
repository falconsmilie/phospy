"""Phosphosite metadata validation and workflow policy enforcement."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypeVar, cast

import numpy as np
import pandas as pd

from phospy.science.sites.identifiers import ParsedSiteToken, try_parse_site_token
from phospy.science.sites.identity import (
    build_phosphosite_identity,
    validate_identity_optional_columns,
)
from phospy.science.sites.site_keys import require_positive_integer_position

if TYPE_CHECKING:
    from phospy.contracts.configs.localisation import LocalisationRequirement

ErrorType = TypeVar("ErrorType", bound=Exception)
_PHOSPHORYLATABLE_RESIDUES = frozenset({"S", "T", "Y"})
_EXAMPLE_LIMIT = 5
_SITE_POSITION_CANDIDATE_COLUMNS = ("site_position", "position")
_SUPPORTED_GAP_SEQUENCE_CHARACTERS = frozenset({"_", "-"})
_CANONICAL_AMINO_ACID_RESIDUES = frozenset(
    {
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    }
)
_SUPPORTED_UNKNOWN_SEQUENCE_CHARACTERS = frozenset({"X"})
_SUPPORTED_BASE_SEQUENCE_SPECIAL_CHARACTERS = frozenset({"_", "-"})
_SUPPORTED_BASE_SEQUENCE_CHARACTERS = frozenset(
    _CANONICAL_AMINO_ACID_RESIDUES
    | _SUPPORTED_UNKNOWN_SEQUENCE_CHARACTERS
    | _SUPPORTED_BASE_SEQUENCE_SPECIAL_CHARACTERS
)
_SUPPORTED_MODIFIED_RESIDUE_SYMBOLS = frozenset({"*", "#", "@", "~"})
_UNKNOWN_SEQUENCE_SOURCE_TOKENS = frozenset({"", "unknown", "none", "na", "n/a"})


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


@dataclass(frozen=True, slots=True)
class SequenceContextContract:
    """Workflow-specific site-sequence context requirements."""

    requires_site_sequence: bool
    requires_centered_site: bool
    required_window_length: int | None
    center_index: int | None
    allowed_residues: frozenset[str]
    allow_terminal_padding: bool
    allow_lowercase: bool
    allow_modified_residue_symbols: bool
    required_center_residues: frozenset[str]
    requires_known_sequence_source: bool = False
    contract_id: str = "sequence_context"

    def __post_init__(self) -> None:
        allowed_residues = _normalise_residue_set(
            self.allowed_residues,
            field_name=f"{self.contract_id}.allowed_residues",
        )
        required_center_residues = _normalise_residue_set(
            self.required_center_residues,
            field_name=f"{self.contract_id}.required_center_residues",
        )
        required_window_length = cast(object, self.required_window_length)
        if required_window_length is not None:
            if isinstance(required_window_length, bool) or not isinstance(
                required_window_length, int
            ):
                raise ValueError(
                    f"{self.contract_id}.required_window_length must be an int or None"
                )
            if required_window_length <= 0:
                raise ValueError(
                    f"{self.contract_id}.required_window_length must be > 0"
                )
        center_index = cast(object, self.center_index)
        if center_index is not None:
            if isinstance(center_index, bool) or not isinstance(center_index, int):
                raise ValueError(f"{self.contract_id}.center_index must be an int")
            if center_index < 0:
                raise ValueError(f"{self.contract_id}.center_index must be >= 0")
            if (
                required_window_length is not None
                and center_index >= required_window_length
            ):
                raise ValueError(
                    f"{self.contract_id}.center_index must be within "
                    "required_window_length"
                )
        if (
            self.requires_centered_site
            and required_window_length is not None
            and center_index is None
        ):
            raise ValueError(
                f"{self.contract_id}.center_index is required when a centred fixed "
                "window is required"
            )
        object.__setattr__(self, "allowed_residues", allowed_residues)
        object.__setattr__(self, "required_center_residues", required_center_residues)


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

        explicit_position: int | None = None
        if site_positions.name is not None:
            raw_position = site_positions.at[site_id]
            try:
                explicit_position = require_positive_integer_position(
                    raw_position,
                    field_name=f"{field_name}.{site_positions.name}",
                    error_type=error_type,
                )
            except error_type:
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
            f"{field_name}.{site_positions.name} must contain positive integer "
            "values when the column is present; "
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
            _index_by_boolean_mask(site_metadata.index, assessment.invalid_mask)
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
            _index_by_boolean_mask(site_metadata.index, assessment.missing_mask)
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
        _index_by_boolean_mask(site_metadata.index, below_threshold)
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
    column = site_metadata[column_name]
    invalid_mask = pd.Series(False, index=pd.Index(column.index), dtype="boolean")
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
    invalid_sites = _site_id_examples(
        _index_by_boolean_mask(site_metadata.index, invalid_mask)
    )
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
    values = site_metadata[column_name]
    for site_id, raw_value in values.items():
        if not isinstance(raw_value, str):
            invalid_rows.append(f"{site_id!r}:{raw_value!r}")
            continue
        sequence = raw_value.strip().upper()
        if sequence == "":
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:blank_sequence")
            continue
        if len(sequence) < 3:
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:sequence_too_short")
            continue
        unsupported_characters = sorted(
            {
                character
                for character in sequence
                if character not in _SUPPORTED_BASE_SEQUENCE_CHARACTERS
            }
        )
        if unsupported_characters:
            invalid_rows.append(
                f"{site_id!r}:{raw_value!r}:unsupported_characters={''.join(unsupported_characters)!r}"
            )
            continue
        if not any(
            character in _CANONICAL_AMINO_ACID_RESIDUES for character in sequence
        ):
            invalid_rows.append(f"{site_id!r}:{raw_value!r}:no_residue_letters")
            continue
    if not invalid_rows:
        return
    raise error_type(
        f"{field_name}.{column_name} must be plausible amino-acid context strings "
        "(allowed residues: ACDEFGHIKLMNPQRSTVWY; allowed unknown: X; "
        "allowed gap placeholders: '_' and '-'); " + _summarise_examples(invalid_rows)
    )


def enforce_site_sequence_context_contract(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    workflow_name: str,
    contract: SequenceContextContract,
    error_type: type[ErrorType],
    scoring_mode: str | None = None,
    site_column: str = "site",
    site_sequence_column: str = "site_sequence",
    residue_column: str = "residue",
    sequence_source_column: str = "site_sequence_source",
    sequence_source_by_site: Mapping[Any, object] | None = None,
    allow_unknown_site_residue: bool = False,
) -> None:
    """Enforce workflow/method-specific sequence context requirements."""

    if (
        not contract.requires_site_sequence
        and site_sequence_column not in site_metadata
    ):
        return
    if site_sequence_column not in site_metadata.columns:
        site_examples = _site_id_examples(site_metadata.index)
        raise error_type(
            _sequence_contract_error_prefix(
                field_name=field_name,
                workflow_name=workflow_name,
                scoring_mode=scoring_mode,
                contract=contract,
                site_sequence_column=site_sequence_column,
            )
            + f"missing required column={field_name}.{site_sequence_column}; "
            f"affected_rows={int(site_metadata.shape[0])}; "
            f"example_site_ids={site_examples}"
        )

    invalid_rows: list[str] = []
    for site_id in site_metadata.index.tolist():
        raw_sequence = site_metadata.at[site_id, site_sequence_column]
        row_failure = _validate_sequence_context_row(
            site_metadata=site_metadata,
            site_id=site_id,
            raw_sequence=raw_sequence,
            field_name=field_name,
            contract=contract,
            site_column=site_column,
            residue_column=residue_column,
            sequence_source_column=sequence_source_column,
            sequence_source_by_site=sequence_source_by_site,
            allow_unknown_site_residue=allow_unknown_site_residue,
        )
        if row_failure is not None:
            invalid_rows.append(row_failure)

    if invalid_rows:
        raise error_type(
            _sequence_contract_error_prefix(
                field_name=field_name,
                workflow_name=workflow_name,
                scoring_mode=scoring_mode,
                contract=contract,
                site_sequence_column=site_sequence_column,
            )
            + _summarise_examples(invalid_rows)
        )

    if not contract.allow_modified_residue_symbols:
        validate_site_sequence_column(
            site_metadata=site_metadata,
            field_name=field_name,
            error_type=error_type,
            column_name=site_sequence_column,
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

    site_ids = site_metadata.index.tolist()
    gene_symbols = site_metadata["gene_symbol"].tolist()
    sites = site_metadata["site"].tolist()
    protein_ids = (
        site_metadata["protein_id"].tolist()
        if "protein_id" in site_metadata.columns
        else None
    )
    protein_accessions = (
        site_metadata["protein_accession"].tolist()
        if "protein_accession" in site_metadata.columns
        else None
    )
    for row_position, site_id in enumerate(site_ids):
        _ = build_phosphosite_identity(
            display_id=site_id,
            gene_symbol=gene_symbols[row_position],
            site=sites[row_position],
            allow_opaque_site_values=allow_opaque_site_values,
            protein_id=(None if protein_ids is None else protein_ids[row_position]),
            protein_accession=(
                None if protein_accessions is None else protein_accessions[row_position]
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
    values = site_metadata[column_name]
    values_index = pd.Index(values.index)
    missing_mask = values.isna()
    blank_string_mask = values.map(
        lambda value: isinstance(value, str) and value.strip() == ""
    )
    bool_mask = values.map(lambda value: isinstance(value, (bool, np.bool_)))
    parse_exempt_mask = missing_mask | blank_string_mask | bool_mask
    numeric_values = pd.to_numeric(values.mask(parse_exempt_mask), errors="coerce")
    finite_mask = pd.Series(
        np.isfinite(numeric_values.to_numpy(dtype=float, copy=False, na_value=np.nan)),
        index=values_index,
    )
    valid_numeric_mask = (
        ~parse_exempt_mask
        & numeric_values.notna()
        & finite_mask
        & numeric_values.ge(0.0)
        & numeric_values.le(1.0)
    )
    missing_mask = pd.Series(
        (missing_mask | blank_string_mask).to_numpy(dtype=bool, copy=False),
        index=values_index,
        dtype="boolean",
    )
    invalid_mask = pd.Series(
        ((~missing_mask.astype(bool)) & (~valid_numeric_mask)).to_numpy(
            dtype=bool,
            copy=False,
        ),
        index=values_index,
        dtype="boolean",
    )
    normalized = pd.Series(pd.NA, index=values_index, dtype="Float64")
    if bool(valid_numeric_mask.any()):
        normalized.loc[valid_numeric_mask] = numeric_values.loc[
            valid_numeric_mask
        ].astype(float)
    invalid_examples: list[str] = []
    invalid_positions = np.flatnonzero(invalid_mask.to_numpy(dtype=bool, copy=False))
    for position in invalid_positions[:_EXAMPLE_LIMIT]:
        site_id = values.index[int(position)]
        raw_value = values.at[site_id]
        parsed = _parse_localisation_probability(raw_value)
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
            series = site_metadata[column_name].copy(deep=True)
            series.name = column_name
            return series
    series = pd.Series(pd.NA, index=pd.Index(site_metadata.index), dtype="object")
    series.name = None
    return series


def _validate_sequence_context_row(
    *,
    site_metadata: pd.DataFrame,
    site_id: object,
    raw_sequence: object,
    field_name: str,
    contract: SequenceContextContract,
    site_column: str,
    residue_column: str,
    sequence_source_column: str,
    sequence_source_by_site: Mapping[Any, object] | None,
    allow_unknown_site_residue: bool,
) -> str | None:
    expected_length = contract.required_window_length
    expected_center_index = contract.center_index
    expected_center_residues = _format_residue_set(contract.required_center_residues)
    if _is_missing(raw_sequence) or not isinstance(raw_sequence, str):
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=None,
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason=f"missing or blank {field_name}.site_sequence values",
        )
    if raw_sequence != raw_sequence.strip():
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=len(raw_sequence),
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason="leading/trailing whitespace is not allowed",
        )
    if raw_sequence == "":
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=0,
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason=f"missing or blank {field_name}.site_sequence values",
        )
    if not contract.allow_lowercase and any(
        character.islower() for character in raw_sequence
    ):
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=len(raw_sequence),
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason="lowercase characters are not allowed by this sequence contract",
        )
    sequence = raw_sequence.upper()
    sequence_length = len(sequence)
    if expected_length is not None and sequence_length != expected_length:
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason=(
                "sequence length is wrong for the selected workflow/scoring mode; "
                f"observed_length={sequence_length}"
            ),
        )
    if expected_length is None and contract.requires_centered_site:
        if sequence_length % 2 == 0:
            return _sequence_context_row_message(
                site_id=site_id,
                sequence=raw_sequence,
                observed_length=sequence_length,
                expected_length=expected_length,
                expected_center_index=expected_center_index,
                expected_center_residues=expected_center_residues,
                reason="site_sequence must be odd length for centred context",
            )

    unsupported_characters = _unsupported_sequence_characters(
        sequence=sequence,
        contract=contract,
    )
    if unsupported_characters:
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason=(
                "sequence alphabet is invalid; unsupported_characters="
                f"{''.join(unsupported_characters)!r}"
            ),
        )
    padding_failure = _validate_padding_policy(sequence=sequence, contract=contract)
    if padding_failure is not None:
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason=padding_failure,
        )
    if not any(character in _CANONICAL_AMINO_ACID_RESIDUES for character in sequence):
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason="sequence contains no residue letters",
        )
    if contract.requires_known_sequence_source:
        source = _resolve_sequence_source(
            site_metadata=site_metadata,
            site_id=site_id,
            sequence_source_column=sequence_source_column,
            sequence_source_by_site=sequence_source_by_site,
        )
        if source is None or source.lower() in _UNKNOWN_SEQUENCE_SOURCE_TOKENS:
            return _sequence_context_row_message(
                site_id=site_id,
                sequence=raw_sequence,
                observed_length=sequence_length,
                expected_length=expected_length,
                expected_center_index=expected_center_index,
                expected_center_residues=expected_center_residues,
                reason="sequence source is unknown",
            )
    if not contract.requires_centered_site:
        return None
    center_index = _resolve_contract_center_index(
        sequence_length=sequence_length,
        contract=contract,
    )
    if center_index is None:
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=expected_center_index,
            expected_center_residues=expected_center_residues,
            reason="required centered sequence is not centered",
        )
    if center_index >= sequence_length:
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=center_index,
            expected_center_residues=expected_center_residues,
            reason="center index does not exist in sequence",
        )
    center_residue = sequence[center_index]
    if center_residue not in contract.required_center_residues:
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=center_index,
            expected_center_residues=expected_center_residues,
            reason=(
                f"center residue {center_residue!r} is not allowed; "
                f"expected_center_residues={expected_center_residues}"
            ),
        )
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
            return None
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=center_index,
            expected_center_residues=expected_center_residues,
            reason="cannot resolve expected phosphosite residue from site/residue metadata",
        )
    if center_residue != expected_residue:
        return _sequence_context_row_message(
            site_id=site_id,
            sequence=raw_sequence,
            observed_length=sequence_length,
            expected_length=expected_length,
            expected_center_index=center_index,
            expected_center_residues=expected_center_residues,
            reason=(
                "centre residue must match the site token residue from "
                f"{field_name}.{site_column}; expected={expected_residue!r}; "
                f"observed={center_residue!r}"
            ),
        )
    return None


def _sequence_context_row_message(
    *,
    site_id: object,
    sequence: object,
    observed_length: int | None,
    expected_length: int | None,
    expected_center_index: int | None,
    expected_center_residues: str,
    reason: str,
) -> str:
    observed = "None" if observed_length is None else str(observed_length)
    expected = "None" if expected_length is None else str(expected_length)
    center_index = (
        "None" if expected_center_index is None else str(expected_center_index)
    )
    return (
        f"site_key={str(site_id)!r}: sequence={sequence!r}: "
        f"observed_length={observed}: expected_length={expected}: "
        f"expected_center_index={center_index}: "
        f"expected_center_residues={expected_center_residues}: {reason}"
    )


def _sequence_contract_error_prefix(
    *,
    field_name: str,
    workflow_name: str,
    scoring_mode: str | None,
    contract: SequenceContextContract,
    site_sequence_column: str,
) -> str:
    centered_label = (
        "centred sequence context"
        if contract.requires_centered_site
        else "sequence context"
    )
    mode_label = "None" if scoring_mode is None else str(scoring_mode)
    return (
        f"{workflow_name} requires {centered_label} in "
        f"{field_name}.{site_sequence_column}; workflow-specific sequence context "
        f"contract failed (contract={contract.contract_id}, "
        f"scoring_mode={mode_label}); "
    )


def _resolve_contract_center_index(
    *,
    sequence_length: int,
    contract: SequenceContextContract,
) -> int | None:
    if contract.center_index is not None:
        return int(contract.center_index)
    if sequence_length <= 0 or sequence_length % 2 == 0:
        return None
    return sequence_length // 2


def _unsupported_sequence_characters(
    *,
    sequence: str,
    contract: SequenceContextContract,
) -> list[str]:
    allowed = set(contract.allowed_residues)
    unsupported = {
        character
        for character in sequence
        if character not in allowed
        and character not in _SUPPORTED_BASE_SEQUENCE_SPECIAL_CHARACTERS
        and not (
            contract.allow_modified_residue_symbols
            and character in _SUPPORTED_MODIFIED_RESIDUE_SYMBOLS
        )
    }
    return sorted(unsupported)


def _validate_padding_policy(
    *,
    sequence: str,
    contract: SequenceContextContract,
) -> str | None:
    padding_positions = [
        position
        for position, character in enumerate(sequence)
        if character in _SUPPORTED_BASE_SEQUENCE_SPECIAL_CHARACTERS
        and character not in contract.allowed_residues
    ]
    if not padding_positions:
        return None
    if not contract.allow_terminal_padding:
        return "terminal padding is not allowed by this sequence contract"
    if not _padding_positions_are_terminal(sequence, padding_positions):
        return "padding is present but is not terminal padding"
    return None


def _padding_positions_are_terminal(
    sequence: str,
    padding_positions: list[int],
) -> bool:
    padding_lookup = set(padding_positions)
    non_padding_positions = [
        position for position in range(len(sequence)) if position not in padding_lookup
    ]
    if not non_padding_positions:
        return False
    first_residue = min(non_padding_positions)
    last_residue = max(non_padding_positions)
    return all(
        position < first_residue or position > last_residue
        for position in padding_positions
    )


def _resolve_sequence_source(
    *,
    site_metadata: pd.DataFrame,
    site_id: object,
    sequence_source_column: str,
    sequence_source_by_site: Mapping[Any, object] | None,
) -> str | None:
    if sequence_source_by_site is not None:
        if site_id in sequence_source_by_site:
            return _optional_source_text(sequence_source_by_site[site_id])
        site_key = str(site_id)
        if site_key in sequence_source_by_site:
            return _optional_source_text(sequence_source_by_site[site_key])
    if sequence_source_column not in site_metadata.columns:
        return None
    return _optional_source_text(site_metadata.at[site_id, sequence_source_column])


def _optional_source_text(value: object) -> str | None:
    if _is_missing(value):
        return None
    if not isinstance(value, str):
        return None
    text = value.strip()
    if text == "":
        return None
    return text


def _normalise_residue_set(values: object, *, field_name: str) -> frozenset[str]:
    if isinstance(values, str):
        iterable: tuple[object, ...] = tuple(values)
    else:
        try:
            iterable = tuple(cast(Iterable[object], values))
        except TypeError as exc:
            raise ValueError(
                f"{field_name} must be an iterable of residue tokens"
            ) from exc
    normalised: set[str] = set()
    for value in iterable:
        if not isinstance(value, str):
            raise ValueError(f"{field_name} must contain string residue tokens")
        token = value.strip().upper()
        if len(token) != 1:
            raise ValueError(f"{field_name} must contain one-character tokens")
        normalised.add(token)
    return frozenset(normalised)


def _format_residue_set(values: frozenset[str]) -> str:
    if not values:
        return "(none)"
    return "/".join(sorted(values))


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


def _index_by_boolean_mask(index: pd.Index, mask: pd.Series) -> pd.Index:
    labels = index.tolist()
    mask_values = cast(list[object], mask.tolist())
    selected = [
        label
        for label, include in zip(labels, mask_values, strict=True)
        if bool(include)
    ]
    return pd.Index(selected)


def _summarise_examples(values: list[str], *, limit: int = _EXAMPLE_LIMIT) -> str:
    if not values:
        return "(none)"
    preview = ", ".join(values[:limit])
    suffix = "" if len(values) <= limit else f", +{len(values) - limit} more"
    return f"[{preview}{suffix}]"


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value = cast(object, value)
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value = cast(object, value)
        return str(temporal_value) == "NaT"
    return False


__all__ = [
    "LocalisationProbabilityAssessment",
    "SequenceContextContract",
    "assess_localisation_confidence_column",
    "assess_localisation_probability_column",
    "enforce_centred_site_sequence_context",
    "enforce_site_sequence_context_contract",
    "enforce_site_identity_rows",
    "enforce_required_non_empty_string_column",
    "enforce_localisation_requirement",
    "validate_localisation_confidence_column",
    "validate_site_sequence_column",
    "validate_localisation_probability_column",
    "validate_site_identity_metadata",
]
