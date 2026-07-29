"""Science-owned phosphosite sequence-context contracts and enforcement."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, TypeVar, cast

import numpy as np
import pandas as pd

from phospy.science.sites.identifiers import ParsedSiteToken, try_parse_site_token
from phospy.science.sites.metadata_validation import validate_site_sequence_column

ErrorType = TypeVar("ErrorType", bound=Exception)
_EXAMPLE_LIMIT = 5
_PHOSPHORYLATABLE_RESIDUES = frozenset({"S", "T", "Y"})
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
_SUPPORTED_BASE_SEQUENCE_SPECIAL_CHARACTERS = frozenset({"_", "-"})
_SUPPORTED_MODIFIED_RESIDUE_SYMBOLS = frozenset({"*", "#", "@", "~"})
_UNKNOWN_SEQUENCE_SOURCE_TOKENS = frozenset({"", "unknown", "none", "na", "n/a"})


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
        and site_sequence_column not in site_metadata.columns
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


def _resolve_expected_residue(
    parsed_site: ParsedSiteToken | None,
    explicit_residue: str | None,
) -> str | None:
    if explicit_residue is not None:
        return explicit_residue
    if parsed_site is not None:
        return parsed_site.residue
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


WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT = SequenceContextContract(
    requires_site_sequence=True,
    requires_centered_site=True,
    required_window_length=None,
    center_index=None,
    allowed_residues=frozenset("ACDEFGHIKLMNPQRSTVWYX_-"),
    allow_terminal_padding=True,
    allow_lowercase=True,
    allow_modified_residue_symbols=False,
    required_center_residues=frozenset({"S", "T", "Y"}),
    contract_id="workflow_centered_phosphosite_sequence_context",
)


__all__ = [
    "SequenceContextContract",
    "WORKFLOW_CENTERED_SEQUENCE_CONTEXT_CONTRACT",
    "enforce_centred_site_sequence_context",
    "enforce_site_sequence_context_contract",
]
