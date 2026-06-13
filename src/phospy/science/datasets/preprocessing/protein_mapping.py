"""Protein-aware phosphosite-to-total-protein mapping resolution."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from numbers import Real
from typing import Any

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.policies import PolicyEnum

_INDEX_KEY = "__index__"
_DEFAULT_PROTEIN_IDENTIFIER_COLUMNS = (
    "protein_accession",
    "protein_id",
    "protein_group_id",
)
_DEFAULT_DISPLAY_LABEL_COLUMNS = ("display_id", "protein_label")
_DEFAULT_GENE_SYMBOL_COLUMNS = ("gene_symbol",)
_DEFAULT_MULTI_VALUE_DELIMITERS = (";", "|", ",")


class ProteinMappingStatus(PolicyEnum):
    """Per-site status for protein-aware mapping resolution."""

    MATCHED = "matched"
    MISSING_SITE_PROTEIN_IDENTIFIER = "missing_site_protein_identifier"
    MISSING_TOTAL_PROTEIN_ROW = "missing_total_protein_row"
    AMBIGUOUS_SITE_PROTEIN_MAPPING = "ambiguous_site_protein_mapping"
    AMBIGUOUS_TOTAL_PROTEIN_MAPPING = "ambiguous_total_protein_mapping"


@dataclass(frozen=True, slots=True)
class ProteinMappingConfig:
    """Configuration for resolving phosphosite rows to total-protein rows.

    Gene symbols are intentionally excluded unless `allow_gene_symbol_matching`
    is set. Display-label fallback is also opt-in because explicit protein
    identifiers should drive protein-aware preparation.
    """

    site_metadata_key: str = _INDEX_KEY
    total_protein_key: str = _INDEX_KEY
    protein_identifier_columns: Sequence[str] = _DEFAULT_PROTEIN_IDENTIFIER_COLUMNS
    allow_display_label_fallback: bool = False
    display_label_columns: Sequence[str] = _DEFAULT_DISPLAY_LABEL_COLUMNS
    allow_gene_symbol_matching: bool = False
    gene_symbol_columns: Sequence[str] = _DEFAULT_GENE_SYMBOL_COLUMNS
    multi_value_delimiters: Sequence[str] = _DEFAULT_MULTI_VALUE_DELIMITERS

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "site_metadata_key",
            _require_non_empty_string(
                self.site_metadata_key,
                field_name="protein mapping site_metadata_key",
            ),
        )
        object.__setattr__(
            self,
            "total_protein_key",
            _require_non_empty_string(
                self.total_protein_key,
                field_name="protein mapping total_protein_key",
            ),
        )
        object.__setattr__(
            self,
            "protein_identifier_columns",
            _normalize_string_tuple(
                self.protein_identifier_columns,
                field_name="protein mapping protein_identifier_columns",
            ),
        )
        object.__setattr__(
            self,
            "display_label_columns",
            _normalize_string_tuple(
                self.display_label_columns,
                field_name="protein mapping display_label_columns",
            ),
        )
        object.__setattr__(
            self,
            "gene_symbol_columns",
            _normalize_string_tuple(
                self.gene_symbol_columns,
                field_name="protein mapping gene_symbol_columns",
            ),
        )
        object.__setattr__(
            self,
            "multi_value_delimiters",
            _normalize_string_tuple(
                self.multi_value_delimiters,
                field_name="protein mapping multi_value_delimiters",
            ),
        )


@dataclass(frozen=True, slots=True)
class ProteinMappingRecord:
    """Resolved mapping record for one phosphosite row key."""

    site_key: str
    protein_identifier: str | None
    total_protein_row_key: str | None
    status: ProteinMappingStatus
    protein_identifier_source: str | None = None
    candidate_protein_identifiers: tuple[str, ...] = ()
    candidate_total_protein_row_keys: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ProteinMappingResult:
    """Complete resolver output for protein-aware preparation."""

    records: tuple[ProteinMappingRecord, ...]
    site_to_protein_identifier: dict[str, str]
    site_to_total_protein_row_key: dict[str, str]
    status_by_site: dict[str, ProteinMappingStatus]


@dataclass(frozen=True, slots=True)
class _CandidateProteinIdentifier:
    identifiers: tuple[str, ...]
    source: str | None


class ProteinMappingResolver:
    """Resolve phosphosite rows to explicit protein and total-protein row keys."""

    def run(
        self,
        site_metadata: pd.DataFrame,
        phospho_matrix_index: pd.Index,
        total_protein_matrix_index: pd.Index,
        config: ProteinMappingConfig | None = None,
    ) -> ProteinMappingResult:
        resolved_config = config or ProteinMappingConfig()
        metadata_by_site = _build_site_metadata_lookup(
            site_metadata=site_metadata,
            site_metadata_key=resolved_config.site_metadata_key,
        )
        total_rows_by_key = _build_total_row_lookup(
            total_protein_matrix_index=total_protein_matrix_index,
            total_protein_key=resolved_config.total_protein_key,
        )

        records: list[ProteinMappingRecord] = []
        site_to_protein_identifier: dict[str, str] = {}
        site_to_total_protein_row_key: dict[str, str] = {}
        status_by_site: dict[str, ProteinMappingStatus] = {}

        for site_key in _normalize_index_values(phospho_matrix_index):
            metadata_rows = metadata_by_site.get(site_key, ())
            candidate = _resolve_candidate_protein_identifier(
                metadata_rows=metadata_rows,
                config=resolved_config,
            )
            record = _resolve_site_mapping_record(
                site_key=site_key,
                candidate=candidate,
                total_rows_by_key=total_rows_by_key,
            )
            records.append(record)
            status_by_site[site_key] = record.status
            if (
                record.protein_identifier is not None
                and record.status
                is not ProteinMappingStatus.AMBIGUOUS_SITE_PROTEIN_MAPPING
            ):
                site_to_protein_identifier[site_key] = record.protein_identifier
            if record.status is ProteinMappingStatus.MATCHED:
                if record.total_protein_row_key is None:  # pragma: no cover
                    raise AssertionError("matched protein mapping lacks total row key")
                site_to_total_protein_row_key[site_key] = record.total_protein_row_key

        return ProteinMappingResult(
            records=tuple(records),
            site_to_protein_identifier=site_to_protein_identifier,
            site_to_total_protein_row_key=site_to_total_protein_row_key,
            status_by_site=status_by_site,
        )


def _resolve_site_mapping_record(
    *,
    site_key: str,
    candidate: _CandidateProteinIdentifier,
    total_rows_by_key: dict[str, tuple[str, ...]],
) -> ProteinMappingRecord:
    identifiers = candidate.identifiers
    if not identifiers:
        return ProteinMappingRecord(
            site_key=site_key,
            protein_identifier=None,
            total_protein_row_key=None,
            status=ProteinMappingStatus.MISSING_SITE_PROTEIN_IDENTIFIER,
            protein_identifier_source=candidate.source,
            candidate_protein_identifiers=(),
            candidate_total_protein_row_keys=(),
        )
    if len(identifiers) > 1:
        return ProteinMappingRecord(
            site_key=site_key,
            protein_identifier=None,
            total_protein_row_key=None,
            status=ProteinMappingStatus.AMBIGUOUS_SITE_PROTEIN_MAPPING,
            protein_identifier_source=candidate.source,
            candidate_protein_identifiers=identifiers,
            candidate_total_protein_row_keys=(),
        )

    protein_identifier = next(iter(identifiers))
    candidate_total_rows = total_rows_by_key.get(protein_identifier, ())
    if not candidate_total_rows:
        return ProteinMappingRecord(
            site_key=site_key,
            protein_identifier=protein_identifier,
            total_protein_row_key=None,
            status=ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
            protein_identifier_source=candidate.source,
            candidate_protein_identifiers=identifiers,
            candidate_total_protein_row_keys=(),
        )
    if len(candidate_total_rows) > 1:
        return ProteinMappingRecord(
            site_key=site_key,
            protein_identifier=protein_identifier,
            total_protein_row_key=None,
            status=ProteinMappingStatus.AMBIGUOUS_TOTAL_PROTEIN_MAPPING,
            protein_identifier_source=candidate.source,
            candidate_protein_identifiers=identifiers,
            candidate_total_protein_row_keys=candidate_total_rows,
        )
    total_protein_row_key = next(iter(candidate_total_rows))
    return ProteinMappingRecord(
        site_key=site_key,
        protein_identifier=protein_identifier,
        total_protein_row_key=total_protein_row_key,
        status=ProteinMappingStatus.MATCHED,
        protein_identifier_source=candidate.source,
        candidate_protein_identifiers=identifiers,
        candidate_total_protein_row_keys=candidate_total_rows,
    )


def _resolve_candidate_protein_identifier(
    *,
    metadata_rows: Sequence[dict[str, Any]],
    config: ProteinMappingConfig,
) -> _CandidateProteinIdentifier:
    explicit = _collect_identifiers_from_columns(
        metadata_rows=metadata_rows,
        columns=config.protein_identifier_columns,
        config=config,
        source_kind="protein_identifier",
    )
    if explicit.identifiers:
        return explicit

    if config.allow_display_label_fallback:
        display_label = _collect_identifiers_from_columns(
            metadata_rows=metadata_rows,
            columns=config.display_label_columns,
            config=config,
            source_kind="display_label",
        )
        if display_label.identifiers:
            return display_label

    if config.allow_gene_symbol_matching:
        return _collect_identifiers_from_columns(
            metadata_rows=metadata_rows,
            columns=config.gene_symbol_columns,
            config=config,
            source_kind="gene_symbol",
        )
    return _CandidateProteinIdentifier(identifiers=(), source=None)


def _collect_identifiers_from_columns(
    *,
    metadata_rows: Sequence[dict[str, Any]],
    columns: Sequence[str],
    config: ProteinMappingConfig,
    source_kind: str,
) -> _CandidateProteinIdentifier:
    values: list[str] = []
    first_source: str | None = None
    for row in metadata_rows:
        for column in columns:
            if _is_gene_symbol_column(column) and not config.allow_gene_symbol_matching:
                continue
            if column not in row:
                continue
            column_values = _identifier_tokens(
                row[column],
                delimiters=config.multi_value_delimiters,
            )
            if not column_values:
                continue
            if first_source is None:
                first_source = f"{source_kind}:{column}"
            values.extend(column_values)
    return _CandidateProteinIdentifier(
        identifiers=_dedupe_preserving_order(values),
        source=first_source,
    )


def _build_site_metadata_lookup(
    *,
    site_metadata: pd.DataFrame,
    site_metadata_key: str,
) -> dict[str, tuple[dict[str, Any], ...]]:
    row_groups: dict[str, list[dict[str, Any]]] = {}
    if site_metadata_key == _INDEX_KEY:
        site_keys = site_metadata.index.tolist()
    elif _index_name_matches(site_metadata.index, site_metadata_key):
        site_keys = site_metadata.index.tolist()
    else:
        if site_metadata_key not in site_metadata.columns:
            raise PhosPyInputError(
                "protein mapping site_metadata_key "
                f"{site_metadata_key!r} is not present in site_metadata"
            )
        site_keys = site_metadata.loc[:, site_metadata_key].tolist()

    records = [
        {str(key): value for key, value in record.items()}
        for record in site_metadata.to_dict(orient="records")
    ]

    for site_key, record in zip(
        site_keys,
        records,
        strict=True,
    ):
        normalized_key = _normalize_identifier(site_key)
        if normalized_key is None:
            continue
        row_groups.setdefault(normalized_key, []).append(record)
    return {key: tuple(rows) for key, rows in row_groups.items()}


def _build_total_row_lookup(
    *,
    total_protein_matrix_index: pd.Index,
    total_protein_key: str,
) -> dict[str, tuple[str, ...]]:
    if total_protein_key != _INDEX_KEY and not _index_name_matches(
        total_protein_matrix_index,
        total_protein_key,
    ):
        raise PhosPyInputError(
            "protein mapping total_protein_key "
            f"{total_protein_key!r} cannot be resolved from the total protein "
            "matrix index"
        )
    row_groups: dict[str, list[str]] = {}
    for value in total_protein_matrix_index.tolist():
        normalized_key = _normalize_identifier(value)
        if normalized_key is None:
            continue
        row_groups.setdefault(normalized_key, []).append(str(value).strip())
    return {key: tuple(rows) for key, rows in row_groups.items()}


def _normalize_index_values(index: pd.Index) -> tuple[str, ...]:
    values: list[str] = []
    for value in index.tolist():
        normalized = _normalize_identifier(value)
        if normalized is None:
            values.append("")
        else:
            values.append(normalized)
    return tuple(values)


def _identifier_tokens(
    value: object,
    *,
    delimiters: Sequence[str],
) -> tuple[str, ...]:
    if _is_missing_scalar(value):
        return ()
    if isinstance(value, str):
        return tuple(
            token
            for token in _split_identifier_string(value, delimiters=delimiters)
            if token
        )
    if isinstance(value, (set, frozenset)):
        iterable_value: Iterable[object] = sorted(value, key=lambda item: str(item))
    elif isinstance(value, Iterable):
        iterable_value = value
    else:
        normalized = _normalize_identifier(value)
        return () if normalized is None else (normalized,)

    tokens: list[str] = []
    for item in iterable_value:
        tokens.extend(_identifier_tokens(item, delimiters=delimiters))
    return tuple(tokens)


def _split_identifier_string(
    value: str,
    *,
    delimiters: Sequence[str],
) -> tuple[str, ...]:
    tokens = (value,)
    for delimiter in delimiters:
        next_tokens: list[str] = []
        for token in tokens:
            next_tokens.extend(token.split(delimiter))
        tokens = tuple(next_tokens)
    return tuple(token.strip() for token in tokens if token.strip())


def _normalize_identifier(value: object) -> str | None:
    if _is_missing_scalar(value):
        return None
    normalized = str(value).strip()
    if normalized == "":
        return None
    return normalized


def _is_missing_scalar(value: object) -> bool:
    if isinstance(value, (list, tuple, set, frozenset, dict)):
        return False
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, Real) and not isinstance(value, bool):
        return bool(pd.isna(float(value)))
    return False


def _dedupe_preserving_order(values: Sequence[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    deduped: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        deduped.append(value)
    return tuple(deduped)


def _normalize_string_tuple(
    values: Sequence[str], *, field_name: str
) -> tuple[str, ...]:
    if isinstance(values, str):
        values = (values,)
    normalized: list[str] = []
    for value in values:
        resolved = _require_non_empty_string(value, field_name=field_name)
        if resolved not in normalized:
            normalized.append(resolved)
    return tuple(normalized)


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def _index_name_matches(index: pd.Index, key: str) -> bool:
    if index.name is None:
        return False
    return str(index.name).strip() == key


def _is_gene_symbol_column(column: str) -> bool:
    return column.strip().lower() in {"gene", "gene_name", "gene_symbol", "symbol"}


__all__ = [
    "ProteinMappingConfig",
    "ProteinMappingRecord",
    "ProteinMappingResolver",
    "ProteinMappingResult",
    "ProteinMappingStatus",
]
