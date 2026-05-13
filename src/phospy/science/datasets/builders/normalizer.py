"""Internal convention normalisation for dataset builder inputs."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.science.sites.identifiers import (
    SiteIdentifierNormalisationRecord,
    SiteIdentifierNormalisationReport,
    build_site_identifier_normalisation_report,
    canonicalize_site_components_series,
    canonicalize_site_index,
    parse_canonical_site_identifier,
)

_GENE_SYMBOL_ALIASES = ("gene_symbol", "gene_name")
_PROTEIN_ID_ALIASES = ("protein_id",)
_SITE_ALIASES = ("site",)
_SITE_SEQUENCE_ALIASES = ("site_sequence", "centralized_sequence")
_LOCALISATION_CONFIDENCE_ALIASES = (
    "localisation_confidence",
    "localization_confidence",
    "localisation_probability",
    "localization_probability",
)


@dataclass(frozen=True, slots=True)
class NormalizedDatasetInputs:
    """Normalised tables ready for sequence derivation and execution."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    site_identifier_normalisation: SiteIdentifierNormalisationReport | None = None


class DatasetConventionNormalizer:
    """Apply narrow, documented shaping rules for supported inputs.

    The supported index derivation convention can populate `gene_symbol` and `site`
    only. It never infers `protein_id`.
    """

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
    ) -> NormalizedDatasetInputs:
        site_identifier_records: list[SiteIdentifierNormalisationRecord] = []
        normalized_phospho = phospho.copy(deep=True)
        normalized_phospho.index = _normalize_supported_site_index_if_present(
            normalized_phospho.index,
            field_name="dataset build request phospho.index",
            site_identifier_records=site_identifier_records,
        )
        normalized_phospho.columns = _normalize_index_labels(
            normalized_phospho.columns,
            field_name="dataset build request phospho.columns",
            policy=_SAMPLE_LABEL_INDEX_POLICY,
        )

        normalized_site_metadata = self._normalize_site_metadata(
            site_metadata.copy(deep=True),
            phospho_index=normalized_phospho.index,
            site_identifier_records=site_identifier_records,
        )

        normalized_sample_metadata = self._normalize_sample_metadata(
            None if sample_metadata is None else sample_metadata.copy(deep=True),
            phospho_columns=normalized_phospho.columns,
        )
        normalized_total = self._normalize_total(
            None if total is None else total.copy(deep=True),
            phospho_columns=normalized_phospho.columns,
        )
        return NormalizedDatasetInputs(
            phospho=normalized_phospho,
            site_metadata=normalized_site_metadata,
            sample_metadata=normalized_sample_metadata,
            total=normalized_total,
            site_identifier_normalisation=build_site_identifier_normalisation_report(
                site_identifier_records
            ),
        )

    def _normalize_site_metadata(
        self,
        site_metadata: pd.DataFrame,
        *,
        phospho_index: pd.Index,
        site_identifier_records: list[SiteIdentifierNormalisationRecord],
    ) -> pd.DataFrame:
        normalized = site_metadata
        normalized = self._normalize_site_metadata_index(
            normalized, site_identifier_records=site_identifier_records
        )
        normalized = self._normalize_site_metadata_columns(normalized)
        normalized = self._derive_site_fields_from_index(normalized)
        normalized = self._normalize_site_metadata_site_identity_fields(normalized)
        if (
            not normalized.index.equals(phospho_index)
            and normalized.index.isin(phospho_index).all()
            and phospho_index.isin(normalized.index).all()
        ):
            normalized = normalized.reindex(phospho_index)
        return normalized

    @staticmethod
    def _normalize_site_metadata_index(
        site_metadata: pd.DataFrame,
        *,
        site_identifier_records: list[SiteIdentifierNormalisationRecord],
    ) -> pd.DataFrame:
        normalized = site_metadata
        index_name = (
            str(normalized.index.name).strip().lower() if normalized.index.name else ""
        )
        has_site_id_column = "site_id" in normalized.columns
        if has_site_id_column and isinstance(normalized.index, pd.RangeIndex):
            normalized_site_id = _canonicalize_site_index_with_label_validation(
                pd.Index(normalized.loc[:, "site_id"], name="site_id"),
                field_name="dataset build request site_metadata.site_id",
                site_identifier_records=site_identifier_records,
                index_name="site_id",
            )
            normalized["site_id"] = normalized_site_id.tolist()
            normalized = normalized.set_index("site_id", drop=True)
        elif has_site_id_column:
            column_as_site_id = _canonicalize_site_index_with_label_validation(
                pd.Index(normalized.loc[:, "site_id"], name="site_id"),
                field_name="dataset build request site_metadata.site_id",
                site_identifier_records=site_identifier_records,
            )
            if not _has_site_like_tokens(normalized.index):
                normalized = normalized.copy()
                normalized["site_id"] = column_as_site_id.tolist()
                normalized = normalized.set_index("site_id", drop=True)
                normalized.index = pd.Index(normalized.index.tolist(), name="site_id")
                return normalized
            index_as_site_id = _canonicalize_site_index_with_label_validation(
                normalized.index,
                field_name="dataset build request site_metadata.index",
                site_identifier_records=site_identifier_records,
            )
            if not index_as_site_id.equals(column_as_site_id):
                raise UnsupportedInputFormatError(
                    "dataset build request site_metadata has conflicting site "
                    "identifiers between index and 'site_id' column. use one source "
                    "of site IDs or make both exactly match"
                )
            if index_name != "site_id":
                normalized.index = pd.Index(
                    index_as_site_id.tolist(),
                    name="site_id",
                )
        normalized.index = _normalize_supported_site_index_if_present(
            normalized.index,
            field_name="dataset build request site_metadata.index",
            site_identifier_records=site_identifier_records,
        )
        return normalized

    @staticmethod
    def _normalize_site_metadata_columns(site_metadata: pd.DataFrame) -> pd.DataFrame:
        _reject_unsupported_historical_aliases(site_metadata)
        rename_map: dict[str, str] = {}
        gene_column = _resolve_alias(
            site_metadata.columns,
            target_column="gene_symbol",
            aliases=_GENE_SYMBOL_ALIASES,
        )
        if gene_column is not None and gene_column != "gene_symbol":
            rename_map[gene_column] = "gene_symbol"
        protein_column = _resolve_alias(
            site_metadata.columns,
            target_column="protein_id",
            aliases=_PROTEIN_ID_ALIASES,
        )
        if protein_column is not None and protein_column != "protein_id":
            rename_map[protein_column] = "protein_id"
        site_column = _resolve_alias(
            site_metadata.columns,
            target_column="site",
            aliases=_SITE_ALIASES,
        )
        if site_column is not None and site_column != "site":
            rename_map[site_column] = "site"
        site_sequence_column = _resolve_alias(
            site_metadata.columns,
            target_column="site_sequence",
            aliases=_SITE_SEQUENCE_ALIASES,
        )
        if site_sequence_column is not None and site_sequence_column != "site_sequence":
            rename_map[site_sequence_column] = "site_sequence"
        localisation_confidence_column = _resolve_alias(
            site_metadata.columns,
            target_column="localisation_confidence",
            aliases=_LOCALISATION_CONFIDENCE_ALIASES,
        )
        if (
            localisation_confidence_column is not None
            and localisation_confidence_column != "localisation_confidence"
        ):
            rename_map[localisation_confidence_column] = "localisation_confidence"
        if not rename_map:
            return site_metadata
        return site_metadata.rename(columns=rename_map)

    @staticmethod
    def _normalize_site_metadata_site_identity_fields(
        site_metadata: pd.DataFrame,
    ) -> pd.DataFrame:
        if (
            "gene_symbol" not in site_metadata.columns
            or "site" not in site_metadata.columns
        ):
            return site_metadata
        gene_symbol = site_metadata.loc[:, "gene_symbol"]
        site = site_metadata.loc[:, "site"]
        gene_token = gene_symbol.astype("string").str.strip()
        site_token = site.astype("string").str.strip()
        has_blank_or_missing = (
            gene_symbol.isna()
            | site.isna()
            | gene_token.isna()
            | site_token.isna()
            | (gene_token == "")
            | (site_token == "")
        )
        if bool(has_blank_or_missing.any()):
            return site_metadata
        canonical_site_ids = canonicalize_site_components_series(
            gene_symbol=gene_symbol,
            site=site,
            field_name="dataset build request site_metadata.gene_symbol/site",
            error_type=UnsupportedInputFormatError,
            output_name="site_id",
        )
        genes: list[str] = []
        sites: list[str] = []
        for canonical_site_id in canonical_site_ids.tolist():
            gene_symbol, site = parse_canonical_site_identifier(
                canonical_site_id,
                field_name="dataset build request site_metadata.gene_symbol/site",
                error_type=UnsupportedInputFormatError,
            )
            genes.append(gene_symbol)
            sites.append(site)
        site_metadata["gene_symbol"] = genes
        site_metadata["site"] = sites
        return site_metadata

    @staticmethod
    def _derive_site_fields_from_index(site_metadata: pd.DataFrame) -> pd.DataFrame:
        normalized = site_metadata
        needs_gene_symbol = "gene_symbol" not in normalized.columns
        needs_site = "site" not in normalized.columns
        if not needs_gene_symbol and not needs_site:
            return normalized

        parsed = _parse_site_metadata_index_convention(normalized.index)
        if parsed is None:
            missing = []
            if needs_gene_symbol:
                missing.append("gene_symbol")
            if needs_site:
                missing.append("site")
            missing_columns = ", ".join(missing)
            raise UnsupportedInputFormatError(
                "dataset build request site_metadata is missing required metadata "
                f"columns ({missing_columns}). Supported conventions: provide "
                "explicit columns, or use site_metadata.index values formatted as "
                "'<gene_symbol>;<site>;' (example: 'MAPK14;Y182;'). This "
                "convention derives only gene_symbol/site; it does not derive "
                "protein_id required by the supported signalome workflow lane."
            )

        genes, sites = parsed
        if needs_gene_symbol:
            normalized.loc[:, "gene_symbol"] = genes.to_numpy(dtype=object, copy=False)
        if needs_site:
            normalized.loc[:, "site"] = sites.to_numpy(dtype=object, copy=False)
        return normalized

    @staticmethod
    def _normalize_sample_metadata(
        sample_metadata: pd.DataFrame | None,
        *,
        phospho_columns: pd.Index,
    ) -> pd.DataFrame | None:
        if sample_metadata is None:
            return None
        normalized = sample_metadata
        normalized.index = _normalize_index_labels(
            normalized.index,
            field_name="dataset build request sample_metadata.index",
            policy=_SAMPLE_LABEL_INDEX_POLICY,
        )
        if (
            not normalized.index.equals(phospho_columns)
            and normalized.index.isin(phospho_columns).all()
            and phospho_columns.isin(normalized.index).all()
        ):
            return normalized.reindex(phospho_columns)
        return normalized

    @staticmethod
    def _normalize_total(
        total: pd.DataFrame | None,
        *,
        phospho_columns: pd.Index,
    ) -> pd.DataFrame | None:
        if total is None:
            return None
        normalized = total
        normalized.index = _normalize_index_labels(
            normalized.index,
            field_name="dataset build request total.index",
            policy=_SAMPLE_LABEL_INDEX_POLICY,
        )
        normalized.columns = _normalize_index_labels(
            normalized.columns,
            field_name="dataset build request total.columns",
            policy=_SAMPLE_LABEL_INDEX_POLICY,
        )
        if (
            not normalized.columns.equals(phospho_columns)
            and normalized.columns.isin(phospho_columns).all()
            and phospho_columns.isin(normalized.columns).all()
        ):
            normalized = normalized.reindex(columns=phospho_columns)
        return normalized


@dataclass(frozen=True, slots=True)
class _IndexLabelNormalizationPolicy:
    allow_non_string_labels: bool = True
    detect_duplicate_labels_after_normalisation: bool = True


@dataclass(frozen=True, slots=True)
class _NormalizedIndexLabels:
    index: pd.Index
    raw_labels: tuple[str, ...]


_SITE_IDENTIFIER_INDEX_POLICY = _IndexLabelNormalizationPolicy(
    allow_non_string_labels=True,
    detect_duplicate_labels_after_normalisation=True,
)
_SAMPLE_LABEL_INDEX_POLICY = _IndexLabelNormalizationPolicy(
    allow_non_string_labels=True,
    detect_duplicate_labels_after_normalisation=True,
)


def _normalize_index_labels(
    index: pd.Index,
    *,
    field_name: str,
    policy: _IndexLabelNormalizationPolicy,
) -> pd.Index:
    normalized = _validate_and_normalize_index_labels(
        index,
        field_name=field_name,
        policy=policy,
    )
    return normalized.index


def _resolve_alias(
    columns: pd.Index,
    target_column: str,
    aliases: tuple[str, ...],
) -> str | None:
    present = _normalized_column_lookup(columns)
    matches: list[str] = []
    for candidate in aliases:
        matches.extend(present.get(candidate.lower(), ()))
    unique_matches = list(dict.fromkeys(matches))
    if len(unique_matches) > 1:
        present_preview = ", ".join(unique_matches)
        accepted = ", ".join(f"'{alias}'" for alias in aliases)
        raise UnsupportedInputFormatError(
            "dataset build request site_metadata has ambiguous columns for "
            f"'{target_column}': {present_preview}. Use exactly one supported "
            f"column name: {accepted}."
        )
    if unique_matches:
        return unique_matches[0]
    return None


def _reject_unsupported_historical_aliases(site_metadata: pd.DataFrame) -> None:
    columns = site_metadata.columns
    present = _normalized_column_lookup(columns)
    unsupported_aliases: dict[str, str] = {
        "sequence": "site_sequence",
        "protein": "protein_id",
        "gene": "gene_symbol",
        "residue": "site",
        "phosphosite": "site",
    }
    for alias_name, canonical_name in unsupported_aliases.items():
        if alias_name not in present:
            continue
        # `residue` can be valid auxiliary metadata when canonical `site` is already
        # present. It remains unsupported as a substitute for `site`.
        if alias_name == "residue" and "site" in present:
            continue
        raise UnsupportedInputFormatError(
            f"dataset build request site_metadata column '{alias_name}' is "
            "unsupported for strict convention mapping. Rename it to "
            f"'{canonical_name}' and provide exactly one '{canonical_name}' column."
        )


def _normalized_column_lookup(columns: pd.Index) -> dict[str, list[str]]:
    present: dict[str, list[str]] = {}
    for column in columns:
        key = str(column).strip().lower()
        present.setdefault(key, []).append(str(column))
    return present


def _parse_site_metadata_index_convention(
    index: pd.Index,
) -> tuple[pd.Series, pd.Series] | None:
    if index.empty:
        return (
            pd.Series(index=index.copy(), dtype="object"),
            pd.Series(index=index.copy(), dtype="object"),
        )
    genes: list[str] = []
    sites: list[str] = []
    for raw_site_id in index.tolist():
        try:
            gene_symbol, site = parse_canonical_site_identifier(
                raw_site_id,
                field_name="dataset build request site_metadata.index",
                error_type=UnsupportedInputFormatError,
            )
        except UnsupportedInputFormatError:
            return None
        genes.append(gene_symbol)
        sites.append(site)
    return (
        pd.Series(genes, index=index.copy(), dtype="object"),
        pd.Series(sites, index=index.copy(), dtype="object"),
    )


def _normalize_supported_site_index_if_present(
    index: pd.Index,
    *,
    field_name: str,
    site_identifier_records: list[SiteIdentifierNormalisationRecord],
) -> pd.Index:
    normalized = _validate_and_normalize_index_labels(
        index,
        field_name=field_name,
        policy=_SITE_IDENTIFIER_INDEX_POLICY,
    )
    if not _has_site_like_tokens(normalized.index):
        return normalized.index
    canonical = canonicalize_site_index(
        normalized.index,
        field_name=field_name,
        error_type=UnsupportedInputFormatError,
    )
    site_identifier_records.extend(
        _site_identifier_normalisation_changes(
            raw_labels=normalized.raw_labels,
            normalized_labels=canonical,
            field_name=field_name,
        )
    )
    return canonical


def _canonicalize_site_index_with_label_validation(
    index: pd.Index,
    *,
    field_name: str,
    site_identifier_records: list[SiteIdentifierNormalisationRecord],
    index_name: str | None = None,
) -> pd.Index:
    normalized = _validate_and_normalize_index_labels(
        index,
        field_name=field_name,
        policy=_SITE_IDENTIFIER_INDEX_POLICY,
    )
    canonical = canonicalize_site_index(
        normalized.index,
        field_name=field_name,
        error_type=UnsupportedInputFormatError,
        index_name=index_name,
    )
    site_identifier_records.extend(
        _site_identifier_normalisation_changes(
            raw_labels=normalized.raw_labels,
            normalized_labels=canonical,
            field_name=field_name,
        )
    )
    return canonical


def _validate_and_normalize_index_labels(
    index: pd.Index,
    *,
    field_name: str,
    policy: _IndexLabelNormalizationPolicy,
) -> _NormalizedIndexLabels:
    raw_objects = index.tolist()
    if not raw_objects:
        return _NormalizedIndexLabels(index=index.copy(), raw_labels=())

    raw_labels: list[str] = []
    normalized_labels: list[str] = []
    missing_positions: list[int] = []
    blank_positions: list[int] = []
    non_string_positions: list[int] = []

    for position, value in enumerate(raw_objects):
        if _is_missing_label(value):
            missing_positions.append(position)
            continue
        if not isinstance(value, str) and not policy.allow_non_string_labels:
            non_string_positions.append(position)
            continue
        raw_label = str(value)
        normalized_label = raw_label.strip()
        if normalized_label == "":
            blank_positions.append(position)
            continue
        raw_labels.append(raw_label)
        normalized_labels.append(normalized_label)

    if missing_positions:
        raise UnsupportedInputFormatError(
            f"{field_name} must not contain missing labels; found missing labels at "
            f"positions: {_position_preview(missing_positions)}"
        )
    if non_string_positions:
        raise UnsupportedInputFormatError(
            f"{field_name} must contain string labels; found non-string labels at "
            f"positions: {_position_preview(non_string_positions)}"
        )
    if blank_positions:
        raise UnsupportedInputFormatError(
            f"{field_name} must contain non-blank labels; found blank labels at "
            f"positions: {_position_preview(blank_positions)}"
        )

    normalized_index = pd.Index(normalized_labels, name=index.name)
    if policy.detect_duplicate_labels_after_normalisation:
        _raise_if_duplicate_labels_introduced_by_normalisation(
            raw_labels=raw_labels,
            normalized_labels=normalized_labels,
            field_name=field_name,
        )
    return _NormalizedIndexLabels(
        index=normalized_index,
        raw_labels=tuple(raw_labels),
    )


def _raise_if_duplicate_labels_introduced_by_normalisation(
    *,
    raw_labels: list[str],
    normalized_labels: list[str],
    field_name: str,
) -> None:
    normalized_index = pd.Index(normalized_labels)
    if normalized_index.is_unique:
        return
    duplicate_labels = list(
        dict.fromkeys(normalized_index[normalized_index.duplicated()])
    )
    introduced_by_normalisation: list[str] = []
    for label in duplicate_labels:
        raw_variants = {
            raw_label
            for raw_label, normalized_label in zip(
                raw_labels, normalized_labels, strict=False
            )
            if normalized_label == label
        }
        if raw_variants != {label}:
            introduced_by_normalisation.append(label)
    if not introduced_by_normalisation:
        return
    preview = ", ".join(repr(label) for label in introduced_by_normalisation[:5])
    suffix = "" if len(introduced_by_normalisation) <= 5 else " ..."
    raise UnsupportedInputFormatError(
        f"{field_name} contains duplicate labels introduced by normalization: "
        f"{preview}{suffix}. Provide unique labels after trimming whitespace."
    )


def _site_identifier_normalisation_changes(
    *,
    raw_labels: tuple[str, ...],
    normalized_labels: pd.Index,
    field_name: str,
) -> tuple[SiteIdentifierNormalisationRecord, ...]:
    records: list[SiteIdentifierNormalisationRecord] = []
    for position, (raw_label, normalized_label) in enumerate(
        zip(raw_labels, normalized_labels.tolist(), strict=False)
    ):
        if raw_label == normalized_label:
            continue
        records.append(
            SiteIdentifierNormalisationRecord(
                field_name=field_name,
                row_position=position,
                original_value=raw_label,
                normalised_value=str(normalized_label),
            )
        )
    return tuple(records)


def _position_preview(positions: list[int]) -> str:
    preview = ", ".join(str(position) for position in positions[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _is_missing_label(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _has_site_like_tokens(index: pd.Index) -> bool:
    for value in index.tolist():
        if ";" in str(value):
            return True
    return False
