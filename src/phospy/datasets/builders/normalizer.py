"""Internal convention normalisation for dataset builder inputs."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.site_ids import canonicalize_site_index

_GENE_SYMBOL_ALIASES = ("gene_symbol", "gene_name")
_PROTEIN_ID_ALIASES = ("protein_id",)
_SITE_ALIASES = ("site",)
_SITE_SEQUENCE_ALIASES = ("site_sequence", "centralized_sequence")


@dataclass(frozen=True, slots=True)
class NormalizedDatasetInputs:
    """Normalised tables ready for sequence derivation and execution."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None


class DatasetConventionNormalizer:
    """Apply narrow, documented shaping rules for supported inputs."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
    ) -> NormalizedDatasetInputs:
        normalized_phospho = phospho
        normalized_phospho.index = canonicalize_site_index(
            normalized_phospho.index,
            field_name="dataset build request phospho.index",
            error_type=UnsupportedInputFormatError,
        )
        normalized_phospho.columns = _normalized_string_index(
            normalized_phospho.columns
        )

        normalized_site_metadata = self._normalize_site_metadata(
            site_metadata,
            phospho_index=normalized_phospho.index,
        )

        normalized_sample_metadata = self._normalize_sample_metadata(
            sample_metadata,
            phospho_columns=normalized_phospho.columns,
        )
        normalized_total = self._normalize_total(
            total,
            phospho_columns=normalized_phospho.columns,
        )
        return NormalizedDatasetInputs(
            phospho=normalized_phospho,
            site_metadata=normalized_site_metadata,
            sample_metadata=normalized_sample_metadata,
            total=normalized_total,
        )

    def _normalize_site_metadata(
        self,
        site_metadata: pd.DataFrame,
        *,
        phospho_index: pd.Index,
    ) -> pd.DataFrame:
        normalized = site_metadata
        normalized = self._normalize_site_metadata_index(normalized)
        normalized = self._normalize_site_metadata_columns(normalized)
        normalized = self._derive_site_fields_from_index(normalized)
        if (
            not normalized.index.equals(phospho_index)
            and normalized.index.isin(phospho_index).all()
            and phospho_index.isin(normalized.index).all()
        ):
            normalized = normalized.reindex(phospho_index)
        return normalized

    @staticmethod
    def _normalize_site_metadata_index(site_metadata: pd.DataFrame) -> pd.DataFrame:
        normalized = site_metadata
        index_name = (
            str(normalized.index.name).strip().lower() if normalized.index.name else ""
        )
        has_site_id_column = "site_id" in normalized.columns
        if has_site_id_column and isinstance(normalized.index, pd.RangeIndex):
            normalized = normalized.set_index("site_id", drop=True)
        elif has_site_id_column:
            index_as_site_id = canonicalize_site_index(
                normalized.index,
                field_name="dataset build request site_metadata.index",
                error_type=UnsupportedInputFormatError,
            )
            column_as_site_id = canonicalize_site_index(
                pd.Index(normalized.loc[:, "site_id"], name="site_id"),
                field_name="dataset build request site_metadata.site_id",
                error_type=UnsupportedInputFormatError,
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
        normalized.index = canonicalize_site_index(
            normalized.index,
            field_name="dataset build request site_metadata.index",
            error_type=UnsupportedInputFormatError,
        )
        return normalized

    @staticmethod
    def _normalize_site_metadata_columns(site_metadata: pd.DataFrame) -> pd.DataFrame:
        _reject_unsupported_legacy_aliases(site_metadata.columns)
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
        if not rename_map:
            return site_metadata
        return site_metadata.rename(columns=rename_map)

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
                "'<gene_symbol>;<site>;' (example: 'MAPK14;Y182;')."
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
        normalized.index = _normalized_string_index(normalized.index)
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
        normalized.index = _normalized_string_index(normalized.index)
        normalized.columns = _normalized_string_index(normalized.columns)
        if (
            not normalized.columns.equals(phospho_columns)
            and normalized.columns.isin(phospho_columns).all()
            and phospho_columns.isin(normalized.columns).all()
        ):
            normalized = normalized.reindex(columns=phospho_columns)
        return normalized


def _normalized_string_index(index: pd.Index) -> pd.Index:
    return pd.Index(index.astype(str).str.strip(), name=index.name)


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


def _reject_unsupported_legacy_aliases(columns: pd.Index) -> None:
    present = _normalized_column_lookup(columns)
    unsupported_aliases: dict[str, str] = {
        "sequence": "site_sequence",
        "protein": "protein_id",
        "gene": "gene_symbol",
        "residue": "site",
        "phosphosite": "site",
        "site_position": "site",
    }
    for legacy_name, canonical_name in unsupported_aliases.items():
        if legacy_name not in present:
            continue
        raise UnsupportedInputFormatError(
            f"dataset build request site_metadata column '{legacy_name}' is "
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
    extracted = (
        index.to_series()
        .astype(str)
        .str.extract(r"^\s*(?P<gene>[^;]+?)\s*;\s*(?P<site>[^;]+?)\s*;\s*$")
    )
    if extracted.isna().any(axis=None):
        return None
    genes = extracted.loc[:, "gene"].astype(str).str.strip()
    sites = extracted.loc[:, "site"].astype(str).str.strip()
    if (genes == "").any() or (sites == "").any():
        return None
    return genes, sites
