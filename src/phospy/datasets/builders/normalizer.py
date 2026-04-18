"""Internal convention normalisation for dataset builder inputs."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.site_ids import canonicalize_site_index

_GENE_SYMBOL_ALIASES = ("gene_symbol", "gene", "gene_name", "protein", "protein_id")
_SITE_ALIASES = ("site", "residue", "phosphosite", "site_position")
_SITE_SEQUENCE_ALIASES = ("site_sequence", "sequence", "centralized_sequence")


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
        normalized_phospho = phospho.copy(deep=True)
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
        normalized = site_metadata.copy(deep=True)
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
        normalized = site_metadata.copy(deep=True)
        index_name = (
            str(normalized.index.name).strip().lower() if normalized.index.name else ""
        )
        if "site_id" in normalized.columns and (
            isinstance(normalized.index, pd.RangeIndex) or index_name != "site_id"
        ):
            normalized = normalized.set_index("site_id", drop=True)
        normalized.index = canonicalize_site_index(
            normalized.index,
            field_name="dataset build request site_metadata.index",
            error_type=UnsupportedInputFormatError,
        )
        return normalized

    @staticmethod
    def _normalize_site_metadata_columns(site_metadata: pd.DataFrame) -> pd.DataFrame:
        rename_map: dict[str, str] = {}
        gene_column = _resolve_alias(site_metadata.columns, _GENE_SYMBOL_ALIASES)
        if gene_column is not None and gene_column != "gene_symbol":
            rename_map[gene_column] = "gene_symbol"
        site_column = _resolve_alias(site_metadata.columns, _SITE_ALIASES)
        if site_column is not None and site_column != "site":
            rename_map[site_column] = "site"
        site_sequence_column = _resolve_alias(
            site_metadata.columns,
            _SITE_SEQUENCE_ALIASES,
        )
        if site_sequence_column is not None and site_sequence_column != "site_sequence":
            rename_map[site_sequence_column] = "site_sequence"
        if not rename_map:
            return site_metadata
        return site_metadata.rename(columns=rename_map)

    @staticmethod
    def _derive_site_fields_from_index(site_metadata: pd.DataFrame) -> pd.DataFrame:
        normalized = site_metadata.copy(deep=True)
        split = (
            normalized.index.to_series()
            .astype(str)
            .str.strip()
            .str.split(";", expand=True)
        )
        if "gene_symbol" not in normalized.columns and split.shape[1] >= 1:
            genes = split.loc[:, 0].astype(str).str.strip()
            if (genes != "").all():
                normalized.loc[:, "gene_symbol"] = genes.to_numpy(
                    dtype=object, copy=False
                )
        if "site" not in normalized.columns and split.shape[1] >= 2:
            sites = split.loc[:, 1].astype(str).str.strip()
            if (sites != "").all():
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
        normalized = sample_metadata.copy(deep=True)
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
        normalized = total.copy(deep=True)
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


def _resolve_alias(columns: pd.Index, aliases: tuple[str, ...]) -> str | None:
    present = {str(column).strip().lower(): str(column) for column in columns}
    for candidate in aliases:
        match = present.get(candidate.lower())
        if match is not None:
            return match
    return None
