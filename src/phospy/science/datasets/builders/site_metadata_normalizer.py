"""Site-metadata focused collaborators for dataset convention normalisation."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.science.datasets.builders.normalization_reporter import (
    DatasetConventionNormalisationReporter,
)
from phospy.science.sites.identifiers import (
    SiteIdentifierNormalisationRecord,
    canonicalize_site_components_series,
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


class SiteMetadataIndexNormalizer:
    """Normalize site-metadata row identity/index conventions."""

    def __init__(
        self,
        *,
        reporter: DatasetConventionNormalisationReporter | None = None,
    ) -> None:
        self._reporter = reporter or DatasetConventionNormalisationReporter()

    def run(
        self,
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
            normalized_site_id = (
                self._reporter.canonicalize_site_index_with_label_validation(
                    pd.Index(normalized.loc[:, "site_id"], name="site_id"),
                    field_name="dataset build request site_metadata.site_id",
                    site_identifier_records=site_identifier_records,
                    index_name="site_id",
                )
            )
            normalized["site_id"] = normalized_site_id.tolist()
            normalized = normalized.set_index("site_id", drop=True)
        elif has_site_id_column:
            column_as_site_id = (
                self._reporter.canonicalize_site_index_with_label_validation(
                    pd.Index(normalized.loc[:, "site_id"], name="site_id"),
                    field_name="dataset build request site_metadata.site_id",
                    site_identifier_records=site_identifier_records,
                )
            )
            if not _has_site_like_tokens(normalized.index):
                normalized = normalized.copy()
                normalized["site_id"] = column_as_site_id.tolist()
                normalized = normalized.set_index("site_id", drop=True)
                normalized.index = pd.Index(normalized.index.tolist(), name="site_id")
                return normalized
            index_as_site_id = (
                self._reporter.canonicalize_site_index_with_label_validation(
                    normalized.index,
                    field_name="dataset build request site_metadata.index",
                    site_identifier_records=site_identifier_records,
                )
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
        normalized.index = self._reporter.normalize_supported_site_index_if_present(
            normalized.index,
            field_name="dataset build request site_metadata.index",
            site_identifier_records=site_identifier_records,
        )
        return normalized


class SiteMetadataColumnAliasResolver:
    """Resolve strict site-metadata alias mapping onto standard column names."""

    def run(self, site_metadata: pd.DataFrame) -> pd.DataFrame:
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


class SiteIdentityFieldNormalizer:
    """Normalize explicit gene/site fields when both are present and usable."""

    def run(
        self,
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
            parsed_gene_symbol, parsed_site = parse_canonical_site_identifier(
                canonical_site_id,
                field_name="dataset build request site_metadata.gene_symbol/site",
                error_type=UnsupportedInputFormatError,
            )
            genes.append(parsed_gene_symbol)
            sites.append(parsed_site)
        site_metadata["gene_symbol"] = genes
        site_metadata["site"] = sites
        return site_metadata


class SiteMetadataNormalizer:
    """Coordinate site-metadata index/alias/identity normalisation."""

    def __init__(
        self,
        *,
        index_normalizer: SiteMetadataIndexNormalizer | None = None,
        column_alias_resolver: SiteMetadataColumnAliasResolver | None = None,
        identity_field_normalizer: SiteIdentityFieldNormalizer | None = None,
    ) -> None:
        self._index_normalizer = index_normalizer or SiteMetadataIndexNormalizer()
        self._column_alias_resolver = (
            column_alias_resolver or SiteMetadataColumnAliasResolver()
        )
        self._identity_field_normalizer = (
            identity_field_normalizer or SiteIdentityFieldNormalizer()
        )

    def run(
        self,
        site_metadata: pd.DataFrame,
        *,
        phospho_index: pd.Index,
        site_identifier_records: list[SiteIdentifierNormalisationRecord],
    ) -> pd.DataFrame:
        normalized = self._index_normalizer.run(
            site_metadata,
            site_identifier_records=site_identifier_records,
        )
        normalized = self._column_alias_resolver.run(normalized)
        normalized = self._derive_site_fields_from_index(normalized)
        normalized = self._identity_field_normalizer.run(normalized)
        if (
            not normalized.index.equals(phospho_index)
            and normalized.index.isin(phospho_index).all()
            and phospho_index.isin(normalized.index).all()
        ):
            normalized = normalized.reindex(phospho_index)
        return normalized

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
            normalized["gene_symbol"] = genes.astype("string")
        if needs_site:
            normalized["site"] = sites.astype("string")
        return normalized


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


def _has_site_like_tokens(index: pd.Index) -> bool:
    for value in index.tolist():
        if ";" in str(value):
            return True
    return False
