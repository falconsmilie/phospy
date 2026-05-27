"""Builder-owned phosphosite display/site-key derivation collaborator."""

from __future__ import annotations

import hashlib

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.frames.ownership import own_dataframe
from phospy.science.references.models import Organism
from phospy.science.sites.identifiers import canonicalize_site_components_series
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)


class DatasetSiteIdentityDeriver:
    """Derive display and protein-scoped phosphosite identity fields."""

    def run(
        self,
        *,
        site_metadata: pd.DataFrame,
        organism: Organism | None,
        allow_gene_symbol_fallback: bool = False,
        fallback_organism: str | None = None,
        allow_non_strict_site_fallback: bool = False,
        copy_site_metadata: bool = True,
    ) -> pd.DataFrame:
        normalized = (
            own_dataframe(
                site_metadata,
                field_name="dataset build request site_metadata",
                error_type=UnsupportedInputFormatError,
            )
            if copy_site_metadata
            else site_metadata
        )
        resolved_fallback_organism = _optional_text(
            fallback_organism,
            field_name=("dataset builder site identity derivation fallback organism"),
        )
        _require_columns(
            normalized,
            field_name="dataset build request site_metadata",
            required_columns=("gene_symbol", "site"),
        )
        display_id = canonicalize_site_components_series(
            gene_symbol=normalized.loc[:, "gene_symbol"],
            site=normalized.loc[:, "site"],
            field_name="dataset build request site_metadata.gene_symbol/site",
            error_type=UnsupportedInputFormatError,
            output_name="display_id",
        )
        normalized.loc[:, "display_id"] = display_id.to_numpy(dtype=object, copy=False)

        site_keys: list[str] = []
        for row_position, row_id in enumerate(normalized.index.tolist()):
            row_field = f"dataset build request site_metadata[{row_id!r}]"
            resolved_organism = _resolve_organism(
                site_metadata=normalized,
                row_position=row_position,
                request_organism=organism,
                field_name=row_field,
                fallback_organism=resolved_fallback_organism,
            )
            protein_source, protein_identifier = _resolve_protein_identifier(
                site_metadata=normalized,
                row_position=row_position,
                field_name=row_field,
                allow_gene_symbol_fallback=allow_gene_symbol_fallback,
            )
            protein_namespace = _resolve_protein_namespace(
                site_metadata=normalized,
                row_position=row_position,
                fallback_namespace=protein_source,
                field_name=row_field,
            )
            isoform_id = _optional_text(
                _series_value_or_none(
                    normalized,
                    row_position=row_position,
                    column="isoform_id",
                ),
                field_name=f"{row_field}.isoform_id",
            )
            parsed_site = _parse_site_for_key(
                value=_series_value_or_none(
                    normalized,
                    row_position=row_position,
                    column="site",
                ),
                field_name=f"{row_field}.site",
                allow_non_strict_site_fallback=allow_non_strict_site_fallback,
            )
            site_key_isoform_id = _merge_site_key_isoform_id(
                isoform_id=isoform_id,
                opaque_site_token=parsed_site.opaque_site_token,
            )
            key = build_protein_scoped_site_key(
                organism=resolved_organism,
                protein_namespace=protein_namespace,
                protein_identifier=protein_identifier,
                residue=parsed_site.residue,
                position=parsed_site.position,
                isoform_id=site_key_isoform_id,
                field_name=f"{row_field}.site_key",
                error_type=UnsupportedInputFormatError,
            )
            site_keys.append(encode_site_key(key))
        normalized.loc[:, "site_key"] = site_keys
        return normalized


def _resolve_organism(
    *,
    site_metadata: pd.DataFrame,
    row_position: int,
    request_organism: Organism | None,
    field_name: str,
    fallback_organism: str | None,
) -> str:
    column_value = _series_value_or_none(
        site_metadata,
        row_position=row_position,
        column="organism",
    )
    explicit = _optional_text(
        column_value,
        field_name=f"{field_name}.organism",
        allow_enum=True,
    )
    if explicit is not None:
        return explicit
    if request_organism is not None:
        return str(request_organism.value).strip()
    if fallback_organism is not None:
        return fallback_organism
    raise UnsupportedInputFormatError(
        f"{field_name} requires organism to derive site_key. Provide non-empty "
        "site_metadata.organism or dataset build request organism."
    )


def _resolve_protein_identifier(
    *,
    site_metadata: pd.DataFrame,
    row_position: int,
    field_name: str,
    allow_gene_symbol_fallback: bool,
) -> tuple[str, str]:
    protein_accession = _optional_text(
        _series_value_or_none(
            site_metadata,
            row_position=row_position,
            column="protein_accession",
        ),
        field_name=f"{field_name}.protein_accession",
    )
    if protein_accession is not None:
        return "protein_accession", protein_accession
    protein_id = _optional_text(
        _series_value_or_none(
            site_metadata,
            row_position=row_position,
            column="protein_id",
        ),
        field_name=f"{field_name}.protein_id",
    )
    if protein_id is not None:
        return "protein_id", protein_id
    if allow_gene_symbol_fallback:
        gene_symbol = _optional_text(
            _series_value_or_none(
                site_metadata,
                row_position=row_position,
                column="gene_symbol",
            ),
            field_name=f"{field_name}.gene_symbol",
        )
        if gene_symbol is not None:
            return "gene_symbol", gene_symbol
    raise UnsupportedInputFormatError(
        f"{field_name} requires protein_accession or protein_id to derive site_key"
    )


def _resolve_protein_namespace(
    *,
    site_metadata: pd.DataFrame,
    row_position: int,
    fallback_namespace: str,
    field_name: str,
) -> str:
    explicit = _optional_text(
        _series_value_or_none(
            site_metadata,
            row_position=row_position,
            column="protein_namespace",
        ),
        field_name=f"{field_name}.protein_namespace",
    )
    if explicit is not None:
        return explicit
    return fallback_namespace


class _ResolvedSiteKeySiteComponents:
    __slots__ = ("residue", "position", "opaque_site_token")

    def __init__(
        self,
        *,
        residue: str,
        position: int,
        opaque_site_token: str | None,
    ) -> None:
        self.residue = residue
        self.position = position
        self.opaque_site_token = opaque_site_token


def _parse_site_for_key(
    *,
    field_name: str,
    value: object,
    allow_non_strict_site_fallback: bool,
) -> _ResolvedSiteKeySiteComponents:
    from phospy.science.sites.identifiers import try_parse_site_token

    parsed = try_parse_site_token(value)
    if parsed is not None:
        return _ResolvedSiteKeySiteComponents(
            residue=parsed.residue,
            position=parsed.position,
            opaque_site_token=None,
        )
    if not allow_non_strict_site_fallback:
        raise UnsupportedInputFormatError(
            f"{field_name} must use strict 'S/T/Y<position>' site tokens to derive "
            f"site_key; got {value!r}"
        )
    opaque_site_token = _optional_text(value, field_name=field_name)
    if opaque_site_token is None:
        raise UnsupportedInputFormatError(
            f"{field_name} must contain a non-empty site token to derive site_key"
        )
    digest = hashlib.sha256(opaque_site_token.upper().encode("utf-8")).hexdigest()
    synthetic_position = (int(digest[:12], 16) % 2_147_483_647) + 1
    synthetic_residue = _resolve_fallback_residue(opaque_site_token)
    return _ResolvedSiteKeySiteComponents(
        residue=synthetic_residue,
        position=synthetic_position,
        opaque_site_token=opaque_site_token.upper(),
    )


def _resolve_fallback_residue(site_token: str) -> str:
    upper = site_token.upper()
    for token in upper:
        if token in {"S", "T", "Y"}:
            return token
    return "S"


def _merge_site_key_isoform_id(
    *,
    isoform_id: str | None,
    opaque_site_token: str | None,
) -> str | None:
    if opaque_site_token is None:
        return isoform_id
    marker = f"opaque_site_token:{opaque_site_token}"
    if isoform_id is None:
        return marker
    return f"{isoform_id}|{marker}"


def _require_columns(
    site_metadata: pd.DataFrame,
    *,
    field_name: str,
    required_columns: tuple[str, ...],
) -> None:
    missing = [name for name in required_columns if name not in site_metadata.columns]
    if not missing:
        return
    missing_csv = ", ".join(missing)
    raise UnsupportedInputFormatError(
        f"{field_name} is missing required metadata columns: {missing_csv}"
    )


def _series_value_or_none(
    frame: pd.DataFrame,
    *,
    row_position: int,
    column: str,
) -> object:
    if column not in frame.columns:
        return None
    column_position = frame.columns.tolist().index(column)
    return frame.iat[row_position, column_position]


def _optional_text(
    value: object,
    *,
    field_name: str,
    allow_enum: bool = False,
) -> str | None:
    if _is_missing(value):
        return None
    if allow_enum and isinstance(value, Organism):
        return str(value.value).strip()
    if not isinstance(value, str):
        raise UnsupportedInputFormatError(
            f"{field_name} must be a string when provided"
        )
    token = value.strip()
    if token == "":
        return None
    return token


def _is_missing(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])
