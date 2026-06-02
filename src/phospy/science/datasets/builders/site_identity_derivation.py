"""Builder-owned phosphosite display/site-key derivation collaborator."""

from __future__ import annotations

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
        normalized["display_id"] = display_id.astype(str).tolist()

        site_keys: list[str] = []
        for row_position, row_id in enumerate(normalized.index.tolist()):
            row_field = f"dataset build request site_metadata[{row_id!r}]"
            resolved_organism = _resolve_organism(
                site_metadata=normalized,
                row_position=row_position,
                request_organism=organism,
                field_name=row_field,
            )
            protein_namespace, protein_identifier = _resolve_protein_identity(
                site_metadata=normalized,
                row_position=row_position,
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
            )
            key = build_protein_scoped_site_key(
                organism=resolved_organism,
                protein_namespace=protein_namespace,
                protein_identifier=protein_identifier,
                residue=parsed_site.residue,
                position=parsed_site.position,
                isoform_id=isoform_id,
                field_name=f"{row_field}.site_key",
                error_type=UnsupportedInputFormatError,
            )
            site_keys.append(encode_site_key(key))
        normalized["site_key"] = site_keys
        return normalized


def _resolve_organism(
    *,
    site_metadata: pd.DataFrame,
    row_position: int,
    request_organism: Organism | None,
    field_name: str,
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
    raise UnsupportedInputFormatError(
        f"{field_name} organism is required to derive site_key. Provide non-empty "
        "site_metadata.organism or dataset build request organism."
    )


def _resolve_protein_identity(
    *,
    site_metadata: pd.DataFrame,
    row_position: int,
    field_name: str,
) -> tuple[str, str]:
    explicit_identifier = _optional_text(
        _series_value_or_none(
            site_metadata,
            row_position=row_position,
            column="protein_identifier",
        ),
        field_name=f"{field_name}.protein_identifier",
    )
    if explicit_identifier is not None:
        explicit_namespace = _optional_text(
            _series_value_or_none(
                site_metadata,
                row_position=row_position,
                column="protein_namespace",
            ),
            field_name=f"{field_name}.protein_namespace",
        )
        if explicit_namespace is None:
            raise UnsupportedInputFormatError(
                f"{field_name} protein context is required to derive site_key; "
                "site_metadata.protein_identifier requires non-empty "
                "site_metadata.protein_namespace."
            )
        return explicit_namespace, explicit_identifier

    protein_accession = _optional_text(
        _series_value_or_none(
            site_metadata,
            row_position=row_position,
            column="protein_accession",
        ),
        field_name=f"{field_name}.protein_accession",
    )
    if protein_accession is not None:
        return (
            _resolve_protein_namespace(
                site_metadata=site_metadata,
                row_position=row_position,
                fallback_namespace="protein_accession",
                field_name=field_name,
            ),
            protein_accession,
        )
    protein_id = _optional_text(
        _series_value_or_none(
            site_metadata,
            row_position=row_position,
            column="protein_id",
        ),
        field_name=f"{field_name}.protein_id",
    )
    if protein_id is not None:
        return (
            _resolve_protein_namespace(
                site_metadata=site_metadata,
                row_position=row_position,
                fallback_namespace="protein_id",
                field_name=field_name,
            ),
            protein_id,
        )
    gene_symbol = _optional_text(
        _series_value_or_none(
            site_metadata,
            row_position=row_position,
            column="gene_symbol",
        ),
        field_name=f"{field_name}.gene_symbol",
    )
    gene_context_note = (
        " gene_symbol is display metadata and is not protein context."
        if gene_symbol is not None
        else ""
    )
    raise UnsupportedInputFormatError(
        f"{field_name} protein context is required to derive site_key; provide "
        "site_metadata.protein_accession, site_metadata.protein_id, or "
        "site_metadata.protein_identifier with site_metadata.protein_namespace."
        f"{gene_context_note}"
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
    __slots__ = ("residue", "position")

    def __init__(
        self,
        *,
        residue: str,
        position: int,
    ) -> None:
        self.residue = residue
        self.position = position


def _parse_site_for_key(
    *,
    field_name: str,
    value: object,
) -> _ResolvedSiteKeySiteComponents:
    from phospy.science.sites.identifiers import try_parse_site_token

    parsed = try_parse_site_token(value)
    if parsed is not None:
        return _ResolvedSiteKeySiteComponents(
            residue=parsed.residue,
            position=parsed.position,
        )
    raise UnsupportedInputFormatError(
        f"{field_name} requires strict residue/position site token ('S/T/Y' "
        "followed by a positive integer) to derive site_key; "
        f"got {value!r}"
    )


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
