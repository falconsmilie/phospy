"""Request and internal state models for FragPipe imports."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from phospy.io.readers.fragpipe.constants import _DEFAULT_INTENSITY_PREFIXES
from phospy.validation.datasets.fragpipe import (
    FRAGPIPE_FLAG_POLICY_REMOVE,
    FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE,
)


@dataclass(frozen=True, slots=True)
class FragPipeColumnMapping:
    """Optional source-column overrides for FragPipe/PTMProphet imports."""

    protein_accession: str | None = None
    gene_symbol: str | None = None
    peptide_sequence: str | None = None
    modified_peptide_sequence: str | None = None
    ptmprophet_probabilities: str | None = None
    protein_start: str | None = None
    site: str | None = None
    site_sequence: str | None = None
    intensity_columns: Mapping[str, str] | Sequence[str] | None = None
    contaminant: str | None = None
    decoy: str | None = None
    row_id: str | None = None
    unique_feature_id: str | None = None


@dataclass(frozen=True, slots=True)
class FragPipePTMProphetImportRequest:
    """Request for importing FragPipe/Philosopher/PTMProphet phosphosite output."""

    source: object
    column_mapping: FragPipeColumnMapping = field(default_factory=FragPipeColumnMapping)
    contaminant_policy: str = FRAGPIPE_FLAG_POLICY_REMOVE
    decoy_policy: str = FRAGPIPE_FLAG_POLICY_REMOVE
    ptmprophet_position_reference: str = FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE
    intensity_column_prefixes: Sequence[str] = _DEFAULT_INTENSITY_PREFIXES
    source_name: str = "fragpipe_ptmprophet"


@dataclass(frozen=True, slots=True)
class _ResolvedFragPipeColumns:
    protein_accession: str
    gene_symbol: str
    peptide_sequence: str
    modified_peptide_sequence: str
    ptmprophet_probabilities: str
    protein_start: str | None
    site: str | None
    site_sequence: str | None
    intensity_columns: dict[str, str]
    contaminant: str | None
    decoy: str | None
    row_id: str | None
    unique_feature_id: str | None


@dataclass(frozen=True, slots=True)
class _LocalisationCandidate:
    residue: str
    position: int | None
    probability: float | None

    @property
    def has_position(self) -> bool:
        return self.position is not None


@dataclass(frozen=True, slots=True)
class _ProteinSiteCandidate:
    residue: str
    protein_position: int
    probability: float | None

    @property
    def token(self) -> str:
        return f"{self.residue}{self.protein_position}"


@dataclass(frozen=True, slots=True)
class _SiteCall:
    site_tokens: tuple[str, ...]
    peptide_site_string: str
    localisation_confidence: object
    candidate_sites: str
    site_probabilities: str
    ambiguous: bool
    phospho_site_count: int


__all__ = [
    "FragPipeColumnMapping",
    "FragPipePTMProphetImportRequest",
]
