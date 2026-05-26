"""Internal models for dataset site-sequence resolution collaborators."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.science.datasets.preprocessing.policy_models import (
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
)
from phospy.science.sequences import (
    FastaProteinSequenceRepository,
    PhosphositeSequenceResolutionRequest,
)


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionContext:
    """Stage-wide immutable context for site-sequence resolution execution."""

    repository: FastaProteinSequenceRepository
    mode: SiteSequenceResolutionMode
    conflict_policy: SiteSequenceConflictPolicy
    flank_size: int
    accession_column: str
    site_column: str
    fasta_source_path: str
    fasta_source_label: str
    fasta_sha256: str


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionRowRequest:
    """Row-level resolution request and optional pre-resolution skip details."""

    row_index: int
    row_id: object
    row_key: str
    site_id: str
    has_existing: bool
    existing_site_sequence: str | None
    resolver_request: PhosphositeSequenceResolutionRequest | None
    skip_status: str | None
    skip_reason: str | None
    skip_action: str | None


@dataclass(frozen=True, slots=True)
class SiteSequenceConflictResolution:
    """Resolved action for existing-vs-FASTA sequence conflicts."""

    status: str
    reason: str
    action: str
    resolved_site_sequence: str
    should_replace_existing: bool
    is_conflict: bool
    is_error: bool


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionDiagnostics:
    """Structured diagnostics payload consumed by stage result assembly."""

    payload: Mapping[str, object]
    has_conflict_errors: bool
