"""Conflict-resolution policy collaborator for site-sequence stage."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.policy_models import (
    SiteSequenceConflictPolicy,
)
from phospy.science.datasets.preprocessing.site_sequence.models import (
    SiteSequenceConflictResolution,
)


class SiteSequenceConflictResolver:
    """Resolve existing-vs-FASTA sequence outcomes under conflict policy."""

    _CONFLICT_REASON = "existing site_sequence conflicts with FASTA-derived sequence"

    def resolve(
        self,
        *,
        existing_sequence: str,
        fasta_sequence: str,
        conflict_policy: SiteSequenceConflictPolicy,
    ) -> SiteSequenceConflictResolution:
        if existing_sequence == fasta_sequence:
            return SiteSequenceConflictResolution(
                status="resolved",
                reason="existing site_sequence validated against FASTA",
                action="validate_existing",
                resolved_site_sequence=existing_sequence,
                should_replace_existing=False,
                is_conflict=False,
                is_error=False,
            )

        if conflict_policy is SiteSequenceConflictPolicy.ERROR:
            return SiteSequenceConflictResolution(
                status="existing_sequence_conflict",
                reason=self._CONFLICT_REASON,
                action="error",
                resolved_site_sequence=existing_sequence,
                should_replace_existing=False,
                is_conflict=True,
                is_error=True,
            )
        if conflict_policy is SiteSequenceConflictPolicy.REPLACE_EXISTING:
            return SiteSequenceConflictResolution(
                status="existing_sequence_conflict",
                reason=self._CONFLICT_REASON,
                action="replace_existing",
                resolved_site_sequence=fasta_sequence,
                should_replace_existing=True,
                is_conflict=True,
                is_error=False,
            )
        if conflict_policy is SiteSequenceConflictPolicy.PRESERVE_EXISTING:
            return SiteSequenceConflictResolution(
                status="existing_sequence_conflict",
                reason=self._CONFLICT_REASON,
                action="preserve_existing",
                resolved_site_sequence=existing_sequence,
                should_replace_existing=False,
                is_conflict=True,
                is_error=False,
            )
        raise PhosPyInputError(
            "dataset preprocessing stage 'site_sequence_resolution' received "
            f"unsupported conflict policy: {conflict_policy!r}"
        )


__all__ = ["SiteSequenceConflictResolver"]
