"""Reference-loading collaborator for dataset site-sequence resolution."""

from __future__ import annotations

from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.datasets.preprocessing.policy_models import (
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
)
from phospy.science.datasets.preprocessing.site_sequence.models import (
    SiteSequenceResolutionContext,
)
from phospy.science.sequences import FastaProteinSequenceRepository


class SiteSequenceReferenceLoader:
    """Load FASTA repository and resolve execution context for a stage run."""

    def load(self, plan: PreprocessingPlan) -> SiteSequenceResolutionContext:
        if plan.site_sequence_resolution_fasta_path is None:
            raise ValueError(
                "site_sequence_resolution_fasta_path is required when stage is enabled"
            )
        repository = FastaProteinSequenceRepository.from_path(
            plan.site_sequence_resolution_fasta_path,
            source_label="dataset.site_sequence_resolution",
        )
        mode = plan.site_sequence_resolution_mode
        conflict_policy = plan.site_sequence_resolution_conflict_policy
        if (
            mode is SiteSequenceResolutionMode.REPLACE_EXISTING
            and conflict_policy is SiteSequenceConflictPolicy.PRESERVE_EXISTING
        ):
            # Backward-compatible behavior for legacy plans that used mode-only
            # replacement semantics before explicit conflict-policy support.
            conflict_policy = SiteSequenceConflictPolicy.REPLACE_EXISTING
        return SiteSequenceResolutionContext(
            repository=repository,
            mode=mode,
            conflict_policy=conflict_policy,
            flank_size=int(plan.site_sequence_resolution_flank_size),
            accession_column=plan.site_sequence_resolution_accession_column,
            site_column=plan.site_sequence_resolution_site_column,
            fasta_source_path=repository.metadata.source_path,
            fasta_source_label=repository.metadata.source_label,
            fasta_sha256=repository.metadata.sha256,
        )


__all__ = ["SiteSequenceReferenceLoader"]
