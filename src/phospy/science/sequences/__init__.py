"""Sequence-domain local protein repository interfaces."""

from phospy.science.sequences.models import (
    FastaSourceMetadata,
    ProteinSequenceLookupResult,
    ProteinSequenceRecord,
)
from phospy.science.sequences.repository import FastaProteinSequenceRepository
from phospy.science.sequences.resolver import (
    PhosphositeSequenceResolutionRequest,
    PhosphositeSequenceResolutionResult,
    PhosphositeSequenceResolver,
)

__all__ = [
    "FastaProteinSequenceRepository",
    "FastaSourceMetadata",
    "PhosphositeSequenceResolutionRequest",
    "PhosphositeSequenceResolutionResult",
    "PhosphositeSequenceResolver",
    "ProteinSequenceLookupResult",
    "ProteinSequenceRecord",
]
