"""Sequence-domain local protein repository interfaces."""

from phospy.sequences.models import (
    FastaSourceMetadata,
    ProteinSequenceLookupResult,
    ProteinSequenceRecord,
)
from phospy.sequences.repository import FastaProteinSequenceRepository

__all__ = [
    "FastaProteinSequenceRepository",
    "FastaSourceMetadata",
    "ProteinSequenceLookupResult",
    "ProteinSequenceRecord",
]
