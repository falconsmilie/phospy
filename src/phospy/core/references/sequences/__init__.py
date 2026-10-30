"""Protein-sequence parsing and repository contracts."""

from phospy.core.references.sequences.models import (
    FastaSourceMetadata,
    ProteinSequenceLookupResult,
    ProteinSequenceRecord,
)
from phospy.core.references.sequences.repository import FastaProteinSequenceRepository

__all__ = [
    "FastaProteinSequenceRepository",
    "FastaSourceMetadata",
    "ProteinSequenceLookupResult",
    "ProteinSequenceRecord",
]
