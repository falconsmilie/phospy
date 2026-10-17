"""Protein-sequence parsing and repository contracts."""

from phospy.references.sequences.models import (
    FastaSourceMetadata,
    ProteinSequenceLookupResult,
    ProteinSequenceRecord,
)
from phospy.references.sequences.repository import FastaProteinSequenceRepository

__all__ = [
    "FastaProteinSequenceRepository",
    "FastaSourceMetadata",
    "ProteinSequenceLookupResult",
    "ProteinSequenceRecord",
]
