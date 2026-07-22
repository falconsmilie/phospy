"""Compatibility re-exports for reference table schemas."""

from phospy.science.tables.references import (
    KinaseSubstrateReference,
    ProteinAccessionReference,
    SiteSequenceReference,
)
from phospy.science.tables.references import (
    _build_index_classification_records as _build_index_classification_records,
)
from phospy.science.tables.references import (
    _build_pair_classification_records as _build_pair_classification_records,
)
from phospy.science.tables.references import (
    _build_protein_accession_classification_records as _build_protein_accession_classification_records,
)
from phospy.science.tables.references import (
    _classify_duplicate_and_conflicting_index_records as _classify_duplicate_and_conflicting_index_records,
)
from phospy.science.tables.references import (
    _classify_duplicate_and_conflicting_pair_records as _classify_duplicate_and_conflicting_pair_records,
)
from phospy.science.tables.references import (
    _classify_duplicate_and_conflicting_protein_accession_records as _classify_duplicate_and_conflicting_protein_accession_records,
)
from phospy.science.tables.references import (
    _group_has_conflicting_payload_rows as _group_has_conflicting_payload_rows,
)
from phospy.science.tables.references import (
    _raise_with_identifier_normalisation_report as _raise_with_identifier_normalisation_report,
)

__all__ = [
    "KinaseSubstrateReference",
    "ProteinAccessionReference",
    "SiteSequenceReference",
]
