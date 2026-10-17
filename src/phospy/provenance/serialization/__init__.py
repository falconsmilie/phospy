"""Generic provenance table serialization."""

from __future__ import annotations

from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload as table_fingerprint_from_payload,
)
from phospy.provenance.serialization.tables import (
    table_fingerprint_to_payload as table_fingerprint_to_payload,
)

__all__ = [
    "table_fingerprint_from_payload",
    "table_fingerprint_to_payload",
]
