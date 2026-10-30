"""Provenance-domain exceptions."""

from phospy.core.errors.input import PhosPyInputError


class PhosPyProvenanceError(PhosPyInputError):
    """Provenance construction or fingerprinting failed."""


class ProvenanceFingerprintError(PhosPyProvenanceError):
    """A provenance fingerprint cannot be constructed safely."""
