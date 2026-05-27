"""Domain helpers for phosphosite identity parsing and validation."""

from phospy.science.sites.identifiers import (
    SITE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION,
    ParsedSiteToken,
    SiteIdentifierNormalisationRecord,
    SiteIdentifierNormalisationReport,
    build_site_identifier_normalisation_report,
    canonicalize_site_components,
    canonicalize_site_components_series,
    canonicalize_site_identifier,
    canonicalize_site_index,
    canonicalize_site_series,
    parse_canonical_site_identifier,
    try_parse_site_token,
)
from phospy.science.sites.identity import (
    PHOSPHOSITE_IDENTITY_OPTIONAL_COLUMNS,
    PHOSPHOSITE_PROTEIN_CONTEXT_COLUMNS,
    PhosphositeIdentity,
    build_phosphosite_identity,
    validate_identity_optional_columns,
    validate_no_conflicting_identity_collisions,
)
from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
)
from phospy.science.sites.validation import (
    require_canonical_site_index,
    require_canonical_site_series,
    require_no_mixed_site_key_isoform_scope,
    require_site_identity_coherence,
)

__all__ = [
    "ParsedSiteToken",
    "SiteIdentifierNormalisationRecord",
    "SiteIdentifierNormalisationReport",
    "SITE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION",
    "build_site_identifier_normalisation_report",
    "canonicalize_site_components",
    "canonicalize_site_components_series",
    "canonicalize_site_identifier",
    "canonicalize_site_index",
    "canonicalize_site_series",
    "parse_canonical_site_identifier",
    "try_parse_site_token",
    "PHOSPHOSITE_IDENTITY_OPTIONAL_COLUMNS",
    "PHOSPHOSITE_PROTEIN_CONTEXT_COLUMNS",
    "PhosphositeIdentity",
    "build_phosphosite_identity",
    "validate_identity_optional_columns",
    "validate_no_conflicting_identity_collisions",
    "ProteinScopedPhosphositeKey",
    "build_protein_scoped_site_key",
    "encode_site_key",
    "decode_site_key",
    "require_no_mixed_site_key_isoform_scope",
    "require_canonical_site_index",
    "require_canonical_site_series",
    "require_site_identity_coherence",
]
