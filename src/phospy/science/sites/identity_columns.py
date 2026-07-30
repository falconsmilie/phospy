"""Compatibility route for science-owned phosphosite identity column guards."""

from phospy.science.sites.identity_contracts import (
    enforce_display_id_column,
    enforce_site_key_column,
)

__all__ = ["enforce_display_id_column", "enforce_site_key_column"]
